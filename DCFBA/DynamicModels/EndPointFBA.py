import cbmpy
from numpy import Inf
from cbmpy.CBModel import Reaction, Species
from ..Models.CommunityModel import CommunityModel
from ..Helpers.BuildEndPointModel import build_time_model
from .DynamicModelBase import DynamicModelBase
from ..Models import KineticsStruct


class EndPointFBA(DynamicModelBase):
    """EndPointFBA class

    This class provides the blueprint and functionality to perform EndPointFBA
    on a CommunityModel.

    Attributes:
        m_model (CommunityModel): The community model on which EndPointFBA is
            performed
        m_times (list[str]): list of time point ids
        m_kinetics (KineticsStruct): A KineticsStruct object holding
            information on reaction kinetics
    """

    m_model: CommunityModel
    m_times: list[str]
    m_kinetics: KineticsStruct
    m_dt: float

    def __init__(
        self,
        community_model: CommunityModel,
        n: int,
        initial_biomasses: dict[str, float],
        initial_concentrations: dict[str, float] = {},
        dt: float = 0.1,
        kinetics: KineticsStruct = None,
    ) -> None:
        """Initializes the EndPointFBA class

        Args:
            community_model (CommunityModel): The community model to perform
                EndPointFBA on.
            n (int): Number of time steps for the model.
            initial_biomasses (dict[str, float]): Dictionary mapping model ID
                to initial biomass concentrations.
            initial_concentrations (dict[str, float], optional): Initial
                concentrations of the metabolites. Defaults to an empty
                dictionary. If no metabolite concentrations are specified,
                the lower bound of each exchange reaction is used as the
                starting value.
            kinetics (KineticsStruct, optional): Kinetic information for the
                model. Defaults to None.
            dt (float, optional): Size of the time step. Defaults to 0.1.
        """
        width = len(str(n))
        self.m_times = [f"time{i:0{width}d}" for i in range(n)]
        self.m_dt = dt

        self.m_model = build_time_model(community_model, self.m_times)
        self.m_kinetics = kinetics
        self.set_objective()
        self.set_constraints(n, initial_biomasses, dt)

        self.set_initial_concentrations(
            initial_biomasses, initial_concentrations
        )

    def set_objective(self):
        """Creates the community biomass reaction and sets it
        to be the objective of the model"""

        self.m_model.createReaction("X_comm", silent=True)
        out: Reaction = self.m_model.getReaction("X_comm")
        out.is_exchange = True
        out.setUpperBound(Inf)
        out.setLowerBound(0)
        out.createReagent("BM_c_" + self.m_times[-1], -1)

        self.m_model.createObjectiveFunction("X_comm")

        self.m_model.setActiveObjective("X_comm_objective")

    def simulate(
        self,
    ):
        return cbmpy.doFBA(self.m_model, quiet=False)

    def set_constraints(
        self, n: int, initial_biomasses: dict[str, float], dt: float
    ):
        """
        Configures the constraints for the EndPointFBA model. Unlike using
        upper and lower bounds, it adjusts constraints for each time point's
        reaction based on biomass, dt, and the initial bound.

        Args:
            n (int): Number of time points.
            initial_biomasses (dict[str, float]): Dictionary mapping model ID
                to initial biomass concentrations.
            dt (float): Time step size.
        """

        rids_t0 = self.m_model.getReactionIds(self.m_times[0])
        for rid in rids_t0:
            reaction = self.m_model.getReaction(rid)
            mid = self.m_model.identify_model_from_reaction(rid)
            # Skip exchange reactions and time to time reactions
            if reaction.is_exchange or mid == "":
                continue
            biomass = initial_biomasses[mid]

            reaction.setLowerBound(reaction.getLowerBound() * biomass * dt)
            reaction.setUpperBound(reaction.getUpperBound() * biomass * dt)

        # From time 1 to last time point
        for i in range(1, n):
            rids = self.m_model.getReactionIds(self.m_times[i])
            for rid in rids:
                reaction: Reaction = self.m_model.getReaction(rid)
                mid = self.m_model.identify_model_from_reaction(rid)

                # If the reaction is exchange, a linking layer reaction skip
                # (BM_) is for the biomass linking reactions
                if reaction.is_exchange or mid == "" or rid.startswith("BM_"):
                    continue

                # Bm of model id
                r_x_t = f"BM_{mid}_{self.m_times[i-1]}_{self.m_times[i]}"

                # r_x_t = biomass_rid + times[i - 1]
                ub = reaction.getUpperBound()
                self.m_model.addUserConstraint(
                    f"{rid}_ub",
                    [
                        [1, rid],
                        [-1 * dt * ub, r_x_t],
                    ],
                    "<=",
                    0.0,
                )
                reaction.setUpperBound(cbmpy.INF)

    def set_initial_concentrations(
        self,
        initial_biomasses: dict[str, float],
        initial_concentrations: dict[str, float],
    ):
        """
        Sets the exchange reactions to the initial concentrations of the
        metabolites and biomasses.

        Args:
            initial_biomasses (dict[str, float]): Dictionary mapping model ID
                to initial biomass concentrations.
            initial_concentrations (dict[str, float]): Dictionary mapping
                metabolite IDs to their initial concentrations.
        """
        for key, value in initial_concentrations.items():
            # get species and it's corresponding exchange reaction
            species: Species = self.m_model.getSpecies(
                key + "_" + self.m_times[0]
            )
            rids = species.getReagentOf()
            for rid in rids:
                reaction: Reaction = self.m_model.getReaction(rid)
                if reaction.is_exchange:
                    reaction.setLowerBound(-value)

        for key, value in initial_biomasses.items():
            self.m_model.setReactionBounds(
                f"BM_{key}_exchange", -value, -value
            )

    def constrain_rates(self, epsilon=0.1):
        """
        Constrains the difference in reaction rates between t_n-1 and t_n.

        See the docs for further explanation

        Args:
            epsilon (float, optional): Maximum allowed rate difference between
                successive time points. Defaults to 0.1.
        """
        old_rids = set(
            [id.split("_time")[0] for id in self.m_model.getReactionIds()]
        )
        for rid in old_rids:
            reaction: Reaction = self.m_model.getReaction(
                rid + "_" + self.m_times[0]
            )

            if reaction is not None and (not reaction.is_exchange):
                for i, time in enumerate(self.m_times[:-1]):
                    id_t0 = f"{rid}_{time}"
                    id_t1 = f"{rid}_{self.m_times[i+1]}"
                    self.m_model.addUserConstraint(
                        "R_constraint_pos" + id_t0,
                        [[1, id_t0], [-1, id_t1]],
                        "<=",
                        epsilon,
                    )

                    self.m_model.addUserConstraint(
                        "R_constraint_neg" + id_t0,
                        [[1, id_t0], [-1, id_t1]],
                        ">=",
                        -epsilon,
                    )

    def set_kinetics(self, kin: KineticsStruct):
        """
        Set the kinetics information for the model.

        Args:
            kin (KineticsStruct): Object containing detailed kinetics
                information for reactions in the model.
        """
        self.m_kinetics = kin

    def mm_approximation(self):
        """Using the KineticsStruct you can set constraints on the reaction
        rates with an approximation of the Michaelis-Menten curve.

        See documentation for further clarification

        Raises:
            Exception: Raised when no limiting substrate is defined in the
                kinetics object for a specific reaction.
        """
        if self.m_kinetics is not None:
            for rid, _ in self.m_kinetics.get_kinetics().items():
                sid, km, vmax = self.m_kinetics.get_reactions_kinetics(rid)
                if sid == "":
                    raise Exception(
                        "No limiting substrate was set in the kinetics object for: "
                        + rid
                    )
                low_line = vmax / km
                high_line = (vmax / 2) / km

                for i in range(0, len(self.m_times) - 1):
                    # Linking reaction is the concentration of S for the
                    # timepoint
                    linking_reaction_id = (
                        f"{sid}_{self.m_times[i]}_{self.m_times[i+1]}"
                    )
                    t_rid = rid + "_" + self.m_times[i + 1]
                    self.m_model.addUserConstraint(
                        f"mm_low_{linking_reaction_id}",
                        [
                            [1, t_rid],
                            [-low_line * self.m_dt, linking_reaction_id],
                        ],
                        ">=",
                        0.0,
                    )

                    self.m_model.addUserConstraint(
                        f"mm_high_{linking_reaction_id}",
                        [
                            [1, t_rid],
                            [-high_line * self.m_dt, linking_reaction_id],
                        ],
                        "<=",
                        0.0,
                    )
