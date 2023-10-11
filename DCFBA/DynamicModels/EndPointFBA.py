import cbmpy
import numpy
import re
from cbmpy.CBModel import Reaction, Species

from ..Exceptions import SpeciesNotFound
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
        # start_time = time.time()
        self.m_model = build_time_model(community_model, self.m_times)
        # print(
        #     f"--- build_time_model: {n} ------ {(time.time() - start_time)} seconds ---"
        # )

        self.m_kinetics = kinetics

        # start_time = time.time()
        self.set_objective()
        # print(
        #     f"--- set Objective: {n} ------ {(time.time() - start_time)} seconds ---"
        # )

        # start_time = time.time()
        self.set_constraints(community_model, n, dt)
        # print(
        #     f"--- set_constraints: {n} ------ {(time.time() - start_time)} seconds ---"
        # )
        # start_time = time.time()
        self.set_initial_concentrations(
            initial_biomasses, initial_concentrations
        )
        # print(
        #     f"--- set_initial_concentrations: {n} ------ {(time.time() - start_time)} seconds ---"
        # )

    def set_objective(self):
        """Creates the community biomass reaction and sets it
        to be the objective of the model"""

        self.m_model.createReaction("X_comm", silent=True)
        out: Reaction = self.m_model.getReaction("X_comm")
        out.is_exchange = True
        out.setUpperBound(cbmpy.INF)
        out.setLowerBound(0)
        out.createReagent("BM_c_" + self.m_times[-1], -1)

        self.m_model.createObjectiveFunction("X_comm")

        self.m_model.setActiveObjective("X_comm_objective")

    def simulate(
        self,
    ):
        return cbmpy.doFBA(self.m_model, quiet=False)

    def set_constraints(
        self,
        initial_model: CommunityModel,
        n: int,
        dt: float,
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
        self.m_model.__FBC_VERSION__ = 3

        # Lookup time of set is on average O(1)
        rids_lb_to_check = set()
        rids_ub_to_check = set()
        for reaction in initial_model.reactions:
            if reaction.is_exchange:
                continue
            lb = reaction.getLowerBound()
            ub = reaction.getUpperBound()
            if not (lb == 0.0 or lb == cbmpy.INF or lb == cbmpy.NINF):
                rids_lb_to_check.add(reaction.getId())
            if not (ub == 0.0 or ub == cbmpy.INF or ub == cbmpy.NINF):
                rids_ub_to_check.add(reaction.getId())

        combined_set = rids_lb_to_check | rids_ub_to_check
        # Reactions at time zero
        for rid in combined_set:
            new_rid = rid + "_" + self.m_times[0]
            reaction = self.m_model.getReaction(new_rid)
            mid = self.m_model.identify_model_from_reaction(rid)

            if rid in rids_lb_to_check:
                udc = self.m_model.createUserDefinedConstraint(
                    f"{new_rid}_lb",
                    0.0,
                    numpy.Inf,
                    components=[
                        (1, new_rid, "linear"),
                        (
                            1 * dt * reaction.getLowerBound(),
                            f"BM_{mid}_exchange",
                            "linear",
                        ),
                    ],
                )
                self.m_model.addUserDefinedConstraint(udc)
                reaction.setLowerBound(cbmpy.NINF)
            if rid in rids_ub_to_check:
                udc = self.m_model.createUserDefinedConstraint(
                    f"{new_rid}_ub",
                    numpy.NINF,
                    0.0,
                    components=[
                        (1, new_rid, "linear"),
                        (
                            1 * dt * reaction.getUpperBound(),
                            f"BM_{mid}_exchange",
                            "linear",
                        ),
                    ],
                )

                self.m_model.addUserDefinedConstraint(udc)
                reaction.setUpperBound(cbmpy.INF)

        for rid in combined_set:
            reaction: Reaction = initial_model.getReaction(rid)
            lb = reaction.getLowerBound()
            ub = reaction.getUpperBound()
            mid = self.m_model.identify_model_from_reaction(rid)

            for i in range(1, n):
                new_rid = rid + "_" + self.m_times[i]
                reactionN: Reaction = self.m_model.getReaction(new_rid)

                # Amount of biomass at time n
                r_x_t = f"BM_{mid}_{self.m_times[i-1]}_{self.m_times[i]}"

                if rid in rids_lb_to_check:
                    udc = self.m_model.createUserDefinedConstraint(
                        f"{new_rid}_lb",
                        0.0,
                        numpy.Inf,
                        components=[
                            (1, new_rid, "linear"),
                            (-1 * dt * lb, r_x_t, "linear"),
                        ],
                    )
                    self.m_model.addUserDefinedConstraint(udc)
                    reactionN.setLowerBound(cbmpy.NINF)
                if rid in rids_ub_to_check:
                    udc = self.m_model.createUserDefinedConstraint(
                        f"{new_rid}_ub",
                        numpy.NINF,
                        0.0,
                        components=[
                            (1, new_rid, "linear"),
                            (-1 * dt * ub, r_x_t, "linear"),
                        ],
                    )

                    self.m_model.addUserDefinedConstraint(udc)
                    reactionN.setUpperBound(cbmpy.INF)

    # For cbmpy < 0.9.0
    # def set_constraints(
    #     self,
    #     initial_model: CommunityModel,
    #     n: int,
    #     initial_biomasses: dict[str, float],
    #     dt: float,
    # ):
    #     """
    #     Configures the constraints for the EndPointFBA model. Unlike using
    #     upper and lower bounds, it adjusts constraints for each time point's
    #     reaction based on biomass, dt, and the initial bound.

    #     Args:
    #         n (int): Number of time points.
    #         initial_biomasses (dict[str, float]): Dictionary mapping model ID
    #             to initial biomass concentrations.
    #         dt (float): Time step size.
    #     """
    #     # Lookup time of set is on average O(1)
    #     rids_lb_to_check = set()
    #     rids_ub_to_check = set()
    #     for reaction in initial_model.reactions:
    #         if reaction.is_exchange:
    #             continue
    #         lb = reaction.getLowerBound()
    #         ub = reaction.getUpperBound()
    #         if not (lb == 0.0 or lb == cbmpy.INF or lb == cbmpy.NINF):
    #             rids_lb_to_check.add(reaction.getId())
    #         if not (ub == 0.0 or ub == cbmpy.INF or ub == cbmpy.NINF):
    #             rids_ub_to_check.add(reaction.getId())

    #     combined_set = rids_lb_to_check | rids_ub_to_check
    #     # Reactions at time zero
    #     for rid in combined_set:
    #         new_rid = rid + "_" + self.m_times[0]
    #         reaction = self.m_model.getReaction(new_rid)
    #         mid = self.m_model.identify_model_from_reaction(rid)
    #         biomass = initial_biomasses[mid]
    #         if rid in rids_lb_to_check:
    #             reaction.setLowerBound(reaction.getLowerBound() * biomass * dt)
    #         if rid in rids_ub_to_check:
    #             reaction.setUpperBound(reaction.getUpperBound() * biomass * dt)

    #     for rid in combined_set:
    #         reaction: Reaction = initial_model.getReaction(rid)
    #         lb = reaction.getLowerBound()
    #         ub = reaction.getUpperBound()
    #         mid = self.m_model.identify_model_from_reaction(rid)

    #         for i in range(1, n):
    #             new_rid = rid + "_" + self.m_times[i]
    #             reactionN: Reaction = self.m_model.getReaction(new_rid)

    #             r_x_t = f"BM_{mid}_{self.m_times[i-1]}_{self.m_times[i]}"

    #             if rid in rids_lb_to_check:
    #                 self.m_model.addUserConstraint(
    #                     f"{new_rid}_lb",
    #                     [
    #                         [1, new_rid],
    #                         [-1 * dt * lb, r_x_t],
    #                     ],
    #                     ">=",
    #                     0.0,
    #                 )
    #                 reactionN.setLowerBound(cbmpy.NINF)
    #             if rid in rids_ub_to_check:
    #                 self.m_model.addUserConstraint(
    #                     f"{new_rid}_ub",
    #                     [
    #                         [1, new_rid],
    #                         [-1 * dt * ub, r_x_t],
    #                     ],
    #                     "<=",
    #                     0.0,
    #                 )
    #                 reactionN.setUpperBound(cbmpy.INF)

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
        Raises:
            SpeciesNotFound: If the species defined in the initial_concentrations
                dictionaries keys is not in the model raise an exception

        """
        for key, value in initial_concentrations.items():
            sid = key + "_" + self.m_times[0]
            if sid not in self.m_model.getSpeciesIds():
                raise SpeciesNotFound(
                    "The species id defined as  \
                                      initial concentrations was not found in the model"
                )

            # get species and it's corresponding exchange reaction
            species: Species = self.m_model.getSpecies(sid)
            rids = species.isReagentOf()
            for rid in rids:
                reaction: Reaction = self.m_model.getReaction(rid)
                if reaction.is_exchange:
                    reaction.setLowerBound(-value)
                    reaction.setUpperBound(-value)

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

    # TODO fix this one
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

    def balanced_growth(self, Xin, Xm):
        self.m_model.setReactionBounds("X_comm", Xm, Xm)
        additional_components = []
        # At least some of each microbe needs to get into the system
        for mid, _ in self.m_model.get_model_biomass_ids().items():
            self.m_model.setReactionBounds(
                f"BM_{mid}_exchange", cbmpy.NINF, 0.0
            )

        for mid, _ in self.m_model.get_model_biomass_ids().items():
            self.m_model.createReaction(
                f"Phi_{mid}",
                f"Phi, fraction of {mid}",
                create_default_bounds=False,
            )
            additional_components.append((1.0, f"Phi_{mid}", "linear"))
            udc = self.m_model.createUserDefinedConstraint(
                f"biomass_fraction_{mid}_{self.m_times[0]}",
                0.0,
                0.0,
                components=[
                    (-1.0, f"BM_{mid}_exchange", "linear"),
                    (-1.0 * Xin, f"Phi_{mid}", "linear"),
                ],
            )

            self.m_model.addUserDefinedConstraint(udc)

            udc = self.m_model.createUserDefinedConstraint(
                f"biomass_fraction_{mid}_{self.m_times[-1]}",
                0.0,
                0.0,
                components=[
                    (1.0, f"BM_{mid}_exchange_final", "linear"),
                    (-1.0 * (Xm + Xin), f"Phi_{mid}", "linear"),
                ],
            )

            self.m_model.addUserDefinedConstraint(udc)

        udc = self.m_model.createUserDefinedConstraint(
            "Phi_add_to_one", 1.0, 1.0, components=additional_components
        )

        self.m_model.addUserDefinedConstraint(udc)

    # TODO In construction
    def remove_balanced_growth_constraints(self, initial_biomasses={}):
        """Restore the EndPointFBA model to before balanced growth constraints
            were added

        Args:
            initial_biomasses (dict, optional): _description_. Defaults to {}.
        """
        if initial_biomasses:
            for mid, value in initial_biomasses.items():
                self.m_model.setReactionBounds(
                    f"BM_{mid}_exchange", -value, -value
                )
        for mid, _ in self.m_model.get_model_biomass_ids().items():
            self.m_model.deleteReactionAndBounds(f"Phi_{mid}")
            self.m_model.__popGlobalId__(
                f"biomass_fraction_{mid}_{self.m_times[0]}"
            )
            self.m_model.__popGlobalId__(
                f"biomass_fraction_{mid}_{self.m_times[-1]}"
            )

    def set_qp(self, solution: float, epsilon=0.01) -> None:
        obj = self.m_model.getActiveObjective()
        obj.setOperation("minimize")
        obj.deleteAllFluxObjectives()
        QP = []

        for i, tid in enumerate(self.m_times[:-1]):
            rids = self.m_model.getReactionIds(tid)
            for rid in rids:
                reaction = self.m_model.getReaction(rid)

                # R_ are all the original reqactions
                # We dont want the exchange reactions
                if rid.startswith("R_") and not reaction.is_exchange:
                    rid = re.match(r"(.*?)_(time\d*)", rid).group(1)
                    rid_t_t1 = f"{rid}_{self.m_times[i+1]}"
                    rid_t_t0 = f"{rid}_{self.m_times[i]}"

                    if i == 0:
                        QP.append([1.0 * 2, rid_t_t0, rid_t_t0, str(i)])
                    else:
                        QP.append([2.0 * 2, rid_t_t0, rid_t_t0, str(i)])
                    if i == len(self.m_times[:-1]) - 1:
                        QP.append([1.0 * 2, rid_t_t1, rid_t_t1, str(i)])

                    QP.append([-2.0, rid_t_t0, rid_t_t1, str(i)])

        obj.createQuadraticFluxObjectives(QP)

        self.m_model.getReaction("X_comm").setLowerBound(solution - epsilon)
        self.m_model.getReaction("X_comm").setUpperBound(solution + epsilon)

    def set_subset_qp(
        self, solution: float, reactions: list[str], epsilon=0.01
    ) -> None:
        obj = self.m_model.getActiveObjective()
        obj.setOperation("minimize")
        obj.deleteAllFluxObjectives()
        QP = []
        all_reactions = [
            f"{r}_{mid}"
            for mid in self.m_model.get_model_ids()
            for r in reactions
        ]

        for i, _ in enumerate(self.m_times[:-1]):
            for rid in all_reactions:
                rid_t_t1 = f"{rid}_{self.m_times[i+1]}"
                rid_t_t0 = f"{rid}_{self.m_times[i]}"

                if i == 0:
                    QP.append([1.0 * 2, rid_t_t0, rid_t_t0, str(i)])
                else:
                    QP.append([2.0 * 2, rid_t_t0, rid_t_t0, str(i)])
                if i == len(self.m_times[:-1]) - 1:
                    QP.append([1.0 * 2, rid_t_t1, rid_t_t1, str(i)])

                QP.append([-2.0, rid_t_t0, rid_t_t1, str(i)])

        obj.createQuadraticFluxObjectives(QP)

        self.m_model.getReaction("X_comm").setLowerBound(solution - epsilon)
        self.m_model.getReaction("X_comm").setUpperBound(solution + epsilon)
