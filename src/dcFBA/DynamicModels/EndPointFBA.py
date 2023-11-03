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
    on a CommunityModel. Inherits from the DynamicModelBase class.
    """

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
        super().__init__()

        width = len(str(n))
        self._times = [f"time{i:0{width}d}" for i in range(n)]
        self._dt = dt

        self._model = build_time_model(community_model, self.times)

        self._kinetics = kinetics

        self._set_objective()

        self._set_constraints(community_model, n, dt)

        self._set_initial_concentrations(
            initial_biomasses, initial_concentrations
        )

    @property
    def model(self) -> CommunityModel:
        """Returns the community model used in the simulation."""
        return self._model

    @property
    def dt(self) -> CommunityModel:
        """Returns the size of the time step."""
        return self._dt

    def _set_fluxes(self) -> None:
        """Private method to set the fluxes from the model solution."""

        solution_vector = self.model.getSolutionVector(names=True)
        self._fluxes = dict(zip(solution_vector[1], solution_vector[0]))

    def _set_biomasses(self) -> None:
        """Private method to set the concentrations of Biomasses over time."""
        mids = self.model.get_model_ids()
        temp_biomasses: dict[str, list[float]] = {}
        for mid in mids:
            temp_biomasses[mid] = [-1 * self.fluxes[f"BM_{mid}_exchange"]]

        for mid in self.model.get_model_ids():
            for i in range(len(self.times[:-1])):
                temp_biomasses[mid].append(
                    self.fluxes[f"BM_{mid}_{self.times[i]}_{self.times[i+1]}"]
                )
            temp_biomasses[mid].append(self.fluxes[f"BM_{mid}_exchange_final"])

        self._biomasses = dict(temp_biomasses)
        del temp_biomasses

    def _set_metabolites(self) -> None:
        """Private method to set the metabolite concentrations."""

        pattern = rf"^(.*?)_{self.times[0]}"
        regex = re.compile(pattern)
        temp_metabolites: dict[str, list[float]] = {}

        # Iterate through exchange reaction IDs
        for eid in self.model.getExchangeReactionIds():
            # Get the species ID for the reaction
            species_id: str = self.model.getReaction(eid).getSpeciesIds()[0]

            # Search for the pattern in the species ID
            match = regex.search(species_id)
            if match is None:
                continue

            # Extract the old species ID from the match
            old_species_id = match.group(1)

            # Check if the old species ID has the desired prefix and update the metabolites dictionary if not
            if not old_species_id.startswith("BM_c"):
                temp_metabolites[old_species_id] = [-1 * self.fluxes[eid]]

        for sid in temp_metabolites.keys():
            for i in range(len(self.times[:-1])):
                temp_metabolites[sid].append(
                    self.fluxes[f"{sid}_{self.times[i]}_{self.times[i+1]}"]
                )
            temp_metabolites[sid].append(self.fluxes[f"{sid}_exchange_final"])

        keys_to_delete = [
            k for k, v in temp_metabolites.items() if sum(v) == 0
        ]

        for key in keys_to_delete:
            del temp_metabolites[key]

        self._metabolites = dict(temp_metabolites)

        del temp_metabolites

    def get_flux_values(self, rid: str) -> list[float]:
        """Returns the flux values for each time point given a reaction ID.

        Args:
            rid (str): Reaction id of the original model

        Returns:
            list[float]: aggregated flux values for each time point
        """
        fluxes: list[float] = []
        for tid in self.times:
            full_id = f"{rid}_{tid}"
            fluxes.append(self.fluxes[full_id])
        return fluxes

    def get_fluxes_values(self, rids: list[str]) -> dict[str, list[float]]:
        """Returns the flux values for a list of reaction IDs for each time point.

        Args:
            rids (list[str]): List of reaction ids for which you want the
                flux vales

        Returns:
            dict[str, list[float]]: dictionary containing reaction id
                and the flux values foreach time-point
        """
        fluxes: dict[str, list[float]] = {}
        for rid in rids:
            fluxes[rid] = self.get_flux_values(rid)
        return fluxes

    def get_specific_flux_values(self, rid: str) -> list[float]:
        """Returns specific flux values for a given reaction ID.
        Specific flux is defined by the aggregated flux divided
        by the time-step size times the biomass

        Args:
            rid (str): reaction id of the original model

        Returns:
            list[float]: specific flux values
        """
        values = self.get_flux_values(rid)
        mid = self.model.identify_model_from_reaction(rid)

        return [
            v / (self.dt * self.biomasses[mid][i])
            for i, v in enumerate(values)
        ]

    def get_community_growth_rate(self) -> list[float]:
        """Calculates and returns the community growth rate over time."

        Returns:
            list[float]: Community growth rate
        """
        total_flux = [0] * len(self.times)
        total_mass = [0] * len(self.times)

        for mid, bid in self.model.get_model_biomass_ids().items():
            total_flux = numpy.add(total_flux, self.get_flux_values(bid))
            total_mass = numpy.add(total_mass, self.biomasses[mid][:-1])

        # Multiply by dt
        numerator = total_mass * self.dt
        return numpy.divide(total_flux, numerator).tolist()

    def get_relative_abundance(self) -> dict[str, list[float]]:
        """Calculates and returns the relative abundance, or percentage of
        each species in the community for each time point.

        Returns:
            dict[str, list[float]]: model_id : relative abundance over time
        """
        mids = self.model.get_model_ids()
        total = [0] * (len(self.times) + 1)

        for mid in mids:
            total = numpy.add(total, self.biomasses[mid])

        return {mid: numpy.divide(self.biomasses[mid], total) for mid in mids}

    def simulate(self, sparse=False) -> None:
        """Performs FBA (Flux Balance Analysis) on the EndPointFBA matrix.

        Args:
            sparse (False): Set to true if you want to use a sparse matrix
                Sparse matrix decreases the amount of memory required
        """

        matrix_type = "scipy_csr" if sparse else "numpy"
        self.model.buildStoichMatrix(matrix_type=matrix_type)
        solution = cbmpy.doFBA(self.model, quiet=False, build_n=False)

        self._set_fluxes()
        self._set_biomasses()
        self._set_metabolites()

        return solution

    def _set_objective(self) -> None:
        """Creates the community biomass reaction and sets it
        to be the objective of the model"""

        self.model.createReaction("X_comm", silent=True)
        out: Reaction = self.model.getReaction("X_comm")
        out.is_exchange = True
        out.setUpperBound(cbmpy.INF)
        out.setLowerBound(0)
        out.createReagent("BM_c_" + self.times[-1], -1)

        self.model.createObjectiveFunction("X_comm")

        self.model.setActiveObjective("X_comm_objective")

    def fva(self, selected_reactions=None) -> None:
        self.model.buildStoichMatrix()
        return cbmpy.doFVA(self.model, selected_reactions=selected_reactions)

    def _set_constraints(
        self,
        initial_model: CommunityModel,
        n: int,
        dt: float,
    ) -> None:
        """
        Private method configures the constraints for the EndPointFBA model.
        Unlike using upper and lower bounds, it adjusts constraints for each
        time point's reaction based on biomass, dt, and the initial bound.

        Args:
            n (int): Number of time points.
            initial_biomasses (dict[str, float]): Dictionary mapping model ID
                to initial biomass concentrations.
            dt (float): Time step size.
        """

        # TODO can be removed once cbmpy version 0.9.0 is online
        self.model.__FBC_VERSION__ = 3

        # Lookup time of set is on average O(1)
        rids_lb_to_check = set()
        rids_ub_to_check = set()
        for reaction in initial_model.reactions:
            if reaction.is_exchange:
                continue
            lb = reaction.getLowerBound()
            ub = reaction.getUpperBound()
            if not (lb == 0.0 or lb == cbmpy.INF or lb == numpy.NINF):
                rids_lb_to_check.add(reaction.getId())
            if not (ub == 0.0 or ub == cbmpy.INF or ub == numpy.NINF):
                rids_ub_to_check.add(reaction.getId())

        combined_set = rids_lb_to_check | rids_ub_to_check
        # Reactions at time zero
        for rid in combined_set:
            # TODO under construction
            # if self.kinetics and self.kinetics.exists(reaction.getId()):
            #     self.mm_approximation(reaction.getId())
            #     continue
            new_rid = rid + "_" + self.times[0]
            reaction = self.model.getReaction(new_rid)
            mid = self.model.identify_model_from_reaction(rid)

            if rid in rids_lb_to_check:
                udc = self.model.createUserDefinedConstraint(
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
                self.model.addUserDefinedConstraint(udc)
                reaction.setLowerBound(numpy.NINF)
            if rid in rids_ub_to_check:
                udc = self.model.createUserDefinedConstraint(
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

                self.model.addUserDefinedConstraint(udc)
                reaction.setUpperBound(cbmpy.INF)

        for rid in combined_set:
            reaction: Reaction = initial_model.getReaction(rid)
            # TODO not production ready yet
            # if self.kinetics and self.kinetics.exists(reaction.getId()):
            #     # just continue here, all bounds where set on the previous step
            #     continue

            lb = reaction.getLowerBound()
            ub = reaction.getUpperBound()
            mid = self.model.identify_model_from_reaction(rid)

            for i in range(1, n):
                new_rid = rid + "_" + self.times[i]
                reactionN: Reaction = self.model.getReaction(new_rid)

                # Amount of biomass at time n
                r_x_t = f"BM_{mid}_{self.times[i-1]}_{self.times[i]}"

                if rid in rids_lb_to_check:
                    udc = self.model.createUserDefinedConstraint(
                        f"{new_rid}_lb",
                        0.0,
                        numpy.Inf,
                        components=[
                            (1, new_rid, "linear"),
                            (-1 * dt * lb, r_x_t, "linear"),
                        ],
                    )
                    self.model.addUserDefinedConstraint(udc)
                    reactionN.setLowerBound(numpy.NINF)
                if rid in rids_ub_to_check:
                    udc = self.model.createUserDefinedConstraint(
                        f"{new_rid}_ub",
                        numpy.NINF,
                        0.0,
                        components=[
                            (1, new_rid, "linear"),
                            (-1 * dt * ub, r_x_t, "linear"),
                        ],
                    )

                    self.model.addUserDefinedConstraint(udc)
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

    def _set_initial_concentrations(
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
            sid = key + "_" + self.times[0]
            if sid not in self.model.getSpeciesIds():
                raise SpeciesNotFound(
                    "The species id defined as  \
                                      initial concentrations was not found in the model"
                )

            # get species and it's corresponding exchange reaction
            species: Species = self.model.getSpecies(sid)
            rids = species.isReagentOf()
            for rid in rids:
                reaction: Reaction = self.model.getReaction(rid)
                if reaction.is_exchange:
                    reaction.setLowerBound(-value)
                    reaction.setUpperBound(-value)

        for key, value in initial_biomasses.items():
            self.model.setReactionBounds(f"BM_{key}_exchange", -value, -value)

    def constrain_rates(self, epsilon=0.1):
        """
        Constrains the difference in reaction rates between t_n-1 and t_n.

        See the docs for further explanation

        Args:
            epsilon (float, optional): Maximum allowed rate difference between
                successive time points. Defaults to 0.1.
        """
        old_rids = set(
            [id.split("_time")[0] for id in self.model.getReactionIds()]
        )
        for rid in old_rids:
            reaction: Reaction = self.model.getReaction(
                rid + "_" + self.times[0]
            )

            if reaction is not None and (not reaction.is_exchange):
                for i, time in enumerate(self.times[:-1]):
                    id_t0 = f"{rid}_{time}"
                    id_t1 = f"{rid}_{self.times[i+1]}"

                    udc = self.model.createUserDefinedConstraint(
                        "R_constraint_pos" + id_t0,
                        numpy.NINF,
                        epsilon,
                        components=[
                            (1, id_t0, "linear"),
                            (-1, id_t1, "linear"),
                        ],
                    )

                    self.model.addUserDefinedConstraint(udc)

                    udc = self.model.createUserDefinedConstraint(
                        "R_constraint_neg" + id_t0,
                        -epsilon,
                        numpy.Inf,
                        components=[
                            (1, id_t0, "linear"),
                            (-1, id_t1, "linear"),
                        ],
                    )

                    self.model.addUserDefinedConstraint(udc)

    def mm_approximation(self, rid: str):
        """
        Approximates the Michaelis-Menten curve for a given reaction using two
        linear lines.

        When the EndPointFBA model is initialized with a `KineticsStruct` object, this method
        approximates the Michaelis-Menten (MM) kinetics of a given reaction by using two linear lines
        instead of directly setting an upper and lower bound.

        See official documentation for a comprehensive explanation of this approximation method.

        Parameters:
            rid (str): The reaction ID for which the MM approximation is to be applied.

        Raises:
            Exception: If no limiting substrate is defined in the `kinetics` object for the specified reaction.
        """
        print("WARNING not production ready")
        sid, km, vmax = self.kinetics.get_reactions_kinetics(rid)
        if sid == "":
            raise Exception(
                "No limiting substrate was set in the kinetics object for: "
                + rid
            )
        low_line = vmax / km
        high_line = (vmax / 2) / km

        for i in range(0, len(self.times) - 1):
            # Linking reaction is the concentration of Substrate for the
            # timepoint
            linking_reaction_id = f"{sid}_{self.times[i]}_{self.times[i+1]}"
            t_rid = rid + "_" + self.times[i + 1]

            udc = self.model.createUserDefinedConstraint(
                f"mm_low_{t_rid}",
                0.0,
                numpy.Inf,
                components=[
                    (1, t_rid, "linear"),
                    (
                        -low_line * self.dt,
                        linking_reaction_id,
                        "linear",
                    ),
                ],
            )

            self.model.addUserDefinedConstraint(udc)

            udc = self.model.createUserDefinedConstraint(
                f"mm_high_{t_rid}",
                numpy.NINF,
                0.0,
                components=[
                    (1, t_rid, "linear"),
                    (
                        -high_line * self.dt,
                        linking_reaction_id,
                        "linear",
                    ),
                ],
            )

            self.model.addUserDefinedConstraint(udc)

    def balanced_growth(self, Xin: float, Xm: float) -> None:
        """Set balanced growth constraint

        Args:
            Xin (float): Total community biomass at the first time point
            Xm (float): Total community biomass at the final time point
        """
        self.model.setReactionBounds("X_comm", Xm, Xm)
        additional_components = []

        for mid, _ in self.model.get_model_biomass_ids().items():
            self.model.setReactionBounds(f"BM_{mid}_exchange", numpy.NINF, 0.0)

        for mid, _ in self.model.get_model_biomass_ids().items():
            self.model.createReaction(
                f"Phi_{mid}",
                f"Phi, fraction of {mid}",
                create_default_bounds=False,
            )
            additional_components.append((1.0, f"Phi_{mid}", "linear"))
            udc = self.model.createUserDefinedConstraint(
                f"biomass_fraction_{mid}_{self.times[0]}",
                0.0,
                0.0,
                components=[
                    (-1.0, f"BM_{mid}_exchange", "linear"),
                    (-1.0 * Xin, f"Phi_{mid}", "linear"),
                ],
            )

            self.model.addUserDefinedConstraint(udc)

            udc = self.model.createUserDefinedConstraint(
                f"biomass_fraction_{mid}_{self.times[-1]}",
                0.0,
                0.0,
                components=[
                    (1.0, f"BM_{mid}_exchange_final", "linear"),
                    (-1.0 * (Xm + Xin), f"Phi_{mid}", "linear"),
                ],
            )

            self.model.addUserDefinedConstraint(udc)

        udc = self.model.createUserDefinedConstraint(
            "Phi_add_to_one", 1.0, 1.0, components=additional_components
        )

        self.model.addUserDefinedConstraint(udc)

    # TODO In construction
    def remove_balanced_growth_constraints(self, initial_biomasses={}):
        """Restore the EndPointFBA model to before balanced growth constraints
            were added

        Args:
            initial_biomasses (dict, optional): _description_. Defaults to {}.
        """
        print("WARNING, not ready")
        if initial_biomasses:
            for mid, value in initial_biomasses.items():
                self.model.setReactionBounds(
                    f"BM_{mid}_exchange", -value, -value
                )
        for mid, _ in self.model.get_model_biomass_ids().items():
            self.model.deleteReactionAndBounds(f"Phi_{mid}")
            self.model.__popGlobalId__(
                f"biomass_fraction_{mid}_{self.times[0]}"
            )
            self.model.__popGlobalId__(
                f"biomass_fraction_{mid}_{self.times[-1]}"
            )

    def set_qp(self, solution: float, epsilon=0.01) -> None:
        """Sets the quadratic objective to minimize
        all consecutive fluxes.

        Args:
           solution (float): Final community biomass

            epsilon (float, optional): How much the solution can differ from
                the final amount of biomass. Defaults to 0.01.
        """
        obj = self.model.getActiveObjective()
        obj.setOperation("minimize")
        obj.deleteAllFluxObjectives()
        QP = []

        for i, tid in enumerate(self.times[:-1]):
            rids = self.model.getReactionIds(tid)
            for rid in rids:
                reaction = self.model.getReaction(rid)

                # R_ are all the original reqactions
                # We dont want the exchange reactions
                if rid.startswith("R_") and not reaction.is_exchange:
                    rid = re.match(r"(.*?)_(time\d*)", rid).group(1)
                    rid_t_t1 = f"{rid}_{self.times[i+1]}"
                    rid_t_t0 = f"{rid}_{self.times[i]}"

                    if i == 0:
                        QP.append([1.0 * 2, rid_t_t0, rid_t_t0, str(i)])
                    else:
                        QP.append([2.0 * 2, rid_t_t0, rid_t_t0, str(i)])
                    if i == len(self.times[:-1]) - 1:
                        QP.append([1.0 * 2, rid_t_t1, rid_t_t1, str(i)])

                    QP.append([-2.0, rid_t_t0, rid_t_t1, str(i)])

        obj.createQuadraticFluxObjectives(QP)

        self.model.getReaction("X_comm").setLowerBound(solution - epsilon)
        self.model.getReaction("X_comm").setUpperBound(solution + epsilon)

    def set_subset_qp(
        self, solution: float, reactions: list[str], epsilon=0.01
    ) -> None:
        """QP for specified reaction ids

        Args:
            solution (float): Final community biomass
            reactions (list[str]): Reactions for which the consecutive
                fluxes are minimized
            epsilon (float, optional): How much the solution can differ from
                the final amount of biomass. Defaults to 0.01.
        """
        obj = self.model.getActiveObjective()
        obj.setOperation("minimize")
        obj.deleteAllFluxObjectives()
        QP = []
        all_reactions = [
            f"{r}_{mid}"
            for mid in self.model.get_model_ids()
            for r in reactions
        ]

        for i, _ in enumerate(self.times[:-1]):
            for rid in all_reactions:
                rid_t_t1 = f"{rid}_{self.times[i+1]}"
                rid_t_t0 = f"{rid}_{self.times[i]}"

                if i == 0:
                    QP.append([1.0 * 2, rid_t_t0, rid_t_t0, str(i)])
                else:
                    QP.append([2.0 * 2, rid_t_t0, rid_t_t0, str(i)])
                if i == len(self.times[:-1]) - 1:
                    QP.append([1.0 * 2, rid_t_t1, rid_t_t1, str(i)])

                QP.append([-2.0, rid_t_t0, rid_t_t1, str(i)])

        obj.createQuadraticFluxObjectives(QP)

        self.model.getReaction("X_comm").setLowerBound(solution - epsilon)
        self.model.getReaction("X_comm").setUpperBound(solution + epsilon)
