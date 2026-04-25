import cbmpy
import numpy
import math
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
    def dt(self) -> float:
        """Returns the size of the time step."""
        return self._dt

    def _set_fluxes(self) -> None:
        """Private method to set the fluxes from the model solution."""

        # solution_vector = self.model.getSolutionVector(names=True)
        # self._fluxes = dict(zip(solution_vector[1], solution_vector[0]))

        # Ensure N exists (scipy_csr: safer for big matrices)
        if self.model.N is None:
            self.model.buildStoichMatrix(matrix_type="scipy_csr")

        rids = self.model.N.col
        values = [self.model.getReaction(r).getValue() for r in rids]

        self._fluxes = dict(zip(rids, values))

    def _set_biomasses(self) -> None:
        """Private method to set the concentrations of Biomasses over time."""
        mids = self.model.custom_model_identifiers
        temp_biomasses: dict[str, list[float]] = {}
        for mid in mids:
            temp_biomasses[mid] = [-1 * self.fluxes[f"BM_{mid}_exchange"]]

        for mid in self.model.custom_model_identifiers:
            for i in range(len(self.times[:-1])):
                temp_biomasses[mid].append(
                    self.fluxes[f"BM_{mid}_{self.times[i]}_{self.times[i+1]}"]
                )
            temp_biomasses[mid].append(self.fluxes[f"BM_{mid}_exchange_final"])

        self._biomasses = dict(temp_biomasses)
        del temp_biomasses

    def _set_metabolites(self, keep_zero_metabolites = False, metabolites_to_keep = None) -> None:
        """
        Private method to set the metabolite concentrations.
        
        Args:
            keep_zero_metabolites (bool, optional): if True, also metabolites that are always zero are kept; default False.
            metabolites_to_keep (list, optional) = list of species ID of the metabolites to keep even if they are always zero
                                                    (in case keep_zero_metabolites == False). Default is None.
        """

        if metabolites_to_keep is None:
            metabolites_to_keep = []

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

        # delete metabolites with concentration always == 0.0 
        if not keep_zero_metabolites:
            keys_to_delete = [
                k for k, v in temp_metabolites.items() if sum(v) == 0
            ]

            for key in keys_to_delete:
                if key not in metabolites_to_keep:
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
        total_flux = numpy.zeros(len(self.times))
        total_mass = numpy.zeros(len(self.times))

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
        mids = self.model.custom_model_identifiers
        total = [0] * (len(self.times) + 1)

        for mid in mids:
            total = numpy.add(total, self.biomasses[mid])

        return {mid: numpy.divide(self.biomasses[mid], total) for mid in mids}

    def simulate(self, sparse=False) -> float:
        """Performs FBA (Flux Balance Analysis) on the EndPointFBA matrix.

        Args:
            sparse (False): Set to true if you want to use a sparse matrix
                Sparse matrix decreases the amount of memory required

         Returns:
            float: Final community flux value,
                NOTE that this is the total created community biomass.
                Since there is no initial community biomass in the FBA.
                To find the total community biomass you have to add
                the initial amount of biomass.
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

    # TODO fix this
    def mm_approximation(self, rid: str):
        """
        .. deprecated:: 0.2
        Use :meth:`mm_approximation_kinetic` instead.

        DEPRECATED: Use `mm_approximation_kinetic()` instead.

        Approximates the Michaelis-Menten curve for a given reaction using two
        linear constraints based on kinetics information.

        When the EndPointFBA model is initialized with a `KineticsStruct` object, this method
        approximates the Michaelis-Menten (MM) kinetics of a given reaction by using two linear constraints.

        See official documentation for a comprehensive explanation of this approximation method.

        Parameters:
            rid (str): The reaction ID for which the MM approximation is to be applied.

        Raises:
            Exception: If no limiting substrate is defined in the `kinetics` object for the specified reaction.
        """
        print("WARNING mm_approximation() is deprecated. "
        "Use mm_approximation_kinetic() instead.")
        #print("WARNING not production ready")
        sid, km, vmax = self.kinetics.get_reactions_kinetics(rid)
        if sid == "":
            raise Exception(
                "No limiting substrate was set in the kinetics object for: "
                + rid
            )
        low_line = vmax / km  # option 1
        high_line = (vmax / 2) / km  # Option 2

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

    def set_species_concentration_dependent_constraints(
        self,
        rid: str,
        sid: str,
        points: list[tuple]
    ) -> tuple[float, float]:
        """
        Add a new set of constraints for the EndPointFBA model.
        Defines a linear relationship between concentration of metabolite `sid`
        and flux through reaction `rid`, approximating kinetics (e.g., MM).

        A new upper bound for reaction 'rid' is added at every time point.
        The new constraint is defined as the line passing through the two points,
        in a plane with concentration of metabolite sid as x-axis 
        and flux through reaction rid as y-axis.   
        (For example, it can be used to approximate the slope of a Michealis-Menten kinetic.)
        See the documentation for more information and examples. 

        Args:
            rid (str): reaction ID of the flux to constrain (as it appears in the community model)
            sid (str): species ID of the metabolite whose concentration is affecting the flux
            points(list[tuple]): list of 2 tuple (c,v) where:
                                - 'c' is a concentrations of metabolites sid;
                                - 'v' is the maximum (aggregated) flux for reaction rid
                                    when concentration of sid is 'c'.
        
        Returns:
            tuple: A tuple containing the slope and offset of the constraint line.
        """
        print("WARNING not production ready")

        # error handling if more than 2 tuples in 'points' list
        if len(points) != 2:
            raise ValueError("The 'points' list must contain exactly 2 tuples.")

        (c1, v1), (c2, v2) = points

        # Compute slope and intercept of the constraint line
        m = float((v2 - v1) / (c2 - c1)) # slope
        q = float((v1 * c2 - v2 * c1) / (c2 - c1)) # offset

        # --- time 0
        # at time 0, the name is not matching the naming scheme, but 
        # only exhange reaction associated with the metabolite  
        for sid0 in self.model.getReactionIdsAssociatedWithSpecies(f"{sid}_{self.times[0]}"):
            if self.model.getReaction(sid0[0]).is_exchange:
                cid = sid0[0]
                break

        vid = rid + "_" + self.times[0]

        udc = self.model.createUserDefinedConstraint(
            f"{vid}_ub_{sid}",
            numpy.NINF,
            q,
            components=[(1, vid, "linear"), (m, cid, "linear")],
        )
        self.model.addUserDefinedConstraint(udc)

        # --- times 1..n-1
        # for all other times, the names match the naming scheme
        for i in range(len(self.times[:-1])):
            cid = f"{sid}_{self.times[i]}_{self.times[i + 1]}"
            vid = rid + "_" + self.times[i + 1]

            udc = self.model.createUserDefinedConstraint(
                f"{vid}_ub_{sid}",
                numpy.NINF,
                q,
                components=[(1, vid, "linear"), (-m, cid, "linear")],
            )
            self.model.addUserDefinedConstraint(udc)

        return m, q
    
    # TODO fix this
    def mm_approximated_kinetic(
        self,
        rid: str,
        conservative: bool = False,
        kinetics_params: tuple[str, float, float] | None = None,
    ) -> None:
        """
        Approximates the Michaelis-Menten curve for a given reaction at every time-point
        using two linear constraints based on kinetics parameters.

        
        This method applies two linear constraints to approximate the MM curve:
        1) a linear constraint v <= m*[S] that approximates
        the MM curve at low substrate concentrations, where the curve is 
        approximately linear with slope vmax/km.
        2) an horizontal constraint v <= vmax that represents the asymptotic behavior
        of the MM curve at high substrate concentrations.

        This method can use the `KineticsStruct` object associated with the EndPointFBA model
        to retrieve the limiting substrate, Km, and Vmax for the specified reaction.

        Warning: This method does not set the horizontal asymptote at vmax.
        Ensure the reaction upper bound is set to vmax in the community model 
        before building the EndPointFBA model.

        Parameters:
            rid (str): The reaction ID for which the MM approximation is to be applied.
            conservative (bool): If True, uses a more conservative slope (vmax/2km),
                                i.e., passing through (km, vmax/2). 
                                If False, uses the true low-substrate slope (vmax/km),
                                i.e., tangent to the MM curve at (0, 0).
                                Default is False.
            kinetics_params (tuple[str, float, float], optional):
                Optional tuple (sid, km, vmax). If provided, these values are used
                instead of retrieving them from `self.kinetics`.

        Raises:
            Exception: If no limiting substrate is defined (either in `kinetics_params`
                    or in the `kinetics` object for the specified reaction).

        Returns:
            None
        """
        print("WARNING not production ready")

        # --- Get kinetic parameters ---
        if kinetics_params is not None:
            # User-supplied kinetics (sid, km, vmax)
            sid, km, vmax = kinetics_params
        else:
            # Default: fetch from kinetics structure
            sid, km, vmax = self.kinetics.get_reactions_kinetics(rid)

        if not sid:
            raise Exception(f"No limiting substrate set in kinetics object or parameters for: {rid}")

        # --- Linear approximation for low substrate regime ---
        if conservative:
            # Line through (km, vmax/2)
            points = [(0, 0), (km, vmax / 2)]
        else:
            # Tangent at (0,0) with slope = vmax/km
            points = [(0, 0), (km, vmax)]

        # --- Apply the constraint ---
        self.set_species_concentration_dependent_constraints(rid, sid, points)

        # TODO --- Set upper bound to vmax (horizontal asymptote) ---
        print(
            "WARNING: setting upper bound to vmax not implemented yet. "
            "Ensure the reaction upper bound is set before building the EndPointFBA model."
        )

    def balanced_growth(self, Xin: float, Xm: float) -> None:
        """Set balanced growth constraint

        Args:
            Xin (float): Total community biomass at the first time point
            Xm (float): Total community biomass at the final time point
        """
        # There is no initial X_c, so subtract the initial biomass from the
        # final biomass
        community_flux = Xm - Xin
        self.model.setReactionBounds("X_comm", community_flux, community_flux)

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
                    (-1.0 * Xm, f"Phi_{mid}", "linear"),
                ],
            )

            self.model.addUserDefinedConstraint(udc)

        udc = self.model.createUserDefinedConstraint(
            "Phi_add_to_one", 1.0, 1.0, components=additional_components
        )

        self.model.addUserDefinedConstraint(udc)


    def binary_search_balanced_growth(self, Xin, Xfin, tolerance: float = 1e-6):
        """
        Perform a binary search to determine the maximum feasible balanced growth
        between initial (Xin) and final (Xfin) community biomass amounts.

        This method assumes that `balanced_growth()` has either already been called
        or that it can be safely called once to initialize the balanced growth
        constraints. It then updates those constraints iteratively to find the
        largest feasible biomass production.

        Args:
            Xin (float): Initial total community biomass.
            Xfin (float): Maximum total community biomass to test. Typically the
                result of an EndPointFBA run without balanced growth.
            tolerance (float, optional): Stopping tolerance for the binary search.
                Defaults to 1e-6.

        Note:
            The search will only explore values within [Xin, Xfin].
            If the model is feasible for Xfin, that value is returned directly.

        Returns:
            float: The highest feasible final community biomass (rounded to match
                the number of digits in the tolerance).
        """
        # set output precision based on tolerance
        #precision = len(str(int((1 / tolerance) - 1)))
        precision = abs(int(round(math.log10(1 / tolerance))))

        # --- Step 1: ensure balanced growth constraints exist

        try:
            # tries to add the balanced growth constraints
            self.balanced_growth(Xin, Xfin)
        except AssertionError:
            # if constraints already exist, update them to current Xin/Xfin
            self.modify_balanced_growth_constraints(Xin, Xfin)

        # --- Step 2: test initial feasibility

        # perform EndPointFBA with imput constraints
        res = self.simulate()

        # if infeasible, run binary search to find the highest feasible biomass production
        if math.isnan(res):
            community_growth = Xfin - Xin
            low, high = 0.0, community_growth # lower, upper bound for binary search

            # binary search loop
            while high - low > tolerance:
                mid = (low + high) / 2
                # update balanced growth constraints with new tentative Xfin
                self.modify_balanced_growth_constraints(Xin, Xin + mid)
                res = self.simulate()

                # if infeasible, all higher values will be excluded from the search
                if math.isnan(res):
                    high = mid
                # if feasible, all lower values will be excluded from the search
                else:
                    low = mid

            # last feasible value
            res = low

        # --- Step 3: finalize model with last feasible constraints
        self.modify_balanced_growth_constraints(Xin, Xin + res)
        res = self.simulate()

        # --- Step 4: round to tolerance precision
        res = round(res, precision)
        return res

    def modify_balanced_growth_constraints(self, Xin: float, Xfin: float):
        """
        Update balanced growth constraints for new initial (Xin)
        and final (Xfin) total biomass values.

        Args:
            Xin (float): Initial total community biomass.
            Xfin (float): Final total community biomass.
        """
        # update community biomass production constraint
        community_growth = Xfin - Xin
        self.model.setReactionBounds("X_comm", community_growth, community_growth)

        # update coefficients in existing UDCs
        for mid, _ in self.model.get_model_biomass_ids().items():
            # initial biomass fraction constraint
            udc_init = self.model.getObject(f"biomass_fraction_{mid}_{self.times[0]}")
            phi_init = udc_init.getConstraintComponent(
                f"udcc_biomass_fraction_{mid}_{self.times[0]}_Phi_{mid}"
            )
            phi_init.setCoefficient(-1.0 * Xin)

            # final biomass fraction constraint
            udc_final = self.model.getObject(f"biomass_fraction_{mid}_{self.times[-1]}")
            phi_final = udc_final.getConstraintComponent(
                f"udcc_biomass_fraction_{mid}_{self.times[-1]}_Phi_{mid}"
            )
            phi_final.setCoefficient(-1.0 * Xfin)


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

    # TODO maybe qp for other objectives than final biomass (last two lines)
    # TODO make QP using LP matrix directly

    # set QP on aggregate fluxes (of all reactions with ID starting with 'R_')
    def set_qp(self, solution: float, epsilon=0.0) -> None:
        """Sets the quadratic objective to minimize
        all consecutive fluxes.

        Args:
           solution (float): Solution for the final X_comm flux

            epsilon (float, optional): How much the solution can differ from
                the final amount of biomass. Defaults to 0.0.
        """
        obj = self.model.getActiveObjective()
        obj.setOperation("minimize")
        obj.deleteAllFluxObjectives()
        QP = []

        for i, tid in enumerate(self.times[:-1]):
            rids = self.model.getReactionIds(tid)
            for rid in rids:
                reaction = self.model.getReaction(rid)

                # R_ are all the original reactions
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
        self, solution: float, reactions: list[str], epsilon=0.0, all_models=False
    ) -> None:
        """
        Set up a quadratic programming (QP) objective to minimize flux changes
        for a specified subset of reactions across consecutive time points.

        This method is used to smooth flux transitions between time steps,penalizing 
        the squared difference between reaction fluxes in consecutive time intervals.

        Example:
            Apply QP on all biomass reactions, fix them, and apply QP on specific fluxes:

            ```python
            bm_reactions = [bid for _, bid in ep.model.get_model_biomass_ids().items()]
            ep.set_subset_qp(solution, bm_reactions)
            ep.simulate()  # run QP on biomasses

            biomasses = ep.biomasses
            ep.set_qp_specific_fluxes(solution, biomasses)
            ep.simulate()  # run QP on specific fluxes
            ```

            Note: the `solution` value should remain the one from the original
            endpoint problem (EP), not from the QP simulation.

        Args:
            solution (float): Target flux value for the final community biomass reaction (`X_comm`).
            reactions (list[str]): Reaction IDs for which consecutive fluxes
                are minimized.
            epsilon (float, optional): Allowed deviation from the final biomass
                flux (`X_comm`). Defaults to 0.0.
            all_models (bool, optional): 
                - If True, use the same `reactions` list for all models in the community.  
                - If False, `reactions` should already include model IDs 
                (e.g. `"R_rid_mid"`). Defaults to False.
        """

        obj = self.model.getActiveObjective()
        obj.setOperation("minimize")
        obj.deleteAllFluxObjectives()
        QP = []

        # select reactions
        if all_models:
            # same reactions across all model IDs
            all_reactions = [
                f"{r}_{mid}"
                for mid in self.model.get_model_ids()
                for r in reactions
            ]
        else:
            # reaction IDs directly
            all_reactions = reactions

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

    def set_qp_specific_fluxes(
        self, solution: float, biomasses: dict[str, list[float]], epsilon: float = 0.0
    ) -> None:
        """
        QP for intracellular reactions using the specific flux (= aggregated flux/biomass).
        Exchange and biomass reactions are excluded from the QP.

        Set up a quadratic programming (QP) objective to minimize temporal
        changes in *specific fluxes* (flux normalized by biomass) for all
        internal reactions across consecutive time points.

        The objective penalizes squared differences in specific fluxes:
            minimize Σ_t Σ_r [ (v_r,t+1 / B_t+1) - (v_r,t / B_t) ]²

        where:
            - v_r,t: flux of reaction r at time t
            - B_t: biomass of the corresponding model at time t

        Args:
            solution (float): Target flux value for the community biomass
                reaction (`X_comm`) at the final time step.
            biomasses (dict[str, list[float]]): Biomass values for each model
                at each time point (keys: model IDs; values: list of biomass
                values at consecutive time points).
            epsilon (float, optional): Allowed deviation from the final
                community biomass flux. Defaults to 0.0.

        Notes:
            - This method is typically run **after** a QP on biomass reactions
                (see `set_subset_qp`) to ensure smooth transitions in internal
                fluxes.
            - Example usage:

                ```python
                bm_reactions = [bid for _, bid in ep.model.get_model_biomass_ids().items()]
                ep.set_subset_qp(solution, bm_reactions)
                ep.simulate()  # run QP on biomasses

                biomasses = ep.biomasses
                ep.set_qp_specific_fluxes(solution, biomasses)
                ep.simulate()  # run QP on specific fluxes
                ```
        """
        obj = self.model.getActiveObjective()
        obj.setOperation("minimize")
        obj.deleteAllFluxObjectives()
        QP = []

        biomass_reactions_ids = [
            bid for _, bid in self.model.get_model_biomass_ids().items()
        ]

        for i, tid in enumerate(self.times[:-1]):
            rids = self.model.getReactionIds(tid)
            for rid in rids:
                reaction = self.model.getReaction(rid)

                # R_ are all the original reactions
                # We dont want the exchange reactions
                if rid.startswith("R_") and not reaction.is_exchange:
                    rid = re.match(r"(.*?)_(time\d*)", rid).group(1)
                    if rid not in biomass_reactions_ids: 
                        rid_t_t1 = f"{rid}_{self.times[i+1]}"
                        rid_t_t0 = f"{rid}_{self.times[i]}"
                        
                        # fix all biomass reactions to their value
                        # usa identify_model_from_reaction?
                        for mid, _ in self.model.get_model_biomass_ids().items():
                            if mid in rid: 
                                bm_t_t1 = biomasses[mid][i+1]
                                bm_t_t0 = biomasses[mid][i]
                                break

                        if i == 0:
                            bm_coeff = bm_t_t0 * bm_t_t0
                            QP.append([1.0 * 2 / bm_coeff, rid_t_t0, rid_t_t0, str(i)])
                        else:
                            bm_coeff = bm_t_t0 * bm_t_t0
                            QP.append([2.0 * 2 / bm_coeff, rid_t_t0, rid_t_t0, str(i)])
                        if i == len(self.times[:-1]) - 1:
                            bm_coeff = bm_t_t1 * bm_t_t1
                            QP.append([1.0 * 2 / bm_coeff, rid_t_t1, rid_t_t1, str(i)])

                        bm_coeff = bm_t_t0 * bm_t_t1
                        QP.append([-2.0 / bm_coeff, rid_t_t0, rid_t_t1, str(i)])

        obj.createQuadraticFluxObjectives(QP)

        self.model.getReaction("X_comm").setLowerBound(solution - epsilon)
        self.model.getReaction("X_comm").setUpperBound(solution + epsilon)

    def set_subset_qp_specific_fluxes(
        self,
        solution: float,
        reactions: list[str],
        biomasses: dict[str, list[float]],
        epsilon: float = 0.0,
        all_models: bool = False,
    ) -> None:
        """
        QP for asubset of intracellular reactions using the specific flux (= aggregated flux/biomass).
        Exchange and biomass reactions are excluded from the QP.

        Set up a quadratic programming (QP) objective to minimize temporal changes
        in *specific fluxes* (flux normalized by biomass) for a specified subset of
        reactions across consecutive time points.

        The QP penalizes squared differences between consecutive *specific fluxes*:
            minimize Σ_t Σ_r [ (v_r,t+1 / B_t+1) - (v_r,t / B_t) ]²

        Args:
            solution (float): Target flux value for the community biomass
                reaction (`X_comm`) at the final time step.
            reactions (list[str]): Reaction IDs to include in the QP objective.
            biomasses (dict[str, list[float]]): Biomass values for each model
                (keys: model IDs; values: list of biomass values at consecutive time points).
            epsilon (float, optional): Allowed deviation from the final community
                biomass flux. Defaults to 0.0.
            all_models (bool, optional):
                - If True, use the same `reactions` list for all models in the community.
                - If False, `reactions` should already include model IDs
                (e.g., "R_rid_mid"). Defaults to False.

        Notes:
            - Internal and exchange reactions are not distinguished here; the
            provided `reactions` determine which fluxes are optimized.
            - Typically run after a biomass QP (`set_subset_qp`) for smoother internal fluxes.

        Example:
            ```python
            bm_reactions = [bid for _, bid in ep.model.get_model_biomass_ids().items()]
            ep.set_subset_qp(solution, bm_reactions)
            ep.simulate()  # run QP on biomasses

            biomasses = ep.biomasses
            subset_reactions = ["R_GLYCOLYSIS_mid1", "R_TCA_mid2"]
            ep.set_subset_qp_specific_fluxes(solution, subset_reactions, biomasses)
            ep.simulate()  # run QP on specific fluxes for the subset
            ```
        """

        obj = self.model.getActiveObjective()
        obj.setOperation("minimize")
        obj.deleteAllFluxObjectives()
        QP = []

        # Select reactions
        if all_models:
            all_reactions = [
                f"{r}_{mid}"
                for mid in self.model.get_model_ids()
                for r in reactions
            ]
        else:
            all_reactions = reactions

        biomass_reactions_ids = [
            bid for _, bid in self.model.get_model_biomass_ids().items()
        ]

        for i, _ in enumerate(self.times[:-1]):
            for rid in all_reactions:
                # skip biomass reactions if accidentally included
                if any(bid in rid for bid in biomass_reactions_ids):
                    continue

                rid_t_t1 = f"{rid}_{self.times[i+1]}"
                rid_t_t0 = f"{rid}_{self.times[i]}"

                # fix all biomass reactions to their value
                # TODO usa identify_model_from_reaction?
                for mid, _ in self.model.get_model_biomass_ids().items():
                    if mid in rid:
                        bm_t_t1 = biomasses[mid][i+1]
                        bm_t_t0 = biomasses[mid][i]
                        break

                if i == 0:
                    bm_coeff = bm_t_t0 * bm_t_t0
                    QP.append([1.0 * 2 / bm_coeff, rid_t_t0, rid_t_t0, str(i)])
                else:
                    bm_coeff = bm_t_t0 * bm_t_t0
                    QP.append([2.0 * 2 / bm_coeff, rid_t_t0, rid_t_t0, str(i)])
                if i == len(self.times[:-1]) - 1:
                    bm_coeff = bm_t_t1 * bm_t_t1
                    QP.append([1.0 * 2 / bm_coeff, rid_t_t1, rid_t_t1, str(i)])

                bm_coeff = bm_t_t0 * bm_t_t1
                QP.append([-2.0 / bm_coeff, rid_t_t0, rid_t_t1, str(i)])

        obj.createQuadraticFluxObjectives(QP)

        self.model.getReaction("X_comm").setLowerBound(solution - epsilon)
        self.model.getReaction("X_comm").setUpperBound(solution + epsilon)
