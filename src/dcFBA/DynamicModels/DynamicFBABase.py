import cbmpy
import math
from cbmpy.CBModel import Reaction
from ..Models import KineticsStruct, CommunityModel
from .StaticOptimizationModel import StaticOptimizationModelBase


class DynamicFBABase(StaticOptimizationModelBase):
    """
    Base class for SingleDynamicFBA and JointFBA.

    Implements the core code for simulating both models. The distinction
    between the two is the presence of a Community Biomass in Dynamic JointFBA.
    """

    def __init__(
        self,
        model: CommunityModel,
        biomasses: list[float],
        initial_concentrations: dict[str, float] = {},
        kinetics: KineticsStruct = KineticsStruct({}),
    ):
        """
        Initialize the DynamicFBABase class.

        Args:
            model (CommunityModel): The community model for simulation.
            biomasses (list[float]): List of initial biomass concentrations for each organism.
            initial_concentrations (dict[str, float], optional): Initial metabolite concentrations. Defaults to an empty dictionary.
            kinetics (KineticsStruct, optional): Kinetic parameters for reactions. Defaults to an empty KineticsStruct object.
        """
        super().__init__()
        self._model = model.clone()

        model_biomasses = model.get_model_biomass_ids()

        initial_biomasses = [[x] for x in biomasses]

        self._biomasses = dict(zip(model_biomasses.keys(), initial_biomasses))

        self._kinetics = kinetics
        self._set_initial_concentrations(self.model, initial_concentrations)

        for rid in self.model.getReactionIds():
            reaction: Reaction = self.model.getReaction(rid)
            self._initial_bounds[rid] = [
                reaction.getLowerBound(),
                reaction.getUpperBound(),
            ]

    @property
    def model(self) -> CommunityModel:
        return self._model

    def get_flux_values(self, rid: str) -> list[float]:
        """
        Get the flux values for a specific reaction over time.

        Args:
            rid (str): The reaction ID.

        Returns:
            list[float]: List of flux values for the specified reaction over time.
        """

        return list(map(lambda d: d[rid], self.fluxes))

    def get_fluxes_values(self, rids: list[str]) -> dict[str, list[float]]:
        """
        Get the flux values for multiple reactions over time.

        Args:
            rids (list[str]): List of reaction IDs.

        Returns:
            dict[str, list[float]]: Dictionary mapping reaction IDs to lists of flux values over time.
        """
        fluxes = {}
        for rid in rids:
            fluxes[rid] = self.get_flux_values(rid)
        return fluxes

    def get_specific_flux_values(self, rid: str) -> list[float]:
        """
        Get the specific flux values for a given reaction.

        Args:
            rid (str): The reaction ID.

        Returns:
            list[float]: List of specific flux values for the specified
                reaction over time. Returns an empty list if the reaction is
                an exchange reaction.
        """

        if rid in self.model.getExchangeReactionIds():
            print("Exchange has no specific flux")
            return []
        ls = self.get_flux_values(rid)
        mid = self.model.identify_model_from_reaction(rid)
        return [ls[t] / v for t, v in enumerate(self.biomasses[mid][:-1])]

    def simulate(
        self,
        dt: float,
        n: int = 10000,
        epsilon=1e-6,
        kinetics_func=None,
        deviate=None,
    ) -> None:
        """
        Perform a dynamic joint FBA simulation.
        """

        used_time = [0]
        dt_hat = -1
        dt_save = dt
        fluxes = []
        run_condition = 0

        for _ in range(1, n):
            if dt_hat != -1:
                dt = dt_hat
                dt_hat = -1
            else:
                dt = dt_save

            if deviate is not None:
                run_condition += deviate(
                    self,
                    used_time,
                    run_condition,
                )

            self.update_reaction_bounds(kinetics_func)
            self.update_exchanges(dt)

            solution = cbmpy.doFBA(self.model, quiet=False)

            if math.isnan(solution) or solution <= epsilon or dt < epsilon:
                break

            FBAsol = self.model.getSolutionVector(names=True)
            FBAsol = dict(zip(FBAsol[1], FBAsol[0]))

            used_time.append(used_time[-1] + dt)

            fluxes.append(FBAsol)

            self.update_concentrations(FBAsol, dt)

            for _, rid in self.model.get_model_biomass_ids().items():
                mid = self.model.identify_model_from_reaction(rid)
                Xt = self.biomasses[mid][-1] + FBAsol[rid] * dt
                self.biomasses[mid].append(Xt)

        self._fluxes = fluxes
        self._times = used_time

    def update_exchanges(self, dt: float) -> None:
        """
        Update exchange reaction lower bounds based on the metabolite
        concentration of the current time step.

        Args:
            dt (float): The time step for the simulation.
        """

        for rid in self.model.getExchangeReactionIds():
            reaction: Reaction = self.model.getReaction(rid)
            # Exchanges only have one species:
            sid = reaction.getSpeciesIds()[0]
            # TODO DISCUSS  (1/dt)
            # How I explain it: We normalize the exchange flux for how much the exchange can take up
            # in 1 unit of time. In all other formulas we multiple by dt, making sure that the flux gets scaled to
            # what it can take up in dt time.
            reaction.setLowerBound(
                min(0, -self.metabolites[sid][-1] * (1 / dt))
            )

    def update_concentrations(
        self, FBAsol: dict[str, float], dt: float
    ) -> None:
        """
        Update metabolite concentrations after an FBA simulation step.

        Args:
            FBAsol (dict): The solution vector from the FBA.
            dt (float): The time step for the simulation.
        """

        for e in self.model.getExchangeReactionIds():
            exchange: Reaction = self.model.getReaction(e)

            sid = exchange.getSpeciesIds()[0]
            if sid not in self.model.m_single_model_biomass_reaction_ids:
                self.metabolites[sid].append(
                    self.metabolites[sid][-1] + FBAsol[e] * dt
                )

    def update_reaction_bounds(self, kinetics_func) -> None:
        """
        Update the reaction bounds for the simulation based on either provided
        kinetics or standard bounds.

        Args:
            kinetics_func (function): A custom function to adjust reaction
                bounds based on kinetics. Uses standard bounds if None.
        """
        for rid in self.model.getReactionIds():
            # TODO fix this
            if kinetics_func is not None:
                kinetics_func(
                    rid,
                    self.model,
                    self.metabolites,
                    self.metabolites,
                )
                continue
            reaction: Reaction = self.model.getReaction(rid)

            # Don't change the exchange reactions bounds
            if not reaction.is_exchange:
                mid_for_reaction = self.model.identify_model_from_reaction(rid)
                # Organism specific biomass at last time point
                X_k_t = self.biomasses[mid_for_reaction][-1]
                reaction.setLowerBound(self.initial_bounds[rid][0] * X_k_t)

                if self.kinetics.exists(rid):
                    self.mm_kinetics(reaction, X_k_t, self.kinetics)

                else:
                    reaction.setUpperBound(self.initial_bounds[rid][1] * X_k_t)

    # def reset_dt(self, species_id: str, FBAsol) -> float:
    #     """
    #     Adjust the time step if a species becomes infeasible during the simulation.

    #     Args:
    #         species_id (str): The ID of the infeasible species.
    #         FBAsol (dict): The solution vector from the FBA.

    #     Returns:
    #         float: The recalculated time step.
    #     """

    #     # Remove last metabolite concentration
    #     self.m_metabolite_concentrations = {
    #         key: lst[:-1]
    #         for key, lst in self.m_metabolite_concentrations.items()
    #     }

    #     # Maybe some metabolite was exported/created, by an organism
    #     # Unfortunately to check what is create we would
    #     # have to know the new dt which we don'...
    #     # So we accept a small error.

    #     total = 0
    #     for (
    #         rid,
    #         species_ids,
    #     ) in self.m_transporters.get_importers(True):
    #         if species_id in species_ids:
    #             total += FBAsol[rid]

    #     return self.m_metabolite_concentrations[species_id][-1] / total
