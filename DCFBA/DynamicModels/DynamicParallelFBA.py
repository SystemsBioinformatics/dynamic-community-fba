import cbmpy
import math
from cbmpy.CBModel import Model, Reaction
from .StaticOptimizationModel import StaticOptimizationModelBase
from ..Models.Kinetics import KineticsStruct
from ..Models.Transporters import Transporters


class DynamicParallelFBA(StaticOptimizationModelBase):
    """A class representing a dynamic parallel FBA simulation.

    This class extends the TimeStepDynamicFBABase and provides functionality
    for performing dynamic FBA simulations on multiple models in parallel.

    Attributes:
    m_models (list[Model]): List of metabolic models to simulate in parallel.
    m_initial_bounds (dict[str, dict[str, tuple[float, float]]]): Original reaction
        bounds for each model, indexed by model ID and reaction ID.
    m_model_transporters (dict[str, Transporters]): Mapping of model ID to associated
        transporters.
    m_model_kinetics (dict[str, KineticsStruct]): Mapping of model ID to associated
        kinetic parameters for reactions.

    """

    m_models: list[Model]
    m_initial_bounds: dict[str, dict[str, tuple[float, float]]] = {}
    m_model_transporters: dict[str, Transporters] = {}
    m_model_kinetics: dict[str, KineticsStruct] = {}

    def __init__(
        self,
        models: list[Model],
        biomasses: list[float],
        initial_concentrations: dict[str, float] = {},
        kinetics: dict[str, KineticsStruct] = {},
    ) -> None:
        """
        Initialize the `DynamicParallelFBA` instance.

        Args:
            models (list[Model]): Metabolic models to simulate.
            biomasses (list[float]): Initial biomass concentrations for each model.
            initial_concentrations (dict[str, float], optional): Mapping of metabolite
                IDs to their initial concentrations. Defaults to an empty dictionary.
                If empty, exchange reaction lower bounds will be used as initial concentrations
            kinetics (dict[str, KineticsStruct], optional): Dictionary pairing model
                IDs with their respective kinetics structures. Defaults to an empty dictionary.
        """

        self.m_models = models
        self.m_model_kinetics = kinetics

        for i in range(0, len(models)):
            model = models[i]
            mid = model.getId()
            self.set_initial_concentrations(model, initial_concentrations)
            self.m_biomass_concentrations[mid] = [biomasses[i]]
            self.m_model_transporters[mid] = Transporters(model)

            for rid in model.getReactionIds():
                reaction: Reaction = model.getReaction(rid)
                if mid not in self.m_initial_bounds.keys():
                    self.m_initial_bounds[mid] = {
                        rid: [
                            reaction.getLowerBound(),
                            reaction.getUpperBound(),
                        ]
                    }
                else:
                    self.m_initial_bounds[mid][rid] = [
                        reaction.getLowerBound(),
                        reaction.getUpperBound(),
                    ]

    def simulate(
        self,
        dt: float,
        n: int = 10000,
        epsilon: float = 0.01,
        kinetics_func=None,
        deviate=None,
    ) -> list[list[float], dict[str, list[float]], dict[str, list[float]]]:
        """
        Execute a dynamic parallel FBA simulation.

        Args:
            dt (float): Time step interval for the simulation.
            epsilon (float, optional): Convergence threshold. Defaults to 0.001.
            kinetics_func (function, optional): Custom function to adjust reaction bounds
                based on kinetics. Defaults to None.
            deviate (function, optional): Function to add perturbations during the
                simulation. Defaults to None.

        Returns:
            list[list[float], dict[str, list[float]], dict[str, list[float]]]: A list containing
            the simulation results in the order of time utilized, metabolite concentrations,
            and biomass concentrations.
        """

        used_time = [0]
        dt_hat = -1
        dt_save = dt
        run_condition = 0
        step = 1

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

            used_time.append(used_time[-1] + dt)
            fluxes = {}
            for model in self.m_models:
                solution = cbmpy.doFBA(model, quiet=True)

                # TODO continue when the others can still grow??????? but if one is nan stop?
                if math.isnan(solution) or solution < epsilon or dt < epsilon:
                    return [
                        used_time[:-1],
                        self.m_metabolite_concentrations,
                        self.m_biomass_concentrations,
                    ]

                FBAsol = model.getSolutionVector(names=True)
                FBAsol = dict(zip(FBAsol[1], FBAsol[0]))

                fluxes[model.getId()] = FBAsol

                self.update_concentrations(model, FBAsol, step, dt)
                self.update_biomasses(model, dt, FBAsol, step)

            species_id = self.check_solution_feasibility()

            if species_id != "":
                dt_hat = self.reset_dt(species_id, fluxes)
                step -= 1
                used_time = used_time[:-1]

            step += 1
        return [
            used_time[:-1],
            self.m_metabolite_concentrations,
            self.m_biomass_concentrations,
        ]

    def update_exchanges(self, dt):
        """
        Adjust exchange reactions based on the present metabolite concentrations

        Sets the lower bounds for exchange reactions to prevent metabolite import beyond
        the current concentration.
        """

        for model in self.m_models:
            for rid in model.getExchangeReactionIds():
                reaction: Reaction = model.getReaction(rid)
                sid = reaction.getSpeciesIds()[0]
                reaction.setLowerBound(
                    -self.m_metabolite_concentrations[sid][-1] * (1 / dt)
                )

    def update_reaction_bounds(self, kinetics_func) -> None:
        """Update the reaction bounds for all models.

        Args:
            kinetics_func (function): Custom function to modify reaction bounds based on
                kinetic data.

        """
        models = self.m_models

        if kinetics_func is not None:
            # TODO Give somethings to the kinetics_func, think about this
            kinetics_func()
            return

        for model in models:
            mid = model.getId()
            for rid in model.getReactionIds():
                reaction: Reaction = model.getReaction(rid)
                if not reaction.is_exchange:
                    # Organism specific biomass
                    X_k_t = self.m_biomass_concentrations[mid][-1]
                    reaction.setLowerBound(
                        self.m_initial_bounds[mid][rid][0] * X_k_t
                    )

                    if (mid in self.m_model_kinetics.keys()) and (
                        self.m_model_kinetics[mid].exists(rid)
                    ):
                        self.mm_kinetics(
                            reaction,
                            X_k_t,
                            self.m_model_transporters[mid],
                            self.m_model_kinetics[mid],
                        )

                    else:
                        reaction.setUpperBound(
                            self.m_initial_bounds[mid][rid][1] * X_k_t
                        )

    def update_concentrations(self, model: Model, FBAsol, step, dt) -> None:
        """
        Modify metabolite concentrations following an FBA step.

        Args:
            model (Model): Current metabolic model being processed.
            FBAsol (dict): Solution vector resulting from the FBA.
            step (int): Current simulation step.
            dt (float): Time step interval for the simulation.
        """

        mid = model.getId()

        if len(next(iter(self.m_metabolite_concentrations.values()))) == step:
            for key in self.m_metabolite_concentrations.keys():
                self.m_metabolite_concentrations[key].append(
                    self.m_metabolite_concentrations[key][-1]
                )
        # TODO Use exchanges?
        for rid, species_ids in self.m_model_transporters[mid].get_exporters(
            True
        ):
            for sid in species_ids:
                self.m_metabolite_concentrations[sid][-1] = round(
                    self.m_metabolite_concentrations[sid][-1]
                    + FBAsol[rid] * dt,
                    8,
                )

        for rid, species_ids in self.m_model_transporters[mid].get_importers(
            True
        ):
            for sid in species_ids:
                self.m_metabolite_concentrations[sid][-1] = round(
                    self.m_metabolite_concentrations[sid][-1]
                    - FBAsol[rid] * dt,
                    8,
                )

    def update_biomasses(self, model: Model, dt, FBAsol, step) -> None:
        """
        Modify biomass concentrations after completing an FBA step.

        Args:
            model (Model): Model for which the biomass concentrations are updated.
            dt (float): Time step interval for the simulation.
            FBAsol (dict): Solution vector resulting from the FBA.
            step (int): Current simulation step.
        """

        mid = model.getId()
        Xt = (
            self.m_biomass_concentrations[mid][step - 1]
            + FBAsol[model.getActiveObjectiveReactionIds()[0]] * dt
        )
        self.m_biomass_concentrations[mid].append(Xt)

    def reset_dt(self, species_id, FBAsol) -> float:
        """
        Recompute a smaller time step if the obtained result was infeasible.

        Args:
            species_id (str): ID of the species for which the time step is recalculated.
            FBAsol (dict): Solution vector from the FBA.

        Returns:
            float: Adjusted time step.
        """

        self.m_metabolite_concentrations = {
            key: lst[:-1]
            for key, lst in self.m_metabolite_concentrations.items()
        }

        self.m_biomass_concentrations = {
            key: lst[:-1] for key, lst in self.m_biomass_concentrations.items()
        }
        total = 0
        for model in self.m_models:
            mid = model.getId()
            for rid, species_ids in self.m_model_transporters[
                mid
            ].get_importers(True):
                if species_id in species_ids:
                    total += FBAsol[mid][rid]

        return self.m_metabolite_concentrations[species_id][-1] / total
