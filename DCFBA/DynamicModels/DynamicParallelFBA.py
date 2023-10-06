import cbmpy
import math
from cbmpy.CBModel import Model, Reaction
from .StaticOptimizationModel import StaticOptimizationModelBase
from ..Models.Kinetics import KineticsStruct


class DynamicParallelFBA(StaticOptimizationModelBase):
    """A class representing a dynamic parallel FBA simulation.

    This class extends the TimeStepDynamicFBABase and provides functionality
    for performing dynamic FBA simulations on multiple models in parallel.

    Attributes:
    m_models (list[Model]): List of metabolic models to simulate in parallel.
    m_initial_bounds (dict[str, dict[str, tuple[float, float]]]): Original
        reaction bounds for each model, indexed by model ID and reaction ID.
    m_model_kinetics (dict[str, KineticsStruct]): Mapping of model ID to
        associated kinetic parameters for reactions.
    """

    m_models: dict[str, Model]
    m_initial_bounds: dict[str, dict[str, tuple[float, float]]] = {}
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

        self.m_models = {m.getId(): m.clone() for m in models}
        self.m_model_kinetics = kinetics

        for i, model in enumerate(self.m_models.values()):
            mid = model.getId()
            self.set_initial_concentrations(model, initial_concentrations)
            self.m_biomass_concentrations[mid] = [biomasses[i]]

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
        final_fluxes = {m: [] for m in self.m_models.keys()}
        used_times = [0]

        dt_save = dt
        dt_hat = -1

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
                    used_times,
                    run_condition,
                )

            self.update_reaction_bounds(kinetics_func)
            self.update_exchanges(dt)

            temp_fluxes = {m: {} for m in self.m_models.keys()}

            # Perform FBA foreach model
            for mid, model in self.m_models.items():
                solution = cbmpy.doFBA(model)
                if math.isnan(solution) or solution <= epsilon:
                    continue

                FBAsol = model.getSolutionVector(names=True)
                FBAsol = dict(zip(FBAsol[1], FBAsol[0]))
                temp_fluxes[mid] = FBAsol

            all_empty = all(
                not child_dict for child_dict in temp_fluxes.values()
            )

            if all_empty:
                return [
                    used_times,
                    self.m_metabolite_concentrations,
                    self.m_biomass_concentrations,
                    final_fluxes,
                ]

            self.update_concentrations(temp_fluxes, dt)
            species_id = self.check_solution_feasibility()

            if species_id:
                dt_hat = self.reset_dt(species_id, temp_fluxes)
                if dt_hat != -1 and dt_hat <= epsilon:
                    return [
                        used_times,
                        self.m_metabolite_concentrations,
                        self.m_biomass_concentrations,
                        final_fluxes,
                    ]
                else:
                    continue

            self.update_biomasses(temp_fluxes, dt)

            for mid, fbasol in temp_fluxes.items():
                final_fluxes[mid].append(fbasol)

            used_times.append(used_times[-1] + dt)
        return [
            used_times,
            self.m_metabolite_concentrations,
            self.m_biomass_concentrations,
            final_fluxes,
        ]

    def update_reaction_bounds(self, kinetics_func) -> None:
        """Update the reaction bounds for all models.

        Args:
            kinetics_func (function): Custom function to modify reaction bounds based on
                kinetic data.

        """

        if kinetics_func is not None:
            # TODO Give somethings to the kinetics_func, think about this
            kinetics_func()
            return

        for mid, model in self.m_models.items():
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
                            self.m_model_kinetics[mid],
                        )

                    else:
                        reaction.setUpperBound(
                            self.m_initial_bounds[mid][rid][1] * X_k_t
                        )

    def update_exchanges(self, dt):
        """
        Adjust exchange reactions based on the present metabolite concentrations

        Sets the lower bounds for exchange reactions to prevent metabolite import beyond
        the current concentration.
        """

        for model in self.m_models.values():
            for rid in model.getExchangeReactionIds():
                reaction: Reaction = model.getReaction(rid)
                sid = reaction.getSpeciesIds()[0]
                reaction.setLowerBound(
                    -self.m_metabolite_concentrations[sid][-1] * (1 / dt)
                )

    def update_concentrations(
        self, fluxes: dict[str, dict[str, float]], dt: float
    ) -> None:
        """
        Modify metabolite concentrations following an FBA step.

        Args:
            model (Model): Current metabolic model being processed.
            FBAsol (dict): Solution vector resulting from the FBA.
            step (int): Current simulation step.
            dt (float): Time step interval for the simulation.
        """
        for key in self.m_metabolite_concentrations.keys():
            self.m_metabolite_concentrations[key].append(
                self.m_metabolite_concentrations[key][-1]
            )
        for mid, fbasol in fluxes.items():
            # Check if dict is empty
            if fbasol:
                model = self.m_models[mid]
                for eid in model.getExchangeReactionIds():
                    reaction: Reaction = model.getReaction(eid)
                    sid = reaction.getSpeciesIds()[
                        0
                    ]  # exchanges should only have 1 reactant
                    self.m_metabolite_concentrations[sid][-1] += round(
                        (fbasol[eid] * dt), 5
                    )

    def update_biomasses(self, fluxes, dt) -> None:
        """
        Modify biomass concentrations after completing an FBA step.

        Args:
            model (Model): Model for which the biomass concentrations are updated.
            dt (float): Time step interval for the simulation.
            FBAsol (dict): Solution vector resulting from the FBA.
            step (int): Current simulation step.
        """
        for mid, fbasol in fluxes.items():
            # Check if dict is empty
            if fbasol:
                model = self.m_models[mid]
                Xt = (
                    self.m_biomass_concentrations[mid][-1]
                    + fbasol[model.getActiveObjectiveReactionIds()[0]] * dt
                )
            else:
                Xt = self.m_biomass_concentrations[mid][-1]

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

        # Remove all metabolites that were updated
        self.m_metabolite_concentrations = {
            key: lst[:-1]
            for key, lst in self.m_metabolite_concentrations.items()
        }

        total = 0
        for mid, model in self.m_models.items():
            for eid in model.getExchangeReactionIds():
                sid = model.getReaction(eid).getSpeciesIds()[0]
                if sid == species_id:
                    total += FBAsol[mid][eid] * -1
        return self.m_metabolite_concentrations[species_id][-1] / total
