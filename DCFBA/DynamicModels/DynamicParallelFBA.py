import cbmpy
import math
from cbmpy.CBModel import Model, Reaction
from .TimeStepDynamicFBABase import TimeStepDynamicFBABase
from ..Models.Kinetics import Kinetics
from ..Models.Transporters import Transporters


class DynamicParallelFBA(TimeStepDynamicFBABase):
    """A class representing a dynamic parallel FBA simulation.

    This class extends the TimeStepDynamicFBABase and provides functionality
    for performing dynamic FBA simulations on multiple models in parallel.

    Attributes:
        m_models (list[Model]): List of models to run in parallel.
        m_initial_bounds (dict[str, dict[str, tuple[float, float]]]):
            Initial reaction bounds for each model.
        m_model_transporters (dict[str, Transporters]): Transporters
             associated with each model.
        m_model_kinetics (dict[str, Kinetics]): Kinetic parameters for
            reactions in each model.

    """

    m_models: list[Model]
    m_initial_bounds: dict[str, dict[str, tuple[float, float]]] = {}
    m_model_transporters: dict[str, Transporters] = {}
    m_model_kinetics: dict[str, Kinetics] = {}

    def __init__(
        self,
        models: list[Model],
        biomasses: list[float],
        initial_concentrations: dict[str, float] = {},
        kinetics: dict[str, Kinetics] = {},
    ) -> None:
        """Initialize the DynamicParallelFBA class.

        Args:
            models (list[Model]): List of models you want to run in parallel.
            biomasses (list[float]): List of initial biomass concentrations.
                Use the same order as the list of models.
            initial_concentrations (dict[str, float], optional): Initial
                metabolite concentrations. Defaults to an empty dictionary.
            kinetics (dict[str, Kinetics], optional): A dictionary with model
                IDs and their respective kinetics objects.
                Defaults to an empty dictionary.

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
        epsilon: float = 0.001,
        kinetics_func=None,
        deviate=None,
    ) -> list[list[float], dict[str, list[float]], dict[str, list[float]]]:
        """Perform a dynamic parallel FBA simulation.

        Args:
            dt float: The time step for the simulation.
            epsilon (float, optional): The convergence threshold.
                Defaults to 0.001.
            kinetics_func (_type_, optional): A custom function to update
                reaction bounds based on kinetics. Defaults to None.
            deviate (function, optional): A custom function to introduce
                perturbations during the simulation.
                Defaults to None.
            deviation_time (float, optional): The time at which to introduce
                the perturbation.
                Defaults to 0.

        Returns:
            list[list[float], dict[str, list[float]], dict[str, list[float]]]: A list containing the simulation results -
                [used_time, metabolite_concentrations, biomass_concentrations].

        """

        step = 1
        used_time = [0]
        dt_hat = math.nan
        dt_save = dt
        run_condition = 0

        while True:
            if not math.isnan(dt_hat):
                dt = dt_hat
                dt_hat = math.nan
            else:
                dt = dt_save

            if deviate is not None:
                run_condition += deviate(
                    self,
                    used_time,
                    run_condition,
                )

            self.update_reaction_bounds(kinetics_func)

            self.constrain_exchanges()

            used_time.append(used_time[-1] + dt)

            for model in self.m_models:
                solution = cbmpy.doFBA(model)
                FBAsol = model.getSolutionVector(names=True)
                FBAsol = dict(zip(FBAsol[1], FBAsol[0]))
                if math.isnan(solution) or solution < epsilon or dt < epsilon:
                    print(dt)
                    return [
                        used_time[:-1],
                        self.m_metabolite_concentrations,
                        self.m_biomass_concentrations,
                    ]

                self.update_concentrations(model, FBAsol, step, dt)
                self.update_biomasses(model, dt, FBAsol, step)

            species_id = self.check_solution_feasibility()

            if species_id != "":
                dt_hat = self.reset_dt(species_id, FBAsol)
                step -= 1
                used_time = used_time[:-1]

            step += 1

    def constrain_exchanges(self):
        """Constrain exchange reactions based on current metabolite
         concentrations and biomass.

        This method sets lower bounds for exchange reactions to prevent
        metabolite import beyond the current concentration and biomass
        of the associated organism.

        """
        for model in self.m_models:
            for rid in model.getExchangeReactionIds():
                reaction: Reaction = model.getReaction(rid)
                sid = reaction.getSpeciesIds()[0]
                reaction.setLowerBound(
                    (
                        (-1 * self.m_metabolite_concentrations[sid][-1])
                        * self.m_biomass_concentrations[model.getId()][-1]
                    )
                )

    def update_reaction_bounds(self, kinetics_func) -> None:
        """Update the reaction bounds for all models.

        Args:
            kinetics_func (function): A custom function to update reaction
                bounds based on kinetics.

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
                    # If the reaction is an importer we need to check
                    # if there is any substrate they can import
                    if self.m_model_transporters[mid].is_importer(rid):
                        self.update_importer_bounds(mid, reaction, X_k_t)

                    elif (mid in self.m_model_kinetics.keys()) and (
                        self.m_model_kinetics[mid].Exists(rid)
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

    def update_importer_bounds(
        self, mid: str, reaction: Reaction, X: float
    ) -> None:
        """Update bounds for importer reactions.

        Args:
            mid (str): The ID of the model for which to update importer bounds.
            reaction (Reaction): The reaction to update bounds for.
            X (float): Biomass concentration for the organism associated with
                the reaction.

        """
        sids = self.m_model_transporters[mid].get_importers_species(
            reaction.getId()
        )
        importers_reagent_concentrations = [
            self.m_metabolite_concentrations[id][-1] for id in sids
        ]

        if all(x > 0 for x in importers_reagent_concentrations):
            if mid in self.m_model_kinetics.keys() and self.m_model_kinetics[
                mid
            ].Exists(reaction.getId()):
                self.mm_kinetics(
                    reaction,
                    X,
                    self.m_model_transporters[mid],
                    self.m_model_kinetics[mid],
                )
            else:
                reaction.setUpperBound(
                    self.m_initial_bounds[mid][reaction.getId()][1] * X
                )

        else:
            reaction.setUpperBound(0)

    def update_concentrations(self, model: Model, FBAsol, step, dt) -> None:
        """Update metabolite concentrations after an FBA step.

        Args:
            model (Model): The current model
            dt (float): The time step for the simulation.
            FBAsol (dict): The solution vector from the FBA.
            step (int): The current simulation step.

        """

        mid = model.getId()

        if len(next(iter(self.m_metabolite_concentrations.values()))) == step:
            for key in self.m_metabolite_concentrations.keys():
                self.m_metabolite_concentrations[key].append(
                    self.m_metabolite_concentrations[key][-1]
                )

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
        """Update biomass concentrations after an FBA step.

        Args:
            model (Model): The model for which to update biomass
                concentrations.
            dt (float): The time step for the simulation.
            FBAsol (dict): The solution vector from the FBA.
            step (int): The current simulation step.

        """
        mid = model.getId()
        Xt = (
            self.m_biomass_concentrations[mid][step - 1]
            + FBAsol[model.getActiveObjectiveReactionIds()[0]] * dt
        )
        self.m_biomass_concentrations[mid].append(Xt)

    def reset_dt(self, species_id, FBAsol) -> float:
        """Recalculate a smaller time step if the result was infeasible.

        Args:
            species_id (str): The ID of the species to reset the time step for.
            FBAsol (dict): The solution vector from the FBA.

        Returns:
            float: The updated time step.

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
                    total += FBAsol[rid]

        return self.m_metabolite_concentrations[species_id][-1] / total
