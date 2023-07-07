import cbmpy
import math
from cbmpy.CBModel import Model, Reaction, Species

from endPointFBA.DynamicFBABase import DynamicFBABase
from endPointFBA.Models.Kinetics import Kinetics
from endPointFBA.Models.Transporters import Transporters


class DynamicParallelFBA(DynamicFBABase):
    m_models: list[Model]
    m_initial_bounds: dict[str, dict[str, tuple[float, float]]] = {}

    # Here the dict represents a model_id: Transporter
    m_model_transporters: dict[str, Transporters] = {}

    def __init__(
        self,
        models: list[Model],
        biomasses: dict[str, float],
        initial_concentrations: dict[str, float],
        kinetics: dict[str, Kinetics] = Kinetics({}),
    ) -> None:
        self.m_models = models
        self.m_kinetics = kinetics

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
        dt,
        epsilon=0.001,
        user_func=None,
    ):
        step = 1
        used_time = [0]
        dt_hat = -1
        dt_save = dt

        self.update_reaction_bounds(user_func)

        while True:
            if dt_hat != -1:
                dt = dt_hat
                dt_hat = -1
            else:
                dt = dt_save

            used_time.append(used_time[-1] + dt)
            for model in self.m_models:
                mid = model.getId()

                solution = cbmpy.doFBA(model)
                FBAsol = model.getSolutionVector(names=True)
                FBAsol = dict(zip(FBAsol[1], FBAsol[0]))

                if math.isnan(solution) or solution == 0 or solution < epsilon:
                    return [
                        used_time[:-1],
                        self.m_metabolite_concentrations,
                        self.m_biomass_concentrations,
                    ]

                self.update_concentrations(mid, dt, FBAsol, step)

                Xt = (
                    self.m_biomass_concentrations[mid][step - 1]
                    + FBAsol[model.getActiveObjectiveReactionIds()[0]] * dt
                )
                self.m_biomass_concentrations[mid].append(Xt)

            species_id = self.check_solution_feasibility()

            if species_id != "":
                # TODO this code needs some propper clean up
                # Fix with DynamicJointFBA

                # remove all last concentrations
                self.m_metabolite_concentrations = {
                    key: lst[:-1]
                    for key, lst in self.m_metabolite_concentrations.items()
                }

                self.m_biomass_concentrations = {
                    key: lst[:-1]
                    for key, lst in self.m_biomass_concentrations.items()
                }

                total = 0
                for model in self.m_models:
                    mid = model.getId()
                    for rid, species_ids in self.m_model_transporters[
                        mid
                    ].get_importers(True):
                        if species_id in species_ids:
                            total += FBAsol[rid]

                dt_hat = (
                    self.m_metabolite_concentrations[species_id][-1] / total
                )

                step -= 1
                used_time = used_time[:-1]

            # After all models updated the external metabolites
            # update all bounds

            self.update_reaction_bounds(user_func)

            step += 1

    def update_reaction_bounds(self, user_func):
        models = self.m_models
        if user_func is not None:
            # TODO Give somethings to the user_func, think about this
            user_func()
            return

        for model in models:
            mid = model.getId()
            for rid in model.getReactionIds():
                reaction: Reaction = model.getReaction(rid)
                # Don't change the exchange reactions bounds
                if not reaction.is_exchange:
                    # Organism specific biomass
                    X_k_t = self.m_biomass_concentrations[mid][-1]
                    reaction.setLowerBound(
                        self.m_initial_bounds[mid][rid][0] * X_k_t
                    )
                    # If the reaction is an importer we need to check
                    # if there is substrate they can import
                    if self.m_model_transporters[mid].is_importer(rid):
                        self.update_importer_bounds(mid, reaction, X_k_t)

                    elif (
                        self.m_kinetics is not None
                        and self.m_kinetics.Exists(rid)
                    ):
                        self.update_kinetics(
                            reaction, X_k_t, self.m_model_transporters[mid]
                        )

                    else:
                        reaction.setUpperBound(
                            self.m_initial_bounds[mid][rid][1] * X_k_t
                        )

    def update_importer_bounds(self, mid: str, reaction: Reaction, X: float):
        sids = self.m_model_transporters[mid].get_importers_species(
            reaction.getId()
        )
        importers_reagent_concentrations = [
            self.m_metabolite_concentrations[id][-1] for id in sids
        ]

        if all(x > 0 for x in importers_reagent_concentrations):
            if self.m_kinetics.Exists(reaction.getId()):
                self.update_kinetics(reaction, X)
            else:
                reaction.setUpperBound(
                    self.m_initial_bounds[mid][reaction.getId()][1] * X
                )

        else:
            reaction.setUpperBound(0)

    def update_concentrations(self, mid, dt, FBAsol, step):
        # first check what was exported
        if len(next(iter(self.m_metabolite_concentrations.values()))) == step:
            for key in self.m_metabolite_concentrations.keys():
                self.m_metabolite_concentrations[key].append(
                    self.m_metabolite_concentrations[key][-1]
                )

        for rid, species_ids in self.m_model_transporters[mid].get_exporters(
            True
        ):
            for sid in species_ids:
                self.m_metabolite_concentrations[sid][-1] += FBAsol[rid] * dt

        for rid, species_ids in self.m_model_transporters[mid].get_importers(
            True
        ):
            for sid in species_ids:
                self.m_metabolite_concentrations[sid][-1] -= FBAsol[rid] * dt