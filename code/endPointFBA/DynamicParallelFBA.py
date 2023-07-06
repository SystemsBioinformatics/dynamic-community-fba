import cbmpy
import math
from cbmpy.CBModel import Model, Reaction, Species

from endPointFBA.DynamicFBABase import DynamicFBABase


class ParallelDynamicFBA(DynamicFBABase):
    m_initial_bounds = dict[str, dict[str, tuple[float, float]]]
    m_importers = dict[str, list[str]]
    m_exporters = dict[str, list[str]]

    def simulate(
        self,
        dt,
        models: list[Model],
        biomasses: list[float],
        initial_concentrations: dict[str, float] = {},
        kinetics={},
        epsilon=0.001,
        user_func=None,
    ):
        self.m_kinetics = kinetics
        for i in range(0, models):
            model = models[i]
            mid = model.getId()
            self.set_initial_concentrations(model, initial_concentrations)
            self.m_biomass_concentrations[mid] = [biomasses[i]]
            self.m_importers[mid] = self.get_importers(model)
            self.m_exporters[mid] = self.get_exporters(model)
            for rid in model.getReactionIds():
                reaction: Reaction = model.getReaction(rid)

                self.m_initial_bounds[mid][rid] = [
                    reaction.getLowerBound(),
                    reaction.getUpperBound(),
                ]

        step = 1
        used_time = [0]
        dt_hat = -1
        dt_save = dt

        while True:
            if dt_hat != -1:
                dt = dt_save - dt_hat
                dt_hat = -1
            else:
                dt = dt_save
            used_time.append(used_time[-1] + dt)

            for model in models:
                solution = cbmpy.doFBA(model)
                if math.isnan(solution) or solution == 0 or solution < epsilon:
                    return [
                        used_time[:-1],
                        self.m_metabolite_concentrations,
                        self.m_biomass_concentrations,
                    ]
                FBAsol = model.getSolutionVector(names=True)
                FBAsol = dict(zip(FBAsol[1], FBAsol[0]))
                self.update_concentrations(FBAsol, dt, step)

                mid = model.getId()
                Xt = (
                    self.m_biomass_concentrations[mid][step - 1]
                    + FBAsol[model.getActiveObjectiveReactionIds()[0]] * dt
                )
                self.m_biomass_concentrations[mid].append(Xt)

                self.update_bounds(model, user_func)

    def update_bounds(self, model: Model, user_func):
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
                if rid in self.m_importers[mid].keys():
                    self.update_importer_bounds(reaction, X_k_t)
                else:
                    reaction.setUpperBound(
                        self.m_initial_bounds[rid][1] * X_k_t
                    )

    def update_importer_bounds(self):
        return super().update_importer_bounds()

    def update_concentrations(self):
        return super().update_concentrations()
