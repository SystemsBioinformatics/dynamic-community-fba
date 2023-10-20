from cbmpy.CBModel import Model, Reaction
from .StaticOptimizationModel import StaticOptimizationModelBase
from ..Models.Kinetics import KineticsStruct


class DynamicParallelFBA(StaticOptimizationModelBase):
    """A class representing a dynamic parallel FBA simulation.

    This class extends the TimeStepDynamicFBABase and provides functionality
    for performing dynamic FBA simulations on multiple models in parallel.
    """

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
        super().__init__()

        self._models = {m.getId(): m.clone() for m in models}
        self._kinetics = kinetics

        for i, model in enumerate(self._models.values()):
            mid = model.getId()
            self.set_initial_concentrations(model, initial_concentrations)
            self._biomasses[mid] = [biomasses[i]]

            for rid in model.getReactionIds():
                reaction: Reaction = model.getReaction(rid)
                if mid not in self.initial_bounds.keys():
                    self._initial_bounds[mid] = {
                        rid: [
                            reaction.getLowerBound(),
                            reaction.getUpperBound(),
                        ]
                    }
                else:
                    self._initial_bounds[mid][rid] = [
                        reaction.getLowerBound(),
                        reaction.getUpperBound(),
                    ]

    @property
    def models(self) -> dict[str, Model]:
        return self._models

    def get_flux_values(self, mid: str, rid: str) -> list[float]:
        return list(map(lambda d: d[rid], self.fluxes[mid]))

    def get_fluxes_values(self, mid, rids) -> dict[str, float]:
        fluxes = {}
        for rid in rids:
            fluxes[rid] = self.get_flux_values(mid, rid)
        return fluxes

    def get_specific_flux_values(self, mid: str, rid: str) -> list[float]:
        for model in self.models.values():
            if rid in model.getExchangeReactionIds():
                print("Exchange reactions do not have a specific flux")
                return None

        ls = self.get_flux_values(mid, rid)

        return [ls[t] / v for t, v in enumerate(self.biomasses[mid][:-1])]

    def get_community_growth_rate(self):
        community_flux = [0] * (len(self.times) - 1)
        total_mass = [0] * (len(self.times) - 1)
        for mid, model in self.models.items():
            rid = model.getActiveObjectiveReactionIds()[0]
            values = self.get_flux_values(mid, rid)
            for i, value in enumerate(values):
                community_flux[i] += value
                total_mass[i] += self.biomasses[mid][i]

        return [flux / total_mass[i] for i, flux in enumerate(community_flux)]

    def get_relative_abundance(self) -> dict[str, list[float]]:
        ans = {}
        total_mass = [0] * (len(self.times))
        for mid in self.models.keys():
            for i in range(0, len(self.times)):
                total_mass[i] += self.biomasses[mid][i]

        for mid in self.models.keys():
            ans[mid] = [
                v / total_mass[i] for i, v in enumerate(self.biomasses[mid])
            ]

        return total_mass

    def simulate(
        self,
        dt: float,
        n: int = 10000,
        epsilon: float = 0.01,
        kinetics_func=None,
        deviate=None,
    ) -> None:
        final_fluxes = {m: [] for m in self.models.keys()}
        used_times = [0]

        dt_save = dt
        dt_hat = -1

        run_condition = 0

        for _ in range(1, n):
            if dt_hat != -1:
                dt = dt_save

            self.update_reaction_bounds(kinetics_func)
            self.update_exchanges(dt)
            temp_fluxes = {m: {} for m in self.models.keys()}
            if deviate is not None:
                run_condition += deviate(
                    self,
                    used_times,
                    run_condition,
                )

            # Perform FBA foreach model
            step = 0
            for mid, model in self.models.items():
                solution = cbmpy.doFBA(model, quiet=True)

                if math.isnan(solution):
                    print(f"model: {mid} had an infeasible solution")
                    self._times = used_times
                    self.fluxes = final_fluxes
                    return

                if solution <= epsilon:
                    step += 1

                FBAsol = model.getSolutionVector(names=True)
                FBAsol = dict(zip(FBAsol[1], FBAsol[0]))
                temp_fluxes[mid] = FBAsol

            if step == len(self.models.keys()):
                print("All models had solution value of 0")
                self._times = used_times
                self.fluxes = final_fluxes
                return

            self.update_concentrations(temp_fluxes, dt)
            species_id = self.check_solution_feasibility()

            if species_id:
                dt_hat = self.reset_dt(species_id, temp_fluxes)
                dt = dt_hat
                if dt_hat <= epsilon:
                    print("dt was smaller than provided value of epsilon")
                    break

                self.update_concentrations(temp_fluxes, dt)
                species_id = self.check_solution_feasibility()

            self.update_biomasses(temp_fluxes, dt)

            for mid, fbasol in temp_fluxes.items():
                final_fluxes[mid].append(fbasol)

            used_times.append(used_times[-1] + dt)

        print(f"All {n} simulation were performed")
        self._times = used_times
        self._fluxes = final_fluxes

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

        for mid, model in self.models.items():
            for rid in model.getReactionIds():
                reaction: Reaction = model.getReaction(rid)
                if not reaction.is_exchange:
                    # Organism specific biomass
                    X_k_t = self.biomasses[mid][-1]
                    lb, ub = self.initial_bounds[mid][rid]
                    reaction.setLowerBound(lb * X_k_t)

                    if (mid in self.kinetics.keys()) and (
                        self.kinetics[mid].exists(rid)
                    ):
                        self.mm_kinetics(
                            reaction,
                            X_k_t,
                            self.kinetics[mid],
                        )

                    else:
                        reaction.setUpperBound(ub * X_k_t)

    def update_exchanges(self, dt):
        """
        Adjust exchange reactions based on the present metabolite concentrations

        Sets the lower bounds for exchange reactions to prevent metabolite import beyond
        the current concentration.
        """

        for model in self.models.values():
            for rid in model.getExchangeReactionIds():
                reaction: Reaction = model.getReaction(rid)
                sid = reaction.getSpeciesIds()[0]
                reaction.setLowerBound(-self.metabolites[sid][-1] * (1 / dt))

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
        for key in self.metabolites.keys():
            self.metabolites[key].append(self.metabolites[key][-1])
        for mid, fbasol in fluxes.items():
            # Check if dict is empty #TODO can be removed?
            if fbasol:
                model = self.models[mid]
                for eid in model.getExchangeReactionIds():
                    reaction: Reaction = model.getReaction(eid)
                    sid = reaction.getSpeciesIds()[
                        0
                    ]  # exchanges should only have 1 reactant
                    self.metabolites[sid][-1] += round((fbasol[eid] * dt), 5)

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
                model = self.models[mid]

                Xt = (
                    self.biomasses[mid][-1]
                    + fbasol[model.getActiveObjectiveReactionIds()[0]] * dt
                )

            else:
                Xt = self.biomasses[mid][-1]

            self.biomasses[mid].append(Xt)

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
        self._metabolites = {
            key: lst[:-1] for key, lst in self.metabolites.items()
        }

        total = 0
        for mid, model in self.models.items():
            for eid in model.getExchangeReactionIds():
                sid = model.getReaction(eid).getSpeciesIds()[0]
                if sid == species_id:
                    total += FBAsol[mid][eid] * -1

        return self.metabolites[species_id][-1] / total
