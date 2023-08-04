import cbmpy
import math
from cbmpy.CBModel import Reaction
from ..Models import KineticsStruct, Transporters, CommunityModel
from .TimeStepDynamicModel import TimeStepDynamicModel


class DynamicFBABase(TimeStepDynamicModel):
    """Base class for SingleDynamicFBA and JointFBA
    Here the code for simulating both is implemented
    The two models only differ in the fact that there is
    a Community Biomass in Dynamic JointFBA
    """

    m_model: CommunityModel
    m_transporters: Transporters
    m_initial_bounds: dict[str, tuple[float, float]] = {}

    def __init__(
        self,
        model: CommunityModel,
        biomasses: list[float],
        initial_concentrations: dict[str, float] = {},
        kinetics: KineticsStruct = KineticsStruct({}),
    ):
        """Initialize the DynamicFBABase class.

        Args:
            model (CommunityModel): The community model to simulate.
            biomasses (list[float]): Initial biomass concentrations for each
                organism.
            initial_concentrations (dict[str, float], optional): Initial
                metabolite concentrations. Defaults to an empty dictionary.
            kinetics (Kinetics, optional): Kinetic parameters for reactions.
                Defaults to an empty Kinetics object.
        """

        self.m_transporters = Transporters(model)

        model_biomasses = model.get_model_biomass_ids()
        self.m_model = model

        # Set X_C to be exporter since it increases over time
        # self.set_community_biomass_reaction()
        # self.m_transporters.add_exporter("X_comm", ["X_c"])

        initial_biomasses = [[x] for x in biomasses]

        self.m_biomass_concentrations = dict(
            zip(model_biomasses.keys(), initial_biomasses)
        )

        self.m_kinetics = kinetics
        self.set_initial_concentrations(self.m_model, initial_concentrations)

        for rid in self.m_model.getReactionIds():
            reaction: Reaction = self.m_model.getReaction(rid)
            self.m_initial_bounds[rid] = [
                reaction.getLowerBound(),
                reaction.getUpperBound(),
            ]

    def simulate(
        self,
        dt: float,
        n: int = 10000,
        epsilon=0.01,
        kinetics_func=None,
        deviate=None,
    ):
        """Perform a dynamic joint FBA simulation.

        Returns:
            list: A list containing the simulation results -
                [used_time, metabolite_concentrations, biomass_concentrations, fluxes].
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
            self.update_exchanges()

            solution = cbmpy.doFBA(self.m_model, quiet=False)
            FBAsol = self.m_model.getSolutionVector(names=True)
            FBAsol = dict(zip(FBAsol[1], FBAsol[0]))

            if math.isnan(solution) or solution <= epsilon or dt < epsilon:
                break

            FBAsol = self.m_model.getSolutionVector(names=True)
            FBAsol = dict(zip(FBAsol[1], FBAsol[0]))

            used_time.append(used_time[-1] + dt)

            fluxes.append(FBAsol)

            self.update_concentrations(FBAsol, dt)

            species_id = self.check_solution_feasibility()

            if species_id != "":
                dt_hat = self.reset_dt(species_id, FBAsol)
                used_time = used_time[:-1]
                continue

            for _, rid in self.m_model.get_model_biomass_ids().items():
                mid = self.m_model.identify_model_from_reaction(rid)
                Xt = self.m_biomass_concentrations[mid][-1] + FBAsol[rid] * dt
                self.m_biomass_concentrations[mid].append(Xt)

        return [
            used_time,
            self.m_metabolite_concentrations,
            self.m_biomass_concentrations,
            fluxes,
        ]

    def update_exchanges(self):
        for rid in self.m_model.getExchangeReactionIds():
            reaction: Reaction = self.m_model.getReaction(rid)
            # Exchanges only have one species:
            sid = reaction.getSpeciesIds()[0]
            reaction.setLowerBound(-self.m_metabolite_concentrations[sid][-1])

    def update_concentrations(self, FBAsol, dt):
        """Update metabolite concentrations after an FBA step.

        Args:
            FBAsol (dict): The solution vector from the FBA.
            dt (float): The time step for the simulation.

        """

        # Update external metabolites

        # TODO just use exchange reactions?

        for key in self.m_metabolite_concentrations.keys():
            self.m_metabolite_concentrations[key].append(
                self.m_metabolite_concentrations[key][-1]
            )

        for rid, species_ids in self.m_transporters.get_exporters(True):
            for sid in species_ids:
                self.m_metabolite_concentrations[sid][-1] = round(
                    self.m_metabolite_concentrations[sid][-1]
                    + FBAsol[rid] * dt,
                    8,
                )

        for rid, species_ids in self.m_transporters.get_importers(True):
            for sid in species_ids:
                self.m_metabolite_concentrations[sid][-1] = round(
                    self.m_metabolite_concentrations[sid][-1]
                    - FBAsol[rid] * dt,
                    8,
                )

    def update_reaction_bounds(self, kinetics_func) -> None:
        """Update the reaction bounds for the simulation.

        Args:
            kinetics_func (function): A custom function to update reaction
                bounds based on kinetics. If None, standard bounds are used.

        """
        for rid in self.m_model.getReactionIds():
            # TODO discuss this
            if kinetics_func is not None:
                kinetics_func(
                    rid,
                    self.m_model,
                    self.m_biomass_concentrations,
                    self.m_biomass_concentrations,
                )
                continue
            reaction: Reaction = self.m_model.getReaction(rid)

            # Don't change the exchange reactions bounds
            if not reaction.is_exchange:
                mid_for_reaction = self.m_model.identify_model_from_reaction(
                    rid
                )
                # Organism specific biomass at last time point
                X_k_t = self.m_biomass_concentrations[mid_for_reaction][-1]
                reaction.setLowerBound(self.m_initial_bounds[rid][0] * X_k_t)

                if self.m_kinetics.exists(rid):
                    self.mm_kinetics(
                        reaction, X_k_t, self.m_transporters, self.m_kinetics
                    )

                else:
                    reaction.setUpperBound(
                        self.m_initial_bounds[rid][1] * X_k_t
                    )

    def reset_dt(self, species_id: str, FBAsol) -> float:
        """Recalculate the time step if a species becomes infeasible.

        Args:
            species_id (str): The ID of the species to reset the time step for.
            FBAsol (dict): The solution vector from the FBA.

        Returns:
            float: The updated time step.

        """

        # Remove last metabolite concentration
        self.m_metabolite_concentrations = {
            key: lst[:-1]
            for key, lst in self.m_metabolite_concentrations.items()
        }

        # Maybe some metabolite was exported/created, by an organism
        # Unfortunately to check what is create we would
        # have to know the new dt which we don'...
        # So we accept a small error.

        total = 0
        for (
            rid,
            species_ids,
        ) in self.m_transporters.get_importers(True):
            if species_id in species_ids:
                total += FBAsol[rid]

        return self.m_metabolite_concentrations[species_id][-1] / total
