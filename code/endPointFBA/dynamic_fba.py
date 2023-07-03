from typing import Tuple
import cbmpy
import math
from cbmpy.CBModel import Model, Reaction, Species, Reagent
from endPointFBA.CommunityModel import CommunityModel
from .KineticModel import KineticModel


class DynamicJointFBA:
    # Maybe save dt here as well??

    m_kinetic_model: KineticModel
    m_model: CommunityModel
    m_biomass_reactions: list[str]
    m_importers: dict[str, list[str]]
    m_exporters: dict[str, list[str]]
    m_biomass_concentrations: dict[str, list[float]] = {}
    m_metabolite_concentrations: dict[str, list[float]] = {}
    m_initial_bounds: dict[str, Tuple[float, float]] = {}

    # Usec to set everything in place for the simulation
    def __init__(
        self,
        km: KineticModel,
        initial_biomasses: list[float],
        initial_conditions: dict[str, float] = {},
    ):
        self.m_kinetic_model = km
        self.m_model: CommunityModel = km.get_model().clone()
        model_biomasses = self.m_model.get_model_biomass_ids()

        self.m_biomass_reactions = list(model_biomasses.values())
        self.m_biomass_concentrations = dict(
            zip(model_biomasses.keys(), [initial_biomasses])
        )

        # Use the exchange reaction lower bounds as initial concentrations
        # If people want to change it they have to change this use inital
        # conditions or change the lower bounds
        # This was implemented in such that you dont need to manually
        # add all external metabolites in the initial conditions
        for exchange in self.m_model.getExchangeReactionIds():
            reaction: Reaction = self.m_model.getReaction(exchange)
            species_id: str = reaction.reagents[0].getSpecies()
            if species_id in initial_conditions.keys():
                self.m_metabolite_concentrations[species_id] = [
                    initial_conditions[species_id]
                ]
            else:
                self.m_metabolite_concentrations[species_id] = [
                    -reaction.getLowerBound()
                ]
        self.m_exporters, self.m_importers = self.get_transporters()
        for rid in self.m_model.getReactionIds():
            reaction: Reaction = self.m_model.getReaction(rid)
            self.m_initial_bounds[rid] = [
                reaction.getLowerBound(),
                reaction.getUpperBound(),
            ]

    def simulate(
        self,
        dt: float,
        user_func=None,
    ):
        # First update all upper and lower bounds to be rates
        # (Corrected for biomass),
        self.update_bounds(user_func)

        used_time = [0]
        step = 1
        while True:
            solution = cbmpy.doFBA(self.m_model)

            if math.isnan(solution) or solution == 0 or solution < 0.001:
                break

            used_time.append(used_time[-1] + dt)

            FBAsol = self.m_model.getSolutionVector(names=True)
            FBAsol = dict(zip(FBAsol[1], FBAsol[0]))

            species_id = self.update_concentrations(FBAsol, dt)

            if species_id != "":
                total = 0
                for rid, species_ids in self.m_importers.items():
                    if species_id in species_ids:
                        total += FBAsol[rid]
                dt_hat = (
                    self.m_metabolite_concentrations[species_id][-2] / total
                )

                self.m_metabolite_concentrations = {
                    key: lst[:-1]
                    for key, lst in self.m_metabolite_concentrations.items()
                }

                # than recalculate
                self.update_concentrations(
                    FBAsol,
                    dt_hat,
                )

                used_time[-1] = used_time[-2] + dt_hat

            # update biomass
            for rid in self.m_biomass_reactions:
                mid = self.m_model.identify_model_from_reaction(rid)
                Xt = self.m_biomass_concentrations[mid][-1] + FBAsol[rid] * dt
                self.m_biomass_concentrations[mid].append(Xt)

            # update lower and upper bounds...
            self.update_bounds(user_func)

            step += 1
        return [
            used_time,
            self.m_metabolite_concentrations,
            self.m_biomass_concentrations,
        ]

    def update_concentrations(self, FBAsol, dt):
        # Update external metabolites
        lowest_concentration = 1e10
        lowest_species_id = ""
        # first check what was exported
        for key in self.m_metabolite_concentrations.keys():
            self.m_metabolite_concentrations[key].append(
                self.m_metabolite_concentrations[key][-1]
            )

        for rid, species_ids in self.m_exporters.items():
            for sid in species_ids:
                self.m_metabolite_concentrations[sid][-1] += FBAsol[rid] * dt
        for rid, species_ids in self.m_importers.items():
            for sid in species_ids:
                self.m_metabolite_concentrations[sid][-1] -= FBAsol[rid] * dt
                new_concentration = self.m_metabolite_concentrations[sid][-1]
                if (
                    new_concentration < 0
                    and new_concentration < lowest_concentration
                ):
                    lowest_concentration = new_concentration
                    lowest_species_id = sid

        return lowest_species_id

    def update_bounds(self, user_func) -> None:
        for rid in self.m_model.getReactionIds():
            reaction: Reaction = self.m_model.getReaction(rid)
            # Don't change the exchange reactions bounds
            if not reaction.is_exchange:
                mid_for_reaction = self.m_model.identify_model_from_reaction(
                    rid
                )
                # Organism specific biomass
                X_k_t = self.m_biomass_concentrations[mid_for_reaction][-1]
                reaction.setLowerBound(self.m_initial_bounds[rid][0] * X_k_t)

                # If the reaction is an importer we need to check
                # if there is substrate they can import
                if rid in self.m_importers.keys():
                    self.update_importer_bounds(reaction, X_k_t)
                else:
                    reaction.setUpperBound(
                        self.m_initial_bounds[rid][1] * X_k_t
                    )

    def update_importer_bounds(self, reaction: Reaction, X: float):
        km = self.m_kinetic_model.get_reaction_km(reaction.getId())
        vmax = self.m_kinetic_model.get_reaction_vmax(reaction.getId())
        sids = self.m_importers[reaction.getId()]
        importers_reagent_concentrations = [
            self.m_metabolite_concentrations[id][-1] for id in sids
        ]

        if all(x > 0 for x in importers_reagent_concentrations):
            reaction.setUpperBound(
                self.m_initial_bounds[reaction.getId()][1] * X
            )
            # TODO fix tis
            if km is not None:
                S = min(importers_reagent_concentrations)
                reaction.setUpperBound(vmax * (S / (km + S)) * X)

        else:
            reaction.setUpperBound(0)

    def get_transporters(self):
        """The dynamic FBA needs to update all external metabolites
            therefore we need to know which reactions uptake metabolites
            from the external space and which are being excreted
            we define importers and exporters
            Both are a {rid : [species]}
        Args:
            model (Model): _description_

        Returns:
            _type_: _description_
        """

        importers: dict[str, list[str]] = {}
        exporters: dict[str, list[str]] = {}
        for rid in self.m_model.getReactionIds():
            reaction: Reaction = self.m_model.getReaction(rid)
            if not reaction.is_exchange:
                for reagent in reaction.reagents:
                    sid: str = reagent.getSpecies()
                    species: Species = self.m_model.getSpecies(sid)
                    if species.getCompartmentId() == "e":
                        if reagent.coefficient == -1:
                            # Reaction imports external metabolites
                            importers[rid] = importers.get(rid, []) + [sid]
                        elif reagent.coefficient == 1:
                            exporters[rid] = exporters.get(rid, []) + [sid]

        return [exporters, importers]
