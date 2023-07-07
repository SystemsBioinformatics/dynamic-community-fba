from cbmpy.CBModel import Model, Reaction, Species
from endPointFBA.Models.Kinetics import Kinetics
from endPointFBA.Models.Transporters import Transporters


class DynamicFBABase:
    m_kinetics: Kinetics
    m_biomass_concentrations: dict[str, list[float]] = {}
    m_metabolite_concentrations: dict[str, list[float]] = {}

    def set_initial_concentrations(
        self, model: Model, initial_concentrations: dict[str, float]
    ):
        # Use the exchange reaction lower bounds as initial concentrations
        # If people want to change it they have to change this use initial
        # conditions or change the lower bounds
        # This was implemented in such that you dont need to manually
        # add all external metabolites in the initial conditions

        for exchange in model.getExchangeReactionIds():
            reaction: Reaction = model.getReaction(exchange)
            species_id: str = reaction.reagents[0].getSpecies()
            if species_id not in self.m_metabolite_concentrations.keys():
                if species_id in initial_concentrations.keys():
                    self.m_metabolite_concentrations[species_id] = [
                        initial_concentrations[species_id]
                    ]
                else:
                    self.m_metabolite_concentrations[species_id] = [
                        -reaction.getLowerBound()
                    ]

    def simulate(self, dt, epsilon=0.001) -> None:
        pass

    def update_reaction_bounds(self):
        pass

    def update_importer_bounds(self):
        pass

    def update_concentrations(self):
        pass

    def check_solution_feasibility(self) -> str:
        low = 1e10
        name = ""
        for key, value in self.m_metabolite_concentrations.items():
            if value[-1] < 0 and value[-1] < low:
                low = value[-1]
                name = key

        return name

    def importers_species_concentration(
        self, rid: str, transporters: Transporters
    ):
        sids = transporters.get_importers_species(rid)
        return [self.m_metabolite_concentrations[id][-1] for id in sids]

    def update_kinetics(
        self, reaction: Reaction, X: float, transporters: Transporters
    ):
        rid: str = reaction.getId()
        km, vmax = self.m_kinetics.get_kinetics(
            rid,
        )
        S = min(self.importers_species_concentration(rid, transporters))
        reaction.setUpperBound(vmax * (S / (km + S)) * X)
