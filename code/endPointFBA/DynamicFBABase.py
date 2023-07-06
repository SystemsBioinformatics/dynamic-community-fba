from cbmpy.CBModel import Model, Reaction, Species


class DynamicFBABase:
    m_biomass_concentrations: dict[str, list[float]]
    m_metabolite_concentrations: dict[str, list[float]] = {}
    m_kinetics: dict[str, tuple[float, float]]

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

    def update_bounds(self):
        pass

    def update_importer_bounds(self):
        pass

    def update_concentrations(self):
        pass

    def get_exporters(self, model) -> dict[str, list[str]]:
        return self.get_transporters(model)[0]

    def get_importers(self, model) -> dict[str, list[str]]:
        return self.get_transporters(model)[1]

    def get_transporters(self, model: Model) -> list[dict[str, list[str]]]:
        """It is convenient to access the importers and exporters quickly
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
        for rid in model.getReactionIds():
            reaction: Reaction = model.getReaction(rid)
            if not reaction.is_exchange:
                for reagent in reaction.reagents:
                    sid: str = reagent.getSpecies()
                    species: Species = model.getSpecies(sid)
                    if species.getCompartmentId() == "e":
                        if reagent.coefficient == -1:
                            # Reaction imports external metabolites
                            importers[rid] = importers.get(rid, []) + [sid]
                        elif reagent.coefficient == 1:
                            exporters[rid] = exporters.get(rid, []) + [sid]

        return [exporters, importers]
