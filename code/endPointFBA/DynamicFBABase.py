from cbmpy.CBModel import Model, Reaction, Species


class DynamicFBABase:
    def simulate(self, dt, epsilon=0.001) -> None:
        pass

    def update_bounds(self):
        pass

    def update_importer_bounds(self):
        pass

    def get_exporters(self, model) -> list[str]:
        return self.get_transporters(model)[0]

    def get_importers(self, model) -> list[str]:
        return self.get_transporters(model)[1]

    def get_transporters(self, model: Model) -> list[str]:
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
