from cbmpy.CBModel import Model, Species, Reaction
import endPointFBA.helpers.build_community_matrix as cm
from endPointFBA.Exceptions.NotInCombinedModel import NotInCombinedModel


class CommunityModel(Model):
    m_identifiers: list[str]
    m_single_model_identifiers: list[str]
    m_single_model_biomass_reaction_ids: list[str] = []

    def __init__(
        self,
        models: list[Model],
        biomass_reaction_ids: list[str],
        ids: list[str] = [],
        combined_model_id: str = "combined_model",
    ) -> None:
        super().__init__(combined_model_id)

        self.createCompartment("e", "extracellular space")
        duplicate_species = cm.create_duplicate_species_dict(models)

        self.m_single_model_identifiers = [model.id for model in models]

        if len(ids) > 0 and len(ids) < len(models):
            raise Exception("Too few ids were provided")
        if len(ids) == 0:
            self.m_identifiers = self.m_single_model_identifiers.copy()
        else:
            self.m_identifiers = ids

        # Save biomass reaction id of old model, make sure
        # the first reaction is the biomass reaction
        for i in range(0, len(biomass_reaction_ids)):
            bm = cm.create_new_id(
                biomass_reaction_ids[i], self.m_identifiers[i]
            )

            self.m_single_model_biomass_reaction_ids.append(bm)

        for i in range(0, len(models)):
            model = models[i]
            new_id = self.m_identifiers[i]

            cm.merge_compartments(model, self, new_id)
            cm.merge_species(duplicate_species, model, new_id)
            cm.merge_reactions(model, self, new_id)

    def __str__(self) -> str:
        return (
            f"Model: {self.getId()} was build from "
            f"{[id for id in self.m_single_model_identifiers]}"
        )

    def add_model_to_community(
        self, model: Model, biomass_reaction: str, new_id: str = None
    ):
        """Adds a model to the CommunityModel

        Args:
            model (Model): The model that needs to be added
            biomass_reaction (str): The reaction id of the biomass reaction of
            the new model
            new_id (str, optional): The user set identifier. Defaults to None.
            If set to None the model.id will be used
        """
        if new_id is None:
            new_id = model.getId()

        duplicate_species = cm.create_duplicate_species_dict([self, model])

        cm.merge_compartments(model, self, new_id)
        cm.merge_species(duplicate_species, model, new_id)
        cm.merge_reactions(model, self, new_id)

        self.m_identifiers.append(new_id)
        self.m_single_model_identifiers.append(model.id)
        self.m_single_model_biomass_reaction_ids.append(biomass_reaction)

    def get_model_specific_reactions(self, mid: str) -> list[Reaction]:
        """Returns a list of reaction ids

        Args:
            mid (str): _description_

        Raises:
            NotInCombinedModel: _description_

        Returns:
            list[Reaction]: _description_
        """
        if mid not in self.m_identifiers:
            raise NotInCombinedModel(
                "The model id provided was not found in the combined model"
            )
        ans = []
        for reaction_id in self.m_model.getReactionIds():
            if mid in reaction_id:
                ans.append(reaction_id)
        return ans

    def get_model_specific_species(self, mid: str) -> list[Species]:
        if mid not in self.m_identifiers:
            raise NotInCombinedModel(
                "The model id provided was not found in the combined model"
            )
        ans = []
        for species_id in self.m_model.getSpeciesIds():
            species: Species = self.m_model.getSpecies(species_id)
            if mid in species.getCompartmentId():
                ans.append(species_id)
        return ans

    def get_reaction_bigg_ids(self, mid="") -> list[str]:
        """Get the reaction bigg ids of all reactions

        Args:
            mid (str, optional): If a model id is provided only reactions from
            the specific model are returned
            Defaults to "".

        Raises:
            NotInCombinedModel: the id provided was not in the combined model

        Returns:
            list[str]: list containing the bigg ids
        """
        if mid != "":
            if mid not in self.m_identifiers:
                raise NotInCombinedModel(
                    "The model id provided was not found in the combined model"
                )
            reaction_ids = self.get_model_specific_reactions(mid)
            reaction_ids = [rid.replace(f"_{mid}", "") for rid in reaction_ids]
        else:
            reaction_ids = self.m_model.getReactionIds()
            for appended_id in self.m_identifiers:
                reaction_ids = [
                    rid.replace(f"_{appended_id}", "") for rid in reaction_ids
                ]
        reaction_ids = [rid.replace("R_", "") for rid in reaction_ids]

        return reaction_ids

    def get_species_bigg_ids(self, mid="") -> list[str]:
        """Get the species bigg ids of all species

        Args:
            mid (str, optional): When provided only the ids of a specific
            model are returned.
            Defaults to "".

        Raises:
            NotInCombinedModel: the id provided was not in the combined model

        Returns:
            list[str]: list containing the bigg ids
        """
        if mid != "":
            if mid not in self.m_identifiers:
                raise NotInCombinedModel(
                    "The model id provided was not found in the combined model"
                )
            species_ids = self.get_model_specific_species(mid)
            species_ids = [sid.replace(f"_{mid}", "") for sid in species_ids]
        else:
            species_ids = self.m_model.getSpeciesIds()
            for appended_id in self.m_identifiers:
                species_ids = [
                    sid.replace(f"_{appended_id}", "") for sid in species_ids
                ]
        species_ids = [sid.replace("R_", "") for sid in species_ids]

        return species_ids

    def identify_model_from_reaction(self, rid: str) -> str:
        """Given a reaction id get the single model this reaction belonged to

        Args:
            rid (str): reaction id of the kinetic model

        Returns:
            str: id of the old model
        """
        for old_id in self.m_identifiers:
            if rid.endswith(old_id):
                return old_id
        return ""

    def identify_biomass_reaction_for_model(self, mid: str) -> list[str]:
        """Given a model id return the biomass reaction

        Args:
            mid (str): _description_

        Returns:
            str: _description_
        """
        if (
            len(self.m_single_model_biomass_reaction_ids) > 0
            and mid in self.m_identifiers
        ):
            return self.m_single_model_biomass_reaction_ids[
                self.m_identifiers.index(mid)
            ]

        return ""

    def identify_biomass_of_model_from_reaction_id(self, rid) -> str:
        model_id = self.identify_model_from_reaction(rid)
        return self.identify_biomass_reaction_for_model(model_id)

    def get_model_biomass_ids(self) -> dict[str, str]:
        return dict(
            zip(self.m_identifiers, self.m_single_model_biomass_reaction_ids)
        )

    def get_exporters(self) -> list[str]:
        return self.get_transporters()[0]

    def get_importers(self) -> list[str]:
        return self.get_transporters()[1]

    def get_transporters(self) -> list[str]:
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
        for rid in self.getReactionIds():
            reaction: Reaction = self.getReaction(rid)
            if not reaction.is_exchange:
                for reagent in reaction.reagents:
                    sid: str = reagent.getSpecies()
                    species: Species = self.getSpecies(sid)
                    if species.getCompartmentId() == "e":
                        if reagent.coefficient == -1:
                            # Reaction imports external metabolites
                            importers[rid] = importers.get(rid, []) + [sid]
                        elif reagent.coefficient == 1:
                            exporters[rid] = exporters.get(rid, []) + [sid]

        return [exporters, importers]
