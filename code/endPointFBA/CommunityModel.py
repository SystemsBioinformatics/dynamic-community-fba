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
        ids: list[str] = [],
        biomass_reaction_ids: list[str] = [],
    ) -> None:
        super().__init__("combined_model")

        self.createCompartment("e", "extracellular space")
        duplicate_species = cm.create_duplicate_species_dict(models)

        self.m_single_model_identifiers = [model.id for model in models]

        if len(ids) > 0 and len(ids) < len(models):
            raise Exception("Too few ids were provided")
        if len(ids) == 0:
            self.m_identifiers = self.m_single_model_identifiers.copy()
        else:
            self.m_identifiers = ids

        for i in range(0, len(models)):
            model = models[i]
            new_id = self.m_identifiers[i]
            "Save biomass reaction id of old model"
            for rid in biomass_reaction_ids:
                if rid in model.getReactionIds():
                    self.m_single_model_biomass_reaction_ids.append(
                        cm.create_new_id(rid, new_id)
                    )
                    break
            cm.merge_compartments(model, self, new_id)
            cm.merge_species(duplicate_species, model, new_id)
            cm.merge_reactions(model, self, new_id)

    def __str__(self) -> str:
        return (
            f"Model: {self.getId()} was build from "
            f"{[id for id in self.m_single_model_identifiers]}"
        )

    def add_model_to_community(self, model: Model, new_id=None):
        if new_id is None:
            new_id = model.getId()

        duplicate_species = cm.create_duplicate_species_dict([self, model])

        cm.merge_compartments(model, self, new_id)
        cm.merge_species(duplicate_species, model, new_id)
        cm.merge_reactions(model, self, new_id)

        self.m_identifiers.append(new_id)
        self.m_single_model_identifiers.append(model.id)

    def get_model_specific_reactions(self, mid: str) -> list[Reaction]:
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

    def identify_biomass_reaction_for_model(self, mid: str) -> str:
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

    def get_model_biomass_ids(self):
        return dict(
            zip(self.m_identifiers, self.m_single_model_biomass_reaction_ids)
        )
