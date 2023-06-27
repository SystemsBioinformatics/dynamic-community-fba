from cbmpy.CBModel import Model, Species, Reaction
import endPointFBA.helpers.build_community_matrix as cm
from endPointFBA.Exceptions.NotInCombinedModel import NotInCombinedModel


class CombinedModel(Model):
    m_identifiers: list[str]
    m_old_model_identifiers: list[str]

    def __init__(
        self,
        models: list[Model],
        ids: list[str] = [],
        objective_function: str = "",
    ) -> None:
        super().__init__("combined_model")

        self.createCompartment("e", "extracellular space")
        duplicate_species = cm.create_duplicate_species_dict(models)

        self.m_old_model_identifiers = [model.id for model in models]

        if len(ids) > 0 and len(ids) < len(models):
            raise Exception("Too few ids were provided")
        if len(ids) == 0:
            self.m_identifiers = self.m_old_model_identifiers.copy()
        else:
            self.m_identifiers = ids

        for i in range(0, len(models)):
            model = models[i]
            new_id = self.m_identifiers[i]
            cm.merge_compartments(model, self, new_id)
            cm.merge_species(duplicate_species, model, new_id)
            cm.merge_reactions(model, self, new_id)
        if len(objective_function) > 0:
            self.createObjectiveFunction(objective_function)

    def __str__(self) -> str:
        return (
            f"Model: {self.getId()} was build from "
            f"{[id for id in self.m_old_model_identifiers]}"
        )

    def add(self, model: Model, new_id=None):
        if new_id is None:
            new_id = model.getId()

        duplicate_species = cm.create_duplicate_species_dict([self, model])

        cm.merge_compartments(model, self, new_id)
        cm.merge_species(duplicate_species, model, new_id)
        cm.merge_reactions(model, self, new_id)

        self.m_identifiers.append(new_id)
        self.m_old_model_identifiers.append(model.id)

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

    def get_model(self) -> Model:
        return self.m_Model

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
