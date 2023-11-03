import copy
from cbmpy.CBModel import Model, Species
from ..Helpers import BuildCommunityMatrix as cm
from ..Exceptions import NotInCombinedModel


class CommunityModel(Model):
    """
    A CommunityModel represents a combined model built from multiple individual models
    to simulate a community of organisms.

    Attributes:
        m_identifiers (list[str]): List of user-specified identifiers for the
            individual models.
        m_single_model_ids (list[str]): List of IDs of the individual models.
        m_single_model_biomass_reaction_ids (list[str]): List of biomass
            reaction IDs of the individual models.
    """

    m_identifiers: list[str]
    m_single_model_ids: list[str]
    m_single_model_biomass_reaction_ids: list[str] = []

    def __init__(
        self,
        models: list[Model],
        biomass_reaction_ids: list[str],
        ids: list[str] = [],
        combined_model_id: str = "combined_model",
    ) -> None:
        """
        Initialize a CommunityModel instance.

        Args:
            models (list[Model]): List of individual models to combine.
            biomass_reaction_ids (list[str]): List of biomass reaction IDs for
                each individual model.
            ids (list[str], optional): List of user-specified identifiers for
                the individual models.
                Defaults to an empty list.
            combined_model_id (str, optional): The ID for the combined CommunityModel.
                Defaults to "combined_model".

        Raises:
            Exception: If too few IDs are provided compared to the number of
                models.

        """
        super().__init__(combined_model_id)

        self.m_single_model_ids = [model.id for model in models]

        if len(ids) == 0:
            self.m_identifiers = self.m_single_model_ids
        else:
            self.m_identifiers = copy.deepcopy(ids)

        cm.check_ids(self.m_identifiers, models)

        self.createCompartment("e", "extracellular space")

        duplicate_species = cm.create_duplicate_species_dict(models)

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
            cm.merge_genes(model, self, new_id)
            cm.merge_compartments(model, self, new_id)
            cm.merge_species(duplicate_species, model, new_id)
            cm.merge_reactions(model, self, new_id)
            cm.setGeneProteinAssociations(model, self, new_id)

        self.__check_gene_activity__ = any(
            [m.__check_gene_activity__ for m in models]
        )

        # TODO Old User Constraints arent added to the new community model,
        # should we implement this?

    def __str__(self) -> str:
        """
        Return a string representation of the CommunityModel.

        Returns:
            str: A string representation of the CommunityModel.

        """
        return (
            f"Model: {self.getId()} was build from "
            f"{[id for id in self.m_single_model_ids]}"
        )

    def get_model_ids(self):
        return self.m_identifiers

    def clone(self):
        """
        Create a deep copy of the CommunityModel instance.

        Returns:
            CommunityModel: A new deep copy of the CommunityModel instance.
        """

        new_instance = super().clone()
        new_instance.m_identifiers = copy.deepcopy(self.m_identifiers)
        new_instance.m_single_model_biomass_reaction_ids = copy.deepcopy(
            self.m_single_model_biomass_reaction_ids
        )
        new_instance.m_single_model_ids = copy.deepcopy(
            self.m_single_model_ids
        )

        return new_instance

    # def clone1(self):
    #     return copy.deepcopy(self)

    # TODO maybe implement __eq__() method

    def add_model_to_community(
        self, model: Model, biomass_reaction: str, new_id: str = None
    ) -> None:
        """
        Adds a model to the CommunityModel.

        Args:
            model (Model): The model to be added.
            biomass_reaction (str): The reaction ID of the biomass reaction of
                the new model.
            new_id (str, optional): The user-set identifier for the model.
                Defaults to None. If set to None, the model ID will be used.

        """
        if new_id is None:
            new_id = model.getId()

        duplicate_species = cm.create_duplicate_species_dict([self, model])

        cm.merge_compartments(model, self, new_id)
        cm.merge_species(duplicate_species, model, new_id)
        cm.merge_reactions(model, self, new_id)

        self.m_identifiers.append(new_id)
        self.m_single_model_ids.append(model.id)
        self.m_single_model_biomass_reaction_ids.append(biomass_reaction)

    def remove_model_from_community(self, mid: str) -> None:
        """
        Remove a model from the CommunityModel.

        Args:
            mid (str): The identifier of the model to be removed.

        Raises:
            Exception: If the provided model ID is not in the CommunityModel.

        """
        if mid in self.m_single_model_ids:
            index = self.m_single_model_ids.index(mid)
        elif mid in self.m_identifiers:
            index = self.m_identifiers.index(mid)
        else:
            raise Exception("Model not in community")

        mid = self.m_identifiers[index]

        for rid in self.getReactionIds():
            if mid in rid:
                self.deleteReactionAndBounds(rid)

        for sid in self.getSpeciesIds():
            species: Species = self.getSpecies(sid)
            if mid in species.getCompartmentId():
                self.deleteSpecies(sid)

        self.m_identifiers.remove(mid)
        del self.m_single_model_ids[index]

        del self.m_single_model_biomass_reaction_ids[index]

    def get_model_specific_reactions(self, mid: str) -> list[str]:
        """
        Returns a list of reaction IDs specific to the given model.

        Args:
            mid (str): The identifier of the model.

        Raises:
            NotInCombinedModel: If the provided model ID is not in the
            CommunityModel.

        Returns:
            list[Reaction]: A list of reaction IDs specific to the given model.

        """

        if mid not in self.m_identifiers:
            raise NotInCombinedModel(
                "The model id provided was not found in the combined model"
            )
        ans = []
        for rid in self.getReactionIds():
            if mid in rid:
                ans.append(rid)
        return ans

    def get_model_specific_species(self, mid: str) -> list[str]:
        """
        Returns a list of species IDs specific to the given model.

        Args:
            mid (str): The identifier of the model.

        Raises:
            NotInCombinedModel: If the provided model ID is not in the
            CommunityModel.

        Returns:
            list[Species]: A list of species IDs specific to the given model.

        """
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
        """Get the reaction BIGG IDs of all reactions

        Args:
            mid (str, optional): If a model id is provided only reactions from
            the specific model are returned
            Defaults to "".

        Raises:
            NotInCombinedModel: the id provided was not in the combined model

        Returns:
            list[str]: list containing the BIGG ids
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
            if f"_{old_id}" in rid:
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
        if model_id == "":
            return ""
        return self.identify_biomass_reaction_for_model(model_id)

    def get_model_biomass_ids(self) -> dict[str, str]:
        return dict(
            zip(self.m_identifiers, self.m_single_model_biomass_reaction_ids)
        )
