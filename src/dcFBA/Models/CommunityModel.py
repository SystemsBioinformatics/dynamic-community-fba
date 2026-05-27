# TODO copy old user constraints to the community model (with cbmpy > 0.9?)
# TODO write the new properties of community model so a sbml file and import the community model

import copy
from cbmpy.CBModel import Model, Species, Reaction
from ..Helpers import BuildCommunityMatrix as cm
from ..Exceptions import NotInCombinedModel
from ..Helpers.DynamicSetup import mark_dynamic_species as _mark_dynamic_species_helper

class CommunityModel(Model):
    """
    A CommunityModel represents a combined model built from multiple
    individual models to simulate a community of organisms.

    Attributes:
        custom_model_identifiers (list[str]):
            User-defined identifiers associated with merged models.

        single_model_ids (list[str]):
            Original IDs of merged models.

        single_model_biomass_reaction_ids (list[str]):
            Biomass reaction IDs associated with each merged model.

        merge_extracellular (bool):
            Whether extracellular metabolites/reactions are shared.

        rename_all_species (bool):
            Whether all species are renamed with model-specific suffixes.

        pool_compartment (str | None):
            Shared extracellular pool compartment used when
            `merge_extracellular=False`.
    """

    def __init__(
        self,
        models: list[Model],
        biomass_reaction_ids: list[str],
        ids: list[str] | None = None,
        combined_model_id: str = "combined_model",
        rename_all_species: bool = False,
        merge_extracellular: bool = True,
        exchange_bound_policy: str = "ignore",
        copy_models: bool = False,
        pool_compartment: str | None = None,
    ) -> None:
        """
        Initialize a CommunityModel instance.

        Args:
            models (list[Model]):
                Individual CBModels to combine.

            biomass_reaction_ids (list[str]):
                Biomass reaction IDs corresponding to each model.

            ids (list[str] | None, optional):
                Model-specific identifiers used during merging.
                If None, original model IDs are used.
                Defaults to None.

            combined_model_id (str, optional):
                Identifier of the combined community model.
                Defaults to "combined_model".

            rename_all_species (bool, optional):
                If True, rename all eligible species using
                model-specific suffixes.
                Defaults to False.

            merge_extracellular (bool, optional):
                If True, extracellular metabolites/reactions
                are shared across models.
                Defaults to True.

            copy_models (bool, optional):
                If True, clone input models before modification.
                Defaults to False.

            pool_compartment (str | None, optional):
                Shared extracellular pool compartment used when
                `merge_extracellular=False`.
                Defaults to None.
        """

        super().__init__(combined_model_id)

        # Store model metadata
        self._single_model_ids = [model.getId() for model in models]

        if not ids:
            self._custom_model_identifiers = copy.deepcopy(
                self._single_model_ids
            )
        else:
            self._custom_model_identifiers = copy.deepcopy(ids)

        # # Store merge configuration
        # self._rename_all_species = rename_all_species
        # self._merge_extracellular = merge_extracellular
        # self._exchange_bound_policy = exchange_bound_policy

        # if pool_compartment is not None:
        #     self._pool_compartment = pool_compartment

        # # Validate IDs
        # cm.check_ids(self.custom_model_identifiers, models)

        # Build biomass reaction mapping
        self._single_model_biomass_reaction_ids: list[str] = []

        if len(biomass_reaction_ids) != len(models):
            raise ValueError(
                "The number of biomass_reaction_ids "
                "must match the number of models."
            )

        for biomass_rid, model_id in zip(
            biomass_reaction_ids,
            self.custom_model_identifiers,
        ):

            biomass_id = cm.create_new_id(
                biomass_rid,
                model_id,
            )

            self._single_model_biomass_reaction_ids.append(
                biomass_id
            )

        # Combine model stoichiometric structures
        cm.combine_models(
            models=models,
            new_ids=self.custom_model_identifiers,
            modelId=combined_model_id,
            rename_all_species=rename_all_species,
            merge_extracellular=merge_extracellular,
            exchange_bound_policy=exchange_bound_policy,
            copy_models=copy_models,
            target_model=self,
            pool_compartment=pool_compartment,
        )

        # Preserve gene activity checks
        self.__check_gene_activity__ = any(
            getattr(model, "__check_gene_activity__", False)
            for model in models
        )


    def __str__(self) -> str:
        """
        Return a string representation of the CommunityModel.
        """

        return (
            f"Model: {self.getId()} was built from "
            f"{self.single_model_ids}"
        )


    @property
    def custom_model_identifiers(self) -> list[str]:
        return self._custom_model_identifiers

    @property
    def single_model_ids(self) -> list[str]:
        return self._single_model_ids

    @property
    def single_model_biomass_reaction_ids(self) -> list[str]:
        return self._single_model_biomass_reaction_ids

    @property
    def merge_extracellular(self) -> bool:
        return self._merge_extracellular

    @property
    def exchange_bound_policy(self) -> str:
        return self._exchange_bound_policy

    @property
    def rename_all_species(self) -> bool:
        return self._rename_all_species

    @property
    def pool_compartment(self) -> str | None:
        return getattr(self, "_pool_compartment", None)


    def clone(self) -> "CommunityModel":
        """
        Create a deep copy of the CommunityModel instance.

        Returns:
            CommunityModel:
                Deep copy of the CommunityModel.
        """

        new_instance = super().clone()

        new_instance._custom_model_identifiers = copy.deepcopy(
            self.custom_model_identifiers
        )

        new_instance._single_model_biomass_reaction_ids = (
            copy.deepcopy(
                self.single_model_biomass_reaction_ids
            )
        )

        new_instance._single_model_ids = copy.deepcopy(
            self.single_model_ids
        )

        new_instance._rename_all_species = (
            self.rename_all_species
        )

        new_instance._merge_extracellular = (
            self.merge_extracellular
        )

        new_instance._exchange_bound_policy = (
            self.exchange_bound_policy
        )

        if hasattr(self, "_pool_compartment"):
            new_instance._pool_compartment = (
                self._pool_compartment
            )

        return new_instance

    # TODO maybe implement __eq__() method

    def add_model_to_community(
        self,
        model: Model,
        biomass_reaction_id: str,
        new_id: str | None = None,
        copy_model: bool = False,
    ) -> None:
        """
        Add a model to the CommunityModel.

        Args:
            model (Model):
                Model to be added.

            biomass_reaction_id (str):
                Biomass reaction ID of the added model.

            new_id (str | None, optional):
                User-set model-specific identifier used during merging.
                If None, the original model ID is used.

            copy_model (bool, optional):
                If True, clone the model before merging.
                Defaults to False.
        """

        if new_id is None:
            new_id = model.getId()

        if new_id in self.custom_model_identifiers:
            raise ValueError(
                f"Model identifier '{new_id}' already exists "
                "in the community model."
            )

        # Merge directly into the current community
        cm.combine_models(
            models=[model],
            new_ids=[new_id],
            target_model=self,
            rename_all_species=self.rename_all_species,
            merge_extracellular=self.merge_extracellular,
            exchange_bound_policy=self.exchange_bound_policy,
            copy_models=copy_model,
            pool_compartment=self.pool_compartment,
        )

        # Update metadata bookkeeping
        self._single_model_ids.append(model.getId())
        self._custom_model_identifiers.append(new_id)

        self._single_model_biomass_reaction_ids.append(
            cm.create_new_id(
                biomass_reaction_id,
                new_id,
            )
        )

        self.__check_gene_activity__ = (
            self.__check_gene_activity__
            or getattr(model, "__check_gene_activity__", False)
        )


    def remove_model_from_community(
        self,
        mid: str,
        strategy: str = "cleanup",
    ) -> None:
        """
        Remove a model from the CommunityModel.

        Args:
            mid (str):
                Identifier of the model to remove.
                Both original model IDs and custom model identifiers
                are accepted.

            strategy (str, optional):
                Removal strategy.

                Available strategies:
                    - "cleanup":
                        Remove model-specific reactions and then
                        delete orphan/non-reacting species.

                    - "legacy":
                        Remove reactions and species based on
                        identifier/compartment matching.

                Defaults to "cleanup".

        Raises:
            ValueError:
                If the model is not present in the community
                or if the strategy is invalid.
        """

        if mid in self.single_model_ids:
            index = self.single_model_ids.index(mid)

        elif mid in self.custom_model_identifiers:
            index = self.custom_model_identifiers.index(mid)

        else:
            raise ValueError(
                f"Model '{mid}' not found in community."
            )

        mid = self.custom_model_identifiers[index]

        if strategy == "cleanup":

            # Remove all model-specific reactions
            reactions_to_delete = [
                rid
                for rid in self.getReactionIds()
                if rid.endswith(f"_{mid}")
            ]

            for rid in reactions_to_delete:
                self.deleteReactionAndBounds(rid)

            # Remove orphan/non-reacting species
            self.deleteNonReactingSpecies(simulate=False)

        elif strategy == "legacy":

            # Remove reactions associated with the model
            for rid in list(self.getReactionIds()):

                reaction: Reaction = self.getReaction(rid)

                if (
                    mid in rid
                    or mid in reaction.getCompartmentId()
                ):
                    self.deleteReactionAndBounds(rid)

            # Remove species associated with the model
            for sid in list(self.getSpeciesIds()):

                species: Species = self.getSpecies(sid)

                if (
                    mid in sid
                    or mid in species.getCompartmentId()
                ):
                    self.deleteSpecies(sid)

        else:
            raise ValueError(
                f"Unknown removal strategy '{strategy}'. "
                "Available strategies are: "
                "'cleanup', 'legacy'."
            )
        
        # TODO:
        # optionally remove orphan genes / GPRs

        # Update bookkeeping
        del self._single_model_ids[index]
        del self._custom_model_identifiers[index]
        del self._single_model_biomass_reaction_ids[index]


    def get_model_specific_reactions(
        self,
        mid: str,
    ) -> list[str]:
        """
        Return reaction IDs specific to a given model.

        Args:
            mid (str):
                Model identifier.

        Raises:
            NotInCombinedModel:
                If the model identifier is not present
                in the CommunityModel.

        Returns:
            list[str]:
                Model-specific reaction IDs.
        """

        if mid not in self.custom_model_identifiers:
            raise NotInCombinedModel(
                "The model id provided was not found "
                "in the combined model"
            )

        return [
            rid
            for rid in self.getReactionIds()
            if rid.endswith(f"_{mid}")
        ]


    def get_model_specific_species(
        self,
        mid: str,
    ) -> list[str]:
        """
        Return species IDs associated with a given model.

        Species ownership is inferred from participation
        in model-specific reactions.

        Args:
            mid (str):
                Model identifier.

        Raises:
            NotInCombinedModel:
                If the model identifier is not present.

        Returns:
            list[str]:
                Species associated with the model.
        """

        if mid not in self.custom_model_identifiers:
            raise NotInCombinedModel(
                "The model id provided was not found "
                "in the combined model"
            )

        species_ids = set()

        for rid in self.get_model_specific_reactions(mid):

            reaction: Reaction = self.getReaction(rid)

            species_ids.update(
                reaction.getSpeciesIds()
            )

        return sorted(species_ids)


    def get_reaction_bigg_ids(
        self,
        mid: str = "",
    ) -> list[str]:
        """
        Return BIGG-style reaction IDs.

        Args:
            mid (str, optional):
                If provided, only reactions associated
                with the specified model are returned.

        Raises:
            NotInCombinedModel:
                If the model identifier is invalid.

        Returns:
            list[str]:
                BIGG-style reaction IDs.
        """

        if mid:

            if mid not in self.custom_model_identifiers:
                raise NotInCombinedModel(
                    "The model id provided was not found "
                    "in the combined model"
                )

            reaction_ids = self.get_model_specific_reactions(mid)

        else:
            reaction_ids = self.getReactionIds()

        cleaned_ids = []

        for rid in reaction_ids:

            cleaned = rid

            for appended_id in self.custom_model_identifiers:

                suffix = f"_{appended_id}"

                if cleaned.endswith(suffix):
                    cleaned = cleaned[:-len(suffix)]
                    break

            if cleaned.startswith("R_"):
                cleaned = cleaned[2:]

            cleaned_ids.append(cleaned)

        return cleaned_ids


    def get_species_bigg_ids(
        self,
        mid: str = "",
    ) -> list[str]:
        """
        Return BIGG-style species IDs.

        Args:
            mid (str, optional):
                If provided, only species associated
                with the specified model are returned.

        Raises:
            NotInCombinedModel:
                If the model identifier is invalid.

        Returns:
            list[str]:
                BIGG-style species IDs.

        Notes:
            Under shared extracellular mode, species may
            participate in multiple models simultaneously.
        """

        if mid:

            if mid not in self.custom_model_identifiers:
                raise NotInCombinedModel(
                    "The model id provided was not found "
                    "in the combined model"
                )

            species_ids = self.get_model_specific_species(mid)

        else:
            species_ids = self.getSpeciesIds()

        cleaned_ids = []

        for sid in species_ids:

            cleaned = sid

            for appended_id in self.custom_model_identifiers:

                suffix = f"_{appended_id}"

                if cleaned.endswith(suffix):
                    cleaned = cleaned[:-len(suffix)]
                    break

            if cleaned.startswith("M_"):
                cleaned = cleaned[2:]

            cleaned_ids.append(cleaned)

        return cleaned_ids


    def identify_model_from_reaction(
        self,
        rid: str,
    ) -> str:
        """
        Given a reaction ID, identify the originating model.

        Args:
            rid (str):
                Reaction ID.

        Returns:
            str:
                Model identifier associated with the reaction.
                Returns an empty string for shared/global reactions.
        """

        for mid in self.custom_model_identifiers:

            if rid.endswith(f"_{mid}"):
                return mid

        return ""


    def identify_biomass_reaction_for_model(
        self,
        mid: str,
    ) -> str:
        """
        Return the biomass reaction associated with a model.

        Args:
            mid (str):
                Model identifier.

        Returns:
            str:
                Biomass reaction ID.
                Returns an empty string if unavailable.
        """

        if (
            self.single_model_biomass_reaction_ids
            and mid in self.custom_model_identifiers
        ):

            return self.single_model_biomass_reaction_ids[
                self.custom_model_identifiers.index(mid)
            ]

        return ""


    def identify_biomass_of_model_from_reaction_id(
        self,
        rid: str,
    ) -> str:
        """
        Identify the biomass reaction associated with
        the model owning a reaction.

        Args:
            rid (str):
                Reaction ID.

        Returns:
            str:
                Biomass reaction ID.
                Returns an empty string if no model
                association exists.
        """

        model_id = self.identify_model_from_reaction(rid)

        if model_id == "":
            return ""

        return self.identify_biomass_reaction_for_model(
            model_id
        )


    def get_model_biomass_ids(self) -> dict[str, str]:
        """
        Return the mapping between model identifiers
        and biomass reaction IDs.

        Returns:
            dict[str, str]:
                Mapping:
                    model_id -> biomass_reaction_id
        """

        return dict(
            zip(
                self.custom_model_identifiers,
                self.single_model_biomass_reaction_ids,
            )
        )

    def mark_dynamic_species(
            self, 
            auto: bool = True, 
            include: list[str] = None, 
            exclude: list[str] = None
        ) -> None:
            """
            Marks species in the community model as either dynamic or non-dynamic 
            for Dynamic FBA simulations (e.g., EndPointFBA, DynamicJointFBA).
            
            In dynamic simulations, extracellular metabolites usually accumulate or deplete 
            over time. However, certain metabolites (like O2, CO2, or H2O) are often 
            treated as infinite sinks/sources or assumed to be in a quasi-steady state 
            with the environment (instantly replenished/consumed). 
            
            This method is a wrapper around the standalone helper `mark_dynamic_species`, 
            which sets the custom `dcFBA_dynamic` attribute on the model's Species:
            - Species explicitly in `include` will be marked True (will accumulate over time).
            - Species explicitly in `exclude` will be marked False (quasi-steady state; 
                their concentrations will not be tracked/linked over time).
            - If `auto=True`, it uses a heuristic to automatically mark valid extracellular 
                or pooled species as dynamic if they participate in both an exchange and a 
                non-exchange reaction.
                
            Args:
                auto (bool): If True, applies auto-detection heuristic to extracellular/pool species.
                include (list[str]): List of species IDs to explicitly mark as dynamic.
                exclude (list[str]): List of species IDs to explicitly exclude from being dynamic.
            """
            # 1. Determine the correct extracellular compartments for this specific community model
            extracellular_comps = ["e", "extracellular"]
            pool_comp = getattr(self, "pool_compartment", None)
            if pool_comp:
                extracellular_comps.append(pool_comp)

            # 2. Pass self (which is a cbmpy Model) and the arguments to the standalone helper
            _mark_dynamic_species_helper(
                model=self,
                auto=auto,
                include=include,
                exclude=exclude,
                extracellular_compartments=extracellular_comps
            )