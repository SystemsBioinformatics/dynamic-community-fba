"""Helper functions for building the CommunityModel 
"""

from cbmpy.CBModel import (
    Model,
    Gene,
    GeneProteinAssociation,
    Compartment,
    Reaction,
    Reagent,
    Species,
)
import cbmpy
import warnings
import numpy as np

# Add extracellular ids here
default_extracellular_compartments = ["e", "extracellular"]


def load_n_models(files: list[str]) -> list[Model]:
    """
    Load multiple CBModels from files and return a list of loaded models.

    Args:
        files (list[str]): A list of file paths of SBML model files.

    Returns:
        list[Model]: A list of loaded CBModels.

    """

    models: list[Model] = []
    for file in files:
        m: Model = cbmpy.loadModel(file)
        if m is not None:
            models.append(m)
        else:
            print(
                f"Model {file} could not be loaded. \n Are you sure it is in a"
                " correct sbml model format?"
            )
            continue
    return models


def combine_models(
    models: list[Model],
    new_ids: list[str] | None = None,
    modelId: str = "combined_model",
    objective_function: str | None = None,
    rename_all_species: bool = False,
    merge_extracellular: bool = True,
    exchange_bound_policy: str = "ignore",
    copy_models: bool = False,
    target_model: Model | None = None,
    pool_compartment: str | None = None,
) -> Model:
    """
    Combine multiple CBModels into a single community model.

    Models are merged by:
    - copying genes, compartments, species, and reactions;
    - renaming model-specific entities when needed;
    - optionally merging extracellular metabolites into a shared pool.

    Args:
        models (list[Model]):
            List of CBModels to combine.

        new_ids (list[str] | None, optional):
            Model-specific suffixes appended to renamed entities
            (species, reactions, genes, compartments).
            If None, model IDs are used.
            Defaults to None. 
            For backward compatibility, empty list behaves like None.

        modelId (str, optional):
            Identifier of the combined community model.
            Defaults to "combined_model".

        objective_function (str, optional):
            Reaction ID to use as objective function.
            Defaults to None (no objective is set).

        rename_all_species (bool, optional):
            If True, rename all eligible species using model-specific
            suffixes. If False, rename only duplicated species.
            Defaults to False.

        merge_extracellular (bool, optional):
            If True, extracellular metabolites and exchange reactions
            are shared across models.
            If False, each model keeps its own extracellular species,
            and a shared extracellular pool is created to connect them.
            Defaults to True.

        exchange_bound_policy (str, optional):
            Policy used to merge bounds of shared exchange
            reactions when `merge_extracellular=True`.

            Supported policies are:

            - "first":
                Keep bounds from the first encountered model.

            - "union":
                Use the most permissive bounds across models.

            - "intersection":
                Use the most restrictive compatible bounds
                across models.

            - "ignore":
                Ignore model-specific exchange bounds and
                assign default unconstrained bounds to shared
                exchange reactions.

            This parameter only affects shared exchange
            reactions created when
            `merge_extracellular=True`.

            Defaults to "ignore".

        copy_models (bool, optional):
            If True, clone input models before modification.
            This preserves the original models at the cost of increased
            memory usage.
            Defaults to False.

        target_model (Model | None, optional):
            Existing community model to which new models
            are appended.
            If provided, compatibility checks are performed
            against the existing merge configuration.
            Defaults to None.

        pool_compartment (str | None, optional):
            Identifier of the shared extracellular pool
            compartment used when
            `merge_extracellular=False`.
            If None, a default pool compartment named
            `"pool"` is created automatically.
            Defaults to None.

    Returns:
        Model:
            Combined community CBModel.

    WARNING:
        By default, this function modifies input models in-place:
        - species may be renamed;
        - species may be deleted/recreated;
        - reactions may be rewired;
        - compartments may be renamed.

        To preserve original models:
        - save models beforehand;
        - manually clone models;
        - or set copy_models=True.
    """

    # validate model-specific suffixes
    if not new_ids:
        new_ids = [model.getId() for model in models]

    check_ids(new_ids, models)

    # validate selected exchange policy 
    valid_exchange_policies = {
        "first",
        "union",
        "intersection",
        "ignore",
    }

    if exchange_bound_policy not in valid_exchange_policies:
        raise ValueError(
            f"Unknown exchange_bound_policy "
            f"'{exchange_bound_policy}'. "
            f"Valid policies are: "
            f"{sorted(valid_exchange_policies)}"
        )
    
    if (
        merge_extracellular
        and exchange_bound_policy == "first"
    ):
        warnings.warn(
            (
                "Shared extracellular exchange reactions "
                "inherit bounds from the first encountered "
                "model.\n"
                "Consider setting "
                "'exchange_bound_policy' explicitly."
            ),
            UserWarning,
        )

    # setup combined model
    if copy_models:
        models = [m.clone() for m in models]

    if target_model is not None:

        if (
            not merge_extracellular
            and pool_compartment is None
            and hasattr(target_model, "_pool_compartment")
        ):
            pool_compartment = target_model._pool_compartment

        check_target_model_compatibility(
            target_model,
            rename_all_species,
            merge_extracellular,
            pool_compartment,
            exchange_bound_policy,
        )

    combined_model = target_model or Model(modelId)

    if target_model is None:
        combined_model.setId(modelId)
    
    if merge_extracellular and "e" not in combined_model.getCompartmentIds():
        combined_model.createCompartment("e", "extracellular space")

    # Store merge configuration for future compatibility checks
    combined_model._merge_extracellular = merge_extracellular
    combined_model._rename_all_species = rename_all_species
    combined_model._exchange_bound_policy = exchange_bound_policy

    _append_models_to_target(
        models=models,
        target_model=combined_model,
        new_ids=new_ids,
        rename_all_species=rename_all_species,
        merge_extracellular=merge_extracellular,
        exchange_bound_policy=exchange_bound_policy,
        pool_compartment=pool_compartment,
    )

    if objective_function:
        combined_model.createObjectiveFunction(objective_function)

    return combined_model


def check_target_model_compatibility(
    target_model: Model,
    rename_all_species: bool,
    merge_extracellular: bool,
    pool_compartment: str | None,
    exchange_bound_policy: str,
) -> None:
    """
    """
    if(
        hasattr(target_model, "_rename_all_species")
        and target_model._rename_all_species != rename_all_species        
    ):
        raise ValueError(
            f"The target model has 'rename_all_species' set to {target_model._rename_all_species}\n"
            f"while the current function has it set to {rename_all_species}.\n"
            "Please check manually."
        )
    
    if(
        hasattr(target_model, "_merge_extracellular")
        and target_model._merge_extracellular != merge_extracellular           
    ):
        raise ValueError(
            f"The target model has 'merge_extracellular' set to {target_model._merge_extracellular}\n"
            f"while the current function has it set to {merge_extracellular}.\n"
            "Please check manually."
        )

    if(
        hasattr(target_model, "_pool_compartment")
        and target_model._pool_compartment != pool_compartment           
    ):
        raise ValueError(
            f"The target model has 'pool_compartment' set to {target_model._pool_compartment}\n"
            f"while the current function has it set to {pool_compartment}.\n"
            "Please check manually."
        )
    
    if (
        hasattr(target_model, "_exchange_bound_policy")
        and target_model._exchange_bound_policy
        != exchange_bound_policy
    ):
        raise ValueError(
            f"The target model has 'exchange_bound_policy' set to {target_model._exchange_bound_policy}\n"
            f"while the current function has it set to {exchange_bound_policy}.\n"
            "Please check manually."
        )
             
 
def _append_models_to_target(
    models: list[Model],
    target_model: Model,
    new_ids: list[str],
    rename_all_species: bool = False,
    merge_extracellular: bool = True,
    exchange_bound_policy: str = "ignore",
    pool_compartment: str | None = None,
) -> None:
    """
    Append one or more CBModels into an existing target model.

    Models are merged by:
    - copying genes, compartments, species, and reactions;
    - renaming model-specific entities when needed;
    - optionally merging extracellular metabolites into
      a shared extracellular space or pool.

    This function performs the low-level merge operations and assumes that:
    - the target model already exists;
    - compatibility checks were already performed;
    - merge configuration metadata was already initialized.

    Args:
        models (list[Model]):
            List of CBModels to append.

        target_model (Model):
            Existing model receiving the merged content.

        new_ids (list[str]):
            Model-specific suffixes appended to renamed entities
            (species, reactions, genes, compartments).

        rename_all_species (bool, optional):
            If True, rename all eligible species using model-specific
            suffixes. If False, rename only duplicated species.
            Defaults to False.

        merge_extracellular (bool, optional):
            If True, extracellular metabolites and exchange reactions
            are shared across models.

            If False, each model keeps its own extracellular species,
            and a shared extracellular pool is used to connect them.
            Defaults to True.

        exchange_bound_policy (str, optional):
            Policy used to merge bounds of shared exchange
            reactions when `merge_extracellular=True`.

            Ignored when `merge_extracellular=False`.

            Defaults to "ignore".

        pool_compartment (str | None, optional):
            Shared pool compartment used when
            `merge_extracellular=False`.

            If None, the default pool compartment is used.
            Defaults to None.

    Returns:
        None

    WARNING:
        This function modifies both:
        - input models;
        - the target model.

        In particular:
        - species may be renamed;
        - species may be deleted/recreated;
        - reactions may be rewired;
        - compartments may be renamed.

        To preserve original models:
        - save models beforehand;
        - or manually clone models.
    """

    # Duplicate species must be identified before models are
    # modified in-place by merge_species()
    duplicate_species = (
        None
        if rename_all_species
        else create_duplicate_species_dict(models)
    )

    for model, new_id in zip(models, new_ids):

        merge_genes(model, target_model, new_id)

        merge_compartments(
            model,
            target_model,
            new_id,
            merge_extracellular=merge_extracellular,
        )

        merge_species(
            model,
            new_id,
            duplicate_species=duplicate_species,
            rename_all_species=rename_all_species,
            merge_extracellular=merge_extracellular,
        )

        merge_reactions(
            model,
            target_model,
            new_id,
            merge_extracellular=merge_extracellular,
            exchange_bound_policy=exchange_bound_policy,
        )

        setGeneProteinAssociations(
            model,
            target_model,
            new_id,
        )

    if not merge_extracellular:

        connect_external_pool(
            target_model,
            new_ids,
            pool_compartment,
        )


def create_duplicate_species_dict(models: list[Model]) -> dict[str, int]:
    """
    Create a dictionary of duplicate species and their occurrence count from
    multiple models. This dictionary will be used later on the know which
    species to rename

    Args:
        models (list[Model]): A list of CBModels.

    Returns:
        dict[str, int]: A dictionary mapping duplicate species IDs to their
        occurrence count.

    """

    combined_species_list = [
        item for sublist in models for item in sublist.getSpeciesIds()
    ]

    species_dict = {}
    for m in combined_species_list:
        if m in species_dict.keys():
            species_dict[m] += 1
        else:
            species_dict[m] = 1

    return {k: v for k, v in species_dict.items() if v >= 2}


def merge_genes(model: Model, combined_model: Model, new_id: str) -> None:
    """Copy the genes from the sub-models to the new model

    Args:
        model (Model): Model from which genes are copied
        combined_model (Model): The new model
        new_id (str): Gene suffix to be used for the new model
    """
    for gid in model.getGeneIds():
        gene: Gene = model.getGene(gid)
        new_gene = gene.clone()
        new_gene.setId(create_new_id(gene.id, new_id))
        combined_model.addGene(new_gene)


def setGeneProteinAssociations(
    model: Model, combined_model: Model, new_id: str
) -> None:
    """Copy the gene protein associations from the sub models
    to the combined model

    Args:
        model (Model): sub-model from which the gene protein associations are copied
        combined_model (Model): The new community model
        new_id (str): The id suffix to be used

    Raises:
        Exception: If the Reaction is not found in the combined model
            throw an error
        Exception: If the gene was not found in the combined model
            throw an error
    """
    new_dict_ids = {
        create_new_id(rid, new_id): list(
            map(lambda gid: create_new_id(gid, new_id), gene_id)
        )
        for rid, gene_id in model.getAllGeneProteinAssociations().items()
    }

    for rid, gls in new_dict_ids.items():
        if rid in combined_model.getReactionIds():
            gid = "{}_assoc".format(rid)
            gpr = GeneProteinAssociation(gid, rid)
            combined_model.addGPRAssociation(gpr)
            for gene_id in gls:
                if gene_id in combined_model.getGeneIds():
                    gpr.createAssociationAndGeneRefsFromString(gene_id)
                else:
                    raise Exception(
                        "something went wrong with setting the gene_id in the model"
                    )
        else:
            raise Exception("Reaction not recognized")


def merge_compartments(
    model: Model,
    combined_model: Model,
    new_id: str,
    merge_extracellular: bool = True,
    extracellular_compartments: list[str] = default_extracellular_compartments,
) -> None:
    """
    Copy compartments from a source model into the combined
    community model.

    If merge_extracellular=True, extracellular compartments
    are skipped because models share a common extracellular
    space.

    If merge_extracellular=False, extracellular compartments
    are copied and renamed like all other compartments.

    Args:
        model (Model):
            Source CBModel.

        combined_model (Model):
            Target community CBModel.

        new_id (str):
            Model-specific suffix appended to compartment IDs.

        merge_extracellular (bool, optional):
            Whether extracellular compartments are shared across
            models. Defaults to True.

        extracellular_compartments (list[str], optional):
            Compartment IDs considered extracellular.
    """

    compartment: Compartment

    for compartment in model.compartments:

        is_external = compartment.id in extracellular_compartments

        if is_external and merge_extracellular:
            continue

        combined_model.createCompartment(
            create_new_id(compartment.id, new_id),
            compartment.name,
            size=compartment.size,
            dimensions=compartment.dimensions,
        )


def set_exchange_reactions(model: Model) -> list[str]:
    """
    Identify and validate exchange reactions in a model.

    Exchange reactions are detected using:
    - reaction IDs starting with "R_EX";
    - reactions already marked as exchange reactions;
    - SBO term "SBO:0000627".

    Reactions are considered valid exchange reactions only if:
    - they contain exactly one reagent;
    - the stoichiometric coefficient is negative.

    Invalid exchange reactions are flagged with warnings and
    their `is_exchange` attribute is reset to False.

    Args:
        model (Model):
            Source CBModel.

    Returns:
        list[str]:
            IDs of anomalous exchange reactions.
    """

    mod_exchange_reactions = model.getExchangeReactionIds()
    anomalous_exchanges: list[str] = []

    for reaction_id in model.getReactionIds():

        reaction: Reaction = model.getReaction(reaction_id)

        if (
            reaction_id.startswith("R_EX")
            or reaction_id in mod_exchange_reactions
            or reaction.getSBOterm() == "SBO:0000627"
        ):

            reaction.is_exchange = True

            stoich = reaction.getStoichiometry()

            # Exchange reactions should involve exactly one reagent
            if len(stoich) != 1:

                anomalous_exchanges.append(reaction_id)
                reaction.is_exchange = False

                warnings.warn(
                    (
                        f"Exchange reaction {reaction_id} "
                        "contains more than one reagent.\n"
                        "Setting 'is_exchange' to False.\n"
                        "Please check manually."
                    ),
                    UserWarning,
                )

                continue

            coef, _ = stoich[0]

            # Exchange reactions should consume the exchanged species
            if coef >= 0:

                anomalous_exchanges.append(reaction_id)
                reaction.is_exchange = False

                warnings.warn(
                    (
                        f"Exchange reaction {reaction_id} "
                        "has a non-negative coefficient.\n"
                        "Setting 'is_exchange' to False.\n"
                        "Please check manually."
                    ),
                    UserWarning,
                )

    return anomalous_exchanges


def merge_reactions(
    model: Model,
    combined_model: Model,
    new_id: str,
    merge_extracellular: bool = True,
    exchange_bound_policy: str = "ignore",
) -> list[str]:
    """
    Copy reactions from a source model into the combined
    community model.

    If merge_extracellular=True:
    - exchange reactions are shared across models;
    - internal reactions are renamed using model-specific suffixes.

    If merge_extracellular=False:
    - all reactions are renamed using model-specific suffixes.

    Args:
        model (Model):
            Source CBModel.

        combined_model (Model):
            Target community CBModel.

        new_id (str):
            Model-specific suffix appended to reaction IDs.

        merge_extracellular (bool, optional):
            Whether extracellular exchange reactions are shared
            across models. Defaults to True.

        exchange_bound_policy (str, optional):
            Policy used to merge bounds of shared exchange
            reactions when `merge_extracellular=True`.

            Ignored when `merge_extracellular=False`.

            Defaults to "ignore".

    Returns:
        list[str]:
            Exchange reaction IDs present in the combined model.
    """

    # Validate and annotate exchange reactions in-place
    set_exchange_reactions(model)

    cmod_exchange_reactions = set(combined_model.getExchangeReactionIds())

    for reaction_id in model.getReactionIds():

        reaction: Reaction = model.getReaction(reaction_id)

        if merge_extracellular:

            # Shared extracellular exchange reactions
            if reaction.is_exchange:

                if reaction_id not in cmod_exchange_reactions:

                    copy_reaction(
                        model,
                        combined_model,
                        reaction.id,
                    )

                    cmod_exchange_reactions.add(reaction.id)

                else:

                    existing_reaction = combined_model.getReaction(
                        reaction_id
                    )

                    merge_exchange_bounds(
                        existing_reaction,
                        reaction,
                        exchange_bound_policy,
                    )

            else:

                # Internal reactions are model-specific
                copy_reaction(
                    model,
                    combined_model,
                    reaction.id,
                    altrid=create_new_id(reaction_id, new_id),
                )

        else:

            # All reactions remain model-specific
            altrid = create_new_id(reaction_id, new_id)

            copy_reaction(
                model,
                combined_model,
                reaction.id,
                altrid=altrid,
            )

            if reaction.is_exchange:
                cmod_exchange_reactions.add(altrid)

    return list(cmod_exchange_reactions)


def merge_exchange_bounds(
    existing_reaction: Reaction,
    incoming_reaction: Reaction,
    policy: str,
) -> None:
    """
    Merge bounds of shared exchange reactions according
    to the selected policy.

    Policies
    --------
    first:
        Keep bounds from the first encountered model.

    union:
        Use the most permissive bounds across models:
        - lower bound = minimum lower bound
        - upper bound = maximum upper bound

    intersection:
        Use the most restrictive compatible bounds:
        - lower bound = maximum lower bound
        - upper bound = minimum upper bound

    ignore:
        Ignore model-specific exchange bounds and assign
        unconstrained bounds to the shared exchange
        reaction, which now represent environmental 
        availability. Organism uptake capabilities
        must be now represented at transport level.

    Args:
        existing_reaction (Reaction):
            Exchange reaction already present in the
            combined model.

        incoming_reaction (Reaction):
            Exchange reaction from the incoming model.

        policy (str):
            Bound merge policy.
    """

    valid_policies = {
        "first",
        "union",
        "intersection",
        "ignore",
    }

    if policy not in valid_policies:
        raise ValueError(
            f"Unknown exchange_bound_policy '{policy}'. "
            f"Valid policies are: {valid_policies}"
        )

    if policy == "first":
        return

    existing_lb = existing_reaction.getLowerBound()
    existing_ub = existing_reaction.getUpperBound()

    incoming_lb = incoming_reaction.getLowerBound()
    incoming_ub = incoming_reaction.getUpperBound()

    if policy == "union":

        new_lb = min(existing_lb, incoming_lb)
        new_ub = max(existing_ub, incoming_ub)

    elif policy == "intersection":

        new_lb = max(existing_lb, incoming_lb)
        new_ub = min(existing_ub, incoming_ub)

        if new_lb > new_ub:
            raise ValueError(
                "Incompatible exchange bounds detected "
                f"while merging reaction "
                f"'{existing_reaction.id}'.\n"
                f"Intersection produced invalid bounds:\n"
                f"lower_bound={new_lb} > upper_bound={new_ub}"
            )

    elif policy == "ignore":

        new_lb = -np.inf
        new_ub = np.inf

    existing_reaction.setLowerBound(new_lb)
    existing_reaction.setUpperBound(new_ub)


def merge_species(
    model: Model,
    new_id: str,
    duplicate_species: dict[str, int] | None = None,
    extracellular_compartments: list[str] = default_extracellular_compartments,
    rename_all_species: bool = False,
    merge_extracellular: bool = True,
) -> list[str]:
    """
    Rename and update species during community model merging.

    Species are renamed to avoid collisions between models
    and to associate them to model-specific compartments.

    Internal species are always processed.

    Extracellular species are:
    - shared unchanged if merge_extracellular=True;
    - compartment-namespaced if merge_extracellular=False.

    Args:
        model (Model):
            Source CBModel.

        new_id (str):
            Model-specific suffix appended to renamed species
            and compartments.

        duplicate_species (dict[str, int] | None, optional):
            Dictionary of duplicated species IDs across models.
            Used when rename_all_species=False. 
            Defaults to None.

        extracellular_compartments (list[str], optional):
            Compartment IDs considered extracellular.

        rename_all_species (bool, optional):
            If True, rename all species.
            If False, rename only duplicated species.
            Defaults to False.

        merge_extracellular (bool, optional):
            Whether extracellular species are shared across
            models.
            Defaults to False.

    Returns:
        list[str]:
            Extracellular species IDs encountered in the model.

    Notes:
        This function modifies the source model in-place:
        - species may be renamed;
        - species may be deleted/recreated;
        - reagent references may be rewired;
        - compartment IDs may be renamed.
    """

    if duplicate_species is None:
        duplicate_species = {}

    ls_species_ids = model.getSpeciesIds()

    external_metabolites: list[str] = []

    for species_id in ls_species_ids:

        species: Species = model.getSpecies(species_id)

        is_external = species.getCompartmentId() in extracellular_compartments

        # External species are skipped only when extracellular
        # metabolites are globally shared across models
        if not is_external or not merge_extracellular:

            if rename_all_species:

                species_id = copy_species_and_reagents(model, species, new_id)

            elif species_id in duplicate_species:

                species_id = copy_species_and_reagents(model, species, new_id)

            else:

                species.setCompartmentId(
                    create_new_id(species.compartment, new_id)
                )

        if is_external:
            external_metabolites.append(species_id)

    return external_metabolites


def copy_species_and_reagents(
    model: Model,
    species: Species,
    new_id: str,
) -> str:
    """
    Clone a species with a model-specific ID and update all
    associated reactions to use the new species.

    The original species is removed after all reagent references
    have been rewired.

    Args:
        model (Model):
            Source CBModel.

        species (Species):
            Species to duplicate and rename.

        new_id (str):
            Model-specific suffix appended to species and
            compartment IDs.

    Returns:
        str:
            ID of the newly created species.

    Notes:
        This function modifies the model in-place:
        - creates a new species;
        - rewires reaction reagents;
        - deletes the original species.
    """

    new_species_id = create_new_id(species.getId(), new_id)
    new_species = species.clone()
    new_species.setId(new_species_id)
    new_species.setName(species.name)
    new_species.setCompartmentId(create_new_id(species.compartment, new_id))

    model.addSpecies(new_species)

    # Rewire all reactions to use the renamed species
    for reaction_id in species.isReagentOf():
        reaction: Reaction = model.getReaction(reaction_id)
        reagent: Reagent = reaction.getReagentWithSpeciesRef(species.id)
        reaction.createReagent(new_species_id, reagent.coefficient)
        reaction.deleteReagentWithSpeciesRef(species.id)

    model.deleteSpecies(species.getId())

    return new_species_id


def copy_reaction(m_src: Model, m_targ: Model, rid, altrid=None):
    """
    Copy a reaction from a source model to a target model, if the required
    species exist in the target then they are mapped as reagents, otherwise
    new metabolites are added as boundary species.

     - *m_src* the source model
     - *m_targ* the target model
     - *rid* the reaction id to copy
     - *altrid* if the reaction name exists in the target,
     try use this one instead

    """
    out = {}
    if m_src.getReaction(rid) is None:
        print(
            'ERROR: reaction with id "{}" does not exist in source model'.format(
                rid
            )
        )
        out = None
    if (
        out is not None
        and m_targ.getReaction(rid) is not None
        and altrid is None
    ):
        print(
            'ERROR: reaction with id "{}" already exists and no alternative id was provided.'.format(
                rid
            )
        )

        out = None

    elif altrid is not None:
        if m_targ.getReaction(altrid) is not None:
            print(
                'ERROR: alternative reaction with id "{}" already exists in target model'.format(
                    rid
                )
            )
            out = None
    if out is None:
        return None
    old_reaction: Reaction = m_src.getReaction(rid)
    R: Reaction = old_reaction.clone()
    if altrid is not None:
        R.setId(altrid)
        for re in R.reagents:
            re.setId("{}_{}".format(altrid, re.getSpecies()))

    tSpecies = m_targ.getSpeciesIds()
    out["new_species"] = []
    out["existing_species"] = []
    out["reagents"] = []
    out["unmapped_species_reactions"] = {}
    for s in R.getSpeciesIds():
        if s not in tSpecies:
            S: Species = m_src.getSpecies(s).clone()
            # S.setBoundary()
            m_targ.addSpecies(S)
            out["new_species"].append(s)
            reag = m_src.getSpecies(s).isReagentOf()
            reag.remove(rid)
            reag.sort()
            out["unmapped_species_reactions"][s] = reag
        else:
            out["existing_species"].append(s)
        out["reagents"].append(s)

    if R.getId().startswith("R_EX") or R.getSBOterm() == "SBO:0000627":
        R.is_exchange = True

    m_targ.addReaction(R, create_default_bounds=False, silent=True)
    m_targ.setReactionBounds(
        R.id, old_reaction.getLowerBound(), old_reaction.getUpperBound()
    )

    del R

    return out


def create_new_id(old_id: str, new_id: str) -> str:
    """
    Function to build new id strings for models
    If empty string is provided we keep the original id for that model

    Args:
        old_id (str): The old identifier
        new_id (str): the new id, if "" than use old id otherwise
            append the old id with the new id

    Returns:
        str: _description_
    """
    if new_id == "":
        return f"{old_id}"

    return f"{old_id}_{new_id}"


def check_ids(new_ids: list[str], models: list[Model]) -> None:
    """
    Validate model-specific suffixes used during model merging.

    Args:
        new_ids (list[str]):
            Suffixes appended to model-specific entities.

        models (list[Model]):
            Models to combine.

    Raises:
        ValueError:
            If identifiers are not unique or may collide with
            existing reaction IDs.
    """

    if len(new_ids) != len(models):
        raise ValueError(
            "The number of new_ids must match the number of models."
        )

    if len(set(new_ids)) != len(new_ids):
        raise ValueError(
            "Model identifiers must be unique."
            )

    rids = {
        rid
        for model in models
        for rid in model.getReactionIds()
    }

    for rid in rids:
        for new_id in new_ids:

            if new_id and f"_{new_id}" in rid:
                raise ValueError(
                    "The provided model identifiers may collide "
                    "with existing reaction IDs. "
                    "Please provide alternative identifiers."
                )


def create_shared_pool(model: Model, cid: str = 'pool') -> None:
    """
    Create the shared extracellular pool compartment.

    Args:
        model (Model):
            Community CBModel.

        cid (str, optional):
            Shared pool compartment ID. Defaults to 'pool'.

    Raises:
        ValueError:
            If the pool compartment already exists.
    """

    if cid in model.getCompartmentIds():
        raise ValueError(f"A compartment with ID '{cid}' already exists")

    model.createCompartment(cid, "shared extracellular pool")


def connect_external_pool(
    combined_model: Model,
    new_ids: list[str],
    pool_compartment: str | None = None,
    extracellular_compartments: list[str] = default_extracellular_compartments,
) -> None:
    """
    Connect model-specific extracellular metabolites through
    a shared extracellular pool.

    Exchange reactions are converted into shuttle reactions
    linking organism-specific extracellular metabolites to
    shared pool metabolites.

    Shared pool exchange reactions are created automatically.

    Notes:

        - organism-specific exchange reactions retain their
          original bounds and represent organism transport
          capabilities;

        - pool exchange reactions represent shared medium
          availability and are initialized using default
          unconstrained bounds.

        Users may later constrain pool exchange reactions
        manually to simulate restricted environmental
        resources shared across organisms.

    Args:
        combined_model (Model):
            Community CBModel.

        new_ids (list[str]):
            Model-specific suffixes used during merging.

        extracellular_compartments (list[str], optional):
            Compartment IDs considered extracellular.
    """
    # setup pool compartment
    if pool_compartment is None:
        pool_compartment = "pool"

    if pool_compartment not in combined_model.getCompartmentIds():
        create_shared_pool(combined_model, pool_compartment)

    combined_model._pool_compartment = pool_compartment

    # Freeze exchange reaction list before modifying reactions
    for rid in list(combined_model.getExchangeReactionIds()):

        r = combined_model.getReaction(rid)
        stoich = r.getStoichiometry()

        if len(stoich) != 1:
            warnings.warn(
                (
                    f"Exchange reaction {rid} contains "
                    "more than one reagent. Skipping."
                ),
                UserWarning,
            )
            continue
        
        coef, met = stoich[0]

        if coef >= 0:
            warnings.warn(
                (
                    f"Exchange reaction {rid} has a "
                    "non-negative coefficient. Skipping."
                ),
                UserWarning,
            )
            continue

        # Recover canonical metabolite name by removing:
        # - model-specific suffixes;
        # - extracellular compartment suffixes.
        base_species_id = met

        for new_id in new_ids:
            if base_species_id.endswith(f"_{new_id}"):
                base_species_id = base_species_id.removesuffix(f"_{new_id}")
                break

        for cid in extracellular_compartments:
            if base_species_id.endswith(f"_{cid}"):
                base_species_id = base_species_id.removesuffix(f"_{cid}")
                break

        pool_id = f"{base_species_id}_pool"

        # Create shared pool metabolite
        if pool_id not in combined_model.getSpeciesIds():
            s = combined_model.getSpecies(met).clone()
            s.setId(pool_id)
            s.setCompartmentId(pool_compartment)
            combined_model.addSpecies(s)

        rexid = f"R_EX_{pool_id}"

        # Create exchange reaction for pool metabolite
        if rexid not in combined_model.getReactionIds():
            rex = Reaction(rexid)
            rex.createReagent(pool_id, -1)
            rex.is_exchange = True
            combined_model.addReaction(rex, create_default_bounds=True)

        # Convert organism-specific exchange into shuttle
        # reaction connected to the shared pool metabolite
        if pool_id not in r.getSpeciesIds():
            r.createReagent(pool_id, 1)

        r.is_exchange = False
