"""Helper functions for building the CommunityModel 
"""

from cbmpy.CBModel import Model, Species, Reaction, Reagent, Compartment
import cbmpy


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
    new_ids: list[str] = [],
    modelId="combined_model",
    objective_function: str = "",
) -> Model:
    """
    Combine N models into a single model.
    For further information see Models.CombinedModel

    Args:
        models (list[Model]): A list of CBModels to combine.
        new_ids (list(str), optional): A list of user specified ids that
        should be used to identify the different species used
        objective_function (str, optional): String of the rid that should be
        set as objective function

    Returns:
        Model: The combined CBModel.

    """
    combined_model = Model(modelId)
    combined_model.createCompartment("e", "extracellular space")
    duplicate_species = create_duplicate_species_dict(models)

    if (len(new_ids) > 0) and (len(new_ids) < len(models)):
        raise Exception("Too few ids were provided")

    if len(new_ids) == 0:
        for model in models:
            new_ids.append(model.id)

    for i in range(0, len(models)):
        model = models[i]
        new_id = new_ids[i]
        merge_compartments(model, combined_model, new_id)
        merge_species(duplicate_species, model, new_id)
        merge_reactions(model, combined_model, new_id)

    if len(objective_function) > 0:
        combined_model.createObjectiveFunction(objective_function)

    return combined_model


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


def merge_compartments(model: Model, combined_model: Model, new_id):
    """
    Merge compartments from a model into a combined model.

    Args:
        model (Model): The source CBModel.
        combined_model (Model): The combined CBModel.

    Returns:
        None
    """

    compartment: Compartment

    for compartment in model.compartments:
        if compartment.id != "e":
            combined_model.createCompartment(
                create_new_id(compartment.id, new_id),
                f"Compartmnet {compartment.name}",
                size=compartment.size,
                dimensions=compartment.dimensions,
            )


def merge_reactions(model: Model, combined_model: Model, new_id: str) -> None:
    """
    Merge reactions from a model into the combined model.

    Args:
        model (Model): The source CBModel.
        combined_model (Model): The combined CBModel.

    Returns:
        None

    """

    exchange_reactions = model.getExchangeReactionIds()
    for reaction_id in model.getReactionIds():
        is_exchange_reaction: bool = False

        # TODO discuss this
        if reaction_id.startswith("R_EX"):
            is_exchange_reaction = True
        if reaction_id in exchange_reactions:
            is_exchange_reaction = True

        reaction: Reaction = model.getReaction(reaction_id)
        # if it is an exchange reaction and it is not yet in the model add it
        if (
            is_exchange_reaction
            and reaction_id not in combined_model.getExchangeReactionIds()
        ):
            copyReaction(
                model,
                combined_model,
                reaction.id,
            )
        # Otherwise add the reaction and add the identifier, we give all
        # reactions a new identifier for consistency
        elif not is_exchange_reaction:
            copyReaction(
                model,
                combined_model,
                reaction.id,
                altrid=create_new_id(reaction_id, new_id),
            )


def merge_species(
    duplicate_species: dict[str, int], model: Model, new_id
) -> None:
    """
    Merge species from a model into a combined model.
    Rules based:
    1) If a species occurs in two or more models in the cytosolic space
    we append the id of that species with the model id, such that we
    distinguish between the two species/metabolites in the organisms
    2) If the species on the other lives in the extracellular space
    we don't merge the species, the two model now use the same extracellular
    space species
    3) If the species of a model is not a duplicate we still change it's
    compartment for clarity. Now the user know where the species is located

    Args:
        models (list[Model]): The list of the models that need to be in the
                              combined model

        combined_model (Model): The new combined model

    Returns:
        None
    """

    ls_species_ids = model.getSpeciesIds()
    for species_id in ls_species_ids:
        species: Species = model.getSpecies(species_id)
        if species.compartment != "e":
            if species_id in duplicate_species.keys():
                copy_species_and_reagents(model, species, new_id)
                model.deleteSpecies(species_id)
            else:
                species.setCompartmentId(
                    create_new_id(species.compartment, new_id)
                )


def copy_species_and_reagents(model: Model, species: Species, new_id):
    """
    Handle duplicate species reagents when merging models.
    If the species exists in one of the base models, add new instance of that
    species in the old model with a new id. Afterwards update
    all the reactions in which the old_species occurs such that the
    reactions uses the new_species, and delete the old reagent

    Args:
        model (Model): The CBModel.

        species (Species): The species to handle.

    Returns:
        None

    """

    new_species_id = create_new_id(species.getId(), new_id)
    new_species = species.clone()
    new_species.setId(new_species_id)
    new_species.setName(species.name)
    new_species.setCompartmentId(create_new_id(species.compartment, new_id))

    model.addSpecies(new_species)

    for reaction_id in species.isReagentOf():
        reaction: Reaction = model.getReaction(reaction_id)
        reagent: Reagent = reaction.getReagentWithSpeciesRef(species.id)
        reaction.createReagent(new_species_id, reagent.coefficient)
        reaction.deleteReagentWithSpeciesRef(species.id)


def copyReaction(m_src: Model, m_targ: Model, rid, altrid=None):
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
        print(
            'INFO: reaction with id "{}" exists in target model trying alternate id'.format(
                rid
            )
        )

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
    # TODO ask Brett
    if R.getId().startswith("R_EX"):
        R.is_exchange = True
    m_targ.addReaction(R, create_default_bounds=True, silent=True)
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
