from cbmpy.CBModel import Model, Species, Reaction, Reagent, Compartment
import cbmpy
import pandas as pd


def load_n_models(files: list[str]) -> list[Model]:
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


def combine_models(models: list[Model]) -> Model:
    # First just rename all species/reagants

    combined_model = Model("Combined_model")
    combined_model.createCompartment("e", "extracellular space")
    duplicate_species = create_duplicate_species_dict(models)
    for model in models:
        fix_duplicate_species(duplicate_species, model)
        merge_compartments(model, combined_model)
        merge_reactions(model, combined_model)

    return combined_model


def create_duplicate_species_dict(models: list[Model]) -> dict[str, int]:
    combined_species_list = [
        item for sublist in models for item in sublist.species
    ]

    species_dict = {}
    for m in combined_species_list:
        if m.id in species_dict.keys():
            species_dict[m.id] += 1
        else:
            species_dict[m.id] = 1

    # Find which species are in more than 1 model
    return {k: v for k, v in species_dict.items() if v >= 2}


# Reanem compartment ids to be organism specific?
def merge_compartments(model: Model, combined_model: Model):
    compartment: Compartment

    for compartment in model.compartments:
        if compartment.id != "e":
            combined_model.createCompartment(
                f"{compartment.id}_{model.id}",
                f"Compartmnet {compartment.name} of model: {model.name}",
                size=compartment.size,
                dimensions=compartment.dimensions,
            )


def merge_reactions(model: Model, combined_model: Model):
    reaction: Reaction
    for reaction in model.reactions:
        # Shoudl be fixed in the package
        res = copyReaction(model, combined_model, reaction.id)
        # If the response is None the reaction was probably already in the model.
        # Check if it is an exchange reaction or if it is specfic to the cellular space
        # of the organisms, and fix accordingly
        if res is None and (not reaction.id.startswith("R_EX")):
            copyReaction(
                model,
                combined_model,
                reaction.id,
                altrid=f"{reaction.id}_{model.id}",
            )
        else:
            continue


def fix_duplicate_species(duplicate_species: dict[str, int], model: Model):
    """If a species occurs in two or more models we append the id of that
    species with the model id, such that we make a distinguishing between the
    two species/metabolites in the organisms, if the species on the other
    hand is present in the extracellular space we can leave the species

    After renaming the species all the linked reagents also need to be changed
    see handle_duplicate_species_reagents

    Args:
        models (list[Model]): The list of the models that need to be in the
                              combined model
        combined_model (Model): The new combined model

    Returns:
        Model: The new combined model
    """

    # Loop over the duplicate species, if the species is in the old_model
    # and if a species is not in the external compartment we have to change it
    # such that we know from which organism it is
    # Next we need to change all reactions in which the old species
    # was a reagent to link to the new wrapped species.

    for species_id in duplicate_species.keys():
        if species_id in model.getSpeciesIds():
            species: Species = model.getSpecies(species_id)
            if species.compartment != "e":
                handle_duplicate_species_reagents(model, species)
                model.deleteSpecies(species_id)
            # TODO what if species has two exchange for the same metabolite
            # TODO reaction with different ids?
        #  else:
        #     for model in models:


def handle_duplicate_species_reagents(model: Model, species: Species) -> Model:
    """_summary_

    Args:
        model (Model): _description_
        combined_model (Model): _description_
        species (Species): _description_

    Returns:
        Model: _description_
    """
    # If the species exists in one of the base models, add new
    # instance of that species in the old model with a new id
    # Afterwards update all the links of the reactions of the model
    # to this new species
    new_species_id = f"{species.id}_{model.id}"
    new_species = species.clone()
    new_species.setId(new_species_id)
    new_species.setName(f"{species.name} of {model.name}")
    new_species.setCompartmentId(f"{species.compartment}_{model.id}")

    model.addSpecies(new_species)

    for reaction_id in species.isReagentOf():
        reaction: Reaction = model.getReaction(reaction_id)
        reagent: Reagent = reaction.getReagentWithSpeciesRef(species.id)
        reaction.createReagent(new_species_id, reagent.coefficient)
        reaction.deleteReagentWithSpeciesRef(species.id)


def copyReaction(m_src, m_targ, rid, altrid=None):
    """
    Copy a reaction from a source model to a target model, if the required species exist in the target
    then they are mapped as reagents, otherwise new metabolites are added as boundary species.

     - *m_src* the source model
     - *m_targ* the target model
     - *rid* the reaction id to copy
     - *altrid* if the reaction name exists in the target, try use this one instead

    """
    out = {}
    targ_exists = False
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
        else:
            targ_exists = True
    if out is None:
        return None

    R = m_src.getReaction(rid).clone()

    if targ_exists and altrid is not None:
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
    m_targ.addReaction(R, create_default_bounds=True, silent=True)
    del R

    return out


def main():
    ls_models_combined = [
        "../Models/toy_organism_A.xml",
        "../Models/toy_organism_B.xml",
        "../Models/toy_organism_A.txt",
    ]

    models = load_n_models(ls_models_combined)

    combined_model = combine_models(models)
    combined_model.createReaction(
        "BM_A_and_B", "Biomass_A_and_B_combined", False
    )

    # I had double exchange reactions for every metabolite
    # because the reaction ids were different but they exchanged the same thing
    combined_model.deleteReactionAndBounds("R_EX_BM_ext_A")
    combined_model.deleteReactionAndBounds("R_EX_BM_ext_B")
    combined_model.deleteReactionAndBounds("R_EX_S_ext_A")
    combined_model.deleteReactionAndBounds("R_EX_B_ext_A")
    combined_model.deleteReactionAndBounds("R_EX_A_ext_B")

    combined_model.getReaction("BM_A_and_B").is_exchange = True

    combined_model.createReactionReagent("BM_A_and_B", "M_BM_ext_A", -1)
    combined_model.createReactionReagent("BM_A_and_B", "M_BM_ext_B", -3)
    combined_model.setReactionLowerBound("BM_A_and_B", 0)

    combined_model.createObjectiveFunction("BM_A_and_B")
    print(combined_model.getActiveObjectiveReactionIds())
    reaction: Reaction
    for reaction in combined_model.reactions:
        reaction.setLowerBound(-1000.0)
        reaction.setUpperBound(1000.0)
    print(cbmpy.doFBA(combined_model))
    FBAsol = combined_model.getSolutionVector(names=True)
    FBAsol = pd.DataFrame(zip(FBAsol[1], FBAsol[0]), columns=["ID", "flux"])

    print(FBAsol)
    # cbmpy.saveModel(combined_model, "../Models/test_community.xml")
    # cbmpy.writeCOBRASBML(combined_model, "../Models/test_community_cobra.xml")


# main()


list_real_models = ["../Models/strep.xml", "../Models/iML1515.xml"]
models = load_n_models(list_real_models)

combined_model = combine_models(models)
print(combined_model.compartments)
print(combined_model.species[3243].compartment)

# print(combined_model.reactions)
