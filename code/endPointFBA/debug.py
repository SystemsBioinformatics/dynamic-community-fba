import pandas as pd
import cbmpy
from cbmpy.CBModel import Reaction
from build_community_matrix import (
    load_n_models,
    combine_models,
    create_duplicate_species_dict,
)
import sys

# sys.path is a list of absolute path strings
sys.path.append("code/tests/helpers")  # <-- relative path
import model_a
import model_b


def main():
    # ls_models_combined = [
    #     "data/toy_model/toy_organism_A.xml",
    #     "data/toy_model/toy_organism_B.xml",
    #     "data/toy_model/toy_organism_A.txt",
    # ]

    # models = load_n_models(ls_models_combined)
    models = [model_a.build_model_A(), model_b.build_model_B()]

    combined_model = combine_models(models)
    combined_model.createReaction(
        "BM_A_and_B", "Biomass_A_and_B_combined", False
    )

    # print(combined_model.getSpeciesIds())

    # I had double exchange reactions for every metabolite
    # because the reaction ids were different but they exchanged the same thing

    # combined_model.getReaction("BM_A_and_B").is_exchange = True

    # combined_model.createReactionReagent("BM_A_and_B", "M_BM_ext_A", -1)
    # combined_model.createReactionReagent("BM_A_and_B", "M_BM_ext_B", -3)
    # combined_model.setReactionLowerBound("BM_A_and_B", 0)


# main()


list_real_models = [
    "data/bigg_models/strep.xml",
    "data/bigg_models/iML1515.xml",
]

models = load_n_models(list_real_models)
print(models[0].getExchangeReactionIds())
print(models[1].getExchangeReactionIds())

combined_model = combine_models(models)
print(combined_model.getCompartmentIds())
print(combined_model.getReactionIds())
print(combined_model.getSpeciesIds())


print(len(models[0].getReactionIds()))
print(len(models[1].getReactionIds()))
print(len(combined_model.getReactionIds()))

print(len(models[0].getSpeciesIds()))
print(len(models[1].getSpeciesIds()))
print(len(combined_model.getSpeciesIds()))


print("\n---------------------\n")

models = load_n_models(list_real_models)
print(create_duplicate_species_dict(models))
print(len(create_duplicate_species_dict(models)))

for model in models:
    for species_id in model.getSpeciesIds():
        if (
            species_id not in combined_model.getSpeciesIds()
            and f"{species_id}_{model.id}"
            not in combined_model.getSpeciesIds()
        ):
            print(model.id)
            print(species_id)

print(combined_model.getSpecies("M_mnl1p_c_iRZ476").id)
