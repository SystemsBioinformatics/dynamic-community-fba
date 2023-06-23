import cbmpy
from .build_community_matrix import combine_models

from cbmpy.CBModel import Model, Reaction


def jointFBA(model: Model, biomass_reaction_ids: list[str]):
    model.createSpecies("X_c", False, "The community biomass")

    for bm_id in biomass_reaction_ids:
        reaction: Reaction = model.getReaction(bm_id)
        reaction.createReagent("X_c", 1)

    model.createReaction("Xcomm")
    out: Reaction = model.getReaction("Xcomm")
    out.is_exchange = True
    out.setUpperBound(1000)
    out.setLowerBound(0)
    out.createReagent("X_c", -1)

    model.createObjectiveFunction("Xcomm")

    model.setActiveObjective("Xcomm_objective")

    return model


model1 = cbmpy.loadModel("data/bigg_models/e_coli_core.xml")
# reaction: Reaction = model1.getReaction("R_GLCpts")
# reaction.setUpperBound(0)
model2 = cbmpy.loadModel("data/bigg_models/e_coli_core.xml")
model = combine_models([model1, model2], ["ecoli_1", "ecoli_2"])

sol = cbmpy.doFBA(model2)
FBAsol = model2.getSolutionVector(names=True)
FBAsol = dict(zip(FBAsol[1], FBAsol[0]))

print(sol)

sol2 = jointFBA(
    model,
    ["R_BIOMASS_Ecoli_core_w_GAM", "R_BIOMASS_Ecoli_core_w_GAM_ecoli_2"],
)


for key, value in FBAsol.items():
    if sol2[key] != value:
        print("combined: ", sol2[key])
        print(key, value)

    cbmpy.doFBA(model)
    FBAsol = model.getSolutionVector(names=True)
    FBAsol = dict(zip(FBAsol[1], FBAsol[0]))
