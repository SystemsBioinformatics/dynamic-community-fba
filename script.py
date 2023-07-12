import cbmpy
import cobra
from cobra.io import save_json_model, read_sbml_model

from DCFBA.DynamicModels import EndPointFBA
from DCFBA.Models import CommunityModel
from DCFBA.ToyModels import model_a, model_b, model_c
from DCFBA.DynamicModels import DynamicJointFBA

# model_a = model_a.build_model_A()
# model_b = model_b.build_model_B()

# community_model: CommunityModel = CommunityModel(
#     [model_a, model_b], ["R_BM_A", "R_BM_B"]
# )
# efb = EndPointFBA()


# model_c_1 = model_c.build_model_C()
# model_c_1.getReaction("R_1").setUpperBound(10)

# community_model = CommunityModel([model_c_1], ["R_BM_C"])
# efb = EndPointFBA(community_model, 3, {"Organism_C": 1}, 1)

# efb.m_model.createReaction(
#     "final_biomass", "FInal time point to biomass", reversible=False
# )

# reaction: cbmpy.CBModel.Reaction = efb.m_model.getReaction("final_biomass")
# reaction.createReagent("BM_e_C_time3", -1)
# reaction.setUpperBound(1e10)
# reaction.setLowerBound(-1)
# reaction.is_exchange = False
# efb.m_model.createObjectiveFunction("final_biomass")

# reaction = efb.m_model.getReaction("BM_e_C_exchange")
# reaction.setLowerBound(-1)

# r: cbmpy.CBModel.Reaction
# for r in efb.m_model.reactions:
#     print("----------")
#     print(r.getId())
#     print("-")
#     for id in r.getSpeciesIds():
#         print(r.getReagentWithSpeciesRef(id).coefficient, id)
#     print([r.getLowerBound(), r.getUpperBound()])
#     print()
#     print("----------")

# cbmpy.doFBA(efb.m_model)
# FBAsol = efb.m_model.getSolutionVector(names=True)
# FBAsol = dict(zip(FBAsol[1], FBAsol[0]))

# print(FBAsol)


# TODO Ecoli core below

textbook = cbmpy.loadModel("data/bigg_models/e_coli_core.xml")
textbook.getReaction("R_GLCpts").setUpperBound(10)
textbook.getReaction("R_EX_glc__D_e").setLowerBound(-1000)

community_model = CommunityModel(
    [textbook, textbook.clone()],
    ["R_BIOMASS_Ecoli_core_w_GAM", "R_BIOMASS_Ecoli_core_w_GAM"],
    ["ecoli_1", "ecoli_2"],
)

efb = EndPointFBA(community_model, 2, {"ecoli_1": 1, "ecoli_2": 1}, 1)


cbmpy.doFBA(efb.m_model)

FBAsol = efb.m_model.getSolutionVector(names=True)
FBAsol = dict(zip(FBAsol[1], FBAsol[0]))
print(FBAsol)
# for reaction in efb.m_model.reactions:
#     print(reaction.getId())
