import cbmpy
import matplotlib.pyplot as plt
from cbmpy.CBModel import Model, Reaction
from DCFBA.ToyModels import model_a, model_b
from DCFBA.Models.CommunityModel import CommunityModel
from DCFBA.Models.Kinetics import Kinetics
from DCFBA.DynamicModels import EndPointFBA

m_a: Model = model_a.build_toy_model_fba_A()

m_a.getReaction("R_1").setUpperBound(10)
m_a.getReaction("R_4").setUpperBound(3)
m_a.getReaction("R_6").setUpperBound(1)

# Need to delete biomass model, since it is created in the endPoint model
reaction: Reaction = m_a.getReaction("R_BM_A").deleteReagentWithSpeciesRef(
    "BM_e_A"
)

# m_a.getReaction("R_import_B").setUpperBound(1)

m_a.getReaction("R_1").setLowerBound(0)
m_a.getReaction("R_4").setLowerBound(0)
m_a.getReaction("R_6").setLowerBound(0)

m_b: Model = model_b.build_toy_model_fba_B()
m_b.getReaction("R_1").setUpperBound(10)
m_b.getReaction("R_3").setUpperBound(1)
m_b.getReaction("R_5").setUpperBound(1)
# Need to delete biomass model, since it is created in the endPoint model
reaction: Reaction = m_b.getReaction("R_BM_B").deleteReagentWithSpeciesRef(
    "BM_e_B"
)


m_b.getReaction("R_1").setLowerBound(0)
m_b.getReaction("R_3").setLowerBound(0)
m_b.getReaction("R_5").setLowerBound(0)


kin = Kinetics(
    {
        "R_1_modelA": ["", 10, 10],
        "R_4_modelA": ["B_e", 5, 3],
        "R_6_modelA": ["B_e", 3, 1],
        # B
        "R_1_modelB": ["", 10, 10],
        "R_3_modelB": ["A_e", 2, 1],
    }
)


community_model = CommunityModel(
    [m_a, m_b], ["R_BM_A", "R_BM_B"], ["modelA", "modelB"]
)

# Are set by the model as often required by other GSMM
community_model.deleteReactionAndBounds("BM_e_A_exchange")
community_model.deleteReactionAndBounds("BM_e_B_exchange")

# community_model.getReaction("B_exchange").setLowerBound(-100)
# community_model.getReaction("A_exchange").setLowerBound(-100)


# community_model.getReaction("S_exchange").setLowerBound(-100)

n = 21
dj = EndPointFBA(
    community_model,
    n,
    {"modelA": 1, "modelB": 1},
    {"S_e": 100, "A_e": 0.0, "B_e": 0.0},
    0.1,
)

#    {"S_e": 100, "A_e": 0.0, "B_e": 0.0},


solution = dj.simulate()


FBAsol = dj.m_model.getSolutionVector(names=True)
FBAsol = dict(zip(FBAsol[1], FBAsol[0]))
# Get biomassess over time
# print(solution)
# print(FBAsol)


biomassA = [1]
biomassB = [1]

for key, value in FBAsol.items():
    if key.startswith("BM_modelA_time"):
        biomassA.append(value)
    if key.startswith("BM_modelB_time"):
        biomassB.append(value)


plt.plot(list(range(0, n)), biomassA, color="blue", label="Model a")
plt.plot(list(range(0, n)), biomassB, color="orange", label="model B")

plt.xlabel("Time")
plt.ylabel("Concentration")
plt.legend()
plt.show()


# TODO uncomment for ecolit model

# model1 = cbmpy.loadModel("data/bigg_models/e_coli_core.xml")
# model1.getReaction("R_EX_glc__D_e").setLowerBound(-1000)
# model1.getReaction("R_GLCpts").setUpperBound(10)

# model2: Model = model1.clone()
# model2.setId("model2")

# cm = CommunityModel(
#     [model1],
#     ["R_BIOMASS_Ecoli_core_w_GAM"],
#     ["ecoli_1"],
# )

# ep = EndPointFBA(cm, 10, {"ecoli_1": 1}, 1)

# # print(ep.m_model.getReaction("R_GLCpts_ecoli_1_time00").getSpeciesIds())

# for e in ep.m_model.getExchangeReactionIds():
#     print(e)
#     reaction: Reaction = ep.m_model.getReaction(e)
#     print(reaction.getSpeciesIds())

# solution = ep.simulate()

# FBAsol = ep.m_model.getSolutionVector(names=True)
# FBAsol = dict(zip(FBAsol[1], FBAsol[0]))

# print(FBAsol)
