import matplotlib.pyplot as plt
from cbmpy.CBModel import Model, Reaction
from dcFBA.ToyModels import model_a, model_b
from dcFBA.Models.CommunityModel import CommunityModel
from dcFBA.Models.Kinetics import KineticsStruct
from dcFBA.DynamicModels import EndPointFBA
from dcFBA.Helpers.PlotsEndPointFBA import plot_biomasses, plot_metabolites
from dcFBA.Helpers.OptimalSearch import (
    balance_search_clean,
    balanced_search_quick,
)
import math

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


kin = KineticsStruct(
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

# n, high = search(
#     community_model,
#     {"modelA": 1.0, "modelB": 2.0},
#     {"S_e": 100, "A_e": 0.0, "B_e": 0.0},
#     0.1,
# )


# print(
#     balance_search(
#         community_model,
#         21,
#         {"S_e": 100, "A_e": 0.0, "B_e": 0.0},
#         0.1,
#         3,
#         12.7777,
#         epsilon=0.0001,
#     )
# )

n = 21
ep = EndPointFBA(
    community_model,
    n,
    {"modelA": 1, "modelB": 2},
    {"S_e": 100, "A_e": 0.0, "B_e": 0.0},
    1,
)

solution1 = ep.simulate()
n = 21
ep = EndPointFBA(
    community_model,
    n,
    {"modelA": 1, "modelB": 2},
    {"S_e": 100, "A_e": 0.0, "B_e": 0.0},
    0.1,
)

e = balance_search_clean(
    community_model,
    n,
    {"S_e": 100, "A_e": 0.0, "B_e": 0.0},
    0.1,
    3,
    12.777,
    0.001,
)
input(e)
e = balanced_search_quick(ep, 3, 12.777, 0.001)


n = 21
ep = EndPointFBA(
    community_model,
    n,
    {"modelA": 1, "modelB": 2},
    {"S_e": 100, "A_e": 0.0, "B_e": 0.0},
    0.1,
)

ep.balanced_growth(3, 12.777 * e)


solution2 = ep.simulate()
print(solution1)
print(solution2)

print(e)
plot_biomasses(ep)

# plot_metabolites(ep, {"S_e": 100, "A_e": 0.0, "B_e": 0.0})
#  0.9957275390625
# 0.9921875
# print(
#     balance_search(
#         community_model,
#         n,
#         {"S_e": 100, "A_e": 0.0, "B_e": 0.0},
#         0.1,
#         3,
#         12.7777,
#     )
# )


# if not math.isnan(solution):
# FBAsol = ep.m_model.getSolutionVector(names=True)
# FBAsol = dict(zip(FBAsol[1], FBAsol[0]))
# print(FBAsol["Phi_modelA"])
# print(FBAsol["Phi_modelB"])
# print(FBAsol["BM_modelA_exchange_final"])
# print(FBAsol["BM_modelB_exchange_final"])
# print("X_comm:" + str(FBAsol["X_comm"]))
# print()
# print(
#     FBAsol["BM_modelA_exchange_final"]
#     / (
#         FBAsol["BM_modelB_exchange_final"]
#         + FBAsol["BM_modelA_exchange_final"]
#     )
# )

# plot_biomasses(ep)
# plot_metabolites(ep, {"S_e": 100, "A_e": 0.0, "B_e": 0.0})
