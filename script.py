import matplotlib.pyplot as plt
from cbmpy.CBModel import Model, Reaction, Objective
from dcFBA.ToyModels import model_a, model_b
from dcFBA.Models.CommunityModel import CommunityModel
from dcFBA.Models.Kinetics import KineticsStruct
from dcFBA.DynamicModels import EndPointFBA
from dcFBA.Helpers.PlotsEndPointFBA import plot_biomasses, plot_metabolites

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


community_model = CommunityModel(
    [m_a, m_b], ["R_BM_A", "R_BM_B"], ["modelA", "modelB"]
)

# Are set by the model as often required by other GSMM
community_model.deleteReactionAndBounds("BM_e_A_exchange")
community_model.deleteReactionAndBounds("BM_e_B_exchange")

print(len(community_model.getReactionIds()))
print(len(community_model.getExchangeReactionIds()))

print(len(community_model.getSpeciesIds()))

n = 4
ep = EndPointFBA(
    community_model,
    n,
    {"modelA": 1.0, "modelB": 2.0},
    {"S_e": 100, "A_e": 0.0, "B_e": 0.0},
    0.1,
)
print(len(ep.m_model.getReactionIds()))
print(len(community_model.getExchangeReactionIds()))

print(len(ep.m_model.getSpeciesIds()))

raise Exception
# ep.balanced_growth(3.0, 12.72)
# ep.set_qp(12.72, epsilon=0.01)
# obj: Objective = ep.m_model.getActiveObjective()
# i = 0

# solution = ep.simulate()
# FBAsol = ep.m_model.getSolutionVector(names=True)
# FBAsol = dict(zip(FBAsol[1], FBAsol[0]))

# print(FBAsol["BM_modelA_exchange_final"])
# print(FBAsol["BM_modelB_exchange_final"])
# print(FBAsol["X_comm"])

# plot_biomasses(dj)
# plot_metabolites(dj, {"S_e": 100, "A_e": 0.0, "B_e": 0.0})
