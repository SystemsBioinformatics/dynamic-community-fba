from cbmpy.CBModel import Model, Reaction
from DCFBA.ToyModels import model_a, model_b, model_d
from DCFBA.Models import CommunityModel
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


m_b.getReaction("R_1").setLowerBound(0)
m_b.getReaction("R_3").setLowerBound(0)
m_b.getReaction("R_5").setLowerBound(0)

# Need to delete biomass model, since it is created in the endPoint model
reaction: Reaction = m_b.getReaction("R_BM_B").deleteReagentWithSpeciesRef(
    "BM_e_B"
)

m_d: Model = model_d.build_toy_model_fba_D()
m_d.getReaction("R_1").setUpperBound(10)
m_d.getReaction("R_3").setUpperBound(1)
m_d.getReaction("R_5").setUpperBound(1)

m_d.getReaction("R_1").setLowerBound(0)
m_d.getReaction("R_3").setLowerBound(0)
m_d.getReaction("R_5").setLowerBound(0)

# Need to delete biomass model, since it is created in the endPoint model
reaction: Reaction = m_d.getReaction("R_BM_D").deleteReagentWithSpeciesRef(
    "BM_e_D"
)


community_model = CommunityModel(
    [m_a, m_b, m_d],
    ["R_BM_A", "R_BM_B", "R_BM_D"],
    ["modelA", "modelB", "modelD"],
)

# Are set by the model as often required by other GSMM
community_model.deleteReactionAndBounds("BM_e_A_exchange")
community_model.deleteReactionAndBounds("BM_e_B_exchange")
community_model.deleteReactionAndBounds("BM_e_D_exchange")

n = 80
width = len(str(n))
# TODO times should not hold the under score
times = [f"time{i:0{width}d}" for i in range(n)]

m = build_time_model(community_model, times)


dj = EndPointFBA(
    community_model,
    n,
    {"modelA": 1.0, "modelB": 1.0, "modelD": 1.0},
    {"S_e": 100, "A_e": 0.0, "B_e": 0.0},
    0.1,
)
