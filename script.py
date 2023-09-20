import matplotlib.pyplot as plt
import cbmpy
from cbmpy.CBModel import Model, Reaction
from dcFBA.Models.CommunityModel import CommunityModel
from dcFBA.DynamicModels import EndPointFBA
import numpy as np
from dcFBA.Helpers.PlotsEndPointFBA import plot_biomasses, plot_metabolites

model1: Model = cbmpy.loadModel(
    "models/bigg_models/e_coli_core.xml"
)  # load e_coli core
model2: Model = cbmpy.loadModel(
    "models/bigg_models/e_coli_core.xml"
)  # load e_coli core

model1.getReaction("R_GLCpts").setUpperBound(6)
model2.getReaction("R_GLCpts").setUpperBound(4)

# model1.getReaction("R_PGK").setLowerBound(cbmpy.NINF)
# model2.getReaction("R_PGK").setLowerBound(cbmpy.NINF)

# model1.getReaction("R_PGM").setUpperBound(cbmpy.INF)
# model2.getReaction("R_PGM").setUpperBound(cbmpy.INF)

community_model = CommunityModel(
    [model1, model2],
    ["R_BIOMASS_Ecoli_core_w_GAM", "R_BIOMASS_Ecoli_core_w_GAM"],
    ["modelA", "modelB"],
)


ep = EndPointFBA(
    community_model,
    3,
    {"modelA": 1.0, "modelB": 1.0},
    {"M_glc__D_e": 1000},
    dt=0.1,
)
print("========")
print(ep.m_model.getReaction("R_EX_glc__D_e").getLowerBound())
print("========")

solution = ep.simulate()

FBAsol = ep.m_model.getSolutionVector(names=True)
FBAsol = dict(zip(FBAsol[1], FBAsol[0]))

print(FBAsol["M_glc__D_e_time0_time1"])
print(FBAsol["M_glc__D_e_time1_time2"])
print(FBAsol["M_glc__D_e_exchange_final"])

print(solution)
plot_biomasses(ep)
plot_metabolites(ep, {"M_glc__D_e": 1000})
