import matplotlib.pyplot as plt
import cbmpy
from cbmpy.CBModel import Model, Reaction
from dcFBA.ToyModels import model_a, model_b
from dcFBA.Models.CommunityModel import CommunityModel
from dcFBA.DynamicModels.DynamicJointFBA import DynamicJointFBA
from dcFBA.Models.Kinetics import KineticsStruct

model1: Model = cbmpy.loadModel(
    "models/bigg_models/e_coli_core.xml"
)  # load e_coli core
model2: Model = cbmpy.loadModel(
    "models/bigg_models/e_coli_core.xml"
)  # load e_coli core

model1.getReaction("R_GLCpts").setUpperBound(6)
model2.getReaction("R_GLCpts").setUpperBound(4)

community_model = CommunityModel(
    [model1, model2],
    ["R_BIOMASS_Ecoli_core_w_GAM", "R_BIOMASS_Ecoli_core_w_GAM"],
    ["modelA", "modelB"],
)


dj = DynamicJointFBA(community_model, [1, 1], {"M_glc__D_e": 20})
T, metabolites, biomasses, fluxes = dj.simulate(0.1)
index = 0

plt.plot(T, metabolites["M_glc__D_e"], color="blue", label="Metabolite s")
plt.xlabel("Time")
plt.ylabel("Concentration")
plt.legend()
plt.show()


plt.plot(T, biomasses["modelA"], color="blue", label="modelA")
plt.plot(T, biomasses["modelB"], color="orange", label="modelB")
plt.xlabel("Time")
plt.ylabel("Concentration")
plt.legend()
plt.show()
