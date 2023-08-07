import matplotlib.pyplot as plt
import cbmpy
from cbmpy.CBModel import Model, Species, Reaction
from DCFBA.Models import CommunityModel
from DCFBA.DynamicModels import DynamicSingleFBA


model1: Model = cbmpy.loadModel("models/bigg_models/e_coli_core.xml")
model2: Model = cbmpy.loadModel("models/bigg_models/strep_therm.xml")


solution = cbmpy.doFBA(model2)
FBAsol = model2.getSolutionVector(names=True)
FBAsol = dict(zip(FBAsol[1], FBAsol[0]))

for key, value in FBAsol.keys():
    if key.startswith("R_EX"):
        print(key, value)

raise Exception
model1.getReaction("R_GLCpts").setUpperBound(10)


model2.getReaction("R_GLCt2").setUpperBound(1)
# model2.getReaction("R_EX_LCTSGALex").setLowerBound(0)
# R_EX_gal_e
# model2.getReaction("R_EX_LCTSGALex").setLowerBound(-10)
# print(model2.getReaction("R_EX_LCTSGALex").getUpperBound())

biomass_reaction_model_1: str = "R_BIOMASS_Ecoli_core_w_GAM"
biomass_reaction_model_2: str = "R_biomass_STR"

{"M_glc__D_e": 100, "M_lcts_e": 100}
dj = DynamicSingleFBA(
    model2, biomass_reaction_model_2, 0.1, {"M_glc__D_e": 1000, "M_gal_e": 25}
)

T, metabolites, biomasses, fluxes = dj.simulate(0.1)

plt.plot(T, metabolites["M_glc__D_e"], color="blue", label="[Glucose]")
plt.xlabel("Time")
plt.ylabel("Concentration")
plt.legend()
plt.show()

# plt.plot(T, biomasses["ecoli"], color="blue", label="modelA")
plt.plot(T, biomasses[""], color="orange", label="modelB")
plt.xlabel("Time")
plt.ylabel("Concentration")
plt.legend()
plt.show()
