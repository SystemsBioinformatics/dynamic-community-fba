import matplotlib.pyplot as plt
import cbmpy
from cbmpy.CBModel import Model
from DCFBA.Models import CommunityModel
from DCFBA.DynamicModels import DynamicParallelFBA


model1: Model = cbmpy.loadModel("models/bigg_models/e_coli_core.xml")
model2: Model = cbmpy.loadModel("models/bigg_models/strep_therm.xml")

for rid in model1.getReactionIds():
    if rid.startswith("R_EX"):
        reaction = model1.getReaction(rid)
        reaction.is_exchange = True

for rid in model2.getReactionIds():
    if rid.startswith("R_EX"):
        reaction = model2.getReaction(rid)
        reaction.is_exchange = True

model1.getReaction("R_GLCpts").setUpperBound(10)
model2.getReaction("R_GLCpts").setUpperBound(6)

parallel_fba = DynamicParallelFBA(
    [model2, model1],
    [1.0, 1.0],
    {"M_glc__D_e": 100, "M_gal_e": 0, "M_lcts_e": 100},
)

T, metabolites, biomasses = parallel_fba.simulate(0.1)

plt.plot(T, metabolites["M_glc__D_e"][:-1], color="blue", label="[Glucose]")
plt.plot(T, metabolites["M_lcts_e"][:-1], color="orange", label="[Lactose]")

plt.xlabel("Time")
plt.ylabel("Concentration")
plt.legend()
plt.show()

# plt.plot(T, biomasses["ecoli"], color="blue", label="modelA")
plt.plot(T, biomasses[model1.getId()], color="orange", label="ecoli")
plt.plot(T, biomasses[model2.getId()][:-1], color="blue", label="strep")

plt.xlabel("Time")
plt.ylabel("Concentration")
plt.legend()
plt.show()
