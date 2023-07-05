import cbmpy
import matplotlib.pyplot as plt
from endPointFBA.CommunityModel import CommunityModel
from endPointFBA.DynamicJointFBA import DynamicJointFBA
from toy_models import model_c
from cbmpy.CBModel import Model, Reaction
import numpy as np

model1: Model = cbmpy.loadModel("data/bigg_models/e_coli_core.xml")
model_2 = model1.clone()
print(model1.getReaction("R_ATPM").setLowerBound(0))

model2 = model_c.build_model_C()
model2.getReaction("R_1").setUpperBound(1)
model2.getReaction("R_1").setLowerBound(-1)

model2.getReaction("S_exchange").setLowerBound(-10000000)

# combined_model = CommunityModel(
#     [model2],
#     ["2"],
#     ["R_BM_C"],
#     "R_BM_C_2",
# )

# km_1: KineticModel = KineticModel(combined_model, {})

# dynamic_fba = DynamicJointFBA(km_1, [1], {"S_e": 100})

# dynamic_fba.simulate(1)

model3 = model_c.build_model_C()
model3.getReaction("R_1").setUpperBound(10)
model3.getReaction("R_1").setLowerBound(0)

model3.getReaction("S_exchange").setLowerBound(-100)

combined_model = CommunityModel(
    [model1],
    ["R_BIOMASS_Ecoli_core_w_GAM"],
    ["ecoli_1"],
)

reaction: Reaction = combined_model.getReaction("R_EX_glc__D_e")
reaction.setLowerBound(-100000)

reaction = combined_model.getReaction("R_GLCpts_ecoli_1")
reaction.setUpperBound(10)


dynamic_fba = DynamicJointFBA(
    combined_model,
    [0.1],
    {"M_glc__D_e": 10},
    {},
)

ts, metabolites, biomasses = dynamic_fba.simulate(0.1)

print(ts)
print(metabolites["M_glc__D_e"])
ax = plt.subplot(111)
ax.plot(ts, biomasses["ecoli_1"])
ax2 = plt.twinx(ax)
ax2.plot(ts, metabolites["M_glc__D_e"], color="r")

ax.set_ylabel("Biomass ecoli 1", color="b")
ax2.set_ylabel("X_comm", color="r")

# ax3 = plt.twinx(ax)
# ax3.plot(ts, biomasses["ecoli_2"], color="y")
# ax3.set_ylabel("Biomass ecoli 2", color="y")

plt.show()
