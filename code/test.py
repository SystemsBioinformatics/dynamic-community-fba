import cbmpy
import matplotlib.pyplot as plt
from endPointFBA.CommunityModel import CommunityModel
from endPointFBA.jointFBA import create_joint_fba_model
from endPointFBA.KineticModel import KineticModel
from endPointFBA.dynamic_fba import DynamicJointFBA
from toy_models import model_c
from cbmpy.CBModel import Model, Reaction
import numpy as np

model1: Model = cbmpy.loadModel("data/bigg_models/e_coli_core.xml")
for reaction_id in model1.getReactionIds():
    if reaction_id.startswith("R_EX"):
        reaction = model1.getReaction(reaction_id)
        reaction.is_exchange = True

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
    ["ecoli_1"],
    ["R_BIOMASS_Ecoli_core_w_GAM"],
)

combined_model.getReaction("R_EX_mal__L_e").setLowerBound(-1000000)

combined_model.createObjectiveFunction("R_BIOMASS_Ecoli_core_w_GAM_ecoli_1")
reaction = combined_model.getReaction("R_EX_glc__D_e")
reaction.setLowerBound(-10000000)

combined_model.getReaction("R_ATPM_ecoli_1").setLowerBound(0)
reaction = combined_model.getReaction("R_GLCpts_ecoli_1")
reaction.setUpperBound(10)

combined_model.getReaction("R_MALt2_2_ecoli_1").setUpperBound(1)

print(reaction.getEquation())

# combined_model = CommunityModel(
#     [model2, model3], ["A1", "B2"], ["R_BM_C", "R_BM_C"]
# )

# combined_model.createReaction("BM_test", "testing biomass function")

# joint_model = create_joint_fba_model(
#     combined_model, ["R_BM_C_A1", "R_BM_C_B2"]
# )

# ts, substrates, biomasses = dynamic_joint_FBA(
#     km_1, 0.1, [0.1, 0.1], {"S_e": 10}
# )

# "R_GLCpts_ecoli_1": [10, 5]
km_1: KineticModel = KineticModel(combined_model, {})

dynamic_fba = DynamicJointFBA(
    km_1, [0.1], {"M_glc__D_e": 0, "M_mal__L_e": 100}
)
ts, metabolites, biomasses = dynamic_fba.simulate(0.1)

y2 = biomasses["ecoli_1"]
y3 = metabolites["M_glc__D_e"]
y4 = metabolites["M_mal__L_e"]

y3 = y3
ts = ts

plt.plot(ts, y2, color="blue", label="biomass")
plt.plot(ts, y3, color="green", label="glucose")
plt.plot(ts, y4, color="red", label="maltase")

# Adding labels and title
plt.xlabel("X-axis")
plt.ylabel("Y-axis")
plt.title("Three Lines Plot")

# Adding legend
plt.legend()

# Displaying the plot
plt.show()
