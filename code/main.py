import cbmpy
from cbmpy.CBModel import Model
from endPointFBA.CommunityModel import CommunityModel
from endPointFBA.DynamicJointFBA import DynamicJointFBA

model1: Model = cbmpy.loadModel("data/bigg_models/e_coli_core.xml")
model_2 = model1.clone()

model1.getReaction("R_GLCpts").setUpperBound(10)
combined_model = CommunityModel(
    [model1],
    ["R_BIOMASS_Ecoli_core_w_GAM"],
    ["ecoli_1"],
)  # Create a CommunityModel of two  E. coli strains competing for resources


# Create the joint FBA object with initial biomasses and the initial
# concentration of glucose
dynamic_fba = DynamicJointFBA(
    combined_model, [0.1], {"M_glc__D_e": 10}, {"R_GLCpts_ecoli_1": [5, 10]}
)

import matplotlib.pyplot as plt

ts, metabolites, biomasses = dynamic_fba.simulate(0.1)
print(metabolites["M_glc__D_e"])
print(metabolites["X_c"])

ax = plt.subplot(161)
ax.plot(ts, metabolites["M_glc__D_e"])
ax2 = plt.twinx(ax)
ax2.plot(ts, metabolites["X_c"], color="r")

ax.set_ylabel("glucose", color="b")
ax2.set_ylabel("X_comm", color="r")

# ax3 = plt.twinx(ax)
# ax3.plot(ts, biomasses["ecoli_2"], color="y")
# ax3.set_ylabel("Biomass ecoli 2", color="y")

plt.show()

# # Perform FBA on the new joint FBA model object
# solution = cbmpy.doFBA(dynamic_fba.get_joint_model())
# print(solution)
