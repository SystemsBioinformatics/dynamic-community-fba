"""Just a small script file to easily run
"""

import cbmpy

from cbmpy.CBModel import Model
from DCFBA.Models.CommunityModel import CommunityModel
from DCFBA.DynamicModels import DynamicJointFBA

model1: Model = cbmpy.loadModel("data/bigg_models/e_coli_core.xml")
model1.getReaction("R_GLCpts").setUpperBound(10)

model_2 = model1.clone()

combined_model = CommunityModel(
    [model1, model_2],
    ["R_BIOMASS_Ecoli_core_w_GAM", "R_BIOMASS_Ecoli_core_w_GAM"],
    ["ecoli_1", "ecoli_2"],
)  # Create a CommunityModel of two  E. coli strains competing for resources


# Create the joint FBA object with initial biomasses and the initial
# concentration of glucose
dynamic_fba = DynamicJointFBA(
    combined_model,
    [0.1, 0.1],
    {"M_glc__D_e": 10},
)

#    Kinetics({"R_GLCpts_ecoli_1": [5, 10], "R_GLCpts_ecoli_2": [5, 10]}),
ts, metabolites, biomasses = dynamic_fba.simulate(0.1)
print(ts)
# print(metabolites["M_glc__D_e"])
# print()
# print(biomasses["ecoli_1"])
# print()
# print(biomasses["ecoli_2"])


# ax = plt.subplot(111)
# ax.plot(ts, metabolites["M_glc__D_e"])
# ax2 = plt.twinx(ax)
# ax2.plot(ts, biomasses["ecoli_2"], color="r")

# ax.set_ylabel("glucose", color="b")
# ax2.set_ylabel("ecoli_2", color="r")

# ax3 = plt.twinx(ax)
# ax3.plot(ts, biomasses["ecoli_2"], color="y")
# ax3.set_ylabel("Biomass ecoli 2", color="y")

# plt.show()

# # Perform FBA on the new joint FBA model object
# solution = cbmpy.doFBA(dynamic_fba.get_joint_model())
# print(solution)
