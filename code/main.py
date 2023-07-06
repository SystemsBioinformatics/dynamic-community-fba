import cbmpy

from toy_models.model_c import build_model_C
from endPointFBA.DynamicParallelFBA import DynamicParallelFBA

# model1 = build_model_C()
# model1.getReaction("R_1").setUpperBound(1)
# model2 = build_model_C()
# model2.setId("Organism_C_2")
# model2.getReaction("R_1").setUpperBound(1)

# print(model2.getId())


model1 = cbmpy.loadModel("data/bigg_models/e_coli_core.xml")
for rid in model1.getReactionIds():
    reaction = model1.getReaction(rid)
    if rid.startswith("R_EX"):
        reaction.is_exchange = True

model1.getReaction("R_GLCpts").setUpperBound(10)

model2 = model1.clone()
model2.setId("ecoli_2")


test: DynamicParallelFBA = DynamicParallelFBA(
    [model1, model2], [0.1, 0.1], {"M_glc__D_e": 10}
)

ts, metabolites, biomasses = test.simulate(0.1)
print(ts)
print(metabolites["M_glc__D_e"])
print(biomasses)


# import cbmpy
# from cbmpy.CBModel import Model
# from endPointFBA.CommunityModel import CommunityModel
# from endPointFBA.DynamicJointFBA import DynamicJointFBA

# model1: Model = cbmpy.loadModel("data/bigg_models/e_coli_core.xml")
# model_2 = model1.clone()

# model1.getReaction("R_GLCpts").setUpperBound(10)
# combined_model = CommunityModel(
#     [model1],
#     ["R_BIOMASS_Ecoli_core_w_GAM"],
#     ["ecoli_1"],
# )  # Create a CommunityModel of two  E. coli strains competing for resources


# # Create the joint FBA object with initial biomasses and the initial
# # concentration of glucose
# dynamic_fba = DynamicJointFBA(
#     combined_model, [0.1], {"M_glc__D_e": 10}, {"R_GLCpts_ecoli_1": [5, 10]}
# )

import matplotlib.pyplot as plt

# ts, metabolites, biomasses = dynamic_fba.simulate(0.1)
# print(metabolites["M_glc__D_e"])
# print(metabolites["X_c"])

ax = plt.subplot(111)
ax.plot(ts, metabolites["M_glc__D_e"])
ax2 = plt.twinx(ax)
ax2.plot(ts, biomasses["e_coli_core"], color="r")

ax.set_ylabel("glucose", color="b")
ax2.set_ylabel("X_comm", color="r")

# ax3 = plt.twinx(ax)
# ax3.plot(ts, biomasses["ecoli_2"], color="y")
# ax3.set_ylabel("Biomass ecoli 2", color="y")

plt.show()

# # Perform FBA on the new joint FBA model object
# solution = cbmpy.doFBA(dynamic_fba.get_joint_model())
# print(solution)
