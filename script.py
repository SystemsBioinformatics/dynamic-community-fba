import matplotlib.pyplot as plt
from cbmpy.CBModel import Model
from DCFBA.ToyModels import model_a, model_b
from DCFBA.Models.CommunityModel import CommunityModel
from DCFBA.DynamicModels.DynamicJointFBA import DynamicJointFBA

m_a: Model = model_a.build_joint_fba_model_A()
m_a.getReaction("R_1").setUpperBound(10)
m_a.getReaction("R_3").setUpperBound(3)
m_a.getReaction("R_5").setUpperBound(1)

m_b: Model = model_b.build_joint_fba_model_B()
m_b.getReaction("R_1").setUpperBound(10)
m_b.getReaction("R_3").setUpperBound(2)
# m_b.getReaction("R_5").setUpperBound(1)
# m_b.getReaction("R_5").setLowerBound(0)


community_model = CommunityModel(
    [m_a, m_b], ["R_BM_A", "R_BM_B"], ["modelA", "modelB"]
)
community_model.getReaction("B_exchange").setLowerBound(-0.001)
community_model.getReaction("A_exchange").setLowerBound(-0.001)
community_model.getReaction("S_exchange").setLowerBound(-100)

# community_model.addUserConstraint(
#     "B_constraint",
#     [
#         [1, "R_import_B_modelA"],
#         [-3, "R_6_modelA"],
#         [-1, "R_4_modelA"],
#     ],
#     ">=",
#     0,
# )

# community_model.addUserConstraint(
#     "A_constraint",
#     [
#         [1, "R_import_A_modelB"],
#         [-3, "R_3_modelB"],
#         [-1, "R_4_modelB"],
#     ],
#     ">=",
#     0,
# )

# community_model.addUserConstraint(
#     "A_exchange_constrained",
#     [
#         [1, "R_export_A_modelA"],
#         [-1, "R_import_A_modelB"],
#         [-1, "A_exchange"],
#     ],
#     ">=",
#     0,
# )

# community_model.addUserConstraint(
#     "B_exchange_constrained",
#     [
#         [1, "R_export_B_modelB"],
#         [-1, "B_exchange"],
#         [-1, "R_import_B_modelA"],
#     ],
#     ">=",
#     0,
# )


# # print(community_model.getReaction("R_1_modelA").getUpperBound())
# print(community_model.getReaction("R_1_modelB").getUpperBound())
# print(community_model.getReaction("R_5_modelB").getUpperBound())

dj = DynamicJointFBA(community_model, [0.1, 0.1], {"S_e": 100})
T, metabolites, biomasses, fluxes = dj.simulate(0.1)

rv1 = list(map(lambda d: d["R_1_modelA"], fluxes))
rv2 = list(map(lambda d: d["R_2_modelA"], fluxes))
rv3 = list(map(lambda d: d["R_3_modelA"], fluxes))
rv4 = list(map(lambda d: d["R_4_modelA"], fluxes))
rv5 = list(map(lambda d: d["R_5_modelA"], fluxes))
rv6 = list(map(lambda d: d["R_6_modelA"], fluxes))


# index = 1
# label_colors = [
#     "blue",
#     "green",
#     "red",
#     "cyan",
#     "magenta",
#     "yellow",
#     "black",
#     "white",
#     "orange",
#     "purple",
# ]

# for r in [rv1, rv2, rv3, rv4, rv5, rv6]:
#     plt.plot(T, r, color=label_colors[index], label=f"r{index}")
#     index += 1

plt.plot(T, metabolites["A_e"], color="blue", label="model a")
# plt.plot(T, biomasses["modelB"], color="orange", label="model b")

# Adding labels and title
plt.xlabel("Time")
plt.ylabel("Value")
plt.title("Plotting Two Lines")

# Adding legend
plt.legend()

# Displaying the plot
plt.show()
