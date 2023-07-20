import matplotlib.pyplot as plt
from cbmpy.CBModel import Model
from ..ToyModels import model_a, model_b
from ..Models.CommunityModel import CommunityModel
from ..DynamicModels.DynamicJointFBA import DynamicJointFBA
from ..Models.Kinetics import Kinetics

m_a: Model = model_a.build_toy_model_fba_A()
m_a.getReaction("R_1").setUpperBound(10)
m_a.getReaction("R_4").setUpperBound(3)
m_a.getReaction("R_6").setUpperBound(1)

# m_a.getReaction("R_import_B").setUpperBound(1)

m_a.getReaction("R_1").setLowerBound(0)
m_a.getReaction("R_4").setLowerBound(0)
m_a.getReaction("R_6").setLowerBound(0)

m_b: Model = model_b.build_toy_model_fba_B()
m_b.getReaction("R_1").setUpperBound(10)
m_b.getReaction("R_3").setUpperBound(1)
m_b.getReaction("R_5").setUpperBound(1)

m_b.getReaction("R_1").setLowerBound(0)
m_b.getReaction("R_3").setLowerBound(0)
m_b.getReaction("R_5").setLowerBound(0)

# m_b.getReaction("R_import_A").setUpperBound(3)


# TODO discuss kinetics.
# thoughts:
#   The concentration of an internal metabolite is zero
#   if we would use the flux * coefficient instead of the concentration
#   It would not work. SInce than we set the upper bound to be the MM
#   calculated with the flux of the previous concentrations
#   and maybe this run many more Y (eg) is made because more
#   S is imorted so than we use an not accurate upper bound
kin = Kinetics(
    {
        "R_1_modelA": ["", 10, 10],
        "R_4_modelA": ["B_e", 5, 3],
        "R_6_modelA": ["B_e", 3, 1],
        # B
        "R_1_modelB": ["", 10, 10],
        "R_3_modelB": ["A_e", 2, 1],
    }
)

community_model = CommunityModel(
    [m_a, m_b], ["R_BM_A", "R_BM_B"], ["modelA", "modelB"]
)

community_model.getReaction("B_exchange").setLowerBound(-0.001)
community_model.getReaction("A_exchange").setLowerBound(-0.001)
community_model.getReaction("S_exchange").setLowerBound(-100)


# # print(community_model.getReaction("R_1_modelA").getUpperBound())
# print(community_model.getReaction("R_1_modelB").getUpperBound())
# print(community_model.getReaction("R_5_modelB").getUpperBound())

dj = DynamicJointFBA(community_model, [1, 1], {"S_e": 100})
T, metabolites, biomasses, fluxes = dj.simulate(0.1)

rv1 = list(map(lambda d: d["R_1_modelA"], fluxes))
rv2 = list(map(lambda d: d["R_2_modelA"], fluxes))
rv3 = list(map(lambda d: d["R_3_modelA"], fluxes))
rv4 = list(map(lambda d: d["R_4_modelA"], fluxes))
rv5 = list(map(lambda d: d["R_5_modelA"], fluxes))
rv6 = list(map(lambda d: d["R_6_modelA"], fluxes))


# index = 1
# label_colors = [
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

plt.show()


# Adding labels and title
plt.xlabel("Time")
plt.ylabel("Value")
plt.title("Plotting Two Lines")
plt.plot(T, metabolites["S_e"], color="blue", label="Se")
plt.show()
# Adding legend
plt.legend()

plt.plot(T, biomasses["modelA"], color="blue", label="modelA")
plt.plot(T, biomasses["modelB"], color="orange", label="modelB")

# Displaying the plot
plt.show()
