import matplotlib.pyplot as plt
from cbmpy.CBModel import Model, Reaction
from DCFBA.ToyModels import model_a, model_b
from DCFBA.Models.CommunityModel import CommunityModel
from DCFBA.DynamicModels.DynamicJointFBA import DynamicJointFBA
from DCFBA.Models.Kinetics import KineticsStruct

m_a: Model = model_a.build_toy_model_fba_A()
m_a.getReaction("R_1").setUpperBound(10)
m_a.getReaction("R_4").setUpperBound(3)
m_a.getReaction("R_6").setUpperBound(1)


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


community_model = CommunityModel(
    [m_a, m_b], ["R_BM_A", "R_BM_B"], ["modelA", "modelB"]
)


community_model.getReaction("S_e_exchange").setLowerBound(-100000)
community_model.getReaction("A_e_exchange").setLowerBound(-100000)
community_model.getReaction("B_e_exchange").setLowerBound(-100000)

dj = DynamicJointFBA(
    community_model, [1, 1], {"S_e": 100, "A_e": 0.01, "B_e": 0.01}
)
T, metabolites, biomasses, fluxes = dj.simulate(0.01)
index = 0


plt.plot(T, metabolites["S_e"], color="blue", label="Metabolite s")
plt.xlabel("Time")
plt.ylabel("Concentration")
plt.legend()
plt.show()


plt.plot(T, metabolites["A_e"], color="blue", label="Metabolite a")
plt.xlabel("Time")
plt.ylabel("Concentration")
plt.legend()
plt.show()


# Adding labels and title

plt.plot(T, metabolites["B_e"], color="blue", label="Metabolite B")
plt.xlabel("Time")
plt.ylabel("Concentration")
# Adding legend
plt.legend()
plt.show()

plt.plot(T, biomasses["modelA"], color="blue", label="modelA")
plt.plot(T, biomasses["modelB"], color="orange", label="modelB")
plt.xlabel("Time")
plt.ylabel("Concentration")
plt.legend()
plt.show()
