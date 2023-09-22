import matplotlib.pyplot as plt
from cbmpy.CBModel import Model, Reaction
from dcFBA.ToyModels import model_a, model_b
from dcFBA.Models.CommunityModel import CommunityModel
from dcFBA.DynamicModels.DynamicJointFBA import DynamicJointFBA
from dcFBA.Models.Kinetics import KineticsStruct

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

kin = KineticsStruct(
    {
        "R_1_modelA": ["S_e", 10, 10],
        "R_4_modelA": ["B_e", 5, 3],
        "R_6_modelA": ["B_e", 3, 1],
        # B
        "R_1_modelB": ["S_e", 10, 10],
        "R_3_modelB": ["A_e", 2, 1],
    }
)


community_model = CommunityModel(
    [m_a, m_b], ["R_BM_A", "R_BM_B"], ["modelA", "modelB"]
)


dj = DynamicJointFBA(
    community_model,
    [1, 1],
    {"S_e": 100, "A_e": 0.01, "B_e": 0.01},
    kinetics=kin,
)
T, metabolites, biomasses, fluxes = dj.simulate(0.1)
index = 0

# Obtain fluxes over time
rv1 = list(map(lambda d: d["R_1_modelA"], fluxes))
rv2 = list(map(lambda d: d["R_2_modelA"], fluxes))


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
