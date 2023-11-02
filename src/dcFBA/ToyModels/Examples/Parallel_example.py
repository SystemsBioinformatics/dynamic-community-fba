import matplotlib.pyplot as plt
from cbmpy.CBModel import Model
from dcFBA.ToyModels import model_a, model_b
from dcFBA.DynamicModels.DynamicParallelFBA import DynamicParallelFBA
from dcFBA.Models.KineticsStruct import KineticsStruct


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


parallel = DynamicParallelFBA(
    [m_a, m_b], [1, 1], {"S_e": 100, "A_e": 1, "B_e": 2}
)
parallel.simulate(0.01)

T = parallel.get_time_points()
metabolites = parallel.get_metabolites()
biomasses = parallel.get_biomasses()

plt.plot(T, metabolites["A_e"], color="blue", label="A")
plt.plot(T, metabolites["B_e"], color="orange", label="B")
plt.legend()
plt.show()

plt.plot(T, metabolites["S_e"], color="orange", label="S_e")
plt.show()


plt.plot(T, biomasses["Organism_A"], color="blue", label="A")
plt.plot(T, biomasses["Organism_B"], color="orange", label="B")

# Adding labels and title
plt.xlabel("Time")
plt.ylabel("Value")
plt.title("Plotting Two Lines")
plt.legend()
plt.show()
