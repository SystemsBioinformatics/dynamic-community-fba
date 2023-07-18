import matplotlib.pyplot as plt
from cbmpy.CBModel import Model
from DCFBA.ToyModels import model_a, model_b
from DCFBA.DynamicModels.DynamicParallelFBA import DynamicParallelFBA
from DCFBA.Models.Kinetics import Kinetics


m_a: Model = model_a.build_toy_model_fba_A()
m_a.getReaction("R_1").setUpperBound(10)
m_a.getReaction("R_4").setUpperBound(3)
m_a.getReaction("R_6").setUpperBound(1)

m_a.getReaction("R_1").setLowerBound(0)
m_a.getReaction("R_4").setLowerBound(0)
m_a.getReaction("R_6").setLowerBound(0)

m_a.getReaction("B_exchange").setLowerBound(-2)
m_a.getReaction("A_exchange").setLowerBound(0)


m_b: Model = model_b.build_toy_model_fba_B()
m_b.getReaction("R_1").setUpperBound(10)
m_b.getReaction("R_3").setUpperBound(1)
m_b.getReaction("R_5").setUpperBound(1)

m_b.getReaction("R_1").setLowerBound(0)
m_b.getReaction("R_3").setLowerBound(0)
m_b.getReaction("R_5").setLowerBound(0)

m_b.getReaction("B_exchange").setLowerBound(0)
m_b.getReaction("A_exchange").setLowerBound(-1)

kin = {
    "Organism_A": Kinetics(
        {
            "R_1": ["S_e", 10, 10],
            "R_4": ["B_e", 5.0, 3.0],
            "R_6": ["B_e", 3, 1],
        }
    ),
    "Organsim_B": Kinetics(
        {
            "R_1": ["S_e", 10, 10],
            "R_3": ["A_e", 2.0, 1.0],
        }
    ),
}

parallel = DynamicParallelFBA(
    [m_a, m_b], [1, 1], {"S_e": 100, "A_e": 1, "B_e": 2}, kinetics=kin
)


T, metabolites, biomasses = parallel.simulate(0.1, 0.05)

print(metabolites)
print(biomasses)

print(len(metabolites["S_e"]))


plt.plot(T, metabolites["A_e"], color="blue", label="A")
plt.plot(T, metabolites["B_e"], color="orange", label="B")
plt.legend()

plt.show()


plt.plot(T, biomasses["Organism_A"], color="blue", label="A")
plt.plot(T, biomasses["Organism_B"], color="orange", label="B")

# Adding labels and title
plt.xlabel("Time")
plt.ylabel("Value")
plt.title("Plotting Two Lines")

# Adding legend
plt.legend()

# Displaying the plot
plt.show()
