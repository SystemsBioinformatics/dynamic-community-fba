import matplotlib.pyplot as plt
from cbmpy.CBModel import Model
from dcFBA.ToyModels import model_a, model_b
from dcFBA.DynamicModels.DynamicParallelFBA import DynamicParallelFBA
from dcFBA.Models.Kinetics import KineticsStruct

# TODO discuss with Francesco; If the upperBounds are multiplied by 0.1 than I get the same results
# But I do not think this should be done? Since we multiply by the delta in the update concentration method
# And if we multiply by it in the upper bounds AND for the calculation of the new concentrations we multiply by (dt)^2
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

kin = {
    "Organism_A": KineticsStruct(
        {
            "R_1": ["S_e", 10, 10],
            "R_4": ["B_e", 5.0, 3.0],
            "R_6": ["B_e", 3, 1],
        }
    ),
    "Organsim_B": KineticsStruct(
        {
            "R_1": ["S_e", 10, 10],
            "R_3": ["A_e", 2.0, 1.0],
        }
    ),
}

parallel = DynamicParallelFBA(
    [m_a, m_b], [1, 1], {"S_e": 100, "A_e": 1, "B_e": 2}
)


T, metabolites, biomasses = parallel.simulate(0.1, 300, 0.05)

print(metabolites)
print(biomasses)
import matplotlib.pyplot as plt
from cbmpy.CBModel import Model
from DCFBA.ToyModels import model_a, model_b
from DCFBA.DynamicModels.DynamicParallelFBA import DynamicParallelFBA
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


kin = {
    "Organism_A": KineticsStruct(
        {
            "R_1": ["S_e", 10, 10],
            "R_4": ["B_e", 5.0, 3.0],
            "R_6": ["B_e", 3, 1],
        }
    ),
    "Organsim_B": KineticsStruct(
        {
            "R_1": ["S_e", 10, 10],
            "R_3": ["A_e", 2.0, 1.0],
        }
    ),
}

parallel = DynamicParallelFBA(
    [m_a, m_b], [1, 1], {"S_e": 100, "A_e": 1, "B_e": 2}
)


T, metabolites, biomasses = parallel.simulate(0.1, 169, 0.1)

print(metabolites)
print(biomasses)

print(len(metabolites["S_e"]))

print(biomasses["Organism_A"][-1])
print(biomasses["Organism_B"][-1])

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

# Adding legend
plt.legend()

# Displaying the plot
plt.show()

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
