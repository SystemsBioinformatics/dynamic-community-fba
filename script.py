import cbmpy
from DCFBA.DynamicModels import DynamicSingleFBA
from DCFBA.Models import KineticsStruct
import matplotlib.pyplot as plt
from DCFBA.ToyModels.model_c import build_model_C

model = build_model_C()
model.getReaction("R_1").setUpperBound(1)

db = DynamicSingleFBA(model, "R_BM_C", 1, {"S_e": 10})

T, metabolites, biomasses, fluxes = db.simulate(dt=0.1)

print(metabolites["S_e"])
print(biomasses)
print(T)
rv6 = list(map(lambda d: d["S_exchange"], fluxes))
print(rv6)

# Plot the results
ax = plt.subplot(111)
ax.plot(T, biomasses[""])
ax2 = plt.twinx(ax)
ax2.plot(T, metabolites["S_e"], color="r")

ax.set_ylabel("Biomass", color="b")
ax2.set_ylabel("S", color="r")
plt.show()
