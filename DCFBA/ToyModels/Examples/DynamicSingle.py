import cbmpy
from DCFBA.DynamicModels import DynamicSingleFBA
from DCFBA.Models import KineticsStruct
import matplotlib.pyplot as plt

model = cbmpy.loadModel("data/bigg_models/e_coli_core.xml")

initial_biomass = 0.1
initial_concentrations = {"M_glc__D_e": 10}
kin = KineticsStruct({"R_GLCpts": ("M_glc__D_e", 5, 10)})

ds = DynamicSingleFBA(
    model,
    "R_BIOMASS_Ecoli_core_w_GAM",
    initial_biomass,
    initial_concentrations,
    kinetics=kin,
)

T, metabolites, biomassess, _ = ds.simulate(0.15)


ax = plt.subplot(111)
ax.plot(T, biomassess[""])
ax2 = plt.twinx(ax)
ax2.plot(T, metabolites["M_glc__D_e"], color="r")

ax.set_ylabel("Biomass", color="b")
ax2.set_ylabel("Glucose", color="r")
plt.show()
