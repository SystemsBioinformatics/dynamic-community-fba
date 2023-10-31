from cbmpy.CBModel import Model
from dcFBA.DynamicModels import DynamicSingleFBA
from dcFBA.DefaultModels import read_default_model
from dcFBA.Models import KineticsStruct
import matplotlib.pyplot as plt

model: Model = read_default_model("e_coli_core")

# Set bounds on the glucose import
model.getReaction("R_GLCpts").setUpperBound(10)

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

ds.simulate(0.15)

T = ds.get_time_points()
metabolites = ds.get_metabolites()
biomass = ds.get_biomass()

ax = plt.subplot(111)
ax.plot(T, biomass)
ax2 = plt.twinx(ax)
ax2.plot(T, metabolites["M_glc__D_e"], color="r")

ax.set_ylabel("Biomass", color="b")
ax2.set_ylabel("Glucose", color="r")
plt.show()
