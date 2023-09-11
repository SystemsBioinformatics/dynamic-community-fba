import matplotlib.pyplot as plt
import cbmpy
from cbmpy.CBModel import Model
from DCFBA.Models import CommunityModel
from DCFBA.DynamicModels import DynamicJointFBA


model1: Model = cbmpy.loadModel(
    "models/bigg_models/e_coli_core.xml"
)  # load e_coli core
model2: Model = cbmpy.loadModel(
    "models/bigg_models/strep_therm.xml"
)  # load the Streptococcus model

# Set the import bounds for glucose in both models
model1.getReaction("R_GLCpts").setUpperBound(10)
model2.getReaction("R_GLCpts").setUpperBound(6)

# The biomass reactions ids
biomass_reaction_model_1: str = "R_BIOMASS_Ecoli_core_w_GAM"
biomass_reaction_model_2: str = "R_biomass_STR"

cm = CommunityModel(
    [model1, model2],
    [biomass_reaction_model_1, biomass_reaction_model_2],
    ["ecoli", "strep"],
)  # Define the community model

dynamic_fba = DynamicJointFBA(
    cm, [1.0, 1.0], {"M_glc__D_e": 100, "M_gal_e": 0, "M_lcts_e": 100}
)  # Create a DynamicJointFBA object, set the initial concentrations of glucose and lactose to 100

T, metabolites, biomasses, fluxes = dynamic_fba.simulate(0.1)

plt.plot(T, metabolites["M_glc__D_e"], color="blue", label="[Glucose]")
plt.plot(T, metabolites["M_lcts_e"], color="orange", label="[Lactose]")

plt.xlabel("Time")
plt.ylabel("Concentration")
plt.legend()
plt.show()
