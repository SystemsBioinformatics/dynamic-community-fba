import matplotlib.pyplot as plt
from cbmpy.CBModel import Model
from dcFBA.Models import CommunityModel
from dcFBA.DynamicModels import DynamicJointFBA
from dcFBA.DefaultModels import read_default_model

model1: Model = read_default_model("e_coli_core")
model2: Model = read_default_model("strep_therm")

# Set the import bounds for glucose in both models
model1.getReaction("R_GLCpts").setUpperBound(10)
model2.getReaction("R_GLCpts").setUpperBound(6)

# The biomass reactions ids
biomass_reaction_model_1: str = "R_BIOMASS_Ecoli_core_w_GAM"
biomass_reaction_model_2: str = "R_biomass_STR"


# Restrict Lactose uptake
# model2.getReaction("R_LCTSGALex").setLowerBound(-24)
# model2.getReaction("R_LCTSt6").setUpperBound(30)


cm = CommunityModel(
    [model1, model2],
    [biomass_reaction_model_1, biomass_reaction_model_2],
    ["ecoli", "strep"],
)  # Define the community model

dynamic_joint_fba = DynamicJointFBA(
    cm,
    [1.0, 1.0],
    {
        "M_glc__D_e": 100,
        "M_succ_e": 0,
        "M_glu__L_e": 0.0,
        "M_gln__L_e": 0.0,
        "M_lcts_e": 100,
    },
)  # Create a DynamicJointFBA object, set the initial concentrations of glucose and lactose to 100

dynamic_joint_fba.simulate(0.1)
# NMATRIX_TYPE = 'scipy_csr'

biomasses = dynamic_joint_fba.get_biomasses()
metabolites = dynamic_joint_fba.get_metabolites()
time_points = dynamic_joint_fba.get_time_points()
fluxes = dynamic_joint_fba.get_fluxes()

plt.plot(
    time_points, metabolites["M_glc__D_e"], color="blue", label="[Glucose]"
)
plt.plot(
    time_points, metabolites["M_lcts_e"], color="orange", label="[Lactose]"
)

plt.xlabel("Time")
plt.ylabel("Concentration")
plt.legend()
plt.show()


plt.plot(time_points, biomasses["ecoli"], color="orange", label="ecoli")
plt.plot(time_points, biomasses["strep"], color="blue", label="strep")

plt.xlabel("Time")
plt.ylabel("Concentration")
plt.legend()
plt.show()
