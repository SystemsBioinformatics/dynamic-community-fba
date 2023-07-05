import cbmpy
from cbmpy.CBModel import Model
from endPointFBA.CommunityModel import CommunityModel
from endPointFBA.DynamicJointFBA import DynamicJointFBA

model1: Model = cbmpy.loadModel("data/bigg_models/e_coli_core.xml")
model_2 = model1.clone()

combined_model = CommunityModel(
    [model1, model_2],
    ["R_BIOMASS_Ecoli_core_w_GAM", "R_BIOMASS_Ecoli_core_w_GAM"],
    ["ecoli_1", "ecoli_2"],
)  # Create a CommunityModel of two  E. coli strains competing for resources


# Create the joint FBA object with initial biomasses and the initial concentration of glucose
dynamic_fba = DynamicJointFBA(
    combined_model,
    [0.1, 0.1],
    {"M_glc__D_e": 10},
)

# Perform FBA on the new joint FBA model object
solution = cbmpy.doFBA(dynamic_fba.get_joint_model())
print(solution)

dynamic_fba.s
