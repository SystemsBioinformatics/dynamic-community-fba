import cbmpy
from cbmpy.CBModel import Model
from endPointFBA import CommunityModel

model1: Model = cbmpy.loadModel("data/bigg_models/e_coli_core.xml")
model2: Model = cbmpy.loadModel("data/bigg_models/strep_therm.xml")


biomass_reaction_model_1 = "R_BIOMASS_Ecoli_core_w_GAM"
biomass_reaction_model_2 = "R_biomass_STR"
community_model: CommunityModel = CommunityModel(
    [model1, model2],
    [
        biomass_reaction_model_1,
        biomass_reaction_model_2,
    ],
)


print(
    model1.getActiveObjectiveReactionIds(),
    model2.getActiveObjectiveReactionIds(),
)
