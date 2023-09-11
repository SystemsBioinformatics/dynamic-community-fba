import cbmpy
from cbmpy.CBModel import Model

lacto: Model = cbmpy.loadModel("models/bigg_models/LBUL_v13.xml")
strep: Model = cbmpy.loadModel("models/bigg_models/strep_therm.xml")

print(lacto.getActiveObjectiveReactionIds())
print(strep.getActiveObjectiveReactionIds())

print(lacto.getReaction("R_biomass_LBUL").getEquation())
print()
print(strep.getReaction("R_biomass_STR").getEquation())
