import cbmpy
import matplotlib.pyplot as plt
from cbmpy.CBModel import Model, Reaction
from DCFBA.ToyModels import model_a, model_b
from DCFBA.Models.CommunityModel import CommunityModel
from DCFBA.Models.Kinetics import KineticsStruct
from DCFBA.DynamicModels import EndPointFBA
from DCFBA.Helpers.PlotsEndPointFBA import plot_biomasses, plot_metabolites

lacto_bac: Model = cbmpy.loadModel("models/bigg_models/LBUL_v13.xml")
strep: Model = cbmpy.loadModel(
    "models/bigg_models/strep_therm.xml"
)  # Weird exchanges: R_EX_glc__D_e, R_EX_lac__L_e, R_EX_lcts_e

strep.setReactionBounds("R_EX_gal_e", 0, 24.04 * 0.2)
print(cbmpy.doFBA(strep))

for e in strep.getReactionIds("R_EX"):
    reaction: Reaction = strep.getReaction(e)
    if reaction.getUpperBound() == 0:
        reaction.setUpperBound(1000)

print(cbmpy.doFBA(strep))

print(strep.getReactionIdsAssociatedWithSpecies("M_gal_e"))

community_model = CommunityModel([strep], ["R_biomass_STR"], ["strep"])


# community_model.getReaction("R_EX_fol_e").setLowerBound(0),

n = 2
dj = EndPointFBA(
    community_model,
    n,
    {"strep": 1.0},
    dt=0.1,
)
# print(dj.m_model.getReactionIdsAssociatedWithSpecies("M_gal_e_time0"))
# print(dj.m_model.getReactionIdsAssociatedWithSpecies("M_gal_e_time1"))
# print(dj.m_model.getReactionIdsAssociatedWithSpecies("M_gal_e_time2"))


solution = dj.simulate()
print(solution)
plot_biomasses(dj)
plot_metabolites(dj, {"M_fol_e": 0, "M_gal_e": 0.0})
