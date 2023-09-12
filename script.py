import cbmpy
from cbmpy.CBModel import Model, Reaction
from DCFBA.DynamicModels import EndPointFBA
from DCFBA.Helpers.OptimalTimeSearch import search
from DCFBA.Models import CommunityModel

ecoli_core = cbmpy.loadModel("models/bigg_models/e_coli_core.xml")
print(cbmpy.doFBA(ecoli_core))

community_model = CommunityModel(
    [ecoli_core], ["R_BIOMASS_Ecoli_core_w_GAM"], ["ecoli"]
)

# high, obj = search(community_model, {"ecoli": 1.0}, {"M_glc__D_e": 10000}, 0.1)

# print(high, obj)

ep = EndPointFBA(community_model, 1, {"ecoli": 1.0}, {}, dt=1.0)

print(ep.simulate())
