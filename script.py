import cbmpy
import matplotlib.pyplot as plt
from cbmpy.CBModel import Model, Reaction
from DCFBA.ToyModels import model_a, model_b
from DCFBA.Models.CommunityModel import CommunityModel
from DCFBA.Models.Kinetics import KineticsStruct
from DCFBA.DynamicModels import EndPointFBA
from DCFBA.Helpers.PlotsEndPointFBA import plot_biomasses, plot_metabolites
from DCFBA.Helpers.OptimalTimeSearch import search


ecoli: Model = cbmpy.loadModel("models/bigg_models/e_coli_core.xml")


print("---------------------")
cm = CommunityModel([ecoli], ["R_BIOMASS_Ecoli_core_w_GAM"], ["ecoli"])

print("\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\")
n, obj = search(cm, {"ecoli": 1.0}, {}, dt=1.0)
print(n, obj)
print("\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\")


print("==========1==========")
ep = EndPointFBA(cm, 1, {"ecoli": 1.0}, initial_concentrations={}, dt=0.1)
print(ep.simulate())

# FBAsol2 = ep.m_model.getSolutionVector(names=True)
# FBAsol2 = dict(zip(FBAsol2[1], FBAsol2[0]))
# print(FBAsol2)
# for key, value FBAsol.
# plot_biomasses(ep)
