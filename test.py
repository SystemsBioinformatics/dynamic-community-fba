import cbmpy
from cbmpy.CBModel import Model

ecoli_core: Model = cbmpy.loadModel("models/bigg_models/e_coli_core.xml")
ecoli_core.getReaction("R_EX_glc__D_e").setLowerBound(-0.48)

sol = cbmpy.doFBA(ecoli_core)
print(sol)
