import cbmpy

model = cbmpy.loadModel("data/bigg_models/e_coli_core.xml")

reaction = model.getReaction("R_EX_glc__D_e")
reaction.setLowerBound(0)

model.getReaction("R_EX_mal__L_e").setLowerBound(-1)

cbmpy.doFBA(model)
