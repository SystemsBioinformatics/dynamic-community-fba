import cbmpy
from cbmpy.CBModel import Model, Reaction
import pandas as pd
from dcFBA.Models import CommunityModel
import cbmpy.CBCPLEX as cbc

iAF1260: Model = cbmpy.loadModel("models/bigg_models/iAF1260.xml")
reaction: Reaction
for reaction in iAF1260.reactions:
    if reaction.getSBOterm() == "SBO:0000627":
        reaction.is_exchange = True

for reaction in iAF1260.reactions:
    if reaction.getLowerBound() == -999999.0:
        reaction.setLowerBound(cbmpy.NINF)
    if reaction.getUpperBound() == 999999.0:
        reaction.setUpperBound(cbmpy.INF)

iAF1260.getReaction("R_EX_lys__L_e").setLowerBound(0)
iAF1260.getReaction("R_EX_leu__L_e").setLowerBound(0)

iAF1260_N_L: Model = iAF1260.clone()
leucine = "M_leu__L_c"
knock_out_gene_leucine = "G_b0074"
leucine_gene_knock_out_associated_reaction = "R_IPPS"
leucine_creating_reaction = "R_LEUTAi"

# Knock out gene
iAF1260_N_L.getGene(knock_out_gene_leucine).setInactive()
iAF1260_N_L.getReaction(
    leucine_gene_knock_out_associated_reaction
).setUpperBound(0)
iAF1260_N_L.getReaction(
    leucine_gene_knock_out_associated_reaction
).setLowerBound(0)

# Negative Lysine (K)
iAF1260_N_K = iAF1260.clone()
lysine = "M_lys__L_c"

knock_out_gene_lysine = "G_b2838"
associated_reaction = "R_DAPDC"

# Make sure no lysine enters the system
iAF1260_N_K.getGene(knock_out_gene_lysine).setInactive()
iAF1260_N_K.getReaction(associated_reaction).setUpperBound(0)
iAF1260_N_K.getReaction(associated_reaction).setLowerBound(0)

community_model = CommunityModel(
    [iAF1260_N_L, iAF1260_N_K],
    ["R_BIOMASS_Ec_iAF1260_core_59p81M", "R_BIOMASS_Ec_iAF1260_core_59p81M"],
    ["dleu", "dlys"],
)
community_model.createSpecies(
    "X_c", False, "The community biomass", compartment="e"
)

for _, biomass_id in community_model.get_model_biomass_ids().items():
    reaction: Reaction = community_model.getReaction(biomass_id)
    reaction.createReagent("X_c", 1)

community_model.createReaction("X_comm")
out: Reaction = community_model.getReaction("X_comm")
out.is_exchange = True
out.setUpperBound(cbmpy.INF)
out.setLowerBound(0)
out.createReagent("X_c", -1)

community_model.setObjectiveFlux("X_comm")
print(cbmpy.doFBA(community_model))
# community_model.buildStoichMatrix()

raise Exception
# community_model.getReaction("R_EX_glc__D_e").setLowerBound(-24)
# fvaResultComm = cbc.cplx_FluxVariabilityAnalysis(community_model)
solution = cbmpy.doFBAMinSum(community_model)
print(solution)

FBAsol = community_model.getSolutionVector(names=True)
FBAsol = dict(zip(FBAsol[1], FBAsol[0]))

print(FBAsol["R_LYStex_dleu"])
print(FBAsol["R_LEUtex_dleu"])

print()
print(FBAsol["R_LYStex_dlys"])
print(FBAsol["R_LEUtex_dlys"])
community_model.setObjectiveFlux()
