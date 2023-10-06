import cbmpy
from cbmpy.CBModel import Model, Reaction
from dcFBA.Models import CommunityModel
from dcFBA.DynamicModels import DynamicJointFBA, EndPointFBA
import pandas as pd
from dcFBA.Helpers.PlotsEndPointFBA import plot_biomasses, plot_metabolites
import matplotlib.pyplot as plt
import time
from dcFBA.Helpers.OptimalSearch import balanced_search_quick, time_search
import json

iAF1260: Model = cbmpy.loadModel("models/bigg_models/iAF1260.xml")

iAF1260.getReaction("R_EX_lys__L_e").setLowerBound(0)
iAF1260.getReaction("R_EX_leu__L_e").setLowerBound(0)


reaction: Reaction
for reaction in iAF1260.reactions:
    if reaction.getSBOterm() == "SBO:0000627":
        reaction.is_exchange = True

for reaction in iAF1260.reactions:
    if reaction.getLowerBound() == -999999.0:
        reaction.setLowerBound(cbmpy.NINF)
    if reaction.getUpperBound() == 999999.0:
        reaction.setUpperBound(cbmpy.INF)

leucine_knock_out = iAF1260.clone()
lysine_knock_out = iAF1260.clone()

leucine_knock_out.getReaction("R_IPPS").setUpperBound(0)
lysine_knock_out.getReaction("R_DAPDC").setUpperBound(0)

# Set creation of the metabolites to zero
leucine_knock_out.getReaction("R_IPPS").setUpperBound(0)
lysine_knock_out.getReaction("R_DAPDC").setUpperBound(0)


# uptake simulation from paper:

leucine_knock_out.getReaction("R_LYStex").setLowerBound(-1000)
leucine_knock_out.getReaction("R_LYStex").setUpperBound(1000)

lysine_knock_out.getReaction("R_LEUtex").setLowerBound(-1000)
lysine_knock_out.getReaction("R_LEUtex").setUpperBound(1000)


# Restrict the uptake of glucose
leucine_knock_out.getReaction("R_GLCtex_copy1").setUpperBound(10)
leucine_knock_out.getReaction("R_GLCtex_copy2").setUpperBound(0)
lysine_knock_out.getReaction("R_GLCtex_copy1").setUpperBound(10)
lysine_knock_out.getReaction("R_GLCtex_copy2").setUpperBound(0)


# R_FE3tex settings from paper
leucine_knock_out.getReaction("R_FE3tex").setUpperBound(0)
lysine_knock_out.getReaction("R_FE3tex").setUpperBound(0)

# Build the community model
community_model: CommunityModel = CommunityModel(
    [leucine_knock_out, lysine_knock_out],
    ["R_BIOMASS_Ec_iAF1260_core_59p81M", "R_BIOMASS_Ec_iAF1260_core_59p81M"],
    ["dleu", "dlys"],
)

# R_O2tex oxygen
# community_model.getReaction("R_O2tex_dleu").setLowerBound(-15)
# community_model.getReaction("R_O2tex_dlys").setLowerBound(-15)

community_model.getReaction("R_EX_o2_e").setLowerBound(
    -30
)  # Vou;d also be set to 2* 18.5
community_model.getReaction("R_EX_cbl1_e").setLowerBound(
    -0.02
)  # 2 * initial value of 0.01

# df = pd.read_csv(
#     "/Users/stevenwijnen/Bioinformatica/year2/SysBio_lab/major_research_project_sysbio/results/code/ecoli/fva_zero_reactions_iAF1260_community.csv"
# )

# fva_list = first_column_list = df.iloc[:, 0].tolist()

# for rid in fva_list:
#     community_model.deleteReactionAndBounds(rid)

# community_model.deleteNonReactingSpecies()

# First try if accumulation of the metabolites is enough to stop
# the corss feeding
community_model.getReaction("R_EX_leu__L_e").setLowerBound(0)
community_model.getReaction("R_EX_lys__L_e").setLowerBound(0)

community_model.getReaction("R_EX_leu__L_e").setUpperBound(0)
community_model.getReaction("R_EX_lys__L_e").setUpperBound(0)

medium_co_culture = {
    "M_glc__D_e": 11.96,
    "M_leu__L_e": 0,
    "M_lys__L_e": 0,
}  # Glucose from the experiment paper Zhang, Reed


# ep = EndPointFBA(
#     community_model,
#     3,
#     {"dleu": 0.0027, "dlys": 0.0027},
#     medium_co_culture,
#     dt=0.5,
# )
# e = balanced_search_quick(ep, 0.0027 * 2, 0.083)


# a, b = time_search(
#     community_model,
#     {"dleu": 0.0027, "dlys": 0.0027},
#     medium_co_culture,
#     0.5,
#     (0.083, 20),
# )
# print(a, b)
start_time = time.time()
ep = EndPointFBA(
    community_model,
    20,
    {"dleu": 0.0027, "dlys": 0.0027},
    medium_co_culture,
    dt=0.2,
)
print(
    f"Elapsed time to ====build==== model n50 dt0.2: {time.time() - start_time}"
)


ep.balanced_growth(0.0027 * 2, 0.083)
print(ep.simulate())
plot_biomasses(ep)

FBAsol = ep.m_model.getSolutionVector(names=True)
FBAsol = dict(zip(FBAsol[1], FBAsol[0]))
with open("FBAsol_n20dt02_balanced.json", "w") as json_file:
    json.dump(FBAsol, json_file)  # 'indent=4' for pretty-printing


# with open("FBAsol_after_search.json", "w") as json_file:
#     json.dump(FBAsol, json_file)  # 'indent=4' for pretty-printing
