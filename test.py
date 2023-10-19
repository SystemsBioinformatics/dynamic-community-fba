import cbmpy
from cbmpy.CBModel import Model, Reaction
from dcFBA.Models import CommunityModel
from dcFBA.DynamicModels import EndPointFBA
import pandas as pd
import time
from dcFBA.Helpers.OptimalSearch import balanced_search_quick, time_search
import json
from dcFBA.Models import KineticsStruct
from dcFBA.DynamicModels import DynamicJointFBA

import matplotlib.pyplot as plt

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
# leucine_knock_out.getReaction("R_LEUtex").setLowerBound(-1)
# leucine_knock_out.getReaction("R_LEUtex").setUpperBound(1)

leucine_knock_out.getReaction("R_LYStex").setLowerBound(-1000)
leucine_knock_out.getReaction("R_LYStex").setUpperBound(1000)


# lysine_knock_out.getReaction("R_LYStex").setLowerBound(-1)
# lysine_knock_out.getReaction("R_LYStex").setUpperBound(1)

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


medium_co_culture = {
    "M_glc__D_e": 11.96,
    "M_leu__L_e": 0,
    "M_lys__L_e": 0,
}  # Glucose from the experiment paper Zhang, Reed


community_model.getReaction("R_IPPS_dleu").setUpperBound(0)
community_model.getReaction("R_DAPDC_dlys").setUpperBound(0)

# Restrict the uptake of glucose
community_model.getReaction("R_GLCtex_copy1_dleu").setUpperBound(10)
community_model.getReaction("R_GLCtex_copy2_dleu").setUpperBound(0)
community_model.getReaction("R_GLCtex_copy1_dlys").setUpperBound(10)
community_model.getReaction("R_GLCtex_copy2_dlys").setUpperBound(0)


# R_O2tex oxygen
community_model.getReaction("R_O2tex_dleu").setLowerBound(-15)
community_model.getReaction("R_O2tex_dlys").setLowerBound(-15)

community_model.getReaction("R_EX_o2_e").setLowerBound(
    -30
)  # could also be set to 2* 18.5
community_model.getReaction("R_EX_cbl1_e").setLowerBound(
    -0.02
)  # 2 * initial value of 0.01


# R_FE3tex
community_model.getReaction("R_FE3tex_dleu").setUpperBound(0)
community_model.getReaction("R_FE3tex_dlys").setUpperBound(0)

# In the paper they forced that all Leucine and Lysine was consumed by the strains
# I.E. there can not be outflow of these metabolites or inflow
community_model.getReaction("R_EX_leu__L_e").setUpperBound(0)
community_model.getReaction("R_EX_lys__L_e").setUpperBound(0)

community_model.getReaction("R_EX_leu__L_e").setLowerBound(0)
community_model.getReaction("R_EX_lys__L_e").setLowerBound(0)

community_model.getReaction("R_LEUtex_dleu").setLowerBound(-1000)
community_model.getReaction("R_LEUtex_dleu").setUpperBound(1000)
community_model.getReaction("R_LYStex_dleu").setLowerBound(-0.056)
community_model.getReaction("R_LYStex_dleu").setUpperBound(-0.056)

community_model.getReaction("R_LYStex_dlys").setLowerBound(-1000)
community_model.getReaction("R_LYStex_dlys").setUpperBound(1000)
community_model.getReaction("R_LEUtex_dlys").setLowerBound(-0.086)
community_model.getReaction("R_LEUtex_dlys").setUpperBound(-0.086)


dj = DynamicJointFBA(
    community_model,
    [0.0027, 0.0027],
    {"M_glc__D_e": 11.96, "M_leu__L_e": 0, "M_lys__L_e": 0},
)


dj.simulate(0.1, n=50, epsilon=0.00001)


pps = dj.get_community_growth_rate()
print(pps)
plt.plot(dj.get_time_points()[:-1], pps)
plt.show()
