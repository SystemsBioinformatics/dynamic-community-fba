import cbmpy
from cbmpy.CBModel import Model, Reaction
from dcFBA.Models import CommunityModel
from dcFBA.DynamicModels import DynamicJointFBA
import pandas as pd
from dcFBA.Helpers.PlotsEndPointFBA import plot_biomasses, plot_metabolites
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

community_model.getReaction("R_EX_o2_e").setLowerBound(
    -30
)  # Vou;d also be set to 2* 18.5
community_model.getReaction("R_EX_cbl1_e").setLowerBound(
    -0.02
)  # 2 * initial value of 0.01


# First try if accumulation of the metabolites is enough to stop
# the corss feeding
community_model.getReaction("R_EX_leu__L_e").setLowerBound(0)
community_model.getReaction("R_EX_lys__L_e").setLowerBound(0)

community_model.getReaction("R_EX_leu__L_e").setUpperBound(0)
community_model.getReaction("R_EX_lys__L_e").setUpperBound(0)

dj_uptake = DynamicJointFBA(
    community_model,
    [0.0027, 0.0027],
    {"M_glc__D_e": 11.96, "M_leu__L_e": 0, "M_lys__L_e": 0},
)


def deviate_func(sim, used_time, run_condition):
    if (
        sim.m_biomass_concentrations["dleu"][-1]
        + sim.m_biomass_concentrations["dlys"][-1]
        >= 0.083
    ):
        # Stop the simulation by setting community reaction to zero, solution will be zero or nan
        sim.m_model.getReaction("X_comm").setUpperBound(0)

    return 0


d = deviate_func
(
    T_uptake,
    metabolites_uptake,
    biomasses_uptake,
    fluxes_uptake,
) = dj_uptake.simulate(0.1, epsilon=0.00001, deviate=d, n=160)

plt.plot(T_uptake, biomasses_uptake["dleu"], color="blue", label="dleusine")
plt.plot(T_uptake, biomasses_uptake["dlys"], color="orange", label="dlysine")
plt.xlabel("Time")
plt.ylabel("Concentration")
plt.legend()
plt.show()
