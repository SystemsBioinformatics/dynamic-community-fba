import cbmpy
import matplotlib.pyplot as plt

from dcFBA.Models import CommunityModel
from dcFBA.DynamicModels import DynamicParallelFBA
from cbmpy.CBModel import Model, Reaction


iAF1260: Model = cbmpy.loadModel("models/bigg_models/iAF1260.xml")
reaction: Reaction = iAF1260.getReaction("R_LEUTAi")

# Set exchange reactions
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


def get_leucine_knock_out_model():
    # Negative leucine (L)
    iAF1260_N_L: Model = iAF1260.clone()
    leucine = "M_leu__L_c"
    knock_out_gene_leucine = "G_b0074"
    leucine_gene_knock_out_associated_reaction = "R_IPPS"
    leucine_creating_reaction = "R_LEUTAi"

    # Knock out gene
    iAF1260_N_L.getGene(knock_out_gene_leucine).setInactive()

    return iAF1260_N_L


def get_lysine_knock_out_model():
    # Negative Lysine (K)
    iAF1260_N_K = iAF1260.clone()
    lysine = "M_lys__L_c"

    knock_out_gene_lysine = "G_b2838"
    associated_reaction = "R_DAPDC"

    # Make sure no lysine enters the system
    iAF1260_N_K.getGene(knock_out_gene_lysine).setInactive()

    return iAF1260_N_K


leucine_knock_out = get_leucine_knock_out_model()
lysine_knock_out = get_lysine_knock_out_model()

# Set creation of the metabolites to zero
leucine_knock_out.getReaction("R_IPPS").setUpperBound(0)
lysine_knock_out.getReaction("R_DAPDC").setUpperBound(0)


leucine_knock_out.getReaction("R_LYStex").setLowerBound(-1000)
# Only way we can make sure the leucine knockout does not take up lysine
leucine_knock_out.getReaction("R_LYStex").setUpperBound(0)


lysine_knock_out.getReaction("R_LEUtex").setLowerBound(-1000)
# Only way we can make sure the lysine knockout does not take up leucine
lysine_knock_out.getReaction("R_LEUtex").setUpperBound(0)


# #R_FE3tex settings from paper
leucine_knock_out.getReaction("R_FE3tex").setUpperBound(0)
lysine_knock_out.getReaction("R_FE3tex").setUpperBound(0)

# # constrain release
# leucine_knock_out.getReaction("R_LEUtex").setLowerBound(0.08)
# leucine_knock_out.getReaction("R_LEUtex").setUpperBound(0.08)
# leucine_knock_out.getReaction("R_LYStex").setLowerBound(-1000)
# leucine_knock_out.getReaction("R_LYStex").setUpperBound(0)

# # constrain release
# lysine_knock_out.getReaction("R_LYStex").setLowerBound(0.06)
# lysine_knock_out.getReaction("R_LYStex").setUpperBound(0.06)
# lysine_knock_out.getReaction("R_LEUtex").setLowerBound(-1000)
# lysine_knock_out.getReaction("R_LEUtex").setUpperBound(0)


# Restrict the release of glucose
leucine_knock_out.getReaction("R_GLCtex_copy1").setUpperBound(10)
leucine_knock_out.getReaction("R_GLCtex_copy2").setUpperBound(0)
lysine_knock_out.getReaction("R_GLCtex_copy1").setUpperBound(10)
lysine_knock_out.getReaction("R_GLCtex_copy2").setUpperBound(0)


# R_FE3tex settings from paper
leucine_knock_out.getReaction("R_FE3tex").setUpperBound(0)
lysine_knock_out.getReaction("R_FE3tex").setUpperBound(0)

# leucine_knock_out.getReaction("R_EX_leu__L_e").setUpperBound(0)
# leucine_knock_out.getReaction("R_EX_lys__L_e").setUpperBound(0)

# lysine_knock_out.getReaction("R_EX_leu__L_e").setUpperBound(0)
# lysine_knock_out.getReaction("R_EX_lys__L_e").setUpperBound(0)


leucine_knock_out.setId("dleu")
lysine_knock_out.setId("dlys")


dpFBA = DynamicParallelFBA(
    [leucine_knock_out, lysine_knock_out],
    [0.0027, 0.0027],
    {
        "M_glc__D_e": 11.96,
        "M_leu__L_e": 0.001,
        "M_lys__L_e": 0.001,
    },
)


def deviate_func(sim, used_time, run_condition):
    if (
        sim.m_biomass_concentrations["dleu"][-1]
        + sim.m_biomass_concentrations["dlys"][-1]
        >= 0.083
    ):
        # Stop the simulation by setting community reaction to zero, solution will be zero or nan
        sim.m_models[0].getReaction(
            "R_BIOMASS_Ec_iAF1260_core_59p81M"
        ).setUpperBound(0)
        sim.m_models[1].getReaction(
            "R_BIOMASS_Ec_iAF1260_core_59p81M"
        ).setUpperBound(0)
    return 0


T, metabolites, biomasses, fluxes = dpFBA.simulate(
    0.1, n=10, epsilon=0.000001, deviate=deviate_func
)

plt.plot(T, biomasses["dleu"], color="blue", label="A")
plt.plot(T, biomasses["dlys"], color="orange", label="B")

plt.show()

plt.plot(T, metabolites["M_lys__L_e"], color="blue", label="A")
plt.plot(T, metabolites["M_leu__L_e"], color="orange", label="B")
plt.legend()
plt.show()
