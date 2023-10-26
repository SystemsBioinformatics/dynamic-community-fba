import cbmpy
from cbmpy.CBModel import Model, Reaction
from dcFBA.DynamicModels import EndPointFBA
from dcFBA.Models import CommunityModel


iAF1260: Model = cbmpy.loadModel("./models/bigg_models/iAF1260.xml")

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
leucine_knock_out.getReaction("R_LYStex").setUpperBound(1000)


lysine_knock_out.getReaction("R_LEUtex").setLowerBound(-1000)
lysine_knock_out.getReaction("R_LEUtex").setUpperBound(1000)


# Restrict the release of glucose
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

community_model.getReaction("R_EX_o2_e").setLowerBound(
    -30
)  # Vou;d also be set to 2* 18.5
community_model.getReaction("R_EX_cbl1_e").setLowerBound(
    -0.02
)  # 2 * initial value of 0.01


ep = EndPointFBA(
    community_model,
    8,
    {"dleu": 0.0027, "dlys": 0.0027},
    {"M_glc__D_e": 11.96, "M_leu__L_e": 0, "M_lys__L_e": 0},
    0.8,
)


Xin = 0.0027 * 2
ep.balanced_growth(Xin, 0.083 - Xin)
