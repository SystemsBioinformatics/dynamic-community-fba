import cbmpy
from cbmpy.CBModel import Model, Reaction
from dcFBA.Models import CommunityModel
from dcFBA.DynamicModels import EndPointFBA
from dcFBA.Helpers import ScanUnusedReactions
import time
import pandas as pd

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

iAF1260_2 = iAF1260.clone()

iAF1260.getReaction("R_GLCtex_copy1").setUpperBound(8)
iAF1260.getReaction("R_GLCtex_copy2").setUpperBound(0)
iAF1260_2.getReaction("R_GLCtex_copy1").setUpperBound(8)
iAF1260_2.getReaction("R_GLCtex_copy2").setUpperBound(0)

# Restrict the uptake of leucine in the knockout:
iAF1260.getReaction("R_LEUtex").setUpperBound(3)
# Restrict the lysine uptake in lysine knock out
iAF1260_2.getReaction("R_LYStex").setUpperBound(3)


community_model = CommunityModel(
    [iAF1260, iAF1260_2],
    ["R_BIOMASS_Ec_iAF1260_core_59p81M", "R_BIOMASS_Ec_iAF1260_core_59p81M"],
    ["dleu", "dlys"],
)
community_model.createObjectiveFunction

community_model.getReaction("R_IPPS_dleu").setUpperBound(0)
community_model.getReaction("R_DAPDC_dlys").setUpperBound(0)


# start_time = time.time()
# EndPointFBA(
#     community_model,
#     2,
#     {"dleu": 1.0, "dlys": 1.0},
#     {"M_glc__D_e": 24, "M_leu__L_e": 0, "M_lys__L_e": 0},
#     dt=0.1,
# )
# print(f"--- total time: {2} ------ {(time.time() - start_time)} seconds ---")
# print()

df = pd.read_csv(
    "/Users/stevenwijnen/Bioinformatica/year2/SysBio_lab/results/code/ecoli/fva_zero_reactions_iAF1260_community.csv"
)
fva_list = first_column_list = df.iloc[:, 0].tolist()
for rid in fva_list:
    community_model.deleteReactionAndBounds(rid)

community_model.deleteNonReactingSpecies()

start_time = time.time()
ep = EndPointFBA(
    community_model,
    20,
    {"dleu": 1.0, "dlys": 1.0},
    {"M_glc__D_e": 24, "M_leu__L_e": 0, "M_lys__L_e": 0},
    dt=0.1,
)
print(f"--- total time: {3} ------ {(time.time() - start_time)} seconds ---")
print(ep.simulate())

# start_time = time.time()
# EndPointFBA(
#     community_model,
#     4,
#     {"dleu": 1.0, "dlys": 1.0},
#     {"M_glc__D_e": 24, "M_leu__L_e": 0, "M_lys__L_e": 0},
#     dt=0.1,
# )

# print(f"--- total time: {4} ------ {(time.time() - start_time)} seconds ---")
# print()

# start_time = time.time()
# EndPointFBA(
#     community_model,
#     5,
#     {"dleu": 1.0, "dlys": 1.0},
#     {"M_glc__D_e": 24, "M_leu__L_e": 0, "M_lys__L_e": 0},
#     dt=0.1,
# )
# print(f"--- total time: {5} ------ {(time.time() - start_time)} seconds ---")
# print()

# start_time = time.time()
# EndPointFBA(
#     community_model,
#     6,
#     {"dleu": 1.0, "dlys": 1.0},
#     {"M_glc__D_e": 24, "M_leu__L_e": 0, "M_lys__L_e": 0},
#     dt=0.1,
# )
# print(f"--- total time: {6} ------ {(time.time() - start_time)} seconds ---")


# ScanUnusedReactions.scan_unused_reactions(
#     community_model, {"M_glc__D_e": 24, "M_leu__L_e": 0, "M_lys__L_e": 0}, True
# )

# start_time = time.time()
# ep = EndPointFBA(
#     community_model,
#     7,
#     {"dleu": 1.0, "dlys": 1.0},
#     {"M_glc__D_e": 24, "M_leu__L_e": 0, "M_lys__L_e": 0},
#     dt=0.1,
# )

# print(
#     "--- EndPointModel 7 time removed reactions------ %s seconds ---"
#     % (time.time() - start_time)
# )

# print(ep.simulate())

# community_model.deleteNonReactingSpecies()
# start_time = time.time()
# ep = EndPointFBA(
#     community_model,
#     7,
#     {"dleu": 1.0, "dlys": 1.0},
#     {"M_glc__D_e": 24, "M_leu__L_e": 0, "M_lys__L_e": 0},
#     dt=0.1,
# )
# print(
#     "--- EndPointModel 7 time no redunant reactions no non reacting species ------ %s seconds ---"
#     % (time.time() - start_time)
# )
# print(ep.simulate())
