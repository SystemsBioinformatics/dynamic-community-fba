import cbmpy
import matplotlib.pyplot as plt
from endPointFBA.dynamic_fba import dynamic_fba
from endPointFBA.CombinedModel import CombinedModel
from endPointFBA.KineticModel import KineticModel

from cbmpy.CBModel import Model
import numpy as np

model1: Model = cbmpy.loadModel("data/bigg_models/e_coli_core.xml")


for reaction_id in model1.getReactionIds():
    if reaction_id.startswith("R_EX"):
        reaction = model1.getReaction(reaction_id)
        reaction.is_exchange = True

model1.getReaction("R_ATPM").setLowerBound(0)
km_1: KineticModel = KineticModel(model1, {"R_GLCpts": [10, 5]})

ts = np.linspace(0, 15, 100)
y = dynamic_fba(km_1, "R_BIOMASS_Ecoli_core_w_GAM", ts, 0.1)

# Assuming y is the DataFrame containing the data for y1, y2, and y3
y1 = y["M_gln__L_e"]
y2 = y["biomass"]
y3 = y["M_glc__D_e"]
ts = ts[0 : len(y1)]

print(y3)
print(y2)
fig, ax1 = plt.subplots()
# Plotting y2 (biomass) on the first y-axis
ax1.plot(ts, y2, color="b")
ax1.set_xlabel("Time")
ax1.set_ylabel("Biomass", color="b")
ax1.tick_params(axis="y", labelcolor="b")

# Creating the third y-axis for y3 (M_glc__D_e)
ax3 = ax1.twinx()
ax3.spines["right"].set_position(("outward", 60))
ax3.plot(ts, y3, color="g")
ax3.set_ylabel("M_glc__D_e", color="g")
ax3.tick_params(axis="y", labelcolor="g")

plt.show()
