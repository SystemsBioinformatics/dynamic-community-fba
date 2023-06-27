import cbmpy
import numpy as np
import matplotlib.pyplot as plt

from endPointFBA.dynamic_fba import dynamic_fba
from endPointFBA.jointFBA import create_joint_fba_model
from endPointFBA.KineticModel import KineticModel
from endPointFBA.CombinedModel import CombinedModel
from cbmpy.CBModel import Model
from endPointFBA.helpers.build_community_matrix import combine_models

model1 = cbmpy.loadModel("data/bigg_models/e_coli_core.xml")
model2: Model = cbmpy.loadModel("data/bigg_models/e_coli_core.xml")

combined_model = combine_models(
    [model1, model2],
    ["e_coli_1", "e_coli_2"],
    "R_BIOMASS_Ecoli_core_w_GAM",
)

cbmpy.saveModel(combined_model, "data/combined_e_coli_core.xml")
raise Exception
joint_model = create_joint_fba_model(
    combined_model,
    ["R_BIOMASS_Ecoli_core_w_GAM", "R_BIOMASS_Ecoli_core_w_GAM_ecoli_2"],
)

km: KineticModel = KineticModel(
    joint_model, {"R_GLCpts": [10, 5], "R_GLCpts_ecoli_2": [10, 5]}
)

ts = np.linspace(0, 15, 100)
y = dynamic_fba(km, "R_BIOMASS_Ecoli_core_w_GAM", ts, 0.1)

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

# Creating the second y-axis for y1 (M_gln__L_e)
# ax2 = ax1.twinx()
# ax2.plot(ts, y1, color="r")
# ax2.set_ylabel("M_gln__L_e", color="r")
# ax2.tick_params(axis="y", labelcolor="r")

# Creating the third y-axis for y3 (M_glc__D_e)
ax3 = ax1.twinx()
ax3.spines["right"].set_position(("outward", 60))
ax3.plot(ts, y3, color="g")
ax3.set_ylabel("M_glc__D_e", color="g")
ax3.tick_params(axis="y", labelcolor="g")

plt.show()
