import cbmpy
import numpy as np
import matplotlib.pyplot as plt

from endPointFBA.dynamic_fba import dynamic_fba
from endPointFBA.jointFBA import create_joint_fba_model
from endPointFBA.KineticModel import KineticModel
from endPointFBA.CombinedModel import CombinedModel
from cbmpy.CBModel import Model

model1: Model = cbmpy.loadModel("data/bigg_models/e_coli_core.xml")
model2: Model = cbmpy.loadModel("data/bigg_models/e_coli_core.xml")
r1 = model1.getReaction("R_ATPM")
r2 = model2.getReaction("R_ATPM")

# r1.setLowerBound(0)
# r2.setLowerBound(0)
combined_model = CombinedModel(
    [model1, model2],
    ["e_coli_1", "e_coli_2"],
    "R_BIOMASS_Ecoli_core_w_GAM_e_coli_1",
)


km_1: KineticModel = KineticModel(
    combined_model, {"R_GLCpts_e_coli_1": [10, 5]}
)


ts = np.linspace(0, 15, 100)
y = dynamic_fba(km_1, "R_BIOMASS_Ecoli_core_w_GAM_e_coli_1", ts, 0.1)

# Assuming y is the DataFrame containing the data for y1, y2, and y3
y1 = y["M_gln__L_e"]
y2 = y["biomass"]
y3 = y["M_glc__D_e"]
ts = ts[0 : len(y1)]

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


joint_model = create_joint_fba_model(
    combined_model,
    [
        "R_BIOMASS_Ecoli_core_w_GAM_e_coli_1",
        "R_BIOMASS_Ecoli_core_w_GAM_e_coli_2",
    ],
    create_new_model=False,
)


km_2: KineticModel = KineticModel(joint_model, {"R_GLCpts_e_coli_1": [10, 5]})

y = dynamic_fba(km_2, "Xcomm", ts, 0.1)

fig2, ax2 = plt.subplots()

y1 = y["M_gln__L_e"]
y2 = y["biomass"]
y3 = y["M_glc__D_e"]
ts = ts[0 : len(y1)]

# Plotting y2 (biomass) on the first y-axis
ax2.plot(ts, y2, color="b")
ax2.set_xlabel("Time")
ax2.set_ylabel("Biomass", color="b")
ax2.tick_params(axis="y", labelcolor="b")

# Creating the third y-axis for y3 (M_glc__D_e)
a4 = ax2.twinx()
a4.spines["right"].set_position(("outward", 60))
a4.plot(ts, y3, color="g")
a4.set_ylabel("M_glc__D_e", color="g")
a4.tick_params(axis="y", labelcolor="g")


plt.show()
