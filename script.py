import cbmpy
import numpy
from cbmpy.CBModel import Model
from DCFBA.Models import CommunityModel
from DCFBA.DynamicModels import DynamicJointFBA

# Build the model and set all bounds
model1: Model = cbmpy.loadModel("data/bigg_models/e_coli_core.xml")
model1.getReaction("R_GLUt2r").setUpperBound(11)
model1.getReaction("R_PYRt2").setUpperBound(6)
model_2 = model1.clone()

model_2.getReaction("R_PYRt2").setUpperBound(10)
model_2.getReaction("R_GLUt2r").setUpperBound(8)

combined_model = CommunityModel(
    [model1, model_2],
    ["R_BIOMASS_Ecoli_core_w_GAM", "R_BIOMASS_Ecoli_core_w_GAM"],
    ["ecoli_1", "ecoli_2"],
)

combined_model.getReaction("R_EX_pyr_e").setLowerBound(-numpy.inf)

dynamic_fba = DynamicJointFBA(
    combined_model,
    [0.1, 0.1],
    {"M_pyr_e": 10, "M_glu__L_e": 0, "M_glc__D_e": 0},
)


def deviate_func(
    DFBA: DynamicJointFBA,
    used_time,
    condition: int,
) -> int:
    if (
        DFBA.m_metabolite_concentrations["M_pyr_e"][-1] <= 5.0
        and condition < 1
    ):
        combined_model.getReaction("R_EX_glu__L_e").setLowerBound(-numpy.inf)

        DFBA.m_metabolite_concentrations["M_glu__L_e"][-1] = 30

        return 1

    return condition


T, metabolites, biomasses, fluxes = dynamic_fba.simulate(
    0.1, deviate=deviate_func
)

import matplotlib.pyplot as plt


plt.figure(1)
plt.plot(T, metabolites["M_glu__L_e"], color="blue", label="[glu__L]")
plt.plot(T, metabolites["M_pyr_e"], color="orange", label="[pyr]")

plt.xlabel("Time")
plt.ylabel("Concentration")
plt.legend()

plt.figure(2)
plt.plot(T, biomasses["ecoli_1"], color="blue", label="Biomass model 1")
plt.plot(T, biomasses["ecoli_2"], color="orange", label="Biomass model 2")

plt.xlabel("Time")
plt.ylabel("Concentration")
plt.legend()

plt.show()
