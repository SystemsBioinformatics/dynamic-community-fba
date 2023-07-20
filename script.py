import cbmpy
import numpy
from cbmpy.CBModel import Model
from DCFBA.Models import CommunityModel
from DCFBA.DynamicModels import DynamicJointFBA

model1: Model = cbmpy.loadModel("data/bigg_models/e_coli_core.xml")
model1.getReaction("R_GLCpts").setUpperBound(10)

model_2 = model1.clone()
# Ecoli 2 imports glucose slower
model_2.getReaction("R_GLCpts").setUpperBound(8)

combined_model = CommunityModel(
    [model1, model_2],
    ["R_BIOMASS_Ecoli_core_w_GAM", "R_BIOMASS_Ecoli_core_w_GAM"],
    ["ecoli_1", "ecoli_2"],
)  # Create a CommunityModel of two  E. coli strains competing for resources


# Set the lower bound so the model does not eat everything at onec
combined_model.getReaction("R_PYRt2_ecoli_1").setUpperBound(3)
combined_model.getReaction("R_PYRt2_ecoli_2").setUpperBound(10)

combined_model.getReaction("R_EX_glc__D_e").setLowerBound(-numpy.inf)
combined_model.getReaction("R_EX_pyr_e").setLowerBound(-numpy.inf)

for exchange in combined_model.getExchangeReactionIds():
    reaction = combined_model.getReaction(exchange)
    reaction.setUpperBound(numpy.inf)

dynamic_fba = DynamicJointFBA(
    combined_model,
    [0.1, 0.1],
    {"M_pyr_e": 0, "M_glc__D_e": 10},
)


def deviate_func(
    DFBA: DynamicJointFBA,
    used_time,
    condition: int,
) -> None:
    if (
        DFBA.m_metabolite_concentrations["M_glc__D_e"][-1] <= 5
        and condition < 1
    ):
        # ACtivate exchange
        combined_model.getReaction("R_EX_pyr_e").setLowerBound(-numpy.inf)

        DFBA.m_metabolite_concentrations["M_pyr_e"][-1] += 10

        return 1

    return condition


T, metabolites, biomasses, fluxes = dynamic_fba.simulate(
    0.1, deviate=deviate_func
)


import matplotlib.pyplot as plt


plt.plot(T, metabolites["M_glc__D_e"], color="blue", label="[Glucose]")
plt.plot(T, metabolites["M_pyr_e"], color="orange", label="[L-Glutamate]")

plt.xlabel("Time")
plt.ylabel("Concentration")
plt.legend()
plt.show()


plt.plot(T, biomasses["ecoli_1"], color="blue", label="Biomass model 1")
plt.plot(T, biomasses["ecoli_2"], color="orange", label="Biomass model 2")

# Adding labels and title
plt.xlabel("Time")
plt.ylabel("Concentration")


# Adding legend
plt.legend()

# Displaying the plot
plt.show()
