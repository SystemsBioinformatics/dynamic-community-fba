import cbmpy
import matplotlib.pyplot as plt
from cbmpy.CBModel import Model, Reaction
from DCFBA.ToyModels import model_a, model_b
from DCFBA.Models.CommunityModel import CommunityModel
from DCFBA.Models.Kinetics import KineticsStruct
from DCFBA.DynamicModels import DynamicJointFBA

model = cbmpy.loadModel("data/bigg_models/e_coli_core.xml")
model.getReaction("R_GLCpts").setUpperBound(10)

community_model = CommunityModel(
    [model], ["R_BIOMASS_Ecoli_core_w_GAM"], ["ecoli_1"]
)

kin = KineticsStruct({"R_GLCpts_ecoli_1": ["", 5, 10]})
dj = DynamicJointFBA(community_model, [0.1], {"M_glc__D_e": 10}, kinetics=kin)


def deviate_func(
    DFBA: DynamicJointFBA,
    used_time,
    condition: int,
) -> int:
    if (
        DFBA.m_metabolite_concentrations["M_glc__D_e"][-1] <= 5.0
        and condition < 1
    ):
        DFBA.m_metabolite_concentrations["M_glu__L_e"][-1] = 30

        return 1

    return condition


T, metabolites, biomasses, _ = dj.simulate(0.1, deviate=deviate_func)


ax = plt.subplot(111)
ax.plot(T, biomasses["ecoli_1"])
ax2 = plt.twinx(ax)
ax2.plot(T, metabolites["M_glu__L_e"], color="r")

ax.set_ylabel("Biomass", color="b")
ax2.set_ylabel("Glucose", color="r")

plt.show()
