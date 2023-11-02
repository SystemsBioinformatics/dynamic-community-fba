import matplotlib.pyplot as plt
from cbmpy.CBModel import Model, Reaction
from dcFBA.ToyModels import model_a, model_b
from dcFBA.Models.CommunityModel import CommunityModel
from dcFBA.Models.KineticsStruct import KineticsStruct
from dcFBA.DynamicModels import EndPointFBA
from dcFBA.Helpers.PlotsEndPointFBA import plot_biomasses, plot_metabolites

m_a: Model = model_a.build_toy_model_fba_A()

m_a.getReaction("R_1").setUpperBound(10)
m_a.getReaction("R_4").setUpperBound(3)
m_a.getReaction("R_6").setUpperBound(1)

# Need to delete biomass model, since it is created in the endPoint model
reaction: Reaction = m_a.getReaction("R_BM_A").deleteReagentWithSpeciesRef(
    "BM_e_A"
)

# m_a.getReaction("R_import_B").setUpperBound(1)

m_a.getReaction("R_1").setLowerBound(0)
m_a.getReaction("R_4").setLowerBound(0)
m_a.getReaction("R_6").setLowerBound(0)

m_b: Model = model_b.build_toy_model_fba_B()
m_b.getReaction("R_1").setUpperBound(10)
m_b.getReaction("R_3").setUpperBound(1)
m_b.getReaction("R_5").setUpperBound(1)
# Need to delete biomass model, since it is created in the endPoint model
reaction: Reaction = m_b.getReaction("R_BM_B").deleteReagentWithSpeciesRef(
    "BM_e_B"
)


m_b.getReaction("R_1").setLowerBound(0)
m_b.getReaction("R_3").setLowerBound(0)
m_b.getReaction("R_5").setLowerBound(0)


kin = KineticsStruct(
    {
        "R_1_modelA": ["S_e", 10, 10],
        "R_4_modelA": ["B_e", 5, 3],
        "R_6_modelA": ["B_e", 3, 1],
        # B
        "R_1_modelB": ["S_e", 10, 10],
        "R_3_modelB": ["A_e", 2, 1],
    }
)


community_model = CommunityModel(
    [m_a, m_b], ["R_BM_A", "R_BM_B"], ["modelA", "modelB"]
)

# Are set by the model as often required by other GSMM
community_model.deleteReactionAndBounds("BM_e_A_exchange")
community_model.deleteReactionAndBounds("BM_e_B_exchange")

# community_model.getReaction("B_exchange").setLowerBound(-100)
# community_model.getReaction("A_exchange").setLowerBound(-100)


# community_model.getReaction("S_exchange").setLowerBound(-100)

n = 21
ep = EndPointFBA(
    community_model,
    n,
    {"modelA": 1.0, "modelB": 2.0},
    {"S_e": 100, "A_e": 0.0, "B_e": 0.0},
    0.1,
)


solution = ep.simulate()


biomasses = ep.get_biomasses()
metabolites = ep.get_metabolites()
fluxes = ep.get_fluxes()
time = ep.get_time_points()

plt.plot(
    range(len(time) + 1),
    metabolites["S_e"],
    color="blue",
    label="Metabolite s",
)
plt.xlabel("Time")
plt.ylabel("Concentration")
plt.legend()
plt.show()


plt.plot(
    range(len(time) + 1),
    metabolites["A_e"],
    color="blue",
    label="Metabolite a",
)
plt.xlabel("Time")
plt.ylabel("Concentration")
plt.legend()
plt.show()


# Adding labels and title

plt.plot(
    range(len(time) + 1),
    metabolites["B_e"],
    color="blue",
    label="Metabolite B",
)
plt.xlabel("Time")
plt.ylabel("Concentration")
# Adding legend
plt.legend()
plt.show()

plt.plot(
    range(len(time) + 1), biomasses["modelA"], color="blue", label="modelA"
)
plt.plot(
    range(len(time) + 1), biomasses["modelB"], color="orange", label="modelB"
)
plt.xlabel("Time")
plt.ylabel("Concentration")
plt.legend()
plt.show()
