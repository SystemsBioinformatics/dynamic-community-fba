import cbmpy
import re
import matplotlib.pyplot as plt
from cbmpy.CBModel import Model, Reaction, Objective
from DCFBA.ToyModels import model_a, model_b
from DCFBA.Models.CommunityModel import CommunityModel
from DCFBA.Models.Kinetics import KineticsStruct
from DCFBA.DynamicModels import EndPointFBA
from DCFBA.Helpers.PlotsEndPointFBA import (
    plot_biomasses,
    plot_metabolites,
    plot_fluxes,
)
from testingQP import cplex_constructPorbfromFBA
import cbmpy.CBCPLEX as wcplex

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
        "R_1_modelA": ["", 10, 10],
        "R_4_modelA": ["B_e", 5, 3],
        "R_6_modelA": ["B_e", 3, 1],
        # B
        "R_1_modelB": ["", 10, 10],
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

n = 27
dj = EndPointFBA(
    community_model,
    n,
    {"modelA": 1.0, "modelB": 1.0},
    {"S_e": 100, "A_e": 0.0, "B_e": 0.0},
    0.1,
)

# solution = dj.simulate()
# plot_fluxes(dj, ["R_1_modelA", "R_1_modelB", "R_2_modelA"])

# plot_biomasses(dj)
# plot_metabolites(dj, {"S_e": 100, "A_e": 0.0, "B_e": 0.0})

obj: Objective = dj.m_model.getActiveObjective()
obj.setOperation("min")
obj.QPObjective = []

for i, tid in enumerate(dj.m_times[:-1]):
    rids = dj.m_model.getReactionIds(tid)
    for rid in rids:
        reaction = dj.m_model.getReaction(rid)

        # R_ are all the original reqactions
        # We dont want the exchange reactions
        if rid.startswith("R_") and not reaction.is_exchange:
            rid = re.match(r"(.*?)_(time\d*)", rid).group(1)
            rid_t_t1 = f"{rid}_{dj.m_times[i+1]}"
            rid_t_t0 = f"{rid}_{dj.m_times[i]}"

            # Min sum (rid_time_t+1 - rid_time_t )^2
            obj.QPObjective.append(([rid_t_t0, rid_t_t0], 1.0))
            obj.QPObjective.append(([rid_t_t0, rid_t_t1], -2.0))
            obj.QPObjective.append(([rid_t_t1, rid_t_t1], 1.0))

# print(obj.QPObjective)

# solution = dj.simulate()
dj.m_model.getReaction("X_comm").setLowerBound(12.77)
dj.m_model.getReaction("X_comm").setUpperBound(12.778)

dj.m_model.buildStoichMatrix()
prob = cplex_constructPorbfromFBA(dj.m_model)
prob.solve()
# wcplex.cplx_setFBAsolutionToModel(dj.m_model, prob)
sol, objname, objval = wcplex.cplx_getOptimalSolution(prob)
(
    dj.m_model.objectives[dj.m_model.activeObjIdx].solution,
    dj.m_model.objectives[dj.m_model.activeObjIdx].value,
) = (sol, objval)

for r in dj.m_model.reactions:
    rid = r.getId()
    if rid in sol:
        r.value = sol[rid]
    else:
        r.value = None

print("Solution status = ", prob.solution.get_status(), ":", end=" ")

# x = prob.solution.get_values()

# numcols = prob.variables.get_num()
# for j in range(numcols):
#     print("Column %s:  Value = %10f" % (prob.variables.get_names(j), x[j]))

print(dj.m_model.getSolutionVector(names=True))

plot_fluxes(dj, ["R_1_modelA", "R_1_modelB", "R_2_modelA"])


# plot_biomasses(dj, n)
# plot_metabolites(dj, {"S_e": 100, "A_e": 0.0, "B_e": 0.0}, n)
