from cobra.io import load_model

model = load_model("textbook")
model.solver = "cplex"
sum_two = model.problem.Constraint(
    model.reactions.FBA.flux_expression + model.reactions.NH4t.flux_expression,
    lb=2,
    ub=2,
)
model.add_cons_vars(sum_two)

quadratic_objective = model.problem.Objective(
    0.5 * model.reactions.NH4t.flux_expression**2
    + 0.5 * model.reactions.FBA.flux_expression**2
    - model.reactions.FBA.flux_expression,
    direction="min",
)
model.objective = quadratic_objective
solution = model.optimize(objective_sense=None)

print(solution.fluxes["NH4t"], solution.fluxes["FBA"])


print(model.reactions.NH4t.lower_bound)
print(model.reactions.FBA.lower_bound)
print(solution)
