from cobra.io import load_model
import cplex


def cobra_example():
    model = load_model("textbook")

    model.solver = "cplex"
    sum_two = model.problem.Constraint(
        model.reactions.FBA.flux_expression
        + model.reactions.NH4t.flux_expression,
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


def setproblemdata(p):
    p.objective.set_sense(p.objective.sense.maximize)

    p.linear_constraints.add(rhs=[20.0, 30.0], senses="LL")

    obj = [1.0, 2.0, 3.0]
    ub = [40.0, cplex.infinity, cplex.infinity]
    cols = [[[0, 1], [-1.0, 1.0]], [[0, 1], [1.0, -3.0]], [[0, 1], [1.0, 1.0]]]

    p.variables.add(
        obj=obj, ub=ub, columns=cols, names=["one", "two", "three"]
    )

    qmat = [
        [[0, 1], [-33.0, 6.0]],
        [[0, 1, 2], [6.0, -22.0, 11.5]],
        [[1, 2], [11.5, -11.0]],
    ]

    p.objective.set_quadratic(qmat)
    Q = cplex.SparseTriple(
        ind1=["one", "two", "three"], ind2=[0, 1, 2], val=[1.0] * 3
    )
    p.quadratic_constraints.add(rhs=1.0, quad_expr=Q, name="Q")


# Quadratic constrianed
def qcpex1():
    p = cplex.Cplex()
    setproblemdata(p)
    print(p.variables)
    raise Exception
    p.solve()

    # solution.get_status() returns an integer code
    print("Solution status = ", p.solution.get_status(), ":", end=" ")
    # the following line prints the corresponding string
    print(p.solution.status[p.solution.get_status()])
    print("Solution value  = ", p.solution.get_objective_value())

    numcols = p.variables.get_num()

    for j in range(numcols):
        print("Column %d:  Value = %10f" % (j, p.solution.get_values(j)))

    print(p.solution.get_linear_slacks(0))

    print()
    print("rhs    = ", p.quadratic_constraints.get_rhs("Q"))
    print("sense  = ", p.quadratic_constraints.get_senses(0))
    # range
    print("[lin]  = ", p.quadratic_constraints.get_linear_components(0, 0))
    # list
    print("[quad] = ", p.quadratic_constraints.get_quadratic_components([0]))
    print("[name] = ", p.quadratic_constraints.get_names())  # all as a list

    print()
    print("Objective function")
    print(" qnnz  = ", p.objective.get_num_quadratic_nonzeros())
    print(" quad  = ", p.objective.get_quadratic())
    print(" lin   = ", p.objective.get_linear())
    print(" sense = ", p.objective.get_sense())


if __name__ == "__main__":
    qcpex1()
