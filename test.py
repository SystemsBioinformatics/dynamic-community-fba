import cplex

c = cplex.Cplex()

# target = c.parameters.optimalitytarget.values
# c.parameters.optimalitytarget.set(target.optimal_global)
c.objective.set_sense(c.objective.sense.minimize)

n = 3
c.variables.add(
    names=["x" + str(i) for i in range(n)],
    types=[c.variables.type.continuous] * n,
    # ub=[10, 20],
    # lb=[-10, -10],
    # obj=[0.0, -1.0],
)

qmat = []

qmat = [
    [["x0", "x1", "x2"], [2.0, -2.0, 0.0]],
    [["x0", "x1", "x2"], [-2.0, 4.0, -2.0]],
    [["x0", "x1", "x2"], [0.0, -2.0, 2.0]],
]


c.objective.set_quadratic(qmat)

print(c.objective.get_quadratic())
# c.objective.set_quadratic_coefficients([("x1", "x1", -1)])

c.solve()

print()
print("Objective function")
print(" qnnz  = ", c.objective.get_num_quadratic_nonzeros())
print(" quad  = ", c.objective.get_quadratic())
print(" lin   = ", c.objective.get_linear())
print(" sense = ", c.objective.get_sense())

print(c.solution.status[c.solution.get_status()])
print("Solution value  = ", c.solution.get_objective_value())

x = c.solution.get_values()
print(x)
print("x0 =", x[0])
print("x1 =", x[1])


# print("x3 =", x[2])
# print("x4 =", x[3])
