import cbmpy
from cbmpy.CBModel import Model, Objective

core = cbmpy.loadModel("./models/bigg_models/e_coli_core.xml")


model: Model = cbmpy.loadModel("models/bigg_models/e_coli_core.xml")
model.__FBC_VERSION__ = 3

udc = model.createUserDefinedConstraint(
    f"crazyyy_lb",
    2.0,
    2.0,
    components=[
        (1, "R_FBA", "linear"),
        (1, "R_NH4t", "linear"),
    ],
)
model.addUserDefinedConstraint(udc)

# model.createUser("sum2", [[1, "R_FBA"], [1, "R_NH4t"]], "=", 2)
model.createObjectiveFunction("R_FBA", -1, "minimize")

obj: Objective = model.getActiveObjective()

QP = [[1.0, "R_NH4t", "R_NH4t"], [1.0, "R_FBA", "R_FBA"]]
obj.createQuadraticFluxObjectives(QP)

sol = cbmpy.doFBA(model)
print(sol)

FBAsol = model.getSolutionVector(names=True)
FBAsol = dict(zip(FBAsol[1], FBAsol[0]))
print(FBAsol["R_FBA"])
print(FBAsol["R_NH4t"])
