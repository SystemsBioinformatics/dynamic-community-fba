import cbmpy
from cbmpy.CBModel import Model

model: Model = cbmpy.readSBML2FBA("cbmpy_test_ecoli.xml")
print(model)
