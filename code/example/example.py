import cbmpy


def add_one(number: int) -> int:
    return number + 1


def multiply_by_two(number: int) -> int:
    return number * 2


def cbmpy_example_model():
    cmod = cbmpy.readSBML3FBC("cbmpy_test_core")
    tt = cbmpy.doFBA(cmod)
    return tt
