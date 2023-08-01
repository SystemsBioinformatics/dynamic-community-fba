# Binary search on answer

import numpy as np
from ..DynamicModels import EndPointFBA
from ..Models import CommunityModel

# TODO improve if remove nodes and add nodes is implemented in EndPointFBA
# class


# Find the smallest N for which the objective function is maximal
def search(
    cm: CommunityModel,
    initial_biomasses: dict[str, float],
    initial_concentrations: dict[str, float] = {},
    dt=0.1,
):
    n = 100  # np.random(0, 10000)

    low = 1
    obj, high = find_upper_bound(
        n, cm, initial_biomasses, initial_concentrations, dt
    )

    while low < high:
        n = (low + high) // 2
        print("trying " + str(n) + ".....")
        ep = EndPointFBA(cm, n, initial_biomasses, initial_concentrations, dt)
        value = ep.simulate()

        if value >= obj:
            high = n
        elif value < obj:
            low = n + 1

    return high


def find_upper_bound(
    n,
    cm: CommunityModel,
    initial_biomasses: dict[str, float],
    initial_concentrations: dict[str, float],
    dt,
):
    prev = 0
    ep = EndPointFBA(cm, n, initial_biomasses, initial_concentrations, dt)
    value = ep.simulate()
    while prev < value:
        prev = value
        n = 2 * n
        ep = EndPointFBA(cm, n, initial_biomasses, initial_concentrations, dt)
        value = ep.simulate()

    return [value, n // 2]
