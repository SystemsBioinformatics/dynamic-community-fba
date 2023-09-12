# Binary search on answer

import numpy as np
from ..DynamicModels import EndPointFBA
from ..Models import CommunityModel
import cbmpy

# Remember which numbers were visited
visited: dict[int, float] = {}

# TODO if nan!!!!!!!!!!!!!!!!


# Find the smallest N for which the objective function is maximal
def search(
    cm: CommunityModel,
    initial_biomasses: dict[str, float],
    initial_concentrations: dict[str, float] = {},
    dt=0.1,
):
    low = 1
    high = find_upper_bound(cm, initial_biomasses, initial_concentrations, dt)
    obj = 0

    while low < high:
        n = (low + high) // 2

        if n in visited.keys():
            value = visited[n]
        else:
            ep = EndPointFBA(
                cm, n, initial_biomasses, initial_concentrations, dt=dt
            )
            value = ep.simulate()
            visited[n] = value

        if round(value, 5) >= obj:
            high = n
            obj = value
        elif round(value, 5) < obj:
            low = n + 1

    return [high, obj]


def find_upper_bound(
    cm: CommunityModel,
    initial_biomasses: dict[str, float],
    initial_concentrations: dict[str, float],
    dt,
):
    n = 1
    prev_value = 0
    while True:
        n *= 2  # Double the value of n
        ep = EndPointFBA(
            cm, n, initial_biomasses, initial_concentrations, dt=dt
        )
        current_value = ep.simulate()
        # Check if current value is NaN or if it doesn't increase from the previous value
        if np.isnan(current_value) or current_value <= prev_value:
            return n // 2

        visited[n] = current_value
        prev_value = current_value
