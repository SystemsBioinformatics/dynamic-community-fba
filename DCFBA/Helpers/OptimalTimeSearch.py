# Binary search on answer

import numpy as np
from ..DynamicModels import EndPointFBA
from ..Models import CommunityModel

# Remember which numbers were visited
visited: dict[int, float] = {}


# Find the smallest N for which the objective function is maximal
def search(
    cm: CommunityModel,
    initial_biomasses: dict[str, float],
    initial_concentrations: dict[str, float] = {},
    dt=0.1,
):
    n = np.random.randint(10, 100)

    low = 1
    obj, high = find_max_solution(
        n, cm, initial_biomasses, initial_concentrations, dt
    )

    while low < high:
        n = (low + high) // 2
        print("trying " + str(n) + ".....")

        if n in visited.keys():
            value = visited[n]
        else:
            ep = EndPointFBA(
                cm, n, initial_biomasses, initial_concentrations, dt
            )
            value = ep.simulate()
            visited[n] = value

        if round(value, 5) >= obj:
            high = n
        elif round(value, 5) < obj:
            low = n + 1

    return high


def find_max_solution(
    n,
    cm: CommunityModel,
    initial_biomasses: dict[str, float],
    initial_concentrations: dict[str, float],
    dt,
):
    prev = 0
    ep = EndPointFBA(cm, n, initial_biomasses, initial_concentrations, dt)
    value = ep.simulate()
    while round(prev, 5) < round(value, 5):
        prev = value
        n = 2 * n
        ep = EndPointFBA(cm, n, initial_biomasses, initial_concentrations, dt)
        value = ep.simulate()
        visited[n] = value

    return [value, n // 2]
