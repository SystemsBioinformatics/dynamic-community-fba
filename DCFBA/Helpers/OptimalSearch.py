# Binary search on answer

import numpy as np
from ..DynamicModels import EndPointFBA
from ..Models import CommunityModel
import cbmpy

# Remember which numbers were visited
visited: dict[int, float] = {}


def time_search(
    cm: CommunityModel,
    initial_biomasses: dict[str, float],
    initial_concentrations: dict[str, float] = {},
    dt=0.1,
):
    low = 1
    high = find_upper_bound(cm, initial_biomasses, initial_concentrations, dt)
    input(high)

    ep = EndPointFBA(
        cm, high, initial_biomasses, initial_concentrations, dt=dt
    )
    obj = ep.simulate()

    while low < high:
        n = (low + high) // 2
        print(f"Trying {n} ...")
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


# TODO when there is a remove UserDefinedConstraint fix this
def balance_search(
    cm: CommunityModel,
    n,
    initial_concentrations: dict[str, float],
    dt,
    X_initial,
    objective,
    epsilon=0.01,
):
    low = 1
    high = objective

    while high - low > epsilon:
        mid = (low + high) / 2
        print(f"Trying {mid} ...")
        ep = EndPointFBA(
            cm,
            n,
            {},
            initial_concentrations,
            dt,
        )
        ep.balanced_growth(X_initial, objective / mid)
        solution = ep.simulate()

        if not np.isnan(solution):  # If solution is not NaN
            high = mid
        else:
            low = mid

    return high  # Return the lowest n value for which the solution is not NaN


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
