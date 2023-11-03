# Binary search on answer

import numpy as np
from ..DynamicModels import EndPointFBA
from ..Models import CommunityModel

# Remember which numbers were visited
visited: dict[int, float] = {}


def time_search(
    cm: CommunityModel,
    initial_biomasses: dict[str, float],
    initial_concentrations: dict[str, float] = {},
    dt=0.1,
    set_values: tuple[float, float] = None,
) -> list[float, float]:
    """Finds the lowest number of time points given initial values and a dt


    Args:
        cm (CommunityModel): _description_
        initial_biomasses (dict[str, float]): _description_
        initial_concentrations (dict[str, float], optional): _description_. Defaults to {}.
        dt (float, optional): _description_. Defaults to 0.1.
        set_values (tuple[float, float], optional): Set to identify when a
            specific value is attained. [value, N]. Where value is the value
            you want to reach and N the initial guess of N.
            Defaults to None


    Returns:
        ist[float, float]: number of time points and the reached objective
            value
    """
    low = 1
    if set_values is None:
        high = find_upper_bound(
            cm, initial_biomasses, initial_concentrations, dt
        )
        ep = EndPointFBA(
            cm, high, initial_biomasses, initial_concentrations, dt=dt
        )
        obj = ep.simulate()

    else:
        obj = set_values[0]
        high = set_values[1]

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
            if not set_values:
                obj = value
        elif round(value, 5) < obj:
            low = n + 1

    if set_values and value < set_values[0]:
        print("WARNING: Set objective can not be reached")

    return [high, visited[high]]


# TODO when there is a remove UserDefinedConstraint fix this
def balance_search_clean(
    community_model,
    n,
    initial_concentrations,
    dt,
    X_initial,
    objective,
    epsilon=0.01,
):
    low = 0
    high = 1

    while high - low > epsilon:
        mid_point = (low + high) / 2
        print(f"Trying {mid_point} ...")
        ep = EndPointFBA(community_model, n, {}, initial_concentrations, dt)
        ep.balanced_growth(X_initial, objective * mid_point)
        solution = ep.simulate()

        if not np.isnan(solution):  # If solution is not NaN
            low = mid_point
        else:
            high = mid_point

    return low  # Return the n closest to 1 for which the solution is not NaN


# TODO Clean this code once deleteUserdefinedConstraint is implemented in cbmpy
def balanced_search_quick(ep: EndPointFBA, X_initial, objective, epsilon=0.01):
    """Run this one if you know what you are doing, quick and dirty solution"""
    low = 0
    high = 1
    ep.balanced_growth(X_initial, objective)
    solution = ep.simulate()
    if not np.isnan(solution):
        return 1.0
    while high - low > epsilon:
        mid_point = (low + high) / 2
        print(f"Trying {mid_point} ...")
        for mid, _ in ep.m_model.get_model_biomass_ids().items():
            # Should be replaced with remove user defined constraint once this is
            #  implemented
            old_udc = ep.m_model.getObject(
                f"biomass_fraction_{mid}_{ep.m_times[-1]}"
            )

            # for cid in old_udc.getConstraintComponentIDs():
            #     ep.m_model.unRegisterObjectInGlobalStore(cid)
            old_udc.constraint_components = []

            ep.m_model.unRegisterObjectInGlobalStore(old_udc.getId())

            value = objective * mid_point
            ep.m_model.setReactionBounds(
                "X_comm", value - 0.001, value + 0.001
            )
            udc = ep.m_model.createUserDefinedConstraint(
                f"biomass_fraction_{mid}_{ep.m_times[-1]}",
                0.0,
                0.0,
                components=[
                    (1.0, f"BM_{mid}_exchange_final", "linear"),
                    (
                        -1.0 * (value + X_initial),
                        f"Phi_{mid}",
                        "linear",
                    ),
                ],
            )

            ep.m_model.addUserDefinedConstraint(udc)
        solution = ep.simulate()

        if not np.isnan(solution):  # If solution is not NaN
            low = mid_point
        else:
            high = mid_point

    return low  # Return the n closest to 1 for which the solution is not NaN


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
