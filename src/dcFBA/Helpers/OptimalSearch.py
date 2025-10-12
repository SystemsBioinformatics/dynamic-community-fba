# Binary search on answer
import numpy as np
from ..DynamicModels import EndPointFBA
from ..Models import CommunityModel

# Remember which numbers were visited
# visited: dict[int, float] = {}
visited: dict[tuple, float] = {}


def time_search(
    cm: CommunityModel,
    initial_biomasses: dict[str, float],
    initial_concentrations: dict[str, float] = None,
    dt=0.1,
    set_values: tuple[float, int] = None,
    use_cache: bool = True,
) -> tuple[float, float]:
    """Finds the lowest number of time points given initial values and a dt

    Args:
        cm (CommunityModel): Genome-scale metabolic model of the community
        initial_biomasses (dict[str, float]): Initial biomass value for each cell population
        initial_concentrations (dict[str, float], optional): Initial metabolite concentrations. Defaults to {}.
        dt (float, optional): Time step. Defaults to 0.1.
        set_values (tuple[float, float], optional): Set to identify when a specific value is attained. 
            [value, N]. Where value is the objective value you want to reach 
            and N the initial guess of minimum number of time points to achieve it.
            Defaults to None
        use_cache (bool, optional): Whether to use global cache. Defaults to True.

    Returns:
        tuple[float, float]: number of time points and the reached objective value
    """

    if initial_concentrations is None:
        initial_concentrations = {}

    def simulate_with_cache(n: int) -> float:
        """Helper to simulate with optional caching"""
        cache_key = _make_cache_key(cm, n, initial_biomasses, initial_concentrations, dt)
        
        if use_cache and cache_key in visited:
            return visited[cache_key]
        
        ep = EndPointFBA(cm, n, initial_biomasses, initial_concentrations, dt=dt)
        value = ep.simulate()
        
        if use_cache:
            visited[cache_key] = value
        
        return value

    low = 1
    if set_values is None:
        high = find_upper_bound(cm, initial_biomasses, initial_concentrations, dt, 
                                simulate_with_cache)
        obj = simulate_with_cache(high)
    else:
        obj = set_values[0]
        high = set_values[1]
        simulate_with_cache(high)  # Ensure it's cached

    while low < high:
        n = (low + high) // 2
        print(f"Trying {n} ...")
        value = simulate_with_cache(n)

        if round(value, 5) >= round(obj, 5):
            high = n
            if not set_values:
                obj = value
        else:
            low = n + 1

    final_value = simulate_with_cache(high) #

    # if set_values and value < set_values[0]:
    if set_values and final_value < set_values[0]:
        print("WARNING: Set objective can not be reached")

    return (high, final_value)


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
            old_udc = ep.m_model.getObject(f"biomass_fraction_{mid}_{ep.m_times[-1]}")

            # for cid in old_udc.getConstraintComponentIDs():
            #     ep.m_model.unRegisterObjectInGlobalStore(cid)
            old_udc.constraint_components = []

            ep.m_model.unRegisterObjectInGlobalStore(old_udc.getId())

            value = objective * mid_point
            ep.m_model.setReactionBounds("X_comm", value - 0.001, value + 0.001)
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
    dt: float,
    simulate_func,
) -> int:
    """Find upper bound by doubling until convergence"""
    n = 1
    prev_value = 0
    while True:
        n *= 2  # Double the value of n
        # ep = EndPointFBA(cm, n, initial_biomasses, initial_concentrations, dt=dt)
        # current_value = ep.simulate()
        current_value = simulate_func(n)

        # Check if current value is NaN or if it doesn't increase from the previous value
        if np.isnan(current_value) or current_value <= prev_value:
            return n // 2

        # visited[n] = current_value
        prev_value = current_value


def clearvisited():
    """Clear the global simulation cache. Useful between different experiments."""
    visited.clear()


def _make_cache_key(
    cm: CommunityModel,
    n: int,
    initial_biomasses: dict[str, float],
    initial_concentrations: dict[str, float],
    dt: float,
) -> tuple:
    """Create a hashable cache key from simulation parameters
    
    Uses model ID + composition to ensure different communities don't share cache.
    """
    # Create identifier from model ID and its composition
    cm_identifier = (
        cm.getId(),
        tuple(cm.single_model_ids),  # The actual models that make up this community
        tuple(cm.single_model_biomass_reaction_ids),  # Their biomass reactions
    )
    
    biomass_tuple = tuple(sorted(initial_biomasses.items()))
    conc_tuple = tuple(sorted(initial_concentrations.items()))
    return (cm_identifier, n, biomass_tuple, conc_tuple, dt)
