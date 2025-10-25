"""
Utility functions for EndPointFBA and general dcFBA postprocessing.

Example usage:
    from dcFBA.Helpers import (
        set_metabolites,
        calculate_growth_rates,
        calculate_biomass_yields,
        calculate_biomass_fractions,
    )

    set_metabolites(ep, keep_zero_metabolites=False)
    growth_rates = calculate_growth_rates(ep.biomasses, dt=0.5)
"""

import math
import numpy as np
from dcFBA.DynamicModels.EndPointFBA import EndPointFBA

#################################################################################
# set_metabolites

def set_metabolites(ep, keep_zero_metabolites=True, metabolites_to_keep=None) -> None:
    """
    Public function to set the metabolite concentrations.
    Wrapper for the private method _set_metabolites of EndPointFBA.

    Args:
        ep (EndPointFBA): EndPointFBA model object.
        keep_zero_metabolites (bool, optional): If True, metabolites that are always zero are kept. Default True.
        metabolites_to_keep (list, optional): List of metabolite IDs to keep even if they are always zero
                                              (only used if keep_zero_metabolites == False). Default None.
    """
    if metabolites_to_keep is None:
        metabolites_to_keep = []

    ep._set_metabolites(keep_zero_metabolites, metabolites_to_keep)



#################################################################################
# round_values

def round_values(dct, precision=1e-6):
    """
    Rounds all float values in a dictionary to the specified precision.

    Args:
        dct (dict): The input dictionary with float values.
        precision (float): The desired precision for rounding (default is 1e-6).

    Returns:
        dict: A new dictionary with rounded float values.
    """
    
    digits = len(str(int((1 / precision) - 1)))
    rounded_dct = {}
    for key, value in dct.items():
        if isinstance(value, list):
            # Round each value in the list
            rounded_dct[key] = [round(v, digits) for v in value]
        else:
            # Round the single value
            rounded_dct[key] = round(value, digits)
    return rounded_dct

#################################################################################
# calculate_growth_rates

def calculate_growth_rates(biomasses, dt):
    """
    Calculates growth rates for each model based on biomass values.

    Args:
        biomasses (dict): A dictionary where keys are model IDs and values are lists of biomass values.
        dt (float): Time interval between biomass measurements.

    Returns:
        dict: A dictionary where keys are model IDs and values are lists of growth rates.
    """

    # Initialize a dictionary where to store the growth rates
    growth_rates = {}

    # Calculate growth rates, for each model at each timpoint
    for model_id, biomass_values in biomasses.items():
        growth_rates[model_id] = [
            math.log(biomass_values[t + 1] / biomass_values[t]) / dt
            for t in range(len(biomass_values) - 1)
        ]
    
    return growth_rates

#################################################################################
# calculate_biomass_yields

def calculate_biomass_yields(growth_rates, uptake_rates):
    """
    Calculates biomass yields given growth rates and uptake rates.

    Args:
        growth_rates (dict): A dictionary where keys are model IDs and values are growth rates.
        uptake_rates (dict): A dictionary where keys are model IDs and values are uptake rates 
                             of the selected substrate(s).

    Returns:
        dict: A dictionary where keys are model IDs and values are biomass yields.
    """

    # Initialize a dictionary where to store the biomass yields
    biomass_yields = {}
    
    # Calculate biomass yields on the selected substrate, for each model at each timpoint
    for model_id, _ in growth_rates.items():
        biomass_yields[model_id] = np.array(growth_rates[model_id]) / np.array(uptake_rates[model_id])
        # set yields to 0.0 when uptake is 0.0
        biomass_yields[model_id] = np.nan_to_num(biomass_yields[model_id], nan=0.0)
    
    return biomass_yields

#################################################################################
# calculate_biomass_fractions

def calculate_biomass_fractions(biomasses):
    """
    Calculates biomass fractions given a dictionary of biomasses.

    Args:
        biomasses (dict): A dictionary where keys are model IDs and values are lists of biomasses.

    Returns:
        dict: A dictionary where keys are model IDs and values are lists of biomass fractions.
    """
    biomass_fractions = {}
    
    # Calculate the sum of biomasses for each time point
    total_biomasses = [
        sum(biomasses[m][i] for m in biomasses)
        for i in range(len(next(iter(biomasses.values()))))
    ]
    
    # Iterate over each model in the input dictionary
    for model_id, model_biomasses in biomasses.items():
        biomass_fractions[model_id] = [
            model_biomasses[i] / total_biomasses[i] for i in range(len(model_biomasses))
        ]
    
    return biomass_fractions


