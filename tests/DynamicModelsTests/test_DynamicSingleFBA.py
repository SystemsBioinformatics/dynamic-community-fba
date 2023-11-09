import pytest
from dcFBA.DynamicModels import DynamicSingleFBA
from dcFBA import DefaultModels
from dcFBA.Models import KineticsStruct


@pytest.fixture
def model_ecoli_core():
    return DefaultModels.read_default_model("e_coli_core")


def test_simulation_is_correct(model_ecoli_core):
    model_ecoli_core.getReaction("R_GLCpts").setUpperBound(10)

    ds = DynamicSingleFBA(
        model_ecoli_core,
        "R_BIOMASS_Ecoli_core_w_GAM",
        0.1,
        {"M_glc__D_e": 10},
    )

    ds.simulate(0.15)

    T = ds.get_time_points()
    metabolites = ds.get_metabolites()
    biomass = ds.get_biomass()

    correct_number_of_time_points = len(T) == 20
    correct_final_time_point = round(T[-1], 2) == 2.85

    correct_final_biomass = round(biomass[-1], 2) == 0.97

    assert (
        correct_number_of_time_points
        and correct_final_time_point
        and correct_final_biomass
        and metabolites != {}
    )


def test_kinetics(model_ecoli_core):
    kin = KineticsStruct({"R_GLCpts": ("M_glc__D_e", 5, 10)})

    ds = DynamicSingleFBA(
        model_ecoli_core.clone(),
        "R_BIOMASS_Ecoli_core_w_GAM",
        0.1,
        {"M_glc__D_e": 10},
        kinetics=kin,
    )

    ds.simulate(0.15)

    T = ds.get_time_points()
    metabolites = ds.get_metabolites()
    biomass = ds.get_biomass()
    correct_number_of_time_points = len(T) == 40
    correct_final_time_point = round(T[-1], 2) == 5.85

    correct_final_biomass = round(biomass[-1], 2) == 0.88

    assert (
        correct_number_of_time_points
        and correct_final_time_point
        and correct_final_biomass
        and metabolites != {}
    )
