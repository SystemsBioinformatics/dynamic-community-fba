import pytest
from dcFBA import DefaultModels
from dcFBA.DynamicModels import DynamicJointFBA
from dcFBA.Models import CommunityModel


@pytest.fixture(scope="module")
def model_ecoli_core():
    return DefaultModels.read_default_model("e_coli_core")


@pytest.fixture(scope="module")
def model_strep_therm():
    return DefaultModels.read_default_model("strep_therm")


@pytest.fixture(scope="module")
def dynamic_joint_fba(model_ecoli_core, model_strep_therm):
    model_ecoli_core.getReaction("R_GLCpts").setUpperBound(10)
    model_strep_therm.getReaction("R_GLCpts").setUpperBound(6)

    cm = CommunityModel(
        [model_ecoli_core, model_strep_therm],
        ["R_BIOMASS_Ecoli_core_w_GAM", "R_biomass_STR"],
        ["ecoli", "strep"],
    )

    dynamic_joint_fba = DynamicJointFBA(
        cm,
        [1.0, 1.0],
        {
            "M_glc__D_e": 100,
            "M_succ_e": 0,
            "M_glu__L_e": 0.0,
            "M_gln__L_e": 0.0,
            "M_lcts_e": 100,
        },
    )

    dynamic_joint_fba.simulate(0.1)

    return dynamic_joint_fba


def test_final_biomass_concentrations(dynamic_joint_fba):
    biomasses = dynamic_joint_fba.get_biomasses()
    assert (
        round(biomasses["ecoli"][-1], 3) == 34.806
        and round(biomasses["strep"][-1], 3) == 37.642
    )


def test_simulation_length(dynamic_joint_fba):
    time_points = dynamic_joint_fba.get_time_points()

    assert len(time_points) == 24


# Not testing intermediate fluxes or metabolite concentrations as I think this
# might vary if you would switch the linear solver. Just checking they are
# Not empty
def test_dynamic_joint_fba(dynamic_joint_fba):
    metabolites = dynamic_joint_fba.get_metabolites()
    fluxes = dynamic_joint_fba.get_fluxes()

    assert metabolites != {} and len(fluxes) != 0


def test_community_growth_rate(dynamic_joint_fba):
    comm = dynamic_joint_fba.get_community_growth_rate()

    assert round(comm[-1], 4) == 0.2914


def test_if_relative_abundance_is_calculated_correct(dynamic_joint_fba):
    abundances = dynamic_joint_fba.get_relative_abundance()

    total = abundances["ecoli"][-1] + abundances["strep"][-1]
    assert (
        round(abundances["ecoli"][-1], 7) == 0.4804275
        and round(abundances["strep"][-1], 7) == 0.5195725
        and round(total, 7) == 1.0
    )
