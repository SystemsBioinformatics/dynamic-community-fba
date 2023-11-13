import pytest
from dcFBA import DefaultModels
from dcFBA.DynamicModels import DynamicParallelFBA


@pytest.fixture(scope="module")
def model_ecoli_core():
    return DefaultModels.read_default_model("e_coli_core")


@pytest.fixture(scope="module")
def model_strep_therm():
    return DefaultModels.read_default_model("strep_therm")


@pytest.fixture(scope="module")
def dynamic_parallel_fba(model_ecoli_core, model_strep_therm):
    for rid in model_ecoli_core.getReactionIds():
        if rid.startswith("R_EX"):
            reaction = model_ecoli_core.getReaction(rid)
            reaction.is_exchange = True

    for rid in model_strep_therm.getReactionIds():
        if rid.startswith("R_EX"):
            reaction = model_strep_therm.getReaction(rid)
            reaction.is_exchange = True

    # Give the two models different glucose uptake rates
    model_ecoli_core.getReaction("R_GLCpts").setUpperBound(10)
    model_strep_therm.getReaction("R_GLCpts").setUpperBound(6)

    # Restrict Lactose uptake
    model_strep_therm.getReaction("R_LCTSGALex").setLowerBound(0)
    model_strep_therm.getReaction("R_LCTSt6").setUpperBound(30)

    model_ecoli_core.setId("ecoli")
    model_strep_therm.setId("strep")

    parallel_fba = DynamicParallelFBA(
        [model_ecoli_core, model_strep_therm],
        [1.0, 1.0],
        {
            "M_glc__D_e": 100,
            "M_succ_e": 0,
            "M_glu__L_e": 0.0,
            "M_gln__L_e": 0.0,
            "M_lcts_e": 100,
        },
    )

    parallel_fba.simulate(0.1)

    return parallel_fba


def test_final_biomass_concentrations(dynamic_parallel_fba):
    biomasses = dynamic_parallel_fba.get_biomasses()
    assert (
        round(biomasses["ecoli"][-1], 3) == 2.998
        and round(biomasses["strep"][-1], 3) == 5.986
    )


def test_simulation_length(dynamic_parallel_fba):
    time_points = dynamic_parallel_fba.get_time_points()

    assert len(time_points) == 14


# # Not testing intermediate fluxes or metabolite concentrations as I think this
# # might vary if you would switch the linear solver. Just checking they are
# # Not empty
def test_dynamic_parallel_fba(dynamic_parallel_fba):
    metabolites = dynamic_parallel_fba.get_metabolites()
    fluxes = dynamic_parallel_fba.get_fluxes()

    assert metabolites != {} and len(fluxes) != 0


def test_community_growth_rate(dynamic_parallel_fba):
    comm = dynamic_parallel_fba.get_community_growth_rate()

    assert round(comm[-1], 4) == 0.7535


def test_if_relative_abundance_is_calculated_correct(dynamic_parallel_fba):
    abundances = dynamic_parallel_fba.get_relative_abundance()

    total = abundances["ecoli"][-1] + abundances["strep"][-1]

    assert (
        round(abundances["ecoli"][-1], 4) == 0.3337
        and round(abundances["strep"][-1], 4) == 0.6663
        and round(total, 7) == 1.0
    )
