import pytest
from dcFBA.DynamicModels import EndPointFBA
from dcFBA.Models import CommunityModel
from dcFBA.ToyModels import model_a, model_b


@pytest.fixture(scope="module")
def model_A():
    m_a = model_a.build_toy_model_fba_A()
    m_a.getReaction("R_1").setUpperBound(10)
    m_a.getReaction("R_4").setUpperBound(3)
    m_a.getReaction("R_6").setUpperBound(1)

    m_a.getReaction("R_1").setLowerBound(0)
    m_a.getReaction("R_4").setLowerBound(0)
    m_a.getReaction("R_6").setLowerBound(0)

    m_a.getReaction("R_BM_A").deleteReagentWithSpeciesRef("BM_e_A")

    return m_a


@pytest.fixture(scope="module")
def model_B():
    m_b = model_b.build_toy_model_fba_B()
    m_b.getReaction("R_1").setUpperBound(10)
    m_b.getReaction("R_3").setUpperBound(1)
    m_b.getReaction("R_5").setUpperBound(1)

    m_b.getReaction("R_1").setLowerBound(0)
    m_b.getReaction("R_3").setLowerBound(0)
    m_b.getReaction("R_5").setLowerBound(0)

    # Delete the original biomass model from model B
    m_b.getReaction("R_BM_B").deleteReagentWithSpeciesRef("BM_e_B")

    return m_b


@pytest.fixture(scope="module")
def EndPointFBA_simulation(model_A, model_B):
    community_model = CommunityModel(
        [model_A, model_B], ["R_BM_A", "R_BM_B"], ["modelA", "modelB"]
    )

    ep = EndPointFBA(
        community_model,
        25,
        {"modelA": 1.0, "modelB": 2.0},
        {"S_e": 100, "A_e": 0.0, "B_e": 0.0},
        dt=0.1,
    )

    ep.simulate()

    return ep


def test_final_biomass_concentrations(EndPointFBA_simulation):
    biomasses = EndPointFBA_simulation.get_biomasses()
    assert (
        round(biomasses["modelA"][-1], 3) == 2.667
        and round(biomasses["modelB"][-1], 3) == 13.111
    )


# Not testing intermediate fluxes or metabolite concentrations as I think this
# might vary if you would switch the linear solver. Just checking they are
# Not empty
def test_EndPointFBA_simulation(EndPointFBA_simulation):
    metabolites = EndPointFBA_simulation.get_metabolites()
    fluxes = EndPointFBA_simulation.get_fluxes()

    assert metabolites != {} and len(fluxes) != 0


def test_community_growth_rate(EndPointFBA_simulation):
    comm = EndPointFBA_simulation.get_community_growth_rate()

    assert round(comm[-1], 4) == 0.9122


def test_if_relative_abundance_is_calculated_correct(EndPointFBA_simulation):
    abundances = EndPointFBA_simulation.get_relative_abundance()

    total = abundances["modelA"][-1] + abundances["modelB"][-1]
    assert (
        round(abundances["modelA"][-1], 7) == 0.1690141
        and round(abundances["modelB"][-1], 7) == 0.8309859
        and round(total, 7) == 1.0
    )


def test_balanced_growth(model_A, model_B):
    community_model = CommunityModel(
        [model_A, model_B], ["R_BM_A", "R_BM_B"], ["modelA", "modelB"]
    )

    ep = EndPointFBA(
        community_model,
        25,
        {"modelA": 1.0, "modelB": 2.0},
        {"S_e": 100, "A_e": 0.0, "B_e": 0.0},
        dt=0.1,
    )

    ep.balanced_growth(3, 12.76)
    ep.simulate()

    abundances = ep.get_relative_abundance()

    total_t0 = abundances["modelA"][0] + abundances["modelB"][0]
    total_tN = abundances["modelA"][-1] + abundances["modelB"][-1]

    phi_modelA_t0 = abundances["modelA"][0] / total_t0
    phi_modelA_tn = abundances["modelA"][-1] / total_t0

    phi_modelB_t0 = abundances["modelB"][0] / total_t0
    phi_modelB_tn = abundances["modelB"][-1] / total_t0

    assert (
        round(total_t0, 7) == 1.0
        and round(total_tN, 7) == 1.0
        and round(phi_modelA_t0, 7) == round(phi_modelA_tn, 7)
        and round(phi_modelB_t0, 7) == round(phi_modelB_tn, 7)
    )


def test_qp(model_A, model_B):
    community_model = CommunityModel(
        [model_A, model_B], ["R_BM_A", "R_BM_B"], ["modelA", "modelB"]
    )

    ep = EndPointFBA(
        community_model,
        25,
        {"modelA": 1.0, "modelB": 2.0},
        {"S_e": 100, "A_e": 0.0, "B_e": 0.0},
        dt=0.1,
    )

    ep.set_qp(12.7777)
    solution = ep.simulate()

    assert round(solution, 3) == 0.426
