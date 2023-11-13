import pytest
from dcFBA.Helpers.OptimalSearch import time_search, balance_search_clean
from dcFBA.ToyModels import model_a, model_b
from dcFBA.Models import CommunityModel


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


def test_optimal_time_search(model_A, model_B):
    community_model = CommunityModel(
        [model_A, model_B], ["R_BM_A", "R_BM_B"], ["modelA", "modelB"]
    )

    n, value = time_search(
        community_model,
        {"modelA": 1.0, "modelB": 2.0},
        {"S_e": 100, "A_e": 0.0, "B_e": 0.0},
        dt=0.1,
    )

    assert n == 21 and round(value, 3) == 12.778


def test_user_defined_objective(model_A, model_B):
    community_model = CommunityModel(
        [model_A, model_B], ["R_BM_A", "R_BM_B"], ["modelA", "modelB"]
    )

    n, value = time_search(
        community_model,
        {"modelA": 1.0, "modelB": 2.0},
        {"S_e": 100, "A_e": 0.0, "B_e": 0.0},
        dt=0.1,
        set_values=[7.5, 21],
    )

    assert n == 14 and round(value, 3) == 7.891


def test_balance_search(model_A, model_B):
    community_model = CommunityModel(
        [model_A, model_B], ["R_BM_A", "R_BM_B"], ["modelA", "modelB"]
    )
    value = balance_search_clean(
        community_model,
        25,
        {"S_e": 100, "A_e": 0.0, "B_e": 0.0},
        0.1,
        3,
        12.77778,
    )

    assert round((value * 12.77778), 4) == 12.678
