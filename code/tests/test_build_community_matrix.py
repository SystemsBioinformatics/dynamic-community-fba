from .helpers import model_a, model_b, combined_model
from cbmpy.CBModel import Model, Reaction
import pytest
import cbmpy
import endPointFBA.build_community_matrix as bcm
import pandas as pd


@pytest.fixture
def model_A():
    return model_a.build_model_A()


@pytest.fixture
def model_B():
    return model_b.build_model_B()


@pytest.fixture
def defined_combined_model():
    return combined_model.build_combined_model()


def test_if_n_models_are_merged_correctly(
    model_A: Model, model_B: Model, defined_combined_model: Model
):
    """For the community matrix to be correct the compartments, species and
    reactions have to be copied accordingly and follow pur template rules for
    a right community matrix

    Args:
        model_A (Model): Dummy model A
        model_B (Model): Dummy model B
        defined_combined_model (Model): Preprogrammed correct community model

    Asserts true if all three demands are correct
    """
    combined_model = bcm.combine_models([model_A, model_B])

    compartments_are_copied_correctly = sorted(
        combined_model.getCompartmentIds()
    ) == sorted(defined_combined_model.getCompartmentIds())

    species_are_copied_correctly = sorted(
        combined_model.getSpeciesIds()
    ) == sorted(defined_combined_model.getSpeciesIds())

    reaction_are_copied_correctly = sorted(
        combined_model.getReactionIds()
    ) == sorted(defined_combined_model.getReactionIds())

    assert all(
        [
            compartments_are_copied_correctly,
            species_are_copied_correctly,
            reaction_are_copied_correctly,
        ]
    )


def test_if_community_matrix_performance_FBA(model_A, model_B):
    combined_model = bcm.combine_models([model_A, model_B])

    # Create objective function
    combined_model.createReaction("BM_A_and_B", "Biomass reaction both", False)
    combined_model.createReactionReagent("BM_A_and_B", "BM_e_A", -1)
    combined_model.createReactionReagent("BM_A_and_B", "BM_e_B", -3)
    combined_model.setReactionLowerBound("BM_A_and_B", 0)
    combined_model.getReaction("BM_A_and_B").is_exchange = True

    combined_model.createObjectiveFunction("BM_A_and_B")

    assert round(cbmpy.doFBA(combined_model), 1) == 57.1


def test_if_a_none_type_is_skipped():
    pass


def test_if_dict_contains_all_duplicate_species():
    pass


def test_merge_compartments():
    pass


def test_merge_reactions():
    pass


def fix_duplicate_species():
    pass


def test_handle_duplicate_species():
    pass
