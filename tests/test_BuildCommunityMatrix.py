import pytest
import cbmpy

from DCFBA.Helpers import BuildCommunityMatrix as bcm
from DCFBA.ToyModels import model_a, model_b, model_c, combined_model


@pytest.fixture
def model_A():
    return model_a.build_model_A()


@pytest.fixture
def model_B():
    return model_b.build_model_B()


@pytest.fixture
def model_C():
    return model_c.build_model_C()


@pytest.fixture
def defined_combined_model():
    return combined_model.build_combined_model()


def test_if_n_models_are_merged_correctly(
    model_A: cbmpy.CBModel.Model,
    model_B: cbmpy.CBModel.Model,
    model_C: cbmpy.CBModel.Model,
    defined_combined_model: cbmpy.CBModel.Model,
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

    combined_model = bcm.combine_models([model_A, model_B, model_C])

    compartments_are_copied_correctly = sorted(
        combined_model.getCompartmentIds()
    ) == sorted(defined_combined_model.getCompartmentIds())

    species_are_copied_correctly = sorted(
        combined_model.getSpeciesIds()
    ) == sorted(defined_combined_model.getSpeciesIds())

    reaction_are_copied_correctly = sorted(
        combined_model.getReactionIds()
    ) == sorted(defined_combined_model.getReactionIds())

    print(combined_model.getReactionIds())
    print(defined_combined_model.getReactionIds())
    assert all(
        [
            compartments_are_copied_correctly,
            species_are_copied_correctly,
            reaction_are_copied_correctly,
        ]
    )


def test_if_community_matrix_performance_FBA(model_A, model_B, model_C):
    combined_model = bcm.combine_models([model_A, model_B, model_C])

    # Create objective function
    combined_model.createReaction("BM_A_B_C", "Biomass reaction both", False)
    combined_model.createReactionReagent("BM_A_B_C", "BM_e_A", -1)
    combined_model.createReactionReagent("BM_A_B_C", "BM_e_B", -3)
    combined_model.createReactionReagent("BM_A_B_C", "BM_e_C", -1)

    combined_model.setReactionLowerBound("BM_A_B_C", 0)
    combined_model.getReaction("BM_A_B_C").is_exchange = True

    combined_model.createObjectiveFunction("BM_A_B_C")

    assert round(cbmpy.doFBA(combined_model), 1) == 55.7


def test_if_dict_contains_all_duplicate_species(model_A, model_B, model_C):
    dic: dict[str, int] = bcm.create_duplicate_species_dict(
        [model_A, model_B, model_C]
    )

    assert dic == {
        "S_e": 3,
        "S_c": 3,
        "B_e": 2,
        "Z_c": 2,
        "X_c": 2,
        "A_e": 2,
        "B_c": 2,
    }


def test_if_species_with_no_reactions_is_not_added(model_A, model_C):
    combined_model = bcm.combine_models([model_A, model_C])

    dummy_species = model_C.getSpecies("Dummy_species")

    assert combined_model.getSpecies(dummy_species.id) is None


# TODO test if Exchange reactions are being correctly handeld

# TODO test if all ids are set correctly if there are ids provided by the user
