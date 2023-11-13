import pytest
import cbmpy

from dcFBA.Helpers import BuildCommunityMatrix as bcm
from dcFBA import DefaultModels
from dcFBA.ToyModels import model_a, model_b, model_c


@pytest.fixture
def model_ecoli_core():
    return DefaultModels.read_default_model("e_coli_core")


@pytest.fixture
def model_strep_therm():
    return DefaultModels.read_default_model("strep_therm")


@pytest.fixture
def model_A():
    return model_a.build_model_A()


@pytest.fixture
def model_B():
    return model_b.build_model_B()


@pytest.fixture
def model_C():
    return model_c.build_model_C()


def test_check_ids_raises_error_on_too_few_ids(
    model_ecoli_core, model_strep_therm
):
    with pytest.raises(Exception) as e_info:
        bcm.check_ids(["test1"], [model_ecoli_core, model_strep_therm])

    assert str(e_info.value) == "Too few ids were provided"


def test_check_ids_raises_error_on_same_ids(
    model_ecoli_core, model_strep_therm
):
    with pytest.raises(Exception) as e_info:
        bcm.check_ids(
            ["test1", "test1"], [model_ecoli_core, model_strep_therm]
        )

    assert str(e_info.value) == "Model identifiers should be unique!"


def test_check_ids_raises_error_on_id_that_exists_as_reaction_id_in_the_model(
    model_ecoli_core, model_strep_therm
):
    with pytest.raises(Exception) as e_info:
        bcm.check_ids(["GLCp", "test1"], [model_ecoli_core, model_strep_therm])

    assert (
        str(e_info.value)
        == "The provided ids are not unique in the model.\
                                Please provide an alternative model \
                                identifier."
    )


def test_check_ids_finishes(model_ecoli_core, model_strep_therm):
    nn = bcm.check_ids(
        ["modelA", "modelB"], [model_ecoli_core, model_strep_therm]
    )

    assert nn is None


def test_if_a_reaction_is_copied_to_the_new_model(
    model_ecoli_core, model_strep_therm
):
    bcm.copy_reaction(model_strep_therm, model_ecoli_core, "R_ASNN")
    old_reaction: cbmpy.CBModel.Reaction = model_strep_therm.getReaction(
        "R_ASNN"
    )
    old_reagents = old_reaction.getReagentObjIds()
    reaction_is_copied: bool = "R_ASNN" in model_ecoli_core.getReactionIds()

    new_reaction = model_ecoli_core.getReaction("R_ASNN")
    new_reagents = new_reaction.getReagentObjIds()

    all_reagents_are_copied: bool = all(
        [reag in new_reagents for reag in old_reagents]
    )

    old_species = old_reaction.getSpeciesIds()
    new_species = new_reaction.getSpeciesIds()

    all_species_are_coped: bool = all(
        [spe in new_species for spe in old_species]
    )

    assert (
        reaction_is_copied
        and all_reagents_are_copied
        and all_species_are_coped
    )


# Good example of a function (copy_species_and_reagents) that does not adhere to
# the single responsibility principle, hence this test is terrible to write
# TODO split the function in copy_species and copy_reagents
def test_copy_species_and_reagents(
    model_ecoli_core: cbmpy.CBModel.Model,
) -> None:
    sid = "M_g6p_c"
    species: cbmpy.CBModel.Species = model_ecoli_core.getSpecies(sid)
    reactions: list[str] = species.isReagentOf()

    old_reagents = []

    for rid in reactions:
        reaction: cbmpy.CBModel.Reaction = model_ecoli_core.getReaction(rid)
        old_reagents.append(reaction.getReagentWithSpeciesRef(sid))

    bcm.copy_species_and_reagents(model_ecoli_core, species, "ecoli16")

    if model_ecoli_core.getSpecies(sid) is not None:
        assert False

    for rid in reactions:
        reaction = model_ecoli_core.getReaction(rid)
        if (
            sid in reaction.getSpeciesIds()
            or f"{sid}_ecoli16" not in reaction.getSpeciesIds()
        ):
            assert False
        reagent: cbmpy.CBModel.Reagent = reaction.getReagentWithSpeciesRef(
            f"{sid}_ecoli16"
        )
        if reagent in old_reagents or reagent is None:
            assert False
    # Check if the old species is deleted
    if model_ecoli_core.getSpecies(sid) is not None:
        assert False


def test_merge_compartments(model_ecoli_core, model_strep_therm):
    bcm.merge_compartments(model_ecoli_core, model_strep_therm, "ecoli16")
    compartment_ids = model_strep_therm.getCompartmentIds()

    assert compartment_ids == ["e", "c", "c_ecoli16"]


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


# def test_if_three_models_are_merged_correctly(
#     model_ecoli_core: cbmpy.CBModel.Model,
#     model_strep_therm: cbmpy.CBModel.Model,
# ):
#     """For the community matrix to be correct the compartments, species and
#     reactions have to be copied accordingly and follow pur template rules for
#     a right community matrix

#     Args:
#         model_A (Model): Dummy model A
#         model_B (Model): Dummy model B
#         defined_combined_model (Model): Preprogrammed correct community model

#     Asserts true if all three demands are correct
#     """

#     combined_model = bcm.combine_models([model_ecoli_core, model_strep_therm])

#     compartments_are_copied_correctly = sorted(
#         combined_model.getCompartmentIds()
#     ) == sorted(defined_combined_model.getCompartmentIds())

#     species_are_copied_correctly = sorted(
#         combined_model.getSpeciesIds()
#     ) == sorted(defined_combined_model.getSpeciesIds())

#     reaction_are_copied_correctly = sorted(
#         combined_model.getReactionIds()
#     ) == sorted(defined_combined_model.getReactionIds())

#     print(combined_model.getReactionIds())
#     print(defined_combined_model.getReactionIds())
#     assert all(
#         [
#             compartments_are_copied_correctly,
#             species_are_copied_correctly,
#             reaction_are_copied_correctly,
#         ]
#     )


# def test_if_community_matrix_performance_FBA(model_A, model_B, model_C):
#     combined_model = bcm.combine_models([model_A, model_B, model_C])

#     # Create objective function
#     combined_model.createReaction("BM_A_B_C", "Biomass reaction both", False)
#     combined_model.createReactionReagent("BM_A_B_C", "BM_e_A", -1)
#     combined_model.createReactionReagent("BM_A_B_C", "BM_e_B", -3)
#     combined_model.createReactionReagent("BM_A_B_C", "BM_e_C", -1)

#     combined_model.setReactionLowerBound("BM_A_B_C", 0)
#     combined_model.getReaction("BM_A_B_C").is_exchange = True

#     combined_model.createObjectiveFunction("BM_A_B_C")

#     assert round(cbmpy.doFBA(combined_model), 1) == 55.7
