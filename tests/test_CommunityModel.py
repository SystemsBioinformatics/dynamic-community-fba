import pytest
import cbmpy
from DCFBA.Models import CommunityModel
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


@pytest.fixture
def e_coli_core():
    return cbmpy.loadModel("data/bigg_models/e_coli_core.xml")


@pytest.fixture
def correctly_combined_e_coli_core():
    return cbmpy.loadModel("data/combined_e_coli_core.xml")


def test_if_combined_model_is_created(
    model_A: cbmpy.CBModel.Model,
    model_B: cbmpy.CBModel.Model,
    model_C: cbmpy.CBModel.Model,
    defined_combined_model: cbmpy.CBModel.Model,
):
    model = CommunityModel(
        [model_A, model_B, model_C], ["R_BM_A", "R_BM_B", "R_BM_C"]
    )

    compartments_are_copied_correctly = sorted(
        model.getCompartmentIds()
    ) == sorted(defined_combined_model.getCompartmentIds())

    species_are_copied_correctly = sorted(model.getSpeciesIds()) == sorted(
        defined_combined_model.getSpeciesIds()
    )

    reaction_are_copied_correctly = sorted(model.getReactionIds()) == sorted(
        defined_combined_model.getReactionIds()
    )

    model_identifiers_set_correctly = model.m_identifiers == [
        model_A.getId(),
        model_B.getId(),
        model_C.getId(),
    ]

    old_model_ids_set_correctly = model.m_single_model_ids == [
        model_A.getId(),
        model_B.getId(),
        model_C.getId(),
    ]

    assert all(
        [
            compartments_are_copied_correctly,
            species_are_copied_correctly,
            reaction_are_copied_correctly,
            model_identifiers_set_correctly,
            old_model_ids_set_correctly,
        ]
    )


def test_if_ids_are_set(e_coli_core):
    model = CommunityModel(
        [e_coli_core, e_coli_core.clone()],
        ["R_BIOMASS_Ecoli_core_w_GAM", "R_BIOMASS_Ecoli_core_w_GAM"],
        ["e_coli_1", "e_coli_2"],
    )

    assert all(
        [
            model.m_identifiers == ["e_coli_1", "e_coli_2"],
            model.m_single_model_ids
            == [e_coli_core.getId(), e_coli_core.getId()],
        ]
    )


def test_if_a_single_model_can_be_added(
    e_coli_core, correctly_combined_e_coli_core
):
    model = CommunityModel(
        [e_coli_core], ["R_BIOMASS_Ecoli_core_w_GAM"], ["e_coli_1"]
    )

    model.add_model_to_community(
        e_coli_core, "R_BIOMASS_Ecoli_core_w_GAM", "e_coli_2"
    )

    compartments_are_copied_correctly = sorted(
        model.getCompartmentIds()
    ) == sorted(correctly_combined_e_coli_core.getCompartmentIds())

    reaction_are_copied_correctly = sorted(model.getReactionIds()) == sorted(
        correctly_combined_e_coli_core.getReactionIds()
    )

    model_identifiers_set_correctly = model.m_identifiers == [
        "e_coli_1",
        "e_coli_2",
    ]

    old_model_ids_set_correctly = model.m_single_model_ids == [
        e_coli_core.getId(),
        e_coli_core.getId(),
    ]

    assert all(
        [
            compartments_are_copied_correctly,
            reaction_are_copied_correctly,
            model_identifiers_set_correctly,
            old_model_ids_set_correctly,
        ]
    )


def test_if_single_model_is_deleted_correctly(e_coli_core):
    single_model = CommunityModel(
        [e_coli_core], ["R_BIOMASS_Ecoli_core_w_GAM"], ["e_coli_1"]
    )

    model = CommunityModel(
        [e_coli_core, e_coli_core.clone()],
        ["R_BIOMASS_Ecoli_core_w_GAM", "R_BIOMASS_Ecoli_core_w_GAM"],
        ["e_coli_1", "e_coli_2"],
    )

    model.remove_model_from_community("e_coli_2")

    model_identifiers_set_correctly = (
        single_model.m_identifiers == model.m_identifiers
    )

    old_model_ids_set_correctly = (
        model.m_single_model_ids == single_model.m_single_model_ids
    )

    biomasses_are_set_correctly = (
        model.m_single_model_biomass_reaction_ids
        == single_model.m_single_model_biomass_reaction_ids
    )

    reaction_are_deleted_correctly = sorted(model.getReactionIds()) == sorted(
        single_model.getReactionIds()
    )

    species_are_deleted_correctly = len(model.getSpeciesIds()) == len(
        single_model.getSpeciesIds()
    )

    assert all(
        [
            model_identifiers_set_correctly,
            old_model_ids_set_correctly,
            biomasses_are_set_correctly,
            reaction_are_deleted_correctly,
            species_are_deleted_correctly,
        ]
    )
