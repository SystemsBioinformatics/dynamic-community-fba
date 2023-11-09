import pytest
import cbmpy
from dcFBA.Models import CommunityModel
from dcFBA.ToyModels import model_a, model_b, model_c, combined_model
from dcFBA import DefaultModels


@pytest.fixture
def model_ecoli_core():
    return DefaultModels.read_default_model("e_coli_core")


@pytest.fixture
def model_strep_therm():
    return DefaultModels.read_default_model("strep_therm")


@pytest.fixture
def two_correctly_merged_ecoli():
    return cbmpy.loadModel("tests/data/two_ecoli_core.xml")


@pytest.fixture
def two_correctly_merged_ecoli_and_one_strep_therm():
    return cbmpy.loadModel(
        "tests/data//two_ecoli_core_and_one_strep_therm.xml"
    )


def test_if_ids_are_set(model_ecoli_core):
    model = CommunityModel(
        [model_ecoli_core, model_ecoli_core.clone()],
        ["R_BIOMASS_Ecoli_core_w_GAM", "R_BIOMASS_Ecoli_core_w_GAM"],
        ["e_coli_1", "e_coli_2"],
    )

    assert all(
        [
            model.m_identifiers == ["e_coli_1", "e_coli_2"],
            model.m_single_model_ids
            == [model_ecoli_core.getId(), model_ecoli_core.getId()],
        ]
    )


def test_if_a_single_model_can_be_added(
    model_ecoli_core, two_correctly_merged_ecoli
):
    model = CommunityModel(
        [model_ecoli_core.clone()],
        ["R_BIOMASS_Ecoli_core_w_GAM"],
        ["e_coli_1"],
    )

    model.add_model_to_community(
        model_ecoli_core.clone(), "R_BIOMASS_Ecoli_core_w_GAM", "e_coli_2"
    )

    compartments_are_copied_correctly = sorted(
        model.getCompartmentIds()
    ) == sorted(two_correctly_merged_ecoli.getCompartmentIds())

    reaction_are_copied_correctly = sorted(model.getReactionIds()) == sorted(
        two_correctly_merged_ecoli.getReactionIds()
    )

    model_identifiers_set_correctly = model.m_identifiers == [
        "e_coli_1",
        "e_coli_2",
    ]

    old_model_ids_set_correctly = model.m_single_model_ids == [
        model_ecoli_core.getId(),
        model_ecoli_core.getId(),
    ]

    del model

    assert all(
        [
            compartments_are_copied_correctly,
            reaction_are_copied_correctly,
            model_identifiers_set_correctly,
            old_model_ids_set_correctly,
        ]
    )


def test_community_model_of_three(
    model_ecoli_core,
    model_strep_therm,
    two_correctly_merged_ecoli_and_one_strep_therm,
):
    cm = CommunityModel(
        [
            model_ecoli_core,
            model_ecoli_core.clone(),
            model_strep_therm,
        ],
        [
            "R_BIOMASS_Ecoli_core_w_GAM",
            "R_BIOMASS_Ecoli_core_w_GAM",
            "R_biomass_STR",
        ],
        ["e_coli_1", "e_coli_2", "strepje"],
    )

    compartments_are_copied_correctly = sorted(
        cm.getCompartmentIds()
    ) == sorted(
        two_correctly_merged_ecoli_and_one_strep_therm.getCompartmentIds()
    )

    species_are_copied_correctly = sorted(cm.getSpeciesIds()) == sorted(
        two_correctly_merged_ecoli_and_one_strep_therm.getSpeciesIds()
    )

    reaction_are_copied_correctly = sorted(cm.getReactionIds()) == sorted(
        two_correctly_merged_ecoli_and_one_strep_therm.getReactionIds()
    )

    model_identifiers_set_correctly = cm.m_identifiers == [
        "e_coli_1",
        "e_coli_2",
        "strepje",
    ]

    old_model_ids_set_correctly = cm.m_single_model_ids == [
        model_ecoli_core.getId(),
        model_ecoli_core.getId(),
        model_strep_therm.getId(),
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


# Basically test joint FBA
def test_community_model_performs_fba(model_ecoli_core, model_strep_therm):
    community_model = CommunityModel(
        [model_ecoli_core, model_ecoli_core.clone(), model_strep_therm],
        [
            "R_BIOMASS_Ecoli_core_w_GAM",
            "R_BIOMASS_Ecoli_core_w_GAM",
            "R_biomass_STR",
        ],
        ["e_coli_1", "e_coli_2", "strepje"],
    )

    community_model.createSpecies(
        "X_c", False, "The community biomass", compartment="e"
    )

    for _, biomass_id in community_model.get_model_biomass_ids().items():
        reaction = community_model.getReaction(biomass_id)
        reaction.createReagent("X_c", 1)

    community_model.createReaction("X_comm")
    out = community_model.getReaction("X_comm")
    out.is_exchange = True
    out.setUpperBound(cbmpy.INF)
    out.setLowerBound(0)
    out.createReagent("X_c", -1)

    community_model.createObjectiveFunction("X_comm")
    community_model.setActiveObjective("X_comm_objective")

    sol = cbmpy.doFBA(community_model)
    assert round(sol, 10) == 25.6978994573


def test_if_single_model_is_deleted_correctly(model_ecoli_core):
    single_model = CommunityModel(
        [model_ecoli_core],
        ["R_BIOMASS_Ecoli_core_w_GAM"],
        ["e_coli_1"],
    )

    model = CommunityModel(
        [model_ecoli_core.clone(), model_ecoli_core.clone()],
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
