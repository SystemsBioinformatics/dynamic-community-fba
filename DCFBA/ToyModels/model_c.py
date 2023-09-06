import cbmpy


def build_model_C():
    model = cbmpy.CBModel.Model("Organism_C")

    model.createSpecies(
        sid="S_e",
        chemFormula="C11H21N2O7PRS",
        name="S external",
        compartment="e",
    )

    model.createSpecies(sid="S_c", name="S cytosol", compartment="c")

    model.createSpecies("F_c", name="F cytosol", compartment="c")

    model.createSpecies(
        "BM_c_C", name="Biomass cytosol species B", compartment="c"
    )

    model.createSpecies(
        "BM_e_C", name="Biomass external species B", compartment="e"
    )

    # Dummy species that won't occur in a reaction to test
    # behavior of the method

    model.createSpecies(
        "Dummy_species",
        name="This species does not occur in any reaction",
        compartment="c",
    )

    model.createCompartment("c", "cytosol")
    model.createCompartment("e", "external")

    model.createReaction(
        "R_1", "Import S external to the cytosol", reversible=True
    )
    model.createReactionReagent("R_1", "S_e", -1)
    model.createReactionReagent("R_1", "S_c", 1)
    model.setReactionBounds("R_1", -1000.0, 1000.0)

    model.createReaction("R_13", "Create F from S_c", reversible=True)
    model.createReactionReagent("R_13", "S_c", -1)
    model.createReactionReagent("R_13", "F_c", 1)
    model.setReactionBounds("R_13", -1000.0, 1000.0)

    model.createReaction("R_14", "Create BM_c_C from F", reversible=True)
    model.createReactionReagent("R_14", "F_c", -1)
    model.createReactionReagent("R_14", "BM_c_C", 1)
    model.setReactionBounds("R_14", -1000.0, 1000.0)

    # Biomass creation
    model.createReaction(
        "R_BM_C",
        "Biomass reaction of species C",
        reversible=True,
    )
    model.createReactionReagent("R_BM_C", "BM_c_C", -1)
    model.createReactionReagent("R_BM_C", "BM_e_C", 1)
    model.setReactionBounds("R_BM_C", 0, 1000.0)

    # EXchange reactions:
    model.createReaction("S_exchange", reversible=True)
    model.createReactionReagent("S_exchange", "S_e", -1)
    model.setReactionBounds("S_exchange", lower=-1000.0, upper=1000.0)
    model.getReaction("S_exchange").is_exchange = True

    # BM sink
    model.createReaction("BM_e_C_exchange", reversible=False)
    model.createReactionReagent("BM_e_C_exchange", "BM_e_C", -1)
    model.setReactionBounds("BM_e_C_exchange", 0, 1000)
    model.getReaction("BM_e_C_exchange").is_exchange = True

    cbmpy.CBModel.FluxObjective(
        "objective",
        "BM_e_C_exchange",
        1,
    )

    model.createObjectiveFunction("test_objective")

    model.setObjectiveFlux("BM_e_C_exchange")

    return model
