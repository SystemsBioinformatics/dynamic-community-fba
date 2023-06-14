import cbmpy


def build_model_B():
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

    model.createCompartment("c", "cytosol")
    model.createCompartment("e", "external")

    model.createReaction(
        "R_1", "Import S external to the cytosol", reversible=True
    )
    model.createReactionReagent("R_1", "S_e", -1)
    model.createReactionReagent("R_1", "S_c", 1)
    model.setReactionBounds("R_1", -1000.0, 1000.0)

    model.createReaction("R_13", "Create X from S_c", reversible=True)
    model.createReactionReagent("R_2", "S_c", -4)
    model.createReactionReagent("R_2", "X_c", 3)
    model.setReactionBounds("R_2", -1000.0, 1000.0)

    model.createReaction("R_14", "Create BM_c_B from X", reversible=True)
    model.createReactionReagent("R_9", "X_c", -7)
    model.createReactionReagent("R_9", "A_c", -3)
    model.createReactionReagent("R_9", "E_c", 1)
    model.createReactionReagent("R_9", "BM_c_B", 1)
    model.setReactionBounds("R_9", -1000.0, 1000.0)

    # Biomass creation
    model.createReaction(
        "R_BM_B",
        "Biomass reaction of species B",
        reversible=True,
    )
    model.createReactionReagent("R_BM_B", "BM_c_B", -1)
    model.createReactionReagent("R_BM_B", "BM_e_B", 1)
    model.setReactionBounds("R_BM_B", -1000.0, 1000.0)

    # EXchange reactions:
    model.createReaction("S_exchange", reversible=True)
    model.createReactionReagent("S_exchange", "S_e", 1)
    model.setReactionBounds("S_exchange", lower=-1000.0, upper=1000.0)
    model.getReaction("S_exchange").is_exchange = True

    model.createReaction(
        "B_exchange", reversible=True
    )  # Q2.0 - should 'reversible' be True or False?
    model.createReactionReagent("B_exchange", "B_e", 1)
    model.setReactionBounds("B_exchange", lower=-1000.0, upper=1000.0)
    model.getReaction("B_exchange").is_exchange = True

    model.createReaction("A_exchange", reversible=True)
    model.createReactionReagent("A_exchange", "A_e", -1)
    model.setReactionBounds("A_exchange", -1000.0, 1000)
    model.getReaction("A_exchange").is_exchange = True

    # E can't be imported
    model.createReaction("E_exchange", reversible=False)
    model.createReactionReagent("E_exchange", "E_e", -1)
    model.setReactionBounds("E_exchange", 0, 1000)
    model.getReaction("E_exchange").is_exchange = True

    # BM sink
    model.createReaction("BM_e_B_exchange", reversible=False)
    model.createReactionReagent("BM_e_B_exchange", "BM_e_B", -1)
    model.setReactionBounds("BM_e_B_exchange", 0, 1000)
    model.getReaction("BM_e_B_exchange").is_exchange = True

    cbmpy.CBModel.FluxObjective(
        "objective",
        "BM_e_B_exchange",
        1,
    )

    model.createObjectiveFunction("test_objective")

    model.setObjectiveFlux("BM_e_B_exchange")

    return model
