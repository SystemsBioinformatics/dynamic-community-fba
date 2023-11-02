import cbmpy


def build_model_B():
    model = cbmpy.CBModel.Model("Organism_B")

    model.createSpecies(
        sid="S_e",
        chemFormula="C11H21N2O7PRS",
        name="S external",
        compartment="e",
    )

    model.createSpecies(sid="S_c", name="S cytosol", compartment="c")

    model.createSpecies("Z_c", name="Z cytosol", compartment="c")

    model.createSpecies("E_c", name="D cytosol", compartment="c")

    model.createSpecies("E_e", name="D external", compartment="e")

    model.createSpecies("X_c", name="X cytosol", compartment="c")

    model.createSpecies("A_c", name="A cytosol", compartment="c")

    model.createSpecies("A_e", name="A external", compartment="e")

    model.createSpecies("B_c", name="B cytosol", compartment="c")

    model.createSpecies("B_e", name="B external", compartment="e")

    model.createSpecies(
        "BM_c_B", name="Biomass cytosol species B", compartment="c"
    )

    model.createSpecies(
        "BM_e_B", name="Biomass external species B", compartment="e"
    )

    model.createCompartment("c", "cytosol")
    model.createCompartment("e", "external")

    model.createReaction(
        "R_1", "Import S external to the cytosol", reversible=True
    )
    model.createReactionReagent("R_1", "S_e", -1)
    model.createReactionReagent("R_1", "S_c", 1)
    model.setReactionBounds("R_1", -1000.0, 1000.0)

    model.createReaction("R_2", "Create X from S_c", reversible=True)
    model.createReactionReagent("R_2", "S_c", -4)
    model.createReactionReagent("R_2", "X_c", 3)
    model.setReactionBounds("R_2", -1000.0, 1000.0)

    model.createReaction("R_9", "Create BM_c_B from X", reversible=True)
    model.createReactionReagent("R_9", "X_c", -7)
    model.createReactionReagent("R_9", "A_c", -3)
    model.createReactionReagent("R_9", "E_c", 1)
    model.createReactionReagent("R_9", "BM_c_B", 1)
    model.setReactionBounds("R_9", -1000.0, 1000.0)

    model.createReaction("R_10", "Import A_c from A_e", reversible=True)
    model.createReactionReagent("R_10", "A_e", -1)
    model.createReactionReagent("R_10", "A_c", 1)
    model.setReactionBounds("R_10", -1000.0, 1000.0)

    model.createReaction("R_5", "Create Z from S_c", reversible=True)
    model.createReactionReagent("R_5", "S_c", -2)
    model.createReactionReagent("R_5", "Z_c", 1)
    model.setReactionBounds("R_5", -1000.0, 1000.0)

    model.createReaction("R_6", "Create BM_c_B and B from Z", reversible=True)
    model.createReactionReagent("R_6", "Z_c", -9)
    model.createReactionReagent("R_6", "B_c", 1)
    model.createReactionReagent("R_6", "BM_c_B", 1)
    model.setReactionBounds("R_6", -1000.0, 1000.0)

    # Import B
    model.createReaction("R_11", "Export B_c to B_e", reversible=True)
    model.createReactionReagent("R_11", "B_c", -1)
    model.createReactionReagent("R_11", "B_e", 1)
    model.setReactionBounds("R_11", -1000.0, 1000.0)

    model.createReaction("R_12", "Export E_c to E_e", reversible=True)
    model.createReactionReagent("R_12", "E_c", -1)
    model.createReactionReagent("R_12", "E_e", 1)
    model.setReactionBounds("R_12", -1000.0, 1000.0)

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
    model.createReactionReagent("S_exchange", "S_e", -1)
    model.setReactionBounds("S_exchange", lower=-1000.0, upper=1000.0)
    model.getReaction("S_exchange").is_exchange = True

    model.createReaction(
        "B_exchange", reversible=True
    )  # Q2.0 - should 'reversible' be True or False?
    model.createReactionReagent("B_exchange", "B_e", -1)
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


def build_toy_model_fba_B():
    model = cbmpy.CBModel.Model("Organism_B")

    model.createSpecies(
        sid="S_e",
        chemFormula="C11H21N2O7PRS",
        name="S external",
        compartment="e",
    )
    model.createSpecies("A_c", name="A cytosol", compartment="c")
    model.createSpecies("B_c", name="B cytosol", compartment="c")
    model.createSpecies(sid="S_c", name="S cytosol", compartment="c")

    model.createSpecies("X_c", name="X cytosol", compartment="c")
    model.createSpecies("Z_c", name="Z cytosol", compartment="c")

    model.createSpecies("BM_c_B", name="Biomass c species B", compartment="c")

    model.createSpecies("A_e", name="A external", compartment="e")

    model.createSpecies("B_e", name="B external", compartment="e")

    model.createSpecies(
        "BM_e_B", name="Biomass external species B", compartment="e"
    )

    model.createCompartment("c", "cytosol")
    model.createCompartment("e", "external")

    model.createReaction(
        "R_1", "Import S external to the cytosol", reversible=True
    )
    model.createReactionReagent("R_1", "S_e", -1)
    model.createReactionReagent("R_1", "S_c", 1)
    model.setReactionBounds("R_1", -1000.0, 1000.0)

    model.createReaction("R_2", "Create X from S_c", reversible=True)
    model.createReactionReagent("R_2", "S_c", -4)
    model.createReactionReagent("R_2", "X_c", 3)
    model.setReactionBounds("R_2", -1000.0, 1000.0)

    model.createReaction("R_3", "Create BM_c_B from X and A", reversible=False)
    model.createReactionReagent("R_3", "X_c", -7)
    model.createReactionReagent("R_3", "A_c", -3)
    model.createReactionReagent("R_3", "BM_c_B", 1)
    model.setReactionBounds("R_3", 0, 1000.0)

    model.createReaction("R_4", "Create Z from A and S", reversible=True)
    model.createReactionReagent("R_4", "A_c", -1)
    model.createReactionReagent("R_4", "S_c", -2)
    model.createReactionReagent("R_4", "Z_c", 1)
    model.setReactionBounds("R_4", -1000.0, 1000.0)

    model.createReaction("R_5", "Create BM_c_B and B from Z", reversible=False)
    model.createReactionReagent("R_5", "Z_c", -3)
    model.createReactionReagent("R_5", "B_c", 1)
    model.createReactionReagent("R_5", "BM_c_B", 1)
    model.setReactionBounds("R_5", -0.0, 1000.0)

    # Import B
    model.createReaction("R_export_B", "Export B_c to B_e", reversible=True)
    model.createReactionReagent("R_export_B", "B_c", -1)
    model.createReactionReagent("R_export_B", "B_e", 1)
    model.setReactionBounds("R_export_B", -1000.0, 1000.0)

    model.createReaction("R_import_A", "Import A_c from A_e", reversible=True)
    model.createReactionReagent("R_import_A", "A_e", -1)
    model.createReactionReagent("R_import_A", "A_c", 1)
    model.setReactionBounds("R_import_A", -1000.0, 1000.0)

    model.createReaction("R_BM_B", "Biomass organism B", reversible=True)
    model.createReactionReagent("R_BM_B", "BM_c_B", -1)
    model.createReactionReagent("R_BM_B", "BM_e_B", 1)
    model.setReactionBounds("R_BM_B", -1000.0, 1000.0)

    # EXchange reactions:
    model.createReaction("S_e_exchange", reversible=True)
    model.createReactionReagent("S_e_exchange", "S_e", -1)
    model.setReactionBounds("S_e_exchange", lower=-1000.0, upper=1000.0)
    model.getReaction("S_e_exchange").is_exchange = True

    model.createReaction(
        "B_e_exchange", reversible=True
    )  # Q2.0 - should 'reversible' be True or False?
    model.createReactionReagent("B_e_exchange", "B_e", -1)
    model.setReactionBounds("B_e_exchange", lower=-1000.0, upper=1000.0)
    model.getReaction("B_e_exchange").is_exchange = True

    model.createReaction("A_e_exchange", reversible=True)
    model.createReactionReagent("A_e_exchange", "A_e", -1)
    model.setReactionBounds("A_e_exchange", -1000.0, 1000)
    model.getReaction("A_e_exchange").is_exchange = True

    # BM exchange
    model.createReaction("BM_e_B_exchange", reversible=False)
    model.createReactionReagent("BM_e_B_exchange", "BM_e_B", -1)
    model.setReactionBounds("BM_e_B_exchange", 0, 1000)
    model.getReaction("BM_e_B_exchange").is_exchange = True

    model.createObjectiveFunction("BM_e_B_exchange")

    model.setActiveObjective("BM_e_B_exchange_objective")

    return model
