import cbmpy


# This is model C from the notebooks
def build_toy_model_fba_D():
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

    model.createSpecies("BM_c_D", name="Biomass c species D", compartment="c")

    model.createSpecies("A_e", name="A external", compartment="e")

    model.createSpecies("B_e", name="B external", compartment="e")

    model.createSpecies(
        "BM_e_D", name="Biomass external species B", compartment="e"
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
    model.createReactionReagent("R_3", "A_c", -4)
    model.createReactionReagent("R_3", "X_c", -9)
    model.createReactionReagent("R_3", "B_c", 1)
    model.createReactionReagent("R_3", "BM_c_D", 1)
    model.setReactionBounds("R_3", 0, 1000.0)

    model.createReaction("R_4", "Create Z from A and S", reversible=True)
    model.createReactionReagent("R_4", "S_c", -2)
    model.createReactionReagent("R_4", "Z_c", 1)
    model.setReactionBounds("R_4", -1000.0, 1000.0)

    model.createReaction("R_5", "Create BM_c_B and B from Z", reversible=False)
    model.createReactionReagent("R_5", "Z_c", -9)
    model.createReactionReagent("R_5", "B_c", -1)
    model.createReactionReagent("R_5", "BM_c_D", 1)
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

    model.createReaction("R_BM_D", "Biomass organism D", reversible=True)
    model.createReactionReagent("R_BM_D", "BM_c_D", -1)
    model.createReactionReagent("R_BM_D", "BM_e_D", 1)
    model.setReactionBounds("R_BM_D", -1000.0, 1000.0)

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
    model.createReaction("BM_e_D_exchange", reversible=False)
    model.createReactionReagent("BM_e_D_exchange", "BM_e_D", -1)
    model.setReactionBounds("BM_e_D_exchange", 0, 1000)
    model.getReaction("BM_e_D_exchange").is_exchange = True

    model.createObjectiveFunction("BM_e_D_exchange")

    model.setActiveObjective("BM_e_D_exchange_objective")

    return model
