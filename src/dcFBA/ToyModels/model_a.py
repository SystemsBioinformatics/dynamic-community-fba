import cbmpy


def build_model_A():
    model = cbmpy.CBModel.Model("Organism_A")

    model.createSpecies(
        sid="S_e",
        chemFormula="C11H21N2O7PRS",
        name="S external",
        compartment="e",
    )

    model.createSpecies(sid="S_c", name="S cytosol", compartment="c")

    model.createSpecies("Z_c", name="Z cytosol", compartment="c")

    model.createSpecies("D_c", name="D cytosol", compartment="c")

    model.createSpecies("D_e", name="D external", compartment="e")

    model.createSpecies("X_c", name="X cytosol", compartment="c")

    model.createSpecies("Y_c", name="Y cytosol", compartment="c")

    model.createSpecies("A_e", name="A external", compartment="e")

    model.createSpecies("B_c", name="B cytosol", compartment="c")

    model.createSpecies("B_e", name="B external", compartment="e")

    model.createSpecies(
        "BM_c_A", name="Biomass cytosol species A", compartment="c"
    )

    model.createSpecies(
        "BM_e_A", name="Biomass external species A", compartment="e"
    )

    model.createCompartment("c", "cytosol")
    model.createCompartment("e", "external")

    model.createReaction(
        "R_1",
        "Import S external to the cytosol",
        reversible=True,
        create_default_bounds=False,
    )
    model.createReactionReagent("R_1", "S_e", -1)
    model.createReactionReagent("R_1", "S_c", 1)
    model.setReactionBounds("R_1", -1000.0, 1000.0)

    model.createReaction(
        "R_2",
        "Create X from S_c",
        create_default_bounds=False,
        reversible=True,
    )
    model.createReactionReagent("R_2", "S_c", -2)
    model.createReactionReagent("R_2", "X_c", 1)
    model.setReactionBounds("R_2", -1000.0, 1000.0)

    model.createReaction(
        "R_3",
        "Create Y_c and A_e from S_c",
        create_default_bounds=False,
        reversible=True,
    )
    model.createReactionReagent("R_3", "X_c", -1)
    model.createReactionReagent("R_3", "Y_c", 1)
    model.createReactionReagent("R_3", "A_e", 2)
    model.setReactionBounds("R_3", -1000.0, 1000.0)

    model.createReaction(
        "R_4",
        "Create Biomass from Y_c and B_c",
        reversible=True,
        create_default_bounds=False,
    )
    model.createReactionReagent("R_4", "Y_c", -10)
    model.createReactionReagent("R_4", "B_c", -1)
    model.createReactionReagent("R_4", "BM_c_A", 1)
    model.setReactionBounds("R_4", -1000.0, 1000.0)

    model.createReaction(
        "R_5",
        "Create Z from S_int",
        reversible=True,
        create_default_bounds=False,
    )
    model.createReactionReagent("R_5", "S_c", -1)
    model.createReactionReagent("R_5", "Z_c", 1)
    model.createReactionReagent("R_5", "D_c", 1)

    model.setReactionBounds("R_5", -1000.0, 1000.0)

    model.createReaction(
        "R_6",
        "Create BM_c_A from Z and B",
        reversible=True,
        create_default_bounds=False,
    )
    model.createReactionReagent("R_6", "Z_c", -8)
    model.createReactionReagent("R_6", "B_c", -3)
    model.createReactionReagent("R_6", "BM_c_A", 1)
    model.setReactionBounds("R_6", -1000.0, 1000.0)

    # Import B
    model.createReaction(
        "R_7",
        "Import B_e to B_c",
        reversible=True,
        create_default_bounds=False,
    )
    model.createReactionReagent("R_7", "B_e", -1)
    model.createReactionReagent("R_7", "B_c", 1)
    model.setReactionBounds("R_7", -1000.0, 1000.0)

    model.createReaction(
        "R_8",
        "Export D_c to D_e",
        reversible=True,
        create_default_bounds=False,
    )
    model.createReactionReagent("R_8", "D_c", -1)
    model.createReactionReagent("R_8", "D_e", 1)
    model.setReactionBounds("R_8", -1000.0, 1000.0)

    # Biomass creation
    model.createReaction(
        "R_BM_A",
        "Biomass reaction of species A",
        reversible=True,
        create_default_bounds=False,
    )
    model.createReactionReagent("R_BM_A", "BM_c_A", -1)
    model.createReactionReagent("R_BM_A", "BM_e_A", 1)
    model.setReactionBounds("R_BM_A", -1000.0, 1000.0)

    # EXchange reactions:
    model.createReaction(
        "S_exchange", reversible=True, create_default_bounds=False
    )
    model.createReactionReagent("S_exchange", "S_e", -1)
    model.setReactionBounds("S_exchange", lower=-1000.0, upper=1000.0)
    model.getReaction("S_exchange").is_exchange = True

    model.createReaction(
        "B_exchange",
        reversible=True,
        create_default_bounds=False,
    )
    model.createReactionReagent("B_exchange", "B_e", -1)
    model.setReactionBounds("B_exchange", lower=-1000.0, upper=1000.0)
    model.getReaction("B_exchange").is_exchange = True

    model.createReaction(
        "A_exchange", reversible=True, create_default_bounds=False
    )
    model.createReactionReagent("A_exchange", "A_e", -1)
    model.setReactionBounds("A_exchange", -1000.0, 1000)
    model.getReaction("A_exchange").is_exchange = True

    # D can't be imported
    model.createReaction(
        "D_exchange", reversible=False, create_default_bounds=False
    )
    model.createReactionReagent("D_exchange", "D_e", -1)
    model.setReactionBounds("D_exchange", 0, 1000)
    model.getReaction("D_exchange").is_exchange = True

    # BM sink
    model.createReaction(
        "BM_e_A_exchange", reversible=False, create_default_bounds=False
    )
    model.createReactionReagent("BM_e_A_exchange", "BM_e_A", -1)
    model.setReactionBounds("BM_e_A_exchange", 0, 1000)
    model.getReaction("BM_e_A_exchange").is_exchange = True

    cbmpy.CBModel.FluxObjective(
        "objective",
        "BM_e_A_exchange",
        1,
    )

    model.createObjectiveFunction("test_objective")

    model.setObjectiveFlux("BM_e_A_exchange")

    return model


def build_toy_model_fba_A():
    model = cbmpy.CBModel.Model("Organism_A")

    model.createSpecies(
        sid="S_e",
        chemFormula="C11H21N2O7PRS",
        name="S external",
        compartment="e",
    )

    model.createSpecies("A_c", name="A cytosol", compartment="c")
    model.createSpecies("B_c", name="B cytosol", compartment="c")
    model.createSpecies("S_c", name="S cytosol", compartment="c")
    model.createSpecies("X_c", name="X cytosol", compartment="c")
    model.createSpecies("Y_c", name="Y cytosol", compartment="c")
    model.createSpecies("Z_c", name="Z cytosol", compartment="c")

    model.createSpecies(
        "BM_c_A", name="Biomass cytosol species A", compartment="c"
    )

    model.createSpecies("A_e", name="A external", compartment="e")

    model.createSpecies("B_e", name="B external", compartment="e")

    model.createSpecies(
        "BM_e_A", name="Biomass external species A", compartment="e"
    )

    model.createCompartment("c", "cytosol")
    model.createCompartment("e", "external")

    model.createReaction(
        "R_1",
        "Import S external to the cytosol",
        reversible=True,
        create_default_bounds=False,
        silent=True,
    )

    model.createReactionReagent("R_1", "S_e", -1)
    model.createReactionReagent("R_1", "S_c", 1)
    model.createReactionBounds("R_1", -1000.0, 1000.0)

    model.createReaction(
        "R_2",
        "Create X from S_c",
        reversible=True,
        create_default_bounds=False,
        silent=True,
    )
    model.createReactionReagent("R_2", "S_c", -2)
    model.createReactionReagent("R_2", "X_c", 1)
    model.setReactionBounds("R_2", -1000.0, 1000.0)

    model.createReaction(
        "R_3",
        "Create Y_c and A_e from S_c",
        reversible=True,
        create_default_bounds=False,
        silent=True,
    )
    model.createReactionReagent("R_3", "X_c", -1)
    model.createReactionReagent("R_3", "Y_c", 1)
    model.createReactionReagent("R_3", "A_c", 2)
    model.setReactionBounds("R_3", -1000.0, 1000.0)

    model.createReaction(
        "R_4",
        "Create Biomass from Y_c and B_c",
        reversible=False,
        create_default_bounds=False,
        silent=True,
    )
    model.createReactionReagent("R_4", "Y_c", -10)
    model.createReactionReagent("R_4", "B_c", -1)
    model.createReactionReagent("R_4", "BM_c_A", 1)
    model.setReactionBounds("R_4", 0.0, 1000.0)

    model.createReaction(
        "R_5",
        "Create Z from S_int",
        reversible=True,
        create_default_bounds=False,
        silent=True,
    )
    model.createReactionReagent("R_5", "S_c", -1)
    model.createReactionReagent("R_5", "Z_c", 1)
    model.setReactionBounds("R_5", -1000.0, 1000.0)

    model.createReaction(
        "R_6",
        "Create BM_c_A from Z and B",
        reversible=False,
        create_default_bounds=False,
        silent=True,
    )
    model.createReactionReagent("R_6", "Z_c", -8)
    model.createReactionReagent("R_6", "B_c", -3)
    model.createReactionReagent("R_6", "BM_c_A", 1)
    model.setReactionBounds("R_6", 0.0, 1000.0)

    # Import B
    model.createReaction(
        "R_import_B",
        "Import B_e to B_c",
        reversible=True,
        create_default_bounds=False,
        silent=True,
    )
    model.createReactionReagent("R_import_B", "B_e", -1)
    model.createReactionReagent("R_import_B", "B_c", 1)
    model.setReactionBounds("R_import_B", -1000.0, 1000.0)

    # Import A
    model.createReaction(
        "R_export_A",
        "Export A_c to A_e",
        reversible=True,
        create_default_bounds=False,
        silent=True,
    )
    model.createReactionReagent("R_export_A", "A_c", -1)
    model.createReactionReagent("R_export_A", "A_e", 1)
    model.setReactionBounds("R_export_A", -1000.0, 1000.0)

    # Import B
    model.createReaction(
        "R_BM_A",
        "Biomass reaction",
        reversible=True,
        create_default_bounds=False,
        silent=True,
    )
    model.createReactionReagent("R_BM_A", "BM_c_A", -1)
    model.createReactionReagent("R_BM_A", "BM_e_A", 1)
    model.setReactionBounds("R_BM_A", -1000.0, 1000.0)

    # EXchange reactions:
    model.createReaction(
        "S_e_exchange",
        reversible=True,
        create_default_bounds=False,
        silent=True,
    )
    model.createReactionReagent("S_e_exchange", "S_e", -1)
    model.setReactionBounds("S_e_exchange", lower=-1000.0, upper=1000.0)
    model.getReaction("S_e_exchange").is_exchange = True

    model.createReaction(
        "B_e_exchange",
        reversible=True,
        create_default_bounds=False,
        silent=True,
    )  # Q2.0 - should 'reversible' be True or False?
    model.createReactionReagent("B_e_exchange", "B_e", -1)
    model.setReactionBounds("B_e_exchange", lower=-1000.0, upper=1000.0)
    model.getReaction("B_e_exchange").is_exchange = True

    model.createReaction(
        "A_e_exchange",
        reversible=True,
        create_default_bounds=False,
        silent=True,
    )
    model.createReactionReagent("A_e_exchange", "A_e", -1, silent=True)
    model.setReactionBounds("A_e_exchange", -1000.0, 1000)
    model.getReaction("A_e_exchange").is_exchange = True

    # BM sink
    model.createReaction(
        "BM_e_A_exchange",
        reversible=False,
        create_default_bounds=False,
        silent=True,
    )
    model.createReactionReagent("BM_e_A_exchange", "BM_e_A", -1, silent=True)
    model.setReactionBounds("BM_e_A_exchange", 0, 1000)
    model.getReaction("BM_e_A_exchange").is_exchange = True

    model.createObjectiveFunction("BM_e_A_exchange")

    model.setActiveObjective("BM_e_A_exchange_objective")

    return model
