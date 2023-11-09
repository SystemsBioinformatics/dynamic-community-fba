import cbmpy


def build_combined_model():
    model = cbmpy.CBModel.Model("combined_model")

    model.createSpecies(
        sid="S_e",
        chemFormula="C11H21N2O7PRS",
        name="S external",
        compartment="e",
    )

    model.createSpecies(
        sid="S_c_Organism_A", name="S cytosol", compartment="c_Organism_A"
    )
    model.createSpecies(
        sid="S_c_Organism_B", name="S cytosol", compartment="c_Organism_B"
    )

    model.createSpecies(
        "Z_c_Organism_A", name="Z cytosol", compartment="c_Organism_A"
    )
    model.createSpecies(
        "Z_c_Organism_B", name="Z cytosol", compartment="c_Organism_B"
    )

    model.createSpecies("D_c", name="D cytosol", compartment="c_Organism_A")
    model.createSpecies("D_e", name="D external", compartment="e")

    model.createSpecies("E_c", name="D cytosol", compartment="c_Organism_B")
    model.createSpecies("E_e", name="D external", compartment="e")

    model.createSpecies(
        "X_c_Organism_A", name="X cytosol", compartment="c_Organism_A"
    )
    model.createSpecies(
        "X_c_Organism_B", name="X cytosol", compartment="c_Organism_B"
    )

    model.createSpecies("Y_c", name="Y cytosol", compartment="c_Organism_A")

    model.createSpecies("A_c", name="A cytosol", compartment="c_Organism_B")
    model.createSpecies("A_e", name="A external", compartment="e")

    model.createSpecies(
        "B_c_Organism_A", name="B cytosol", compartment="c_Organism_A"
    )
    model.createSpecies(
        "B_c_Organism_B", name="B cytosol", compartment="c_Organism_B"
    )

    model.createSpecies("B_e", name="B external", compartment="e")

    model.createSpecies(
        "BM_c_A", name="Biomass cytosol species A", compartment="c_Organism_A"
    )

    model.createSpecies(
        "BM_e_A", name="Biomass external species A", compartment="e"
    )

    model.createSpecies(
        "BM_c_B", name="Biomass cytosol species A", compartment="c_Organism_B"
    )

    model.createSpecies(
        "BM_e_B", name="Biomass external species A", compartment="e"
    )

    # Model C species
    model.createSpecies(
        "S_c_Organism_C", name="S_c organism C", compartment="c"
    )

    model.createSpecies("F_c", name="F of organism C", compartment="c")

    model.createSpecies(
        "BM_c_C", name="Biomass cytosol species C", compartment="c"
    )

    model.createSpecies(
        "BM_e_C", name="Biomass external species C", compartment="e"
    )

    # Compartments
    model.createCompartment("c_Organism_A", "cytosol")
    model.createCompartment("c_Organism_B", "cytosol")
    model.createCompartment("c_Organism_C", "cytosol")

    model.createCompartment("e", "external")

    # Model A Reactions
    model.createReaction(
        "R_1_Organism_A", "Import S external to the cytosol", reversible=True
    )
    model.createReactionReagent("R_1_Organism_A", "S_e", -1)
    model.createReactionReagent("R_1_Organism_A", "S_c_Organism_A", 1)
    model.setReactionBounds("R_1_Organism_A", -1000.0, 1000.0)

    model.createReaction(
        "R_2_Organism_A", "Create X from S_c", reversible=True
    )
    model.createReactionReagent("R_2_Organism_A", "S_c_Organism_A", -2)
    model.createReactionReagent("R_2_Organism_A", "X_c_Organism_A", 1)
    model.setReactionBounds("R_2_Organism_A", -1000.0, 1000.0)

    model.createReaction(
        "R_3_Organism_A", "Create Y_c and A_e from S_c", reversible=True
    )
    model.createReactionReagent("R_3_Organism_A", "X_c_Organism_A", -1)
    model.createReactionReagent("R_3_Organism_A", "Y_c", 1)
    model.createReactionReagent("R_3_Organism_A", "A_e", 2)
    model.setReactionBounds("R_3_Organism_A", -1000.0, 1000.0)

    model.createReaction(
        "R_4_Organism_A", "Create Biomass from Y_c and B_c", reversible=True
    )
    model.createReactionReagent("R_4_Organism_A", "Y_c", -10)
    model.createReactionReagent("R_4_Organism_A", "B_c_Organism_A", -1)
    model.createReactionReagent("R_4_Organism_A", "BM_c_A", 1)
    model.setReactionBounds("R_4_Organism_A", -1000.0, 1000.0)

    model.createReaction(
        "R_5_Organism_A", "Create Z from S_int", reversible=True
    )
    model.createReactionReagent("R_5_Organism_A", "S_c_Organism_A", -1)
    model.createReactionReagent("R_5_Organism_A", "Z_c_Organism_A", 1)
    model.createReactionReagent("R_5_Organism_A", "D_c", 1)

    model.setReactionBounds("R_5_Organism_A", -1000.0, 1000.0)

    model.createReaction(
        "R_6_Organism_A", "Create BM_c_A from Z and B", reversible=True
    )
    model.createReactionReagent("R_6_Organism_A", "Z_c_Organism_A", -8)
    model.createReactionReagent("R_6_Organism_A", "B_c_Organism_A", -3)
    model.createReactionReagent("R_6_Organism_A", "BM_c_A", 1)
    model.setReactionBounds("R_6_Organism_A", -1000.0, 1000.0)

    # Import B
    model.createReaction(
        "R_7_Organism_A", "Import B_e to B_c", reversible=True
    )
    model.createReactionReagent("R_7_Organism_A", "B_e", -1)
    model.createReactionReagent("R_7_Organism_A", "B_c_Organism_A", 1)
    model.setReactionBounds("R_7_Organism_A", -1000.0, 1000.0)

    model.createReaction(
        "R_8_Organism_A", "Export D_c to D_e", reversible=True
    )
    model.createReactionReagent("R_8_Organism_A", "D_c", -1)
    model.createReactionReagent("R_8_Organism_A", "D_e", 1)
    model.setReactionBounds("R_8_Organism_A", -1000.0, 1000.0)

    # Biomass creation
    model.createReaction(
        "R_BM_A_Organism_A",
        "Biomass reaction of species A",
        reversible=True,
    )
    model.createReactionReagent("R_BM_A_Organism_A", "BM_c_A", -1)
    model.createReactionReagent("R_BM_A_Organism_A", "BM_e_A", 1)
    model.setReactionBounds("R_BM_A_Organism_A", -1000.0, 1000.0)

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
    model.getReaction("S_exchange").is_exchange = True

    # D can't be imported
    model.createReaction("D_exchange", reversible=False)
    model.createReactionReagent("D_exchange", "D_e", -1)
    model.setReactionBounds("D_exchange", 0, 1000)
    model.getReaction("D_exchange").is_exchange = True

    # BM sink
    model.createReaction("BM_e_A_exchange", reversible=False)
    model.createReactionReagent("BM_e_A_exchange", "BM_e_A", -1)
    model.setReactionBounds("BM_e_A_exchange", 0, 1000)
    model.getReaction("BM_e_A_exchange").is_exchange = True

    # Model B reactions

    model.createReaction(
        "R_1_Organism_B", "Import S external to the cytosol", reversible=True
    )
    model.createReactionReagent("R_1_Organism_B", "S_e", -1)
    model.createReactionReagent("R_1_Organism_B", "S_c_Organism_B", 1)
    model.setReactionBounds("R_1_Organism_B", -1000.0, 1000.0)

    model.createReaction(
        "R_2_Organism_B", "Create X from S_c", reversible=True
    )
    model.createReactionReagent("R_2_Organism_B", "S_c_Organism_B", -4)
    model.createReactionReagent("R_2_Organism_B", "X_c_Organism_B", 3)
    model.setReactionBounds("R_2_Organism_B", -1000.0, 1000.0)

    model.createReaction(
        "R_9_Organism_B", "Create BM_c_B from X", reversible=True
    )
    model.createReactionReagent("R_9_Organism_B", "X_c_Organism_B", -7)
    model.createReactionReagent("R_9_Organism_B", "A_c", -3)
    model.createReactionReagent("R_9_Organism_B", "E_c", 1)
    model.createReactionReagent("R_9_Organism_B", "BM_c_B", 1)
    model.setReactionBounds("R_9_Organism_B", -1000.0, 1000.0)

    model.createReaction(
        "R_10_Organism_B", "Import A_c from A_e", reversible=True
    )
    model.createReactionReagent("R_10_Organism_B", "A_e", -1)
    model.createReactionReagent("R_10_Organism_B", "A_c", 1)
    model.setReactionBounds("R_10_Organism_B", -1000.0, 1000.0)

    model.createReaction(
        "R_5_Organism_B", "Create Z from S_c", reversible=True
    )
    model.createReactionReagent("R_5_Organism_B", "S_c_Organism_B", -2)
    model.createReactionReagent("R_5_Organism_B", "Z_c_Organism_B", 1)
    model.setReactionBounds("R_5_Organism_B", -1000.0, 1000.0)

    model.createReaction(
        "R_6_Organism_B", "Create BM_c_B and B from Z", reversible=True
    )
    model.createReactionReagent("R_6_Organism_B", "Z_c_Organism_B", -9)
    model.createReactionReagent("R_6_Organism_B", "B_c_Organism_B", 1)
    model.createReactionReagent("R_6_Organism_B", "BM_c_B", 1)
    model.setReactionBounds("R_6_Organism_B", -1000.0, 1000.0)

    # Import B
    model.createReaction(
        "R_11_Organism_B", "Export B_c to B_e", reversible=True
    )
    model.createReactionReagent("R_11_Organism_B", "B_c_Organism_B", -1)
    model.createReactionReagent("R_11_Organism_B", "B_e", 1)
    model.setReactionBounds("R_11_Organism_B", -1000.0, 1000.0)

    model.createReaction(
        "R_12_Organism_B", "Export E_c to E_e", reversible=True
    )
    model.createReactionReagent("R_12_Organism_B", "E_c", -1)
    model.createReactionReagent("R_12_Organism_B", "E_e", 1)
    model.setReactionBounds("R_12_Organism_B", -1000.0, 1000.0)

    # Biomass creation
    model.createReaction(
        "R_BM_B_Organism_B",
        "Biomass reaction of species B",
        reversible=True,
    )
    model.createReactionReagent("R_BM_B_Organism_B", "BM_c_B", -1)
    model.createReactionReagent("R_BM_B_Organism_B", "BM_e_B", 1)
    model.setReactionBounds("R_BM_B_Organism_B", -1000.0, 1000.0)

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

    # Model C reactions
    model.createReaction("R_1_Organism_C", reversible=True)
    model.createReactionReagent("R_1_Organism_C", "S_e", -1)
    model.createReactionReagent("R_1_Organism_C", "S_c_Organism_C", 1)
    model.setReactionBounds("R_1_Organism_C", -1000, 1000)

    model.createReaction("R_13_Organism_C", reversible=True)
    model.createReactionReagent("R_13_Organism_C", "S_c_Organism_C", -2)
    model.createReactionReagent("R_13_Organism_C", "F_c", 3)
    model.setReactionBounds("R_13_Organism_C", -1000, 1000)

    model.createReaction("R_14_Organism_C", reversible=True)
    model.createReactionReagent("R_14_Organism_C", "F_c", -2)
    model.createReactionReagent("R_14_Organism_C", "BM_c_C", 1)
    model.setReactionBounds("R_14_Organism_C", -1000, 1000)

    model.createReaction("R_BM_C_Organism_C", reversible=True)
    model.createReactionReagent("R_BM_C_Organism_C", "BM_c_C", -1)
    model.createReactionReagent("R_BM_C_Organism_C", "BM_e_C", 1)
    model.setReactionBounds("R_BM_C_Organism_C", -1000, 1000)

    model.createReaction("BM_e_C_exchange", reversible=False)
    model.createReactionReagent("BM_e_C_exchange", "BM_e_C", -1)
    model.setReactionBounds("BM_e_C_exchange", 0, 1000)
    model.getReaction("BM_e_C_exchange").is_exchange = True

    return model
