import cbmpy
import math
from cbmpy.CBModel import Model, Reaction, Species, Reagent
from .KineticModel import KineticModel


def dynamic_fba(
    kinetic_model: KineticModel,
    biomass_reaction_id: str,
    time_steps: list[float],
    initial_biomass: float = 0.1,
):
    model = kinetic_model.get_model()
    update_conditions = {}

    # Use the exchange reaction lower bounds as initial concentrations
    # If people want to change it they have to change the model
    # before input, we will use the initial

    for exchange in model.getExchangeReactionIds():
        reaction: Reaction = model.getReaction(exchange)
        species_id: str = reaction.reagents[0].getSpecies()
        update_conditions[species_id] = [-reaction.getLowerBound()]

    update_conditions["biomass"] = [initial_biomass]

    species_reaction_dict: dict[
        Species, list[Reaction]
    ] = create_species_reactions_dict(model)

    import_reactions: list[Reaction] = get_import_reactions(model)

    t_old = time_steps[0]

    for t in time_steps[1:]:
        dt = t - t_old
        Xt = update_conditions["biomass"][-1]
        # Calculate new upper bounds for all the import reaction
        for reaction in import_reactions:
            species_concentrations = []
            reagent: Reagent
            for reagent in reaction.reagents:
                species = model.getSpecies(reagent.getSpecies())
                if species.getCompartmentId() == "e":
                    species_concentrations.append(
                        update_conditions[reagent.getSpecies()][-1]
                    )
            if all(c > 0 for c in species_concentrations):
                reaction_upper_bound = calc_reaction_upper_bound(
                    kinetic_model,
                    reaction.id,
                    reaction.getUpperBound(),
                    min(species_concentrations),
                    Xt,
                    dt,
                )
                reaction.setUpperBound(reaction_upper_bound)

            else:
                reaction.setUpperBound(0)

        solution = cbmpy.doFBA(model)

        if math.isnan(solution):
            break

        FBAsol = model.getSolutionVector(names=True)
        FBAsol = dict(zip(FBAsol[1], FBAsol[0]))

        # calculate new concentrations
        for species_id, ls in species_reaction_dict.items():
            new_concentration = update_conditions[species_id][-1]
            reaction: Reaction
            for reaction_id in ls:
                new_concentration = new_concentration - (
                    FBAsol[reaction_id] * Xt * dt
                )
            # update in dictonary
            update_conditions[species_id].append(new_concentration)
        # Update Biomass

        update_conditions["biomass"].append(
            Xt + (FBAsol[biomass_reaction_id] * Xt * dt)
        )

        t_old = t

    return update_conditions


def get_import_reactions(model: Model) -> list[Reaction]:
    import_reactions: list[Reaction] = []
    for reaction in model.reactions:
        if not reaction.is_exchange:
            for reagent in reaction.reagents:
                species_id: str = reagent.getSpecies()
                species = model.getSpecies(species_id)
                if (
                    species.getCompartmentId() == "e"
                    and reagent.coefficient == -1
                ):
                    import_reactions.append(reaction)
                    break
    return import_reactions


def create_species_reactions_dict(
    model: Model,
) -> dict[Species, list[Reaction]]:
    """Helper function such that we can easily access all reactions that
    use/create a certain species by it's species id.

    Since we are only interested in the external metabolites for dynamic FBA
    We create a dict of all external metabolites and their reactions

    Args:
        model (Model): A cbmpy.CBModel
    """

    dict = {}
    species: Species
    for species in model.species:
        if species.getCompartmentId() == "e":
            dict[species.id] = []
            for reaction in species.isReagentOf():
                if reaction not in model.getExchangeReactionIds():
                    dict[species.id].append(reaction)

    my_dict = {key: value for key, value in dict.items() if len(value) > 0}

    return my_dict


def calc_reaction_upper_bound(kinetic_model: KineticModel, rid, ub, S, X, dt):
    km = kinetic_model.get_reaction_km(rid)
    vmax = kinetic_model.get_reaction_vmax(rid)
    if km is None:
        # Calculate Without kinetics
        v_hat = S / (X * dt)
        if ub <= v_hat:
            return ub
        else:
            return v_hat
    else:
        return vmax * (S / (km + S))


from .build_community_matrix import combine_models
from numpy import linspace

model = combine_models([cbmpy.loadModel("data/bigg_models/e_coli_core.xml")])
import_reactions = []

model.createObjectiveFunction("R_BIOMASS_Ecoli_core_w_GAM")

# {"R_GLCpts": [10, 5]}
kinetic_model: KineticModel = KineticModel(model, {})

ts = linspace(0, 15, 100)
y = dynamic_fba(kinetic_model, "R_BIOMASS_Ecoli_core_w_GAM", ts, 0.1)


# y2 = dynamic_fba(
#     kinetic_model,
#     "R_BIOMASS_Ecoli_core_w_GAM",
#     ts,
#     {
#         "M_co2_e": [1],
#         "M_h_e": [1],
#         "M_h2o_e": [1],
#         "M_nh4_e": [1],
#         "M_o2_e": [1],
#         "M_pi_e": [1],
#         "M_glc__D_e": [10],
#         "biomass": [0.1],
#     },
# )


import matplotlib.pyplot as plt

# Assuming y is the DataFrame containing the data for y1, y2, and y3
y1 = y["M_gln__L_e"]
y2 = y["biomass"]
y3 = y["M_glc__D_e"]
ts = ts[0 : len(y1)]

print(y3)
print(y2)

fig, ax1 = plt.subplots()

# Plotting y2 (biomass) on the first y-axis
ax1.plot(ts, y2, color="b")
ax1.set_xlabel("Time")
ax1.set_ylabel("Biomass", color="b")
ax1.tick_params(axis="y", labelcolor="b")

# Creating the second y-axis for y1 (M_gln__L_e)
# ax2 = ax1.twinx()
# ax2.plot(ts, y1, color="r")
# ax2.set_ylabel("M_gln__L_e", color="r")
# ax2.tick_params(axis="y", labelcolor="r")

# Creating the third y-axis for y3 (M_glc__D_e)
ax3 = ax1.twinx()
ax3.spines["right"].set_position(("outward", 60))
ax3.plot(ts, y3, color="g")
ax3.set_ylabel("M_glc__D_e", color="g")
ax3.tick_params(axis="y", labelcolor="g")

plt.show()
