import cbmpy
import math
from cbmpy.CBModel import Model, Reaction, Species, Reagent
from .KineticModel import KineticModel


def dynamic_fba(
    kinetic_model: KineticModel,
    biomass_reaction_id: str,
    time_steps: list[float],
    initial_biomass: float = 0.1,
    user_func=None,
):
    model = kinetic_model.get_model()
    update_conditions: dict[str, list[float]] = {}

    # Use the exchange reaction lower bounds as initial concentrations
    # If people want to change it they have to change the model
    # before input, we will use the initial

    for exchange in model.getExchangeReactionIds():
        # print(exchange)
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
                    user_func,
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
    for reaction_id in model.getReactionIds():
        reaction: Reaction = model.getReaction(reaction_id)
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


def calc_reaction_upper_bound(
    kinetic_model: KineticModel, rid, ub, S, X, dt, user_func=None
):
    km = kinetic_model.get_reaction_km(rid)
    vmax = kinetic_model.get_reaction_vmax(rid)
    if user_func is not None:
        return user_func(kinetic_model, rid, ub, S, X, dt)
    elif km is None:
        return ub

        # # Calculate Without kinetics
        # # TODO discuss if this should be turned on or off with Francesco
        # v_hat = S / (X * dt)
        # if ub <= v_hat:
        #     return ub
        # else:
        #     return v_hat
    else:
        return vmax * (S / (km + S))
