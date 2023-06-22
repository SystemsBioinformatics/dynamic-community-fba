import cbmpy
import math
from cbmpy.CBModel import Model, Reaction, Species, Reagent
from .KineticModel import KineticModel


def dynamic_fba(
    kinetic_model: KineticModel,
    biomass_reaction_id: str,
    time_steps: list[float],
    initial_conditions: dict[str, list[float]],
):
    if "biomass" not in initial_conditions.keys():
        raise Exception("Biomass must be defined in the initial condition")

    model = kinetic_model.get_model()
    update_conditions = initial_conditions

    # If we have Vmax of a reaction set this to be the upper bound
    # Else just use the model defined upper bound
    for rid, tt in kinetic_model.get_model_kinetics().items():
        reaction: Reaction = model.getReaction(rid)
        reaction.setUpperBound(tt[0])

    species_reaction_dict: dict[
        Species, list[Reaction]
    ] = create_species_reactions_dict(model)

    for species in model.species:
        if species.id not in update_conditions.keys():
            update_conditions[species.id] = [0]

    t_old = time_steps[0]

    for t in time_steps[1:]:
        dt = t - t_old

        Xt = update_conditions["biomass"][-1]
        # Calculate new upper bounds for import reaction
        for reaction in model.reactions:
            if reaction.id not in model.getExchangeReactionIds():
                reagents: list[Reagent] = reaction.reagents
                slowest_reaction = 1e10
                for reagent in reagents:
                    species_id: str = reagent.getSpecies()
                    species = model.getSpecies(species_id)

                    species_concentration = update_conditions[species.id][-1]

                    # Check if the reaction is an importer (imports external
                    # metabolite) something to import set the upper bound
                    if (
                        species.getCompartmentId() == "e"
                        and reagent.coefficient == -1
                        and species_concentration > 0
                    ):
                        reaction_upper_bound = calc_reaction_upper_bound(
                            kinetic_model,
                            reaction.id,
                            reaction.getUpperBound(),
                            update_conditions[species.id][-1],
                            Xt,
                            dt,
                        )
                        # TODO This code is for the fact that if a importer has
                        # two substrates in the external compartment we choose
                        # The speed of the bottleneck species
                        if (
                            reaction_upper_bound < slowest_reaction
                            and reaction_upper_bound != 0
                        ):
                            reaction.setUpperBound(reaction_upper_bound)
                            slowest_reaction = reaction_upper_bound

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
