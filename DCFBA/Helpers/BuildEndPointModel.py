"""Helper function for the creation of the EndPointFBA model
For further information see DynamicModels.EndPointFBA
"""

import re
from cbmpy.CBModel import Model, Compartment, Reaction, Reagent, Species
from ..Models.CommunityModel import CommunityModel


def build_time_model(
    initial_model: CommunityModel, times: list[str]
) -> CommunityModel:
    """
    Build a time-dependent CommunityModel based on the initial CommunityModel.

    Args:
        initial_model (CommunityModel): The initial CommunityModel to be used as a base.
        times (list): List of the ids which are used to to identify the
            reactions and species of the different time instances

    Returns:
        CommunityModel: The final time-dependent CommunityModel.

    """

    final_model = CommunityModel(
        [],
        [],
        [],
        "Timed_community_model",
    )

    final_model.m_identifiers = initial_model.m_identifiers
    final_model.m_single_model_biomass_reaction_ids = (
        initial_model.m_single_model_biomass_reaction_ids
    )
    final_model.m_single_model_ids = initial_model.m_single_model_ids

    # Exchanges only exchange through the first time point
    for exchange in initial_model.getExchangeReactionIds():
        reaction: Reaction = initial_model.getReaction(exchange)

        # Exchange only has one species
        sid = reaction.getSpeciesIds()[0]
        species: Species = initial_model.getSpecies(sid)
        final_model.createReaction(exchange, reaction.name, reversible=True)
        new_reaction = final_model.getReaction(exchange)
        new_reaction.setLowerBound(reaction.getLowerBound())
        new_reaction.setUpperBound(reaction.getUpperBound())
        new_reaction.is_exchange = True

        new_species = species.clone()
        new_species.setId(f"{sid}_time0")

        new_species.setCompartmentId(f"{species.getCompartmentId()}_time0")
        new_reaction.createReagent(new_species.getId(), -1)

        # Also create an exchange in the final time point. Such that
        # metabolites can flow out of the system
        add_final_exchange(final_model, new_reaction, times[-1])

    for time_id in times:
        copy_compartments(initial_model, final_model, time_id)
        copy_reactions(initial_model, final_model, time_id)
        copy_species_and_reagents(initial_model, final_model, time_id)

    create_time_links(final_model, times)
    return final_model


def copy_compartments(
    initial_model: Model, final_model: Model, time_id: str
) -> None:
    """
    Copy compartments from the initial model to the final model with the
    specified time_id.

    Args:
        initial_model (Model): The source CBModel.
        final_model (Model): The target CBModel
        time_id (str): The time identifier to be appended to the compartment IDs.

    Returns:
        None

    """
    for cid in initial_model.getCompartmentIds():
        c: Compartment = initial_model.getCompartment(cid)
        final_model.createCompartment(
            f"{cid}{time_id}", c.name, c.size, c.dimensions, c.volume
        )


def copy_reactions(
    initial_model: Model, final_model: Model, time_id: str
) -> None:
    """
    Copy reactions from the initial model to the final model with the
    specified time_id.

    Args:
        initial_model (Model): The source CBModel.
        final_model (Model): The target CBModel where reactions will be copied.
        time_id (str): The time identifier to be appended to the reaction IDs.

    Returns:
        None

    """

    for rid in initial_model.getReactionIds():
        reaction: Reaction = initial_model.getReaction(rid)
        if not reaction.is_exchange:
            new_id = rid + time_id
            final_model.createReaction(
                new_id, reaction.name, reaction.reversible
            )
            new_reaction: Reaction = final_model.getReaction(new_id)

            new_reaction.setLowerBound(reaction.getLowerBound())
            new_reaction.setUpperBound(reaction.getUpperBound())


def copy_species_and_reagents(
    initial_model: Model, final_model: Model, time_id
) -> None:
    """
    Copy species and reagents from the initial model to the final model with
    the specified time_id.

    Args:
        initial_model (Model): The source CBModel.
        final_model (Model): The target CBModel where species and reagents
            will be copied.
        time_id (str): The time identifier to be appended to the species and
            reagent IDs.

    Returns:
        None

    """
    for sid in initial_model.getSpeciesIds():
        species: Species = initial_model.getSpecies(sid)
        new_id = sid + time_id
        new_species = species.clone()
        new_species.setId(new_id)
        new_species.setCompartmentId(species.getCompartmentId() + time_id)

        final_model.addSpecies(new_species)

        for rid in species.isReagentOf():
            old_reaction: Reaction = initial_model.getReaction(rid)
            if not old_reaction.is_exchange:
                new_reaction: Reaction = final_model.getReaction(rid + time_id)
                reagent: Reagent = old_reaction.getReagentWithSpeciesRef(sid)

                new_reaction.createReagent(
                    new_species.getId(), reagent.coefficient
                )


def create_time_links(
    final_model: CommunityModel, time_ids: list[str]
) -> None:
    """
    Create links between species from different time points
    in the CommunityModel.
    Links represent reactions foreach reaction,time point
        => species_time0 -> species_time1

    Args:
        final_model (CommunityModel): The final time-dependent CommunityModel.
        time_ids (list): Ids of the different time points

    Returns:
        None

    """
    index = 0
    t = time_ids[index]

    while t != time_ids[-1]:
        for sid in final_model.getSpeciesIds():
            old_id = re.match(r"(.*?)_time\d+", sid).group(1)

            species: Species = final_model.getSpecies(sid)
            if (
                species.getCompartmentId() == f"e{t}"
                or old_id in final_model.m_single_model_biomass_reaction_ids
            ):
                rid = f"{sid}{time_ids[index+1]}"
                print(rid)
                linking_reaction = final_model.createReaction(
                    rid,
                    reversible=False,
                )

                linking_reaction = final_model.getReaction(rid)
                linking_reaction.createReagent(sid, -1)
                linking_reaction.createReagent(
                    f"{old_id}{time_ids[index+1]}", 1
                )
                linking_reaction.setUpperBound(1e10)

        index += 1
        t = time_ids[index]

    print(final_model.getReactionIds("ecoli_1"))
    raise Exception()


def add_final_exchange(
    final_model: CommunityModel, exchange_reaction: Reaction, time_id: str
) -> None:
    """
    Add a final exchange reaction to the final CommunityModel for a given time_id.

    Args:
        final_model (CommunityModel): The final time-dependent CommunityModel.
        exchange_reaction (Reaction): The exchange reaction from the initial model.
        time_id (str): The time identifier for the final exchange reaction.

    Returns:
        None

    """
    # Exchanges only have one species
    sid = exchange_reaction.getSpeciesIds()[0]
    print(sid)
    print(time_id)
    id = exchange_reaction.getId() + "_final"

    # Irreversible, no new species will be imported in the final time step
    final_model.createReaction(
        id, "final exchange " + exchange_reaction.id, False
    )

    final_exchange = final_model.getReaction(id)

    final_exchange.setUpperBound(exchange_reaction.getUpperBound())

    final_exchange.createReagent(re.sub(r"_time\d*", time_id, sid), -1)
