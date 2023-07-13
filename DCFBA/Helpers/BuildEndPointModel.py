import re
from cbmpy.CBModel import Model, Compartment, Reaction, Reagent, Species
from ..Models.CommunityModel import CommunityModel


def build_time_model(initial_model: CommunityModel, times):
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

        add_final_exchange(final_model, new_reaction, times[-1])

    for time_id in times:
        copy_compartments(initial_model, final_model, time_id)
        copy_reactions(initial_model, final_model, time_id)
        copy_species_and_reagents(initial_model, final_model, time_id)

    create_time_links(final_model, times)
    return final_model


def copy_compartments(initial_model: Model, final_model: Model, time_id: str):
    for cid in initial_model.getCompartmentIds():
        c: Compartment = initial_model.getCompartment(cid)
        final_model.createCompartment(
            f"{cid}{time_id}", c.name, c.size, c.dimensions, c.volume
        )


def copy_reactions(initial_model: Model, final_model: Model, time_id: str):
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
):
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


def create_time_links(final_model: CommunityModel, times):
    index = 0
    t = times[index]

    while t != times[-1]:
        for sid in final_model.getSpeciesIds():
            old_id = re.match(r"(.*?)_time\d+", sid).group(1)

            species: Species = final_model.getSpecies(sid)
            if (
                species.getCompartmentId() == f"e{t}"
                or old_id in final_model.m_single_model_biomass_reaction_ids
            ):
                rid = f"{sid}{times[index+1]}"
                linking_reaction = final_model.createReaction(
                    rid,
                    reversible=False,
                )

                linking_reaction = final_model.getReaction(rid)
                linking_reaction.createReagent(sid, -1)
                linking_reaction.createReagent(f"{old_id}{times[index+1]}", 1)
                linking_reaction.setUpperBound(1e10)

        index += 1
        t = times[index]


def add_final_exchange(
    final_model: CommunityModel, exchange_reaction: Reaction, time_id: str
):
    # Exchanges only have one species
    sid = exchange_reaction.getSpeciesIds()[0]
    print(sid)
    print(time_id)
    id = exchange_reaction.getId() + "_final"

    # Irreversible, we don't want to import new species in the final step
    final_model.createReaction(
        id, "final exchange " + exchange_reaction.id, False
    )

    final_exchange = final_model.getReaction(id)

    final_exchange.setUpperBound(exchange_reaction.getUpperBound())

    final_exchange.createReagent(re.sub(r"_time\d*", time_id, sid), -1)
