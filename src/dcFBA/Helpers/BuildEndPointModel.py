"""Helper function for the creation of the EndPointFBA model
For further information see DynamicModels.EndPointFBA
"""
import numpy
import re
from cbmpy.CBModel import (
    Model,
    Compartment,
    Reaction,
    Reagent,
    Species,
    FluxBound,
)

from ..Models.CommunityModel import CommunityModel

import weakref


def build_time_model(cm: CommunityModel, times: list[str]) -> CommunityModel:
    """
    Build a time-dependent CommunityModel based on the initial CommunityModel.

    Args:
        initial_model (CommunityModel): The initial CommunityModel to be used as a base.
        times (list): List of the ids which are used to to identify the
            reactions and species of the different time instances

    Returns:
        CommunityModel: The final time-dependent CommunityModel.

    """
    initial_model: CommunityModel = cm.clone()

    # strip the intial model gene protein associations for the EndPointModel
    # TODO SetUpperBOund for inactive genes!!!
    initial_model.gpr = None
    initial_model.genes = None
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
    add_biomass_species(initial_model)

    set_exchanges(initial_model, final_model, times)

    add_time_points(initial_model, final_model, times)

    return final_model


def add_time_points(src_model, target_model, times):
    for i, time in enumerate(times):
        add_time_point(src_model, target_model, time)

        # Check if it's not the last element to avoid index out of range error
        if i < len(times) - 1:
            add_time_link(target_model, times[i], times[i + 1])


def add_time_point(
    src_model: CommunityModel, target_model: CommunityModel, time_id
):
    # start_time = time.time()
    add_time_compartments(src_model, target_model, time_id)

    # start_time = time.time()
    add_reactions(src_model, target_model, time_id)

    # start_time = time.time()

    copy_species_and_reagents(src_model, target_model, time_id)


def set_exchanges(
    initial_model: CommunityModel, final_model: CommunityModel, times
):
    # Exchanges only exchange through the first time point
    for exchange in initial_model.getExchangeReactionIds():
        reaction: Reaction = initial_model.getReaction(exchange)

        # Exchange only has one species
        sid = reaction.getSpeciesIds()[0]
        species: Species = initial_model.getSpecies(sid)
        final_model.createReaction(
            exchange, reaction.name, reversible=True, silent=True
        )
        new_reaction = final_model.getReaction(exchange)
        new_reaction.setLowerBound(reaction.getLowerBound())
        new_reaction.setUpperBound(reaction.getUpperBound())
        new_reaction.is_exchange = True

        new_species = species.clone()
        new_species.setId(f"{sid}_{times[0]}")

        new_species.setCompartmentId(
            f"{species.getCompartmentId()}_{times[0]}"
        )
        new_reaction.createReagent(new_species.getId(), -1)

        # Also create an exchange in the final time point. Such that
        # metabolites can flow out of the system
        add_final_exchange(final_model, new_reaction, times[-1])


def add_time_compartments(
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


def add_reactions(
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
            new_id = rid + "_" + time_id
            # start_time = time.time()
            new_reaction: Reaction = Reaction(
                new_id, reaction.name, reaction.reversible
            )
            final_model.addReaction(
                new_reaction,
                create_default_bounds=False,
                silent=True,
            )

            # Dirty work around, should be fixed in cbmpy 0.9.0
            # Than just use: final_model.createReactionBounds(new_id, reaction.getLowerBound(), reaction.getUpperBound())
            boundId = "%s_%s_bnd" % (new_id, "lower")
            flux = FluxBound(
                boundId, new_id, "greaterEqual", reaction.getLowerBound()
            )
            final_model.__pushGlobalId__(flux.getId(), flux)

            flux.__objref__ = weakref.ref(final_model)
            final_model.flux_bounds.append(flux)
            boundId = "%s_%s_bnd" % (new_id, "upper")
            flux = FluxBound(
                boundId, new_id, "lessEqual", reaction.getUpperBound()
            )
            final_model.__pushGlobalId__(flux.getId(), flux)

            flux.__objref__ = weakref.ref(final_model)
            final_model.flux_bounds.append(flux)


def add_biomass_species(initial_model: CommunityModel) -> None:
    initial_model.createSpecies(
        "BM_c", False, "The community biomass", compartment="e"
    )

    for mid, biomass_id in initial_model.get_model_biomass_ids().items():
        reaction: Reaction = initial_model.getReaction(biomass_id)
        # Create one model biomass and one community biomass
        initial_model.createSpecies(
            f"BM_{mid}", False, f"Biomass of {mid}", compartment="e"
        )
        reaction.createReagent(f"BM_{mid}", 1)
        # Add community biomass
        reaction.createReagent("BM_c", 1)

        exchange_reaction = Reaction(
            f"BM_{mid}_exchange", f"Exchange of biomass {mid}", reversible=True
        )
        exchange_reaction.is_exchange = True
        exchange_reaction.createReagent(f"BM_{mid}", -1)

        initial_model.addReaction(exchange_reaction, False, silent=True)
        exchange_reaction.setLowerBound(0)
        exchange_reaction.setUpperBound(numpy.inf)


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
        new_id = sid + "_" + time_id
        new_species = species.clone()
        new_species.setId(new_id)
        new_species.setCompartmentId(
            species.getCompartmentId() + "_" + time_id
        )

        final_model.addSpecies(new_species)

        for rid in species.isReagentOf():
            old_reaction: Reaction = initial_model.getReaction(rid)
            if not old_reaction.is_exchange:
                new_reaction: Reaction = final_model.getReaction(
                    rid + "_" + time_id
                )
                reagent: Reagent = old_reaction.getReagentWithSpeciesRef(sid)

                new_reaction.createReagent(
                    new_species.getId(), reagent.coefficient
                )


def add_time_link(model: CommunityModel, time0, time1):
    for sid in model.getSpeciesIds():
        if "final" not in sid:
            old_id = re.match(r"(.*?)_time\d+", sid).group(1)

            species: Species = model.getSpecies(sid)
            if species.getCompartmentId() == f"e_{time0}":
                # TODO maybe prefix with LINK_?
                rid = f"{sid}_{time1}"
                linking_reaction = model.createReaction(
                    rid, reversible=False, silent=True
                )

                linking_reaction: Reaction = model.getReaction(rid)
                linking_reaction.createReagent(sid, -1)
                linking_reaction.createReagent(f"{old_id}_{time1}", 1)
                linking_reaction.setLowerBound(0)
                linking_reaction.setUpperBound(numpy.inf)


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
    old_id = re.match(r"(.*?)_time\d+", sid).group(1)
    id = old_id + "_exchange_final"

    # Irreversible, no new species will be imported in the final time step
    final_model.createReaction(
        id,
        "final exchange " + exchange_reaction.id,
        reversible=False,
        silent=True,
    )

    final_exchange = final_model.getReaction(id)
    final_exchange.setLowerBound(0)
    # TODO discuss Do we want to set this this to Inf or to the upperbound?
    # final_exchange.setUpperBound(max(exchange_reaction.getUpperBound(), 0))
    final_exchange.setUpperBound(numpy.inf)

    final_exchange.createReagent(f"{old_id}_{time_id}", -1)
    final_exchange.is_exchange = True
