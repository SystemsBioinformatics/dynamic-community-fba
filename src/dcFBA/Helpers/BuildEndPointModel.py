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

from .DynamicSetup import ensure_dynamic_flags, is_dynamic_species

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

    # ensure the initial model has the dynamic metabolites flagged
    # TODO: double check that flags are preserved with cloning
    ensure_dynamic_flags(initial_model)

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

    final_model._custom_model_identifiers = list(
        initial_model.custom_model_identifiers
    )
    final_model._single_model_biomass_reaction_ids = list(
        initial_model.single_model_biomass_reaction_ids
    )

    final_model._single_model_ids = list(initial_model.single_model_ids)

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
    initial_model: CommunityModel,
    final_model: CommunityModel,
    times: list[str],
) -> None:
    """
    Create the initial exchange pool and final sink reactions for
    dynamically tracked metabolites.

    Dynamic metabolites:
    - receive an initial exchange reaction at the first time point
    - receive a final irreversible sink exchange at the last time point
    - accumulate/deplete through LINK reactions between time points

    Quasi-steady-state metabolites:
    - do not receive dynamic pool exchanges
    - do not accumulate over time
    - are instead duplicated as ordinary exchange reactions at every
      time point in add_reactions()

    Dynamic behavior is determined by is_dynamic_species().

    Args:
        initial_model (CommunityModel):
            The source community model.

        final_model (CommunityModel):
            The time-expanded model being constructed.

        times (list[str]):
            Ordered list of time identifiers.
    """

    # Exchanges only exchange through the first time point
    for exchange in initial_model.getExchangeReactionIds():
        reaction: Reaction = initial_model.getReaction(exchange)

        # Exchange only has one species
        species_ids = reaction.getSpeciesIds()
        assert len(species_ids) == 1, (
            f"Exchange reaction {exchange} must contain exactly one species."
        )
        sid = species_ids[0]
        species: Species = initial_model.getSpecies(sid)

        # Skip quasi-steady-state species
        # If it is not dynamic, we treat it as a regular reaction (handled in add_reactions)
        if not is_dynamic_species(species, initial_model):
            continue

        final_model.createReaction(
            exchange,
            reaction.name,
            reversible=True,
            silent=True,
        )

        new_reaction = final_model.getReaction(exchange)
        new_reaction.setLowerBound(reaction.getLowerBound())
        new_reaction.setUpperBound(reaction.getUpperBound())
        new_reaction.is_exchange = True

        # new_species = species.clone()
        new_species_id =f"{sid}_{times[0]}"
        # new_species.setId(new_species_id)
        # new_species.setCompartmentId(
        #     f"{species.getCompartmentId()}_{times[0]}"
        # )

        # if new_species_id not in final_model.getSpeciesIds():
        #     final_model.addSpecies(new_species)

        new_reaction.createReagent(new_species_id, -1)

        # Create irreversible sink at final time point
        add_final_exchange(
            final_model,
            new_reaction,
            times[-1],
        )


def add_time_compartments(
    initial_model: Model, final_model: Model, time_id: str
) -> None:
    """
    Copy compartments from the initial model to the final model with the
    specified time_id.

    Args:
        initial_model (Model): The source CBModel.
        final_model (Model): The target CBModel
        time_id (str): The time identifier to be appended to the compartment
            IDs.

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

        # Determine if this is a dynamic exchange
        is_dynamic_exchange = False

        if getattr(reaction, "is_exchange", False):

            species_ids = reaction.getSpeciesIds()

            assert len(species_ids) == 1, (
                f"Exchange reaction {rid} must contain exactly one species."
            )

            sid = species_ids[0]
            species = initial_model.getSpecies(sid)

            is_dynamic_exchange = is_dynamic_species(
                species,
                initial_model,
            )

        # If it's a standard internal reaction OR a quasi-steady state exchange,
        # we duplicate it at every time step
        if not is_dynamic_exchange:
            new_id = f"{rid}_{time_id}"
            
            new_reaction: Reaction = Reaction(
                new_id, 
                reaction.name, 
                reaction.reversible
            )
            
            # Retain the exchange flag for quasi-steady state exchanges
            if getattr(reaction, "is_exchange", False):
                new_reaction.is_exchange = True

            final_model.addReaction(
                new_reaction, 
                create_default_bounds=False, 
                silent=True,
            )
            
            final_model.createReactionBounds(
                new_id,
                reaction.getLowerBound(),
                reaction.getUpperBound(),
            )


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

        # if new_id not in final_model.getSpeciesIds():
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


def add_time_link(model: CommunityModel, time0: str, time1: str) -> None:
    """
    Creates linking reactions that carry dynamic metabolites from one time point to the next.

    Dynamic ans quasi-steady-state species are determined using is_dynamic_species().
    Species classified as dynamic:
    - receive LINK reactions
    - accumulate over time
    while those clasisfied quasi-steady-state are instead treated as independent exchange fluxes 
    at each time point.
    
    Advanced users can override this behavior by setting the `dcFBA_dynamic` attribute on a 
    Species object (e.g., using DynamicSetup.mark_dynamic_species()). 
    - If `dcFBA_dynamic` is True, a linking reaction is created regardless of compartment.
    - If `dcFBA_dynamic` is False, the species is skipped (useful for putting species like 
      O2 in a quasi-steady state where they do not accumulate).

    Args:
        model (CommunityModel): The time-expanded community model being built.
        time0 (str): The current time step identifier (e.g., 'time0').
        time1 (str): The next time step identifier (e.g., 'time1').
    """
    
    for sid in model.getSpeciesIds():

        # Final sink species are not propagated forward
        if "final" in sid:
            continue

        match = re.match(r"(.*?)_time\d+", sid)
        if not match:
            continue
        
        old_id = match.group(1)
        species: Species = model.getSpecies(sid)

        # Skip quasi-steady-state species
        if not is_dynamic_species(species, model):
            continue

        # Create the reaction that moves the metabolite from time0 to time1
        if sid.endswith(time0):
            rid = f"LINK_{sid}_{time1}"
            model.createReaction(
                rid, reversible=False, silent=True
            )

            linking_reaction = model.getReaction(rid)
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


