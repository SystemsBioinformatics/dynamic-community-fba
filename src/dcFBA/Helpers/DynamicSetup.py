
"""Helper functions to prepare standard cbmpy models for Dynamic FBA."""

from cbmpy.CBModel import Model, Species


def mark_dynamic_species(
    model: Model,
    auto: bool = True,
    include: list[str] | None = None,
    exclude: list[str] | None = None,
    extracellular_compartments: list[str] | None = None,
) -> None:
    """
    Mark species in a CBMPy model as dynamic or quasi-steady-state
    for Dynamic FBA simulations, configuring the custom `dcFBA_dynamic` 
    attribute.

    In dynamic simulations, extracellular metabolites usually accumulate or deplete 
    over time. However, certain metabolites (like O2, CO2, or H2O) are often 
    treated as infinite sinks/sources or assumed to be in a quasi-steady state 
    with the environment (instantly replenished/consumed).

    Dynamic species (s.dcFBA_dynamic == True):
    - accumulate/deplete over time
    - in EndPointFBA, receive inter-timepoint linking reactions
    - require concentration tracking

    Quasi-steady-state species (s.dcFBA_dynamic == False):
    - are assumed to have an infinite/constant supply
    - do not accumulate/deplete over time
    - behave as instantaneous exchanges with the environment at each timepoint. 

    Priority of assignment:
    1. Species explicitly listed in ``include`` are marked dynamic.
    2. Species explicitly listed in ``exclude`` are marked non-dynamic.
    3. If auto=True, remaining species are classified using autodetect_dynamic_species().

    Automatic detection marks species as dynamic if:
    - they belong to a valid extracellular/shared compartment
    - they participate in both:
        - an exchange reaction
        - a non-exchange reaction

    If the model defines a shared pool compartment, only species
    in the pool compartment are automatically considered dynamic.
    Otherwise, species in extracellular compartments are considered.

    Args:
        model (Model): The CBMPy model to modify.
        auto (bool): If True, apply automatic dynamic-species detection.
        include (list[str] | None): Species IDs explicitly marked as dynamic.
        exclude (list[str] | None): Species IDs explicitly marked as quasi-steady-state.
        extracellular_compartments (list[str] | None): List of compartment IDs considered extracellular.

    How to use it:
        from dcFBA.Helpers.DynamicSetup import mark_dynamic_species

        # Prep models
        for model in my_models:
            mark_dynamic_species(model, exclude=["O2_e", "H2O_e"])

        # Pass prepped models to the parallel solver
        dpFBA = DynamicParallelFBA(my_models, biomasses, ...)
    """

    if include is None:
        include = []

    if exclude is None:
        exclude = []

    for sid in model.getSpeciesIds():
        species = model.getSpecies(sid)

        # Explicit inclusion
        if sid in include:
            species.dcFBA_dynamic = True
            continue

        # Explicit exclusion
        if sid in exclude:
            species.dcFBA_dynamic = False
            continue

        # Automatic detection
        if auto:
            autodetect_dynamic_species(
                species,
                model,
                extracellular_compartments,
            )


def autodetect_dynamic_species(
    species: Species | str,
    model: Model,
    extracellular_compartments: list[str] | None = None,
) -> bool:
    """
    Automatically determine whether a species should be treated
    dynamically in Dynamic FBA simulations.

    Dynamic species must:
    - belong to a valid extracellular/shared compartment
    - participate in both:
        - an exchange reaction
        - a non-exchange reaction

    If the model defines a shared pool compartment, only species
    in the pool compartment are automatically considered dynamic.
    Otherwise, species in extracellular compartments are considered.

    The result is written to the species attribute
    ``species.dcFBA_dynamic``.

    Args:
        species (Species | str):
            Species object or species ID.

        model (Model):
            CBMPy model containing the species.

        extracellular_compartments (list[str] | None):
            Compartments considered extracellular when no shared
            pool compartment exists.

    Returns:
        bool:
            True if the species is classified as dynamic,
            False otherwise.
    """

    if isinstance(species, str):
        species = model.getSpecies(species)

    compartment_id = species.getCompartmentId()

    # Shared pool takes precedence
    pool_compartment = getattr(model, "pool_compartment", None)

    if pool_compartment is not None:
        valid_compartments = {pool_compartment}

    else:
        if extracellular_compartments is None:
            extracellular_compartments = ["e", "extracellular"]

        valid_compartments = set(extracellular_compartments)

    # Default assumption
    is_dynamic = False

    # Only extracellular/shared species can become dynamic
    if compartment_id in valid_compartments:

        # Fetch all reactions this species participates in
        rxns = [model.getReaction(rid) for rid in species.isReagentOf()]
        
        # Determine if it has both exchange and non-exchange roles
        has_exchange = any(getattr(r, "is_exchange", False) for r in rxns)
        has_non_exchange = any(not getattr(r, "is_exchange", False) for r in rxns)

        is_dynamic = (has_exchange and has_non_exchange)

    species.dcFBA_dynamic = is_dynamic

    return is_dynamic

def is_dynamic_species(
    species: Species,
    model: Model,
) -> bool:
    """
    Determine whether a species should be treated dynamically.

    Priority:
    1. Explicit `dcFBA_dynamic` attribute
    2. Automatic heuristic fallback

    Args:
        species (Species):
            Species to evaluate.

        model (Model):
            Model containing the species.

    Returns:
        bool:
            True if dynamic.
    """

    if not hasattr(species, "dcFBA_dynamic"):
        raise AttributeError(
            f"Species '{species.getId()}' has no dcFBA_dynamic attribute. "
            "Run mark_dynamic_species() to assign it."
        )

    return species.dcFBA_dynamic


def ensure_dynamic_flags(model: Model) -> None:
    """
    Ensure all species possess a dcFBA_dynamic attribute.

    Missing flags are assigned using the automatic heuristic.
    Existing explicit assignments are preserved.
    """

    for sid in model.getSpeciesIds():
        species = model.getSpecies(sid)

        if not hasattr(species, "dcFBA_dynamic"):
            autodetect_dynamic_species(species, model)