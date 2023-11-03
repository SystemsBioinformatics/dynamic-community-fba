"""Helper function which loops over all reactions in a Model or CommunityModel. Sets the reaction as the objective function and checks, given a medium if the reaction can take place
"""
import math
import cbmpy
import numpy
from cbmpy.CBModel import Model, Reaction


def scan_unused_reactions(model: Model, medium: dict[str, float]) -> list[str]:
    """Generate list of unused reactions

    Args:
        model (Model): a GSMM
        medium (dict[str, float]): The medium the model is in

    Returns:
        list[str]: All the unused reactions given the medium
    """
    for eid in model.getExchangeReactionIds():
        reaction: Reaction = model.getReaction(eid)
        sid = reaction.getSpeciesIds()[0]
        if sid in medium.keys():
            # Set all the exchange reactions of the medium
            reaction.setLowerBound(-medium[sid])

    unused_reactions = []

    for rid in model.getReactionIds():
        reaction = model.getReaction(rid)

        # Check forward flux possibility
        old_ub = reaction.getUpperBound()
        if old_ub == cbmpy.INF or old_ub == numpy.Inf:
            reaction.setUpperBound(1)

        model.createObjectiveFunction(rid, osense="maximize")
        model.setActiveObjective(f"{rid}_objective")

        solution = cbmpy.doFBA(model, quiet=True)

        model.deleteObjective(f"{rid}_objective")
        reaction.setUpperBound(old_ub)

        if not (math.isnan(solution) or solution == 0):
            continue

        # Check reverse flux possibility
        old_lb = reaction.getLowerBound()
        if old_lb == cbmpy.NINF or old_lb == numpy.NINF:
            reaction.setLowerBound(-1)

        model.createObjectiveFunction(rid, osense="minimize")
        model.setActiveObjective(f"{rid}_objective")

        solution = cbmpy.doFBA(model, quiet=True)

        model.deleteObjective(f"{rid}_objective")

        reaction.setLowerBound(old_lb)

        if not (math.isnan(solution) or solution == 0):
            continue

        unused_reactions.append(rid)
    return unused_reactions


def reduce_model(model: Model, medium: dict[str, float]) -> None:
    """Remove unused reactions and metabolites from the model"""
    rids: list[str] = scan_unused_reactions(model, medium)

    for rid in rids:
        model.deleteReactionAndBounds(rid)

    model.deleteNonReactingSpecies()
