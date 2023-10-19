"""Helper function which loops over all reactions in a Model or CommunityModel. Sets the reaction as the objective function and checks, given a medium if the reaction can take place
"""
import math
import cbmpy
from cbmpy.CBModel import Model, Reaction


def scan_unused_reactions(
    model: Model, medium: dict[str, float], inplace=False
) -> list[str]:
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
        if reaction.getUpperBound() == cbmpy.INF:
            reaction.setUpperBound(1)

        model.createObjectiveFunction(rid, osense="maximize")
        model.setActiveObjective(f"{rid}_objective")

        solution = cbmpy.doFBA(model, quiet=True)

        model.deleteObjective(f"{rid}_objective")
        reaction.setUpperBound(cbmpy.INF)

        if not (math.isnan(solution) or solution == 0):
            continue

        # Check reverse flux possibility
        if reaction.getLowerBound() == cbmpy.NINF:
            reaction.setLowerBound(-1)

        model.createObjectiveFunction(rid, osense="minimize")
        model.setActiveObjective(f"{rid}_objective")

        solution = cbmpy.doFBA(model, quiet=True)

        model.deleteObjective(f"{rid}_objective")

        reaction.setLowerBound(cbmpy.NINF)

        if not (math.isnan(solution) or solution == 0):
            continue

        if inplace:
            model.deleteReactionAndBounds(rid)
        unused_reactions.append(rid)
    return unused_reactions
