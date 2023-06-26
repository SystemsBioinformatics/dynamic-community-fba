from cbmpy.CBModel import Model, Reaction


def create_joint_fba_model(
    model: Model, biomass_reaction_ids: list[str]
) -> Model:
    model.createSpecies("X_c", False, "The community biomass")

    for bm_id in biomass_reaction_ids:
        reaction: Reaction = model.getReaction(bm_id)
        reaction.createReagent("X_c", 1)

    model.createReaction("Xcomm")
    out: Reaction = model.getReaction("Xcomm")
    out.is_exchange = True
    out.setUpperBound(1000)
    out.setLowerBound(0)
    out.createReagent("X_c", -1)

    model.createObjectiveFunction("Xcomm")

    model.setActiveObjective("Xcomm_objective")

    return model
