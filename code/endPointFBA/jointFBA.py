from cbmpy.CBModel import Model, Reaction


def create_joint_fba_model(
    model: Model, biomass_reaction_ids: list[str]
) -> Model:
    new_model = model.clone()
    new_model.createSpecies("X_c", False, "The community biomass")

    for bm_id in biomass_reaction_ids:
        reaction: Reaction = new_model.getReaction(bm_id)
        reaction.createReagent("X_c", 1)

    new_model.createReaction("Xcomm")
    out: Reaction = new_model.getReaction("Xcomm")
    out.is_exchange = True
    out.setUpperBound(1000)
    out.setLowerBound(0)
    out.createReagent("X_c", -1)

    new_model.createObjectiveFunction("Xcomm")

    new_model.setActiveObjective("Xcomm_objective")

    return new_model
