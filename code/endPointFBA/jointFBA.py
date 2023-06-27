from cbmpy.CBModel import Reaction
from endPointFBA.CombinedModel import CombinedModel


def create_joint_fba_model(
    model: CombinedModel,
    biomass_reaction_ids: list[str],
    create_new_model=False,
) -> CombinedModel:
    if create_new_model:
        joint_model = model.clone()
    else:
        joint_model = model
    joint_model.createSpecies("X_c", False, "The community biomass")

    for bm_id in biomass_reaction_ids:
        reaction: Reaction = joint_model.getReaction(bm_id)
        reaction.createReagent("X_c", 1)

    joint_model.createReaction("Xcomm")
    out: Reaction = joint_model.getReaction("Xcomm")
    out.is_exchange = True
    out.setUpperBound(1000)
    out.setLowerBound(0)
    out.createReagent("X_c", -1)

    joint_model.createObjectiveFunction("Xcomm")

    joint_model.setActiveObjective("Xcomm_objective")

    return joint_model
