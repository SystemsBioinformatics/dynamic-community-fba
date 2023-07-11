from cbmpy.CBModel import Reaction
from ..Models.CommunityModel import CommunityModel
from ..Helpers.BuildEndPointModel import build_time_model


class EndPointFBA:
    m_model = CommunityModel

    def __init__(
        self,
        community_model: CommunityModel,
        n: int,
        initial_biomasses: dict[str, float],
        dt: float = 0.1,
    ) -> None:
        self.m_model = build_time_model(community_model, n)

        self.set_constraints(n, initial_biomasses, dt)

    def set_constraints(self, n: int, initial_biomasses, dt: float):
        width = len(str(n))
        times = [f"_time{i:0{width}d}" for i in range(n + 1)]

        rids_t0 = self.m_model.getReactionIds(times[0])
        for rid in rids_t0:
            reaction = self.m_model.getReaction(rid)
            mid = self.m_model.identify_model_from_reaction(rid)
            # Skip exchange reactions and time to time reactions
            if reaction.is_exchange or mid == "":
                continue
            biomass = initial_biomasses[mid]
            reaction.setLowerBound(reaction.getLowerBound() * biomass * dt)
            reaction.setUpperBound(reaction.getUpperBound() * biomass * dt)

        # From time 1 to last time point
        for i in range(1, n + 1):
            rids = self.m_model.getReactionIds(times[i])
            for rid in rids:
                reaction: Reaction = self.m_model.getReaction(rid)
                biomass_rid = (
                    self.m_model.identify_biomass_of_model_from_reaction_id(
                        rid
                    )
                )

                if reaction.is_exchange or biomass_rid == "":
                    continue

                r_x_t = biomass_rid + times[i - 1]
                ub = reaction.getUpperBound()
                self.m_model.addUserConstraint(
                    f"{rid}_ub",
                    [
                        [1, rid],
                        [-1 * dt * ub, r_x_t],
                    ],
                    "<=",
                    0.0,
                )
                reaction.setUpperBound(1e10)
