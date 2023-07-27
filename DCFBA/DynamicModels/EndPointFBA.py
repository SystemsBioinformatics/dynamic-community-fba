import cbmpy
from cbmpy.CBModel import Reaction, Species
from ..Models.CommunityModel import CommunityModel
from ..Helpers.BuildEndPointModel import build_time_model
from .DynamicFBABase import DynamicFBABase


class EndPointFBA(DynamicFBABase):
    m_model = CommunityModel
    m_times: list[str]

    def __init__(
        self,
        community_model: CommunityModel,
        n: int,
        initial_biomasses: dict[str, float],
        initial_concentrations: dict[str, float] = {},
        dt: float = 0.1,
    ) -> None:
        width = len(str(n))
        # TODO times should not hold the under score
        m_times = [f"time{i:0{width}d}" for i in range(n)]
        self.m_model = build_time_model(community_model, m_times)

        self.set_constraints(n, initial_biomasses, dt, m_times)
        self.set_initial_concentrations(
            initial_biomasses, initial_concentrations, m_times
        )

    def simulate(
        self,
    ):
        return cbmpy.doFBA(self.m_model, quiet=True)

    def set_constraints(
        self, n: int, initial_biomasses, dt: float, times: list[str]
    ):
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
        for i in range(1, n):
            rids = self.m_model.getReactionIds(times[i])
            for rid in rids:
                reaction: Reaction = self.m_model.getReaction(rid)
                mid = self.m_model.identify_model_from_reaction(rid)

                # If the reaction is exchange, a linking layer reaction skip
                # (BM_) is for the biomass linking reactions
                if reaction.is_exchange or mid == "" or rid.startswith("BM_"):
                    continue

                # Bm of model id
                r_x_t = f"BM_{mid}_{times[i-1]}_{times[i]}"

                # r_x_t = biomass_rid + times[i - 1]
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

    def set_initial_concentrations(
        self,
        initial_biomasses: dict[str, float],
        initial_concentrations: dict[str, float],
        times,
    ):
        for key, value in initial_concentrations.items():
            # get species and it's corresponding exchange reaction
            species: Species = self.m_model.getSpecies(key + "_" + times[0])
            rids = species.getReagentOf()
            for rid in rids:
                reaction: Reaction = self.m_model.getReaction(rid)
                if reaction.is_exchange:
                    reaction.setLowerBound(-value)

        for key, value in initial_biomasses.items():
            self.m_model.setReactionBounds(
                f"BM_{key}_exchange", -value, -value
            )

    def remove_n_time_points(self, n):
        pass
