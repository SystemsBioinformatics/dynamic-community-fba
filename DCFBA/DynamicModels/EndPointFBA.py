import cbmpy
import numpy
from cbmpy.CBModel import Reaction
from ..Models.CommunityModel import CommunityModel
from ..Helpers.BuildEndPointModel import build_time_model
from .DynamicFBABase import DynamicFBABase


class EndPointFBA(DynamicFBABase):
    m_model = CommunityModel

    def __init__(
        self,
        community_model: CommunityModel,
        n: int,
        initial_biomasses: dict[str, float],
        dt: float = 0.1,
    ) -> None:
        width = len(str(n))
        times = [f"_time{i:0{width}d}" for i in range(n)]

        self.add_biomass_species(community_model)
        self.m_model = build_time_model(community_model, times)

        self.set_constraints(n, initial_biomasses, dt, times)
        self.set_objective(times[-1])

    def simulate(
        self,
    ):
        return cbmpy.doFBA(self.m_model)

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
                reaction.setUpperBound(numpy.inf)

    def add_biomass_species(self, model: CommunityModel):
        model.createSpecies(
            "X_c", False, "The community biomass", compartment="e"
        )

        for _, biomass_id in model.get_model_biomass_ids().items():
            reaction: Reaction = model.getReaction(biomass_id)
            reaction.createReagent("X_c", 1)

    def set_objective(self, time_id):
        self.m_model.createReaction("X_comm")
        out: Reaction = self.m_model.getReaction("X_comm")
        out.is_exchange = True
        out.setUpperBound(1e10)
        out.setLowerBound(0)
        out.createReagent("X_c" + time_id, -1)

        self.m_model.createObjectiveFunction("X_comm")

        self.m_model.setActiveObjective("X_comm_objective")
