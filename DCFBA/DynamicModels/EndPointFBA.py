import cbmpy
import re
from numpy import Inf
from cbmpy.CBModel import Reaction, Species
from ..Models.CommunityModel import CommunityModel
from ..Helpers.BuildEndPointModel import build_time_model
from .DynamicFBABase import DynamicFBABase
from ..Models import KineticsStruct


class EndPointFBA(DynamicFBABase):
    m_model = CommunityModel
    m_times: list[str]
    m_kinetics: KineticsStruct

    # m_original_reactions: list[str]

    def __init__(
        self,
        community_model: CommunityModel,
        n: int,
        initial_biomasses: dict[str, float],
        initial_concentrations: dict[str, float] = {},
        kinetics: KineticsStruct = None,
        dt: float = 0.1,
    ) -> None:
        width = len(str(n))
        self.m_times = [f"time{i:0{width}d}" for i in range(n)]
        self.m_model = build_time_model(community_model, self.m_times)
        self.m_kinetics = kinetics
        self.set_objective()
        self.set_constraints(n, initial_biomasses, dt)
        self.set_initial_concentrations(
            initial_biomasses, initial_concentrations
        )

    def set_objective(self):
        self.m_model.createReaction("X_comm", silent=True)
        out: Reaction = self.m_model.getReaction("X_comm")
        out.is_exchange = True
        out.setUpperBound(Inf)
        out.setLowerBound(0)
        out.createReagent("BM_c_" + self.m_times[-1], -1)

        self.m_model.createObjectiveFunction("X_comm")

        self.m_model.setActiveObjective("X_comm_objective")

    def simulate(
        self,
    ):
        return cbmpy.doFBA(self.m_model, quiet=True)

    def set_constraints(self, n: int, initial_biomasses, dt: float):
        rids_t0 = self.m_model.getReactionIds(self.m_times[0])
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
            rids = self.m_model.getReactionIds(self.m_times[i])
            for rid in rids:
                reaction: Reaction = self.m_model.getReaction(rid)
                mid = self.m_model.identify_model_from_reaction(rid)

                # If the reaction is exchange, a linking layer reaction skip
                # (BM_) is for the biomass linking reactions
                if reaction.is_exchange or mid == "" or rid.startswith("BM_"):
                    continue

                # Bm of model id
                r_x_t = f"BM_{mid}_{self.m_times[i-1]}_{self.m_times[i]}"

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
    ):
        for key, value in initial_concentrations.items():
            # get species and it's corresponding exchange reaction
            species: Species = self.m_model.getSpecies(
                key + "_" + self.m_times[0]
            )
            rids = species.getReagentOf()
            for rid in rids:
                reaction: Reaction = self.m_model.getReaction(rid)
                if reaction.is_exchange:
                    reaction.setLowerBound(-value)

        for key, value in initial_biomasses.items():
            self.m_model.setReactionBounds(
                f"BM_{key}_exchange", -value, -value
            )

    def constrain_rates(self, epsilon=0.1):
        old_rids = set(
            [id.split("_time")[0] for id in self.m_model.getReactionIds()]
        )
        for rid in old_rids:
            reaction: Reaction = self.m_model.getReaction(
                rid + "_" + self.m_times[0]
            )

            if reaction is not None and (not reaction.is_exchange):
                for i, time in enumerate(self.m_times[:-1]):
                    id_t0 = f"{rid}_{time}"
                    id_t1 = f"{rid}_{self.m_times[i+1]}"
                    self.m_model.addUserConstraint(
                        "R_constraint_pos" + id_t0,
                        [[1, id_t0], [-1, id_t1]],
                        "<=",
                        epsilon,
                    )

                    self.m_model.addUserConstraint(
                        "R_constraint_neg" + id_t0,
                        [[1, id_t0], [-1, id_t1]],
                        ">=",
                        -epsilon,
                    )

    def set_kinetics(self, kin: KineticsStruct):
        self.m_kinetics = kin

    def mm_approximation(self, dt):
        if self.m_kinetics is not None:
            for rid, _ in self.m_kinetics.get_kinetics().items():
                sid, km, vmax = self.m_kinetics.get_reactions_kinetics(rid)
                if sid == "":
                    raise Exception(
                        "No limiting substrate was set in the kinetics object for: "
                        + rid
                    )
                low_line = vmax / km
                high_line = (vmax / 2) / km

                for i in range(0, len(self.m_times) - 1):
                    # Linking reaction is the concentration of S for the
                    # timepoint
                    linking_reaction_id = (
                        f"{sid}_{self.m_times[i]}_{self.m_times[i+1]}"
                    )
                    t_rid = rid + "_" + self.m_times[i + 1]
                    print(t_rid)
                    self.m_model.addUserConstraint(
                        f"mm_low_{linking_reaction_id}",
                        [[1, t_rid], [-low_line * dt, linking_reaction_id]],
                        ">=",
                        0.0,
                    )

                    self.m_model.addUserConstraint(
                        f"mm_high_{linking_reaction_id}",
                        [[1, t_rid], [-high_line * dt, linking_reaction_id]],
                        "<=",
                        0.0,
                    )

    def set_qp(self, solution, epsilon=0.01):
        obj = self.m_model.getActiveObjective()
        obj.setOperation("minimize")
        obj.QPObjective = []

        for i, tid in enumerate(self.m_times[:-1]):
            rids = self.m_model.getReactionIds(tid)
            for rid in rids:
                reaction = self.m_model.getReaction(rid)

                # R_ are all the original reqactions
                # We dont want the exchange reactions
                if rid.startswith("R_") and not reaction.is_exchange:
                    rid = re.match(r"(.*?)_(time\d*)", rid).group(1)
                    rid_t_t1 = f"{rid}_{self.m_times[i+1]}"
                    rid_t_t0 = f"{rid}_{self.m_times[i]}"

                    if i == 0:
                        obj.QPObjective.append(((rid_t_t0, rid_t_t0), 1.0))
                    else:
                        obj.QPObjective.append(((rid_t_t0, rid_t_t0), 2.0))

                    obj.QPObjective.append(((rid_t_t0, rid_t_t1), -2.0))

                    obj.QPObjective.append(((rid_t_t1, rid_t_t1), 1.0))

        self.m_model.getReaction("X_comm").setLowerBound(solution - epsilon)
        self.m_model.getReaction("X_comm").setUpperBound(solution + epsilon)
