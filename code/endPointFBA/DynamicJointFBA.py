import cbmpy
import math
from cbmpy.CBModel import Model, Reaction
from endPointFBA.CommunityModel import CommunityModel
from endPointFBA.DynamicFBABase import DynamicFBABase
from endPointFBA.Models.Kinetics import Kinetics
from endPointFBA.Models.Transporters import Transporters


class DynamicJointFBA(DynamicFBABase):
    m_model: CommunityModel
    m_transporters: Transporters
    m_initial_bounds: dict[str, tuple[float, float]] = {}

    # use to set everything in place for the simulation
    def __init__(
        self,
        model: CommunityModel,
        biomasses: list[float],
        initial_concentrations: dict[str, float] = {},
        kinetics: Kinetics = Kinetics({}),
    ):
        self.m_transporters = Transporters(model)

        model_biomasses = model.get_model_biomass_ids()

        self.set_community_biomass_reaction(model)

        initial_biomasses = [[x] for x in biomasses]

        self.m_biomass_concentrations = dict(
            zip(model_biomasses.keys(), initial_biomasses)
        )

        self.m_kinetics = kinetics
        self.set_initial_concentrations(model, initial_concentrations)

        for rid in self.m_model.getReactionIds():
            reaction: Reaction = self.m_model.getReaction(rid)
            self.m_initial_bounds[rid] = [
                reaction.getLowerBound(),
                reaction.getUpperBound(),
            ]

    def set_community_biomass_reaction(self, model: CommunityModel):
        # TODO implement clone funciton in CommunirtModel
        joint_model = model
        joint_model.createSpecies("X_c", False, "The community biomass")

        for _, biomass_id in model.get_model_biomass_ids().items():
            reaction: Reaction = joint_model.getReaction(biomass_id)
            reaction.createReagent("X_c", 1)

        joint_model.createReaction("X_comm")
        out: Reaction = joint_model.getReaction("X_comm")
        out.is_exchange = True
        out.setUpperBound(1000)
        out.setLowerBound(0)
        out.createReagent("X_c", -1)

        joint_model.createObjectiveFunction("X_comm")

        joint_model.setActiveObjective("X_comm_objective")

        # Set X_C to be exporter since it increases over time
        self.m_transporters.add_exporter("X_comm", ["X_c"])
        self.m_model = joint_model

        # we return the joint model such that you can also perform a regular
        # FBA

        return joint_model

    def get_joint_model(self) -> CommunityModel:
        return self.m_model

    # TODO add a continue function for time point t
    def simulate(
        self,
        dt: float,
        epsilon=0.001,
        user_func=None,
    ):
        # First update all upper and lower bounds to be rates
        # (Corrected for biomass),
        self.update_reaction_bounds(user_func)

        used_time = [0]
        dt_hat = -1
        dt_save = dt

        while True:
            solution = cbmpy.doFBA(self.m_model)

            if math.isnan(solution) or solution == 0 or solution < epsilon:
                break

            if dt_hat != -1:
                dt = dt_save - dt_hat
                dt_hat = -1
            else:
                dt = dt_save

            used_time.append(used_time[-1] + dt)

            FBAsol = self.m_model.getSolutionVector(names=True)
            FBAsol = dict(zip(FBAsol[1], FBAsol[0]))

            self.update_concentrations(FBAsol, dt)

            species_id = self.check_solution_feasibility()

            if species_id != "":
                # TODO This code needs some proper clean up
                # Same as for DynamicParallelFBA

                # This time point is not be feasible remove all last
                # concentrations

                self.m_metabolite_concentrations = {
                    key: lst[:-1]
                    for key, lst in self.m_metabolite_concentrations.items()
                }

                # Maybe some metabolite was exported/created, by an organism
                # Unfortunately to check what is create we would
                # have to know the new dt which we don'...
                # So we accept a small error, but in real life we
                # also would not know if the exporters would have
                # Created new metabolites before the importers ran out

                total = 0
                for (
                    rid,
                    species_ids,
                ) in self.m_transporters.get_importers(True):
                    if species_id in species_ids:
                        total += FBAsol[rid]
                dt_hat = (
                    self.m_metabolite_concentrations[species_id][-1] / total
                )

                # than recalculate
                self.update_concentrations(
                    FBAsol,
                    dt_hat,
                )

                used_time[-1] = used_time[-2] + dt_hat

            # update biomass
            for _, rid in self.m_model.get_model_biomass_ids().items():
                mid = self.m_model.identify_model_from_reaction(rid)
                Xt = self.m_biomass_concentrations[mid][-1] + FBAsol[rid] * dt
                self.m_biomass_concentrations[mid].append(Xt)

            # update lower and upper bounds...
            self.update_reaction_bounds(user_func)

        return [
            used_time,
            self.m_metabolite_concentrations,
            self.m_biomass_concentrations,
        ]

    def update_concentrations(self, FBAsol, dt):
        # Update external metabolites
        for key in self.m_metabolite_concentrations.keys():
            self.m_metabolite_concentrations[key].append(
                self.m_metabolite_concentrations[key][-1]
            )
        # first check what was exported (created) because new metabolite can
        # be created
        for rid, species_ids in self.m_transporters.get_exporters(True):
            for sid in species_ids:
                self.m_metabolite_concentrations[sid][-1] += FBAsol[rid] * dt

        # Next check the importers (consumption)
        for rid, species_ids in self.m_transporters.get_importers(True):
            for sid in species_ids:
                self.m_metabolite_concentrations[sid][-1] -= FBAsol[rid] * dt

    def update_reaction_bounds(self, user_func) -> None:
        if user_func is not None:
            # TODO Give somethings to the user_func, think about this
            user_func()
            return
        for rid in self.m_model.getReactionIds():
            reaction: Reaction = self.m_model.getReaction(rid)
            # Don't change the exchange reactions bounds
            if not reaction.is_exchange:
                mid_for_reaction = self.m_model.identify_model_from_reaction(
                    rid
                )
                # Organism specific biomass
                X_k_t = self.m_biomass_concentrations[mid_for_reaction][-1]
                reaction.setLowerBound(self.m_initial_bounds[rid][0] * X_k_t)

                # If the reaction is an importer we need to check
                # if there is substrate they can import
                if self.m_transporters.is_importer(rid):
                    self.update_importer_bounds(reaction, X_k_t)

                elif self.m_kinetics.Exists(rid):
                    self.update_kinetics(reaction, X_k_t, self.m_transporters)

                else:
                    reaction.setUpperBound(
                        self.m_initial_bounds[rid][1] * X_k_t
                    )

    def update_importer_bounds(self, reaction: Reaction, X: float):
        if all(
            x > 0
            for x in self.importers_species_concentration(
                reaction.getId(), self.m_transporters
            )
        ):
            if self.m_kinetics.Exists(reaction.getId()):
                self.update_kinetics(reaction, X, self.m_transporters)
            else:
                reaction.setUpperBound(
                    self.m_initial_bounds[reaction.getId()][1] * X
                )
        else:
            reaction.setUpperBound(0)
