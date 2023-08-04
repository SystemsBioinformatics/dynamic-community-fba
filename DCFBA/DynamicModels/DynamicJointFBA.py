import cbmpy
import math
from .DynamicFBABase import DynamicFBABase
from cbmpy.CBModel import Reaction
from ..Models import KineticsStruct, Transporters, CommunityModel


class DynamicJointFBA(DynamicFBABase):
    """A class representing a dynamic joint Flux Balance Analysis simulation.

    This class extends the TimeStepDynamicFBABase and provides functionality
    for performing dynamic joint FBA simulations on a community model.

    """

    m_model: CommunityModel
    m_transporters: Transporters
    m_initial_bounds: dict[str, tuple[float, float]] = {}

    def __init__(
        self,
        model: CommunityModel,
        biomasses: list[float],
        initial_concentrations: dict[str, float] = {},
        kinetics: KineticsStruct = KineticsStruct({}),
    ):
        super().__init__(model, biomasses, initial_concentrations, kinetics)

        # Set X_C to be exporter since it increases over time
        self.set_community_biomass_reaction()
        self.m_transporters.add_exporter("X_comm", ["X_c"])
        self.m_metabolite_concentrations["X_c"] = [sum(biomasses)]

    def get_joint_model(self) -> CommunityModel:
        """Get the community model appended with the Community
        Biomass function.

        Returns:
            CommunityModel: The community model used for the simulation.

        """
        return self.m_model

    def set_community_biomass_reaction(self):
        """Create and set up the community biomass reaction.

        This method creates the community biomass species 'X_c' which
        will be created in each model's biomass reaction. Next the
        community biomass exchange reaction (X_comm) is created and set to be
        the objective function of the model.
        """
        self.m_model.createSpecies(
            "X_c", False, "The community biomass", compartment="e"
        )

        for _, biomass_id in self.m_model.get_model_biomass_ids().items():
            reaction: Reaction = self.m_model.getReaction(biomass_id)
            reaction.createReagent("X_c", 1)

        self.m_model.createReaction("X_comm")
        out: Reaction = self.m_model.getReaction("X_comm")
        out.is_exchange = True
        out.setUpperBound(1000)
        out.setLowerBound(0)
        out.createReagent("X_c", -1)

        self.m_model.createObjectiveFunction("X_comm")

        self.m_model.setActiveObjective("X_comm_objective")
