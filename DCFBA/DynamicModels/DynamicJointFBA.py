import cbmpy
from .DynamicFBABase import DynamicFBABase
from cbmpy.CBModel import Reaction
from ..Models import KineticsStruct, CommunityModel


class DynamicJointFBA(DynamicFBABase):
    """
    A class to perform dynamic joint Flux Balance Analysis (FBA) on a community model.

    This class facilitates dynamic joint FBA simulations by providing methods to set up
    and manage the community model, including the creation and management of community biomass reactions.

    Attributes:
        m_model (CommunityModel): The community metabolic model being simulated.
        m_initial_bounds (dict[str, tuple[float, float]], optional): Initial bounds for the model reactions.
            Defaults to an empty dictionary.
    """

    def __init__(
        self,
        model: CommunityModel,
        biomasses: list[float],
        initial_concentrations: dict[str, float] = {},
        kinetics: KineticsStruct = KineticsStruct({}),
    ):
        """
        Initialize the DynamicJointFBA class.

        Args:
            model (CommunityModel): The community metabolic model to simulate.
            biomasses (list[float]): List of initial biomass concentrations for each model in the community.
            initial_concentrations (dict[str, float], optional): Initial concentrations for metabolites.
                Defaults to an empty dictionary.
            kinetics (KineticsStruct, optional): Kinetic parameters for reactions in the model.
                Defaults to an empty KineticsStruct.
        """

        super().__init__(model, biomasses, initial_concentrations, kinetics)

        # Set X_C to be exporter since it increases over time
        self.set_community_biomass_reaction()
        self._metabolites["X_c"] = [sum(biomasses)]

    def get_joint_model(self) -> CommunityModel:
        """
        Retrieve the community model, which includes the Community Biomass function.

        Returns:
            CommunityModel: The underlying community model used for simulation.
        """

        return self.m_model

    def set_community_biomass_reaction(self) -> None:
        """
        Set up the community biomass reaction.

        This method establishes the community biomass species 'X_c' and associates it with
        each model's biomass reaction. Additionally, a community biomass exchange reaction (X_comm)
        is created and designated as the model's objective function.
        """
        self.model.createSpecies(
            "X_c", False, "The community biomass", compartment="e"
        )

        for _, biomass_id in self.model.get_model_biomass_ids().items():
            reaction: Reaction = self.model.getReaction(biomass_id)
            reaction.createReagent("X_c", 1)

        self.model.createReaction("X_comm")
        out: Reaction = self.model.getReaction("X_comm")
        out.is_exchange = True
        out.setUpperBound(cbmpy.INF)
        out.setLowerBound(0)
        out.createReagent("X_c", -1)

        self.model.createObjectiveFunction("X_comm")

        self.model.setActiveObjective("X_comm_objective")
