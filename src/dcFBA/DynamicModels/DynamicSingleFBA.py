from .DynamicFBABase import DynamicFBABase
from ..Models import CommunityModel, KineticsStruct
from cbmpy.CBModel import Model


class DynamicSingleFBA(DynamicFBABase):
    """
    DynamicSingleFBA Class.

    This class provides functionality to perform a dynamic Flux Balance Analysis (FBA)
    on a single model. It inherits from the DynamicFBABase and internally converts
    the single model to a community model representation for consistent handling of
    dynamic simulations.
    """

    def __init__(
        self,
        model: Model,
        biomass_id: str,
        initial_biomass: float,
        initial_concentrations: dict[str, float] = {},
        kinetics: KineticsStruct = KineticsStruct({}),
    ):
        """Initializes the DynamicSingleFBA class.

        Args:
            model (Model): The metabolic model you want to simulate.
            biomass_id (str): Identifier for the biomass reaction in the model.
            initial_biomass (float): Initial biomass concentration for the model.
            initial_concentrations (dict[str, float], optional): Initial concentrations
                for metabolites. Defaults to an empty dictionary.
            kinetics (KineticsStruct, optional): Kinetic parameters for reactions in
                the model. Defaults to an empty KineticsStruct.
        """
        cm = CommunityModel([model], [biomass_id], [""], model.getId())

        if biomass_id not in cm.getReactionIds():
            raise Exception(f"Biomass id: {biomass_id} not a reaction")
        cm.createObjectiveFunction(biomass_id)

        cm.setActiveObjective(f"{biomass_id}_objective")

        super().__init__(
            cm, [initial_biomass], initial_concentrations, kinetics
        )

    def get_biomass(self) -> list[float]:
        """Get the total biomass over time

        Returns:
            list[float]: Biomass concentration over time
        """
        biomasses = super().get_biomasses()
        return biomasses[""]
