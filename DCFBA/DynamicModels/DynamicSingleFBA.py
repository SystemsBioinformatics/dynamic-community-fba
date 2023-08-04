from .DynamicFBABase import DynamicFBABase
from ..Models import CommunityModel, KineticsStruct
from cbmpy.CBModel import Model


class DynamicSingleFBA(DynamicFBABase):
    def __init__(
        self,
        model: Model,
        biomass_id: str,
        initial_biomass: float,
        initial_concentrations: dict[str, float] = {},
        kinetics: KineticsStruct = KineticsStruct({}),
    ):
        cm = CommunityModel([model], [biomass_id], [""], model.getId())

        cm.createObjectiveFunction(biomass_id)

        super().__init__(
            cm, [initial_biomass], initial_concentrations, kinetics
        )
