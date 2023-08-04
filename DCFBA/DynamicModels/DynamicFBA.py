import TimeStepDynamicFBABase
import DynamicJointFBA
from ..Models import KineticsStruct
from ..Models import CommunityModel

x


class DynamicFBA(TimeStepDynamicFBABase):
    def __init__(
        self,
        model,
        biomass_id: str,
        biomasses: list[float],
        initial_concentrations: dict[str, float] = {},
        kinetics: KineticsStruct = KineticsStruct({}),
    ):
        cm = CommunityModel([model], [biomass_id], [""], model.getId())
        dj = DynamicJointFBA()

    def simulate(self, dt):
        pass
