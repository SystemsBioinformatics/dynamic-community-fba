from cbmpy.CBModel import Model

"""
 Here we define a kinetic model which is simply stated just a cbmpy model with
 additional kinetic information.
"""


class KineticModel(Model):
    m_model: Model
    m_kinetics: dict[str, tuple[float, float]]

    def __init__(
        self, model: Model, kinetics: dict[str, tuple[float, float]]
    ) -> None:
        self.m_model = model
        self.m_kinetics = kinetics

    def get_model_kinetics(self) -> dict[str, tuple[float, float]]:
        return self.m_kinetics

    def get_model(self) -> Model:
        return self.m_model

    def get_reaction_kinetics(self, rid) -> tuple[float, float]:
        try:
            ans = self.m_kinetics[rid]
        # No such reaction
        except KeyError:
            ans = None
        return ans

    def get_reaction_vmax(self, rid) -> float:
        kinetics = self.get_reaction_kinetics(rid)
        if kinetics is not None:
            return kinetics[0]
        return None

    def get_reaction_km(self, rid) -> float:
        kinetics = self.get_reaction_kinetics(rid)
        if kinetics is not None:
            return kinetics[1]
        return None

    def set_model(self, model: Model) -> None:
        self.m_model = model

    def set_kinetics(self, kinetics: dict[str, tuple[float, float]]):
        self.m_kinetics = kinetics

    def set_reaction_kinetics(self, rid: str, kinetics: tuple[float, float]):
        self.m_kinetics[rid] = kinetics
