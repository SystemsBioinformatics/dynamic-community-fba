class Kinetics:
    """
    Kinetics

    This class represents a collection of reaction kinetics data. It stores
    the kinetics information as a dictionary, where the reaction identifier
    (rid) is the key, and the associated tuple of vmax and km values is the
    value.

    Kinetics reaction_id => tuple(limiting_species, km. vmax)

    Usage:
    - Add kinetics data for a reaction using the Add() method.
    - Remove kinetics data for a reaction using the Remove() method.
    - Check if kinetics data exists for a reaction using the Exists() method.
    - Retrieve km value for a reaction using the get_km() method.
    - Retrieve vmax value for a reaction using the get_vmax() method.
    - Retrieve the tuple of vmax and km values for a reaction using the
    get_kinetics() method.
    """

    # Kinetics tuple(limiting_species, km. vmax)
    Kinetics: dict[str, tuple[str, float, float]] = {}

    def __init__(
        self, kinetics: dict[str, tuple[str, float, float]] = {}
    ) -> None:
        self.Kinetics = kinetics

    def Add(
        self, rid: str, km: float, vmax: float, limiting_substrate: str
    ) -> None:
        self.Kinetics[rid] = [limiting_substrate, km, vmax]

    def Remove(self, rid: str) -> None:
        del self.Kinetics[rid]

    def Exists(self, rid: str) -> bool:
        if rid in self.Kinetics.keys():
            return True

    def get_km(self, rid: str) -> float:
        return self.Kinetics[rid][1]

    def get_vmax(self, rid) -> float:
        return self.Kinetics[rid][2]

    def get_kinetics(self, rid) -> tuple[float, float]:
        return self.Kinetics[rid]

    def get_limiting_substrate(self, rid) -> str:
        return self.Kinetics[rid][0]

    # def has_limiting_substrate(self, rid: str) -> bool:
    #     return self.Kinetics[rid][0] != ""
