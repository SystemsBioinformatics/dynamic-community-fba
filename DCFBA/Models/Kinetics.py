class KineticsStruct:
    """
    This class represents a collection of reaction kinetics data. It stores
    the kinetics information as a dictionary, where the reaction identifier
    (rid) is the key, and the associated tuple of vmax and km values is the
    value.

    kinetics reaction_id => tuple(limiting_species_id, km, vmax)

    Usage:
        - Add kinetics data for a reaction using the Add() method.
        - Remove kinetics data for a reaction using the Remove() method.
        - Check if kinetics data exists for a reaction using the Exists() method.
        - Retrieve km value for a reaction using the get_km() method.
        - Retrieve vmax value for a reaction using the get_vmax() method.
        - Retrieve the tuple of vmax and km values for a reaction using the get_kinetics() method.
    """

    # kinetics rid: tuple(limiting_species, km. vmax)
    kinetics: dict[str, tuple[str, float, float]] = {}

    def __init__(
        self, kinetics: dict[str, tuple[str, float, float]] = {}
    ) -> None:
        self.kinetics = kinetics

    def add(
        self, rid: str, km: float, vmax: float, limiting_substrate: str
    ) -> None:
        self.kinetics[rid] = [limiting_substrate, km, vmax]

    def remove(self, rid: str) -> None:
        del self.kinetics[rid]

    def exists(self, rid: str) -> bool:
        """Checks if kinetic data is available for reaction

        Args:
            rid (str): _description_

        Returns:
            bool: _description_
        """
        if rid in self.kinetics.keys():
            return True

    def get_km(self, rid: str) -> float:
        return self.kinetics[rid][1]

    def get_vmax(self, rid) -> float:
        return self.kinetics[rid][2]

    def get_reactions_kinetics(self, rid) -> tuple[float, float]:
        return self.kinetics[rid]

    def get_limiting_substrate(self, rid) -> str:
        return self.kinetics[rid][0]

    def get_kinetics(self) -> dict[str, tuple[str, float, float]]:
        return self.kinetics

    # def has_limiting_substrate(self, rid: str) -> bool:
    #     return self.kinetics[rid][0] != ""
