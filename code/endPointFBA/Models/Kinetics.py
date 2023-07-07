class Kinetics:
    Kinetics: dict[str, tuple[float, float]] = {}

    def __init__(self, kinetics=dict[str, tuple[float, float]]) -> None:
        self.Kinetics = kinetics

    def Add(self, rid: str, vmax: float, km: float) -> None:
        self.Kinetics[rid] = [km, vmax]

    def Remove(self, rid: str) -> None:
        del self.Kinetics[rid]

    def Exists(self, rid: str) -> bool:
        if rid in self.Kinetics.keys():
            return True

    def get_km(self, rid: str) -> float:
        return self.Kinetics[rid][0]

    def get_vmax(self, rid) -> float:
        return self.Kinetics[rid][1]

    def get_kinetics(self, rid) -> tuple[float, float]:
        return self.Kinetics[rid]
