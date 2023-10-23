from ..Models.Kinetics import KineticsStruct


class DynamicModelBase:
    """Base class for the Dynamic FBA classes"""

    def __init__(self) -> None:
        self._biomasses = {}
        self._metabolites = {}
        self._fluxes = {}
        self._times = []
        self._kinetics: KineticsStruct = None

    def simulate(self):
        pass

    @property
    def fluxes(self):
        return self._fluxes

    @property
    def metabolites(self):
        return self._metabolites

    @property
    def biomasses(self):
        return self._biomasses

    @property
    def kinetics(self) -> KineticsStruct:
        return self._kinetics

    @property
    def times(self):
        return self._times

    def get_biomasses(self):
        return self.biomasses

    def get_metabolites(self):
        return self.metabolites

    def get_fluxes(self):
        return self.fluxes

    def get_time_points(self):
        return self.times

    def get_flux_values(self) -> list[float]:
        """Get the flux values for each time point
            given a reaction ID

        Args:
            rid (str): _description_
        """
        pass

    def get_fluxes_values(self) -> dict[str, float]:
        """Get the flux values for each time point
            given a reaction ID

        Args:
            rid (str): _description_
        """
        pass

    def get_specific_flux_values(self) -> list[float]:
        pass

    def get_community_growth_rate(self):
        """Get the community growth rate over time"""
        pass

    def get_relative_abundance(self):
        """Get the percentage of each species for each time point"""
        pass
