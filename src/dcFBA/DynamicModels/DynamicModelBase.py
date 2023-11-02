from ..Models.KineticsStruct import KineticsStruct


class DynamicModelBase:
    """Base class for the Dynamic FBA (Flux Balance Analysis) models.

    Attributes:
        _biomasses (dict): Dictionary containing biomass concentrations.
        _metabolites (dict): Dictionary containing metabolite concentrations.
        _fluxes (dict): Dictionary containing flux values.
        _times (list): List of time points used in the simulation.
        _kinetics (KineticsStruct): Structure holding kinetic information.
    """

    def __init__(self) -> None:
        """Initialize the DynamicModelBase class with default values."""
        self._biomasses = {}
        self._metabolites = {}
        self._fluxes = {}
        self._times = []
        self._kinetics: KineticsStruct = None

    def simulate(self):
        pass

    @property
    def fluxes(self):
        """Get the fluxes dictionary."""
        return self._fluxes

    @property
    def metabolites(self):
        """Get the metabolites dictionary."""
        return self._metabolites

    @property
    def biomasses(self):
        """Get the biomasses dictionary."""
        return self._biomasses

    @property
    def kinetics(self) -> KineticsStruct:
        """Get the kinetic information."""
        return self._kinetics

    @property
    def times(self):
        """Get the list of time-points used"""
        return self._times

    def get_biomasses(self):
        """Return the biomasses dictionary.

        Returns:
            dict: Dictionary containing biomass information.
        """
        return self.biomasses

    def get_metabolites(self):
        """Return the metabolites dictionary.

        Returns:
            dict: Dictionary containing metabolite information.
        """

        return self.metabolites

    def get_fluxes(self) -> list[float]:
        """Return the fluxes dictionary.

        Returns:
            dict: Dictionary containing flux information.
        """
        return self.fluxes

    def get_time_points(self) -> list:
        """Return the list of time points for the simulation.

        Returns:
            list: List of time points.
        """
        return self.times

    def get_flux_values(self) -> list[float]:
        pass

    def get_fluxes_values(self) -> dict[str, float]:
        """Get the flux values for each time point given a reaction ID.

        Returns:
            dict[str, float]: Dictionary of flux values for each time point.
        """
        pass

    def get_specific_flux_values(self) -> list[float]:
        """Get specific flux values.

        Returns:
            list[float]: Specific flux values.
        """

        pass

    def get_community_growth_rate(self) -> list[float]:
        """Get the community growth rate over time.

        Returns:
            list[float]: Community growth rate values over time.
        """
        pass

    def get_relative_abundance(self) -> list[float]:
        """Get the percentage of each species for each time point.

        Returns:
            list[float]: list containing the relative abundance of each species for each time point.
        """
        pass
