from cbmpy.CBModel import Model, Reaction
from ..Models.KineticsStruct import KineticsStruct
from ..Models.Transporters import Transporters
from .DynamicModelBase import DynamicModelBase
from ..Exceptions import NoLimitingSubstrateFound


class StaticOptimizationModelBase(DynamicModelBase):

    """Base class providing a framework for the static optimization approaches
    (SOA)
    Base class for methods using time steps to track biomass and metabolite
    concentrations over time.

    Attributes:
        _initial_bounds (dict): Containing the original upper and lower bound
            for all reactions

    """

    def __init__(self) -> None:
        super().__init__()
        self._initial_bounds = {}

    def _set_initial_concentrations(
        self, model: Model, initial_concentrations: dict[str, float]
    ) -> None:
        """Sets the initial concentrations using either the provided
        concentrations or the lower bounds of the model's exchange reactions.

        Args:
            model (Model): A cbmpy.CBModel instance.
            initial_concentrations (dict[str, float]): Dictionary containing
                the initial concentrations as key-value pairs for different
                metabolites.

        """

        for exchange in model.getExchangeReactionIds():
            reaction: Reaction = model.getReaction(exchange)
            species_id: str = reaction.reagents[0].getSpecies()

            if species_id in initial_concentrations.keys():
                self.metabolites[species_id] = [
                    initial_concentrations[species_id]
                ]
            else:
                # Max such that positive lower bounds in exchanges are not set
                # code to make compatible with ParallelFBA if no species concentrairon
                # is defined in the initial_concentrations than we keep track of the lowest
                # exchange reaction of all models provided
                # TODO check this
                concentration = (
                    0
                    if species_id not in self.metabolites.keys()
                    else self.metabolites[species_id][-1]
                )
                self.metabolites[species_id] = [
                    max(-reaction.getLowerBound(), concentration)
                ]

        # This is in place for bad GSMM's which do not set an exchange for each
        # extracellular metabolite
        for species in model.species:
            if (
                species.getCompartmentId() == "e"
                and species.getId() not in self.metabolites.keys()
            ):
                self.metabolites[species.getId()] = [0]

    @property
    def initial_bounds(self):
        return self._initial_bounds

    def simulate(
        self,
        dt: float,
        n: int = 10000,
        epsilon=0.001,
        kinetics_func=None,
        deviate=None,
    ):
        """Placeholder for dynamic FBA simulation over specified time intervals.
        This method should be overridden in subclasses to implement specific simulation logic.

        Args:
            dt (float): The time step for simulation.
            n (int, optional): number of simulation
                Defaults to 10000
            epsilon (float, optional): The tolerance value. When the solution
                or time step is less than epsilon, the simulation stops.
                Defaults to 0.001.
            kinetics_func (function, optional): A user-defined function to
                calculate kinetics. Defaults to None.
            deviate (function, optional): A function to apply model changes during the simulation.
                Should accept: the model, biomass concentrations, metabolite concentrations, and dt as parameters.
                Defaults to None.
        """
        pass

    def update_reaction_bounds(self, kinetics_func) -> None:
        """Updates all reaction bounds using the new concentrations.

        Args:
            kinetics_func (function): A user-defined function that calculates
                the new upper bound for a reaction. It should accept:
                - cbmpy.CBModel, a string (reaction_id),
                - dictionary of biomass concentrations,
                - a dictionary of metabolite concentrations.
        """
        pass

    def update_concentrations(self) -> None:
        """Update all metabolite concentrations"""
        pass

    def update_exchanges(self, dt: float) -> None:
        """Updates the lower and upper bound of the
        exchange reaction, according to the amount of external metabolite
        present

        Args:
            dt (float): time step size
        """
        pass

    def update_biomasses(self) -> None:
        """Update all biomass concentration"""
        pass

    def check_solution_feasibility(self) -> str:
        """Checks if the current solution has any metabolite concentration
        below zero.

        Returns:
            str: The species ID with the lowest negative concentration,
                or an empty string if all concentrations are positive.
        """

        low = 1e10
        sid = ""
        for key, value in self.metabolites.items():
            if round(value[-1], 8) < 0 and value[-1] < low:
                low = value[-1]
                sid = key

        return sid

    def importers_species_concentration(
        self, rid: str, transporters: Transporters
    ) -> list[float]:
        """Given an importer reaction, returns the species concentrations that
        this importer imports.

        Args:
            rid (str): The reaction id of the importer.
            transporters (Transporters): A Transporters object.

        Returns:
            list[float]: List of species concentrations.
        """
        sids = transporters.get_importers_species(rid)
        return [self.metabolites[id][-1] for id in sids]

    def mm_kinetics(
        self,
        reaction: Reaction,
        X: float,
        kinetics: KineticsStruct,
    ) -> None:
        """Computes the Michaelis-Menten kinetics for a given reaction.

         If the reaction is an importer and no limiting substrate was provided,
         the function uses the minimum concentration of all substrates.


         Args:
             reaction (Reaction): A reaction object.
             X (float): The biomass of the species the reaction belongs to.
             transporters (Transporters): A Transporters object of the model.
             kinetics (Kinetics): A Kinetics object for kinetic parameter
             retrieval.

        Raises:
             NoLimitingSubstrateFound: Raised if the reaction isn't an import reaction and
             no limiting substrate was provided or if the provided substrate isn't external.
        """

        rid: str = reaction.getId()
        sid, km, vmax = kinetics.get_reactions_kinetics(
            rid,
        )

        if sid in self.metabolites.keys():
            S = self.metabolites[sid][-1]
        else:
            raise NoLimitingSubstrateFound(
                "The limiting substrate was not an external species"
            )

        reaction.setUpperBound(vmax * (S / (km + S)) * X)
