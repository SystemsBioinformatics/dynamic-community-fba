from cbmpy.CBModel import Model, Reaction
from ..Models.Kinetics import KineticsStruct
from ..Models.Transporters import Transporters
from .DynamicModelBase import DynamicModelBase
from ..Exceptions import NoLimitingSubstrateFound


class StaticOptimizationModelBase(DynamicModelBase):

    """Base class providing a framework for the static optimization approaches (SOA)
    methods; using time steps to track biomass and metabolite concentrations over time.
    """

    m_biomass_concentrations: dict[str, list[float]] = {}
    m_metabolite_concentrations: dict[str, list[float]] = {}
    m_initial_bounds = {}
    m_kinetics: KineticsStruct

    def set_initial_concentrations(
        self, model: Model, initial_concentrations: dict[str, float]
    ) -> None:
        """Sets initial concentrations of metabolites based on provided values
        or reaction bounds.
        This method sets the initial concentrations using either the provided
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
                self.m_metabolite_concentrations[species_id] = [
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
                    if species_id
                    not in self.m_metabolite_concentrations.keys()
                    else self.m_metabolite_concentrations[species_id][-1]
                )
                self.m_metabolite_concentrations[species_id] = [
                    max(-reaction.getLowerBound(), concentration)
                ]

        # This is in place for bad GSMM's which do not set an exchange for each
        # extracellular metabolite
        for species in model.species:
            if (
                species.getCompartmentId() == "e"
                and species.getId()
                not in self.m_metabolite_concentrations.keys()
            ):
                self.m_metabolite_concentrations[species.getId()] = [0]

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
            deviation_time (int, optional): The time step when the deviation
                function should be called.
                Defaults to 0.

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
        pass

    def update_biomasses(self) -> None:
        """Update all biomass concentration"""
        pass

    def reset_dt(self) -> float:
        """Resets the simulation and calculates dt_hat (smaller time step)."""
        pass

    def check_solution_feasibility(self) -> str:
        """Checks if the current solution has any metabolite concentration below zero.

        Returns:
            str: The metabolite ID with the lowest negative concentration, or an empty string if all concentrations are positive.
        """

        low = 1e10
        name = ""
        for key, value in self.m_metabolite_concentrations.items():
            if value[-1] < 0 and value[-1] < low:
                low = value[-1]
                name = key

        return name

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
        return [self.m_metabolite_concentrations[id][-1] for id in sids]

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

        if sid in self.m_metabolite_concentrations.keys():
            S = self.m_metabolite_concentrations[sid][-1]
        else:
            raise NoLimitingSubstrateFound(
                "The limiting substrate was not an external species"
            )

        reaction.setUpperBound(vmax * (S / (km + S)) * X)
