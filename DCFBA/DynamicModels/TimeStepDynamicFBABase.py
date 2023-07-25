from cbmpy.CBModel import Model, Reaction
from ..Models.Kinetics import Kinetics
from ..Models.Transporters import Transporters
from .DynamicFBABase import DynamicFBABase
from ..Exceptions import NoLimitingSubstrateFound


class TimeStepDynamicFBABase:
    """Base class for handling dynamic flux balance analysis simulations
    using time steps to track biomass and metabolite concentrations over time.

    """

    m_biomass_concentrations: dict[str, list[float]] = {}
    m_metabolite_concentrations: dict[str, list[float]] = {}

    def set_initial_concentrations(
        self, model: Model, initial_concentrations: dict[str, float]
    ) -> None:
        """Sets the initial concentrations of the metabolites.

        The method uses the exchange reaction lower bounds as initial
        concentrations. If you want to change the initial concentrations,
        you can modify either the initial concentrations or the lower bounds of
        exchange reactions.

        Args:
            model (Model): A cbmpy.CBModel instance.
            initial_concentrations (dict[str, float]): Dictionary containing
                the initial concentrations as key-value pairs for different
                metabolites.
        """

        for exchange in model.getExchangeReactionIds():
            reaction: Reaction = model.getReaction(exchange)
            species_id: str = reaction.reagents[0].getSpecies()
            if species_id not in self.m_metabolite_concentrations.keys():
                if species_id in initial_concentrations.keys():
                    self.m_metabolite_concentrations[species_id] = [
                        initial_concentrations[species_id]
                    ]
                    # TODO set species exchange to this value
                else:
                    self.m_metabolite_concentrations[species_id] = [
                        -reaction.getLowerBound()
                    ]

    def simulate(
        self,
        dt: float,
        n: int = 10000,
        epsilon=0.001,
        kinetics_func=None,
        deviate=None,
    ):
        """Simulates the dynamic flux balance analysis (dFBA) over specified time steps.

        This method should be implemented in the subclasses with their own rules for simulation.

        Args:
            dt (float): The time step for simulation.
            n (int, optional): number of simulation
                Defaults to 10000
            epsilon (float, optional): The tolerance value. When the solution
                or time step is less than epsilon, the simulation stops.
                Defaults to 0.001.
            kinetics_func (function, optional): A user-defined function to
                calculate kinetics. Defaults to None.
            deviate (function, optional): A user-defined function to apply
                changes during the simulation. It should accept the model,
                a dictionary of biomass concentrations, a dictionary of
                metabolite concentrations, and dt.
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

    def update_importer_bounds(self) -> None:
        """Update the upper bounds of the importer reactions"""
        pass

    def update_concentrations(self) -> None:
        """Update all metabolite concentrations"""
        pass

    def update_biomasses(self) -> None:
        """Update all biomass concentration"""
        pass

    def reset_dt(self) -> float:
        """Resets the simulation and calculates dt_hat (smaller time step)."""
        pass

    def check_solution_feasibility(self) -> str:
        """Checks the feasibility of the current solution.

        If the current solution isn't feasible due to a lack of metabolite
        concentrations.

        Returns:
            str: The name of the metabolite with the lowest concentration
                below zero.
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
        transporters: Transporters,
        kinetics: Kinetics,
    ) -> None:
        """Calculates the Michaelis Menten Kinetics for a given reaction.

        Args:
            reaction (Reaction): A reaction object.
            X (float): The biomass of the species the reaction belongs to.
            transporters (Transporters): A Transporters object of the model.
            kinetics (Kinetics): A Kinetics object for kinetic parameter
            retrieval.
        Raises:
            NoLimitingSubstrateFound: If the reaction is not an importer and no
            limiting substrate was supplied, or if the limiting substrate was
            not an external species.
        """
        rid: str = reaction.getId()
        sid, km, vmax = kinetics.get_kinetics(
            rid,
        )

        # If the reaction is an importer and no limiting substrate
        # was supplied use the min of all substrates
        # this is only possible for importers
        # If the reaciton is not an import reaction and no limiting substrate
        # Was set Raise an error
        if transporters.is_importer(rid) and (
            sid not in self.m_metabolite_concentrations.keys()
        ):
            S = min(self.importers_species_concentration(rid, transporters))
        elif sid in self.m_metabolite_concentrations.keys():
            S = self.m_metabolite_concentrations[sid][-1]
        else:
            raise NoLimitingSubstrateFound(
                "The limiting substrate was not an external species"
            )
        reaction.setUpperBound(vmax * (S / (km + S)) * X)
