from cbmpy.CBModel import Model, Reaction
from ..Models.Kinetics import Kinetics
from ..Models.Transporters import Transporters
from .DynamicFBABase import DynamicFBABase
from ..Exceptions import NoLimitingSubstrateFound


class TimeStepDynamicFBABase:
    """The base class for Dynamic FBA"""

    m_biomass_concentrations: dict[str, list[float]] = {}
    m_metabolite_concentrations: dict[str, list[float]] = {}

    def set_initial_concentrations(
        self, model: Model, initial_concentrations: dict[str, float]
    ):
        """Set the initial concentrations of the metabolites

        Uses the exchange reaction lower bounds as initial concentrations
        If people want to change it they have to change this use initial
        conditions or change the lower bounds
        This was implemented in such that you dont need to manually
        add all external metabolites in the initial conditions

        Args:
            model (Model): A cbmpy.CBModel
            initial_concentrations (dict[str, float]): dictionary containing
            the initial concentrations as key value pairs
        """

        for exchange in model.getExchangeReactionIds():
            reaction: Reaction = model.getReaction(exchange)
            species_id: str = reaction.reagents[0].getSpecies()
            if species_id not in self.m_metabolite_concentrations.keys():
                if species_id in initial_concentrations.keys():
                    self.m_metabolite_concentrations[species_id] = [
                        initial_concentrations[species_id]
                    ]
                else:
                    self.m_metabolite_concentrations[species_id] = [
                        -reaction.getLowerBound()
                    ]

    def simulate(
        self,
        dt: float,
        epsilon=0.001,
        kinetics_func=None,
        deviate=None,
        deviation_time=0,
    ):
        """The simulation methods, calculates the fluxes for each dt
        as specified in the function call.

        Each specific dynamic FBA implementation performs it's own simulation
        with their own set of rules

        Args:
            dt (float): The time step

            epsilon (float, optional): When the solution < epsilon stop the

            simulation. Defaults to 0.001.

            kinetics_func (_type_, optional): A user function that can be used.
            to calculate kinetics.  Defaults to None.

            deviate (_type_, optional): If you want to apply changes during
            simulation you can provide a function which should accept
            the model
            the dictionary of biomass concentrations
            the dictionary of metabolite concentrations
            dt.
            Defaults to None.

            deviation_time (int, optional): the time step the deviation
            function should be called. Defaults to 0.
        """
        pass

    def update_reaction_bounds(self, kinetics_func) -> None:
        """Update all reaction bounds using the new
        concentrations.

        Args:
            kinetics_func (_type_): A user defined function
            that calculates the new upper bound for a reaction.
            The function should accept a cbmpy.CBmodel, a
            string (reaction_id) a dictionary of biomasses concentrations
            and a dictionary of metabolite
            concentrations
        """
        pass

    def update_importer_bounds(self) -> None:
        """Update the upper bounds of the importer reaction"""
        pass

    def update_concentrations(self) -> None:
        """Update all metabolite concentrations"""
        pass

    def update_biomasses(self) -> None:
        """Update all Biomass concentration"""
        pass

    def reset_dt(self) -> float:
        """Reset the simulation and calculate dt_hat"""
        pass

    def check_solution_feasibility(self) -> str:
        """If the current solution isn't feasible due to
        a lack of metabolite concentrations return the
        metabolite which dropped furthest below zero
        the concentration of this metabolite will be used
        to calculate dt_hat.
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
        """Given an importer reaction get the
        species concentrations that this importer imports

        Args:
            rid (str): the reaction id of the importer

            transporters (Transporters): A Transporters object

        Returns:
            list[float]: list of species ids
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
        """Given a reaction calculate the Michaelis Menten
        Kinetics


        Args:
            reaction (Reaction): A reaction objet

            X (float): The biomass of the species the reaction belongs to

            transporters (Transporters): A Transporters object of the model
        """
        rid: str = reaction.getId()
        sid, km, vmax = kinetics.get_kinetics(
            rid,
        )

        # If the reaction is an importer and no limiting substrate
        # was supplied use the min of all substrates
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
