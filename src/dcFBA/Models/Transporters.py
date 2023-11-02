from cbmpy.CBModel import Model, Reaction, Species


class Transporters:
    """
    Facilitate access and operations on importers, exporters, and transporters of a metabolic model.

    This class is designed to manage and interact with importers (reactions that uptake metabolites from external space)
    and exporters (reactions that excrete metabolites to external space). Both importers and exporters are stored
    as dictionaries where reaction identifiers (rid) are keys and the associated species are values.

    Attributes:
        m_importers (dict[str, list[str]]): Dictionary storing the import reactions and their corresponding species.
        m_exporters (dict[str, list[str]]): Dictionary storing the export reactions and their corresponding species.

    """

    m_importers: dict[str, list[str]] = {}
    m_exporters: dict[str, list[str]] = {}

    def __init__(self, model: Model) -> None:
        """Initialize the Transporters class based on a metabolic model."""

        importers = {}
        exporters = {}
        for rid in model.getReactionIds():
            reaction: Reaction = model.getReaction(rid)
            if not reaction.is_exchange:
                for reagent in reaction.reagents:
                    sid: str = reagent.getSpecies()
                    species: Species = model.getSpecies(sid)
                    if species.getCompartmentId() == "e":
                        if reagent.coefficient == -1:
                            # Reaction imports external metabolites
                            importers[rid] = importers.get(rid, []) + [sid]
                        elif reagent.coefficient == 1:
                            exporters[rid] = exporters.get(rid, []) + [sid]
        self.m_importers = importers
        self.m_exporters = exporters

    def add_importer(self, rid: str, species: list[str]) -> None:
        """Add a new importer reaction along with its associated species."""

        self.m_importers[rid] = species

    def add_exporter(self, rid: str, species: list[str]) -> list[str]:
        """Add a new exporter reaction along with its associated species."""
        self.m_exporters[rid] = species

    def get_transporters(self) -> dict[str, list[str]]:
        """
        Retrieve combined dictionary of both importers and exporters.

        Returns:
            dict[str, list[str]]: Combined dictionary of importers and exporters.
        """

        return {**self.m_importers, **self.m_exporters}

    def get_importers(self, iterable=False) -> dict[str, list[str]]:
        """
        Retrieve the import reactions and their associated species.

        Args:
            iterable (bool): If True, return as dictionary items. Default is False.

        Returns:
            dict[str, list[str]]: Dictionary of import reactions and species.
        """
        if iterable:
            return self.m_importers.items()

        return self.m_importers

    def get_exporters(self, iterable=False) -> dict[str, list[str]]:
        """
        Retrieve the export reactions and their associated species.

        Args:
            iterable (bool): If True, return as dictionary items. Default is False.

        Returns:
            dict[str, list[str]]: Dictionary of export reactions and species.
        """
        if iterable:
            return self.m_exporters.items()

        return self.m_exporters

    def get_reaction_species(self, rid) -> list[str] | None:
        """Given a reaction id return the species belonging to it

        Args:
            rid (_type_): reaction id

        Returns:
            list[str] | None: the species or None if rid was not in the model
        """
        if rid in self.m_exporters.keys():
            return self.m_exporters[rid]
        if rid in self.m_importers.keys():
            return self.m_importers[rid]

        return None

    def get_importers_reactions(self) -> list[str]:
        return list(self.m_importers.keys())

    def get_exports_reactions(self) -> list[str]:
        return list(self.m_exporters.keys())

    def get_importers_species(self, rid: str) -> list[str]:
        return self.m_importers[rid]

    def get_exporter_species(self, rid: str) -> list[str]:
        return self.m_exporters[rid]

    def is_importer(self, rid: str) -> bool:
        return rid in self.m_importers.keys()

    def is_exporter(self, rid: str) -> bool:
        return rid in self.m_exporters.keys()
