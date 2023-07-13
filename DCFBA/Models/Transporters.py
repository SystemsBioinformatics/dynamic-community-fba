from cbmpy.CBModel import Model, Reaction, Species


class Transporters:
    """Accessing Importers, Exporters, and Transporters

    The purpose of this class is to provide easy access to importers,
    exporters, and transporters within a system. Importers represent reactions
    that uptake metabolites from the external space, while exporters represent
    reactions that excrete metabolites to the external space. Both importers
    and exporters are stored as dictionaries, where the
    reaction identifiers (rid) serve as keys, and the associated species are
    stored as values.

    This class helps users identify, work with, and analyze importers,
    exporters, and transporters, which are frequently accessed during various
    computational tasks.

    Usage:
    Add new importers or exporters using the add_importer() and add_exporter()
    Retrieve the combined importers and exporters dictionary using get_transporters()
    Retrieve importers using the get_importers() method,
    and optionally iterate over the results using the iterable flag.
    Retrieve exporters using the get_exporters() method, and optionally
    iterate over the results using the iterable flag.
    Obtain the species associated with a specific reaction using the
    get_reaction_species() method.
    Check if a reaction is an importer or exporter using the is_importer() and
    is_exporter() methods.
    Access specific lists of importers or exporters using the provided methods.
    """

    m_importers: dict[str, list[str]] = {}
    m_exporters: dict[str, list[str]] = {}

    def __init__(self, model: Model) -> None:
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
        self.m_importers[rid] = species

    def add_exporter(self, rid: str, species: list[str]) -> list[str]:
        self.m_exporters[rid] = species

    def get_transporters(self) -> dict[str, list[str]]:
        """Get the combined importers and exporters

        Returns:
            dict[str, list[str]]: dict with reactions and all species
        """
        return {**self.m_importers, **self.m_exporters}

    def get_importers(self, iterable=False) -> dict[str, list[str]]:
        """Get all the import reactions of the model

        Returns:
            dict[str, list[str]]: Return the entire dictionary
        """
        if iterable:
            return self.m_importers.items()

        return self.m_importers

    def get_exporters(self, iterable=False) -> dict[str, list[str]]:
        """Get the exporters of the model

        Returns:
            dict[str, list[str]]:: The dictionary containing the export
            reactions
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
