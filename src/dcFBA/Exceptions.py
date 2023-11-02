class NotInCombinedModel(Exception):
    def __init__(self, message):
        super().__init__(message)


class NoLimitingSubstrateFound(Exception):
    def __init__(self, message):
        super().__init__(message)


class SpeciesNotFound(Exception):
    def __init__(self, message):
        super().__init__(message)
