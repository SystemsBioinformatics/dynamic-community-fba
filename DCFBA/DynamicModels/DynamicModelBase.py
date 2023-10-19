class DynamicModelBase:
    """Base class for the Dynamic FBA classes"""

    def simulate(self):
        pass

    # TODO fix this
    # IDEA of the refactor.
    # place all properties on nthis base class and than
    # make private setters, also in epFBA, such that many rdeundant code can be removed
    # eg in epFBA result you can just ask for the metabolites which got stripped of their time ids
    # most of the code in djFBA en dpFBA will be the same however also more easy acces to fluxes etc
    # @property
    # def fluxes(self):
    #     return self._fluxes

    # @property
    # def metabolites(self):
    #     return._metabolites

    # def get_flux_values(self, rid):
    #     pass

    # def get_community_growth_rate(self):
    #     pass

    # def get_relative_abundance(self):
    #     pass
