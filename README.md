# Dynamic Community FBA

Welcome to Dynamic Community FBA (dcFBA): The Python package that makes modeling Microbial Communities dynamically a breeze!

This project is licensed under the MIT License - please refer to the [LICENSE](LICENSE) file for more details.

## About

Dynamic Community FBA (dcFBA) is a versatile tool designed for modeling microbial communities as single organisms using Genome Scale Metabolic Models (GSMMs). This package builds upon the solid foundation of [cbmpy](https://github.com/SystemsBioinformatics/cbmpy) and seamlessly integrates with SBML and COBRApy models. dcFBA empowers users with three distinct dynamic modeling methods:

1. **Dynamic Joint FBA** - Incrementally updates the concentrations of biomass and metabolites within the combined stoichiometric matrix of the provided models.

2. **Dynamic Parallel FBA** - Simultaneously updates the concentrations of biomass and metabolites while performing FBA on individual models.

3. **EndPointFBA** - Duplicates the CommunityMatrix N times and performs FBA on the community's time-dependent stoichiometric matrix.

For a comprehensive understanding of these methods and their underlying mathematics, please consult [^1].

Whether you're exploring parasitic interactions or investigating costly cross-feeding behaviors in microbial communities, dcFBA offers an elegant and efficient solution.

## Installation

### Prerequisites

Before installing dcFBA, make sure you have the following prerequisites in place:

- Python >= 3.8, < 3.11
- [cbmpy](https://github.com/SystemsBioinformatics/cbmpy)
- [cplex](https://www.ibm.com/products/ilog-cplex-optimization-studio)

### Installation Steps

To install dcFBA, follow these simple steps:

1. Install using pip:

    ```bash
    pip install dcFBA
    ```

## Usage

For basic usage examples and detailed documentation on both [cbmpy](https://pythonhosted.org/cbmpy/modules_doc.html) and [dcFBA](https://dynamic-community-fba.readthedocs.io/en/latest), please refer to their respective documentation pages.

[^1]: [Cite paper.]
