# Dynamic Community FBA

Introducing Dynamic Community FBA (dcFBA): The first ever Python Package to effortlessly model Microbial Communities dynamically!

This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details.

## About

Dynamic Community FBA (dcFBA) can be used to model both microbial communities as single organisms using Genome Scale Metabolic models (GSMMs). The package is an extension of Built upon the foundations of [cbmpy](https://github.com/SystemsBioinformatics/cbmpy), 
and is compatible with SBML and COBRApy models. 
It is designed to perform three distinct types of dynamic modelling methods: 

1. **Dynamic Joint FBA** - Incrementally updates the concentrations of biomass and metabolites within the combined stoichiometric matrix of the supplied models.
2. **Dynamic Parallel FBA** -Incrementally updates the concentrations of biomass and metabolites, while performing FBA on the separate models 
3. **EndPointFBA** - Duplicate the CommunityMatrix N times and performs FBA on the community time dependant stoichiometric matrix 
^*For further details and mathematical description see [^1]

Whether you're modelling or parasitic or costly cross-feeding behaviors in microbial communities, dcFBA provides a simple and effective solution.


## Install

### Prerequisites

- Python 3 >= 3.8, < 3.11
- [cbmpy](https://github.com/SystemsBioinformatics/cbmpy)
- [cplex](https://www.ibm.com/products/ilog-cplex-optimization-studio)

### Installation Steps

1. Install using pip:

    ```bash
    pip install dcFBA
    ```

   
## Usage

For basic usage of [cbmpy](https://pythonhosted.org/cbmpy/modules_doc.html) and [dcFBA](https://dynamic-community-fba.readthedocs.io/en/latest/2_getting_started/home.html) we refer to the respective documentation pages.
[^1]: [Your citation for Dynamic Joint FBA goes here.]