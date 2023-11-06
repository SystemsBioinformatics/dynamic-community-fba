4. The community model 
=======================

4.1. Definition
---------------

Here we will introduce the concept of the community model, which plays a vital role in upcoming modeling techniques. 
Strictly put, the community model represents the combined stoichiometry matrices of N Genome-Scale Metabolic Models (GSMMs) 
and thus can again be used to perform FBA.

To streamline the handling of multiple models, we have developed the ``CommunityModel`` class. The class offers a structured representation of combined metabolic networks, integrating stoichiometric 
information from individual GSMMs/SBML models [1]_ [2]_. By merging the information of each model, we can examine the intricate world microbial communities, their dynamics, and phenotypic potentials.

The ``CommunityModel`` class is a direct subclass of the ``cbmpy.CBModel`` class. By building upon the foundational ``cbmpy.CBModel``, the ``CommunityModel`` inherits all of its capabilities, ensuring compatibility and ease of use. 
All functionalities available in the ``cbmpy.CBModel`` class are retained in the ``CommunityModel``. Additionally, the latter can be effortlessly exported to an SBML file, enabling its analysis in various tools.

In this chapter, we elaborate on the principles and guidelines essential for constructing an accurate community model. By doing so, we aim to maximize its utility and prevent inadvertent errors. 
To further illustrate the adaptability of the "CommunityModel" class, we provide practical, hands-on examples.


4.2. Creating and modifying the community model
-----------------------------------------------

To set up a ``CommunityModel``, you must provide a list of N GSMMs along with their corresponding biomass reaction IDs:

.. code-block:: python

    from dcFBA.Models import CommunityModel
    from dcFBA.DefaultModels import read_default_model
    from cbmpy.CBModel import Model

    model1:Model = read_default_model("e_coli_core")
    model2:Model = read_default_model("strep_therm")

    biomass_reaction_model_1 = (
        "R_BIOMASS_Ecoli_core_w_GAM"  # Biomass reaction id of ecoli
    )

    biomass_reaction_model_2 = "R_biomass_STR"  # Biomass reaction id of strep.
    community_model: CommunityModel = CommunityModel(
        [model1, model2],
        [
            biomass_reaction_model_1,
            biomass_reaction_model_2,
        ],
    )


If you don't know the biomass reaction ID of the model, there are a couple of ways to find it. 
First of all, you can look it up in the SBML model. 
Alternatively, you can assess if the model's active objective 
function is set to the biomass reaction with :code:`model.getActiveObjectiveReactionIds()`. 
This returns a list of objective IDs where you can look for the biomass reaction.

In our example, we use the GSMMs of: *E. coli core metabolism* and *Streptococcus thermophilus*. 
If desired, you can provide an optional list of alternative IDs to refer to the single models.
This can be done by supplying a list of strings as the third argument to the function. 
If no alternative IDs are provided, the community matrix will default to using the original model identifiers.

.. warning::
    If you want to create a community model of N identical models, it is mandatory to specify the alternative ids.
    Without providing alternative IDs, both models would have the same ID during the community model creation process. 
    As a result, the code would be unable to distinguish between the two models, leading to undesired behavior.


4.3. Rules 
----------

During the initialization of the community model a new ``cbmpy.CBModel`` object is created. In the new object all compartments,
genes, reactions, reagents and species are copied from the provided models. To build a functional community model, we made the following considerations for each section of the old models:

Compartments:
*************

In our community model, we utilize a compartmentalized design to accurately simulate the natural structure and behavior of microbial communities.
Therefore we handle the ``Compartments`` as follows: 

**Unified External Environment:** 
All organisms in the model share a common external environment.

**Compartmentalization:** 
While organisms share the same external environment, we consider each individual organism to be a separate compartment in the model. Furthermore we retain the 
compartments that were already present in the original models.

**Design Specifics:**
The `external` or `e` compartment from individual models is singular and shared in the community model. This ensures that all external species from all models coexist. When compartments from individual models are copied into the community model, they are appended with a suffix based on the provided ID. This procedure helps in tracking and identifying which organism each compartment corresponds to.

Mathematically, the total number of compartments in the community model is given by:

.. math::

   \Sigma_{i=1}^{n} (c_i-1) + 1

Where :math:`c_i` represents the number of compartments in each individual model \(i\).

Reactions:
**********

For the ``Reactions`` we consider the following: 

**Reaction Duplication:** 
All reactions, are duplicated within the community model. To maintain traceability, each reaction ID is supplemented with the originating model ID.

**Handling Exchange Reactions:** 
A distinctive approach is employed for exchange reactions. When identical exchange reactions are present in multiple models,
only a single instance is retained in the combined community model.


Reagents and Species
********************

Lastly the ``Reagents`` and their corresponding ``Species``:

**External Metabolites:** 
Given that the community shares all external metabolites, only a single instance of each external metabolite is retained in the combined model.

**Species Duplication and Identification:** 
All other Species are copied to the community model. However, to maintain clarity and avoid confusion, 
a species receives a distinguishing suffix only when it's present in more than one original model.

.. warning:: 
    It is crucial to verify that the identical reactions and species within different models have consistent IDs before 
    creating the community model. This is particularly significant for exchange reactions and species localized in the 
    extracellular space. If these IDs are not uniform, despite referring to the same reactions or species, the CommunityModel 
    class cannot determine their equivalence accurately.

    Please ensure that the corresponding IDs for these reactions and species are compatible to guarantee the proper 
    functioning of the community model.


.. [1] Khandelwal, R. A., Olivier, B. G., Röling, W. F. M., Teusink, B., and Bruggeman, F. J. (2013). Community flux balance analysis for microbial consortiaat balanced growth. PLoS ONE, 8(5):e64567.
.. [2] Stolyar, S., Dien, S. V., Hillesland, K. L., Pinel, N., Lie, T. J., Leigh, J. A., and Stahl, D. A. (2007). Metabolic modeling of a mutualistic microbial community. Molecular Systems Biology, 3(1).
