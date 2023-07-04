4. The community model 
=======================

4.1. Definition
---------------

Here we will introduce the concept of the community model, which plays a vital role in upcoming modeling techniques. 
In simple terms, the community model represents the combined stoichiometry matrices of N Genome-Scale Metabolic Models (GSMMs) 
for performing Flux Balance Analysis (FBA).

To streamline the handling of multiple models, we have developed the ``CommunityModel`` class. The class offers a structured representation of combined metabolic networks, integrating stoichiometric 
information from individual GSMMs. This facilitates in-depth analysis of complex microbial communities, their dynamics, 
and metabolic potentials.

The documentation delves into the underlying principles and rules governing the community model's construction, ensuring 
effective utilization and preventing erroneous model creation. Additionally, we provide practical usage examples to further 
illustrate the versatility of the ``CommunityModel`` class.

4.2. Creating and modifying the community model
-----------------------------------------------

To initialize a ``CommunityModel`` you have to give a list of N GSMMs as well as there biomass reaction ids:

.. code-block:: python

    import cbmpy
    from endPointFBA import CommunityModel

    model1 = cbmpy.loadModel("data/bigg_models/e_coli_core.xml")
    model2 = cbmpy.loadModel("data/bigg_models/strep_therm.xml")

    biomass_reaction_model_1 = "R_BIOMASS_Ecoli_core_w_GAM"
    biomass_reaction_model_2 = "R_biomass_STR"
    community_model: CommunityModel = CommunityModel(
        [model1, model2],
        [
            biomass_reaction_model_1,
            biomass_reaction_model_2,
        ],
    )

If you are not familiar with the biomass reaction ID, there are a couple of ways to identify it. 
First, you can look it up in the SBML model. Alternatively, you can check if the active objective function of the model is 
set to the biomass reaction by using the following code: :code:`model.getActiveObjectiveReactionIds()`. 
This will display a list of objective IDs, and you can verify if the biomass reaction is included.

In this example, we have utilized two models included in the package: *E. coli core metabolism*
and *Streptococcus thermophilus*. If desired, you can provide an optional list of alternative IDs to refer to the single models.
This can be done by supplying a list of strings as the third argument to the function. 
If no alternative IDs are provided, the community matrix will default to using the model identifiers.
Lastly, you have the option to assign an ID of your own choosing to the new combined model. If not provided the default ID will be used.

It is important to note that the ``CommunityModel`` class is derived from the base ``cbmpy.CBModel`` class, meaning that the 
``CommunityModel`` is still an instance of ``cbmpy.CBModel``. With the newly initialized object we can obtain even more information.
Check the API for all functionality. 

.. warning::
    If you want to create a community model of N identical models, it is mandatory to specify the alternative ids.
    Without providing alternative IDs, both models would have the same ID during the community model creation process. 
    As a result, the code would be unable to distinguish between the two models, leading to undesired behavior.


4.3. Rules 
----------

Maybe move this to an advanced section?

During the initialization of the community model a new ``cbmpy.CBModel`` object is created. In the new object all compartments,
reactions, reagents and species are copied from the provided model. The following considerations were made for each:

Compartments:
*************

The new model comprises a total of :math:`\Sigma_{i=1}^{n} (c_i-1) + 1` compartments, where c represents the number 
of compartments in each model i.

In the new model, all compartments from the individual models are copied, except for the `external` or `e` compartment. 
The external compartment is reserved to be added at the end. All other compartment ids get a prefix of the id of the model
they belonged to such that we know which compartment corresponds to which organism.

This design ensures that there is only one external compartment in which all species and reactions from all models coexist.

By following this approach, the new model achieves compartmental organization while consolidating all species and reactions 
within a unified external compartment. 

Reactions:
**********

Just like the compartments, all reactions are duplicated, and each reaction ID is augmented with its corresponding model ID. This 
holds up for all reactions except for the exchange reactions. If two models share an exchange reaction only one is saved in the new 
combined model such that we have more control over all. 

Reagents and species
********************

Before the reactions can be copied to the new model there is a check for which species occur in more than one model.
For these species a new species in the original model is created and all reactions associated to this species have there reagents changed. 
By doing so the initial models already have everything set correctly into place to have the reactions copied. 

In contrast with the compartment IDs, and the reaction IDs the species IDs are not changed by default. But only if the species id 
occurs in two or more models. This is done since we can already quickly lookup to which original model the species belonged by checking 
the compartment which it is in.




.. note:: 
    It is crucial to verify that the identical reactions and species within different models have consistent IDs before 
    creating the community model. This is particularly significant for exchange reactions and species localized in the 
    extracellular space. If these IDs are not uniform, despite referring to the same reactions or species, the CommunityModel 
    class cannot determine their equivalence accurately.

    Please ensure that the corresponding IDs for these reactions and species are harmonized to guarantee the proper 
    functioning of the CommunityModel.