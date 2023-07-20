5. Dynamic Joint FBA 
====================

The multiple metabolic models of different or the same organism that were combined in the ``CommunityModel`` as described
in the previous chapter can now be used in dynamic joint FBA. The community model incorporates the metabolic reactions and 
interactions between the organisms, allowing for the study of their collective behavior and the emergent properties of the 
community as a whole.

The joint FBA approach enables the investigation of metabolic exchanges, such as the exchange of nutrients or byproducts, 
between organisms within the community. By simulating the community-level metabolic interactions, researchers can gain 
insights into the dependencies, cooperation, competition, and overall dynamics of the organisms in the community.

Joint FBA
---------

After defining the ``CommunityModel`` it is easy to refine it in such a way that we can perform a Joint FBA
The only thing left to do is append the biomass reaction of each individual model to create the so called `Community biomass`.
We can now define the reaction :literal:`X_c ->` where  ``X_c`` is the total sum of the individual models' biomasses.

To perform joint FBA run the following: 

.. code-block:: python
   
    import cbmpy
    from cbmpy.CBModel import Model
    from DCFBA.Models import CommunityModel
    from DCFBA.DynamicModels import DynamicJointFBA

    model1: Model = cbmpy.loadModel("data/bigg_models/e_coli_core.xml")
    model1.getReaction("R_GLCpts").setUpperBound(10)

    model_2 = model1.clone()
    
    # Ecoli 2 imports glucose slower
    model_2.getReaction("R_GLCpts").setUpperBound(8)

    combined_model = CommunityModel(
        [model1, model_2],
        ["R_BIOMASS_Ecoli_core_w_GAM", "R_BIOMASS_Ecoli_core_w_GAM"],
        ["ecoli_1", "ecoli_2"],
    )  # Create a CommunityModel of two  E. coli strains competing for resources


    # Create the joint FBA object with initial biomasses and the initial concentration of glucose
    dynamic_fba = DynamicJointFBA(
        combined_model,
        [0.1, 0.1],
        {"M_glc__D_e": 10},
    )

    # Perform FBA on the new joint FBA model object
    solution = cbmpy.doFBA(dynamic_fba.get_joint_model())
    print(solution)



Making it dynamic
-----------------

Dynamic FBA is an extension of the traditional FBA approach that incorporates the element of 
time. In dynamic FBA, a specific time step `dt` is selected, and the concentrations of external
metabolites and biomass concentrations are calculated at each time point.
This enables the modeling and analysis of dynamic processes such as metabolic fluxes, 
nutrient uptake, and product secretion over time. 

The same technique can be applied for using the previously described Joint FBA which we call
`Dynamic Joint FBA`

To perform the joint FBA over time using the ``DynamicJointFBA`` model:

.. code-block:: python

    T, metabolites, biomasses, fluxes = dynamic_fba.simulate(0.1)

The ``simulate`` method returns a tuple with four elements. First, it provides a list of time points for the simulation. 
Second, it returns a ``Dictionary`` containing the species_ids and their corresponding concentrations at each time point. Third, we get a 
``Dictionary`` containing the biomasses of the models, accessed through their IDs. Lastly we get the fluxes of all reactions for each time point.

You can now easily plot the species concentration over time:

.. code-block:: python

    plt.plot(T, metabolites["M_glc__D_e"], color="blue", label="[Glucose]")
    
    plt.xlabel("Time")
    plt.ylabel("Concentration")
    plt.legend()
    plt.show()
  
And the biomasses of both species

.. code-block:: python
    
    plt.plot(T, biomasses["ecoli_1"], color="blue", label="Biomass model 1")
    plt.plot(T, biomasses["ecoli_2"], color="orange", label="Biomass model 2")

    
    plt.xlabel("Time")
    plt.ylabel("Concentration")
    plt.legend()
    plt.show()

.. tip::

    If you create a ``DynamicJointFBA`` object with a ``CommunityModel`` build from just one organism and call the simulate function you
    perform just regular dynamic FBA!


