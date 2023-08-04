5. Dynamic Joint FBA 
====================

The multiple metabolic models of different (or the same) organism that were combined in the ``CommunityModel`` as described
in the previous chapter can now be used in dynamic joint FBA. The community model incorporates the metabolic reactions and 
interactions between the organisms, allowing for the study of their collective behavior and the emergent properties of the 
community as a whole.

The joint FBA approach enables the investigation of metabolic exchanges, such as the exchange of nutrients or byproducts, 
between organisms within the community. By simulating the community-level metabolic interactions, researchers can gain 
insights into the dependencies, cooperation, competition, and overall dynamics of the organisms in the community.


Making it dynamic!
------------------

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

.. .. tip::

..     If you create a ``DynamicJointFBA`` object with a ``CommunityModel`` build from just one organism and call the simulate function you
..     perform just regular dynamic FBA!


