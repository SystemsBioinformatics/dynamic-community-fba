6. Parallel Dynamic FBA
=======================

The second method to model microbial interactions available through this package is Dynamic Parallel FBA. 
The key idea is to perform FBA on individual models while keeping track of a *pool* of external metabolites from which all models grow.
Again we keep track of the overall concentrations of biomass and external metabolites and update this using the set time interval.

Example
-------
To run a simulation with Parallel Dynamic FBA, you first need to define your ``ParallelDynamicFBA`` model.

.. code-block:: python

    import matplotlib.pyplot as plt
    import cbmpy
    from cbmpy.CBModel import Model
    from DCFBA.Models import CommunityModel
    from DCFBA.DynamicModels import DynamicParallelFBA


    model1: Model = cbmpy.loadModel("models/bigg_models/e_coli_core.xml") #load ecoli model
    model2: Model = cbmpy.loadModel("models/bigg_models/strep_therm.xml") #load streptococcus model

    #Exchange reaction are not set correctly, here we fix this
    for rid in model1.getReactionIds():
        if rid.startswith("R_EX"):
            reaction = model1.getReaction(rid)
            reaction.is_exchange = True

    for rid in model2.getReactionIds():
        if rid.startswith("R_EX"):
            reaction = model2.getReaction(rid)
            reaction.is_exchange = True

    #Give the two models different glucose uptake rates
    model1.getReaction("R_GLCpts").setUpperBound(10)
    model2.getReaction("R_GLCpts").setUpperBound(6)

    parallel_fba = DynamicParallelFBA(
        [model1, model2],
        [1.0, 1.0],
        {"M_glc__D_e": 100, "M_gal_e": 0, "M_lcts_e": 100},
    )


The ``DynamicParallelFBA`` class accepts a list of N models. In this example, we use two toy models described in the `toy model section`.
Next, we specify the initial concentrations of both models and the initial concentrations of the external species.

To run the simulation, use the following code:

.. code-block:: python

    T, metabolites, biomasses = parallelModel.simulate(0.1)

The ``simulate`` method returns a tuple with three elements. First, it provides a list of time points for the simulation. 
Next, it returns a ``Dictionary`` containing the species_ids and their corresponding concentrations at each time point. Lastly, we get a 
``Dictionary`` containing the biomasses of the models, accessed through their IDs.

Using simple plotting functions from ``matplotlib.pyplot``, you can visualize the results easily:

.. code-block:: python

    #plot external metabolites
    plt.plot(T, metabolites["M_glc__D_e"], color="blue", label="[Glucose]")
    plt.plot(T, metabolites["M_lcts_e"], color="orange", label="[Lactose]")

    plt.xlabel("Time")
    plt.ylabel("Concentration")
    plt.legend()
    plt.show()

.. image:: ../_static/images/ParallelFBA2_metabolites.png
    :width: 500px
    :align: center
    :alt: Metabolite concentrations


Alternatively, you can plot the biomasses over time:

.. code-block:: python

    #Plot biomasses
    plt.plot(T, biomasses[model1.getId()], color="orange", label="ecoli")
    plt.plot(T, biomasses[model2.getId()], color="blue", label="strep")

    plt.xlabel("Time")
    plt.ylabel("Concentration")
    plt.legend()
    plt.show()


.. image:: ../_static/images/ParallelFBA2_Biomass_concentrations.png
    :width: 500px
    :align: center
    :alt: Biomass concentrations
