3. Parallel Dynamic FBA
=======================

Parallel dynamic FBA is a simple yet highly effective method for modeling microbial interactions. The key idea is to perform FBA on individual models while keeping track of a *pool* of external metabolites and their concentrations over time.

Example
-------
To run a simulation with Parallel Dynamic FBA, you first need to define your ``ParallelDynamicFBA`` model.

.. code-block:: python

    from cbmpy.CBModel import Model
    from DCFBA.ToyModels import model_a, model_b
    from DCFBA.DynamicModels.DynamicParallelFBA import DynamicParallelFBA

    # Load two toy models
    model_a_instance: Model = model_a.build_toy_model_fba_A()
    model_a_instance.getReaction("B_exchange").setLowerBound(-2)
    model_a_instance.getReaction("A_exchange").setLowerBound(0)

    model_b_instance: Model = model_b.build_toy_model_fba_B()
    model_b_instance.getReaction("B_exchange").setLowerBound(0)
    model_b_instance.getReaction("A_exchange").setLowerBound(-1)

    # Create an instance of the model
    parallelModel = DynamicParallelFBA(
        [model_a_instance, model_b_instance], [1, 1], {"S_e": 100, "A_e": 1, "B_e": 2}
    )

The ``DynamicParallelFBA`` class accepts a list of N models. In this example, we use two toy models described in the `toy model section`.
Next, we specify the initial concentrations of both models and the initial concentrations of the external species.

To run the simulation, use the following code:

.. code-block:: python

    T, metabolites, biomasses = parallelModel.simulate(0.1, 0.05)

The ``simulate`` method returns a tuple with three elements. First, it provides a list of time points for the simulation. 
Next, it returns a ``Dictionary`` containing the species_ids and their corresponding concentrations at each time point. Lastly, we get a 
``Dictionary`` containing the biomasses of the models, accessed through their IDs.

Using simple plotting functions from ``matplotlib.pyplot``, you can visualize the results easily:

.. code-block:: python

    # Plotting the external metabolite concentrations of species A and B
    plt.plot(T, metabolites["A_e"], color="blue", label="Species A")
    plt.plot(T, metabolites["B_e"], color="orange", label="Species B")
    plt.legend()

    plt.show()

.. image:: ../_static/images/ParallelFBA_metabolites.png
    :width: 500px
    :align: center
    :alt: Metabolite concentrations


Alternatively, you can plot the biomasses over time:

.. code-block:: python

    plt.plot(T, biomasses["Organism_A"], color="blue", label="Model A")
    plt.plot(T, biomasses["Organism_B"], color="orange", label="Model B")

    plt.legend()

    plt.show()

.. image:: ../_static/images/ParallelFBA_Biomass_concentrations.png
    :width: 500px
    :align: center
    :alt: Biomass concentrations
