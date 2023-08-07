7. Advanced Topics
==================

This is the Advanced Topics section of our Python documentation. Here, we explore how you can manipulate the simulation during it's run and we explore the
"Kinetics Object" which can be used for Michaelis Menten kinetics. Next we discuss the "Kinetics function." 
These tools enable you to add more biological information in the modelling.

Let's dive in!

Manipulating the simulation
---------------------------

Each simulation for objects of the ``TimeStepDynamicFBABase`` base class can be given an extra function enabling the user to manipulate the simulation on the go.
By manipulating you can think of adding or removing external metabolites from the system. But following some simple rules nearly everything can be changed during simulation.

Here we will give an example where we grow two strains of *E coli.* on Pyruvate and add L-Glutamate when the concentration of Pyruvate is half it's initial concentration.

First we build a ``DynamicJointFBA`` model:

.. code-block:: python

    import cbmpy
    import numpy
    from cbmpy.CBModel import Model
    from DCFBA.Models import CommunityModel
    from DCFBA.DynamicModels import DynamicJointFBA

    # Build the model and set all bounds
    model1: Model = cbmpy.loadModel("data/bigg_models/e_coli_core.xml")
    model1.getReaction("R_GLUt2r").setUpperBound(11)
    model1.getReaction("R_PYRt2").setUpperBound(6)
    model_2 = model1.clone()

    model_2.getReaction("R_PYRt2").setUpperBound(10)
    model_2.getReaction("R_GLUt2r").setUpperBound(8)

    combined_model = CommunityModel(
        [model1, model_2],
        ["R_BIOMASS_Ecoli_core_w_GAM", "R_BIOMASS_Ecoli_core_w_GAM"],
        ["ecoli_1", "ecoli_2"],
    )

    combined_model.getReaction("R_EX_pyr_e").setLowerBound(-numpy.inf)

    #Set initial concentrations
    dynamic_fba = DynamicJointFBA(
        combined_model,
        [0.1, 0.1],
        {"M_pyr_e": 10, "M_glu__L_e": 0, "M_glc__D_e": 0},
    )
    
Next we will define our deviation function. Remember that the function, no matter if you use it or not, needs to accept the ``DynamicJointFBA`` model, the ``time steps`` of the simulation 
and a ``condition`` which we will explain in a minute.

.. code-block:: python

    def deviate_func(DFBA: DynamicJointFBA, used_time, condition) -> int:
    if (
        DFBA.m_metabolite_concentrations["M_pyr_e"][-1] <= 5.0
        and condition < 1
    ):
        combined_model.getReaction("R_EX_glu__L_e").setLowerBound(-numpy.inf)

        DFBA.m_metabolite_concentrations["M_glu__L_e"][-1] = 30

        return 1

    return condition

Here we define the ``deviate_func`` as follows, run the function when the concentration of the dynamic FBA object runs under 5.0 and run it only once. If the concentration is 
not below oir equal to 5.0 return the condition. The ``condition`` parameter is by default set to 0 at the start of the simulation and is always passed to the deviate function.
This is done such that for example you can rerun the function for consecutive time steps. The value returned by the function will always be added to the global condition variable. 
By doing so we can for example add Glucose to the system for N consecutive runs by modifying the conditional statement. 
You can also skip on using the condition and define your own conditions when the function should run. 

We can now pass the function to the simulation method and plot the results:

.. code-block:: python 

    import matplotlib.pyplot as plt

    T, metabolites, biomasses, _ = dynamic_fba.simulate(
        0.1, deviate=deviate_func
    )

    import matplotlib.pyplot as plt


    plt.figure(1)
    plt.plot(T, metabolites["M_glu__L_e"], color="blue", label="[glu__L]")
    plt.plot(T, metabolites["M_pyr_e"], color="orange", label="[pyr]")

    plt.xlabel("Time")
    plt.ylabel("Concentration")
    plt.legend()

    plt.figure(2)
    plt.plot(T, biomasses["ecoli_1"], color="blue", label="Biomass model 1")
    plt.plot(T, biomasses["ecoli_2"], color="orange", label="Biomass model 2")

    plt.xlabel("Time")
    plt.ylabel("Concentration")
    plt.legend()

    plt.show()

.. image:: ../_static/images/Deviate_function_metabolites.png
    :width: 500px
    :align: center
    :alt: Metabolite concentrations   

.. image:: ../_static/images/deviation_function_biomasses.png
    :width: 500px
    :align: center
    :alt: Metabolite concentrations 


The Kinetics Object for Michaelis Menten
----------------------------------------



The Kinetics function
---------------------