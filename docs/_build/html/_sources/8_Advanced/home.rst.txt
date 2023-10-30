8. Advanced Topics
==================

This is the Advanced Topics section of our Python documentation. Here, we explore how you can manipulate the simulation during it's run and we explore the
``KineticStruct`` which can be used for Michaelis Menten kinetics. Next we discuss the ``Kinetics function``` 
These tools enable you to add more biological information in the modelling.

Let's dive in!

Manipulating the simulation
---------------------------

Each simulation for objects of the ``TimeStepDynamicFBABase`` base class can be given an extra function enabling the user to manipulate the simulation on the go.
By manipulating you can think of adding or removing external metabolites from the system. But following some simple rules nearly everything can be changed during simulation.

Here we will give an example where we grow two strains of *E coli.* on Pyruvate and add L-Glutamate when the concentration of Pyruvate is half it's initial concentration.

First we build a ``DynamicJointFBA`` model:

.. code-block:: python

    import numpy
    from cbmpy.CBModel import Model
    from dcFBA.Models import CommunityModel
    from dcFBA.DynamicModels import DynamicJointFBA
    from dcFBA.DefaultModels import read_default_model

    # Build the model and set all bounds
    model1: Model = read_default_model("e_coli_core")
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

    # Set initial concentrations
    dynamic_fba = DynamicJointFBA(
        combined_model,
        [0.1, 0.1],
        {"M_pyr_e": 10, "M_glu__L_e": 0, "M_glc__D_e": 0},
    )
    
Next we will define our deviation function. Remember that the function, no matter if you use it or not, needs to accept the ``DynamicJointFBA`` model, the ``time steps`` of the simulation 
and a ``condition`` which we will explain in a minute. Furthermore, you always need to return an ``int``.

.. code-block:: python

    def deviate_func(DFBA: DynamicJointFBA, used_time, condition) -> int:
    if DFBA.metabolites["M_pyr_e"][-1] <= 5.0 and condition < 1:
        combined_model.getReaction("R_EX_glu__L_e").setLowerBound(-numpy.inf)

        DFBA.metabolites["M_glu__L_e"][-1] = 30

        return 1

    return condition


Here we define the ``deviate_func`` as follows: run the function when the concentration of `M_pyr_e` runs under 5.0, run it only once. If the concentration is 
not below or equal to 5.0 return the condition. The ``condition`` parameter is by default set to 0 at the start of the simulation and is always passed to the deviate function.
This is done such that for example you can rerun the function for consecutive time steps. The value returned by the function will always be added to the global condition variable. 
By doing so we can for example add Glucose to the system for N consecutive runs by modifying the conditional statement. 

We can now pass the function to the simulation method and plot the results:

.. code-block:: python 

    import matplotlib.pyplot as plt

    dynamic_fba.simulate(0.1, deviate=deviate_func)


    T = dynamic_fba.get_time_points()
    metabolites = dynamic_fba.get_metabolites()
    biomasses = dynamic_fba.get_biomasses()


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

All methods permit the passing of a ``KineticStruct`` object during initialization. By doing so, the upper bounds of the reactions defined in this object
are not statically set, but are calculated using the Michaelis-Menten equation. 

A ``KineticStruct`` van be defined using the reaction id followed by the limiting substrate and the `Km` and `Vmax`

.. code-block:: python

    from cbmpy.CBModel import Model
    from dcFBA.DynamicModels import DynamicSingleFBA
    from dcFBA.DefaultModels import read_default_model
    from dcFBA.Models import KineticsStruct
    import matplotlib.pyplot as plt

    model: Model = read_default_model("e_coli_core")

    # Set bounds on the glucose import
    model.getReaction("R_GLCpts").setUpperBound(10)

    initial_biomass = 0.1
    initial_concentrations = {"M_glc__D_e": 10}
    #Define kinietcStruct
    kin = KineticsStruct({"R_GLCpts": ("M_glc__D_e", 5, 10)})

    ds = DynamicSingleFBA(
        model,
        "R_BIOMASS_Ecoli_core_w_GAM",
        initial_biomass,
        initial_concentrations,
        kinetics=kin,
    )

    ds.simulate(0.15)

    T = ds.get_time_points()
    metabolites = ds.get_metabolites()
    biomass = ds.get_biomass()

    ax = plt.subplot(111)
    ax.plot(T, biomass)
    ax2 = plt.twinx(ax)
    ax2.plot(T, metabolites["M_glc__D_e"], color="r")

    ax.set_ylabel("Biomass", color="b")
    ax2.set_ylabel("Glucose", color="r")
    plt.show()


.. image:: ../_static/images/kinetics_example.png
    :width: 500px
    :align: center
    :alt: Biomass concentrations


There is one exception, un the case of ``EndPointFBA`` this is done by a MM approximation, that is using the `Km` and `Vmax` of a reaction we derive two linear lines and add these as extra constraints
to the linear model. To set these extra constraints we first need to create a ``KineticStruct`` holding all the kinetic information for the reactions 
for which we know the Kinetics.

.. math::

   \text{lower_bound} = \frac{V_{\text{max}}}{Km}

   \text{upper_bound} = \frac{\frac{V_{\text{max}}}{2}}{Km}



The Kinetics function
---------------------

.. admonition:: Under Construction
   :class: warning

   This section is currently under construction. Check back later for updates.

Unused reactions
----------------

The bottleneck for ``EndPointFBA`` is the amount of RAM required to store the joint stoichiometric matrix. We introduce two methods here to lower the total amount of RAM needed.

1. Scan unused reactions and species
""""""""""""""""""""""""""""""""""""

By initializing the ``EndPointFBA`` model, we define the medium. Based on the specified medium and its associated metabolites, certain reactions might never be activated, and some external metabolites may never be produced. To identify these inactive reactions, invoke the ``dcFBA.Helpers.ScanUnusedReactions`` method. 
This function runs consecutive FBA's , each time  setting the a reaction as an objective. Reactions that remain inactive get omitted from the initial ``CommunityModel``. 
Additionally, species and reagents that are never involved are also pruned. This ensures that these reactions aren't duplicated at each time point, effectively minimizing the total number of reactions and species.


2. Sparse Matrix
""""""""""""""""

When invoking the ``EndPointFBA.simulate()`` method, the stoichiometric matrix is constructed. Setting the `sparse` option to `True` during this call will result in the creation of a ``scipy.sparse.csr_matrix``. This structure is employed to depict a sparse matrix in the compressed sparse row (CSR) format, thereby conserving a significant amount of system RAM.

