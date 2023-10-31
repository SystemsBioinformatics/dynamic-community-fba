5. Dynamic Joint FBA 
====================

Dynamic Joint Flux Balance Analysis (djFBA) is a method to dynamically model microbial communities. The core principle of djFBA is the formation of 
`community biomass`, denoted as ``X_c``. Which is introduced in the ``CommunityModel`` to be the total sum of individual biomasses. The formation of 
this mock species is than set to be the objective function for FBA.

Using this in combination with the static optimization approaches as introduced in the dFBA section. We now have our first method to model microbial communities!

** citation from Dynamic Joint FBA

Example
--------

Here we will give use the example of a combining the  *E. coli core metabolism* model with the *Streptococcus thermophilus* (iRZ476) model.
To perform the Dynamic Joint FBA we first define the ``CommunityModel`` and use this to initialize a ``DynamicJointFBA`` object:

.. code-block:: python

    import matplotlib.pyplot as plt
    from cbmpy.CBModel import Model
    from dcFBA.Models import CommunityModel
    from dcFBA.DynamicModels import DynamicJointFBA
    from dcFBA.DefaultModels import read_default_model

    model1: Model = read_default_model("e_coli_core")
    model2: Model = read_default_model("strep_therm")

    # Set the import bounds for glucose in both models
    model1.getReaction("R_GLCpts").setUpperBound(10)
    model2.getReaction("R_GLCpts").setUpperBound(6)

    # The biomass reactions ids
    biomass_reaction_model_1: str = "R_BIOMASS_Ecoli_core_w_GAM"
    biomass_reaction_model_2: str = "R_biomass_STR"

    cm = CommunityModel(
        [model1, model2],
        [biomass_reaction_model_1, biomass_reaction_model_2],
        ["ecoli", "strep"],
    )  # Define the community model

    dynamic_joint_fba = DynamicJointFBA(
        cm,
        [1.0, 1.0],
        {
            "M_glc__D_e": 100,
            "M_succ_e": 0,
            "M_glu__L_e": 0.0,
            "M_gln__L_e": 0.0,
            "M_lcts_e": 100,
        },
    )  # Create a DynamicJointFBA object, set the initial concentrations of glucose and lactose to 100


Now that the model is in place we can run the simulation using time steps of 0.1:

.. code-block:: python

    dynamic_joint_fba.simulate(0.1)

The simulate method conducts the simulation and stores all values of biomass concentrations, metabolite concentrations, and fluxes at each time point. You can access these values by retrieving the respective objects.

.. code-block:: python

    biomasses = dynamic_joint_fba.get_biomasses()
    metabolites = dynamic_joint_fba.get_metabolites()
    time_points = dynamic_joint_fba.get_time_points()
    fluxes = dynamic_joint_fba.get_fluxes()

Since the fluxes are stored are the aggregated fluxes (the flux multiplied by the total amount of biomass present at that time), we can also restore the specific flux for a reaction:

.. code-block:: python

    specific_flux_values = dynamic_joint_fba.get_specific_flux_values("R_GLCpts_ecoli")


You can now easily plot the species concentration over time:

.. code-block:: python

    plt.plot(
        time_points, metabolites["M_glc__D_e"], color="blue", label="[Glucose]"
    )
    plt.plot(
        time_points, metabolites["M_lcts_e"], color="orange", label="[Lactose]"
    )

    plt.xlabel("Time")
    plt.ylabel("Concentration")
    plt.legend()
    plt.show()
  


.. image:: ../_static/images/Metabolites_DJFBA.png
    :width: 500px
    :align: center
    :alt: Biomass concentrations
     
And the biomasses of both species over time

.. code-block:: python
    
    plt.plot(time_points, biomasses["ecoli"], color="orange", label="ecoli")
    plt.plot(time_points, biomasses["strep"], color="blue", label="strep")

    plt.xlabel("Time")
    plt.ylabel("Concentration")
    plt.legend()
    plt.show()


.. image:: ../_static/images/Biomass_DJFBA.png
    :width: 500px
    :align: center
    :alt: Biomass concentrations
     
