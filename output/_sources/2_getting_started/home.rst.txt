2. CBMPy Quick Start Guide
==========================

Before delving into the modeling of microbial communities, let's establish a solid foundation by exploring the core principles of the `CBMPy` package. 
We will guide you through a series of instructive examples that illuminate key operations for GSMM's. These examples will demonstrate the foundational techniques 
for loading and manipulating models. For these examples we will use the *E. coli core iJR904* which comes with the `CBMPy` package. 
Through these hands-on demonstrations, you'll swiftly grasp the fundamental operations that underpin effective model manipulation.
For a more comprehensive understanding and explanation of all the functionalities available in  `CBMPy`, we recommend you to read the extensive documentation_.

For those who have prior experience with `CPBRApy`, the transition to `CBMPy` will feel remarkably intuitive, as the two share notable similarities in their approach and methodology. 

.. _documentation: https://pythonhosted.org/cbmpy/modules_doc.html


2.1. Loading a model using CBMPy
--------------------------------

To load a model and perform a simple FBA analysis on it type:

.. code-block:: python

   import cbmpy
   from cbmpy.CBModel import Model

   iJR904: Model = cbmpy.readSBML3FBC("cbmpy_test_ecoli") #Load the model
   solution = cbmpy.doFBA(iJR904) #Perform FBA, returns the objective value
   FBAsol = iJR904.getSolutionVector(names=True) #Get all reaction ids with there flux value
   
Here model *iJR904* refers to a ``cbmpy.CBModel`` which represents the GSMM of the loaded organism .


2.2. SBML and COBRA models 
--------------------------

The ``cbmpy.loadModel()`` function is designed to efficiently handle a wide range of models. It seamlessly supports the 
import of models encoded in the standardized `Systems Biology Markup Language (SBML)`_ format, as well as models exported by 
`COBRApy`. This means that you can easily work with different versions of SBML and `COBRApy` models without having to 
specify them explicitly. This flexibility simplifies the model loading process. 

.. note::
    Sometimes the conversion of exchanges, sinks or other boundary conditions are not properly set when exporting or importing 
    a `COBRA` model into `CBMPy` therefore always check if these reactions are set correctly in the loaded model.

.. _Systems Biology Markup Language (SBML): https://sbml.org/


2.3. Saving a model
-------------------
There are two ways to save a cbmpy model. The easiest way is to save your altered model to the latest version of SBML:

.. code-block:: python
    
    from cbmpy.CBModel import Reaction 

    reaction: Reaction = iJR904.getReaction("R_EX_glc__D_e") #Get a reaction from the model
    print(reaction.getLowerBound()) #-10.0
    reaction.setLowerBound(0) #Alter the reaction in the model
    print(reaction.getLowerBound()) #0.0

    cbmpy.saveModel(iJR904, "adjusted_model.xml") #Save the new model to a XML file

You can save the modified model in a file format of your choice using any of the following methods:

.. code-block:: python
    
    cbmpy.writeCOBRASBML(...)
    cbmpy.writeFVAtoCSV(...)
    cbmpy.writeModelToExcel97(...)    
    cbmpy.writeSBML3FBCV2(...)

2.4. Reactions, Reagents and Species
------------------------------------

In `CBMPy`, the ``cbmpy.CBModel`` object stores the model and all it's attributes. When working with the model, 
most modifications will involve manipulating this object. In the previous section, 
we demonstrated how to load the *E. coli core iJR904* model and perform a FBA on it. Now, let's explore some basic alterations that can be made 
to the model.

Reactions
*********

To list all the reactions in the model, or list the reaction containing a certain string you can call the following functions:

.. code-block:: python 
    
    modelRxns = iJR904.getReactionIds() #All the reactions, as a list[str]
    print(modelRxns)

    print(iJR904.getReactionIds('PG'))  #Outputs only reactions with "PG" in their ID

Once you have identified your reaction of interest, you can easily access its key details, including the reagents, upper and lower bounds, and equation, as follows:

.. code-block:: python
    
    from cbmpy.CBModel import Reagent, Species 
    
    reaction: Reaction = iJR904.getReaction("R_PGK")

    reagents: list[Reagent] = reaction.getReagentObjIds()  # Get all reagent ids of the reaction
    print(reagents)

    bounds = [reaction.getLowerBound(), reaction.getUpperBound()] # Get the lower and upper bound
    print(bounds)

    equation = reaction.getEquation() # Get the reactions equation
    print(equation)

Furthermore you can check if a reaction is reversible and if it is an exchange reaction:

.. code-block:: python
    
    print(reaction.is_exchange) #True if the reaction is an exchange reaction

    print(reaction.reversible) #True if the reaction is reversible


You can easily add your own defined reactions to the model using the ``createReaction()`` method, if we for example want to add the 
irreversible reaction: :literal:`ATP + H2O -> ADP + Pi` we can do this with the following code:

.. code-block:: python 

    iJR904.createReaction('ATPsink', reversible = False) # Create a new empty irreversible reaction
   
    # Add the reagents to the reaction, All metabolites already existed in the model so we did not 
    # Need to create them 
    iJR904.createReactionReagent('ATPsink', metabolite = "M_atp_c" , coefficient = -1) 
    iJR904.createReactionReagent('ATPsink', metabolite = "M_adp_c", coefficient =1)
    iJR904.createReactionReagent('ATPsink', metabolite =  "M_h2o_c", coefficient = -1)
    iJR904.createReactionReagent('ATPsink', metabolite = "M_pi_c" , coefficient = 1)


Reagents
********

The ``Reagent`` class  represents a reagent within a reaction, providing essential information about its properties and characteristics. 
Within the class, users can access and manipulate the reagents associated with a specific reaction within the model. The reagent itself 
is linked to a ``Species`` which we will cover shortly. 
You can access a reagent by retrieving it from an instance of the ``Reaction`` class, given the `R_PGK` reaction from the previous example
we can access information about a reagent as follows:

.. code-block:: python

    reagent: Reagent = reaction.getReagent("R_PGK_M_3pg_c")

    reagent.getCoefficient() # Get the reagent's stoichiometric coefficient

    reagent.getCompartmentId() #Get the compartment 

    reagent.getSpecies() # Get the species id corresponding to this reagent 

If a reagent has a negative coefficient it is consumed by the reaction, if the reagent has a positive coefficient it is created by the reaction.

Species
*******

Species represent the metabolites in the system using the ``Species`` object you can easily retrieve details such as the molecular formula, charge, and the compartment of the species.
Furthermore you can list the reactions in which a species is consumed or synthesized

.. code-block:: python 

    species: Species = iJR904.getSpecies("M_pi_c")

    species.getChemFormula() 
    species.getCharge()
    species.getCompartmentId() # Gives the id of the compartment in which the species lives
    species.isReagentOf() # Returns a list of reaction ids in which the species is present




Objective function 
******************

To perform FBA on the model you need to set an objective function. The output of FBA 
will be a flux distribution which minimizes/maximizes this objective function. 

To check what the active objective function of the model is you can write: 

.. code-block:: python 

    objective_ids = iJR904.getActiveObjectiveReactionIds() #Returns the IDs of the reactions which have been set as objective reaction
    
    print(objective_ids)
    #['R_BIOMASS_Ecoli']
    
    objective = iJR904.getActiveObjective()
    print(objective.getOperation())
    #Maximize

    reaction: Reaction = iJR904.getReaction("R_EX_glc__D_e") 
    reaction.setLowerBound(-10) #Reset lower bound
    solution = cbmpy.doFBA(iJR904) #0.922


Calling ``cbmpy.doFBA(iJR904)``  will calculate the fluxes such that the flux through the 
reaction with id `R_BIOMASS_Ecoli` is maximized. 

Next, we'll delve into dynamic modeling of CBMPy models. Once we lay this foundation, we'll journey into the fascinating realm of modeling microbial communities.
