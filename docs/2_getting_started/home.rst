2. Getting Started
==================

With the package a the models for *E. coli core metabolism*, *Streptococcus thermophilus* and *Lactobacillus delbrueckii*, here 
we will show you the basic operations for loading and modifying a cbmpy model. If you've worked with `cobrapy` previously
you will see `cbmpy` is not that different. If this is your first time working with GSMM's and FBA we recommend  you to start with
the extensive cbmpy tutorial (if we have completed this)

2.1. Loading a model using cbmpy
--------------------------------

To load a model and perform a simple FBA analysis on it type:

.. code-block:: python

   import cbmpy
   from cbmpy.CBModel import Model

   model: Model = cbmpy.loadModel("data/bigg_models/e_coli_core.xml")
   cbmpy.doFBA(model)
   FBAsol = model.getSolutionVector(names=True)
   FBAsol = dict(zip(FBAsol[1], FBAsol[0]))
   
Here the model refers to a ``cbmpy.CBModel`` which represents the GSMM of the loaded organism .


2.2. SBML and cobra models 
--------------------------

The ``cbmpy.load_model()`` function is designed to efficiently handle a wide range of models. It seamlessly supports the 
import of models encoded in the standardized `Systems Biology Markup Language (SBML)`_ format, as well as models exported by 
`cobrapy`. This means that you can easily work with different versions of SBML and `cobrapy` models without having to 
specify them explicitly. This flexibility simplifies the model loading process. 

.. note::
    Sometimes the conversion of exchanges, sinks or other boundary conditions are not properly set when exporting or importing 
    a `cobra` model into `cbmpy` therefore always check if these reactions are set correctly in the loaded model.*

.. _Systems Biology Markup Language (SBML): https://sbml.org/


2.3. Saving a model
-------------------
There are two ways to save a cbmpy model. The easiest way is to save your altered model to the latest version of SBML:

.. code-block:: python
    
    reaction = model.getReaction("R_EX_glc__D_e") #Get a reaction from the model
    reaction.setLowerBound(0) #Alter the reaction in the model

    cbmpy.saveModel(model, "adjusted_model.xml") #Save the new model to a XML file

If you don't want to save the model to to the latest version of SBML or you wan't to save it to a `cobrapy` model you can call one of the below functions:

.. code-block:: python
    
    cbmpy.writeCOBRASBML(...)
    cbmpy.writeFVAtoCSV(...)
    cbmpy.writeModelToExcel97(...)    
    cbmpy.writeSBML3FBCV2(...)

2.4. Reactions, Reagents and Species
------------------------------------

In `cbmpy`, the ``cbmpy.CBModel`` object forms the basis of the model. When working with the model, 
most modifications will involve manipulating this object. In the previous section, 
we demonstrated how to load the *E. coli core metabolism*. Now, let's explore some basic alterations that can be made 
to the model. For a more comprehensive understanding of the functionalities available in the ``cbmpy.CBModel`` object and 
other features of cbmpy, we recommend referring to the extensive documentation_.



.. _documentation: https://pythonhosted.org/cbmpy/modules_doc.html

Reactions
*********

To list all the reactions in the model, or list the reaction containing a certain string you can call the following functions:

.. code-block:: python 
    
    modelRxns = model.getReactionIds() #All the reactions, as a list[str]
    print(modelRxns)

    model.getReactionIds('PG') #Outputs only reactions with "PG" in their ID

Once you have identified your reaction of interest, you can easily access its key details, including the reagents, upper and lower bounds, and equation, as follows:

.. code-block:: python
    
    from cbmpy.CBModel import Reaction, Reagent, Species 
    
    reaction: Reaction = model.getReaction("R_PGK")

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


You can easily add your own defined reactions to the model using ``model.createReaction()``, if we for example want to add the 
reaction: :literal:`ATP + H2O -> ADP + Pi` we can do this with the following code:

.. code-block:: python 

    model.createReaction('ATPsink', reversible = False) # Create a new empty irreversible reaction
   
    # Add the reagents to the reaction, All metabolites already existed in the model so we did not 
    # Need to create them 
    model.createReactionReagent('ATPsink', metabolite = "M_atp_c" , coefficient = -1) 
    model.createReactionReagent('ATPsink', metabolite = "M_adp_c", coefficient =1)
    model.createReactionReagent('ATPsink', metabolite =  "M_h2o_c", coefficient = -1)
    model.createReactionReagent('ATPsink', metabolite = "M_pi_c" , coefficient = 1)


Reagents
********

Tohe ``Reagent`` class  represents a reagent within a reaction, providing essential information about its properties and characteristics. 
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

Species represent the metabolites in the system using the ``Species`` object you can easily retrieve details such as the species' molecular formula, charge, and compartment information.
Furthermore you can list the reactions in which a species is consumed or created

.. code-block:: python 

    species: Species = model.getSpecies("M_pi_c")

    species.getChemFormula() 
    species.getCharge()
    species.getCompartmentId() # Gives the id of the compartment in which the species lives
    species.isReagentOf() # Returns a list of reaction ids in which the species is present




Objective function 
******************

To perform FBA on the model you need to set an objective function. This is the reaction for which the maximal or minimal flux will 
be calculated. 
To check what the active objective function of the model is you can write: 

.. code-block:: python 

    objective_ids = model.getActiveObjectiveReactionIds() 
    #['R_BIOMASS_Ecoli_core_w_GAM']
    
    objective = model.getActiveObjective()
    objective.getOperation()
    #Maximize


If you would call the function ``cbmpy.doFBA(model)`` FBA will calculate the fluxes such that the flux through the 
reaction with id `R_BIOMASS_Ecoli_core_w_GAM` is maximal. 


2.5. Transitioning to cbmpy from cobrapy
----------------------------------------

As previously mentioned, any model build using `cobrapy` or any other toolbox for that matter can easily be opened in `cbmpy`. By 
just exporting your model of interest to either SBML format or to a cobra model you can import it as an cbmpy model.



Next we will explore how cbmpy models can be used to model the behaviors of microbial communities. 

