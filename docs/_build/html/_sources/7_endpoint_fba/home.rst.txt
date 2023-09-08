7. EndPointFBA 
==============

The final method implemented in this package is called `EndPointFBA`. 
EndPointFBA is a dynamic simulation method for microbial 
communities provided by this package. 

Unlike other methods provided by the package that run FBA foreach time point, EndPointFBA adopts a distinct approach.
It generates a new `CommunityModel` for each time point. 
These individual models are then interconnected, allowing for the forward flow of metabolites and biomass from one time point to the next without 
the possibility of reverse flow. With this structure in place, a singular FBA can be executed across the interconnected models.

7.1 Example
-----------

We'll use a lightweight toy model for this example to ensure the 
computational load remains manageable. 

The toy model is bundled with the package:

.. code-block:: python 
    
    from DCFBA.ToyModels import model_a, model_b

    m_a: Model = model_a.build_toy_model_fba_A()

    #Setting some initial bounds on reactions
    m_a.getReaction("R_1").setUpperBound(10) 
    m_a.getReaction("R_4").setUpperBound(3)
    m_a.getReaction("R_6").setUpperBound(1)

    m_a.getReaction("R_1").setLowerBound(0)
    m_a.getReaction("R_4").setLowerBound(0)
    m_a.getReaction("R_6").setLowerBound(0)

    # Delete the original biomass model from model A as it is redefined in the endPoint model
    reaction: Reaction = m_a.getReaction("R_BM_A").deleteReagentWithSpeciesRef(
        "BM_e_A"
    ) 

    m_b: Model = model_b.build_toy_model_fba_B()
    
    #Setting some initial bounds on reactions
    m_b.getReaction("R_1").setUpperBound(10)
    m_b.getReaction("R_3").setUpperBound(1)
    m_b.getReaction("R_5").setUpperBound(1)
    
    m_b.getReaction("R_1").setLowerBound(0)
    m_b.getReaction("R_3").setLowerBound(0)
    m_b.getReaction("R_5").setLowerBound(0)

    # Delete the original biomass model from model B
    reaction: Reaction = m_b.getReaction("R_BM_B").deleteReagentWithSpeciesRef(
        "BM_e_B"
    ) 


Next, let's set the initial ``CommunityModel``:

.. code-block:: python 

    from DCFBA.Models import CommunityModel

    community_model = CommunityModel(
        [m_a, m_b], ["R_BM_A", "R_BM_B"], ["modelA", "modelB"]
    )

To run a simulation with ``EndPointFBA``, you'll need to initialize the object by providing the following:

- The ``CommunityModel``
- The desired number of time points
- Initial concentrations for both biomass and metabolites
- The time step size


.. code-block:: python

    from DCFBA.DynamicModels import EndPointFBA

    n = 25
    ep = EndPointFBA(
        community_model,
        n,
        {"modelA": 1.0, "modelB": 2.0},
        {"S_e": 100, "A_e": 0.0, "B_e": 0.0},
        0.1,
    )


    solution = ep.simulate()
    print(solution) 
    #12.778

This provides the community biomass value after 25 intervals of 0.1 time units each.

7.2 Examining the results
-------------------------

From the ``EndPointFBA`` object we can retrieve the ``CommunityModel``, which we can examine as any other ``cbmpy.CBModel`` object.

For instance, to view biomass values of model A at the start of 
the 10th interval and at the last time point:

.. code-block:: python 

    FBAsol = ep.m_model.getSolutionVector(names=True)
    FBAsol = dict(zip(FBAsol[1], FBAsol[0]))

    print(FBAsol["BM_modelA_time09_time10"])
    #1.629

    print(FBAsol["BM_modelA_exchange_final"])
    #2.667

Reactions in the ``EndPointFBA`` model have unique IDs formed by 
appending their original IDs with the time ID (e.g., R_1_timeNN). 
Also, metabolite flow between time points are represented with a 
unique ID combining the metabolite's ID and the start and end 
interval IDs.

Built-in plotting functions allow for visualizing metabolite flows and biomass concentrations over time:

.. code-block:: python 
    
    from DCFBA.Helpers.PlotsEndPointFBA import plot_biomasses, plot_metabolites

    plot_biomasses(ep)
    plot_metabolites(ep, {"S_e": 100, "A_e": 0.0, "B_e": 0.0}) #The method requires that you give the ids of the metabolites you want to plot 



7.3 Advanced constraints
------------------------

EndPointFBA offers additional methods for refining your simulation results, enhancing the biological relevance. Here, we discuss these methods and demonstrate their application.

7.3.1 Optimal time search
"""""""""""""""""""""""""

The ``OptimalTimeSearch`` module helps determine the optimal number of intervals required to reach the maximum objective method value. By default, the chosen number of intervals might either stretch resources too thinly or not utilize them fully, leading to skewed FBA results.

.. code-block:: python
   
   from DCFBA.Helpers.OptimalTimeSearch import search

   optimal_timepoints = search(
       community_model,
       {"modelA": 1.0, "modelB": 2.0},
       {"S_e": 100, "A_e": 0.0, "B_e": 0.0},
       0.1,
   )

   print(optimal_timepoints)  # e.g., 21

.. warning::
   Keep in mind, this search might be time-consuming as it involves iterative model construction and simulation.

Use the determined optimal number of time points to create a new ``EndPointFBA`` 
model and examine the results of the optimal solution.

7.3.2 Restricting reaction bounds
"""""""""""""""""""""""""""""""""

In the context of steady-state microbial community dynamics, drastic changes in reaction rates between two successive time points are unlikely. 
To capture this behavior more accurately, we have implemented two methods to minimize these fluctuations.

The first approach introduces a constraint on the model to ensure that the reaction rate at time ``N`` doesn't vary more than a predefined value, 
:math:`\epsilon`, from its rate at time ``N-1``. 
This method effectively smoothens the fluctuations 
in reaction rates, leading to more biologically plausible results.

.. code-block:: python

    # Activate the additional constraints with epsilon set to 0.1
    ep.constrain_rates(epsilon=0.1) 

    solution = ep.simulate()
    print(solution) # 12.599
    
    FBAsol = ep.m_model.getSolutionVector(names=True)
    FBAsol = dict(zip(FBAsol[1], FBAsol[0]))

     # Print the rates for a particular reaction at successive time points
    print(FBAsol["R_1_modelA_time01"]) #1.062
    print(FBAsol["R_1_modelA_time02"]) #1.128
    print(FBAsol["R_1_modelA_time03"]) #1.975

As demonstrated, the variation in reaction rates between consecutive time points remains within the :math:`\epsilon` boundary of 0.1.

A more sophisticated approach to minimizing fluctuations in reaction rates involves quadratic programming. The idea is to set the `EndPointFBA` model's objective to minimize the sum of squared differences between reaction rates at successive time points. Here's how you can employ this approach:

1. First, simulate the `EndPointFBA` model to obtain the optimal solution value.
2. Use the obtained optimal solution value as a constraint for the `EndPointFBA` model.
3. Set the new objective function to minimize the sum of squared differences between reaction rates at two consecutive time points.

The resulting objective function is expressed as:

.. math::
   \min \sum_{j, n} (r_{j,n} - r_{j,n+1})^2

Where r\ :sub:`j,n`\ denotes the rate of reaction ``j`` at time point ``n``. 

The following code demonstrates how to implement this approach:

.. code-block:: python

   # Set the objective value of community biomass to 12.77
   ep.set_qp(12.77)

   # Simulate the model with the new quadratic objective
   ep.simulate()

Remember, this approach is useful when you want the model to focus on smooth transitions between time points, rather than maximizing any particular objective.

7.3.3 Michaelis-Menten approximation
""""""""""""""""""""""""""""""""""""

The last method to refine the model is by making use of the famous Michaelis-Menten curve. Since deriving the MM curve foreach reaction
is computational infeasible. We can still make use of any kinetic knowledge that we have for a reaction.

Here we make use of the MM approximation, that is using the `Km` and `Vmax` of a reaction we derive two linear lines and add these as extra constraints
to the linear model. To set these extra constraints we first need to create a ``KineticStruct`` holding all the kinetic information for the reactions 
for which we know the Kinetics. The ``KineticStruct`` is created using a dictionary with reaction ids followed by a list of restricting metabolite id, 
Km and Vmax of the reaction.

.. code-block:: python

    from DCFBA.Models import KineticsStruct
    kin = KineticsStruct(
        {
            "R_1_modelA": ["S_e", 10, 10],
            "R_4_modelA": ["B_e", 5, 3],
            "R_6_modelA": ["B_e", 3, 1],
            # B
            "R_1_modelB": ["S_e", 10, 10],
            "R_3_modelB": ["A_e", 2, 1],
        }
    )

Using the kinetic struct we can set the MM approximation as follows:

.. code-block:: python

    n = 21
    ep = EndPointFBA(
        community_model,
        n,
        {"modelA": 1.0, "modelB": 2.0},
        {"S_e": 100, "A_e": 0.0, "B_e": 0.0},
        kin,
        0.1,
    )

    #Activate the extra constraint
    ep.mm_approximation()
    
    solution = ep.simulate()
    print(solution)



Happy modelling!
