1. Installation and Requirements
================================

Prerequisites
-------------
The dynamic-community-fba package utilizes the CBMPy library [#ref_cbmpy]_ for managing constraint-based metabolic models. Similar to 
`COBRApy` [#ref_cobra]_, `CBMPy` is a Python package designed for loading, creating, and manipulating constraint-based models. 
Furthermore, it streamlines the interaction with linear program solvers, such as `CPLEX` [#ref_cplex]_.

To install `CBMPy`, you can use the following command:

.. code-block:: bash

   pip install cbmpy

For more information and detailed documentation on using `CBMPy`, please refer to the `CBMPy GitHub repository`_ and the `CBMPy documentation`_.

Please make sure `CBMPy` is up and running correctly before preceding with the analyses provided by `dcFBA`.
See also: :ref:`cbmpy-guide` section

.. _CBMPy GitHub repository: https://github.com/SystemsBioinformatics/cbmpy
.. _CBMPy documentation: https://pythonhosted.org/cbmpy/modules_doc.html
.. _COBRApy: https://opencobra.github.io/cobrapy/#:~:text=cobrapy%20is%20a%20python%20package,io.


Dynamic Community FBA
---------------------

After the installation of `CBMPy` you can install dynamic community FBA using the following command

.. code-block:: bash

   pip install dcFBA


.. Escher
.. ----------

.. Maybe we can write some easy converter functions for known maps. To display the models
.. have to think about this



.. [#ref_cbmpy] PySCeS Constraint Based Modelling (http://cbmpy.sourceforge.net) Copyright (C) 2010-2023 Brett G. Olivier, Vrije Universiteit Amsterdam, Amsterdam, The Netherlands
.. [#ref_cobra] Ebrahim, A., Lerman, J.A., Palsson, B.O. et al. COBRApy: COnstraints-Based Reconstruction and Analysis for Python. BMC Syst Biol 7, 74 (2013). https://doi.org/10.1186/1752-0509-7-74
.. [#ref_cplex] IBM (2017) IBM ILOG CPLEX 12.7 User’s Manual (IBM ILOG CPLEX Division, Incline Village, NV).1
