1. Installation and Requirements
================================

Prerequisites
-------------
The dynamic-community-fba package relies on the `cbmpy` library [#ref_cbmpy]_ for handling constraint-based metabolic models. Like 
`cobrapy`_ `cbmpy`` is a Python package that simplifies the creation, loading, and manipulation of constraint-based models, allowing interaction with LP-solvers like 'cplex' or 'GLPK'.

To install `cbmpy`, you can use the following command:

.. code-block:: bash

   pip install cbmpy

For more information and detailed documentation on using `cbmpy`, please refer to the `cbmpy GitHub repository`_ and the `cbmpy documentation`_.

.. _cbmpy GitHub repository: https://github.com/SystemsBioinformatics/cbmpy
.. _cbmpy documentation: https://pythonhosted.org/cbmpy/modules_doc.html
.. _cobrapy: https://opencobra.github.io/cobrapy/#:~:text=cobrapy%20is%20a%20python%20package,io.

.. [#ref_cbmpy] PySCeS Constraint Based Modelling (http://cbmpy.sourceforge.net) Copyright (C) 2010-2023 Brett G. Olivier, Vrije Universiteit Amsterdam, Amsterdam, The Netherlands


Dynamic Community FBA
---------------------

After the installation of `cbmpy` you can install dynamic community FBA using the following command

.. code-block:: bash

   pip install NAME

More text if needed about the installation

Escher
----------

Maybe we can write some easy converter functions for known maps. To display the models
have to think about this