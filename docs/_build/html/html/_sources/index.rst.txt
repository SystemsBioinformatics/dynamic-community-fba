.. dynamic-community-fba documentation master file, created by
   sphinx-quickstart on Mon Jul  3 11:45:14 2023.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

Dynamic-community-fba's documentation
================================================

Welcome to the documentation of dynamic-community-fba. The documentation provides you with a comprehensive guidance on using the dynamic-community-fba Python package. The package was designed
to conduct dynamic flux balance analysis (dFBA) in microbial communities and to easily analyze the results. dcFBA aims to accurately model and study the collaborative and 
competitive behaviors of organisms in a shared environment or community. This approach extends traditional FBA by considering the metabolic exchange and cross-feeding between 
the members of the consortia, offering a holistic perspective on microbial ecosystems. By utilizing the Genome Scale Metabolic Models (GSMM's) of organisms of interest, our package  
enables you to analyze and explore the intricate interactions among two or more organisms.

The documentation starts off with a quick introduction to `CBMPy`, providing an overview of the fundamentals for working with GSMMs in Python. 
With this knowledge we will first give an example on how to dynamically model single organisms using Dynamic FBA (dFBA). 
Ultimately, we will explore the three distinct approaches offered by this package for the dynamic modeling microbial communities:

- Dynamic Parallel FBA (dpFBA) [#ref_dynamic]_ [#ref_djfba]_ 
- Dynamic Joint FBA (djFBA) [#ref_dynamic]_ [#ref_dpfba]_ 
- EndPointFBA (epFBA) 
  
One of the key features of dynamic-community-fba is the ability to construct and export a community model, which represents the joint stoichiometry matrix of the provided 
GSMM's models. The documentation provides in-depth explanations on how to build and utilize this matrix to study the dynamics of multi-organism interactions and metabolic 
networks. Furthermore, we describe how you can export the Community Model to the standardized SBML format. 

To ensure a smooth start, the documentation includes a section outlining the prerequisites and installation guide. We describe step-by-step instructions on 
installing the necessary dependencies and setting up the dynamic-community-fba package in their Python environment. Additionally, the documentation highlights the 
compatibility requirements and recommends best practices for a successful installation. By following the documentation, you will gain a comprehensive understanding of 
the dynamic-community-fba package and its capabilities. 

The dynamic-community-fba documentation serves as a valuable resource for researchers, scientists, and students working in the field of systems biology, 
metabolic engineering, and microbial ecology. It empowers effective modelling and analysis of dynamic flux balance analyses, facilitating a deeper understanding of the 
complex interactions between organisms in various biological systems.


.. toctree::
   :maxdepth: 2
   :hidden:

   1_installation/home
   2_getting_started/home
   3_DynamicFBA/home
   4_community_matrix/home
   5_djoint/home
   6_parallel_dfba/home
   7_endpoint_fba/home   
   8_Advanced/home
   9_faq/home
   10_API/DCFBA

.. [#ref_dynamic] Mahadevan R, Edwards JS, Doyle FJ 3rd. Dynamic flux balance analysis of diauxic growth in Escherichia coli. Biophys J. 2002 Sep;83(3):1331-40. doi: 10.1016/S0006-3495(02)73903-9. PMID: 12202358; PMCID: PMC1302231.
.. [#ref_djfba] Stolyar, Sergey, Van Dien, Steve, Hillesland, Kristina Linnea, Pinel, Nicolas, Lie, Thomas J, Leigh, John A, Stahl, David A, (2007) Metabolic modeling of a mutualistic microbial community. Molecular Systems Biology, 3. 92. doi: accession:10.1038/msb4100131
.. [#ref_dpfba] Tzamali, E., Poirazi, P., Tollis, I.G. et al. A computational exploration of bacterial metabolic diversity identifying metabolic interactions and growth-efficient strain communities. BMC Syst Biol 5, 167 (2011). https://doi.org/10.1186/1752-0509-5-167