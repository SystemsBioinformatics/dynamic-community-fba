.. dynamic-community-fba documentation master file, created by
   sphinx-quickstart on Mon Jul  3 11:45:14 2023.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.
.. meta::
    :description lang=en:
         Documentation for dynamic-community-fba, a Python package that facilitates dynamic flux balance analysis (dFBA) 
         of microbial communities, modeling metabolic exchanges and interactions. 
         Perfect for researchers in systems biology, metabolic engineering, and microbial ecology.

    :Keywords: dynamic-community-fba, Python, dynamic flux balance analysis, microbial communities, Genome Scale Metabolic Models (GSMM), Systems biology

Dynamic-community-fba's documentation
================================================

Welcome to the documentation of `dynamic-community-fba`. The documentation provides you with a comprehensive guidance on using the `dynamic-community-fba` Python package. The package was designed
to conduct dynamic flux balance analysis (dFBA) in microbial communities and to easily analyze the results. DcFBA aims to accurately model and study the collaborative and 
competitive behaviors of organisms in a shared environment or community. The extended FBA methodologies incorporate 
dynamic elements by simulating metabolic exchanges and cross-feeding in microbial consortia over time. 
By integrating a temporal dimension, these approaches offer a more nuanced view of microbial ecosystems, capturing the evolving interactions and flux distributions that traditional FBA overlooks. 
This enhancement allows for a deeper understanding of the dynamic nature of microbial communities.
By utilizing the Genome Scale Metabolic Models (GSMM's) of organisms of interest, our package  
enables you to analyze and explore the intricate interactions among two or more organisms.

The documentation starts off with a quick introduction to `CBMPy`, providing an overview of the fundamentals for working with GSMMs in Python. 
With this knowledge we will first give an example on how to dynamically model single organisms using Dynamic FBA (dFBA). 
Ultimately, we will explore the three distinct approaches offered by this package for the dynamic modeling microbial communities:

- Dynamic Parallel FBA (dpFBA) [1]_ [2]_
- Dynamic Joint FBA (djFBA) [1]_ [3]_
- EndPointFBA (epFBA)
  
One of the key features of `dynamic-community-fba` is the ability to construct and export a community model, which represents the joint stoichiometry matrix of the provided 
GSMM's models. The documentation provides in-depth explanations on how to build and utilize this matrix to study the dynamics of multi-organism interactions and metabolic 
networks. Furthermore, we describe how you can export the Community Model to the standardized SBML format. 

To ensure a smooth start, the documentation includes a section outlining the prerequisites and installation guide. We describe step-by-step instructions on 
installing the necessary dependencies and setting up the `dynamic-community-fba` package in their Python environment. Additionally, the documentation highlights the 
compatibility requirements and recommends best practices for a successful installation. By following the documentation, you will gain a comprehensive understanding of 
the `dynamic-community-fba` package and its capabilities. 

This documentation serves as a valuable resource for researchers, scientists, and students working in the field of systems biology, 
metabolic engineering, and microbial ecology. It empowers effective modelling and analysis of dynamic flux balance analyses, facilitating a deeper understanding of the 
complex interactions between organisms in various biological systems.

This documentation is designed to cater to researchers, scientists, and students who are actively involved in the domains of systems biology, metabolic engineering, and microbial ecology. It serves as a comprehensive resource to kickstart your usage of the package.
The primary goal of this package is to streamline the process of dynamic flux balance analysis, enabling more efficient modeling and analysis. 
Thereby we hope it will help in enhancing our comprehension of intricate interactions between organisms within diverse biological systems.
We warmly invite you to delve into this documentation and uncover the ways it can empower your research and studies within these fields.


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

.. [1] Mahadevan R, Edwards JS, Doyle FJ 3rd. Dynamic flux balance analysis of diauxic growth in Escherichia coli. Biophys J. 2002 Sep;83(3):1331-40. doi: 10.1016/S0006-3495(02)73903-9. PMID: 12202358; PMCID: PMC1302231.
.. [2] Stolyar, Sergey, Van Dien, Steve, Hillesland, Kristina Linnea, Pinel, Nicolas, Lie, Thomas J, Leigh, John A, Stahl, David A, (2007) Metabolic modeling of a mutualistic microbial community. Molecular Systems Biology, 3. 92. doi: accession:10.1038/msb4100131
.. [3] Tzamali, E., Poirazi, P., Tollis, I.G. et al. A computational exploration of bacterial metabolic diversity identifying metabolic interactions and growth-efficient strain communities. BMC Syst Biol 5, 167 (2011). https://doi.org/10.1186/1752-0509-5-167
