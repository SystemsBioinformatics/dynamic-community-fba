.. dynamic-community-fba documentation master file, created by
   sphinx-quickstart on Mon Jul  3 11:45:14 2023.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

Dynamic-community-fba's documentation
================================================

Welcome to the documentation of dynamic-community-fba. The documentation provides you with a comprehensive guidance on using the dynamic-community-fba Python package for 
modeling dynamic flux balance analysis (FBA) in microbial consortia. Community Flux Balance Analysis (Community FBA) aims to accurately model and study the collaborative and 
competitive behaviors of multiple organisms in a shared environment. This approach extends traditional FBA by considering the metabolic exchange and cross-feeding between 
community members, offering a holistic perspective on microbial ecosystems. By utilizing the Genome Scale Metabolic Models (GSMM’s) of organisms of interest, our package  
enables you to analyze and explore the intricate interactions among two or more organisms.

The documentation starts off with a quick introduction on how to dynamically model single organisms using Dynamic FBA (dFBA). 
Next, we will cover the three different approaches for modelling microbial communities:

- Parallel FBA [reference to paper]
- Joint dynamic FBA [refference to paper]
- EndPointFBA [refference to paper]

One of the key features of dynamic-community-fba is the ability to construct and export a community model, which represents the joint stoichiometry matrix of the provided 
GSMM’s models. The documentation provides in-depth explanations on how to build and utilize this matrix to study the dynamics of multi-organism interactions and metabolic 
networks. Furthermore, we describe how you can export the Community Model to the standardized SBML format. 

To ensure a smooth start, the documentation includes a section outlining the prerequisites and installation guide. We describe step-by-step instructions on 
installing the necessary dependencies and setting up the dynamic-community-fba package in their Python environment. Additionally, the documentation highlights the 
compatibility requirements and recommends best practices for a successful installation. By following the documentation, you will gain a comprehensive understanding of 
the dynamic-community-fba package and its capabilities. They will learn how to leverage its functionalities to perform Parallel FBA, joint FBA, and endPointFBA analyses, 
enabling them to study the dynamic behavior of metabolic networks in diverse biological systems.

The dynamic-community-fba documentation serves as a valuable resource for researchers, scientists, and students working in the field of systems biology, 
metabolic engineering, and microbial ecology. It empowers you to effectively model and analyze dynamic flux balance analyses, facilitating a deeper understanding of the 
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

