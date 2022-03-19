Current Flow Allostery
============
Documentation: https://current-flow-allostery.readthedocs.io/en/latest/index.html

[//]: # (Badges)
[![GitHub Actions Build Status](https://github.com/REPLACE_WITH_OWNER_ACCOUNT/current_flow_allostery/workflows/CI/badge.svg)](https://github.com/REPLACE_WITH_OWNER_ACCOUNT/current_flow_allostery/actions?query=workflow%3ACI)
[![codecov](https://codecov.io/gh/REPLACE_WITH_OWNER_ACCOUNT/current_flow_allostery/branch/master/graph/badge.svg)](https://codecov.io/gh/REPLACE_WITH_OWNER_ACCOUNT/current_flow_allostery/branch/master)

![Robust Determination of Protein Allosteric Signaling Pathways](./pics/Robust_Determination_of_Protein_Allosteric_Signaling_Pathways.png)

Allosteric network from molecular dynamics (MD) simulations is powerful for understanding how protein function changes upon an allosteric perturbation, such as ligand binding and mutation. The main challenge in mapping out the information propagating between amino acids is due to the large fluctuations in protein dynamics that cause instability of the network topology. To solve this problem, we introduce the current-flow betweenness scheme, which originated from electrical network theory. The current-flow betweenness provides a significant improvement in the convergence of the allosteric networks.

https://pubs.acs.org/doi/10.1021/acs.jctc.8b01197

Install
=======
As the current module is still under development and clean up, the best way to get started is to clone the repository to your local computer.

This clones the repository:
   'git clone git@github.com:charliewen95/current_flow_allostery.git'

Then you have 2 choices:
   1. Create environemnt using "./creating_environment.sh"
   2. Install the python package using "pip install -e ."
      Note: if the code aboe doesn't work, try "python -m pip install -e ."
   3. Now you can try running the code.

Usage 
=====
This package is written in 2 ways:
   1. Specific for bash submission
   2. As a python function

e.g.
import current_flow_allostery
current_flow_allostery.betweenness()

For further and more specific examples of how to use, check out the attached bash submission script. Within you will be able to identify 4 different types of submission:
   1. Parallel bash submission
   2. Parallel python function submission
   3. Single run bash submission
   4. Single run python function submission

Contributors and Acknowledgement
================================
Author : Botello-Smith, W., Luo, Y.

Code clean-up and packaging : Chen-yun Wen

Citation and References
=======================
1. Botello-Smith, W., Luo, Y., Concepts, practices, and interactive tutorial for allosteric network analysis of molecular dynamics simulations. Methods Mol. Biol. 2021, 2302, 311
9:47

2. Botello-Smith, W., Luo, Y.*, Robust determination of protein allosteric signaling pathways, J. Chem. Theory Comput., 2019, 15, 2116-2126.
9:47

3. Botello-Smith, W, Luo, Y.*, Investigating protein-protein allosteric network using current-flow scheme, J. Comput. Chem, 2020, 41, 552-560

Question/Suggestion?
====================
   1. Try to reach Charlie at chenyun.wen@westernu.edu


### Copyright

Copyright (c) 2022, Chen Yun Wen


#### Acknowledgements
 
Project based on the 
[Computational Molecular Science Python Cookiecutter](https://github.com/molssi/cookiecutter-cms) version 1.6.


