DART_LAB Tutorial
=================

Overview
--------

The files in this directory contain PDF tutorial materials on DART, and Matlab exercises. See below for links to the PDF
files and a list of the corresponding matlab scripts.

This tutorial begins at a more introductory level than the materials in the tutorial directory, and includes hands-on
exercises at several points. In a workshop setting, these materials and exercises took about 1.5 days to complete.

DART tutorial presentations
---------------------------

Here are the PDF files for the presentation part of the tutorial:

- :download:`Section 1: <presentation/DART_LAB_Section01.pdf>` Ensemble Data Assimilation Concepts in 1D.
- :download:`Section 2: <presentation/DART_LAB_Section02.pdf>` How Should Observations Impact an Unobserved 
  State Variable? Multivariate Assimilation.
- :download:`Section 3: <presentation/DART_LAB_Section03.pdf>` Inflation and Localization to Improve Performance.
- :download:`Section 4: <presentation/DART_LAB_Section04.pdf>` Nonlinear and Non-Gaussian Extensions.
- :download:`Section 5: <presentation/DART_LAB_Section05.pdf>` Adaptive Inflation.

Matlab hands-on exercises
-------------------------

In the ``matlab`` subdirectory are a set of Matlab scripts and GUI (graphical user interface) programs which are
exercises that go with the tutorial. Each is interactive with settings that can be changed and rerun to explore various
options. A valid `Matlab <http://www.mathworks.com/products/matlab/>`__ license is needed to run these scripts.

The exercises use the following functions:

-  bounded_oned_ensemble
-  gaussian_product
-  oned_cycle
-  oned_ensemble
-  oned_model
-  oned_model_inf
-  run_lorenz_63
-  run_lorenz_96
-  run_lorenz_96_inf
-  twod_ensemble
-  twod_ppi_ensemble

To run these, cd into the DART_LAB/matlab directory, start matlab, and type the names at the prompt.
