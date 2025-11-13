.. _DART_LAB:

DART_LAB Tutorial
=================

DART_LAB is a MATLAB®-based tutorial that introduces the fundamental concepts of 
ensemble data assimilation through a combination of PDF slides and hands-on MATLAB® 
exercises. The full tutorial is designed to take about 1.5 days in a workshop setting.

If you are already familiar with ensemble DA, you can skip straight to the :ref:`DART-tutorial`.

DART_LAB tutorial slides
------------------------

1. `Ensemble Data Assimilation Concepts in 1D. <https://ncar.github.io/dart-tutorial/DART_LAB_Section01.pdf>`_

2. `How Should Observations Impact an Unobserved State Variable? Multivariate Assimilation. <https://ncar.github.io/dart-tutorial/DART_LAB_Section02.pdf>`_

3. `Inflation and Localization to Improve Performance. <https://ncar.github.io/dart-tutorial/DART_LAB_Section03.pdf>`_

4. `Nonlinear and Non-Gaussian Extensions. <https://ncar.github.io/dart-tutorial/DART_LAB_Section04.pdf>`_

5. `Adaptive Inflation. <https://ncar.github.io/dart-tutorial/DART_LAB_Section05.pdf>`_


MATLAB hands-on exercises
-------------------------

In the ``guide/DART_LAB/matlab`` directory are a set of MATLAB scripts and GUI (graphical user interface) programs which are
exercises that go with the tutorial. Each is interactive with settings that can be changed and rerun to explore various
options. A valid `MATLAB <http://www.mathworks.com/products/matlab/>`__ license is needed to run these scripts.

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

To run these, cd into the ``guide/DART_LAB/matlab`` directory, start MATLAB, and type the names at the prompt.
