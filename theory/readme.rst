.. _DART-tutorial:

DART Tutorial
=============

The slides for the Data Assimilation Research Testbed tutorial:  
`DART_LAB part 6: Using the Real DART System <https://ncar.github.io/dart-tutorial/DART_LAB_Section06.pdf>`_,
guide you through running ensemble data-assimilation experiments using the real DART software environment.

Building on assimilation concepts covered in the :ref:`DART_LAB` part 1-5 you will:

* Install and configure DART, gaining an understanding of its directory structure and build system.

* Run the Lorenz-96 model to reproduce the DART_LAB assimilation experiments within the full DART framework.

* Explore configurable parameters such as ensemble size, localization, inflation method, and 
  observation networks — and learn how they influence assimilation performance. 

* Generate synthetic observations via OSSEs (Observing System Simulation Experiments), 
  assimilate them, and diagnose results using MATLAB tools (time-series error/spread plots, 
  rank histograms, bias/variance diagnostics, etc.).

* Examine observation-space diagnostics, including rank histograms, RMSE/bias evolution, 
  and outlier-observation rejection (via a configurable threshold).

* Extend your familiarity to other low-order models (such as the Lorenz 63 system) and prepare 
  for full-scale geoscience applications using DART’s model interfaces.

The diagnostics in the DART tutorial use MATLAB®. To learn how to configure your
environment to use MATLAB and the DART diagnostics, see the documentation for
:doc:`/guide/matlab-observation-space`.


.. raw:: html

   <iframe src="https://ncar.github.io/dart-tutorial/DART_LAB_Section06.pdf" width="100%" height="500px">
   </iframe>