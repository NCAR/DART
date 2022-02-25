9-variable
==========

Overview
--------

The 9-variable model is described in Lorenz (1980). [1]_ Lorenz developed this
primitive-equation model using shallow-water equations as a starting point and
manipulating the divergence equations so that the model exhibits
quasi-geostrophic behavior and transient gravity waves that dissipate with
time. Gent and McWilliams (1982) [2]_ explore the behavior of this model
extensively. For an introduction to shallow-water equations, we recommend
consulting the relevant section of a meteorology textbook such as section 4.5
of Holton and Hakim (2013). [3]_

The model's three *X* variables are at 0, 1/9, and 2/9, three *Y* variables are
at 3/9, 4/9 and 5/9, and three *Z* variables are at 6/9, 7/9, and 8/9 on a
cyclic [0, 1] domain.

In the 9-variable model, DART advances the model, gets the model state and
metadata describing this state. The model can be configured by altering the
``&model_nml`` `namelist`_ in the ``input.nml`` file. The details of the
``&model_nml`` namelist are always model-specific (there are no generic
namelist values). The model time step defaults to 1 hour (3600 seconds) but is
settable by altering the namelist.

The 9-variable model has a ``work/workshop_setup.csh`` script that compiles 
and runs an example. This example is referenced in Sections 7 and 10 of the
:doc:`DART_tutorial <../../../theory/readme>`
and is intended to provide insight into model/assimilation behavior.
The example **may or may not** result in good (*or even decent!*) results!

Namelist
--------

The ``&model_nml`` namelist is read from the ``input.nml`` file. Namelists
start with an ampersand ``&`` and terminate with a slash ``/``. Character
strings that contain a ``/`` must be enclosed in quotes to prevent them from
prematurely terminating the namelist.

.. code-block:: fortran

  &model_nml
     g                 = 8.0,
     deltat            = 0.0833333333333333,
     time_step_days    = 0,
     time_step_seconds = 3600
  /

Description of each namelist entry
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

+-------------------+----------+-------------------------------------+
| Item              | Type     | Description                         |
+===================+==========+=====================================+
| g                 | real(r8) | Model parameter, see comp_dt in     |
|                   |          | code for equations.                 |
+-------------------+----------+-------------------------------------+
| delta_t           | real(r8) | Non-dimensional timestep. This is   |
|                   |          | mapped to the dimensional timestep  |
|                   |          | specified by time_step_days and     |
|                   |          | time_step_seconds.                  |
+-------------------+----------+-------------------------------------+
| time_step_days    | real(r8) | Number of days for dimensional      |
|                   |          | timestep, mapped to delta_t.        |
+-------------------+----------+-------------------------------------+
| time_step_seconds | real(r8) | Number of seconds for dimensional   |
|                   |          | timestep, mapped to delta_t.        |
+-------------------+----------+-------------------------------------+

References
----------

.. [1] Lorenz, Edward N., 1980: Attractor Sets and Quasi-Geostrophic
   Equilibrium. *Journal of the Atmospheric Sciences*, **37**, 1685-1699.
   `doi:10.1175/1520-0469(1980)037\<1685:ASAQGE\>2.0.CO;2
   <https://doi.org/10.1175/1520-0469(1980)037\<1685:ASAQGE\>2.0.CO;2>`__
.. [2] Gent, Peter R., and James C. McWilliams, 1982: Intermediate Model
   Solutions to the Lorenz Equations: Strange Attractors and Other Phenomena.
   *Journal of the Atmospheric Sciences*, **39**, 3-13.
   `doi:10.1175/1520-0469(1982)039\<0003:IMSTTL\>2.0.CO;2
   <https://doi.org/10.1175/1520-0469(1982)039\<0003:IMSTTL\>2.0.CO;2>`__
.. [3] Holton, James R., and Gregory J. Hakim, 2013: *An Introduction to
   Dynamic Meteorology -- Fifth Edition.* Academic Press, 532 pp.
