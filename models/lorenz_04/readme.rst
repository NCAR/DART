Lorenz 05
=========

Naming History
--------------

In earlier versions of DART, this collection of models was referred to as
Lorenz 04. Edward Lorenz provided James A. Hansen these model formulations
before they had been published, since both Lorenz and Hansen were faculty
members at MIT at the time. Hansen developed the DART model interface and
incorporated it into the DART codebase in 2004. Thus, within DART, it was named
Lorenz 04.

The collection of models was published a year later in Lorenz (2005), [1]_
thus, within the wider community, the models are typically referred to as
Lorenz 05. To reflect this fact, the collection of models was renamed within
DART from Lorenz 04 to Lorenz 05 during the Manhattan release.

Overview
--------

Lorenz (2005) provides a fascinating account of the difficulties involved in
designing simple models that exhibit chaotic behavior and realistically
simulate aspects of atmospheric flow. It presents three models of increasing
complexity:

- Model I is a single-scale model, similar to Lorenz (1996), [2]_ intended to
  represent the atmosphere at a specific height and latitude.
- Model II is also a single-scale model, similar to Model I, but with spatial
  continuity in the waves.
- Model III is a two-scale model. It is fundamentally different from the Lorenz
  96 two-scale model because of the spatial continuity and the fact that both
  scales are projected onto a single variable of integration. The scale
  separation is achieved by a spatial filter and is therefore not perfect (i.e.
  there is leakage).

Model II and Model III are implemented in this DART model interface, and the
user is free to choose Model II or III by editing the `namelist`_. For users
interested in Model I, please use Lorenz 96. The slow scale in Model III is
Model II, and thus Model II is a deficient form of Model III.

The Lorenz 05 model has a ``work/workshop_setup.csh`` script that compiles and 
runs an example.  This example may be used anywhere in the
:doc:`DART tutorial <../../theory/readme>` to explore 
multiscale dynamics
and to provide insight into model/assimilation behavior.
The example **may or may not** result in good (*or even decent!*) results!

Model Formulation
~~~~~~~~~~~~~~~~~

For Lorenz 05, DART to advances the model, gets the model state and metadata
describing this state, finds state variables that are close to a given
location, and does spatial interpolation for model state variables.

Namelist
--------

The ``&model_nml`` namelist is read from the ``input.nml`` file. Namelists
start with an ampersand ``&`` and terminate with a slash ``/``. Character
strings that contain a ``/`` must be enclosed in quotes to prevent them from
prematurely terminating the namelist.

.. code-block:: fortran

  &model_nml
     model_size        = 960,
     forcing           = 15.00,
     delta_t           = 0.001,
     space_time_scale  = 10.00,
     coupling          = 3.00,
     K                 = 32,
     smooth_steps      = 12,
     time_step_days    = 0,
     time_step_seconds = 3600,
     model_number      = 3
  /

Description of each namelist entry
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

+-------------------+----------+-------------------------------------+
| Contents          | Type     | Description                         |
+===================+==========+=====================================+
| model_size        | integer  | Number of variables in model        |
+-------------------+----------+-------------------------------------+
| forcing           | real(r8) | Forcing, F, for model               |
+-------------------+----------+-------------------------------------+
| delta_t           | real(r8) | Non-dimensional timestep            |
+-------------------+----------+-------------------------------------+
| space_time_scale  | real(r8) | Determines temporal and spatial     |
|                   |          | relationship between fast and slow  |
|                   |          | variables (model III)               |
+-------------------+----------+-------------------------------------+
| coupling          | real(r8) | Linear coupling between fast and    |
|                   |          | slow variables (model III)          |
+-------------------+----------+-------------------------------------+
| K                 | integer  | Determines the wavenumber of the    |
|                   |          | slow variables (K=1, smooth_steps=0 |
|                   |          | reduces model II to Lorenz 96)      |
+-------------------+----------+-------------------------------------+
| smooth_steps      | integer  | Determines filter length to         |
|                   |          | separate fast and slow scales       |
+-------------------+----------+-------------------------------------+
| time_step_days    | integer  | Arbitrary real time step days       |
+-------------------+----------+-------------------------------------+
| time_step_seconds | integer  | Arbitrary real time step seconds    |
|                   |          | (could choose this for proper       |
|                   |          | scaling)                            |
+-------------------+----------+-------------------------------------+
| model_number      | integer  | 2 = single-scale, 3 = 2-scale.      |
|                   |          | (This follows the notation in the   |
|                   |          | paper.)                             |
+-------------------+----------+-------------------------------------+

References
----------

.. [1] Lorenz, Edward N., 2005: Designing Chaotic Models. *Journal of the Atmospheric Sciences*, **62**, 1574-1587.
.. [2] Lorenz, Edward N., 1996: Predictability: A Problem Partly Solved. *Seminar on Predictability.* **1**, ECMWF, Reading, Berkshire, UK, 1-18.
