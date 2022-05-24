Lorenz 96 2-scale
=================

Overview
--------

The Lorenz 96 2-scale model was first described by Edward Lorenz during a
seminar at the European Centre for Medium-Range Weather Forecasts in the Autumn
of 1995, the proceedings of which were published as Lorenz (1996) [1]_ the
following year, hence the model is commonly referred to as Lorenz 96.

The model state varies on two separate time scales, one for the X dimension and
another in the Y dimension. It is constructed by coupling together two
implementations of the Lorenz 96 single-scale model. The constant *F* term in
Lorenz 96 single-scale model is replaced by a term that couples the two scales
together.

Lorenz 96 2-scale is a widely studied model because the differing timescales
can be viewed as an analog of processes that occur on different time and
spatial scales in the atmosphere such as large-scale flow and localized
convection. The `references`_ contain some of the earlier studies including
Palmer (2001), [2]_ Smith (2001), [3]_ Orrell (2002), [4]_ Orrel (2003), [5]_
Vannitsem and Toth (2002), [6]_ Roulston and Smith (2003), [7]_ and Wilks
(2005). [8]_

The Lorenz 96 2-scale model has a ``work/workshop_setup.csh`` script that 
compiles and runs an example. This example may be explored in the
:doc:`DART tutorial <../../theory/readme>`
and is intended to provide insight into model/assimilation behavior.
The example **may or may not** result in good (*or even decent!*) results!

Development History
~~~~~~~~~~~~~~~~~~~

This DART model interface was developed by Josh Hacker as an adaptation of 
the Lorenz 96 implementation. The 2-scale model is the second model
described in Lorenz (1996).

Quick Start
-----------

To run Lorenz 96 2-scale with its default settings:

1. Ensure you have the correct settings in mkmf.template in
   ``<DARTROOT>/build_templates/mkmf.template``
2. Build the DART executables using the ``quickbuild.sh`` script in the
   ``./work`` directory.
3. Once the executables have been built, the two Perl scripts provided in the
   ``./shell_scripts`` directory, ``spinup_model.pl`` and ``run_expt.pl``, can
   be used to spin up the model and run an experiment.

Namelist
--------

The model also implements the variant of Smith (2001), which can be invoked by
setting ``local_y = .true.`` in the ``&model_nml`` namelist in the
``input.nml`` file.

The ``&model_nml`` namelist is read from the ``input.nml`` file. Namelists
start with an ampersand ``&`` and terminate with a slash ``/``. Character
strings that contain a ``/`` must be enclosed in quotes to prevent them from
prematurely terminating the namelist.

.. code-block:: fortran

  &model_nml
     model_size_x        = 36,
     y_per_x             = 10,
     forcing             = 15.00,
     delta_t             = 0.005,
     coupling_b          = 10.0,
     coupling_c          = 10.0,
     coupling_h          = 1.0,
     local_y             = .false.,
     time_step_days      = 0,
     time_step_seconds   = 3600
     template_file       = 'filter_input.nc'
  /


Description of each namelist entry
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

+-------------------+--------------------+-------------------------------------+
| Item              | Type               | Description                         |
+===================+====================+=====================================+
| model_size_x      | integer            | Number of variables in x-dimension. |
+-------------------+--------------------+-------------------------------------+
| y_per_x           | integer            | Scaling factor for number of        |
|                   |                    | variables in y-dimension compared   |
|                   |                    | to x-dimension.                     |
+-------------------+--------------------+-------------------------------------+
| forcing           | real(r8)           | Forcing, F, for model.              |
+-------------------+--------------------+-------------------------------------+
| delta_t           | real(r8)           | Non-dimensional timestep. This is   |
|                   |                    | mapped to the dimensional timestep  |
|                   |                    | specified by time_step_days and     |
|                   |                    | time_step_seconds.                  |
+-------------------+--------------------+-------------------------------------+
| coupling_b        | real(r8)           |                                     |
+-------------------+--------------------+-------------------------------------+
| coupling_c        | real(r8)           |                                     |
+-------------------+--------------------+-------------------------------------+
| coupling_h        | real(r8)           |                                     |
+-------------------+--------------------+-------------------------------------+
| local_y           | boolean            |                                     |
+-------------------+--------------------+-------------------------------------+
| time_step_days    | integer            | Number of days for dimensional      |
|                   |                    | timestep, mapped to delta_t.        |
+-------------------+--------------------+-------------------------------------+
| time_step_seconds | integer            | Number of seconds for dimensional   |
|                   |                    | timestep, mapped to delta_t.        |
+-------------------+--------------------+-------------------------------------+
| template_file     | character(len=256) | this in script                      |
+-------------------+--------------------+-------------------------------------+

References
~~~~~~~~~~

.. [1] Lorenz, Edward N., 1996: Predictability: A Problem Partly Solved. *Seminar on Predictability*. **1**, ECMWF, Reading, Berkshire, UK, 1-18.

.. [2] Palmer, Timothy N., 2001: A nonlinear dynamical perspective on model error: A proposal for non‐local stochastic‐dynamic parametrization in weather and climate prediction models. *Quarterly Journal of the Royal Meteorological Society*, **127**, 279–304. https://doi.org/10.1002/qj.49712757202

.. [3] Smith, Leonard A., 2001: Disentangling uncertainty and error: On the predictability of nonlinear systems. *Nonlinear dynamics and statistics,* Alistair I. Mees, Editor, Birkhauser, Boston, USA, 31–64.

.. [4] Orrell, David, 2002: Role of the metric in forecast error growth: How chaotic is the weather? *Tellus*, **54A**, 350–362.

.. [5] Orrell, David, 2003: Model error and predictability over different timescales in the Lorenz '96 Systems. *Journal of the Atmospheric Sciences*, **60**, 2219–2228.

.. [6] Vannitsem, Stéphane and Zoltan Toth, 2002: Short-term dynamics of model errors. *Journal of the Atmospheric Sciences*, **59**, 2594–2604.

.. [7] Roulston, Mark S. and Leonard A. Smith, 2003: Combining dynamical and statistical ensembles. *Tellus*, **55A**, 16–30.

.. [8] Wilks, Daniel S., 2005: Effects of stochastic parametrizations in the Lorenz ’96 system. *Quarterly Journal of the Royal Meteorological Society*. **131**. 389-407. https://doi.org/10.1256/qj.04.03
