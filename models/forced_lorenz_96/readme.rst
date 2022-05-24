Forced Lorenz 96
================

Overview
--------

The *forced_lorenz_96* model implements the standard Lorenz (1996) [1]_
equations except that the forcing term, ``F``, is added to the state vector and
is assigned an independent value at each gridpoint. The result is a model that
is twice as big as the standard L96 model. The forcing can be allowed to vary in
time or can be held fixed so that the model looks like the standard L96 but with
a state vector that includes the constant forcing term. An option is also
included to add random noise to the forcing terms as part of the time tendency
computation which can help in assimilation performance. If the random noise
option is turned off (see namelist) the time tendency of the forcing terms is 0.

DART state vector composition:

+----------------------------------------+----------------------------------------+
| **state variables**                    |  **forcing terms**                     |
+========================================+========================================+
| *traditional Lorenz_96 state*          | *"extended" state*                     |
+----------------------------------------+----------------------------------------+
| ``indices 1 - 40``                     | ``indices 41 - 80``                    |
+----------------------------------------+----------------------------------------+

The *forced_lorenz_96* model has a ``work/workshop_setup.csh`` script that 
compiles and runs an example.  This example is referenced in Section 20 of the
`DART_tutorial <https://dart.ucar.edu/pages/Tutorial.html>`__ 
and is intended to provide insight into parameter estimation and model/assimilation 
behavior. 
Be aware that the ``input.nml`` file is modified by the ``workshop_setup.csh`` script.

Quick Start
-----------

To become familiar with the model, try this quick experiment.

#. compile everything in the ``model/forced_lorenz_96/work`` directory.

   .. code-block::
   
      cd $DARTROOT/models/forced_lorenz_96/work
      ./quickbuild.ch

#. make sure the ``input.nml`` looks like the following (there is a lot
   that has been left out for clarity, these are the settings of
   interest for this example):

   .. code-block:: fortran

      &perfect_model_obs_nml
         start_from_restart    = .true.,
         output_restart        = .true.,
         async                 = 0,
         restart_in_file_name  = "perfect_ics",
         obs_seq_in_file_name  = "obs_seq.in",
         obs_seq_out_file_name = "obs_seq.out",
         ...
      /

      &filter_nml
         async                    = 0,
         ens_size                 = 80,
         start_from_restart       = .true.,
         output_restart           = .true.,
         obs_sequence_in_name     = "obs_seq.out",
         obs_sequence_out_name    = "obs_seq.final",
         restart_in_file_name     = "filter_ics",
         restart_out_file_name    = "filter_restart",
         num_output_state_members = 80,
         num_output_obs_members   = 80,
         ...
      /

      &model_nml
         num_state_vars    = 40,
         forcing           = 8.00,
         delta_t           = 0.05,
         time_step_days    = 0,
         time_step_seconds = 3600,
         reset_forcing     = .false.,
         random_forcing_amplitude = 0.10
      /

#. Run ``perfect_model_obs`` to generate ``true_state.nc`` and
   ``obs_seq.out``. The default ``obs_seq.in`` will cause the model to
   advance for 1000 time steps.

   .. code-block::

      ./perfect_model_obs

#. If you have *ncview*, explore the ``true_state.nc``. Notice that the
   State Variable indices from 1-40 are the dynamical part of the model
   and 41-80 are the Forcing variables.

   .. code-block::
   
      ncview true_state.nc

#. Run ``filter`` to generate ``preassim.nc``, ``analysis.nc`` and
   ``obs_seq.final``.

   .. code-block::

      ./filter

#. Launch Matlab and run ``plot_ens_time_series``.

   .. code-block::

      >> plot_ens_time_series
      Input name of prior or posterior diagnostics file for preassim.nc:
      preassim.nc
      OPTIONAL: if you have the true state and want it superimposed,
      provide the name of the input file. If not, enter a dummy filename.
      Input name of True State file for true_state.nc:
      true_state.nc
      Using state state variable IDs 1 13 27
      If these are OK, ;
      If not, please enter array of state variable ID's
      To choose from entire state enter A 25 50 75 (between 1 and 80)
      To choose traditional model state enter S 1 23 40 (between 1 and 40)
      To choose forcing estimates enter F 2 12 22 (between 1 and 40)
      (no intervening syntax required)
      A 20 30 40 60 70 80

   Indices 20, 30, and 40 will be from the dynamical part of the
   lorenz_96 attractor, indices 60, 70, and 80 will be the corresponding
   Forcing values. Here are some images for just indices 20 and 60.
   Click on each image for a high-res version.

Repeat the experiment with *reset_forcing = .true.* when creating the
true state and *reset_forcing = .false.* when assimilating. What
happens?

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
      num_state_vars    = 40,
      forcing           = 8.00,
      delta_t           = 0.05,
      time_step_days    = 0,
      time_step_seconds = 3600,
      reset_forcing     = .false.,
      random_forcing_amplitude = 0.10  
   /

Description of each namelist entry
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

+--------------------------+----------+----------------------------+
| Item                     | Type     | Description                |
+==========================+==========+============================+
| num_state_vars           | integer  | Number of variables in     |
|                          |          | model.                     |
+--------------------------+----------+----------------------------+
| forcing                  | real(r8) | Forcing, F, for model.     |
+--------------------------+----------+----------------------------+
| delta_t                  | real(r8) | Non-dimensional timestep.  |
+--------------------------+----------+----------------------------+
| time_step_days           | real(r8) | Base model time step maps  |
|                          |          | to this much real time.    |
+--------------------------+----------+----------------------------+
| time_step_seconds        | real(r8) | Base model time step maps  |
|                          |          | to this.                   |
+--------------------------+----------+----------------------------+
| reset_forcing            | logical  | If true, all forcing       |
|                          |          | values are held fixed at   |
|                          |          | the value specified for    |
|                          |          | the forcing namelist.      |
+--------------------------+----------+----------------------------+
| random_forcing_amplitude | real(r8) | Standard deviation of the  |
|                          |          | gaussian noise with zero   |
|                          |          | mean that is added to each |
|                          |          | forcing value's time step. |
+--------------------------+----------+----------------------------+

References
----------

.. [1] Lorenz, Edward N., 1996: Predictability: A Problem Partly Solved. *Seminar on Predictability*. **1**, ECMWF, Reading, Berkshire, UK, 1-18.
