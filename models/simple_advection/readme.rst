Simple advection
================

Overview
--------

This simple advection model simulates a wind field using Burger's Equation with
an upstream semi-lagrangian differencing on a periodic one-dimensional domain.
This diffusive numerical scheme is stable and forcing is provided by adding in 
random gaussian noise to each wind grid variable independently at each 
timestep. The domain mean value of the wind is relaxed to a constant fixed 
value set by the namelist parameter ``mean_wind``. The random forcing magnitude is 
set by namelist parameter ``wind_random_amp`` and the damping of the mean wind is 
controlled by parameter ``wind_damping_rate``. An Eulerian option with centered in 
space differencing is also provided and can be used by setting namelist 
parameter ``lagrangian_for_wind`` to ``.false.`` The Eulerian differencing is both 
numerically unstable and subject to shock formation. However, it can sometimes 
be made stable in assimilation mode (see recent work by Majda and
collaborators).

The model state includes a single passive tracer that is advected by the wind
field using semi-lagrangian upstream differencing. The state also includes a
tracer source value at each gridpoint. At each time step, the source is added
into the concentration at each gridpoint. There is also a constant global
destruction of tracer that is controlled by the namelist parameter
destruction_rate. The appropriate percentage of tracer is destroyed at each
gridpoint at each timestep.

The model also includes an associated model for the tracer source rate. At each
gridpoint, there is a value of the time mean source rate and a value of the
phase offset for a diurnal component of the source rate. The diurnal source
rate has an amplitude that is proportional to the source rate (this proportion
is controlled by namelist parameter ``source_diurnal_rel_amp``). At each grid
point, the source is the sum of the source rate plus the appropriate diurnally
varying component. The phase_offset at the gridpoint controls the diurnal
phase. The namelist parameter ``source_phase_noise`` controls the amplitude of
random gaussian noise that is added into the source phase at each time step.
If ``source_phase_noise`` is zero then the phase offset is fixed. Finally, the time
mean source rate is constant in time in the present model version. The time
mean source rate controls the amplitude of the diurnal cycle of the tracer
source.

For the simple advection model, DART advances the model, gets the model state
and metadata describing this state, finds state variables that are close to a
given location, and does spatial interpolation for model state variables.

The simple advection model has a ``work/workshop_setup.csh`` script that compiles 
and runs an example.  This example is referenced in Section 25 of the
:doc:`DART tutorial <../../theory/readme>`.
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
     num_grid_points        = 10,
     grid_spacing_meters    = 100000.0,
     time_step_days         = 0,
     time_step_seconds      = 3600,
     mean_wind              = 20.0,
     wind_random_amp        = 0.00027778,
     wind_damping_rate      = 0.0000027878,
     lagrangian_for_wind    = .true.,
     destruction_rate       = 0.000055556,
     source_random_amp_frac = 0.00001,
     source_damping_rate    = 0.0000027878,
     source_diurnal_rel_amp = 0.05,
     source_phase_noise     = 0.0,
     output_state_vector  = .false.
  /

Description of each namelist entry
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

+------------------------+----------+----------------------------+
| Item                   | Type     | Description                |
+========================+==========+============================+
| num_grid_points        | integer  | Number of grid points in   |
|                        |          | model. State vector size   |
|                        |          | is 5 times this number.    |
+------------------------+----------+----------------------------+
| grid_spacing_meters    | integer  | Grid spacing in meters.    |
+------------------------+----------+----------------------------+
| time_step_days         | real(r8) | Number of days for         |
|                        |          | dimensional timestep,      |
|                        |          | mapped to delta_t.         |
+------------------------+----------+----------------------------+
| time_step_seconds      | real(r8) | Number of seconds for      |
|                        |          | dimensional timestep,      |
|                        |          | mapped to delta_t.         |
+------------------------+----------+----------------------------+
| mean_wind              | real(r8) | Base wind velocity         |
|                        |          | (expected value over time) |
|                        |          | in meters/second.          |
+------------------------+----------+----------------------------+
| wind_random_amp        | real(r8) | Random walk amplitude for  |
|                        |          | wind in                    |
|                        |          | meters/second\ :sup:`2`.   |
+------------------------+----------+----------------------------+
| wind_damping_rate      | real(r8) | Rate of damping towards    |
|                        |          | mean wind value in         |
|                        |          | fraction/second.           |
+------------------------+----------+----------------------------+
| lagrangian_for_wind    | logical  | Can use Lagrangian         |
|                        |          | (stable) or Eulerian       |
|                        |          | (unstable) scheme for      |
|                        |          | wind.                      |
+------------------------+----------+----------------------------+
| destruction_rate       | real(r8) | Tracer destruction rate in |
|                        |          | fraction/second.           |
+------------------------+----------+----------------------------+
| source_random_amp_frac | real(r8) | Random walk amplitude for  |
|                        |          | source as a fraction of    |
|                        |          | mean source (per           |
|                        |          | second)\ :sup:`2`.         |
+------------------------+----------+----------------------------+
| source_damping_rate    | real(r8) | Damping towards mean       |
|                        |          | source rate in             |
|                        |          | fraction/second.           |
+------------------------+----------+----------------------------+
| source_diurnal_rel_amp | real(r8) | Relative amplitude of      |
|                        |          | diurnal cycle of source    |
|                        |          | (dimensionless).           |
+------------------------+----------+----------------------------+
| source_phase_noise     | real(r8) | Amplitude of gaussian      |
|                        |          | noise to be added to       |
|                        |          | source phase offset (per   |
|                        |          | second).                   |
+------------------------+----------+----------------------------+
| output_state_vector    | logical  | Controls the output to     |
|                        |          | netCDF files. If .true.,   |
|                        |          | output the raw dart state  |
|                        |          | vector. If .false. output  |
|                        |          | the prognostic version     |
|                        |          | (gridded data) for easier  |
|                        |          | plotting (recommended).    |
+------------------------+----------+----------------------------+
