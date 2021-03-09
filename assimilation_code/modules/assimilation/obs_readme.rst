MODULE obs_model_mod
====================

Overview
--------

The code in this module computes the assimilation windows, and decides if the model needs to run in order for the data
to be at the appropriate time to assimilate the next available observations. It also has the code to write out the
current states, advance the model (in a variety of ways) and then read back in the updated states.

Other modules used
------------------

::

   types_mod
   utilities_mod
   assim_model_mod
   obs_sequence_mod
   obs_def_mod
   time_manager_mod
   ensemble_manager_mod
   mpi_utilities_mod

Public interfaces
-----------------

=========================== =============
*use obs_model_mod, only :* advance_state
\                           move_ahead
=========================== =============

| 

.. container:: routine

   *call move_ahead(ens_handle, ens_size, seq, last_key_used, window_time, key_bounds, num_obs_in_set, curr_ens_time,
   next_ens_time, trace_messages)*
   ::

      type(ensemble_type),     intent(in)  :: ens_handle
      integer,                 intent(in)  :: ens_size
      type(obs_sequence_type), intent(in)  :: seq
      integer,                 intent(in)  :: last_key_used
      type(time_type),         intent(in)  :: window_time
      integer, dimension(2),   intent(out) :: key_bounds
      integer,                 intent(out) :: num_obs_in_set
      type(time_type),         intent(out) :: curr_ens_time
      type(time_type),         intent(out) :: next_ens_time
      logical, optional,       intent(in)  :: trace_messages

.. container:: indent1

   | Given an observation sequence and an ensemble, determines how to advance the model so that the next set of
     observations can be assimilated. Also returns the first and last keys and the number of observations to be
     assimilated at this time. The algorithm implemented here (one might want to have other variants) first finds the
     time of the next observation that has not been assimilated at a previous time. It also determines the time of the
     ensemble state vectors. It then uses information about the model's time stepping capabilities to determine the time
     to which the model can be advanced that is CLOSEST to the time of the next observation. For now, this algorithm
     assumes that the model's timestep is a constant. A window of width equal to the model timestep is centered around
     the closest model time to the next observation and all observations in this window are added to the set to be
     assimilated.
   | Previous versions of this routine also made the call which actually advanced the model before returning. This is no
     longer true. The routine only determines the time stepping and number of observations. The calling code must then
     call advance_state() if indeed the next observation to be assimilated is not within the current window. This is
     determined by comparing the current ensemble time with the next ensemble time. If equal no advance is needed.
     Otherwise, next ensemble time is the target time for advance_state().

   +--------------------+------------------------------------------------------------------------------------------------+
   | ``ens_handle``     | Identifies the model state ensemble                                                            |
   +--------------------+------------------------------------------------------------------------------------------------+
   | ``ens_size``       | Number of ensemble members                                                                     |
   +--------------------+------------------------------------------------------------------------------------------------+
   | ``seq``            | An observation sequence                                                                        |
   +--------------------+------------------------------------------------------------------------------------------------+
   | ``last_key_used``  | Identifies the last observation from the sequence that has been used                           |
   +--------------------+------------------------------------------------------------------------------------------------+
   | ``window_time``    | Reserved for future use.                                                                       |
   +--------------------+------------------------------------------------------------------------------------------------+
   | ``key_bounds``     | Returned lower and upper bound on observations to be used at this time                         |
   +--------------------+------------------------------------------------------------------------------------------------+
   | ``num_obs_in_set`` | Number of observations to be used at this time                                                 |
   +--------------------+------------------------------------------------------------------------------------------------+
   | ``curr_ens_time``  | The time of the ensemble data passed into this routine.                                        |
   +--------------------+------------------------------------------------------------------------------------------------+
   | ``next_ens_time``  | The time the ensemble data should be advanced to. If equal to curr_ens_time, the model does    |
   |                    | not need to advance to assimilate the next observation.                                        |
   +--------------------+------------------------------------------------------------------------------------------------+
   | ``trace_messages`` | Optional argument. By default, detailed time trace messages are disabled but can be turned on  |
   |                    | by passing this in as .True. . The messages will print the current window times, data time,    |
   |                    | next observation time, next window time, next data time, etc.                                  |
   +--------------------+------------------------------------------------------------------------------------------------+

| 

.. container:: routine

   *call advance_state(ens_handle, ens_size, target_time, async, adv_ens_command, tasks_per_model_advance)*
   ::

      type(ensemble_type), intent(inout) :: ens_handle
      integer, intent(in)                :: ens_size
      type(time_type), intent(in)        :: target_time
      integer, intent(in)                :: async
      character(len=*), intent(in)       :: adv_ens_command
      integer, intent(in)                :: tasks_per_model_advance

.. container:: indent1

   Advances all ensemble size copies of an ensemble stored in ens_handle to the target_time. If async=0 this is done by
   repeated calls to the ``adv_1step()`` subroutine. If async=2, a call to the shell with the command
   ``adv_ens_command`` is used. If async=4, the filter program synchronizes with the MPI job shell script using the
   ``block_task()`` and ``restart_task()`` routines to suspend execution until all model advances have completed. The
   script can start the model advances using MPI and have it execute in parallel in this mode.

   +-----------------------------------------------------------+-----------------------------------------------------------+
   | ``ens_handle``                                            | Structure for holding ensemble information and data       |
   +-----------------------------------------------------------+-----------------------------------------------------------+
   | ``ens_size``                                              | Ensemble size.                                            |
   +-----------------------------------------------------------+-----------------------------------------------------------+
   | ``target_time``                                           | Time to which model is to be advanced.                    |
   +-----------------------------------------------------------+-----------------------------------------------------------+
   | ``async``                                                 | How to advance model:                                     |
   |                                                           |                                                           |
   |                                                           | +-----------------------------------------------------+   |
   |                                                           | | 0 = subroutine adv_1step                            |   |
   |                                                           | +-----------------------------------------------------+   |
   |                                                           | | 2 = shell executes adv_ens_command                  |   |
   |                                                           | +-----------------------------------------------------+   |
   |                                                           | | 4 = MPI job script advances models and syncs with   |   |
   |                                                           | | filter task                                         |   |
   |                                                           | +-----------------------------------------------------+   |
   +-----------------------------------------------------------+-----------------------------------------------------------+
   | ``adv_ens_command``                                       | Command to be issued to shell to advance model if         |
   |                                                           | async=2.                                                  |
   +-----------------------------------------------------------+-----------------------------------------------------------+
   | ``tasks_per_model_advance``                               | Reserved for future use.                                  |
   +-----------------------------------------------------------+-----------------------------------------------------------+

| 

Namelist
--------

This module does not have a namelist.

Files
-----

+------------------------------+--------------------------------------------------------------------------------------+
| filename                     | purpose                                                                              |
+==============================+======================================================================================+
| assim_model_state_ic\ *####* | a binary representation of the state vector prepended by a small header consisting   |
|                              | of the 'advance-to' time and the 'valid-time' of the state vector. The *####*        |
|                              | represents the ensemble member number if ``&ensemble_manager_nml``:                  |
|                              | ``single_restart_file_out = .true.``.                                                |
+------------------------------+--------------------------------------------------------------------------------------+
| assim_model_state_ud\ *####* | a binary representation of the state vector prepended by a small header consisting   |
|                              | of the 'valid-time' of the state vector. This is the 'updated' model state (after    |
|                              | the model has advanced the state to the desired 'advance-to' time).                  |
+------------------------------+--------------------------------------------------------------------------------------+
| filter_control\ *####*       | a text file containing information needed to advance the ensemble members; i.e., the |
|                              | ensemble member number, the input state vector file, the output state vector file -  |
|                              | that sort of thing.                                                                  |
+------------------------------+--------------------------------------------------------------------------------------+

References
----------

-  none

Private components
------------------

N/A
