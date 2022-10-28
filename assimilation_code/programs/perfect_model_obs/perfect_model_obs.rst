program ``perfect_model_obs``
=============================

Overview
--------

Main program for creating synthetic observation sequences given a model for use in filter assimilations. Reads in an
observation sequence file which has only observation definitions and generates synthetic observation values for an
output observation sequence file. The execution of perfect_model_obs is controlled by the input observation sequence
file and the model time-stepping capabilities in a manner analogous to that used by the filter program.

Namelist
--------

This namelist is read from the file ``input.nml``. Namelists start with an ampersand '&' and terminate with a slash '/'.
Character strings that contain a '/' must be enclosed in quotes to prevent them from prematurely terminating the
namelist.

::

   &perfect_model_obs_nml
      single_file_in             = .false.,
      read_input_state_from_file = .false.,
      input_state_files          = "",
      init_time_days             = 0,
      init_time_seconds          = 0,

      single_file_out            = .false.,
      output_state_files         = "",
      write_output_state_to_file = .false.,
      output_interval            = 1,

      async                      = 0,
      adv_ens_command            = "./advance_model.csh",
      tasks_per_model_advance    = 1,

      obs_seq_in_file_name       = "obs_seq.in",
      obs_seq_out_file_name      = "obs_seq.out",
      first_obs_days             = -1,
      first_obs_seconds          = -1,
      last_obs_days              = -1,
      last_obs_seconds           = -1,
      obs_window_days            = -1,
      obs_window_seconds         = -1,

      trace_execution            = .false.,
      output_timestamps          = .false.,
      print_every_nth_obs        = 0,
      output_forward_op_errors   = .false.,
      silence                    = .false.,
   /

| 

.. container::

   +---------------------------------------+---------------------------------------+---------------------------------------+
   | Item                                  | Type                                  | Description                           |
   +=======================================+=======================================+=======================================+
   | read_input_state_from_file            | logical                               | If false, model_mod must provide the  |
   |                                       |                                       | input state.                          |
   +---------------------------------------+---------------------------------------+---------------------------------------+
   | single_file_in                        | logical                               | Get all states from a single file.    |
   +---------------------------------------+---------------------------------------+---------------------------------------+
   | input_state_files                     | character(len=256)                    | A list of files, one per domain. Each |
   |                                       | dimension(MAX_NUM_DOMS)               | file must be a text file containing   |
   |                                       |                                       | the name of the NetCDF file to open.  |
   +---------------------------------------+---------------------------------------+---------------------------------------+
   | write_output_state_to_file            | logical                               | If false, state is not written out.   |
   +---------------------------------------+---------------------------------------+---------------------------------------+
   | single_file_out                       | logical                               | Write all states to a single file.    |
   +---------------------------------------+---------------------------------------+---------------------------------------+
   | output_state_files                    | character(len=256)                    | A list of files, one per domain. Each |
   |                                       | dimension(MAX_NUM_DOMS)               | file must be a text file containing   |
   |                                       |                                       | the names of the NetCDF file to open. |
   +---------------------------------------+---------------------------------------+---------------------------------------+
   | init_time_days                        | integer                               | If negative, don't use. If            |
   |                                       |                                       | non-negative, override the initial    |
   |                                       |                                       | data time read from restart file.     |
   +---------------------------------------+---------------------------------------+---------------------------------------+
   | init_time_seconds                     | integer                               | If negative don't use. If             |
   |                                       |                                       | non-negative, override the initial    |
   |                                       |                                       | data time read from restart file.     |
   +---------------------------------------+---------------------------------------+---------------------------------------+
   | output_interval                       | integer                               | Output state and observation          |
   |                                       |                                       | diagnostics every nth assimilation    |
   |                                       |                                       | time, n is output_interval.           |
   +---------------------------------------+---------------------------------------+---------------------------------------+
   | async                                 | integer                               | Controls method for advancing model:  |
   |                                       |                                       |                                       |
   |                                       |                                       | -  0 = subroutine call                |
   |                                       |                                       | -  2 = shell command, single task     |
   |                                       |                                       |    model                              |
   |                                       |                                       | -  4 = shell command, parallel model  |
   +---------------------------------------+---------------------------------------+---------------------------------------+
   | adv_ens_command                       | character(len=129)                    | Command sent to shell if async == 2   |
   |                                       |                                       | or 4.                                 |
   +---------------------------------------+---------------------------------------+---------------------------------------+
   | tasks_per_model_advance               | integer                               | Number of tasks to use while          |
   |                                       |                                       | advancing the model.                  |
   +---------------------------------------+---------------------------------------+---------------------------------------+
   | obs_seq_in_file_name                  | character(len=256)                    | File name from which to read an       |
   |                                       |                                       | observation sequence.                 |
   +---------------------------------------+---------------------------------------+---------------------------------------+
   | obs_seq_out_file_name                 | character(len=256)                    | File name to which to write output    |
   |                                       |                                       | observation sequence.                 |
   +---------------------------------------+---------------------------------------+---------------------------------------+
   | first_obs_days                        | integer                               | If negative, don't use. If            |
   |                                       |                                       | non-negative, ignore any observations |
   |                                       |                                       | before this time.                     |
   +---------------------------------------+---------------------------------------+---------------------------------------+
   | first_obs_seconds                     | integer                               | If negative, don't use. If            |
   |                                       |                                       | non-negative, ignore any observations |
   |                                       |                                       | before this time.                     |
   +---------------------------------------+---------------------------------------+---------------------------------------+
   | last_obs_days                         | integer                               | If negative, don't use. If            |
   |                                       |                                       | non-negative, ignore any observations |
   |                                       |                                       | after this time.                      |
   +---------------------------------------+---------------------------------------+---------------------------------------+
   | last_obs_seconds                      | integer                               | If negative, don't use. If            |
   |                                       |                                       | non-negative, ignore any observations |
   |                                       |                                       | after this time.                      |
   +---------------------------------------+---------------------------------------+---------------------------------------+
   | obs_window_days                       | integer                               | If negative, don't use. If            |
   |                                       |                                       | non-negative, reserved for future     |
   |                                       |                                       | use.                                  |
   +---------------------------------------+---------------------------------------+---------------------------------------+
   | obs_window_seconds                    | integer                               | If negative, don't use. If            |
   |                                       |                                       | non-negative, reserved for future     |
   |                                       |                                       | use.                                  |
   +---------------------------------------+---------------------------------------+---------------------------------------+
   | trace_execution                       | logical                               | True means output very detailed       |
   |                                       |                                       | messages about what routines are      |
   |                                       |                                       | being called in the main loop. Useful |
   |                                       |                                       | if a job hangs or otherwise doesn't   |
   |                                       |                                       | execute as expected.                  |
   +---------------------------------------+---------------------------------------+---------------------------------------+
   | output_timestamps                     | logical                               | True means output timestamps before   |
   |                                       |                                       | and after the model advance and the   |
   |                                       |                                       | forward observation computation       |
   |                                       |                                       | phases.                               |
   +---------------------------------------+---------------------------------------+---------------------------------------+
   | print_every_nth_obs                   | integer                               | If negative, don't use. If            |
   |                                       |                                       | non-negative, print a message noting  |
   |                                       |                                       | the processing of every Nth           |
   |                                       |                                       | observation.                          |
   +---------------------------------------+---------------------------------------+---------------------------------------+
   | output_forward_op_errors              | logical                               | True means output errors from forward |
   |                                       |                                       | observation operators. This is the    |
   |                                       |                                       | 'istatus' error return code from the  |
   |                                       |                                       | model interpolate routine. An ascii   |
   |                                       |                                       | text file 'forward_op_errors' will be |
   |                                       |                                       | created in the current directory.     |
   |                                       |                                       | Each line will contain an observation |
   |                                       |                                       | key number, and the istatus return    |
   |                                       |                                       | code.                                 |
   +---------------------------------------+---------------------------------------+---------------------------------------+
   | silence                               | logical                               | True means output almost no runtime   |
   |                                       |                                       | messages. Not recommended for general |
   |                                       |                                       | use, but can speed test programs if   |
   |                                       |                                       | the execution time becomes dominated  |
   |                                       |                                       | by the volume of output.              |
   +---------------------------------------+---------------------------------------+---------------------------------------+

| 

Modules used
------------

::

   types_mod
   utilities_mod
   time_manager_mod
   obs_sequence_mod
   obs_def_mod
   obs_model_mod
   assim_model_mod
   mpi_utilities_mod
   random_seq_mod
   ensemble_manager_mod

Files
-----

-  observation sequence input file; name comes from obs_seq_in_file_name
-  observation sequence output file; name comes from obs_seq_out_file_name
-  input state vector file; name comes from restart_in_file_name
-  output state vector file; name comes from restart_out_file_name
-  perfect_model_mod.nml in input.nml

References
----------

-  none
