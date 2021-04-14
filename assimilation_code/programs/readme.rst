Programs
========

This list of programs is separated (imprecisely) into groups which have similar functionality.
Within each group they are sorted (approximately) by the order
in which they might be used andor by their usefulness.

">>" marks the most commonly used and up to date programs.

Setting Up Experiments
-----------------------------------
 
`>> model_mod_check/model_mod_check.f90 <model_mod_check/model_mod_check.html>`_
   "Hackable" program to test some of the more fundamental routines in any ``model_mod``, 
   especially a new one.

`>> preprocess/preprocess.f90 <preprocess/preprocess.html>`_
   Program to insert observation specific code into DART before filter or perfect_model_obs is compiled.
 
`>> fill_inflation_restart/ fill_inflation_restart.f90* <fill_inflation_restart/fill_inflation_restart.html>`_
   Create inflation restart files with constant values taken from ``filter_nml``.

`>> obs_impact_tool/obs_impact_tool.f90 <obs_impact_tool/obs_impact_tool.html>`_
   Construct a table that is read by filter at run-time to localize the
   impact of sets of observation types on sets of state vector quantities.
 
`perturb_single_instance/perturb_single_instance.f90 <perturb_single_instance/perturb_single_instance.html>`_
   Generate an ensemble of perturbed ensemble member restart files.
   (Alternatively, you might perturb the model state using ``model_nml`` variables).
 
`gen_sampling_err_table/gen_sampling_err_table.f90 <gen_sampling_err_table/gen_sampling_err_table.html>`_
   Computes a table of values needed to apply Sampling Error Correction (SEC),
   which corrects covariances based on small sample size statistics.
 
Creating Observation Sequence Files
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

>> obs_utils/create_obs_grid.f90
   Create a set of observations located on a regular grid.  Obs have no data values.

`>> create_fixed_network_seq/create_fixed_network_seq.f90 <create_fixed_network_seq/create_fixed_network_seq.html>`_ 
   Reads observation sequence file information from standard input 
   and replicates it multiple times in a second observation sequence file, at user specified dates. 
 
obs_utils/obs_timejitter.f90
   Randomly perturb the times of the observations in a (usually) set_def.out file.
   Writes the results to (usually) obs_seq.in.

`>> create_obs_sequence/create_obs_sequence.f90 <create_obs_sequence/create_obs_sequence.html>`_
   Creates a short andor synthetic observation sequence file using values read from standard input.
 
Querying Observation Sequence Files
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

obs_utils/obs_info.f90
   Summarize obs types, times, counts found in observation sequence file(s).

>> obs_utils/obs_assim_count.f90
   Prints out a quick table of obs types and counts, overall start and stop times, 
   and metadata strings and counts.  See obs_diag for more.

obs_assim_count/obs_assim_count.f90
   Prints out a quick table of obs types and counts, overall start and stop times, 
   and metadata strings and counts.  This is an older version of obs_utils/obs_assim_count.f90.

`obs_seq_coverage/obs_seq_coverage.f90 <obs_seq_coverage/obs_seq_coverage.html>`_
   Queries a set of observation sequence files to determine which observation locations report
   frequently enough to be useful for a verification study.
 
obs_total_error/obs_total_error.f90
   Prints the total error in the mean and spread from an obs_seq file 
   which has been through both perfect_model_obs and filter, so it has copies
   'truth', 'ensemble mean', and 'ensemble spread'.
   You can get more information by running the obs_diag program.

Changing Observation Sequence Files
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

`obs_keep_a_few/obs_keep_a_few.f90 <obs_keep_a_few/obs_keep_a_few.html>`_
   Creates an output observation sequence file that is shorter than the input obs_seq file.
 
`obs_selection/obs_selection.f90 <obs_selection/obs_selection.html>`_
   Extracts observations out of one or more obs_sequence files
   according to a  list of observation types, times, and locations.
   The list is usually created by :doc:`obs_seq_coverage/obs_seq_coverage`, 
   but can be an observation sequence file,
 
`>> obs_sequence_tool/obs_sequence_tool.f90 <obs_sequence_tool/obs_sequence_tool.html>`_
   Subsets, combines, or alters observations from one or more observation sequence files 
   and optionally writes them into a single output obs_seq file.

`obs_loop/obs_loop.f90 <obs_loop/obs_loop.html>`_
   A template to read in observations from one obs_seq file and write them,
   optionally modified by user supplied code, to another obs_seq file.
 
obs_utils/obs_sort.f90
   Do a complete sort of an obs_seq file by location, observation type, then variance.
   An ancestor of obs_remove_dups.

obs_utils/obs_remove_dups.f90
   Removes duplicate observations from an obs_seq file, which involves a complete sort
   by time, location, observation type, then variance.
 
`>> obs_common_subset/obs_common_subset.f90 <obs_common_subset/obs_common_subset.html>`_
   Select the subset of observations, which were successfully assimilated, 
   from two or more assimilation cases (which used the same obs_seq.out file).
 
`obs_seq_verify/obs_seq_verify.f90 <obs_seq_verify/obs_seq_verify.html>`_
   Reorders the observations from a forecast run of DART into a structure 
   that is amenable for the evaluation of the forecast.
 

obs_utils/obs_data_denial.f90
   THIS IS NOT YET DONE!
   Help implement a data-denial experiment by randomly changing the error variance
   of N of each obs type in an observation sequence file to a huge value.
 
Assimilation Programs
-----------------------------------
 
`>> perfect_model_obs/perfect_model_obs.f90 <perfect_model_obs/perfect_model_obs.html>`_
   Creates synthetic observation sequences from a hindcast model.
 
`>> filter/filter.f90 <filter/filter.html>`_
   Main Fortran program for driving ensemble filter assimilations.

filter/filter.separate_seq.f90
   Like filter.f90, but each task updates its own sequence in obs_space_diagnostics.
   Included here just for future reference.

`>> advance_time/advance_time.f90 <advance_time/advance_time.html>`_
   Provides a shell-scripting-friendly way to increment and decrement calendar dates and times.
 
`integrate_model/integrate_model.f90 <integrate_model/integrate_model.html>`_
   Generic main program which advances a single ensemble member in ``perfect_model_obs`` 
   or the serial or parallel version of the ``filter`` program.

`integrate_model/integrate_model_parallel.f90 <integrate_model/integrate_model.html>`_
   Generic main program which advances a single
   ensemble member in ``perfect_model_obs`` or the serial ``
   or parallel version of the ``filter`` program.``

`wakeup_filter/wakeup_filter.f90 <wakeup_filter/wakeup_filter.html>`_
   For use in the "async=4" case where both the main filter program and the hindcast model are MPI programs. 
   The main MPI job script runs each of the model advances for the ensemble members, 
   and then runs this program to restart the filter program.
   
Evaluating Results
-----------------------------------
 
`obs_diag/oned/obs_diag.f90 <obs_diag/oned/obs_diag.html>`_
   Reads obs_seq.final files from a model with 1D locations calculates statistics, 
   and writes them to NetCDF files for use by Matlab (or other) plotting scripts.
 
`obs_diag/threed_cartesian/obs_diag.f90 <obs_diag/threed_cartesian/obs_diag.html>`_
   Reads obs_seq.final files from a model with a 3D cartesian coordinate system,
   calculates statistics, and writes them to NetCDF files for use by Matlab 
   (or other) plotting scripts.
 
`obs_diag/threed_sphere/obs_diag.f90 <obs_diag/threed_sphere/obs_diag.html>`_
   Reads obs_seq.final files from a model with a 3D spherical coordinate system,
   calculates statistics, and writes them to NetCDF files for use by Matlab 
   (or other) plotting scripts.

obs_diag/threed_sphere/streamflow_obs_diag.f90
   Like the threed_sphere obs_diag, but for identity
   streamflow observations, (needed for wrf-Hydro support).
 
`obs_seq_to_netcdf/obs_seq_to_netcdf.f90 <obs_seq_to_netcdf/obs_seq_to_netcdf.html>`_
   Extracts the observation components from observation sequence files and writes out
   netCDF files that can be used by other applications.
   such as ``diagnostics/matlab/plot_obs_netcdf*``

obs_seq_to_netcdf/radiance_obs/obs_seq_to_netcdf.f90
   Like obs_seq_to_netcdf/obs_seq_to_netcdf, but filters out radiance metadata
   which is not needed by the scripts which use the resulting NetCDF file.

`compare_states/compare_states.f90 <compare_states/compare_states.html>`_
   Compare fields in two NetCDF files and print out the min and max values from each file and of
   the difference between the two files.

`compute_error/compute_error.f90 <compute_error/compute_error.html>`_
   Compute the time-mean ensemble error and spread in the same manner as the DART MATLAB diagnostic
   routine ``plot_total_err``; in state space from true_state.nc and preassim.nc (or analysis.nc).
 
`closest_member_tool/closest_member_tool.f90 <closest_member_tool/closest_member_tool.html>`_
   Prints out a sorted order of which ensemble members are 'closest' to the mean, 
   where 'close' is selectable by namelist option.
 
Historical and Deprecated
-------------------------
 
`system_simulation <system_simulation/system_simulation.html>`_
   A collection of standalone programs for simulating various properties of ensembles.
   Talk to Jeff Anderson about the programs in this directory.

system_simulation/system_simulation.f90 
   This program begins attempts to analyze the value of particular 
   observations. Begin by trying to determine the value of 
   observations with a given correlation to a state variable using an 
   N member ensemble to compute the correlations.

`restart_file_tool/restart_file_tool.f90 <restart_file_tool/restart_file_tool.html>`_
   Deprecated, since in Manhattan all DART initial and restart files are in NetCDF format.
 

