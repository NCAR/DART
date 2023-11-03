
.. _DART programs:

Programs included with DART
===========================

This list of programs is separated into groups which have similar functionality.
Within each group they are sorted  by the order
in which they might be used and/or by how widely they are used.

Setting Up Experiments
-----------------------------------

In many cases, you won't need to use any programs in Setting Up Experiments
except for ``preprocess``, because you're using an existing model interface 
and have the observation sequence files.
In that case, you the programs you're looking for are probably in `Assimilation Programs`_.
 
:doc:`preprocess <preprocess/preprocess>`
   Program to insert observation specific code into DART before filter or perfect_model_obs is compiled.
 
:doc:`fill_inflation_restart <fill_inflation_restart/fill_inflation_restart>`
   Create inflation restart files with constant values taken from ``fill_inflation_restart_nml``.

:doc:`obs_impact_tool <obs_impact_tool/obs_impact_tool>`
   Construct a table that is read by filter at run-time to localize the
   impact of sets of observation types on sets of state vector quantities.
 
:doc:`model_mod_check <model_mod_check/model_mod_check>` 
  Program to test some of the more fundamental routines in any ``model_mod``, especially a for a new model.

:doc:`perturb_single_instance <perturb_single_instance/perturb_single_instance>`
   Generate an ensemble of perturbed ensemble member restart files.
   (Alternatively, you might perturb the model state using ``model_nml`` variables).
 
:doc:`gen_sampling_err_table <gen_sampling_err_table/gen_sampling_err_table>`
   Computes a table of values needed to apply Sampling Error Correction (SEC),
   which corrects covariances based on small sample size statistics.
 
Creating Observation Sequence Files
-----------------------------------

:doc:`create_obs_sequence <create_obs_sequence/create_obs_sequence>`
   Creates a short andor synthetic observation sequence file using values read from standard input.
 
:doc:`create_fixed_network_seq <create_fixed_network_seq/create_fixed_network_seq>` 
   Reads observation sequence file information from standard input 
   and replicates it multiple times in a second observation sequence file, at user specified dates. 
 
obs_utils/create_obs_grid
   Create a set of observations located on a regular grid.  
   Obs have no data values, but they are time ordered.

obs_utils/obs_timejitter
   Randomly perturb the times of the observations in a (usually) set_def.out file.
   Writes the results to (usually) obs_seq.in.

Querying Observation Sequence Files
-----------------------------------

obs_utils/obs_info
   Summarize obs types, times, counts found in observation sequence file(s).

obs_utils/obs_assim_count
   Prints out a quick table of obs types and counts, overall start and stop times, 
   and metadata strings and counts.  See obs_diag for more.
   There is an older version in the obs_assim_count directory.

:doc:`obs_seq_coverage <obs_seq_coverage/obs_seq_coverage>`
   Queries a set of observation sequence files to determine which observation locations report
   frequently enough to be useful for a verification study.
 
obs_total_error
   Prints the total error in the mean and spread from an obs_seq file 
   which has been through both perfect_model_obs and filter, so it has copies
   'truth', 'ensemble mean', and 'ensemble spread'.
   You can get more information by running the obs_diag program.

Changing Observation Sequence Files
-----------------------------------

:doc:`obs_sequence_tool <obs_sequence_tool/obs_sequence_tool>`
   Subsets, combines, or alters observations from one or more observation sequence files 
   and optionally writes them into a single output obs_seq file.

:doc:`obs_loop <obs_loop/obs_loop>`
   A template to read in observations from one obs_seq file and write them,
   optionally modified by user supplied code, to another obs_seq file.
 
obs_utils/obs_sort
   Do a complete sort of an obs_seq file by location, observation type, then variance.
   An ancestor of obs_remove_dups.

obs_utils/obs_remove_dups
   Removes duplicate observations from an obs_seq file, which involves a complete sort
   by time, location, observation type, then variance.
 
:doc:`obs_selection <obs_selection/obs_selection>`
   Extracts observations out of one or more obs_sequence files
   according to a  list of observation types, times, and locations.
   The list is usually created by :doc:`obs_seq_coverage <obs_seq_coverage/obs_seq_coverage>`, 
   but can be an observation sequence file.
 
:doc:`obs_common_subset <obs_common_subset/obs_common_subset>`
   Select the subset of observations, which were successfully assimilated, 
   from two or more assimilation cases (which used the same obs_seq.out file).
 
:doc:`obs_keep_a_few <obs_keep_a_few/obs_keep_a_few>`
   Creates an output observation sequence file that is shorter than the input obs_seq file.
 
:doc:`obs_seq_verify <obs_seq_verify/obs_seq_verify>`
   Reorders the observations from a forecast run of DART into a structure 
   that is amenable for the evaluation of the forecast.
 

obs_utils/obs_data_denial
   THIS IS NOT YET DONE!
   Help implement a data-denial experiment by randomly changing the error variance
   of N of each obs type in an observation sequence file to a huge value.
 
Assimilation Programs
-----------------------------------
 
:doc:`perfect_model_obs <perfect_model_obs/perfect_model_obs>`
   Creates synthetic observation sequences from a hindcast model.
 
:doc:`filter <filter/filter>`
   Main Fortran program for driving ensemble filter assimilations.

:doc:`advance_time <advance_time/advance_time>`
   Provides a shell-scripting-friendly way to increment and decrement calendar dates and times.
 
:doc:`integrate_model <integrate_model/integrate_model>`
   Generic main program which advances a single ensemble member in ``perfect_model_obs`` 
   or the serial or parallel version of the ``filter`` program.

Evaluating Results
-----------------------------------
 
obs_diag 
   Reads obs_seq.final files, calculates statistics, and writes them to NetCDF files 
   for use by Matlab (or other) plotting scripts.
   There are separate versions for models with different coordinate systems:

   - :doc:`1D <obs_diag/oned/obs_diag>`
   - :doc:`3D Cartesian <obs_diag/threed_cartesian/obs_diag>`
   - :doc:`3D spherical <obs_diag/threed_sphere/obs_diag>`
   - 3D spherical with streamflow.
   
obs_seq_to_netcdf
   Extracts the observation components from observation sequence files and writes out
   netCDF files that can be used by other applications.
   such as ``diagnostics/matlab/plot_obs_netcdf*``
   There are two versions; the :doc:`standard version <obs_seq_to_netcdf/obs_seq_to_netcdf>`
   and one which filters out radiance metadata which is not needed by the scripts 
   which use the resulting NetCDF file.

:doc:`compare_states <compare_states/compare_states>`
   Compare fields in two NetCDF files and print out the min and max values from each file and of
   the difference between the two files.

:doc:`compute_error <compute_error/compute_error>`
   Compute the time-mean ensemble error and spread in the same manner as the DART MATLAB diagnostic
   routine ``plot_total_err``; in state space from true_state.nc and preassim.nc (or analysis.nc).
 
:doc:`closest_member_tool <closest_member_tool/closest_member_tool>`
   Prints out a sorted order of which ensemble members are 'closest' to the mean, 
   where the method for computing the 'close' metric is selectable by namelist option.
 
Historical and Deprecated
-------------------------
 
:doc:`system_simulation <system_simulation/system_simulation>`
   A collection of standalone programs for simulating various properties of ensembles.
   Talk to Jeff Anderson about the programs in this directory.

:doc:`wakeup_filter <wakeup_filter/wakeup_filter>`
   For use in the "async=4" case where both the main filter program and the hindcast model are MPI programs. 
   The main MPI job script runs each of the model advances for the ensemble members, 
   and then runs this program to restart the filter program.
   
