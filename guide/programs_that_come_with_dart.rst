
Programs included with DART
===========================


DART comes with a collection of programs that can be used in
various ways when doing data assimilation tasks.  


Assimilation
------------

*filter* 
   `main program <../assimilation_code/programs/filter/filter.rst>`__  which performs data assimilation



OSE or OSSE - simulation experiments
------------------------------------

*perfect_model_obs* 
   generates `synthetic observations from a single model state <../assimilation_code/programs/perfect_model_obs/perfect_model_obs.rst>`__

*create_obs_seq* 
   generates `synthetic observations of any kind <../assimilation_code/programs/create_obs_seq/create_obs_seq.rst>`__

*create_fixed_network_seq* 
   generates `a time sequence of observations <../assimilation_code/programs/create_fixed_network_seq/create_fixed_network_seq.rst>`__



Diagnostics
-----------

*obs_diag* 
   `main diagnostic <../assimilation_code/programs/obs_diag/obs_diag.rst>`__ program which generates summary statistics for bins of obs

*obs_seq_to_netcdf* 
   converts DART's `obs_seq format to netcdf <../assimilation_code/programs/obs_seq_to_netcdf/obs_seq_to_netcdf.rst>`__

*matlab routines* 
   `generate plots <../diagnostics/matlab>`__  based on either the output of obs_diag or netcdf files



Observation manipulation and analysis
-------------------------------------

*obs_sequence_tool* 
   `manipulates obs_seq <../assimilation_code/programs/obs_sequence_tool/obs_sequence_tool.rst>`__ files in many ways

*obs_common_subset* 
   `selects obs <../assimilation_code/programs/obs_common_subset/obs_common_subset.rst>`__ which are common to multiple files

*obs_selection, obs_seq_coverage, obs_seq_verify*
   `analyzes time sequences <../assimilation_code/programs/obs_selection/obs_selection.rst>`__ of repeated assimilations of the same network of obs

*obs_loop* 
   `generic obs_seq program source <../assimilation_code/programs/obs_loop/obs_loop.rst>`__ that can be adapted to do specialized transformations

*observation converters*  
   see the programs in the `observations/obs_converters <../observations/obs_converters/README.rst>`__ directory

*observation utilities*  
   see the programs in the observations/utilities/threed_sphere or ../oned



Model specific programs
-----------------------

Additional programs may be available for specific models.  See the ``mkmf_xxx`` files 
in the ``models/xxx/work`` directories.



Utility programs
----------------

*preprocess* 
   generates `customized Fortran code <../assimilation_code/programs/preprocess/preprocess.rst>`__ based on what Quantities and Obs Types are used

*perturb_single_instance* 
   generates `an ensemble of states <../assimilation_code/programs/perturb_single_instance/perturb_single_instance.rst>`__ based on a single state

*closest_member_tool* 
   computes `which ensemble member is closest <../assimilation_code/programs/closest_member_tool/closest_member_tool.rst>`__ to the ensemble mean

*fill_inflation_restart* 
   generates `initial inflation restart files <../assimilation_code/programs/fill_inflation_restart/fill_inflation_restart.rst>`__ 

*model_mod_check* 
   `utility program for debugging <../assimilation_code/programs/model_mod_check/model_mod_check.rst>`__ the model_mod.f90 interface code

*obs_impact_tool* 
   generates a `custom table of localization values <../assimilation_code/programs/obs_impact_tool/obs_impact_tool.rst>`__ for obs and state quantities

*advance_time* 
   general `time utility <../assimilation_code/programs/advance_time/advance_time.rst>`__



Less commonly used Utilities
----------------------------

*compare_states* 
   `print differences <../assimilation_code/programs/compare_states/compare_states.rst>`__ between two netcdf state files

*compute_error* 
   `print the same error <../assimilation_code/programs/compute_error/compute_error.rst>`__ as the matlab routines in a standalone program

*obs_total_error* 
   `print the same total error <../assimilation_code/programs/obs_total_error/obs_total_error.f90>`__ as the matlab routines in a standalone program

*gen_sampling_err_table* 
   generate `sampling correction error tables <../assimilation_code/programs/gen_sampling_err_table/gen_sampling_err_table.rst>`__ for ensemble sizes not already present

*obs_assim_count* 
   `print number of assimilated obs <../assimilation_code/programs/obs_assim_count/obs_assim_count.rst>`__ of each type from an obs_seq file

*integrate_model* 
   used with subroutine callable models to `test running as a standalone program <../assimilation_code/programs/integrate_model/integrate_model.rst>`__ 


