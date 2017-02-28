# DART software - Copyright UCAR. This open source software is provided
# by UCAR, "as is", without charge, subject to all terms of use at
# http://www.image.ucar.edu/DAReS/DART/DART_download
#
# DART $Id$

++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
+ Jan 13th 2017 :: rma_fixed_filenames merge changes Revision: 10902
++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

Specific namelist changes include:

1. ) Earlier versions of the RMA branch code supported both direct NetCDF 
     reads/writes and the original binary/ascii DART format restart files.  
     As of the next update DART format files are no longer supported.  All 
     I/O is NetCDF only.  If your model does not use NetCDF you will still 
     need a model_to_dart and dart_to_model converter; otherwise all DART 
     programs read the model's NetCDF files directly.  The namelist options 
     related to selecting direct netcdf I/O have been removed.

2. ) Diagnostic and state space data (such as inflation, mean and sd 
     information) that were previously stored in {Prior,Posterior}_Diag.nc 
     are now broken up into multiple files and have fixed filenames. This 
     decreases the IO time for diagnostic output and reduces the number of 
     namelist options.

3. ) There is no longer support for observation space inflation 
     (i.e. inf_flavor = 1).  Contact us at dart@ucar.edu if you have an 
     interest in using this option.

------------------------------------------------------------------------------
Changes to the filter_nml are :
------------------------------------------------------------------------------
* restart_in_file_name      -- has been replaced with input_restart_file_list.
                               The namelist must contain one or more file names, 
                               each of which is a textfile containing a list of N 
                               NetCDF restart files, one per line for each ensemble member.
                               For models with multiple domains (e.g. nested WRF or 
                               CLM) you must specify a listfile for each domain.

* restart_out_file_name      -- has been replaced with output_restart_file_list.
                                Same format as input_restart_file_list.

* inf_in_file_name           -- REMOVED, now have fixed names of the form 
                                input_{prior,posterior}inf_{mean,sd}.nc

* inf_out_file_name          -- REMOVED, now have fixed names of the form 
                                output_{prior,posterior}inf_{mean,sd}.nc.

* inf_diag_filename          -- REMOVED

* inf_output_restart         -- REMOVED, inflation restarts will be written 
                                out if inflation is turned on

* output_inflation           -- REMOVED, inflation diagnostic files will be written 
                                if inflation is turned on

* stages_to_write            -- There is more control over what state data 
                                to write.  Options are at stages : 
                                'input', 'preassim', postassim', 'output'.  
                                Stages preassim and postassim will output 
                                state data originally contained within the 
                                copies of Prior_Diag.nc and Posterior_Diag.nc.
                                See rma_doc/rma.html for details on the 
                                filename conventions. For example, running 
                                filter with prior inflation enabled with 
                                stage 'preassim' enabled will produce files 
                                with names:
                                   preassim_member_####.nc
                                   preassim_{mean,sd}.nc
                                   preassim_priorinf_{mean,sd}.nc

* write_all_stages_at_end    -- important for large models - all output file 
                                I/O is deferred until the end of filter, but 
                                will use more memory to store the data.  More
                                detailed info is in rma_doc/rma.html

* output_restart_mean        -- renamed output_mean

* output_restart             -- renamed output_restarts

* direct_netcdf_{read,write} -- REMOVED, always true

* restart_list_file          -- renamed input_restart_file_list

* single_restart_file_in     -- renamed single_file_in

* single_restart_file_out    -- renamed single_file_out

* add_domain_extension       -- REMOVED

* use_restart_list           -- REMOVED

* overwrite_state_input      -- REMOVED, equivalent functionality can be set 
                                with single_restart_file_in = single_restart_file_out

------------------------------------------------------------------------------
Changes to the perfect_model_obs_nml are :
------------------------------------------------------------------------------
* restart_in_filename        -- renamed restart_in_file_names takes a NetCDF 
                                file. For multiple domains you can specify a 
                                list.

* direct_netcdf_{read,write} -- REMOVED, always true

------------------------------------------------------------------------------
Changes to the state_space_diag_nml are :
------------------------------------------------------------------------------
* single_file               -- REMOVED, diagnostic files are now controlled 
                               in filter_nml with stages_to_write

* make_diagnostic_files     -- REMOVED, no longer produce original
                               Prior_Diag.nc and Posterior_Diag.nc

* netCDF_large_file_support -- REMOVED, always true

------------------------------------------------------------------------------
Changes to the state_vector_io_nml are :
------------------------------------------------------------------------------
* write_binary_restart_files -- REMOVED

------------------------------------------------------------------------------
Changes to the ensemble_manager_nml are :
------------------------------------------------------------------------------
* flag_unneeded_transposes -- REMOVED

------------------------------------------------------------------------------
Changes to the integrate_model_nml are :
------------------------------------------------------------------------------
* advance_restart_format -- REMOVED, only supporting NetCDF format.

------------------------------------------------------------------------------
Scripting with CESM
------------------------------------------------------------------------------
See models/cam-fv/scripts_cesm1_5/assimilate.csh for an example of how to 
handle the new filename conventions.

(To help find things:  input_priorinf_mean output_priorinf_mean )
{in,out}put_{prior,post}inf_{mean,sd}.nc   ARE in use;
    Search for stage_metadata%filenames turned up
    interface set_file_metadata
       module procedure set_explicit_file_metadata
       module procedure set_stage_file_metadata

      ! stage_name is {input,preassim,postassim,output}
      ! base_name  is {mean,sd,{prior,post}inf_{mean,sd}} from filter/filter_mod.f90.
      write(string1,'(A,''.nc'')') trim(stage_name)//'_'//trim(base_name)
      file_info%stage_metadata%filenames(my_copy,1) = trim(string1)

    This shows where inflation file names are defined.
      > grep -I set_file_metadata */*.f90 | grep inf
    filter/filter_mod.f90:
       call set_file_metadata(file_info, PRIOR_INF_MEAN, stage, 'priorinf_mean', 'prior inflation mean')
       call set_file_metadata(file_info, PRIOR_INF_SD,   stage, 'priorinf_sd',   'prior inflation sd')
       call set_file_metadata(file_info, POST_INF_MEAN,  stage, 'postinf_mean',  'posterior inflation mean')
       call set_file_metadata(file_info, POST_INF_SD,    stage, 'postinf_sd',    'posterior inflation sd')

    subroutine set_member_file_metadata(file_info, ens_size, my_copy_start)
       call set_file_metadata(file_info, icopy, stage_name, base_name, desc, offset)

    subroutine set_stage_file_metadata(file_info, copy_number, stage, base_name, desc, offset)
       write(string1,'(A,''.nc'')') trim(stage_name)//'_'//trim(base_name)

    subroutine set_explicit_file_metadata(file_info, cnum, fnames, desc)
       file_info%stage_metadata%filenames(cnum,idom)        = trim(fnames(idom))
       file_info%stage_metadata%file_description(cnum,idom) = trim(string1)

    function construct_file_names(file_info, ens_size, copy, domain)
       write(construct_file_names, '(A, ''_member_'', I4.4, A, ''.nc'')') &
                           trim(file_info%root_name), copy, trim(dom_str)



Also see
   harnesses/filename_harness/files:  ENS_MEAN_COPY       PriorDiag_mean.nc


------------------------------------------------------------------------------
ADDITIONAL NOTES :
------------------------------------------------------------------------------
* currently the closest_member_tool is broken but plans on being fixed soon.
* restart_file_tool and most model_to_dart/dart_to_model programs have been
  deprecated, since DART formated restarts are no longer supported.
* some programs such as model_mod_check have not been fully tested and need
  to be exercised with the new naming conventions.



++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
+ Feb 15th 2017 :: rma_single_file merge changes Revision: 11136 
++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

Filter and PMO can now run with multiple cycles for low order models. The output
for this is only supported with single file output (members, inflation, mean, sd
are all in the same file).

Added matlab support for diagnostics format in lower order models.

------------------------------------------------------------------------------
Changes to the filter_nml are :
------------------------------------------------------------------------------

output_restart          -- RENAMED to output_members
restart_in_file_name    -- RENAMED to input_state_file_list
restart_out_file_name   -- RENAMED to output_state_file_list
single_restart_file_in  -- RENAMED to single_file_in
single_restart_file_out -- RENAMED to single_file_out

input_state_files       -- ADDED - not currently working
output_state_files      -- ADDED - not currently working

has_cycling             -- ADDED for low order models

------------------------------------------------------------------------------
Changes to the perfect_model_obs_nml are :
------------------------------------------------------------------------------

start_from_restart    -- RENAMED read_input_state_from_file
output_restart        -- RENAMED write_output_state_to_file
restart_in_file_name  -- RENAMED input_state_files
restart_out_file_name -- RENAMED output_state_files

single_file_in        -- ADDED for low order models
single_file_out       -- ADDED for low order models
has_cycling           -- ADDED for low order models

++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
+ Feb 15th 2017 :: filter updates  Revision: 11160 
++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

The postassim diagnostics file was being incorrectly written after
posterior inflation was applied.  It is now written immediately after
the assimilation update, and then posterior inflation, if enabled,
is applied.

Sampling Error Correction now reads data from a single netcdf file
for any ensemble size.  To add other sizes, a program can generate
any ensemble size and append it to this file.  The default file is
currently in system_simulation:

system_simulation/work/sampling_error_correction_table.nc

Filter and PMO no longer need the "has_cycling" flag.

------------------------------------------------------------------------------
Changes to the filter_nml are :
------------------------------------------------------------------------------

has_cycling             -- REMOVED for low order models

------------------------------------------------------------------------------
Changes to the perfect_model_obs_nml are :
------------------------------------------------------------------------------

has_cycling           -- REMOVED for low order models


++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
+ Month Day 2017 :: change $Revision$
++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

# <next few lines under version control, do not edit>
# $URL$
# $Revision$
# $Date$
