# DART software - Copyright UCAR. This open source software is provided
# by UCAR, "as is", without charge, subject to all terms of use at
# http://www.image.ucar.edu/DAReS/DART/DART_download
#
# DART $Id$

++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
+ Jan 13th 2017 :: rma_fixed_filenames merge changes                         +
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

------------------------------------------------------------------------------
ADDITIONAL NOTES :
------------------------------------------------------------------------------
* currently the closest_member_tool is broken but plans on being fixed soon.
* restart_file_tool and most model_to_dart/dart_to_model programs have been
  deprecated, since DART formated restarts are no longer supported.
* some programs such as model_mod_check have not been fully tested and need
  to be exercised with the new naming conventions.

# <next few lines under version control, do not edit>
# $URL$
# $Revision$
# $Date$
