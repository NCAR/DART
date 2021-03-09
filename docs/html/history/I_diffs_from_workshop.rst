DART Iceland revisions
======================

The DART Iceland release (23 Nov 2005) is expected to be the last update that makes major modifications to the user
interfaces to DART and the DART interfaces to models or observations. The primary purpose of this release is to improve
the way that DART handles observation sequences. The goal is to provide a more modular fashion for adding new types of
observations with potentially complex forward operators and arbitrary amounts and types of metadata that are needed to
define a particular observation. Limited backward compatibility has been implemented for existing observation sequence
files, however, it is strongly recommended that users switch to the new format as quickly as possible.

Improvements
------------

The changes to the observation sequences impact large portions of the DART code. In addition, Iceland also implements a
number of other improvements and additions to DART and fixes several known bugs from earlier releases. Highlights of the
changes from the workshop and Hawaii releases are:

#. Namelist error detection.
   Enhanced error checking capabilities for DART namelist reads are included. The input.nml must now contain an entry
   for each namelist that is read by a program, even if none of the values in that namelist entry are set. For instance,
   if a program that uses the assim_tools_mod is run, input.nml MUST contain an entry &assim_tools_nml. Any namelist
   values to be set must follow this and the list is terminated by a /. If no entries are to be set in a given namelist,
   it must still be terminated by a line that contains only a slash (/). If a variable that is NOT in a given namelist
   is found in input.nml, a fatal error occurs. Failing to terminate a namelist with slash is also fatal. Tools to
   support this improved namelist error detection are implemented in utilities_mod and can be used to check for errors
   in other applications.
#. Automatic detection of input file format.
   Restart and observation sequence files can be written in binary or ascii as in previous releases. However, file
   formats for reads are now detected automatically. The namelist entries 'read_binary_restart_files' in assim_model_nml
   and 'read_binary_obs_sequence' in obs_sequence_nml have been removed.
#. Error corrected in threed_sphere/location_mod.
   An error in the Hawaii and Workshop releases led to the possibility of erroneous computations of distance in the
   threed_sphere/location_mod. The error only occurred if the namelist entry 'approximate_distance' was set to .true.
   (the default was .false.) and errors were generally small and confined to observations located in polar regions. The
   result could be that localization was not applied properly to polar observations.
#. Support for reduced precision real computation.
   For some large model applications, using reduced precision reals can enhance performance. The definition of
   'real(r8)' in types_mod can now be modified to support different real precision. Some computations in the
   time_manager can fail with reduced precision, so all time computations now use an independent fixed precision.
#. Quality control definition change.
   values written to quality control (qc) fields by perfect_model_obs and filter have been modified. If a prior
   (posterior) forward operator computation fails, qc for that observation is incremented by 1000 (1000000). If the
   prior (posterior) observation lies outside the outlier_threshold, the qc for that observation is incremented by 100
   (400).
#. Added obs_sequence file header.
   Observation sequence files (obs_sequence files) now include a metadata header that includes a list of the names of
   each observation type that might be used in the file along with a mapping to an integer for each name in the file.
   Old format files without a header are automatically detected and a default mapping of integer indices for observation
   types to observation type names is used. This default mapping works for most low-order model, bgrid and CAM
   obs_sequence files. It is strongly recommended that legacy obs_sequence files be converted to the new format by
   inserting a header block with the appropriate names and integer mappings.
#. Interactive creation input for obs_sequence files.
   A standard method for creating obs_sequence files is to redirect a text file to standard input and use the
   interactive creation facility in create_obs_sequence. Old input files for this procedure may no longer work correctly
   and may need to be updated. The new interactive creation interface allows the name of observations to be used, or an
   integer index. However, the integer index is now defined by a table in the obs_kind_mod and may change dynamically.
   Users should understand this procedure.
#. obs_def_nml moved to obs_kind_nml.
   The obs_def_nml in previous releases had entries to choose what observation types to assimilate or evaluate. These
   entries have been moved to obs_kind_nml and obs_def_nml no longer exists.
#. Preprocessor functionality has changed.
   While the preprocess program still needs to be run when setting up DART experiments, it now works differently.
   Previously, the preprocessor read obs_def_mod.F90 and created obs_def_mod.f90 in a fashion very similar to the
   standard CPP preprocessor. The new preprocessor does not look like CPP. It reads input from DEFAULT_obs_kind_mod.F90,
   DEFAULT_obs_def_mod.F90, and from any requested special obs_def modules and creates an obs_kind_mod.f90 and an
   obs_def_mod.f90. Documentation of this new functionality is available in tutorial section 21 and in the html files
   for preprocess, obs_kind_mod, and obs_def_mod.
#. Improved observation space diagnostics interfaces.
   The observation space diagnostics program and associated diagnostic tools are now all located in the diagnostics
   directory. The interfaces and controls have been modified and are described in html documentation and in the tutorial
   section 18. Consistent interfaces for low-order models and large three-dimensional atmospheric models have been
   provided.

Changes
-------

A summary list of changes occurring in each DART directory follows:

+---------------------------------+-----------------------------------------------------------------------------------+
| assim_model                     | Uses digits12 for netcdf time computations allowing reduced precision models;     |
|                                 | automatic detection of restart file input format, read_binary_restart removed     |
|                                 | from namelist, namelist error checking.                                           |
+---------------------------------+-----------------------------------------------------------------------------------+
| assim_tools                     | Namelist error checking added. Order in which computation of ensemble prior and   |
|                                 | observation to state distance are computed has been changed in attempt to         |
|                                 | increase efficiency (should have only round-off level impact on results). Added   |
|                                 | ability to print out regional regression coefficients when using the parallel     |
|                                 | implementation.                                                                   |
+---------------------------------+-----------------------------------------------------------------------------------+
| common: types_mod.f90           | Added digits12 real precision kind for doing times with required precision when   |
|                                 | rest of assimilation is done with reduced real precision, added gravity to        |
|                                 | constants list.                                                                   |
+---------------------------------+-----------------------------------------------------------------------------------+
| converters                      | These WRF specific routines have been significantly modified to use the updated   |
|                                 | release of WRF. See the files and the WRF specific documentation for details.     |
+---------------------------------+-----------------------------------------------------------------------------------+
| cov_cutoff                      | added optional argument localization_override to comp_cov_factor to override the  |
|                                 | namelist selection of the localization_type for a single call. Changed to more    |
|                                 | computationally efficient method for computing value of Gaspari Cohn function.    |
|                                 | Namelist error checking.                                                          |
+---------------------------------+-----------------------------------------------------------------------------------+
| diagnostics                     | Observation space diagnostics are now included in this directory. There are       |
|                                 | directories for oned models and threed_sphere models, each contains an            |
|                                 | obs_diag.f90 program along with its namelist and documentation. The matlab        |
|                                 | directory contains matlab scripts for plotting the results of observation space   |
|                                 | diagnostics.                                                                      |
+---------------------------------+-----------------------------------------------------------------------------------+
| ensemble_manager                | Includes commented block needed to write out ensemble mean for WRF boundary       |
|                                 | forcing computations. Namelist error checking.                                    |
+---------------------------------+-----------------------------------------------------------------------------------+
| filter                          | Incorporated new namelist error checking, modified calls to read_obs_seq_header   |
|                                 | to support automatic file format detection, changed to new qc values (see summary |
|                                 | above). Namelist error checking.                                                  |
+---------------------------------+-----------------------------------------------------------------------------------+
| integrate_model                 | Namelist error checking.                                                          |
+---------------------------------+-----------------------------------------------------------------------------------+
| location/threed_sphere          | Added 5 VERTIS***\* variables for describing vertical location kinds. Corrected   |
|                                 | error in table lookup for approximate computation of cos and sin by doubling      |
|                                 | range of lookup table. Added public logical functions vert_is_undef and           |
|                                 | vert_is_surface. Improved menu for interactive definition of locations. Namelist  |
|                                 | error checking.                                                                   |
+---------------------------------+-----------------------------------------------------------------------------------+
| matlab                          | Minor modifications to several scripts.                                           |
+---------------------------------+-----------------------------------------------------------------------------------+
| mkmf                            | Templates cleaned up and templates for additional platforms added.                |
+---------------------------------+-----------------------------------------------------------------------------------+
| models                          | All with namelists have namelist error detection.                                 |
+---------------------------------+-----------------------------------------------------------------------------------+
| models/bgrid_solo               | Use new generic kind definitions to decide how to interpolate observations.       |
+---------------------------------+-----------------------------------------------------------------------------------+
| models/lorenz_04                | Added nc_read_model_vars to read in netcdf file format.                           |
+---------------------------------+-----------------------------------------------------------------------------------+
| ncep_obs                        | The code from NCEP to read bufr files has been added to the directory. This is    |
|                                 | not technically part of DART but is required as a first phase for BUFR file       |
|                                 | translation. Program create_real_obs has been generalized to read in portions of  |
|                                 | days if required and to use the new obs_kind and obs_def modules and the          |
|                                 | obs_def_reanalysis_bufr_mod.f90 to include much more detailed descriptions of the |
|                                 | types of observations. The obs_diag programs have been moved to the diagnostics   |
|                                 | directory. The matlab diagnostic routines have also been moved to the diagnostics |
|                                 | directory and generalized.                                                        |
+---------------------------------+-----------------------------------------------------------------------------------+
| DEFAULT_obs_def_mod.F90         | Replaces obs_def_mod.f90, preprocessed to create obs_def_mod.f90. No longer has a |
|                                 | namelist (previous namelist moved to obs_kind_nml). Function get_obs_name returns |
|                                 | the name string for an observation kind given the integer kind index. Routine     |
|                                 | set_obs_def_key added to set the value of the integer key associated with an      |
|                                 | obs_def_type. Provides default mapping for obs_sequence files in the old format   |
|                                 | that do not have a header table mapping indices to obs_kind strings.              |
+---------------------------------+-----------------------------------------------------------------------------------+
| obs_def_dew_point_mod.f90       | New module for doing dew point forward operators.                                 |
+---------------------------------+-----------------------------------------------------------------------------------+
| obs_def_metar_mod.f90           | New module for doing surface observation forward operators.                       |
+---------------------------------+-----------------------------------------------------------------------------------+
| obs_def_radar_mod.f90           | Revised version of radar forward operator module that works with                  |
|                                 | DEFAULT_obs_def_mod.F90.                                                          |
+---------------------------------+-----------------------------------------------------------------------------------+
| obs_def_1d_state_mod.f90        | Computes forward operators for interpolation and integrals of low-order models    |
|                                 | with a single state variable type on a cyclic domain.                             |
+---------------------------------+-----------------------------------------------------------------------------------+
| obs_def_reanalysis_bufr_mod.f90 | Computes forward operators for all types of observations available in the         |
|                                 | reanalysis BUFR files.                                                            |
+---------------------------------+-----------------------------------------------------------------------------------+
| DEFAULT_obs_kind_mod.F90        | Replaces obs_kind_mod.f90, preprocessed to create obs_kind_mod.f90. Includes new  |
|                                 | 'generic' kind definitions list with associated integers. Each observation kind   |
|                                 | must be associated with one of these generic kinds. Now has namelist to define    |
|                                 | what observation kinds are being assimilated or evaluated plus new namelist error |
|                                 | checking. Provides new interfaces to get information about obs_kind:              |
|                                 | get_obs_kind_name returns the observation name string given a kind index;         |
|                                 | get_obs_kind_index does the inverse, assimilate_this_obs_kind and                 |
|                                 | evaluate_this_obs_kind return true if this observation index is one that is to be |
|                                 | used in this way; get_obs_kind_var_type returns the generic kind associated with  |
|                                 | an observation type, get_kind_from_menu offers interactive creation capability.   |
+---------------------------------+-----------------------------------------------------------------------------------+
| obs_sequence_mod.f90            | obs_sequence files now have a header that maps from obs_kind indices to a string  |
|                                 | that uniquely identifies the observation kind. Automatic detection of             |
|                                 | obs_sequence file formats and old format without header. Automatic namelist error |
|                                 | detection and removal of read_binary_obs_sequence from namelist. Removal of code  |
|                                 | for WRF radar observations.                                                       |
+---------------------------------+-----------------------------------------------------------------------------------+
| perfect_model_obs.f90           | Uses revised calls to read_obs_seq_header and read_obs_seq to use automatic file  |
|                                 | format detection. Automatic namelist error detection.                             |
+---------------------------------+-----------------------------------------------------------------------------------+
| preprocess.f90                  | Now preprocesses the DEFAULT_obs_kind_mod.F90 and DEFAULT_obs_def_mod.F90 and     |
|                                 | inputs information from ``obs_def_???_mod.f90`` files such as                     |
|                                 | obs_def_reanalysis_bufr_mod. Looks for fixed format text strings in the input     |
|                                 | files to determine what sections of code to extract and where to insert them in   |
|                                 | the DEFAULT files. Namelist includes the names of the two input DEFAULT files,    |
|                                 | the names of the output preprocessed files (normally obs_def_mod.f90 and          |
|                                 | obs_kind_mod.f90 in the appropriate directories) and a list of all the            |
|                                 | ``obs_def_???_mod.f90`` files that are to be incorporated.                        |
+---------------------------------+-----------------------------------------------------------------------------------+
| reg_factor_mod.f90              | Automatic namelist error detection.                                               |
+---------------------------------+-----------------------------------------------------------------------------------+
| shell_scripts                   | Several new scripts for managing files and cleaning up DART directories have been |
|                                 | added. Significant modifications have been made to the platform specific scripts  |
|                                 | advance_ens, assim_filter, and filter_server. Versions for additional platforms   |
|                                 | have been added.                                                                  |
+---------------------------------+-----------------------------------------------------------------------------------+
| time_manager_mod.f90            | Use of digits12 precision for real computations allows reduced precision to be    |
|                                 | used for rest of dart. Optional error return added to read_time to support        |
|                                 | automatic file detection for dart state vector files.                             |
+---------------------------------+-----------------------------------------------------------------------------------+
| tutorial                        | The workshop tutorial scripts have been updated to correct several errors and to  |
|                                 | be consistent with the preprocessing changes. Section 21 has been added to        |
|                                 | describe the new obs_sequence implementation.                                     |
+---------------------------------+-----------------------------------------------------------------------------------+
| utilities_mod.f90               | Namelist error detection added.                                                   |
+---------------------------------+-----------------------------------------------------------------------------------+

Future enhancements / work
--------------------------

-  Extend PBL_1d support for all matlab scripts.
   currently only supported by the observation-space diagnostics and a crude implementation for 'plot_total_err'.
-  Unify the machine-specific scripts to handle PBS, LSF and interactive submission in one script.
-  Incorporate support for 'null_model'.
   A useful exercise to test many facets of the DART code without a chaotic model. Should provide capability to perform
   regression testing of DART infrasturcture.
-  Improve netcdf error messages.
   Will incorporate an additional argument to the 'check' routine to append a string to the netCDF error library string.
