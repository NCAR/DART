
![DARTlogo](docs/images/Dartboard7.png)

----

DART software - Copyright UCAR. This open source software is provided     
by UCAR, "as is", without charge, subject to all terms of use at     
http://www.image.ucar.edu/DAReS/DART/DART_download     

----

This file documents the most user-visible changes to the DART code.
It is not intended to document every change, but instead is intended
to inform people what features are now available or have been removed.
Detailed changes are always available through the version control framework.

DART now uses Git for version control but has preserved the revision history
from when subversion (and CVS before that) was used. The previous revision numbers
can be related to git hashes by searching the output of `git log`

```
0[1011] machine:dartGIT % git log > full_git_log.txt
```

A reminder that since many files were moved or renamed, the best way to get the complete
log is to use `git log --follow` for information on individual files.

## The changes are now listed with the most recent at the top.

------------------------------------------------------------------------------
## Oct 29 2020 :: radiance support, MPAS, obs converters           Tag: v9.9.0

- Use RTTOV (Radiative Transfer for TOVS) routines to support radiance assimilation.
    - [Introduction to DART support for RTTOV](https://dart.ucar.edu/pages/Radiance_support.html)
    - WRF, MPAS, and CAM-FV model interfaces now support radiance assimilation.
    - Added GOES 16-19 ABI converter

- *_NOTE_*: The `build_templates/mkmf.template` file has been removed from version control. You must now explicitly copy the best example `mkmf.template` into place before compiling. If there is no `mkmf.template` when you try to build, an error message is displayed.

- MPAS regional configurations now supported.

- Converted CHANGELOG to a markdown document, put newest content at top.

- Converted many HTML documents to markdown
   - renamed `observations/obs_converters/observations.html` to `observations/obs_converters/README.md` for example.

- [Updated Publications](https://dart.ucar.edu/pages/Publications.html)

- declare hexadecimal constants according to the Fortran standard.

- GSI2DART converter updated - Thanks to **Craig Schwartz** & **Jamie Bresch**.
    
- The WRF-DART tutorial has been rewritten as `models/wrf/tutorial/README.md`

- Hydro-DART (AKA wrf-hydro/DART) has been updated to be Manhattan-compliant.
    - also support masked bucket
    - added perturbed forcing capability

- The support for POP and CESM2 has been implemented and documented.

- `obs_diag` now correctly handles the special case when the observation is properly assimilated or evaluated but the posterior forward operator fails. The posterior DART QC in the `obs_diag_output.nc` should be a '2', not a '4'. The prior DART QC value in obs_diag_output.nc can still be a 7 if need be.

- `obs_def_tower_mod.f90` was refactored into `obs_def_tower_mod.f90` and `obs_def_land_mod.f90`.

- WRF-Chem/DART documentation and datasets have been updated for Manhattan.
  [WRF-Chem contact information.](https://dart.ucar.edu/pages/Models.html#wrf-chem)

- Fixed bug in `obs_seq_to_netcdf` to correctly append to existing netCDF files.
    
- Support absolute humidity observations - Thanks to **Michael Ying**.

- `DEFAULT_obs_kind_mod.F90` has many added quantities.

- new observation converters including (but not limited to):
   - absolute humidity
   - streamflow observations from the Mexican water agency
   - streamflow observations from the USGS
   - total water storage observations from GRACE
   - radiance observations from GOES

- the following forward operator modules are either new or modified:
   - (M) `observations/forward_operators/DEFAULT_obs_def_mod.F90`
   - (M) `observations/forward_operators/obs_def_GRACE_mod.f90`
   - (A) `observations/forward_operators/obs_def_abs_humidity_mod.f90`
   - (M) `observations/forward_operators/obs_def_altimeter_mod.f90`
   - (A) `observations/forward_operators/obs_def_land_mod.f90`
   - (A) `observations/forward_operators/obs_def_mesonet_mod.f90`
   - (M) `observations/forward_operators/obs_def_oxygen_ion_density_mod.f90`
   - (M) `observations/forward_operators/obs_def_reanalysis_bufr_mod.f90`
   - (M) `observations/forward_operators/obs_def_rel_humidity_mod.f90`
   - (A) `observations/forward_operators/obs_def_rttov_mod.f90`
   - (A) `observations/forward_operators/obs_def_streamflow_mod.f90`
   - (M) `observations/forward_operators/obs_def_tower_mod.f90`
   - (M) `observations/forward_operators/obs_def_upper_atm_mod.f90`
   - (A) `observations/forward_operators/rttov_sensor_db.csv`

- `fill_inflation_restart` now correctly creates inflation values for all variables in the DART state, regardless of the setting of the `no update` metadata.
    
- GITM is now fully Manhattan compliant.

- fix bug in madis rawin converter

- avoid computing posterior inflation if using the 'relaxation to prior spread' inflation option -- Thanks to **Craig Schwartz**.

- add additional reporting options to the `obs_assim_count` utility


------------------------------------------------------------------------------
## Nov 20 2019 :: FESOM,NOAH-MP model support, better testing       Tag: v9.8.0

- first release entirely from GIT

- fixed bug in `fill_inflation_restart` tool which used the prior inflation mean
  and sd for both prior and posterior inflation files.  now correctly uses the
  posterior mean/sd if requested.

- fixed a typo in the location test script that prevented it from running

- additional functionality in the quad interpolation code, now supports grids
  which start at 90 (north) and end at -90 (south).

- if possible, send shorter MPI messages.  improves performance on some platforms
  and MPI implementations.

- add explicit call to `initalize_utilities()` where it was missing in a couple of
  the WRF utility routines.

- added an example of how to use a namelist to the `text_to_obs.f90` observation
  converter program.

- Removing the clamping messages in `clamp_variable()` of clamped values

- changed argument names using reserved keywords.
  - `state_vector_io_mod:read_state() 'time' to 'model_time'`
  - `random_seq_mod:random_gamma() 'shape' to 'rshape', 'scale' to 'rscale'`
  - `random_seq_mod:random_inverse_gamma() 'shape' to 'rshape', 'scale' to 'rscale'`
  - `obs_def_mod:init_obs_def() 'kind' to 'obkind', 'time' to 'obtime'`
  - `obs_utilities_mod: 'start' to 'varstart', 'count' to 'varcount'`

- The **FESOM** model is now Manhattan-ready. Thanks to **Ali Aydogdu**

- The **noah** model is now Manhattan-ready and may be used with NOAH-MP.

- bugfixed references to the `documentation` directory that was renamed
  `docs` to comply with GitHub Pages.

- improved `test_dart.csh` functionality.


------------------------------------------------------------------------------
## Apr 30 2019 :: cam-fv refactor, posteriors optional, QC 8    Revision: 13138

- The CAM Finite Volume (**cam-fv**)  `model_mod.f90` has undergone substantial
  refactoring to improve simplicity and remove code for unsupported CAM variants
  while also supporting WACCM and WACCM-X.  Namelist changes will be required.

- **cam-fv** setup and scripting support added for CESM 2.1, including advanced
  archiving and compression

- fix for WRF's wind direction vectors when using the Polar Stereographic
  map projection.  Thanks to **Kevin Manning** for the fix.

- Add filter namelist option to avoid calling the posterior forward operators
  and to not create those copies in the `obs_seq.final` file.

- Use less memory if writing ensemble member values into the `obs_seq.final` file.

- added a DART QC of 8 for failed vertical conversions

- updated Matlab scripts support QC=8 and no posterior in obs sequence files.

- sampling error correction table now has all ensemble sizes between 3 and 200

- `closest_member_tool` can be compiled with other MPI targets

- `COSMIC_ELECTRON_DENSITY` has been moved from `obs_def_gps_mod.f90` to
  `obs_def_upper_atm_mod.f90`, which has new quantities for
  `ION_O_MIXING_RATIO` and `ATOMIC_H_MIXING_RATIO`

- `obs_converters/gps/convert_cosmic_ionosphere.f90` has a test dataset

- support for NAG compiler

- fixed Intel compiler bug in `lorenz_96` comparing long integers to integer loop indices

- `get_maxdist()` now a required routine all location modules

- Default routines now create a time variable as `time(time)` to allow multiple
  files to be concatenated along the unlimited dimension more easily.  Also
  conforms to the netCDF convention for coordinate dimensions.

- `obs_impact_tool` handles a continuum of values, not just discrete 0 or 1.

- `fill_inflation_restart` now produces files with names consistent with filter defaults.

- expanded functionality in `xyz_location_mod.f90`

- Removed 'slow' sorting routines from `sort_mod.f90`

- replacing some repeated native netCDF library calls with routines from
  the `netcdf_utilities_mod.f90`

- Updated dewpoint equation to avoid dividing by zero given a very unlikely
  scenario (r12832)

- More efficient implementation of adaptive inflation

- *Yongfei Zhang* and *Cecilia Bitz* added improvements to the CICE model and
  observation converters and forward operators. These changes also use the
  locations of the 'new' glade filesystem. They used CESM tag: cesm2_0_alpha06n

- Worked with Yongfei Zhang to remove prototype codes and more completely
  document observation converters and data sources for cice assimilation.

- removed `allow_missing_in_clm` flag from the `&assim_tools_nml` namelist in
  the CICE work directory.  The flag moved to a different namelist and the
  CICE model doesn't care about it.

- increased the maximum number of input files to `obs_diag` from 100 to 10000.

- Updated the `developer_tests` to include more cases.

- Updated `oned/obs_diag.f90` to support `obs_seq.out` files.

- Better error and informational messages in various routines.


------------------------------------------------------------------------------
## Aug 03 2018 :: performance fix for distributed mean          Revision: 12758

- Important performance fix if model does vertical conversion for localization.
  Results were not wrong but performance was poor if `distribute_mean = .true.`
  was selected in the `&assim_tools_nml` namelist.

  Now distributing the mean runs in close to the non-distributed time and uses
  much less memory for large models. This only impacts models which do a vertical
  conversion of either the observation or state vertical coordinate for localization
  AND which set `&assim_tools_nml :: distribute_mean = .true.` to use less memory.

  When using a distributed mean `convert_all_obs_verticals_first = .true.` should
  be set.  If your observations will impact most of the model state, then
  `convert_all_state_verticals_first = .true.` can also be set.


------------------------------------------------------------------------------
## Jun 18 2018 :: CAM/CESM 2.0, DART QC 8, closest_member_tool  Revision: 12682

- Support for **cam-fv** assimilations in the CESM 2.0 release.  See
  documentation in `models/cam-fv/doc/README_cam-fv` for details.

- `obs_diag` and matlab scripts updated to report statistics on DART QC 8,
  observation failed vertical conversion

- Updates to fix minor problems with the new WRF scripts

- Added the `inf_sd_max_change` namelist item to all `input.nml` files for
  the enhanced inflation option

- Revival of the `closest_member_tool`, which now runs in parallel on
  all ensemble members at one time.  This tool can be used as a template
  for any other tools which need to process something for all ensemble
  members in parallel.

- Revival of the `fill_inflation_restart` tool as a Fortran 90 program.
  Using `ncap2` is still possible, but if the correct version is not
  installed or available this tool can be used.

- Added more functions to the `netcdf_utilities_mod.f90`


------------------------------------------------------------------------------
## May 21 2018 :: enhanced inflation option, scripting          Revision: 12591

- Enhanced inflation algorithm added. See the `filter_mod.html` for new
documentation on this option.

- Updated WRF scripts for the Manhattan release.

- `obs_diag` reports statistics on DART QC 8, observation failed vertical
conversion.  Matlab scripts also updated to support QC 8.

- New parallel conversion scripts for GPS Radio Occultation observations and
  NCEP prepbufr conversions.

- Further updates to documentation files to change KIND to QTY or Quantity.

- Documented required changes when moving from the Lanai/Classic release to
Manhattan in `documentation/html/Manhattan_diffs_from_Lanai.html`

- Expanded the routines in the `netcdf_utilities_mod.f90`

- Add an ensemble handle parameter to the 6 ensemble manager routines
where it was missing.

- The `advance_time` program can read/generate CESM format time strings
(YYYY-MM-DD-SSSSS).

- Fixed a bug in the netcdf read routines that under certain circumstances
could report an array was using the unlimited dimension incorrectly.

- Removed the option to try to bitwise reproduce Lanai results; due to the
number of changes this is no longer possible.

- Minor bug fixes to the (seldom used) perturb routines in the **WRF**
and **mpas_atm** `model_mod.f90` files.  (used to add gaussian noise to a
single model state to generate an ensemble; this is never the recommended
method of starting a new experiment but the code remains for testing
purposes.)

- Several remaining model-specific `model_mod_check` programs were
removed in favor of a single common program source file.

- Keep `filter_mod.dopplerfold.f90` in sync with `filter_mod.f90`,
and `assim_tools_mod.pf.f90` in sync with `assim_tools_mod.f90`.

- Removed makefiles for the obsolete `trans_time` program.


------------------------------------------------------------------------------
## Mar 01 2018 :: ROMS, MMC, PMO, mpas_atm debug, etc           Revision: 12419

- Fix a debug message in the **mpas_atm** model which might have caused a buffer
overflow crash when formatting a message for a larger ensemble size.

- Update the **ROMS** shell scripts to support PBS, SLURM, as well as LSF.
Update the ROMS model_mod html documentation.

- Update the default **cam-fv** `input.nml` to have more realistic values
for the highest observation assimilated, and for where the ramp starts
that decreases the increments at the model top. If running with a higher
model top than the default check these items carefully.

- Fixed variable type for `time` variables we create in diagnostic files

- Miscellaneous minor Bug fixes:
   - Print format wider for fractional levels in `threed_sphere` locations
   - Fixed a deallocate call at program shutdown time
   - Fixed an indexing problem computing **cam-fv** U_WIND observations if the
     observation used HEIGHT as the vertical coordinate (very unusual).
   - Fixed grid creation bug in a test program used with `model_mod_check`.
     Now uses correct spacing for grids in the x,y coordinates.
   - Fixed an allocate problem in a test interpolate routine.

- Add surface pressure to the default state list in the **wrf** `work/input.nml`

- `developer_tests/test_dart.csh` can run PMO for more models. required
updates to the `work/input.nml` in several directories 
(wrf, cm1, POP, mpas_atm) to match the current namelist.

- several `model_mod_check` programs were combined into a single version
that allows for selection of individual tests.  many of the input.nml
`models/xxx/work/input.nml` files have either had a `&model_mod_check_nml`
section added or updated to match the updated interface.

- the DART QTYs are now available via the state structure in the **wrf**
and **clm** `model_mod`s.

- support the NAG compiler better.  (contact dart@ucar.edu for more
help if you want to use this compiler.  some hand work is still needed.)

- streamlined the debug output from the `state_structure_info()` call to
avoid replicating information that was the same for all variables.

- minor formatting change to the dart log file output for the list of
observation types being assimilated, evaluated, and using precomputed
forward operators.

- fixed an uninitialized variable in the BGRID model code in a routine
that isn't normally used.

- Updated the `threed_sphere` location module documentation with some usage
notes about issues commonly encountered.

- Fixed an incorrect test when printing out a log message describing if the
inflation would be variance-adaptive or not.

- Change the location of the POP MDT reference file to be relative to the
current run directory and not an absolute file location on cheyenne.

- Make the ROMS, CM1, and POP model_mod log namelist information to the
namelist log file and not the main DART log file.

- Updated several html documentation files, including the `template/model_mod.html`
which describes the current model_mod required interfaces.

- Updated the instructions for the GSI to DART obs converter to suggest some
needed compiler flags in certain cases.

- Updated the location module test programs.


------------------------------------------------------------------------------
## Dec 01 2017 :: ROMS scripting, debugging aids                Revision: 12166  

- Added an option to the ROMS model scripting to advance the model ensemble
members in parallel using a job array.

- Updated the DART_LAB Matlab GUIs to log a history of the settings and results.

- Added a debug option to the filter namelist, `write_obs_every_cycle`,
to output the full `obs_seq.final` during each cycle of filter.  
(Very slow - use only when debugging a filter crash.)

- Allow the test grid in `model_mod_check` to cross the prime meridian for testing
longitude interpolation in grids that cross the 360/0 line.


------------------------------------------------------------------------------
## Nov 22 2017 :: minor updates for DA challenge files          Revision: 12144

- added `obs_seq.in.power` to the Lorenz 96 directory

- added new obs types to the workshop version of the `input.nml` assimilation list


------------------------------------------------------------------------------
## Nov 21 2017 :: 1D obs_diag fix, 1D power forward operator    Revision: 12138

- fixed a bad URL reference in tutorial section 18

- fixed a crash with the 1D version of the observation diagnostics program
  when including identity observations.

- all models with a `workshop_setup.csh` now build the same set of programs.
  (some/most did not build obs_diag - which is used in the tutorial)

- added a 1D obs-to-a-power forward operator.

- updates to the matlab plotting routines for NetCDF observation formats

- World Ocean Database (WOD) converter supports partial year conversions
  and 2013 file formats.


------------------------------------------------------------------------------
## Oct 17 2017 :: mpas_atm bug fix, various other updates.     Revision: 12002

- Fixed a bug in the **mpas_atm** `model_mod` that affected surface observations,
  in particular altimeter obs.  also fixed a bug in the vertical conversion
  if using 'scale height' as the vertical localization type.

- Fixed a bug in the **cam-fv** `model_mod` which might have excluded observations
  with a vertical coordinate of height (meters) which were in fact below the
  equivalent highest_obs_pressure_Pa namelist setting.  also fixed a possible
  memory leak.

- Added two new modules: `options_mod.f90` and `obs_def_utilities_mod.f90`
  this was required so we didn't have circular dependencies in our modules
  as we reused common code in more places.
  We have updated all the `path_names*` files which are in the repository.
  if you have your own path_names files you may need to add these new modules
  to your path lists.
  - `assimilation_code/modules/utilities/options_mod.f90`
  - `observations/forward_operators/obs_def_utilities_mod.f90`

- Removed `QTY_SURFACE_TEMPERATURE` from the default obs quantities list
  and added `QTY_2M_SPECIFIC_HUMIDITY`.  `QTY_2M_TEMPERATURE` exists for
  atmospheric models, and `QTY_SKIN_TEMPERATURE` and `QTY_SOIL_TEMPERATURE`
  exist for other models.  if you were using `QTY_SURFACE_TEMPERATURE`
  please replace it with the corresponding other temperature quantity.

- Updated and improved the observation converter for ionospheric observations
  from the COSMIC GPS satellite.

- Updated the **cam-fv** scripts for cesm2_0_beta05.

- Updated the Matlab diagnostics documentation.  'help DART' or 'doc DART'
  will give an overview of the available Matlab diagnostics shipped with the
  dart distribution.

- Added the observation type `COSMIC_ELECTRON_DENSITY` to the `obs_def_upper_atm_mod`

- `dart_to_clm` and `clm_to_dart` were resurrected to correctly handle conversions
  for the SWE (snow water equivalent) field.

- Updated the channel and column location modules to be compatible with
  the current required interfaces.

- Updated the `model_mod_check.f90` program (most often used when porting
  DART to a new model).  there is now more control over exactly which
  tests are being run.  updated the nml and html documentation files to
  match the current code and describe the tests in more detail.

- Fixed a misleading status message in the `obs_sequence_tool` when all obs
  are excluded by the min/max lon/lat box namelist items.  the incorrect
  message blamed it on observation height instead of the bounding box.

- Added some additional debugging options to the mpi utilities module.
  if you have problems that appear to be MPI related, contact us for
  more help in enabling them.

- Improved some error messages in `location_io_mod` and `state_structure_mod`


------------------------------------------------------------------------------
## Aug 2 2017 :: single filenames, random distributions, bug fixes.  Revision: 11864

- added code to support listing input and output filenames directly in the
  namelist instead of having to go through an indirect text file.  most useful
  for programs that take a single input and output file, but works for all cases.

- bug fix in `location_io_mod.f90` that affected `obs_seq_to_netcdf` (error in adding
  vertical location types to output file).

- fix to `convert_gpsro_bufr.f90` converter (GPS obs from BUFR files) that failed
  if r8 defined to be r4.

- added draws from gamma, inverse gamma, and exponential distributions to the
  random sequence module.

- various updates to the **cam** scripts to work more smoothly with the most
  recent CIME changes and DART Manhattan updates.

- added `QTY_CWP_PATH` and `QTY_CWP_PATH_ZERO` to the default quantities list for
  the `obs_def_cwp_mod.f90` forward operator.

- improved some error messages in the diagnostic matlab scripts


------------------------------------------------------------------------------
## July 18 2017 :: bug fixes, documentation updates.  Revision: 11830

- fixed bug in `obs_impact_tool` when generating the run-time table.  specifying
  a generic quantity resulted in selecting the wrong specific obs types.

- fixed a bug that would not allow filter to start from a single ensemble member
  if `single_file_in = .true.`

- updates to HTML documentation especially for types/quantities (replacing kinds)

- updates to `input.nml` namelists, code comments, and shell scripts where
  names changed from `restart` to `state` for input and output files.


------------------------------------------------------------------------------
## July 7th 2017 :: cam-fv, mpas_atm scripts, single file i/o.  Revision: 11807

- **mpas_atm**: scripts completely revised for the Manhattan release.
  Many thanks to **Soyoung Ha** and **Ryan Torn** for the contributed code.

- **cam-fv**:  scripts and `model_mod.f90` updated for cesm2_0_beta05.

Single File I/O:

- Now we are able to run `single_file_in` and `single_file_out` with MPI.

- `single_file_io_mod.f90` has been removed and its functionality has been moved
  to `direct_netcdf_mod.f90`.

- `single_file_io_mod.f90` has been removed from all of the `path_names_*` 
   files in the repository. (Remove it from any private `path_names_*` files.)


------------------------------------------------------------------------------
## June 27rd 2017 :: CICE 5, model_mod_check, tutorial.  Revision: 11770

- Updated support for CICE5.

- Updated support for `model_mod_check` - now compatible with netCDF input files,
  input is through [input,output]_state_files namelist variable (variables renamed).

- Ensured consistency between low-order namelists and the updated DART tutorial.
  Updated documentation of many namelists. More to come.

- `location_mod`: namelist variable `maintain_original_vert` was deprecated,
  it is now removed.  You must remove it from your existing namelists or
  DART will error out immediately.

- `obs_diag`: namelist variables `rat_cri` and `input_qc_threshold` have been
  deprecated for years, they have been removed. You must remove them from
  your existing namelists or obs_diag will error out immediately.


------------------------------------------------------------------------------
## Jun 2nd 2017 :: tutorial, DART_LAB, and various updates.  Revision: 11696

- bring the DART tutorial pdf slides up to date with the current release.

- include new GUIs with adaptive inflation options in DART_LAB:
  - `oned_model_inf.m`
  - `run_lorenz_96_inf.m`

- added the **lorenz_96_2scale** model - additional kinds of
  `QTY_SMALL_SCALE_STATE` and `QTY_LARGE_SCALE_STATE` added as required.

- add useful attributes to the variables in the diagnostic files

- updates and minor bug fixes to the matlab diagnostic scripts

- updates to the default input.nmls for models

- updates to the **cam-fv** shell scripts to work with the CESM2.0 framework

- updates to the **cam-fv** `model_mod` for support of `cam-chem` variables
  Added more QUANTITIES/KINDS for chemistry species.
  Removed support for 'stand-alone' **cam** and **cam-se** (**cam-se** will be a separate 'model').

- major bug fix in the **simple_advection** `model_mod`: Fixed an error with
  the layout of the state vector.

- `obs_def_radar_mod`: Fixed a serious bug in the fall velocity forward operator.
  If the fall speed field is not in the state the test for a bad istatus from
  the interpolate() call was looking at the wrong variable and returning ok
  even if interpolate() had set bad values.

- bug fix in the **wrf** model_mod for fields which have a vertical stagger

- fix to the makefiles for the GSI2DART observation converter

- added additional netcdf and location utility routines

- various fixes to documentation and test code

- renamed `QTY_RAW_STATE_VARIABLE` to `QTY_STATE_VARIABLE`  (RAW is redundant)

- `direct_netcdf_mod`: Renamed `limit_mem` to `buffer_state_io`.
  `buffer_state_io` is now a logical that states if a variable that tells
  DART it it should read and write variables all at once or variable-by-variable.


------------------------------------------------------------------------------
## May 5th 2017 :: major changes to model_mod interfaces.  Revision: 11615

A long-awaited overhaul of the model_mod interfaces. All models which are
in our subversion repository and are supported in the Manhattan release
have been updated to match the new interfaces.  If you have model_mods with
extensive changes, our recommendation is to diff your changes with the version
you checked out and insert those changes into the new version.  The changes for
this update are unfortunately extensive.

The detailed list of changes:

`model_mod::get_state_meta_data()` is no longer passed an ensemble_handle as the
first argument.  it should not do vertical coordinate conversion.  that will be
done as a separate step by `convert_vertical_state()`

`model_mod::vert_convert` is replaced by `convert_vertical_state()` and `convert_vertical_obs()`
Any vertical conversion code that was in `get_state_meta_data` should be moved
to `convert_vertical_state()` which has access to the state vector index, so the
code should move easily.

`model_mod::query_vert_localization_coord` is no longer a required interface
`model_mod::get_close_maxdist_init` is not longer a required interface
`model_mod::get_close_obs_init` is not longer a required interface

`model_mod::get_close_obs` has a different calling convention and is split into
`get_close_obs()` and `get_close_state()`.  the close obs routine is passed both the
obs types and quantities, and the close state routine is passed both the
state quantities and the state index, for ease in vertical conversion if needed.

`model_mod::nc_write_model_vars()` is deprecated for now; it may return in a
slightly different form in the future.  

`model_mod::nc_write_model_atts()` is now a subroutine with different arguments.
it should now only write any global attributes wanted, and possibly some grid
information.  it should NOT write any of the state variables; those will be
written by DART routines.

`model_mod::get_model_size()` needs to return an `i8` (a long integer) for the size.

A new module `default_model_mod` supplies default routines for any required
interfaces that don't need to be specialized for this model.

A new module `netcdf_utilities_mod` can do some simple netcdf functions for
you and we plan to add many more over the next couple months.

`model_mod::get_model_time_step` has been replaced by `shortest_time_between_assimilations()`
since in fact it has always controlled the minimum time filter would request a model advance
and never had anything to do with the internal time step of the dynamics of the model.

We have removed `output_state_vector` from the namelist of all model_mods since
we no longer output a single 1d vector.  all i/o is now in netcdf format.

Models now have more control over when vertical conversion happens - on demand
as needed, or all up front before assimilation.

Models that were doing vertical conversion in `get_state_meta_data` should set:
```
&assim_tools_nml
   convert_all_state_verticals_first = .true.
   convert_all_obs_verticals_first = .true.

Models which were not should set:
   convert_all_state_verticals_first = .false.
   convert_all_obs_verticals_first = .true.
```

The `location_mod::vert_is_xxx()` routines have become a single `is_vertical(loc, "string")` where
string is one of: "PRESSURE", "HEIGHT", "SURFACE", "LEVEL", "UNDEFINED", "SCALE_HEIGHT"

Models doing vertical localization should add a call to `set_vertical_localization_coord()`
in their `static_init_model()` routine to tell dart what vertical coordinate system they
are expecting to convert to for vert localization

Most `path_names_xxx` files have been updated to add additional modules.  compare against
what is checked out to see the differences.

Some of the internal changes include pulling common code from the locations
modules into a `location_io_mod` which contains common functions for creating
and writing 'location' variables for any location type.

`QTY_RAW_STATE_VARIABLE` is redundant and was shortened to `QTY_STATE_VARIABLE`

Many utility programs use the `template/model_mod.f90` because they do not
depend on any model-specific functions.  this file was also updated to
match the new interfaces.

The `obs_impact` facility is enabled in the `assim_tools` namelist. you can
use the `obs_impact_tool` to construct a table which prevents one class of
observations from impacting another class of state.

Sampling Error Correction now reads the values it needs from a single
netcdf file found in `assimilation_code/programs/gen_sampling_err_table/work`.
Copy it to the same directory as where filter is running.  All ensemble
sizes which were previously in `final_full.XX` files are included, and there
is a tool to generate and append to the file any other ensemble size required.

------------------------------------------------------------------------------
## April 27th 2017 :: diagnostic file changes.  Revision: 11545

Two additional Diagnostic Files (forecast and analysis) in Filter
which can be set with the namelist option (stages_to_write)

- **input** writes out mean and sd if requested.
   - For low order models, mean and sd are only inserted into restart 
     files with a single time step.

- **forecast**
   - contains the forecast and potentially the mean and sd for the,
     this is mostly important for lower order models which cycle

- **preassim** before assimilation
     - No Inflation:       same as forecast
     - Prior Inf:          the inflated ensemble and damped prior inf
     - Post Inf:           same as forecast
     - Prior and Post Inf: the inflated ensemble and damped prior inf

- **postassim** after assimilation (before posterior infation)
     - No Inflation:       same as analysis
     - Prior Inf:          same as analysis
     - Post Inf:           assimilated ensemble and damped posterior inflation
     - Prior and Post Inf: assimilated ensemble and damped posterior inflation

- **analysis** after assimilation and before potentially update posterior inflation ensemble and updated prior inf
     - No Inflation:       assimilated ensemble
     - Prior Inf:          assimilated ensemble and updated prior inf
     - Post Inf:           post inflated ensemble and updated posterior inflation
     - Prior and Post Inf: post inflated ensemble and updated prior inf and posterior inflation

- **output**
   - a single time step of the output ensemble and potentially updated prior inf and posterior inflation


------------------------------------------------------------------------------
## Feb 15th 2017 :: filter updates.   Revision: 11160

The postassim diagnostics file was being incorrectly written after
posterior inflation was applied.  It is now written immediately after
the assimilation update, and then posterior inflation, if enabled,
is applied.

Sampling Error Correction now reads data from a single netcdf file
for any ensemble size.  To add other sizes, a program can generate
any ensemble size and append it to this file.  The default file is
currently in `system_simulation`:

`system_simulation/work/sampling_error_correction_table.nc`

Filter and PMO no longer need the "has_cycling" flag.

#### Changes to the filter_nml are :

- `has_cycling` REMOVED for low order models

#### Changes to the perfect_model_obs_nml are :

- `has_cycling` REMOVED for low order models


------------------------------------------------------------------------------
## Feb 15th 2017 :: rma_single_file merge changes.   Revision: 11136

Filter and PMO can now run with multiple cycles for low order models. The output
for this is only supported with single file output (members, inflation, mean, sd
are all in the same file).

Added matlab support for diagnostics format in lower order models.

#### Changes to the filter_nml are :

- `output_restart`          RENAMED to `output_members`
- `restart_in_file_name`    RENAMED to `input_state_file_list`
- `restart_out_file_name`   RENAMED to `output_state_file_list`
- `single_restart_file_in`  RENAMED to `single_file_in`
- `single_restart_file_out` RENAMED to `single_file_out`

- `input_state_files`       ADDED - not currently working
- `output_state_files`      ADDED - not currently working

- `has_cycling`             ADDED for low order models

#### Changes to the perfect_model_obs_nml are :

- `start_from_restart`    RENAMED `read_input_state_from_file`
- `output_restart`        RENAMED `write_output_state_to_file`
- `restart_in_file_name`  RENAMED `input_state_files`
- `restart_out_file_name` RENAMED `output_state_files`
- `single_file_in`        ADDED for low order models
- `single_file_out`       ADDED for low order models
- `has_cycling`           ADDED for low order models


------------------------------------------------------------------------------
## Jan 13th 2017 :: rma_fixed_filenames merge changes.  Revision: 10902

Specific namelist changes include:

1. Earlier versions of the RMA branch code supported both direct NetCDF
     reads/writes and the original binary/ascii DART format restart files.  
     As of the next update DART format files are no longer supported.  All
     I/O is NetCDF only.  If your model does not use NetCDF you will still
     need a model_to_dart and dart_to_model converter; otherwise all DART
     programs read the model's NetCDF files directly.  The namelist options
     related to selecting direct netcdf I/O have been removed.

1. Diagnostic and state space data (such as inflation, mean and sd
     information) that were previously stored in {Prior,Posterior}_Diag.nc
     are now broken up into multiple files and have fixed filenames. This
     decreases the IO time for diagnostic output and reduces the number of
     namelist options.

1. There is no longer support for observation space inflation
     (i.e. inf_flavor = 1).  Contact us at dart@ucar.edu if you have an
     interest in using this option.

#### Changes to the filter_nml are :

- `restart_in_file_name` has been replaced with `input_restart_file_list`.
   The namelist must contain one or more file names,
   each of which is a textfile containing a list of N
   NetCDF restart files, one per line for each ensemble member.
   For models with multiple domains (e.g. nested WRF or CLM)
   you must specify a listfile for each domain.

- `restart_out_file_name` has been replaced with `output_restart_file_list`.
   Same format as `input_restart_file_list`.

- `inf_in_file_name` REMOVED, now have fixed names of the form
   input_{prior,posterior}inf_{mean,sd}.nc

- `inf_out_file_name` REMOVED, now have fixed names of the form
   output_{prior,posterior}inf_{mean,sd}.nc.

- `inf_diag_filename` REMOVED

- `inf_output_restart` REMOVED, inflation restarts will be written
   out if inflation is turned on

- `output_inflation` REMOVED, inflation diagnostic files will be written
   if inflation is turned on

- `stages_to_write` There is more control over what state data
   to write.  Options are at stages : 'input', 'preassim', postassim', 'output'.  
   Stages preassim and postassim will output state data originally contained 
   within the copies of `Prior_Diag.nc` and `Posterior_Diag.nc`.
   See rma_doc/rma.html for details on the filename conventions. 
   For example, running filter with prior inflation enabled with stage 
   'preassim' enabled will produce files with names:
   - preassim_member_####.nc
   - preassim_{mean,sd}.nc
   - preassim_priorinf_{mean,sd}.nc

- `write_all_stages_at_end` important for large models - all output file 
   I/O is deferred until the end of filter, but will use more memory to store 
   the data.  More detailed info is in rma_doc/rma.html

- `output_restart_mean`        renamed output_mean

- `output_restart`             renamed output_restarts

- `direct_netcdf_{read,write}` REMOVED, always true

- `restart_list_file`          renamed input_restart_file_list

- `single_restart_file_in`     renamed single_file_in

- `single_restart_file_out`    renamed single_file_out

- `add_domain_extension`       REMOVED

- `use_restart_list`           REMOVED

- `overwrite_state_input`      REMOVED, equivalent functionality 
   can be set with `single_restart_file_in = single_restart_file_out`

#### Changes to the perfect_model_obs_nml are :

- `restart_in_filename` renamed `restart_in_file_names` takes a NetCDF 
   file. For multiple domains you can specify a list.

- `direct_netcdf_{read,write}` REMOVED, always true

#### Changes to the state_space_diag_nml are :

- `single_file` REMOVED, diagnostic files are now controlled 
   in `filter_nml` with `stages_to_write`

- `make_diagnostic_files` REMOVED, no longer produce 
   original `Prior_Diag.nc` and `Posterior_Diag.nc`

- `netCDF_large_file_support` REMOVED, always true

#### Changes to the state_vector_io_nml are :

- `write_binary_restart_files` REMOVED

#### Changes to the ensemble_manager_nml are :

- `flag_unneeded_transposes` -- REMOVED

#### Changes to the integrate_model_nml are :

- `advance_restart_format` -- REMOVED, only supporting NetCDF format.

#### Scripting with CESM

See `models/cam-fv/scripts_cesm1_5/assimilate.csh` for an example of how to
handle the new filename conventions.

```
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
```

#### ADDITIONAL NOTES :

1. currently the closest_member_tool is broken but plans on being fixed soon.
1. restart_file_tool and most model_to_dart/dart_to_model programs have been
  deprecated, since DART formated restarts are no longer supported.
1. some programs such as model_mod_check have not been fully tested and need
  to be exercised with the new naming conventions.

------------------------------------------------------------------------------
## ancient history 

To see previous history, it is probably best to use 

- `git log --follow` 
- `git diff --name-status XXXX YYYY` where XXXX and YYYY are commits, branches, ...

or something along those lines.



