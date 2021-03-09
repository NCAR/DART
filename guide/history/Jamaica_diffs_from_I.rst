DART "Jamaica release" difference document
==========================================

| This document describes the changes from the "Iceland" or "Post-Iceland" DART releases to the new **jamaica** release.
| **A word about MPI**: It is very important to remember that all the MPI routines are in one module --
  ``mpi_utilities_mod.f90``. If you are NOT using MPI - there is a simple equivalent module named
  ``null_mpi_utilities_mod.f90`` which is the default.

These files have changed (or are added)
---------------------------------------

| 1. ``adaptive_inflate_mod.f90``: This module was almost completely rewritten although fundamental operations remain
  unchanged. The only change to have direct impact on users is that there is no longer a namelist for this module
  (inflation is now controlled through the filter_nml) and inflation restarts are always in a single file, no matter
  what the choice of multiple restart files for state restarts is.
| There have been a number of structural and algorithmic changes. First, the module is now somewhat more object-oriented
  with an **adaptive_inflate_type** introduced that contains all the information about a particular type of inflation.
  The **adaptive_inflate_type** includes: the inflation_flavor (observation, spatially fixed state space, or
  spatially-varying state space) the unit open for inflation diagnostic output, logicals that control whether an
  inflation restart is read or written, and a logical that controls whether the inflation is done with a deterministic
  or stochastic algorithm. The file names for the input and output restarts and the diagnostic file are also included.
  The mean value and standard deviation of the inflation for the spatially fixed state space are in the type. Also a
  lower bound on the standard deviation and lower and upper bounds on the mean inflation value are included. All calls
  to apply or adjust inflation values now require an inflate_handle of type adaptive_inflate_type.
| Several algorithmic improvements were made to attempt to increase the numerical stability of the inflation algorithms.
  First, an improved quadratic solver was added. The rest of the inflation code was restructured to reduce the loss of
  precision in computing the mean of an updated inflation. A minor code change removed the possibility of inflation
  values of exactly 1.0 introducing inflation through round-off on some compilers. Numerical stability for extreme
  outliers in ensembles was greatly improved.
| There is evidence that inflation computations can continue to be sensitive to numerical round-off when run with (r4)
  precision instead of (r8). There may be a need to switch to enforced higher precision for challenging applications
  that require the use of (r4) as the default precision (like some WRF applications).

| 

| 2. ``assim_model_mod.f90``: There were several changes to the public interfaces. First, get_close_states and
  get_num_close_states no longer exist. This is part of moving the default ability to find close points out of model_mod
  and into the location mod. Second, an additional optional argument, override_write_format, was added to function
  open_restart_write. This allows a caller to force a write in either binary or ascii independent of the namelist
  setting.
| A number of non-algorithmic changes were made internal to the code. These include using the utilities_mod routines to
  check arguments for netcdf calls and addition of output control for the mpi implementation. More precision was added
  in diagnostic output of days and seconds to avoid a problem with WRF.

| 

| 3. ``assim_tools_mod.f90``: Large portions of this module were entirely rewritten. The parallel assimilation algorithm
  is implemented in the subroutine filter_assim which is now the ONLY public interface into assim_tools. The new
  parallel assimilation algorithm and its associated inflation computations are documented in the manuscript by Anderson
  and Collins that is to appear in the Journal of Atmospheric and Oceanic Technology and in section 22 of the tutorial.
  The manual definition of 'regions' no longer exists. assim_region.csh no longer exists.
| The namelist has also seen major changes. Namelist entries do_parallel, num_domains and parallel_command no longer
  exist. Namelist entry num_close_threshold has been renamed as adaptive_localization_threshold and has a default value
  of -1 which means no adaptive localization is in use. Namelist entry print_every_nth_obs is new and requests that a
  message be printed after every nth observation is assimilated in order to monitor progress of large test jobs.
| Algorithmically, the ``obs_increment_????`` routines are unchanged but the ensemble mean and variance are passed in as
  arguments now to avoid a redundant computation. The update_from_obs_increment routine also has the mean and variance
  for the observed variable input to avoid extra computation while the correl argument has been made optional and the
  correlation is only computed if it is requested in the call.
| An additional observation space increment option has been added along with the supporting subroutine
  obs_increment_hybrid subroutine. This interface is not documented at present and should not be used without consulting
  the DART development team. Additional support for a class of particle filters is provided by subroutine
  update_ens_from_weights. This routine is also undocumented and should not be used without consulting the DART
  development team.

| 

4. ``types_mod.f90``: Types i4, r4, c4 and c8 were added to the public list and additional comment documentation added.

| 

5. ``cov_cutoff_mod.f90``: The interface to comp_cov_factor was modified to add four new arguments, obs_loc, obs_kind,
target_loc, and target_kind. These arguments are not used in the default implementation but are made available for
applications in which users want to implement more complicated cutoff functions that can be a function of the location
and type of both the observation and the state variable.

| 

6. ``diagnostics/oned/obs_diag.f90``: The default value of rat_cri in the namelist was changed from 3.0 to 4.0 to avoid
throwing out too many observations as outside the expected range. Observation-space diagnostics have been extended to
handle the fact that inflation can happen in the prior or posterior (by keying on the new dart_qc variable). The new
dart_qc values were used to determine which observations were assimilated for posterior and prior estimates and a
histogram of the qc values is output. The namelist variable qc_threshold was renamed input_qc_threshold to further
distinguish the fact there are now two types of QC values: dart_qc and whatever comes with the observation to begin
with.

| 

7. ``diagnostics/threed_sphere/obs_diag``: Observation-space diagnostics have been extended to handle which observations
were assimilated for posterior and/or prior estimates (by keying on the new dart_qc variable). The namelist variable
qc_threshold was renamed input_qc_threshold to further distinguish the fact there are now two types of QC values:
dart_qc and whatever comes with the observation to begin with.

| 

| 8. ``ensemble_manager_mod.f90``: This module has been almost entirely rewritten for the mpi implementation. Its
  continues to manage ensemble storage but adds a general transpose capability that is central the the parallel
  implementation.
| The namelist retains the single_restart_file_in and single_restart_file_out but drops the in_core entry and adds
  perturbation_amplitude which controls the standard deviation of perturbations added to the state when starting from a
  single state estimate and generating a default ensemble. See the html documentation for details on the new interfaces.

| 

| 9. ``filter.f90``: There are major changes throughout filter for the mpi implementation. Foremost is the fact that the
  async=3 option no longer exists, which removes the need for filter_server.csh; nor do manually-defined 'regions'
  exist, which removes the need for assim_region.csh.
| Changes to namelist: for full details see the html documentation for filter and the tutorial documentation for
  adaptive inflation. The namelist can be divided into two parts, a new part that controls inflation (previously done in
  the now defunct adaptive_inflate.nml) and an old part that controls other aspects of the assimilation. In the old
  part, the namelist entries output_state_ens_mean, output_state_ens_spread, output_obs_ens_mean and
  output_obs_ens_spread have been removed. Mean and spread are now always output. New entries first_obs_days,
  first_obs_seconds, last_obs_days and last_obs_seconds have been added. These specify the time of the first and last
  observation in the obs_sequence file that are to be used. The default values of -1 indicate that these are to be
  ignored (see html documentation). The entry input_qc_threshold has been added: observations with an associated qc
  field in the input obs_sequence that is larger than the threshold are not assimilated. The entry
  output_forward_op_errors (default false) has been added: it outputs a detailed file containing a list of all failed
  forward observation operators in the assimilation when true. The entry output_timestamps has been added: when true it
  generates diagnostic ouput bounding each model advance call in filter.
| The inflation portion of the namelist is divided into two columns, the first controlling prior inflation and the
  second controlling posterior inflation. Details of these controls can be found in the html documentation. They are
  related to the defunct namelist values that were found in adaptive_inflate_mod.nml in the post-I release which
  supported only prior inflation.
| A number of other internal details were changed. The main program is now a single statement that calls a subroutine,
  filter_main. This allows for better memory management and avoids lots of shared 'module' storage from the main
  program. In order to support a fully modular smoother capability, the subroutine filter_state_space_diagnostics has
  been moved to smoother_mod.f90 and smoother_mod is used by filter. The arguments for all previously existing internal
  subroutines have been modified and are described in the html documentation.

| 

10. ``integrate_model.f90``: The namelist integrate_model.nml no longer exists. The old namelist had entries for the
input and output files, but they are now hard-coded to 'temp_ic' and 'temp_ud'. The target time for the advance is now
taken directly from the input file and the old namelist entries for target_time_days and target_time_seconds are no
longer needed. It is essential that integrate_model be compiled WITHOUT the parallel mpi_utilities.

| 

| 11. ``oned/location_mod.f90``: location modules are now responsible for providing interfaces to (efficiently) find a
  set of locations that are close to a given location. The new and modified public interfaces are: get_close_obs,
  get_close_obs_destroy, get_close_maxdist_init, and get_close_obs_init. In addition, a new type, get_close_type, is
  defined to store information that helps to do this type of search efficiently given a fixed set of locations to check.
  The oned location_mod doesn't have a need to do this efficiently, so these new interfaces are generally just stubs.
  The old get_close_obs still works as before to find locations that are within a certain distance of a given location.
| The get_dist interface has been changed to include two new arguments, kind1 and kind2. These are the kinds associated
  with the two locations. These arguments are not used in the default implementation of get_dist, but are made available
  for users who want to define distances using not only the location but also the kinds.

| 

| 12. ``threed_sphere/location_mod.f90``: The location module is now primarily responsible for the efficient search for
  close states / obs. A series of new interfaces have been added, along with namelist modifications, to support finding
  a subset of locations that are close to a single given location. This can be used both for get_close_obs and
  get_close_state computations and replaces the get_close_state interfaces that were in the model_mod.
| A new type, the get_close_type, is defined in a partially object-oriented fashion. For the threed_sphere, the
  algorithm works by partitioning the surface of the sphere using a longitude/latitude equally-spaced grid. This grid
  divides the sphere's surface into a set of nlon by nlat boxes. The first step in the efficient search computes the
  minimum distance between points in boxes that are separated by a given number of boxes in latitude and in longitude.
  This is accomplished by the new interface get_close_maxdist_init. This routine also accepts a cutoff radius and keeps
  a list of all box offsets from a box at a given latitude that are possibly within the radius. The second step takes a
  list of locations and places them into the appropriate boxes. This is performed by new interface get_close_obs_init.
  Finally, interface get_close_obs finds all observations that are less than the cutoff distance from a single input
  location and returns their indices in the original list of locations, along with the distance between them and the
  single base location if requested. An interface, get_close_obs_destroy, is also provided to destroy an instance of the
  get_close_type.
| Three new namelist entries control the number of boxes used in the search, nlon and nlat, and allow for detailed
  diagnostic output of the performance of the close search, output_box_info. The public interface print_get_close_type
  is also provided for debug and diagnostic use.
| The get_close_obs algorithm only partitions a subset of the sphere's surface into boxes if the input set of locations
  is confined to a small region of the surface. The algorithm works most efficiently when the average number of
  locations in a box is small compared to the total number of locations being searched but large compared to 1.
  Additional guidance in tuning the nlon and nlat control over the number of boxes is available from the DART
  development team.
| The get_dist interface has been changed to include two new arguments, kind1 and kind2. These are the kinds associated
  with the two locations. These arguments are not used in the default implementation of get_dist, but are made available
  for users who want to define distances using not only the location but also the kinds.

| 

| 13. ``model_mod.f90``: Unfortunately, there are minor changes to the model_mod public interfaces required to work with
  the switch to using the location_mod to find close locations. The public interface model_get_close_states is no longer
  required. Three new interfaces: get_close_maxdist_init, get_close_obs_init, and get_close_obs are required to use the
  new location module. In a minimal implementation, these three interfaces can just be satisfied by using the interfaces
  of the same name in the location module. The models/template/model_mod.f90 demonstrates this minimal implementation.
  Large models can implement their own modified algorithms for the get_close interfaces if desired for efficiency or
  correctness purposes. Examples of this can be seen in the model_mod.f90 for cam or wrf.
| An additional new interface, ens_mean_for_model has also been added. This routine is used to pass the ensemble mean
  state vector into model_mod by the filter before each assimilation step. This allows the model_mod to save this
  ensemble mean state if it is needed for computing forward operators as in some large atmospheric models (see cam). For
  low-order models, this interface can be a stub as shown in the template/model_mod.f90.

| 

14. ``PBL_1d/create_real_network.f90`` originated from create_fixed_network. It uses module_wrf to get obs from smos
file, with file, date, and interval controlled via the wrf1d namelist. Note that an obs_def is still required to control
which obs are actually written out. Normally, this would be created with create_obs_sequence. This would be run in place
of both create_fixed_network and perfect_model_obs.

| 

| 15. ``cam/model_mod.f90``: model_mod can now automatically handle the eulerian and finite volume CAMs (and should
  handle the Semi-Lagrangian core too), both single threaded and MPI-enabled. The latter enables efficient use of more
  processors than there are ensemble members. This situation is becoming more common, since DART can now assimilate
  using smaller ensembles, and massively parallel machines are becoming more common. This new mode of running jobs
  (async=4) replaces the async=3 and requires only 2 scripts instead of 4.
| The multi-core capability required reorganizing the state vector, so new filter_ic files will be necessary. These can
  be created, from the CAM initial files for the relevant dynamical core, using program trans_sv_pv_time0.
| The namelist has changed: state_names_pert has been replaced by pert_names, pert_sd and pert_base_vals as described in
  the cam/model_mod.html page. highest_obs_pressure_mb, which prevents obs above this height to be assimilated, has been
  joined by 2 other parameters. Observations on heights (or model levels) will have their vertical location converted to
  pressure, and be restricted by highest_obs_pressure_mb. highest_state_pressure_mb damps the influence of all obs on
  state variables above this height. max_obs_lat_degree restricts obs to latitudes less than this parameter, which is
  needed by some GPS observation sets.
| If more fields from the CAM initial files are to be added to the state vector, it may be necessary to add more 'TYPE's
  in model_mod, and more 'KIND's in DART/obs_kind/DEFAULT_obs_kind_mod.F90, and possibly a new obs_def_ZZZ_mod.f90.
| There is a new program, trans_pv_sv_pert0.f90, which can be useful in parameterization studies. It takes a model
  parameter, which has been added to the CAM initial files, and gives it a spread of values among the filter_ic files
  that it creates.

| 

16. ``wrf/model_mod.f90``: several researchers had their own subtly-different versions of WRF model_mod.f90. These
versions have been integrated (assimilated? ;) into one version. The new version performs vertical localization, support
for soil parameters, and a host of other features. Hui Liu (DAReS), Altug Aksoy, Yongsheng Chen, and David Dowell of
NCAR's MMM group are extensively using this model.

| 

17. ``DEFAULT_obs_def_mod.F90``: A new public interface, get_obs_def_key, was added to return the integer key given an
obs_def_type. The interface to read_obs_def has an additional intent(inout) argument, obs_val. This is NOT used in the
default implementation, but is required for the implementation of certain special observations like radar reflectivity.

17a. ``obs_def_radar_mod.f90``: added nyquist velocity metadata field to radial velocity.

17b. ``obs_def_QuikSCAT_mod.f90``: New module for the SeaWinds instrument (on the QuikSCAT satellite) data as available
in the NCEP BUFR files.

17c. ``obs_def_reanalysis_bufr_mod.f90``: Added land surface and ATOVS temperature and moisture values.

17d. ``obs_def_GWD_mod.f90``: module to define 'observations' of gravity wave drag that are needed to perform parameter
estimation studies.

| 

18. ``DEFAULT_obs_kind_mod.F90``: Added in several new raw variable types including KIND_CLOUD_LIQUID_WATER,
KIND_CLOUD_ICE, KIND_CONDENSATION_HEATING, KIND_VAPOR_MIXING_RATIO, KIND_ICE_NUMBER_CONCENTRATION,
KIND_GEOPOTENTIAL_HEIGHT, KIND_SOIL_MOISTURE, KIND_GRAV_WAVE_DRAG_EFFIC, and KIND_GRAV_WAVE_STRESS_FRACTION.

| 

19. ``obs_model_mod.f90``: There were major internal changes to this routine which implements the shell interface
advance of models that is used in the new async = 2 and async = 4 advance options. The public interfaces also changed
with the old get_close_states being removed and the advance_state interface being made public. The advance_state
interface appears the same except that the intent(in) argument model_size no longer exists. The advance_state interface
allows the model to be advanced to a particular target_time without going through the move_ahead interface that uses the
observation sequence to drive the advance.

| 

20. ``merge_obs_seq.f90``: This routine is now MUCH faster for both insertions and simple 'appends' and can now handle
multiple input files. Conversion between ASCII and binary formats (even for a single file) is now supported. Sorting by
time and the removal of unused blocks of observations is also possible.

| 

| 21. ``obs_sequence_mod.f90``: The obs_sequence_mod presents abstractions for both the obs_sequence and the obs_type.
| For the observation sequence, public interfaces delete_seq_head, delete_seq_tail, get_next_obs_from_key and
  get_prev_obs_from_key were added. Given a time and a sequence, the delete_seq_tail deletes all observations later than
  the time from the sequence and delete_seq_head deletes all observations earlier than the time. Given a sequence and an
  integer key, get_next_obs_from_key returns the next observation after the one with 'key' and get_prev_obs_from_key
  returns the previous observation. A bug in get_obs_time_range that could occur when the entire time range was after
  the end of the sequence was corrected. A bug that occurred when the only observation in a sequence was deleted with
  delete_obs_from_seq was corrected.
| For the obs_type section, public interfaces replace_obs_values and replace_qc were added. These replace the values of
  either the observations or the qc fields given a sequence and a key to an observation in that sequence. The interface
  read_obs had an optional argument, max_obs, added. It allows error checking to make sure that the maximum storage
  space in the sequence is not exceeded during a read. The value of the observation is now passed as an argument to
  read_obs_def where it can be used to make observed value dependent modifications to the definition. This is only used
  at present by the doppler velocity obs_def_mod when it does unfolding of aliased doppler velocities.

| 

22. ``perfect_model_obs.f90``: There were major internal changes to be consistent with the new ensemble_manager_mod and
to use a one-line main program that calls a subroutine to avoid lots of shared storage. The namelist has 4 additional
arguments, first_obs_days, first_obs_seconds, last_obs_days and last_obs_seconds. These specify times before which and
after which observations in the input obs_sequence should be ignored. The default value is -1 which implies that all
observations are to be used.

| 

23. ``random_nr_mod.f90``: Converted to use digits12 for real computations to avoid possible change in sequences when
reduced precision is used for the r8 kind defined in types_mod.f90.

| 

24. ``random_seq_mod.f90``: Interface init_random_seq was modified to accept an additional optional argument, seed,
which is the seed for the sequence if present.

| 

| 25. ``utilities_mod.f90``: Several modules had duplicate netCDF error checking routines; these have been consolidated
  into an nc_check() routine in the utilities module. A new set_output() routine can control which tasks in a multi-task
  MPI job write output messages (warnings and errors are written from any task). The default is for task 0 to write and
  all others not to. A routine do_output() returns .true. if this task should write messages. This is true by default in
  a single process job, so user code can always safely write: if (do_output()) write(**,**) 'informative message' In an
  MPI job only task 0 will return true and only one copy of the message will appear in the log file or on standard
  output.
| In an MPI job messages written via the error_handler() will prefix the message with the task number. The
  initialize_utilities() routine now takes an alternative log filename which overrides the default in the input.nml
  namelist; this allows utility programs to select their own separate log files and avoid conflicts with other DART
  programs. The MPI initialization and finalize routines call the utility init and finalize routines internally, so
  programs which use the MPI utilities no longer need to initialize the utilities separately.

| 

| 26. ``mpi_utilities_mod.f90``: A new module which isolates all calls to the MPI libraries to this one module. Includes
  interfaces for sending and receiving arrays of data, broadcasts, barriers for synchronization, data reduction (e.g.
  global sum), and routines for identifying the local task number and total number of tasks. Also contains a block and
  restart routine for use with the async=4 mode of interacting with a parallel MPI model advance. Programs using this
  module must generally be compiled with either an MPI wrapper script (usually called mpif90) or with the proper command
  line flags. Some MPI installations use an include file to define the MPI parameters, others use an F90 module. If the
  mpi_utilities_mod does not compile as distributed search the source code of this module for the string 'BUILD TIP' for
  more detailed suggestions on getting it to compile.
| When using MPI the call to initialize_mpi_utilities() must be made as close to the start of the execution of the
  program as possible, and the call to finalize_mpi_utilities() as close to the end of execution as possible. Some
  implementations of the MPICH library (which is common on Linux clusters) require that MPI be initialized before any
  I/O is done, and other implementations (SGI in particular) will not allow I/O after MPI is finalized. These routines
  call the normal utilities init and finalize routines internally, so at the user level only the mpi versions need to be
  called.

| 

27. ``null_mpi_utilities_mod.f90``: A module which has all the same entry points as the mpi_utilities_mod but does not
require the MPI library. A program which compiles with this module instead of the real MPI utilities module can only be
run with a single task since it cannot do real parallel communication, but does not require the MPI libraries to compile
or link. This is the default module -- you cannot simultaneously use both the mpi_utilities_mod and the
null_mpi_utilities_mod.

| 

28. ``mkmf/mkmf``: The mkmf program takes a new -w argument. If specified, the resulting makefile will call 'wrappers'
for the fortran compiler and loader. The default compiler and loader are $(FC) and $(LD); with the -w flag they will
become $(MPIFC) and $(MPILD). In the mkmf.template file you can then define both the MPI wrappers (generally 'mpif90')
and the regular F90 compiler.

| 

29. ``mkmf.template.*``: The mkmf.template files have been consolidated where possible and heavily commented to reflect
nuances of implementations on different comiler/OS combinations. The habit of appending individual platform names (which
led to file creep) is hopefully broken.

| 

30. ``input.nml``: All the default input.nml namelists are now easily 'diff'ed against the corresponding
input.nml.*_template files. This enhances the ability to determine what values are different than the default values.

| 

| 31. ``DART/shell_scripts/DiffCVS_SVN``: is a new script that identifies differences in two parallel instantiations of
  source files. I used it so many times during the migration from CVS to SVN that I decided to add it to the DART
  project. This script compares all the source files in two directories ... after removing the trivial differences like
  the form of the copyright tags (its a new year, you know) and the fact that SVN has different revision numbers than
  CVS. If you set an environment variable XDIFF to a graphical comparator like 'xdiff' or 'xxdiff' it is used to display
  the differences between the files. With no arguments, a usage note is printed.
| The hope is that you can use it see on your 'old' sandboxes in the directories where you have modified code and
  compare that to the new code to see if there are any conflicts. The directories do not have to be under CVS or SVN
  control, by the way.

| 

32. ``models/PBL_1d/src/*``: These source files are 'directly' from other developers and had file extensions that
required the use of special compilation flags. There is now a script PBL_1d/shell_scripts/ChangeExtensions.csh that not
only changes the file extensions (to something the compilers understand i.e. F90) it also modifies the ``path_names_*``
files appropriately. The original files had an extension of .F even though they used F90-like formatting. .F is
generally understood to mean the contents of the file abide by the F77 fixed-format syntax ... columns 2-5 are for line
numbers, column 7 is a line-continuation flag ... etc. Now if we can only get them to not rely on the 64bit real
autopromotion ...

| 

| 33. ``DART/matlab``: The matlab scripts have experienced no major overhauls, but do handle a few more models than
  previously. The next release is scheduled to have full matlab support for CAM, and WRF. The ReadASCIIObsSeq.m function
  has a couple (backwards compatible) tweaks to accomodate some changes in R2006a, R2006b.
| The biggest change is the addition of ``plot_observation_locations.m``, which facilitates exploring observation
  locations (by type) for any 3D (i.e. real-world) observation sequence. This ability is derived by running ``obs_diag``
  with a namelist variable set such that a matlab-readable dataset is written out.

| 

34. ``models/ikeda``: There is a whole new model - courtesy of Greg Lawson of CalTech. A nice 2-variable system that
does not require any fancy time-stepping routines. Thanks Greg!

| 

35. ``obs_def_dew_point_mod.f90``: implements a more robust method (based on Bolton's Approximation) for computing dew
point.

| 

These files are obsolete - and have been deleted
------------------------------------------------

1. ``assim_tools/assim_region.f90``: new algorithms and methodologies fundamentally replace the need to manually define
and then assimilate regions. Think of all those intermediate files that are not needed!

2. ``filter_server*`` are no longer needed because async == 3 is no longer supported (again - replaced by MPI).

3. scripts that advance 'ensembles' are all gone - again because of the MPI implementation. The only script now needed
is 'advance_model.csh'.

4. ``smoother/smoother.f90``: The standalone smoother program has become a module and the functionality is now part of
the filter program.
