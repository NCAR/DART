---
title: Documentation
layout: default
---

# DART Manhattan Release Notes

## Dart Overview

The Data Assimilation Research Testbed (DART) is designed to facilitate
the combination of assimilation algorithms, models, and real (or
synthetic) observations to allow increased understanding of all three.
The DART programs are highly portable, having been compiled with many
Fortran 90 compilers and run on linux compute-servers, linux clusters,
OSX laptops/desktops, SGI Altix clusters, supercomputers running AIX,
and more. Read the
[Customizations](Getting_Started.html#customizations)
section for help in building on new platforms.

DART employs a modular programming approach to apply an Ensemble Kalman
Filter which adjusts model values toward a state that is more consistent
with information from a set of observations. Models may be swapped in
and out, as can different algorithms in the Ensemble Kalman Filter. The
method requires running multiple instances of a model to generate an
ensemble of states. A forward operator appropriate for the type of
observation being assimilated is applied to each of the states to
generate the model's estimate of the observation. Comparing these
estimates and their uncertainty to the observation and its uncertainty
ultimately results in the adjustments to the model states. See the
[DART_LAB](dart_lab.html) demos or read more
[in the DART tutorial](Tutorial.html).

DART diagnostic output can be written that contains the model state
before and after the adjustment, along with the ensemble mean and
standard deviation, and prior or posterior inflation values if inflation
is enabled. There is also a text file, `obs_seq.final`, with the model
estimates of the observations. There is a suite of MATLAB® functions
that facilitate exploration of the results, but the netCDF files are
inherently portable and contain all the necessary metadata to interpret
the contents with other analysis programs such as NCL, R, etc.

To get started running with Lorenz 63 model refer to the 
[Lorenz 63 section](Getting_Started.html#Lorenz63) of the *Getting Started* document.

<span id="CurrentUsers"></span>

\[[top](#)\]

-----

## About the Manhattan layout change ...

If you have been updating from the rma\_trunk branch of the DART
subversion repository you will notice that the code tree has been
simplified to be more intuitive for users. The new top level directory
structure looks like :

  - *README*
  - *COPYRIGHT*
  - *assimilation\_code*
  - *build\_templates*
  - *diagnostics*
  - *docs*
  - *models*
  - *observations*

We suggest that current users checkout a fresh version of Manhattan in a
new location. To see which files need to be moved, run 'svn status' on
your original checked out version. Anything with an M or ? in the first
column needs to be moved to the new location in the new tree. Please
[contact DART](mailto:dart@ucar.edu) if you have any issues migrating
your existing code to the new tree structure.

There is a list of non-backwards compatible changes 
([see below](#Nonbackward)), and a list of new options and functions.

The Manhattan release will continue to be updated for the next few
months as we continue to add features. Checking out the Manhattan
release branch and running 'svn update' from time to time is the
recommended way to update your DART tree.

<span id="Nonbackward"></span>

\[[top](#)\]

-----

## Non-backwards Compatible Changes

Unlike previous releases of DART, this version contains more
non-backwards compatible changes than usual. Please examine the
following list carefully. We do suggest you check out the Manhattan
release into a new location and migrate any local changes from previous
versions as a second step.

Changes in the Manhattan release (15 May 2015) which are *not* backwards
compatible with the Lanai release (13 Dec 2013):

1.  We no longer require model data to be converted to DART format
    restart files. We directly read and write netCDF format only. To
    specify the input and output files for filter, there are new
    namelist items in the &filter_nml namelist:
    `input_state_file_list` and `output_state_file_list` .

2.  The information formerly in `Prior_Diag.nc` and
    `Posterior_Diag.nc` has been moved. If you are reading and writing
    ensemble members from different files, the state information, the
    ensemble mean and standard deviation, and the inflation mean and
    standard deviation will all be read and written to separate files:  

      - *\[stage\]_member\_\#\#\#\#.nc*
      - *\[stage\]_mean.nc*
      - *\[stage\]_sd.nc*
      - *\[stage\]_priorinf_{mean,sd}.nc* (if prior inflation is turned on)
      - *\[stage\]_postinf_{mean,sd}.nc* (if posterior inflation is turned on)

    If you are reading and writing ensemble members from a single file,
    all this information will now be in a single netCDF file but will be
    stored in different variables inside that file:  

      - *\[var\].nc*
      - *\[var\]_mean.nc*
      - *\[var\]_sd.nc*
      - *\[var\]_priorinf_{mean,sd}.nc* (if prior inflation is turned on)
      - *\[var\]_postinf_{mean,sd}.nc* (if posterior inflation is turned on)
      
    We also now have options for writing files at four stages of the
    assimilation cycle: `'input', 'preassim', 'postassim', 'output'`.
    This is set in the \&filter_nml namelist with ```stages_to_write```.

3.  New `model_mod.f90` required routines:  

      - *vert_convert()*
      - *query_vert_localization_coord()*
      - *pert_model_copies()*
      - *read_model_time()*
      - *write_model_time()*
    
    There are default version of these available to use if you have no
    special requirements.

4.  Several of the `model_mod.f90` argument lists have changed
    
      - *model_interpolate()* now takes in the 'state_handle' as an
        argument rather than a state vector array. It also return an
        array of 'expected_obs' and 'istatus' for each of the ensemble
        members
      - *get_state_meta_data()* also requires the 'state_handle' as
        an argument rather than a state vector array.
      - *nc_write_model_atts()* has an additional argument
        `model_mod_writes_state_variables`. If `TRUE` then the model_mod
        is expected to write out the state variables, if `FALSE` DART will
        write out the state variable (this is the prefered method for
        adding new models, it requires less code from the model developer)

5.  There are several namelist changes mainly in the &filter_nml and
    &perfect_model_mod which are outlined in detail in
    [Manhattan_diffs_from_Lanai](Manhattan_diffs_from_Lanai.html)

6.  All modules have been moved to *DART/assimilation_code/modules/*
    directory. And similarly all of the programs have moved to
    *DART/assimilation_code/programs/*

7.  The location modules which were stored in *locations* have moved to
    *DART/assimilation_code/location* directory

8.  The observation converters which were stored in *observations* have
    moved to *DART/observations/obs_converters* directory

9.  The forward operators have moved from
    *obs_def/obs_def_\*_mod.f90* to
    *observations/forward_operators*

10. The tutorial files have moved to *DART/docs/tutorial* directory

11. The *fill_inflation_restart* program now creates
    inflation files in netCDF format. You can still create inflation files
    using *ncap2*. Here is an example using version 4.4.2 or
    later of the NCO tools:
    
    ``` 
      ncap2 -s "T=1.0;U=1.0;V=1.0" wrfinput_d01 prior_inf.nc'
      ncap2 -s "T=0.6;U=0.6;V=0.6" wrfinput_d01 prior_sd.nc'
    ```

12. The default flags in the mkmf_template.XXX files have been updated
    to be more consistent with current compiler versions.

13. If you enable the sampling error correction option, the required
    data is now read from a single netCDF file which supports multiple
    ensemble sizes. A program is provided to compute additional ensemble
    sizes if they are not in the default file.

14. Our use of TYPES and KINDS has been very confusing in the past. In
    Manhattan we have tried to make it clearer which things in DART are
    generic quantities (QTY) - temperature, pressure, etc - and which
    things are specific types of observations - Radiosonde_temperature,
    Argo_salinity etc.
    
    Below is a mapping between old and new routine names here for
    reference. We have made these changes to all files distributed with
    DART. If you have lots of code developed outside of the subversion
    repository, please contact [DART](mailto:dart@ucar.edu) for a sed
    script to help automate the changes.
    
    Public routines, existing name on left, replacement on right:
    
    ~~~ 
        assimilate_this_obs_kind()     =>     assimilate_this_type_of_obs(type_index)
        evaluate_this_obs_kind()       =>       evaluate_this_type_of_obs(type_index)
        use_ext_prior_this_obs_kind()  =>  use_ext_prior_this_type_of_obs(type_index)
        
        get_num_obs_kinds()            =>  get_num_types_of_obs()
        get_num_raw_obs_kinds()        =>  get_num_quantities()
        
        get_obs_kind_index()           => get_index_for_type_of_obs(type_name)
        get_obs_kind_name()            => get_name_for_type_of_obs(type_index)
        
        get_raw_obs_kind_index()       =>  get_index_for_quantity(qty_name)
        get_raw_obs_kind_name()        =>  get_name_for_quantity(qty_index)
        
        get_obs_kind_var_type()        =>  get_quantity_for_type_of_obs(type_index)
        
        get_obs_kind()                 =>  get_obs_def_type_of_obs(obs_def)
        set_obs_def_kind()             =>  set_obs_def_type_of_obs(obs_def)
        
        get_kind_from_menu()           =>  get_type_of_obs_from_menu()
        
        read_obs_kind()                =>   read_type_of_obs_table(file_unit, file_format)
        write_obs_kind()               =>  write_type_of_obs_table(file_unit, file_format)
        
        maps obs_seq nums to specific type nums, only used in read_obs_seq:
        map_def_index()                => map_type_of_obs_table()
        
        removed this.  apparently unused, and simply calls get_obs_kind_name():
        get_obs_name()
        
        apparently unused anywhere, removed:
        add_wind_names()
        do_obs_form_pair()
    ~~~
    
    Public integer parameter constants and routine formal argument
    names, old on left, new on right:
    
    ~~~ 
       KIND_ => QTY_
       kind  => quantity
       
       TYPE_ => TYPE_
       type  => type_of_obs
       
       integer parameters:
       max_obs_generic  =>  max_defined_quantities  (not currently public, stays private)
       max_obs_kinds    =>  max_defined_types_of_obs 
    ~~~

15. For smaller models we support single file input and output. These
    files contain all of the member information, mean, standard
    deviation and inflation values for all of the state variables. This
    can be run with cycling and all time steps will be appended to the
    file.
    
    For *perfect_model_obs* we provide a `perfect_input.cdl` file
    which contains a single ensemble member which will be considered the
    'truth' and observations will be generated based on those values.
    The output will contain all of the cycling timesteps all of the
    state variables.
    
    For *filter* we provide a `filter_input.cdl` file which contains
    all of the state member variables and potentially inflation mean and
    standard deviation values. The output will contain all of the
    cycling timesteps all of the state variables. Additionally you have
    the option to write out different stages during the assimilation in
    the \&filter_nml 'stages_to_write' mentioned above.
    
    To generate a netCDF file from a .cdl file run:
   
    ``` 
    ncgen -o perfect_input.nc perfect_input.cdl
    ncgen -o filter_input.nc filter_input.cdl
    ```

<span id="NewFeatures"></span>

\[[top](#)\]

-----

## New Features

  - DART now reads and writes netCDF files for the model state
    information. If your model uses netCDF file format, you no longer
    need *model_to_dart* or *dart_to_model* to translate to a DART
    format file. If your model does not use netCDF, you can adapt your
    *model_to_dart* and *dart_to_model* executables to read and write a
    netCDF file for DART to use.  
      
    The read/write code is part of the core DART routines so no code is
    needed in the model_mod model-specific module. There is a new
    routine [add_domain()](state_structure.html) that a
    model_mod::static_init_model() can user to define which netCDF
    variables should be part of the model state, and what DART quantity
    (formerly kind) they correspond to.
  - DART no longer limits the size of a model state to the size of a
    single MPI task's memory. The state is read in variable by variable
    and distributed across all MPI tasks, so the memory use is much
    smaller than previous versions of DART. One-sided MPI communication
    is used during the computation of forward operator values to get
    required parts of the state from other tasks.
  - Many of the DART namelists have been simplified, and some items have
    moved to a more specific namelist.
  - Observation sequence files can include externally computed forward
    operator values which can be used in the assimilation instead of
    calling a forward operator inside DART.
  - The DART directory structure has been reorganized to make it easier
    to identify the various software tools, modules, documentation and
    tutorials supplied with the system.
  - The MATLAB® diagnostic routines have been updated to not require the
    MEXNC toolbox. These routines use the built-in netCDF support that
    comes with MATLAB®.
  - There is a new Particle Filter type. Please contact us if you are
    interested in using it.
  - DART can now take subsets of observation types and restrict them
    from impacting certain quantities in the state during the
    assimilation. A tool to simplify constructing the table of
    interactions is provided (obs_impact_tool).
  - State Structure
      - Contains information about dimensions and size of variables in
        your state. There is a number of accessor functions to get
        variable information such as `get_variable_size()`. See the
        [state_structure.html](state_structure.html) for more details.
  - The POP model_mod now can interpolate Sea Surface Anomaly
    observations.

<span id="SupportedModels"></span>

\[[top](#)\]

-----

## Supported Models

Currently we support the models listed below. There are several new
models that have been added that are not on the Lanai Release including
CM1, CICE, and ROMS. Any previously supported models not on this list
are still supported in DART
[classic](http://www.image.ucar.edu/DAReS/DART/classic/index.html)

  - **9var**
      - DART interface documentation for the
        [9var](models/9var/model_mod.html) model.
  - **bgrid_solo**
      - DART interface documentation for the
       [bgrid solo](models/bgrid_solo/model_mod.html) model.
  - **CAM-FV, CAM-CHEM, WACCM, WACCM-X**
      - DART interface documentation for the
       [CAM finite volume](models/cam-fv/model_mod.html) family of global atmospheric models.
      - Documentation for the
       [CAM model](http://www.cesm.ucar.edu/models/atm-cam/).
  - **CICE (NEW)**
      - DART interface documentation for the
        [CICE](models/cice/model_mod.html) model.
      - Documentation for the
        [CICE model](http://www.cesm.ucar.edu/models/ccsm4.0/cice/).
  - **CM1 (NEW)**
      - DART interface documentation for the
        [CM1 cloud-resolving model](models/cm1/model_mod.html).
      - Documentation for the
        [CM1 model](http://www2.mmm.ucar.edu/people/bryan/cm1/).
  - **forced_lorenz_96**
      - DART interface documentation for the
        [forced lorenz_96](models/forced_lorenz_96/model_mod.html) model.
  - **lorenz_63**
      - DART interface documentation for the
        [lorenz_96](models/lorenz_63/model_mod.html) model.
  - **lorenz_84**
      - DART interface documentation for the
        [lorenz_84](models/lorenz_84/model_mod.html) model.
  - **lorenz_96**
      - DART interface documentation for the
        [lorenz_96](models/lorenz_96/model_mod.html) model.
  - **lorenz_04**
      - DART interface documentation for the
        [lorenz_04](models/lorenz_04/model_mod.html) model.
  - **MPAS_ATM** (netCDF overwrite not supported for update_u_from_reconstruct = .true. )
      - DART interface documentation for the
        [MPAS atmosphere](models/mpas_atm/model_mod.html) model.
      - Documentation for the
        [MPAS model](https://mpas-dev.github.io/atmosphere/atmosphere.html).
  - **POP**
      - DART interface documentation for the
        [POP](models/POP/model_mod.html) global ocean model.
      - Documentation for the
        [POP model](http://www.cesm.ucar.edu/models/ccsm2.0/pop/).
  - **ROMS (NEW)**
      - DART interface documentation for the
        [ROMS](models/ROMS/model_mod.html) regional ocean model.
      - Documentation for the [ROMS model](https://www.myroms.org/).
  - **simple_advection**
      - DART interface documentation for the
        [simple advection](models/simple_advection/model_mod.html) model.
  - **WRF**
      - DART interface documentation for the
        [WRF](models/wrf/model_mod.html) regional forecast model.
      - Documentation for the [WRF model](http://www.wrf-model.org/index.php).
  - **NOAH-MP**
      - Documentation for [(NOAH-MP LSM)](https://ral.ucar.edu/solutions/products/noah-multiparameterization-land-surface-model-noah-mp-lsm)
  - **FESOM**
      - DART interface documentation for the
        [FESOM model interface](models/FESOM/Readme.md)
  - **GITM**
      - Documentation for the The Global Ionosphere Thermosphere Model
        [(GITM)](https://ccmc.gsfc.nasa.gov/models/modelinfo.php?model=GITM)


The *DART/models/template* directory contains sample files for adding a
new model. See the 
[Adding a Model](Models.html#adding_a_model)
section of the DART web pages for more help on adding a new model.

<span id="ChangedModels"></span>

\[[top](#)\]

-----

## Changed Models

  - WRF
      - Allow advanced microphysics schemes (needed interpolation for 7
        new kinds)
      - Interpolation in the vertical is now done in log(p) instead of
        linear pressure space. log(p) is the default, but a compile-time
        variable can restore the linear interpolation.
      - Added support in the namelist to avoid writing updated fields
        back into the wrf netCDF files. The fields are still updated
        during the assimilation but the updated data is not written back
        to the wrfinput file during the dart_to_wrf step.
      - Fixed an obscure bug in the vertical convert routine of the wrf
        model_mod that would occasionally fail to convert an obs. This
        would make tiny differences in the output as the number of mpi
        tasks change. No quantitative differences in the results but
        they were not bitwise compatible before and they are again now.
  - CAM
      - DART/CAM now runs under the CESM framework, so all options
        available with the framework can be used.
      - Support for the SE core (HOMME) has been developed but is NOT
        part of this release. Please contact the [DART Development
        Group](mailto:dart@ucar.edu) if you have an interest in this
        configuration of CAM.
  - Simple Advection Model
      - Fixed a bug where the random number generator was being used
        before being called with an initial seed.

<span id="NewFOs"></span>

\[[top](#)\]

-----

## New Observation Types/Forward Operators

  - Many new observation types related to land and atmospheric chemistry
    have been added. See the
    [obs_kind_mod.f90](../../assimilation_code/modules/observations/DEFAULT_obs_kind_mod.F90)
    for a list of the generic quantities now available.
  - New forward operator for Sea Ice (cice) ice thickness observations.
    See the
    [obs_def_cice_mod.f90](../../observations/forward_operators/obs_def_cice_mod.f90)
    file for details.
  - New forward operator for Carbon Monoxide (CO) Nadir observations.
    See the
    [obs_def_CO_Nadir_mod.f90](../../observations/forward_operators/obs_def_CO_Nadir_mod.f90)
    file for details.
  - New forward operator for Total Cloud Water in a column observations.
    See the
    [obs_def_cwp_mod.f90](../../observations/forward_operators/obs_def_cwp_mod.f90)
    file for details.

<span id="NewObs"></span>

\[[top](#)\]

-----

## New Observation Types/Sources

  - AVISO  
    Added an observation converter for Sea Surface Height Anomaly
    observations. Documentation in
    [convert_aviso.f90](../../observations/obs_converters/AVISO/convert_aviso.f90)
    (source).
  - cice  
    Added an obs_sequence converter for Sea Ice observations.
    Documentation in
    [cice_to_obs.html](../../observations/obs_converters/cice/cice_to_obs.html).
  - GPSPW  
    Added an obs_sequence converter for GPS precipitable water
    observations. Documentation in
    [convert_gpspw.f90](../../observations/obs_converters/GPSPW/convert_gpspw.f90)
    (source).
  - MODIS  
    Added an obs_sequence converter for MODIS FPAR (Fraction of
    Photosynthetically Active Radiation) and LAI (Leaf Area Index)
    obseverations. Documentation in
    [MOD15A2_to_obs.html](../../observations/obs_converters/MODIS/MOD15A2_to_obs.html).
  - ok_mesonet  
    Added an obs_sequence converter for the Oklahoma Mesonet
    observations. Documentation in
    [ok_mesonet.html](../../observations/obs_converters/ok_mesonet/ok_mesonet.html).
  - ROMS  
    Added an obs_sequence converter for ROMS ocean data. This converter
    includes externally computed forward operators output from the ROMS
    model using FGAT (First Guess At Time) during the model run.
    Documentation in
    [convert_roms_obs.f90](../../observations/obs_converters/ROMS/convert_roms_obs.f90)
    (source).
  - SSUSI  
    Added an obs_sequence converter for wind profiler observations.
    Documentation in
    [convert_f16_edr_dsk.html](../../observations/obs_converters/SSUSI/convert_f16_edr_dsk.html).
  - tropical_cyclone  
    Added an obs_sequence converter for ASCII format tropical cyclone
    track observations. Documentation in
    [tc_to_obs.html](../../observations/obs_converters/tropical_cyclone/tc_to_obs.html).

<span id="NewDiagnostics"></span>

\[[top](#)\]

-----

## New Diagnostics and Documentation

**Better Web Pages.** We've put a lot of effort into expanding our
documentation. For example, please check out
[the MATLAB diagnostics section](Diagnostics.html#mat_obs)
or the pages outlining the
[observation sequence file contents](Observations.html#obs_seq_overview).  

  - The MATLAB® diagnostic routines have been updated to remove the
    dependency on third-party toolboxes. These routines use the built-in
    netCDF support that comes with basic MATLAB® (no other toolboxes needed).

  - The MATLAB® routines have been updated to work consistently across
    many MATLAB® versions (which is no small trick considering how many non-backward-compatible
    changes they have made!).
  
But there's always more to add.
**Please use the GitHub Wiki or Issues facilities to let us know where we are lacking or what doesn't work.**

<span id="NewUtilities"></span>

\[[top](#)\]

-----

## New Utilities

This section describes updates and changes to the tutorial materials,
scripting, setup, and build information since the Lanai release.

  - *obs_impact_tool* please refer to
    [Website](https://www.image.ucar.edu/DAReS/DART/Manhattan/assimilation_code/programs/obs_impact_tool/obs_impact_tool.html)
    or 
    [local file](../../assimilation_code/programs/obs_impact_tool/obs_impact_tool.html)
  - *gen_sampling_error_table* now computes sampling error correction
    tables for any ensemble size.
  - *compute_error*
    [Website](https://www.image.ucar.edu/DAReS/DART/Manhattan/assimilation_code/programs/compute_error/compute_error.html)
    or
    [local file](../../assimilation_code/programs/compute_error/compute_error.html)

<span id="KnownProblems"></span>

\[[top](#)\]

-----

## Known Problems

  - The move to provide DART through GitHub offers us the opportunity
    to revise the way we distribute documentation. We are still exploring
    options for the best way to provide content that is accessible both
    offline and through GitHub Pages. Please bear with us, but please do
    point out broken links and offer suggestions on ways we can improve.

  - There are many changes in this release and more updates are expected
    to come soon. We are not aware of any obvious bugs, but if you
    encounter any unexpected behavior please contact us. Please watch
    the dart-users email list for announcements of updates to the
    release code, and be prepared to issue a GitHub Pull Request from
    the DART repository to your fork, or a git 'fetch' from the DART repository
    to get updated files.

\[[top](#)\]

-----
