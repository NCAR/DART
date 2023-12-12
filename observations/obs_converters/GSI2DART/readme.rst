GSI2DART
========

Overview
--------

The GSI2DART converter was contributed by **Craig Schwartz** and **Jamie Bresch** of the 
Mesoscale & Microscale Meteorology Lab at NSF NCAR. *Thanks Craig and Jamie!* 

This converter is designed to convert observation files created by the Gridpoint 
Statistical Interpolation (GSI) system maintained by the National Oceanic and 
Atmospheric Administration (NOAA) into DART observation sequence files.
The files created by GSI are 'BIG_ENDIAN' and have filenames such as:

- diag_amsua_metop-a_ges.ensmean
- diag_amsua_metop-a_ges.mem001
- diag_amsua_metop-a_ges.mem002
- diag_amsua_n18_ges.ensmean
- diag_amsua_n18_ges.mem001
- diag_amsua_n18_ges.mem002
- diag_amsua_n19_ges.ensmean
- diag_amsua_n19_ges.mem001
- diag_amsua_n19_ges.mem002
- diag_conv_ges.ensmean
- diag_conv_ges.mem001
- diag_conv_ges.mem002

The DART converter uses routines from the GSI system that use the Message Passing 
Interface (MPI) to process observations in parallel (even when converting a small 
amount of observations) so MPI is required to execute this observation converter.

Due to these prerequisites, we provide a detailed description of this directory to 
guide the user.

This directory contains copies of several source code files from GSI. 
The GSI source code is available via a Github repository managed by NOAA's 
Environmental Modeling Center (EMC):

https://github.com/NOAA-EMC/GSI

To differentiate between the sets of code, we refer to the root directory of the 
NOAA-EMC repository as ``GSI`` and refer to the root directory of this observation 
converter as ``GSI2DART``.

``GSI2DART/enkf`` copies seven files from ``GSI/src`` mostly without modification:

1. ``GSI2DART/enkf/constants.f90`` from ``GSI/src/gsi/constants.f90``
2. ``GSI2DART/enkf/kinds.F90`` from ``GSI/src/gsi/kinds.F90``
3. ``GSI2DART/enkf/mpi_readobs.f90`` from ``GSI/src/enkf/mpi_readobs.f90``
4. ``GSI2DART/enkf/readconvobs.f90`` from ``GSI/src/enkf/readconvobs.f90``
5. ``GSI2DART/enkf/read_diag.f90`` from ``GSI/src/gsi/read_diag.f90``
6. ``GSI2DART/enkf/readozobs.f90`` from ``GSI/enkf/readozobs.f90``
7. ``GSI2DART/enkf/readsatobs.f90`` from ``GSI/enkf/readsatobs.f90``

Note that within ``GSI`` the source file ``kinds.F90`` has an upper-case ``F90`` 
suffix. Within the ``GSI2DART`` observation converter, it gets preprocessed 
into ``mykinds.f90`` with a lower-case ``f90`` suffix. Case-insensitive filesystems 
should be banned ... until then, it is more robust to implement some name change 
during preprocessing. 

The following three files had their open() statements modified to read 
'BIG_ENDIAN' files without the need to compile EVERYTHING with 
the ``-convert big_endian`` compiler option. Using the DART open_file() 
routine also provides some nice error handling.

- original: ``open(iunit,form="unformatted",file=obsfile,iostat=ios)``
- modified: ``iunit = open_file(obsfile,form='unformatted',action='read',convert='BIG_ENDIAN')``

1. ``GSI2DART/enkf/readconvobs.f90``
2. ``GSI2DART/enkf/readozobs.f90``
3. ``GSI2DART/enkf/readsatobs.f90``

DART Modifications
------------------

Within GSI2DART
~~~~~~~~~~~~~~~

The source files within ``GSI2DART`` are:

1. ``gsi_to_dart.f90``: the main program.
2. ``dart_obs_seq_mod.f90``: the DART obs_seq output subroutine.
3. ``params.f90``: the same module name as ``GSI/src/enkf/params.f90`` but with different content. This version is used to avoid modifying ``GSI2DART/enkf/read*.f90``.
4. ``radinfo.f90``: the same module name as ``GSI/src/gsi/radinfo.f90`` but with different content. This version is used to avoid modifying ``GSI2DART/enkf/read*.f90``.
5. ``mpisetup.f90``: the same module name as ``GSI/src/enkf/mpisetup.f90`` but with different content. This version is used to avoid dependency on ``GSI``.

Elsewhere in the repository
~~~~~~~~~~~~~~~~~~~~~~~~~~~

This observation converter required modifying two files and adding a module for 
radiance observation types.

- Modified ``../../forward_operators/DEFAULT_obs_def_mod.F90``
- Modified ``../../DEFAULT_obs_kind_mod.F90``
- Added ``../../forward_operators/obs_def_radiance_mod.f90`` which has radiance observation types


Additional files and directories
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

1. ``satinfo`` is a file read by ``radinfo.f90`` and must exist in the ``GSI2DART/work`` directory.
2. ``datapath`` specifies the directory containing the data to be converted -- it is specified in the ``gsi_to_dart_nml`` namelist in ``GSI2DART/work/input.nml``.
3. ``submit.csh`` is contained in ``GSI2DART/work/`` -- it runs the gsi_to_dart converter once it has been compiled. Again, since GSI requires MPI, multiple processors must be requested to run the gsi_to_dart executable.

Issues
------

1. The converter requires an ensemble size greater than one and will MPI_Abort() 
if only one ensemble member is requested.

The following are issues previously recorded in the README:

1. Radiance and surface pressure bias correction
2. Surface pressure altimeter adjustment?
3. Specific humidity obs are transformed to relative humidity.  What to do? [Just run EnSRF with psuedo_rh=.false. and assimilate RH obs]
4. DART must use W and PH as control variables [okay, EnSRF can do this too (nvars=6 for WRF-ARW)]
5. Does DART not do vertical localization for surface obs?

.. code-block:: fortran

   ! If which_vert has no vertical definition for either location do only horizontal
   if(loc1%which_vert == VERTISUNDEF .or. loc2%which_vert == VERTISUNDEF) comp_h_only = .true.
   ! If both verts are surface, do only horizontal
   if(loc1%which_vert == VERTISSURFACE .and. loc2%which_vert == VERTISSURFACE) comp_h_only = .true.

Running with 32 bit reals
~~~~~~~~~~~~~~~~~~~~~~~~~

The converter has been tested with 64-bit reals as well as 32-bit reals 
(i.e. r8=r4 and -D_REAL_4). The answers are different only at the roundoff level.

This requires changes in two places:

1. ``DART/assimilation_code/modules/utilities/types_mod.f90`` change required:  r8 = r4
2. ``GSI2DART/work/quickbuild.sh`` change required: ``-D_REAL4_``

If these are not set in a compatible fashion, you will fail to compile with the
following error (or something similar):

.. code-block:: bash

   ../../../../observations/obs_converters/GSI2DART/dart_obs_seq_mod.f90(213): error #6284:
   There is no matching specific function for this generic function reference.   [SET_LOCATION]
   location = set_location(lon, lat, vloc, which_vert)
   -----------------^
