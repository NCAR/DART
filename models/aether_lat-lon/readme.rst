Aether Rectangular Grid Interface
=================================

Overview
--------

The Aether ("eether") space weather model can be implemented 
on a logically rectangular grid "lat-lon", or on a cubed-sphere grid.
This is the interface to the lat-lon version.
The model code is available on 
`GitHub <https://github.com/AetherModel/Aether>`_ .

Aether writes history and restart files.
The restart fields are divided among 2 types of files: neutrals and ions.
They are further divided into "blocks", which are subdomains of the globe.
The numbering of blocks starts in the southwest corner of the lat-lon grid 
and goes east first, then to the west end of the next row north,
and ends in the northeast corner. 
Each block has a halo around it filled with field values from neighboring blocks.
All of these need to be combined to make a single state vector for filter.
There's a unique set of these files for each member.
The restart file names reflect this information ::  

  {neutrals,ions}_mMMMM_gBBBB.nc
  MMMM = ensemble member (0-based)
  BBBB = block number (0-based)

The restart files do not have grid information in them. 
Grid information must be read from ::

  grid_gBBBB.nc

Programs ``aether_to_dart`` and ``dart_to_aether`` read the same namelist; 
``transform_state_nml``.
The fields chosen to be part of the model state are specified in 'variables'.
``Aether_to_dart`` will read the specified fields, from all the restarts
for a member plus grid files, and repackage them into an ensemble state vector file
(filter_input.nc).  Filter_input.nc has a single domain and no halos.
The field names will be transformed into CF-compliant names in filter_input.nc.

``Filter`` will read a list of variables from ``model_nml`` (not ``transform_state_nml``),
then read the ensemble of filter_input.nc files, assimilate, 
and write an ensemble of filter_output.nc files.

``Dart_to_aether`` will convert the fields' names to the CF-compliant filter names,
find those names in filter_output.nc, extract the updated field data, 
and overwrite those fields in the appropriate Aether restart files.

Namelists
---------

- The namelists are read from the file ``input.nml``. 
- Namelists start with an ampersand '&' and terminate with a slash '/'.
- Character strings that contain a '/' must be enclosed in quotes 
  to prevent them from prematurely terminating the namelist.

transform_state_nml
...................

   aether_restart_dirname 
      The directory where the Aether restart files reside, 
      and will be transformed (the "run" directory).

   nblocks_lon, nblocks_lat, nblocks_lev 
      Number of Aether domain "blocks" in the longitudinal, latitudinal, 
      and vertical directions.  Vertical is always 1 (2024-2).
      The total number of blocks (nblocks_lon x nblocks_lat x nblocks_lev)
      is defined by the number of processors used by Aether.

   variables
      The Aether fields to be transformed into a model state are specified
      in the 'variables' namelist variable in transform_state_nml.
      The following information must be provided for each field
      
         1) Aether field name
         2) which file contains the field ("neutrals" or "ions")
      
      Aether field names are not CF-compliant and are translated 
      to CF-compliant forms by aether_to_dart.  

      In ``transform_state_nml`` there is no association of DART "quantities" 
      (QTY\_\*) with fields.  
      A subset of the transformed variables to be included in the model state 
      is specified in :ref:`model_nml:variables<model_nml>`, using the CF-compliant names.
      That is where the associations with QTYs are made. 
      See the :ref:`QTY<QTY>` section, below.
      
      The neutrals restart files contain the following fields.
      The most important fields are **noted in bold text**
      
        |  **Temperature**, **velocity_east**, **velocity_north**, 
        |  velocity_up, N, O2, N2, NO, He, N_2D, N_2P, H, O_1D, CO2
      
      Similarly for the ions restart files
      
        |  **O+**, **O+_2D**, **O+_2P**, **O2+**, **N2+**, NO+, N+, He+,
        |  Temperature_bulk_ion, Temperature_electron

      In addition, there are 7 (independent) fields associated with *each* ion density
      ::
      
         - Temperature\ \(O+\)
         - velocity_parallel_east\ \(O+\)
         - velocity_parallel_north\ \(O+\)
         - velocity_parallel_up\ \(O+\)
         - velocity_perp_east\ \(O+\)
         - velocity_perp_north\ \(O+\)
         - velocity_perp_up\ \(O+\)

.. WARNING:: 
   As of this writing (2024-1-30) the electron density and solar radiation
   parameter "f10.7" are not available through the restart files, 
   even though electron temperature is.
   They may be available in the history files.
      

.. _model_nml:

model_nml
.........

template_file  
   = 'filter_input_0001.nc' is the default

variables
   Each field to be included in the state vector requires 5 descriptors:
   
      1) field name (transformed to CF-compliant)
      #) DART "quantity" to be associated with the field
      #) min value
      #) max value
      #) update the field in the restart file? {UPDATE,NO_COPY_BACK}

   The field names listed in 'variables' must be the *transformed* names,
   as found in the filter_input.nc files (see :ref:`Usage`).  
   In general the transformation does the following
   
      - Remove all '\\', '(', and ')'
      - Replace blanks with underscores
      - Replace '+' with 'pos' and '-' with 'neg'
      - For ions, move the ion name from the end to the beginning.
   
   For example 'velocity_parallel_east\\ \\(O+_2D\\)' becomes 'Opos_2D_velocity_parallel_east'.
   
.. _QTY:

   The DART QTY associated with each field is an open question,
   depending on the forward operators required for the available observations
   and on the scientific objective.   The default choices are not necessarily correct
   for your assimilation.  For the fields identified as most important
   in early Aether assimilation experiments, these are the defaults:

==============   ====================
variables        quantity (kind)
==============   ====================
Temperature      QTY_TEMPERATURE
velocity_east    QTY_U_WIND_COMPONENT
velocity_north   QTY_V_WIND_COMPONENT
Opos             QTY_DENSITY_ION_OP
O2pos            QTY_DENSITY_ION_O2P
N2pos            QTY_DENSITY_ION_N2P
O2pos_2D         QTY_DENSITY_ION_O2DP
O2pos_2P         QTY_DENSITY_ION_O2PP
==============   ====================
      
   Some fields could have one of several QTYs associated with them.  
   For example, the field 'Opos_velocity_parallel_up'
   could potentially have these existing QTYs associated with it::

   - QTY_VELOCITY_W 
   - QTY_VELOCITY_W_ION 
   - QTY_VERTICAL_VELOCITY

   It's possible that several fields could have the same QTY.
   A third possibility is that the experiment may require the creation of a new QTY.
   The example above may require something like QTY_VEL_PARALLEL_VERT_OP.

.. WARNING:: 
   The size of these parameters may be limited to 31 characters (``types_mod.f90``)

time_step_days, time_step_seconds
   = 0, 3600  The hindcast period between assimilations.

.. _Usage:

Usage
-----

The workflow and scripting for fully cycling assimilation
(ensemble hindcast, then assimilation, repeat as needed)
has not been defined yet for Aether (2024-2),
but we expect that all of the DART executables will be in a directory
which is defined in the script.
So the script will be able to run the programs using a full pathname.
In addition, all of the Aether restart files will be in a "run" directory,
which has plenty of space for the data.
The DART executables will be run in this directory using their full pathnames.

To run a more limited test (no assimilation),
which is just the transformation of files for a member (0) 
use the following csh commands, or equivalents in your preferred languange.  
These build the ``aether_to_dart`` and ``dart_to_aether`` executables
in $DART/models/aether_lat-lon/work directory.
Also in that directory, edit input.nml to set ``transform_state_nml:`` ``aether_restart_dirname``
to be the full pathname of the directory where the Aether restart and grid files are.

::

> set exec_dir = $DART/models/aether_lat-lon/work
> cd $exec_dir
> ./quick_build.sh
> cd {aether_restart_dirname}
> mkdir Orig
> cp *m0000* Orig/
> cp ${exec_dir}/input.nml .
> ${exec_dir}/aether_to_dart  0
> cp filter_input_0001.nc filter_output_0001.nc
> ${exec_dir}/dart_to_aether  0

| Compare the modified Aether restart files with those in Orig.
| The filter\_ files will contain the CF-compliant field names 
  which must be used in ``model_nml:variables``.

.. NOTE::
   Some halo parts may have no data in them because Aether currently (2024-2) 
   does not use those regions.
.. WARNING::
   The restart files have dimensions ordered such that common viewing tools 
   (e.g. ncview) may display the pictures transposed from what is expected.

