Aether Rectangular Grid Interface
=================================

Overview
--------

The `Aether` ("eether") space weather model can be implemented 
on a logically rectangular grid "lat-lon", 
or on an the cubed-sphere grid (see ../aether_cubed_shere).
This is the interface to the lat-lon version.

.. Aether: https://aetherdocumentation.readthedocs.io/en/latest/

Aether writes history and restart files, with some overlap of the fields (?).
The restart fields are divided among 2 types of files: neutrals and ions.
They are further divided into "blocks", which are subdomains of the globe.
Blocks start in the southwest corner of the lat/lon grid and go east first, 
then to the west end of the next row north and end in the northeast corner. 
Each block has a halo around it filled with field values from neighboring blocks.
All of these need to be combined to make a single state vector for filter.
There's a unique set of these files for each member.
The restart file names reflect this information:  

|   {neutrals,ions}_mMMMM_gBBBB.nc
|   MMMM = ensemble member (0-based)
|   BBBB = block number (0-based)

The restart files do not have grid information in them, which must be read from
   grid_gBBBB.nc

Aether_to_dart and dart_to_aether read the same namelist; transform_state_nml.
The fields chosen to be part of the model state are specified in ``variables``.
Program aether_to_dart will read the specified fields, from all the restarts
for a member plus grid files, and repackage them into an ensemble state vector file
(filter_input.nc), which has a single domain and no halos.
The field names will be transformed into CF-compliant names in filter_input.nc.

Filter will read the ensemble of filter_input.nc files, assimilate, 
and write an ensemble of filter_output.nc files.

Dart_to_aether will convert the fields' names to the CF-compliant filter names,
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
      and vertical directions.  (vertical is always 1 as of 2024-2)

   variables
      The Aether fields to be included in the model state are specified
      in the ``variables`` namelist variable in transform_state_nml.
      The following information must be provided for each field
      
         1) Aether field name
         2) which file contains the field ("neutrals" or "ions")
      
      Aether field names are not CF-compliant and are translated 
      to CF-compliant forms by aether_to_dart.
      The suggested DART quantity to associate with some fields are listed
      in ./aether_to_dart.nml.
      
      The neutrals restart files contain the following fields.
      The most important fields are **highlighted**
      
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
   As of this writing (2024-1-30) the electron density is not available 
   through the restart files, even though electron temperature is.
   It can be written to the history files.
      

model_nml
.........

filter_io_filename  
   = 'if other than filter_input_0001.nc'

variables
   Each field to be included in the state vector requires 5 descriptors:
   
      1) field name (transformed to CF-compliant)
      #) DART "quantity" to be associated with the field
      #) min value
      #) max value
      #) update the field in the restart file? {UPDATE,NO_COPY_BACK}

   The field names listed in ``variables`` must be the *transformed* names,
   as found in the filter_input.nc files (see :ref:`Usage`).  
   In general the transformation does the following
   
      - Remove all '\\', '(', and ')'
      - Replace blanks with underscores
      - Replace '+' with 'pos' and '-' with 'neg'
      - For ions, move the ion name from the end to the beginning.
   
   For example 'velocity_parallel_east\\ \\(O+_2D\\)' becomes 'Opos_2D_velocity_parallel_east'.
   
time_step_days, time_step_seconds
   = 0, 3600  The hindcast period between assimilations.

.. _Usage:

Usage
-----

To test the transformation of files for member 0:
::

> cd {aether_restart_dirname}
> mkdir Orig
> cp *m0000* Orig
> ./aether_to_dart  0
> cp filter_input_0001.nc filter_output_0001.nc
> ./dart_to_aether  0

| The filter\_ files now contain the CF-compliant field names which must be used in model_nml:variables.
| Compare the modified Aether restart files with those in Orig.
| (Some halo regions may have no data in them because Aether currently (2024-2) does not use those regions.)

