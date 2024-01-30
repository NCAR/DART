Aether Rectangular Grid Interface
=================================

Overview
--------

The Aether ("eether") space weather model (TODO: reference) can be implemented 
on a logically rectangular grid"lon-lat", 
or on an the cubed-sphere grid (see ../aether_cubed_shere).
This is the interface to the lon-lat version.

Aether writes history and restart files, with some overlap of the fields (?).
The restart fields are divided among 2 types of files: neutrals and ions.
They are further divided into "blocks", which are subdomains of the globe.
All of these need to be combined to make a single state vector for filter.
There's a unique set of these files for each member.
The restart file names reflect this information:  
   {neutrals,ions}_mMMMM_gBBBB.nc
   MMMM = ensemble member (0-based)
   BBBB = block number (0-based)
These files do not have grid information in them, which must be read from
   grid_gBBBB.nc

Program aether_to_dart will read a selection of fields from all the restart 
and grid files for a member and repackage them into an ensemble state vector 
(filter_input.nc).

Filter will read the ensemble of filter_input.nc files, assimilate, 
and write an ensemble of filter_output.nc files.

Program dart_to_aether will extract the updated field data from them
and overwrite those fields in the Aether restart files.

Namelists
---------

- The namelists are read from the file ``input.nml``. 
- Namelists start with an ampersand '&' and terminate with a slash '/'.
- Character strings that contain a '/' must be enclosed in quotes 
  to prevent them from prematurely terminating the namelist.

aether_to_dart_nml
.....................

The Aether fields to be included in the model state are specified
in the ``variables`` namelist variable.
The following information must be provided for each field:

1) Aether field name
#  DART "quantity" to be associated with the field
#  max value
#  min value
#  which file contains the field ("neutrals" or "ions")
#  whether the field should be updated in the assimilation

Aether field names are not CF-compliant and are translated 
to CF-compliant forms by aether_to_dart.
The suggested DART quantity to associate with some fields are listed
in ./aether_to_dart.nml.

The neutrals restart files contain the following fields.
The most important fields are **highlighted**.::
   **Temperature**, **velocity_east**, **velocity_north**, 
   velocity_up, N, O2, N2, NO, He, N_2D, N_2P, H, O_1D, CO2

Similarly for the ions restart files: ::
   **O+**, **O+_2D**, **O+_2P**, **O2+**, **N2+**, NO+, N+, He+,
     Temperature_bulk_ion, Temperature_electron
   **NOTE** As of this writing (2024-1-30) the electron density is not available 
   through the restart files, even though electron temperature is.
   It can be written to the history files.

In addition, there are 7 (independent) fields associated with *each* ion density:

- Temperature\ \(O+\)
- velocity_parallel_east\ \(O+\)
- velocity_parallel_north\ \(O+\)
- velocity_parallel_up\ \(O+\)
- velocity_perp_east\ \(O+\)
- velocity_perp_north\ \(O+\)
- velocity_perp_up\ \(O+\)


dart_to_aether_nml
.....................

The ``variables`` in this namelist must match the list in aether_to_dart_nml.
Dart_to_aether_nml will convert these fields names to the CF-compliant filter names,
find those names in filter_output.nc, and transfer the updated fields
from filter_output.nc to the Aether appropriate restart files.


model_nml
.........

The fields listed in ``variables`` must be the translated names,
as found in the filter_input.nc files.  
In general the transformation does the following:

- Remove all '\', '(', and ')'
- Replace blanks with underscores
- Replace '+' with 'pos' and '-' with 'neg'
- For ions, move the ion name from the end to the beginning.

For example 'velocity_parallel_east\ \(O+_2D\)' becomes
'Opos_2D_velocity_parallel_east'
::

   &model_nml 
    /

| 

