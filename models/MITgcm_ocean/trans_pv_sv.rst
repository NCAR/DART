PROGRAM ``trans_pv_sv``
=======================

| ``trans_pv_sv`` is responsible for converting the ocean model 'snapshot' files to a DART 'initial conditions' file. In
  order to do that, the valid time for the snapshot files must be calculated from several pieces of information: the
  filename contains a timestep index, the ``data``\ ``&PARM03`` namelist contains information about the amount of time
  per timestep, and the ``data.cal``\ ``&CAL_NML`` namelist contains the start time. Additionally, the grid
  characteristics must be read from ``data``\ ``&PARM04``. Consequently, the files ``data``, and ``data.cal`` as well as
  the general ``input.nml`` are needed in addition to the snapshot files.
| This program has a number of options that are driven from namelists and **one** piece of input read from STDIN: the
  integer representing the timestep index of the snapshot file set.

Usage
-----

| The output filename is hardwired to that expected by ``filter``. This example creates an output file named
  ``assim_model_state_ud`` from the following files in the local directory:
| ``S.0000000096.data``
| ``T.0000000096.data``
| ``U.0000000096.data``
| ``V.0000000096.data``
| ``Eta.0000000096.data``

.. container:: unix

   ./trans_pv_sv < 96

| 

Modules used
------------

::

   types_mod
   utilities_mod
   model_mod
   assim_model_mod
   time_manager_mod

Namelist
--------

This program has no namelist of its own, but some of the underlying modules require namelists. To avoid duplication and,
possibly, some inconsistency in the documentation, only a list of the required namelists is provided here, with a
hyperlink to the full documentation for each namelist.

+----------------------------------------------------------+----------------------------------------------------------+
| Namelist                                                 | Primary Purpose                                          |
+==========================================================+==========================================================+
| `utilities_nml <../../assimilatio                        | set the termination level and file name for the run-time |
| n_code/modules/utilities/utilities_mod.html#Namelist>`__ | log                                                      |
+----------------------------------------------------------+----------------------------------------------------------+
| `assim_model_mod_nml <../../assimilation_cod             | write DART restart files in binary or ASCII              |
| e/modules/assimilation/assim_model_mod.html#Namelist>`__ |                                                          |
+----------------------------------------------------------+----------------------------------------------------------+
| `model_nml <model_mod.html#Namelist>`__                  | write netCDF files with prognostic variables             |
+----------------------------------------------------------+----------------------------------------------------------+
| `CAL_NML <model_mod.html#namelist_cal_nml>`__            | determine start time of the ocean model                  |
+----------------------------------------------------------+----------------------------------------------------------+
| `PARM03 <model_mod.html#namelist_parm03>`__              | the amount of time per model timestep for deciphering    |
|                                                          | snapshot filenames                                       |
+----------------------------------------------------------+----------------------------------------------------------+
| `PARM04 <model_mod.html#namelist_parm04>`__              | ocean model grid parameters                              |
+----------------------------------------------------------+----------------------------------------------------------+

Files
-----

-  input namelist files: ``data, data.cal, input.nml``
-  input snapshot files: ``[S,T,U,V,Eta].nnnnnnnnnn.[data[,.meta]]``
-  output initial conditions file: ``assim_model_state_ud``

References
----------

-  none

Private components
------------------

N/A
