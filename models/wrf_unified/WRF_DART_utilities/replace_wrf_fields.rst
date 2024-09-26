PROGRAM ``replace_wrf_fields``
==============================

Overview
--------

Program to copy various fields from one WRF netCDF file to another.

There are many existing utilities to process netCDF files, i.e. the NCO operators and NCL scripts, which have more
functionality than this program. The only purpose for having this one is that it is a standalone program with no
prerequisites or dependencies other than the netCDF libraries. If you already have other tools available they can do the
same functions that this program does.

This program copies the given data fields from the input file to the output file, failing if their sizes, shapes, or
data types do not match exactly. The expected use is to copy fields which are updated by the WRF program but are not
part of the DART state vector, for example, sea surface temperature or soil fields. After DART has updated the WRF
restart ``wrfinput_d01`` file, this program can be used to update other fields in the file before running the model.

Namelist
--------

This namelist is read from the file ``input.nml``. Namelists start with an ampersand '&' and terminate with a slash '/'.
Character strings that contain a '/' must be enclosed in quotes to prevent them from prematurely terminating the
namelist.

::

   &replace_wrf_fields_nml
      fieldnames = 'SST',
      fieldlist_file = '',
      fail_on_missing_field = .true.
      debug = .false.,
      /

| 

.. container::

   +-----------------------+------------------------+-------------------------------------------------------------------+
   | Item                  | Type                   | Description                                                       |
   +=======================+========================+===================================================================+
   | fieldnames            | character(len=129) (:) | An array of ASCII field names to be copied from the input netCDF  |
   |                       |                        | file to the output netCDF file. The names must match exactly, and |
   |                       |                        | the size and shape of the data must be the same in the input and  |
   |                       |                        | output files for the data to be copied. If the field names are    |
   |                       |                        | set here, the fieldlist_file item must be ' '.                    |
   +-----------------------+------------------------+-------------------------------------------------------------------+
   | fieldlist_file        | character(len=129)     | An alternative to an explicit list of field names to copy. This   |
   |                       |                        | is a single string, the name of a file which contains a single    |
   |                       |                        | field name, one per line. If this option is set, the fieldnames   |
   |                       |                        | namelist item must be ' '.                                        |
   +-----------------------+------------------------+-------------------------------------------------------------------+
   | fail_on_missing_field | logical                | If any fields in the input list are not found in either the input |
   |                       |                        | or output netcdf files, fail if this is set to true. If false, a  |
   |                       |                        | warning message will be printed but execution will continue.      |
   +-----------------------+------------------------+-------------------------------------------------------------------+
   | debug                 | logical                | If true, print out debugging messages about which fields are      |
   |                       |                        | found in the input and output files.                              |
   +-----------------------+------------------------+-------------------------------------------------------------------+

| 

Modules used
------------

::

   types_mod
   utilities_mod
   parse_args_mod

Files
-----

-  input namelist ; ``input.nml``
-  Input - output WRF state netCDF files; ``wrfinput_d01, wrfinput_d02, ...``
-  fieldlist_file (if specified in namelist)

File formats
~~~~~~~~~~~~

This utility works on any pair of netCDF files, doing a simple read and copy from one to the other.

Error codes and conditions
~~~~~~~~~~~~~~~~~~~~~~~~~~

+--------------------+---------------------------------------------------------+----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------+
|       Routine      |                         Message                         |                                                                                            Comment                                                                                           |
+====================+=========================================================+==============================================================================================================================================================================================+
| replace_wrf_fields | Usage: echo infile.nc outfile.nc | ./replace_wrf_fields | The program did not read 2 filenames from the console.                                                                                                                                       |
+--------------------+---------------------------------------------------------+----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------+
| replace_wrf_fields | cannot specify both fieldnames and fieldlist_file       | In the namelist you must either specify an explicit list of fieldnames to copy between the files, or give a single filename which contains the list of field names. You cannot specify both. |
+--------------------+---------------------------------------------------------+----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------+
| replace_wrf_fields | field not found in input/output file                    | If 'fail_on_missing_field' is true in the namelist and a field is not found in either the input or output file.                                                                              |
+--------------------+---------------------------------------------------------+----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------+
| replace_wrf_fields | field does not match                                    | If the input and output files have different sizes, number of dimensions, or data types, the program cannot copy the data.                                                                   |
+--------------------+---------------------------------------------------------+----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------+

References
----------

-  none
