PROGRAM ``gen_sampling_err_table``
==================================

Overview
--------

Utility program which computes a table of values needed to apply Sampling Error Correction (SEC) during assimilation.
These values are used to correct covariances based on small sample size statistics. See reference below.

To enable the sampling error correction algorithm in ``filter``, 
set the namelist item `&assim_tools_nml : sampling_error_correction 
<../../modules/assimilation/assim_tools_mod.html#Namelist>`__ to *.true.*, 
and copy the netCDF file gen_sampling_err_table/work/sampling_error_correction_table.nc 
into the run directory.

The name of the SEC table is always ``sampling_error_correction_table.nc``. This is a NetCDF format file. 
If this file already exists in the current directory any tables for new ensemble sizes will be appended 
to the existing file. If the file does not exist a new file will be created by this tool. 
The resulting file should be copied into the current working directory when ``filter`` is run.

A file with almost 200 common ensemble sizes is distributed with the system. 
Any new ensemble sizes can be generated on demand.
Be aware that the computation can be time consuming. The job may need to be submitted to a batch system if many new
ensemble sizes are being generated, or start the job on a laptop and leave it to run overnight.

The file contains a "sparse array" of ensemble sizes. 
Only sizes which have an existing table are stored in the file so
large ensemble sizes do not require a large jump in the size of the output file.

The structure of sampling_error_correction_table.nc is shown in this example.

.. container::

   ::

      0[1095] desktop:system_simulation/work % ncdump -v ens_sizes *nc
      netcdf sampling_error_correction_table {
      dimensions:
         bins = 200 ;
         ens_sizes = UNLIMITED ; // (40 currently)
      variables:
         int count(ens_sizes, bins) ;
                 count:description = "number of samples in each bin" ;
         double true_corr_mean(ens_sizes, bins) ;
         double alpha(ens_sizes, bins) ;
                 alpha:description = "sampling error correction factors" ;
         int ens_sizes(ens_sizes) ;
                 ens_sizes:description = "ensemble size used for calculation" ;
 
         // global attributes:
                 :num_samples = 100000000 ;
                 :title = "Sampling Error Corrections for fixed ensemble sizes." ;
                 :reference = "Anderson, J., 2012: Localization and Sampling Error Correction 
                  in Ensemble Kalman Filter Data Assimilation. Mon. Wea. Rev., 140, 2359-2371, 
                  doi: 10.1175/MWR-D-11-00013.1." ;
                 :version = "$Id: gen_sampling_err_table.f90 13014 2019-03-14 21:59:07Z nancy@ucar.edu $" ;
 
       data:

       ens_sizes = 5,  6,  7,  8,  9, 10, 12, 14, 15, 16, 18, 20, 22, 24, 28, 30, 32, 36, 40, 44,
                  48, 49, 50, 52, 56, 60, 64, 70, 72, 80, 84, 88, 90, 96, 100, 120, 140, 160, 180, 200
      }

This program uses the random number generator to compute the correction factors. The generator is seeded with the
ensemble size so repeated runs of the program will generate the same values for the tables.

Namelist
--------

This namelist is read from the file ``input.nml``. Namelists start with an ampersand '&' and terminate with a slash '/'.
Character strings that contain a '/' must be enclosed in quotes to prevent them from prematurely terminating the
namelist.

::

   &gen_sampling_error_table_nml
      ens_sizes = -1
      debug = .false.
      /

Description of each namelist entry
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

ens_sizes
   *type:* integer(200)

   List of ensemble sizes to compute Sampling Error Correction tables for. These do not need to be in any particular
   order. Duplicates will be removed and any sizes which already have tables computed in the output file will be
   skipped. The file which comes with the system already has tables computed for these ensemble sizes:

   ::

      ens_sizes = 3,  4,  5,  6,  7,  8,  9, 10, 
         11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 
         21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 
         31, 32, 33, 34, 35, 36, 37, 38, 39, 40, 
         41, 42, 43, 44, 45, 46, 47, 48, 49, 50, 
         51, 52, 53, 54, 55, 56, 57, 58, 59, 60, 
         61, 62, 63, 64, 65, 66, 67, 68, 69, 70, 
         71, 72, 73, 74, 75, 76, 77, 78, 79, 80, 
         81, 82, 83, 84, 85, 86, 87, 88, 89, 90, 
         91, 92, 93, 94, 95, 96, 97, 98, 99, 100, 
         101, 102, 103, 104, 105, 106, 107, 108, 109, 110, 
         111, 112, 113, 114, 115, 116, 117, 118, 119, 120, 
         121, 122, 123, 124, 125, 126, 127, 128, 129, 130, 
         131, 132, 133, 134, 135, 136, 137, 138, 139, 140, 
         141, 142, 143, 144, 145, 146, 147, 148, 149, 150, 
         151, 152, 153, 154, 155, 156, 157, 158, 159, 160, 
         161, 162, 163, 164, 165, 166, 167, 168, 169, 170, 
         171, 172, 173, 174, 175, 176, 177, 178, 179, 180, 
         181, 182, 183, 184, 185, 186, 187, 188, 189, 190, 
         191, 192, 193, 194, 195, 196, 197, 198, 199, 200


debug
   *type:* logical

   If true print out debugging info.

Examples
--------

To add tables for ensemble sizes 220 and 256 run the program with this namelist:

.. container::

   ::

      &gen_sampling_error_table_nml
         ens_sizes = 220, 256,
         debug = .false.
         /

Modules used
------------

::

   types_mod
   utilities_mod
   random_seq_mod
   netcdf

Files
-----

-  output file is always ``sampling_error_correction_table.nc`` If one exists new ensemble sizes will be appended. If
   it doesn't exist a new file will be created. This is a NetCDF format file.

References
----------

-  Ref: Anderson, J., 2012: Localization and Sampling Error Correction in Ensemble Kalman Filter Data Assimilation. Mon.
   Wea. Rev., 140, 2359-2371, doi: 10.1175/MWR-D-11-00013.1.
