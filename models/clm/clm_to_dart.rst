PROGRAM ``clm_to_dart``
=======================

Overview
--------

``clm_to_dart`` is used to replace indeterminate values in the CLM restart file
with the *_FillValue* value specified as the variable attribute. This is necessary
to ensure that all ensemble members have the proper *_FillValue* applied to any
variable that may be part of the DART state. Only CLM *restart* files have variables
that require this treatment.

The issue arises because all CLM columns have a fixed number of snow layers but only
some of the layers may actually contain snow. CLM internally tracks the number of
active snow layers in the **SNLSNO** variable and the indeterminate values in the
unused snow layers are of no consequence. The internal logic in DART requires these
indeterminate values to be replaced by the DART **MISSING_R8** code, which is done
automatically if the variable has a *_FillValue*.

``clm_to_dart`` **overwrites** the input CLM restart file.
As such, the intended use of ``clm_to_dart`` is to copy the CLM restart file and 
then feed the copy into ``clm_to_dart``.  This preserves the original CLM restart 
file which will ultimately be updated with the posterior values 
**with the DART MISSING_R8 values being replaced with whatever is in the 
original CLM restart file**.

The **frac_sno** variable (snow cover fraction) is used to determine when
there is a reliable trace amount of snow and the original values in the snow layer
closest to the ground are left alone. This was determined to be a more useful
indicator than **H2OSNO**, as evidenced by the values from a CLM5.0.34 run at low resolution:

::

    variable      snow layer   column      value              H2OSNO(column)          frac_sno(column)        SNOW_DEPTH(column)
    T_SOISNO          12         488   273.150000000000       0.000000000000000E+000  0.000000000000000E+000  0.000000000000000E+000
    T_SOISNO          12         503   273.150000000000       0.000000000000000E+000  0.000000000000000E+000  0.000000000000000E+000
    T_SOISNO          12         522   273.150000000000       0.000000000000000E+000  0.000000000000000E+000  0.000000000000000E+000
    T_SOISNO          12         523   273.150000000000       0.000000000000000E+000  0.000000000000000E+000  0.000000000000000E+000
    T_SOISNO          12         524   273.150000000000       0.000000000000000E+000  0.000000000000000E+000  0.000000000000000E+000
    T_SOISNO          12         525   273.150000000000       0.000000000000000E+000  0.000000000000000E+000  0.000000000000000E+000
    T_SOISNO          12         526   273.150000000000       0.000000000000000E+000  0.000000000000000E+000  0.000000000000000E+000
    T_SOISNO          12         527   273.150000000000       0.000000000000000E+000  0.000000000000000E+000  0.000000000000000E+000
    T_SOISNO          12         528   273.150000000000       0.000000000000000E+000  0.000000000000000E+000  0.000000000000000E+000
    T_SOISNO          12         533   273.046188426878       0.000000000000000E+000  0.000000000000000E+000  0.000000000000000E+000
    T_SOISNO          12         538   273.150000000000       0.000000000000000E+000  0.000000000000000E+000  0.000000000000000E+000
    T_SOISNO          12         574   272.772918365707       0.000000000000000E+000  0.000000000000000E+000  0.000000000000000E+000
    T_SOISNO          12         601   271.544922903262       0.000000000000000E+000  0.000000000000000E+000  0.000000000000000E+000
    T_SOISNO          12         603   271.546454129895       0.000000000000000E+000  0.000000000000000E+000  0.000000000000000E+000
    T_SOISNO          12         607   271.577321775855       0.000000000000000E+000  0.000000000000000E+000  0.000000000000000E+000
    T_SOISNO          12         628   272.985865594965       0.000000000000000E+000  0.000000000000000E+000  0.000000000000000E+000
    T_SOISNO          12         762   272.454672000344       0.000000000000000E+000  0.000000000000000E+000  0.000000000000000E+000
    T_SOISNO          12         763   271.342543883932       0.000000000000000E+000  0.000000000000000E+000  0.000000000000000E+000
    T_SOISNO          12         769   271.349325345831       0.000000000000000E+000  0.000000000000000E+000  0.000000000000000E+000
    T_SOISNO          12        2386   273.040404688590       0.000000000000000E+000  0.000000000000000E+000  0.000000000000000E+000
    T_SOISNO          12        2479   273.149911535050        1.36948617665057       0.117959985216703       4.628451295940120E-002
    T_SOISNO          12        2760   273.085830012049       0.000000000000000E+000  0.000000000000000E+000  0.000000000000000E+000
    T_SOISNO          12        2761   273.044553466693       0.000000000000000E+000  0.000000000000000E+000  0.000000000000000E+000
    T_SOISNO          12        2762   272.913329164589       0.000000000000000E+000  0.000000000000000E+000  0.000000000000000E+000
    T_SOISNO          12        2763   273.121598850622       0.000000000000000E+000  0.000000000000000E+000  0.000000000000000E+000
    T_SOISNO          12        2768   272.898113844949       0.000000000000000E+000  0.000000000000000E+000  0.000000000000000E+000
    T_SOISNO          12        2769   273.044999576594       0.000000000000000E+000  0.000000000000000E+000  0.000000000000000E+000
    T_SOISNO          12        2885   273.149791334259       6.691868878551743E-007  6.691868781327770E-008  0.633442724952117
    T_SOISNO          12        2887   273.149641652539       2.249255099351152E-020  0.000000000000000E+000  0.000000000000000E+000
    T_SOISNO          12        2889   272.532944664285       1.518962027158998E-020  0.000000000000000E+000  0.000000000000000E+000
    T_SOISNO          12        2891   273.149453320764       5.422236468627066E-020  0.000000000000000E+000  0.000000000000000E+000
    T_SOISNO          12        2893   273.149803417158       1.420168305844218E-020  0.000000000000000E+000  0.000000000000000E+000
    T_SOISNO          12        2895   273.150000000000       3.552624797218375E-020  0.000000000000000E+000  0.000000000000000E+000
    T_SOISNO          12        2896   273.149994048124       9.687394969495016E-006  9.687393148771761E-007  0.208724250498646
    T_SOISNO          12        2899   273.149704033900       6.058000457650322E-022  0.000000000000000E+000  0.000000000000000E+000
    T_SOISNO          12        2901   273.149575781764       1.067059557323987E-021  0.000000000000000E+000  0.000000000000000E+000
    T_SOISNO          12        2903   273.149476385003       1.111026747434362E-021  0.000000000000000E+000  0.000000000000000E+000
    T_SOISNO          12        2905   273.149355249455       2.554257750508619E-022  0.000000000000000E+000  0.000000000000000E+000
    T_SOISNO          12        2907   273.149587195164       1.908144655210132E-022  0.000000000000000E+000  0.000000000000000E+000
    T_SOISNO          12        2997   273.150000000000       0.537874174219609       1.291252285193245E-002  2.177153167805032E-003
    T_SOISNO          12        3146   273.150000000000       0.000000000000000E+000  0.000000000000000E+000  0.000000000000000E+000
    T_SOISNO          12        3329   273.150000000000       0.000000000000000E+000  0.000000000000000E+000  0.000000000000000E+000
    T_SOISNO          12        3409   273.150000000000       0.000000000000000E+000  0.000000000000000E+000  0.000000000000000E+000
    T_SOISNO          12        3410   273.150000000000       0.000000000000000E+000  0.000000000000000E+000  0.000000000000000E+000
    ...


Usage
-----

- copy a CLM restart file to 'clm_restart.nc' - **warning** - this file *will* be modified.

- set the verbosity flag in the *input.nml*

- run ``clm_to_dart``



Namelist
--------

Namelists start with an ampersand '&' and terminate with a slash '/'. 
Character strings that contain a '/' must be enclosed in quotes to prevent 
them from prematurely terminating the namelist. These are the defaults:

::

   &clm_to_dart_nml
     clm_restart_file  = 'clm_restart.nc'
     verbose           = 0
    /


.. container::

   +------------------+--------------------+-----------------------------------------------------------------+
   | Item             | Type               | Description                                                     |
   +==================+====================+=================================================================+
   | clm_restart_file | character(len=256) | Path name of the CLM restart file to be overwritten.            |
   |                  |                    | to be preprocessed. The default is                              |
   +------------------+--------------------+-----------------------------------------------------------------+
   | verbose          | integer            | Flag to control how much run-time output is created.            |
   |                  |                    | 0   is very little output                                       |
   |                  |                    | 1   reports which variables are being updated                   |
   |                  |                    | 2   reports all 2d variables with columns as one dimension      |
   |                  |                    | 3   reports the variables and columns that have traces of snow. |
   |                  |                    |     Warning, this can generate a lot of output.                 |
   +------------------+--------------------+-----------------------------------------------------------------+


Modules used
------------

::

   types_mod
   utilities_mod
   netcdf_utilities_mod
   time_manager_mod
   null_mpi_utilities_mod


Files
-----

-  clm_restart_file, specified by namelist;
-  namelistfile; ``input.nml``

References
----------

-  none, but https://escomp.github.io/ctsm-docs/versions/master/html/tech_note/Snow_Hydrology/CLM50_Tech_Note_Snow_Hydrology.html
         is very relevant.
