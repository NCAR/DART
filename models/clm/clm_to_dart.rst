PROGRAM ``clm_to_dart``
=======================

Overview
--------

``clm_to_dart`` is used to replace indeterminate values in the CLM restart file
with the *_FillValue* value specified by the variable attribute. This is necessary
to ensure that all ensemble members have the proper *_FillValue* applied to any
variable that may be part of the DART state. Only CLM restart files have variables
that require this treatment.

``clm_to_dart`` reads every variable in the input file and determines if the variable 
has snow layers by looking for a *levtot* or *levsno* dimension in the variable. 
The *SNLSNO* variable is used to detemine which layers are unused. However, when there 
is only a trace of snow, *SNLSNO* may indicate there are no snow layers in use and yet 
still have some valid values in the snow layer closest to the ground. This situation 
is confirmed by checking for non-zero values of snowcover fraction (*frac_sno*) as it 
seems to be a conservative predictor of the presence of snow. Therefore, in general, 
the indeterminate value for unused snow layer is preserved during the assimilation step,
**unless there is a trace of snow** in which case it may be adjusted.
 
See the `Discussion of Indeterminate Values`_ section for an example.

Usage
-----

The issue arises because all CLM columns have a fixed number of snow layers but only
some of the layers may actually contain snow. CLM internally tracks the number of
active snow layers in the **SNLSNO** variable and the indeterminate values in the
unused snow layers are of no consequence. To prevent these indeterminate values
from being operated on by ``filter``, DART requires these
indeterminate values to be replaced by the DART **MISSING_R8** code, which is done
automatically if the variable has a *_FillValue*.

``clm_to_dart`` **overwrites** the input file.
Consequently, the intended use is to **copy** the CLM restart file and feed the copy 
into ``clm_to_dart``. This strategy allows the copies to be read into DART and updated 
with the posterior values. The companion :doc:`dart_to_clm` uses the *_FillValue* 
values as a mask to only update the original CLM restart file with valid posterior 
values. Where the DART posterior file has *_FillValue* values, the original CLM 
values are left unchanged.

The *SNLSNO* variable is required. Any file that does not have *SNLSNO* will be 
unchanged and the program will report the error and terminate. This is why only
CLM restart files need the treatment. CLM history files do not have *SNLSNO*.

It is a useful exercise to dump the CLM restart file and the output of ``clm_to_dart``
and compare the results (for select variables). You have to use your own 
(low-resolution!) restart file.

.. container:: unix

   :: 

      cp T_SOISNO clm5_tag34.clm2.r.0003-01-01-00000.nc clm_restart.nc
      ncdump -f F -v T_SOISNO clm_restart.nc > T_SOISNO_org.txt
      ./clm_to_dart
      ncdump -f F -v T_SOISNO clm_restart.nc > T_SOISNO_new.txt
      diff T_SOISNO_org.txt T_SOISNO_new.txt


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

   ================== ==================== ========================================================================== 
   Item               Type                 Description                                                     
   ================== ==================== ========================================================================== 
   clm_restart_file   character(len=256)   Path name of the CLM restart file to be overwritten.
   verbose            integer              | Flag to control how much run-time output is created.
                                           | 0   is very little output.
                                           | 1   reports which variables are being updated.
                                           | 2   reports all 2d variables with columns as one dimension.
                                           | 3   reports all variables and columns which contain indeterminate values
                                           | *Warning*, 3 can generate a lot of output.
   ================== ==================== ==========================================================================


Modules used
------------

::

   assimilation_code/modules/utilities/netcdf_utilities_mod.f90
   assimilation_code/modules/utilities/null_mpi_utilities_mod.f90
   assimilation_code/modules/utilities/time_manager_mod.f90
   assimilation_code/modules/utilities/types_mod.f90
   assimilation_code/modules/utilities/utilities_mod.f90


Files
-----

-  clm_restart_file, specified by namelist;
-  namelistfile; ``input.nml``


Discussion of Indeterminate Values
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

To explore the most robust way to replace the indeterminate values with the 
DART-required *_FillValue*, a low-resolution run (10 deg x 15 deg) was performed 
using the **f10_f10_musgs** tag on the **release-clm5.0.34** branch.
The ``clm_to_dart_nml`` verbosity was set to 3 and the ``clm_to_dart``  run-time output was captured 
and is summarized below. The output has been pruned to improve clarity and the 
column headers have been added. Only the output for the *T_SOISNO* variable is 
shown as it clearly illustrates the situation.

The **frac_sno** variable (snow cover fraction) is used to determine when
there is a reliable trace amount of snow and the original values in the snow layer
closest to the ground are left alone. This was determined to be a more useful
indicator than **H2OSNO**, as demonstrated in the example below.  Please note that
if the **frac_sno** > 0 this indicates a trace of snow exists for that column/layer,
and the indeterminate value for *T_SOISNO* is left alone.  On the other hand, 
if **frac_sno** < 0 this indicates that particular column/layer is truly empty,
and the indeterminate value will be overwritten with the *_FillValue* attribute value,
such that DART will not adjust it.

.. note::

 The program clm_to_dart identifies indeterminate values for any variable that
 contains snow layers. These variables are identified if they have dimensions such as
 *column*, *levtot* or *levsno*. These include variables related to water mass (e.g. 
 **H2OSOI_ICE**, **H2OSOI_LIQ**), thickness (e.g. **DZSNO**, **ZSNO**, **ZISNO**) 
 and other snow layer characteristics that influence albedo (e.g. **Snw_rds**, 
 **Mss_bcpho**, **Mss_dst1**).

In the example below, the execution of ``clm_to_dart`` will overwrite the indeterminate
values within the clm restart file with the *_FillValue* attribute value to prevent DART
from performing an update to that value.  

:: 

  --------------------------------------
  Starting ... at YYYY MM DD HH MM SS = 
                  2021  6 14 16 57  4
  Program clm_to_dart
  --------------------------------------
 
   set_nml_output Echo NML values to log file only
  minval of SNLSNO is          -10
  There are          520  variables in the file.
   inquire variable number          204 varname T_SOISNO dimensions are levtot column
   clm_to_dart  variable #          204  is "T_SOISNO" - replacing indeterminate values
 
  VARIABLE  LEVEL  COLUMN   VALUE               H2OSNO                    frac_sno                  SNOW_DEPTH
  T_SOISNO  12        488   273.14999999999998  0.0000000000000000        0.0000000000000000        0.0000000000000000     
  T_SOISNO  12        503   273.14999999999998  0.0000000000000000        0.0000000000000000        0.0000000000000000     
  T_SOISNO  12        522   273.14999999999998  0.0000000000000000        0.0000000000000000        0.0000000000000000     
  T_SOISNO  12        523   273.14999999999998  0.0000000000000000        0.0000000000000000        0.0000000000000000     
  T_SOISNO  12        524   273.14999999999998  0.0000000000000000        0.0000000000000000        0.0000000000000000     
  T_SOISNO  12        525   273.14999999999998  0.0000000000000000        0.0000000000000000        0.0000000000000000     
  T_SOISNO  12        526   273.14999999999998  0.0000000000000000        0.0000000000000000        0.0000000000000000     
  T_SOISNO  12        527   273.14999999999998  0.0000000000000000        0.0000000000000000        0.0000000000000000     
  T_SOISNO  12        528   273.14999999999998  0.0000000000000000        0.0000000000000000        0.0000000000000000     
  T_SOISNO  12        533   273.04618842687802  0.0000000000000000        0.0000000000000000        0.0000000000000000     
  T_SOISNO  12        538   273.14999999999998  0.0000000000000000        0.0000000000000000        0.0000000000000000     
  T_SOISNO  12        574   272.77291836570703  0.0000000000000000        0.0000000000000000        0.0000000000000000     
  T_SOISNO  12        601   271.54492290326198  0.0000000000000000        0.0000000000000000        0.0000000000000000     
  T_SOISNO  12        603   271.54645412989498  0.0000000000000000        0.0000000000000000        0.0000000000000000     
  T_SOISNO  12        607   271.57732177585501  0.0000000000000000        0.0000000000000000        0.0000000000000000     
  T_SOISNO  12        628   272.98586559496499  0.0000000000000000        0.0000000000000000        0.0000000000000000     
  T_SOISNO  12        762   272.45467200034398  0.0000000000000000        0.0000000000000000        0.0000000000000000     
  T_SOISNO  12        763   271.34254388393202  0.0000000000000000        0.0000000000000000        0.0000000000000000     
  T_SOISNO  12        769   271.34932534583101  0.0000000000000000        0.0000000000000000        0.0000000000000000     
  T_SOISNO  12       2386   273.04040468859000  0.0000000000000000        0.0000000000000000        0.0000000000000000     
  T_SOISNO  12       2479   273.14991153505002  1.3694861766505699       0.11795998521670301        4.6284512959401197E-002
  T_SOISNO  12       2760   273.08583001204897  0.0000000000000000        0.0000000000000000        0.0000000000000000     
  T_SOISNO  12       2761   273.04455346669300  0.0000000000000000        0.0000000000000000        0.0000000000000000     
  T_SOISNO  12       2762   272.91332916458902  0.0000000000000000        0.0000000000000000        0.0000000000000000     
  T_SOISNO  12       2763   273.12159885062198  0.0000000000000000        0.0000000000000000        0.0000000000000000     
  T_SOISNO  12       2768   272.89811384494902  0.0000000000000000        0.0000000000000000        0.0000000000000000     
  T_SOISNO  12       2769   273.04499957659402  0.0000000000000000        0.0000000000000000        0.0000000000000000     
  T_SOISNO  12       2885   273.14979133425902  6.6918688785517395E-007   6.6918687813277700E-008  0.63344272495211695     
  T_SOISNO  12       2887   273.14964165253900  2.2492550993511499E-020   0.0000000000000000        0.0000000000000000     
  T_SOISNO  12       2889   272.53294466428503  1.5189620271589999E-020   0.0000000000000000        0.0000000000000000     
  T_SOISNO  12       2891   273.14945332076400  5.4222364686270699E-020   0.0000000000000000        0.0000000000000000     
  T_SOISNO  12       2893   273.14980341715801  1.4201683058442200E-020   0.0000000000000000        0.0000000000000000     
  T_SOISNO  12       2895   273.14999999999998  3.5526247972183701E-020   0.0000000000000000        0.0000000000000000     
  T_SOISNO  12       2896   273.14999404812397  9.6873949694950192E-006   9.6873931487717592E-007  0.20872425049864601     
  T_SOISNO  12       2899   273.14970403389998  6.0580004576503196E-022   0.0000000000000000        0.0000000000000000     
  T_SOISNO  12       2901   273.14957578176399  1.0670595573239900E-021   0.0000000000000000        0.0000000000000000     
  T_SOISNO  12       2903   273.14947638500303  1.1110267474343600E-021   0.0000000000000000        0.0000000000000000     
  T_SOISNO  12       2905   273.14935524945503  2.5542577505086200E-022   0.0000000000000000        0.0000000000000000     
  T_SOISNO  12       2907   273.14958719516397  1.9081446552101299E-022   0.0000000000000000        0.0000000000000000     
  T_SOISNO  12       2997   273.14999999999998  0.53787417421960904       1.2912522851932499E-002   2.1771531678050301E-003
  T_SOISNO  12       3146   273.14999999999998  0.0000000000000000        0.0000000000000000        0.0000000000000000     
  T_SOISNO  12       3329   273.14999999999998  0.0000000000000000        0.0000000000000000        0.0000000000000000     
  T_SOISNO  12       3409   273.14999999999998  0.0000000000000000        0.0000000000000000        0.0000000000000000     
  T_SOISNO  12       3410   273.14999999999998  0.0000000000000000        0.0000000000000000        0.0000000000000000     
  T_SOISNO  12       3411   273.14999999999998  1.0623461506454501E-002   1.0618761027172401E-003  0.10534145620281500     
  T_SOISNO  12       3412   273.14999999999998  1.0717589342154700E-002   1.0712805155617101E-003  0.10527378255254000     
  T_SOISNO  12       3413   273.14999999999998  1.1963762509646500E-002   1.1957631326947600E-003  0.11023440903539100     
 
  --------------------------------------
  Finished ... at YYYY MM DD HH MM SS = 
                  2021  6 14 16 57  4
  Program clm_to_dart
  --------------------------------------


References
----------

-  none, but https://escomp.github.io/ctsm-docs/versions/master/html/tech_note/Snow_Hydrology/CLM50_Tech_Note_Snow_Hydrology.html is very relevant.
