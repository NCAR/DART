############################
forward operator test README
############################

Contents
========

#. `Overview`_
#. `rttov_test.f90`_
#. `rttov_unit_tests.f90`_
#. `make_COS_input`_
#. `make_assim_list`_
#. `Terms of Use`_

Overview
========

The ``developer_tests/forward_operators`` directory contains the testing
framework for three kinds of tests:

- `rttov_test.f90` which tests basic functionality of the RTTOV radiance forward operators for AMSUA and AIRS observations (MW and IR, respectively)

  - create a dummy rttov_sensor_db_file to read
  - check the IR instrument ID
  - check the MW instrument ID
  - exercise the IR (direct) forward operator
  - exercise the MW (scatt) forward operator

- `rttov_unit_tests.f90` performs a host of unit tests. 

  - If run with a TERMLEVEL of 3, all tests will be completed even if previous tests are not successful.
  
- `make_COS_input` and `make_assim_list` creates input for create_obs_sequence (COS) and the appropriate namelist settings to test the forward operator code.  

Different sets of observations are grouped into separate files based on certain criteria - are they atmospheric observations, oceanic ... do they require special metadata, etc. The following files are intended to be supplied as input to `make_COS_input` and will result in a text file that will generate an observation sequence file when used as input to `create_obs_sequence`. 

- all_atm_obs_types
- all_commoncode_atm_obs_types
- all_f90s
- all_fwdop_atm_obs_types
- all_obs_types
- forward_op_code
- no_special_forward_op_code

See the `make_COS_input`_ section for more detail. 

rttov_test.f90
==============

This test requires several coefficient files that are not part of the default
set provided by the RTTOV 12.3 distribution. Specifically:

- `rtcoef_eos_2_amsua.dat`
- `rtcoef_eos_2_airs.H5`
- `mietable_eos_amsua.dat` (same file as `mietable_noaa_amsua.dat`)

These coefficient files may be downloaded by using the `rtcoef_rttov12/rttov_coef_download.sh`
script provided in the RTTOV distribution.

rttov_unit_tests.f90
====================

These unit tests are best run with a TERMLEVEL of 3, which allows DART to
continue past errors that would otherwise be fatal.
If any of the unit tests are unable to start, the error code from
*rttov_unit_tests* is 102.  This is to give an error for *test_dart.csh* to detect.


.. list-table::


  * - **Test**
    - **Pass** 
    - **Fail** 
  * - metadata growth
    - metadata arrays grow correctly as observations are added
    - incorrect metadata length
  * - metadata content  
    - metadata arrays contain correct data
    - incorrect data
  

make_COS_input
==============

*make_COS_input* takes one filename as an argument and creates a text file that
can be used as input for *create_obs_sequence*. The output text file has
a name based on the input filename. For example:

.. code-block::

  <prompt> ./make_COS_input forward_op_code
   ready to run create_obs_sequence < forward_op_code_COS.in

*create_obs_sequence* must be created with the `preprocess_nml`
settings to support the observation definitions required by the input file.

make_assim_list
===============

*make_assim_list* is a follow-on step to *make_COS_input* and simply creates
the text for the `input.nml:filter_nml:assimilate_these_obs` variable.

.. code-block::

  <prompt> forward_operators > ./make_assim_list forward_op_code
   created forward_op_code_obskind.nml
   add this section to your &obs_kind_nml in input.nml
   <prompt> head -n 10 forward_op_code_obskind.nml
   assimilate_these_obs_types =
   'ACARS_DEWPOINT,',
   'ACARS_RELATIVE_HUMIDITY,',
   'AIRCRAFT_DEWPOINT,',
   'AIRCRAFT_RELATIVE_HUMIDITY,',
   'AIREP_DEWPOINT,',
   'AIRS_DEWPOINT,',
   'AIRS_RELATIVE_HUMIDITY,',
   'AMDAR_DEWPOINT,',
   'AMSR_TOTAL_PRECIPITABLE_WATER,',


Terms of Use
============

|Copyright| University Corporation for Atmospheric Research

Licensed under the `Apache License, Version 2.0
<http://www.apache.org/licenses/LICENSE-2.0>`__. Unless required by applicable
law or agreed to in writing, software distributed under this license is
distributed on an "as is" basis, without warranties or conditions of any kind,
either express or implied.

.. |Copyright| unicode:: 0xA9 .. copyright sign
