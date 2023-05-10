program ``model_mod_check``
===========================

Overview
--------

``model_mod_check`` tests some of the more fundamental routines in any ``model_mod``. This is intended to be used when
adding a new model to DART - test the pieces as they are written. As such, this program is meant to be hacked up and
customized to your own purpose. Right now, it reads in model netCDF file(s) - one per domain/nest/whatever - and writes
out files, queries the metdata, etc. It also exercises ``static_init_model()``, which is the first routine to get right
...

Namelist
--------

This namelist is read from the file ``input.nml``. Namelists start with an ampersand '&' and terminate with a slash '/'.
Character strings that contain a '/' must be enclosed in quotes to prevent them from prematurely terminating the
namelist.

::

   &model_mod_check_nml
      num_ens               = 1
      single_file           = .FALSE.
      input_state_files     = 'null'
      output_state_files    = 'null'
      all_metadata_file     = 'metadata.txt'

      test1thru             = 7
      run_tests             = -1

      x_ind                 = -1
      loc_of_interest       = -1.0, -1.0, -1.0
      quantity_of_interest  = 'NONE'

      interp_test_dlon      = 10.0
      interp_test_dlat      = 10.0
      interp_test_dvert     = 10.0

      interp_test_lonrange  = 0.0, 120.0
      interp_test_latrange  = 0.0, 120.0
      interp_test_vertrange = 0.0, 100.0

      interp_test_dx        = -888888.0
      interp_test_dy        = -888888.0
      interp_test_dz        = -888888.0

      interp_test_xrange    = -888888.0, -888888.0
      interp_test_yrange    = -888888.0, -888888.0
      interp_test_zrange    = -888888.0, -888888.0

      interp_test_vertcoord = 'VERTISHEIGHT'
      verbose               = .FALSE.
      /

| 

.. container::

   +------------------------+------------------------+----------------------------------------------+
   | Item                   | Type                   | Description                                  |
   +========================+========================+==============================================+
   | num_ens                | integer                | Provided for future use. Must be 1.          |
   |                        |                        | Ultimately, The number of ensemble           |
   |                        |                        | members you would like to read in for        |
   |                        |                        | testing.                                     |
   +------------------------+------------------------+----------------------------------------------+
   | single_file            | logical                | If .TRUE. all members are stored in a        |
   |                        |                        | single restart file.                         |
   +------------------------+------------------------+----------------------------------------------+
   | input_state_files(:)   | character(len=256)     | The name(s) of the NetCDF file(s)            |
   |                        |                        | containing the model states, one per         |
   |                        |                        | domain. If num_ens > 1 and not               |
   |                        |                        | single_file, specify a filename for          |
   |                        |                        | each ensemble member (num_ens). If           |
   |                        |                        | you have both multiple ensemble              |
   |                        |                        | members in separate files AND                |
   |                        |                        | multiple domains, specify all the            |
   |                        |                        | ensemble member filenames for domain         |
   |                        |                        | 1, then all the ensemble member              |
   |                        |                        | filenames for domain 2, etc.                 |
   +------------------------+------------------------+----------------------------------------------+
   | output_state_files(:)  | character(len=256)     | The name(s) of the output NetCDF             |
   |                        |                        | file(s) for testing IO, one per              |
   |                        |                        | domain. If num_ens > 1 and not               |
   |                        |                        | single_file, specify a filename for          |
   |                        |                        | each ensemble member (num_ens). If           |
   |                        |                        | you have both multiple ensemble              |
   |                        |                        | members in separate files AND                |
   |                        |                        | multiple domains, specify all the            |
   |                        |                        | ensemble member filenames for domain         |
   |                        |                        | 1, then all the ensemble member              |
   |                        |                        | filenames for domain 2, etc.                 |
   +------------------------+------------------------+----------------------------------------------+
   | all_metadata_file      | character(len=256)     | Test 6 produces an exhaustive list of        |
   |                        |                        | metadata for EVERY element in the            |
   |                        |                        | DART state vector. The metadata get          |
   |                        |                        | written to this file name.                   |
   +------------------------+------------------------+----------------------------------------------+
   | x_ind                  | integer(i8)            | An integer index into the DART state         |
   |                        |                        | vector. This will be used to test the        |
   |                        |                        | metadata routines. Answers questions         |
   |                        |                        | about location, what variable type is        |
   |                        |                        | stored there, etc.                           |
   +------------------------+------------------------+----------------------------------------------+
   | loc_of_interest        | real(r8), dimension(3) | The lat/lon/level for a                      |
   |                        |                        | **particular** location. Used in Test        |
   |                        |                        | 4, the single-point interpolation            |
   |                        |                        | test. Indirectly tests the routine to        |
   |                        |                        | find the closest gridpoint.                  |
   +------------------------+------------------------+----------------------------------------------+
   | quantity_of_interest   | character(len=32)      | Specifies the QUANTITY of the model          |
   |                        |                        | state to use in Tests 4, 5, and 7.           |
   +------------------------+------------------------+----------------------------------------------+
   | interp_test_dlon       | real(r8)               | The distance (measured in degrees) on        |
   |                        |                        | the longitude interpolation grid.            |
   |                        |                        | Ignored if interpolating with                |
   |                        |                        | cartesian coordinates. Used in Test 5.       |
   +------------------------+------------------------+----------------------------------------------+
   | interp_test_dlat       | real(r8)               | The distance (measured in degrees) on        |
   |                        |                        | the latitude interpolation grid.             |
   |                        |                        | Ignored if interpolating with                |
   |                        |                        | cartesian coordinates. Used in Test 5.       |
   +------------------------+------------------------+----------------------------------------------+
   | interp_test_dvert      | real(r8)               | The distance (measured in                    |
   |                        |                        | interp_vertcoord) on the vertical            |
   |                        |                        | interpolation grid. Ignored if               |
   |                        |                        | interpolating with cartesian                 |
   |                        |                        | coordinates. Used in Test 5.                 |
   +------------------------+------------------------+----------------------------------------------+
   | interp_test_lonrange   | real(r8)               | The range of y to be tested with             |
   |                        |                        | model_interpolate, with spacing              |
   |                        |                        | ``interp_test_dlon``. Ignored if             |
   |                        |                        | interpolating with cartesian                 |
   |                        |                        | coordinates. Used in Test 5.                 |
   +------------------------+------------------------+----------------------------------------------+
   | interp_test_latrange   | real(r8)               | The range of y to be tested with             |
   |                        |                        | model_interpolate, with spacing              |
   |                        |                        | ``interp_test_dlat``. Ignored if             |
   |                        |                        | interpolating with cartesian                 |
   |                        |                        | coordinates. Used in Test 5.                 |
   +------------------------+------------------------+----------------------------------------------+
   | interp_test_vertrange  | real(r8)               | The range in the vertical direction          |
   |                        |                        | to be tested with model_interpolate,         |
   |                        |                        | with spacing ``interp_test_dvert``.          |
   |                        |                        | Ignored if interpolating with                |
   |                        |                        | cartesian coordinates. Used in Test 5.       |
   +------------------------+------------------------+----------------------------------------------+
   | interp_test_dx         | real(r8)               | The interval on the x axis of the            |
   |                        |                        | interpolation grid. This is used in          |
   |                        |                        | Test 5 for models with                       |
   |                        |                        | threed_cartesian coordinates.                |
   +------------------------+------------------------+----------------------------------------------+
   | interp_test_dy         | real(r8)               | The interval on the y axis of the            |
   |                        |                        | interpolation grid. This is used in          |
   |                        |                        | Test 5 for models with                       |
   |                        |                        | threed_cartesian coordinates.                |
   +------------------------+------------------------+----------------------------------------------+
   | interp_test_dz         | real(r8)               | The interval on the z axis of the            |
   |                        |                        | interpolation grid. This is used in          |
   |                        |                        | Test 5 for models with                       |
   |                        |                        | threed_cartesian coordinates.                |
   +------------------------+------------------------+----------------------------------------------+
   | interp_test_xrange     | real(r8)               | The range of x to be tested with             |
   |                        |                        | model_interpolate in Test 5, with            |
   |                        |                        | spacing ``interp_test_dx``.                  |
   +------------------------+------------------------+----------------------------------------------+
   | interp_test_yrange     | real(r8)               | The range of y to be tested with             |
   |                        |                        | model_interpolate in Test 5, with            |
   |                        |                        | spacing ``interp_test_dy``.                  |
   +------------------------+------------------------+----------------------------------------------+
   | interp_test_zrange     | real(r8)               | The range in the vertical direction          |
   |                        |                        | to be tested with model_interpolate          |
   |                        |                        | in Test 5, with spacing                      |
   |                        |                        | ``interp_test_dz``.                          |
   +------------------------+------------------------+----------------------------------------------+
   | interp_test_vertcoord  | character(len=32)      | Specifies the vertical coordinate            |
   |                        |                        | system to use during the                     |
   |                        |                        | interpolation tests. Valid values            |
   |                        |                        | are: 'VERTISHEIGHT',                         |
   |                        |                        | 'VERTISPRESSURE',                            |
   |                        |                        | 'VERTISLEVEL', and                           |
   |                        |                        | 'VERTISSCALEHEIGHT'.                         |
   +------------------------+------------------------+----------------------------------------------+
   | test1thru              | integer                | If ``test1thru > 0``, specifies the          |
   |                        |                        | last test to be performed. All tests         |
   |                        |                        | get performed sequentially. If               |
   |                        |                        | ``test1thru < 0``, ``run_tests`` is          |
   |                        |                        | used to specify the tests to perform.        |
   |                        |                        |                                              |
   |                        |                        | +------+--------------------------------+    |
   |                        |                        | | test | summary                        |    |
   |                        |                        | +======+================================+    |
   |                        |                        | | 0    | Mandatory. Tests               |    |
   |                        |                        | |      | ``static_init_model()``        |    |
   |                        |                        | |      | by calling                     |    |
   |                        |                        | |      | ``static_init_assim_model()``. |    |
   |                        |                        | |      | Reads ``input.nml``            |    |
   |                        |                        | |      | ``&model_nml``                 |    |
   |                        |                        | +------+--------------------------------+    |
   |                        |                        | | 1    | Tests                          |    |
   |                        |                        | |      | ``get_model_size()`` and       |    |
   |                        |                        | |      | reports on the makeup of       |    |
   |                        |                        | |      | the DART state vector.         |    |
   |                        |                        | +------+--------------------------------+    |
   |                        |                        | | 2    | Reads and writes a             |    |
   |                        |                        | |      | restart file.                  |    |
   |                        |                        | +------+--------------------------------+    |
   |                        |                        | | 3    | Tests                          |    |
   |                        |                        | |      | ``get_state_meta_data()``      |    |
   |                        |                        | |      | for a single index into        |    |
   |                        |                        | |      | the DART state. Helps          |    |
   |                        |                        | |      | determine if the state         |    |
   |                        |                        | |      | vector is constructed          |    |
   |                        |                        | |      | correctly.                     |    |
   |                        |                        | +------+--------------------------------+    |
   |                        |                        | | 4    | Tests                          |    |
   |                        |                        | |      | ``model_interpolate()``        |    |
   |                        |                        | |      | for a single point.            |    |
   |                        |                        | +------+--------------------------------+    |
   |                        |                        | | 5    | Tests                          |    |
   |                        |                        | |      | ``model_interpolate()``        |    |
   |                        |                        | |      | for a range of                 |    |
   |                        |                        | |      | interpolation points.          |    |
   |                        |                        | +------+--------------------------------+    |
   |                        |                        | | 6    | Long, expensive test to        |    |
   |                        |                        | |      | return the metadata for        |    |
   |                        |                        | |      | every element of the           |    |
   |                        |                        | |      | state vector. May be           |    |
   |                        |                        | |      | useful to decide on            |    |
   |                        |                        | |      | known locations for            |    |
   |                        |                        | |      | subsequent testing.            |    |
   |                        |                        | +------+--------------------------------+    |
   |                        |                        | | 7    | Find the closest               |    |
   |                        |                        | |      | gridpoint to a known           |    |
   |                        |                        | |      | location.                      |    |
   |                        |                        | +------+--------------------------------+    |
   +------------------------+------------------------+----------------------------------------------+
   | run_tests(:)           | integer                | Specifies a list of tests to be              |
   |                        |                        | performed. Same test numbers as              |
   |                        |                        | described in test1thru. There are            |
   |                        |                        | some dependencies. Tests 4 and 5             |
   |                        |                        | require a valid model state - which          |
   |                        |                        | is read by Test 2. If a required test        |
   |                        |                        | is not specified, the required test          |
   |                        |                        | is enabled and run. A value of -1            |
   |                        |                        | means that ``test1thru`` will be             |
   |                        |                        | used.                                        |
   +------------------------+------------------------+----------------------------------------------+
   | verbose                | logical                | Print extra info about the                   |
   |                        |                        | ``model_mod_check`` run. This is only        |
   |                        |                        | used for more reporting during Test 5.       |
   |                        |                        | Be warned - it will generate                 |
   |                        |                        | several lines of output for each             |
   |                        |                        | point in the test!                           |
   +------------------------+------------------------+----------------------------------------------+

A more typical namelist for a single ensemble member for a model with an outer grid and a single nested grid is shown
below.

::

   &model_mod_check_nml
      input_state_files     = 'dart_vector1.nc','dart_vector2.nc'
      output_state_files    = 'check_me1.nc', 'check_me2.nc'
      all_metadata_file     = 'metadata.txt'
      verbose               = .TRUE.
      test1thru             = 5
      run_tests             = -1
      loc_of_interest       = 243.72386169, 52.78578186, 10.0
      x_ind                 = 12666739
      quantity_of_interest  = 'QTY_POTENTIAL_TEMPERATURE'
      interp_test_lonrange  = 144.0, 326.0
      interp_test_dlon      = 1.0
      interp_test_latrange  = -5.0, 80.0
      interp_test_dlat      = 1.0
      interp_test_vertrange = 100.0, 11000.0
      interp_test_dvert     = 200.0
      interp_test_vertcoord = 'VERTISHEIGHT'
     /


Files
-----

-  ``input.nml`` is used for ``model_mod_check_nml``
-  The ``"input_state_files"`` can either be a single file containing multiple restart files, or a single NetCDF restart
   file. One file per domain.
-  The ``"output_state_files"`` is the output netCDF files from Test 2. Check the attributes, values, etc.
-  ``check_me_interptest.nc`` and ``check_me_interptest.m`` are the result of Test 5.
-  ``"all_metadata_file"`` is the run-time output of Test 6.

Usage
-----

Unlike other DART components, you are expected to modify ``model_mod_check.f90`` to suit your needs as you develop your
``model_mod``. The code is roughly divided into the following categories:

#. Check the geometry information,
#. Read/write a restart file,
#. Check the construction of the state vector ... i.e. the metadata,
#. Interpolate at a single point,
#. Interpolate for a range of points.

Test 0. mandatory
~~~~~~~~~~~~~~~~~

The first test in ``model_mod_check`` reads the namelist and runs ``static_init_model`` - which generally sets the
geometry of the grid, the number of state variables and their shape, etc. Virtually everything requires knowledge of the
grid and state vector, so this block cannot be skipped.

Test 1. checking the geometry information
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The first test in ``model_mod_check`` exercises a basic required interface ``get_model_size()``. This also generates a
report on the geometry of the grid, the number of state variables and their shape, etc. as well as the total number of
elements in the DART state vector.

Test 2. read/writing a restart file
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

This directly reads and write state variables from the model netCDF file. This is a nice sanity check to make sure that
the DART state vector is being read in properly.

Test 3. check the construction of the state vector
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

It is critical to return the correct metadata for any given index into the DART state vector. This code block tests the
two most common features of the metadata. As a bonus, this routine is also quite useful to determine EXACTLY where to
place your first test observation. If you test precisely at a grid location, you should be able to really get a handle
on debugging your ``model_interpolate()`` routine.

Test 4. test interpolation on a single point
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

This tests your model's interpolation routine on a single point and returns the interpolated value. This requires that
Test 2 works - it needs a valid model state with data. Test 2 is automatically run if this test is selected.

Test 5. test interpolation on a range of values
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

This tests your model's interpolation routine on a range of values returns the interpolated grid in
``check_me_interptest.nc`` and ``check_me_interptest.m`` which can be read in Matlab and used to visualize the result.
This requires that Test 2 works - it needs a valid model state with data. Test 2 is automatically run if this test is
selected.

Test 6. exhaustively test the construction of the state vector
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

This can be a long test, depending on the size of your state vector. This returns the same data as in Test 3 - but *for
every element* in the state vector. The metadata are written to a file specified by ``all_metadata_file`` and
``check_me_interptest.m`` which can be read in Matlab and used to visualize the result.

Test 7. find the closest gridpoint to a test location
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

This is a good test to verify that *get_state_meta_data()* and the grid information are correct. Typically, one would
put in a location that is actually **on** the grid and see if the correct gridpoint index is returned. Repeat the test
with slightly different locations until the next gridpoint is closer. Repeat ...

References
----------

-  none
