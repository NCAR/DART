PROGRAM ``littler_tf_dart``
===========================

Overview
--------

Programs to convert littler data files into DART observation sequence files, and vice versa. The capability of the
program is limited to wind and temperature from radiosondes.

The littler data files do not contain observation errors. The observation errors are in a separate file called
``obserr.txt``. The littler file generated here has to be preprocessed by the program ``3dvar_obs.exe`` before beeing
ingested in the WRF 3D-Var system.

Modules used
------------

::

   types_mod
   obs_sequence_mod
   obs_def_mod
   obs_kind_mod
   location/threed_sphere/location_mod
   time_manager_mod
   utilities_mod

Modules indirectly used
-----------------------

::

   assim_model_mod
   models/wrf/model_mod
   models/wrf/module_map_utils
   random_seq_mod

| 

Namelist
--------

The program does not have its own namelist. However, an ``input.nml`` file is required for the modules used by the
program.

Files
-----

-  input namelist ; ``input.nml``
-  Input - output observation files; ``obs_seq.out`` and ``little-r.dat``
-  Input - output littler observation error files ; ``obserr.txt``

File formats
~~~~~~~~~~~~

If there are no observation error at a particular pressure level, the default value of -1 is written in ``obserr.txt``.

References
----------

-  `3DVAR GROUP PAGE <http://www.mmm.ucar.edu/wrf/WG4/>`__

Private components
------------------

| 

.. container:: routine

   *call set_str_date(timestring, dart_time)*
   ::

      type(time_type),   intent(in)  ::  dart_time 
      character(len=20), intent(out) ::  timestring 

.. container:: indent1

   Given a dart_time (seconds, days), returns date as bbbbbbyyyymmddhhmmss, where b is a blank space.

| 

.. container:: routine

   *call set_dart_time(tstring, dart_time)*
   ::

      character(len=20), intent(in)  ::  tstring 
      type(time_type),   intent(out) ::  dart_time 

.. container:: indent1

   Given a date as bbbbbbyyyymmddhhmmss, where b is a blank space, returns the dart_time (seconds, days).

| 

.. container:: routine

   *call StoreObsErr(obs_err_var, pres, plevel, nlev, obs_err_std)*
   ::

      integer,  intent(in)    ::  nlev, pres 
      real(r8), intent(in)    ::  obs_err_var 
      integer,  intent(in)    ::  plevel(nlev) 
      real(r8), intent(inout) ::  obs_err_std(nlev) 

.. container:: indent1

   If the incoming pres corresponds exactly to a pressure level in plevel, then transfers the incoming obs_err_var into
   the array obs_err_std at the corresponding level.

| 

.. container:: routine

   *level_index = GetClosestLevel(ilev, vlev, nlev)*
   ::

      integer,  intent(in) ::  nlev, ilev 
      integer,  intent(in) ::  vlev(nlev) 

.. container:: indent1

   Returns the index of the closest level in vlev to the incoming ilev.

| 

.. container:: routine

   *call READ_OBSERR(filein, platform, sensor_name, err, nlevels)*
   ::

      CHARACTER (LEN=80), intent(in)  ::  filein 
      CHARACTER (LEN=80), intent(in)  ::  platform 
      CHARACTER (LEN=80), intent(in   ::  sensor_name 
      INTEGER,            intent(in)  ::  nlevels 
      REAL(r8),           intent(out) ::  err(nlevels) 

.. container:: indent1

   Read observational error on pressure levels (in hPa) from the incoming filein and store the result in the array err.
   It is assumed that filein has the same format as WRF 3D-Var ``obserr.txt`` file. It reads observational error for a
   specific platform (e.g. RAOBS) and a specific sensor (e.g. WIND SENSOR ERRORS).

| 

.. container:: routine

   *f_obstype = obstype(line)*
   ::

      CHARACTER (LEN= 80), intent(in) ::  line 

.. container:: indent1

   Read in a line the string present after keyword 'BOGUS', which should be the sensor name.

| 

.. container:: routine

   *f_sensor = sensor(line)*
   ::

      CHARACTER (LEN= 80), intent(in) ::  line 

.. container:: indent1

   Read in a line the string present after numbers, which should be the platform name.

| 

.. container:: routine

   *val = intplin(x,xx,yy)*
   ::

      INTEGER,  DIMENSION (:), intent(in) ::  xx 
      REAL(r8), DIMENSION (:), intent(in) ::  yy 
      REAL(r8),                intent(in) ::  x 

.. container:: indent1

   Do a linear interpolation.

| 

.. container:: routine

   *val = intplog(x,xx,yy)*
   ::

      INTEGER,  DIMENSION (:), intent(in) ::  xx 
      REAL(r8), DIMENSION (:), intent(in) ::  yy 
      REAL(r8),                intent(in) ::  x 

.. container:: indent1

   Do a log-linear interpolation.

| 

.. container:: routine

   *index = locate(x,xx)*
   ::

      INTEGER, DIMENSION (:), intent(in) ::  xx 
      REAL(r8),               intent(in) ::  x 

.. container:: indent1

   Return the index in xx such that xx(index) < x < xx(index+1).

| 
