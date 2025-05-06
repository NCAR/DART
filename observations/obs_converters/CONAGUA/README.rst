CONAGUA
=======

Overview: 
--------- 
CONAGUA stands for Comision Nacional del Agua (National Water Commision).
It is Mexico's central authority for water management.
Established in 1989, CONAGUA is responsible for management, preservation,
and distribution of water resources in Mexico.

Data Source:
------------
The streamflow observations from CONAGUA are naturally in a Microsoft
database format. These can be obtained from the
`CONAGUA Webpage <https://www.gob.mx/conagua>`_.
Mirce, a PhD student at the time from UNAM,  converts these
one-at-a-time to a "csv" format. The filenames have a gauge
identifier in them. There is also another file
that has the lat/lon of the gauge. The streamflow 
is in cubic meters per second (cms). Each of the column 
headers in the daily observation files are as follows:

``pk_ano = Year``, ``pk_mes = Month``, ``ngasto_d01, d02 ... = Streamflow in days``. 

Observation Converter:
----------------------
The obs converter is a program called ``CONAGUA_convert_streamflow`` and 
meta_data_filehas a namelist by the name ``&CONAGUA_convert_streamflow``.
Namelists start with an ampersand '&' and terminate with a slash '/'.

  .. code-block:: fortran 
  
        &CONAGUA_convert_streamflow_nml
         meta_data_file         = 'LaSierra/Observations/LaSierra_strm_pts.csv'
         data_file_list         = 'LaSierra/Observations/daily_files_list.txt'
         obs_out_file           = 'obs_seq.out'
         obs_fraction_for_error = 0.05
         obs_min_err_sd         = 0.5
         debug                  = .true. 
        /

This namelist provides control over the kind of observations (streamflow)
to extract from the file in addition to their uncertainties.

+----------------------------+--------------------+-------------------------------------------------------------+
| Contents                   | Type               | Description                                                 |
+============================+====================+=============================================================+
| ``meta_data_file``         | character(len=256) | Pathname to the data file                                   |
+----------------------------+--------------------+-------------------------------------------------------------+
| ``data_file_list``         | character(len=256) | List of data files if processing several days at once       |
+----------------------------+--------------------+-------------------------------------------------------------+
| ``obs_out_file``           | character(len=256) | Name of output DART-style file: ``obs_seq.out``             |
+----------------------------+--------------------+-------------------------------------------------------------+
| ``obs_fraction_for_error`` | real               | Factor to parametrize the streamflow observation error:     |
|                            |                    | obs error variance = obs value * ``obs_fraction_for_error`` |
+----------------------------+--------------------+-------------------------------------------------------------+
| ``obs_min_err_sd``         | character(len=256) | Lower bound for the observation error variance              |
+----------------------------+--------------------+-------------------------------------------------------------+
| ``debug``                  | logical            | A switch to make the converter more verbose                 |
+----------------------------+--------------------+-------------------------------------------------------------+ 
