BATS
====

Overview:
---------
BATS stands for the `Bermuda Atlantic Time-series Study <https://bats.bios.asu.edu/>`_. 
BATS has collected data on the physical, biological, and chemical properties of
the ocean every month since 1988. BATS was established to uncover mysteries of the 
deep ocean by analyzing important hydrographic and biological parameters 
throughout the water column.

Data Source: 
------------
Water column data from BATS (roughly 31N, 64W) can be downloaded from https://bats.bios.asu.edu/bats-data/
The data is stored in in different forms. This converter operates on the ASCII formatted file, 
often named ``bats_bottle.txt``

The data is huge extending from Oct 1988 to present day. It consists of a long list of information
such as 
`Depth, Oxygen, CO2, Nitrate, Phosphate, Silicate, Alkalinity, Organic Carbon, Bacteria, ...`

Observation Converter:
----------------------
The obs converter is a program called ``bats_to_obs`` and has a namelist by the name ``&bats_to_obs_nml`` 
Namelists start with an ampersand '&' and terminate with a slash '/'.

  .. code-block:: fortran 
  
        &bats_to_obs_nml
         text_input_file       = "../bats_bottle.txt"
         max_lines             = 68000
         read_starting_at_line = 61
         date_firstcol         = 14
         hourminute_firstcol   = 35
         lat_cols              = 42, 47
         lon_cols              = 51, 56
         vert_cols             = 64, 69
         scalar_obs_cols       = 113, 119,
                                 137, 143,
                                 145, 151,
                                 153, 159,
                                 170, 176,
                                 178, 184
         obs_uncertainties     = 0.2,
                                 0.2,
                                 0.2,
                                 0.2,
                                 0.2,
                                 0.2
         obs_out_dir           = '../obs_seq_files',
         debug                 = .true.
        /

This namelist provides control over the kind of observations to extract from the file in addition to their uncertainties. 
In its current form, the observations that are extracted from the data file are: 

``BATS_OXYGEN``, ``BATS_INORGANIC_CARBON``, ``BATS_ALKALINITY``, ``BATS_NITRATE``, ``BATS_PHOSPHATE``, ``BATS_SILICATE``

+-------------------------------------+--------------------+------------------------------------------------------------+
| Contents                            | Type               | Description                                                |
+=====================================+====================+============================================================+
| ``text_input_file``                 | character(len=256) | Pathname to the data file: ``bats_bottle.txt``             |
+-------------------------------------+--------------------+------------------------------------------------------------+
| ``max_lines``                       | integer            | Upper bound on the number of lines in the file that record | 
|                                     |                    | observations.                                              |
+-------------------------------------+--------------------+------------------------------------------------------------+
| ``read_starting_at_line``           | integer            | Skip the information in the header of the file.            |
+-------------------------------------+--------------------+------------------------------------------------------------+ 
| ``date_firstcol``                   | integer            | First column of the YYYYMMDD date code at each line.       |
+-------------------------------------+--------------------+------------------------------------------------------------+ 
| ``hourminute_firstcol``             | integer            | First column of the HHMM time stamp at each line.          |
+-------------------------------------+--------------------+------------------------------------------------------------+ 
| ``lat_cols``                        | integer(2)         | First and last columns where latitude is recorded.         | 
+-------------------------------------+--------------------+------------------------------------------------------------+
| ``lon_cols``                        | integer(2)         | First and last columns where longitude is recorded.        | 
+-------------------------------------+--------------------+------------------------------------------------------------+
| ``vert_cols``                       | integer(2)         | First and last columns where depth is recorded.            | 
+-------------------------------------+--------------------+------------------------------------------------------------+
| ``scalar_obs_cols``                 | integer(:, 2)      | ith row of this table should list the first and last       |
|                                     |                    | columns where the value of the ith observation variable    |
|                                     |                    | is recorded. Ordering of observation variables is defined  |
|                                     |                    | by the OTYPE_ORDERING parameter in bats_to_obs.f90.        |
+-------------------------------------+--------------------+------------------------------------------------------------+
| ``obs_uncertainties``               | real(:)            | ith entry of this list gives the uncertainty associated    |
|                                     |                    | with the ith observation variable.                         |
|                                     |                    |                                                            |
|                                     |                    | The observation error variance is defind as the square of  |
|                                     |                    | the product of the ``obs_uncertainties`` and the           |
|                                     |                    | observation value.                                         |
+-------------------------------------+--------------------+------------------------------------------------------------+
| ``obs_out_dir``                     | character(len=256) | Pathname to obs_seq files resulting from the converter.    |
+-------------------------------------+--------------------+------------------------------------------------------------+
| ``debug``                           | logical            | A switch that makes the converetr prints useful            | 
|                                     |                    | information as it runs.                                    |
+-------------------------------------+--------------------+------------------------------------------------------------+

Climatology:
~~~~~~~~~~~~
On top of assimilating real-time data, we often observe the quasi-cyclostationary behavior of the biogeochemical system 
over the period of one year, and we update MARBL parameters by comparing this observed climatology to a climatology 
predicted by MARBL. This usually involves running different forms of the ensmeble smoother where the model is re-run 
using the updated parameters over long periods of times. 

To access the observed climatology at BATS, the script ``bats_climatology.py`` can be used to generate the climatology
by averaging the data over time. The program ``bats_to_clim_obs`` can then be executed to generate DART-style 
observation sequence files using the climatological data. This code also supports 
`Multiple Data Assimilation (MDA) <https://www.sciencedirect.com/science/article/abs/pii/S0098300412000994>`_ in 
which the observations are assimilated multiple times with inflated observation error variance. 
