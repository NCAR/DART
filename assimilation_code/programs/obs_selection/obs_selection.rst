program ``obs_selection``
=========================

Overview
--------

This specialized tool selects a subset of input observations from an observation sequence file. For a more general
purpose observation sequence file tool, see the :doc:`../obs_sequence_tool/obs_sequence_tool`. This tool takes a
selected list of observation types, times, and locations, and extracts only the matching observations out of one or more
obs_sequence files. The tool which creates the input selection file is usually
:doc:`../obs_seq_coverage/obs_seq_coverage`. Alternatively, the selection file can be a full observation sequence file,
in which case the types, times, and locations of those observations are used as the selection criteria.

This tool processes each observation sequence file listed in the input namelist ``filename_seq`` or
``filename_seq_list``. If the observation type, time and location matches an entry in the selection file, it is copied
through to the output. Otherwise it is ignored.

The actions of the ``obs_selection`` program are controlled by a Fortran namelist, read from a file named ``input.nml``
in the current directory. A detailed description of each namelist item is described in the namelist section of this
document. The names used in this discussion refer to these namelist items.

Namelist
--------

This namelist is read from the file ``input.nml``. Namelists start with an ampersand '&' and terminate with a slash '/'.
Character strings that contain a '/' must be enclosed in quotes to prevent them from prematurely terminating the
namelist.

::

   &obs_selection_nml
      filename_seq          = ''
      filename_seq_list     = ''
      filename_out          = 'obs_seq.processed'
      num_input_files       = 0
      selections_file       = 'obsdef_mask.txt'
      selections_is_obs_seq = .false.
      latlon_tolerance      = 0.000001
      match_vertical        = .false.
      surface_tolerance     = 0.0001
      pressure_tolerance    = 0.001
      height_tolerance      = 0.0001
      scaleheight_tolerance = 0.001
      level_tolerance       = 0.00001
      print_only            = .false.
      partial_write         = .false.
      print_timestamps      = .false.
      calendar              = 'Gregorian'
     /

| 

.. container::

   +-----------------------+------------------------------------+-------------------------------------------------------+
   | Item                  | Type                               | Description                                           |
   +=======================+====================================+=======================================================+
   | filename_seq          | character(len=256), dimension(500) | The array of names of the observation sequence files  |
   |                       |                                    | to process, up to a max count of 500 files. (Specify  |
   |                       |                                    | only the actual number of input files. It is not      |
   |                       |                                    | necessary to specify 500 entries.)                    |
   +-----------------------+------------------------------------+-------------------------------------------------------+
   | filename_seq_list     | character(len=256)                 | An alternative way to specify the list of input       |
   |                       |                                    | files. The name of a text file which contains, one    |
   |                       |                                    | per line, the names of the observation sequence files |
   |                       |                                    | to process. You can only specify one of filename_seq  |
   |                       |                                    | OR filename_seq_list, not both.                       |
   +-----------------------+------------------------------------+-------------------------------------------------------+
   | num_input_files       | integer                            | Optional. The number of observation sequence files to |
   |                       |                                    | process. Maximum of 500. If 0, the length is set by   |
   |                       |                                    | the number of input files given. If non-zero, must    |
   |                       |                                    | match the given input file list length. (Can be used  |
   |                       |                                    | to verify the right number of input files were        |
   |                       |                                    | processed.)                                           |
   +-----------------------+------------------------------------+-------------------------------------------------------+
   | filename_out          | character(len=256)                 | The name of the resulting output observation sequence |
   |                       |                                    | file. There is only a single output file from this    |
   |                       |                                    | tool. If the input specifies multiple obs_seq input   |
   |                       |                                    | files, the results are concatinated into a single     |
   |                       |                                    | output file.                                          |
   +-----------------------+------------------------------------+-------------------------------------------------------+
   | selections_file       | character(len=256)                 | The name of the input file containing the mask of     |
   |                       |                                    | observation definitions (the textfile output of       |
   |                       |                                    | :doc:`../obs_seq_coverage/obs_seq_coverage`).         |
   |                       |                                    | Alternatively, this can be the name of a full         |
   |                       |                                    | observation sequence file. In this case, the types,   |
   |                       |                                    | times, and locations are extracted from this file and |
   |                       |                                    | then used in the same manner as a mask file from the  |
   |                       |                                    | coverage tool.                                        |
   +-----------------------+------------------------------------+-------------------------------------------------------+
   | selections_is_obs_seq | logical                            | If .TRUE. the filename given for the                  |
   |                       |                                    | "selections_file" is a full obs_sequence file and not |
   |                       |                                    | a text file from the coverage tool.                   |
   +-----------------------+------------------------------------+-------------------------------------------------------+
   | latlon_tolerance      | real(r8)                           | Specified in degrees. For observations to match in    |
   |                       |                                    | the horizontal the difference in degrees for each of  |
   |                       |                                    | latitude and longitude must be less than this         |
   |                       |                                    | threshold. If less than or equal to 0, the values     |
   |                       |                                    | must match exactly.                                   |
   +-----------------------+------------------------------------+-------------------------------------------------------+
   | match_vertical        | logical                            | If .TRUE. the locations of the observations in the    |
   |                       |                                    | input files have to match the selection list not only |
   |                       |                                    | the horizontal but also in the vertical.              |
   +-----------------------+------------------------------------+-------------------------------------------------------+
   | surface_tolerance     | real(r8)                           | Specified in meters. If "match_vertical" is .FALSE.   |
   |                       |                                    | this value is ignored. If "match_vertical" is .TRUE., |
   |                       |                                    | this applies to observations with a vertical type of  |
   |                       |                                    | VERTISSURFACE. For observations which match in the    |
   |                       |                                    | horizontal, the vertical surface elevation difference |
   |                       |                                    | must be less than this to be considered the same.     |
   +-----------------------+------------------------------------+-------------------------------------------------------+
   | pressure_tolerance    | real(r8)                           | Specified in pascals. If "match_vertical" is .FALSE.  |
   |                       |                                    | this value is ignored. If "match_vertical" is .TRUE., |
   |                       |                                    | this applies to observations with a vertical type of  |
   |                       |                                    | VERTISPRESSURE. For observations which match in the   |
   |                       |                                    | horizontal, the vertical difference must be less than |
   |                       |                                    | this to be considered the same.                       |
   +-----------------------+------------------------------------+-------------------------------------------------------+
   | height_tolerance      | real(r8)                           | Specified in meters. If "match_vertical" is .FALSE.   |
   |                       |                                    | this value is ignored. If "match_vertical" is .TRUE., |
   |                       |                                    | this applies to observations with a vertical type of  |
   |                       |                                    | VERTISHEIGHT. For observations which match in the     |
   |                       |                                    | horizontal, the vertical difference must be less than |
   |                       |                                    | this to be considered the same.                       |
   +-----------------------+------------------------------------+-------------------------------------------------------+
   | scaleheight_tolerance | real(r8)                           | Specified in unitless values. If "match_vertical" is  |
   |                       |                                    | .FALSE. this value is ignored. If "match_vertical" is |
   |                       |                                    | .TRUE., this applies to observations with a vertical  |
   |                       |                                    | type of VERTISSCALEHEIGHT. For observations which     |
   |                       |                                    | match in the horizontal, the vertical difference must |
   |                       |                                    | be less than this to be considered the same.          |
   +-----------------------+------------------------------------+-------------------------------------------------------+
   | level_tolerance       | real(r8)                           | Specified in fractional model levels. If              |
   |                       |                                    | "match_vertical" is .FALSE. this value is ignored. If |
   |                       |                                    | "match_vertical" is .TRUE., this applies to           |
   |                       |                                    | observations with a vertical type of VERTISLEVEL. For |
   |                       |                                    | observations which match in the horizontal, the       |
   |                       |                                    | vertical difference must be less than this to be      |
   |                       |                                    | considered the same. Note that some models only       |
   |                       |                                    | support integer level values, but others support      |
   |                       |                                    | fractional levels. The vertical value in an           |
   |                       |                                    | observation is a floating point/real value, so        |
   |                       |                                    | fractional levels are possible to specify for an      |
   |                       |                                    | observation.                                          |
   +-----------------------+------------------------------------+-------------------------------------------------------+
   | print_only            | logical                            | If .TRUE. do not create an output file, but print a   |
   |                       |                                    | summary of the number and types of each observation   |
   |                       |                                    | in each input file, and then the number of            |
   |                       |                                    | observations and types which would have been created  |
   |                       |                                    | in an output file.                                    |
   +-----------------------+------------------------------------+-------------------------------------------------------+
   | partial_write         | logical                            | Generally only used for debugging problems. After     |
   |                       |                                    | each input obs_seq file is processed, this flag, if   |
   |                       |                                    | .TRUE., causes the code to write out the partial      |
   |                       |                                    | results to the output file. The default is to process |
   |                       |                                    | all input files (if more than a single file is        |
   |                       |                                    | specified) and write the output file only at the end  |
   |                       |                                    | of the processing.                                    |
   +-----------------------+------------------------------------+-------------------------------------------------------+
   | print_timestamps      | logical                            | Generally only used for debugging very slow execution |
   |                       |                                    | runs. This flag, if .TRUE., causes the code to output |
   |                       |                                    | timestamps (wall clock time) at various locations     |
   |                       |                                    | during the processing phases. It may help isolate     |
   |                       |                                    | where particularly slow execution times are           |
   |                       |                                    | occurring. For very large input files, or long lists  |
   |                       |                                    | of input files, it can also help to estimate what the |
   |                       |                                    | eventual run time of the job will be.                 |
   +-----------------------+------------------------------------+-------------------------------------------------------+
   | calendar              | character(len=32)                  | Set to the name of the calendar; only controls the    |
   |                       |                                    | printed output for the dates of the first and last    |
   |                       |                                    | observations in the file. Set this to "no_calendar"   |
   |                       |                                    | if the observations are not using any calendar.       |
   +-----------------------+------------------------------------+-------------------------------------------------------+

| 

Building
--------

Most ``$DART/models/*/work`` directories contain files needed to build this tool along with the other executable
programs. It is also possible to build this tool in the ``$DART/observations/utilities`` directory. In either case the
``preprocess`` program must be built and run first to define what set of observation types will be supported. See the
:doc:`../../../assimilation_code/programs/preprocess/preprocess` for more details on how to define the list and run it.
The ``&preprocess_nml`` namelist in the ``input.nml`` file must contain files with definitions for the combined set of
all observation types which will be encountered over all input obs_seq files.

Usually the directories where executables are built will include a ``quickbuild.sh`` script which builds and runs
preprocess and then builds the rest of the executables.

Modules used
------------

::

   types_mod
   utilities_mod
   time_manager_mod
   obs_def_mod
   obs_sequence_mod

Files
-----

-  ``input.nml``
-  The input files specified in the ``filename_seq`` namelist variable.
-  The output file specified in the ``filename_out`` namelist variable.

References
----------

-  none
