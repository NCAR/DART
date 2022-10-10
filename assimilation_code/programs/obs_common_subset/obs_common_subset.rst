program ``obs_common_subset``
=============================

Overview
--------

This specialized tool allows you to select subsets of observations from two or more observation sequence files output
from ``filter``. It creates a new set of output observation sequence files containing only the observations which were
successfully assimilated in all experiments.

Experiments using the same input observation sequence file but with different configurations (e.g. different inflation
values, different localization radii, etc) can assimilate different numbers of the available observations. In that case
there will be differences in the diagnostic plots which are not directly relatable to the differences in the quality of
the assimilation. If this tool is run on the ``obs_seq.final`` files from all the experiments and then the diagnostics
are generated, only the observations which were assimilated in all experiments will contribute to the summary
statistics. A more direct comparison can be made and improvements can be correctly attributed to the differences in the
experimental parameters.

This tool is intended to be used when comparing the results from a group of related experiments in which **the exact
same input observation sequence file** is used for all runs. The tool cannot process observation sequence files which
differ in anything other than whether an observation was successfully assimilated/evaluated or not. Note that it is fine
to add or remove observation types from the ``assimilate_these_obs_types`` or ``evaluate_these_obs_types`` namelist
items for different experiments. The output observation sequence files will still contain an identical list of
observations, with some marked with a DART QC indicating 'not assimilated because of namelist control'.

See the :doc:`"two_experiment" diagnostic plots in <../../../guide/matlab-observation-space>`
documentation for Matlab scripts
supplied with DART to directly compare the observation diagnostic output from multiple experiments
(it does more than two, the script has a poor name).

This is one of a set of tools which operate on observation sequence files. For a more general purpose tool see the
:doc:`../obs_sequence_tool/obs_sequence_tool`, and for a more flexible selection tool see the obs_selection_tool.

Creating an input filelist
^^^^^^^^^^^^^^^^^^^^^^^^^^

One of the inputs to this tool is a list of filenames to compare. The filenames can be directly in the namelist file, or
they can be in a set of separate text files. The latter may be easier when there are more than just a few files to
compare.

For experiments where there are multiple job steps, and so multiple output observation sequence files per experiment,
the input to this tool would then be a list of lists of filenames. Each set of names must be put into a text file with
each filename on a separate line.

If each experiment was run in a different set of directories, and if a list of observation sequence filenames was made
with the ``ls`` command:

::

   > ls exp1/*/obs_seq.final > exp1list
   > cat exp1list
   exp1/step1/obs_seq.final
   exp1/step2/obs_seq.final
   exp1/step3/obs_seq.final
   exp1/step4/obs_seq.final
   > ls exp2/*/obs_seq.final > exp2list
   > cat exp2list
   exp2/step1/obs_seq.final
   exp2/step2/obs_seq.final
   exp2/step3/obs_seq.final
   exp2/step4/obs_seq.final
   > ls exp3/*/obs_seq.final > exp3list
   > cat exp2list
   exp3/step1/obs_seq.final
   exp3/step2/obs_seq.final
   exp3/step3/obs_seq.final
   exp3/step4/obs_seq.final

Then the namelist entries would be:

::

    filename_seq = ''
    filename_seq_list = 'exp1list', 'exp2list', exp3list'
    num_to_compare_at_once = 3

Namelist
--------

This namelist is read from the file ``input.nml``. Namelists start with an ampersand '&' and terminate with a slash '/'.
Character strings that contain a '/' must be enclosed in quotes to prevent them from prematurely terminating the
namelist.

::

   &obs_common_subset_nml
    num_to_compare_at_once = 2,
    filename_seq           = '',
    filename_seq_list      = '',
    filename_out_suffix    = '.common' ,
    print_every            = 10000,
    dart_qc_threshold      = 3,
    calendar               = 'Gregorian',
    print_only             = .false.,
    eval_and_assim_can_match = .false.,
   /

| 

.. container::

   +--------------------------+-------------------------------------+---------------------------------------------------+
   | Item                     | Type                                | Description                                       |
   +==========================+=====================================+===================================================+
   | num_to_compare_at_once   | integer                             | Number of observation sequence files to compare   |
   |                          |                                     | together at a time. Most commonly the value is 2, |
   |                          |                                     | but can be any number. If more than this number   |
   |                          |                                     | of files are listed as inputs, the tool will loop |
   |                          |                                     | over the list N files at a time.                  |
   +--------------------------+-------------------------------------+---------------------------------------------------+
   | filename_seq             | character(len=256), dimension(5000) | The array of names of the observation sequence    |
   |                          |                                     | files to process. If more than N files (where N   |
   |                          |                                     | is num_to_compare_at_once) are listed, they       |
   |                          |                                     | should be ordered so the first N files are        |
   |                          |                                     | compared together, followed by the next set of N  |
   |                          |                                     | files, etc. You can only specify one of           |
   |                          |                                     | filename_seq OR filename_seq_list, not both.      |
   +--------------------------+-------------------------------------+---------------------------------------------------+
   | filename_seq_list        | character(len=256), dimension(100)  | An alternative way to specify the list of input   |
   |                          |                                     | observation sequence files. Give a list of N      |
   |                          |                                     | filenames which contain, one per line, the names  |
   |                          |                                     | of the observation sequence files to process.     |
   |                          |                                     | There should be N files specified (where N is     |
   |                          |                                     | num_to_compare_at_once), and the first            |
   |                          |                                     | observation sequence filename listed in each file |
   |                          |                                     | will be compared together, then the second, until |
   |                          |                                     | the lists are exhausted. You can only specify one |
   |                          |                                     | of filename_seq OR filename_seq_list, not both.   |
   +--------------------------+-------------------------------------+---------------------------------------------------+
   | filename_out_suffix      | character(len=32)                   | A string to be appended to each of the input      |
   |                          |                                     | observation sequence file names to create the     |
   |                          |                                     | output filenames.                                 |
   +--------------------------+-------------------------------------+---------------------------------------------------+
   | print_every              | integer                             | To indicate progress, a count of the successfully |
   |                          |                                     | processed observations is printed every Nth set   |
   |                          |                                     | of obs. To decrease the output volume set this to |
   |                          |                                     | a larger number. To disable this output           |
   |                          |                                     | completely set this to -1.                        |
   +--------------------------+-------------------------------------+---------------------------------------------------+
   | dart_qc_threshold        | integer                             | Observations with a DART QC value larger than     |
   |                          |                                     | this threshold will be discarded. Note that this  |
   |                          |                                     | is the QC value set by ``filter`` to indicate the |
   |                          |                                     | outcome of trying to assimilate an observation.   |
   |                          |                                     | This is not related to the incoming data QC. For  |
   |                          |                                     | an observation which was successfully assimilated |
   |                          |                                     | or evaluated in both the Prior and Posterior this |
   |                          |                                     | should be set to 1. To also include observations  |
   |                          |                                     | which were successfully processed in the Prior    |
   |                          |                                     | but not the Posterior, set to 3. To ignore the    |
   |                          |                                     | magnitude of the DART QC values and keep          |
   |                          |                                     | observations only if the DART QCs match, set this |
   |                          |                                     | to any value higher than 7.                       |
   +--------------------------+-------------------------------------+---------------------------------------------------+
   | calendar                 | character(len=32)                   | Set to the name of the calendar; only controls    |
   |                          |                                     | the printed output for the dates of the first and |
   |                          |                                     | last observations in the file. Set this to        |
   |                          |                                     | "no_calendar" if the observations are not using   |
   |                          |                                     | any calendar.                                     |
   +--------------------------+-------------------------------------+---------------------------------------------------+
   | print_only               | logical                             | If .TRUE. do not create the output files, but     |
   |                          |                                     | print a summary of the number and types of each   |
   |                          |                                     | observation in each of the input and output       |
   |                          |                                     | files.                                            |
   +--------------------------+-------------------------------------+---------------------------------------------------+
   | eval_and_assim_can_match | logical                             | Normally .FALSE. . If .TRUE. then observations    |
   |                          |                                     | which were either successfully evaluated OR       |
   |                          |                                     | assimilated will match and are kept.              |
   +--------------------------+-------------------------------------+---------------------------------------------------+

| 

Building
--------

Most ``$DART/models/*/work`` directories will build the tool along with other executable programs. It is also possible
to build the tool in the ``$DART/observations/utilities`` directory. The ``preprocess`` program must be built and run
first, to define what set of observation types will be supported. See the
:doc:`../../../assimilation_code/programs/preprocess/preprocess` for more details on how to define the list and run it.
The combined list of all observation types which will be encountered over all input files must be in the preprocess
input list. The other important choice when building the tool is to include a compatible locations module. For the
low-order models, the ``oned`` module should be used; for real-world observations, the ``threed_sphere`` module should
be used.

Generally the directories where executables are built will include a "quickbuild.sh" script which will build and run
preprocess and then build the rest of the executables. The "input.nml" namelists will need to be edited to include all
the required observation types first.

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
-  The input files specified in the ``filename_seq`` or ``filename_seq_list`` namelist variable.
-  The output files are specified by appending the string from the ``filename_out_suffix`` namelist item to the input
   filenames.

References
----------

-  none
