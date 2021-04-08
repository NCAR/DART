MODULE utilities_mod
====================

Overview
--------

Provides a number of tools used by most DART modules including tools for file IO, diagnostic tools for registering
modules and recording namelist arguments, and an error handler.

Namelist
--------

This namelist is read from the file ``input.nml``. Namelists start with an ampersand '&' and terminate with a slash '/'.
Character strings that contain a '/' must be enclosed in quotes to prevent them from prematurely terminating the
namelist.

::

   &utilities_nml
      TERMLEVEL      = 2,
      logfilename    = 'dart_log.out',
      nmlfilename    = 'dart_log.nml',
      module_details = .true.,
      print_debug    = .false.,
      write_nml      = 'file'
   /

| 

The namelist controls how the logging, namelist, messages, and general utility routines behave.

.. container::

   +----------------+--------------------+------------------------------------------------------------------------------+
   | Item           | Type               | Description                                                                  |
   +================+====================+==============================================================================+
   | TERMLEVEL      | integer            | Level at which calls to error manager terminate program. The default setting |
   |                |                    | is warnings and errors terminate the program. Setting this to 2 (E_ERR)      |
   |                |                    | means only errors terminate. Setting this to 3 means even errors do not      |
   |                |                    | cause an exit (which is not a good idea).                                    |
   +----------------+--------------------+------------------------------------------------------------------------------+
   | logfilename    | character(len=256) | File to which the log messages are written.                                  |
   +----------------+--------------------+------------------------------------------------------------------------------+
   | nmlfilename    | character(len=256) | File to which the namelist output is written. Can be the same name as the    |
   |                |                    | logfile.                                                                     |
   +----------------+--------------------+------------------------------------------------------------------------------+
   | module_details | logical            | Each source code module can write out the repository version number and      |
   |                |                    | filename to the logfile. Verbose, but useful for knowing what version of the |
   |                |                    | code was used during the run.                                                |
   +----------------+--------------------+------------------------------------------------------------------------------+
   | print_debug    | logical            | Setting this to .true. causes additional debug messages to print. These can  |
   |                |                    | be very verbose and by default are turned off.                               |
   +----------------+--------------------+------------------------------------------------------------------------------+
   | write_nml      | character(len=32)  | String which controls where to write the namelist values that are being used |
   |                |                    | for this execution. Valid values are: 'none', 'file', 'terminal', 'both'.    |
   |                |                    | 'none' turns off this write. 'file' writes a copy only to the                |
   |                |                    | ``nmlfilename``. Writes are always in append mode, so the most recent        |
   |                |                    | information will be at the end of an existing file. 'terminal' will write to |
   |                |                    | the job's standard output. 'both' will write both to the nml file and the    |
   |                |                    | standard output unit.                                                        |
   +----------------+--------------------+------------------------------------------------------------------------------+

| 

Other modules used
------------------

::

   types_mod
   netCDF

Public interfaces
-----------------

======================= =====================
*use utilities, only :* file_exist
\                       get_unit
\                       open_file
\                       close_file
\                       timestamp
\                       register_module
\                       error_handler
\                       to_upper
\                       nc_check
\                       logfileunit
\                       nmlfileunit
\                       initialize_utilities
\                       finalize_utilities
\                       dump_unit_attributes
\                       find_namelist_in_file
\                       check_namelist_read
\                       find_textfile_dims
\                       file_to_text
\                       is_longitude_between
\                       get_next_filename
\                       set_filename_list
\                       set_tasknum
\                       set_output
\                       do_output
\                       E_DBG, DEBUG
\                       E_MSG, MESSAGE
\                       E_WARN, WARNING
\                       E_ERR, FATAL
======================= =====================

A note about documentation style. Optional arguments are enclosed in brackets *[like this]*.

| 

.. container:: routine

   *var = file_exist(file_name)*
   ::

      logical                      :: file_exist
      character(len=*), intent(in) :: file_name

.. container:: indent1

   Returns true if file_name exists in the working directory, else false.

   ============= ==============================================
   ``var``       True if file_name exists in working directory.
   ``file_name`` Name of file to look for.
   ============= ==============================================

| 

.. container:: routine

   *var = get_unit()*
   ::

      integer :: get_unit

.. container:: indent1

   Returns an unused unit number for IO.

   ======= ======================
   ``var`` An unused unit number.
   ======= ======================

| 

.. container:: routine

   *var = open_file(fname [, form, action])*
   ::

      integer                                :: open_file
      character(len=*), intent(in)           :: fname
      character(len=*), optional, intent(in) :: form
      character(len=*), optional, intent(in) :: action

.. container:: indent1

   Returns a unit number that is opened to the file fname. If form is not present or if form is "formatted" or
   "FORMATTED", file is opened for formatted IO. Otherwise, it is unformatted. The action string is the standard action
   string for Fortran IO (see F90 language description).

   ========= ======================================================================================================
   ``var``   Unit number opened to file fname.
   ``fname`` Name of file to be opened.
   *form*    Format: 'formatted' or 'FORMATTED' give formatted, anything else is unformatted. Default is formatted.
   *action*  Standard fortran string description of requested file open action.
   ========= ======================================================================================================

| 

.. container:: routine

   *call timestamp([string1, string2, string3,] pos)*
   ::

      character(len=*), optional, intent(in) :: string1
      character(len=*), optional, intent(in) :: string2
      character(len=*), optional, intent(in) :: string3
      character(len=*), intent(in)           :: pos

.. container:: indent1

   Prints the message 'Time is YYYY MM DD HH MM SS' to the logfile along with three optional message strings. If the pos
   argument is 'end', the message printed is 'Finished... at YYYY MM DD HH MM SS' and the logfile is closed.

   ========= ====================================
   *string1* An optional message to be printed.
   *string2* An optional message to be printed.
   *string3* An optional message to be printed.
   ``pos``   If 'end' terminates log_file output.
   ========= ====================================

| 

.. container:: routine

   *call close_file(iunit)*
   ::

      integer, intent(in) :: iunit

.. container:: indent1

   Closes the given unit number. If the unit is not open, nothing happens.

   ========= =======================
   ``iunit`` File unit to be closed.
   ========= =======================

| 

.. container:: routine

   *call register_module(src, rev, rdate)*
   ::

      character(len=*), intent(in) :: src
      character(len=*), optional, intent(in) :: rev
      character(len=*), optional, intent(in) :: rdate

.. container:: indent1

   Writes the source name to both the logfileunit and to standard out. The rev and revdate are deprecated as they are
   unsupported by git.

   ========= =================
   ``src``   source file name.
   ``rev``   ignored
   ``rdate`` ignored
   ========= =================

| 

.. container:: routine

   *call error_handler(level, routine, text, src, rev, rdate [, aut, text2, text3])*
   ::

      integer, intent(in)                    :: level
      character(len=*), intent(in)           :: routine
      character(len=*), intent(in)           :: text
      character(len=*), intent(in)           :: src
      character(len=*), intent(in)           :: rev
      character(len=*), intent(in)           :: rdate
      character(len=*), optional, intent(in) :: aut
      character(len=*), optional, intent(in) :: text2
      character(len=*), optional, intent(in) :: text3

.. container:: indent1

   Prints an error message to standard out and to the logfileunit. The message contains the routine name, an error
   message, the source file, revision and revision date, and optionally the author. The level of severity is message,
   debug, warning, or error. If the level is greater than or equal to the TERMLEVEL (set in the namelist), execution is
   terminated. The default TERMLEVEL only stops for ERRORS.

   =========== ===============================================================================
   ``level``   Error severity (message, debug, warning, error). See below for specific ations.
   ``routine`` Name of routine generating error.
   ``text``    Error message.
   ``src``     Source file containing routine generating message.
   ``rev``     Revision number of source file.
   ``rdate``   Revision date of source file.
   *aut*       Author of routine.
   *text2*     If specified, the second line of text for the error message.
   *text3*     If specified, the third line of text for the error message.
   =========== ===============================================================================

| 

.. container:: routine

   *call find_namelist_in_file(namelist_file_name, nml_name, iunit, [,write_to_logfile_in])*
   ::

      character(len=*),  intent(in)          :: namelist_file_name
      character(len=*),  intent(in)          :: nml_name
      integer,           intent(out)         :: iunit
      logical, optional, intent(in)          :: write_to_logfile_in

.. container:: indent1

   Opens the file namelist_file_name if it exists on unit iunit. A fatal error occurs if the file does not exist (DART
   requires an input.nml to be available, even if it contains no values). Searches through the file for a line
   containing ONLY the string &nml_name (for instance &filter_nml if nml_name is "filter_nml"). If this line is found,
   the file is rewound and the routine returns. Otherwise, a fatal error message is issued.

   +-----------------------+---------------------------------------------------------------------------------------------+
   | ``namelist``          | Name of file assumed to hold the namelist.                                                  |
   +-----------------------+---------------------------------------------------------------------------------------------+
   | ``nml_name``          | Name of the namelist to be searched for in the file, for instance, filter_nml.              |
   +-----------------------+---------------------------------------------------------------------------------------------+
   | ``iunit``             | Channel number on which file is opened.                                                     |
   +-----------------------+---------------------------------------------------------------------------------------------+
   | *write_to_logfile_in* | When the namelist for the utilities module is read, the logfile has not yet been open       |
   |                       | because its name is in the namelist. If errors are found, have to write to standard out.    |
   |                       | So, when utilities module calls this internally, this optional argument is set to false.    |
   |                       | For all other applications, it is normally not used (default is false).                     |
   +-----------------------+---------------------------------------------------------------------------------------------+

| 

.. container:: routine

   *call check_namelist_read(iunit, iostat_in, nml_name, [, write_to_logfile_in])*
   ::

      integer, intent(in)                    :: iunit
      integer, intent(in)                    :: iostat_in
      character(len=*), intent(in)           :: nml_name
      logical, optional, intent(in)          :: write_to_logfile_in

.. container:: indent1

   Once a namelist has been read from an opened namelist file, this routine checks for possible errors in the read. If
   the namelist read was successful, the file opened on iunit is closed and the routine returns. If iostat is not zero,
   an attempt is made to rewind the file on iunit and read the last line that was successfully read. If this can be
   done, this last line is printed with the preamble "INVALID NAMELIST ENTRY". If the attempt to read the line after
   rewinding fails, it is assumed that the original read (before the call to this subroutine) failed by reaching the end
   of the file. An error message stating that the namelist started but was never terminated is issued.

   +-----------------------+---------------------------------------------------------------------------------------------+
   | ``iunit``             | Channel number on which file is opened.                                                     |
   +-----------------------+---------------------------------------------------------------------------------------------+
   | ``iostat_in``         | Error status return from an attempted read of a namelist from this file.                    |
   +-----------------------+---------------------------------------------------------------------------------------------+
   | ``nml_name``          | The name of the namelist that is being read (for instance filter_nml).                      |
   +-----------------------+---------------------------------------------------------------------------------------------+
   | *write_to_logfile_in* | When the namelist for the utilities module is read, the logfile has not yet been open       |
   |                       | because its name is in the namelist. If errors are found, have to write to standard out.    |
   |                       | So, when utilities module calls this internally, this optional argument is set to false.    |
   |                       | For all other applications, it is normally not used (default is false).                     |
   +-----------------------+---------------------------------------------------------------------------------------------+

| 

.. container:: routine

   *call find_textfile_dims (fname, nlines, linelen)*
   ::

      character(len=*), intent (IN)  :: fname
      integer,          intent (OUT) :: nlines
      integer,          intent (OUT) :: linelen

.. container:: indent1

   Determines the number of lines and maximum line length of an ASCII text file.

   =========== ==========================================
   ``fname``   input, character string file name
   ``nlines``  output, number of lines in the file
   ``linelen`` output, length of longest line in the file
   =========== ==========================================

| 

.. container:: routine

   *call file_to_text (fname, textblock)*
   ::

      character(len=*),               intent (IN)  :: fname
      character(len=*), dimension(:), intent (OUT) :: textblock

.. container:: indent1

   Opens the given filename and reads ASCII text lines into a character array.

   ============= ===========================================
   ``fname``     input, character string file name
   ``textblock`` output, character array of text in the file
   ============= ===========================================

| 

.. container:: routine

   *var = is_longitude_between(lon, minlon, maxlon [, doradians])*
   ::

      real(r8), intent(in)           :: lon
      real(r8), intent(in)           :: minlon
      real(r8), intent(in)           :: maxlon
      logical,  intent(in), optional :: doradians
      logical                        :: is_longitude_between

.. container:: indent1

   Uniform way to test longitude ranges, in degrees, on a globe. Returns true if lon is between min and max, starting at
   min and going EAST until reaching max. Wraps across 0 longitude. If min equals max, all points are inside. Includes
   endpoints. If optional arg doradians is true, do computation in radians between 0 and 2*PI instead of default 360.
   There is no rejection of input values based on range; they are all converted to a known range by calling modulo()
   first.

   +-------------+-------------------------------------------------------------------------------------------------------+
   | ``var``     | True if lon is between min and max.                                                                   |
   +-------------+-------------------------------------------------------------------------------------------------------+
   | ``lon``     | Location to test.                                                                                     |
   +-------------+-------------------------------------------------------------------------------------------------------+
   | ``minlon``  | Minimum longitude. Region will start here and go east.                                                |
   +-------------+-------------------------------------------------------------------------------------------------------+
   | ``maxlon``  | Maximum longitude. Region will end here.                                                              |
   +-------------+-------------------------------------------------------------------------------------------------------+
   | *doradians* | Optional argument. Default computations are in degrees. If this argument is specified and is .true.,  |
   |             | do the computation in radians, and wrap across the globe at 2 \* PI. All inputs must then be          |
   |             | specified in radians.                                                                                 |
   +-------------+-------------------------------------------------------------------------------------------------------+

| 

.. container:: routine

   *var = get_next_filename( listname, lineindex )*
   ::

      character(len=*),  intent(in) :: listname
      integer,           intent(in) :: lineindex
      character(len=128)            :: get_next_filename

.. container:: indent1

   Returns the specified line of a text file, given a filename and a line number. It returns an empty string when the
   line number is larger than the number of lines in a file.

   Intended as an easy way to process a list of files. Use a command like 'ls > out' to create a file containing the
   list, in order, of files to be processed. Then call this function with an increasing index number until the return
   value is empty.

   +---------------+-----------------------------------------------------------------------------------------------------+
   | ``var``       | An ascii string, up to 128 characters long, containing the contents of line ``lineindex`` of the    |
   |               | input file.                                                                                         |
   +---------------+-----------------------------------------------------------------------------------------------------+
   | ``listname``  | The filename to open and read lines from.                                                           |
   +---------------+-----------------------------------------------------------------------------------------------------+
   | ``lineindex`` | Integer line number, starting at 1. If larger than the number of lines in the file, the empty       |
   |               | string '' will be returned.                                                                         |
   +---------------+-----------------------------------------------------------------------------------------------------+

| 

.. container:: routine

   *var = set_filename_list( name_array, listname, caller_name )*
   ::

      character(len=*),  intent(inout) :: name_array
      character(len=*),  intent(in)    :: listname
      character(len=*),  intent(in)    :: caller_name
      integer                          :: var

.. container:: indent1

   Returns the count of filenames specified. Verifies that one of either the name_array or the listname was specified
   but not both. If the input was a listname copy the names into the name_array so when this routine returns all the
   filenames are in name_array(). Verifies that no more than the allowed number of names was specified if the input was
   a listname file.

   +-----------------+---------------------------------------------------------------------------------------------------+
   | ``var``         | The count of input files specified.                                                               |
   +-----------------+---------------------------------------------------------------------------------------------------+
   | ``name_array``  | Array of input filename strings. Either this item or the listname must be specified, but not      |
   |                 | both.                                                                                             |
   +-----------------+---------------------------------------------------------------------------------------------------+
   | ``listname``    | The filename to open and read filenames from, one per line. Either this item or the name_array    |
   |                 | must be specified but not both.                                                                   |
   +-----------------+---------------------------------------------------------------------------------------------------+
   | ``caller_name`` | Calling subroutine name, used for error messages.                                                 |
   +-----------------+---------------------------------------------------------------------------------------------------+

| 

.. container:: routine

   *call to_upper(string)*
   ::

      character(len=*), intent (INOUT) :: string

.. container:: indent1

   Converts the character string to UPPERCASE - in place. The input string **is** modified.

   ========== ====================
   ``string`` any character string
   ========== ====================

| 

.. container:: routine

   *call nc_check(istatus, subr_name [, context])*
   ::

      integer, intent(in)                    :: istatus
      character(len=*), intent(in)           :: subr_name
      character(len=*), optional, intent(in) :: context

.. container:: indent1

   Check the return code from a netcdf call. If no error, return without taking any action. If an error is indicated (in
   the ``istatus`` argument) then call the error handler with the subroutine name and any additional context information
   (e.g. which file or which variable was being processed at the time of the error). All errors are currently hardcoded
   to be ``FATAL`` and this routine will not return.

   This routine calls a netCDF library routine to construct the text error message corresponding to the error code in
   the first argument. An example use of this routine is:
   ::

      call nc_check(nf90_create(path = trim(ncFileID%fname), cmode = nf90_share, ncid = ncFileID%ncid), &
                   'init_diag_output', 'create '//trim(ncFileID%fname))

   +---------------+-----------------------------------------------------------------------------------------------------+
   | ``istatus``   | The return value from any netCDF call.                                                              |
   +---------------+-----------------------------------------------------------------------------------------------------+
   | ``subr_name`` | String name of the current subroutine, used in case of error.                                       |
   +---------------+-----------------------------------------------------------------------------------------------------+
   | *context*     | Additional text to be used in the error message, for example to indicate which file or which        |
   |               | variable is being processed.                                                                        |
   +---------------+-----------------------------------------------------------------------------------------------------+

| 

.. container:: routine

   *call set_tasknum(tasknum)*
   ::

      integer, intent(in)               :: tasknum

.. container:: indent1

   Intended to be used in the MPI multi-task case. Sets the local task number, which is then prepended to subsequent
   messages.

   +-------------+-------------------------------------------------------------------------------------------------------+
   | ``tasknum`` | Task number returned from MPI_Comm_Rank(). MPI task numbers are 0 based, so for a 4-task job these    |
   |             | numbers are 0-3.                                                                                      |
   +-------------+-------------------------------------------------------------------------------------------------------+

| 

.. container:: routine

   *call set_output(doflag)*
   ::

      logical, intent(in)               :: doflag

.. container:: indent1

   Set the status of output. Can be set on a per-task basis if you are running with multiple tasks. If set to false only
   warnings and fatal errors will write to the log. The default in the multi-task case is controlled by the MPI module
   initialization code, which sets task 0 to .TRUE. and all other tasks to .FALSE.

   +------------+--------------------------------------------------------------------------------------------------------+
   | ``doflag`` | Sets, on a per-task basis, whether messages are to be written to the logfile or standard output.       |
   |            | Warnings and errors are always output.                                                                 |
   +------------+--------------------------------------------------------------------------------------------------------+

| 

.. container:: routine

   *var = do_output()*
   ::

      logical                      :: do_output

.. container:: indent1

   Returns true if this task should write to the log, false otherwise. Set by the ``set_output()`` routine. Defaults to
   true for the single task case. Can be used in code like so:

   ::

      if (do_output()) then
       write(*,*) 'At this point in the code'
      endif

   ======= ======================================
   ``var`` True if this task should write output.
   ======= ======================================

| 

.. container:: routine

   *call initialize_utilities( [progname] [, alternatename] )*
   ::

      character(len=*), intent(in), optional :: progname
      character(len=*), intent(in), optional :: alternatename

.. container:: indent1

   Reads the namelist and opens the logfile. Records the values of the namelist and registers this module.

   +-----------------+---------------------------------------------------------------------------------------------------+
   | *progname*      | If given, use in the timestamp message in the log file to say which program is being started.     |
   +-----------------+---------------------------------------------------------------------------------------------------+
   | *alternatename* | If given, log filename to use instead of the value in the namelist. This permits, for example,    |
   |                 | different programs sharing the same input.nml file to have different logs. If not given here and  |
   |                 | no value is specified in the namelist, this defaults to dart_log.out                              |
   +-----------------+---------------------------------------------------------------------------------------------------+

| 

.. container:: routine

   *call finalize_utilities()*

.. container:: indent1

   Closes the logfile; using utilities after this call is a bad idea.

| 

.. container:: routine

   *call dump_unit_attributes(iunit)*
   ::

      integer, intent(in) :: iunit

.. container:: indent1

   Writes all information about the status of the IO unit to the error handler with error level message.

   ========= ==========================================
   ``iunit`` Unit about which information is requested.
   ========= ==========================================

| 

.. container:: routine

   ::

      integer :: E_DBG, DEBUG
      integer :: E_MSG, MESSAGE
      integer :: E_WARN, WARNING
      integer :: E_ERR, FATAL

.. container:: indent1

   +---+-----------------------------------------------------------------------------------------------------------------+
   |   | Severity levels to be passed to error handler. Levels are debug, message, warning and fatal. The namelist       |
   |   | parameter TERMLEVEL can be used to control at which level program termination should occur.                     |
   +---+-----------------------------------------------------------------------------------------------------------------+

| 

.. container:: routine

   ::

      integer :: logfileunit

.. container:: indent1

   =============== ==========================================
   ``logfileunit`` Unit opened to file for diagnostic output.
   =============== ==========================================

| 

.. container:: routine

   ::

      integer :: nmlfileunit

.. container:: indent1

   +-----------------+---------------------------------------------------------------------------------------------------+
   | ``nmlfileunit`` | Unit opened to file for diagnostic output of namelist files. Defaults to same as ``logfileunit``. |
   |                 | Provides the flexibility to log namelists to a separate file, reducing the clutter in the log     |
   |                 | files and perhaps increasing readability.                                                         |
   +-----------------+---------------------------------------------------------------------------------------------------+

| 

Files
-----

-  assim_model_mod.nml in input.nml
-  logfile, name specified in namelist

References
----------

-  none

Error codes and conditions
--------------------------

+-----------------------+--------------------------------------------------+--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------+
|        Routine        |                      Message                     |                                                                                                                                                                Comment                                                                                                                                                               |
+=======================+==================================================+======================================================================================================================================================================================================================================================================================================================================+
| get_unit              | No available units                               | Unable to open enough IO channels                                                                                                                                                                                                                                                                                                    |
+-----------------------+--------------------------------------------------+--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------+
| check_nml_error       | while reading namelist _____                     | Fatal error reading namelist. This could be caused by having an entry in the namelist input file that is not in the namelist, by having illegal values for namelist variables, or by a variety of other compiler dependent problems.                                                                                                 |
+-----------------------+--------------------------------------------------+--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------+
| find_namelist_in_file | Namelist entry &____ must exist in namelist_nml. | There must be an entry for the required namelist, for instance &filter_nml, in the input.nml namelist file. Even if no values are to be changed from the default, an entry like &filter_nml followed by a line containing only / is required.                                                                                        |
+-----------------------+--------------------------------------------------+--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------+
| find_namelist_in_file | Namelist input file: input.nml must exist        | The namelist input file (usually input.nml) must exist.                                                                                                                                                                                                                                                                              |
+-----------------------+--------------------------------------------------+--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------+
| check_namelist_read   | INVALID NAMELIST ENTRY: ___ in namelist ____     | While reading the namelist, either a bad entry was found or an end of file was encountered. The most confusing case is when a namelist is being read successfully but is not appropriately terminated with a /. The line printed out by the error message will be the start of the next namelist in the input.nml file in this case. |
+-----------------------+--------------------------------------------------+--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------+

Private components
------------------

N/A
