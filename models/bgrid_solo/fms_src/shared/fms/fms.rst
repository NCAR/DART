Module fms_mod
==============

Overview
--------

The fms module provides routines that are commonly used by most FMS modules.

.. container::

   Here is a summary of the functions performed by routines in the fms module.
   1. Output module version numbers to a common (``log``) file using a common format.
   2. Open specific types of files common to many FMS modules. These include namelist files, restart files, and 32-bit
   IEEE data files. There also is a matching interface to close the files. If other file types are needed the
   ``mpp_open`` and ``mpp_close`` interfaces in module
   ` <http://www.gfdl.noaa.gov/fms-cgi-bin/cvsweb.cgi/FMS/shared/mpp/models/bgrid_solo/fms_src/shared/mpp/mpp_io.html>`__
   must be used.
   3. Read and write distributed data to simple native unformatted files. This type of file (called a restart file) is
   used to checkpoint model integrations for a subsequent restart of the run.
   4. Time sections of code (using a wrapper for
   `mpp_mod <http://www.gfdl.noaa.gov/fms-cgi-bin/cvsweb.cgi/FMS/shared/mpp/models/bgrid_solo/fms_src/shared/mpp/mpp.html>`__).
   5. For convenience there are several routines published from the
   ` <http://www.gfdl.noaa.gov/fms-cgi-bin/cvsweb.cgi/FMS/shared/mpp/models/bgrid_solo/fms_src/shared/mpp/mpp.html>`__
   module. These are routines for getting processor numbers, commonly used I/O unit numbers, and error handling.

| 

Other modules used
------------------

.. container::

   ::

      mpp_mod
      mpp_domains_mod
           mpp_io_mod
           fms_io_mod

Public interface
----------------

.. container::

   fms_init:
      Initializes the FMS module and also calls the initialization routines for all modules in the MPP package. Will be
      called automatically if the user does not call it.
   fms_end:
      Calls the termination routines for all modules in the MPP package.
   file_exist:
      Checks the existence of a given file name.
   error_mesg:
      Print notes, warnings and error messages; terminates program for warning and error messages. (use error levels
      NOTE,WARNING,FATAL, see example below)
   check_nml_error:
      Checks the iostat argument that is returned after reading a namelist and determines if the error code is valid.
   write_version_number:
      Prints to the log file (or a specified unit) the (cvs) version id string and (cvs) tag name.
   mpp_clock_init:
      Returns an identifier for performance timing a section of code (similar to mpp_clock_id).
   lowercase:
      Convert character strings to all lower case.
   uppercase:
      Convert character strings to all upper case.
   string_array_index:
      match the input character string to a string in an array/list of character strings
   monotonic_array:
      Determines if a real input array has monotonically increasing or decreasing values.

| 

Public data
-----------

.. container::

   None.

Public routines
---------------

a. .. rubric:: Fms_init
      :name: fms_init

   ::

      call fms_init ( )

   **DESCRIPTION**
      Initialization routine for the fms module. It also calls initialization routines for the mpp, mpp_domains, and
      mpp_io modules. Although this routine will be called automatically by other fms_mod routines, users should
      explicitly call fms_init. If this routine is called more than once it will return silently. There are no
      arguments.

b. .. rubric:: Fms_end
      :name: fms_end

   ::

      call fms_end ( )

   **DESCRIPTION**
      Termination routine for the fms module. It also calls destructor routines for the mpp, mpp_domains, and mpp_io
      modules. If this routine is called more than once it will return silently. There are no arguments.

c. .. rubric:: File_exist
      :name: file_exist

   ::

       
      file_exist ( file_name )

   **DESCRIPTION**
      Checks the existence of the given file name. If the file_name string has zero length or the first character is
      blank return a false result.
   **INPUT**
      +-----------------------------------------------------------+-----------------------------------------------------------+
      | ``file_name``                                             | A file name (or path name) that is checked for existence. |
      |                                                           | [character]                                               |
      +-----------------------------------------------------------+-----------------------------------------------------------+

   **OUTPUT**
      +-----------------------------------------------------------+-----------------------------------------------------------+
      |                                                           | This function returns a logical result. If file_name      |
      |                                                           | exists the result is true, otherwise false is returned.   |
      |                                                           | If the length of character string "file_name" is zero or  |
      |                                                           | the first character is blank, then the returned value     |
      |                                                           | will be false. When reading a file, this function is      |
      |                                                           | often used in conjunction with routine open_file.         |
      |                                                           | []                                                        |
      +-----------------------------------------------------------+-----------------------------------------------------------+

d. .. rubric:: Error_mesg
      :name: error_mesg

   ::

      call error_mesg ( routine, message, level )

   **DESCRIPTION**
      Print notes, warnings and error messages; and terminates the program for error messages. This routine is a wrapper
      around mpp_error, and is provided for backward compatibility. This module also publishes mpp_error, **users should
      try to use the mpp_error interface**.
   **INPUT**
      +-----------------------------------------------------------+-----------------------------------------------------------+
      | ``routine``                                               | Routine name where the warning or error has occurred.     |
      |                                                           | [character]                                               |
      +-----------------------------------------------------------+-----------------------------------------------------------+
      | ``message``                                               | Warning or error message to be printed.                   |
      |                                                           | [character]                                               |
      +-----------------------------------------------------------+-----------------------------------------------------------+
      | ``level``                                                 | Level of severity; set to NOTE, WARNING, or FATAL         |
      |                                                           | Termination always occurs for FATAL, never for NOTE, and  |
      |                                                           | is settable for WARNING (see namelist).                   |
      |                                                           | [integer]                                                 |
      +-----------------------------------------------------------+-----------------------------------------------------------+

   **NOTE**
      Examples:

      ::

                 use fms_mod, only: error_mesg, FATAL, NOTE
                 call error_mesg ('fms_mod', 'initialization not called', FATAL)
                 call error_mesg ('fms_mod', 'fms_mod message', NOTE)

e. .. rubric:: Check_nml_error
      :name: check_nml_error

   ::

       
      check_nml_error ( iostat, nml_name )

   **DESCRIPTION**
      The FMS allows multiple namelist records to reside in the same file. Use this interface to check the iostat
      argument that is returned after reading a record from the namelist file. If an invalid iostat value is detected
      this routine will produce a fatal error. See the NOTE below.
   **INPUT**
      +-----------------------------------------------------------+-----------------------------------------------------------+
      | ``iostat``                                                | The iostat value returned when reading a namelist record. |
      |                                                           | [integer]                                                 |
      +-----------------------------------------------------------+-----------------------------------------------------------+
      | ``nml_name``                                              | The name of the namelist. This name will be printed if an |
      |                                                           | error is encountered, otherwise the name is not used.     |
      |                                                           | [character]                                               |
      +-----------------------------------------------------------+-----------------------------------------------------------+

   **OUTPUT**
      +-----------------------------------------------------------+-----------------------------------------------------------+
      |                                                           | This function returns the input iostat value (integer) if |
      |                                                           | it is an allowable error code. If the iostat error code   |
      |                                                           | is not allowable, an error message is printed and the     |
      |                                                           | program terminated.                                       |
      |                                                           | [integer]                                                 |
      +-----------------------------------------------------------+-----------------------------------------------------------+

   **NOTE**
      | Some compilers will return non-zero iostat values when reading through files with multiple namelist. This
        routine will try skip these errors and only terminate for true namelist errors.
      | Examples
      | The following example checks if a file exists, reads a namelist input from that file, and checks for errors in
        that namelist. When the correct namelist is read and it has no errors the routine check_nml_error will return
        zero and the while loop will exit. This code segment should be used to read namelist files.

      ::

                   integer :: unit, ierr, io

                   if ( file_exist('input.nml') ) then
                       unit = open_namelist_file ( )
                       ierr=1
                       do while (ierr /= 0)
                         read  (unit, nml=moist_processes_nml, iostat=io, end=10)
                         ierr = check_nml_error(io,'moist_processes_nml')
                       enddo
                 10    call close_file (unit)
                   endif

f. .. rubric:: Write_version_number
      :name: write_version_number

   ::

      call write_version_number ( version [, tag, unit] )

   **DESCRIPTION**
      Prints to the log file (stdlog) or a specified unit the (cvs) version id string and (cvs) tag name.
   **INPUT**
      +-----------------------------------------------------------+-----------------------------------------------------------+
      | ``version``                                               | string that contains routine name and version number.     |
      |                                                           | [character(len=*)]                                        |
      +-----------------------------------------------------------+-----------------------------------------------------------+
      | ``tag``                                                   | The tag/name string, this is usually the Name string      |
      |                                                           | returned by CVS when checking out the code.               |
      |                                                           | [character(len=*)]                                        |
      +-----------------------------------------------------------+-----------------------------------------------------------+
      | ``unit``                                                  | The Fortran unit number of an open formatted file. If     |
      |                                                           | this unit number is not supplied the log file unit number |
      |                                                           | is used (stdlog).                                         |
      |                                                           | [integer]                                                 |
      +-----------------------------------------------------------+-----------------------------------------------------------+

g. .. rubric:: Mpp_clock_init
      :name: mpp_clock_init

   ::

      id = mpp_clock_init ( name, level [, flags] )

   **DESCRIPTION**
      Returns an identifier for performance timing sections of code. Should be used in conjunction with mpp_clock_begin
      and mpp_clock_end. For more details see the documentation for the MPP module and look at the example below.
   **INPUT**
      +-----------------------------------------------------------+-----------------------------------------------------------+
      | ``name``                                                  | A unique name string given to the code segment to be      |
      |                                                           | timed. The length should not exceed 32 characters.        |
      |                                                           | [character]                                               |
      +-----------------------------------------------------------+-----------------------------------------------------------+
      | ``level``                                                 | Level of timing. When level > timing_level, which is set  |
      |                                                           | by namelist &fms_nml, an identifier of zero is returned.  |
      |                                                           | This will turn off performance timing for the code        |
      |                                                           | section.                                                  |
      |                                                           | [integer]                                                 |
      +-----------------------------------------------------------+-----------------------------------------------------------+
      | ``flags``                                                 | Use the flags published via the mpp_mod to control        |
      |                                                           | whether synchronization or extra detail is desired.       |
      |                                                           | (flags = MPP_CLOCK_SYNC, MPP_CLOCK_DETAILED)              |
      |                                                           | [integer]                                                 |
      +-----------------------------------------------------------+-----------------------------------------------------------+

   **OUTPUT**
      +-----------------------------------------------------------+-----------------------------------------------------------+
      | ``id``                                                    | The identification index returned by mpp_clocks_id. A     |
      |                                                           | zero value is returned (turning clocks off) when input    |
      |                                                           | argument level > namelist variable timing_level.          |
      |                                                           | [integer]                                                 |
      +-----------------------------------------------------------+-----------------------------------------------------------+

   **NOTE**
      | 1.The MPP_CLOCK_SYNC flag should be used whenever possible. This flag causes mpp_sync to be called at the begin
        of a code segment, resulting in more accurate performance timings. **Do not use the MPP_CLOCK_SYNC flag for code
        sections that may not be called on all processors.**
      | 2.There is some amount of coordination required throughout an entire program for consistency of the "timing
        levels". As a guideline the following levels may be used, with higher levels added as desired to specific
        component models.

      ::

                         level 
                               example code section
                          1 
                               main program
                          2 
                               components models
                          3 
                               atmosphere dynamics or physics

      | Examples:
      | The mpp_clock_init interface should be used in conjunction with the mpp_mod interfaces mpp_clock_begin and
        mpp_clock_end. For example:

      ::

         use fms_mod, only: mpp_clock_init, mpp_clock_begin, &
                                      mpp_clock_end. MPP_CLOCK_SYNC
                   integer :: id_mycode
                   integer :: timing_level = 5

                   id_mycode = mpp_clock_init ('mycode loop', timing_level, &
                                               flags=MPP_CLOCK_SYNC)
                   call mpp_clock_begin (id_mycode)
                                 :
                                 :
                    ~~ this code will be timed ~~ 
                                 :
                                 :
                   call mpp_clock_end (id_mycode)

h. .. rubric:: Lowercase
      :name: lowercase

   ::

      string = lowercase ( cs )

   **DESCRIPTION**
      Converts a character string to all lower case letters. The characters "A-Z" are converted to "a-z", all other
      characters are left unchanged.
   **INPUT**
      +-----------------------------------------------------------+-----------------------------------------------------------+
      | ``cs``                                                    | Character string that may contain upper case letters.     |
      |                                                           | [character(len=*), scalar]                                |
      +-----------------------------------------------------------+-----------------------------------------------------------+

   **OUTPUT**
      +-----------------------------------------------------------+-----------------------------------------------------------+
      | ``string``                                                | Character string that contains all lower case letters.    |
      |                                                           | The length of this string must be the same as the input   |
      |                                                           | string.                                                   |
      |                                                           | [character(len=len(cs)), scalar]                          |
      +-----------------------------------------------------------+-----------------------------------------------------------+

i. .. rubric:: Uppercase
      :name: uppercase

   ::

      string = uppercase ( cs )

   **DESCRIPTION**
      Converts a character string to all upper case letters. The characters "a-z" are converted to "A-Z", all other
      characters are left unchanged.
   **INPUT**
      +-----------------------------------------------------------+-----------------------------------------------------------+
      | ``cs``                                                    | Character string that may contain lower case letters.     |
      |                                                           | [character(len=*), scalar]                                |
      +-----------------------------------------------------------+-----------------------------------------------------------+

   **OUTPUT**
      +-----------------------------------------------------------+-----------------------------------------------------------+
      | ``string``                                                | Character string that contains all upper case letters.    |
      |                                                           | The length of this string must be the same as the input   |
      |                                                           | string.                                                   |
      |                                                           | [character(len=len(cs)), scalar]                          |
      +-----------------------------------------------------------+-----------------------------------------------------------+

j. .. rubric:: String_array_index
      :name: string_array_index

   ::

       
      string_array_index ( string, string_array [, index] )

   **DESCRIPTION**
      Tries to find a match for a character string in a list of character strings. The match is case sensitive and
      disregards blank characters to the right of the string.
   **INPUT**
      +-----------------------------------------------------------+-----------------------------------------------------------+
      | ``string``                                                | Character string of arbitrary length.                     |
      |                                                           | [character(len=*), scalar]                                |
      +-----------------------------------------------------------+-----------------------------------------------------------+
      | ``string_array``                                          | Array/list of character strings.                          |
      |                                                           | [character(len=*), dimension(:)]                          |
      +-----------------------------------------------------------+-----------------------------------------------------------+

   **OUTPUT**
      +-----------------------------------------------------------+-----------------------------------------------------------+
      | ``index``                                                 | The index of string_array where the first match was       |
      |                                                           | found. If no match was found then index = 0.              |
      |                                                           | []                                                        |
      +-----------------------------------------------------------+-----------------------------------------------------------+
      | ``found``                                                 | If an exact match was found then TRUE is returned,        |
      |                                                           | otherwise FALSE is returned.                              |
      |                                                           | [logical]                                                 |
      +-----------------------------------------------------------+-----------------------------------------------------------+

   **NOTE**
      Examples

      ::

                string = "def"
                string_array = (/ "abcd", "def ", "fghi" /)
                string_array_index ( string, string_array, index )
                Returns: TRUE, index = 2

k. .. rubric:: Monotonic_array
      :name: monotonic_array

   ::

       
      monotonic_array ( array [, direction] )

   **DESCRIPTION**
      Determines if the real input array has monotonically increasing or decreasing values.
   **INPUT**
      +-----------------------------------------------------------+-----------------------------------------------------------+
      | ``array``                                                 | An array of real values. If the size(array) < 2 this      |
      |                                                           | function assumes the array is not monotonic, no fatal     |
      |                                                           | error will occur.                                         |
      |                                                           | [real, dimension(:)]                                      |
      +-----------------------------------------------------------+-----------------------------------------------------------+

   **OUTPUT**
      +-----------------------------------------------------------+-----------------------------------------------------------+
      | ``direction``                                             | If the input array is: >> monotonic (small to large) then |
      |                                                           | direction = +1. >> monotonic (large to small) then        |
      |                                                           | direction = -1. >> not monotonic then direction = 0.      |
      |                                                           | [integer]                                                 |
      +-----------------------------------------------------------+-----------------------------------------------------------+
      |                                                           | If the input array of real values either increases or     |
      |                                                           | decreases monotonically then TRUE is returned, otherwise  |
      |                                                           | FALSE is returned.                                        |
      |                                                           | [logical]                                                 |
      +-----------------------------------------------------------+-----------------------------------------------------------+

Namelist
--------

.. container::

   **&fms_nml**

   .. container::

      ``timing_level``
      The level of performance timing. If calls to the performance timing routines have been inserted into the code then
      code sections with a level <= timing_level will be timed. The resulting output will be printed to STDOUT. See the
      MPP module or mpp_clock_init for more details.
      [integer, default: 0]
      ``read_all_pe``
      Read global data on all processors extracting local part needed (TRUE) or read global data on PE0 and broadcast to
      all PEs (FALSE).
      [logical, default: true]
      ``warning_level``
      Sets the termination condition for the WARNING flag to interfaces error_mesg/mpp_error. set warning_level =
      'fatal' (program crashes for warning messages) or 'warning' (prints warning message and continues).
      [character, default: 'warning']
      ``iospec_ieee32``
      iospec flag used with the open_ieee32_file interface.
      [character, default: '-F f77,cachea:48:1']
      ``stack_size``
      The size in words of the MPP user stack. If stack_size > 0, the following MPP routine is called: call
      mpp_set_stack_size (stack_size). If stack_size = 0 (default) then the default size set by mpp_mod is used.
      [integer, default: 0]
      ``domains_stack_size``
      The size in words of the MPP_DOMAINS user stack. If domains_stack_size > 0, the following MPP_DOMAINS routine is
      called: call mpp_domains_set_stack_size (domains_stack_size). If domains_stack_size = 0 (default) then the default
      size set by mpp_domains_mod is used.
      [integer, default: 0]

| 

Data sets
---------

.. container::

   None.

Error messages
--------------

.. container::

   **FATAL in fms_init**
      invalid entry for namelist variable warning_level
      The namelist variable warning_level must be either 'fatal' or 'warning' (case-insensitive).
   **FATAL in file_exist**
      set_domain not called
      Before calling write_data you must first call set_domain with domain2d data type associated with the distributed
      data you are writing.
   **FATAL in check_nml_error**
      while reading namelist ...., iostat = ####
      There was an error message reading the namelist specified. Carefully examine all namelist variables for
      misspellings of type mismatches (e.g., integer vs. real).

References
----------

.. container::

   None.

| 

Compiler specifics
------------------

.. container::

   None.

| 

Precompiler options
-------------------

.. container::

   None.

| 

Loader options
--------------

.. container::

   None.

Test PROGRAM
------------

.. container::

   None.

| 

Notes
-----

.. container::

   1) If the **MPP** or **MPP_DOMAINS** stack size is exceeded the program will terminate after printing the required
   size.
   2) When running on a very small number of processors or for high resolution models the default domains_stack_size
   will probably be insufficient.

| 
