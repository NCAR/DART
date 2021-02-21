<!--

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!                                                                   !!
!!                   GNU General Public License                      !!
!!                                                                   !!
!! This file is part of the Flexible Modeling System (FMS).          !!
!!                                                                   !!
!! FMS is free software; you can redistribute it and/or modify       !!
!! it and are expected to follow the terms of the GNU General Public !!
!! License as published by the Free Software Foundation.             !!
!!                                                                   !!
!! FMS is distributed in the hope that it will be useful,            !!
!! but WITHOUT ANY WARRANTY; without even the implied warranty of    !!
!! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the     !!
!! GNU General Public License for more details.                      !!
!!                                                                   !!
!! You should have received a copy of the GNU General Public License !!
!! along with FMS; if not, write to:                                 !!
!!          Free Software Foundation, Inc.                           !!
!!          59 Temple Place, Suite 330                               !!
!!          Boston, MA  02111-1307  USA                              !!
!! or see:                                                           !!
!!          http://www.gnu.org/licenses/gpl.txt                      !!
!!                                                                   !!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

-->
<!DOCTYPE html>
<html>
<head>
<meta name="generator" content=
"HTML Tidy for HTML5 for Apple macOS version 5.6.0">
<title>module utilities</title>
</head>
<body bgcolor="#AABBCC" text="#332211">
<div align="center"><font size="1"><a href="#INTERFACE">PUBLIC
INTERFACE</a> / <a href="#DATA_TYPES">DATA</a> / <a href=
"#ROUTINES">ROUTINES</a> / <a href="#NAMELIST">NAMELIST</a> /
<a href="#CHANGES">CHANGES</a> / <a href="#ERRORS">ERRORS</a> /
<a href="#NOTES">NOTES</a></font><br>
<br></div>
<hr>
<h1>Module Utilities</h1>
<a name="OVERVIEW" id="OVERVIEW">
<hr>
<h2>OVERVIEW</h2>
<!-- BEGIN OVERVIEW -->
<pre>

     A collection of routines for performing commonly used
     (predominately I/O) functions.
     
</pre></a><!-- END OVERVIEW -->
<!--------------------------------------------------------------------->
 <a name="DESCRIPTION" id="DESCRIPTION"><!-- BEGIN DESCRIPTION -->
<pre>

     The utilities module contains routines for opening files,
     error handling, and outputting version numbers.

</pre></a><!-- END DESCRIPTION -->
<!--------------------------------------------------------------------->
 <a name="MODULES_USED" id="MODULES_USED">
<hr>
<h2>OTHER MODULES USED</h2>
<!-- BEGIN MODULES_USED -->
<pre>

             mpp_mod
     mpp_domains_mod
          mpp_io_mod

</pre></a><!-- END MODULES_USED -->
<!--------------------------------------------------------------------->
 <a name="INTERFACE" id="INTERFACE">
<hr>
<h2>PUBLIC INTERFACE</h2>
<!-- BEGIN INTERFACE -->
<pre>

use utilities_mod [, only: utilities_init, utilities_end,
                           error_mesg, NOTE, WARNING, FATAL,
                           file_exist, open_file, close_file,
                           read_data, write_data,
                           check_nml_error, print_version_number,
                           get_my_pe, get_num_pes, sync_all_pes,
                           get_root_pe,
                           set_domain, get_domain_decomp,
                           mpp_sum, mpp_min, mpp_max,
                           set_system_clock, check_system_clock
                           check_sum

   utilities_init   Initializes the utilities module and also the mpp_io 
                    module (which initializes all mpp mods).
                    Will be called automatically if the user does not call it.

   utilities_end    Calls mpp_io_exit and check_system_clock.

   error_mesg       Print notes, warnings and error messages;
                    terminates program for warning and error messages.
                    (use error levels NOTE,WARNING,FATAL, see example below)

   file_exist       Function that returns whether a requested file name exists.

   open_file        Opens a given file name for I/O and returns a unit number.
                    If the file is already open the unit number is returned.

   close_file       Closes the requested unit number for I/O.

   read_data        Reads distributed data.

   write_data       Writes distributed data.

   check_nml_error  Checks the iostat returned when reading namelist input
                    and determines if the error code is valid, if not
                    the program is terminated.

   write_version_number    Prints out the version/id string and tag name
                           to a specified unit (presumably a log file).

   set_domain          Sets a pointer for a private utilities_mod
                       domain2d data type. This domain2d type is
                       used for subsequent calls by read_data and
                       write data.

   get_domain_decomp   Returns domain information from the domain2d
                       type set by the most recent call to set_domain.

   set_system_clock    Sets (or resets) a real-time clock.
                       This is automatically called by utilities_init.

   check_system_clock  Outputs the elapsed real-time since the last
                       call to set_system_clock (or utilities_init).

   lowercase           Convert character strings to all lower case

   uppercase           Convert character strings to all upper case

   mpp_clock_init      Sets up a identifier for performance timing
                         (similar to mpp_clock_id)


   check_sum           Returns a global (from all PEs) integer check sum.


   <b>--- Routines with exact functionality to mpp_mod routines ----</b>

   get_my_pe           Returns the pe number. (mpp_pe)

   get_root_pe         Returns the root pe number. (mpp_root_pe)

   get_num_pes         Returns the total number pe's. (mpp_npes)

   sync_all_pes        Synchronizes pe's. (mpp_sync)

   mpp_sum             Compute summation. (mpp_sum)
   mpp_min             Find minimum. (mpp_min)
   mpp_max             Find maximum. (mpp_max)

   check_sum           Performs integer check sums (mpp_chksum)

   mpp_clock_begin     Initiates performance timing (mpp_clock_begin)
   mpp_clock_end       Terminates performance timing (mpp_clock_end)


</pre></a><!-- END INTERFACE -->
<!--------------------------------------------------------------------->
 <a name="DATA_TYPES" id="DATA_TYPES">
<hr>
<h2>PUBLIC DATA</h2>
<!-- BEGIN DATA_TYPES -->
<pre>

     None.

</pre></a><!-- END DATA_TYPES -->
<!--------------------------------------------------------------------->
 <a name="ROUTINES" id="ROUTINES">
<hr>
<h2>PUBLIC ROUTINES</h2>
<!-- BEGIN ROUTINES -->
<pre>

<b>file_exist</b> (file_name)

   Input

       file_name    A file name (or path name) that is checked for existence.
                       [character]

   Returns

       This function returns a logical result. 
       If file_name exists the result is true, otherwise false is returned.
       When reading a file, this function is often used in conjunction with
       routine open_file.

---------------------------------------------------------------------

<b>call error_mesg</b> ( routine, message, level )

   Input

      routine   Routine name where the warning or error has occurred.
                   [character]

      message   Warning or error message to be printed to standard output.
                   [character]

      level     Level of severity; <b>set to NOTE, WARNING, or FATAL</b>
                Termination always occurs for FATAL, never for NOTE,
                and is setable for WARNING (see namelist).
                   [integer]

  Examples:

  use utilities_mod, only: error_mesg, FATAL, NOTE

  call error_mesg ('utilities_mod', 'initialization not called', FATAL)

  call error_mesg ('utilities_mod', 'utilities_mod message', NOTE)

---------------------------------------------------------------------

check_nml_error ( iostat, nml_name )

   Input

      iostat    The iostat value returned when reading a namelist record.
                  [integer]

      nml_name  The name of the namelist.
                  [character]

   Returns

      This function returns the input iostat value (integer) if it is an
      allowable error code.  If the iostat error code is not allowable,
      an error message is printed to standard output and the program
      terminated.

   Note

      See a <a href=
"ioerrors.txt">table of possible iostat error codes</a>
      you may encounter on the Cray T90 or an SGI workstation.

   Example

      The following example check if a file exists, reads a namelist input
      from that file, and checks for errors in that namelist.

      if ( file_exist('input.nml') ) then
          unit = open_file (file='input.nml')
          ierr=1
          do while (ierr /= 0)
            read  (unit, nml=moist_processes_nml, iostat=io, end=10)
            ierr = check_nml_error(io,'moist_processes_nml')
          enddo
    10    close (unit)
      endif

---------------------------------------------------------------------

unit = <b>open_file</b> ( file, form, action, access, threading, recl )

   Input

       file       A file name (or path name).
                    [character]

   Input (optional)

       form       The file format, use "formatted", "unformatted", "ascii",
                  "native", "ieee32", or "netcdf".  Formatted and ascii are
                  the same, unformated and native are the same. 
                  For ieee32 files, the iospec option for mpp_open is set
                  to the value in the utilities module namelist.
                    [character, default: "formatted"]

       action     I/O action to be performed, use "read", "write", or "append".
                  In this version there is no default, <b>this argument must be
                  supplied when using the open_file interface</b>.
                    [character, default: none ????]

       access     Direct or sequential file access, use "direct" or "sequential".
                    [character, default: "sequential"]

       threading  Single or multi-threading, use "single" or "multi".
                  For single threading, the no_header flag for mpp_open
                  is set to TRUE.
                    [character, default: "single"]

       recl       Record length for direct access files.
                    [integer, default: none]

   Returns

       This function opens the requested file and returns a fortran unit
       number (integer). If the file is already open, the unit assigned to
       that file is returned.

   Example

       When reading a file, this function is often used in conjunction with
       routine file_exist.

       if ( file_exist (`your_file') ) then
            unit = open_file ( file='your_file', ..... )
       endif

---------------------------------------------------------------------

<b>call close_file</b> ( unit, status )

---------------------------------------------------------------------

<b>call write_version_number</b> ( unit, version, tagname )

   Input

      unit       The fortran unit number of an open formatted file.
                 This file is usually a log file.
                    [integer]

      tagname    The tag/name string usually returned by CVS
                 when checking out the code.
                    [character]

      version    The version/id string usually returned by CVS
                 when checking out code.
                    [character]

---------------------------------------------------------------------

<b>call set_domain</b> ( Domain2 )

   Input

      Domain2   Domain decomposition information.
                 [type(domain2d), target]

This interface must be called prior to calling interfaces
get_domain_decomp, read_data, or write_data.

---------------------------------------------------------------------

<b>call get_domain_decomp</b> ( x, y )

   Input

      x, y    Arrays containing the current (since the last call to
              set_domain) global and local start and end indices,
              for the x- and y-axis, respectivly.  The array values are
              (/ global_start, global_end, local_start, local_end /).
                 [integer, dimension(4)]

---------------------------------------------------------------------

<b>call read_data</b> ( unit, data, end )

   Input

      unit    The fortran unit number of an open unformatted file.
              This value is returned by a previous call to open_file.
                [integer]

   Output

      data    Distributed data to be read.
              The following data formats are allowed:
                [real   , dimension(isd:,jsd:)    ]
                [logical, dimension(isd:,jsd:)    ]
                [integer, dimension(isd:,jsd:)    ]
                [complex, dimension(isd:,jsd:)    ]
                [real   , dimension(isd:,jsd:,:)  ]
                [complex, dimension(isd:,jsd:,:)  ]
                [real   , dimension(isd:,jsd:,:,:)]
                [complex, dimension(isd:,jsd:,:,:)]

<b>NOTE:</b>  You must call set_domain before calling read_data.

---------------------------------------------------------------------

<b>call write_data</b> ( unit, data, end )

   Input

      unit       The fortran unit number of an open unformatted file.
                 This value is returned by a previous call to open_file.
                    [integer]

      data       Distributed data to be written.
                 The following data formats are allowed:
                    [real   , dimension(isd:,jsd:)    ]
                    [logical, dimension(isd:,jsd:)    ]
                    [integer, dimension(isd:,jsd:)    ]
                    [complex, dimension(isd:,jsd:)    ]
                    [real   , dimension(isd:,jsd:,:)  ]
                    [complex, dimension(isd:,jsd:,:)  ]
                    [real   , dimension(isd:,jsd:,:,:)]
                    [complex, dimension(isd:,jsd:,:,:)]

<b>NOTE:</b>  You must call set_domain before calling write_data.

---------------------------------------------------------------------
                                                            
<b>call set_system_clock</b>

There are no arguments.

---------------------------------------------------------------------
                                                            
<b>call check_system_clock</b> ( string )

   Optional Input

      string   An optional label that is printed before any timing numbers.

---------------------------------------------------------------------

<b>call utilities_init</b>

There are no arguments.

---------------------------------------------------------------------

<b>call utilities_end</b>

There are no arguments.

---------------------------------------------------------------------

string = <b>lowercase</b> ( cs )

   INPUT

      cs   character string [character(len=*)]

   RETURNS

      string   Character string that contain only lower case.
               The characters "A-Z" are converted to "a-z".

---------------------------------------------------------------------

string = <b>uppercase</b> ( cs )

   INPUT

      cs   character string [character(len=*)]

   RETURNS

      string   Character string that contain only upper case.
               The characters "a-z" are converted to "A-Z".

---------------------------------------------------------------------

id = <b>mpp_clock_init</b> ( module, routine, level )

   INPUT

      module    Module name to be timed. [character]

      routine   Routine name and label for section of code to be timed.
                   [character]

      level     Level of timing.  When this level equals "timing_level"
                defined by namelist &amp;utilities_nml a non-zero identification
                integer will be returned. 
                   [integer]

   RETURNS

      id        The identification index returned by mpp_clocks_id.
                A zero value is returned (turning clocks off) when input
                argument level does not equal namelist variable timing_level.

</pre></a><!-- END ROUTINES -->
<!--------------------------------------------------------------------->
 <a name="NAMELIST" id="NAMELIST">
<hr>
<h2>NAMELIST</h2>
<!-- BEGIN NAMELIST -->
<pre>

<b>&amp;utilities_nml</b>

  timing_level          Level/degree of detailed performance timing.
                        If timing_level = 0 then no timing is done.
                        The timing information will be saved to a separate
                        file for each processor.

  read_all_pe           Read global data on all processors extracting local
                        part needed (TRUE) or read global data on PE0 and
                        broadcast to all PEs (FALSE).
                           [logical, default: read_all_pe=.true.]

  warning_level         Sets the termination condition for the WARNING flag
                        to the interface error_mesg.
                        set warning_level = 'fatal' (program crashs) or
                        'warning' (print message and continues).
                           [character, default: warning_level='fatal']

  iospec_ieee32         iospec flag for all IEEE32 (32 bit) files.
                           [character, default: iospec_ieee32='-F f77,cachea:48:1']

   stack_size          The size in words of the MPP user stack.
                       If stack_size &gt; 0, the following MPP routine is called:
                       call mpp_set_stack_size (stack_size).
                       If stack_size = 0 (default) then the default size set
                       by mpp_mod is used.
                           [integer, default: stack_size = 0]

   domains_stack_size  The size in words of the MPP_DOMAINS user stack.
                       If domains_stack_size &gt; 0, the following
                       MPP_DOMAINS routine is called:
                       call mpp_domains_set_stack_size (domains_stack_size).
                       If domains_stack_size = 0 (default) then the default
                       size set by mpp_domains_mod is used.
                           [integer, default: domains_stack_size = 0]

NOTE:

    If the MPP or MPP_DOMAINS stack size is exceeded the program will
    terminate after printing the required size.


</pre></a><!-- END NAMELIST -->
<!--------------------------------------------------------------------->
 <a name="CHANGES" id="CHANGES">
<hr>
<h2>CHANGE HISTORY</h2>
<!-- BEGIN CHANGES -->
<pre>
<a href=".doc.log#utilities.f90">Revision history</a>

<b>changes</b> (5/1/2000)

   1) Added namelist variables read_all_pe, iospec_ieee32, and one_level_restarts.
   2) Added complex data overloads for interfaces read_data and write_data.


<b>changes</b> (11/23/1999)

   1) added a check_sum interface; uses the type domain2d from the
      previous call set_domain to get the global data onto a PE.
   2) Fixed some bugs in the namelist error checking, bugs still exist.
   3) Timing interface check_system_clock now prints an optional string.
   4) Subroutine utilities_end now calls mpp_io_exit.

</pre></a><!-- END CHANGES -->
<!--------------------------------------------------------------------->
 <a name="ERRORS" id="ERRORS">
<hr>
<h2>ERROR MESSAGES</h2>
<!-- BEGIN ERRORS -->
<pre>

<b>FATAL ERROR in check_nml_error</b>

    <b>while reading namelist ...., iostat = ####</b>
        There was an error message reading the namelist specified.
        Carefully examine all namelist variables for misspellings
        of type mismatches (e.g., integer vs. real).  Checking the
        iostat number with the table of iostat errors may also help.

<b>FATAL ERROR in read_data in utilities_mod</b>

    <b>set_domain not called</b>

<b>FATAL ERROR in write_data in utilities_mod</b>

    <b>set_domain not called</b>

<b>FATAL ERROR in utilities_init</b>

    <b>invalid entry for namelist variable warning_level</b>

</pre></a><!-- END ERRORS -->
<!--------------------------------------------------------------------->
 <a name="BUGS" id="BUGS">
<hr>
<h2>KNOWN BUGS</h2>
<!-- BEGIN BUGS -->
<pre>

   Namelist error checking may not work correctly with some compilers.

   Users should beware when mixing straight reads and read_data calls.
   If straight reads follow read_data calls and namelist variable
   read_all_pe = .false. (not the default), then the code will fail.
   It is safest if straight reads precede the read_data calls.

</pre></a><!-- END BUGS -->
<!--------------------------------------------------------------------->
 <a name="NOTES" id="NOTES">
<hr>
<h2>NOTES</h2>
<!-- BEGIN NOTES -->
<pre>

   Version numbers are printed out from the CVS keyword "$Id$".
   The tag that the code was checked out with is also included
   by using the keyword "$Name$". The following code can serve as
   a guideline for how to generate and output this information.

   character(len=128) :: version = '$Id$'
   character(len=128) :: tag = '$Name$'

   if (get_my_pe() == 0) &
        write (unit,'(/,80("="),/(a))') trim(version), trim(tag)

   This output format may now done with the routine
   <b>write_version_number</b>.

</pre></a><!-- END NOTES -->
<!--------------------------------------------------------------------->
 <a name="PLANS" id="PLANS">
<hr>
<h2>FUTURE PLANS</h2>
<!-- BEGIN PLANS -->
<pre>

   1) The name of this module will be changed to fms_utilities.f90
      (or something like that).
   2) MPP routines will be usable through this module.

</pre></a><!-- END PLANS -->
<!--------------------------------------------------------------------->
<hr>
</body>
</html>
