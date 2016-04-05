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

module fms_mod

! <CONTACT EMAIL="bw@gfdl.noaa.gov">
!   Bruce Wyman
! </CONTACT>

! <HISTORY SRC="http://www.gfdl.noaa.gov/fms-cgi-bin/cvsweb.cgi/FMS/"/>

! <OVERVIEW>
!   The fms module provides routines that are commonly used
!   by most FMS modules.
! </OVERVIEW>

! <DESCRIPTION>
!   Here is a summary of the functions performed by routines
!     in the fms module.
!
!   1. Output module version numbers to a common (<TT>log</TT>) file
!     using a common format.<BR/>
!   2. Open specific types of files common to many FMS modules.
!     These include namelist files, restart files, and 32-bit IEEE
!     data files. There also is a matching interface to close the files.
!     If other file types are needed the <TT>mpp_open</TT> and <TT>mpp_close</TT>
!     interfaces in module <LINK SRC="http://www.gfdl.noaa.gov/fms-cgi-bin/cvsweb.cgi/FMS/shared/mpp/mpp_io.html"><TT>mpp_io</TT></LINK>must be used.<BR/>
!    3. Read and write distributed data to simple native unformatted files.
!     This type of file (called a restart file) is used to checkpoint
!     model integrations for a subsequent restart of the run.<BR/>
!    4. Time sections of code (using a wrapper for <LINK SRC="http://www.gfdl.noaa.gov/fms-cgi-bin/cvsweb.cgi/FMS/shared/mpp/mpp.html">mpp_mod</LINK>).<BR/>
!    5. For convenience there are several routines published from
!     the <LINK SRC="http://www.gfdl.noaa.gov/fms-cgi-bin/cvsweb.cgi/FMS/shared/mpp/mpp.html"><TT>mpp</TT></LINK> module. These are routines for getting processor
!     numbers, commonly used I/O unit numbers, and error handling.
! </DESCRIPTION>

!-----------------------------------------------------------------------
!
!         A collection of commonly used routines.
!
!  The routines are primarily I/O related, however, there also
!  exists several simple miscellaneous utility routines.
!
!-----------------------------------------------------------------------
!
!  file_exist         Checks the existence of the given file name.
!
!  check_nml_error    Checks the iostat argument that is returned after
!                     reading a namelist and determines if the error
!                     code is valid.
!
!  write_version_number  Prints to the log file (or a specified unit)
!                        the (cvs) version id string and (cvs) tag name.
!
!  error_mesg          Print notes, warnings and error messages, 
!                      terminates program for error messages.
!                      (use error levels NOTE,WARNING,FATAL)
!
!  open_namelist_file  Opens namelist file for reading only.
!
!  open_restart_file   Opens a file that will be used for reading or writing
!                      restart files with native unformatted data.
!
!  open_ieee32_file    Opens a file that will be used for reading or writing
!                      unformatted 32-bit ieee data.
!
!  close_file          Closes a file that was opened using 
!                      open_namelist_file, open_restart_file, or
!                      open_ieee32_file.
!
!  set_domain          Call this routine to internally store in fms_mod the
!                      domain2d data type prior to calling the distributed
!                      data I/O routines read_data and write_data.
!
!  read_data           Reads distributed data from a single threaded file.
!
!  write_data          Writes distributed data to a single threaded file.
!
!  fms_init            Initializes the fms module and also the
!                      mpp_io module (which initializes all mpp mods).
!                      Will be called automatically if the user does
!                      not call it.
!
!  fms_end             Calls mpp exit routines.
!
!  lowercase           Convert character strings to all lower case
!
!  uppercase           Convert character strings to all upper case
!
!  monotonic_array     Determines if the real input array has
!                      monotonically increasing or decreasing values.
!
!  string_array_index  Match the input character string to a string
!                      in an array/list of character strings.
!
!  mpp_clock_init      Sets up a identifier for performance timing
!                        (similar to mpp_clock_id)
!
!-----------------------------------------------------------------------
!
!-----------------------------------------------------------------------

use types_mod, only : r8
use fms_io_mod, only : fms_io_init, fms_io_exit, &
                       open_ieee32_file, open_namelist_file, open_restart_file, close_file
use utilities_mod, only : get_unit, find_namelist_in_file, check_namelist_read

implicit none
private

! routines for initialization and termination of module
public :: fms_init, fms_end

! routines for opening/closing specific types of file
public :: open_namelist_file, open_restart_file, &
          close_file, open_ieee32_file

! miscellaneous i/o routines
public :: file_exist,      &
          write_version_number, error_mesg

! miscellaneous utilities (non i/o)
public :: lowercase, uppercase,                &
          string_array_index, monotonic_array

! public mpp interfaces
public ::  NOTE, WARNING, FATAL, stdlog
           

!------ namelist interface -------
!------ adjustable severity level for warnings ------

  integer           :: timing_level  = 0
  logical           :: read_all_pe   = .true.
  character(len=8)  :: warning_level = 'warning'
  character(len=64) :: iospec_ieee32 = '-N ieee_32'
  integer           :: stack_size = 0
  integer           :: domains_stack_size = 0

integer, parameter :: NOTE = 0, WARNING = 1, FATAL = 2

  character(len=32) :: configfile='logfile.out'
integer :: log_unit

!------ namelist interface -------

! <NAMELIST NAME="fms_nml">
!   <DATA NAME="timing_level"  TYPE="integer"  DEFAULT="0">
!     The level of performance timing. If calls to the performance timing routines 
!     have been inserted into the code then code sections with a level &lt;= 
!     timing_level will be timed. The resulting output will be printed to STDOUT.
!     See the MPP module or <LINK SRC="#mpp_clock_init">mpp_clock_init</LINK> for 
!     more details.
!   </DATA>
!   <DATA NAME="read_all_pe"  TYPE="logical"  DEFAULT="true">
!     Read global data on all processors extracting local part needed (TRUE) or
!     read global data on PE0 and broadcast to all PEs (FALSE).
!   </DATA>
!   <DATA NAME="warning_level"  TYPE="character"  DEFAULT="'warning'">
!     Sets the termination condition for the WARNING flag to interfaces
!     error_mesg/mpp_error. set warning_level = 'fatal' (program crashes for
!     warning messages) or 'warning' (prints warning message and continues).
!   </DATA>
!   <DATA NAME="iospec_ieee32"  TYPE="character"  DEFAULT="'-F f77,cachea:48:1'">
!     iospec flag used with the open_ieee32_file interface.
!   </DATA>
!   <DATA NAME="stack_size"  TYPE="integer"  DEFAULT="0">
!     The size in words of the MPP user stack. If stack_size > 0, the following
!     MPP routine is called: call mpp_set_stack_size (stack_size). If stack_size
!     = 0 (default) then the default size set by mpp_mod is used.
!   </DATA>
!   <DATA NAME="domains_stack_size" TYPE="integer"  DEFAULT="0">
!     The size in words of the MPP_DOMAINS user stack. If
!     domains_stack_size > 0, the following MPP_DOMAINS routine is called:
!     call mpp_domains_set_stack_size (domains_stack_size). If
!     domains_stack_size = 0 (default) then the default size set by
!     mpp_domains_mod is used. 
!   </DATA>
! </NAMELIST>

  namelist /fms_nml/  timing_level, read_all_pe,    &
                      warning_level, iospec_ieee32, &
                      stack_size, domains_stack_size


!   ---- private data for check_nml_error ----

!   integer, private :: num_nml_error_codes, nml_error_codes(20)
!   logical, private :: do_nml_error_init = .true.
!   private  nml_error_init


!  ---- version number -----

  character(len=128) :: version = '$Revision$'
  character(len=128) :: tagname = '$Id$'

  logical :: module_is_initialized = .FALSE.


contains

!#######################################################################

! <SUBROUTINE NAME="fms_init">

!   <OVERVIEW>
!     Initializes the FMS module and also calls the initialization routines for all
!     modules in the MPP package. Will be called automatically if the user does
!     not call it. 
!   </OVERVIEW>
!   <DESCRIPTION>
!      Initialization routine for the fms module. It also calls initialization routines
!      for the mpp, mpp_domains, and mpp_io modules. Although this routine
!      will be called automatically by other fms_mod routines, users should
!      explicitly call fms_init. If this routine is called more than once it will
!      return silently. There are no arguments.
!   </DESCRIPTION>
!   <TEMPLATE>
!     call fms_init ( )
!   </TEMPLATE>


!   <ERROR MSG="invalid entry for namelist variable warning_level" STATUS="FATAL">
!     The namelist variable warning_level must be either 'fatal' or 'warning'
!     (case-insensitive). 
!   </ERROR>

! initializes the fms module/package
! also calls mpp initialization routines and reads fms namelist

subroutine fms_init ( )

 integer :: unit, ierr, io, iunit

    if (module_is_initialized) return    ! return silently if already called
    module_is_initialized = .true.
!---- initialize mpp routines ----

!!!    call mpp_init
!!!    call mpp_domains_init
    call fms_io_init

!---- read namelist input ----

!    call nml_error_init  ! first initialize namelist iostat error codes

! 8 June, 2006 Switching to better namelist reading; original block should go away
call find_namelist_in_file("input.nml", "fms_nml", iunit)
read(iunit, nml=fms_nml, iostat = io)
call check_namelist_read(iunit, io, "fms_nml")

!!!    if (file_exist('input.nml')) then
!!!       unit = open_namelist_file ( )
!!!       ierr=1; do while (ierr /= 0)
!!!          read  (unit, nml=fms_nml, iostat=io, end=10)
!!!          ierr = check_nml_error(io,'fms_nml')  ! also initializes nml error codes
!!!       enddo
!!! 10    close (unit)
!!!    endif

!---- define mpp stack sizes if non-zero -----

!!!    if (        stack_size > 0) call         mpp_set_stack_size (        stack_size)
!!!    if (domains_stack_size > 0) call mpp_domains_set_stack_size (domains_stack_size)

!---- set severity level for warnings ----

!    if ( lowercase(trim(warning_level)) == 'fatal' ) then
!            call mpp_set_warn_level ( FATAL )
!    else if ( lowercase(trim(warning_level)) == 'warning' ) then
!            call mpp_set_warn_level ( WARNING )
!    else
!            call error_mesg ( 'fms_init',  &
!            'invalid entry for namelist variable warning_level', FATAL )
!    endif

!--- write version info and namelist to logfile ---

    call write_version_number (version, tagname)
      write (stdlog(), nml=fms_nml)
!!!      write (stdlog(),*) 'nml_error_codes=', nml_error_codes(1:num_nml_error_codes)


end subroutine fms_init
! </SUBROUTINE>

!#######################################################################


! <SUBROUTINE NAME="fms_end">

!   <OVERVIEW>
!     Calls the termination routines for all modules in the MPP package.
!   </OVERVIEW>
!   <DESCRIPTION>
!     Termination routine for the fms module. It also calls destructor routines
!      for the mpp, mpp_domains, and mpp_io modules. If this routine is called
!      more than once it will return silently. There are no arguments. 
!   </DESCRIPTION>
!   <TEMPLATE>
!     call fms_end ( )
!   </TEMPLATE>

! terminates the fms module/package
! also calls mpp destructor routines

subroutine fms_end ( )

    if (.not.module_is_initialized) return  ! return silently
    call fms_io_exit
!!!    call mpp_domains_exit
!!!    call mpp_exit
    module_is_initialized =.FALSE.

end subroutine fms_end
! </SUBROUTINE>

!#######################################################################
! check the existence of the given file name
! if the file_name string has zero length or the
! first character is blank return a false result
! <FUNCTION NAME="file_exist">

!   <OVERVIEW>
!     Checks the existence of a given file name.
!   </OVERVIEW>
!   <DESCRIPTION>
!     Checks the existence of the given file name.
!     If the file_name string has zero length or the
!     first character is blank return a false result.
!   </DESCRIPTION>
!   <TEMPLATE>
!     file_exist ( file_name )
!   </TEMPLATE>

!   <IN NAME="file_name"  TYPE="character" >
!     A file name (or path name) that is checked for existence.
!   </IN>
!   <OUT NAME=""  TYPE="" >
!     This function returns a logical result.  If file_name exists the result 
!     is true, otherwise false is returned.
!     If the length of character string "file_name" is zero or the first
!     character is blank, then the returned value will be false.
!     When reading a file, this function is often used in conjunction with
!     routine open_file.
!   </OUT>
!   <ERROR MSG="set_domain not called" STATUS="FATAL">
!     Before calling write_data you must first call set_domain with domain2d data 
!     type associated with the distributed data you are writing.
!   </ERROR>

 function file_exist (file_name)
  character(len=*), intent(in) :: file_name
  logical  file_exist

   file_exist = .false.
   if (len_trim(file_name) == 0) return
   if (file_name(1:1) == ' ')    return

   inquire (file=trim(file_name), exist=file_exist)

 end function file_exist
! </FUNCTION>

!#######################################################################
! <SUBROUTINE NAME="error_mesg">

!   <OVERVIEW>
!     Print notes, warnings and error messages; terminates program for warning 
!     and error messages. (use error levels NOTE,WARNING,FATAL, see example below)
!   </OVERVIEW>
!   <DESCRIPTION>
!     Print notes, warnings and error messages; and terminates the program for 
!     error messages. This routine is a wrapper around mpp_error, and is provided 
!     for backward compatibility. This module also publishes mpp_error,
!      <B>users should try to use the mpp_error interface</B>. 
!   </DESCRIPTION>
!   <TEMPLATE>
!     call error_mesg ( routine, message, level )
!   </TEMPLATE>

!   <IN NAME="routine"  TYPE="character" >
!     Routine name where the warning or error has occurred.
!   </IN>
!   <IN NAME="message"  TYPE="character" >
!     Warning or error message to be printed.
!   </IN>
!   <IN NAME="level"  TYPE="integer" >
!     Level of severity; set to NOTE, WARNING, or FATAL Termination always occurs 
!     for FATAL, never for NOTE, and is settable for WARNING (see namelist).
!   </IN>
!   <NOTE>
!
!     Examples:
!     <PRE>
!        use fms_mod, only: error_mesg, FATAL, NOTE

!        call error_mesg ('fms_mod', 'initialization not called', FATAL)
!        call error_mesg ('fms_mod', 'fms_mod message', NOTE)
!     </PRE>
!   </NOTE>
! wrapper for the mpp error handler
! users should try to use the mpp_error interface

 subroutine error_mesg (routine, message, level)
  character(len=*), intent(in) :: routine, message
  integer,          intent(in) :: level

character(len=128) :: text
!  input:
!      routine   name of the calling routine (character string)
!      message   message written to output   (character string)
!      level     set to NOTE, MESSAGE, or FATAL (integer)

    if (.not.module_is_initialized) call fms_init ( )

      select case( level)
      case(NOTE)
          text = 'NOTE'         !just FYI
      case(WARNING)
          text = 'WARNING'      !probable error
      case(FATAL)
          text = 'FATAL'        !fatal error
      case default
          text = 'WARNING: non-existent level (must be NOTE|WARNING|FATAL)'
      end select
      text = trim(text)//': '//trim(routine) //': ' //trim(message)

      select case( level )
      case(NOTE)
          write( *,'(a)' )trim(text)
      case default
          write( 0,'(/a/)' )trim(text)
          if( level.EQ.FATAL)then
              stop
          end if                                                                            
      end select                                                                            

 end subroutine error_mesg
! </SUBROUTINE>

!#######################################################################
! <FUNCTION NAME="check_nml_error">

!   <OVERVIEW>
!     Checks the iostat argument that is returned after reading a namelist 
!     and determines if the error code is valid. 
!   </OVERVIEW>
!   <DESCRIPTION>
!     The FMS allows multiple namelist records to reside in the same file. 
!     Use this interface to check the iostat argument that is returned after 
!     reading a record from the namelist file. If an invalid iostat value 
!     is detected this routine will produce a fatal error. See the NOTE below.
!   </DESCRIPTION>
!   <TEMPLATE>
!     check_nml_error ( iostat, nml_name )
!   </TEMPLATE>

!   <IN NAME="iostat"  TYPE="integer" >
!     The iostat value returned when reading a namelist record.
!   </IN>
!   <IN NAME="nml_name"  TYPE="character" >
!     The name of the namelist. This name will be printed if an error is 
!     encountered, otherwise the name is not used.
!   </IN>
!   <OUT NAME=""  TYPE="integer" >
!     This function returns the input iostat value (integer) if it is an 
!     allowable error code. If the iostat error code is not
!     allowable, an error message is printed and the program terminated.
!   </OUT>
!   <NOTE>
!     Some compilers will return non-zero iostat values when reading through 
!     files with multiple namelist. This routine
!     will try skip these errors and only terminate for true namelist errors.
!
!     Examples
!
!       The following example checks if a file exists, reads a namelist input 
!       from that file, and checks for errors in that
!       namelist. When the correct namelist is read and it has no errors the 
!       routine check_nml_error will return zero and the while loop will exit. 
!       This code segment should be used to read namelist files. 
!       <PRE>
!          integer :: unit, ierr, io
!
!          if ( file_exist('input.nml') ) then
!              unit = open_namelist_file ( )
!              ierr=1
!              do while (ierr /= 0)
!                read  (unit, nml=moist_processes_nml, iostat=io, end=10)
!                ierr = check_nml_error(io,'moist_processes_nml')
!              enddo
!        10    call close_file (unit)
!          endif
!       </PRE>
!   </NOTE>

!   <ERROR MSG="while reading namelist ...., iostat = ####" STATUS="FATAL">
!     There was an error message reading the namelist specified. Carefully 
!     examine all namelist variables for
!     misspellings of type mismatches (e.g., integer vs. real).
!   </ERROR>

! used to check the iostat argument that is
! returned after reading a namelist
! see the online documentation for how this routine might be used

! function check_nml_error (iostat, nml_name) result (error_code)
!
!  integer,          intent(in) :: iostat
!  character(len=*), intent(in) :: nml_name
!  integer   error_code, i
!  character(len=128) :: err_str
!
!   if (.not.module_is_initialized) call fms_init ( )
!
!   error_code = iostat
!
!   do i = 1, num_nml_error_codes
!        if (error_code == nml_error_codes(i)) return
!   enddo
!
!!  ------ fatal namelist error -------
!!  ------ only on root pe ----------------
!       write (err_str,*) 'while reading namelist ',  &
!                         trim(nml_name), ', iostat = ',error_code
!       call error_mesg ('check_nml_error in fms_mod', err_str, FATAL)
!       call error_mesg ('check_nml_error in fms_mod', err_str, FATAL)
!
!end function check_nml_error
! </FUNCTION>

!-----------------------------------------------------------------------
!   private routine for initializing allowable error codes

!subroutine nml_error_init
!
!! some compilers return non-zero iostat values while
!! reading through files with multiple namelist records
!! this routines "attempts" to identify the iostat values associated
!! with records not belonging to the requested namelist
!
!   integer  unit, io, ir
!   real    ::  a=1.
!   integer ::  b=1
!   logical ::  c=.true.
!   character(len=8) ::  d='testing'
!   namelist /b_nml/  a,b,c,d
!
!      nml_error_codes(1) = 0
!
!!     ---- create dummy namelist file that resembles actual ----
!!     ---- (each pe has own copy) ----
!      call mpp_open (unit, '_read_error.nml', form=MPP_ASCII,  &
!                     action=MPP_OVERWR, access=MPP_SEQUENTIAL, &
!                     threading=MPP_MULTI)
!!     ---- due to namelist bug this will not always work ---
!      write (unit, 10)
!  10  format ('    ', &
!             /' &a_nml  a=1.  /',    &
!             /'#------------------', &
!             /' &b_nml  a=5., b=0, c=.false., d=''test'',  &end')
!      call mpp_close (unit)
!
!!     ---- read namelist files and save error codes ----
!      call mpp_open (unit, '_read_error.nml', form=MPP_ASCII,  &
!                     action=MPP_RDONLY, access=MPP_SEQUENTIAL, &
!                     threading=MPP_MULTI)
!      ir=1; io=1; do
!         read  (unit, nml=b_nml, iostat=io, end=20)
!         if (io == 0) exit
!         ir=ir+1; nml_error_codes(ir)=io
!      enddo
!  20  call mpp_close (unit, action=MPP_DELETE)
!
!      num_nml_error_codes = ir
!      do_nml_error_init = .false.
!
!end subroutine nml_error_init

!#######################################################################
! <SUBROUTINE NAME="write_version_number">

!   <OVERVIEW>
!     Prints to the log file (or a specified unit) the (cvs) version id string and
!     (cvs) tag name.
!   </OVERVIEW>
!   <DESCRIPTION>
!     Prints to the log file (stdlog) or a specified unit the (cvs) version id string
!      and (cvs) tag name.
!   </DESCRIPTION>
!   <TEMPLATE>
!    call write_version_number ( version [, tag, unit] )
!   </TEMPLATE>

!   <IN NAME="version" TYPE="character(len=*)">
!    string that contains routine name and version number.
!   </IN>
!   <IN NAME="tag" TYPE="character(len=*)">
!    The tag/name string, this is usually the Name string
!    returned by CVS when checking out the code.
!   </IN>
!   <IN NAME="unit" TYPE="integer">
!    The Fortran unit number of an open formatted file. If this unit number 
!    is not supplied the log file unit number is used (stdlog). 
!   </IN>
! prints module version number to the log file of specified unit number

 subroutine write_version_number (version, tag, unit)

!   in:  version = string that contains routine name and version number
!
!   optional in:
!        tag = cvs tag name that code was checked out with
!        unit    = alternate unit number to direct output  
!                  (default: unit=stdlog)

   character(len=*), intent(in) :: version
   character(len=*), intent(in), optional :: tag 
   integer,          intent(in), optional :: unit 

   integer :: logunit 

   if (.not.module_is_initialized) call fms_init ( )

     logunit = stdlog()
     if (present(unit)) then
         logunit = unit
     endif   

     if (present(tag)) then
         write (logunit,'(/,80("="),/(a))') trim(version), trim(tag)
     else    
         write (logunit,'(/,80("="),/(a))') trim(version)
     endif   

 end subroutine write_version_number
! </SUBROUTINE>



!#######################################################################
!  functions for changing the case of character strings
!#######################################################################


! <FUNCTION NAME="lowercase">

!   <OVERVIEW>
!     Convert character strings to all lower case.
!   </OVERVIEW>
!   <DESCRIPTION>
!     Converts a character string to all lower case letters. The characters "A-Z"
!      are converted to "a-z", all other characters are left unchanged.
!   </DESCRIPTION>
!   <TEMPLATE>
!     string = lowercase ( cs )
!   </TEMPLATE>

!   <IN NAME="cs"  TYPE="character(len=*), scalar" >
!     Character string that may contain upper case letters.
!   </IN>
!   <OUT NAME="string"  TYPE="character(len=len(cs)), scalar" >
!     Character string that contains all lower case letters. The
!     length of this string must be the same as the input string.
!   </OUT>

!   change to all lower case

 function lowercase (cs) 
 character(len=*), intent(in) :: cs
 character(len=len(cs))       :: lowercase 
 character :: ca(len(cs)) 

 integer, parameter :: co=iachar('a')-iachar('A') ! case offset
    
    !ca = transfer(cs,"x",len(cs)) 
    !where (ca >= "A" .and. ca <= "Z") ca = achar(iachar(ca)+co) 
    !lowercase = transfer(ca,cs) 
    
    ! the transfer function is not implemented in the gfortran compiler.
    ! make this a no-op!
    lowercase = cs

 end function lowercase 
! </FUNCTION>

!#######################################################################


! <FUNCTION NAME="uppercase">

!   <OVERVIEW>
!     Convert character strings to all upper case.
!   </OVERVIEW>
!   <DESCRIPTION>
!     Converts a character string to all upper case letters. The characters "a-z"
!      are converted to "A-Z", all other characters are left unchanged. 
!   </DESCRIPTION>
!   <TEMPLATE>
!     string = uppercase ( cs )
!   </TEMPLATE>

!   <IN NAME="cs"  TYPE="character(len=*), scalar" >
!     Character string that may contain lower case letters.
!   </IN>
!   <OUT NAME="string"  TYPE="character(len=len(cs)), scalar" >
!     Character string that contains all upper case letters. The
!             length of this string must be the same as the input string.
!   </OUT>
!   change to all upper case

 function uppercase (cs) 
 character(len=*), intent(in) :: cs
 character(len=len(cs))       :: uppercase 
 character :: ca(len(cs)) 

 integer, parameter :: co=iachar('A')-iachar('a') ! case offset
    
    !ca = transfer(cs,"x",len(cs)) 
    !where (ca >= "a" .and. ca <= "z") ca = achar(iachar(ca)+co) 
    !uppercase = transfer(ca,cs) 

    ! the transfer function is not implemented in the gfortran compiler.
    ! make this a no-op!
    uppercase = cs
    
 end function uppercase 
! </FUNCTION>

!#######################################################################


! <FUNCTION NAME="string_array_index">

!   <OVERVIEW>
!     match the input character string to a string
!     in an array/list of character strings
!   </OVERVIEW>
!   <DESCRIPTION>
!      Tries to find a match for a character string in a list of character strings.
!      The match is case sensitive and disregards blank characters to the right of
!      the string. 
!   </DESCRIPTION>
!   <TEMPLATE>
!      string_array_index ( string, string_array [, index] )
!   </TEMPLATE>

!   <IN NAME="string"  TYPE="character(len=*), scalar" >
!     Character string of arbitrary length.
!   </IN>
!   <IN NAME="string_array"  TYPE="character(len=*)" DIM="(:)">
!     Array/list of character strings.
!   </IN>
!   <OUT NAME="index"  TYPE="" >
!     The index of string_array where the first match was found. If
!            no match was found then index = 0.
!   </OUT>
!   <OUT NAME="found"  TYPE="logical" >
!     If an exact match was found then TRUE is returned, otherwise FALSE is returned.
!   </OUT>
!   <NOTE>
!     Examples
!      <PRE>
!       string = "def"
!       string_array = (/ "abcd", "def ", "fghi" /)

!       string_array_index ( string, string_array, index )

!       Returns: TRUE, index = 2
!      </PRE>
!   </NOTE>
! match the input character string to a string
! in an array/list of character strings

function string_array_index ( string, string_array, index ) result (found)
character(len=*),  intent(in)  :: string, string_array(:)
integer, optional, intent(out) :: index
logical :: found
integer :: i

! initialize this function to false
! loop thru string_array and exit when a match is found

  found = .false.
  if (present(index)) index = 0

  do i = 1, size(string_array)
    ! found a string match ?
    if ( trim(string) == trim(string_array(i)) ) then
         found = .true.
         if (present(index)) index = i
         exit
    endif
  enddo

end function string_array_index
! </FUNCTION>

!#######################################################################

! <FUNCTION NAME="monotonic_array">

!   <OVERVIEW>
!     Determines if a real input array has monotonically increasing or
!     decreasing values.
!   </OVERVIEW>
!   <DESCRIPTION>
!     Determines if the real input array has monotonically increasing or
!     decreasing values.
!   </DESCRIPTION>
!   <TEMPLATE>
!     monotonic_array ( array [, direction] )
!   </TEMPLATE>

!   <IN NAME="array"  TYPE="real" DIM="(:)">
!     An array of real values. If the size(array) < 2 this function
!     assumes the array is not monotonic, no fatal error will occur.
!   </IN>
!   <OUT NAME="direction"  TYPE="integer" >
!     If the input array is:
!                >> monotonic (small to large) then direction = +1.
!                >> monotonic (large to small) then direction = -1.
!                >> not monotonic then direction = 0. 
!   </OUT>
!   <OUT NAME=""  TYPE="logical" >
!     If the input array of real values either increases or decreases monotonically
!      then TRUE is returned, otherwise FALSE is returned. 
!   </OUT>
! determines if the real input array has
! monotonically increasing or decreasing values

function monotonic_array ( array, direction )
real(r8),    intent(in)            :: array(:)
integer, intent(out), optional :: direction
logical :: monotonic_array
integer :: i

! initialize
  monotonic_array = .false.
  if (present(direction)) direction = 0

! array too short
  if ( size(array) < 2 ) return

! ascending
  if ( array(1) < array(size(array)) ) then
     do i = 2, size(array)
       if (array(i-1) < array(i)) cycle
       return
     enddo
     monotonic_array = .true.
     if (present(direction)) direction = +1

! descending
  else
     do i = 2, size(array)
       if (array(i-1) > array(i)) cycle
       return
     enddo
     monotonic_array = .true.
     if (present(direction)) direction = -1
  endif

end function monotonic_array
! </FUNCTION>



!---------------------------------------------------------------
    function stdlog()
      integer :: stdlog
      logical :: opened
          inquire( file=trim(configfile), opened=opened )
          if( opened )then
              !!!call FLUSH(log_unit)
          else
              log_unit=get_unit()
!!!              open( unit=log_unit, status='OLD', file=trim(configfile), position='APPEND', err=10 )
              open( unit=log_unit, file=trim(configfile), position='APPEND', err=10 )
          end if
          stdlog = log_unit 
      return
   10 call error_mesg( 'fms module', 'STDLOG: unable to open '//trim(configfile)//'.' , FATAL)
    end function stdlog

!#######################################################################

end module fms_mod
! <INFO>
!   <BUG>              
!     Namelist error checking may not work correctly with some compilers.
!
!     Users should beware when mixing Fortran reads and read_data calls. If a
!     Fortran read follows read_data and namelist variable read_all_pe = FALSE
!     (not the default), then the code will fail. It is safest if Fortran reads 
!     precede calls to read_data.
!   </BUG>
!   <ERROR MSG="unexpected EOF" STATUS="FATAL">
!     An unexpected end-of-file was encountered in a read_data call.
!     You may want to use the optional end argument to detect the EOF. 
!   </ERROR>
!   <NOTE>
!     1) If the <B>MPP</B> or <B>MPP_DOMAINS</B> stack size is exceeded the
!     program will terminate after printing the required size. 
!   
!     2) When running on a very small number of processors or for high
!     resolution models the default domains_stack_size will
!     probably be insufficient. 
!   </NOTE>
!   <FUTURE>           
!     NetCDF facilities for reading and writing restart files and (IEEE32) 
!       data files.
!    </FUTURE>
!    <FUTURE>
!     May possible split the FMS module into two modules. 
!
!      i.general utilities (FMS_MOD) <BR/>
!     ii.I/O utilities (FMS_IO_MOD) 
!    </FUTURE>
! </INFO>
