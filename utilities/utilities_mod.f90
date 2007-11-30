! Data Assimilation Research Testbed -- DART
! Copyright 2004-2007, Data Assimilation Research Section
! University Corporation for Atmospheric Research
! Licensed under the GPL -- www.gpl.org/licenses/gpl.html

module utilities_mod

! <next few lines under version control, do not edit>
! $URL$
! $Id$
! $Revision$
! $Date$

!-----------------------------------------------------------------------
!
!   A collection of simple useful programs.
!
!      file_exist       Function that returns if a given
!                       file name exists
!
!      get_unit         Function that returns an available 
!                       Fortran unit number
!
!      error_handler    Print warning and error messages, 
!                       terminates program for error messages.
!
!      open_file        Opens a given file name for i/o and returns
!                       a unit number.  If the file is already open
!                       the unit number is returned.
!
!      close_file       Closes the given unit_number. If the file is 
!                       already closed, nothing happens.
!
!      print_version_number    Prints out a routine name and
!                              version number to a specified unit
!
!      set_output       Set the status of printing.  Can be set on a per-task
!                       basis if you are running with multiple tasks.
!                       If set to false only warnings and fatal errors will 
!                       write to the log.
!
!      do_output        Logical function which returns whether informational
!                       messages should be output.  Controlled by the setting
!                       made from set_output.  Useful for messages which cannot
!                       go through the normal error handler (e.g. namelists).
!
!      set_tasknum      Only called for an MPI job with multiple tasks.
!                       Sets the 'multi-task' flag and records the local task
!                       number for later error and info messages.
!
!      nc_check         Check netcdf return codes, and if not ok, extract
!                       the netcdf error string and pass that to the error
!                       handler routine.  Two optional strings allow the caller
!                       to provide the subroutine name and some context.
!
!      write_time       Writes a timestamp in a standard format.
!
! nsc start 31jan07
!   idea - add some unit number routine here?
!   you can extract the filename associated with a fortran unit number
!   with the inquire function on the unit.  this seems useful for
!   automatically generating filenames in messages.  here is an example
!   of how that code works:
!
!character(len=128) :: filename
!logical :: is_named
!integer :: rc
!
!inquire(ncFileID, named=is_named, name=filename, iostat=rc)
!print *, 'is_named =', is_named, 'name = ', trim(filename)
!if ((rc /= 0) .or. (.not. is_named)) filename = 'unknown file'
!
! nsc end 31jan07
!
!-----------------------------------------------------------------------

use types_mod, only : r8
use netcdf

implicit none
private

!   ---- private data for check_nml_error ----

integer, private :: num_nml_error_codes, nml_error_codes(5)
logical, private :: do_output_flag = .true.
logical, private :: single_task = .true.
integer, private :: task_number = 0

integer, parameter :: E_DBG = -1,   E_MSG = 0,  E_WARN = 1, E_ERR = 2
integer, parameter :: DEBUG = -1, MESSAGE = 0, WARNING = 1, FATAL = 2

public :: file_exist, get_unit, open_file, timestamp, set_tasknum, &
       close_file, register_module, error_handler, logfileunit, nc_check, &
       initialize_utilities, finalize_utilities, dump_unit_attributes, &
       find_namelist_in_file, check_namelist_read, do_output, set_output, &
       E_DBG, E_MSG, E_WARN, E_ERR, & 
       DEBUG, MESSAGE, WARNING, FATAL

! version controlled file description for error handling, do not edit
character(len=128), parameter :: &
   source   = "$URL$", &
   revision = "$Revision$", &
   revdate  = "$Date$"

logical, save :: module_initialized = .false.
integer, save :: logfileunit = -1

!----------------------------------------------------------------
! Namelist input with default values
integer  :: TERMLEVEL = E_ERR     ! E_ERR All warnings/errors are assumed fatal.
character(len=129) :: logfilename = 'dart_log.out'
logical  :: module_details = .true.  ! print svn details about each module

namelist /utilities_nml/ TERMLEVEL, logfilename, module_details

contains

!#######################################################################

   subroutine initialize_utilities(progname, alternatename, output_flag)
   character(len=*), intent(in), optional :: progname
   character(len=*), intent(in), optional :: alternatename
   logical, intent(in), optional          :: output_flag
   ! integer :: logfileunit -- public module variable
   integer :: iunit, io

   character(len=129) :: lname


      if ( module_initialized ) then ! nothing to do

         ! write(*,*)'Module initialized ... carry on.'

         return

      else ! initialize the module
         
         module_initialized = .true.

         if (present(output_flag)) do_output_flag = output_flag

         ! Since the logfile is not open yet, the error terminations
         ! must be handled differently than all other cases.
         ! The routines that normally write to the logfile cannot
         ! be used just yet. If we cannot open a logfile, we
         ! always abort execution at this step.

         if ( present(progname) ) then
            if (do_output_flag) write(*,*)'Starting program ',trim(adjustl(progname))
         endif

         if (do_output_flag) write(*,*)'Initializing the utilities module.'

         ! Read the namelist entry
         call find_namelist_in_file("input.nml", "utilities_nml", iunit, .false.)
         read(iunit, nml = utilities_nml, iostat = io)
         call check_namelist_read(iunit, io, "utilities_nml", .false.)

         ! Open the log file with the name from the namelist 
         logfileunit = nextunit()
         if ( logfileunit < 0 ) then
            write(*,*)'   unable to get a unit to use for the logfile.'
            write(*,*)'   stopping.'
            stop 99
         endif

         if (present(alternatename)) then
            lname = alternatename
         else
            lname = logfilename
         endif

         if (do_output_flag) write(*,*)'Trying to log to unit ', logfileunit
         if (do_output_flag) write(*,*)'Trying to open file ', trim(adjustl(lname))

         open(logfileunit, file=trim(adjustl(lname)), form='formatted', &
                           position='append', iostat = io )
         if ( io /= 0 ) then
            write(*,*)'FATAL ERROR in initialize_utilities'
            write(*,*)'  ',trim(source)
            write(*,*)'  ',trim(revision)
            write(*,*)'  ',trim(revdate)
            write(*,*)'   unable to open the logfile.'
            write(*,*)'   the intended file name was <',trim(lname),'>'
            write(*,*)'   stopping.'
            stop 99
         endif

         ! Log the run-time 

         if (do_output_flag) then
         if ( present(progname) ) then
               call write_time (logfileunit, label='Starting ', &
                                string1='Program '//trim(progname))
               call write_time (             label='Starting ', &
                                string1='Program '//trim(progname))
         else
               call write_time (logfileunit, label='Starting ')
               call write_time (             label='Starting ')
         endif 
         endif

         ! Check to make sure termlevel is set to a reasonable value
         call checkTermLevel

         ! Echo the module information using normal mechanism
         call register_module(source, revision, revdate)

         ! Echo the namelist values for this module using normal mechanism
         if (do_output_flag) write(logfileunit, nml=utilities_nml)
         if (do_output_flag) write(     *     , nml=utilities_nml)

      endif

   contains

      function nextunit() result(iunit)
         integer :: iunit

         logical :: open
         integer :: i

         iunit = -1
         UnitLoop : do i = 10, 80
            inquire (i, opened=open)
            if (.not. open) then
               iunit = i
               exit UnitLoop
            endif
         enddo UnitLoop
         if ( iunit < 0 ) then 
            write(*,*)'FATAL ERROR in initialize_utilities'
            write(*,*)'  ',trim(source)
            write(*,*)'  ',trim(revision)
            write(*,*)'  ',trim(revdate)
         endif
      end function nextunit

      subroutine checktermlevel
         select case (TERMLEVEL)
             case (E_MSG)
                ! do nothing
             case (E_WARN)
                ! do nothing
             case (E_ERR)
                ! do nothing
             case default
                print *, ' MESSAGE from initialize_utilities'
                print *, ' namelist input of TERMLEVEL is ',TERMLEVEL
                print *, ' possible values are ',E_MSG, E_WARN, E_ERR
                if (TERMLEVEL < E_WARN ) TERMLEVEL = E_WARN
                if (TERMLEVEL > E_ERR  ) TERMLEVEL = E_ERR
                print *, ' using ',TERMLEVEL
         end select
      end subroutine checktermlevel

   end subroutine initialize_utilities


   subroutine finalize_utilities
   ! integer :: logfileunit -- private module variable

      close(logfileunit)

   end subroutine finalize_utilities


!#######################################################################

   subroutine register_module(src, rev, rdate)
   character(len=*), intent(in) :: src, rev, rdate

      if ( .not. module_initialized ) call initialize_utilities
      if ( .not. do_output_flag) return
      if ( .not. module_details) return


      write(logfileunit,*)
      write(logfileunit,*)'Registering module :'
      write(logfileunit,*)trim(src)
      write(logfileunit,*)trim(rev)
      write(logfileunit,*)trim(rdate)
      write(logfileunit,*)'Registration complete.'
      write(logfileunit,*)

      write(     *     ,*)
      write(     *     ,*)'Registering module :'
      write(     *     ,*)trim(src)
      write(     *     ,*)trim(rev)
      write(     *     ,*)trim(rdate)
      write(     *     ,*)'Registration complete.'
      write(     *     ,*)

   end subroutine register_module

!#######################################################################

   subroutine timestamp(string1,string2,string3,pos)
      ! values(1) year
      ! values(2) month
      ! values(3) day
      ! values(4) minutes diff from UTC
      ! values(5) hour
      ! values(6) minutes
      ! values(7) seconds
      ! values(8) milliseconds

      character(len=*), optional, intent(in) :: string1
      character(len=*), optional, intent(in) :: string2
      character(len=*), optional, intent(in) :: string3
      character(len=*),           intent(in) :: pos

      integer, dimension(8) :: values

      if ( .not. module_initialized ) call initialize_utilities
      if ( .not. do_output_flag) return

      if (trim(adjustl(pos)) == 'end') then
         call write_time (logfileunit, label='Finished ', &
                          string1=string1, string2=string2, string3=string3)
         call write_time (             label='Finished ', &
                          string1=string1, string2=string2, string3=string3)

         call finalize_utilities()
      else
         call write_time (logfileunit, & 
                          string1=string1, string2=string2, string3=string3)
         call write_time (string1=string1, string2=string2, string3=string3)
          
      endif

   end subroutine timestamp

!#######################################################################

   function file_exist (file_name)

      character(len=*), intent(in) :: file_name
      logical  file_exist

      if ( .not. module_initialized ) call initialize_utilities

      inquire (file=file_name(1:len_trim(file_name)), exist=file_exist)

   end function file_exist

!#######################################################################

   function get_unit () result (iunit)

      integer :: i, iunit
      logical :: open

      if ( .not. module_initialized ) call initialize_utilities

! ---- get available unit ----

      iunit = -1
      do i = 10, 80
         inquire (i, opened=open)
         if (.not.open) Then
            iunit = i
            exit
         endif
      enddo

      if (iunit == -1) call error_handler(E_ERR,'get_unit', &
             'No available units.', source, revision, revdate)

   end function get_unit

!#######################################################################

   subroutine dump_unit_attributes(iunit) 
!--------------------------------------------------------------------------------
!  subroutine dump_unit_attributes(iunit) 
!
!  Useful for dumping all the attributes for a file 'unit'
!  A debugging routine, really. TJH Oct 2004

      integer, intent(in) :: iunit

      logical :: exists, connected, named_file
      character(len=129) :: file_name, srname, str1
      character(len=32)  :: ynu     ! YES, NO, UNDEFINED ... among others
      integer :: ios, reclen, nextrecnum

      if ( .not. module_initialized ) call initialize_utilities

      srname = "dump_unit_attributes"

! --- start querying 

      write(str1,*)'for unit ',iunit 
      call error_handler(E_MSG, srname, str1, source, revision, revdate)

      inquire(iunit, opened = connected, iostat=ios)
      if ( connected .and. (ios == 0) ) &
         call error_handler(E_MSG, srname, ' connected', source, revision, revdate)

      inquire(iunit, named = named_file, iostat=ios)
      if ( named_file .and. (ios == 0) ) &
         call error_handler(E_MSG, srname, ' file is named.', source, revision, revdate)

      inquire(iunit, name = file_name, iostat=ios)
      if ( ios == 0 ) then
         write(str1,*)'file name is ' // trim(adjustl(file_name))
         call error_handler(E_MSG, srname, str1, source, revision, revdate)
      endif

      inquire(iunit, exist = exists, iostat=ios)
      if ( exists .and. (ios == 0) ) &
         call error_handler(E_MSG, srname, ' file exists', source, revision, revdate)

      inquire(iunit, recl = reclen, iostat=ios)
      if ( ios == 0 ) then
         write(str1,*)'record length is ', reclen
         call error_handler(E_MSG, srname, str1, source, revision, revdate)
      endif

      inquire(iunit, nextrec = nextrecnum, iostat=ios)
      if ( ios == 0 ) then
         write(str1,*)'next record is ', nextrecnum
         call error_handler(E_MSG, srname, str1, source, revision, revdate)
      endif

      inquire(iunit, access = ynu, iostat=ios)
      if ( ios == 0 ) then
         write(str1,*)'access_type is ', ynu
         call error_handler(E_MSG, srname, str1, source, revision, revdate)
      endif

      inquire(iunit, sequential = ynu, iostat=ios)
      if ( ios == 0 ) then
         write(str1,*)'is file sequential ', ynu
         call error_handler(E_MSG, srname, str1, source, revision, revdate)
      endif

      inquire(iunit, direct = ynu, iostat=ios)
      if ( ios == 0 ) then
         write(str1,*)'is file direct ', ynu
         call error_handler(E_MSG, srname, str1, source, revision, revdate)
      endif

      inquire(iunit, form = ynu, iostat=ios)
      if ( ios == 0 ) then
         write(str1,*)'file format ', ynu
         call error_handler(E_MSG, srname, str1, source, revision, revdate)
      endif

      inquire(iunit, action = ynu, iostat=ios)
      if ( ios == 0 ) then
         write(str1,*)'action ', ynu
         call error_handler(E_MSG, srname, str1, source, revision, revdate)
      endif

      inquire(iunit, read = ynu, iostat=ios)
      if ( ios == 0 ) then
         write(str1,*)'read ', ynu
         call error_handler(E_MSG, srname, str1, source, revision, revdate)
      endif

      inquire(iunit, write = ynu, iostat=ios)
      if ( ios == 0 ) then
         write(str1,*)'write ', ynu
         call error_handler(E_MSG, srname, str1, source, revision, revdate)
      endif

      inquire(iunit, readwrite = ynu, iostat=ios)
      if ( ios == 0 ) then
         write(str1,*)'readwrite ', ynu
         call error_handler(E_MSG, srname, str1, source, revision, revdate)
      endif

      inquire(iunit, blank = ynu, iostat=ios)
      if ( ios == 0 ) then
         write(str1,*)'blank ', ynu
         call error_handler(E_MSG, srname, str1, source, revision, revdate)
      endif

      inquire(iunit, position = ynu, iostat=ios)
      if ( ios == 0 ) then
         write(str1,*)'position ', ynu
         call error_handler(E_MSG, srname, str1, source, revision, revdate)
      endif

      inquire(iunit, delim = ynu, iostat=ios)
      if ( ios == 0 ) then
         write(str1,*)'delim ', ynu
         call error_handler(E_MSG, srname, str1, source, revision, revdate)
      endif

      inquire(iunit, pad = ynu, iostat=ios)
      if ( ios == 0 ) then
         write(str1,*)'pad ', ynu
         call error_handler(E_MSG, srname, str1, source, revision, revdate)
      endif

   end subroutine dump_unit_attributes

!#######################################################################

   subroutine error_mesg (routine, message, level)

!             ------------------------------------
!             |                                  |
!             |    a very simple error handler   |
!             |                                  |
!             ------------------------------------
!
!  input:
!      routine   name of the calling routine (character string)
!      message   message written to standard output (character string)
!      level     if not equal to zero then the program terminates
!
          character(len=*), intent(in) :: routine, message
          integer,          intent(in) :: level

          select case (iabs(level))
             case (0)
                if ( .not. do_output_flag) return
                print *, ' MESSAGE from ',routine(1:len_trim(routine))
                print *, ' ',message(1:len_trim(message))
             case (1)
                print *, ' WARNING in ',routine(1:len_trim(routine))
                print *, ' ',message(1:len_trim(message))
             case default
                print *, ' ERROR in ',routine(1:len_trim(routine))
                print *, ' ',message(1:len_trim(message))
                stop 99
          end select

!         --------------------------------------------

   end subroutine error_mesg

!#######################################################################

  subroutine error_handler(level, routine, text, src, rev, rdate, aut )
!----------------------------------------------------------------------
! subroutine error_handler(level, routine, text, src, rev, rdate, aut )
!
! logs warning/error 
implicit none

integer, intent(in) :: level
character(len = *), intent(in) :: routine, text, src, rev, rdate
character(len = *), intent(in), OPTIONAL :: aut

character(len = 8) :: taskstr

if ( .not. module_initialized ) call initialize_utilities

! current choice is to log all errors and warnings regardless
! of setting of output flag.  messages only print from those
! tasks which are allowed to print.

select case(level)
   case (E_MSG)

      if ( .not. do_output_flag) return
      if ( single_task ) then
        write(     *     , *) trim(routine),' ', trim(text)
        write(logfileunit, *) trim(routine),' ', trim(text)
      else
        if (task_number < 10) then
            write(taskstr, '(a,i1)' ) "PE ", task_number
        else if (task_number < 100) then
            write(taskstr, '(a,i2)' ) "PE ", task_number
        else
            write(taskstr, '(a,i5)' ) "PE ", task_number
        endif
        write(     *     , *) trim(taskstr),': ',trim(routine),' ', trim(text)
        write(logfileunit, *) trim(taskstr),': ',trim(routine),' ', trim(text)
      endif

   case (E_WARN)

      write(     *     , *) 'WARNING FROM:'
      if ( .not. single_task ) &
      write(     *     , *) '   task id       : ', task_number
      write(     *     , *) '   routine       : ', trim(routine)
      write(     *     , *) '   source file   : ', trim(src)
      write(     *     , *) '   file revision : ', trim(rev)
      write(     *     , *) '   revision date : ', trim(rdate)
      if(present(aut)) &
      write(     *     , *) '   last editor   : ', trim(aut)
      write(     *     , *) '   message       : ', trim(text)

      write(logfileunit, *) 'WARNING FROM:'
      if ( .not. single_task ) &
      write(logfileunit, *) '   task id       : ', task_number
      write(logfileunit, *) '   routine       : ', trim(routine)
      write(logfileunit, *) '   source file   : ', trim(src)
      write(logfileunit, *) '   file revision : ', trim(rev)
      write(logfileunit, *) '   revision date : ', trim(rdate)
      if(present(aut)) &
      write(logfileunit, *) '   last editor   : ', trim(aut)
      write(logfileunit, *) '   message       : ', trim(text)

   case(E_ERR)

      write(     *     , *) 'ERROR FROM:'
      if ( .not. single_task ) &
      write(     *     , *) '   task id       : ', task_number
      write(     *     , *) '   routine       : ', trim(routine)
      write(     *     , *) '   source file   : ', trim(src)
      write(     *     , *) '   file revision : ', trim(rev)
      write(     *     , *) '   revision date : ', trim(rdate)
      if(present(aut)) &
      write(     *     , *) '   last editor   : ', trim(aut)
      write(     *     , *) '   message       : ', trim(text)

      write(logfileunit, *) 'ERROR FROM:'
      if ( .not. single_task ) &
      write(logfileunit, *) '   task id       : ', task_number
      write(logfileunit, *) '   routine       : ', trim(routine)
      write(logfileunit, *) '   source file   : ', trim(src)
      write(logfileunit, *) '   file revision : ', trim(rev)
      write(logfileunit, *) '   revision date : ', trim(rdate)
      if(present(aut)) &
      write(logfileunit, *) '   last editor   : ', trim(aut)
      write(logfileunit, *) '   message       : ', trim(text)

end select

! TERMLEVEL gets set in the namelist
if( level >= TERMLEVEL ) call exit( 99 ) 

end subroutine error_handler

!#######################################################################

   function open_file (fname, form, action) result (iunit)

   character(len=*), intent(in) :: fname
   character(len=*), intent(in), optional :: form, action
   integer  :: iunit

   integer           :: nc
   logical           :: open
   character(len=11) :: format
   character(len=6)  :: location

   if ( .not. module_initialized ) call initialize_utilities

      inquire (file=fname(1:len_trim(fname)), opened=open, number=iunit,  &
               form=format)

     if (open) then
! ---------- check format ??? ---------
! ---- (skip this and let fortran i/o catch bug) -----

    !    if (present(form)) then
    !        nc = min(11,len(form))
    !        if (format == 'UNFORMATTED') then
    !             if (form(1:nc) /= 'unformatted' .and.  &
    !                 form(1:nc) /= 'UNFORMATTED')       &
    !                 call error_mesg ('open_file in utilities_mod', &
    !                                  'invalid form argument', 2)
    !        else if (format(1:9) == 'FORMATTED') then
    !             if (form(1:nc) /= 'formatted' .and.  &
    !                 form(1:nc) /= 'FORMATTED')       &
    !                 call error_mesg ('open_file in utilities_mod', &
    !                                  'invalid form argument', 2)
    !        else
    !             call error_mesg ('open_file in utilities_mod', &
    !                       'unexpected format returned by inquire', 2)
    !        endif
    !    endif
         
     else
! ---------- open file ----------

         format   = 'formatted  '
         location = 'asis  '

         if (present(form)) then
             nc = min(11,len(form))
             format(1:nc) = form(1:nc)
         endif

         if (present(action)) then
             nc = min(6,len(action))
             location(1:nc) = action(1:nc)
             if(location /= 'append' .and. location /= 'APPEND') &
                location = 'rewind'
         endif

         iunit = get_unit()

         if (format == 'formatted  ' .or. format == 'FORMATTED  ') then
             open (iunit, file=fname(1:len_trim(fname)),      &
                         form=format(1:len_trim(format)),  &
                     position=location(1:len_trim(location)),  &
                        delim='apostrophe')
         else
             open (iunit, file=fname(1:len_trim(fname)),      &
                         form=format(1:len_trim(format)),  &
                     position=location(1:len_trim(location)) )
         endif
     endif


   end function open_file

!#######################################################################

   subroutine print_version_number (iunit, routine, version)

! *** prints routine name and version number to a log file ***
!
!    in:  iunit    = unit number to direct output
!         routine = routine name (character, max len = 20)
!         version = version name or number (character, max len = 8)

   integer,          intent(in) :: iunit
   character(len=*), intent(in) :: routine, version

   integer           :: n
   character(len=20) :: name
   character(len=8)  :: vers

   if ( .not. module_initialized ) call initialize_utilities
   if ( .not. do_output_flag) return

     n = min(len(routine),20); name = adjustl(routine(1:n))
     n = min(len(version), 8); vers = adjustl(version(1:n))

     if (iunit > 0) then
         write (iunit,10) name, vers
     else
         write (*,10) name, vers
     endif

  10 format (/,60('-'),  &
             /,10x, 'ROUTINE = ',a20, '  VERSION = ', a8, &
             /,60('-'))

! 10 format (/,1x, 12('>'), 1x, 'ROUTINE = ',a20, '  VERSION = ', a8, &
!              1x, 12('<'),/)

   end subroutine print_version_number

!#######################################################################

   subroutine write_time (unit, label, string1, string2, string3, tz)

! ***  Write the current time to a log file or standard output ***
!
!    in: unit number (default is * if not specified)
!    in: label (default is  "Time is" if not specified)
!    in: string1,2,3 (no defaults)

   integer,          optional, intent(in) :: unit
   character(len=*), optional, intent(in) :: label
   character(len=*), optional, intent(in) :: string1
   character(len=*), optional, intent(in) :: string2
   character(len=*), optional, intent(in) :: string3
   logical,          optional, intent(in) :: tz


   integer :: lunit
   character(len= 8) :: cdate
   character(len=10) :: ctime
   character(len= 5) :: zone
   integer, dimension(8) :: values

   if (present(unit)) then
      lunit = unit
   else
      lunit = 6   ! this should be *
   endif

   call DATE_AND_TIME(cdate, ctime, zone, values)

   ! give up if no good values were returned
   if (.not. any(values /= -HUGE(0)) ) return 

   ! ok to proceed
   write(lunit,*)
   write(lunit,*)'--------------------------------------'
   if ( present(label) ) then
      write(lunit,*) label // '... at YYYY MM DD HH MM SS = '
   else
      write(lunit,*) 'Time is  ... at YYYY MM DD HH MM SS = '
   endif 
   write(lunit,'(17x,i4,5(1x,i2))') values(1), values(2), &
                     values(3),  values(5), values(6), values(7)

   if(present(string1)) write(lunit,*)trim(string1)
   if(present(string2)) write(lunit,*)trim(string2)
   if(present(string3)) write(lunit,*)trim(string3)

   if (present(tz)) then
      if ( values(4) /= -HUGE(0) .and. tz) &
         write(lunit,*)'time zone offset is ',values(4),' minutes.'
   endif

   write(lunit,*)'--------------------------------------'
   write(lunit,*)

   end subroutine write_time

!#######################################################################

   subroutine set_output (doflag)

! *** set whether output is written to a log file or simply ignored ***
!
!    in:  doflag  = whether to output log information or not

   logical, intent(in) :: doflag

!! THIS ONE IS DIFFERENT.  Set the flag FIRST before doing the
!! standard initialization, so if you are turning off writing
!! for some tasks you do not get output you are trying to avoid.

   do_output_flag = doflag

   if ( .not. module_initialized ) call initialize_utilities

   end subroutine set_output


!#######################################################################

   function do_output ()

! *** return whether output should be written from this task ***
!

   logical :: do_output

   if ( .not. module_initialized ) call initialize_utilities

   do_output = do_output_flag

   end function do_output

!#######################################################################


   subroutine set_tasknum (tasknum)

! *** for multiple-task jobs, set the task number for error msgs ***
!
!    in:  tasknum  = task number, 0 to N-1

   integer, intent(in) :: tasknum

   if ( .not. module_initialized ) call initialize_utilities

   single_task = .false. 
   task_number = tasknum

   end subroutine set_tasknum


!#######################################################################


subroutine close_file(iunit)
!-----------------------------------------------------------------------
!
! Closes the given unit_number. If the file is already closed, 
! nothing happens. Pretty dramatic, eh?
!

integer, intent(in) :: iunit

integer :: ios
logical :: open

if ( .not. module_initialized ) call initialize_utilities

inquire (unit=iunit, opened=open, iostat=ios)
if ( ios /= 0 ) then
   print *,'Dagnabbit. Cannot inquire about unit # ',iunit
   print *,'Error status was ',ios
   print *,'Hoping for the best and continuing.'
endif

if (open) close(iunit)

end subroutine close_file

!#######################################################################



subroutine find_namelist_in_file(namelist_file_name, nml_name, iunit, &
   write_to_logfile_in)
!-----------------------------------------------------------------------
! 
! Opens namelist_file_name if it exists on unit iunit, error if it
! doesn't exist.
! Searches file for a line containing ONLY the string
! &nml_name, for instance &filter_nml. If found, rewinds the file and
! returns true. Otherwise, error message and terminates
!

character(len = *), intent(in) :: namelist_file_name
character(len = *), intent(in) :: nml_name
integer, intent(out)           :: iunit
logical, intent(in), optional :: write_to_logfile_in

character(len = 169) :: nml_string, test_string, err_string
integer              :: io
logical              :: write_to_logfile

! Decide if there is a logfile or not
write_to_logfile = .true.
if(present(write_to_logfile_in)) write_to_logfile = write_to_logfile_in

! Check for file existence; no file is an error
if(file_exist(trim(namelist_file_name))) then
   ! Open the file
   iunit = open_file(trim(namelist_file_name), action = 'read')

   ! Now look for the start of a namelist with &nml_name
   test_string = '&' // trim(adjustl(nml_name))

   ! Read each line until end of file is found
   do
      read(iunit, '(A)', iostat = io) nml_string
      if(io /= 0) then
         ! No values for this namelist; error
         write(err_string, *) 'Namelist entry &', nml_name, ' must exist in ', namelist_file_name
         ! Can't write to logfile if it hasn't yet been opened
         if(write_to_logfile) then
            call error_handler(E_ERR, 'find_namelist_in_file', err_string, &
               source, revision, revdate)
         else
            write(*, *) 'FATAL ERROR before logfile initialization in utilities_mod'
            write(*, *) 'Error is in subroutine find_namelist_in_file'
            write(*, *) err_string
            write(*,*)'  ',trim(source)
            write(*,*)'  ',trim(revision)
            write(*,*)'  ',trim(revdate)
            stop 
         endif
      else
         if(trim(adjustl(nml_string)) == trim(adjustl(test_string))) then
            rewind(iunit)
            return
         endif
      endif
   end do
else
   ! No namelist_file_name file is an error
   write(err_string, *) 'Namelist input file: ', namelist_file_name, ' must exist.'
   if(write_to_logfile) then
      call error_handler(E_ERR, 'find_namelist_in_file', err_string, &
         source, revision, revdate)
   else
      write(*, *) 'FATAL ERROR before logfile initialization in utilities_mod'
      write(*, *) 'Error is in subroutine find_namelist_in_file'
      write(*, *) err_string
      write(*,*)'  ',trim(source)
      write(*,*)'  ',trim(revision)
      write(*,*)'  ',trim(revdate)
      stop 
   endif
endif

end subroutine find_namelist_in_file


!#######################################################################



subroutine check_namelist_read(iunit, iostat_in, nml_name, &
   write_to_logfile_in)
!-----------------------------------------------------------------------
! 
! Confirms that a namelist read was successful. If it failed
! produces an error message and stops execution.
!

integer,            intent(in) :: iunit, iostat_in
character(len = *), intent(in) :: nml_name
logical, intent(in), optional :: write_to_logfile_in

character(len=159) :: nml_string, err_string
integer            :: io
logical            :: write_to_logfile

! Decide if there is a logfile or not
write_to_logfile = .true.
if(present(write_to_logfile_in)) write_to_logfile = write_to_logfile_in

if(iostat_in == 0) then
   ! If the namelist read was successful, just close the file
   call close_file(iunit)
else
   ! If it wasn't successful, print the line on which it failed  
   BACKSPACE iunit
   read(iunit, '(A)', iostat = io) nml_string
   ! A failure in this read means that the namelist started but never terminated
   ! Result was falling off the end, so backspace followed by read fails
   if(io /= 0) then
      write(err_string, *) 'Namelist ', trim(nml_name), ' started but never terminated'
      if(write_to_logfile) then
         call error_handler(E_ERR, 'check_namelist_read', err_string, &
            source, revision, revdate)
      else
         write(*, *) 'FATAL ERROR before logfile initialization in utilities_mod'
         write(*, *) 'Error is in subroutine check_namelist_read'
         write(*, *) err_string
         write(*,*)'  ',trim(source)
         write(*,*)'  ',trim(revision)
         write(*,*)'  ',trim(revdate)
         stop 
      endif
   else
      ! Didn't fall off end so bad entry in the middle of namelist
      ! TEMP HELP FOR USERS; remove after next release
      if ((nml_name(1:10) == 'filter_nml') .and. (index(nml_string,'inf_start_from_restart') > 0)) then
         write(err_string, *) 'inf_start_from_restart obsolete'
         call error_handler(E_MSG, 'filter_nml: ', err_string, "", "", "")
         write(err_string, *) 'use inf_initial_from_restart and inf_sd_initial_from_restart'
         call error_handler(E_MSG, 'filter_nml: ', err_string, "", "", "")
      endif 
      write(err_string, *) 'INVALID NAMELIST ENTRY: ', trim(nml_string), ' in namelist ', trim(nml_name)
      if(write_to_logfile) then
         call error_handler(E_ERR, 'check_namelist_read', err_string, &
            source, revision, revdate)
      else
         write(*, *) 'FATAL ERROR before logfile initialization in utilities_mod'
         write(*, *) 'Error is in subroutine check_namelist_read'
         write(*, *) err_string
         write(*,*)'  ',trim(source)
         write(*,*)'  ',trim(revision)
         write(*,*)'  ',trim(revdate)
         stop 
      endif
   endif
endif

end subroutine check_namelist_read

!#######################################################################

   subroutine nc_check(istatus, subr_name, context)
      integer, intent (in)                   :: istatus
      character(len=*), intent(in)           :: subr_name
      character(len=*), intent(in), optional :: context
  
      character(len=129) :: error_msg
  
      ! if no error, nothing to do here.  we are done.
      if( istatus == nf90_noerr) return


      ! something wrong.  construct an error string and call the handler.

      ! context is optional, but is very useful if specified.
      ! if context + error code > 129, the assignment will truncate.
      if (present(context) ) then
          error_msg = trim(context) // ': ' // trim(nf90_strerror(istatus))
      else
          error_msg = nf90_strerror(istatus)
      endif

      ! this does not return 
      call error_mesg(subr_name, error_msg, FATAL)
  

   end subroutine nc_check

!=======================================================================
! End of utilities_mod
!=======================================================================

end module utilities_mod

