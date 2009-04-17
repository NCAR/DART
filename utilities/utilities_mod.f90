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
!   A collection of simple useful routines:
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
!      logfileunit      Global integer unit numbers for the log file and
!      nmlfileunit      for the namelist file (which defaults to the same as log)
!
!      to_upper         converts a character string to uppercase
!
!      find_textfile_dims    finds number of lines and max line length in a 
!                            text file. Used so we can record the namelist files 
!                            in the netcdf output files.
!
!      file_to_text     converts the contents of a (hopefully small) file to
!                       a single text variable ... to record in the
!                       netcdf output files.
!
!      is_longitude_between  checks whether a given longitude is between
!                            the two given limits, starting at the first and
!                            going EAST until reaching the second.  the end
!                            points are included. if min=max, all points are
!                            considered inside.  there is no rejection of the
!                            input values based on range; they are all converted
!                            to [0-360) by calling modulo() before starting.
!                            default is degrees, but there is an optional
!                            argument to select radians instead.
!                           
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
!inquire(unitnum, named=is_named, name=filename, iostat=rc)
!print *, 'is_named =', is_named, 'name = ', trim(filename)
!if ((rc /= 0) .or. (.not. is_named)) filename = 'unknown file'
!
! nsc end 31jan07
!
!-----------------------------------------------------------------------

use types_mod, only : r8, PI
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

real(r8), parameter :: TWOPI = PI * 2.0_r8

public :: file_exist, get_unit, open_file, close_file, timestamp, &
       register_module, error_handler, to_upper, next_file, &
       nc_check, logfileunit, nmlfileunit, &
       find_textfile_dims, file_to_text, is_longitude_between, &
       initialize_utilities, finalize_utilities, dump_unit_attributes, &
       find_namelist_in_file, check_namelist_read, &
       set_tasknum, set_output, do_output,  &
       E_DBG, E_MSG, E_WARN, E_ERR, & 
       DEBUG, MESSAGE, WARNING, FATAL

! this routine is either in the null_mpi_utilities_mod.f90, or in
! the mpi_utilities_mod.f90 file, but it is not a module subroutine.
! the mpi files use this module, and you cannot have circular module
! references.  the point of this call is that in the mpi multi-task
! case, you want to call MPI_Abort() to kill the other tasks associated
! with this job when you exit.  in the non-mpi case, it just calls exit.
interface
 subroutine exit_all(exitval)
  integer :: exitval
 end subroutine exit_all
end interface

! version controlled file description for error handling, do not edit
character(len=128), parameter :: &
   source   = "$URL$", &
   revision = "$Revision$", &
   revdate  = "$Date$"

logical, save :: module_initialized = .false.
integer, save :: logfileunit = -1
integer, save :: nmlfileunit = -1


character(len = 169) :: msgstring

!----------------------------------------------------------------
! Namelist input with default values
integer  :: TERMLEVEL = E_ERR     ! E_ERR All warnings/errors are assumed fatal.
character(len=129) :: logfilename = 'dart_log.out'
character(len=129) :: nmlfilename = 'dart_log.nml'
logical  :: module_details = .true.  ! print svn details about each module
logical  :: print_debug    = .false. ! print messages labeled DBG

namelist /utilities_nml/ TERMLEVEL, logfilename, module_details, &
                         nmlfilename, print_debug

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
            call exit_all(77)
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
            call exit_all(77)
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

         ! If nmlfilename != logfilename, open it.  otherwise set nmlfileunit
         ! to be same as logunit.
         if (trim(adjustl(nmlfilename)) /= trim(adjustl(lname))) then
            if (do_output_flag) &
             write(*,*)'Trying to open namelist log ', trim(adjustl(nmlfilename))
    
            nmlfileunit = nextunit()
            if (nmlfileunit < 0) &
               call error_handler(E_ERR,'initialize_utilities', &
                 'Cannot get unit for nm log file', source, revision, revdate)

            open(nmlfileunit, file=trim(adjustl(nmlfilename)), form='formatted', &
                 position='append', iostat = io )
            if ( io /= 0 ) then
               call error_handler(E_ERR,'initialize_utilities', &
                   'Cannot open nm log file', source, revision, revdate)
            endif
    
         else
           nmlfileunit = logfileunit
         endif

         ! Echo the namelist values for this module using normal mechanism
         ! including a separator line for this run.
         if (do_output_flag) then
            if (nmlfileunit /= logfileunit) then
               if ( present(progname) ) then
                  write(nmlfileunit, *) '!Starting Program '//trim(progname)
               else
                  write(nmlfileunit, *) '!Starting Program '
               endif 
            endif
            write(nmlfileunit, nml=utilities_nml)
            write(    *     , nml=utilities_nml)
         endif

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

      if (do_output_flag) then
         if (nmlfileunit /= logfileunit) then
            write(nmlfileunit, *) '!Ending Program '
         endif 
      endif

      close(logfileunit)
      if (nmlfileunit /= logfileunit) close(nmlfileunit)

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

      character(len=*), optional, intent(in) :: string1
      character(len=*), optional, intent(in) :: string2
      character(len=*), optional, intent(in) :: string3
      character(len=*),           intent(in) :: pos

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
                call exit_all(99)
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
character(len = *), intent(in) :: routine, text
character(len = *), intent(in), optional :: src, rev, rdate, aut

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

   case (E_DBG)
      if (print_debug) then

         ! what about do_output_flag?  want messages from all procs or just PE0?
         if ( single_task ) then
           write(     *     , *) 'DEBUG FROM: ', trim(routine),' ', trim(text)
           write(logfileunit, *) 'DEBUG FROM: ', trim(routine),' ', trim(text)
         else
           if (task_number < 10) then
               write(taskstr, '(a,i1)' ) "PE ", task_number
           else if (task_number < 100) then
               write(taskstr, '(a,i2)' ) "PE ", task_number
           else
               write(taskstr, '(a,i5)' ) "PE ", task_number
           endif
           write(     *     , *) trim(taskstr),': DEBUG FROM: ',trim(routine),' ', trim(text)
           write(logfileunit, *) trim(taskstr),': DEBUG FROM: ',trim(routine),' ', trim(text)
         endif
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
if( level >= TERMLEVEL ) call exit_all( 99 ) 

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

character(len = 169) :: nml_string, test_string
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
         write(msgstring, *) 'Namelist entry &', nml_name, ' must exist in ', namelist_file_name
         ! Can't write to logfile if it hasn't yet been opened
         if(write_to_logfile) then
            call error_handler(E_ERR, 'find_namelist_in_file', msgstring, &
               source, revision, revdate)
         else
            write(*, *) 'FATAL ERROR before logfile initialization in utilities_mod'
            write(*, *) 'Error is in subroutine find_namelist_in_file'
            write(*, *) msgstring
            write(*,*)'  ',trim(source)
            write(*,*)'  ',trim(revision)
            write(*,*)'  ',trim(revdate)
            call exit_all(88) 
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
   write(msgstring, *) 'Namelist input file: ', namelist_file_name, ' must exist.'
   if(write_to_logfile) then
      call error_handler(E_ERR, 'find_namelist_in_file', msgstring, &
         source, revision, revdate)
   else
      write(*, *) 'FATAL ERROR before logfile initialization in utilities_mod'
      write(*, *) 'Error is in subroutine find_namelist_in_file'
      write(*, *) msgstring
      write(*,*)'  ',trim(source)
      write(*,*)'  ',trim(revision)
      write(*,*)'  ',trim(revdate)
      call exit_all(88) 
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

character(len=159) :: nml_string
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
      write(msgstring, *) 'Namelist ', trim(nml_name), ' started but never terminated'
      if(write_to_logfile) then
         call error_handler(E_ERR, 'check_namelist_read', msgstring, &
            source, revision, revdate)
      else
         write(*, *) 'FATAL ERROR before logfile initialization in utilities_mod'
         write(*, *) 'Error is in subroutine check_namelist_read'
         write(*, *) msgstring
         write(*,*)'  ',trim(source)
         write(*,*)'  ',trim(revision)
         write(*,*)'  ',trim(revdate)
         call exit_all(66) 
      endif
   else
      ! Didn't fall off end so bad entry in the middle of namelist
      ! TEMP HELP FOR USERS; remove after next release
      if (len(nml_name) >= 10) then
         if ((nml_name(1:10) == 'filter_nml') .and. (index(nml_string,'inf_start_from_restart') > 0)) then
            write(msgstring, *) 'inf_start_from_restart obsolete'
            call error_handler(E_MSG, 'filter_nml: ', msgstring, "", "", "")
            write(msgstring, *) 'use inf_initial_from_restart and inf_sd_initial_from_restart'
            call error_handler(E_MSG, 'filter_nml: ', msgstring, "", "", "")
         endif 
      endif 
      write(msgstring, *) 'INVALID NAMELIST ENTRY: ', trim(nml_string), ' in namelist ', trim(nml_name)
      if(write_to_logfile) then
         call error_handler(E_ERR, 'check_namelist_read', msgstring, &
            source, revision, revdate)
      else
         write(*, *) 'FATAL ERROR before logfile initialization in utilities_mod'
         write(*, *) 'Error is in subroutine check_namelist_read'
         write(*, *) msgstring
         write(*,*)'  ',trim(source)
         write(*,*)'  ',trim(revision)
         write(*,*)'  ',trim(revdate)
         call exit_all(66) 
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


!#######################################################################

subroutine to_upper( string )
! Converts 'string' to uppercase
character(len=*), intent(INOUT) :: string
integer :: ismalla, ibiga, i

ismalla = ichar('a')
ibiga   = ichar('A')

do i = 1,len(string)
   if ((string(i:i) >= 'a') .and. (string(i:i) <= 'z')) then
        string(i:i)  = char( ichar(string(i:i)) + ibiga - ismalla)
   endif
enddo

end subroutine to_upper

!#######################################################################

subroutine find_textfile_dims( fname, nlines, linelen )
! Determines the number of lines and maximum line length
! of the file. Sometimes you need to know this stuff.
character(len=*), intent(IN)  :: fname
integer,          intent(OUT) :: nlines, linelen

integer :: i, mylen, ios, funit

character(len=1024) :: oneline
character(len=129)  :: error_msg

nlines  = 0
linelen = 0
funit   = open_file(fname, form="FORMATTED", action="READ")

READLOOP : do i = 1,100000

   read(funit, '(A)', iostat=ios) oneline
   if (ios < 0) exit READLOOP  ! end of file
   if (ios > 0) then
      write(error_msg,'(A,'' read around line '',i8)')trim(fname),nlines
      call error_handler(E_ERR,'find_textfile_dims', error_msg)
   endif

   nlines = nlines + 1
   mylen  = len_trim(oneline)

   if (mylen > linelen) linelen = mylen

enddo READLOOP

call close_file(funit)

end subroutine find_textfile_dims

!#######################################################################

subroutine file_to_text( fname, textblock )
!
! Reads a text file into a character variable.
! Initially needed to read a namelist file into a variable that could 
! then be inserted into a netCDF file. Due to a quirk in the way Fortran
! and netCDF play together, I have not figured out how to dynamically
! create the minimal character length ... so any line longer than
! the declared length of the textblock variable is truncated.

character(len=*),               intent(IN)  :: fname
character(len=*), dimension(:), intent(OUT) :: textblock

integer :: i, ios, funit
integer :: mynlines, mylinelen, strlen

character(len=512)  :: string

call find_textfile_dims(fname, mynlines, mylinelen)

strlen = len(textblock)

if ( ( mynlines /= size(textblock) ) .or. &
     (mylinelen >      strlen    ) ) then
   write(string,'(A, '' file shape is '',i6,'' by '',i4, &
                    &'' truncating to '',i6,'' by '',i4)') &
   trim(fname),mynlines,mylinelen,size(textblock),strlen
   call error_handler(E_MSG,'file_to_text', trim(string))
endif

funit   = open_file(fname, form="FORMATTED", action="READ")

strlen  = min(mylinelen, strlen)

PARSELOOP : do i = 1,mynlines

   read(funit, '(A)', iostat=ios) string

   write(textblock(i),'(A)') string(1:strlen)

   if ( ios /= 0 ) then
      write(string,'(A,'' read around line '',i8)')trim(fname),i
      call error_handler(E_ERR,'file_to_text', trim(string))
   endif

enddo PARSELOOP

call close_file(funit)

end subroutine file_to_text

!#######################################################################

function is_longitude_between (lon, minlon, maxlon, doradians)

!  uniform way to treat longitude ranges, in degrees, on a globe.
!  returns true if lon is between min and max, starting at min
!  and going EAST until reaching max.  wraps across 0 longitude.
!  if min == max, all points are inside.  includes edges.
!  if optional arg doradians is true, do computation in radians 
!  between 0 and 2*PI instead of 360.

real(r8), intent(in)           :: lon, minlon, maxlon
logical,  intent(in), optional :: doradians
logical :: is_longitude_between

real(r8) :: minl, maxl, lon2, circumf

circumf = 360.0_r8
if (present(doradians)) then
  if (doradians) circumf = TWOPI
endif

minl = modulo(minlon, circumf)
maxl = modulo(maxlon, circumf)

if (minl == maxl) then
   is_longitude_between = .true.  ! entire globe
   return
endif

lon2  = modulo(lon, circumf)

if (minl > maxl) then
   maxl = maxl + circumf
   if (lon2 < minl) lon2 = lon2 + circumf
endif

is_longitude_between = ((lon2 >= minl) .and. (lon2 <= maxl))


end function is_longitude_between 



Function next_file(fname,ifile)
!----------------------------------------------------------------------
! The file name can take one of three forms:
! /absolute/path/to/nirvana/obs_001/obs_seq.final   (absolute path)
! obs_0001/obs_seq.final    (relative path)
! obs_seq.final      (no path ... local)
!
! If there is a '/' in the file name, we grab the portion before the
! slash and look for an underscore. Anything following the underscore
! is presumed to be the portion to increment.
!
! If there is no slash AND ifile is > 1 ... we have already read
! the 'one and only' obs_seq.final file and we return 'done'.
!----------------------------------------------------------------------

character(len=*), intent(in) :: fname
integer,          intent(in) :: ifile

character(len=len(fname)) :: next_file
character(len=len(fname)) :: dir_name

integer,            SAVE :: filenum = 0
integer,            SAVE :: dir_prec = 0
character(len=129), SAVE :: dir_base
character(len=129), SAVE :: filename
character(len=129), SAVE :: dir_ext

integer :: slashindex, splitindex, i, strlen

if (len(fname) > len(dir_base) ) then
   write(msgstring,*)'input filename not guaranteed to fit in local variables'
   call error_handler(E_MSG,'next_file',msgstring, source,revision,revdate)

   write(msgstring,'('' input filename (len='',i3,'') tempvars are (len='',i3,'')'')') &
   len(fname),len(dir_base)

   call error_handler(E_MSG,'next_file',msgstring, source,revision,revdate)
   write(msgstring,*)'increase len of dir_base, filename, dir_ext and recompile'
   call error_handler(E_MSG,'next_file',msgstring, source,revision,revdate)
endif

if (ifile == 1) then ! First time through ... find things.

   ! Start looking (right-to-left) for the 'slash'.
   ! Anything to the right of it must be a filename.
   ! Anything to the left must be the part that gets incremented.

   filename   = adjustl(fname)
   next_file  = trim(filename)
   strlen     = len_trim(filename)
   slashindex = 0

   SlashLoop : do i = strlen,1,-1
   if ( next_file(i:i) == '/' ) then
      slashindex = i
      exit SlashLoop
   endif
   enddo SlashLoop

   if (slashindex > 0) then ! we have a directory structure

      dir_name   = trim(fname(1:slashindex-1))
      filename   = trim(fname(slashindex+1:129))
      strlen     = len_trim(dir_name)
      splitindex = 0

      SplitLoop : do i = strlen,1,-1
      if ( dir_name(i:i) == '_' ) then
         splitindex = i
         exit SplitLoop
      elseif (dir_name(i:i) == '/') then
         ! there is no underscore in the directory node
         ! immediately preceeding the filename
         exit SplitLoop
      endif
      enddo SplitLoop

      if (splitindex <= 0) then
         filenum  = -1 ! indicates no next file
      else
         dir_base   = dir_name(1:splitindex-1)
         dir_ext    = dir_name(splitindex+1:slashindex-1)
         dir_prec   = slashindex - splitindex - 1

         read(dir_ext,*) filenum
      endif

   else ! we have one single file - on the first trip through

      filenum  = -1 ! indicates no next file

   endif 

else

   if (filenum < 0) then
      next_file = 'doneDONEdoneDONE'
   else

      filenum = filenum + 1
      if (dir_prec == 1) then
      write(next_file,'(a,''_'',i1.1,''/'',a)') trim(dir_base),filenum,trim(filename)
      elseif (dir_prec == 2) then
      write(next_file,'(a,''_'',i2.2,''/'',a)') trim(dir_base),filenum,trim(filename)
      elseif (dir_prec == 3) then
      write(next_file,'(a,''_'',i3.3,''/'',a)') trim(dir_base),filenum,trim(filename)
      elseif (dir_prec == 4) then
      write(next_file,'(a,''_'',i4.4,''/'',a)') trim(dir_base),filenum,trim(filename)
      else
      write(next_file,'(a,''_'',i5.5,''/'',a)') trim(dir_base),filenum,trim(filename)
      endif

   endif

endif

end Function next_file


!=======================================================================
! End of utilities_mod
!=======================================================================

end module utilities_mod

