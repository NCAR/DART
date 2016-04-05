! DART software - Copyright 2004 - 2013 UCAR. This open source software is
! provided by UCAR, "as is", without charge, subject to all terms of use at
! http://www.image.ucar.edu/DAReS/DART/DART_download
!
! $Id$

module utilities_mod

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
!      initialize_utilities    Call this routine once before using any
!                              of the other routines in this file.  If you
!                              are using the MPI library, do not call this
!                              routine -- call mpi_initialize_utilities instead.
!                              It will call this routine as part of initializing
!                              the MPI code.
!
!      finalize_utilities      Call this routine at the end of a program to close
!                              the log file and flush pending output.  For an MPI
!                              program, call mpi_finalize_utilities instead.
!
!      dump_unit_attributes    A debug routine that dumps out a long list of
!                              attributes that can be queried from an open file unit.
!
!      print_version_number    Prints out a routine name and
!                              version number to a specified unit
!
!      set_output       Set the status of printing.  Can be set on a per-task
!                       basis if you are running with multiple tasks.
!                       By default all warnings and errors print no matter
!                       which task executes the code; messages only print
!                       from task 0 to avoid N copies of identical messages.
!
!      do_output        Logical function which returns whether informational
!                       messages should be output.  Controlled by the setting
!                       made from set_output.  Useful for messages which cannot
!                       go through the normal error handler (e.g. namelists).
!
!      find_namelist_in_file   Opens the namelist file and verifies the given
!                              namelist name exists.  If found, rewind the file
!                              and return true.  Then a F90 namelist read command
!                              can be used to read in the namelist values.
!                              Otherwise, print an error message and terminate.
!
!      check_namelist_read    Confirms that a namelist read was successful. 
!                             If it failed print an error message and terminate.
!
!      set_nml_output   Set the status of printing namelist data.  By default,
!                       only print to the nml log file.  Can be set to print
!                       to stdout, both, or none.  Argument is a string; valid
!                       values are 'none', 'file', 'terminal', or 'both'.
!
!      do_nml_file      Logical function which returns whether informational
!                       messages should be output to the file.  Controlled
!                       by a call to set_nml_output().
!                     
!      do_nml_term      Logical function which returns whether informational
!                       messages should be output to * (unit 6?).  Controlled
!                       by a call to set_nml_output().
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
!      nmlfileunit      for the namelist file (which defaults to same as log)
!
!      to_upper         converts a character string to uppercase.
!
!      find_textfile_dims    finds number of lines and max line length in a 
!                            text file. Used so we can record the namelist 
!                            file content in the netcdf output files.
!
!      file_to_text     converts the contents of a (hopefully small) file to
!                       a single text variable ... to record in the
!                       netcdf output files.
!
!      get_next_filename     returns the next filename, given the name of
!                            an ascii file which contains a filename per line.
!                            it returns an empty string at end of list.
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
!      ascii_file_format  function that returns true if the string argument
!                         is indicating the requested format is ascii/text.
!                         false means unformatted/binary.  
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

use types_mod, only : r4, r8, digits12, i4, i8, PI
use netcdf

implicit none
private

! module local data

integer, parameter :: E_DBG = -1,   E_MSG = 0,  E_WARN = 1, E_ERR = 2
integer, parameter :: DEBUG = -1, MESSAGE = 0, WARNING = 1, FATAL = 2
integer, parameter :: NML_NONE = 0, NML_FILE = 1, NML_TERMINAL = 2, NML_BOTH = 3

real(r8), parameter :: TWOPI = PI * 2.0_r8

logical :: do_output_flag = .true.
integer :: nml_flag       = NML_FILE
logical :: single_task    = .true.
integer :: task_number    = 0
logical :: module_initialized = .false.
integer :: logfileunit = -1
integer :: nmlfileunit = -1

public :: file_exist, get_unit, open_file, close_file, timestamp,           &
          register_module, error_handler, to_upper, nc_check, next_file,    &
          logfileunit, nmlfileunit, find_textfile_dims, file_to_text,       &
          initialize_utilities, finalize_utilities, dump_unit_attributes,   &
          find_namelist_in_file, check_namelist_read, do_nml_term,          &
          set_tasknum, set_output, do_output, set_nml_output, do_nml_file,  &
          E_DBG, E_MSG, E_WARN, E_ERR, DEBUG, MESSAGE, WARNING, FATAL,      &
          is_longitude_between, get_next_filename, ascii_file_format

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
character(len=256), parameter :: source   = &
   "$URL$"
character(len=32 ), parameter :: revision = "$Revision$"
character(len=128), parameter :: revdate  = "$Date$"

character(len = 169) :: msgstring

!----------------------------------------------------------------
! Namelist input with default values

! E_ERR All warnings/errors are assumed fatal.
integer            :: TERMLEVEL      = E_ERR   

! default log and namelist output filenames
character(len=129) :: logfilename    = 'dart_log.out'
character(len=129) :: nmlfilename    = 'dart_log.nml'

! output each module subversion details
logical            :: module_details = .true.  

! print messages labeled DBG
logical            :: print_debug    = .false. 

! where to write namelist values.
! valid strings:  'none', 'file', 'terminal', 'both'
character(len=32)  :: write_nml      = 'file'  

namelist /utilities_nml/ TERMLEVEL, logfilename, module_details, &
                         nmlfilename, print_debug, write_nml

contains

!#######################################################################

   subroutine initialize_utilities(progname, alternatename, output_flag)
   character(len=*), intent(in), optional :: progname
   character(len=*), intent(in), optional :: alternatename
   logical, intent(in), optional          :: output_flag
   ! integer :: logfileunit -- public module variable
   integer :: iunit, io

   character(len=129) :: lname
   character(len=169) :: string1,string2,string3

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
                           action='write', position='append', iostat = io )
         if ( io /= 0 ) then
            write(*,*)'FATAL ERROR in initialize_utilities'
            write(*,*)'  ',trim(source)
            write(*,*)'  ',trim(revision)
            write(*,*)'  ',trim(revdate)
            write(*,*)'   unable to open the logfile for writing.'
            write(*,*)'   the logfile name is "',trim(lname),'"'
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

         ! Set the defaults for logging the namelist values
         call set_nml_output(write_nml)

         ! If nmlfilename != logfilename, open it.  otherwise set nmlfileunit
         ! to be same as logunit.
         if (do_nml_file()) then
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
         endif

         ! Echo the namelist values for this module using normal mechanism
         ! including a separator line for this run.
         if (do_output_flag) then
            if (do_nml_file() .and. (nmlfileunit /= logfileunit)) then
               if ( present(progname) ) then
                  write(nmlfileunit, *) '!Starting Program '//trim(progname)
               else
                  write(nmlfileunit, *) '!Starting Program '
               endif 
            endif
            if (do_nml_file()) write(nmlfileunit, nml=utilities_nml)
            if (do_nml_term()) write(     *     , nml=utilities_nml)
         endif

         ! Record the values used for variable types:
         if (do_output_flag .and. print_debug) then
        
            write(     *     ,*)  ! a little whitespace is nice
            write(logfileunit,*)  ! a little whitespace is nice

            write(string1,*)'..  digits12 is ',digits12
            write(string2,*)'r8       is ',r8
            write(string3,*)'r4       is ',r4
            call error_handler(E_DBG, 'initialize_utilities', string1, &
                               source, revision, revdate, text2=string2, text3=string3)

            write(string1,*)'..  integer  is ',kind(iunit) ! any integer variable will do
            write(string2,*)'i8       is ',i8
            write(string3,*)'i4       is ',i4
            call error_handler(E_DBG, 'initialize_utilities', string1, &
                               source, revision, revdate, text2=string2, text3=string3)
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


   subroutine finalize_utilities(progname)
   character(len=*), intent(in), optional :: progname
   ! integer :: logfileunit -- private module variable

      ! if called multiple times, just return
      if (.not. module_initialized) return

      if (do_output_flag) then
         if ( present(progname) ) then
            call write_time (logfileunit, label='Finished ', &
                             string1='Program '//trim(progname))
            call write_time (             label='Finished ', &
                             string1='Program '//trim(progname))
         else
            call write_time (logfileunit, label='Finished ')
            call write_time (             label='Finished ')
         endif 

         if (do_nml_file() .and. (nmlfileunit /= logfileunit)) then
            if ( present(progname) ) then
               write(nmlfileunit, *) '!Finished Program '//trim(progname)
            else
               write(nmlfileunit, *) '!Finished Program '
            endif 
         endif 
      endif

      close(logfileunit)
      if ((nmlfileunit /= logfileunit) .and. (nmlfileunit /= -1)) then
         close(nmlfileunit)
      endif

      module_initialized = .false.

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
         call finalize_utilities()
      else if (trim(adjustl(pos)) == 'brief') then
         call write_time (logfileunit, brief=.true., & 
                          string1=string1, string2=string2, string3=string3)
         call write_time (             brief=.true., &
                          string1=string1, string2=string2, string3=string3)
          
      else
         call write_time (logfileunit, & 
                          string1=string1, string2=string2, string3=string3)
         call write_time (string1=string1, string2=string2, string3=string3)
          
      endif

   end subroutine timestamp

!#######################################################################

   function file_exist (file_name)

      character(len=*), intent(in) :: file_name
      logical :: file_exist
      integer :: trimlen

      if ( .not. module_initialized ) call initialize_utilities

      trimlen = len_trim(file_name)

      inquire (file=file_name(1:trimlen), exist=file_exist)

   end function file_exist

!#######################################################################

   function get_unit () result (iunit)

      integer :: i, iunit
      logical :: available

      if ( .not. module_initialized ) call initialize_utilities

! ---- get available unit ----

      iunit = -1
      do i = 10, 80
         inquire (i, opened=available)
         if (.not. available) then
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

  subroutine error_handler(level, routine, text, src, rev, rdate, aut, text2, text3 )
!----------------------------------------------------------------------
! subroutine error_handler(level, routine, text, src, rev, rdate, aut , text2, text3)
!
! logs warning/error 
implicit none

integer, intent(in) :: level
character(len = *), intent(in) :: routine, text
character(len = *), intent(in), optional :: src, rev, rdate, aut, text2, text3

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
        if ( present(text2)) then
           write(     *     , *) trim(routine),' ... ', trim(text2)
           write(logfileunit, *) trim(routine),' ... ', trim(text2)
        endif
        if ( present(text3)) then
           write(     *     , *) trim(routine),' ... ', trim(text3)
           write(logfileunit, *) trim(routine),' ... ', trim(text3)
        endif
      else
        ! FIXME: should they just all use i5? but most common case is only
        ! messages from PE0, so it's tempting not to waste all those columns.
        if (task_number == 0) then
            write(taskstr, '(a)' ) "PE 0"
        else
            write(taskstr, '(a,i5)' ) "PE ", task_number
        endif
        write(     *     , *) trim(taskstr),': ',trim(routine),' ', trim(text)
        write(logfileunit, *) trim(taskstr),': ',trim(routine),' ', trim(text)
        if ( present(text2)) then
           write(     *     , *) trim(taskstr),': ',trim(routine),' ... ', trim(text2)
           write(logfileunit, *) trim(taskstr),': ',trim(routine),' ... ', trim(text2)
        endif
        if ( present(text3)) then
           write(     *     , *) trim(taskstr),': ',trim(routine),' ... ', trim(text3)
           write(logfileunit, *) trim(taskstr),': ',trim(routine),' ... ', trim(text3)
        endif
      endif

   case (E_DBG)
      if (print_debug) then

         ! what about do_output_flag?  want messages from all procs or just PE0?

         if ( single_task ) then
           write(     *     , *) 'DEBUG FROM: ', trim(routine),' ', trim(text)
           write(logfileunit, *) 'DEBUG FROM: ', trim(routine),' ', trim(text)
           if ( present(text2)) then
              write(     *     , *) 'DEBUG FROM: ', trim(routine),' ... ', trim(text2)
              write(logfileunit, *) 'DEBUG FROM: ', trim(routine),' ... ', trim(text2)
           endif
           if ( present(text3)) then
              write(     *     , *) 'DEBUG FROM: ', trim(routine),' ... ', trim(text3)
              write(logfileunit, *) 'DEBUG FROM: ', trim(routine),' ... ', trim(text3)
           endif
         else
           if (task_number == 0) then
               write(taskstr, '(a)' ) "PE 0"
           else
               write(taskstr, '(a,i5)' ) "PE ", task_number
           endif
           write(     *     , *) trim(taskstr),': DEBUG FROM: ',trim(routine),' ', trim(text)
           write(logfileunit, *) trim(taskstr),': DEBUG FROM: ',trim(routine),' ', trim(text)
           if ( present(text2)) then
              write(     *     , *) trim(taskstr),': DEBUG FROM: ',trim(routine),' ... ', trim(text2)
              write(logfileunit, *) trim(taskstr),': DEBUG FROM: ',trim(routine),' ... ', trim(text2)
           endif
           if ( present(text3)) then
              write(     *     , *) trim(taskstr),': DEBUG FROM: ',trim(routine),' ... ', trim(text3)
              write(logfileunit, *) trim(taskstr),': DEBUG FROM: ',trim(routine),' ... ', trim(text3)
           endif
         endif
      endif

   case (E_WARN)

      write(     *     , *) 'WARNING FROM:'
      if ( .not. single_task ) &
      write(     *     , *) ' task id: ', task_number
      write(     *     , *) ' routine: ', trim(routine)
      write(     *     , *) ' message: ', trim(text)
      if ( present(text2)) &
      write(     *     , *) ' message: ... ', trim(text2)
      if ( present(text3)) &
      write(     *     , *) ' message: ... ', trim(text3)
      write(     *     , *) ' '
      write(     *     , *) ' source file: ', trim(src)
      write(     *     , *) ' file revision: ', trim(rev)
      write(     *     , *) ' revision date: ', trim(rdate)
      if(present(aut)) &
      write(     *     , *) ' last editor: ', trim(aut)

      write(logfileunit, *) 'WARNING FROM:'
      if ( .not. single_task ) &
      write(logfileunit, *) ' task id: ', task_number
      write(logfileunit, *) ' routine: ', trim(routine)
      write(logfileunit, *) ' message: ', trim(text)
      if ( present(text2)) &
      write(logfileunit, *) ' message: ... ', trim(text2)
      if ( present(text3)) &
      write(logfileunit, *) ' message: ... ', trim(text3)
      write(logfileunit, *) ' '
      write(logfileunit, *) ' source file: ', trim(src)
      write(logfileunit, *) ' file revision: ', trim(rev)
      write(logfileunit, *) ' revision date: ', trim(rdate)
      if(present(aut)) &
      write(logfileunit, *) ' last editor: ', trim(aut)

   case(E_ERR)

      write(     *     , *) 'ERROR FROM:'
      if ( .not. single_task ) &
      write(     *     , *) ' task id: ', task_number
      write(     *     , *) ' routine: ', trim(routine)
      write(     *     , *) ' message: ', trim(text)
      if ( present(text2)) &
      write(     *     , *) ' message: ... ', trim(text2)
      if ( present(text3)) &
      write(     *     , *) ' message: ... ', trim(text3)
      write(     *     , *) ' '
      write(     *     , *) ' source file: ', trim(src)
      write(     *     , *) ' file revision: ', trim(rev)
      write(     *     , *) ' revision date: ', trim(rdate)
      if(present(aut)) &
      write(     *     , *) ' last editor: ', trim(aut)

      write(logfileunit, *) 'ERROR FROM:'
      if ( .not. single_task ) &
      write(logfileunit, *) ' task id: ', task_number
      write(logfileunit, *) ' routine: ', trim(routine)
      write(logfileunit, *) ' message: ', trim(text)
      if ( present(text2)) &
      write(logfileunit, *) ' message: ... ', trim(text2)
      if ( present(text3)) &
      write(logfileunit, *) ' message: ... ', trim(text3)
      write(logfileunit, *) ' '
      write(logfileunit, *) ' source file: ', trim(src)
      write(logfileunit, *) ' file revision: ', trim(rev)
      write(logfileunit, *) ' revision date: ', trim(rdate)
      if(present(aut)) &
      write(logfileunit, *) ' last editor: ', trim(aut)

end select

! TERMLEVEL gets set in the namelist
if( level >= TERMLEVEL ) call exit_all( 99 ) 

end subroutine error_handler

!#######################################################################

   function open_file (fname, form, action) result (iunit)

   character(len=*), intent(in) :: fname
   character(len=*), intent(in), optional :: form, action
   integer  :: iunit

   integer           :: nc, rc
   logical           :: open
   character(len=11) :: format
   character(len=6)  :: pos
   character(len=9)  :: act
   character(len=7)  :: stat

   if ( .not. module_initialized ) call initialize_utilities

   inquire (file=trim(fname), opened=open, number=iunit,  &
            form=format, iostat=rc)

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

      ! this code used to only set the form and position, not the action.
      ! not specifying 'read' meant that many compilers would create an
      ! empty file instead of returning a read error.  this leads to lots
      ! of confusion.  add an explicit action here.  if the incoming argument
      ! is read, make sure the open() call passes that in as an action.

      format   = 'formatted'
      act      = 'readwrite'
      pos      = 'rewind'
      stat     = 'unknown'

      if (present(form)) then
          nc = min(len(format),len(form))
          format(1:nc) = form(1:nc)
      endif

      if (present(action)) then
          select case(action)

             case ('read', 'READ')
             ! open existing file.  fail if not found.  read from start.
                act  = 'read'
                stat = 'old'
                pos  = 'rewind'

             case ('write', 'WRITE')
             ! create new file/replace existing file.  write at start.
                act  = 'write'
                stat = 'replace'
                pos  = 'rewind'

             case ('append', 'APPEND')
             ! create new/open existing file.  write at end.
                act  = 'readwrite'
                stat = 'unknown'
                pos  = 'append'

             case default
             ! leave defaults specified above, currently
             ! create new/open existing file.  write at start.
                !print *, 'action specified, and is ', action
          end select
      endif

      iunit = get_unit()

      if (format == 'formatted' .or. format == 'FORMATTED') then
          open (iunit, file=trim(fname), form=format,     &
                position=pos, delim='apostrophe',    &
                action=act, status=stat, iostat=rc)
      else
          open (iunit, file=trim(fname), form=format,     &
                position=pos, action=act, status=stat, iostat=rc)
      endif

      if (rc /= 0) then
         write(msgstring,*)'Cannot open file "'//trim(fname)//'" for '//trim(act)
         call error_handler(E_ERR, msgstring, source, revision, revdate)
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

   subroutine write_time (unit, label, string1, string2, string3, tz, brief)

! ***  Write the current time to a log file or standard output ***
!
!    in: unit number (default is * if not specified)
!    in: label (default is  "Time is" if not specified)
!    in: string1,2,3 (no defaults)
!
!  default output is a block of 3-4 lines, with dashed line separators
!  and up to 3 descriptive text strings.
!  if brief specified as true, only string1 printed if given,
!  and time printed on same line in YYYY/MM/DD HH:MM:SS format
!  with the tag 'TIME:' before it.  should be easier to postprocess.

   integer,          optional, intent(in) :: unit
   character(len=*), optional, intent(in) :: label
   character(len=*), optional, intent(in) :: string1
   character(len=*), optional, intent(in) :: string2
   character(len=*), optional, intent(in) :: string3
   logical,          optional, intent(in) :: tz
   logical,          optional, intent(in) :: brief


   integer :: lunit
   character(len= 8) :: cdate
   character(len=10) :: ctime
   character(len= 5) :: zone
   integer, dimension(8) :: values
   logical :: oneline

   if (present(unit)) then
      lunit = unit
   else
      lunit = 6   ! this should be *
   endif

   call DATE_AND_TIME(cdate, ctime, zone, values)

   ! give up if no good values were returned
   if (.not. any(values /= -HUGE(0)) ) return 

   oneline = .false.
   if (present(brief)) oneline = brief

   if (oneline) then
      if (present(string1)) then
         write(lunit,'(A,1X,I4,5(A1,I2.2))') string1//' TIME:', &
                        values(1), '/', values(2), '/', values(3), &
                        ' ', values(5), ':', values(6), ':', values(7)
      else
         write(lunit,'(A,1X,I4,5(A1,I2.2))') 'TIME: ', &
                        values(1), '/', values(2), '/', values(3), &
                        ' ', values(5), ':', values(6), ':', values(7)
      endif
   else
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
   endif

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

   subroutine set_nml_output (nmlstring)

! *** set whether nml output is written to stdout file or only nml file
!
!    in:  doflag  = whether to output nml information to stdout 

   character(len=*), intent(in) :: nmlstring

   if ( .not. module_initialized ) call initialize_utilities

   select case (nmlstring)
      case ('NONE', 'none')
         nml_flag = NML_NONE
         call error_handler(E_MSG, 'set_nml_output', &
                            'No echo of NML values')

      case ('FILE', 'file')
         nml_flag = NML_FILE
         call error_handler(E_MSG, 'set_nml_output', &
                            'Echo NML values to log file only')
  
      case ('TERMINAL', 'terminal')
         nml_flag = NML_TERMINAL
         call error_handler(E_MSG, 'set_nml_output', &
                            'Echo NML values to terminal output only')
   
      case ('BOTH', 'both')
         nml_flag = NML_BOTH
         call error_handler(E_MSG, 'set_nml_output', &
                            'Echo NML values to both log file and terminal')

      case default
         call error_handler(E_ERR, 'set_nml_output', &
           'unrecognized input string: '//trim(nmlstring), &
           source, revision, revdate)
 
   end select

   end subroutine set_nml_output


!#######################################################################

   function do_nml_file ()

! *** return whether nml should be written to nml file
!

   logical :: do_nml_file

   if ( .not. module_initialized ) call initialize_utilities

   if ( .not. do_output()) then
      do_nml_file = .false.
   else
      do_nml_file = (nml_flag == NML_FILE .or. nml_flag == NML_BOTH)
   endif

   end function do_nml_file

!#######################################################################

   function do_nml_term ()

! *** return whether nml should be written to terminal
!

   logical :: do_nml_term

   if ( .not. module_initialized ) call initialize_utilities

   if ( .not. do_output()) then
      do_nml_term = .false.
   else
      do_nml_term = (nml_flag == NML_TERMINAL .or. nml_flag == NML_BOTH)
   endif

   end function do_nml_term

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

character(len = 169) :: nml_string, test_string, string1
integer              :: io
logical              :: write_to_logfile


! Decide if there is a logfile or not
write_to_logfile = .true.
if(present(write_to_logfile_in)) write_to_logfile = write_to_logfile_in

! Check for file existence; no file is an error
if(file_exist(trim(namelist_file_name))) then

   iunit = open_file(trim(namelist_file_name), action = 'read')

   ! Read each line until end of file is found
   ! Look for the start of a namelist with &nml_name
   ! Convert test string to all uppercase ... since that is
   ! what happens if Fortran writes a namelist.

   string1 = adjustl(nml_name)
   call to_upper(string1)             ! works in-place
   test_string = '&' // trim(string1)

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

         string1 = adjustl(nml_string)
         call to_upper(string1)

         if(trim(string1) == trim(test_string)) then
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
            call error_handler(E_MSG, 'filter_nml: ', msgstring)
            write(msgstring, *) 'use inf_initial_from_restart and inf_sd_initial_from_restart'
            call error_handler(E_MSG, 'filter_nml: ', msgstring)
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
character(len=*),  intent(IN)  :: fname
integer,           intent(OUT) :: nlines
integer, optional, intent(OUT) :: linelen

integer :: i, maxlen, mylen, ios, funit

character(len=1024) :: oneline
character(len=129)  :: error_msg

! if there is no file, return -1 for both counts
if (.not. file_exist(fname)) then
  nlines = -1
  if (present(linelen)) linelen = -1
  return
endif

! the file exists, go count things up.
nlines  = 0
maxlen  = 0
funit   = open_file(fname, form="FORMATTED", action="READ")

READLOOP : do i = 1,100000

   read(funit, '(A)', iostat=ios) oneline
   if (ios < 0) exit READLOOP  ! end of file
   if (ios > 0) then
      write(error_msg,'(A,'' read around line '',i8)')trim(fname),nlines
      call error_handler(E_ERR,'find_textfile_dims', error_msg, &
                         source, revision, revdate)
   endif

   nlines = nlines + 1
   mylen  = len_trim(oneline)

   if (mylen > maxlen) maxlen = mylen

enddo READLOOP

call close_file(funit)

if (present(linelen)) linelen = maxlen

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
      call error_handler(E_ERR,'file_to_text', trim(string), &
                         source, revision, revdate)
   endif

enddo PARSELOOP

call close_file(funit)

end subroutine file_to_text

!#######################################################################

function get_next_filename( listname, index )

! Arguments are the name of a file which contains a list of filenames.
! This routine opens the listfile, and returns the index-th one.
!
character(len=*),  intent(in) :: listname
integer,           intent(in) :: index
character(len=128)            :: get_next_filename

integer :: i, ios, funit

character(len=512)  :: string

funit   = open_file(listname, form="FORMATTED", action="READ")

PARSELOOP : do i=1, index

   read(funit, '(A)', iostat=ios) string

   ! reached end of file, return '' as indicator.
   if ( ios /= 0 ) then
      get_next_filename = ''
      call close_file(funit)
      return
   endif

enddo PARSELOOP

! check for length problems
if (len_trim(string) > len(get_next_filename)) then
   call error_handler(E_ERR, 'get_next_filename', &
                      'maximum filename length of 128 exceeded', &
                      source, revision, revdate)   
endif

get_next_filename = adjustl(string(1:len(get_next_filename)))
call close_file(funit)

end function get_next_filename

!#######################################################################

function is_longitude_between (lon, minlon, maxlon, doradians, newlon)

!  uniform way to treat longitude ranges, in degrees, on a globe.
!  returns true if lon is between min and max, starting at min
!  and going EAST until reaching max.  wraps across 0 longitude.
!  if min == max, all points are inside.  includes edges.
!  if optional arg doradians is true, do computation in radians 
!  between 0 and 2*PI instead of 360.   if given, return the
!  'lon' value possibly + 360 (or 2PI) which can be used for averaging
!  or computing on a consistent set of longitude values.  after the
!  computation is done if the answer is > 360 (or 2PI), subtract that
!  value to get back into the 0 to 360 (or 2PI) range.

real(r8), intent(in)            :: lon, minlon, maxlon
logical,  intent(in),  optional :: doradians
real(r8), intent(out), optional :: newlon
logical :: is_longitude_between

real(r8) :: minl, maxl, lon2, circumf

circumf = 360.0_r8
if (present(doradians)) then
  if (doradians) circumf = TWOPI
endif

! ensure the valid region boundaries are between 0 and one circumference
! (must use modulo() and not mod() so negative vals are handled ok)
minl = modulo(minlon, circumf)
maxl = modulo(maxlon, circumf)

! boundary points are included in the valid region so if min=max 
! the 'region' is the entire globe and you can return early.
if (minl == maxl) then
   is_longitude_between = .true. 
   if (present(newlon)) newlon = lon
   return
endif

! ensure the test point is between 0 and one circumference
lon2  = modulo(lon, circumf)

! here's where the magic happens:
! minl will be bigger than maxl if the region of interest crosses the prime 
! meridian (longitude = 0).  in this case add one circumference to the 
! eastern boundary so maxl is guarenteed to be larger than minl (and valid 
! values are now between 0 and 2 circumferences).  
!
! if the test point longitude is west of the minl boundary add one circumference
! to it as well before testing against the bounds.  values that were east of 
! longitude 0 but west of maxl will now be shifted so they are again correctly 
! within the new range; values that were west of the prime meridian but east 
! of minl will stay in range; values west of minl and east of maxl will be 
! correctly shifted out of range.

if (minl > maxl) then
   maxl = maxl + circumf
   if (lon2 < minl) lon2 = lon2 + circumf
endif

is_longitude_between = ((lon2 >= minl) .and. (lon2 <= maxl))

! if requested, return the value that was tested against the bounds, which 
! will always be between 0 and 2 circumferences and monotonically increasing
! from minl to maxl.  if the region of interest doesn't cross longitude 0
! this value will be the same as the input value.  if the region does
! cross longitude 0 this value will be between 0 and 2 circumferences.
! it's appropriate for averaging values together or comparing them against
! other values returned from this routine with a simple greater than or less
! than without further computation for longitude 0.  to convert the values
! back into the range from 0 to one circumference, compare it to the
! circumference and if larger, subtract one circumference from the value.

if (present(newlon)) newlon = lon2

end function is_longitude_between 


!#######################################################################

function next_file(fname,ifile)
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

integer :: slashindex, splitindex, i, strlen, ios

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
    
         read(dir_ext,*,iostat=ios) filenum
         if(ios /= 0) then
            ! Directory has an '_' separating two alphabetic parts
            ! Nothing to increment.
            filenum = -1
         endif
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

end function next_file

!#######################################################################


function ascii_file_format(fform)

!----------------------------------------------------------------------
! Common routine for setting read/write file format.

character(len=*), intent(in), optional :: fform
logical                                :: ascii_file_format

character(len=129) :: lj_fform ! Left Justified version of optional argument 

! Returns .true. for formatted/ascii file, .false. is unformatted/binary
! Defaults (if fform not specified) to formatted/ascii.

if ( .not. module_initialized ) call initialize_utilities

! Default to formatted/ascii.
if ( .not. present(fform)) then
   ascii_file_format = .true.
   return
endif

! Check to make sure we don't put 10lbs of stuff in a 5lb bag

if (len(fform) > len(lj_fform)) then
   write(msgstring,*)'fform is long: increase len of lj_fform to ',&
                     len(fform),' and recompile.'
   call error_handler(E_ERR,'ascii_file_format', msgstring, source, revision, revdate)
endif

lj_fform = adjustl(fform)

SELECT CASE (trim(lj_fform))
   CASE("unf", "UNF", "unformatted", "UNFORMATTED")
      ascii_file_format = .false.
   CASE DEFAULT
      ascii_file_format = .true.
END SELECT

end function ascii_file_format


!=======================================================================
! End of utilities_mod
!=======================================================================

end module utilities_mod

! <next few lines under version control, do not edit>
! $URL$
! $Id$
! $Revision$
! $Date$
