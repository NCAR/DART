! Data Assimilation Research Testbed -- DART
! Copyright 2004, Data Assimilation Initiative, University Corporation for Atmospheric Research
! Licensed under the GPL -- www.gpl.org/licenses/gpl.html

module utilities_mod

! <next four lines automatically updated by CVS, do not edit>
! $Source$ 
! $Revision$ 
! $Date$ 
! $Author$ 
!
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
!      check_nml_error  Checks the iostat returned when reading
!                       namelists and determines if the error code
!                       is valid, if not the program is terminated.
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
!-----------------------------------------------------------------------

use types_mod, only : r8

implicit none
private

!   ---- private data for check_nml_error ----

integer, private :: num_nml_error_codes, nml_error_codes(5)
logical, private :: do_nml_error_init = .true.
private  nml_error_init

integer, parameter ::   E_MSG = 0,  E_WARN = 1, E_ERR = 2
integer, parameter :: MESSAGE = 0, WARNING = 1, FATAL = 2

public file_exist, get_unit, check_nml_error, open_file, timestamp, &
       close_file, register_module, error_handler, initialize_utilities, &
       finalize_utilities, E_MSG, E_WARN, E_ERR, MESSAGE, WARNING, FATAL, logfileunit

! CVS Generated file description for error handling, do not edit
character(len=128) :: &
source   = "$Source$", &
revision = "$Revision$", &
revdate  = "$Date$"

logical, save :: module_initialized = .false.
integer, save :: logfileunit = -1

!----------------------------------------------------------------
! Namelist input with default values
integer  :: TERMLEVEL = E_ERR     ! E_ERR All warnings/errors are assumed fatal.
character(len=129) :: logfilename = 'dart_log.out'

namelist /utilities_nml/TERMLEVEL, logfilename

contains

!#######################################################################

   subroutine initialize_utilities
   ! integer :: logfileunit -- public module variable
   integer :: iunit, io
   logical :: lfile

   character(len= 8) :: date
   character(len=10) :: time
   character(len= 5) :: zone
   integer, dimension(8) :: values


      if ( module_initialized ) then ! nothing to do

         ! write(*,*)'Module initialized ... carry on.'

         return

      else ! initialize the module

         ! Since the logfile is not open yet, the error terminations
         ! must be handled differently than all other cases.
         ! The routines that normally write to the logfile cannot
         ! be used just yet. If we cannot open a logfile, we
         ! always abort execution at this step.

         write(*,*)'Initializing the utilities module.'

         inquire (file = 'input.nml', exist = lfile) ! does a namelist exist

         ! Read the namelist input if it exists
         if( lfile ) then

            ! write(*,*)'Determine unit for namelist'

            ! determine an open unit for the namelist 
            iunit = nextunit()
            if ( iunit < 0 ) then
               write(*,*)'   unable to get a unit to use to read the namelist.'
               write(*,*)'   stopping.'
               stop 99
            endif

            ! If the file exists, it might still not contain the utilities_nml namelist
            open( iunit, file='input.nml')
            read( iunit, nml = utilities_nml, iostat = io)
            if (io /= 0 ) then
               write(*,*)'  No utilities_nml in input.nml, iostat was ',io 
               write(*,*)'  using default values ...'
            endif
            close(iunit)

         endif

         ! Open the log file with the name from the namelist 
         logfileunit = nextunit()
         if ( logfileunit < 0 ) then
            write(*,*)'   unable to get a unit to use for the logfile.'
            write(*,*)'   stopping.'
            stop 99
         endif

         write(*,*)'Trying to read from unit ', logfileunit
         write(*,*)'Trying to open file ', trim(adjustl(logfilename))

         open(logfileunit, file=trim(adjustl(logfilename)), form='formatted', &
                           position='append', iostat = io )
         if ( io /= 0 ) then
            write(*,*)'FATAL ERROR in initialize_utilities'
            write(*,*)'  ',trim(source)
            write(*,*)'  ',trim(revision)
            write(*,*)'  ',trim(revdate)
            write(*,*)'   unable to open the logfile.'
            write(*,*)'   the intended file name was <',trim(logfilename),'>'
            write(*,*)'   stopping.'
            stop 99
         endif

         ! Log the run-time 

         call DATE_AND_TIME(date, time, zone, values)

         write(logfileunit,*)
         write(logfileunit,*)'--------------------------------------'
         write(logfileunit,*)'Running  ... at YYYY MM DD HH MM SS = '
         write(logfileunit,'(17x,i4,5(1x,i2))') values(1), values(2), &
                           values(3),  values(5), values(6), values(7)
         if ( values(4) /= -HUGE(0) ) &
         write(logfileunit,*)'time zone offset is ',values(4),' minutes.'
         write(logfileunit,*)'--------------------------------------'
         write(logfileunit,*)

         ! Check to make sure termlevel is set to a reasonable value
         call checkTermLevel

         module_initialized = .true.

         ! Echo the module information using normal mechanism
         call register_module(source, revision, revdate)

         ! Echo the namelist values for this module using normal mechanism
         write(logfileunit, nml=utilities_nml)
         write(     *     , nml=utilities_nml)

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

      character(len= 8) :: date
      character(len=10) :: time
      character(len= 5) :: zone
      integer, dimension(8) :: values

      if ( .not. module_initialized ) call initialize_utilities

      call DATE_AND_TIME(date, time, zone, values)

      ! check to see if values are valid on this system

      if ( any(values /= -HUGE(0)) ) then ! at least one is valid

         SELECT CASE ( pos )

            CASE ( 'end' )

               write(logfileunit,*)
               write(logfileunit,*)'--------------------------------------'
               write(logfileunit,*)'Finished ... at YYYY MM DD HH MM SS = '
               write(logfileunit,'(17x,i4,5(1x,i2))') values(1), values(2), &
                           values(3),  values(5), values(6), values(7)

               if(present(string1)) write(logfileunit,*)trim(string1)
               if(present(string2)) write(logfileunit,*)trim(string2)
               if(present(string3)) write(logfileunit,*)trim(string3)
               write(logfileunit,*)'--------------------------------------'

               call finalize_utilities
         
            CASE DEFAULT

               write(logfileunit,*)
               write(logfileunit,*)'--------------------------------------'
               write(logfileunit,*)'Time is YYYY MM DD HH MM SS = '
               write(logfileunit,'(9x,i4,5(1x,i2))') values(1), values(2), &
                           values(3),  values(5), values(6), values(7)

               if(present(string1)) write(logfileunit,*)trim(string1)
               if(present(string2)) write(logfileunit,*)trim(string2)
               if(present(string3)) write(logfileunit,*)trim(string3)

               if ( values(4) /= -HUGE(0) ) &
               write(logfileunit,*)'time zone offset is ',values(4),' minutes.'
               write(logfileunit,*)'--------------------------------------'

         END SELECT
          
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

if ( .not. module_initialized ) call initialize_utilities

select case(level)
   case (E_MSG)

      write(     *     , *) trim(routine),' ', trim(text)
      write(logfileunit, *) trim(routine),' ', trim(text)

   case (E_WARN)

      write(     *     , *) 'WARNING FROM:'
      write(     *     , *) '   routine       : ', trim(routine)
      write(     *     , *) '   source file   : ', trim(src)
      write(     *     , *) '   file revision : ', trim(rev)
      write(     *     , *) '   revision date : ', trim(rdate)
      write(     *     , *) '   message       : ', trim(text)

      write(logfileunit, *) 'WARNING FROM:'
      write(logfileunit, *) '   routine       : ', trim(routine)
      write(logfileunit, *) '   source file   : ', trim(src)
      write(logfileunit, *) '   file revision : ', trim(rev)
      write(logfileunit, *) '   revision date : ', trim(rdate)
      if(present(aut)) &
      write(logfileunit, *) '   last editor   : ', trim(aut)
      write(logfileunit, *) '   message       : ', trim(text)

   case(E_ERR)

      write(     *     , *) 'ERROR FROM:'
      write(     *     , *) '   routine       : ', trim(routine)
      write(     *     , *) '   source file   : ', trim(src)
      write(     *     , *) '   file revision : ', trim(rev)
      write(     *     , *) '   revision date : ', trim(rdate)
      write(     *     , *) '   message       : ', trim(text)

      write(logfileunit, *) 'ERROR FROM:'
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

function check_nml_error (iostat, nml_name) result (error_code)

   integer,          intent(in) :: iostat
   character(len=*), intent(in) :: nml_name
   integer   error_code, i
   character(len=128) :: err_str

   if ( .not. module_initialized ) call initialize_utilities

   if (do_nml_error_init) call nml_error_init

   error_code = iostat

   do i = 1, num_nml_error_codes
        if (error_code == nml_error_codes(i)) return
   enddo

!  ------ fatal namelist error -------
   write (err_str,*) 'while reading namelist ',  &
                     nml_name(1:len_trim(nml_name)),  &
                     ', iostat = ',error_code

   call error_handler(E_ERR, 'check_nml_error', &
            err_str, source, revision, revdate)

end function check_nml_error

!-----------------------------------------------------------------------

subroutine nml_error_init

!   private routine for initializing allowable error codes

   integer  iunit, io, ir
   real(r8) ::  b=1.0_r8
   namelist /b_nml/  b

      if ( .not. module_initialized ) call initialize_utilities

      nml_error_codes(1) = 0

!     ---- create dummy namelist file ----
      iunit=get_unit(); open (iunit, file='_read_error.nml')
      write (iunit, 10)
  10  format (' &a_nml  a=1.  /',  &
             /'-------------------',  &
             /' &b_nml  e=5.  &end')
      close (iunit)

!     ---- read namelist file and save error codes ----
      iunit=get_unit(); open (iunit, file='_read_error.nml')
      ir=1; io=1; do
         read  (iunit, nml=b_nml, iostat=io, end=20)
         if (io == 0) exit
         ir=ir+1; nml_error_codes(ir)=io
      enddo
  20  close (iunit, status='delete')

      num_nml_error_codes = ir
!del  print *, 'nml_error_codes=', nml_error_codes(1:ir)
      do_nml_error_init = .false.

end subroutine nml_error_init

!#######################################################################
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

!
!=======================================================================
! End of utilities_mod
!=======================================================================
!
end module utilities_mod
