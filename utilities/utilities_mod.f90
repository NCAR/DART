module utilities_mod
!
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
!      error_mesg       Print warning and error messages, 
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

implicit none

private

!   ---- private data for check_nml_error ----

   integer, private :: num_nml_error_codes, nml_error_codes(5)
   logical, private :: do_nml_error_init = .true.
   private  nml_error_init

integer, parameter :: E_MSG = 0, E_WARN = 1, E_ERR = 2
integer, parameter :: MESSAGE = 0, WARNING = 1, FATAL = 2

public file_exist, get_unit, error_mesg, check_nml_error, open_file, &
       close_file, print_version_number, output_err, &
       E_MSG, E_WARN, E_ERR, MESSAGE, WARNING, FATAL 

! let CVS fill strings ... DO NOT EDIT ...
character(len=128) :: &
   source   = "$Source$", &
   revision = "$Revision$", &
   revdate  = "$Date$"

contains

!#######################################################################

   function file_exist (file_name)

      character(len=*), intent(in) :: file_name
      logical  file_exist

      inquire (file=file_name(1:len_trim(file_name)), exist=file_exist)

   end function file_exist

!#######################################################################

   function get_unit () result (unit)

      integer  i, unit
      logical  open

! ---- get available unit ----

      unit = -1
      do i = 10, 80
         inquire (i, opened=open)
         if (.not.open) Then
            unit = i
            exit
         endif
      enddo

      if (unit == -1) Then
         call error_mesg ('get_unit', 'no available units.', 1)
      endif

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
                stop
             case default
                print *, ' ERROR in ',routine(1:len_trim(routine))
                print *, ' ',message(1:len_trim(message))
                stop
          end select

!         --------------------------------------------

   end subroutine error_mesg

!#######################################################################

subroutine output_err(level, src, rev, date, aut, routine, text)

implicit none

integer, intent(in) :: level
character(len = *), intent(in) :: src, rev, date, aut, routine, text

select case(level)
   case (E_MSG)
      write(*, *) 'MESSAGE FROM:'
   case (E_WARN)
      write(*, *) 'WARNING FROM:'
   case(E_ERR)
      write(*, *) 'ERROR FROM:'
end select

write(*, *) trim(src)
write(*, *) trim(rev)
write(*, *) trim(date)
write(*, *) trim(aut)
write(*, *) 'In routine ', trim(routine)
write(*, *) trim(text)

! Stop for all but message; 
if(level /= E_MSG) stop

end subroutine output_err

!#######################################################################

function check_nml_error (iostat, nml_name) result (error_code)

   integer,          intent(in) :: iostat
   character(len=*), intent(in) :: nml_name
   integer   error_code, i
   character(len=128) :: err_str

   if (do_nml_error_init) call nml_error_init

   error_code = iostat

   do i = 1, num_nml_error_codes
        if (error_code == nml_error_codes(i)) return
   enddo

!  ------ fatal namelist error -------
   write (err_str,*) 'while reading namelist ',  &
                     nml_name(1:len_trim(nml_name)),  &
                     ', iostat = ',error_code

   call error_mesg ('check_nml_error', err_str, 3)

end function check_nml_error

!-----------------------------------------------------------------------

subroutine nml_error_init

!   private routine for initializing allowable error codes

   integer  unit, io, ir
   real ::  b=1.
   namelist /b_nml/  b

      nml_error_codes(1) = 0

!     ---- create dummy namelist file ----
      unit=get_unit(); open (unit, file='_read_error.nml')
      write (unit, 10)
  10  format (' &a_nml  a=1.  /',  &
             /'-------------------',  &
             /' &b_nml  e=5.  &end')
      close (unit)

!     ---- read namelist file and save error codes ----
      unit=get_unit(); open (unit, file='_read_error.nml')
      ir=1; io=1; do
         read  (unit, nml=b_nml, iostat=io, end=20)
         if (io == 0) exit
         ir=ir+1; nml_error_codes(ir)=io
      enddo
  20  close (unit, status='delete')

      num_nml_error_codes = ir
!del  print *, 'nml_error_codes=', nml_error_codes(1:ir)
      do_nml_error_init = .false.

end subroutine nml_error_init

!#######################################################################
!#######################################################################

   function open_file (file, form, action) result (unit)

   character(len=*), intent(in) :: file
   character(len=*), intent(in), optional :: form, action
   integer  :: unit

   integer           :: nc
   logical           :: open
   character(len=11) :: format
   character(len=6)  :: location

      inquire (file=file(1:len_trim(file)), opened=open, number=unit,  &
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

         unit = get_unit()

         if (format == 'formatted  ' .or. format == 'FORMATTED  ') then
             open (unit, file=file(1:len_trim(file)),      &
                         form=format(1:len_trim(format)),  &
                     position=location(1:len_trim(location)),  &
                        delim='apostrophe')
         else
             open (unit, file=file(1:len_trim(file)),      &
                         form=format(1:len_trim(format)),  &
                     position=location(1:len_trim(location)) )
         endif
     endif


   end function open_file

!#######################################################################

   subroutine print_version_number (unit, routine, version)

! *** prints routine name and version number to a log file ***
!
!    in:  unit    = unit number to direct output
!         routine = routine name (character, max len = 20)
!         version = version name or number (character, max len = 8)

   integer,          intent(in) :: unit
   character(len=*), intent(in) :: routine, version

   integer           :: n
   character(len=20) :: name
   character(len=8)  :: vers

     n = min(len(routine),20); name = adjustl(routine(1:n))
     n = min(len(version), 8); vers = adjustl(version(1:n))

     if (unit > 0) then
         write (unit,10) name, vers
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



subroutine close_file(unit)
!-----------------------------------------------------------------------
!
! Closes the given unit_number. If the file is already closed, 
! nothing happens. Pretty dramatic, eh?
!

integer, intent(in) :: unit

integer :: ios
logical :: open

inquire (unit=unit, opened=open, iostat=ios)
if ( ios /= 0 ) then
   print *,'Dagnabbit. Cannot inquire about unit # ',unit
   print *,'Error status was ',ios
   print *,'Hoping for the best and continuing.'
endif

if (open) close(unit)

end subroutine close_file

!
!=======================================================================
! End of utilities_mod
!=======================================================================
!
end module utilities_mod
