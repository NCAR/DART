module transform_state_mod

use types_mod,     only : varnamelength
use utilities_mod, only : open_file, close_file, find_namelist_in_file, &
                          check_namelist_read

implicit none

integer :: io, iunit
integer :: nfaces, nghost_rows, nghost_columns
integer :: nargs, iarg
character(len=varnamelength), dimension(:), allocatable :: args
character(len=varnamelength) :: padded_stringy

namelist /transform_state_nml/ nfaces, nghost_rows, nghost_columns

contains

subroutine read_namelist_and_command_line_arguments()

   call find_namelist_in_file('input.nml', 'transform_state_nml', iunit)
   read(iunit, nml = transform_state_nml, iostat = io)
   call check_namelist_read(iunit, io, 'transform_state_nml')

   print *, nfaces
   print *, nghost_rows
   print *, nghost_columns

   nargs = command_argument_count()
   allocate(args(nargs))

   print *, 'nargs', nargs

   do iarg = 1, nargs
      call get_command_argument(iarg, args(iarg))
      padded_stringy = zfill(args(iarg), 4)
      print *, padded_stringy
   end do

end subroutine read_namelist_and_command_line_arguments

function str2int(string) result(int)

   character(len=varnamelength), intent(in) :: string
   
   integer :: length_of_string
   integer :: int
   
   length_of_string = len_trim(adjustl(string))
   if (length_of_string > 9) then
      print *, 'Error: integer stringing length is greater than 9 digits long.'
      stop
   end if 
   
   read(string,*) int
   
end function str2int

function int2str(int) result(string)

   integer, intent(in) :: int
   character(len=varnamelength) :: string

   write(string,'(I0)') int
   string = trim(string)

end function int2str

function zfill(string, desired_length) result(filled_string)

   character(len=*), intent(in) :: string
   integer, intent(in) :: desired_length

   character(len=varnamelength) :: filled_string
   integer :: length_of_string
   integer :: string_index, difference_of_string_lengths

   filled_string = ''
   length_of_string = len_trim(string)
   difference_of_string_lengths = desired_length - length_of_string

   if (difference_of_string_lengths < 0) then
      print *, 'Error: input string is longer than the desired output string.'
      stop
   else if (difference_of_string_lengths > 0) then
      do string_index = 1, difference_of_string_lengths
         filled_string(string_index:string_index) = '0'
      end do
   end if

   filled_string(difference_of_string_lengths+1:desired_length) = trim(string)

end function zfill

end module transform_state_mod