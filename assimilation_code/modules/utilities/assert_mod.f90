! DART software - Copyright UCAR. This open source software is provided
! by UCAR, "as is", without charge, subject to all terms of use at
! http://www.image.ucar.edu/DAReS/DART/DART_download


!> Aim: collection of assertions for use in test code
module assert_mod

use types_mod, only : r8, r4, i8

implicit none

public

interface assert_equal
   module procedure assert_equal_real
   module procedure assert_equal_double
   module procedure assert_equal_int
   module procedure assert_equal_int8
   module procedure assert_equal_strings
   module procedure assert_equal_int_array
   module procedure assert_equal_logical
   module procedure assert_equal_logical_array
end interface

interface assert_not_equal
   module procedure assert_not_equal_real
   module procedure assert_not_equal_int
end interface

interface assert_greater
   module procedure assert_greater_int
end interface assert_greater

contains

!-------------------------------
! Assert equal
!-------------------------------
subroutine assert_equal_real(a, b, message)

real(r4), intent(in) :: a, b
character(len=*), intent(in) :: message

if (a /= b) print*, 'FAIL: ', trim(message),' assertion ', a, '==', b, 'failed'

end subroutine assert_equal_real

!-------------------------------
subroutine assert_equal_double(a, b, message)

double precision, intent(in) :: a, b
character(len=*), intent(in) :: message

if (a /= b) print*, 'FAIL: ',  trim(message),' assertion ', a, '==', b, 'failed'

end subroutine assert_equal_double

!-------------------------------
subroutine assert_equal_int(a, b, message)

integer, intent(in) :: a, b
character(len=*), intent(in) :: message

if (a /= b) print*, 'FAIL: ',  trim(message),' assertion ', a, '==', b, 'failed'

end subroutine assert_equal_int

!-------------------------------
subroutine assert_equal_int8(a, b, message)

integer(i8),        intent(in) :: a, b
character(len=*), intent(in) :: message

if (a /= b) print*, 'FAIL: ',  trim(message),' assertion ', a, '==', b, 'failed'

end subroutine assert_equal_int8

!-------------------------------
subroutine assert_equal_strings(a, b, message)

character(len=*), intent(in) :: a, b
character(len=*), intent(in) :: message

if (trim(a) == "" .and. trim(b) == "") return ! both strings blank

if (trim(a) /= trim(b)) print*, 'FAIL: ',  trim(message), &
           ' assertion "', trim(a), '" == "', trim(b), '" failed'

end subroutine assert_equal_strings

!-------------------------------
subroutine assert_equal_int_array(a, b, message)

integer, dimension(:), intent(in) :: a, b
character(len=*),      intent(in) :: message

integer :: i

if (size(a) /= size(b)) print*, 'FAIL: ',  trim(message), &
           ' array assertion failed because of unequal lengths'

if (any(a /= b)) then

   print*, 'FAIL: ', trim(message), ' array assertion failed.'

   if (size(a) < 100) then
      do i = 1,size(a)
         write(*,'(''       element('',i3,'') '',i3,'' ?==? '',i3)')i, a(i), b(i)
      enddo
   else
      print*, 'arrays too long to concisely specify where/how failed.'
   endif
endif

end subroutine assert_equal_int_array

!-------------------------------
subroutine assert_equal_logical(a, b, message)

logical,          intent(in) :: a, b
character(len=*), intent(in) :: message

if (a .neqv. b) print*, 'FAIL: ',  trim(message), &
           ' assertion ', a, '==', b, ' failed'

end subroutine assert_equal_logical


!-------------------------------
subroutine assert_equal_logical_array(a, b, message)

logical, dimension(:), intent(in) :: a, b
character(len=*),      intent(in) :: message

integer :: i

if (size(a) /= size(b)) print*, 'FAIL: ',  trim(message), &
           ' array assertion failed because of unequal lengths'

if (any(a .neqv. b)) then

   print*, 'FAIL: ', trim(message), ' array assertion failed.'

   if (size(a) < 100) then
      do i = 1,size(a)
         write(*,'(''       element('',i3,'') '',L1,'' ?==? '',L1)')i, a(i), b(i)
      enddo
   else
      print*, 'arrays too long to concisely specify where/how failed.'
   endif
endif

end subroutine assert_equal_logical_array

!-------------------------------
! Assert greater
!-------------------------------
subroutine assert_greater_int(a, b, message)

integer,          intent(in) :: a, b
character(len=*), intent(in) :: message

if (a <= b) print*, 'FAIL: ',  trim(message),' assertion ', a, ' >', b, 'failed'

end subroutine assert_greater_int

!-------------------------------
! Assert not equal
!-------------------------------
subroutine assert_not_equal_real(a, b, message)

real, intent(in) :: a, b
character(len=*), intent(in) :: message

if (a == b) print*, 'FAIL: ',  trim(message),' assertion ', a, '/=', b, 'failed'

end subroutine assert_not_equal_real

!-------------------------------
subroutine assert_not_equal_int(a, b, message)

integer, intent(in) :: a, b
character(len=*), intent(in) :: message

if (a == b) print*, 'FAIL: ',  trim(message),' assertion ', a, '/=', b, 'failed'

end subroutine assert_not_equal_int


end module assert_mod

