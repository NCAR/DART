! DART software - Copyright 2004 - 2013 UCAR. This open source software is
! provided by UCAR, "as is", without charge, subject to all terms of use at
! http://www.image.ucar.edu/DAReS/DART/DART_download
!
! $Id$


!> Aim: collection of assertions for use in test code
module assert_mod

implicit none

public

interface assert_equal
   module procedure assert_equal_real
   module procedure assert_equal_double
   module procedure assert_equal_int
   module procedure assert_equal_int8
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

real, intent(in) :: a, b
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

integer*8, intent(in) :: a, b
character(len=*), intent(in) :: message

if (a /= b) print*, 'FAIL: ',  trim(message),' assertion ', a, '==', b, 'failed'

end subroutine assert_equal_int8

!-------------------------------
! Assert greater
!-------------------------------
subroutine assert_greater_int(a, b, message)

integer, intent(in) :: a, b
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