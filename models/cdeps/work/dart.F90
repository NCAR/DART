! three required subroutines for ESMX
module dart

implicit none

public :: init, run, finalize

contains

!-----------------
subroutine init()

print*, 'Hello Anh from init'

! example of initializing a dart module
!call assim_tools_init()

end subroutine init

!-----------------
subroutine run

print*, 'Hello Anh from run'

end subroutine run

!-----------------
subroutine finalize

print*, 'Hello Anh from finalize' 

end subroutine finalize

end module dart
