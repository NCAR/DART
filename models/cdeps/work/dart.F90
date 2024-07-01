! three required subroutines for ESMX
module dart

use mpi_utilities_mod, only: initialize_mpi_utilities, &
                             finalize_mpi_utilities, &
                             my_task_id

use dart_comp_nuopc, only : SetServices, &
                            SetVM

implicit none
public :: init, run, finalize

contains

!-----------------
subroutine init()

call initialize_mpi_utilities()
print*, 'Hello Anh from init', my_task_id()

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
call finalize_mpi_utilities()

end subroutine finalize

end module dart
