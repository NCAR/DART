program gswm_test

use types_mod
use gswm_mod,          only : psi, static_model_init
use time_manager_mod
use location_mod
use utilities_mod,     only : get_unit, open_file, close_file, &
                              check_nml_error, file_exist

implicit none


!----------------------------------------------------------------------
! Declare local variables
!----------------------------------------------------------------------

real(r8) :: xi, dec, eta, sfc



!----------------------------------------------------------------------
! Initialize the model
!----------------------------------------------------------------------

call static_model_init

!----------------------------------------------------------------------
! Test the PSI routine ...
!----------------------------------------------------------------------

if ( 1 == 1 ) then

   xi  = 10.0
   dec = 13.0
   eta = 90.0_r8*pi/360.0_r8

   sfc = psi(xi,eta,dec)

   write(*,*)'xi         is ',xi
   write(*,*)'dec        is ',dec
   write(*,*)'eta        is ',eta
   write(*,*)'psi result is ',sfc

endif


!----------------------------------------------------------------------
! Test some other routine ...
!----------------------------------------------------------------------

if ( 0 == 1 ) then


endif

END PROGRAM gswm_test
