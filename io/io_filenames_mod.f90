module io_filenames_mod

!> \defgroup io_filenames io_filenames
!> Aim is to store the io filenames
!>  * Restarts
!>  * Diagnostics
!>  * Inflation files
!>
!> Any module can set the filenames here, then state_vector_io_mod
!> can read from this module to get the filenames.
!> Maybe this is a bit lazy, but I just want a way for the 
!> different modules to set the filenames
!>
!> Diagnostic files could have different netcdf variable ids
!> @{

implicit none

private

! These should probably be set and get functions rather than 
! direct access

public :: restart_stub, &
          prior_diagnostic_file, &
          post_diagnostic_file, &
          prior_mean_inf_file, &
          prior_sd_inf_file, &
          post_mean_inf_file, &
          post_sd_inf_file

! How do people name there restart files?
! What about domains?
character(len=256) :: restart_stub = 'restart'

character(len=256) :: prior_diagnostic_file = 'Prior_diag'
character(len=256) :: post_diagnostic_file = 'Posterior_diag'

character(len=256) :: prior_mean_inf_file = 'prior_inf_ic_old_d'
character(len=256) :: prior_sd_inf_file   = 'prior_inf_ic_sd_old_d'
character(len=256) :: post_mean_inf_file  = 'post_inf_ic_new_d'
character(len=256) :: post_sd_inf_file    = 'post_inf_ic_sd_new_d'

contains


!> @}
end module io_filenames_mod