! Data Assimilation Research Testbed -- DART
! Copyright 2004, Data Assimilation Initiative, University Corporation for Atmospheric Research
! Licensed under the GPL -- www.gpl.org/licenses/gpl.html

module cov_cutoff_mod

! <next four lines automatically updated by CVS, do not edit>
! $Source$ 
! $Revision$ 
! $Date$ 
! $Author$ 
!

use     types_mod, only : r8
use utilities_mod, only : file_exist, open_file, check_nml_error, &
                          close_file, error_handler, FATAL

implicit none
private

public :: comp_cov_factor

! CVS Generated file description for error handling, do not edit
character(len=128) :: &
source   = "$Source$", &
revision = "$Revision$", &
revdate  = "$Date$"


!============================================================================

!---- namelist with default values
logical :: namelist_initialized = .false.

integer :: select_localization = 1
! Value 1 selects default Gaspari-Cohn cutoff
! Value 2 selects boxcar
! Value 3 selects ramped boxcar

namelist / cov_cutoff_nml / select_localization

!============================================================================

contains

!======================================================================



function comp_cov_factor(z_in, c)
!----------------------------------------------------------------------
! function comp_cov_factor(z_in, c)
!
! Computes a covariance cutoff function from Gaspari and Cohn
! QJRMS, 125, 723-757.  (their eqn. 4.10)
!
! z_in is the distance while c is the cutoff distance. 
! For distances greater than 2c, the cov_factor returned goes to 0.

! Other ramping shapes are also available and can be selected by a namelist
! parameter. At present, these include a boxcar with the given halfwidth
! and a ramp in which the weight is set at 1.0 up to the half-width 
! distance and then decreases linearly to 0 at twice the half-width 
! distance.

implicit none

real(r8), intent(in) :: z_in, c
real(r8)             :: comp_cov_factor

real(r8) :: z, r
integer  :: iunit, ierr, io

z = abs(z_in)

!--------------------------------------------------------
! Initialize namelist if not already done
if(.not. namelist_initialized) then

   namelist_initialized = .true.

   if(file_exist('input.nml')) then
      iunit = open_file(file = 'input.nml', action = 'read')
      ierr = 1

      READBLOCK: do while(ierr /= 0)
         read(iunit, nml = cov_cutoff_nml, iostat = io)
         if ( io < 0 ) exit READBLOCK          ! end-of-file
         ierr = check_nml_error(io, 'cov_cutoff_nml')
      enddo READBLOCK

      call close_file(iunit)
   endif
endif
!---------------------------------------------------------

if(select_localization == 3) then ! Ramped localization

   if(z >= 2.0_r8 * c) then
      comp_cov_factor = 0.0_r8
   else if(z >= c .and. z < 2.0_r8 * c) then
      comp_cov_factor = (2.0_r8 * c - z) / c
   else
      comp_cov_factor = 1.0_r8
   endif

else if(select_localization == 2) then ! BOXCAR localization

   if(z < 2 * c) then
      comp_cov_factor = 1.0_r8
   else
      comp_cov_factor = 0.0_r8
   endif

else if(select_localization == 1) then ! Standard Gaspari Cohn localization

   if( z >= c*2.0_r8 ) then

      comp_cov_factor = 0.0_r8

   else if( z >= c .and. z < c*2.0_r8 ) then

      r = z / c
      comp_cov_factor = r**5 / 12.0_r8  -  &
                        r**4 / 2.0_r8   +  &
                        r**3 * 5.0_r8 / 8.0_r8 + &
                        r**2 * 5.0_r8 / 3.0_r8 - &
                        5.0_r8*r + 4.0_r8 - (c * 2.0_r8) / (3.0_r8 * z) 
   else
      r = z / c
      comp_cov_factor = r**5 * (-0.25_r8 ) + &
                        r**4 / 2.0_r8 +              &
                        r**3 * 5.0_r8/8.0_r8 -       &
                        r**2 * 5.0_r8/3.0_r8 + 1.0_r8
   endif

else ! Otherwise namelist parameter is illegal; this is an error

     call error_handler(FATAL,'comp_cov_factor', &
              'Illegal value of "select_localization" in cov_cutoff_mod namelist', &
               source, revision, revdate )

endif

end function comp_cov_factor

end module cov_cutoff_mod
