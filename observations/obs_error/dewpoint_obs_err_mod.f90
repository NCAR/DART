! DART software - Copyright UCAR. This open source software is provided
! by UCAR, "as is", without charge, subject to all terms of use at
! http://www.image.ucar.edu/DAReS/DART/DART_download
!
! $Id$

module dewpoint_obs_err_mod

use        types_mod, only : r8
use       meteor_mod, only : rh_and_temp_to_dewpoint, temp_and_dewpoint_to_rh

implicit none
private

public :: dewpt_error_from_rh_and_temp, rh_error_from_dewpt_and_temp

contains


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!   dewpt_error_from_rh_and_temp - computes estimated uncertainty in
!           a dewpoint observation, when the dewpoint is derived from
!           temperature and relative humidity observations
!
!     created Dec. 2008 David Dowell, NCAR/MMM
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

function dewpt_error_from_rh_and_temp(tmpk, rh_in)

! reference:  Lin and Hubbard 2004, Journal of Applied Meteorology, 821-825

real(r8)             :: dewpt_error_from_rh_and_temp
real(r8), intent(in) :: tmpk                   ! temperature (Kelvin)
real(r8), intent(in) :: rh_in                  ! relative humidity (0.00-1.00)

real(r8), parameter  :: rh_error = 0.05_r8     ! guess for instrument + representativeness error
real(r8), parameter  :: t_error = 1.0_r8       ! guess for instrument + representativeness error

real(r8), parameter  :: rh_min = 0.20_r8
real(r8), parameter  :: rh_max = 1.00_r8
real(r8), parameter  :: delta_rh = 0.01_r8     ! perturbation for finite differencing
real(r8), parameter  :: delta_t  = 0.1_r8      ! perturbation for finite differencing
real(r8)             :: rh
real(r8)             :: rh1
real(r8)             :: rh2
real(r8)             :: t1
real(r8)             :: t2
real(r8)             :: td_deriv_rh
real(r8)             :: td_deriv_t


if ( ( rh_in < 0.00_r8 ) .or. ( rh_in > 1.00_r8 ) ) then
  print*,'dewpt_error_from_rh_and_temp:  bad rh ',rh_in
  stop
end if

! The uncertainty in dewpoint derived from temperature and relative humidity grows without bound
! as the relative humidity approaches 0.  A lower bound on the input relative humidity will be
! applied to keep the calculation bounded and to acknowledge uncertainty in the relative humidity
! measurement.

rh = max(rh_in, rh_min)

! finite-difference approximations of derivatives
rh1 = max(rh - delta_rh, rh_min)
rh2 = min(rh + delta_rh, rh_max)
t1 = tmpk - delta_t
t2 = tmpk + delta_t
td_deriv_rh = ( rh_and_temp_to_dewpoint(rh2, tmpk) - rh_and_temp_to_dewpoint(rh1, tmpk) ) / (rh2-rh1)
td_deriv_t = ( rh_and_temp_to_dewpoint(rh, t2) - rh_and_temp_to_dewpoint(rh, t1) ) / (t2-t1)

dewpt_error_from_rh_and_temp = sqrt ( td_deriv_t*td_deriv_t * t_error*t_error &
                                     +td_deriv_rh*td_deriv_rh * rh_error*rh_error)

return
end function dewpt_error_from_rh_and_temp


! GSR  - Added rh error function following above
function rh_error_from_dewpt_and_temp(tmpk, dewpt)

! reference:  Lin and Hubbard 2004, Journal of Applied Meteorology, 821-825

real(r8)             :: rh_error_from_dewpt_and_temp
real(r8), intent(in) :: tmpk                   ! temperature (Kelvin)
real(r8), intent(in) :: dewpt                  ! dewpoint temperature (Kelvin) 

real(r8), parameter  :: td_error = 1.50_r8     ! guess for instrument + representativeness error
real(r8), parameter  :: t_error = 1.0_r8       ! guess for instrument + representativeness error

real(r8), parameter  :: delta_td = 0.1_r8      ! perturbation for finite differencing
real(r8), parameter  :: delta_t  = 0.1_r8      ! perturbation for finite differencing
real(r8)             :: td1
real(r8)             :: td2
real(r8)             :: t1
real(r8)             :: t2
real(r8)             :: rh_deriv_td
real(r8)             :: rh_deriv_t


if ( ( dewpt > tmpk ) ) then
  print*,'rh_error_from_dewpt_and_temp:  bad dewpt ',dewpt, tmpk
  stop
end if

! finite-difference approximations of derivatives
td1 = dewpt - delta_td
td2 = min (dewpt + delta_td, tmpk + delta_t)
t1 = tmpk - delta_t
t2 = tmpk + delta_t
rh_deriv_td = ( temp_and_dewpoint_to_rh(tmpk, td2) - temp_and_dewpoint_to_rh(tmpk, td1) ) / (td2-td1)
rh_deriv_t = ( temp_and_dewpoint_to_rh(t2, dewpt) - temp_and_dewpoint_to_rh(t1, dewpt) ) / (t2-t1)

rh_error_from_dewpt_and_temp = sqrt ( rh_deriv_t*rh_deriv_t * t_error*t_error &
                                     +rh_deriv_td*rh_deriv_td * td_error*td_error)

return
end function rh_error_from_dewpt_and_temp

end module dewpoint_obs_err_mod

! <next few lines under version control, do not edit>
! $URL$
! $Id$
! $Revision$
! $Date$
