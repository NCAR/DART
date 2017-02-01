! DART software - Copyright UCAR. This open source software is provided
! by UCAR, "as is", without charge, subject to all terms of use at
! http://www.image.ucar.edu/DAReS/DART/DART_download
!
! $Id$

module obs_err_mod

use        types_mod, only : r8, missing_r8

implicit none
private

integer, parameter :: nobs_level = 15
real(r8)           :: obs_prs(nobs_level)

! this array is ordered from top of atm down to surface
data obs_prs/ 10.0_r8,  20.0_r8,  30.0_r8,  50.0_r8,  70.0_r8, & 
             100.0_r8, 150.0_r8, 200.0_r8, 250.0_r8, 300.0_r8, & 
             400.0_r8, 500.0_r8, 700.0_r8, 850.0_r8, 1050.0_r8/


public :: fixed_marine_pres_error,      &
          fixed_marine_rel_hum_error,   &
          fixed_marine_temp_error,      &
          fixed_marine_wind_error,      &
          moving_marine_pres_error,     &
          moving_marine_rel_hum_error,  &
          moving_marine_temp_error,     &
          moving_marine_wind_error

public :: land_pres_error,              &
          land_rel_hum_error,           &
          land_temp_error,              &
          land_wind_error

public :: acars_rel_hum_error,          &
          acars_temp_error,             &
          acars_wind_error

public :: sat_wind_error,               &
          sat_wv_wind_error

public :: rawin_pres_error,             &
          rawin_rel_hum_error,          &
          rawin_temp_error,             &
          rawin_wind_error

public :: drop_pres_error,              &
          drop_rel_hum_error,           &
          drop_temp_error,              &
          drop_wind_error

public :: prof_wind_error        

contains


function acars_rel_hum_error(pres, tmpk, rh)

real(r8), intent(in) :: pres, tmpk, rh   !  (mb)

real(r8) :: acars_rel_hum_error

acars_rel_hum_error = missing_r8

return
end function acars_rel_hum_error


function acars_temp_error(pres)

real(r8), intent(in) :: pres  !  (mb)

integer  :: k0
real(r8) :: obs_err(nobs_level), wght, acars_temp_error

data obs_err/2.2_r8, 2.0_r8, 1.8_r8, 1.6_r8, 1.5_r8, &
             1.4_r8, 1.4_r8, 1.4_r8, 1.3_r8, 1.3_r8, &
             1.2_r8, 1.2_r8, 1.2_r8, 1.3_r8, 1.4_r8/

call find_pressure_level_weight(pres, k0, wght)
acars_temp_error = wght * obs_err(k0) + (1.0_r8 - wght) * obs_err(k0+1)
acars_temp_error = dble(nint(acars_temp_error * 10.0_r8)) * 0.1_r8

return
end function acars_temp_error


function acars_wind_error(pres)

real(r8), intent(in) :: pres  !  (mb)

integer  :: k0
real(r8) :: obs_err(nobs_level), wght, acars_wind_error

data obs_err/4.0_r8, 4.0_r8, 4.0_r8, 4.0_r8, 4.0_r8, &
             4.0_r8, 4.0_r8, 4.0_r8, 4.0_r8, 4.0_r8, &
             4.0_r8, 3.5_r8, 3.0_r8, 2.5_r8, 2.5_r8/

call find_pressure_level_weight(pres, k0, wght)
acars_wind_error = wght * obs_err(k0) + (1.0_r8 - wght) * obs_err(k0+1)
acars_wind_error = dble(nint(acars_wind_error * 10.0_r8)) * 0.1_r8

return
end function acars_wind_error


function land_pres_error(pres)

real(r8), intent(in) :: pres  !  (mb)

integer  :: k0
real(r8) :: obs_err(nobs_level), wght, land_pres_error

data obs_err/114.2_r8, 96.0_r8, 69.8_r8, 59.3_r8, 50.3_r8, & 
              39.4_r8, 32.4_r8, 27.7_r8, 25.4_r8, 18.8_r8, & 
              14.9_r8, 12.1_r8,  8.6_r8,  8.0_r8,  7.0_r8/

call find_pressure_level_weight(pres, k0, wght)
land_pres_error = wght * obs_err(k0) + (1.0_r8 - wght) * obs_err(k0+1)
land_pres_error = land_pres_error * 1.225_r8 * 0.25_r8
land_pres_error = dble(nint(land_pres_error * 10.0_r8)) * 0.1_r8

return
end function land_pres_error


function land_rel_hum_error(pres, tmpk, rh)

real(r8), intent(in) :: pres, tmpk, rh

real(r8) :: land_rel_hum_error

if ( RH < 0.2_r8 ) then
  land_rel_hum_error = 0.23_r8
elseif ( tmpk < 233.0_r8 ) then
  land_rel_hum_error = 0.28_r8
else
  land_rel_hum_error = 0.13_r8
end if

return
end function land_rel_hum_error


function land_temp_error(pres)

real(r8), intent(in) :: pres  !  (mb)

integer :: k0
real(r8) :: obs_err(nobs_level), wght, land_temp_error

data obs_err/2.5_r8, 2.5_r8, 2.5_r8, 2.4_r8, 2.2_r8, &
             2.0_r8, 1.9_r8, 1.8_r8, 1.8_r8, 1.5_r8, &
             1.3_r8, 1.2_r8, 1.3_r8, 1.5_r8, 2.0_r8/

call find_pressure_level_weight(pres, k0, wght)
land_temp_error = wght * obs_err(k0) + (1.0_r8 - wght) * obs_err(k0+1)
land_temp_error = dble(nint(land_temp_error * 10.0_r8)) * 0.1_r8

return
end function land_temp_error


function land_wind_error(pres)

real(r8), intent(in) :: pres  !  (mb)

integer  :: k0
real(r8) :: obs_err(nobs_level), wght, land_wind_error

data obs_err/3.0_r8, 2.5_r8, 2.0_r8, 2.0_r8, 2.0_r8, &
             2.2_r8, 2.4_r8, 3.2_r8, 3.2_r8, 3.8_r8, &
             3.6_r8, 3.4_r8, 3.0_r8, 3.0_r8, 3.0_r8/

call find_pressure_level_weight(pres, k0, wght)
land_wind_error = wght * obs_err(k0) + (1.0_r8 - wght) * obs_err(k0+1)
land_wind_error = dble(nint(land_wind_error * 10.0_r8)) * 0.1_r8

return
end function land_wind_error


function fixed_marine_pres_error(pres)

real(r8), intent(in) :: pres  !  (mb)

real(r8) :: fixed_marine_pres_error

fixed_marine_pres_error = land_pres_error(pres)

return
end function fixed_marine_pres_error


function fixed_marine_rel_hum_error(pres, tmpk, rh)

real(r8), intent(in) :: pres, tmpk, rh

real(r8) :: fixed_marine_rel_hum_error

fixed_marine_rel_hum_error = land_rel_hum_error(pres,tmpk,rh)

return
end function fixed_marine_rel_hum_error


function fixed_marine_temp_error(pres)

real(r8), intent(in) :: pres  !  (mb)

real(r8) :: fixed_marine_temp_error

fixed_marine_temp_error = land_temp_error(pres)

return
end function fixed_marine_temp_error


function fixed_marine_wind_error(pres)

real(r8), intent(in) :: pres  !  (mb)

real(r8) :: fixed_marine_wind_error

fixed_marine_wind_error = land_wind_error(pres)

return
end function fixed_marine_wind_error


function moving_marine_pres_error(pres)

real(r8), intent(in) :: pres  !  (mb)

integer  :: k0
real(r8) :: obs_err(nobs_level), wght, moving_marine_pres_error

data obs_err/228.4_r8, 192.0_r8, 139.5_r8, 118.6_r8, 100.6_r8, & 
              78.8_r8,  64.8_r8,  55.4_r8,  50.8_r8,  37.6_r8, & 
              29.8_r8,  24.2_r8,  17.2_r8,  11.5_r8,  11.5_r8/

call find_pressure_level_weight(pres, k0, wght)
moving_marine_pres_error = wght * obs_err(k0) + (1.0_r8 - wght) * obs_err(k0+1)
moving_marine_pres_error = moving_marine_pres_error * 1.225_r8 * 0.25_r8
moving_marine_pres_error = dble(nint(moving_marine_pres_error * 10.0_r8)) * 0.1_r8

return
end function moving_marine_pres_error


function moving_marine_rel_hum_error(pres, tmpk, rh)

real(r8), intent(in) :: pres, tmpk, rh

real(r8) :: moving_marine_rel_hum_error

moving_marine_rel_hum_error = missing_r8

return
end function moving_marine_rel_hum_error


function moving_marine_temp_error(pres)

real(r8), intent(in) :: pres  !  (mb)

integer  :: k0
real(r8) :: obs_err(nobs_level), wght, moving_marine_temp_error

data obs_err/2.5_r8, 2.5_r8, 2.5_r8, 2.4_r8, 2.2_r8, &
             2.0_r8, 1.9_r8, 1.8_r8, 1.8_r8, 1.5_r8, &
             1.3_r8, 1.2_r8, 1.3_r8, 1.5_r8, 1.8_r8/

call find_pressure_level_weight(pres, k0, wght)
moving_marine_temp_error = wght * obs_err(k0) + (1.0_r8 - wght) * obs_err(k0+1)
moving_marine_temp_error = dble(nint(moving_marine_temp_error * 10.0_r8)) * 0.1_r8

return
end function moving_marine_temp_error


function moving_marine_wind_error(pres)

real(r8), intent(in) :: pres  !  (mb)

integer  :: k0
real(r8) :: obs_err(nobs_level), wght, moving_marine_wind_error

data obs_err/3.0_r8, 2.5_r8, 2.0_r8, 2.0_r8, 2.0_r8, &
             2.2_r8, 2.4_r8, 3.2_r8, 3.2_r8, 3.8_r8, &
             3.6_r8, 3.4_r8, 2.5_r8, 2.4_r8, 2.4_r8/

call find_pressure_level_weight(pres, k0, wght)
moving_marine_wind_error = wght * obs_err(k0) + (1.0_r8 - wght) * obs_err(k0+1)
moving_marine_wind_error = dble(nint(moving_marine_wind_error * 10.0_r8)) * 0.1_r8

return
end function moving_marine_wind_error


function sat_wind_error(pres)

real(r8), intent(in) :: pres  !  (mb)

integer  :: k0
real(r8) :: obs_err(nobs_level), wght, sat_wind_error

data obs_err/5.7_r8, 5.0_r8, 5.0_r8, 5.0_r8, 5.0_r8, &
             5.0_r8, 5.0_r8, 5.0_r8, 5.0_r8, 5.0_r8, &
             4.3_r8, 3.5_r8, 2.0_r8, 2.0_r8, 2.0_r8/

call find_pressure_level_weight(pres, k0, wght)
sat_wind_error = wght * obs_err(k0) + (1.0_r8 - wght) * obs_err(k0+1)
sat_wind_error = dble(nint(sat_wind_error * 10.0_r8)) * 0.1_r8

return
end function sat_wind_error


function sat_wv_wind_error(pres)

real(r8), intent(in) :: pres  !  (mb)

real(r8) :: sat_wv_wind_error

sat_wv_wind_error = sat_wind_error(pres)

return
end function sat_wv_wind_error


function rawin_pres_error(pres)

real(r8), intent(in) :: pres  !  (mb)

integer  :: k0
real(r8) :: obs_err(nobs_level), wght, rawin_pres_error

data obs_err/40.0_r8, 32.0_r8, 25.0_r8, 22.5_r8, 19.5_r8, & 
             18.1_r8, 15.2_r8, 13.2_r8, 11.8_r8, 10.7_r8, & 
              9.8_r8,  8.4_r8,  5.2_r8,  4.4_r8,  4.3_r8/

call find_pressure_level_weight(pres, k0, wght)
rawin_pres_error = wght * obs_err(k0) + (1.0_r8 - wght) * obs_err(k0+1)
rawin_pres_error = rawin_pres_error * 1.225_r8 * 0.25_r8
rawin_pres_error = dble(nint(rawin_pres_error * 10.0_r8)) * 0.1_r8

return
end function rawin_pres_error


function rawin_rel_hum_error(pres, tmpk, rh)

real(r8), intent(in) :: pres, tmpk, rh  !  (mb)

real(r8) :: rawin_rel_hum_error

if ( rh < 0.2_r8 ) then
  rawin_rel_hum_error = 0.23
else if ( tmpk < 233.0_r8 ) then
  rawin_rel_hum_error = 0.28
else
  rawin_rel_hum_error = 0.17
end if

return
end function rawin_rel_hum_error


function rawin_temp_error(pres)

real(r8), intent(in) :: pres  !  (mb)

integer  :: k0
real(r8) :: obs_err(nobs_level), wght, rawin_temp_error

data obs_err/2.5_r8, 2.2_r8, 2.0_r8, 1.9_r8, 1.8_r8, & 
             1.7_r8, 1.6_r8, 1.5_r8, 1.5_r8, 1.4_r8, & 
             1.2_r8, 1.2_r8, 1.3_r8, 1.5_r8, 1.7_r8/

call find_pressure_level_weight(pres, k0, wght)
rawin_temp_error = wght * obs_err(k0) + (1.0_r8 - wght) * obs_err(k0+1)
rawin_temp_error = dble(nint(rawin_temp_error * 10.0_r8)) * 0.1_r8

return
end function rawin_temp_error


function rawin_wind_error(pres)

real(r8), intent(in) :: pres  !  (mb)

integer  :: k0
real(r8) :: obs_err(nobs_level), wght, rawin_wind_error

data obs_err/4.5_r8, 3.6_r8, 3.3_r8, 3.2_r8, 3.2_r8, &
             3.3_r8, 3.4_r8, 3.5_r8, 3.5_r8, 3.7_r8, &
             3.5_r8, 3.0_r8, 2.5_r8, 2.3_r8, 2.3_r8/

call find_pressure_level_weight(pres, k0, wght)
rawin_wind_error = wght * obs_err(k0) + (1.0_r8 - wght) * obs_err(k0+1)
rawin_wind_error = dble(nint(rawin_wind_error * 10.0_r8)) * 0.1_r8

return
end function rawin_wind_error


function drop_pres_error(pres)

real(r8), intent(in) :: pres  !  (mb)

real(r8) :: drop_pres_error

drop_pres_error = rawin_pres_error(pres)

return
end function drop_pres_error


function drop_rel_hum_error(pres, tmpk, rh)

real(r8), intent(in) :: pres, tmpk, rh

real(r8) :: drop_rel_hum_error

drop_rel_hum_error = rawin_rel_hum_error(pres, tmpk, rh)

return
end function drop_rel_hum_error


function drop_temp_error(pres)

real(r8), intent(in) :: pres  !  (mb)

real(r8) :: drop_temp_error

drop_temp_error = rawin_temp_error(pres)

return
end function drop_temp_error


function drop_wind_error(pres)

real(r8), intent(in) :: pres  !  (mb)

real(r8) :: drop_wind_error

drop_wind_error = rawin_wind_error(pres)

return
end function drop_wind_error


function prof_wind_error(pres)

real(r8), intent(in) :: pres  !  (mb)

prof_wind_error = missing_r8

return
end function prof_wind_error


! uses global obs_prs() array.  assumes order of pressure values in
! the array starts at the top and ends at the surface.
subroutine find_pressure_level_weight(pres, zloc, wght)

real(r8), intent(in)  :: pres
integer, intent(out)  :: zloc
real(r8), intent(out) :: wght

integer :: k

if ( pres < obs_prs(1) .or. pres > obs_prs(nobs_level) ) then
  print*,'bad pressure level',pres
  stop
end if

do k = 2, nobs_level

  if ( pres <= obs_prs(k) )  then
    zloc = k - 1
    wght = (log(obs_prs(k)) - log(pres)) / &
           (log(obs_prs(k)) - log(obs_prs(k-1)))
    exit
  end if

end do

return
end subroutine find_pressure_level_weight

end module obs_err_mod

! <next few lines under version control, do not edit>
! $URL$
! $Id$
! $Revision$
! $Date$
