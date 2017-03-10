! DART software - Copyright UCAR. This open source software is provided
! by UCAR, "as is", without charge, subject to all terms of use at
! http://www.image.ucar.edu/DAReS/DART/DART_download
!
! $Id$

module obs_err_mod

use        types_mod, only : r8, missing_r8

implicit none
private

integer,  parameter :: nobs_level = 33
real(r8)           :: obs_prs(nobs_level)

! this array is ordered from top of atm down to surface
data obs_prs/   0.0_r8,    1.0_r8,    2.0_r8,   3.0_r8,   4.0_r8, &
                5.0_r8,   10.0_r8,   20.0_r8,  30.0_r8,  40.0_r8, &
               50.0_r8,   75.0_r8,  100.0_r8, 150.0_r8, 200.0_r8, &
              250.0_r8,  300.0_r8,  350.0_r8, 400.0_r8, 450.0_r8, &
              500.0_r8,  550.0_r8,  600.0_r8, 650.0_r8, 700.0_r8, &
              750.0_r8,  800.0_r8,  850.0_r8, 900.0_r8, 950.0_r8, &
             1000.0_r8, 1050.0_r8, 1100.0_r8/

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

public :: metar_pres_error,              &
          metar_rel_hum_error,           &
          metar_temp_error,              &
          metar_wind_error

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

acars_rel_hum_error = 0.2_r8

return
end function acars_rel_hum_error


function acars_temp_error(pres)

real(r8), intent(in) :: pres  !  (mb)

integer  :: k0
real(r8) :: obs_err(nobs_level), wght, acars_temp_error
! this array is ordered from top of atm down to surface

data obs_err/1.0_r8, 1.0_r8, 1.0_r8, 1.0_r8, 1.0_r8, 1.0_r8, &
             1.0_r8, 1.0_r8, 1.0_r8, 1.0_r8, 1.0_r8, 1.0_r8, &
             1.0_r8, 1.0_r8, 1.0_r8, 1.0_r8, 1.0_r8, 1.0_r8, &
             1.0_r8, 1.0_r8, 1.0_r8, 1.0_r8, 1.0_r8, 1.0_r8, &
             1.0_r8, 1.0_r8, 1.0_r8, 1.11_r8, 1.24_r8, 1.35_r8, &
             1.47_r8, 1.47_r8, 1.47_r8/
!data obs_prs/   0.0_r8,    1.0_r8,    2.0_r8,   3.0_r8,    4.0_r8,   5.0_r8,
!               10.0_r8,   20.0_r8,   30.0_r8,  40.0_r8,   50.0_r8,  75.0_r8,
!              100.0_r8,  150.0_r8,  200.0_r8, 250.0_r8,  300.0_r8, 350.0_r8, 
!              400.0_r8,  450.0_r8,  500.0_r8, 550.0_r8,  600.0_r8, 650.0_r8,  
!              700.0_r8,  750.0_r8,  800.0_r8, 850.0_r8,  900.0_r8, 950.0_r8, &
!             1000.0_r8, 1050.0_r8, 1100.0_r8/


call find_pressure_level_weight(pres, k0, wght)
acars_temp_error = wght * obs_err(k0) + (1.0_r8 - wght) * obs_err(k0+1)
acars_temp_error = dble(nint(acars_temp_error * 10.0_r8)) * 0.1_r8

return
end function acars_temp_error


function acars_wind_error(pres)

real(r8), intent(in) :: pres  !  (mb)

real(r8) :: acars_wind_error

acars_wind_error = 2.5_r8

return
end function acars_wind_error


function land_pres_error(pres)

real(r8), intent(in) :: pres  !  (mb)

real(r8) :: land_pres_error

if ( pres >= 600.0_r8 ) then
  land_pres_error = 1.0_r8
else
  land_pres_error = missing_r8
end if

return
end function land_pres_error


function metar_pres_error(pres)

real(r8), intent(in) :: pres  !  (mb)

real(r8) :: metar_pres_error

if ( pres >= 600.0_r8 ) then
  metar_pres_error = 1.0_r8
else
  metar_pres_error = missing_r8
end if

return
end function metar_pres_error


function land_rel_hum_error(pres, tmpk, rh)

real(r8), intent(in) :: pres, tmpk, rh

real(r8) :: land_rel_hum_error

land_rel_hum_error = 0.2_r8

return
end function land_rel_hum_error


function metar_rel_hum_error(pres, tmpk, rh)

real(r8), intent(in) :: pres, tmpk, rh

real(r8) :: metar_rel_hum_error

metar_rel_hum_error = 0.2_r8

return
end function metar_rel_hum_error


function land_temp_error(pres)

real(r8), intent(in) :: pres  !  (mb)

real(r8) :: land_temp_error

land_temp_error = 2.5_r8

return
end function land_temp_error


function metar_temp_error(pres)

real(r8), intent(in) :: pres  !  (mb)

real(r8) :: metar_temp_error

metar_temp_error = 2.5_r8

return
end function metar_temp_error


function land_wind_error(pres)

real(r8), intent(in) :: pres  !  (mb)

real(r8) :: land_wind_error

land_wind_error = 3.5_r8

return
end function land_wind_error


function metar_wind_error(pres)

real(r8), intent(in) :: pres  !  (mb)

real(r8) :: metar_wind_error

metar_wind_error = 3.5_r8

return
end function metar_wind_error


function fixed_marine_pres_error(pres)

real(r8), intent(in) :: pres  !  (mb)

real(r8) :: fixed_marine_pres_error

if ( pres >= 600.0_r8 ) then
  fixed_marine_pres_error = 1.6_r8
else
  fixed_marine_pres_error = missing_r8
end if

return
end function fixed_marine_pres_error


function fixed_marine_rel_hum_error(pres, tmpk, rh)

real(r8), intent(in) :: pres, tmpk, rh

real(r8) :: fixed_marine_rel_hum_error

fixed_marine_rel_hum_error = 0.2_r8

return
end function fixed_marine_rel_hum_error


function fixed_marine_temp_error(pres)

real(r8), intent(in) :: pres  !  (mb)

real(r8) :: fixed_marine_temp_error

fixed_marine_temp_error = 2.5_r8

return
end function fixed_marine_temp_error


function fixed_marine_wind_error(pres)

real(r8), intent(in) :: pres  !  (mb)

real(r8) :: fixed_marine_wind_error

fixed_marine_wind_error = 2.5_r8

return
end function fixed_marine_wind_error


function moving_marine_pres_error(pres)

real(r8), intent(in) :: pres  !  (mb)

real(r8) :: moving_marine_pres_error

moving_marine_pres_error = fixed_marine_pres_error(pres)

return
end function moving_marine_pres_error


function moving_marine_rel_hum_error(pres, tmpk, rh)

real(r8), intent(in) :: pres, tmpk, rh

real(r8) :: moving_marine_rel_hum_error

moving_marine_rel_hum_error = fixed_marine_rel_hum_error(pres, tmpk, rh)

return
end function moving_marine_rel_hum_error


function moving_marine_temp_error(pres)

real(r8), intent(in) :: pres  !  (mb)

real(r8) :: moving_marine_temp_error

moving_marine_temp_error = fixed_marine_temp_error(pres)

return
end function moving_marine_temp_error


function moving_marine_wind_error(pres)

real(r8), intent(in) :: pres  !  (mb)

real(r8) :: moving_marine_wind_error

moving_marine_wind_error = fixed_marine_wind_error(pres)

return
end function moving_marine_wind_error


function sat_wind_error(pres)

real(r8), intent(in) :: pres  !  (mb)

integer  :: k0
real(r8) :: obs_err(nobs_level), wght, sat_wind_error

data obs_err/5.0_r8, 5.0_r8, 5.0_r8, 5.0_r8, 5.0_r8, 5.0_r8, &
             5.0_r8, 5.0_r8, 5.0_r8, 5.0_r8, 5.0_r8, 5.0_r8, &
             5.0_r8, 5.0_r8, 5.0_r8, 5.0_r8, 4.6_r8, 4.3_r8, &
             4.0_r8, 3.0_r8, 2.1_r8, 2.0_r8, 2.0_r8, 1.9_r8, &
             1.9_r8, 1.8_r8, 1.8_r8, 1.8_r8, 1.8_r8, 1.8_r8, &
             1.8_r8, 1.8_r8, 1.8_r8/

call find_pressure_level_weight(pres, k0, wght)
sat_wind_error = wght * obs_err(k0) + (1.0_r8 - wght) * obs_err(k0+1)
sat_wind_error = dble(nint(sat_wind_error * 10.0_r8)) * 0.1_r8

return
end function sat_wind_error


function sat_wv_wind_error(pres)

real(r8), intent(in) :: pres  !  (mb)

integer  :: k0
real(r8) :: obs_err(nobs_level), wght, sat_wv_wind_error

data obs_err/7.0_r8, 7.0_r8, 7.0_r8, 7.0_r8, 7.0_r8, 7.0_r8, &
             7.0_r8, 7.0_r8, 7.0_r8, 7.0_r8, 7.0_r8, 7.0_r8, &
             7.0_r8, 7.0_r8, 7.0_r8, 7.0_r8, 6.6_r8, 6.3_r8, &
             6.0_r8, 5.0_r8, 4.1_r8, 4.0_r8, 4.0_r8, 3.9_r8, &
             3.9_r8, 3.8_r8, 3.8_r8, 3.8_r8, 3.8_r8, 3.8_r8, &
             3.8_r8, 3.8_r8, 3.8_r8/

call find_pressure_level_weight(pres, k0, wght)
sat_wv_wind_error = wght * obs_err(k0) + (1.0_r8 - wght) * obs_err(k0+1)
sat_wv_wind_error = dble(nint(sat_wv_wind_error * 10.0_r8)) * 0.1_r8

return
end function sat_wv_wind_error


function rawin_pres_error(pres)

real(r8), intent(in) :: pres  !  (mb)

real(r8) :: rawin_pres_error

rawin_pres_error = 2.0_r8

return
end function rawin_pres_error


function rawin_rel_hum_error(pres, tmpk, rh)

real(r8), intent(in) :: pres, tmpk, rh

real(r8) :: rawin_rel_hum_error

rawin_rel_hum_error = 0.2_r8

return
end function rawin_rel_hum_error


function rawin_temp_error(pres)

real(r8), intent(in) :: pres  !  (mb)

integer  :: k0
real(r8) :: obs_err(nobs_level), wght, rawin_temp_error

data obs_err/1.5_r8, 1.5_r8, 1.5_r8, 1.5_r8, 1.5_r8, 1.5_r8, &
             1.5_r8, 1.25_r8, 1.0_r8, 0.95_r8, 0.9_r8, 0.8_r8, &
             0.8_r8, 1.0_r8, 1.2_r8, 1.2_r8, 0.9_r8, 0.8_r8, &
             0.8_r8, 0.8_r8, 0.8_r8, 0.8_r8, 0.8_r8, 0.8_r8, &
             0.8_r8, 0.8_r8, 0.8_r8, 0.8_r8, 0.9_r8, 1.1_r8, &
             1.2_r8, 1.2_r8, 1.2_r8/

call find_pressure_level_weight(pres, k0, wght)
rawin_temp_error = wght * obs_err(k0) + (1.0_r8 - wght) * obs_err(k0+1)
rawin_temp_error = dble(nint(rawin_temp_error * 10.0_r8)) * 0.1_r8

return
end function rawin_temp_error


function rawin_wind_error(pres)

real(r8), intent(in) :: pres  !  (mb)

integer  :: k0
real(r8) :: obs_err(nobs_level), wght, rawin_wind_error

data obs_err/2.1_r8, 2.1_r8, 2.1_r8, 2.1_r8, 2.1_r8, 2.1_r8, &
             2.1_r8, 2.1_r8, 2.1_r8, 2.1_r8, 2.1_r8, 2.1_r8, &
             2.1_r8, 2.4_r8, 2.7_r8, 3.2_r8, 3.0_r8, 2.8_r8, &
             2.6_r8, 2.3_r8, 2.1_r8, 2.0_r8, 1.9_r8, 1.8_r8, &
             1.6_r8, 1.6_r8, 1.6_r8, 1.5_r8, 1.5_r8, 1.5_r8, &
             1.4_r8, 1.4_r8, 1.4_r8/
!data obs_prs/   0.0_r8,    1.0_r8,    2.0_r8,   3.0_r8,    4.0_r8,   5.0_r8,
!               10.0_r8,   20.0_r8,   30.0_r8,  40.0_r8,   50.0_r8,  75.0_r8,
!              100.0_r8,  150.0_r8,  200.0_r8, 250.0_r8,  300.0_r8, 350.0_r8, 
!              400.0_r8,  450.0_r8,  500.0_r8, 550.0_r8,  600.0_r8, 650.0_r8,  
!              700.0_r8,  750.0_r8,  800.0_r8, 850.0_r8,  900.0_r8, 950.0_r8, &
!             1000.0_r8, 1050.0_r8, 1100.0_r8/

call find_pressure_level_weight(pres, k0, wght)
rawin_wind_error = wght * obs_err(k0) + (1.0_r8 - wght) * obs_err(k0+1)
rawin_wind_error = dble(nint(rawin_wind_error * 10.0_r8)) * 0.1_r8

return
end function rawin_wind_error


function drop_pres_error(pres)

real(r8), intent(in) :: pres  !  (mb)

real(r8) :: drop_pres_error

drop_pres_error = 2.0_r8

return
end function drop_pres_error


function drop_rel_hum_error(pres, tmpk, rh)

real(r8), intent(in) :: pres, tmpk, rh   !  (mb)

real(r8) :: drop_rel_hum_error

drop_rel_hum_error = 0.2_r8

return
end function drop_rel_hum_error


function drop_temp_error(pres)

real(r8), intent(in) :: pres  !  (mb)

integer  :: k0
real(r8) :: obs_err(nobs_level), wght, drop_temp_error

data obs_err/1.5_r8, 1.5_r8, 1.5_r8, 1.5_r8, 1.5_r8, 1.5_r8, &
             1.5_r8, 1.25_r8, 1.0_r8, 0.95_r8, 0.9_r8, 0.8_r8, &
             0.8_r8, 1.0_r8, 1.2_r8, 1.2_r8, 0.9_r8, 0.8_r8, &
             0.8_r8, 0.8_r8, 0.8_r8, 0.8_r8, 0.8_r8, 0.8_r8, &
             0.8_r8, 0.8_r8, 0.8_r8, 0.8_r8, 0.9_r8, 1.1_r8, &
             1.2_r8, 1.2_r8, 1.2_r8/

call find_pressure_level_weight(pres, k0, wght)
drop_temp_error = wght * obs_err(k0) + (1.0_r8 - wght) * obs_err(k0+1)
drop_temp_error = dble(nint(drop_temp_error * 10.0_r8)) * 0.1_r8

return
end function drop_temp_error


function drop_wind_error(pres)

real(r8), intent(in) :: pres  !  (mb)

integer  :: k0
real(r8) :: obs_err(nobs_level), wght, drop_wind_error

data obs_err/2.7_r8, 2.7_r8, 2.7_r8, 2.7_r8, 2.7_r8, 2.7_r8, &
             2.7_r8, 2.7_r8, 2.7_r8, 2.7_r8, 2.7_r8, 2.6_r8, &
             2.5_r8, 2.73_r8, 2.95_r8, 3.18_r8, 3.4_r8, 3.25_r8, &
             3.1_r8, 2.95_r8, 2.8_r8, 2.7_r8, 2.6_r8, 2.5_r8, &
             2.4_r8, 2.4_r8, 2.4_r8, 2.4_r8, 2.4_r8, 2.4_r8, &
             2.4_r8, 2.4_r8, 2.4_r8/

call find_pressure_level_weight(pres, k0, wght)
drop_wind_error = wght * obs_err(k0) + (1.0_r8 - wght) * obs_err(k0+1) 
drop_wind_error = dble(nint(drop_wind_error * 10.0_r8)) * 0.1_r8

return
end function drop_wind_error


function prof_wind_error(pres)

real(r8), intent(in) :: pres  !  (mb)

integer  :: k0
real(r8) :: obs_err(nobs_level), wght, prof_wind_error

data obs_err/ 3.7902_r8, 3.7997_r8, 3.8008_r8, 3.8006_r8, 3.7995_r8, &
              3.7924_r8, 3.7654_r8, 3.6912_r8, 3.5278_r8, 3.2342_r8, &
              2.8106_r8, 2.3461_r8, 1.9945_r8, 1.7914_r8, 1.6938_r8, &
              1.6377_r8, 1.5902_r8, 1.5445_r8, 1.5150_r8, 1.5225_r8, &
              1.5700_r8, 1.6504_r8, 1.7469_r8, 1.8341_r8, 1.8947_r8, &
              1.9288_r8, 1.9390_r8, 1.9362_r8, 1.9392_r8, 1.9652_r8, &
              2.0075_r8, 2.0468_r8, 2.0694_r8 /

call find_pressure_level_weight(pres, k0, wght)
prof_wind_error = wght * obs_err(k0) + (1.0_r8 - wght) * obs_err(k0+1)
prof_wind_error = dble(nint(prof_wind_error * 10000.0_r8)) * 0.0001_r8

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
