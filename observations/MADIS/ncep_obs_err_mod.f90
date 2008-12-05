module ncep_obs_err_mod

use        types_mod, only : r8

implicit none
private

real(r8), parameter :: ncep_marine_wind_error    = 2.5_r8, &
                       ncep_marine_temp_error    = 2.5_r8, &
                       ncep_marine_moist_error   = 0.2_r8, &
                       ncep_marine_altim_error   = 1.6_r8

real(r8), parameter :: ncep_land_wind_error      = 3.5_r8, &
                       ncep_land_temp_error      = 2.5_r8, &
                       ncep_land_moist_error     = 0.2_r8, &
                       ncep_land_altim_error     = 1.0_r8

real(r8), parameter :: ncep_acars_wind_error     = 3.0_r8,  &
                       ncep_acars_moist_error    = 0.20_r8

real(r8), parameter :: ncep_rawin_sfcpres_error  = 1.0_r8,  &
                       ncep_rawin_moist_error    = 0.20_r8

real(r8), parameter :: ncep_tc_pos_error         = 0.1_r8,  &
                       ncep_tc_pmin_error        = 2.0_r8,  &
                       ncep_tc_wmax_error        = 3.0_r8

public :: ncep_marine_wind_error,   & 
          ncep_marine_temp_error,   &
          ncep_marine_moist_error,  &
          ncep_marine_altim_error

public :: ncep_land_wind_error,     &
          ncep_land_temp_error,     &
          ncep_land_moist_error,    &
          ncep_land_altim_error

public :: ncep_acars_wind_error,    &
          ncep_acars_temp_error,    &
          ncep_acars_moist_error

public :: ncep_cloud_wind_error,    &
          ncep_wv_wind_error

public :: ncep_rawin_temp_error,    &
          ncep_rawin_wind_error,    &
          ncep_rawin_sfcpres_error, &
          ncep_rawin_moist_error

public :: ncep_tc_pos_error,        &
          ncep_tc_pmin_error,       &
          ncep_tc_wmax_error

contains


function ncep_acars_temp_error(pres)

real(r8), intent(in) :: pres

real(r8) :: ncep_acars_temp_error

if     ( pres <= 800.0_r8  ) then
  ncep_acars_temp_error = 1.00_r8
elseif ( pres <= 850.0_r8  ) then
  ncep_acars_temp_error = 1.12_r8
elseif ( pres <= 900.0_r8  ) then
  ncep_acars_temp_error = 1.24_r8
elseif ( pres <= 950.0_r8  ) then
  ncep_acars_temp_error = 1.35_r8
elseif ( pres <= 1100.0_r8 ) then
  ncep_acars_temp_error = 1.47_r8
else
  print*,'bad pressure level',pres
  stop
end if

return
end function ncep_acars_temp_error


function ncep_cloud_wind_error(pres)

real(r8), intent(in) :: pres  !  (mb)

real(r8) :: ncep_cloud_wind_error

if ( pres <= 250.0_r8 ) then
  ncep_cloud_wind_error = 5.0_r8
  return
else if ( pres <= 300.0_r8 ) then
  ncep_cloud_wind_error = 4.6_r8
  return
else if ( pres <= 350.0_r8 ) then
  ncep_cloud_wind_error = 4.3_r8
  return
else if ( pres <= 400.0_r8 ) then
  ncep_cloud_wind_error = 4.0_r8
  return
else if ( pres <= 450.0_r8 ) then
  ncep_cloud_wind_error = 3.0_r8
  return
else if ( pres <= 500.0_r8 ) then
  ncep_cloud_wind_error = 2.1_r8
  return
else if ( pres <= 600.0_r8 ) then
  ncep_cloud_wind_error = 2.0_r8
  return
else if ( pres <= 700.0_r8 ) then
  ncep_cloud_wind_error = 1.9_r8
  return
else if ( pres <= 1100.0_r8 ) then
  ncep_cloud_wind_error = 1.8_r8
  return
else
  print*,'bad pressure level',pres
  stop
end if

return
end function ncep_cloud_wind_error


function ncep_wv_wind_error(pres)

real(r8), intent(in) :: pres  !  (mb)

real(r8) :: ncep_wv_wind_error

if ( pres <= 250.0_r8 ) then
  ncep_wv_wind_error = 7.0_r8
  return
else if ( pres <= 300.0_r8 ) then
  ncep_wv_wind_error = 6.6_r8
  return
else if ( pres <= 350.0_r8 ) then
  ncep_wv_wind_error = 6.3_r8
  return
else if ( pres <= 400.0_r8 ) then
  ncep_wv_wind_error = 6.0_r8
  return
else if ( pres <= 450.0_r8 ) then
  ncep_wv_wind_error = 5.0_r8
  return
else if ( pres <= 500.0_r8 ) then
  ncep_wv_wind_error = 4.1_r8
  return
else if ( pres <= 600.0_r8 ) then
  ncep_wv_wind_error = 4.0_r8
  return
else if ( pres <= 700.0_r8 ) then
  ncep_wv_wind_error = 3.9_r8
  return
else if ( pres <= 1100.0_r8 ) then
  ncep_wv_wind_error = 3.8_r8
  return
else
  print*,'bad pressure level',pres
  stop
end if

return
end function ncep_wv_wind_error


function ncep_rawin_temp_error(pres)

real(r8), intent(in) :: pres  !  (mb)

real(r8) :: ncep_rawin_temp_error

if ( pres <= 10.0_r8 ) then
  ncep_rawin_temp_error = 1.50_r8
  return
else if ( pres <= 20.0_r8 ) then
  ncep_rawin_temp_error = 1.25_r8
  return
else if ( pres <= 30.0_r8 ) then
  ncep_rawin_temp_error = 1.00_r8
  return
else if ( pres <= 40.0_r8 ) then
  ncep_rawin_temp_error = 0.95_r8
  return
else if ( pres <= 50.0_r8 ) then
  ncep_rawin_temp_error = 0.90_r8
  return
else if ( pres <= 100.0_r8 ) then
  ncep_rawin_temp_error = 0.80_r8
  return
else if ( pres <= 150.0_r8 ) then
  ncep_rawin_temp_error = 1.00_r8
  return
else if ( pres <= 250.0_r8 ) then
  ncep_rawin_temp_error = 1.20_r8
  return
else if ( pres <= 300.0_r8 ) then
  ncep_rawin_temp_error = 0.90_r8
  return
else if ( pres <= 850.0_r8 ) then
  ncep_rawin_temp_error = 0.80_r8
  return
else if ( pres <= 900.0_r8 ) then
  ncep_rawin_temp_error = 0.90_r8
  return
else if ( pres <= 950.0_r8 ) then
  ncep_rawin_temp_error = 1.10_r8
  return
else if ( pres <= 1100.0_r8 ) then
  ncep_rawin_temp_error = 1.20_r8
  return
else
  print*,'bad pressure level',pres
  stop
end if

return
end function ncep_rawin_temp_error


function ncep_rawin_wind_error(pres)

real(r8), intent(in) :: pres  !  (mb)

real(r8) :: ncep_rawin_wind_error

if ( pres <= 100.0_r8 ) then
  ncep_rawin_wind_error = 2.1_r8
  return
else if ( pres <= 150.0_r8 ) then
  ncep_rawin_wind_error = 2.4_r8
  return
else if ( pres <= 200.0_r8 ) then
  ncep_rawin_wind_error = 2.7_r8
  return
else if ( pres <= 250.0_r8 ) then
  ncep_rawin_wind_error = 3.2_r8
  return
else if ( pres <= 300.0_r8 ) then
  ncep_rawin_wind_error = 3.0_r8
  return
else if ( pres <= 350.0_r8 ) then
  ncep_rawin_wind_error = 2.8_r8
  return
else if ( pres <= 400.0_r8 ) then
  ncep_rawin_wind_error = 2.6_r8
  return
else if ( pres <= 450.0_r8 ) then
  ncep_rawin_wind_error = 2.3_r8
  return
else if ( pres <= 500.0_r8 ) then
  ncep_rawin_wind_error = 2.1_r8
  return
else if ( pres <= 550.0_r8 ) then
  ncep_rawin_wind_error = 2.0_r8
  return
else if ( pres <= 600.0_r8 ) then
  ncep_rawin_wind_error = 1.9_r8
  return
else if ( pres <= 650.0_r8 ) then
  ncep_rawin_wind_error = 1.8_r8
  return
else if ( pres <= 800.0_r8 ) then
  ncep_rawin_wind_error = 1.6_r8
  return
else if ( pres <= 950.0_r8 ) then
  ncep_rawin_wind_error = 1.5_r8
  return
else if ( pres <= 1100.0_r8 ) then
  ncep_rawin_wind_error = 1.4_r8
  return
else
  print*,'bad pressure level',pres
  stop
end if

return
end function ncep_rawin_wind_error


end module ncep_obs_err_mod
