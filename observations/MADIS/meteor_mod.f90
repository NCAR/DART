module meteor_mod

use        types_mod, only : r8, deg2rad

implicit none
private

public :: invert_altimeter,   &
          pres_alt_to_pres,   & 
          sat_vapor_pressure, & 
          specific_humidity,  &
          theta_to_temp,      & 
          wind_dirspd_to_uv

real(r8), parameter :: grav    = 9.81_r8,      & ! gravitational constant
                       Cp      = 1004.5_r8,    &
                       Rd      = 287.0_r8,     &
                       Rv      = 461.6_r8,     &
                       Lvap    = 2500000.0_r8, &
                       R_earth = 6370.0_r8,    &
                       Pref    = 100000.0,     & ! reference pressure
                       Pralt   = 101325.0,     & ! altimeter reference P
                       Talt    = 288.15,       & ! attimeter reference T
                       es0C    = 611.0_r8,     & ! vapor pressure at 0 C (Pa)
                       Tfrez   = 273.15_r8,    & ! water freezing point
                       RvRd    = Rv/Rd,        & ! Rd/Rv
                       kappa   = Rd/Cp,        & ! kappa for pot. temp
                       dTdzsta =0.0065          ! standard atmosphere

contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!   invert_altimeter - function that computes the surface pressure
!                      given an altimeter setting and the station
!                      elevation.
!
!    altimeter_setting - altimeter setting (hPa)
!    elevation         - elevation of station (m)
!    invert_altimeter  - surface pressure value (hPa)
!
!     created Dec. 2007 Ryan Torn, NCAR/MMM
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
function invert_altimeter(altimeter_setting, elevation)

use        types_mod, only : r8

implicit none

real(r8), parameter :: k1 = 0.190284_r8
real(r8), parameter :: k2 = 8.4228807E-5_r8

real(r8), intent(in) :: altimeter_setting, elevation

real(r8) :: invert_altimeter !  (hPa)

invert_altimeter = (altimeter_setting ** k1 - k2 * elevation) ** (1 / k1) + 0.3_r8

return
end function invert_altimeter

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!   pres_alt_to_pres - function that computes the pressure level given
!                      the pressure altitude.  Used mostly for ACARS 
!                      observations.
!
!    hght - pressure-height level (m)
!
!     created Oct. 2007 Ryan Torn, NCAR/MMM
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
function pres_alt_to_pres(hght)

real(r8), parameter  :: Po   = 101325.0_r8
real(r8), intent(in) :: hght

real(r8) :: C1, pres_alt_to_pres

C1 = grav / (dTdzsta * Rd)
pres_alt_to_pres = Po * exp( C1 * log( 1 - (dTdzsta * hght) / Talt) )

end function pres_alt_to_pres

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!   sat_vapor_pressure - function that computes the water vapor 
!                        saturation vapor pressure given a temperature.
!
!    tmpk - temperature (K)
!
!     created Oct. 2007 Ryan Torn, NCAR/MMM
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
function sat_vapor_pressure(tmpk)

real(r8), intent(in) :: tmpk

real(r8) :: sat_vapor_pressure

sat_vapor_pressure = es0C * exp((Lvap/Rv)*(1.0_r8/Tfrez - 1.0_r8/tmpk))

return
end function sat_vapor_pressure

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!   specific_humidity - function that computes the specific humidity 
!                       given the vapor pressure and atmospheric dry 
!                       air pressure.
!
!    vapor_pres - vapor pressure (Pa)
!    pres       - atmospheric pressure (Pa)
!
!     created Dec. 2007 Ryan Torn, NCAR/MMM 
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
function specific_humidity(vapor_pres, pres)

real(r8), intent(in) :: vapor_pres, pres

real(r8) :: specific_humidity

specific_humidity = (vapor_pres / (pres-vapor_pres)) / RvRd

return
end function specific_humidity

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!   theta_to_temp - function that computes the temperature given a 
!                   potential temperature and the pressure level.
!
!    thta - potential temperature (K)
!    pres - pressure (Pa)
!
!     created Oct. 2007 Ryan Torn, NCAR/MMM
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
function theta_to_temp(thta, pres)

real(r8), intent(in) :: thta, pres

real(r8) :: theta_to_temp

theta_to_temp = thta * (Pref/pres) ** (-kappa) 

return
end function theta_to_temp

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!   wind_dirspd_to_uv - subroutine that converts a wind direction and 
!                       wind speed to a zonal and meridional wind 
!                       component.
!
!    wdir - wind direction
!    wspd - wind speed (m/s)
!    uwnd - u component of the wind (m/s)
!    vwnd - v component of the wind (m/s)
!
!     created Oct. 2007 Ryan Torn, NCAR/MMM
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine wind_dirspd_to_uv(wdir, wspd, uwnd, vwnd)

real(r8), intent(in)  :: wdir, wspd
real(r8), intent(out) :: uwnd, vwnd

uwnd =  wspd * cos(deg2rad * (90.0_r8 + wdir))
vwnd = -wspd * sin(deg2rad * (90.0_r8 + wdir))

return
end subroutine wind_dirspd_to_uv

end module meteor_mod
