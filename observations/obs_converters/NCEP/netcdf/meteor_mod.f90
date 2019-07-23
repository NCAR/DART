! DART software - Copyright UCAR. This open source software is provided
! by UCAR, "as is", without charge, subject to all terms of use at
! http://www.image.ucar.edu/DAReS/DART/DART_download
!
! $Id$

module meteor_mod

use        types_mod, only : r8, deg2rad

implicit none
private

public :: invert_altimeter,        &
          pres_alt_to_pres,        & 
          sat_vapor_pressure,      & 
          specific_humidity,       &
          theta_to_temp,           & 
          wind_dirspd_to_uv,       &
          rh_and_temp_to_dewpoint, &
          temp_and_dewpoint_to_rh

real(r8), parameter :: grav    = 9.81_r8,      & ! gravitational constant
                       Cp      = 1004.5_r8,    &
                       Rd      = 287.0_r8,     &
                       Rv      = 461.6_r8,     &
                       Lvap    = 2500000.0_r8, &
                       R_earth = 6370.0_r8,    &
                       Pref    = 100000.0,     & ! reference pressure
                       Pralt   = 101325.0,     & ! altimeter reference P
                       Talt    = 288.15,       & ! altimeter reference T
                       es0C    = 611.0_r8,     & ! vapor pressure at 0 C (Pa)
                       Tfrez   = 273.15_r8,    & ! water freezing point (K)
                       RvRd    = Rv/Rd,        & ! Rv/Rd
                       RdRv    = Rd/Rv,        & ! Rd/Rv  added 11/2009
                       kappa   = Rd/Cp,        & ! kappa for pot. temp
                       dTdzsta =0.0065           ! standard atmosphere lapse rate K/m

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
! Clausius-Clapeyron w/ constant Lv in Pa
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

! 11/2009 changed to Emanuel (1994) eqn 4.1.4
 specific_humidity = (RdRv * vapor_pres) / (pres - vapor_pres*(1.0_r8-RdRv))

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

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!   sat_vapor_press_bolton - function that uses Bolton's approximation to
!                            compute saturation vapor pressure given
!                            temperature.
!
!   reference:  Bolton 1980, MWR, 1046-1053
!
!    sat_vapor_press_bolton - saturation vapor pressure (Pa)
!    tmpk                   - temperature (K)
!
!     created Dec. 2008 David Dowell, NCAR/MMM
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
function sat_vapor_press_bolton(tmpk)

real(r8)             :: sat_vapor_press_bolton
real(r8), intent(in) :: tmpk

real(r8)             :: tmpc               ! temperature (Celsius)

tmpc = tmpk - Tfrez
if ( tmpc <= -200.0_r8 ) then
  print*,'sat_vapor_press_bolton:  tmpc too low ',tmpc
  stop
end if
sat_vapor_press_bolton = es0C * exp( 17.67_r8 * tmpc / (tmpc + 243.5_r8) )

return
end function sat_vapor_press_bolton

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!   temp_and_dewpoint_to_rh - function that computes the relative humidity
!                             given temperature and dewpoint
!
!    temp_and_dewpoint_to_rh - relative humidity (0.00 - 1.00)
!    tmpk                    - temperature (Kelvin)
!    dptk                    - dewpoint (Kelvin)
!
!     created Dec. 2008 David Dowell, NCAR/MMM
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
function temp_and_dewpoint_to_rh(tmpk, dptk)

real(r8)             :: temp_and_dewpoint_to_rh
real(r8), intent(in) :: tmpk
real(r8), intent(in) :: dptk

real(r8)             :: e                  ! vapor pressure (Pa)
real(r8)             :: es                 ! saturation vapor pressure (Pa)

e = sat_vapor_press_bolton(dptk)
es = sat_vapor_press_bolton(tmpk)

temp_and_dewpoint_to_rh = e / es

if (temp_and_dewpoint_to_rh > 1.00_r8) then
  print*,'rh = ', temp_and_dewpoint_to_rh, ', resetting to 1.00'
  temp_and_dewpoint_to_rh = 1.00_r8
end if

return
end function temp_and_dewpoint_to_rh

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!   rh_and_temp_to_dewpoint - function that computes the dewpoint
!                             given relative humidity and temperature
!
!    rh_and_temp_to_dewpoint - dewpoint (Kelvin)
!    rh                      - relative humidity (0.00 - 1.00)
!    tmpk                    - temperature (Kelvin)
!
!     created Dec. 2008 David Dowell, NCAR/MMM
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
function rh_and_temp_to_dewpoint(rh, tmpk)

real(r8)             :: rh_and_temp_to_dewpoint
real(r8), intent(in) :: rh
real(r8), intent(in) :: tmpk

real(r8)             :: e                  ! vapor pressure (Pa)
real(r8)             :: es                 ! saturation vapor pressure (Pa)
real(r8)             :: dptc               ! dptc (Celsius)

if ( ( rh <= 0.00_r8 ) .or. ( rh > 1.00_r8 ) ) then
  print*,'rh_and_temp_to_dewpoint:  bad rh ',rh
  stop
end if
if ( rh <= 0.01_r8 ) then
  print*,'rh_and_temp_to_dewpoint: low rh ',rh
end if

es = sat_vapor_press_bolton(tmpk)
e = rh * es

dptc = 243.5_r8 / (17.67_r8 / log(e/es0C) - 1.0_r8)

rh_and_temp_to_dewpoint = dptc + Tfrez

return
end function rh_and_temp_to_dewpoint


end module meteor_mod

! <next few lines under version control, do not edit>
! $URL$
! $Id$
! $Revision$
! $Date$
