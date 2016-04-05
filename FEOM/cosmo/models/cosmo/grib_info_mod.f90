MODULE grib_info_mod

! <next few lines under version control, do not edit>
! $URL$
! $Id$
! $Revision$
! $Date$

!----------------------------------------------------------------------
! purpose: test routines
!----------------------------------------------------------------------

  USE byte_mod, ONLY: concat_bytes1,to_positive
  USE obs_kind_mod, only : KIND_CLOUD_FRACTION, &
                           KIND_CLOUD_ICE, &
                           KIND_CLOUD_LIQUID_WATER, &
                           KIND_PRESSURE, &
                           KIND_PRESSURE_PERTURBATION, &
                           KIND_RELATIVE_HUMIDITY, &
                           KIND_SEA_SURFACE_PRESSURE, &
                           KIND_SOIL_MOISTURE, &
                           KIND_SPECIFIC_HUMIDITY, &
                           KIND_SURFACE_ELEVATION, &
                           KIND_SURFACE_PRESSURE, &
                           KIND_TEMPERATURE, &
                           KIND_U_WIND_COMPONENT, &
                           KIND_V_WIND_COMPONENT, &
                           KIND_VERTICAL_VELOCITY,&
                           KIND_SURFACE_GEOPOTENTIAL

  implicit none

  ! version controlled file description for error handling, do not edit
character(len=128), parameter :: &
   source   = "$URL$", &
   revision = "$Revision$", &
   revdate  = "$Date$"

  TYPE varnames
    CHARACTER(len=256) :: longname
    CHARACTER(len=16)  :: shortname
    CHARACTER(len=32)  :: unit
  END TYPE varnames

  TYPE(varnames) :: gribtab(1:255,1:4,1:22)
  INTEGER        :: dartkinds(1:255,1:4,1:22)

  private :: set_gribtab
  private :: set_dartkinds
  public :: get_varname
  public :: get_dart_kind
  public :: get_level
  public :: get_dims

CONTAINS
  
  SUBROUTINE set_gribtab()  
!    gribtab(2,001)%longname="Pressure"
!    gribtab(2,001)%shortname="PP"
!    gribtab(2,001)%unit="Pascal"
!    gribtab(2,002)%longname="Pressure reduced to MSL"
!    gribtab(2,002)%shortname="PMSL"
!    gribtab(2,002)%unit="Pascal"
!    gribtab(2,003)%longname="Pressure tendency"
!    gribtab(2,003)%shortname="PPT"
!    gribtab(2,003)%unit="Pascal/second"
!    gribtab(2,006)%longname="Geopotential"
!    gribtab(2,006)%shortname="FI"
!    gribtab(2,006)%unit="square meters/square second"
!    gribtab(2,007)%longname="Geopotential height"
!    gribtab(2,007)%shortname="GPH"
!    gribtab(2,007)%unit="Geopotential meters"
!    gribtab(2,008)%longname="Geometric height"
!    gribtab(2,008)%shortname="HH"
!    gribtab(2,008)%unit="meters"
!    gribtab(2,009)%longname="Standard deviation of height"
!    gribtab(2,009)%shortname="HHsd"
!    gribtab(2,009)%unit="meters"
!    gribtab(2,011)%longname="Temperature"
!    gribtab(2,011)%shortname="T"
!    gribtab(2,011)%unit="Kelvin"
!    gribtab(2,012)%longname="Virtual temperature"
!    gribtab(2,012)%shortname="TV"
!    gribtab(2,012)%unit="Kelvin"
!    gribtab(2,013)%longname="Potential temperature"
!    gribtab(2,013)%shortname="THETA"
!    gribtab(2,013)%unit="Kelvin"
!    gribtab(2,014)%longname="Pseudo-adiabatic potential temperature"
!    gribtab(2,014)%shortname="Tpa"
!    gribtab(2,014)%unit="Kelvin"
!    gribtab(2,015)%longname="Maximum temperature"
!    gribtab(2,015)%shortname="Tmax"
!    gribtab(2,015)%unit="Kelvin"
!    gribtab(2,016)%longname="Minimum temperature"
!    gribtab(2,016)%shortname="Tmin"
!    gribtab(2,016)%unit="Kelvin"
!    gribtab(2,017)%longname="Dew point temperature"
!    gribtab(2,017)%shortname="TD"
!    gribtab(2,017)%unit="Kelvin"
!    gribtab(2,018)%longname="Dew point depression (or deficit)"
!    gribtab(2,018)%shortname="TDDP"
!    gribtab(2,018)%unit="Kelvin"
!    gribtab(2,019)%longname="Lapse rate"
!    gribtab(2,019)%shortname="LR"
!    gribtab(2,019)%unit="Kelvin/meter"
!    gribtab(2,020)%longname="Visibility"
!    gribtab(2,020)%shortname="VIS"
!    gribtab(2,020)%unit="meters"
!    gribtab(2,021)%longname="Radar Spectra (1)"
!    gribtab(2,021)%shortname="RSp 1"
!    gribtab(2,021)%unit="-"
!    gribtab(2,022)%longname="Radar Spectra (2)"
!    gribtab(2,022)%shortname="RSp 2"
!    gribtab(2,022)%unit="-"
!    gribtab(2,023)%longname="Radar Spectra (3)"
!    gribtab(2,023)%shortname="RSp 3"
!    gribtab(2,023)%unit="-"
!    gribtab(2,025)%longname="Temperature anomaly"
!    gribtab(2,025)%shortname="Tanom"
!    gribtab(2,025)%unit="Kelvin"
!    gribtab(2,026)%longname="Pressure anomaly"
!    gribtab(2,026)%shortname="Panom"
!    gribtab(2,026)%unit="Pascal"
!    gribtab(2,027)%longname="Geopotential height anomaly"
!    gribtab(2,027)%shortname="GPHanom"
!    gribtab(2,027)%unit="Geopotential meters"
!    gribtab(2,028)%longname="Wave Spectra (1)"
!    gribtab(2,028)%shortname="WSp 1"
!    gribtab(2,028)%unit="-"
!    gribtab(2,029)%longname="Wave Spectra (2)"
!    gribtab(2,029)%shortname="WSp 2"
!    gribtab(2,029)%unit="-"
!    gribtab(2,030)%longname="Wave Spectra (3)"
!    gribtab(2,030)%shortname="WSp 3"
!    gribtab(2,030)%unit="-"
!    gribtab(2,031)%longname="Wind direction"
!    gribtab(2,031)%shortname="DD"
!    gribtab(2,031)%unit="degree"
!    gribtab(2,032)%longname="Wind speed"
!    gribtab(2,032)%shortname="FF"
!    gribtab(2,032)%unit="meters/second"
!    gribtab(2,033)%longname="u-component of wind"
!    gribtab(2,033)%shortname="U"
!    gribtab(2,033)%unit="meters/second"
!    gribtab(2,034)%longname="v-component of wind"
!    gribtab(2,034)%shortname="V"
!    gribtab(2,034)%unit="meters/second"
!    gribtab(2,035)%longname="Stream function"
!    gribtab(2,035)%shortname="PSI"
!    gribtab(2,035)%unit="square meters/second"
!    gribtab(2,036)%longname="Velocity potential"
!    gribtab(2,036)%shortname="VPOT"
!    gribtab(2,036)%unit="square meters/second"
!    gribtab(2,037)%longname="Montgomery stream function"
!    gribtab(2,037)%shortname="MontPSI"
!    gribtab(2,037)%unit="square meters/square second"
!    gribtab(2,038)%longname="Sigma coord. vertical velocity"
!    gribtab(2,038)%shortname="WSIGMA"
!    gribtab(2,038)%unit="seconds/second"
!    gribtab(2,039)%longname="Pressure Vertical velocity"
!    gribtab(2,039)%shortname="OMEGA"
!    gribtab(2,039)%unit="Pascal/second"
!    gribtab(2,040)%longname="Geometric Vertical velocity"
!    gribtab(2,040)%shortname="W"
!    gribtab(2,040)%unit="meters/second"
!    gribtab(2,041)%longname="Absolute vorticity"
!    gribtab(2,041)%shortname="VOR"
!    gribtab(2,041)%unit="/second"
!    gribtab(2,042)%longname="Absolute divergence"
!    gribtab(2,042)%shortname="DIV"
!    gribtab(2,042)%unit="/second"
!    gribtab(2,043)%longname="Relative vorticity"
!    gribtab(2,043)%shortname="VORrel"
!    gribtab(2,043)%unit="/second"
!    gribtab(2,044)%longname="Relative divergence"
!    gribtab(2,044)%shortname="DIVrel"
!    gribtab(2,044)%unit="/second"
!    gribtab(2,045)%longname="Vertical u-component shear"
!    gribtab(2,045)%shortname="U shear"
!    gribtab(2,045)%unit="/second"
!    gribtab(2,046)%longname="Vertical v-component shear"
!    gribtab(2,046)%shortname="V shear"
!    gribtab(2,046)%unit="/second"
!    gribtab(2,047)%longname="Direction of current"
!    gribtab(2,047)%shortname="DirCurr"
!    gribtab(2,047)%unit="degree"
!    gribtab(2,048)%longname="Speed of current"
!    gribtab(2,048)%shortname="VCurr"
!    gribtab(2,048)%unit="meters/second"
!    gribtab(2,049)%longname="u-component of current"
!    gribtab(2,049)%shortname="Ucurr"
!    gribtab(2,049)%unit="meters/second"
!    gribtab(2,050)%longname="v-component of current"
!    gribtab(2,050)%shortname="Vcurr"
!    gribtab(2,050)%unit="meters/second"
!    gribtab(2,051)%longname="Specific humidity"
!    gribtab(2,051)%shortname="QV"
!    gribtab(2,051)%unit="kilogram/kilogram"
!    gribtab(2,052)%longname="Relative humidity"
!    gribtab(2,052)%shortname="RH"
!    gribtab(2,052)%unit="percent"
!    gribtab(2,053)%longname="Humidity mixing ratio"
!    gribtab(2,053)%shortname="MR"
!    gribtab(2,053)%unit="kilogram/kilogram"
!    gribtab(2,054)%longname="Precipitable water"
!    gribtab(2,054)%shortname="PrWater"
!    gribtab(2,054)%unit="kilogram/square meter"
!    gribtab(2,055)%longname="Vapor pressure"
!    gribtab(2,055)%shortname="VPres"
!    gribtab(2,055)%unit="Pascal"
!    gribtab(2,056)%longname="Saturation deficit"
!    gribtab(2,056)%shortname="SatDef"
!    gribtab(2,056)%unit="Pascal"
!    gribtab(2,057)%longname="Evaporation"
!    gribtab(2,057)%shortname="EVA"
!    gribtab(2,057)%unit="kilogram/square meter"
!    gribtab(2,058)%longname="Cloud Ice"
!    gribtab(2,058)%shortname="CI"
!    gribtab(2,058)%unit="kilogram/square meter"
!    gribtab(2,059)%longname="Precipitation rate"
!    gribtab(2,059)%shortname="PR"
!    gribtab(2,059)%unit="kilogram/square meter/s"
!    gribtab(2,060)%longname="Thunderstorm probability"
!    gribtab(2,060)%shortname="ThProb"
!    gribtab(2,060)%unit="percent"
!    gribtab(2,061)%longname="Total precipitation"
!    gribtab(2,061)%shortname="TOTPrec"
!    gribtab(2,061)%unit="kilogram/square meter"
!    gribtab(2,062)%longname="Large scale precipitation"
!    gribtab(2,062)%shortname="LSPrec"
!    gribtab(2,062)%unit="kilogram/square meter"
!    gribtab(2,063)%longname="Convective precipitation"
!    gribtab(2,063)%shortname="CONVPrec"
!    gribtab(2,063)%unit="kilogram/square meter"
!    gribtab(2,064)%longname="Snowfall rate water equivalent"
!    gribtab(2,064)%shortname="SNOWrate"
!    gribtab(2,064)%unit="kilogram/square meter * second"
!    gribtab(2,065)%longname="Water equiv. of accum. snow depth"
!    gribtab(2,065)%shortname="ACCSNOWdep"
!    gribtab(2,065)%unit="kilogram/square meter"
!    gribtab(2,066)%longname="Snow depth"
!    gribtab(2,066)%shortname="SNOWdep"
!    gribtab(2,066)%unit="meters"
!    gribtab(2,067)%longname="Mixed layer depth"
!    gribtab(2,067)%shortname="MLdep"
!    gribtab(2,067)%unit="meters"
!    gribtab(2,068)%longname="Transient thermocline depth"
!    gribtab(2,068)%shortname="transTCLdep"
!    gribtab(2,068)%unit="meters"
!    gribtab(2,069)%longname="Main thermocline depth"
!    gribtab(2,069)%shortname="TCLdep"
!    gribtab(2,069)%unit="meters"
!    gribtab(2,070)%longname="Main thermocline anomaly"
!    gribtab(2,070)%shortname="TCLanom"
!    gribtab(2,070)%unit="meters"
!    gribtab(2,071)%longname="Total cloud cover"
!    gribtab(2,071)%shortname="CLCT"
!    gribtab(2,071)%unit="percent"
!    gribtab(2,072)%longname="Convective cloud cover"
!    gribtab(2,072)%shortname="CLCC"
!    gribtab(2,072)%unit="percent"
!    gribtab(2,073)%longname="Low cloud cover"
!    gribtab(2,073)%shortname="CLCL"
!    gribtab(2,073)%unit="percent"
!    gribtab(2,074)%longname="Medium cloud cover"
!    gribtab(2,074)%shortname="CLCM"
!    gribtab(2,074)%unit="percent"
!    gribtab(2,075)%longname="High cloud cover"
!    gribtab(2,075)%shortname="CLCH"
!    gribtab(2,075)%unit="percent"
!    gribtab(2,076)%longname="Cloud water"
!    gribtab(2,076)%shortname="QC"
!    gribtab(2,076)%unit="kilogram/square meter"
!    gribtab(2,078)%longname="Convective snow"
!    gribtab(2,078)%shortname="ConvSNOW"
!    gribtab(2,078)%unit="kilogram/square meter"
!    gribtab(2,079)%longname="Large scale snow"
!    gribtab(2,079)%shortname="LSSNOW"
!    gribtab(2,079)%unit="kilogram/square meter"
!    gribtab(2,080)%longname="Water Temperature"
!    gribtab(2,080)%shortname="TWater"
!    gribtab(2,080)%unit="Kelvin"
!    gribtab(2,081)%longname="Land-sea mask"
!    gribtab(2,081)%shortname="LandSea"
!    gribtab(2,081)%unit="Fraction"
!    gribtab(2,082)%longname="Deviation of sea level from mean"
!    gribtab(2,082)%shortname="SeaLevDev"
!    gribtab(2,082)%unit="meters"
!    gribtab(2,083)%longname="Surface roughness"
!    gribtab(2,083)%shortname="Z0"
!    gribtab(2,083)%unit="meters"
!    gribtab(2,084)%longname="Albedo"
!    gribtab(2,084)%shortname="Albedo"
!    gribtab(2,084)%unit="percent"
!    gribtab(2,085)%longname="Soil temperature"
!    gribtab(2,085)%shortname="Tsoil"
!    gribtab(2,085)%unit="Kelvin"
!    gribtab(2,086)%longname="Soil moisture content"
!    gribtab(2,086)%shortname="SM"
!    gribtab(2,086)%unit="kilogram/square meter"
!    gribtab(2,087)%longname="Vegetation"
!    gribtab(2,087)%shortname="VEG"
!    gribtab(2,087)%unit="percent"
!    gribtab(2,088)%longname="Salinity"
!    gribtab(2,088)%shortname="SAL"
!    gribtab(2,088)%unit="kilogram/kilogram"
!    gribtab(2,089)%longname="Density"
!    gribtab(2,089)%shortname="RHO"
!    gribtab(2,089)%unit="kilogram/m3"
!    gribtab(2,090)%longname="Water run off"
!    gribtab(2,090)%shortname="RUNOFF"
!    gribtab(2,090)%unit="kilogram/square meter"
!    gribtab(2,091)%longname="Ice concentration"
!    gribtab(2,091)%shortname="ICEconc"
!    gribtab(2,091)%unit="Fraction"
!    gribtab(2,092)%longname="Ice thickness"
!    gribtab(2,092)%shortname="ICEthick"
!    gribtab(2,092)%unit="meters"
!    gribtab(2,093)%longname="Direction of ice drift"
!    gribtab(2,093)%shortname="DirIceDrift"
!    gribtab(2,093)%unit="degree"
!    gribtab(2,094)%longname="Speed of ice drift"
!    gribtab(2,094)%shortname="VIceDrift"
!    gribtab(2,094)%unit="meters/second"
!    gribtab(2,095)%longname="u-component of ice drift"
!    gribtab(2,095)%shortname="UIceDrift"
!    gribtab(2,095)%unit="meters/second"
!    gribtab(2,096)%longname="v-component of ice drift"
!    gribtab(2,096)%shortname="VIceDrift"
!    gribtab(2,096)%unit="meters/second"
!    gribtab(2,097)%longname="Ice growth rate"
!    gribtab(2,097)%shortname="IceGR"
!    gribtab(2,097)%unit="meters/second"
!    gribtab(2,098)%longname="Ice divergence"
!    gribtab(2,098)%shortname="IceDiv"
!    gribtab(2,098)%unit="/second"
!    gribtab(2,099)%longname="Snow melt"
!    gribtab(2,099)%shortname="SNOWmelt"
!    gribtab(2,099)%unit="kilogram/square meter"
!    gribtab(2,100)%longname="Significant height of combined wind"
!    gribtab(2,100)%shortname="HcombWind"
!    gribtab(2,100)%unit="meters"
!    gribtab(2,101)%longname="Direction of wind waves"
!    gribtab(2,101)%shortname="DirWave"
!    gribtab(2,101)%unit="degree"
!    gribtab(2,102)%longname="Significant height of wind waves"
!    gribtab(2,102)%shortname="HWave"
!    gribtab(2,102)%unit="meters"
!    gribtab(2,103)%longname="Mean period of wind waves"
!    gribtab(2,103)%shortname="PerWave"
!    gribtab(2,103)%unit="seconds"
!    gribtab(2,104)%longname="Direction of swell waves"
!    gribtab(2,104)%shortname="DirSwell"
!    gribtab(2,104)%unit="degree"
!    gribtab(2,105)%longname="Significant height of swell waves"
!    gribtab(2,105)%shortname="HSwell"
!    gribtab(2,105)%unit="meters"
!    gribtab(2,106)%longname="Mean period of swell waves"
!    gribtab(2,106)%shortname="PerSwell"
!    gribtab(2,106)%unit="seconds"
!    gribtab(2,107)%longname="Primary wave direction"
!    gribtab(2,107)%shortname="DirPrimWave"
!    gribtab(2,107)%unit="degree"
!    gribtab(2,108)%longname="Primary wave mean period"
!    gribtab(2,108)%shortname="PerPrimWave"
!    gribtab(2,108)%unit="seconds"
!    gribtab(2,109)%longname="Secondary wave direction"
!    gribtab(2,109)%shortname="DirSecWave"
!    gribtab(2,109)%unit="degree"
!    gribtab(2,110)%longname="Secondary wave mean period"
!    gribtab(2,110)%shortname="PerSecWave"
!    gribtab(2,110)%unit="seconds"
!    gribtab(2,111)%longname="Net short-wave radiation (surface)"
!    gribtab(2,111)%shortname="netSWRadSurf"
!    gribtab(2,111)%unit="W/square meter"
!    gribtab(2,112)%longname="Net long wave radiation (surface)"
!    gribtab(2,112)%shortname="netLWRadSurf"
!    gribtab(2,112)%unit="W/square meter"
!    gribtab(2,113)%longname="Net short-wave radiation (top of atmos)"
!    gribtab(2,113)%shortname="netSWRadTop"
!    gribtab(2,113)%unit="W/square meter"
!    gribtab(2,114)%longname="Net long wave radiation (top of atmos)"
!    gribtab(2,114)%shortname="netLWRadTop"
!    gribtab(2,114)%unit="W/square meter"
!    gribtab(2,115)%longname="Long wave radiation"
!    gribtab(2,115)%shortname="LWRad"
!    gribtab(2,115)%unit="W/square meter"
!    gribtab(2,116)%longname="Short wave radiation"
!    gribtab(2,116)%shortname="SWRad"
!    gribtab(2,116)%unit="W/square meter"
!    gribtab(2,117)%longname="Global radiation"
!    gribtab(2,117)%shortname="RadGlob"
!    gribtab(2,117)%unit="W/square meter"
!    gribtab(2,121)%longname="Latent heat net flux"
!    gribtab(2,121)%shortname="netLH"
!    gribtab(2,121)%unit="W/square meter"
!    gribtab(2,122)%longname="Sensible heat net flux"
!    gribtab(2,122)%shortname="netSH"
!    gribtab(2,122)%unit="W/square meter"
!    gribtab(2,123)%longname="Boundary layer dissipation"
!    gribtab(2,123)%shortname="BLdiss"
!    gribtab(2,123)%unit="W/square meter"
!    gribtab(2,124)%longname="Momentum flux, u component"
!    gribtab(2,124)%shortname="UMomFlux"
!    gribtab(2,124)%unit="N/square meter"
!    gribtab(2,125)%longname="Momentum flux, v component"
!    gribtab(2,125)%shortname="VMomFlux"
!    gribtab(2,125)%unit="N/square meter"
!    gribtab(2,126)%longname="Wind mixing energy"
!    gribtab(2,126)%shortname="WindMixEner"
!    gribtab(2,126)%unit="Joule"

    gribtab(:,:,:)%longname=""
    gribtab(:,:,:)%shortname=""
    gribtab(:,:,:)%unit=""

    gribtab(8,1,19)%shortname="HHL"
    gribtab(8,1,19)%longname="geometric height layer boundaries above sea level"
    gribtab(8,1,19)%unit="m"
    gribtab(6,1,1)%shortname="FIS"
    gribtab(6,1,1)%longname="geopotential at surface"
    gribtab(6,1,1)%unit="m2/s2"
    gribtab(8,1,1)%shortname="HSURF"
    gribtab(8,1,1)%longname="geometric height surface"
    gribtab(8,1,1)%unit="m"
    gribtab(81,1,1)%shortname="FR_LAND"
    gribtab(81,1,1)%longname="land fraction"
    gribtab(81,1,1)%unit="1"
    gribtab(46,3,1)%shortname="SSO_STDH"
    gribtab(46,3,1)%longname="standard deviation of sub-scale orography"
    gribtab(46,3,1)%unit="m"
    gribtab(47,3,1)%shortname="SSO_GAMMA"
    gribtab(47,3,1)%longname="horizontal anisotropy sub-scale orography"
    gribtab(47,3,1)%unit="-"
    gribtab(48,3,1)%shortname="SSO_THETA"
    gribtab(48,3,1)%longname="angle,between main axis of SSO and East"
    gribtab(48,3,1)%unit="rad"
    gribtab(49,3,1)%shortname="SSO_SIGMA"
    gribtab(49,3,1)%longname="mean inclination of sub-scale orography"
    gribtab(49,3,1)%unit="1"
    gribtab(57,3,1)%shortname="SOILTYP"
    gribtab(57,3,1)%longname="soil type"
    gribtab(57,3,1)%unit="-"
    gribtab(114,3,1)%shortname="RLAT"
    gribtab(114,3,1)%longname="latitude"
    gribtab(114,3,1)%unit="degree north"
    gribtab(115,3,1)%shortname="RLON"
    gribtab(115,3,1)%longname="longitude"
    gribtab(115,3,1)%unit="degree east"
    gribtab(113,3,1)%shortname="FC"
    gribtab(113,3,1)%longname="Coriolis parameter"
    gribtab(113,3,1)%unit="1/s"
    gribtab(85,1,21)%shortname="T_CL"
    gribtab(85,1,21)%longname="temperature at the lower boundary of the middle soil layer"
    gribtab(85,1,21)%unit="K"
    gribtab(86,1,22)%shortname="W_CL"
    gribtab(86,1,22)%longname="water content at the lower boundary of the lower soil layer"
    gribtab(86,1,22)%unit="kg/m2"
    gribtab(87,1,1)%shortname="PLCOV"
    gribtab(87,1,1)%longname="vegetation cover"
    gribtab(87,1,1)%unit="%"
    gribtab(75,3,1)%shortname="FOR_E"
    gribtab(75,3,1)%longname="vegetation cover evergreen forest"
    gribtab(75,3,1)%unit="1"
    gribtab(76,3,1)%shortname="FOR_D"
    gribtab(76,3,1)%longname="vegetation cover seasonal forest"
    gribtab(76,3,1)%unit="1"
    gribtab(61,3,1)%shortname="LAI"
    gribtab(61,3,1)%longname="leaf area index"
    gribtab(61,3,1)%unit="1"
    gribtab(62,3,1)%shortname="ROOTDP"
    gribtab(62,3,1)%longname="root length"
    gribtab(62,3,1)%unit="m"
    gribtab(55,3,1)%shortname="FR_LAKE"
    gribtab(55,3,1)%longname="lake fraction"
    gribtab(55,3,1)%unit="1"
    gribtab(96,2,1)%shortname="DEPTH_LK"
    gribtab(96,2,1)%longname="lake depth"
    gribtab(96,2,1)%unit="m"
    gribtab(64,3,1)%shortname="HMO3"
    gribtab(64,3,1)%longname="height of ozone maximum"
    gribtab(64,3,1)%unit="Pa"
    gribtab(65,3,1)%shortname="VIO3"
    gribtab(65,3,1)%longname="vertically integrated ozone content"
    gribtab(65,3,1)%unit="Pa(O3)"
    gribtab(33,1,20)%shortname="U"
    gribtab(33,1,20)%longname="zonal wind"
    gribtab(33,1,20)%unit="m/s"
    gribtab(34,1,20)%shortname="V"
    gribtab(34,1,20)%longname="meridional wind"
    gribtab(34,1,20)%unit="m/s"
    gribtab(40,1,19)%shortname="W"
    gribtab(40,1,19)%longname="vertical wind (geometric)"
    gribtab(40,1,19)%unit="m/s"
    gribtab(1,1,20)%shortname="P"
    gribtab(1,1,20)%longname="pressure"
    gribtab(1,1,20)%unit="Pa"
    gribtab(139,2,20)%shortname="PP"
    gribtab(139,2,20)%longname="pressure perturbation"
    gribtab(139,2,20)%unit="Pa"
    gribtab(11,1,20)%shortname="T"
    gribtab(11,1,20)%longname="temperature"
    gribtab(11,1,20)%unit="K"
    gribtab(51,1,20)%shortname="QV"
    gribtab(51,1,20)%longname="specific moisture"
    gribtab(51,1,20)%unit="kg/kg"
    gribtab(31,2,20)%shortname="QC"
    gribtab(31,2,20)%longname="specific cloud water content"
    gribtab(31,2,20)%unit="kg/kg"
    gribtab(33,2,20)%shortname="QI"
    gribtab(33,2,20)%longname="specific cloud ice content"
    gribtab(33,2,20)%unit="kg/kg"
    gribtab(35,2,20)%shortname="QR"
    gribtab(35,2,20)%longname="specific rain water content"
    gribtab(35,2,20)%unit="kg/kg"
    gribtab(36,2,20)%shortname="QS"
    gribtab(36,2,20)%longname="specific snow water content"
    gribtab(36,2,20)%unit="kg/kg"
    gribtab(29,2,20)%shortname="CLC"
    gribtab(29,2,20)%longname="cloud cover"
    gribtab(29,2,20)%unit="%"
    gribtab(152,2,19)%shortname="TKE"
    gribtab(152,2,19)%longname="turbulent kinetic energy"
    gribtab(152,2,19)%unit="m2/s2"
    gribtab(153,2,19)%shortname="TKVM"
    gribtab(153,2,19)%longname="turbulent diffusion coefficient for vertical transport of momentum"
    gribtab(153,2,19)%unit="m2/s"
    gribtab(154,2,19)%shortname="TKVH"
    gribtab(154,2,19)%longname="turbulent diffusion coefficient for vertical transport of heat and moisture"
    gribtab(154,2,19)%unit="m2/s"
    gribtab(1,1,1)%shortname="PS"
    gribtab(1,1,1)%longname="unreduced surface pressure"
    gribtab(1,1,1)%unit="Pa"
    gribtab(2,1,12)%shortname="PMSL"
    gribtab(2,1,12)%longname="surface pressure reduced to sea level"
    gribtab(2,1,12)%unit="Pa"
    gribtab(203,2,1)%shortname="T_SNOW"
    gribtab(203,2,1)%longname="snow temperature"
    gribtab(203,2,1)%unit="K"
    gribtab(85,1,21)%shortname="T_S"
    gribtab(85,1,21)%longname="soil temperature (level 0 = surface layer)"
    gribtab(85,1,21)%unit="K"
    gribtab(11,1,1)%shortname="T_G"
    gribtab(11,1,1)%longname="earth surface temperature"
    gribtab(11,1,1)%unit="K"
    gribtab(85,1,21)%shortname="T_M"
    gribtab(85,1,21)%longname="temperature upper boundary of middle surface layer"
    gribtab(85,1,21)%unit="K"
    gribtab(51,1,1)%shortname="QV_S"
    gribtab(51,1,1)%longname="surface specific moisture"
    gribtab(51,1,1)%unit="kg/kg"
    gribtab(65,1,1)%shortname="W_SNOW"
    gribtab(65,1,1)%longname="water content of snow cover"
    gribtab(65,1,1)%unit="kg/m2"
    gribtab(133,2,1)%shortname="RHO_SNOW"
    gribtab(133,2,1)%longname="snow density"
    gribtab(133,2,1)%unit="kg/m3"
    gribtab(66,1,1)%shortname="H_SNOW"
    gribtab(66,1,1)%longname="height of snow cover"
    gribtab(66,1,1)%unit="m"
    gribtab(200,2,1)%shortname="W_I"
    gribtab(200,2,1)%longname="water content of interception storage"
    gribtab(200,2,1)%unit="kg/m2"
    gribtab(86,1,22)%shortname="W_G"
    gribtab(86,1,22)%longname="water content of soil layer"
    gribtab(86,1,22)%unit="kg/m2"
    gribtab(84,1,1)%shortname="ALB_RAD"
    gribtab(84,1,1)%longname="surface albedo"
    gribtab(84,1,1)%unit="%"
    gribtab(129,2,1)%shortname="FRESHSNW"
    gribtab(129,2,1)%longname="indicator of snow aging (for albedo)"
    gribtab(129,2,1)%unit="1"
    gribtab(111,1,1)%shortname="ASOB_S"
    gribtab(111,1,1)%longname="net short wave radiation at surface"
    gribtab(111,1,1)%unit="W/m2"
    gribtab(112,1,1)%shortname="ATHB_S"
    gribtab(112,1,1)%longname="net long wave radiation at surface"
    gribtab(112,1,1)%unit="W/m2"
    gribtab(5,2,1)%shortname="APAB_S"
    gribtab(5,2,1)%longname="net active radiation of photosynthesis at surface"
    gribtab(5,2,1)%unit="W/m2"
    gribtab(22,2,1)%shortname="ASWDIR_S"
    gribtab(22,2,1)%longname="direct short wave radiation"
    gribtab(22,2,1)%unit="W/m2"
    gribtab(23,2,1)%shortname="ASWDIFD_S"
    gribtab(23,2,1)%longname="diffuse downward short wave radiation at surface"
    gribtab(23,2,1)%unit="W/m2"
    gribtab(24,2,1)%shortname="ASWDIFU_S"
    gribtab(24,2,1)%longname="diffuse upwards short wave radiation at surface"
    gribtab(24,2,1)%unit="W/m2"
    gribtab(113,1,8)%shortname="ASOB_T"
    gribtab(113,1,8)%longname="net short wave radiation at upper model boundary"
    gribtab(113,1,8)%unit="W/m2"
    gribtab(114,1,8)%shortname="ATHB_T"
    gribtab(114,1,8)%longname="net long wave radiation at upper model boundary"
    gribtab(114,1,8)%unit="W/m2"
    gribtab(102,2,1)%shortname="RAIN_GSP"
    gribtab(102,2,1)%longname="grid scale rain (surface)"
    gribtab(102,2,1)%unit="kg/m2"
    gribtab(79,1,1)%shortname="SNOW_GSP"
    gribtab(79,1,1)%longname="grid scale snow (surface)"
    gribtab(79,1,1)%unit="kg/m2"
    gribtab(113,2,1)%shortname="RAIN_CON"
    gribtab(113,2,1)%longname="convective rain (surface)"
    gribtab(113,2,1)%unit="kg/m2"
    gribtab(78,1,1)%shortname="SNOW_CON"
    gribtab(78,1,1)%longname="convective snow (surface)"
    gribtab(78,1,1)%unit="kg/m2"
    gribtab(61,1,1)%shortname="TOT_PREC"
    gribtab(61,1,1)%longname="total recipitation (surface)"
    gribtab(61,1,1)%unit="kg/m2"
    gribtab(90,1,22)%shortname="RUNOFF"
    gribtab(90,1,22)%longname="runoff in soil layer (levels = layer boundaries)"
    gribtab(90,1,22)%unit="kg/m2"
    gribtab(57,1,1)%shortname="AEVAP_S"
    gribtab(57,1,1)%longname="moisture flux at surface"
    gribtab(57,1,1)%unit="kg/m2"
    gribtab(42,2,1)%shortname="TDIV_HUM"
    gribtab(42,2,1)%longname="vertically integrated specific moisture divergence"
    gribtab(42,2,1)%unit="kg/m2"
    gribtab(41,2,1)%shortname="TWATER"
    gribtab(41,2,1)%longname="vertically integrated water"
    gribtab(41,2,1)%unit="kg/m2"
    gribtab(54,1,1)%shortname="TQV"
    gribtab(54,1,1)%longname="vertically integrated water vapor"
    gribtab(54,1,1)%unit="kg/m2"
    gribtab(76,1,1)%shortname="TQC"
    gribtab(76,1,1)%longname="vertically integrated cloud water"
    gribtab(76,1,1)%unit="kg/m2"
    gribtab(58,1,1)%shortname="TQI"
    gribtab(58,1,1)%longname="vertically integrated cloud ice"
    gribtab(58,1,1)%unit="kg/m2"
    gribtab(33,1,15)%shortname="U_2M"
    gribtab(33,1,15)%longname="zonal wind at 10m height"
    gribtab(33,1,15)%unit="m/s"
    gribtab(34,1,15)%shortname="V_2M"
    gribtab(34,1,15)%longname="meridional wind at 10m height"
    gribtab(34,1,15)%unit="m/s"
    gribtab(11,1,15)%shortname="T_2M"
    gribtab(11,1,15)%longname="temperature at 2m height"
    gribtab(11,1,15)%unit="K"
    gribtab(17,1,15)%shortname="TD_2M"
    gribtab(17,1,15)%longname="dew point temperture at 2m height"
    gribtab(17,1,15)%unit="K"
    gribtab(52,1,15)%shortname="RELHUM_2M"
    gribtab(52,1,15)%longname="relative humidity at 2m height"
    gribtab(52,1,15)%unit="%"
    gribtab(16,1,15)%shortname="TMIN_2M"
    gribtab(16,1,15)%longname="minimum temperature at 2m height"
    gribtab(16,1,15)%unit="K"
    gribtab(15,1,15)%shortname="TMAX_2M"
    gribtab(15,1,15)%longname="maximum temperature at 2m height"
    gribtab(15,1,15)%unit="K"
    gribtab(187,2,15)%shortname="VMAX_10M"
    gribtab(187,2,15)%longname="maximum wind speed at 2m height"
    gribtab(187,2,15)%unit="m/s"
    gribtab(71,1,1)%shortname="CLCT"
    gribtab(71,1,1)%longname="total cloud cover"
    gribtab(71,1,1)%unit="%"
    gribtab(75,1,1)%shortname="CLCH"
    gribtab(75,1,1)%longname="high cloud cover (0-400 hPa)"
    gribtab(75,1,1)%unit="%"
    gribtab(74,1,1)%shortname="CLCM"
    gribtab(74,1,1)%longname="middle cloud cover(400-800 hPa)"
    gribtab(74,1,1)%unit="%"
    gribtab(73,1,1)%shortname="CLCL"
    gribtab(73,1,1)%longname="low cloud cover (800 hPa-surface)"
    gribtab(73,1,1)%unit="%"
    gribtab(203,4,1)%shortname="CLDEPTH"
    gribtab(203,4,1)%longname="modificated cloud depth"
    gribtab(203,4,1)%unit="1"
    gribtab(204,4,1)%shortname="CLCT_MOD"
    gribtab(204,4,1)%longname="modificated total cloud cover"
    gribtab(204,4,1)%unit="1"
    gribtab(68,2,2)%shortname="HBAS_CON"
    gribtab(68,2,2)%longname="convective cloud base height above sea level"
    gribtab(68,2,2)%unit="m"
    gribtab(69,2,3)%shortname="HTOP_CON"
    gribtab(69,2,3)%longname="convective cloud top above sea level"
    gribtab(69,2,3)%unit="m"
    gribtab(72,2,1)%shortname="BAS_CON"
    gribtab(72,2,1)%longname="level index of convective cloud base"
    gribtab(72,2,1)%unit="-"
    gribtab(73,2,1)%shortname="TOP_CON"
    gribtab(73,2,1)%longname="level index of convective cloud top"
    gribtab(73,2,1)%unit="-"
    gribtab(82,2,1)%shortname="HTOP_DC"
    gribtab(82,2,1)%longname="upper boundary height of dry convection above sea level"
    gribtab(82,2,1)%unit="m"
    gribtab(84,2,4)%shortname="HZEROCL"
    gribtab(84,2,4)%longname="height of 0 degree Celsius above sea level"
    gribtab(84,2,4)%unit="m"
    gribtab(85,2,4)%shortname="SNOWLMT"
    gribtab(85,2,4)%longname="height of snow line above sea level"
    gribtab(85,2,4)%unit="m"
    gribtab(240,2,1)%shortname="MFLX_CON"
    gribtab(240,2,1)%longname="mass flux at convective cloud base"
    gribtab(240,2,1)%unit="kg/(sm2)"
    gribtab(241,2,1)%shortname="CAPE_CON"
    gribtab(241,2,1)%longname="convective available potential energy"
    gribtab(241,2,1)%unit="J/kg"
    gribtab(243,2,1)%shortname="QCVG_CON"
    gribtab(243,2,1)%longname="moisture convergence below convective cloud base"
    gribtab(243,2,1)%unit="1/s"
    gribtab(147,2,1)%shortname="TKE_CON"
    gribtab(147,2,1)%longname="convective turbulent kinetic energy"
    gribtab(147,2,1)%unit="J/kg"
    gribtab(145,2,1)%shortname="CAPE_ML"
    gribtab(145,2,1)%longname="mixed layer CAPE"
    gribtab(145,2,1)%unit="J/kg"
    gribtab(146,2,1)%shortname="CIN_ML"
    gribtab(146,2,1)%longname="mixed layer convective inhibition"
    gribtab(146,2,1)%unit="J/kg"
    gribtab(124,1,1)%shortname="AUMFL_S"
    gribtab(124,1,1)%longname="U-momentum at surface"
    gribtab(124,1,1)%unit="N/m2"
    gribtab(125,1,1)%shortname="AVMFL_S"
    gribtab(125,1,1)%longname="V-momentum at surface"
    gribtab(125,1,1)%unit="N/m2"
    gribtab(122,1,1)%shortname="ASHFL_S"
    gribtab(122,1,1)%longname="sensible heat flux at surface"
    gribtab(122,1,1)%unit="W/m2"
    gribtab(121,1,1)%shortname="ALHFL"
    gribtab(121,1,1)%longname="S latent heat flux at surface"
    gribtab(121,1,1)%unit="W/m2"
    gribtab(170,2,1)%shortname="TCM"
    gribtab(170,2,1)%longname="turbulent momentum transfer coefficient at surface"
    gribtab(170,2,1)%unit="-"
    gribtab(171,2,1)%shortname="TCH"
    gribtab(171,2,1)%longname="turbulent heat and moisture transfer coefficient at surface"
    gribtab(171,2,1)%unit="-"
    gribtab(83,1,1)%shortname="Z0"
    gribtab(83,1,1)%longname="roughness length"
    gribtab(83,1,1)%unit="m"
    gribtab(215,2,1)%shortname="T_ICE"
    gribtab(215,2,1)%longname="sea ice surface temperature"
    gribtab(215,2,1)%unit="K"
    gribtab(92,1,1)%shortname="H_ICE"
    gribtab(92,1,1)%longname="sea ice thickness"
    gribtab(92,1,1)%unit="m"
    gribtab(193,2,1)%shortname="T_WML_LK"
    gribtab(193,2,1)%longname="lake mixing layer temperture"
    gribtab(193,2,1)%unit="K"
    gribtab(95,2,1)%shortname="H_ML_LK"
    gribtab(95,2,1)%longname="lake mixing layer thickness"
    gribtab(95,2,1)%unit="m"
    gribtab(194,2,1)%shortname="T_MNW_LK"
    gribtab(194,2,1)%longname="mean temperture of lake water column"
    gribtab(194,2,1)%unit="K"
    gribtab(191,2,1)%shortname="T_BOT_LK"
    gribtab(191,2,1)%longname="temperature of lake bed"
    gribtab(191,2,1)%unit="K"
    gribtab(91,2,1)%shortname="C_T_LK"
    gribtab(91,2,1)%longname="form factor of lake model"
    gribtab(91,2,1)%unit="-"
    gribtab(33,1,10)%shortname="U"
    gribtab(33,1,10)%longname="zonal wind"
    gribtab(33,1,10)%unit="m/s"
    gribtab(34,1,10)%shortname="V"
    gribtab(34,1,10)%longname="meridional Wind"
    gribtab(34,1,10)%unit="m/s"
    gribtab(39,1,10)%shortname="OMEGA"
    gribtab(39,1,10)%longname="vertical motion (pressure)"
    gribtab(39,1,10)%unit="Pa/s"
    gribtab(6,1,10)%shortname="FI"
    gribtab(6,1,10)%longname="geopotential"
    gribtab(6,1,10)%unit="m2/s2"
    gribtab(11,1,10)%shortname="T"
    gribtab(11,1,10)%longname="temperature"
    gribtab(11,1,10)%unit="K"
    gribtab(52,1,10)%shortname="RELHUM"
    gribtab(52,1,10)%longname="relative humidity"
    gribtab(52,1,10)%unit="%"
    gribtab(33,1,13)%shortname="U"
    gribtab(33,1,13)%longname="zonal wind"
    gribtab(33,1,13)%unit="m/s"
    gribtab(34,1,13)%shortname="V"
    gribtab(34,1,13)%longname="meridional wind"
    gribtab(34,1,13)%unit="m/s"
    gribtab(40,1,13)%shortname="W"
    gribtab(40,1,13)%longname="vertikal wind (geometric)"
    gribtab(40,1,13)%unit="m/s"
    gribtab(1,1,13)%shortname="P"
    gribtab(1,1,13)%longname="pressure"
    gribtab(1,1,13)%unit="Pa"
    gribtab(11,1,13)%shortname="T"
    gribtab(11,1,13)%longname="temperature"
    gribtab(11,1,13)%unit="K"
    gribtab(52,1,13)%shortname="RELHUM"
    gribtab(52,1,13)%longname="relative humidity"
    gribtab(52,1,13)%unit="%"
    gribtab(197,2,21)%shortname="T_SO"
    gribtab(197,2,21)%longname="multi layer soil temperture"
    gribtab(197,2,21)%unit="K"
    gribtab(198,2,21)%shortname="W_SO"
    gribtab(198,2,21)%longname="multi layer soil water content"
    gribtab(198,2,21)%unit="kg/m2"
    gribtab(199,2,21)%shortname="W"
    gribtab(199,2,21)%longname="SO multi layer soil ice content"
    gribtab(199,2,21)%unit="kg/m2"
    gribtab(173,2,1)%shortname="MH"
    gribtab(173,2,1)%longname="mixing layer height above ground"
    gribtab(173,2,1)%unit="m"
    gribtab(99,4,1)%shortname="WW"
    gribtab(99,4,1)%longname="interpreted weather (WMO-key)"
    gribtab(99,4,1)%unit="-"

  END SUBROUTINE set_gribtab

  SUBROUTINE set_dartkinds()  

    dartkinds(:,:,:)=-1

    dartkinds(6,1,1)=KIND_SURFACE_GEOPOTENTIAL
    dartkinds(8,1,1)=KIND_SURFACE_ELEVATION
    dartkinds(33,1,20)=KIND_U_WIND_COMPONENT
    dartkinds(34,1,20)=KIND_V_WIND_COMPONENT
    dartkinds(40,1,19)=KIND_VERTICAL_VELOCITY
!    dartkinds(1,1,20)=KIND_PRESSURE
    dartkinds(139,2,20)=KIND_PRESSURE_PERTURBATION
    dartkinds(11,1,20)=KIND_TEMPERATURE
    dartkinds(51,1,20)=KIND_SPECIFIC_HUMIDITY
    dartkinds(31,2,20)=KIND_CLOUD_LIQUID_WATER
    dartkinds(33,2,20)=KIND_CLOUD_ICE
    dartkinds(29,2,20)=KIND_CLOUD_FRACTION
    dartkinds(1,1,1)=KIND_SURFACE_PRESSURE
    dartkinds(2,1,12)=KIND_SEA_SURFACE_PRESSURE
    dartkinds(86,1,22)=KIND_SOIL_MOISTURE
!    dartkinds(33,1,10)=KIND_U_WIND_COMPONENT
!    dartkinds(34,1,10)=KIND_V_WIND_COMPONENT
!    dartkinds(39,1,10)=KIND_VERTICAL_VELOCITY
!    dartkinds(6,1,10)=KIND_GEOPOTENTIAL_HEIGHT !!! geopotential (not height)
!    dartkinds(11,1,10)=KIND_TEMPERATURE
    dartkinds(52,1,10)=KIND_RELATIVE_HUMIDITY
!    dartkinds(33,1,13)=KIND_U_WIND_COMPONENT
!    dartkinds(34,1,13)=KIND_V_WIND_COMPONENT
!    dartkinds(1,1,13)=KIND_PRESSURE
!    dartkinds(11,1,13)=KIND_TEMPERATURE
    dartkinds(198,2,21)=KIND_SOIL_MOISTURE

  END SUBROUTINE set_dartkinds

  FUNCTION get_varname(tab,code,ltype)

    CHARACTER(len=256) :: get_varname(3)
    INTEGER,INTENT(in) :: tab
    INTEGER,INTENT(in) :: code
    INTEGER,INTENT(in) :: ltype

    INTEGER            :: ltab,lltype

    CALL set_gribtab()

    ltab=tab
    lltype=ltype

    if (tab==2) ltab=1
    if (tab>=200) ltab=ltab-199
    if (ltype>=100) lltype=lltype-90

    get_varname(1)=gribtab(code,ltab,lltype)%longname
    get_varname(2)=gribtab(code,ltab,lltype)%shortname
    get_varname(3)=gribtab(code,ltab,lltype)%unit

    RETURN

  END FUNCTION get_varname

  FUNCTION get_level(ltype,b)

    INTEGER                    :: get_level(1:3)
    INTEGER,INTENT(in)         :: ltype
    INTEGER(kind=1),INTENT(in) :: b(2)

    get_level(:)=-999.

    SELECT CASE (ltype)
      CASE (100)
        get_level(1)=concat_bytes1(b,2,.TRUE.)
        get_level(3)=1
      CASE (101)
        get_level(1:2)=b(1:2)
        get_level(3)=2
      CASE (102)
        get_level(1)=0
        get_level(3)=1
      CASE (103)
        get_level(1)=concat_bytes1(b,2,.TRUE.)
        get_level(3)=1
      CASE (104)
        get_level(1:2)=b(1:2)
        get_level(3)=2
      CASE (105)
        get_level(1)=concat_bytes1(b,2,.TRUE.)
        get_level(3)=1
      CASE (106)
        get_level(1:2)=b(1:2)
        get_level(3)=2
      CASE (107)
        get_level(1)=concat_bytes1(b,2,.TRUE.)
        get_level(3)=1
      CASE (108)
        get_level(1:2)=b(1:2)
        get_level(3)=2
      CASE (109)
        get_level(1)=concat_bytes1(b,2,.TRUE.)
        get_level(3)=1
      CASE (110)
        get_level(1:2)=b(1:2)
        get_level(3)=2
      CASE (111)
        get_level(1)=concat_bytes1(b,2,.TRUE.)
        get_level(3)=1
      CASE (112)
        get_level(1:2)=b(1:2)
        get_level(3)=2
      CASE (113)
        get_level(1)=concat_bytes1(b,2,.TRUE.)
        get_level(3)=1
      CASE (114)
        get_level(1:2)=b(1:2)
        get_level(3)=2
      CASE (121)
        get_level(1:2)=b(1:2)
        get_level(3)=2
      CASE (125)
        get_level(1)=concat_bytes1(b,2,.TRUE.)
        get_level(3)=1
      CASE (128)
        get_level(1:2)=b(1:2)
        get_level(3)=2
      CASE (141)
        get_level(1:2)=b(1:2)
        get_level(3)=2
      CASE (160)
        get_level(1)=concat_bytes1(b,2,.TRUE.)
        get_level(3)=1
      CASE (200)
        get_level(1)=0
        get_level(3)=1
      CASE (201)
        get_level(1)=0
        get_level(3)=1
    END SELECT

    RETURN
  END FUNCTION get_level

  FUNCTION get_dims(b)

    INTEGER                    :: get_dims(3)
    INTEGER(kind=1),INTENT(in) :: b(10)

    get_dims(1) = concat_bytes1(b(7:8),2,.FALSE.)
    get_dims(2) = concat_bytes1(b(9:10),2,.FALSE.)
    get_dims(3) = to_positive(b(4)) - 5  ! 4 parameters + 1 edge

    RETURN

  END FUNCTION get_dims

  FUNCTION get_dart_kind(tab,code,ltype)

    INTEGER            :: get_dart_kind
    INTEGER,INTENT(in) :: tab
    INTEGER,INTENT(in) :: code
    INTEGER,INTENT(in) :: ltype

    INTEGER            :: ltab,lltype

    CALL set_dartkinds()

    ltab=tab
    lltype=ltype

    if (tab==2) ltab=1
    if (tab>=200) ltab=ltab-199
    if (ltype>=100) lltype=lltype-90

    get_dart_kind=dartkinds(code,ltab,lltype)

    RETURN

  END FUNCTION get_dart_kind

END MODULE grib_info_mod
