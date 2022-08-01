! DART software - Copyright UCAR. This open source software is provided
! by UCAR, "as is", without charge, subject to all terms of use at
! http://www.image.ucar.edu/DAReS/DART/DART_download

! ! in this section, define the quantity of interest.  it must
! ! start with QTY_xxx and be less than 32 characters total.
! ! if there are units, min/max values, pdf, etc you can add one or
! ! more name=value pairs.  at this point there *cannot* be spaces
! ! around the equals sign.  the value can have spaces if it is
! ! surrounded by double quotes.  e.g. desc="assumes dry air density"
! ! trailing comments can be added if they start with !
!
! BEGIN DART PREPROCESS QUANTITY DEFINITIONS
!
!   QTY_U_WIND_COMPONENT
!   QTY_V_WIND_COMPONENT
!   QTY_SURFACE_PRESSURE
!   QTY_TEMPERATURE
!   QTY_SPECIFIC_HUMIDITY
!   QTY_PRESSURE
!   QTY_VERTICAL_VELOCITY
!   QTY_RAINWATER_MIXING_RATIO
!   QTY_DEWPOINT
!   QTY_DENSITY
!   QTY_MEAN_SOURCE
!   QTY_VELOCITY
!   QTY_RADAR_REFLECTIVITY
!   QTY_GRAUPEL_MIXING_RATIO
!   QTY_SNOW_MIXING_RATIO
!   QTY_GPSRO
!   QTY_CLOUD_LIQUID_WATER
!   QTY_CLOUD_ICE
!   QTY_CONDENSATIONAL_HEATING
!   QTY_VAPOR_MIXING_RATIO
!   QTY_ICE_NUMBER_CONCENTRATION
!   QTY_GEOPOTENTIAL_HEIGHT
!   QTY_POTENTIAL_TEMPERATURE
!   QTY_SOIL_MOISTURE
!   QTY_SOURCE_PHASE
!   QTY_SURFACE_ELEVATION
!
! ! kinds for Gravity Wave Drag (CAM - kevin)
!
!   QTY_GRAV_WAVE_DRAG_EFFIC
!   QTY_GRAV_WAVE_STRESS_FRACTION
!
!   QTY_POWER_WEIGHTED_FALL_SPEED
!   QTY_CLOUDWATER_MIXING_RATIO
!   QTY_ICE_MIXING_RATIO
!   QTY_DROPLET_NUMBER_CONCENTR
!   QTY_SNOW_NUMBER_CONCENTR
!   QTY_RAIN_NUMBER_CONCENTR
!   QTY_GRAUPEL_NUMBER_CONCENTR
!   QTY_CLOUD_FRACTION
!   QTY_ICE_FRACTION
!   QTY_RELATIVE_HUMIDITY
!
! ! kinds for generic parameters that aren't going to be
! ! directly observed but are going to be adjusted by the
! ! assimilation.
!
!   QTY_1D_PARAMETER
!   QTY_2D_PARAMETER
!   QTY_3D_PARAMETER
!
! ! kinds for tendencies
!
!   QTY_ALTIMETER_TENDENCY
!
! ! kind for precip water; contrast with
! ! total precip water (also in this file),
! ! which is the total column integrated value.
!
!   QTY_PRECIPITABLE_WATER
!
! ! proposed new kinds for COSMIC GPS/RO obs (currently unused)
!
!   QTY_OCCULTATION_REFRACTIVITY
!   QTY_OCCULTATION_EXCESSPHASE
!
! ! kind for the measuring elevation - contrast this with geopotential height
!
!   QTY_GEOMETRIC_HEIGHT
!
!   QTY_INFRARED_RADIANCE
!   QTY_INFRARED_BRIGHT_TEMP
!   QTY_LANDMASK
!
! ! kind for unstructured grids
!
!   QTY_EDGE_NORMAL_SPEED
!
! ! kind for cloud liquid water path
!
!   QTY_CLW_PATH
!   QTY_CWP_PATH
!   QTY_CWP_PATH_ZERO
!
! ! kind for wind power
!
!   QTY_WIND_TURBINE_POWER
!
! ! kinds for surface fields
!
!   QTY_2M_SPECIFIC_HUMIDITY
!   QTY_2M_TEMPERATURE
!   QTY_10M_U_WIND_COMPONENT
!   QTY_10M_V_WIND_COMPONENT
!
! ! kinds for planetary remote sensing (wglawson, c.lee)
!
!   QTY_SKIN_TEMPERATURE
!   QTY_NADIR_RADIANCE
!   QTY_TRACER_1_MIXING_RATIO
!   QTY_TRACER_2_MIXING_RATIO
!   QTY_TRACER_CONCENTRATION
!   QTY_TRACER_SOURCE
!   QTY_SOIL_TEMPERATURE
!   QTY_SOIL_LIQUID_WATER
!
! ! kinds for NCOMMAS  (Lou W., Ted M.)
!
!   QTY_VERTICAL_VORTICITY
!
! ! more kinds for planetary remote sensing (c.lee)
!
!   QTY_SURFACE_ALBEDO
!   QTY_SURFACE_EMISSIVITY
!   QTY_DUST_OPACITY_7MB
!   QTY_THC
!
! ! kinds for vortex tracking (WRF - yongsheng)
!
!   QTY_VORTEX_LON
!   QTY_VORTEX_LAT
!   QTY_VORTEX_PMIN
!   QTY_VORTEX_WMAX
!
! ! kinds for COAMPS (Tim Whitcomb)
!
!   QTY_EXNER_FUNCTION
!   QTY_TURBULENT_KINETIC_ENERGY
!   QTY_TOTAL_PRECIPITABLE_WATER
!   QTY_VERTLEVEL
!   QTY_MICROWAVE_BRIGHT_TEMP
!
! ! kinds for NAAPS (Walter R. Sessions)
!
!   QTY_INTEGRATED_SULFATE
!   QTY_INTEGRATED_DUST
!   QTY_INTEGRATED_SMOKE
!   QTY_INTEGRATED_SEASALT
!   QTY_INTEGRATED_AOD
!   QTY_SO2
!   QTY_SULFATE
!   QTY_DUST
!   QTY_SMOKE
!   QTY_SEASALT
!
! ! kinds for ZVD (advanced microphysics)
!
!   QTY_HAIL_MIXING_RATIO
!   QTY_HAIL_NUMBER_CONCENTR
!   QTY_GRAUPEL_VOLUME
!   QTY_HAIL_VOLUME
!   QTY_DIFFERENTIAL_REFLECTIVITY
!   QTY_SPECIFIC_DIFFERENTIAL_PHASE
!   QTY_FLASH_RATE_2D
!
! ! kinds for radiance
!
!   QTY_RADIANCE                        ! L1 radiance (mW/cm^-1/sr/m^2)
!   QTY_BRIGHTNESS_TEMPERATURE          ! L1 brightness temperature (K)
!   QTY_BI_DIRECTIONAL_REFLECTANCE      ! L1 bi-directional reflectance (BDRF, unitless)
!   QTY_SURFACE_TYPE                    ! land = 0, sea = 1, seaice = 2
!   QTY_WIND_FETCH                      ! Wind fetch, m
!   QTY_WATER_TYPE                      ! fresh = 0, ocean = 1
!   QTY_FOAM_FRAC                       ! Fraction of foam on ocean surface (0-1)
!   QTY_INSOLUBLE_AER                   ! Insoluble aerosol OPAC aerosol (INSO)
!   QTY_H2O_SOLUBLE_AER                 ! Soluble aerosol OPAC aerosol (WASO)
!   QTY_SOOT                            ! Soot aerosol OPAC aerosol (SOOT)
!   QTY_SEASALT_ACCUM                   ! Sea salt (accumulation mode) OPAC aerosol (SSAM)
!   QTY_SEASALT_COARSE                  ! Sea salt (coarse) OPAC aerosol (SSCM)
!   QTY_MINERAL_NUCLEUS                 ! Mineral (nucleus) OPAC aerosol (MINM)
!   QTY_MINERAL_ACCUM                   ! Mineral (accumulation mode) OPAC aerosol (MIAM)
!   QTY_MINERAL_COARSE                  ! Mineral (coarse mode) OPAC aerosol (MICM)
!   QTY_MINERAL_TRANSPORTED             ! Mineral (transported mode) OPAC aerosol (MITR)
!   QTY_SULPHATED_DROPS                 ! Sulphated droplets OPAC aerosol (SUSO)
!   QTY_VOLCANIC_ASH                    ! Volcanic ash OPAC aerosol (VOLA)
!   QTY_NEW_VOLCANIC_ASH                ! New volcanic ash OPAC aerosol (VAPO)
!   QTY_ASIAN_DUST                      ! Asian dust OPAC aerosol (ASDU)
!   QTY_BLACK_CARBON                    ! Black carbon CAMS aerosol (BCAR)
!   QTY_DUST_BIN1                       ! Dust bin 1 CAMS aerosol (DUS1)
!   QTY_DUST_BIN2                       ! Dust bin 2 CAMS aerosol (DUS2)
!   QTY_DUST_BIN3                       ! Dust bin 3 CAMS aerosol (DUS3)
!   QTY_AMMONIUM_SULPHATE               ! Ammonium sulphate CAMS aerosol (SULP)
!   QTY_SEA_SALT_BIN1                   ! Sea salt bin 1 CAMS aerosol (SSA1)
!   QTY_SEA_SALT_BIN2                   ! Sea salt bin 2 CAMS aerosol (SSA2)
!   QTY_SEA_SALT_BIN3                   ! Sea salt bin 3 CAMS aerosol (SSA3)
!   QTY_HYDROPHILIC_ORGANIC_MATTER      ! Hydrophilic organic matter CAMS aerosol (OMAT)
!   QTY_CLOUDWATER_DE                   ! Cloud liquid water effective diameter (microns)
!   QTY_CLOUD_ICE_DE                    ! Cloud ice effective diameter (microns)
!   QTY_COLUMN_CLOUD_FRAC               ! Simple cloud fraction (0-1)
!   QTY_CLOUD_TOP_PRESSURE
!   QTY_ABSOLUTE_HUMIDITY
!
! END DART PREPROCESS QUANTITY DEFINITIONS

