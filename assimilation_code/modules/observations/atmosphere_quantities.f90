! DART software - Copyright UCAR. This open source software is provided
! by UCAR, "as is", without charge, subject to all terms of use at
! http://www.image.ucar.edu/DAReS/DART/DART_download
!
! $Id$


 
! ! in this section you can have comments (which need a second !) or
! ! lines with two words on them:      QTY_xxx    "m/s"
! ! or lines with 4 words on them:     QTY_xxx    "hPa"   0.0   MISSING_R8
! ! can have trailing comments after:  QTY_xxx  "m"  ! comment
!
! BEGIN DART PREPROCESS QUANTITY DEFINITIONS
!
!     QTY_STATE_VARIABLE                "none"   MISSING_R8   MISSING_R8
!     QTY_U_WIND_COMPONENT              "m/s"    MISSING_R8   MISSING_R8
!     QTY_V_WIND_COMPONENT              "m/s"    MISSING_R8   MISSING_R8
!     QTY_SURFACE_PRESSURE              "hPa"    0.0          MISSING_R8
!     QTY_TEMPERATURE                   "K"      0.0          MISSING_R8
!     QTY_SPECIFIC_HUMIDITY             "none"   MISSING_R8   MISSING_R8
!     QTY_PRESSURE                      "hPa"    0.0          MISSING_R8
!     QTY_VERTICAL_VELOCITY             "m/s"    MISSING_R8   MISSING_R8
!     QTY_RAINWATER_MIXING_RATIO        "none"   MISSING_R8   MISSING_R8
!     QTY_DEWPOINT                      "none"   MISSING_R8   MISSING_R8
!     QTY_DENSITY                       "none"   MISSING_R8   MISSING_R8
!     QTY_VELOCITY                      "m/s"    MISSING_R8   MISSING_R8
!     QTY_RADAR_REFLECTIVITY            "none"   MISSING_R8   MISSING_R8
!     QTY_GRAUPEL_MIXING_RATIO          "none"   MISSING_R8   MISSING_R8
!     QTY_SNOW_MIXING_RATIO             "none"   MISSING_R8   MISSING_R8
!     QTY_GPSRO                         "none"   MISSING_R8   MISSING_R8
!     QTY_CLOUD_LIQUID_WATER            "none"   MISSING_R8   MISSING_R8
!     QTY_CLOUD_ICE                     "none"   MISSING_R8   MISSING_R8
!     QTY_CONDENSATIONAL_HEATING        "none"   MISSING_R8   MISSING_R8
!     QTY_VAPOR_MIXING_RATIO            "none"   MISSING_R8   MISSING_R8
!     QTY_ICE_NUMBER_CONCENTRATION      "none"   MISSING_R8   MISSING_R8
!     QTY_GEOPOTENTIAL_HEIGHT           "m"      MISSING_R8   MISSING_R8
!     QTY_POTENTIAL_TEMPERATURE         "none"   MISSING_R8   MISSING_R8
!     QTY_SOIL_MOISTURE                 "none"   MISSING_R8   MISSING_R8
!     QTY_SURFACE_ELEVATION             "m"      MISSING_R8   MISSING_R8
! 
! ! kinds for Gravity Wave Drag (CAM - kevin)
!     QTY_GRAV_WAVE_DRAG_EFFIC          "none"   MISSING_R8   MISSING_R8
!     QTY_GRAV_WAVE_STRESS_FRACTION     "none"   MISSING_R8   MISSING_R8
! 
!     QTY_POWER_WEIGHTED_FALL_SPEED     "none"   MISSING_R8   MISSING_R8
!     QTY_CLOUDWATER_MIXING_RATIO       "none"   MISSING_R8   MISSING_R8
!     QTY_ICE_MIXING_RATIO              "none"   MISSING_R8   MISSING_R8
!     QTY_DROPLET_NUMBER_CONCENTR       "none"   MISSING_R8   MISSING_R8
!     QTY_SNOW_NUMBER_CONCENTR          "none"   MISSING_R8   MISSING_R8
!     QTY_RAIN_NUMBER_CONCENTR          "none"   MISSING_R8   MISSING_R8
!     QTY_GRAUPEL_NUMBER_CONCENTR       "none"   MISSING_R8   MISSING_R8
!     QTY_CLOUD_FRACTION                "none"   MISSING_R8   MISSING_R8
!     QTY_ICE_FRACTION                  "none"   MISSING_R8   MISSING_R8
!     QTY_RELATIVE_HUMIDITY             "none"   MISSING_R8   MISSING_R8
! 
! ! kinds for generic parameters that aren't going to be
! ! directly observed but are going to be adjusted by the
! ! assimilation.
!     QTY_1D_PARAMETER                  "none"   MISSING_R8   MISSING_R8
!     QTY_2D_PARAMETER                  "none"   MISSING_R8   MISSING_R8
!     QTY_3D_PARAMETER                  "none"   MISSING_R8   MISSING_R8
! 
! ! kinds for tendencies
!     QTY_ALTIMETER_TENDENCY            "none"   MISSING_R8   MISSING_R8
! 
! ! kind for precip water; contrast with
! ! total precip water (also in this file), 
! ! which is the total column integrated value. 
!     QTY_PRECIPITABLE_WATER            "none"   MISSING_R8   MISSING_R8
! 
! ! proposed new kinds for COSMIC GPS/RO obs
! ! (currently unused)
!     QTY_OCCULTATION_REFRACTIVITY      "none"   MISSING_R8   MISSING_R8
!     QTY_OCCULTATION_EXCESSPHASE       "none"   MISSING_R8   MISSING_R8
! 
! ! kind for the other way of measuring elevation
! ! contrast this with geopotential height
!     QTY_GEOMETRIC_HEIGHT              "m"      MISSING_R8   MISSING_R8
! 
!     QTY_INFRARED_RADIANCE             "none"   MISSING_R8   MISSING_R8
!     QTY_INFRARED_BRIGHT_TEMP          "none"   MISSING_R8   MISSING_R8
!     QTY_LANDMASK                      "none"   MISSING_R8   MISSING_R8
! 
! ! kind for unstructured grids
!     QTY_EDGE_NORMAL_SPEED             "m/s"    MISSING_R8   MISSING_R8
! 
! ! kind for cloud liquid water path
!     QTY_CLW_PATH                      "none"   MISSING_R8   MISSING_R8
!     QTY_CWP_PATH                      "none"   MISSING_R8   MISSING_R8
!     QTY_CWP_PATH_ZERO                 "none"   MISSING_R8   MISSING_R8
! 
! 
! ! kind for wind power
!     QTY_WIND_TURBINE_POWER            "none"   MISSING_R8   MISSING_R8
! 
! ! kinds for surface fields
!     QTY_2M_SPECIFIC_HUMIDITY          "none"   MISSING_R8   MISSING_R8
!     QTY_2M_TEMPERATURE                "K"      MISSING_R8   MISSING_R8
!     QTY_10M_U_WIND_COMPONENT          "m/s"    MISSING_R8   MISSING_R8
!     QTY_10M_V_WIND_COMPONENT          "m/s"    MISSING_R8   MISSING_R8
! 
! ! kinds for planetary remote sensing (wglawson, c.lee)
!     QTY_SKIN_TEMPERATURE              "K"      MISSING_R8   MISSING_R8
!     QTY_NADIR_RADIANCE                "none"   MISSING_R8   MISSING_R8
!     QTY_TRACER_1_MIXING_RATIO         "none"   MISSING_R8   MISSING_R8
!     QTY_TRACER_2_MIXING_RATIO         "none"   MISSING_R8   MISSING_R8
!     ! this or QTY_TRACER_CONCENTRATION?
!     QTY_SOIL_TEMPERATURE              "K"      MISSING_R8   MISSING_R8
!     QTY_SOIL_LIQUID_WATER             "none"   MISSING_R8   MISSING_R8
! 
! ! kinds for NCOMMAS  (Lou W., Ted M.)
!     QTY_VERTICAL_VORTICITY            "none"   MISSING_R8   MISSING_R8
! 
! ! more kinds for planetary remote sensing (c.lee)
!     QTY_SURFACE_ALBEDO                "none"   MISSING_R8   MISSING_R8
!     QTY_SURFACE_EMISSIVITY            "none"   MISSING_R8   MISSING_R8
!     QTY_DUST_OPACITY_7MB              "none"   MISSING_R8   MISSING_R8
!     QTY_THC                           "none"   MISSING_R8   MISSING_R8
! 
! ! kinds for vortex tracking (WRF - yongsheng)
!     QTY_VORTEX_LON                    "degrees" MISSING_R8   MISSING_R8
!     QTY_VORTEX_LAT                    "degrees" MISSING_R8   MISSING_R8
!     QTY_VORTEX_PMIN                   "hPa"     MISSING_R8   MISSING_R8
!     QTY_VORTEX_WMAX                   "m/s"     MISSING_R8   MISSING_R8
! 
! ! kinds for COAMPS (Tim Whitcomb)
!     QTY_EXNER_FUNCTION                "none"   MISSING_R8   MISSING_R8
!     QTY_TURBULENT_KINETIC_ENERGY      "none"   MISSING_R8   MISSING_R8
!     QTY_TOTAL_PRECIPITABLE_WATER      "none"   MISSING_R8   MISSING_R8
!     QTY_VERTLEVEL                     "none"   MISSING_R8   MISSING_R8
!     QTY_MICROWAVE_BRIGHT_TEMP         "none"   MISSING_R8   MISSING_R8
! 
! ! kinds for NAAPS (Walter R. Sessions)
!     QTY_INTEGRATED_SULFATE            "none"   MISSING_R8   MISSING_R8
!     QTY_INTEGRATED_DUST               "none"   MISSING_R8   MISSING_R8
!     QTY_INTEGRATED_SMOKE              "none"   MISSING_R8   MISSING_R8
!     QTY_INTEGRATED_SEASALT            "none"   MISSING_R8   MISSING_R8
!     QTY_INTEGRATED_AOD                "none"   MISSING_R8   MISSING_R8
!     QTY_SO2                           "none"   MISSING_R8   MISSING_R8
!     QTY_SULFATE                       "none"   MISSING_R8   MISSING_R8
!     QTY_DUST                          "none"   MISSING_R8   MISSING_R8
!     QTY_SMOKE                         "none"   MISSING_R8   MISSING_R8
!     QTY_SEASALT                       "none"   MISSING_R8   MISSING_R8
! 
! ! kinds for ZVD (advanced microphysics)
!     QTY_HAIL_MIXING_RATIO             "none"   MISSING_R8   MISSING_R8
!     QTY_HAIL_NUMBER_CONCENTR          "none"   MISSING_R8   MISSING_R8
!     QTY_GRAUPEL_VOLUME                "none"   MISSING_R8   MISSING_R8
!     QTY_HAIL_VOLUME                   "none"   MISSING_R8   MISSING_R8
!     QTY_DIFFERENTIAL_REFLECTIVITY     "none"   MISSING_R8   MISSING_R8
!     QTY_SPECIFIC_DIFFERENTIAL_PHASE   "none"   MISSING_R8   MISSING_R8
!     QTY_FLASH_RATE_2D                 "none"   MISSING_R8   MISSING_R8
! 
! ! kinds for radiance 
!     QTY_RADIANCE                      "none"   MISSING_R8   MISSING_R8   ! L1 radiance (mW/cm^-1/sr/m^2)
!     QTY_BRIGHTNESS_TEMPERATURE        "none"   MISSING_R8   MISSING_R8   ! L1 brightness temperature (K)
!     QTY_BI_DIRECTIONAL_REFLECTANCE    "none"   MISSING_R8   MISSING_R8   ! L1 bi-directional reflectance (BDRF, unitless)
!     QTY_SURFACE_TYPE                  "none"   MISSING_R8   MISSING_R8   ! land = 0, sea = 1, seaice = 2
!     QTY_WIND_FETCH                    "none"   MISSING_R8   MISSING_R8   ! Wind fetch, m
!     QTY_WATER_TYPE                    "none"   MISSING_R8   MISSING_R8   ! fresh = 0, ocean = 1
!     QTY_FOAM_FRAC                     "none"   MISSING_R8   MISSING_R8   ! Fraction of foam on ocean surface (0-1)
!     QTY_INSOLUBLE_AER                 "none"   MISSING_R8   MISSING_R8   ! Insoluble aerosol OPAC aerosol (INSO)
!     QTY_H2O_SOLUBLE_AER               "none"   MISSING_R8   MISSING_R8   ! Soluble aerosol OPAC aerosol (WASO)
!     QTY_SOOT                          "none"   MISSING_R8   MISSING_R8   ! Soot aerosol OPAC aerosol (SOOT)
!     QTY_SEASALT_ACCUM                 "none"   MISSING_R8   MISSING_R8   ! Sea salt (accumulation mode) OPAC aerosol (SSAM)
!     QTY_SEASALT_COARSE                "none"   MISSING_R8   MISSING_R8   ! Sea salt (coarse) OPAC aerosol (SSCM)
!     QTY_MINERAL_NUCLEUS               "none"   MISSING_R8   MISSING_R8   ! Mineral (nucleus) OPAC aerosol (MINM)
!     QTY_MINERAL_ACCUM                 "none"   MISSING_R8   MISSING_R8   ! Mineral (accumulation mode) OPAC aerosol (MIAM)
!     QTY_MINERAL_COARSE                "none"   MISSING_R8   MISSING_R8   ! Mineral (coarse mode) OPAC aerosol (MICM)
!     QTY_MINERAL_TRANSPORTED           "none"   MISSING_R8   MISSING_R8   ! Mineral (transported mode) OPAC aerosol (MITR)
!     QTY_SULPHATED_DROPS               "none"   MISSING_R8   MISSING_R8   ! Sulphated droplets OPAC aerosol (SUSO)
!     QTY_VOLCANIC_ASH                  "none"   MISSING_R8   MISSING_R8   ! Volcanic ash OPAC aerosol (VOLA)
!     QTY_NEW_VOLCANIC_ASH              "none"   MISSING_R8   MISSING_R8   ! New volcanic ash OPAC aerosol (VAPO)
!     QTY_ASIAN_DUST                    "none"   MISSING_R8   MISSING_R8   ! Asian dust OPAC aerosol (ASDU)
!     QTY_BLACK_CARBON                  "none"   MISSING_R8   MISSING_R8   ! Black carbon CAMS aerosol (BCAR)
!     QTY_DUST_BIN1                     "none"   MISSING_R8   MISSING_R8   ! Dust bin 1 CAMS aerosol (DUS1)
!     QTY_DUST_BIN2                     "none"   MISSING_R8   MISSING_R8   ! Dust bin 2 CAMS aerosol (DUS2)
!     QTY_DUST_BIN3                     "none"   MISSING_R8   MISSING_R8   ! Dust bin 3 CAMS aerosol (DUS3)
!     QTY_AMMONIUM_SULPHATE             "none"   MISSING_R8   MISSING_R8   ! Ammonium sulphate CAMS aerosol (SULP)
!     QTY_SEA_SALT_BIN1                 "none"   MISSING_R8   MISSING_R8   ! Sea salt bin 1 CAMS aerosol (SSA1)
!     QTY_SEA_SALT_BIN2                 "none"   MISSING_R8   MISSING_R8   ! Sea salt bin 2 CAMS aerosol (SSA2)
!     QTY_SEA_SALT_BIN3                 "none"   MISSING_R8   MISSING_R8   ! Sea salt bin 3 CAMS aerosol (SSA3)
!     QTY_HYDROPHILIC_ORGANIC_MATTER    "none"   MISSING_R8   MISSING_R8   ! Hydrophilic organic matter CAMS aerosol (OMAT)
!     QTY_CLOUDWATER_DE                 "none"   MISSING_R8   MISSING_R8   ! Cloud liquid water effective diameter (microns)
!     QTY_CLOUD_ICE_DE                  "none"   MISSING_R8   MISSING_R8   ! Cloud ice effective diameter (microns)
!     QTY_COLUMN_CLOUD_FRAC             "none"   MISSING_R8   MISSING_R8   ! Simple cloud fraction (0-1)
!     QTY_CLOUD_TOP_PRESSURE            "none"   MISSING_R8   MISSING_R8
! 
! END DART PREPROCESS QUANTITY DEFINITIONS


! <next few lines under version control, do not edit>
! $URL$
! $Id$
! $Revision$
! $Date$
