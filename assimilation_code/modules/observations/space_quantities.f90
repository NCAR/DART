! DART software - Copyright UCAR. This open source software is provided
! by UCAR, "as is", without charge, subject to all terms of use at
! http://www.image.ucar.edu/DAReS/DART/DART_download
!
! $Id$


 
! ! in this section you can have comments (which need a second !) or
! ! lines with a single word on them that begins QTY_    "none"   MISSING_R8   MISSING_R8
! ! can have trailing comments after:  QTY_xxx  ! comment   "none"   MISSING_R8   MISSING_R8
!
! BEGIN DART PREPROCESS QUANTITY DEFINITIONS
!
! ! kinds for TIEgcm
!     QTY_ELECTRON_DENSITY              "none"   MISSING_R8   MISSING_R8
! 
! ! kinds for generic parameters that aren't going to be
! ! directly observed but are going to be adjusted by the
! ! assimilation.
!     QTY_1D_PARAMETER                  "none"   MISSING_R8   MISSING_R8
!     QTY_2D_PARAMETER                  "none"   MISSING_R8   MISSING_R8
!     QTY_3D_PARAMETER                  "none"   MISSING_R8   MISSING_R8
! 
! ! kinds for CHAMP upper atmosphere computations
!     QTY_ATOMIC_OXYGEN_MIXING_RATIO    "none"   MISSING_R8   MISSING_R8
!     QTY_MOLEC_OXYGEN_MIXING_RATIO     "none"   MISSING_R8   MISSING_R8
! 
! ! proposed new kinds for COSMIC GPS/RO obs
! ! (currently unused)
!     QTY_OCCULTATION_REFRACTIVITY      "none"   MISSING_R8   MISSING_R8
!     QTY_OCCULTATION_EXCESSPHASE       "none"   MISSING_R8   MISSING_R8
! 
! ! kind for the other way of measuring elevation
! ! contrast this with geopotential height
!     QTY_GEOMETRIC_HEIGHT              "none"   MISSING_R8   MISSING_R8
! 
! ! kinds for satellite radiances (jason o.)
!     QTY_INFRARED_RADIANCE             "none"   MISSING_R8   MISSING_R8
!     QTY_INFRARED_BRIGHT_TEMP          "none"   MISSING_R8   MISSING_R8
!     QTY_LANDMASK                      "none"   MISSING_R8   MISSING_R8
! 
! ! kinds for planetary remote sensing (wglawson, c.lee)
!     QTY_SKIN_TEMPERATURE              "none"   MISSING_R8   MISSING_R8
!     QTY_NADIR_RADIANCE                "none"   MISSING_R8   MISSING_R8
!     QTY_TRACER_1_MIXING_RATIO         "none"   MISSING_R8   MISSING_R8
!     QTY_TRACER_2_MIXING_RATIO         "none"   MISSING_R8   MISSING_R8
!     ! Is QTY_TRACER_MIXING_RATIO necessary with QTY_TRACER_CONCENTRATION   "none"   MISSING_R8   MISSING_R8
!     !   (= 29) available from the simple advection model?
!     QTY_SOIL_TEMPERATURE              "none"   MISSING_R8   MISSING_R8
!     QTY_SOIL_LIQUID_WATER             "none"   MISSING_R8   MISSING_R8
! 
! ! more kinds for planetary remote sensing (c.lee)
!     QTY_SURFACE_ALBEDO                "none"   MISSING_R8   MISSING_R8
!     QTY_SURFACE_EMISSIVITY            "none"   MISSING_R8   MISSING_R8
!     QTY_DUST_OPACITY_7MB              "none"   MISSING_R8   MISSING_R8
!     QTY_THC                           "none"   MISSING_R8   MISSING_R8
! 
! ! kinds for COAMPS (Tim Whitcomb)
!     QTY_EXNER_FUNCTION                "none"   MISSING_R8   MISSING_R8
!     QTY_TURBULENT_KINETIC_ENERGY      "none"   MISSING_R8   MISSING_R8
!     QTY_TOTAL_PRECIPITABLE_WATER      "none"   MISSING_R8   MISSING_R8
!     QTY_VERTLEVEL                     "none"   MISSING_R8   MISSING_R8
!     QTY_MICROWAVE_BRIGHT_TEMP         "none"   MISSING_R8   MISSING_R8
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
! 
! ! more kinds for TIEGCM Alex Chartier 
!     QTY_VERTICAL_TEC                  "none"   MISSING_R8   MISSING_R8
!     QTY_O_N2_COLUMN_DENSITY_RATIO     "none"   MISSING_R8   MISSING_R8
! 
! 
! ! kinds for GITM (Alexey Morozov)
!   QTY_TEMPERATURE_ELECTRON            "none"   MISSING_R8   MISSING_R8
!   QTY_TEMPERATURE_ION                 "none"   MISSING_R8   MISSING_R8
!   QTY_DENSITY_NEUTRAL_O3P             "none"   MISSING_R8   MISSING_R8
!   QTY_DENSITY_NEUTRAL_O2              "none"   MISSING_R8   MISSING_R8
!   QTY_DENSITY_NEUTRAL_N2              "none"   MISSING_R8   MISSING_R8
!   QTY_DENSITY_NEUTRAL_N4S             "none"   MISSING_R8   MISSING_R8
!   QTY_DENSITY_NEUTRAL_NO              "none"   MISSING_R8   MISSING_R8
!   QTY_DENSITY_NEUTRAL_N2D             "none"   MISSING_R8   MISSING_R8
!   QTY_DENSITY_NEUTRAL_N2P             "none"   MISSING_R8   MISSING_R8
!   QTY_DENSITY_NEUTRAL_H               "none"   MISSING_R8   MISSING_R8
!   QTY_DENSITY_NEUTRAL_HE              "none"   MISSING_R8   MISSING_R8
!   QTY_DENSITY_NEUTRAL_CO2             "none"   MISSING_R8   MISSING_R8
!   QTY_DENSITY_NEUTRAL_O1D             "none"   MISSING_R8   MISSING_R8
!   QTY_DENSITY_ION_O4SP                "none"   MISSING_R8   MISSING_R8
!   QTY_DENSITY_ION_O2P                 "none"   MISSING_R8   MISSING_R8
!   QTY_DENSITY_ION_N2P                 "none"   MISSING_R8   MISSING_R8
!   QTY_DENSITY_ION_NP                  "none"   MISSING_R8   MISSING_R8
!   QTY_DENSITY_ION_NOP                 "none"   MISSING_R8   MISSING_R8
!   QTY_DENSITY_ION_O2DP                "none"   MISSING_R8   MISSING_R8
!   QTY_DENSITY_ION_O2PP                "none"   MISSING_R8   MISSING_R8
!   QTY_DENSITY_ION_HP                  "none"   MISSING_R8   MISSING_R8
!   QTY_DENSITY_ION_HEP                 "none"   MISSING_R8   MISSING_R8
!   QTY_DENSITY_ION_E                   "none"   MISSING_R8   MISSING_R8
!   QTY_VELOCITY_U                      "none"   MISSING_R8   MISSING_R8
!   QTY_VELOCITY_V                      "none"   MISSING_R8   MISSING_R8
!   QTY_VELOCITY_W                      "none"   MISSING_R8   MISSING_R8
!   QTY_VELOCITY_U_ION                  "none"   MISSING_R8   MISSING_R8
!   QTY_VELOCITY_V_ION                  "none"   MISSING_R8   MISSING_R8
!   QTY_VELOCITY_W_ION                  "none"   MISSING_R8   MISSING_R8
!   QTY_VELOCITY_VERTICAL_O3P           "none"   MISSING_R8   MISSING_R8
!   QTY_VELOCITY_VERTICAL_O2            "none"   MISSING_R8   MISSING_R8
!   QTY_VELOCITY_VERTICAL_N2            "none"   MISSING_R8   MISSING_R8
!   QTY_VELOCITY_VERTICAL_N4S           "none"   MISSING_R8   MISSING_R8
!   QTY_VELOCITY_VERTICAL_NO            "none"   MISSING_R8   MISSING_R8
!   QTY_GND_GPS_VTEC                    "none"   MISSING_R8   MISSING_R8
!   QTY_DENSITY_ION_OP                  "none"   MISSING_R8   MISSING_R8
!  
! 
! END DART PREPROCESS QUANTITY DEFINITIONS


! <next few lines under version control, do not edit>
! $URL$
! $Id$
! $Revision$
! $Date$
