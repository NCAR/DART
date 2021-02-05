! DART software - Copyright UCAR. This open source software is provided
! by UCAR, "as is", without charge, subject to all terms of use at
! http://www.image.ucar.edu/DAReS/DART/DART_download

! ! in this section, define the quantity of interest.  it must
! ! start with QTY_xxx and be less than 32 characters total.
! ! if there are units, min/max values, pdf, etc you can add one or
! ! more name=value pairs.  at this point there *cannot* be spaces
! ! around the equals sign.  the value can have spaces if it is
! ! surrounded by double quotes.  e.g. desc="assumes dry air density"
! ! comment lines (like these) can be added if they start with an additional !

! BEGIN DART PREPROCESS QUANTITY DEFINITIONS
!
! ! kinds for TIE-GCM
!   QTY_ELECTRON_DENSITY
!
! ! kinds for generic parameters that aren't going to be
! ! directly observed but are going to be adjusted by the
! ! assimilation.
!   QTY_1D_PARAMETER
!   QTY_2D_PARAMETER
!   QTY_3D_PARAMETER
!
! ! kinds for CHAMP upper atmosphere computations
!   QTY_ATOMIC_OXYGEN_MIXING_RATIO
!   QTY_MOLEC_OXYGEN_MIXING_RATIO
!
! ! proposed new kinds for COSMIC GPS/RO obs
! ! (currently unused)
!   QTY_OCCULTATION_REFRACTIVITY
!   QTY_OCCULTATION_EXCESSPHASE
!
! ! kind for the other way of measuring elevation
! ! contrast this with geopotential height
!   QTY_GEOMETRIC_HEIGHT
!
! ! kinds for satellite radiances (jason o.)
!   QTY_INFRARED_RADIANCE
!   QTY_INFRARED_BRIGHT_TEMP
!   QTY_LANDMASK
!
! ! kinds for planetary remote sensing (wglawson, c.lee)
!   QTY_SKIN_TEMPERATURE
!   QTY_NADIR_RADIANCE
!   QTY_TRACER_1_MIXING_RATIO
!   QTY_TRACER_2_MIXING_RATIO
!   QTY_SOIL_TEMPERATURE
!   QTY_SOIL_LIQUID_WATER
!
! ! more kinds for planetary remote sensing (c.lee)
!   QTY_SURFACE_ALBEDO
!   QTY_SURFACE_EMISSIVITY
!   QTY_DUST_OPACITY_7MB
!   QTY_THC
!
! ! kinds for COAMPS (Tim Whitcomb)
!   QTY_EXNER_FUNCTION
!   QTY_TURBULENT_KINETIC_ENERGY
!   QTY_TOTAL_PRECIPITABLE_WATER
!   QTY_VERTLEVEL
!   QTY_MICROWAVE_BRIGHT_TEMP
!
! ! kinds for ZVD (advanced microphysics)
!   QTY_HAIL_MIXING_RATIO
!   QTY_HAIL_NUMBER_CONCENTR
!   QTY_GRAUPEL_VOLUME
!   QTY_HAIL_VOLUME
!   QTY_DIFFERENTIAL_REFLECTIVITY
!   QTY_SPECIFIC_DIFFERENTIAL_PHASE
!   QTY_FLASH_RATE_2D
!
! ! more kinds for TIEGCM Alex Chartier
!   QTY_VERTICAL_TEC
!   QTY_O_N2_COLUMN_DENSITY_RATIO
!
! ! kinds for GITM (Alexey Morozov)
!   QTY_TEMPERATURE_ELECTRON
!   QTY_TEMPERATURE_ION
!   QTY_DENSITY_NEUTRAL_O3P
!   QTY_DENSITY_NEUTRAL_O2
!   QTY_DENSITY_NEUTRAL_N2
!   QTY_DENSITY_NEUTRAL_N4S
!   QTY_DENSITY_NEUTRAL_NO
!   QTY_DENSITY_NEUTRAL_N2D
!   QTY_DENSITY_NEUTRAL_N2P
!   QTY_DENSITY_NEUTRAL_H
!   QTY_DENSITY_NEUTRAL_HE
!   QTY_DENSITY_NEUTRAL_CO2
!   QTY_DENSITY_NEUTRAL_O1D
!   QTY_DENSITY_ION_O4SP
!   QTY_DENSITY_ION_O2P
!   QTY_DENSITY_ION_N2P
!   QTY_DENSITY_ION_NP
!   QTY_DENSITY_ION_NOP
!   QTY_DENSITY_ION_O2DP
!   QTY_DENSITY_ION_O2PP
!   QTY_DENSITY_ION_OP
!   QTY_DENSITY_ION_HP
!   QTY_DENSITY_ION_HEP
!   QTY_DENSITY_ION_E
!   QTY_VELOCITY_U
!   QTY_VELOCITY_V
!   QTY_VELOCITY_W
!   QTY_VELOCITY_U_ION
!   QTY_VELOCITY_V_ION
!   QTY_VELOCITY_W_ION
!   QTY_VELOCITY_VERTICAL_O3P
!   QTY_VELOCITY_VERTICAL_O2
!   QTY_VELOCITY_VERTICAL_N2
!   QTY_VELOCITY_VERTICAL_N4S
!   QTY_VELOCITY_VERTICAL_NO
!   QTY_GND_GPS_VTEC
!
! END DART PREPROCESS QUANTITY DEFINITIONS

