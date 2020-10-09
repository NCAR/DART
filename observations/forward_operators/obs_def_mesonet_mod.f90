! DART software - Copyright UCAR. This open source software is provided
! by UCAR, "as is", without charge, subject to all terms of use at
! http://www.image.ucar.edu/DAReS/DART/DART_download

! BEGIN DART PREPROCESS KIND LIST
!MESONET_U_WIND_COMPONENT,  QTY_U_WIND_COMPONENT,     COMMON_CODE
!MESONET_V_WIND_COMPONENT,  QTY_V_WIND_COMPONENT,     COMMON_CODE
!MESONET_GEOPOTENTIAL_HGT,  QTY_GEOPOTENTIAL_HEIGHT,  COMMON_CODE
!MESONET_SURFACE_PRESSURE,  QTY_SURFACE_PRESSURE,     COMMON_CODE
!MESONET_TEMPERATURE,       QTY_TEMPERATURE,          COMMON_CODE
!MESONET_SPECIFIC_HUMIDITY, QTY_SPECIFIC_HUMIDITY,    COMMON_CODE
! END DART PREPROCESS KIND LIST

! !!! Note about Specific Humidity observations:
! !!! UNITS in original BUFR are g/kg; This is converted to kg/kg by
! !!! the BUFR to obs_sequence conversion programs making it unnecessary
! !!! to multiply by 1000 at assimilation time.
! !!! PLEASE pay attention to units for specific humidity in models.

