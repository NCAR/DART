! DART software - Copyright UCAR. This open source software is provided
! by UCAR, "as is", without charge, subject to all terms of use at
! http://www.image.ucar.edu/DAReS/DART/DART_download
!


! ! in this section, define the quantity of interest.  it must
! ! start with QTY_xxx and be less than 32 characters total.
! ! if there are units, min/max values, pdf, etc you can add one or
! ! more name=value pairs.  at this point there *cannot* be spaces
! ! around the equals sign.  the value can have spaces if it is
! ! surrounded by double quotes.  e.g. desc="assumes dry air density"
! ! comment lines (like these) can be added if they start with an additional !
 
! BEGIN DART PREPROCESS QUANTITY DEFINITIONS
!
! ! ocean related quantities
!
!     QTY_SALINITY                                  minval=0.0
!     QTY_TEMPERATURE                 units="K"     minval=0.0
!     QTY_POTENTIAL_TEMPERATURE                     
!     QTY_PRESSURE                                  minval=0.0
!     QTY_VELOCITY                    units="m/s" 
!     QTY_U_CURRENT_COMPONENT         units="m/s"
!     QTY_V_CURRENT_COMPONENT         units="m/s"
!     QTY_W_CURRENT_COMPONENT         units="m/s"
!     QTY_DRY_LAND                                  minval=0.0  maxval=1.0
!     QTY_SEA_SURFACE_PRESSURE        units="mb"
!     QTY_SEA_SURFACE_HEIGHT          units="m"
!     QTY_SEA_SURFACE_ANOMALY         units="m"
! 
!
! ! fixme - these are hardcoded in obs_diag and shouldn't be
!     QTY_U_WIND_COMPONENT            units="m/s"
!     QTY_V_WIND_COMPONENT            units="m/s"
! 
! END DART PREPROCESS QUANTITY DEFINITIONS

