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
! ! ocean related quantities
!
!   QTY_SALINITY
!   QTY_TEMPERATURE
!   QTY_POTENTIAL_TEMPERATURE
!   QTY_PRESSURE
!   QTY_VELOCITY
!   QTY_U_CURRENT_COMPONENT
!   QTY_V_CURRENT_COMPONENT
!   QTY_W_CURRENT_COMPONENT
!   QTY_DRY_LAND
!   QTY_SEA_SURFACE_PRESSURE
!   QTY_SEA_SURFACE_HEIGHT
!   QTY_SEA_SURFACE_ANOMALY
!   QTY_NITRATE_CONCENTRATION
!   QTY_PHOSPHATE_CONCENTRATION
!   QTY_DISSOLVED_OXYGEN
!   QTY_PHYTOPLANKTON_BIOMASS
!   QTY_ALKALINITY
!   QTY_DISSOLVED_INORGANIC_CARBON
!   QTY_DISSOLVED_ORGANIC_NITROGEN
!   QTY_DISSOLVED_ORGANIC_P
!   QTY_DISSOLVED_INORGANIC_IRON
!   QTY_SURFACE_CHLOROPHYLL
!   QTY_LAYER_THICKNESS
!
! ! fixme - these units are hardcoded in obs_diag and shouldn't be
!   QTY_U_WIND_COMPONENT
!   QTY_V_WIND_COMPONENT
!
! END DART PREPROCESS QUANTITY DEFINITIONS

