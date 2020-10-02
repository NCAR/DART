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
! ! trailing comments can be added if they start with !
!
! BEGIN DART PREPROCESS QUANTITY DEFINITIONS
!
! ! kinds for CLM - Community Land Model (Tim Hoar)
!     QTY_SNOW_THICKNESS                units="none"
!     QTY_SNOW_WATER                    units="none"
!     QTY_SNOWCOVER_FRAC                units="none"
!     QTY_LIQUID_WATER                  units="none"
!     QTY_ICE                           units="none"
!     QTY_CARBON                        units="none"
!     QTY_SOIL_CARBON                   units="none"
!     QTY_ROOT_CARBON                   units="none"
!     QTY_STEM_CARBON                   units="none"
!     QTY_LEAF_CARBON                   units="none"
!     QTY_LEAF_AREA_INDEX               units="none"
!     QTY_NET_CARBON_FLUX               units="none"
!     QTY_LATENT_HEAT_FLUX              units="none"
!     QTY_SENSIBLE_HEAT_FLUX            units="none"
!     QTY_RADIATION                     units="none"
!     QTY_NET_CARBON_PRODUCTION         units="none"
!     QTY_NITROGEN                      units="none"
!     QTY_SOIL_NITROGEN                 units="none"
!     QTY_ROOT_NITROGEN                 units="none"
!     QTY_STEM_NITROGEN                 units="none"
!     QTY_LEAF_NITROGEN                 units="none"
!     QTY_WATER_TABLE_DEPTH             units="none"
!     QTY_FPAR                          units="none"
!     QTY_TOTAL_WATER_STORAGE           units="none"
! 
! ! kinds for NOAH  (Tim Hoar)
!     QTY_NEUTRON_INTENSITY             units="none"
!     QTY_CANOPY_WATER                  units="none"
!     QTY_GROUND_HEAT_FLUX              units="none"
!  
! ! more land kinds
!   QTY_BRIGHTNESS_TEMPERATURE          units="none"
!   QTY_VEGETATION_TEMPERATURE          units="none"
!   QTY_CANOPY_HEIGHT                   units="none"
!   QTY_FPAR_DIRECT                     units="none"
!   QTY_FPAR_DIFFUSE                    units="none"
!   QTY_FPAR_SUNLIT_DIRECT              units="none"
!   QTY_FPAR_SUNLIT_DIFFUSE             units="none"
!   QTY_FPAR_SHADED_DIRECT              units="none"
!   QTY_FPAR_SHADED_DIFFUSE             units="none"
!   QTY_FPSN                            units="none"
!   QTY_FSIF                            units="none"
! 
! 
! END DART PREPROCESS QUANTITY DEFINITIONS

