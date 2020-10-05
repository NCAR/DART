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
! QTY_ABSORBED_PAR
! QTY_BIOMASS
! QTY_BRIGHTNESS_TEMPERATURE        units="none"
! QTY_CANOPY_HEIGHT                 units="none"
! QTY_CANOPY_TEMPERATURE
! QTY_CANOPY_WATER                  units="none"
! QTY_CARBON                        units="none"
! QTY_CLAY_FRACTION
! QTY_DEAD_ROOT_CARBON
! QTY_DEAD_ROOT_NITROGEN
! QTY_DEAD_STEM_CARBON
! QTY_DEAD_STEM_NITROGEN
! QTY_ER_FLUX
! QTY_FPAR                          units="none"
! QTY_FPAR_DIFFUSE                  units="none"
! QTY_FPAR_DIRECT                   units="none"
! QTY_FPAR_SHADED_DIFFUSE           units="none"
! QTY_FPAR_SHADED_DIRECT            units="none"
! QTY_FPAR_SUNLIT_DIFFUSE           units="none"
! QTY_FPAR_SUNLIT_DIRECT            units="none"
! QTY_FPSN                          units="none"
! QTY_FRACTION_ABSORBED_PAR
! QTY_FRAC_PHOTO_AVAIL_RADIATION
! QTY_FSIF                          units="none"
! QTY_GEOPOTENTIAL_HEIGHT           units="m"         
! QTY_GROSS_PRIMARY_PROD_FLUX
! QTY_GROUND_HEAT_FLUX              units="none"
! QTY_ICE                           units="none"
! QTY_LATENT_HEAT_FLUX              units="none"
! QTY_LEAF_AREA_INDEX               units="none"
! QTY_LEAF_CARBON                   units="none"
! QTY_LEAF_NITROGEN                 units="none"
! QTY_LIQUID_WATER                  units="none"
! QTY_LIVE_ROOT_CARBON
! QTY_LIVE_ROOT_NITROGEN
! QTY_LIVE_STEM_CARBON
! QTY_LIVE_STEM_NITROGEN
! QTY_NET_CARBON_FLUX               units="none"
! QTY_NET_CARBON_PRODUCTION         units="none"
! QTY_NET_PRIMARY_PROD_FLUX
! QTY_NEUTRON_INTENSITY             units="none"
! QTY_NITROGEN                      units="none"
! QTY_PAR_DIFFUSE
! QTY_PAR_DIRECT
! QTY_PHOTO_AVAILABLE_RADIATION
! QTY_RADIATION                     units="none"
! QTY_RADIATION_NEAR_IR_DOWN
! QTY_RADIATION_NEAR_IR_UP
! QTY_RADIATION_VISIBLE_DOWN
! QTY_RADIATION_VISIBLE_UP
! QTY_ROOT_CARBON                   units="none"
! QTY_ROOT_NITROGEN                 units="none"
! QTY_RTM_PARAMETERS_N
! QTY_RTM_PARAMETERS_P
! QTY_RUNOFF_MULTIPLIER
! QTY_SAND_FRACTION
! QTY_SENSIBLE_HEAT_FLUX            units="none"
! QTY_SNOWCOVER_FRAC                units="none"
! QTY_SNOW_DEPTH
! QTY_SNOW_GRAIN_SIZE
! QTY_SNOW_MIXING_RATIO
! QTY_SNOW_NUMBER_CONCENTR
! QTY_SNOW_TEMPERATURE
! QTY_SNOW_THICKNESS                units="none"
! QTY_SNOW_WATER                    units="none"
! QTY_SOIL_CARBON                   units="none"
! QTY_SOIL_ICE
! QTY_SOIL_LIQUID_WATER
! QTY_SOIL_MOISTURE
! QTY_SOIL_NITROGEN                 units="none"
! QTY_SOIL_RESPIRATION_FLUX
! QTY_SOIL_TEMPERATURE              units="K"     desc="is soil temp really in K or C?"
! QTY_SOLAR_INDUCED_FLUORESCENCE
! QTY_STEM_AREA_INDEX
! QTY_STEM_CARBON                   units="none"
! QTY_STEM_NITROGEN                 units="none"
! QTY_STREAM_FLOW
! QTY_STREAM_HEIGHT
! QTY_SURFACE_ALBEDO
! QTY_SURFACE_EMISSIVITY
! QTY_SURFACE_HEAD
! QTY_SURFACE_PRESSURE              units="hPa"    minval=0.0          
! QTY_SURFACE_RUNOFF
! QTY_SURFACE_TYPE
! QTY_TEMPERATURE                   units="K"      minval=0.0          
! QTY_TOTAL_WATER_STORAGE           units="none"
! QTY_UNCONFINED_WATER
! QTY_UNDER_RUNOFF
! QTY_U_WIND_COMPONENT              units="m/s"       
! QTY_VEGETATED_AREA_FRACTION
! QTY_VEGETATION_TEMPERATURE        units="none"
! QTY_V_WIND_COMPONENT              units="m/s"       
! QTY_WATER_TABLE_DEPTH             units="none"
! 
! END DART PREPROCESS QUANTITY DEFINITIONS

