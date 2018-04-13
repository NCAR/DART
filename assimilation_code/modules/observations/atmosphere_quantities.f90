! DART software - Copyright UCAR. This open source software is provided
! by UCAR, "as is", without charge, subject to all terms of use at
! http://www.image.ucar.edu/DAReS/DART/DART_download
!
! $Id$


 
! ! in this section you can have comments (which need a second !) or
! ! lines with a single word on them that begins QTY_ 
! ! can have trailing comments after QTY_xxx  word
!
! BEGIN DART PREPROCESS QUANTITY DEFINITIONS
!
!     QTY_STATE_VARIABLE             
!     QTY_U_WIND_COMPONENT             ! can have trailing comments
!     QTY_V_WIND_COMPONENT! can have no whitespace           
!     QTY_SURFACE_PRESSURE           
!     QTY_TEMPERATURE                
!     QTY_SPECIFIC_HUMIDITY          
!     QTY_PRESSURE                   
!     QTY_VERTICAL_VELOCITY          
!     QTY_RAINWATER_MIXING_RATIO     
!     QTY_DEWPOINT                   
!     QTY_DENSITY                    
!     QTY_VELOCITY                   
!     QTY_RADAR_REFLECTIVITY         
!     QTY_GRAUPEL_MIXING_RATIO       
!     QTY_SNOW_MIXING_RATIO          
!     QTY_GPSRO                      
!     QTY_CLOUD_LIQUID_WATER         
!     QTY_CLOUD_ICE                  
!     QTY_CONDENSATIONAL_HEATING     
!     QTY_VAPOR_MIXING_RATIO         
!     QTY_ICE_NUMBER_CONCENTRATION   
!     QTY_GEOPOTENTIAL_HEIGHT        
!     QTY_POTENTIAL_TEMPERATURE      
!     QTY_SOIL_MOISTURE              
!     QTY_SURFACE_ELEVATION          
! 
! ! kinds for Gravity Wave Drag (CAM - kevin)
!     QTY_GRAV_WAVE_DRAG_EFFIC       
!     QTY_GRAV_WAVE_STRESS_FRACTION  
! 
!     QTY_POWER_WEIGHTED_FALL_SPEED  
!     QTY_CLOUDWATER_MIXING_RATIO    
!     QTY_ICE_MIXING_RATIO           
!     QTY_DROPLET_NUMBER_CONCENTR    
!     QTY_SNOW_NUMBER_CONCENTR       
!     QTY_RAIN_NUMBER_CONCENTR       
!     QTY_GRAUPEL_NUMBER_CONCENTR    
!     QTY_CLOUD_FRACTION             
!     QTY_ICE_FRACTION               
!     QTY_RELATIVE_HUMIDITY          
! 
! ! kinds for generic parameters that aren't going to be
! ! directly observed but are going to be adjusted by the
! ! assimilation.
!     QTY_1D_PARAMETER               
!     QTY_2D_PARAMETER               
!     QTY_3D_PARAMETER               
! 
! ! kinds for tendencies
!     QTY_ALTIMETER_TENDENCY         
! 
! ! kind for precip water; contrast with
! ! total precip water (also in this file), 
! ! which is the total column integrated value. 
!     QTY_PRECIPITABLE_WATER         
! 
! ! proposed new kinds for COSMIC GPS/RO obs
! ! (currently unused)
!     QTY_OCCULTATION_REFRACTIVITY   
!     QTY_OCCULTATION_EXCESSPHASE    
! 
! ! kind for the other way of measuring elevation
! ! contrast this with geopotential height
!     QTY_GEOMETRIC_HEIGHT           
! 
!     QTY_INFRARED_RADIANCE          
!     QTY_INFRARED_BRIGHT_TEMP       
!     QTY_LANDMASK                   
! 
! ! kind for unstructured grids
!     QTY_EDGE_NORMAL_SPEED          
! 
! ! kind for cloud liquid water path
! ! see also kinds 363, 364
!     QTY_CLW_PATH                   
! 
! ! kind for wind power
!     QTY_WIND_TURBINE_POWER         
! 
! ! kinds for surface fields
!     QTY_2M_SPECIFIC_HUMIDITY       
!     QTY_2M_TEMPERATURE             
!     QTY_10M_U_WIND_COMPONENT       
!     QTY_10M_V_WIND_COMPONENT       
! 
! ! kinds for planetary remote sensing (wglawson, c.lee)
!     QTY_SKIN_TEMPERATURE           
!     QTY_NADIR_RADIANCE             
!     QTY_TRACER_1_MIXING_RATIO      
!     QTY_TRACER_2_MIXING_RATIO      
!     ! Is QTY_TRACER_MIXING_RATIO necessary with QTY_TRACER_CONCENTRATION
!     !   (= 29) available from the simple advection model?
!     QTY_SOIL_TEMPERATURE           
!     QTY_SOIL_LIQUID_WATER          
! 
! ! kinds for NCOMMAS  (Lou W., Ted M.)
!     QTY_VERTICAL_VORTICITY         
! 
! ! more kinds for planetary remote sensing (c.lee)
!     QTY_SURFACE_ALBEDO             
!     QTY_SURFACE_EMISSIVITY         
!     QTY_DUST_OPACITY_7MB           
!     QTY_THC                        
! 
! ! kinds for vortex tracking (WRF - yongsheng)
!     QTY_VORTEX_LON                 
!     QTY_VORTEX_LAT                 
!     QTY_VORTEX_PMIN                
!     QTY_VORTEX_WMAX                
! 
! ! kinds for COAMPS (Tim Whitcomb)
!     QTY_EXNER_FUNCTION             
!     QTY_TURBULENT_KINETIC_ENERGY   
!     QTY_TOTAL_PRECIPITABLE_WATER   
!     QTY_VERTLEVEL                  
!     QTY_MICROWAVE_BRIGHT_TEMP      
! 
! ! kinds for NAAPS (Walter R. Sessions)
!     QTY_INTEGRATED_SULFATE         
!     QTY_INTEGRATED_DUST            
!     QTY_INTEGRATED_SMOKE           
!     QTY_INTEGRATED_SEASALT         
!     QTY_INTEGRATED_AOD             
!     QTY_SO2                        
!     QTY_SULFATE                    
!     QTY_DUST                       
!     QTY_SMOKE                      
!     QTY_SEASALT                    
! 
! ! kinds for ZVD (advanced microphysics)
!     QTY_HAIL_MIXING_RATIO          
!     QTY_HAIL_NUMBER_CONCENTR       
!     QTY_GRAUPEL_VOLUME             
!     QTY_HAIL_VOLUME                
!     QTY_DIFFERENTIAL_REFLECTIVITY  
!     QTY_SPECIFIC_DIFFERENTIAL_PHASE
!     QTY_FLASH_RATE_2D              
! 
! 
!   QTY_CWP_PATH                     
!   QTY_CWP_PATH_ZERO                
! 
! END DART PREPROCESS QUANTITY DEFINITIONS


! <next few lines under version control, do not edit>
! $URL$
! $Id$
! $Revision$
! $Date$
