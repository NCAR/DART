! Data Assimilation Research Testbed -- DART
! Copyright 2004-2007, Data Assimilation Research Section
! University Corporation for Atmospheric Research
! Licensed under the GPL -- www.gpl.org/licenses/gpl.html

module model_mod

! <next few lines under version control, do not edit>
! $URL$
! $Id$
! $Revision$
! $Date$

! Assimilation interface for WRF model

!-----------------------------------------------------------------------
!
!     interface for WRF
!
!-----------------------------------------------------------------------
!---------------- m o d u l e   i n f o r m a t i o n ------------------
!-----------------------------------------------------------------------

use         types_mod, only : r8, deg2rad, missing_r8, ps0, earth_radius, &
                              gas_constant, gas_constant_v, gravity

use  time_manager_mod, only : time_type, set_time, set_calendar_type, GREGORIAN

use      location_mod, only : location_type, get_location, set_location, &
                              horiz_dist_only, &
                              LocationDims, LocationName, LocationLName, &
                              query_location, vert_is_undef, vert_is_surface, &
                              vert_is_level, vert_is_pressure, vert_is_height, &
                              VERTISUNDEF, VERTISSURFACE, VERTISLEVEL, VERTISPRESSURE, &
                              VERTISHEIGHT,&
                              get_close_type, get_dist, get_close_maxdist_init, &
                              get_close_obs_init, loc_get_close_obs => get_close_obs

use     utilities_mod, only : file_exist, open_file, close_file, &
                              register_module, error_handler, E_ERR, &
                              E_MSG, logfileunit, find_namelist_in_file, check_namelist_read

use      obs_kind_mod, only : KIND_U_WIND_COMPONENT, KIND_V_WIND_COMPONENT, &
                              KIND_SURFACE_PRESSURE, KIND_TEMPERATURE, &
                              KIND_SPECIFIC_HUMIDITY, &
                              KIND_PRESSURE, KIND_VERTICAL_VELOCITY, &
                              KIND_RAINWATER_MIXING_RATIO, KIND_DENSITY, &
                              KIND_GRAUPEL_MIXING_RATIO, KIND_SNOW_MIXING_RATIO, &
                              KIND_CLOUD_LIQUID_WATER, KIND_CLOUD_ICE, &
                              KIND_CONDENSATIONAL_HEATING, KIND_VAPOR_MIXING_RATIO, &
                              KIND_ICE_NUMBER_CONCENTRATION, KIND_GEOPOTENTIAL_HEIGHT, &
                              KIND_POTENTIAL_TEMPERATURE, KIND_SOIL_MOISTURE, &
                              KIND_VORTEX_LAT, KIND_VORTEX_LON, &
                              KIND_VORTEX_PMIN, KIND_VORTEX_WMAX

use         map_utils, only : proj_info, map_init, map_set, latlon_to_ij, &
                              PROJ_LATLON, PROJ_MERC, PROJ_LC, PROJ_PS, &
                              ij_to_latlon, gridwind_to_truewind
use netcdf
use typesizes

implicit none
private

!  public routines and data for the WRF model

public ::  get_model_size,                    &
           get_state_meta_data,               &
           model_interpolate,                 &
           vert_interpolate,                  &
           get_model_time_step,               &
           static_init_model,                 &
           pert_model_state,                  &
           nc_write_model_atts,               &
           nc_write_model_vars,               &
           get_close_obs,                     &
           ens_mean_for_model,                &
           get_close_maxdist_init,            &
           get_close_obs_init


!  public stubs 

public ::  adv_1step,       &
           end_model,       &
           init_time,       &
           init_conditions

!-----------------------------------------------------------------------
! version controlled file description for error handling, do not edit
character(len=128), parameter :: &
   source   = "$URL$", &
   revision = "$Revision$", &
   revdate  = "$Date$"

!-----------------------------------------------------------------------
! Model namelist parameters with default values.
!
! center_search_half_length:  half length (in meter) of the searching box to locate 
!                             minimum pressure at a grid point
! center_spline_scale: coarse grid to spline interp. fine grid ratio
!-----------------------------------------------------------------------

logical :: output_state_vector  = .false.     ! output prognostic variables
integer :: num_moist_vars       = 3
integer :: num_domains          = 1
integer :: calendar_type        = GREGORIAN
integer :: assimilation_period_seconds = 21600
logical :: surf_obs             = .true.
logical :: soil_data            = .true.
logical :: h_diab               = .false.
character(len = 72) :: adv_mod_command = './wrf.exe'
real (kind=r8) :: center_search_half_length = 500000.0_r8
integer :: center_search_half_size
integer :: center_spline_grid_scale = 10
integer :: vert_localization_coord = VERTISHEIGHT

real(r8), allocatable :: ens_mean(:)

namelist /model_nml/ output_state_vector, num_moist_vars, &
                     num_domains, calendar_type, surf_obs, soil_data, h_diab, &
                     adv_mod_command, assimilation_period_seconds, &
                     vert_localization_coord, &
                     center_search_half_length, center_spline_grid_scale

!-----------------------------------------------------------------------

! Private definition of domain map projection use by WRF

integer, parameter :: map_sphere = 0, map_lambert = 1, map_polar_stereo = 2, map_mercator = 3

! Private definition of model variable types

integer, parameter :: TYPE_U     = 1,   TYPE_V     = 2,   TYPE_W     = 3,  &
                      TYPE_GZ    = 4,   TYPE_T     = 5,   TYPE_MU    = 6,  &
                      TYPE_TSK   = 7,   TYPE_QV    = 8,   TYPE_QC    = 9,  &
                      TYPE_QR    = 10,  TYPE_QI    = 11,  TYPE_QS    = 12, &
                      TYPE_QG    = 13,  TYPE_QNICE = 14,  TYPE_U10   = 15, &
                      TYPE_V10   = 16,  TYPE_T2    = 17,  TYPE_TH2   = 18, &
                      TYPE_Q2    = 19,  TYPE_PS    = 20,  TYPE_TSLB  = 21, &
                      TYPE_SMOIS = 22,  TYPE_SH2O  = 23,  TYPE_HDIAB = 24
integer, parameter :: num_model_var_types = 24

real (kind=r8), PARAMETER    :: kappa = 2.0_r8/7.0_r8 ! gas_constant / cp
real (kind=r8), PARAMETER    :: ts0 = 300.0_r8        ! Base potential temperature for all levels.

!---- private data ----

logical :: surf_var             = .false.

TYPE wrf_static_data_for_dart

   integer  :: bt, bts, sn, sns, we, wes, sls
   real(r8) :: dx, dy, dt, p_top
   integer  :: map_proj
   real(r8) :: cen_lat,cen_lon
   type(proj_info) :: proj

   integer  :: n_moist
   logical  :: surf_obs
   logical  :: soil_data
   integer  :: vert_coord
   real(r8), dimension(:),     pointer :: znu, dn, dnw, zs
   real(r8), dimension(:,:),   pointer :: mub, latitude, longitude, hgt
   real(r8), dimension(:,:),   pointer :: mapfac_m, mapfac_u, mapfac_v
   real(r8), dimension(:,:,:), pointer :: phb

   integer :: number_of_wrf_variables
   integer, dimension(:,:), pointer :: var_index
   integer, dimension(:,:), pointer :: var_size
   integer, dimension(:),   pointer :: var_type
   integer, dimension(:),   pointer :: dart_kind
   integer, dimension(:,:), pointer :: land

   integer, dimension(:,:,:,:), pointer :: dart_ind

end type wrf_static_data_for_dart

type wrf_dom
   type(wrf_static_data_for_dart), pointer :: dom(:)
   integer :: model_size
end type wrf_dom

type(wrf_dom) :: wrf


contains

!#######################################################################

subroutine static_init_model()

! INitializes class data for WRF???

integer :: ncid
integer :: io, iunit

character (len=80)    :: name
character (len=1)     :: idom
logical, parameter    :: debug = .false.
integer               :: var_id, ind, i, j, k, id, dart_index, model_type

integer  :: proj_code
real(r8) :: stdlon,truelat1,truelat2,dt
character(len=129) :: errstring

!----------------------------------------------------------------------

! Register the module
call register_module(source, revision, revdate)

! Begin by reading the namelist input
call find_namelist_in_file("input.nml", "model_nml", iunit)
read(iunit, nml = model_nml, iostat = io)
call check_namelist_read(iunit, io, "model_nml")

! Record the namelist values used for the run ...
call error_handler(E_MSG,'static_init_model','model_nml values are',' ',' ',' ')
write(logfileunit, nml=model_nml)
write(     *     , nml=model_nml)

allocate(wrf%dom(num_domains))

wrf%dom(:)%n_moist = num_moist_vars

if( num_moist_vars > 7) then
   write(*,'(''num_moist_vars = '',i3)')num_moist_vars
   call error_handler(E_ERR,'static_init_model', &
        'num_moist_vars is too large.', source, revision,revdate)
endif

wrf%dom(:)%surf_obs = surf_obs
wrf%dom(:)%soil_data = soil_data

if ( debug ) then
   if ( output_state_vector ) then
      write(*,*)'netcdf file in state vector format'
   else
      write(*,*)'netcdf file in prognostic vector format'
   endif
endif

call set_calendar_type(calendar_type)

! Store vertical localization coordinate
! Only 3 are allowed: level(1), pressure(2), or height(3)
! Everything else is assumed height
if (vert_localization_coord == VERTISLEVEL) then
   wrf%dom(:)%vert_coord = VERTISLEVEL
elseif (vert_localization_coord == VERTISPRESSURE) then
   wrf%dom(:)%vert_coord = VERTISPRESSURE
elseif (vert_localization_coord == VERTISHEIGHT) then
   wrf%dom(:)%vert_coord = VERTISHEIGHT
else
   write(errstring,*)'vert_localization_coord must be one of ', &
                     VERTISLEVEL, VERTISPRESSURE, VERTISHEIGHT
   call error_handler(E_MSG,'static_init_model', errstring, source, revision,revdate)
   write(errstring,*)'vert_localization_coord is ', vert_localization_coord
   call error_handler(E_ERR,'static_init_model', errstring, source, revision,revdate)
endif

dart_index = 1

call read_dt_from_wrf_nml()

do id=1,num_domains

   write( idom , '(I1)') id

   write(*,*) '******************'
   write(*,*) '**  DOMAIN # ',idom,'  **'
   write(*,*) '******************'

   if(file_exist('wrfinput_d0'//idom)) then

      call check( nf90_open('wrfinput_d0'//idom, NF90_NOWRITE, ncid) )

   else

      call error_handler(E_ERR,'static_init_model', &
           'Please put wrfinput_d0'//idom//' in the work directory.', source, revision,revdate)

   endif

   if(debug) write(*,*) ' ncid is ',ncid

! get wrf grid dimensions

   call check( nf90_inq_dimid(ncid, "bottom_top", var_id) )
   call check( nf90_inquire_dimension(ncid, var_id, name, wrf%dom(id)%bt) )

   call check( nf90_inq_dimid(ncid, "bottom_top_stag", var_id) ) ! reuse var_id, no harm
   call check( nf90_inquire_dimension(ncid, var_id, name, wrf%dom(id)%bts) )

   call check( nf90_inq_dimid(ncid, "south_north", var_id) )
   call check( nf90_inquire_dimension(ncid, var_id, name, wrf%dom(id)%sn) )

   call check( nf90_inq_dimid(ncid, "south_north_stag", var_id)) ! reuse var_id, no harm
   call check( nf90_inquire_dimension(ncid, var_id, name, wrf%dom(id)%sns) )

   call check( nf90_inq_dimid(ncid, "west_east", var_id) )
   call check( nf90_inquire_dimension(ncid, var_id, name, wrf%dom(id)%we) )

   call check( nf90_inq_dimid(ncid, "west_east_stag", var_id) )  ! reuse var_id, no harm
   call check( nf90_inquire_dimension(ncid, var_id, name, wrf%dom(id)%wes) )

   call check( nf90_inq_dimid(ncid, "soil_layers_stag", var_id) )  ! reuse var_id, no harm
   call check( nf90_inquire_dimension(ncid, var_id, name, wrf%dom(id)%sls) )

   if(debug) then
      write(*,*) ' dimensions bt, sn, we are ',wrf%dom(id)%bt, &
           wrf%dom(id)%sn, wrf%dom(id)%we
      write(*,*) ' staggered  bt, sn, we are ',wrf%dom(id)%bts, &
           wrf%dom(id)%sns,wrf%dom(id)%wes
   endif

! get meta data and static data we need

   call check( nf90_get_att(ncid, nf90_global, 'DX', wrf%dom(id)%dx) )
   call check( nf90_get_att(ncid, nf90_global, 'DY', wrf%dom(id)%dy) )
   call check( nf90_get_att(ncid, nf90_global, 'DT', dt) )
   print*,'dt from wrfinput is: ',dt
   print*,'Using dt from namelist.input: ',wrf%dom(id)%dt
   if(debug) write(*,*) ' dx, dy, dt are ',wrf%dom(id)%dx, &
        wrf%dom(id)%dy, wrf%dom(id)%dt

   call check( nf90_get_att(ncid, nf90_global, 'MAP_PROJ', wrf%dom(id)%map_proj) )
   if(debug) write(*,*) ' map_proj is ',wrf%dom(id)%map_proj

   call check( nf90_get_att(ncid, nf90_global, 'CEN_LAT', wrf%dom(id)%cen_lat) )
   if(debug) write(*,*) ' cen_lat is ',wrf%dom(id)%cen_lat

   call check( nf90_get_att(ncid, nf90_global, 'CEN_LON', wrf%dom(id)%cen_lon) )

   call check( nf90_get_att(ncid, nf90_global, 'TRUELAT1', truelat1) )
   if(debug) write(*,*) ' truelat1 is ',truelat1

   call check( nf90_get_att(ncid, nf90_global, 'TRUELAT2', truelat2) )
   if(debug) write(*,*) ' truelat2 is ',truelat2

   call check( nf90_get_att(ncid, nf90_global, 'STAND_LON', stdlon) )

   call check( nf90_inq_varid(ncid, "P_TOP", var_id) )
   call check( nf90_get_var(ncid, var_id, wrf%dom(id)%p_top) )

!  get 1D (z) static data defining grid levels

   allocate(wrf%dom(id)%dn(1:wrf%dom(id)%bt))
   call check( nf90_inq_varid(ncid, "DN", var_id) )
   call check( nf90_get_var(ncid, var_id, wrf%dom(id)%dn) )
   if(debug) write(*,*) ' dn ',wrf%dom(id)%dn

   allocate(wrf%dom(id)%znu(1:wrf%dom(id)%bt))
   call check( nf90_inq_varid(ncid, "ZNU", var_id) )
   call check( nf90_get_var(ncid, var_id, wrf%dom(id)%znu) )
   if(debug) write(*,*) ' znu is ',wrf%dom(id)%znu

   allocate(wrf%dom(id)%dnw(1:wrf%dom(id)%bt))
   call check( nf90_inq_varid(ncid, "DNW", var_id) )
   call check( nf90_get_var(ncid, var_id, wrf%dom(id)%dnw) )
   if(debug) write(*,*) ' dnw is ',wrf%dom(id)%dnw

   allocate(wrf%dom(id)%zs(1:wrf%dom(id)%sls))
   call check( nf90_inq_varid(ncid, "ZS", var_id) )
   call check( nf90_get_var(ncid, var_id, wrf%dom(id)%zs) )

!  get 2D (x,y) base state for mu, latitude, longitude

   allocate(wrf%dom(id)%mub(1:wrf%dom(id)%we,1:wrf%dom(id)%sn))
   call check( nf90_inq_varid(ncid, "MUB", var_id) )
   call check( nf90_get_var(ncid, var_id, wrf%dom(id)%mub) )
   if(debug) then
      write(*,*) ' corners of mub '
      write(*,*) wrf%dom(id)%mub(1,1),wrf%dom(id)%mub(wrf%dom(id)%we,1),  &
           wrf%dom(id)%mub(1,wrf%dom(id)%sn),wrf%dom(id)%mub(wrf%dom(id)%we, &
           wrf%dom(id)%sn)
   end if

   allocate(wrf%dom(id)%longitude(1:wrf%dom(id)%we,1:wrf%dom(id)%sn))
   call check( nf90_inq_varid(ncid, "XLONG", var_id) )
   call check( nf90_get_var(ncid, var_id, wrf%dom(id)%longitude) )

   allocate(wrf%dom(id)%latitude(1:wrf%dom(id)%we,1:wrf%dom(id)%sn))
   call check( nf90_inq_varid(ncid, "XLAT", var_id) )
   call check( nf90_get_var(ncid, var_id, wrf%dom(id)%latitude) )

   allocate(wrf%dom(id)%land(1:wrf%dom(id)%we,1:wrf%dom(id)%sn))
   call check( nf90_inq_varid(ncid, "XLAND", var_id) )
   call check( nf90_get_var(ncid, var_id, wrf%dom(id)%land) )
   if(debug) then
      write(*,*) ' corners of land '
      write(*,*) wrf%dom(id)%land(1,1),wrf%dom(id)%land(wrf%dom(id)%we,1),  &
           wrf%dom(id)%land(1,wrf%dom(id)%sn),wrf%dom(id)%land(wrf%dom(id)%we, &
           wrf%dom(id)%sn)
   end if

   if(debug) then
      write(*,*) ' corners of lat '
      write(*,*) wrf%dom(id)%latitude(1,1),wrf%dom(id)%latitude(wrf%dom(id)%we,1),  &
           wrf%dom(id)%latitude(1,wrf%dom(id)%sn), &
           wrf%dom(id)%latitude(wrf%dom(id)%we,wrf%dom(id)%sn)
      write(*,*) ' corners of long '
      write(*,*) wrf%dom(id)%longitude(1,1),wrf%dom(id)%longitude(wrf%dom(id)%we,1),  &
           wrf%dom(id)%longitude(1,wrf%dom(id)%sn), &
           wrf%dom(id)%longitude(wrf%dom(id)%we,wrf%dom(id)%sn)
   end if

   allocate(wrf%dom(id)%mapfac_m(1:wrf%dom(id)%we,1:wrf%dom(id)%sn))
   call check( nf90_inq_varid(ncid, "MAPFAC_M", var_id) )
   call check( nf90_get_var(ncid, var_id, wrf%dom(id)%mapfac_m) )

   allocate(wrf%dom(id)%hgt(1:wrf%dom(id)%we,1:wrf%dom(id)%sn))
   call check( nf90_inq_varid(ncid, "HGT", var_id) )
   call check( nf90_get_var(ncid, var_id, wrf%dom(id)%hgt) )

   allocate(wrf%dom(id)%mapfac_u(1:wrf%dom(id)%wes,1:wrf%dom(id)%sn))
   call check( nf90_inq_varid(ncid, "MAPFAC_U", var_id) )
   call check( nf90_get_var(ncid, var_id, wrf%dom(id)%mapfac_u) )

   allocate(wrf%dom(id)%mapfac_v(1:wrf%dom(id)%we,1:wrf%dom(id)%sns))
   call check( nf90_inq_varid(ncid, "MAPFAC_V", var_id) )
   call check( nf90_get_var(ncid, var_id, wrf%dom(id)%mapfac_v) )

! get 3D base state geopotential

   allocate(wrf%dom(id)%phb(1:wrf%dom(id)%we,1:wrf%dom(id)%sn,1:wrf%dom(id)%bts))
   call check( nf90_inq_varid(ncid, "PHB", var_id) )
   call check( nf90_get_var(ncid, var_id, wrf%dom(id)%phb) )
   if(debug) then
      write(*,*) ' corners of phb '
      write(*,*) wrf%dom(id)%phb(1,1,1),wrf%dom(id)%phb(wrf%dom(id)%we,1,1),  &
           wrf%dom(id)%phb(1,wrf%dom(id)%sn,1),wrf%dom(id)%phb(wrf%dom(id)%we, &
           wrf%dom(id)%sn,1)
      write(*,*) wrf%dom(id)%phb(1,1,wrf%dom(id)%bts), &
           wrf%dom(id)%phb(wrf%dom(id)%we,1,wrf%dom(id)%bts),  &
           wrf%dom(id)%phb(1,wrf%dom(id)%sn,wrf%dom(id)%bts), &
           wrf%dom(id)%phb(wrf%dom(id)%we,wrf%dom(id)%sn,wrf%dom(id)%bts)
   end if

! close data file, we have all we need

   call check( nf90_close(ncid) )

! Initializes the map projection structure to missing values

   call map_init(wrf%dom(id)%proj)

! Populate the map projection structure

   if(wrf%dom(id)%map_proj == map_sphere) then
      proj_code = PROJ_LATLON
   elseif(wrf%dom(id)%map_proj == map_lambert) then
      proj_code = PROJ_LC
   elseif(wrf%dom(id)%map_proj == map_polar_stereo) then
      proj_code = PROJ_PS
   elseif(wrf%dom(id)%map_proj == map_mercator) then
      proj_code = PROJ_MERC
   else
      call error_handler(E_ERR,'static_init_model', &
        'Map projection no supported.', source, revision, revdate)
   endif

   call map_set(proj_code,wrf%dom(id)%latitude(1,1),wrf%dom(id)%longitude(1,1), &
        1.0_r8,1.0_r8,wrf%dom(id)%dx,stdlon,truelat1,truelat2,wrf%dom(id)%proj)

!  build the map into the 1D DART vector for WRF data

   wrf%dom(id)%number_of_wrf_variables = 7 + wrf%dom(id)%n_moist
   if( wrf%dom(id)%surf_obs ) then
      wrf%dom(id)%number_of_wrf_variables = wrf%dom(id)%number_of_wrf_variables + 6
   endif
   if( wrf%dom(id)%soil_data ) then
      wrf%dom(id)%number_of_wrf_variables = wrf%dom(id)%number_of_wrf_variables + 3
   endif
   if( h_diab ) then
      wrf%dom(id)%number_of_wrf_variables = wrf%dom(id)%number_of_wrf_variables + 1
   endif
   allocate(wrf%dom(id)%var_type(wrf%dom(id)%number_of_wrf_variables))

   allocate(wrf%dom(id)%dart_kind(wrf%dom(id)%number_of_wrf_variables))
   wrf%dom(id)%var_type(1)  = TYPE_U
   wrf%dom(id)%dart_kind(1) = KIND_U_WIND_COMPONENT
   wrf%dom(id)%var_type(2)  = TYPE_V
   wrf%dom(id)%dart_kind(2) = KIND_V_WIND_COMPONENT
   wrf%dom(id)%var_type(3)  = TYPE_W
   wrf%dom(id)%dart_kind(3) = KIND_VERTICAL_VELOCITY
   wrf%dom(id)%var_type(4)  = TYPE_GZ
   wrf%dom(id)%dart_kind(4) = KIND_GEOPOTENTIAL_HEIGHT
   wrf%dom(id)%var_type(5)  = TYPE_T
   wrf%dom(id)%dart_kind(5) = KIND_TEMPERATURE
   wrf%dom(id)%var_type(6)  = TYPE_MU
   wrf%dom(id)%dart_kind(6) = KIND_PRESSURE
   wrf%dom(id)%var_type(7)  = TYPE_TSK
   wrf%dom(id)%dart_kind(7) = KIND_TEMPERATURE

   ind = 7
   if( wrf%dom(id)%n_moist >= 1) then
      ind = ind + 1
      wrf%dom(id)%var_type(ind)  = TYPE_QV
      wrf%dom(id)%dart_kind(ind) = KIND_VAPOR_MIXING_RATIO
   end if
   if( wrf%dom(id)%n_moist >= 2) then
      ind = ind + 1
      wrf%dom(id)%var_type(ind)  = TYPE_QC
      wrf%dom(id)%dart_kind(ind) = KIND_CLOUD_LIQUID_WATER
   end if
   if( wrf%dom(id)%n_moist >= 3) then
      ind = ind + 1
      wrf%dom(id)%var_type(ind) = TYPE_QR
      wrf%dom(id)%dart_kind(ind) = KIND_RAINWATER_MIXING_RATIO
   end if
   if( wrf%dom(id)%n_moist >= 4) then
      ind = ind + 1
      wrf%dom(id)%var_type(ind) = TYPE_QI
      wrf%dom(id)%dart_kind(ind) = KIND_CLOUD_ICE
   end if
   if( wrf%dom(id)%n_moist >= 5) then
      ind = ind + 1
      wrf%dom(id)%var_type(ind) = TYPE_QS
      wrf%dom(id)%dart_kind(ind) = KIND_SNOW_MIXING_RATIO
   end if
   if( wrf%dom(id)%n_moist >= 6) then
      ind = ind + 1
      wrf%dom(id)%var_type(ind) = TYPE_QG
      wrf%dom(id)%dart_kind(ind) = KIND_GRAUPEL_MIXING_RATIO
   end if
   if( wrf%dom(id)%n_moist == 7) then
      ind = ind + 1
      wrf%dom(id)%var_type(ind) = TYPE_QNICE
      wrf%dom(id)%dart_kind(ind) = KIND_ICE_NUMBER_CONCENTRATION
   end if
   if( wrf%dom(id)%surf_obs ) then
      ind = ind + 1
      wrf%dom(id)%var_type(ind) = TYPE_U10
      wrf%dom(id)%dart_kind(ind) = KIND_U_WIND_COMPONENT
      ind = ind + 1
      wrf%dom(id)%var_type(ind) = TYPE_V10
      wrf%dom(id)%dart_kind(ind) = KIND_V_WIND_COMPONENT
      ind = ind + 1
      wrf%dom(id)%var_type(ind) = TYPE_T2
      wrf%dom(id)%dart_kind(ind) = KIND_TEMPERATURE
      ind = ind + 1
      wrf%dom(id)%var_type(ind) = TYPE_TH2
      wrf%dom(id)%dart_kind(ind) = KIND_POTENTIAL_TEMPERATURE
      ind = ind + 1
      wrf%dom(id)%var_type(ind) = TYPE_Q2
      wrf%dom(id)%dart_kind(ind) = KIND_SPECIFIC_HUMIDITY
      ind = ind + 1
      wrf%dom(id)%var_type(ind) = TYPE_PS
      wrf%dom(id)%dart_kind(ind) = KIND_PRESSURE
   end if
   if( wrf%dom(id)%soil_data ) then
      ind = ind + 1
      wrf%dom(id)%var_type(ind)  = TYPE_TSLB
      wrf%dom(id)%dart_kind(ind) = KIND_TEMPERATURE
      ind = ind + 1
      wrf%dom(id)%var_type(ind)  = TYPE_SMOIS
      wrf%dom(id)%dart_kind(ind) = KIND_SOIL_MOISTURE
      ind = ind + 1
      wrf%dom(id)%var_type(ind)  = TYPE_SH2O
      wrf%dom(id)%dart_kind(ind) = KIND_SOIL_MOISTURE
   end if
   if( h_diab ) then
      ind = ind + 1
      wrf%dom(id)%var_type(ind) = TYPE_HDIAB
      wrf%dom(id)%dart_kind(ind) = KIND_CONDENSATIONAL_HEATING
   end if

! indices into 1D array
   allocate(wrf%dom(id)%dart_ind(wrf%dom(id)%wes,wrf%dom(id)%sns,wrf%dom(id)%bts,num_model_var_types))
   allocate(wrf%dom(id)%var_index(2,wrf%dom(id)%number_of_wrf_variables))
! dimension of variables
   allocate(wrf%dom(id)%var_size(3,wrf%dom(id)%number_of_wrf_variables))

   wrf%dom(id)%dart_ind = 0

   ind = 1                         ! *** u field ***
   wrf%dom(id)%var_size(1,ind) = wrf%dom(id)%wes
   wrf%dom(id)%var_size(2,ind) = wrf%dom(id)%sn
   wrf%dom(id)%var_size(3,ind) = wrf%dom(id)%bt
   wrf%dom(id)%var_index(1,ind) = dart_index
   do k=1,wrf%dom(id)%var_size(3,ind)
      do j=1,wrf%dom(id)%var_size(2,ind)
         do i=1,wrf%dom(id)%var_size(1,ind)
            wrf%dom(id)%dart_ind(i,j,k,TYPE_U) = dart_index
            dart_index = dart_index + 1
         enddo
      enddo
   enddo
   wrf%dom(id)%var_index(2,ind) = dart_index - 1

   ind = ind + 1                   ! *** v field ***
   wrf%dom(id)%var_size(1,ind) = wrf%dom(id)%we
   wrf%dom(id)%var_size(2,ind) = wrf%dom(id)%sns
   wrf%dom(id)%var_size(3,ind) = wrf%dom(id)%bt
   wrf%dom(id)%var_index(1,ind) = dart_index
   do k=1,wrf%dom(id)%var_size(3,ind)
      do j=1,wrf%dom(id)%var_size(2,ind)
         do i=1,wrf%dom(id)%var_size(1,ind)
            wrf%dom(id)%dart_ind(i,j,k,TYPE_V) = dart_index
            dart_index = dart_index + 1
         enddo
      enddo
   enddo
   wrf%dom(id)%var_index(2,ind) = dart_index - 1

   ind = ind + 1                   ! *** w field ***
   wrf%dom(id)%var_size(1,ind) = wrf%dom(id)%we
   wrf%dom(id)%var_size(2,ind) = wrf%dom(id)%sn
   wrf%dom(id)%var_size(3,ind) = wrf%dom(id)%bts
   wrf%dom(id)%var_index(1,ind) = dart_index
   do k=1,wrf%dom(id)%var_size(3,ind)
      do j=1,wrf%dom(id)%var_size(2,ind)
         do i=1,wrf%dom(id)%var_size(1,ind)
            wrf%dom(id)%dart_ind(i,j,k,TYPE_W) = dart_index
            dart_index = dart_index + 1
         enddo
      enddo
   enddo
   wrf%dom(id)%var_index(2,ind) = dart_index - 1

   ind = ind + 1                   ! *** geopotential field ***
   wrf%dom(id)%var_size(1,ind) = wrf%dom(id)%we
   wrf%dom(id)%var_size(2,ind) = wrf%dom(id)%sn
   wrf%dom(id)%var_size(3,ind) = wrf%dom(id)%bts
   wrf%dom(id)%var_index(1,ind) = dart_index
   do k=1,wrf%dom(id)%var_size(3,ind)
      do j=1,wrf%dom(id)%var_size(2,ind)
         do i=1,wrf%dom(id)%var_size(1,ind)
            wrf%dom(id)%dart_ind(i,j,k,TYPE_GZ) = dart_index
            dart_index = dart_index + 1
         enddo
      enddo
   enddo
   wrf%dom(id)%var_index(2,ind) = dart_index - 1

   ind = ind + 1                   ! *** theta field ***
   wrf%dom(id)%var_size(1,ind) = wrf%dom(id)%we
   wrf%dom(id)%var_size(2,ind) = wrf%dom(id)%sn
   wrf%dom(id)%var_size(3,ind) = wrf%dom(id)%bt
   wrf%dom(id)%var_index(1,ind) = dart_index
   do k=1,wrf%dom(id)%var_size(3,ind)
      do j=1,wrf%dom(id)%var_size(2,ind)
         do i=1,wrf%dom(id)%var_size(1,ind)
            wrf%dom(id)%dart_ind(i,j,k,TYPE_T) = dart_index
            dart_index = dart_index + 1
         enddo
      enddo
   enddo
   wrf%dom(id)%var_index(2,ind) = dart_index - 1

   ind = ind + 1                   ! *** mu field ***
   wrf%dom(id)%var_size(1,ind) = wrf%dom(id)%we
   wrf%dom(id)%var_size(2,ind) = wrf%dom(id)%sn
   wrf%dom(id)%var_size(3,ind) = 1
   wrf%dom(id)%var_index(1,ind) = dart_index
   do k=1,wrf%dom(id)%var_size(3,ind)
      do j=1,wrf%dom(id)%var_size(2,ind)
         do i=1,wrf%dom(id)%var_size(1,ind)
            wrf%dom(id)%dart_ind(i,j,k,TYPE_MU) = dart_index
            dart_index = dart_index + 1
         enddo
      enddo
   enddo
   wrf%dom(id)%var_index(2,ind) = dart_index - 1

   ind = ind + 1                   ! *** tsk field ***
   wrf%dom(id)%var_size(1,ind) = wrf%dom(id)%we
   wrf%dom(id)%var_size(2,ind) = wrf%dom(id)%sn
   wrf%dom(id)%var_size(3,ind) = 1
   wrf%dom(id)%var_index(1,ind) = dart_index
   do k=1,wrf%dom(id)%var_size(3,ind)
      do j=1,wrf%dom(id)%var_size(2,ind)
         do i=1,wrf%dom(id)%var_size(1,ind)
            wrf%dom(id)%dart_ind(i,j,k,TYPE_TSK) = dart_index
            dart_index = dart_index + 1
         enddo
      enddo
   enddo
   wrf%dom(id)%var_index(2,ind) = dart_index - 1

   do model_type = TYPE_QV, TYPE_QV + wrf%dom(id)%n_moist - 1
      ind = ind + 1                   ! *** moisture field ***
      wrf%dom(id)%var_size(1,ind) = wrf%dom(id)%we
      wrf%dom(id)%var_size(2,ind) = wrf%dom(id)%sn
      wrf%dom(id)%var_size(3,ind) = wrf%dom(id)%bt
      wrf%dom(id)%var_index(1,ind) = dart_index
      do k=1,wrf%dom(id)%var_size(3,ind)
         do j=1,wrf%dom(id)%var_size(2,ind)
            do i=1,wrf%dom(id)%var_size(1,ind)
               wrf%dom(id)%dart_ind(i,j,k,model_type) = dart_index
               dart_index = dart_index + 1
            enddo
         enddo
      enddo
      wrf%dom(id)%var_index(2,ind) = dart_index - 1
   enddo

   if(wrf%dom(id)%surf_obs ) then
      do model_type = TYPE_U10, TYPE_PS
         ind = ind + 1                   ! *** Surface variable ***
         wrf%dom(id)%var_size(1,ind) = wrf%dom(id)%we
         wrf%dom(id)%var_size(2,ind) = wrf%dom(id)%sn
         wrf%dom(id)%var_size(3,ind) = 1
         wrf%dom(id)%var_index(1,ind) = dart_index
         do k=1,wrf%dom(id)%var_size(3,ind)
            do j=1,wrf%dom(id)%var_size(2,ind)
               do i=1,wrf%dom(id)%var_size(1,ind)
                  wrf%dom(id)%dart_ind(i,j,k,model_type) = dart_index
                  dart_index = dart_index + 1
               enddo
            enddo
         enddo
         wrf%dom(id)%var_index(2,ind) = dart_index - 1
      enddo
   end if

   if(wrf%dom(id)%soil_data ) then
      ind = ind + 1                   ! *** tslb field ***
      wrf%dom(id)%var_size(1,ind) = wrf%dom(id)%we
      wrf%dom(id)%var_size(2,ind) = wrf%dom(id)%sn
      wrf%dom(id)%var_size(3,ind) = wrf%dom(id)%sls
      wrf%dom(id)%var_index(1,ind) = dart_index
      do k=1,wrf%dom(id)%var_size(3,ind)
         do j=1,wrf%dom(id)%var_size(2,ind)
            do i=1,wrf%dom(id)%var_size(1,ind)
               wrf%dom(id)%dart_ind(i,j,k,TYPE_TSLB) = dart_index
               dart_index = dart_index + 1
            enddo
         enddo
      enddo
      wrf%dom(id)%var_index(2,ind) = dart_index - 1

      ind = ind + 1                   ! *** smois field ***
      wrf%dom(id)%var_size(1,ind) = wrf%dom(id)%we
      wrf%dom(id)%var_size(2,ind) = wrf%dom(id)%sn
      wrf%dom(id)%var_size(3,ind) = wrf%dom(id)%sls
      wrf%dom(id)%var_index(1,ind) = dart_index
      do k=1,wrf%dom(id)%var_size(3,ind)
         do j=1,wrf%dom(id)%var_size(2,ind)
            do i=1,wrf%dom(id)%var_size(1,ind)
               wrf%dom(id)%dart_ind(i,j,k,TYPE_SMOIS) = dart_index
               dart_index = dart_index + 1
            enddo
         enddo
      enddo
      wrf%dom(id)%var_index(2,ind) = dart_index - 1  
      
      ind = ind + 1                   ! *** sh2o field ***
      wrf%dom(id)%var_size(1,ind) = wrf%dom(id)%we
      wrf%dom(id)%var_size(2,ind) = wrf%dom(id)%sn
      wrf%dom(id)%var_size(3,ind) = wrf%dom(id)%sls
      wrf%dom(id)%var_index(1,ind) = dart_index
      do k=1,wrf%dom(id)%var_size(3,ind)
         do j=1,wrf%dom(id)%var_size(2,ind)
            do i=1,wrf%dom(id)%var_size(1,ind)
               wrf%dom(id)%dart_ind(i,j,k,TYPE_SH2O) = dart_index
               dart_index = dart_index + 1
            enddo
         enddo
      enddo
      wrf%dom(id)%var_index(2,ind) = dart_index - 1  
   end if

   if(h_diab ) then
      ind = ind + 1                   ! *** h_diabatic variable ***
      wrf%dom(id)%var_size(1,ind) = wrf%dom(id)%we
      wrf%dom(id)%var_size(2,ind) = wrf%dom(id)%sn
      wrf%dom(id)%var_size(3,ind) = wrf%dom(id)%bt
      wrf%dom(id)%var_index(1,ind) = dart_index
      do k=1,wrf%dom(id)%var_size(3,ind)
         do j=1,wrf%dom(id)%var_size(2,ind)
            do i=1,wrf%dom(id)%var_size(1,ind)
               wrf%dom(id)%dart_ind(i,j,k,TYPE_HDIAB) = dart_index
               dart_index = dart_index + 1
            enddo
         enddo
      enddo
      wrf%dom(id)%var_index(2,ind) = dart_index - 1
   end if

enddo

write(*,*)

wrf%model_size = dart_index - 1
allocate (ens_mean(wrf%model_size))
if(debug) write(*,*) ' wrf model size is ',wrf%model_size

contains

  ! Internal subroutine - checks error status after each netcdf, prints
  !                       text message each time an error code is returned.
  subroutine check(istatus)
    integer, intent (in) :: istatus

    if(istatus /= nf90_noerr) call error_handler(E_ERR, 'static_init_model', &
       trim(nf90_strerror(istatus)), source, revision, revdate)

  end subroutine check

end subroutine static_init_model


!#######################################################################

function get_model_size()

integer :: get_model_size

get_model_size = wrf%model_size

end function get_model_size

!#######################################################################

function get_model_time_step()
!------------------------------------------------------------------------
! function get_model_time_step()
!
! Returns the time step of the model. In the long run should be replaced
! by a more general routine that returns details of a general time-stepping
! capability.
!
! toward that end ... we are now reading a namelist variable for the
! width of the assimilation time window.

type(time_type) :: get_model_time_step
integer :: model_dt, assim_dt

! We need to coordinate the desired assimilation window to be a 
! multiple of the model time step (which has no precision past integer seconds).

model_dt = nint(wrf%dom(1)%dt)

! The integer arithmetic does its magic.
assim_dt = (assimilation_period_seconds / model_dt) * model_dt

get_model_time_step = set_time(assim_dt)

end function get_model_time_step


!#######################################################################


subroutine get_state_meta_data(index_in, location, var_type_out, id_out)

! Given an integer index into the DART state vector structure, returns the
! associated location. This is not a function because the more general
! form of the call has a second intent(out) optional argument kind.
! Maybe a functional form should be added?

integer,             intent(in)  :: index_in
type(location_type), intent(out) :: location
integer, optional,   intent(out) :: var_type_out, id_out

integer  :: var_type, dart_type
integer  :: index, ip, jp, kp
integer  :: nz, ny, nx
logical  :: var_found
real(r8) :: lon, lat, lev

integer :: i, id
logical, parameter :: debug = .false.
character(len=129) :: errstring

if(debug) then
   write(errstring,*)' index_in = ',index_in
   call error_handler(E_MSG,'get_state_meta_data',errstring,' ',' ',' ')
endif

! index_in can be negative if ob is identity ob...
index = abs(index_in)

var_found = .false.

!  first find var_type

if(debug) then
   do id=1,num_domains
      do i=1, wrf%dom(id)%number_of_wrf_variables
         write(errstring,*)' domain, var, var_type(i) = ',id,i,wrf%dom(id)%var_type(i)
         call error_handler(E_MSG,'get_state_meta_data',errstring,' ',' ',' ')
      enddo
   enddo
endif

! first find var_type and domain id
i = 0
id = 1
do while (.not. var_found)
   i = i + 1
   if(i .gt. wrf%dom(id)%number_of_wrf_variables) then
      i = 1
      if (id < num_domains) then
         id = id + 1
      else
         write(errstring,*)' size of vector ',wrf%model_size
         call error_handler(E_MSG,'get_state_meta_data', errstring, ' ', ' ', ' ')
         write(errstring,*)' dart_index ',index_in
         call error_handler(E_ERR,'get_state_meta_data', 'index out of range', &
              source, revision, revdate)
      end if
   end if
   if( (index .ge. wrf%dom(id)%var_index(1,i) ) .and.  &
       (index .le. wrf%dom(id)%var_index(2,i) )       )  then
      var_found = .true.
      var_type  = wrf%dom(id)%var_type(i)
      dart_type = wrf%dom(id)%dart_kind(i)
      index = index - wrf%dom(id)%var_index(1,i) + 1
   end if
end do

!  now find i,j,k location.
!  index has been normalized such that it is relative to
!  array starting at (1,1,1)

nx = wrf%dom(id)%var_size(1,i)
ny = wrf%dom(id)%var_size(2,i)
nz = wrf%dom(id)%var_size(3,i)

kp = 1 + (index-1)/(nx*ny)
jp = 1 + (index - (kp-1)*nx*ny - 1)/nx
ip = index - (kp-1)*nx*ny - (jp-1)*nx

! at this point, (ip,jp,kp) refer to indices in the variable's own grid

if(debug) write(*,*) ' ip, jp, kp for index ',ip,jp,kp,index
if(debug) write(*,*) ' Var type: ',var_type

! first obtain lat/lon from (ip,jp)
call get_wrf_horizontal_location( ip, jp, var_type, id, lon, lat )

! now convert to desired vertical coordinate (defined in the namelist)
if (wrf%dom(id)%vert_coord == VERTISLEVEL) then
   ! here we need level index of mass grid
   if( (var_type == type_w ) .or. (var_type == type_gz) ) then
      lev = real(kp) - 0.5_r8
   else
      lev = real(kp)
   endif
elseif (wrf%dom(id)%vert_coord == VERTISPRESSURE) then
   ! directly convert to pressure
   lev = model_pressure(ip,jp,kp,id,var_type,ens_mean)
elseif (wrf%dom(id)%vert_coord == VERTISHEIGHT) then
   lev = model_height(ip,jp,kp,id,var_type,ens_mean)
endif

if(debug) write(*,*) 'lon, lat, lev: ',lon, lat, lev

! convert to DART location type
location = set_location(lon, lat, lev, wrf%dom(id)%vert_coord)

! return DART variable kind if requested
if(present(var_type_out)) var_type_out = dart_type

! return domain id if requested
if(present(id_out)) id_out = id

end subroutine get_state_meta_data

!#######################################################################


subroutine model_interpolate(x, location, obs_kind, obs_val, istatus)

! This is the main forward operator subroutine for WRF.
! Given an ob (its DART location and kind), the corresponding model
! value is computed at nearest i,j,k. Thus, first i,j,k is obtained
! from ob lon,lat,z and then the state value that corresponds to
! the ob kind is interpolated.

! No location conversions are carried out in this subroutine. See
! get_close_obs, where ob vertical location information is converted
! to the requested vertical coordinate type.

! x:       Full DART state vector relevant to what's being updated
!          in the filter (mean or individual members).
! istatus: Returned 0 if everything is OK, 1 if error occured.
!                  -1 if the station height is lower than the lowest model level 
!                     while the station is located inside the horizontal model domain.

! modified 26 June 2006 to accomodate vortex attributes
! modified 13 December 2006 to accomodate changes for the mpi version

real(r8),            intent(in) :: x(:)
type(location_type), intent(in) :: location
integer,             intent(in) :: obs_kind
real(r8),           intent(out) :: obs_val
integer,            intent(out) :: istatus

logical, parameter  :: debug = .false.
real(r8)            :: xloc, yloc, zloc, xloc_u, yloc_v, xyz_loc(3)
integer             :: i, i_u, j, j_v, k, k2, i1,i2
real(r8)            :: dx,dy,dz,dxm,dym,dzm,dx_u,dxm_u,dy_v,dym_v
real(r8)            :: a1,utrue,vtrue,ugrid,vgrid
integer             :: id

character(len=129) :: errstring

real(r8), dimension(2) :: fld
real(r8), allocatable, dimension(:) :: v_h, v_p

! local vars, used in finding sea-level pressure and vortex center
real(r8), allocatable, dimension(:)   :: t1d, p1d, qv1d, z1d
real(r8), allocatable, dimension(:,:) :: psea, pp, pd
real(r8), allocatable, dimension(:)   :: x1d,y1d,xx1d,yy1d
integer  :: xlen, ylen, xxlen, yylen, ii1, ii2
real(r8) :: clat, clon, cxmin, cymin

! center_track_*** used to define center search area
integer :: center_track_xmin, center_track_ymin, &
           center_track_xmax, center_track_ymax

! local vars, used in calculating density, pressure, height
real(r8)            :: rho1 , rho2 , rho3, rho4
real(r8)            :: pres1, pres2, pres3, pres4, pres


istatus = 0
obs_val = missing_r8


! - - - - - - - - - !

! This part determines model i,j,k from ob lon,lat,z (location).
! This information is only needed if this is not an identity ob.
! For an identity ob, (absolute) obs_kind is directly the index
! into the state vector.

xyz_loc = get_location(location)

if (obs_kind > 0) then
   ! first obtain domain id, and mass points (i,j)
   call get_domain_info(xyz_loc(1),xyz_loc(2),id,xloc,yloc)

   if (id==0) then
      istatus = 1
      return
   endif

   if(debug) then
      i = xloc
      j = yloc
      print*,xyz_loc(2), xyz_loc(1), xloc,yloc
      write(*,*) ' corners of lat '
      write(*,*) wrf%dom(id)%latitude(i,j),wrf%dom(id)%latitude(i+1,j),  &
           wrf%dom(id)%latitude(i,j+1), &
           wrf%dom(id)%latitude(i+1,j+1)
      write(*,*) ' corners of long '
      write(*,*) wrf%dom(id)%longitude(i,j),wrf%dom(id)%longitude(i+1,j),  &
           wrf%dom(id)%longitude(i,j+1), &
           wrf%dom(id)%longitude(i+1,j+1)
   endif

   allocate(v_h(0:wrf%dom(id)%bt), v_p(0:wrf%dom(id)%bt))

   ! get integer (west/south) grid point and distances to neighboring grid points
   ! distances are used as weights to carry out horizontal interpolations
   call toGrid(xloc,i,dx,dxm)
   call toGrid(yloc,j,dy,dym)

   ! Determine corresponding model level for obs location
   ! This depends on the obs vertical coordinate
   if(vert_is_level(location)) then
      ! Ob is by model level
      zloc = xyz_loc(3)

   elseif(vert_is_pressure(location)) then
      ! Ob is by pressure: get corresponding mass level zloc from
      ! computed column pressure profile
      call get_model_pressure_profile(i,j,dx,dy,dxm,dym,wrf%dom(id)%bt,x,id,v_p)
      ! get pressure vertical co-ordinate
      call pres_to_zk(xyz_loc(3), v_p, wrf%dom(id)%bt,zloc)
      if(debug.and.obs_kind /= KIND_SURFACE_PRESSURE) &
                print*,' obs is by pressure and zloc =',zloc
      if(debug) print*,'model pressure profile'
      if(debug) print*,v_p
      
   elseif(vert_is_height(location)) then
      ! Ob is by height: get corresponding mass level zloc from
      ! computed column height profile
      call get_model_height_profile(i,j,dx,dy,dxm,dym,wrf%dom(id)%bt,x,id,v_h)
      if(debug) print*,'model height profile'
      if(debug) print*,v_h
      ! get height vertical co-ordinate
      call height_to_zk(xyz_loc(3), v_h, wrf%dom(id)%bt,zloc)
      if(debug) print*,' obs is by height and zloc =',zloc
   
   elseif(vert_is_surface(location)) then
      zloc = 1.0_r8
      if(debug) print*,' obs is at the surface = ', xyz_loc(3)

   elseif(vert_is_undef(location)) then
      zloc  = missing_r8

   else
      write(errstring,*) 'wrong option for which_vert ', &
                         nint(query_location(location,'which_vert'))
      call error_handler(E_ERR,'model_interpolate', errstring, &
           source, revision, revdate)

   endif

   surf_var = .false.
   if(vert_is_surface(location)) surf_var = .true.
   if(zloc == missing_r8) then
     if(xyz_loc(3) < v_h(0)) then
      zloc = 1.0_r8         ! Lower than the model terrain.
      surf_var = .true.     ! Estimate U,V,T,and Q from the model sfc states.
     else                
      obs_val = missing_r8
      istatus = 1
      deallocate(v_h, v_p)
      return
     endif
   endif

   k = max(1,int(zloc))

endif

! - - - - - - - - - !

! This part is the forward operator -- compute desired model state value for given point.

! Get the desired field to be interpolated
if( obs_kind == KIND_U_WIND_COMPONENT .or. obs_kind == KIND_V_WIND_COMPONENT) then   ! U, V

   if(.not. vert_is_surface(location) .or. .not. surf_var) then

      xloc_u = xloc + 0.5
      yloc_v = yloc + 0.5
      call toGrid(xloc_u,i_u,dx_u,dxm_u)
      call toGrid(yloc_v,j_v,dy_v,dym_v)

      if(  i_u >= 1 .and. i_u < wrf%dom(id)%var_size(1,TYPE_U) .and. &
           j   >= 1 .and. j   < wrf%dom(id)%var_size(2,TYPE_U) .and. &
           i   >= 1 .and. i   < wrf%dom(id)%var_size(1,TYPE_V) .and. &
           j_v >= 1 .and. j_v < wrf%dom(id)%var_size(2,TYPE_V) .and. &
           k   >= 1 .and. k   < wrf%dom(id)%var_size(3,TYPE_U)) then

         do k2=1,2

            ! For memory economy over speed, dart_ind should be removed from
            ! the wrf structure and get_wrf_index reactivated.

!!$            i1 = get_wrf_index(i_u,j  ,k+k2-1,TYPE_U,id)
!!$            i2 = get_wrf_index(i_u,j+1,k+k2-1,TYPE_U,id)

            i1 = wrf%dom(id)%dart_ind(i_u,j  ,k+k2-1,TYPE_U)
            i2 = wrf%dom(id)%dart_ind(i_u,j+1,k+k2-1,TYPE_U)

            ugrid = dym*( dxm_u*x(i1) + dx_u*x(i1+1) ) + dy*( dxm_u*x(i2) + dx_u*x(i2+1) )

!!$            i1 = get_wrf_index(i,j_v  ,k+k2-1,TYPE_V,id)
!!$            i2 = get_wrf_index(i,j_v+1,k+k2-1,TYPE_V,id)

            i1 = wrf%dom(id)%dart_ind(i,j_v  ,k+k2-1,TYPE_V)
            i2 = wrf%dom(id)%dart_ind(i,j_v+1,k+k2-1,TYPE_V)

            vgrid = dym_v*( dxm*x(i1) + dx*x(i1+1) ) + dy_v*( dxm*x(i2) + dx*x(i2+1) )

            call gridwind_to_truewind(xyz_loc(1), wrf%dom(id)%proj, ugrid, vgrid, &
                 utrue, vtrue)

            if( obs_kind == KIND_U_WIND_COMPONENT) then

               fld(k2) = utrue

            else   ! must want v

               fld(k2) = vtrue

            endif

         enddo

      else

         fld(:) = missing_r8

      endif

   else

      if(i >= 1 .and. i < wrf%dom(id)%var_size(1,TYPE_T) .and. &
         j >= 1 .and. j < wrf%dom(id)%var_size(2,TYPE_T) .and. &
         wrf%dom(id)%surf_obs) then

         if( obs_kind == KIND_U_WIND_COMPONENT) then

!!$            i1 = get_wrf_index(i,j,1,TYPE_U10,id)
!!$            i2 = get_wrf_index(i,j+1,1,TYPE_U10,id)
            i1 = wrf%dom(id)%dart_ind(i,j,1,TYPE_U10)
            i2 = wrf%dom(id)%dart_ind(i,j+1,1,TYPE_U10)

         else
!!$            i1 = get_wrf_index(i,j,1,TYPE_V10,id)
!!$            i2 = get_wrf_index(i,j+1,1,TYPE_V10,id)
            i1 = wrf%dom(id)%dart_ind(i,j,1,TYPE_V10)
            i2 = wrf%dom(id)%dart_ind(i,j+1,1,TYPE_V10)

         endif

         fld(1) = dym*( dxm*x(i1) + dx*x(i1+1) ) + dy*( dxm*x(i2) + dx*x(i2+1) )

      else

         fld(1) = missing_r8

      endif

   endif

else if( obs_kind == KIND_TEMPERATURE ) then

   if(.not. vert_is_surface(location) .or. .not. surf_var) then

      if(i >= 1 .and. i < wrf%dom(id)%var_size(1,TYPE_T) .and. &
         j >= 1 .and. j < wrf%dom(id)%var_size(2,TYPE_T) .and. &
         k >= 1 .and. k < wrf%dom(id)%var_size(3,TYPE_T)) then

!!$         i1 = get_wrf_index(i,j  ,k,TYPE_T,id)
!!$         i2 = get_wrf_index(i,j+1,k,TYPE_T,id)
         i1 = wrf%dom(id)%dart_ind(i,j  ,k,TYPE_T)
         i2 = wrf%dom(id)%dart_ind(i,j+1,k,TYPE_T)
         a1 = dym*( dxm*x(i1) + dx*x(i1+1) ) + dy*( dxm*x(i2) + dx*x(i2+1) )
         pres1 = model_pressure_t(i  ,j  ,k,id,x)
         pres2 = model_pressure_t(i+1,j  ,k,id,x)
         pres3 = model_pressure_t(i  ,j+1,k,id,x)
         pres4 = model_pressure_t(i+1,j+1,k,id,x)
         pres = dym*( dxm*pres1 + dx*pres2 ) + dy*( dxm*pres3 + dx*pres4 )
         fld(1) = (ts0 + a1)*(pres/ps0)**kappa

!!$         i1 = get_wrf_index(i,j  ,k+1,TYPE_T,id)
!!$         i2 = get_wrf_index(i,j+1,k+1,TYPE_T,id)
         i1 = wrf%dom(id)%dart_ind(i,j  ,k+1,TYPE_T)
         i2 = wrf%dom(id)%dart_ind(i,j+1,k+1,TYPE_T)
         a1 = dym*( dxm*x(i1) + dx*x(i1+1) ) + dy*( dxm*x(i2) + dx*x(i2+1) )
         pres1 = model_pressure_t(i  ,j  ,k+1,id,x)
         pres2 = model_pressure_t(i+1,j  ,k+1,id,x)
         pres3 = model_pressure_t(i  ,j+1,k+1,id,x)
         pres4 = model_pressure_t(i+1,j+1,k+1,id,x)
         pres = dym*( dxm*pres1 + dx*pres2 ) + dy*( dxm*pres3 + dx*pres4 )
         fld(2) = (ts0 + a1)*(pres/ps0)**kappa

      else

         fld(:) = missing_r8

      endif

   else

      if(i >= 1 .and. i < wrf%dom(id)%var_size(1,TYPE_T) .and. &
         j >= 1 .and. j < wrf%dom(id)%var_size(2,TYPE_T) .and. &
         wrf%dom(id)%surf_obs) then

!!$         i1 = get_wrf_index(i,j,1,TYPE_T2,id)
!!$         i2 = get_wrf_index(i,j+1,1,TYPE_T2,id)
         i1 = wrf%dom(id)%dart_ind(i,j,1,TYPE_T2)
         i2 = wrf%dom(id)%dart_ind(i,j+1,1,TYPE_T2)
         fld(1) = dym*( dxm*x(i1) + dx*x(i1+1) ) + dy*( dxm*x(i2) + dx*x(i2+1) )

      else

         fld(1) = missing_r8

      endif

   endif

else if( obs_kind == KIND_POTENTIAL_TEMPERATURE ) then

   if(.not. vert_is_surface(location) .or. .not. surf_var) then

      if(i >= 1 .and. i < wrf%dom(id)%var_size(1,TYPE_T) .and. &
         j >= 1 .and. j < wrf%dom(id)%var_size(2,TYPE_T) .and. &
         k >= 1 .and. k < wrf%dom(id)%var_size(3,TYPE_T)) then

!        Note:  T is perturbation potential temperature (potential temperature - ts0)
!               TH2 is potential temperature at 2 m

!!$         i1 = get_wrf_index(i,j  ,k,TYPE_T,id)
!!$         i2 = get_wrf_index(i,j+1,k,TYPE_T,id)
         i1 = wrf%dom(id)%dart_ind(i,j  ,k,TYPE_T)
         i2 = wrf%dom(id)%dart_ind(i,j+1,k,TYPE_T)
         fld(1) = ts0 + dym*( dxm*x(i1) + dx*x(i1+1) ) + dy*( dxm*x(i2) + dx*x(i2+1) )

!!$         i1 = get_wrf_index(i,j  ,k+1,TYPE_T,id)
!!$         i2 = get_wrf_index(i,j+1,k+1,TYPE_T,id)
         i1 = wrf%dom(id)%dart_ind(i,j  ,k+1,TYPE_T)
         i2 = wrf%dom(id)%dart_ind(i,j+1,k+1,TYPE_T)
         fld(2) = ts0 + dym*( dxm*x(i1) + dx*x(i1+1) ) + dy*( dxm*x(i2) + dx*x(i2+1) )

      else

         fld(:) = missing_r8

      endif

   else

      if(i >= 1 .and. i < wrf%dom(id)%var_size(1,TYPE_T) .and. &
         j >= 1 .and. j < wrf%dom(id)%var_size(2,TYPE_T) .and. &
         wrf%dom(id)%surf_obs) then

!!$         i1 = get_wrf_index(i,j,1,TYPE_TH2,id)
!!$         i2 = get_wrf_index(i,j+1,1,TYPE_TH2,id)
         i1 = wrf%dom(id)%dart_ind(i,j,1,TYPE_TH2)
         i2 = wrf%dom(id)%dart_ind(i,j+1,1,TYPE_TH2)
         fld(1) = dym*( dxm*x(i1) + dx*x(i1+1) ) + dy*( dxm*x(i2) + dx*x(i2+1) )

      else

         fld(1) = missing_r8

      endif

   endif

else if( obs_kind == KIND_DENSITY ) then

   ! Rho calculated at mass points, and so is like "TYPE_T"
   if(i >= 1 .and. i < wrf%dom(id)%var_size(1,TYPE_T) .and. &
      j >= 1 .and. j < wrf%dom(id)%var_size(2,TYPE_T) .and. &
      k >= 1 .and. k < wrf%dom(id)%var_size(3,TYPE_T)) then

      ! calculate full rho at corners of interp box
      ! and interpolate to desired horizontal location

      rho1 = model_rho_t(i  ,j  ,k,id,x)
      rho2 = model_rho_t(i+1,j  ,k,id,x)
      rho3 = model_rho_t(i  ,j+1,k,id,x)
      rho4 = model_rho_t(i+1,j+1,k,id,x)
      fld(1) = dym*( dxm*rho1 + dx*rho2 ) + dy*( dxm*rho3 + dx*rho4 )

      rho1 = model_rho_t(i  ,j  ,k+1,id,x)
      rho2 = model_rho_t(i+1,j  ,k+1,id,x)
      rho3 = model_rho_t(i  ,j+1,k+1,id,x)
      rho4 = model_rho_t(i+1,j+1,k+1,id,x)
      fld(2) = dym*( dxm*rho1 + dx*rho2 ) + dy*( dxm*rho3 + dx*rho4 )

   else

      fld(:) = missing_r8

   endif

else if( obs_kind == KIND_VERTICAL_VELOCITY ) then

   zloc = zloc + 0.5_r8
   k = max(1,int(zloc))

   if(i >= 1 .and. i < wrf%dom(id)%var_size(1,TYPE_W) .and. &
      j >= 1 .and. j < wrf%dom(id)%var_size(2,TYPE_W) .and. &
      k >= 1 .and. k < wrf%dom(id)%var_size(3,TYPE_W)) then

!!$      i1 = get_wrf_index(i,j  ,k,TYPE_W,id)
!!$      i2 = get_wrf_index(i,j+1,k,TYPE_W,id)
      i1 = wrf%dom(id)%dart_ind(i,j  ,k,TYPE_W)
      i2 = wrf%dom(id)%dart_ind(i,j+1,k,TYPE_W)
      fld(1) = dym*( dxm*x(i1) + dx*x(i1+1) ) + dy*( dxm*x(i2) + dx*x(i2+1) )

!!$      i1 = get_wrf_index(i,j  ,k+1,TYPE_W,id)
!!$      i2 = get_wrf_index(i,j+1,k+1,TYPE_W,id)
      i1 = wrf%dom(id)%dart_ind(i,j  ,k+1,TYPE_W)
      i2 = wrf%dom(id)%dart_ind(i,j+1,k+1,TYPE_W)
      fld(2) = dym*( dxm*x(i1) + dx*x(i1+1) ) + dy*( dxm*x(i2) + dx*x(i2+1) )

   else

      fld(:) = missing_r8

   endif

else if( obs_kind == KIND_SPECIFIC_HUMIDITY ) then

!  Convert water vapor mixing ratio to specific humidity:

   if ( wrf%dom(id)%n_moist >= 1) then

      if(.not. vert_is_surface(location) .or. .not. surf_var) then

         if(i >= 1 .and. i < wrf%dom(id)%var_size(1,TYPE_T) .and. &
            j >= 1 .and. j < wrf%dom(id)%var_size(2,TYPE_T) .and. &
            k >= 1 .and. k < wrf%dom(id)%var_size(3,TYPE_T)) then


!!$            i1 = get_wrf_index(i,j ,k,TYPE_QV,id)
!!$            i2 = get_wrf_index(i,j+1,k,TYPE_QV,id)
            i1 = wrf%dom(id)%dart_ind(i,j ,k,TYPE_QV)
            i2 = wrf%dom(id)%dart_ind(i,j+1,k,TYPE_QV)
            a1 = dym*( dxm*x(i1) + dx*x(i1+1) ) + dy*( dxm*x(i2) + dx*x(i2+1) )
            fld(1) = a1 /(1.0_r8 + a1)

!!$            i1 = get_wrf_index(i,j ,k+1,TYPE_QV,id)
!!$            i2 = get_wrf_index(i,j+1,k+1,TYPE_QV,id)
            i1 = wrf%dom(id)%dart_ind(i,j ,k+1,TYPE_QV)
            i2 = wrf%dom(id)%dart_ind(i,j+1,k+1,TYPE_QV)
            a1 = dym*( dxm*x(i1) + dx*x(i1+1) ) + dy*( dxm*x(i2) + dx*x(i2+1) )
            fld(2) = a1 /(1.0_r8 + a1)

         else

            fld(:) = missing_r8

         endif

      else

         if(i >= 1 .and. i < wrf%dom(id)%var_size(1,TYPE_T) .and. &
            j >= 1 .and. j < wrf%dom(id)%var_size(2,TYPE_T) .and. &
            wrf%dom(id)%surf_obs) then

!!$            i1 = get_wrf_index(i,j,1,TYPE_Q2,id)
!!$            i2 = get_wrf_index(i,j+1,1,TYPE_Q2,id)
            i1 = wrf%dom(id)%dart_ind(i,j,1,TYPE_Q2)
            i2 = wrf%dom(id)%dart_ind(i,j+1,1,TYPE_Q2)
            fld(1) = dym*( dxm*x(i1) + dx*x(i1+1) ) + dy*( dxm*x(i2) + dx*x(i2+1) )
            fld(1) = fld(1) / (1.0_r8 + fld(1))

         else

            fld(1) = missing_r8

         endif

      endif

   else

      fld(:) = 0.0_r8

   endif

else if( obs_kind == KIND_VAPOR_MIXING_RATIO ) then

   if ( wrf%dom(id)%n_moist >= 1) then

      if(.not. vert_is_surface(location) .or. .not. surf_var) then

         if(i >= 1 .and. i < wrf%dom(id)%var_size(1,TYPE_T) .and. &
            j >= 1 .and. j < wrf%dom(id)%var_size(2,TYPE_T) .and. &
            k >= 1 .and. k < wrf%dom(id)%var_size(3,TYPE_T)) then

!!$            i1 = get_wrf_index(i,j  ,k,TYPE_QV,id)
!!$            i2 = get_wrf_index(i,j+1,k,TYPE_QV,id)
            i1 = wrf%dom(id)%dart_ind(i,j  ,k,TYPE_QV)
            i2 = wrf%dom(id)%dart_ind(i,j+1,k,TYPE_QV)
            fld(1) = dym*( dxm*x(i1) + dx*x(i1+1) ) + dy*( dxm*x(i2) + dx*x(i2+1) )

!!$            i1 = get_wrf_index(i,j  ,k+1,TYPE_QV,id)
!!$            i2 = get_wrf_index(i,j+1,k+1,TYPE_QV,id)
            i1 = wrf%dom(id)%dart_ind(i,j  ,k+1,TYPE_QV)
            i2 = wrf%dom(id)%dart_ind(i,j+1,k+1,TYPE_QV)
            fld(2) = dym*( dxm*x(i1) + dx*x(i1+1) ) + dy*( dxm*x(i2) + dx*x(i2+1) )

         else

            fld(:) = missing_r8

         endif

      else

         if(i >= 1 .and. i < wrf%dom(id)%var_size(1,TYPE_T) .and. &
            j >= 1 .and. j < wrf%dom(id)%var_size(2,TYPE_T) .and. &
            wrf%dom(id)%surf_obs) then

!!$            i1 = get_wrf_index(i,j,1,TYPE_Q2,id)
!!$            i2 = get_wrf_index(i,j+1,1,TYPE_Q2,id)
            i1 = wrf%dom(id)%dart_ind(i,j,1,TYPE_Q2)
            i2 = wrf%dom(id)%dart_ind(i,j+1,1,TYPE_Q2)
            fld(1) = dym*( dxm*x(i1) + dx*x(i1+1) ) + dy*( dxm*x(i2) + dx*x(i2+1) )

         else

            fld(1) = missing_r8

         endif

      endif

   else

      fld(:) = 0.0_r8

   endif

else if( obs_kind == KIND_RAINWATER_MIXING_RATIO) then

   if(i >= 1 .and. i < wrf%dom(id)%var_size(1,TYPE_T) .and. &
      j >= 1 .and. j < wrf%dom(id)%var_size(2,TYPE_T) .and. &
      k >= 1 .and. k < wrf%dom(id)%var_size(3,TYPE_T)) then

      if( wrf%dom(id)%n_moist >= 3) then

!!$      i1 = get_wrf_index(i,j  ,k,TYPE_QR,id)
!!$      i2 = get_wrf_index(i,j+1,k,TYPE_QR,id)
         i1 = wrf%dom(id)%dart_ind(i,j  ,k,TYPE_QR)
         i2 = wrf%dom(id)%dart_ind(i,j+1,k,TYPE_QR)
         fld(1) = dym*( dxm*x(i1) + dx*x(i1+1) ) + dy*( dxm*x(i2) + dx*x(i2+1) )

!!$      i1 = get_wrf_index(i,j  ,k+1,TYPE_QR,id)
!!$      i2 = get_wrf_index(i,j+1,k+1,TYPE_QR,id)
         i1 = wrf%dom(id)%dart_ind(i,j  ,k+1,TYPE_QR)
         i2 = wrf%dom(id)%dart_ind(i,j+1,k+1,TYPE_QR)
         fld(2) = dym*( dxm*x(i1) + dx*x(i1+1) ) + dy*( dxm*x(i2) + dx*x(i2+1) )

         fld = max(0.0_r8, fld)

      else

         fld(:) = 0.0_r8

      endif


   else

      fld(:) = missing_r8

   endif

else if( obs_kind == KIND_GRAUPEL_MIXING_RATIO) then

   if(i >= 1 .and. i < wrf%dom(id)%var_size(1,TYPE_T) .and. &
      j >= 1 .and. j < wrf%dom(id)%var_size(2,TYPE_T) .and. &
      k >= 1 .and. k < wrf%dom(id)%var_size(3,TYPE_T)) then

      if( wrf%dom(id)%n_moist >= 6) then

         i1 = wrf%dom(id)%dart_ind(i,j  ,k,TYPE_QG)
         i2 = wrf%dom(id)%dart_ind(i,j+1,k,TYPE_QG)
         fld(1) = dym*( dxm*x(i1) + dx*x(i1+1) ) + dy*( dxm*x(i2) + dx*x(i2+1) )

         i1 = wrf%dom(id)%dart_ind(i,j  ,k+1,TYPE_QG)
         i2 = wrf%dom(id)%dart_ind(i,j+1,k+1,TYPE_QG)
         fld(2) = dym*( dxm*x(i1) + dx*x(i1+1) ) + dy*( dxm*x(i2) + dx*x(i2+1) )
         fld = max(0.0_r8, fld)

      else

         fld(:) = 0.0_r8

      endif


   else

      fld(:) = missing_r8

   endif

else if( obs_kind == KIND_SNOW_MIXING_RATIO) then

   if(i >= 1 .and. i < wrf%dom(id)%var_size(1,TYPE_T) .and. &
      j >= 1 .and. j < wrf%dom(id)%var_size(2,TYPE_T) .and. &
      k >= 1 .and. k < wrf%dom(id)%var_size(3,TYPE_T)) then

      if( wrf%dom(id)%n_moist >= 5) then

         i1 = wrf%dom(id)%dart_ind(i,j  ,k,TYPE_QS)
         i2 = wrf%dom(id)%dart_ind(i,j+1,k,TYPE_QS)
         fld(1) = dym*( dxm*x(i1) + dx*x(i1+1) ) + dy*( dxm*x(i2) + dx*x(i2+1) )

         i1 = wrf%dom(id)%dart_ind(i,j  ,k+1,TYPE_QS)
         i2 = wrf%dom(id)%dart_ind(i,j+1,k+1,TYPE_QS)
         fld(2) = dym*( dxm*x(i1) + dx*x(i1+1) ) + dy*( dxm*x(i2) + dx*x(i2+1) )

         fld = max(0.0_r8, fld)

      else

         fld(:) = 0.0_r8

      endif

   else

      fld(:) = missing_r8

   endif

else if( obs_kind == KIND_PRESSURE) then

   if(i >= 1 .and. i < wrf%dom(id)%var_size(1,TYPE_T) .and. &
      j >= 1 .and. j < wrf%dom(id)%var_size(2,TYPE_T) .and. &
      k >= 1 .and. k < wrf%dom(id)%var_size(3,TYPE_T)) then

      pres1 = model_pressure_t(i  ,j  ,k,id,x)
      pres2 = model_pressure_t(i+1,j  ,k,id,x)
      pres3 = model_pressure_t(i  ,j+1,k,id,x)
      pres4 = model_pressure_t(i+1,j+1,k,id,x)
      fld(1) = dym*( dxm*pres1 + dx*pres2 ) + dy*( dxm*pres3 + dx*pres4 )

      pres1 = model_pressure_t(i  ,j  ,k+1,id,x)
      pres2 = model_pressure_t(i+1,j  ,k+1,id,x)
      pres3 = model_pressure_t(i  ,j+1,k+1,id,x)
      pres4 = model_pressure_t(i+1,j+1,k+1,id,x)
      fld(2) = dym*( dxm*pres1 + dx*pres2 ) + dy*( dxm*pres3 + dx*pres4 )

   else

      fld(:) = missing_r8

   endif

else if( obs_kind == KIND_SURFACE_PRESSURE) then

   if(i >= 1 .and. i < wrf%dom(id)%var_size(1,TYPE_T) .and. &
      j >= 1 .and. j < wrf%dom(id)%var_size(2,TYPE_T) .and. &
      wrf%dom(id)%surf_obs) then

!!$      i1 = get_wrf_index(i,j,1,TYPE_PS,id)
!!$      i2 = get_wrf_index(i,j+1,1,TYPE_PS,id)
      i1 = wrf%dom(id)%dart_ind(i,j,1,TYPE_PS)
      i2 = wrf%dom(id)%dart_ind(i,j+1,1,TYPE_PS)
      if(x(i1) /= 0.0_r8 .and. x(i1+1) /= 0.0_r8 .and. &
           x(i2) /= 0.0_r8 .and. x(i2+1) /= 0.0_r8) then

         fld(1) = dym*( dxm*x(i1) + dx*x(i1+1) ) + dy*( dxm*x(i2) + dx*x(i2+1) )

      else

         fld(1) = missing_r8

      endif

   else

      fld(1) = missing_r8

   endif

!!-------------------------------------------------------------------------
!! find vortex center location and with minimum sea level pressure
!! steps: 1) define a search box centered at the obs. location (i,j)
!!        2) compute sea-level pressure (slp)
!!        3) spline-interpolation of slp within the refined search box 
!!           -- by a factor of center_spline_grid_scale (default=10)
!!        4) find the minimum slp and grid point location 
!!           (default=1/10 precision)
!!        5) return the minimum slp and lat/lon of this grid point
!!-------------------------------------------------------------------------
else if( obs_kind == KIND_VORTEX_LAT .or. &
         obs_kind == KIND_VORTEX_LON .or. &
         obs_kind == KIND_VORTEX_PMIN .or. &
         obs_kind == KIND_VORTEX_WMAX) then

   if(i >= 1 .and. i < wrf%dom(id)%var_size(1,TYPE_T)  .and. &
      j >= 1 .and. j < wrf%dom(id)%var_size(2,TYPE_T)  .and. &
      k >= 1 .and. k < wrf%dom(id)%var_size(3,TYPE_T)) then

!!   define a search box bounded by center_track_***
     center_search_half_size = nint(center_search_half_length/wrf%dom(id)%dx)
     center_track_xmin = max(1,i-center_search_half_size)
     center_track_xmax = min(wrf%dom(id)%var_size(1,TYPE_MU),i+center_search_half_size)
     center_track_ymin = max(1,j-center_search_half_size)
     center_track_ymax = min(wrf%dom(id)%var_size(2,TYPE_MU),j+center_search_half_size)
     if(center_track_xmin<1 .or. center_track_xmax>wrf%dom(id)%var_size(1,TYPE_MU) .or. &
        center_track_ymin<1 .or. center_track_ymax>wrf%dom(id)%var_size(2,TYPE_MU) .or. &
        center_track_xmin >= center_track_xmax .or. &
        center_track_ymin >= center_track_ymax) then
          print*,'i,j,center_search_half_length,center_track_xmin(max),center_track_ymin(max)'
          print*,i,j,center_search_half_length,center_track_xmin,center_track_xmax,center_track_ymin,center_track_ymax
         write(errstring,*)'Wrong setup in center_track_nml'
         call error_handler(E_ERR,'model_interpolate', errstring, source, revision, revdate)
     endif 

!!   define spline interpolation box dimensions
     xlen = center_track_xmax - center_track_xmin + 1
     ylen = center_track_ymax - center_track_ymin + 1
     xxlen = (center_track_xmax - center_track_xmin)*center_spline_grid_scale + 1
     yylen = (center_track_ymax - center_track_ymin)*center_spline_grid_scale + 1
     allocate(p1d(wrf%dom(id)%bt), t1d(wrf%dom(id)%bt))
     allocate(qv1d(wrf%dom(id)%bt), z1d(wrf%dom(id)%bt))
     allocate(psea(xlen,ylen))
     allocate(pd(xlen,ylen))
     allocate(pp(xxlen,yylen))
     allocate(x1d(xlen))
     allocate(y1d(ylen))
     allocate(xx1d(xxlen))
     allocate(yy1d(yylen))

!!   compute sea-level pressure
     do i1 = center_track_xmin, center_track_xmax
     do i2 = center_track_ymin, center_track_ymax
        do k2 = 1,wrf%dom(id)%var_size(3,TYPE_T)
           p1d(k2) = model_pressure_t(i1,i2,k2,id,x)
           t1d(k2) = ts0 + x(wrf%dom(id)%dart_ind(i1,i2,k2,TYPE_T))
           qv1d(k2)= x(wrf%dom(id)%dart_ind(i1,i2,k2,TYPE_QV))
           z1d(k2) = ( x(wrf%dom(id)%dart_ind(i1,i2,k2,TYPE_GZ))+wrf%dom(id)%phb(i1,i2,k2) + &
                      x(wrf%dom(id)%dart_ind(i1,i2,k2+1,TYPE_GZ))+wrf%dom(id)%phb(i1,i2,k2+1) &
                    )*0.5_r8/gravity
        enddo
        call compute_seaprs(wrf%dom(id)%bt, z1d, t1d, p1d, qv1d, &
                          psea(i1-center_track_xmin+1,i2-center_track_ymin+1),debug)
     enddo
     enddo

!!   spline-interpolation
     do i1 = 1,xlen
        x1d(i1) = (i1-1)+center_track_xmin
     enddo
     do i2 = 1,ylen
        y1d(i2) = (i2-1)+center_track_ymin
     enddo
     do ii1 = 1,xxlen
        xx1d(ii1) = center_track_xmin+real(ii1-1,r8)*1_r8/real(center_spline_grid_scale,r8)
     enddo
     do ii2 = 1,yylen
        yy1d(ii2) = center_track_ymin+real(ii2-1,r8)*1_r8/real(center_spline_grid_scale,r8)
     enddo

     call splie2(x1d,y1d,psea,xlen,ylen,pd)

     pres1 = 1.e20
     cxmin = -1
     cymin = -1
     do ii1=1,xxlen
     do ii2=1,yylen
        call splin2(x1d,y1d,psea,pd,xlen,ylen,xx1d(ii1),yy1d(ii2),pp(ii1,ii2))
        if (pres1 .gt. pp(ii1,ii2)) then
           pres1=pp(ii1,ii2)
           cxmin=xx1d(ii1)
           cymin=yy1d(ii2)
        endif
     enddo
     enddo

!!   if too close to the edge of the box, reset to observed center
     if( cxmin-xx1d(1) < 1_r8 .or. xx1d(xxlen)-cxmin < 1_r8 .or.  &
         cymin-yy1d(1) < 1_r8 .or. yy1d(yylen)-cymin < 1_r8 ) then
       cxmin = xloc
       cymin = yloc
       call splin2(x1d,y1d,psea,pd,xlen,ylen,cxmin,cymin,pres1)
     endif

     call ij_to_latlon(wrf%dom(id)%proj, cxmin, cymin, clat, clon)

     if( obs_kind == KIND_VORTEX_LAT ) then
        fld(1) = clat
     else if( obs_kind == KIND_VORTEX_LON ) then
        fld(1) = clon
     else if( obs_kind == KIND_VORTEX_PMIN ) then
        fld(1) = pres1
     else if( obs_kind == KIND_VORTEX_WMAX ) then
        fld(1) = missing_r8
        ! not implemented yet
     endif

     deallocate(p1d, t1d, qv1d, z1d)
     deallocate(psea,pd,pp,x1d,y1d,xx1d,yy1d)


   else

      fld(1) = missing_r8

   endif

elseif (obs_kind < 0) then

   ! identity observation -> -(obs_kind)=DART state vector index
   ! obtain state value directly from index
   fld(1) = x(-1*obs_kind)

else

   write(errstring,*)'Obs kind not recognized for following kind: ',obs_kind
   call error_handler(E_ERR,'model_interpolate', errstring, &
        source, revision, revdate)

end if

! Do vertical interpolation (only for non-surface, non-indetity obs)
if ((.not. vert_is_surface(location)).and.(obs_kind > 0).and. .not. surf_var) then

   call toGrid(zloc, k, dz, dzm)

   if( k >= 1 ) then
      obs_val = dzm*fld(1) + dz*fld(2)
   elseif( k == 0 ) then
      obs_val = fld(1) - (fld(2)-fld(1))*dzm
   else
      obs_val = missing_r8
   endif

else

   obs_val = fld(1)

endif

if(obs_val == missing_r8) istatus = 1

if(debug) print*,' interpolated value= ',obs_val

deallocate(v_h, v_p)


end subroutine model_interpolate


!#######################################################################


subroutine vert_interpolate(x, location, obs_kind, istatus)

! This subroutine converts a given ob/state vertical coordinate to
! the vertical coordinate type requested through the model_mod namelist.

! Notes: (1) obs_kind is only necessary to check whether the ob
!            is an identity ob.
!        (2) This subroutine can convert both obs' and state points'
!            vertical coordinates. Remember that state points get
!            their DART location information from get_state_meta_data
!            which is called by filter_assim during the assimilation
!            process.
!        (3) x is the relevant DART state vector for carrying out
!            interpolations necessary for the vertical coordinate
!            transformations. As the vertical coordinate is only used
!            in distance computations, this is actually the "expected"
!            vertical coordinate, so that computed distance is the
!            "expected" distance. Thus, under normal circumstances,
!            x that is supplied to vert_interpolate should be the
!            ensemble mean. Nevertheless, the subroutine has the
!            functionality to operate on any DART state vector that
!            is supplied to it.

real(r8),            intent(in)    :: x(:)
integer,             intent(in)    :: obs_kind
type(location_type), intent(inout) :: location
integer,             intent(out)   :: istatus

real(r8)            :: xloc, yloc, zloc, xyz_loc(3), zvert
integer             :: id, i, j, k
real(r8)            :: dx,dy,dz,dxm,dym,dzm

character(len=129) :: errstring

real(r8), allocatable, dimension(:) :: v_h, v_p

! local vars, used in calculating pressure and height
real(r8)            :: pres1, pres2, pres3, pres4, pres
real(r8)            :: presa, presb
real(r8)            :: hgt1, hgt2, hgt3, hgt4, hgta, hgtb


istatus = 0

! first off, check if ob is identity ob
if (obs_kind < 0) then
   call get_state_meta_data(obs_kind,location)
   return
endif

xyz_loc = get_location(location)

! first obtain domain id, and mass points (i,j)
call get_domain_info(xyz_loc(1),xyz_loc(2),id,xloc,yloc)

if (id==0) then
   ! Note: need to reset location using the namelist variable directly because
   ! wrf%dom(id)%vert_coord is not defined for id=0
   location = set_location(xyz_loc(1),xyz_loc(2),missing_r8,vert_localization_coord)
   istatus = 1
   return
endif

allocate(v_h(0:wrf%dom(id)%bt), v_p(0:wrf%dom(id)%bt))

! get integer (west/south) grid point and distances to neighboring grid points
! distances are used as weights to carry out horizontal interpolations
call toGrid(xloc,i,dx,dxm)
call toGrid(yloc,j,dy,dym)

! Determine corresponding model level for obs location
! This depends on the obs vertical coordinate
! Obs vertical coordinate will also be converted to the desired
! vertical coordinate as specified by the namelist variable
! "vert_localization_coord" (stored in wrf structure pointer "vert_coord")
if(vert_is_level(location)) then
   ! If obs is by model level: get neighboring mass level indices
   ! and compute weights to zloc
   zloc = xyz_loc(3)
   ! convert obs vert coordinate to desired coordinate type
   if (wrf%dom(id)%vert_coord == VERTISPRESSURE) then
      call toGrid(zloc,k,dz,dzm)
      if ( ((k<1) .or. (k>=wrf%dom(id)%bt)) .or. &
           ((j<1) .or. (j>=wrf%dom(id)%sn)) .or. &
           ((i<1) .or. (i>=wrf%dom(id)%we)) ) then
         zloc  = missing_r8
         zvert = missing_r8
      else
         ! need to compute pressure at all neighboring mass points
         ! and interpolate
         presa = model_pressure_t(i  ,j  ,k  ,id,x)
         presb = model_pressure_t(i  ,j  ,k+1,id,x)
         pres1 = dzm*presa + dz*presb
         presa = model_pressure_t(i+1,j  ,k  ,id,x)
         presb = model_pressure_t(i+1,j  ,k+1,id,x)
         pres2 = dzm*presa + dz*presb
         presa = model_pressure_t(i  ,j+1,k  ,id,x)
         presb = model_pressure_t(i  ,j+1,k+1,id,x)
         pres3 = dzm*presa + dz*presb
         presa = model_pressure_t(i+1,j+1,k  ,id,x)
         presb = model_pressure_t(i+1,j+1,k+1,id,x)
         pres4 = dzm*presa + dz*presb
         zvert = dym*( dxm*pres1 + dx*pres2 ) + dy*( dxm*pres3 + dx*pres4 )
      endif
   elseif (wrf%dom(id)%vert_coord == VERTISHEIGHT) then
      ! need to add half a grid to get to staggered vertical coordinate
      call toGrid(zloc+0.5,k,dz,dzm)
      if ( ((k<1) .or. (k>=wrf%dom(id)%bts)) .or. &
           ((j<1) .or. (j>=wrf%dom(id)%sn )) .or. &
           ((i<1) .or. (i>=wrf%dom(id)%we )) ) then
         zloc  = missing_r8
         zvert = missing_r8
      else
         ! need to compute height at all neighboring vertically staggered points
         ! and interpolate
         hgta = model_height_w(i  ,j  ,k  ,id,x)
         hgtb = model_height_w(i  ,j  ,k+1,id,x)
         hgt1 = dzm*hgta + dz*hgtb
         hgta = model_height_w(i+1,j  ,k  ,id,x)
         hgtb = model_height_w(i+1,j  ,k+1,id,x)
         hgt2 = dzm*hgta + dz*hgtb
         hgta = model_height_w(i  ,j+1,k  ,id,x)
         hgtb = model_height_w(i  ,j+1,k+1,id,x)
         hgt3 = dzm*hgta + dz*hgtb
         hgta = model_height_w(i+1,j+1,k  ,id,x)
         hgtb = model_height_w(i+1,j+1,k+1,id,x)
         hgt4 = dzm*hgta + dz*hgtb
         zvert = dym*( dxm*hgt1 + dx*hgt2 ) + dy*( dxm*hgt3 + dx*hgt4 )
      endif
   else
      if ((zloc<1.0).or.(zloc>=real(wrf%dom(id)%bt))) then
         zloc  = missing_r8
         zvert = missing_r8
      else
         zvert = zloc
      endif
   endif

elseif(vert_is_pressure(location)) then
   ! If obs is by pressure: get corresponding mass level zk,
   ! then get neighboring mass level indices
   ! and compute weights to zloc
   ! get model pressure profile
   call get_model_pressure_profile(i,j,dx,dy,dxm,dym,wrf%dom(id)%bt,x,id,v_p)
   ! get pressure vertical co-ordinate
   call pres_to_zk(xyz_loc(3), v_p, wrf%dom(id)%bt,zloc)
   ! convert obs vert coordinate to desired coordinate type
   if (zloc==missing_r8) then
      zvert = missing_r8
   else
      if (wrf%dom(id)%vert_coord == VERTISLEVEL) then
         zvert = zloc
      elseif (wrf%dom(id)%vert_coord == VERTISHEIGHT) then
         ! adding 0.5 to get to the staggered vertical grid
         ! because height is on staggered vertical grid
         call toGrid(zloc+0.5,k,dz,dzm)
         if ( ((k<1) .or. (k>=wrf%dom(id)%bts)) .or. &
              ((j<1) .or. (j>=wrf%dom(id)%sn )) .or. &
              ((i<1) .or. (i>=wrf%dom(id)%we )) ) then
            zloc  = missing_r8
            zvert = missing_r8
         else
            ! need to compute height at all neighboring vertically staggered points
            ! and interpolate
            hgta = model_height_w(i  ,j  ,k  ,id,x)
            hgtb = model_height_w(i  ,j  ,k+1,id,x)
            hgt1 = dzm*hgta + dz*hgtb
            hgta = model_height_w(i+1,j  ,k  ,id,x)
            hgtb = model_height_w(i+1,j  ,k+1,id,x)
            hgt2 = dzm*hgta + dz*hgtb
            hgta = model_height_w(i  ,j+1,k  ,id,x)
            hgtb = model_height_w(i  ,j+1,k+1,id,x)
            hgt3 = dzm*hgta + dz*hgtb
            hgta = model_height_w(i+1,j+1,k  ,id,x)
            hgtb = model_height_w(i+1,j+1,k+1,id,x)
            hgt4 = dzm*hgta + dz*hgtb
            zvert = dym*( dxm*hgt1 + dx*hgt2 ) + dy*( dxm*hgt3 + dx*hgt4 )
         endif
      else
         ! take pressure directly
         zvert  = xyz_loc(3)
      endif
   endif

elseif(vert_is_height(location)) then
   ! If obs is by height: get corresponding mass level zk,
   ! then get neighboring mass level indices
   ! and compute weights to zloc
   ! get model height profile
   call get_model_height_profile(i,j,dx,dy,dxm,dym,wrf%dom(id)%bt,x,id,v_h)
   ! get height vertical co-ordinate
   call height_to_zk(xyz_loc(3), v_h, wrf%dom(id)%bt,zloc)
   ! convert obs vert coordinate to desired coordinate type
   if (zloc==missing_r8) then
      zvert = missing_r8
   else
      if (wrf%dom(id)%vert_coord == VERTISLEVEL) then
         zvert = zloc
      elseif (wrf%dom(id)%vert_coord == VERTISPRESSURE) then
         call toGrid(zloc,k,dz,dzm)
         if ( ((k<1) .or. (k>=wrf%dom(id)%bt)) .or. &
              ((j<1) .or. (j>=wrf%dom(id)%sn)) .or. &
              ((i<1) .or. (i>=wrf%dom(id)%we)) ) then
            zloc  = missing_r8
            zvert = missing_r8
         else
            ! need to compute pressure at all neighboring mass points
            ! and interpolate
            presa = model_pressure_t(i  ,j  ,k  ,id,x)
            presb = model_pressure_t(i  ,j  ,k+1,id,x)
            pres1 = dzm*presa + dz*presb
            presa = model_pressure_t(i+1,j  ,k  ,id,x)
            presb = model_pressure_t(i+1,j  ,k+1,id,x)
            pres2 = dzm*presa + dz*presb
            presa = model_pressure_t(i  ,j+1,k  ,id,x)
            presb = model_pressure_t(i  ,j+1,k+1,id,x)
            pres3 = dzm*presa + dz*presb
            presa = model_pressure_t(i+1,j+1,k  ,id,x)
            presb = model_pressure_t(i+1,j+1,k+1,id,x)
            pres4 = dzm*presa + dz*presb
            zvert = dym*( dxm*pres1 + dx*pres2 ) + dy*( dxm*pres3 + dx*pres4 )
         endif
      else
         ! take height directly
         zvert  = xyz_loc(3)
      endif
   endif

elseif(vert_is_surface(location) .or. surf_var) then
   zloc = 1.0_r8
   ! convert obs vert coordinate to desired coordinate type
   if (wrf%dom(id)%vert_coord == VERTISLEVEL) then
      zvert = zloc
   elseif (wrf%dom(id)%vert_coord == VERTISPRESSURE) then
      ! need to compute surface pressure at all neighboring mass points
      ! and interpolate
      pres1 = model_pressure_s(i  ,j  ,id,x)
      pres2 = model_pressure_s(i+1,j  ,id,x)
      pres3 = model_pressure_s(i  ,j+1,id,x)
      pres4 = model_pressure_s(i+1,j+1,id,x)
      zvert = dym*( dxm*pres1 + dx*pres2 ) + dy*( dxm*pres3 + dx*pres4 )
   else
      ! a surface ob is assumed to have height as vertical coordinate...
      ! this may need to be revised if this is not always true (in which
      ! case, just need to uncomment below lines to get terrain height
      ! from model)
      zvert = xyz_loc(3)
      !! directly interpolate terrain height at neighboring mass points
      !zvert = dym*( dxm*wrf%dom(id)%hgt(i,j) + &
      !             dx*wrf%dom(id)%hgt(i+1,j) ) + &
      !        dy*( dxm*wrf%dom(id)%hgt(i,j+1) + &
      !             dx*wrf%dom(id)%hgt(i+1,j+1) )
   endif

elseif(vert_is_undef(location)) then
   zloc  = missing_r8
   zvert = missing_r8

else
   write(errstring,*) 'Vertical coordinate not recognized: ',nint(query_location(location,'which_vert'))
   call error_handler(E_ERR,'vert_interpolate', errstring, &
        source, revision, revdate)

endif

deallocate(v_h, v_p)

if(zvert == missing_r8) istatus = 1

! Reset location   
location = set_location(xyz_loc(1),xyz_loc(2),zvert,wrf%dom(id)%vert_coord)

end subroutine vert_interpolate


!#######################################################################


function get_wrf_index( i,j,k,var_type,id )

integer, intent(in) :: i,j,k,var_type,id

integer :: get_wrf_index
integer :: in

character(len=129) :: errstring

do in = 1, wrf%dom(id)%number_of_wrf_variables
   if(var_type == wrf%dom(id)%var_type(in) ) then
      exit
   endif
enddo

if(i >= 1 .and. i <= wrf%dom(id)%var_size(1,in) .and. &
   j >= 1 .and. j <= wrf%dom(id)%var_size(2,in) .and. &
   k >= 1 .and. k <= wrf%dom(id)%var_size(3,in)) then

   get_wrf_index = wrf%dom(id)%dart_ind(i,j,k,var_type)

!!$   get_wrf_index = wrf%dom(id)%var_index(1,in)-1 +   &
!!$        i + wrf%dom(id)%var_size(1,in)*((j-1) + &
!!$        wrf%dom(id)%var_size(2,in)*(k-1))

else

  write(errstring,*)'Indices ',i,j,k,' exceed grid dimensions: ', &
       wrf%dom(id)%var_size(1,in), &
       wrf%dom(id)%var_size(2,in),wrf%dom(id)%var_size(3,in)
  call error_handler(E_ERR,'get_wrf_index', errstring, &
       source, revision, revdate)

endif

end function get_wrf_index


!***********************************************************************


subroutine get_wrf_horizontal_location( i, j, var_type, id, long, lat )

integer,  intent(in)  :: i,j,var_type, id
real(r8), intent(out) :: long, lat

! find lat and long, must
! correct for possible u or v staggering in x, y

if (var_type == type_u) then

   if (i == 1) then
      long = wrf%dom(id)%longitude(1,j) - &
           0.5_r8*(wrf%dom(id)%longitude(2,j)-wrf%dom(id)%longitude(1,j))
      if ( abs(wrf%dom(id)%longitude(2,j) - wrf%dom(id)%longitude(1,j)) > 180.0_r8 ) then
         long = long - 180.0_r8
      endif
      lat = wrf%dom(id)%latitude(1,j) - &
           0.5_r8*(wrf%dom(id)%latitude(2,j)-wrf%dom(id)%latitude(1,j))
   else if (i == wrf%dom(id)%wes) then
      long = wrf%dom(id)%longitude(i-1,j) + &
           0.5_r8*(wrf%dom(id)%longitude(i-1,j)-wrf%dom(id)%longitude(i-2,j))
      if ( abs(wrf%dom(id)%longitude(i-1,j) - wrf%dom(id)%longitude(i-2,j)) > 180.0_r8 ) then
         long = long - 180.0_r8
      endif
      lat = wrf%dom(id)%latitude(i-1,j) + &
           0.5_r8*(wrf%dom(id)%latitude(i-1,j)-wrf%dom(id)%latitude(i-2,j))
   else
      long = 0.5_r8*(wrf%dom(id)%longitude(i,j)+wrf%dom(id)%longitude(i-1,j))
      if ( abs(wrf%dom(id)%longitude(i,j) - wrf%dom(id)%longitude(i-1,j)) > 180.0_r8 ) then
         long = long - 180.0_r8
      endif
      lat = 0.5_r8*(wrf%dom(id)%latitude(i,j) +wrf%dom(id)%latitude(i-1,j))
   end if

elseif (var_type == type_v) then

   if (j == 1) then
      long = wrf%dom(id)%longitude(i,1) - &
           0.5_r8*(wrf%dom(id)%longitude(i,2)-wrf%dom(id)%longitude(i,1))
      if ( abs(wrf%dom(id)%longitude(i,2) - wrf%dom(id)%longitude(i,1)) > 180.0_r8 ) then
         long = long - 180.0_r8
      endif
      lat = wrf%dom(id)%latitude(i,1) - &
           0.5_r8*(wrf%dom(id)%latitude(i,2)-wrf%dom(id)%latitude(i,1))
   else if (j == wrf%dom(id)%sns) then
      long = wrf%dom(id)%longitude(i,j-1) + &
           0.5_r8*(wrf%dom(id)%longitude(i,j-1)-wrf%dom(id)%longitude(i,j-2))
      if ( abs(wrf%dom(id)%longitude(i,j-1) - wrf%dom(id)%longitude(i,j-2)) > 180.0_r8 ) then
         long = long - 180.0_r8
      endif
      lat = wrf%dom(id)%latitude(i,j-1) + &
           0.5_r8*(wrf%dom(id)%latitude(i,j-1)-wrf%dom(id)%latitude(i,j-2))
   else
      long = 0.5_r8*(wrf%dom(id)%longitude(i,j)+wrf%dom(id)%longitude(i,j-1))
      if ( abs(wrf%dom(id)%longitude(i,j) - wrf%dom(id)%longitude(i,j-1)) > 180.0_r8 ) then
         long = long - 180.0_r8
      endif
      lat  = 0.5_r8*(wrf%dom(id)%latitude(i,j) +wrf%dom(id)%latitude(i,j-1))

   end if

else

   long = wrf%dom(id)%longitude(i,j)
   lat  = wrf%dom(id)%latitude(i,j)

end if

do while (long <   0.0_r8)
   long = long + 360.0_r8
end do
do while (long > 360.0_r8)
   long = long - 360.0_r8
end do

end subroutine get_wrf_horizontal_location



!***********************************************************************


function nc_write_model_atts( ncFileID ) result (ierr)
!-----------------------------------------------------------------
! Writes the model-specific attributes to a netCDF file
! A. Caya May 7 2003
! T. Hoar Mar 8 2004 writes prognostic flavor

integer, intent(in)  :: ncFileID      ! netCDF file identifier
integer              :: ierr          ! return value of function

!-----------------------------------------------------------------

integer :: nDimensions, nVariables, nAttributes, unlimitedDimID
integer :: StateVarDimID, StateVarVarID, StateVarID, TimeDimID

integer, dimension(num_domains) :: weDimID, weStagDimID, snDimID, snStagDimID, &
     btDimID, btStagDimID, slSDimID

integer :: MemberDimID, DomDimID
integer :: DXVarID, DYVarID, TRUELAT1VarID, TRUELAT2VarID
integer :: CEN_LATVarID, CEN_LONVarID, MAP_PROJVarID

integer, dimension(num_domains) :: DNVarID, ZNUVarID, DNWVarID, phbVarID, &
     MubVarID, LonVarID, LatVarID, ilevVarID, XlandVarID, &
     MapFacMVarID, MapFacUVarID, MapFacVVarID, hgtVarID

integer :: var_id
integer :: i, id

character(len=129) :: errstring

character(len=8)      :: crdate      ! needed by F90 DATE_AND_TIME intrinsic
character(len=10)     :: crtime      ! needed by F90 DATE_AND_TIME intrinsic
character(len=5)      :: crzone      ! needed by F90 DATE_AND_TIME intrinsic
integer, dimension(8) :: values      ! needed by F90 DATE_AND_TIME intrinsic
character(len=NF90_MAX_NAME) :: str1

character (len=1)     :: idom

!-----------------------------------------------------------------

ierr = 0     ! assume normal termination

!-----------------------------------------------------------------
! make sure ncFileID refers to an open netCDF file, 
! and then put into define mode.
!-----------------------------------------------------------------

call check(nf90_Inquire(ncFileID, nDimensions, nVariables, nAttributes, unlimitedDimID))
call check(nf90_Redef(ncFileID))

!-----------------------------------------------------------------
! We need the dimension ID for the number of copies 
!-----------------------------------------------------------------

call check(nf90_inq_dimid(ncid=ncFileID, name="copy", dimid=MemberDimID))
call check(nf90_inq_dimid(ncid=ncFileID, name="time", dimid=  TimeDimID))

if ( TimeDimID /= unlimitedDimId ) then
   write(errstring,*)'Time Dimension ID ',TimeDimID, &
        ' must match Unlimited Dimension ID ',unlimitedDimID
   call error_handler(E_ERR,'nc_write_model_atts', errstring, source, revision, revdate)
endif

!-----------------------------------------------------------------
! Define the model size, state variable dimension ... whatever ...
!-----------------------------------------------------------------
call check(nf90_def_dim(ncid=ncFileID, name="StateVariable", &
                        len=wrf%model_size, dimid = StateVarDimID))

!-----------------------------------------------------------------
! Write Global Attributes 
!-----------------------------------------------------------------
call DATE_AND_TIME(crdate,crtime,crzone,values)
write(str1,'(''YYYY MM DD HH MM SS = '',i4,5(1x,i2.2))') &
                  values(1), values(2), values(3), values(5), values(6), values(7)

call check(nf90_put_att(ncFileID, NF90_GLOBAL, "creation_date",str1))
call check(nf90_put_att(ncFileID, NF90_GLOBAL, "model","WRF"))
call check(nf90_put_att(ncFileID, NF90_GLOBAL, "model_source",source))
call check(nf90_put_att(ncFileID, NF90_GLOBAL, "model_revision",revision))
call check(nf90_put_att(ncFileID, NF90_GLOBAL, "model_revdate",revdate))

! how about namelist input? might be nice to save ...

!-----------------------------------------------------------------
! Define the dimensions IDs
!-----------------------------------------------------------------

call check(nf90_def_dim(ncid=ncFileID, name="domain",           &
          len = num_domains,  dimid = DomDimID))

do id=1,num_domains
   write( idom , '(I1)') id
   call check(nf90_def_dim(ncid=ncFileID, name="west_east_d0"//idom,        &
        len = wrf%dom(id)%we,  dimid = weDimID(id)))
   call check(nf90_def_dim(ncid=ncFileID, name="west_east_stag_d0"//idom,   &
        len = wrf%dom(id)%wes, dimid = weStagDimID(id)))
   call check(nf90_def_dim(ncid=ncFileID, name="south_north_d0"//idom,      &
        len = wrf%dom(id)%sn,  dimid = snDimID(id)))
   call check(nf90_def_dim(ncid=ncFileID, name="south_north_stag_d0"//idom, &
        len = wrf%dom(id)%sns, dimid = snStagDimID(id)))
   call check(nf90_def_dim(ncid=ncFileID, name="bottom_top_d0"//idom,       &
        len = wrf%dom(id)%bt,  dimid = btDimID(id)))
   call check(nf90_def_dim(ncid=ncFileID, name="bottom_top_stag_d0"//idom,  &
        len = wrf%dom(id)%bts, dimid = btStagDimID(id)))
   call check(nf90_def_dim(ncid=ncFileID, name="soil_layers_stag_d0"//idom,  &
        len = wrf%dom(id)%sls, dimid = slSDimID(id)))
enddo

!-----------------------------------------------------------------
! Create the (empty) Variables and the Attributes
!-----------------------------------------------------------------

!-----------------------------------------------------------------
! Create the (empty) static variables and their attributes
! Commented block is from wrfinput
!-----------------------------------------------------------------

call check(nf90_def_var(ncFileID, name="DX", xtype=nf90_real, &
     dimids= DomDimID, varid=DXVarID) )
call check(nf90_put_att(ncFileID, DXVarID, "long_name", &
     "X HORIZONTAL RESOLUTION"))
call check(nf90_put_att(ncFileID, DXVarID, "units", "m"))

call check(nf90_def_var(ncFileID, name="DY", xtype=nf90_real, &
     dimids= DomDimID, varid=DYVarID) )
call check(nf90_put_att(ncFileID, DYVarID, "long_name", &
     "Y HORIZONTAL RESOLUTION"))
call check(nf90_put_att(ncFileID, DYVarID, "units", "m"))

call check(nf90_def_var(ncFileID, name="TRUELAT1", xtype=nf90_real, &
     dimids= DomDimID, varid=TRUELAT1VarID) )
call check(nf90_put_att(ncFileID, TRUELAT1VarID, "long_name", &
     "first standard parallel"))
call check(nf90_put_att(ncFileID, TRUELAT1VarID, "units", &
     "degrees, negative is south"))

call check(nf90_def_var(ncFileID, name="TRUELAT2", xtype=nf90_real, &
     dimids= DomDimID, varid=TRUELAT2VarID) )
call check(nf90_put_att(ncFileID, TRUELAT2VarID, "long_name", &
     "second standard parallel"))
call check(nf90_put_att(ncFileID, TRUELAT2VarID, "units", &
     "degrees, negative is south"))

call check(nf90_def_var(ncFileID, name="CEN_LAT", xtype=nf90_real, &
     dimids= DomDimID, varid=CEN_LATVarID) )
call check(nf90_put_att(ncFileID, CEN_LATVarID, "long_name", &
     "center latitude"))
call check(nf90_put_att(ncFileID, CEN_LATVarID, "units", &
     "degrees, negative is south"))

call check(nf90_def_var(ncFileID, name="CEN_LON", xtype=nf90_real, &
     dimids= DomDimID, varid=CEN_LONVarID) )
call check(nf90_put_att(ncFileID, CEN_LONVarID, "long_name", &
     "central longitude"))
call check(nf90_put_att(ncFileID, CEN_LONVarID, "units", &
     "degrees, negative is west"))

call check(nf90_def_var(ncFileID, name="MAP_PROJ", xtype=nf90_real, &
     dimids= DomDimID, varid=MAP_PROJVarID) )
call check(nf90_put_att(ncFileID, MAP_PROJVarID, "long_name", &
     "domain map projection"))
call check(nf90_put_att(ncFileID, MAP_PROJVarID, "units", &
     "0=none, 1=Lambert, 2=polar, 3=Mercator"))

do id=1,num_domains
   write( idom , '(I1)') id

   call check(nf90_def_var(ncFileID, name="DN_d0"//idom, xtype=nf90_real, &
        dimids= btDimID(id), varid=DNVarID(id)) )
   call check(nf90_put_att(ncFileID, DNVarID(id), "long_name", &
        "dn values on half (mass) levels"))
   call check(nf90_put_att(ncFileID, DNVarID(id), "units", "dimensionless"))

   call check(nf90_def_var(ncFileID, name="ZNU_d0"//idom, xtype=nf90_real, &
        dimids= btDimID(id), varid=ZNUVarID(id)) )
   call check(nf90_put_att(ncFileID, ZNUVarID(id), "long_name", &
        "eta values on half (mass) levels"))
   call check(nf90_put_att(ncFileID, ZNUVarID(id), "units", "dimensionless"))

   call check(nf90_def_var(ncFileID, name="DNW_d0"//idom, xtype=nf90_real, &
        dimids= btDimID(id), varid=DNWVarID(id)) )
   call check(nf90_put_att(ncFileID, DNWVarID(id), "long_name", &
        "dn values on full (w) levels"))
   call check(nf90_put_att(ncFileID, DNWVarID(id), "units", "dimensionless"))

!
!    float MUB(Time, south_north, west_east) ;
!            MUB:FieldType = 104 ;
!            MUB:MemoryOrder = "XY " ;
!            MUB:stagger = "" ;
   call check(nf90_def_var(ncFileID, name="MUB_d0"//idom, xtype=nf90_real, &
        dimids= (/ weDimID(id), snDimID(id) /), varid=MubVarID(id)) )
   call check(nf90_put_att(ncFileID, MubVarID(id), "long_name", &
        "base state dry air mass in column"))
   call check(nf90_put_att(ncFileID, MubVarID(id), "units", "Pa"))

! Longitudes
!      float XLONG(Time, south_north, west_east) ;
!         XLONG:FieldType = 104 ;
!         XLONG:MemoryOrder = "XY " ;
!         XLONG:stagger = "" ;
   call check(nf90_def_var(ncFileID, name="XLON_d0"//idom, xtype=nf90_real, &
        dimids= (/ weDimID(id), snDimID(id) /), varid=LonVarID(id)) )
   call check(nf90_put_att(ncFileID, LonVarID(id), "long_name", "longitude"))
   call check(nf90_put_att(ncFileID, LonVarID(id), "units", "degrees_east"))
   call check(nf90_put_att(ncFileID, LonVarID(id), "valid_range", &
        (/ -180.0_r8, 180.0_r8 /)))
   call check(nf90_put_att(ncFileID, LonVarID(id), "description", &
        "LONGITUDE, WEST IS NEGATIVE"))

! Latitudes
!      float XLAT(Time, south_north, west_east) ;
!         XLAT:FieldType = 104 ;
!         XLAT:MemoryOrder = "XY " ;
!         XLAT:stagger = "" ;
   call check(nf90_def_var(ncFileID, name="XLAT_d0"//idom, xtype=nf90_real, &
        dimids=(/ weDimID(id), snDimID(id) /), varid=LatVarID(id)) ) 
   call check(nf90_put_att(ncFileID, LatVarID(id), "long_name", "latitude"))
   call check(nf90_put_att(ncFileID, LatVarID(id), "units", "degrees_north"))
   call check(nf90_put_att(ncFileID, LatVarID(id), "valid_range", &
        (/ -90.0_r8, 90.0_r8 /)))
   call check(nf90_put_att(ncFileID, LatVarID(id), "description", &
        "LATITUDE, SOUTH IS NEGATIVE"))

! grid levels
   call check(nf90_def_var(ncFileID, name="level_d0"//idom, xtype=nf90_short, &
        dimids=btDimID(id), varid=ilevVarID(id)) )
   call check(nf90_put_att(ncFileID, ilevVarID(id), "long_name", &
        "placeholder for level"))
   call check(nf90_put_att(ncFileID, ilevVarID(id), "units", &
        "at this point, indexical"))

! Land Mask
!    float XLAND(Time, south_north, west_east) ;
!            XLAND:FieldType = 104 ;
!            XLAND:MemoryOrder = "XY " ;
!            XLAND:units = "NA" ;
!            XLAND:stagger = "" ;
   call check(nf90_def_var(ncFileID, name="XLAND_d0"//idom, xtype=nf90_short, &
        dimids= (/ weDimID(id), snDimID(id) /), varid=XlandVarID(id)) )
   call check(nf90_put_att(ncFileID, XlandVarID(id), "long_name", "land mask"))
   call check(nf90_put_att(ncFileID, XlandVarID(id), "units", "NA"))
   call check(nf90_put_att(ncFileID, XlandVarID(id), "valid_range", (/ 1, 2 /)))
   call check(nf90_put_att(ncFileID, XlandVarID(id), "description", &
        "1 = LAND, 2 = WATER"))

! Map Scale Factor on m-grid
!    float MAPFAC_M(Time, south_north, west_east) ;
!            MAPFAC_M:FieldType = 104 ;
!            MAPFAC_M:MemoryOrder = "XY " ;
!            MAPFAC_M:stagger = "" ;
   call check(nf90_def_var(ncFileID, name="MAPFAC_M_d0"//idom, xtype=nf90_real, &
        dimids= (/ weDimID(id), snDimID(id) /), varid=MapFacMVarID(id)) )
   call check(nf90_put_att(ncFileID, MapFacMVarID(id), "long_name", &
        "Map scale factor on mass grid"))
   call check(nf90_put_att(ncFileID, MapFacMVarID(id), "units", "dimensionless"))

! Map Scale Factor on u-grid
!    float MAPFAC_U(Time, south_north, west_east_stag) ;
!            MAPFAC_U:FieldType = 104 ;
!            MAPFAC_U:MemoryOrder = "XY " ;
!            MAPFAC_U:stagger = "X" ;
   call check(nf90_def_var(ncFileID, name="MAPFAC_U_d0"//idom, xtype=nf90_real, &
        dimids= (/ weStagDimID(id), snDimID(id) /), varid=MapFacUVarID(id)) )
   call check(nf90_put_att(ncFileID, MapFacUVarID(id), "long_name", &
        "Map scale factor on u-grid"))
   call check(nf90_put_att(ncFileID, MapFacUVarID(id), "units", "dimensionless"))

! Map Scale Factor on v-grid
!    float MAPFAC_V(Time, south_north_stag, west_east) ;
!            MAPFAC_V:FieldType = 104 ;
!            MAPFAC_V:MemoryOrder = "XY " ;
!            MAPFAC_V:stagger = "Y" ;
   call check(nf90_def_var(ncFileID, name="MAPFAC_V_d0"//idom, xtype=nf90_real, &
        dimids= (/ weDimID(id), snStagDimID(id) /), varid=MapFacVVarID(id)) )
   call check(nf90_put_att(ncFileID, MapFacVVarID(id), "long_name", &
        "Map scale factor on v-grid"))
   call check(nf90_put_att(ncFileID, MapFacVVarID(id), "units", "dimensionless"))

! PHB
!    float PHB(Time, bottom_top_stag, south_north, west_east) ;
!            PHB:FieldType = 104 ;
!            PHB:MemoryOrder = "XYZ" ;
!            PHB:stagger = "Z" ;
   call check(nf90_def_var(ncFileID, name="PHB_d0"//idom, xtype=nf90_real, &
        dimids= (/ weDimID(id), snDimID(id), btStagDimID(id) /), varid=phbVarId(id)) )
   call check(nf90_put_att(ncFileID, phbVarId(id), "long_name", &
        "base-state geopotential"))
   call check(nf90_put_att(ncFileID, phbVarId(id), "units", "m^2/s^2"))
   call check(nf90_put_att(ncFileID, phbVarId(id), "units_long_name", "m{2} s{-2}"))

   call check(nf90_def_var(ncFileID, name="HGT_d0"//idom, xtype=nf90_real, &
        dimids= (/ weDimID(id), snDimID(id) /), varid=hgtVarId(id)) )
   call check(nf90_put_att(ncFileID, hgtVarId(id), "long_name", "Terrain Height"))
   call check(nf90_put_att(ncFileID, hgtVarId(id), "units", "m"))
   call check(nf90_put_att(ncFileID, hgtVarId(id), "units_long_name", "meters"))

enddo

if ( output_state_vector ) then

   !-----------------------------------------------------------------
   ! Create attributes for the state vector 
   !-----------------------------------------------------------------

   ! Define the state vector coordinate variable

   call check(nf90_def_var(ncid=ncFileID,name="StateVariable", xtype=nf90_int, &
              dimids=StateVarDimID, varid=StateVarVarID))

   call check(nf90_put_att(ncFileID, StateVarVarID, "long_name", &
        "State Variable ID"))
   call check(nf90_put_att(ncFileID, StateVarVarID, "units", &
        "indexical") )
   call check(nf90_put_att(ncFileID, StateVarVarID, "valid_range", &
        (/ 1, wrf%model_size /)))

   ! Define the actual state vector

   call check(nf90_def_var(ncid=ncFileID, name="state", xtype=nf90_real, &
              dimids = (/ StateVarDimID, MemberDimID, unlimitedDimID /), &
              varid=StateVarID))
   call check(nf90_put_att(ncFileID, StateVarID, "long_name", &
        "model state or fcopy"))
   call check(nf90_put_att(ncFileID, StateVarId, "U_units","m/s"))
   call check(nf90_put_att(ncFileID, StateVarId, "V_units","m/s"))
   call check(nf90_put_att(ncFileID, StateVarId, "W_units","m/s"))
   call check(nf90_put_att(ncFileID, StateVarId, "GZ_units","m2/s2"))
   call check(nf90_put_att(ncFileID, StateVarId, "T_units","K"))
   call check(nf90_put_att(ncFileID, StateVarId, "MU_units","Pa"))
   call check(nf90_put_att(ncFileID, StateVarId, "TSK_units","K"))
   if( wrf%dom(num_domains)%n_moist >= 1) then
      call check(nf90_put_att(ncFileID, StateVarId, "QV_units","kg/kg"))
   endif
   if( wrf%dom(num_domains)%n_moist >= 2) then
      call check(nf90_put_att(ncFileID, StateVarId, "QC_units","kg/kg"))
   endif
   if( wrf%dom(num_domains)%n_moist >= 3) then
      call check(nf90_put_att(ncFileID, StateVarId, "QR_units","kg/kg"))
   endif
   if( wrf%dom(num_domains)%n_moist >= 4) then
      call check(nf90_put_att(ncFileID, StateVarId, "QI_units","kg/kg"))
   endif
   if( wrf%dom(num_domains)%n_moist >= 5) then
      call check(nf90_put_att(ncFileID, StateVarId, "QS_units","kg/kg"))
   endif
   if( wrf%dom(num_domains)%n_moist >= 6) then
      call check(nf90_put_att(ncFileID, StateVarId, "QG_units","kg/kg"))
   endif
   if( wrf%dom(num_domains)%n_moist == 7) then
      call check(nf90_put_att(ncFileID, StateVarId, "QNICE_units","kg-1"))
   endif
   if(wrf%dom(num_domains)%surf_obs ) then
      call check(nf90_put_att(ncFileID, StateVarId, "U10_units","m/s"))
      call check(nf90_put_att(ncFileID, StateVarId, "V10_units","m/s"))
      call check(nf90_put_att(ncFileID, StateVarId, "T2_units","K"))
      call check(nf90_put_att(ncFileID, StateVarId, "TH2_units","K"))
      call check(nf90_put_att(ncFileID, StateVarId, "Q2_units","kg/kg"))
      call check(nf90_put_att(ncFileID, StateVarId, "PS_units","Pa"))
   endif
   if(wrf%dom(num_domains)%soil_data ) then
      call check(nf90_put_att(ncFileID, StateVarId, "TSLB_units","K"))
      call check(nf90_put_att(ncFileID, StateVarId, "SMOIS_units","m3/m3"))
      call check(nf90_put_att(ncFileID, StateVarId, "SH2O_units","m3/m3"))
   endif
   if(h_diab ) then
      call check(nf90_put_att(ncFileID, StateVarId, "H_DIAB_units",""))
   endif

   ! Leave define mode so we can actually fill the variables.

   call check(nf90_enddef(ncfileID))

   call check(nf90_put_var(ncFileID, StateVarVarID, &
        (/ (i,i=1,wrf%model_size) /) ))

else

do id=1,num_domains
   write( idom , '(I1)') id

   !----------------------------------------------------------------------------
   ! Create the (empty) Prognostic Variables and their attributes
   !----------------------------------------------------------------------------

   !      float U(Time, bottom_top, south_north, west_east_stag) ;
   !         U:FieldType = 104 ;
   !         U:MemoryOrder = "XYZ" ;
   !         U:stagger = "X" ;
   call check(nf90_def_var(ncid=ncFileID, name="U_d0"//idom, xtype=nf90_real, &
         dimids = (/ weStagDimID(id), snDimId(id), btDimID(id), MemberDimID, &
         unlimitedDimID /), varid  = var_id))
   call check(nf90_put_att(ncFileID, var_id, "long_name", "x-wind component"))
   call check(nf90_put_att(ncFileID, var_id, "units", "m/s"))
   call check(nf90_put_att(ncFileID, var_id, "units_long_name", "m s{-1}"))


   !      float V(Time, bottom_top, south_north_stag, west_east) ;
   !         V:FieldType = 104 ;
   !         V:MemoryOrder = "XYZ" ;
   !         V:stagger = "Y" ;
   call check(nf90_def_var(ncid=ncFileID, name="V_d0"//idom, xtype=nf90_real, &
         dimids = (/ weDimID(id), snStagDimID(id), btDimID(id), MemberDimID, &
         unlimitedDimID /), varid  = var_id))
   call check(nf90_put_att(ncFileID, var_id, "long_name", "y-wind component"))
   call check(nf90_put_att(ncFileID, var_id, "units", "m/s"))
   call check(nf90_put_att(ncFileID, var_id, "units_long_name", "m s{-1}"))


   !      float W(Time, bottom_top_stag, south_north, west_east) ;
   !         W:FieldType = 104 ;
   !         W:MemoryOrder = "XYZ" ;
   !         W:stagger = "Z" ;
   call check(nf90_def_var(ncid=ncFileID, name="W_d0"//idom, xtype=nf90_real, &
         dimids = (/ weDimID(id), snDimID(id), btStagDimID(id), MemberDimID, &
         unlimitedDimID /), varid  = var_id))
   call check(nf90_put_att(ncFileID, var_id, "long_name", "z-wind component"))
   call check(nf90_put_att(ncFileID, var_id, "units", "m/s"))
   call check(nf90_put_att(ncFileID, var_id, "units_long_name", "m s{-1}"))


   !      float PH(Time, bottom_top_stag, south_north, west_east) ;
   !         PH:FieldType = 104 ;
   !         PH:MemoryOrder = "XYZ" ;
   !         PH:stagger = "Z" ;
   call check(nf90_def_var(ncid=ncFileID, name="PH_d0"//idom, xtype=nf90_real, &
         dimids = (/ weDimID(id), snDimID(id), btStagDimID(id), MemberDimID, &
         unlimitedDimID /), varid  = var_id))
   call check(nf90_put_att(ncFileID, var_id, "long_name", &
        "perturbation geopotential"))
   call check(nf90_put_att(ncFileID, var_id, "units", "m^2/s^2"))
   call check(nf90_put_att(ncFileID, var_id, "units_long_name", "m{2} s{-2}"))


   !      float T(Time, bottom_top, south_north, west_east) ;
   !         T:FieldType = 104 ;
   !         T:MemoryOrder = "XYZ" ;
   !         T:units = "K" ;
   !         T:stagger = "" ;
   call check(nf90_def_var(ncid=ncFileID, name="T_d0"//idom, xtype=nf90_real, &
         dimids = (/ weDimID(id), snDimID(id), btDimID(id), MemberDimID, &
         unlimitedDimID /), varid  = var_id))
   call check(nf90_put_att(ncFileID, var_id, "long_name", "temperature"))
   call check(nf90_put_att(ncFileID, var_id, "units", "K"))
   call check(nf90_put_att(ncFileId, var_id, "description", &
        "perturbation potential temperature (theta-t0)"))


   !      float MU(Time, south_north, west_east) ;
   !         MU:FieldType = 104 ;
   !         MU:MemoryOrder = "XY " ;
   !         MU:description = "perturbation dry air mass in column" ;
   !         MU:units = "pascals" ;
   !         MU:stagger = "" ;
   call check(nf90_def_var(ncid=ncFileID, name="MU_d0"//idom, xtype=nf90_real, &
         dimids = (/ weDimID(id), snDimID(id), MemberDimID, &
         unlimitedDimID /), varid  = var_id))
   call check(nf90_put_att(ncFileID, var_id, "long_name", "mu field"))
   call check(nf90_put_att(ncFileID, var_id, "units", "pascals"))
   call check(nf90_put_att(ncFileId, var_id, "description", &
        "perturbation dry air mass in column"))


   !      float TSK(Time, south_north, west_east) ;
   !         TSK:FieldType = 104 ;
   !         TSK:MemoryOrder = "XY " ;
   !         TSK:description = "SURFACE SKIN TEMPERATURE" ;
   !         TSK:units = "K" ;
   !         TSK:stagger = "" ;
   call check(nf90_def_var(ncid=ncFileID, name="TSK_d0"//idom, xtype=nf90_real, &
         dimids = (/ weDimID(id), snDimID(id), MemberDimID, &
         unlimitedDimID /), varid  = var_id))
   call check(nf90_put_att(ncFileID, var_id, "long_name", "tsk field"))
   call check(nf90_put_att(ncFileID, var_id, "units", "K"))
   call check(nf90_put_att(ncFileId, var_id, "description", &
        "SURFACE SKIN TEMPERATURE"))


   !      float QVAPOR(Time, bottom_top, south_north, west_east) ;
   !         QVAPOR:FieldType = 104 ;
   !         QVAPOR:MemoryOrder = "XYZ" ;
   !         QVAPOR:stagger = "" ;
   if( wrf%dom(id)%n_moist >= 1) then
      call check(nf90_def_var(ncid=ncFileID, name="QVAPOR_d0"//idom, xtype=nf90_real, &
           dimids = (/ weDimID(id), snDimID(id), btDimID(id), MemberDimID, &
           unlimitedDimID /), varid  = var_id))
      call check(nf90_put_att(ncFileID, var_id, "units", "kg/kg"))
      call check(nf90_put_att(ncFileId, var_id, "description", &
           "Water vapor mixing ratio"))
   endif


   !      float QCLOUD(Time, bottom_top, south_north, west_east) ;
   !         QCLOUD:FieldType = 104 ;
   !         QCLOUD:MemoryOrder = "XYZ" ;
   !         QCLOUD:stagger = "" ;
   if( wrf%dom(id)%n_moist >= 2) then
      call check(nf90_def_var(ncid=ncFileID, name="QCLOUD_d0"//idom, xtype=nf90_real, &
           dimids = (/ weDimID(id), snDimID(id), btDimID(id), MemberDimID, &
           unlimitedDimID /), varid  = var_id))
      call check(nf90_put_att(ncFileID, var_id, "units", "kg/kg"))
      call check(nf90_put_att(ncFileID, var_id, "description", &
           "Cloud water mixing ratio"))
   endif


   !      float QRAIN(Time, bottom_top, south_north, west_east) ;
   !         QRAIN:FieldType = 104 ;
   !         QRAIN:MemoryOrder = "XYZ" ;
   !         QRAIN:stagger = "" ;
   if( wrf%dom(id)%n_moist >= 3) then
      call check(nf90_def_var(ncid=ncFileID, name="QRAIN_d0"//idom, xtype=nf90_real, &
           dimids = (/ weDimID(id), snDimID(id), btDimID(id), MemberDimID, &
           unlimitedDimID /), varid  = var_id))
      call check(nf90_put_att(ncFileID, var_id, "units", "kg/kg"))
      call check(nf90_put_att(ncFileID, var_id, "description", &
           "Rain water mixing ratio"))
   endif

   if( wrf%dom(id)%n_moist >= 4) then
      call check(nf90_def_var(ncid=ncFileID, name="QICE_d0"//idom, xtype=nf90_real, &
           dimids = (/ weDimID(id), snDimID(id), btDimID(id), MemberDimID, &
           unlimitedDimID /), varid  = var_id))
      call check(nf90_put_att(ncFileID, var_id, "units", "kg/kg"))
      call check(nf90_put_att(ncFileID, var_id, "description", &
           "Ice mixing ratio"))
   endif

   if( wrf%dom(id)%n_moist >= 5) then
      call check(nf90_def_var(ncid=ncFileID, name="QSNOW_d0"//idom, xtype=nf90_real, &
           dimids = (/ weDimID(id), snDimID(id), btDimID(id), MemberDimID, &
           unlimitedDimID /), varid  = var_id))
      call check(nf90_put_att(ncFileID, var_id, "units", "kg/kg"))
      call check(nf90_put_att(ncFileID, var_id, &
           "description", "Snow mixing ratio"))
   endif

   if( wrf%dom(id)%n_moist >= 6) then
      call check(nf90_def_var(ncid=ncFileID, name="QGRAUP_d0"//idom, xtype=nf90_real, &
           dimids = (/ weDimID(id), snDimID(id), btDimID(id), MemberDimID, &
           unlimitedDimID /), varid  = var_id))
      call check(nf90_put_att(ncFileID, var_id, "units", "kg/kg"))
      call check(nf90_put_att(ncFileID, var_id, "description", &
           "Graupel mixing ratio"))
   endif

   if( wrf%dom(id)%n_moist == 7) then
      call check(nf90_def_var(ncid=ncFileID, name="QNICE_d0"//idom, xtype=nf90_real, &
           dimids = (/ weDimID(id), snDimID(id), btDimID(id), MemberDimID, &
           unlimitedDimID /), varid  = var_id))
      call check(nf90_put_att(ncFileID, var_id, "units", "kg-1"))
      call check(nf90_put_att(ncFileID, var_id, "description", &
           "Ice Number concentration"))
   endif

   if(wrf%dom(id)%surf_obs ) then

      call check(nf90_def_var(ncid=ncFileID, name="U10_d0"//idom, xtype=nf90_real, &
           dimids = (/ weDimID(id), snDimID(id), MemberDimID, unlimitedDimID /), &
           varid  = var_id))
      call check(nf90_put_att(ncFileID, var_id, "units", "m/s"))
      call check(nf90_put_att(ncFileID, var_id, "description", "U at 10 m"))

      call check(nf90_def_var(ncid=ncFileID, name="V10_d0"//idom, xtype=nf90_real, &
           dimids = (/ weDimID(id), snDimID(id), MemberDimID, unlimitedDimID /), &
           varid  = var_id))
      call check(nf90_put_att(ncFileID, var_id, "units", "m/s"))
      call check(nf90_put_att(ncFileID, var_id, "description", "V at 10 m"))

      call check(nf90_def_var(ncid=ncFileID, name="T2_d0"//idom, xtype=nf90_real, &
           dimids = (/ weDimID(id), snDimID(id), MemberDimID, unlimitedDimID /), &
           varid  = var_id))
      call check(nf90_put_att(ncFileID, var_id, "units", "K"))
      call check(nf90_put_att(ncFileID, var_id, "description", "TEMP at 2 m"))

      call check(nf90_def_var(ncid=ncFileID, name="TH2_d0"//idom, xtype=nf90_real, &
           dimids = (/ weDimID(id), snDimID(id), MemberDimID, unlimitedDimID /), &
           varid  = var_id))
      call check(nf90_put_att(ncFileID, var_id, "units", "K"))
      call check(nf90_put_att(ncFileID, var_id, "description", "POT TEMP at 2 m"))

      call check(nf90_def_var(ncid=ncFileID, name="Q2_d0"//idom, xtype=nf90_real, &
           dimids = (/ weDimID(id), snDimID(id), MemberDimID, unlimitedDimID /), &
           varid  = var_id))
      call check(nf90_put_att(ncFileID, var_id, "units", "kg/kg"))
      call check(nf90_put_att(ncFileID, var_id, "description", "QV at 2 m"))

      call check(nf90_def_var(ncid=ncFileID, name="PS_d0"//idom, xtype=nf90_real, &
           dimids = (/ weDimID(id), snDimID(id), MemberDimID, unlimitedDimID /), &
           varid  = var_id))
      call check(nf90_put_att(ncFileID, var_id, "units", "Pa"))
      call check(nf90_put_att(ncFileID, var_id, "description", "Total surface pressure"))

   end if

   if(wrf%dom(id)%soil_data) then   
      !      float TSLB(Time, soil_layers_stag, south_north, west_east) ;
      !         TSLB:FieldType = 104 ;
      !         TSLB:MemoryOrder = "XYZ" ;
      !         TSLB:description = "SOIL TEMPERATURE" ;
      !         TSLB:units = "K" ;
      !         TSLB:stagger = "Z" ;
      call check(nf90_def_var(ncid=ncFileID, name="TSLB_d0"//idom, xtype=nf90_real, &
           dimids = (/ weDimID(id), snDimID(id), slSDimID(id), MemberDimID, &
           unlimitedDimID /), varid  = var_id))
      call check(nf90_put_att(ncFileID, var_id, "long_name", "soil temperature"))
      call check(nf90_put_att(ncFileID, var_id, "units", "K"))
      call check(nf90_put_att(ncFileId, var_id, "description", &
           "SOIL TEMPERATURE"))
           
      call check(nf90_def_var(ncid=ncFileID, name="SMOIS_d0"//idom, xtype=nf90_real, &
            dimids = (/ weDimID(id), snDimID(id), slSDimID(id), MemberDimID, &
            unlimitedDimID /), varid  = var_id))
      call check(nf90_put_att(ncFileID, var_id, "long_name", "soil moisture"))
      call check(nf90_put_att(ncFileID, var_id, "units", "m3/m3"))
      call check(nf90_put_att(ncFileId, var_id, "description", &
           "SOIL MOISTURE"))

      call check(nf90_def_var(ncid=ncFileID, name="SH2O_d0"//idom, xtype=nf90_real, &
            dimids = (/ weDimID(id), snDimID(id), slSDimID(id), MemberDimID, &
            unlimitedDimID /), varid  = var_id))
      call check(nf90_put_att(ncFileID, var_id, "long_name", "soil liquid water"))
      call check(nf90_put_att(ncFileID, var_id, "units", "m3/m3"))
      call check(nf90_put_att(ncFileId, var_id, "description", &
           "SOIL LIQUID WATER"))

   endif

   if(h_diab ) then
      !    float H_DIABATIC(Time, bottom_top, south_north, west_east) ;
      !            H_DIABATIC:FieldType = 104 ;
      !            H_DIABATIC:MemoryOrder = "XYZ" ;
      !            H_DIABATIC:description = "PREVIOUS TIMESTEP CONDENSATIONAL HEATING" ;
      !            H_DIABATIC:units = "" ;
      !            H_DIABATIC:stagger = "" ;
      call check(nf90_def_var(ncid=ncFileID, name="H_DIAB_d0"//idom, xtype=nf90_real, &
            dimids = (/ weDimID(id), snDimID(id), btDimID(id), MemberDimID, &
            unlimitedDimID /), varid  = var_id))
      call check(nf90_put_att(ncFileID, var_id, "long_name", "diabatic heating"))
      call check(nf90_put_att(ncFileID, var_id, "units", ""))
      call check(nf90_put_att(ncFileID, var_id, "FieldType", 104))
      call check(nf90_put_att(ncFileID, var_id, "MemoryOrder", "XYZ"))
      call check(nf90_put_att(ncFileID, var_id, "stagger", ""))
      call check(nf90_put_att(ncFileId, var_id, "description", &
           "previous timestep condensational heating"))
   endif

enddo

endif

!-----------------------------------------------------------------
! Fill the variables we can
!-----------------------------------------------------------------
call check(nf90_enddef(ncfileID))

call check(nf90_put_var(ncFileID,       DXVarID, wrf%dom(1:num_domains)%dx        ))
call check(nf90_put_var(ncFileID,       DYVarID, wrf%dom(1:num_domains)%dy        ))
call check(nf90_put_var(ncFileID, TRUELAT1VarID, wrf%dom(1:num_domains)%proj%truelat1  ))
call check(nf90_put_var(ncFileID, TRUELAT2VarID, wrf%dom(1:num_domains)%proj%truelat2  ))
call check(nf90_put_var(ncFileID,  CEN_LATVarID, wrf%dom(1:num_domains)%cen_lat   ))
call check(nf90_put_var(ncFileID,  CEN_LONVarID, wrf%dom(1:num_domains)%cen_lon   ))
call check(nf90_put_var(ncFileID, MAP_PROJVarID, wrf%dom(1:num_domains)%map_proj  ))

do id=1,num_domains

! defining grid levels
   call check(nf90_put_var(ncFileID,       DNVarID(id), wrf%dom(id)%dn        ))
   call check(nf90_put_var(ncFileID,      ZNUVarID(id), wrf%dom(id)%znu       ))
   call check(nf90_put_var(ncFileID,      DNWVarID(id), wrf%dom(id)%dnw       ))

! defining horizontal
   call check(nf90_put_var(ncFileID,      mubVarID(id), wrf%dom(id)%mub       ))
   call check(nf90_put_var(ncFileID,      LonVarID(id), wrf%dom(id)%longitude ))
   call check(nf90_put_var(ncFileID,      LatVarID(id), wrf%dom(id)%latitude  ))
   call check(nf90_put_var(ncFileID,     ilevVarID(id), (/ (i,i=1,wrf%dom(id)%bt) /) ))
   call check(nf90_put_var(ncFileID,    XlandVarID(id), wrf%dom(id)%land      ))
   call check(nf90_put_var(ncFileID,  MapFacMVarID(id), wrf%dom(id)%mapfac_m  ))
   call check(nf90_put_var(ncFileID,  MapFacUVarID(id), wrf%dom(id)%mapfac_u  ))
   call check(nf90_put_var(ncFileID,  MapFacVVarID(id), wrf%dom(id)%mapfac_v  ))
   call check(nf90_put_var(ncFileID,      phbVarID(id), wrf%dom(id)%phb       ))
   call check(nf90_put_var(ncFileID,      hgtVarID(id), wrf%dom(id)%hgt       ))

enddo

!-----------------------------------------------------------------
! Flush the buffer and leave netCDF file open
!-----------------------------------------------------------------

call check(nf90_sync(ncFileID))

write (*,*)'nc_write_model_atts: netCDF file ',ncFileID,' is synched ...'

contains

  ! Internal subroutine - checks error status after each netcdf, prints
  !                       text message each time an error code is returned.
  subroutine check(istatus)
    integer, intent ( in) :: istatus

    if(istatus /= nf90_noerr) call error_handler(E_ERR, 'nc_write_model_atts', &
       trim(nf90_strerror(istatus)), source, revision, revdate)

  end subroutine check

end function nc_write_model_atts



function nc_write_model_vars( ncFileID, statevec, copyindex, timeindex ) result (ierr)
!-----------------------------------------------------------------
! Writes the model-specific variables to a netCDF file
! TJH 25 June 2003
!
! TJH 29 July 2003 -- for the moment, all errors are fatal, so the
! return code is always '0 == normal', since the fatal errors stop execution.


integer,                intent(in) :: ncFileID      ! netCDF file identifier
real(r8), dimension(:), intent(in) :: statevec
integer,                intent(in) :: copyindex
integer,                intent(in) :: timeindex
integer                            :: ierr          ! return value of function

!-----------------------------------------------------------------

logical, parameter :: debug = .false.  
integer :: nDimensions, nVariables, nAttributes, unlimitedDimID
integer :: StateVarID, VarID, id
integer :: i,j
real(r8), allocatable, dimension(:,:)   :: temp2d
real(r8), allocatable, dimension(:,:,:) :: temp3d
character(len=10) :: varname
character(len=1) :: idom

ierr = 0     ! assume normal termination

!-----------------------------------------------------------------
! make sure ncFileID refers to an open netCDF file, 
! then get all the Variable ID's we need.
!-----------------------------------------------------------------

call check(nf90_Inquire(ncFileID, nDimensions, nVariables, nAttributes, unlimitedDimID))

if ( output_state_vector ) then

   call check(NF90_inq_varid(ncFileID, "state", StateVarID) )
   call check(NF90_put_var(ncFileID, StateVarID, statevec,  &
                start=(/ 1, copyindex, timeindex /)))

else

j = 0

do id=1,num_domains

   write( idom , '(I1)') id

   !----------------------------------------------------------------------------
   ! Fill the variables, the order is CRITICAL  ...   U,V,W,GZ,T,MU,TSK,QV,QC,QR,...
   !----------------------------------------------------------------------------

   !----------------------------------------------------------------------------
   varname = 'U_d0'//idom
   !----------------------------------------------------------------------------
   call check(NF90_inq_varid(ncFileID, trim(adjustl(varname)), VarID))
   i       = j + 1
   j       = i + wrf%dom(id)%wes * wrf%dom(id)%sn * wrf%dom(id)%bt - 1 
   if (debug) write(*,'(a10,'' = statevec('',i7,'':'',i7,'') with dims '',3(1x,i3))') &
              trim(adjustl(varname)),i,j,wrf%dom(id)%wes,wrf%dom(id)%sn,wrf%dom(id)%bt 
   allocate ( temp3d(wrf%dom(id)%wes, wrf%dom(id)%sn, wrf%dom(id)%bt) )
   temp3d  = reshape(statevec(i:j), (/ wrf%dom(id)%wes, wrf%dom(id)%sn, wrf%dom(id)%bt /) ) 
   call check(nf90_put_var( ncFileID, VarID, temp3d, &
                            start=(/ 1, 1, 1, copyindex, timeindex /) ))
   deallocate(temp3d)


   !----------------------------------------------------------------------------
   varname = 'V_d0'//idom
   !----------------------------------------------------------------------------
   call check(NF90_inq_varid(ncFileID, trim(adjustl(varname)), VarID))
   i       = j + 1
   j       = i + wrf%dom(id)%we * wrf%dom(id)%sns * wrf%dom(id)%bt - 1
   if (debug) write(*,'(a10,'' = statevec('',i7,'':'',i7,'') with dims '',3(1x,i3))') &
              trim(adjustl(varname)),i,j,wrf%dom(id)%we,wrf%dom(id)%sns,wrf%dom(id)%bt
   allocate ( temp3d(wrf%dom(id)%we, wrf%dom(id)%sns, wrf%dom(id)%bt) )
   temp3d  = reshape(statevec(i:j), (/ wrf%dom(id)%we, wrf%dom(id)%sns, wrf%dom(id)%bt /) ) 
   call check(nf90_put_var( ncFileID, VarID, temp3d, &
                            start=(/ 1, 1, 1, copyindex, timeindex /) ))
   deallocate(temp3d)


   !----------------------------------------------------------------------------
   varname = 'W_d0'//idom
   !----------------------------------------------------------------------------
   call check(NF90_inq_varid(ncFileID, trim(adjustl(varname)), VarID))
   i       = j + 1
   j       = i + wrf%dom(id)%we * wrf%dom(id)%sn * wrf%dom(id)%bts - 1
   if (debug) write(*,'(a10,'' = statevec('',i7,'':'',i7,'') with dims '',3(1x,i3))') &
              trim(adjustl(varname)),i,j,wrf%dom(id)%we,wrf%dom(id)%sn,wrf%dom(id)%bts
   allocate ( temp3d(wrf%dom(id)%we, wrf%dom(id)%sn, wrf%dom(id)%bts) )
   temp3d  = reshape(statevec(i:j), (/ wrf%dom(id)%we, wrf%dom(id)%sn, wrf%dom(id)%bts /) ) 
   call check(nf90_put_var( ncFileID, VarID, temp3d, &
                            start=(/ 1, 1, 1, copyindex, timeindex /) ))


   !----------------------------------------------------------------------------
   varname = 'PH_d0'//idom       ! AKA "GZ"
   !----------------------------------------------------------------------------
   call check(NF90_inq_varid(ncFileID, trim(adjustl(varname)), VarID))
   i       = j + 1
   j       = i + wrf%dom(id)%we * wrf%dom(id)%sn * wrf%dom(id)%bts - 1
   if (debug) write(*,'(a10,'' = statevec('',i7,'':'',i7,'') with dims '',3(1x,i3))') &
              trim(adjustl(varname)),i,j,wrf%dom(id)%we,wrf%dom(id)%sn,wrf%dom(id)%bts
   temp3d  = reshape(statevec(i:j), (/ wrf%dom(id)%we, wrf%dom(id)%sn, wrf%dom(id)%bts /) ) 
   call check(nf90_put_var( ncFileID, VarID, temp3d, &
                            start=(/ 1, 1, 1, copyindex, timeindex /) ))
   deallocate(temp3d)


   !----------------------------------------------------------------------------
   varname = 'T_d0'//idom
   !----------------------------------------------------------------------------
   call check(NF90_inq_varid(ncFileID, trim(adjustl(varname)), VarID))
   i       = j + 1
   j       = i + wrf%dom(id)%we * wrf%dom(id)%sn * wrf%dom(id)%bt - 1
   if (debug) write(*,'(a10,'' = statevec('',i7,'':'',i7,'') with dims '',3(1x,i3))') &
              trim(adjustl(varname)),i,j,wrf%dom(id)%we,wrf%dom(id)%sn,wrf%dom(id)%bt
   allocate ( temp3d(wrf%dom(id)%we, wrf%dom(id)%sn, wrf%dom(id)%bt) )
   temp3d  = reshape(statevec(i:j), (/ wrf%dom(id)%we, wrf%dom(id)%sn, wrf%dom(id)%bt /) ) 
   call check(nf90_put_var( ncFileID, VarID, temp3d, &
                            start=(/ 1, 1, 1, copyindex, timeindex /) ))
   deallocate(temp3d)


   !----------------------------------------------------------------------------
   varname = 'MU_d0'//idom
   !----------------------------------------------------------------------------
   call check(NF90_inq_varid(ncFileID, trim(adjustl(varname)), VarID))
   i       = j + 1
   j       = i + wrf%dom(id)%we * wrf%dom(id)%sn - 1
   if (debug) write(*,'(a10,'' = statevec('',i7,'':'',i7,'') with dims '',3(1x,i3))') &
              trim(adjustl(varname)),i,j,wrf%dom(id)%we,wrf%dom(id)%sn
   allocate ( temp2d(wrf%dom(id)%we, wrf%dom(id)%sn) )
   temp2d  = reshape(statevec(i:j), (/ wrf%dom(id)%we, wrf%dom(id)%sn /) ) 
   call check(nf90_put_var( ncFileID, VarID, temp2d, &
                            start=(/ 1, 1, copyindex, timeindex /) ))


   !----------------------------------------------------------------------------
   varname = 'TSK_d0'//idom
   !----------------------------------------------------------------------------
   call check(NF90_inq_varid(ncFileID, trim(adjustl(varname)), VarID))
   i       = j + 1
   j       = i + wrf%dom(id)%we * wrf%dom(id)%sn - 1
   if (debug) write(*,'(a10,'' = statevec('',i7,'':'',i7,'') with dims '',3(1x,i3))') &
              trim(adjustl(varname)),i,j,wrf%dom(id)%we,wrf%dom(id)%sn
   temp2d  = reshape(statevec(i:j), (/ wrf%dom(id)%we, wrf%dom(id)%sn /) ) 
   call check(nf90_put_var( ncFileID, VarID, temp2d, &
                            start=(/ 1, 1, copyindex, timeindex /) ))


   if( wrf%dom(id)%n_moist >= 1) then
      !----------------------------------------------------------------------------
      varname = 'QVAPOR_d0'//idom
      !----------------------------------------------------------------------------
      call check(NF90_inq_varid(ncFileID, trim(adjustl(varname)), VarID))
      i       = j + 1
      j       = i + wrf%dom(id)%we * wrf%dom(id)%sn * wrf%dom(id)%bt - 1
      if (debug) write(*,'(a10,'' = statevec('',i7,'':'',i7,'') with dims '',3(1x,i3))') &
           trim(adjustl(varname)),i,j,wrf%dom(id)%we,wrf%dom(id)%sn,wrf%dom(id)%bt
      allocate ( temp3d(wrf%dom(id)%we, wrf%dom(id)%sn, wrf%dom(id)%bt) )
      temp3d  = reshape(statevec(i:j), (/ wrf%dom(id)%we, wrf%dom(id)%sn, wrf%dom(id)%bt /) ) 
      call check(nf90_put_var( ncFileID, VarID, temp3d, &
           start=(/ 1, 1, 1, copyindex, timeindex /) ))
   endif


   if( wrf%dom(id)%n_moist >= 2) then
      !----------------------------------------------------------------------------
      varname = 'QCLOUD_d0'//idom
      !----------------------------------------------------------------------------
      call check(NF90_inq_varid(ncFileID, trim(adjustl(varname)), VarID))
      i       = j + 1
      j       = i + wrf%dom(id)%we * wrf%dom(id)%sn * wrf%dom(id)%bt - 1
      if (debug) write(*,'(a10,'' = statevec('',i7,'':'',i7,'') with dims '',3(1x,i3))') &
           trim(adjustl(varname)),i,j,wrf%dom(id)%we,wrf%dom(id)%sn,wrf%dom(id)%bt
      temp3d  = reshape(statevec(i:j), (/ wrf%dom(id)%we, wrf%dom(id)%sn, wrf%dom(id)%bt /) ) 
      call check(nf90_put_var( ncFileID, VarID, temp3d, &
           start=(/ 1, 1, 1, copyindex, timeindex /) ))
   endif


   if( wrf%dom(id)%n_moist >= 3) then
      !----------------------------------------------------------------------------
      varname = 'QRAIN_d0'//idom
      !----------------------------------------------------------------------------
      call check(NF90_inq_varid(ncFileID, trim(adjustl(varname)), VarID))
      i       = j + 1
      j       = i + wrf%dom(id)%we * wrf%dom(id)%sn * wrf%dom(id)%bt - 1
      if (debug) write(*,'(a10,'' = statevec('',i7,'':'',i7,'') with dims '',3(1x,i3))') &
           trim(adjustl(varname)),i,j,wrf%dom(id)%we,wrf%dom(id)%sn,wrf%dom(id)%bt
      temp3d  = reshape(statevec(i:j), (/ wrf%dom(id)%we, wrf%dom(id)%sn, wrf%dom(id)%bt /) ) 
      call check(nf90_put_var( ncFileID, VarID, temp3d, &
           start=(/ 1, 1, 1, copyindex, timeindex /) ))
   endif


   if( wrf%dom(id)%n_moist >= 4) then
      !----------------------------------------------------------------------------
      varname = 'QICE_d0'//idom
      !----------------------------------------------------------------------------
      call check(NF90_inq_varid(ncFileID, trim(adjustl(varname)), VarID))
      i       = j + 1
      j       = i + wrf%dom(id)%we * wrf%dom(id)%sn * wrf%dom(id)%bt - 1
      if (debug) write(*,'(a10,'' = statevec('',i7,'':'',i7,'') with dims '',3(1x,i3))') &
           trim(adjustl(varname)),i,j,wrf%dom(id)%we,wrf%dom(id)%sn,wrf%dom(id)%bt
      temp3d  = reshape(statevec(i:j), (/ wrf%dom(id)%we, wrf%dom(id)%sn, wrf%dom(id)%bt /) ) 
      call check(nf90_put_var( ncFileID, VarID, temp3d, &
           start=(/ 1, 1, 1, copyindex, timeindex /) ))
   endif
   if( wrf%dom(id)%n_moist >= 5) then
      !----------------------------------------------------------------------------
      varname = 'QSNOW_d0'//idom
      !----------------------------------------------------------------------------
      call check(NF90_inq_varid(ncFileID, trim(adjustl(varname)), VarID))
      i       = j + 1
      j       = i + wrf%dom(id)%we * wrf%dom(id)%sn * wrf%dom(id)%bt - 1
      if (debug) write(*,'(a10,'' = statevec('',i7,'':'',i7,'') with dims '',3(1x,i3))') &
           trim(adjustl(varname)),i,j,wrf%dom(id)%we,wrf%dom(id)%sn,wrf%dom(id)%bt
      temp3d  = reshape(statevec(i:j), (/ wrf%dom(id)%we, wrf%dom(id)%sn, wrf%dom(id)%bt /) ) 
      call check(nf90_put_var( ncFileID, VarID, temp3d, &
           start=(/ 1, 1, 1, copyindex, timeindex /) ))
   endif
   if( wrf%dom(id)%n_moist >= 6) then
      !----------------------------------------------------------------------------
      varname = 'QGRAUP_d0'//idom
      !----------------------------------------------------------------------------
      call check(NF90_inq_varid(ncFileID, trim(adjustl(varname)), VarID))
      i       = j + 1
      j       = i + wrf%dom(id)%we * wrf%dom(id)%sn * wrf%dom(id)%bt - 1
      if (debug) write(*,'(a10,'' = statevec('',i7,'':'',i7,'') with dims '',3(1x,i3))') &
           trim(adjustl(varname)),i,j,wrf%dom(id)%we,wrf%dom(id)%sn,wrf%dom(id)%bt
      temp3d  = reshape(statevec(i:j), (/ wrf%dom(id)%we, wrf%dom(id)%sn, wrf%dom(id)%bt /) ) 
      call check(nf90_put_var( ncFileID, VarID, temp3d, &
           start=(/ 1, 1, 1, copyindex, timeindex /) ))
   endif
   if( wrf%dom(id)%n_moist == 7) then
      !----------------------------------------------------------------------------
      varname = 'QNICE_d0'//idom
      !----------------------------------------------------------------------------
      call check(NF90_inq_varid(ncFileID, trim(adjustl(varname)), VarID))
      i       = j + 1
      j       = i + wrf%dom(id)%we * wrf%dom(id)%sn * wrf%dom(id)%bt - 1
      if (debug) write(*,'(a10,'' = statevec('',i7,'':'',i7,'') with dims '',3(1x,i3))') &
           trim(adjustl(varname)),i,j,wrf%dom(id)%we,wrf%dom(id)%sn,wrf%dom(id)%bt
      temp3d  = reshape(statevec(i:j), (/ wrf%dom(id)%we, wrf%dom(id)%sn, wrf%dom(id)%bt /) ) 
      call check(nf90_put_var( ncFileID, VarID, temp3d, &
           start=(/ 1, 1, 1, copyindex, timeindex /) ))
   endif

   deallocate(temp3d)

   if ( wrf%dom(id)%n_moist > 7 ) then
      write(*,'(''wrf%dom(id)%n_moist = '',i3)')wrf%dom(id)%n_moist
      call error_handler(E_ERR,'nc_write_model_vars', &
               'num_moist_vars is too large.', source, revision, revdate)
   endif

   if(wrf%dom(id)%surf_obs ) then

      !----------------------------------------------------------------------------
      varname = 'U10_d0'//idom
      !----------------------------------------------------------------------------
      call check(NF90_inq_varid(ncFileID, trim(adjustl(varname)), VarID))
      i       = j + 1
      j       = i + wrf%dom(id)%we * wrf%dom(id)%sn - 1
      if (debug) write(*,'(a10,'' = statevec('',i7,'':'',i7,'') with dims '',3(1x,i3))') &
           trim(adjustl(varname)),i,j,wrf%dom(id)%we,wrf%dom(id)%sn
      temp2d  = reshape(statevec(i:j), (/ wrf%dom(id)%we, wrf%dom(id)%sn /) ) 
      call check(nf90_put_var( ncFileID, VarID, temp2d, &
           start=(/ 1, 1, copyindex, timeindex /) ))

      !----------------------------------------------------------------------------
      varname = 'V10_d0'//idom
      !----------------------------------------------------------------------------
      call check(NF90_inq_varid(ncFileID, trim(adjustl(varname)), VarID))
      i       = j + 1
      j       = i + wrf%dom(id)%we * wrf%dom(id)%sn - 1
      if (debug) write(*,'(a10,'' = statevec('',i7,'':'',i7,'') with dims '',3(1x,i3))') &
           trim(adjustl(varname)),i,j,wrf%dom(id)%we,wrf%dom(id)%sn
      temp2d  = reshape(statevec(i:j), (/ wrf%dom(id)%we, wrf%dom(id)%sn /) ) 
      call check(nf90_put_var( ncFileID, VarID, temp2d, &
           start=(/ 1, 1, copyindex, timeindex /) ))

      !----------------------------------------------------------------------------
      varname = 'T2_d0'//idom
      !----------------------------------------------------------------------------
      call check(NF90_inq_varid(ncFileID, trim(adjustl(varname)), VarID))
      i       = j + 1
      j       = i + wrf%dom(id)%we * wrf%dom(id)%sn - 1
      if (debug) write(*,'(a10,'' = statevec('',i7,'':'',i7,'') with dims '',3(1x,i3))') &
           trim(adjustl(varname)),i,j,wrf%dom(id)%we,wrf%dom(id)%sn
      temp2d  = reshape(statevec(i:j), (/ wrf%dom(id)%we, wrf%dom(id)%sn /) ) 
      call check(nf90_put_var( ncFileID, VarID, temp2d, &
           start=(/ 1, 1, copyindex, timeindex /) ))

      !----------------------------------------------------------------------------
      varname = 'TH2_d0'//idom
      !----------------------------------------------------------------------------
      call check(NF90_inq_varid(ncFileID, trim(adjustl(varname)), VarID))
      i       = j + 1
      j       = i + wrf%dom(id)%we * wrf%dom(id)%sn - 1
      if (debug) write(*,'(a10,'' = statevec('',i7,'':'',i7,'') with dims '',3(1x,i3))') &
           trim(adjustl(varname)),i,j,wrf%dom(id)%we,wrf%dom(id)%sn
      temp2d  = reshape(statevec(i:j), (/ wrf%dom(id)%we, wrf%dom(id)%sn /) ) 
      call check(nf90_put_var( ncFileID, VarID, temp2d, &
           start=(/ 1, 1, copyindex, timeindex /) ))

      !----------------------------------------------------------------------------
      varname = 'Q2_d0'//idom
      !----------------------------------------------------------------------------
      call check(NF90_inq_varid(ncFileID, trim(adjustl(varname)), VarID))
      i       = j + 1
      j       = i + wrf%dom(id)%we * wrf%dom(id)%sn - 1
      if (debug) write(*,'(a10,'' = statevec('',i7,'':'',i7,'') with dims '',3(1x,i3))') &
           trim(adjustl(varname)),i,j,wrf%dom(id)%we,wrf%dom(id)%sn
      temp2d  = reshape(statevec(i:j), (/ wrf%dom(id)%we, wrf%dom(id)%sn /) ) 
      call check(nf90_put_var( ncFileID, VarID, temp2d, &
           start=(/ 1, 1, copyindex, timeindex /) ))

      !----------------------------------------------------------------------------
      varname = 'PS_d0'//idom
      !----------------------------------------------------------------------------
      call check(NF90_inq_varid(ncFileID, trim(adjustl(varname)), VarID))
      i       = j + 1
      j       = i + wrf%dom(id)%we * wrf%dom(id)%sn - 1
      if (debug) write(*,'(a10,'' = statevec('',i7,'':'',i7,'') with dims '',3(1x,i3))') &
           trim(adjustl(varname)),i,j,wrf%dom(id)%we,wrf%dom(id)%sn
      temp2d  = reshape(statevec(i:j), (/ wrf%dom(id)%we, wrf%dom(id)%sn /) ) 
      call check(nf90_put_var( ncFileID, VarID, temp2d, &
           start=(/ 1, 1, copyindex, timeindex /) ))

   endif

   if(wrf%dom(id)%soil_data ) then   

      !----------------------------------------------------------------------------
      varname = 'TSLB_d0'//idom
      !----------------------------------------------------------------------------
      call check(NF90_inq_varid(ncFileID, trim(adjustl(varname)), VarID))
      i       = j + 1
      j       = i + wrf%dom(id)%we * wrf%dom(id)%sn * wrf%dom(id)%sls - 1
      if (debug) write(*,'(a10,'' = statevec('',i7,'':'',i7,'') with dims '',3(1x,i3))') &
                 trim(adjustl(varname)),i,j,wrf%dom(id)%we,wrf%dom(id)%sn,wrf%dom(id)%sls
      allocate ( temp3d(wrf%dom(id)%we, wrf%dom(id)%sn, wrf%dom(id)%sls) )
      temp3d  = reshape(statevec(i:j), (/ wrf%dom(id)%we, wrf%dom(id)%sn, wrf%dom(id)%sls /) ) 
      call check(nf90_put_var( ncFileID, VarID, temp3d, &
                               start=(/ 1, 1, 1, copyindex, timeindex /) ))

      !----------------------------------------------------------------------------
      varname = 'SMOIS_d0'//idom
      !----------------------------------------------------------------------------
      call check(NF90_inq_varid(ncFileID, trim(adjustl(varname)), VarID))
      i       = j + 1
      j       = i + wrf%dom(id)%we * wrf%dom(id)%sn * wrf%dom(id)%sls - 1
      if (debug) write(*,'(a10,'' = statevec('',i7,'':'',i7,'') with dims '',3(1x,i3))') &
                 trim(adjustl(varname)),i,j,wrf%dom(id)%we,wrf%dom(id)%sn,wrf%dom(id)%sls
      temp3d  = reshape(statevec(i:j), (/ wrf%dom(id)%we, wrf%dom(id)%sn, wrf%dom(id)%sls /) ) 
      call check(nf90_put_var( ncFileID, VarID, temp3d, &
                               start=(/ 1, 1, 1, copyindex, timeindex /) ))

      !----------------------------------------------------------------------------
      varname = 'SH2O_d0'//idom
      !----------------------------------------------------------------------------
      call check(NF90_inq_varid(ncFileID, trim(adjustl(varname)), VarID))
      i       = j + 1
      j       = i + wrf%dom(id)%we * wrf%dom(id)%sn * wrf%dom(id)%sls - 1
      if (debug) write(*,'(a10,'' = statevec('',i7,'':'',i7,'') with dims '',3(1x,i3))') &
                 trim(adjustl(varname)),i,j,wrf%dom(id)%we,wrf%dom(id)%sn,wrf%dom(id)%sls
      temp3d  = reshape(statevec(i:j), (/ wrf%dom(id)%we, wrf%dom(id)%sn, wrf%dom(id)%sls /) ) 
      call check(nf90_put_var( ncFileID, VarID, temp3d, &
                               start=(/ 1, 1, 1, copyindex, timeindex /) ))

      deallocate(temp3d)

   endif

   deallocate(temp2d)

   if( h_diab) then

      !----------------------------------------------------------------------------
      varname = 'H_DIAB_d0'//idom
      !----------------------------------------------------------------------------
      call check(NF90_inq_varid(ncFileID, trim(adjustl(varname)), VarID))
      i       = j + 1
      j       = i + wrf%dom(id)%we * wrf%dom(id)%sn * wrf%dom(id)%bt - 1
      if (debug) write(*,'(a10,'' = statevec('',i7,'':'',i7,'') with dims '',3(1x,i3))') &
                 trim(adjustl(varname)),i,j,wrf%dom(id)%we,wrf%dom(id)%sn,wrf%dom(id)%bt
      allocate ( temp3d(wrf%dom(id)%we, wrf%dom(id)%sn, wrf%dom(id)%bt) )
      temp3d  = reshape(statevec(i:j), (/ wrf%dom(id)%we, wrf%dom(id)%sn, wrf%dom(id)%bt /) ) 
      call check(nf90_put_var( ncFileID, VarID, temp3d, &
                               start=(/ 1, 1, 1, copyindex, timeindex /) ))
      deallocate(temp3d)

   end if

enddo

endif

!-----------------------------------------------------------------
! Flush the buffer and leave netCDF file open
!-----------------------------------------------------------------

write (*,*)'Finished filling variables ...'
call check(nf90_sync(ncFileID))
write (*,*)'netCDF file is synched ...'

contains

  ! Internal subroutine - checks error status after each netcdf, prints
  !                       text message each time an error code is returned.
  subroutine check(istatus)
    integer, intent ( in) :: istatus

    if(istatus /= nf90_noerr) call error_handler(E_ERR, 'nc_write_model_vars', &
         trim(nf90_strerror(istatus)), source, revision, revdate)

  end subroutine check

end function nc_write_model_vars

!-------------------------------

!  public stubs

!**********************************************

subroutine adv_1step(x, Time)

! Does single time-step advance with vector state as
! input and output.

  real(r8), intent(inout) :: x(:)

! Time is needed for more general models like this; need to add in to 
! low-order models
  type(time_type), intent(in) :: Time

end subroutine adv_1step

!**********************************************

subroutine end_model()
end subroutine end_model

!**********************************************

subroutine init_time(i_time)
! For now returns value of Time_init which is set in initialization routines.

  type(time_type), intent(out) :: i_time

!Where should initial time come from here?
! WARNING: CURRENTLY SET TO 0
  i_time = set_time(0, 0)

end subroutine init_time

!**********************************************

subroutine init_conditions(x)
! Reads in restart initial conditions and converts to vector

! Following changed to intent(inout) for ifc compiler;should be like this
  real(r8), intent(inout) :: x(:)

end subroutine init_conditions



!#######################################################################

subroutine toGrid (x, j, dx, dxm)

!  Transfer obs. x to grid j and calculate its
!  distance to grid j and j+1

  real(r8), intent(in)  :: x
  real(r8), intent(out) :: dx, dxm
  integer,  intent(out) :: j

  j = int (x)

  dx = x - real (j)

  dxm= 1.0_r8 - dx

end subroutine toGrid

!#######################################################################

subroutine pres_to_zk(pres, mdl_v, n3, zk)

! Calculate the model level "zk" on half (mass) levels,
! corresponding to pressure "pres".

  integer,  intent(in)  :: n3
  real(r8), intent(in)  :: pres
  real(r8), intent(in)  :: mdl_v(0:n3)
  real(r8), intent(out) :: zk

  integer  :: k

  zk = missing_r8

  if (pres > mdl_v(0) .or. pres < mdl_v(n3)) return

  do k = 0,n3-1
     if(pres <= mdl_v(k) .and. pres >= mdl_v(k+1)) then
        zk = real(k) + (mdl_v(k) - pres)/(mdl_v(k) - mdl_v(k+1))
        exit
     endif
  enddo

end subroutine pres_to_zk

!#######################################################################

subroutine height_to_zk(obs_v, mdl_v, n3, zk)

! Calculate the model level "zk" on half (mass) levels,
! corresponding to height "obs_v".

  real(r8), intent(in)  :: obs_v
  integer,  intent(in)  :: n3
  real(r8), intent(in)  :: mdl_v(0:n3)
  real(r8), intent(out) :: zk

  integer   :: k

  zk = missing_r8

  if (obs_v < mdl_v(0) .or. obs_v > mdl_v(n3)) return

  do k = 0,n3-1
     if(obs_v >= mdl_v(k) .and. obs_v <= mdl_v(k+1)) then
        zk = real(k) + (mdl_v(k) - obs_v)/(mdl_v(k) - mdl_v(k+1))
        exit
     endif
  enddo

end subroutine height_to_zk

!#######################################################

subroutine get_model_pressure_profile(i,j,dx,dy,dxm,dym,n,x,id,v_p)

! Calculate the full model pressure profile on half (mass) levels,
! horizontally interpolated at the observation location.

integer,  intent(in)  :: i,j,n,id
real(r8), intent(in)  :: dx,dy,dxm,dym
real(r8), intent(in)  :: x(:)
real(r8), intent(out) :: v_p(0:n)

integer  :: i1,i2,k
real(r8) :: pres1, pres2, pres3, pres4

if(i >= 1 .and. i < wrf%dom(id)%var_size(1,TYPE_T) .and. &
   j >= 1 .and. j < wrf%dom(id)%var_size(2,TYPE_T)) then

   do k=1,n
      pres1 = model_pressure_t(i  ,j  ,k,id,x)
      pres2 = model_pressure_t(i+1,j  ,k,id,x)
      pres3 = model_pressure_t(i  ,j+1,k,id,x)
      pres4 = model_pressure_t(i+1,j+1,k,id,x)
      v_p(k) = dym*( dxm*pres1 + dx*pres2 ) + dy*( dxm*pres3 + dx*pres4 )
   enddo

   if(wrf%dom(id)%surf_obs ) then

!!$      i1 = get_wrf_index(i,j,1,TYPE_PS,id)
!!$      i2 = get_wrf_index(i,j+1,1,TYPE_PS,id)
      i1 = wrf%dom(id)%dart_ind(i,j,1,TYPE_PS)
      i2 = wrf%dom(id)%dart_ind(i,j+1,1,TYPE_PS)
      if(x(i1) /= 0.0_r8 .and. x(i1+1) /= 0.0_r8 .and. &
           x(i2) /= 0.0_r8 .and. x(i2+1) /= 0.0_r8) then

         v_p(0) = dym*( dxm*x(i1) + dx*x(i1+1) ) + dy*( dxm*x(i2) + dx*x(i2+1) )

      else

         pres1 = model_pressure_t(i  ,j  ,2,id,x)
         pres2 = model_pressure_t(i+1,j  ,2,id,x)
         pres3 = model_pressure_t(i  ,j+1,2,id,x)
         pres4 = model_pressure_t(i+1,j+1,2,id,x)
         v_p(0) = (3.0_r8*v_p(1) - &
              dym*( dxm*pres1 + dx*pres2 ) - dy*( dxm*pres3 + dx*pres4 ))/2.0_r8

      endif

   else

      pres1 = model_pressure_t(i  ,j  ,2,id,x)
      pres2 = model_pressure_t(i+1,j  ,2,id,x)
      pres3 = model_pressure_t(i  ,j+1,2,id,x)
      pres4 = model_pressure_t(i+1,j+1,2,id,x)
      v_p(0) = (3.0_r8*v_p(1) - &
           dym*( dxm*pres1 + dx*pres2 ) - dy*( dxm*pres3 + dx*pres4 ))/2.0_r8

   endif

else

   v_p(:) = missing_r8

endif

end subroutine get_model_pressure_profile

!#######################################################

function model_pressure(i,j,k,id,var_type,x)

! Calculate the pressure at grid point (i,j,k), domain id.
! The grid is defined according to var_type.

integer,  intent(in)  :: i,j,k,id,var_type
real(r8), intent(in)  :: x(:)
real(r8)              :: model_pressure

integer  :: ips, imu
real(r8) :: pres1, pres2

model_pressure = missing_r8

if( (var_type == type_w) .or. (var_type == type_gz) ) then

   if( k == 1 ) then

      pres1 = model_pressure_t(i,j,k,id,x)
      pres2 = model_pressure_t(i,j,k+1,id,x)
      model_pressure = (3.0_r8*pres1 - pres2)/2.0_r8

   elseif( k == wrf%dom(id)%var_size(3,TYPE_W) ) then

      pres1 = model_pressure_t(i,j,k-1,id,x)
      pres2 = model_pressure_t(i,j,k-2,id,x)
      model_pressure = (3.0_r8*pres1 - pres2)/2.0_r8

   else

      pres1 = model_pressure_t(i,j,k,id,x)
      pres2 = model_pressure_t(i,j,k-1,id,x)
      model_pressure = (pres1 + pres2)/2.0_r8

   endif

elseif( var_type == type_u ) then

   if( i == wrf%dom(id)%var_size(1,TYPE_U) ) then

      pres1 = model_pressure_t(i-1,j,k,id,x)
      pres2 = model_pressure_t(i-2,j,k,id,x)
      model_pressure = (3.0_r8*pres1 - pres2)/2.0_r8

   elseif( i == 1 ) then

      pres1 = model_pressure_t(i,j,k,id,x)
      pres2 = model_pressure_t(i+1,j,k,id,x)
      model_pressure = (3.0_r8*pres1 - pres2)/2.0_r8

   else

      pres1 = model_pressure_t(i,j,k,id,x)
      pres2 = model_pressure_t(i-1,j,k,id,x)
      model_pressure = (pres1 + pres2)/2.0_r8

   endif

elseif( var_type == type_v ) then

   if( j == wrf%dom(id)%var_size(2,TYPE_V) ) then

      pres1 = model_pressure_t(i,j-1,k,id,x)
      pres2 = model_pressure_t(i,j-2,k,id,x)
      model_pressure = (3.0_r8*pres1 - pres2)/2.0_r8

   elseif( j == 1 ) then

      pres1 = model_pressure_t(i,j,k,id,x)
      pres2 = model_pressure_t(i,j+1,k,id,x)
      model_pressure = (3.0_r8*pres1 - pres2)/2.0_r8

   else

      pres1 = model_pressure_t(i,j,k,id,x)
      pres2 = model_pressure_t(i,j-1,k,id,x)
      model_pressure = (pres1 + pres2)/2.0_r8

   endif

elseif( var_type == type_mu    .or. var_type == type_tslb .or. &
        var_type == type_ps    .or. var_type == type_u10  .or. &
        var_type == type_v10   .or. var_type == type_t2   .or. &
        var_type == type_th2   .or.                            &
        var_type == type_q2    .or. var_type == type_tsk  .or. &
        var_type == type_smois .or. var_type == type_sh2o) then

   model_pressure = model_pressure_s(i,j,id,x)

else

   model_pressure = model_pressure_t(i,j,k,id,x)

endif

end function model_pressure

!#######################################################

function model_pressure_t(i,j,k,id,x)

! Calculate total pressure on mass point (half (mass) levels, T-point).

integer,  intent(in)  :: i,j,k,id
real(r8), intent(in)  :: x(:)
real(r8)              :: model_pressure_t

real (kind=r8), PARAMETER    :: rd_over_rv = gas_constant / gas_constant_v
real (kind=r8), PARAMETER    :: cpovcv = 1.4_r8        ! cp / (cp - gas_constant)

integer  :: iqv,it
real(r8) :: qvf1,rho

model_pressure_t = missing_r8

! Adapted the code from WRF module_big_step_utilities_em.F ----
!         subroutine calc_p_rho_phi      Y.-R. Guo (10/20/2004)

! Simplification: alb*mub = (phb(i,j,k+1) - phb(i,j,k))/dnw(k)

!!$iqv = get_wrf_index(i,j,k,TYPE_QV,id)
!!$it  = get_wrf_index(i,j,k,TYPE_T,id)
iqv = wrf%dom(id)%dart_ind(i,j,k,TYPE_QV)
it  = wrf%dom(id)%dart_ind(i,j,k,TYPE_T)

qvf1 = 1.0_r8 + x(iqv) / rd_over_rv

rho = model_rho_t(i,j,k,id,x)

! .. total pressure:
model_pressure_t = ps0 * ( (gas_constant*(ts0+x(it))*qvf1) / &
     (ps0/rho) )**cpovcv

end function model_pressure_t

!#######################################################

function model_pressure_s(i,j,id,x)

! compute pressure at surface at mass point

integer,  intent(in)  :: i,j,id
real(r8), intent(in)  :: x(:)
real(r8)              :: model_pressure_s

integer  :: ips, imu

if(wrf%dom(id)%surf_obs ) then
   ips = wrf%dom(id)%dart_ind(i,j,1,TYPE_PS)
   model_pressure_s = x(ips)

else
   imu = wrf%dom(id)%dart_ind(i,j,1,TYPE_MU)
   model_pressure_s = wrf%dom(id)%p_top + wrf%dom(id)%mub(i,j) + x(imu)

endif


end function model_pressure_s

!#######################################################

function model_rho_t(i,j,k,id,x)

! Calculate the total density on mass point (half (mass) levels, T-point).

integer,  intent(in)  :: i,j,k,id
real(r8), intent(in)  :: x(:)
real(r8)              :: model_rho_t

integer  :: imu,iph,iphp1
real(r8) :: ph_e

model_rho_t = missing_r8

! Adapted the code from WRF module_big_step_utilities_em.F ----
!         subroutine calc_p_rho_phi      Y.-R. Guo (10/20/2004)

! Simplification: alb*mub = (phb(i,j,k+1) - phb(i,j,k))/dnw(k)

!!$imu = get_wrf_index(i,j,1,TYPE_MU,id)
!!$iph = get_wrf_index(i,j,k,TYPE_GZ,id)
!!$iphp1 = get_wrf_index(i,j,k+1,TYPE_GZ,id)
imu = wrf%dom(id)%dart_ind(i,j,1,TYPE_MU)
iph = wrf%dom(id)%dart_ind(i,j,k,TYPE_GZ)
iphp1 = wrf%dom(id)%dart_ind(i,j,k+1,TYPE_GZ)

ph_e = ( (x(iphp1) + wrf%dom(id)%phb(i,j,k+1)) - (x(iph) + wrf%dom(id)%phb(i,j,k)) ) &
     /wrf%dom(id)%dnw(k)

! now calculate rho = - mu / dphi/deta

model_rho_t = - (wrf%dom(id)%mub(i,j)+x(imu)) / ph_e

end function model_rho_t

!#######################################################

subroutine get_model_height_profile(i,j,dx,dy,dxm,dym,n,x,id,v_h)

! Calculate the model height profile on half (mass) levels,
! horizontally interpolated at the observation location.

integer,  intent(in)  :: i,j,n,id
real(r8), intent(in)  :: dx,dy,dxm,dym
real(r8), intent(in)  :: x(:)
real(r8), intent(out) :: v_h(0:n)

real(r8)  :: fll(n+1)
integer   :: i1,i2,k

if(i >= 1 .and. i < wrf%dom(id)%var_size(1,TYPE_GZ) .and. &
   j >= 1 .and. j < wrf%dom(id)%var_size(2,TYPE_GZ)) then

   do k = 1, wrf%dom(id)%var_size(3,TYPE_GZ)
!!$      i1 = get_wrf_index(i,j,k,TYPE_GZ,id)
!!$      i2 = get_wrf_index(i,j+1,k,TYPE_GZ,id)
      i1 = wrf%dom(id)%dart_ind(i,j,k,TYPE_GZ)
      i2 = wrf%dom(id)%dart_ind(i,j+1,k,TYPE_GZ)
      fll(k) = (dym*( dxm*(wrf%dom(id)%phb(i,j,k)+x(i1)) + &
           dx*(wrf%dom(id)%phb(i+1,j,k)+x(i1+1))) + &
           dy*(dxm*(wrf%dom(id)%phb(i,j+1,k)+x(i2)) + &
           dx*(wrf%dom(id)%phb(i+1,j+1,k)+x(i2+1)) ))/gravity
   end do

   do k=1,n
      v_h(k) = 0.5_r8*(fll(k) + fll(k+1) )
   end do

   v_h(0) = dym*( dxm*wrf%dom(id)%hgt(i  ,j  ) + &
                   dx*wrf%dom(id)%hgt(i+1,j  ) ) + &
             dy*( dxm*wrf%dom(id)%hgt(i  ,j+1) + &
                   dx*wrf%dom(id)%hgt(i+1,j+1) )

else

   print*,'Not able the get height_profile'
   print*,i,j,dx,dy,dxm,dym,n,id,wrf%dom(id)%var_size(1,TYPE_GZ), &
        wrf%dom(id)%var_size(2,TYPE_GZ)

   v_h(:) =  missing_r8

endif

end subroutine get_model_height_profile

!#######################################################

function model_height(i,j,k,id,var_type,x)

integer,  intent(in)  :: i,j,k,id,var_type
real(r8), intent(in)  :: x(:)
real(r8)              :: model_height

integer   :: i1, i2, i3, i4

model_height = missing_r8

if( (var_type == type_w) .or. (var_type == type_gz) ) then

!!$   i1 = get_wrf_index(i,j,k,TYPE_GZ,id)
   i1 = wrf%dom(id)%dart_ind(i,j,k,TYPE_GZ)
   model_height = (wrf%dom(id)%phb(i,j,k)+x(i1))/gravity

elseif( var_type == type_u ) then

   if( i == wrf%dom(id)%var_size(1,TYPE_U) ) then

!!$      i1 = get_wrf_index(i-1,j,k  ,TYPE_GZ,id)
!!$      i2 = get_wrf_index(i-1,j,k+1,TYPE_GZ,id)
      i1 = wrf%dom(id)%dart_ind(i-1,j,k  ,TYPE_GZ)
      i2 = wrf%dom(id)%dart_ind(i-1,j,k+1,TYPE_GZ)

      model_height = ( 3.0_r8*(wrf%dom(id)%phb(i-1,j,k  )+x(i1)) &
                      +3.0_r8*(wrf%dom(id)%phb(i-1,j,k+1)+x(i2)) &
                             -(wrf%dom(id)%phb(i-2,j,k  )+x(i1-1)) &
                             -(wrf%dom(id)%phb(i-2,j,k+1)+x(i2-1)) )/(4.0_r8*gravity)

   elseif( i == 1 ) then

!!$      i1 = get_wrf_index(i,j,k  ,TYPE_GZ,id)
!!$      i2 = get_wrf_index(i,j,k+1,TYPE_GZ,id)
      i1 = wrf%dom(id)%dart_ind(i,j,k  ,TYPE_GZ)
      i2 = wrf%dom(id)%dart_ind(i,j,k+1,TYPE_GZ)

      model_height = ( 3.0_r8*(wrf%dom(id)%phb(i  ,j,k  )+x(i1)) &
                      +3.0_r8*(wrf%dom(id)%phb(i  ,j,k+1)+x(i2)) &
                             -(wrf%dom(id)%phb(i+1,j,k  )+x(i1+1)) &
                             -(wrf%dom(id)%phb(i+1,j,k+1)+x(i2+1)) )/(4.0_r8*gravity)

   else

!!$      i1 = get_wrf_index(i,j,k  ,TYPE_GZ,id)
!!$      i2 = get_wrf_index(i,j,k+1,TYPE_GZ,id)
      i1 = wrf%dom(id)%dart_ind(i,j,k  ,TYPE_GZ)
      i2 = wrf%dom(id)%dart_ind(i,j,k+1,TYPE_GZ)

      model_height = ( (wrf%dom(id)%phb(i  ,j,k  )+x(i1)) &
                      +(wrf%dom(id)%phb(i  ,j,k+1)+x(i2)) &
                      +(wrf%dom(id)%phb(i-1,j,k  )+x(i1-1)) &
                      +(wrf%dom(id)%phb(i-1,j,k+1)+x(i2-1)) )/(4.0_r8*gravity)

   endif

elseif( var_type == type_v ) then

   if( j == wrf%dom(id)%var_size(2,TYPE_V) ) then

!!$      i1 = get_wrf_index(i,j-1,k  ,TYPE_GZ,id)
!!$      i2 = get_wrf_index(i,j-1,k+1,TYPE_GZ,id)
!!$      i3 = get_wrf_index(i,j-2,k  ,TYPE_GZ,id)
!!$      i4 = get_wrf_index(i,j-2,k+1,TYPE_GZ,id)
      i1 = wrf%dom(id)%dart_ind(i,j-1,k  ,TYPE_GZ)
      i2 = wrf%dom(id)%dart_ind(i,j-1,k+1,TYPE_GZ)
      i3 = wrf%dom(id)%dart_ind(i,j-2,k  ,TYPE_GZ)
      i4 = wrf%dom(id)%dart_ind(i,j-2,k+1,TYPE_GZ)

      model_height = ( 3.0_r8*(wrf%dom(id)%phb(i,j-1,k  )+x(i1)) &
                      +3.0_r8*(wrf%dom(id)%phb(i,j-1,k+1)+x(i2)) &
                             -(wrf%dom(id)%phb(i,j-2,k  )+x(i3)) &
                             -(wrf%dom(id)%phb(i,j-2,k+1)+x(i4)) )/(4.0_r8*gravity)

   elseif( j == 1 ) then

!!$      i1 = get_wrf_index(i,j  ,k  ,TYPE_GZ,id)
!!$      i2 = get_wrf_index(i,j  ,k+1,TYPE_GZ,id)
!!$      i3 = get_wrf_index(i,j+1,k  ,TYPE_GZ,id)
!!$      i4 = get_wrf_index(i,j+1,k+1,TYPE_GZ,id)
      i1 = wrf%dom(id)%dart_ind(i,j  ,k  ,TYPE_GZ)
      i2 = wrf%dom(id)%dart_ind(i,j  ,k+1,TYPE_GZ)
      i3 = wrf%dom(id)%dart_ind(i,j+1,k  ,TYPE_GZ)
      i4 = wrf%dom(id)%dart_ind(i,j+1,k+1,TYPE_GZ)

      model_height = ( 3.0_r8*(wrf%dom(id)%phb(i,j  ,k  )+x(i1)) &
                      +3.0_r8*(wrf%dom(id)%phb(i,j  ,k+1)+x(i2)) &
                             -(wrf%dom(id)%phb(i,j+1,k  )+x(i3)) &
                             -(wrf%dom(id)%phb(i,j+1,k+1)+x(i4)) )/(4.0_r8*gravity)

   else

!!$      i1 = get_wrf_index(i,j  ,k  ,TYPE_GZ,id)
!!$      i2 = get_wrf_index(i,j  ,k+1,TYPE_GZ,id)
!!$      i3 = get_wrf_index(i,j-1,k  ,TYPE_GZ,id)
!!$      i4 = get_wrf_index(i,j-1,k+1,TYPE_GZ,id)
      i1 = wrf%dom(id)%dart_ind(i,j  ,k  ,TYPE_GZ)
      i2 = wrf%dom(id)%dart_ind(i,j  ,k+1,TYPE_GZ)
      i3 = wrf%dom(id)%dart_ind(i,j-1,k  ,TYPE_GZ)
      i4 = wrf%dom(id)%dart_ind(i,j-1,k+1,TYPE_GZ)

      model_height = ( (wrf%dom(id)%phb(i,j  ,k  )+x(i1)) &
                      +(wrf%dom(id)%phb(i,j  ,k+1)+x(i2)) &
                      +(wrf%dom(id)%phb(i,j-1,k  )+x(i3)) &
                      +(wrf%dom(id)%phb(i,j-1,k+1)+x(i4)) )/(4.0_r8*gravity)

   endif

elseif( var_type == type_mu .or. var_type == type_ps .or. &
        var_type == type_tsk) then

   model_height = wrf%dom(id)%hgt(i,j)

elseif( var_type == type_tslb .or. var_type == type_smois .or. &
        var_type == type_sh2o ) then

   model_height = wrf%dom(id)%hgt(i,j) - wrf%dom(id)%zs(k)

elseif( var_type == type_u10 .or. var_type == type_v10 ) then

   model_height = wrf%dom(id)%hgt(i,j) + 10.0_r8

elseif( var_type == type_t2 .or. var_type == type_th2 .or. var_type == type_q2 ) then

   model_height = wrf%dom(id)%hgt(i,j) + 2.0_r8

else

!!$   i1 = get_wrf_index(i,j,k  ,TYPE_GZ,id)
!!$   i2 = get_wrf_index(i,j,k+1,TYPE_GZ,id)
   i1 = wrf%dom(id)%dart_ind(i,j,k  ,TYPE_GZ)
   i2 = wrf%dom(id)%dart_ind(i,j,k+1,TYPE_GZ)

   model_height = ( (wrf%dom(id)%phb(i,j,k  )+x(i1)) &
                   +(wrf%dom(id)%phb(i,j,k+1)+x(i2)) )/(2.0_r8*gravity)

endif

end function model_height

!#######################################################

function model_height_w(i,j,k,id,x)

! return total height at staggered vertical coordinate
! and horizontal mass coordinates

integer,  intent(in)  :: i,j,k,id
real(r8), intent(in)  :: x(:)
real(r8)              :: model_height_w

integer   :: i1

i1 = wrf%dom(id)%dart_ind(i,j,k,TYPE_GZ)

model_height_w = (wrf%dom(id)%phb(i,j,k) + x(i1))/gravity

end function model_height_w

!#######################################################


subroutine pert_model_state(state, pert_state, interf_provided)

! Perturbs a model state for generating initial ensembles
! Returning interf_provided means go ahead and do this with uniform
! small independent perturbations.

real(r8), intent(in)  :: state(:)
real(r8), intent(out) :: pert_state(:)
logical,  intent(out) :: interf_provided

interf_provided = .false.
pert_state = state

end subroutine pert_model_state

!#######################################################

subroutine read_dt_from_wrf_nml()

real(r8) :: dt

integer :: time_step, time_step_fract_num, time_step_fract_den
integer :: max_dom, feedback, smooth_option
integer, dimension(3) :: s_we, e_we, s_sn, e_sn, s_vert, e_vert
integer, dimension(3) :: dx, dy, ztop, grid_id, parent_id
integer, dimension(3) :: i_parent_start, j_parent_start, parent_grid_ratio
integer, dimension(3) :: parent_time_step_ratio
integer :: io, iunit, id

namelist /domains/ time_step, time_step_fract_num, time_step_fract_den
namelist /domains/ max_dom
namelist /domains/ s_we, e_we, s_sn, e_sn, s_vert, e_vert
namelist /domains/ dx, dy, ztop, grid_id, parent_id
namelist /domains/ i_parent_start, j_parent_start, parent_grid_ratio
namelist /domains/ parent_time_step_ratio
namelist /domains/ feedback, smooth_option


! Begin by reading the namelist input
call find_namelist_in_file("namelist.input", "domains", iunit)
read(iunit, nml = domains, iostat = io)
call check_namelist_read(iunit, io, "domains")

! Record the namelist values used for the run ...
call error_handler(E_MSG,'read_dt_from_wrf_nml','domains namelist values are',' ',' ',' ')
write(logfileunit, nml=domains)
write(     *     , nml=domains)

if (max_dom /= num_domains) then

   write(*,*) 'max_dom in namelist.input = ',max_dom
   write(*,*) 'num_domains in input.nml  = ',num_domains
   call error_handler(E_ERR,'read_dt_from_wrf_nml', &
        'Make them consistent.', source, revision,revdate)

endif

if (time_step_fract_den /= 0) then
   dt = real(time_step) + real(time_step_fract_num) / real(time_step_fract_den)
else
   dt = real(time_step)
endif

do id=1,num_domains
   wrf%dom(id)%dt = dt / real(parent_time_step_ratio(id))
enddo

end subroutine read_dt_from_wrf_nml



subroutine compute_seaprs ( nz, z, t, p , q ,          &
                            sea_level_pressure, debug)
!-------------------------------------------------------------------------
! compute_seaprs    Estimate sea level pressure.
!
! This routines has been taken "as is" from wrf_user_fortran_util_0.f
!
! This routine assumes
!    index order is (i,j,k)
!    wrf staggering
!    units: pressure (Pa), temperature(K), height (m), mixing ratio (kg kg{-1})
!    availability of 3d p, t, and qv; 2d terrain; 1d half-level zeta string
!    output units of SLP are Pa, but you should divide that by 100 for the
!          weather weenies.
!    virtual effects are included
!
! Dave
!
! cys: change to 1d
! TJH: verified intent() qualifiers, declaration syntax, uses error_handler

      IMPLICIT NONE
      INTEGER,  intent(in)    :: nz
      REAL(r8), intent(in)    :: z(nz), p(nz), q(nz)
      REAL(r8), intent(inout) :: t(nz)
      REAL(r8), intent(out)   :: sea_level_pressure
      LOGICAL,  intent(in)    :: debug

      INTEGER  :: level
      REAL(r8) :: t_surf, t_sea_level

!     Some required physical constants:

      REAL(r8) :: R, G, GAMMA
      PARAMETER (R=287.04_r8, G=9.81_r8, GAMMA=0.0065_r8)

!     Specific constants for assumptions made in this routine:

      REAL(r8) :: TC, PCONST
      PARAMETER (TC=273.16_r8 + 17.5_r8, PCONST = 10000.0_r8)

      LOGICAL  :: ridiculous_mm5_test
      PARAMETER  (ridiculous_mm5_test = .TRUE.)
!     PARAMETER  (ridiculous_mm5_test = .false.)

!     Local variables:

      character(len=129) :: errstring
      INTEGER :: k
      INTEGER :: klo, khi

      REAL(r8) :: plo, phi, tlo, thi, zlo, zhi
      REAL(r8) :: p_at_pconst, t_at_pconst, z_at_pconst
      REAL(r8) :: z_half_lowest

      REAL(r8), PARAMETER :: cp           = 7.0_r8*R/2.0_r8
      REAL(r8), PARAMETER :: rcp          = R/cp
      REAL(r8), PARAMETER :: p1000mb      = 100000.0_r8

      LOGICAL ::  l1 , l2 , l3, found

!     Find least zeta level that is PCONST Pa above the surface.  We later use this
!     level to extrapolate a surface pressure and temperature, which is supposed
!     to reduce the effect of the diurnal heating cycle in the pressure field.

      t = t*(p/p1000mb)**rcp

      level = -1

      k = 1
      found = .false.
      do while( (.not. found) .and. (k.le.nz))
         IF ( p(k) .LT. p(1)-PCONST ) THEN
            level = k
            found = .true.
         END IF
         k = k+1
      END DO

      IF ( level .EQ. -1 ) THEN
         PRINT '(A,I4,A)','Troubles finding level ',   &
               NINT(PCONST)/100,' above ground.'
         print*, 'p=',p
         print*, 't=',t
         print*, 'z=',z
         print*, 'q=',q
         write(errstring,*)'Error_in_finding_100_hPa_up'
         call error_handler(E_ERR,'compute_seaprs',errstring,' ',' ',' ')
      END IF


!     Get temperature PCONST Pa above surface.  Use this to extrapolate
!     the temperature at the surface and down to sea level.

      klo = MAX ( level - 1 , 1      )
      khi = MIN ( klo + 1        , nz - 1 )

      IF ( klo .EQ. khi ) THEN
         PRINT '(A)','Trapping levels are weird.'
         PRINT '(A,I3,A,I3,A)','klo = ',klo,', khi = ',khi, &
               ': and they should not be equal.'
         write(errstring,*)'Error_trapping_levels'
         call error_handler(E_ERR,'compute_seaprs',errstring,' ',' ',' ')
      END IF

      plo = p(klo)
      phi = p(khi)
      tlo = t(klo)*(1. + 0.608 * q(klo) )
      thi = t(khi)*(1. + 0.608 * q(khi) )
!     zlo = zetahalf(klo)/ztop*(ztop-terrain(i,j))+terrain(i,j)
!     zhi = zetahalf(khi)/ztop*(ztop-terrain(i,j))+terrain(i,j)
      zlo = z(klo)
      zhi = z(khi)

      p_at_pconst = p(1) - pconst
      t_at_pconst = thi-(thi-tlo)*LOG(p_at_pconst/phi)*LOG(plo/phi)
      z_at_pconst = zhi-(zhi-zlo)*LOG(p_at_pconst/phi)*LOG(plo/phi)

      t_surf = t_at_pconst*(p(1)/p_at_pconst)**(gamma*R/g)
      t_sea_level = t_at_pconst+gamma*z_at_pconst


!     If we follow a traditional computation, there is a correction to the sea level
!     temperature if both the surface and sea level temnperatures are *too* hot.

      IF ( ridiculous_mm5_test ) THEN
         l1 = t_sea_level .LT. TC
         l2 = t_surf      .LE. TC
         l3 = .NOT. l1
         IF ( l2 .AND. l3 ) THEN
            t_sea_level = TC
         ELSE
            t_sea_level = TC - 0.005*(t_surf-TC)**2
         END IF
      END IF

!     The grand finale: ta da!

!     z_half_lowest=zetahalf(1)/ztop*(ztop-terrain(i,j))+terrain(i,j)
      z_half_lowest=z(1)
      sea_level_pressure = p(1) *              &
                           EXP((2.*g*z_half_lowest)/   &
                           (R*(t_sea_level+t_surf)))

!        sea_level_pressure(i,j) = sea_level_pressure(i,j)*0.01

    if (debug) then
      print *,'slp=',sea_level_pressure
    endif
!      print *,'t=',t(10:15,10:15,1),t(10:15,2,1),t(10:15,3,1)
!      print *,'z=',z(10:15,1,1),z(10:15,2,1),z(10:15,3,1)
!      print *,'p=',p(10:15,1,1),p(10:15,2,1),p(10:15,3,1)
!      print *,'slp=',sea_level_pressure(10:15,10:15),     &
!         sea_level_pressure(10:15,10:15),sea_level_pressure(20,10:15)

end subroutine compute_seaprs


      
SUBROUTINE splin2(x1a,x2a,ya,y2a,m,n,x1,x2,y)
      INTEGER m,n,NN
      REAL(r8), intent(in) :: x1,x2,x1a(m),x2a(n),y2a(m,n),ya(m,n)
      real(r8), intent(out) :: y
      PARAMETER (NN=100)
!     USES spline,splint
      INTEGER j,k
      REAL(r8) y2tmp(NN),ytmp(NN),yytmp(NN)
      do 12 j=1,m
        do 11 k=1,n
          ytmp(k)=ya(j,k)
          y2tmp(k)=y2a(j,k)
11      continue
        call splint(x2a,ytmp,y2tmp,n,x2,yytmp(j))
12    continue
      call spline(x1a,yytmp,m,1.e30_r8,1.e30_r8,y2tmp)
      call splint(x1a,yytmp,y2tmp,m,x1,y)
      return
END subroutine splin2

SUBROUTINE splie2(x1a,x2a,ya,m,n,y2a)
      INTEGER m,n,NN
      REAL(r8), intent(in) :: x1a(m),x2a(n),ya(m,n)
      REAL(r8), intent(out) :: y2a(m,n)
      PARAMETER (NN=100)
!     USES spline
      INTEGER j,k
      REAL(r8) y2tmp(NN),ytmp(NN)
      do 13 j=1,m
        do 11 k=1,n
          ytmp(k)=ya(j,k)
11      continue
        call spline(x2a,ytmp,n,1.e30_r8,1.e30_r8,y2tmp)
        do 12 k=1,n
          y2a(j,k)=y2tmp(k)
12      continue
13    continue
      return
END subroutine splie2

SUBROUTINE spline(x,y,n,yp1,ypn,y2)
      INTEGER n,NMAX
      REAL(r8), intent(in) :: yp1,ypn,x(n),y(n)
      REAL(r8), intent(out) :: y2(n)
      PARAMETER (NMAX=500)
      INTEGER i,k
      REAL(r8) p,qn,sig,un,u(NMAX)
      if (yp1.gt..99e30) then
        y2(1)=0.
        u(1)=0.
      else
        y2(1)=-0.5
        u(1)=(3./(x(2)-x(1)))*((y(2)-y(1))/(x(2)-x(1))-yp1)
      endif
      do 11 i=2,n-1
        sig=(x(i)-x(i-1))/(x(i+1)-x(i-1))
        p=sig*y2(i-1)+2.
        y2(i)=(sig-1.)/p
        u(i)=(6.*((y(i+1)-y(i))/(x(i+ &
      1)-x(i))-(y(i)-y(i-1))/(x(i)-x(i-1)))/(x(i+1)-x(i-1))-sig* &
      u(i-1))/p
11    continue
      if (ypn.gt..99e30) then
        qn=0.
        un=0.
      else
        qn=0.5
        un=(3./(x(n)-x(n-1)))*(ypn-(y(n)-y(n-1))/(x(n)-x(n-1)))
      endif
      y2(n)=(un-qn*u(n-1))/(qn*y2(n-1)+1.)
      do 12 k=n-1,1,-1
        y2(k)=y2(k)*y2(k+1)+u(k)
12    continue
      return
END subroutine spline


SUBROUTINE splint(xa,ya,y2a,n,x,y)
      INTEGER n
      REAL(r8),intent(in) :: x,xa(n),y2a(n),ya(n)
      REAL(r8),intent(out) :: y
      INTEGER k,khi,klo
      REAL(r8) a,b,h
      klo=1
      khi=n
1     if (khi-klo.gt.1) then
        k=(khi+klo)/2
        if(xa(k).gt.x)then
          khi=k
        else
          klo=k
        endif
      goto 1
      endif
      h=xa(khi)-xa(klo)
      if (h.eq.0.) pause 'bad xa input in splint'
      a=(xa(khi)-x)/h
      b=(x-xa(klo))/h
      y=a*ya(klo)+b*ya(khi)+((a**3-a)*y2a(klo)+(b**3-b)*y2a(khi))*(h**2)/6.
      return
END subroutine splint


!#######################################################################


subroutine ens_mean_for_model(filter_ens_mean)

! Not used in low-order models
! Stores provided ensemble mean within the module for later use

real(r8), intent(in) :: filter_ens_mean(:)

ens_mean = filter_ens_mean

end subroutine ens_mean_for_model


!#######################################################################

subroutine get_domain_info(obslon,obslat,id,iloc,jloc)

real(r8), intent(in)  :: obslon, obslat
integer, intent(out)  :: id
real(r8), intent(out) :: iloc, jloc

logical               :: dom_found

! given arbitrary lat and lon values, returns closest domain id and
! horizontal mass point grid points (xloc,yloc)

dom_found = .false.

id = num_domains
do while (.not. dom_found)

   ! Checking for exact equality on real variable types is generally a bad idea.

   if( (wrf%dom(id)%proj%hemi ==  1.0_r8 .and. obslat == -90.0_r8) .or. &
       (wrf%dom(id)%proj%hemi == -1.0_r8 .and. obslat ==  90.0_r8) .or. &
       (wrf%dom(id)%proj%code == PROJ_MERC .and. abs(obslat) >= 90.0_r8) ) then

   else
      call latlon_to_ij(wrf%dom(id)%proj,obslat,obslon,iloc,jloc)
      if ( (iloc >= 1 .and. iloc <= wrf%dom(id)%we .and. &
           jloc >= 1 .and. jloc <= wrf%dom(id)%sn) ) then
         dom_found = .true.
      endif

   endif

   if (.not. dom_found) then
      id = id - 1
      if (id == 0) return
   endif

end do

end subroutine get_domain_info

!#######################################################################

subroutine get_close_obs(gc, base_obs_loc, base_obs_kind, obs_loc, obs_kind, &
                            num_close, close_ind, dist)

! Given a DART ob (referred to as "base") and a set of obs priors or state variables
! (obs_loc, obs_kind), returns the subset of close ones to the "base" ob, their
! indices, and their distances to the "base" ob...

! For vertical distance computations, general philosophy is to convert all vertical
! coordinates to a common coordinate. This coordinate type is defined in the namelist
! with the variable "vert_localization_coord".

! Vertical conversion is carried out by the subroutine vert_interpolate.

! Note that both base_obs_loc and obs_loc are intent(inout), meaning that these
! locations are possibly modified here and returned as such to the calling routine.
! The calling routine is always filter_assim and these arrays are local arrays
! within filter_assim. In other words, these modifications will only matter within
! filter_assim, but will not propagate backwards to filter.
      
implicit none

type(get_close_type), intent(in)     :: gc
type(location_type),  intent(inout)  :: base_obs_loc, obs_loc(:)
integer,              intent(in)     :: base_obs_kind, obs_kind(:)
integer,              intent(out)    :: num_close, close_ind(:)
real(r8),             intent(out)    :: dist(:)

integer                :: t_ind, istatus1, istatus2, k
integer                :: base_which, local_obs_which
real(r8), dimension(3) :: base_array, local_obs_array
type(location_type)    :: local_obs_loc


! Initialize variables to missing status
num_close = 0
close_ind = -99
dist      = 1.0e9

istatus1 = 0
istatus2 = 0

! Convert base_obs vertical coordinate to requested vertical coordinate if necessary

base_array = get_location(base_obs_loc) 
base_which = nint(query_location(base_obs_loc))

if (.not. horiz_dist_only) then
   if (base_which /= wrf%dom(1)%vert_coord) then
      call vert_interpolate(ens_mean, base_obs_loc, base_obs_kind, istatus1)
   elseif (base_array(3) == missing_r8) then
      istatus1 = 1
   end if
endif

if (istatus1 == 0) then

   ! Get all the potentially close obs but no dist (optional argument dist(:) is not present)
   ! This way, we are decreasing the number of distance computations that will follow.
   ! This is a horizontal-distance operation and we don't need to have the relevant vertical
   ! coordinate information yet (for obs_loc).
   call loc_get_close_obs(gc, base_obs_loc, base_obs_kind, obs_loc, obs_kind, &
                          num_close, close_ind)

   ! Loop over potentially close subset of obs priors or state variables
   do k = 1, num_close

      t_ind = close_ind(k)
      local_obs_loc   = obs_loc(t_ind)
      local_obs_which = nint(query_location(local_obs_loc))

      ! Convert local_obs vertical coordinate to requested vertical coordinate if necessary.
      ! This should only be necessary for obs priors, as state location information already
      ! contains the correct vertical coordinate (filter_assim's call to get_state_meta_data).
      if (.not. horiz_dist_only) then
         if (local_obs_which /= wrf%dom(1)%vert_coord) then
            call vert_interpolate(ens_mean, local_obs_loc, obs_kind(t_ind), istatus2)
            ! Store the "new" location into the original full local array
            obs_loc(t_ind) = local_obs_loc
         endif
      endif

      ! Compute distance - set distance to a very large value if vert coordinate is missing
      ! or vert_interpolate returned error (istatus2=1)
      local_obs_array = get_location(local_obs_loc)
      if (((.not. horiz_dist_only).and.(local_obs_array(3) == missing_r8)).or.(istatus2 == 1)) then
         dist(k) = 1.0e9        
      else
         dist(k) = get_dist(base_obs_loc, local_obs_loc, base_obs_kind, obs_kind(t_ind))
      end if

   end do
endif

end subroutine get_close_obs



end module model_mod
