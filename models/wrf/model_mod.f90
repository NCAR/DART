! DART software - Copyright UCAR. This open source software is provided
! by UCAR, "as is", without charge, subject to all terms of use at
! http://www.image.ucar.edu/DAReS/DART/DART_download
!

module model_mod

use        types_mod, only : r8, i8, MISSING_R8, digits12, &
                             gas_constant, gas_constant_v, ps0, gravity

use time_manager_mod, only : time_type, set_time, GREGORIAN, set_date, &
                             set_calendar_type

use     location_mod, only : location_type, get_close_type, &
                             set_location, set_location_missing, &
                             set_vertical_localization_coord, set_vertical, &
                             VERTISHEIGHT, VERTISLEVEL, VERTISPRESSURE, &
                             VERTISSURFACE, VERTISUNDEF, VERTISSCALEHEIGHT, &
                             loc_get_close => get_close, get_location, &
                             query_location, is_vertical,  vertical_localization_on, &
                             get_dist

use    utilities_mod, only : register_module, error_handler, &
                             E_ERR, E_MSG, &
                             nmlfileunit, do_output, do_nml_file, do_nml_term,  &
                             find_namelist_in_file, check_namelist_read, &
                             to_upper

use netcdf_utilities_mod, only : nc_add_global_attribute, nc_synchronize_file, &
                                 nc_add_global_creation_time, &
                                 nc_begin_define_mode, nc_end_define_mode, &
                                 NF90_MAX_NAME, nc_get_variable_size, &
                                 nc_get_variable, nc_close_file, nc_check, &
                                 nc_open_file_readonly, nc_get_variable_size, &
                                 nc_get_global_attribute, nc_get_dimension_size

use state_structure_mod, only : add_domain, get_domain_size, get_model_variable_indices, &
                                get_dim_name, get_num_dims, get_dart_vector_index, &
                                get_varid_from_kind, get_varid_from_varname

use distributed_state_mod, only : get_state_array, get_state

use obs_kind_mod, only : get_index_for_quantity, &
                         QTY_U_WIND_COMPONENT, &
                         QTY_v_WIND_COMPONENT, &
                         QTY_10M_U_WIND_COMPONENT, &
                         QTY_10M_V_WIND_COMPONENT, &
                         QTY_DENSITY, &
                         QTY_GEOPOTENTIAL_HEIGHT, &
                         QTY_PRESSURE, &
                         QTY_SURFACE_TYPE, &
                         QTY_SURFACE_ELEVATION, &
                         QTY_LANDMASK, &
                         QTY_SURFACE_PRESSURE, &
                         QTY_VAPOR_MIXING_RATIO, &
                         QTY_TEMPERATURE, &
                         QTY_POTENTIAL_TEMPERATURE, &
                         QTY_DENSITY, &
                         QTY_VERTICAL_VELOCITY, &
                         QTY_SPECIFIC_HUMIDITY, &
                         QTY_VAPOR_MIXING_RATIO, &
                         QTY_SURFACE_PRESSURE, &
                         QTY_VORTEX_LAT, &
                         QTY_VORTEX_LON, &
                         QTY_VORTEX_PMIN,QTY_VORTEX_WMAX, &
                         QTY_SKIN_TEMPERATURE, &
                         QTY_SURFACE_TYPE, &
                         QTY_2M_TEMPERATURE, &
                         QTY_2M_SPECIFIC_HUMIDITY, &
                         QTY_RAINWATER_MIXING_RATIO,    &
                         QTY_GRAUPEL_MIXING_RATIO,      &
                         QTY_HAIL_MIXING_RATIO,         &
                         QTY_SNOW_MIXING_RATIO,         &
                         QTY_ICE_MIXING_RATIO,          &
                         QTY_CLOUDWATER_MIXING_RATIO,   &
                         QTY_DROPLET_NUMBER_CONCENTR,   &
                         QTY_ICE_NUMBER_CONCENTRATION,  &
                         QTY_SNOW_NUMBER_CONCENTR,      &
                         QTY_RAIN_NUMBER_CONCENTR,      &
                         QTY_GRAUPEL_NUMBER_CONCENTR,   &
                         QTY_HAIL_NUMBER_CONCENTR, &
                         QTY_SOIL_TEMPERATURE, &
                         QTY_SOIL_MOISTURE, &
                         QTY_SOIL_LIQUID_WATER

use ensemble_manager_mod, only : ensemble_type

use default_model_mod, only : write_model_time, &
                              init_time => fail_init_time, &
                              init_conditions => fail_init_conditions, &
                              adv_1step

use         map_utils, only : latlon_to_ij, &
                              proj_info, &
                              map_set, &
                              map_init, &
                              gridwind_to_truewind, &
                              PROJ_LATLON, &
                              PROJ_LC, &
                              PROJ_PS, &
                              PROJ_PS_WGS84, &
                              PROJ_MERC, &
                              PROJ_GAUSS, &
                              PROJ_CYL, &
                              PROJ_CASSINI, &
                              PROJ_ALBERS_NAD83, &
                              PROJ_ROTLL

use netcdf ! no get_char in netcdf_utilities_mod

implicit none
private

! routines required by DART code - will be called from filter and other
! DART executables. 
public :: get_model_size,         &
          get_state_meta_data,    &
          model_interpolate,      &
          end_model,              &
          static_init_model,      &
          nc_write_model_atts,    &
          get_close_obs,          &
          get_close_state,        &
          pert_model_copies,      &
          convert_vertical_obs,   &
          convert_vertical_state, &
          read_model_time,        &
          adv_1step,              &
          init_time,              &
          init_conditions,        &
          shortest_time_between_assimilations, &
          write_model_time

! for wrf_dart_obs_preprocess.f90
public :: get_domain_info


! module variables
character(len=256), parameter :: source   = "wrf/model_mod.f90"
logical :: module_initialized = .false.

integer, parameter :: MAX_STATE_VARIABLES = 100
integer, parameter :: NUM_STATE_TABLE_COLUMNS = 4
integer, parameter :: NUM_BOUNDS_TABLE_COLUMNS = 4
 
integer, allocatable :: wrf_dom(:) ! This needs a better name, it is the id from add_domain
                                   ! for each wrf_domain added to the state

!-- Namelist with default values --
logical :: default_state_variables = .true.
character(len=NF90_MAX_NAME) :: wrf_state_variables(MAX_STATE_VARIABLES*NUM_STATE_TABLE_COLUMNS) = 'NULL'
character(len=NF90_MAX_NAME) :: wrf_state_bounds(num_bounds_table_columns,max_state_variables) = 'NULL'
integer :: num_domains = 1
integer :: calendar_type        = GREGORIAN
integer :: assimilation_period_seconds = 21600
! Max height a surface obs can be away from the actual model surface
! and still be accepted (in meters)
real (kind=r8) :: sfc_elev_max_diff  = -1.0_r8   ! could be something like 200.0_r8

real (kind=r8) :: center_search_half_length = 500000.0_r8
real(r8) :: circulation_pres_level = 80000.0_r8
real(r8) :: circulation_radius     = 108000.0_r8
integer :: center_spline_grid_scale = 10
integer :: vert_localization_coord = VERTISHEIGHT

! Allow observations above the surface but below the lowest sigma level.
logical :: allow_obs_below_vol = .false.

! Do the interpolation of pressure values only after taking the log (.true.)
! vs doing a linear interpolation directly in pressure units (.false.)
logical :: log_vert_interp  = .true.
logical :: log_horz_interpM = .false.
logical :: log_horz_interpQ = .false.

! wrf options, apply to domain 1 only
logical :: polar = .false.
logical :: periodic_x = .false.
logical :: periodic_y = .false.

logical :: allow_perturbed_ics = .false.
!-------------------------------

logical, parameter :: restrict_polar = .false. !HK what is this for? Hardcoded in original code

namelist /model_nml/ &
default_state_variables, &
wrf_state_variables, &
wrf_state_bounds, &
num_domains, &
calendar_type, &
assimilation_period_seconds, &
sfc_elev_max_diff, &
center_search_half_length, &
circulation_pres_level, &
circulation_radius, &
center_spline_grid_scale, &
vert_localization_coord, &
allow_perturbed_ics, &
allow_obs_below_vol, &
allow_obs_below_vol, &
log_horz_interpM, &
log_horz_interpQ

type grid_ll
  integer :: map_proj ! MAP_PROJ in wrf netcdf file
  type(proj_info) :: proj ! wrf map projection structure
  real(r8) :: dx
  real(r8) :: truelat1, truelat2, stand_lon
  real(r8), dimension(:,:),   allocatable :: latitude, latitude_u, latitude_v
  real(r8), dimension(:,:),   allocatable :: longitude, longitude_u, longitude_v
  integer :: we, sn ! west-east, south-north number of grid points
  integer :: wes, sns ! west-east staggered, south-north staggered number of grid points
  integer :: bt ! bottom-top number of grid points
  integer :: bts ! staggered bottom-top number of grid points
  real(r8) :: dt ! time step

  ! wrf options, apply to domain 1 only.
  logical :: polar = .false.
  logical :: periodic_x = .false.
  logical :: periodic_y = .false.

end type grid_ll

type static_data
   real(r8), allocatable :: phb(:,:,:) ! base-state geopotential
   real(r8), allocatable :: mub(:,:)   ! base state dry air mass in column
   real(r8), allocatable :: hgt(:,:)   ! Terrain Height
   real(r8), allocatable :: dnw(:)     ! d(eta) values between full (w) level
   real(r8), allocatable :: land(:,:)  ! land mask (1 for land, 2 for water)
   real(r8), allocatable :: zs(:)      ! depths of center of soil layers
   real(r8)              :: p_top      ! Pressure top of the model
end type static_data

! need grid for each domain
type(grid_ll), allocatable :: grid(:)
type(static_data), allocatable :: stat_dat(:)

! Physical constants
real(r8), parameter :: rd_over_rv = gas_constant / gas_constant_v
real(r8), parameter :: cpovcv = 1.4_r8        ! cp / (cp - gas_constant)
real(r8), parameter :: ts0 = 300.0_r8         ! Base potential temperature for all levels.
real(r8), parameter :: kappa = 2.0_r8/7.0_r8  ! gas_constant / cp

contains

!------------------------------------------------------------------
subroutine static_init_model()

integer  :: iunit, io

character(len=NF90_MAX_NAME) :: varname(MAX_STATE_VARIABLES)
integer :: state_qty(MAX_STATE_VARIABLES)
logical :: update_var(MAX_STATE_VARIABLES)
character(len=9) :: in_domain(MAX_STATE_VARIABLES) ! assumes <=9 or 999

integer :: nfields
logical, allocatable :: domain_mask(:)
integer :: i, field ! loop indices
character (len=1)     :: idom ! assumes <=9

module_initialized = .true.

call find_namelist_in_file("input.nml", "model_nml", iunit)
read(iunit, nml = model_nml, iostat = io)
call check_namelist_read(iunit, io, "model_nml")

! Record the namelist values used for the run 
if (do_nml_file()) write(nmlfileunit, nml=model_nml)
if (do_nml_term()) write(     *     , nml=model_nml)

call set_calendar_type(calendar_type)

allocate(wrf_dom(num_domains), grid(num_domains), stat_dat(num_domains))

call verify_state_variables(nfields, varname, state_qty, update_var, in_domain)

allocate(domain_mask(nfields))

do i = 1, num_domains

  do field = 1, nfields
     domain_mask(field) = variable_is_on_domain(in_domain(field), i)
  end do

  write( idom , '(I1)') i
  wrf_dom(i) = add_domain('wrfinput_d0'//idom, &
                          num_vars=count(domain_mask), &
                          var_names = pack(varname(1:nfields), domain_mask), &
                          kind_list = pack(state_qty(1:nfields), domain_mask), &
                          !clamp_vals  = &
                          update_list = pack(update_var(1:nfields), domain_mask) )
   
enddo

call read_grid()
call read_static_data()

call set_vertical_localization_coord(vert_localization_coord)
 
deallocate(domain_mask)

end subroutine static_init_model

!------------------------------------------------------------------
! Returns the number of items in the state vector as an i8 integer.
function get_model_size()

integer(i8) :: get_model_size
integer :: i

if ( .not. module_initialized ) call static_init_model

get_model_size = 0
do i = 1, num_domains
   get_model_size = get_model_size + get_domain_size(wrf_dom(i))
enddo

end function get_model_size

!------------------------------------------------------------------
subroutine model_interpolate(state_handle, ens_size, location, qty_in, expected_obs, istatus)


type(ensemble_type), intent(in) :: state_handle
integer,             intent(in) :: ens_size
type(location_type), intent(in) :: location
integer,             intent(in) :: qty_in
real(r8),           intent(out) :: expected_obs(ens_size) ! array of interpolated values
integer,            intent(out) :: istatus(ens_size)

integer, parameter :: FAILED_BOUNDS_CHECK = 44
integer, parameter :: CANNOT_INTERPOLATE_QTY = 55
integer, parameter :: POLAR_RESTRICTED = 10 ! polar observation while restrict_polar = .true.
integer, parameter :: NOT_IN_ANY_DOMAIN = 11
integer, parameter :: VERTICAL_LOCATION_FAIL = 66
real(r8) :: lon_lat_vert(3)
real(r8) :: xloc, yloc ! WRF i,j in the grid
integer  :: i, j ! grid
real(r8) :: dx, dxm, dy, dym ! grid fractions
integer :: ll(2), ul(2), lr(2), ur(2) !(x,y) of four corners
integer :: rc ! return code getCorners
integer :: id
integer :: k(ens_size)      ! level
integer :: which_vert       ! vertical coordinate of the observation
real(r8) :: zloc(ens_size)  ! vertical location of the obs for each ens member
real(r8) :: fld_k1(ens_size), fld_k2(ens_size) ! value at level k and k+1
real(r8) :: fld1(ens_size), fld2(ens_size) ! value at level k and k+1 non-negative if required
logical :: fail
integer :: qty

if ( .not. module_initialized ) call static_init_model

expected_obs(:) = MISSING_R8

istatus(:) = 1

which_vert = nint(query_location(location))
lon_lat_vert = get_location(location)
call get_domain_info(lon_lat_vert(1),lon_lat_vert(2),id,xloc,yloc) ! mass points

if (id == 0) then
   istatus(:) = NOT_IN_ANY_DOMAIN
   return
endif

if (.not. able_to_interpolate_qty(id, qty_in) ) then
   istatus(:) = CANNOT_INTERPOLATE_QTY
   return
endif

qty = update_qty_if_location_is_surface(qty_in, location)

! horizontal location mass point
call toGrid(xloc,i,dx,dxm)
call toGrid(yloc,j,dy,dym)

if ( .not. within_bounds_horizontal(i, j, id, qty) ) then
   istatus(:) = FAILED_BOUNDS_CHECK
   return
endif

call getCorners(i, j, id, qty, ll, ul, lr, ur, rc)

! vertical location
call get_level_below_obs(which_vert, id, lon_lat_vert, ens_size, state_handle, ll, ul, lr, ur, dx, dy, dxm, dym, k, zloc, fail)
if (fail) then
   istatus(:) = VERTICAL_LOCATION_FAIL
   return
endif


select case (qty)
   case (QTY_U_WIND_COMPONENT, QTY_V_WIND_COMPONENT )
      fld_k1 = wind_interpolate(ens_size, state_handle, qty, id, k,   xloc, yloc, i, j, dxm, dx, dy, dym, lon_lat_vert(1))
      fld_k2 = wind_interpolate(ens_size, state_handle, qty, id, k+1, xloc, yloc, i, j, dxm, dx, dy, dym, lon_lat_vert(1))
   case (QTY_TEMPERATURE)
      fld_k1 = temperature_interpolate(ens_size, state_handle, id, ll, ul, lr, ur, k, dxm, dx, dy, dym)
      fld_k2 = temperature_interpolate(ens_size, state_handle, id, ll, ul, lr, ur, k+1, dxm, dx, dy, dym)
   case (QTY_POTENTIAL_TEMPERATURE)
      fld_k1 = simple_interpolation(ens_size, state_handle, qty, id, ll, ul, lr, ur, k, dxm, dx, dy, dym) + ts0
      fld_k2 = simple_interpolation(ens_size, state_handle, qty, id, ll, ul, lr, ur, k+1, dxm, dx, dy, dym) + ts0
   case (QTY_DENSITY)
      fld_k1 = density_interpolate(ens_size, state_handle, id, ll, ul, lr, ur, k, dxm, dx, dy, dym)
      fld_k2 = density_interpolate(ens_size, state_handle, id, ll, ul, lr, ur, k+1, dxm, dx, dy, dym)
   case (QTY_VERTICAL_VELOCITY)
      zloc(:) = zloc(:) + 0.5_r8 ! Adjust zloc for staggered
      k(:) = max(1,int(zloc(:)))  ! Adjust corresponding level k
      fld_k1(:) = simple_interpolation(ens_size, state_handle, qty, id, ll, ul, lr, ur, k, dxm, dx, dy, dym) 
      fld_k2(:) = simple_interpolation(ens_size, state_handle, qty, id, ll, ul, lr, ur, k+1, dxm, dx, dy, dym) 
   case (QTY_SPECIFIC_HUMIDITY)
      fld_k1 = specific_humidity_interpolate(ens_size, state_handle, qty, id, ll, ul, lr, ur, k, dxm, dx, dy, dym)
      fld_k2 = specific_humidity_interpolate(ens_size, state_handle, qty, id, ll, ul, lr, ur, k+1, dxm, dx, dy, dym)
   case (QTY_VAPOR_MIXING_RATIO)
      fld_k1(:) = surface_type_interpolate(ens_size, id, ll, ul, lr, ur, dxm, dx, dy, dym)
   case (QTY_PRESSURE)
      fld_k1 = pressure_interpolate(ens_size, state_handle, id, ll, ul, lr, ur, k, dxm, dx, dy, dym)
      fld_k2 = pressure_interpolate(ens_size, state_handle, id, ll, ul, lr, ur, k+1, dxm, dx, dy, dym)
   case (QTY_SURFACE_PRESSURE)
      fld_k1 = surface_pressure_interpolate(ens_size, state_handle, id, ll, ul, lr, ur, dxm, dx, dy, dym)
   case (QTY_VORTEX_LAT, QTY_VORTEX_LON, QTY_VORTEX_PMIN, QTY_VORTEX_WMAX)
      call vortex()
   case (QTY_GEOPOTENTIAL_HEIGHT)
      zloc(:) = zloc(:) + 0.5_r8 ! Adjust zloc for staggered
      k(:) = max(1,int(zloc(:)))  ! Adjust corresponding level k
      fld_k1 = geopotential_height_interpolate(ens_size, state_handle, qty, id, ll, ul, lr, ur, k, dxm, dx, dy, dym)
      fld_k2 = geopotential_height_interpolate(ens_size, state_handle, qty, id, ll, ul, lr, ur, k+1, dxm, dx, dy, dym)
   case (QTY_SURFACE_ELEVATION)
      fld_k1 = surface_elevation_interpolate(ens_size, id, ll, ul, lr, ur, dxm, dx, dy, dym)
   case (QTY_SURFACE_TYPE)
      fld_k1(:) = surface_type_interpolate(ens_size, id, ll, ul, lr, ur, dxm, dx, dy, dym)
   case (QTY_SKIN_TEMPERATURE, QTY_10M_U_WIND_COMPONENT, QTY_10M_V_WIND_COMPONENT, QTY_2M_TEMPERATURE, QTY_2M_SPECIFIC_HUMIDITY)
      fld_k1(:) = simple_interpolation(ens_size, state_handle, qty, id, ll, ul, lr, ur, k, dxm, dx, dy, dym)
   case default ! simple interpolation
      ! HK todo 2D variables
      fld_k1(:) = simple_interpolation(ens_size, state_handle, qty, id, ll, ul, lr, ur, k, dxm, dx, dy, dym)
      fld_k2(:) = simple_interpolation(ens_size, state_handle, qty, id, ll, ul, lr, ur, k+1, dxm, dx, dy, dym)
end select

if (surface_qty(qty)) then
  expected_obs(:) = fld_k1(:)
else
   fld1  = force_non_negative_if_required(ens_size, qty, fld_k1)
   fld2  = force_non_negative_if_required(ens_size, qty, fld_k2)
   expected_obs(:) = vertical_interpolation(ens_size, zloc, fld1, fld2)
endif

istatus(:) = 0

end subroutine model_interpolate

!------------------------------------------------------------------
! Returns the smallest increment in time that the model is capable 
! of advancing the state in a given implementation, or the shortest
! time you want the model to advance between assimilations.
function shortest_time_between_assimilations()

type(time_type) :: shortest_time_between_assimilations
integer :: model_dt

if ( .not. module_initialized ) call static_init_model

model_dt = nint(grid(1)%dt) ! model time step in seconds lowest res domain

if (assimilation_period_seconds < model_dt ) then
   shortest_time_between_assimilations = set_time(model_dt)
else
   shortest_time_between_assimilations = set_time(assimilation_period_seconds)
endif

end function shortest_time_between_assimilations

!------------------------------------------------------------------
! Given an integer index into the state vector, returns the
! associated location and optionally the physical quantity.
subroutine get_state_meta_data(index_in, location, qty_out)

integer(i8),         intent(in)  :: index_in
type(location_type), intent(out) :: location
integer, optional,   intent(out) :: qty_out

integer :: i, j, k, id, var_id, state_id, qty

if ( .not. module_initialized ) call static_init_model

! wrf domain may not equal state_id
call get_model_variable_indices(index_in, i, j, k, var_id=var_id, dom_id=state_id, kind_index=qty)

location = convert_indices_to_lon_lat_lev(i, j, k, var_id, state_id)

! return DART variable qty if requested
if(present(qty_out)) qty_out = qty

end subroutine get_state_meta_data

!------------------------------------------------------------------
! observations have a type and qty
!  observation type not taken in to account for wrf get close calculations
subroutine get_close_obs(gc, base_loc, base_type, locs, loc_qtys, loc_types, &
                         num_close, close_ind, dist, state_handle)

type(get_close_type),          intent(in)    :: gc            ! handle to a get_close structure
integer,                       intent(in)    :: base_type     ! observation TYPE
type(location_type),           intent(inout) :: base_loc      ! location of interest
type(location_type),           intent(inout) :: locs(:)       ! obs/state locations
integer,                       intent(in)    :: loc_qtys(:)   ! QTYS for obs
integer,                       intent(in)    :: loc_types(:)  ! types for obs
integer,                       intent(out)   :: num_close     ! how many are close
integer,                       intent(out)   :: close_ind(:)  ! indices into the locs array
real(r8),            optional, intent(out)   :: dist(:)       ! distances in radians
type(ensemble_type), optional, intent(in)    :: state_handle

character(len=*), parameter :: routine = 'get_close_obs'
integer :: istatus(1), loc_qtys_ar(1), loc_types_ar(1), i, t_ind
type(location_type) :: loc_ar(1)

if (vertical_localization_on()) then
   loc_ar(1) = base_loc
   ! dummy qty 1 not used in convert_vertical_obs
   call convert_vertical_obs(state_handle, 1, loc_ar, (/1/), (/base_type/), vert_localization_coord, istatus)
   if (istatus(1) /= 0) then
      num_close = 0
      return
   endif
   base_loc = loc_ar(1)
endif

call loc_get_close(gc, base_loc, base_type, locs, loc_qtys, &
                          num_close, close_ind)

if (.not. present(dist)) return

do i = 1, num_close
   t_ind = close_ind(i)
   loc_ar(1) = locs(t_ind)
   loc_qtys_ar(1) = loc_qtys(t_ind)   ! HK not used in convert_vertical_obs todo: use dummy instead?
   loc_types_ar(1) = loc_types(t_ind) ! HK not used in convert_vertical_obs todo: use dummy instead?

   istatus(1) = 0 
   if (vertical_localization_on()) then
     call convert_vertical_obs(state_handle, 1, loc_ar, loc_qtys_ar, loc_types_ar, vert_localization_coord, istatus)
   endif
   if (istatus(1) == 0) then
      dist(i) = get_dist(base_loc, loc_ar(1), base_type, loc_qtys(t_ind))
   else
      dist(i) = 1.0e9
   endif
enddo


end subroutine get_close_obs

!------------------------------------------------------------------
! state only has qty
subroutine get_close_state(gc, base_loc, base_type, locs, loc_qtys, loc_indx, &
                           num_close, close_ind, dist, state_handle)

type(get_close_type),          intent(in)     :: gc
type(location_type),           intent(inout)  :: base_loc, locs(:)
integer,                       intent(in)     :: base_type, loc_qtys(:)
integer(i8),                   intent(in)     :: loc_indx(:)
integer,                       intent(out)    :: num_close, close_ind(:)
real(r8),            optional, intent(out)    :: dist(:)
type(ensemble_type), optional, intent(in)     :: state_handle

character(len=*), parameter :: routine = 'get_close_state'
integer :: istatus(1), loc_qtys_ar(1), i, t_ind
integer(i8) :: loc_indx_ar(1)
type(location_type) :: loc_ar(1)

loc_ar(1) = base_loc
! dummy qty 1 not used in convert_vertical_obs
call convert_vertical_obs(state_handle, 1, loc_ar, (/1/), (/base_type/), vert_localization_coord, istatus)
if (istatus(1) /= 0) then
   num_close = 0
   return
endif
base_loc = loc_ar(1)

call loc_get_close(gc, base_loc, base_type, locs, loc_qtys, num_close, close_ind)

if (.not. present(dist)) return

do i = 1, num_close
   t_ind = close_ind(i)
   loc_ar(1) = locs(t_ind)
   loc_qtys_ar(1) = loc_qtys(t_ind)   ! HK not used in convert_vertical_state todo: use dummy instead?
   loc_indx_ar(1) = loc_indx(t_ind)

   if (vertical_localization_on()) then
       if (nint(query_location(loc_ar(1))) /= vert_localization_coord) then
          ! convert_vertical_state always returns istatus = 0
          call convert_vertical_state(state_handle, 1, loc_ar, loc_qtys_ar, loc_indx_ar, vert_localization_coord, istatus(1))
       endif
   endif
       
   dist(i) = get_dist(base_loc, loc_ar(1), base_type, loc_qtys(t_ind))
 
enddo

end subroutine get_close_state


!------------------------------------------------------------------
! write any additional attributes to netcdf files
! HK todo nc_write_model_atts
subroutine nc_write_model_atts(ncid, domain_id)

integer, intent(in) :: ncid      ! netCDF file identifier
integer, intent(in) :: domain_id

if ( .not. module_initialized ) call static_init_model

! put file into define mode.

call nc_begin_define_mode(ncid)

call nc_add_global_creation_time(ncid)

call nc_add_global_attribute(ncid, "model_source", source )
call nc_add_global_attribute(ncid, "model", "template")

call nc_end_define_mode(ncid)

! Flush the buffer and leave netCDF file open
call nc_synchronize_file(ncid)

end subroutine nc_write_model_atts

!------------------------------------------------------------------
subroutine pert_model_copies(ens_handle, ens_size,  dummy_pert_amp, interf_provided)

type(ensemble_type), intent(inout) :: ens_handle
integer,             intent(in)    :: ens_size
real(r8),            intent(in)    :: dummy_pert_amp ! not used
logical,             intent(out)   :: interf_provided

if (.not. allow_perturbed_ics) then
call error_handler(E_ERR,'pert_model_copies', &
                     'starting WRF model from a single vector requires additional steps', &
                      text2='see comments in wrf/model_mod.f90::pert_model_copies()')
endif

print*, 'not done'

end subroutine pert_model_copies

!------------------------------------------------------------------
function read_model_time(filename)

character(len=*),  intent(in) :: filename
type(time_type)               :: read_model_time

character(len=*), parameter :: routine = 'read_model_time'
integer :: ncid, ret, var_id
integer :: dim_size(2) ! wrf netcdf: char Times(Time, DateStrLen)
character(len=19) :: timestring ! e.g. 2007-04-26_00:00:00
integer           :: year, month, day, hour, minute, second

ncid = nc_open_file_readonly(filename, routine)

call nc_get_variable_size(ncid, 'Times', dim_size, routine)

ret = nf90_inq_varid(ncid, "Times", var_id)
call nc_check(ret, routine, 'inq_varid Times')

! last slice of Time dimension
ret = nf90_get_var(ncid, var_id, timestring, start = (/ 1, dim_size(2) /))
call nc_check(ret, routine, 'get_var Times')

call get_wrf_date(timestring, year, month, day, hour, minute, second)
read_model_time = set_date(year, month, day, hour, minute, second)

call  nc_close_file(ncid, routine)

end function read_model_time

!------------------------------------------------------------------
subroutine end_model()

deallocate(wrf_dom, grid, stat_dat)

end subroutine end_model

!------------------------------------------------------------------
! end of public routines
!------------------------------------------------------------------
function force_non_negative_if_required(n, qty, fld)

integer,  intent(in) :: n
integer,  intent(in) :: qty
real(r8), intent(in) :: fld(n)

real(r8) :: force_non_negative_if_required(n)

select case (qty)
   case (QTY_RAINWATER_MIXING_RATIO,    &
         QTY_GRAUPEL_MIXING_RATIO,      &
         QTY_HAIL_MIXING_RATIO,         &
         QTY_SNOW_MIXING_RATIO,         &
         QTY_ICE_MIXING_RATIO,          &
         QTY_CLOUDWATER_MIXING_RATIO,   &
         QTY_DROPLET_NUMBER_CONCENTR,   &
         QTY_ICE_NUMBER_CONCENTRATION,  &
         QTY_SNOW_NUMBER_CONCENTR,      &
         QTY_RAIN_NUMBER_CONCENTR,      &
         QTY_GRAUPEL_NUMBER_CONCENTR,   &
         QTY_HAIL_NUMBER_CONCENTR )
      force_non_negative_if_required = max(0.0_r8, fld)
   case default
      force_non_negative_if_required = fld
end select

end function force_non_negative_if_required

!------------------------------------------------------------------
! 2D surface variables, or soil variables (z = soil_layers_stag)
function surface_qty(qty)

integer, intent(in) :: qty
logical :: surface_qty

!var_id get_varid_from_kind(wrf_dom(id), qty)
! HK todo soil_layers_sta
!do i = 1, get_num_dims(var_id)
!if () then
!get_dim_name
!endif

select case (qty)

   case (QTY_2M_TEMPERATURE, &
         QTY_2M_SPECIFIC_HUMIDITY, &
         QTY_10M_U_WIND_COMPONENT, &
         QTY_10M_V_WIND_COMPONENT, &
         QTY_SURFACE_TYPE, &
         QTY_SKIN_TEMPERATURE)
      surface_qty = .true.
   case default
      surface_qty = .false.

end select

end function surface_qty

!------------------------------------------------------------------
! WRF has separate model qtys for surface variables
function update_qty_if_location_is_surface(qty_in, location) result(qty)

integer, intent(in) :: qty_in
type(location_type), intent(in) :: location
integer :: qty

if (.not. is_vertical(location,"SURFACE")) then
   qty = qty_in
   return
endif

select case (qty)

   case (QTY_U_WIND_COMPONENT); qty = QTY_10M_U_WIND_COMPONENT
   case (QTY_V_WIND_COMPONENT); qty = QTY_10M_V_WIND_COMPONENT
   case (QTY_POTENTIAL_TEMPERATURE); qty = QTY_2M_TEMPERATURE
   case (QTY_SPECIFIC_HUMIDITY); qty = QTY_2M_SPECIFIC_HUMIDITY
   case (QTY_VAPOR_MIXING_RATIO); qty = QTY_2M_SPECIFIC_HUMIDITY  ! Vapor Mixing Ratio (QV, Q2)
   case (QTY_PRESSURE); qty = QTY_SURFACE_PRESSURE

end select

end function update_qty_if_location_is_surface

!------------------------------------------------------------------
function simple_interpolation(ens_size, state_handle, qty, id, ll, ul, lr, ur, k, dxm, dx, dy, dym)

integer,             intent(in) :: ens_size
type(ensemble_type), intent(in) :: state_handle
integer,             intent(in) :: qty
integer,             intent(in) :: id
integer,             intent(in) :: ll(2), ul(2), lr(2), ur(2) ! (x,y) at  four corners
integer,             intent(in) :: k(ens_size) ! k may be different across the ensemble
real(r8),            intent(in) :: dxm, dx, dy, dym

real(r8) ::  simple_interpolation(ens_size)

integer :: e ! loop variable
! lower left, upper left, lower right, upper right
integer(i8), dimension(ens_size)  :: ill, iul, ilr, iur  ! dart index at four corners
real(r8),    dimension(ens_size)  :: x_ill, x_iul, x_ilr, x_iur ! state value at four corners
integer :: var_id


var_id = get_varid_from_kind(wrf_dom(id), qty)

do e = 1, ens_size
                              !   x,     y,     z,    domain,      variable
   ill(e) = get_dart_vector_index(ll(1), ll(2), k(e), wrf_dom(id), var_id)
   iul(e) = get_dart_vector_index(ul(1), ul(2), k(e), wrf_dom(id), var_id)
   ilr(e) = get_dart_vector_index(lr(1), lr(2), k(e), wrf_dom(id), var_id)
   iur(e) = get_dart_vector_index(ur(1), ur(2), k(e), wrf_dom(id), var_id)

enddo

call get_state_array(x_ill, ill, state_handle)
call get_state_array(x_iul, iul, state_handle)
call get_state_array(x_ilr, ilr, state_handle)
call get_state_array(x_iur, iur, state_handle)

simple_interpolation = dym*( dxm*x_ill(:) + dx*x_ilr(:) ) + dy*( dxm*x_iul(:) + dx*x_iur(:) )

end function simple_interpolation

!------------------------------------------------------------------
function vertical_interpolation(ens_size, zloc, fld1, fld2)

integer,  intent(in) :: ens_size
real(r8), intent(in) :: zloc(ens_size)
real(r8), intent(in) :: fld1(ens_size), fld2(ens_size)

real(r8) :: vertical_interpolation(ens_size)

real(r8) :: dz, dzm
integer  :: z ! level
integer :: e

do e = 1, ens_size
   call toGrid(zloc(e), z, dz, dzm)
   ! comment from original code:
   ! If you get here and zloc < 1.0, then z will be 0, and
   ! we should extrapolate.  fld(1,:) and fld(2,:) where computed
   ! at levels 1 and 2.

   if (z >= 1) then
      ! Linearly interpolate between levels
      vertical_interpolation(e) = dzm*fld1(e) + dz*fld2(e)
   else
      ! Extrapolate below first level.
      vertical_interpolation(e) = fld1(e) - (fld2(e)-fld1(e))*dzm
   endif
enddo

end function vertical_interpolation

!------------------------------------------------------------------
subroutine get_wrf_date(tstring, year, month, day, hour, minute, second)

character(len=19), intent(in)  :: tstring ! YYYY-MM-DD_hh:mm:ss
integer,           intent(out) :: year, month, day, hour, minute, second

read(tstring( 1: 4),'(i4)') year
read(tstring( 6: 7),'(i2)') month
read(tstring( 9:10),'(i2)') day
read(tstring(12:13),'(i2)') hour
read(tstring(15:16),'(i2)') minute
read(tstring(18:19),'(i2)') second

end subroutine get_wrf_date

!------------------------------------------------------------------
! Verify that the namelist was filled in correctly, and check
! that there are valid entries for the dart_kind.
subroutine verify_state_variables(nvar, varname, qty, update, in_domain)

integer,           intent(out) :: nvar
character(len=NF90_MAX_NAME), intent(out) :: varname(MAX_STATE_VARIABLES)
integer,           intent(out) :: qty(MAX_STATE_VARIABLES)
logical,           intent(out) :: update(MAX_STATE_VARIABLES)
character(len=9),  intent(out) :: in_domain(MAX_STATE_VARIABLES) ! assumes <=9 or 999

integer :: i
character(len=NF90_MAX_NAME) :: qty_str, update_str
character(len=256) :: string1, string2

if ( .not. module_initialized ) call static_init_model

nvar = 0
varloop: do i = 1, MAX_STATE_VARIABLES

   if ( wrf_state_variables(NUM_STATE_TABLE_COLUMNS*i -3) == 'NULL ' ) exit varloop ! Found end of list. !HK do you need to test all columns for ' '?

   varname(i)   = trim(wrf_state_variables(NUM_STATE_TABLE_COLUMNS*i -3))
   qty_str      = trim(wrf_state_variables(NUM_STATE_TABLE_COLUMNS*i -2))
   update_str   = trim(wrf_state_variables(NUM_STATE_TABLE_COLUMNS*i -1))
   in_domain(i) = trim(wrf_state_variables(NUM_STATE_TABLE_COLUMNS*i   ))

   call to_upper(update_str)
   

   ! Make sure DART qty is valid
   qty(i) = get_index_for_quantity(qty_str)
   if( qty(i)  < 0 ) then
      write(string1,'(''there is no obs_kind <'',a,''> in obs_kind_mod.f90'')') trim(qty_str)
      call error_handler(E_ERR,'verify_state_variables',string1)
   endif
   
   ! Force QTY_TEMPERATURE to  QTY_POTENTIAL_TEMPERATURE
   if (qty(i) == QTY_TEMPERATURE) qty(i) = QTY_POTENTIAL_TEMPERATURE
   ! Make sure the update variable has a valid name
   select case (update_str)
      case ('UPDATE')
         update(i) = .true.
      case ('NO_COPY_BACK')
         update(i) = .false.
      case default
         write(string1,'(A)')  'only UPDATE or NO_COPY_BACK supported in model_state_variable namelist'
         write(string2,'(6A)') 'you provided : ', trim(varname(i)), ', ', trim(qty_str), ', ', trim(update_str)
         call error_handler(E_ERR,'verify_state_variables',string1, text2=string2)
   end select
   
   nvar = nvar + 1

enddo varloop

end subroutine verify_state_variables

!------------------------------------------------------------------
! This is assuming that the number of domains <=9
! do while loop could be replaced with intrinsic scan
logical function variable_is_on_domain(domain_id_string, id)

integer,           intent(in) :: id
character(len=*),  intent(in) :: domain_id_string

integer                       :: domain_int, i

variable_is_on_domain = .false.

! if '999' then counts all domains
if ( trim(domain_id_string) == '999' ) then
   variable_is_on_domain = .true.
else
i = 1
   do while ( domain_id_string(i:i) /= ' ' )
   read(domain_id_string(i:i),'(i1)') domain_int
   if ( domain_int == id ) variable_is_on_domain = .true.
   i = i+1
   enddo
endif

end function variable_is_on_domain

!------------------------------------------------------------------
subroutine read_grid()

integer :: ncid, i
character (len=1) :: idom ! assumes <=9
character(len=*), parameter :: routine = 'read_grid'

integer :: dim_size(3) ! (west_east, south_north, Time)

! This is assuming the unlimited dimension is length 1
! Should we be reading the latest time slice instead?

do i = 1, num_domains

    write( idom , '(I1)') i
    ncid = nc_open_file_readonly('wrfinput_d0'//idom, routine)

    call nc_get_variable_size(ncid, 'XLONG', dim_size)
    allocate(grid(i)%longitude(dim_size(1), dim_size(2)))
    call nc_get_variable(ncid, 'XLONG', grid(i)%longitude, routine)
    grid(i)%we = dim_size(1); grid(i)%sn = dim_size(2)

    call nc_get_variable_size(ncid, 'XLONG_U', dim_size)
    allocate(grid(i)%longitude_u(dim_size(1), dim_size(2)))
    call nc_get_variable(ncid, 'XLONG_U', grid(i)%longitude_u, routine)
    grid(i)%wes = dim_size(1)
    
    call nc_get_variable_size(ncid, 'XLONG_V', dim_size)
    allocate(grid(i)%longitude_v(dim_size(1), dim_size(2)))
    call nc_get_variable(ncid, 'XLONG_V', grid(i)%longitude_v, routine)
    grid(i)%sns = dim_size(2)   
 
    call nc_get_variable_size(ncid, 'XLAT', dim_size)
    allocate(grid(i)%latitude(dim_size(1), dim_size(2)))
    call nc_get_variable(ncid, 'XLAT', grid(i)%latitude, routine)
    
    call nc_get_variable_size(ncid, 'XLAT_U', dim_size)
    allocate(grid(i)%latitude_u(dim_size(1), dim_size(2)))
    call nc_get_variable(ncid, 'XLAT_U', grid(i)%latitude_u, routine)
    
    call nc_get_variable_size(ncid, 'XLAT_V', dim_size)
    allocate(grid(i)%latitude_v(dim_size(1), dim_size(2)))
    call nc_get_variable(ncid, 'XLAT_V', grid(i)%latitude_v, routine)

    call nc_get_global_attribute(ncid, 'MAP_PROJ', grid(i)%map_proj)
    call nc_get_global_attribute(ncid, 'DX', grid(i)%dx)
    call nc_get_global_attribute(ncid, 'TRUELAT1', grid(i)%truelat1)
    call nc_get_global_attribute(ncid, 'TRUELAT2', grid(i)%truelat2)
    call nc_get_global_attribute(ncid, 'STAND_LON', grid(i)%stand_lon)
    call nc_get_global_attribute(ncid, 'DT', grid(i)%dt)

    grid(i)%bt = nc_get_dimension_size(ncid, 'bottom_top', routine)
    grid(i)%bts = nc_get_dimension_size(ncid, 'bottom_top_stag', routine)

    call nc_close_file(ncid, routine)

    call setup_map_projection(i)

    if (i == 1) then
       grid(i)%periodic_x = periodic_x
       grid(i)%periodic_y = periodic_y
       grid(i)%polar = polar
    else
       grid(i)%periodic_x = .false.
       grid(i)%periodic_y = .false.
       grid(i)%polar = .false.
    endif

enddo

end subroutine read_grid

!------------------------------------------------------------------
subroutine read_static_data()

integer :: ncid, i
character (len=1) :: idom ! assumes <=9
character(len=*), parameter :: routine = 'read_static_data'

integer :: dim_size(4) ! (west_east, south_north, bottom_top{_stag}, Time)

do i = 1, num_domains

   write( idom , '(I1)') i
   ncid = nc_open_file_readonly('wrfinput_d0'//idom, routine)

   call nc_get_variable_size(ncid, 'PHB', dim_size)
   allocate(stat_dat(i)%phb(dim_size(1), dim_size(2), dim_size(3)))
   call nc_get_variable(ncid, 'PHB', stat_dat(i)%phb, routine)
 
   call nc_get_variable_size(ncid, 'MUB', dim_size)
   allocate(stat_dat(i)%mub(dim_size(1), dim_size(2)))
   call nc_get_variable(ncid, 'MUB', stat_dat(i)%mub, routine)

   call nc_get_variable_size(ncid, 'HGT', dim_size)
   allocate(stat_dat(i)%hgt(dim_size(1), dim_size(2)))
   call nc_get_variable(ncid, 'HGT', stat_dat(i)%hgt, routine)

   call nc_get_variable_size(ncid, 'DNW', dim_size)
   allocate(stat_dat(i)%dnw(dim_size(1)))
   call nc_get_variable(ncid, 'DNW', stat_dat(i)%dnw, routine)

   call nc_get_variable_size(ncid, 'XLAND', dim_size)
   allocate(stat_dat(i)%land(dim_size(1), dim_size(2)))
   call nc_get_variable(ncid, 'XLAND', stat_dat(i)%land, routine)

   call nc_get_variable_size(ncid, 'ZS', dim_size)
   allocate(stat_dat(i)%zs(dim_size(1))) ! soil_layers_stag
   call nc_get_variable(ncid, 'ZS', stat_dat(i)%zs, routine)

   call nc_get_variable(ncid, 'P_TOP', stat_dat(i)%p_top, routine)


   call nc_close_file(ncid, routine)

end do

end subroutine read_static_data

!------------------------------------------------------------------
pure function compute_geometric_height(geopot, lat)

real(r8), intent(in)  :: geopot
real(r8), intent(in)  :: lat
real(r8)              :: compute_geometric_height


real(digits12) :: pi2, latr
real(digits12) :: semi_major_axis, semi_minor_axis, grav_polar, grav_equator
real(digits12) :: earth_omega, grav_constant, flattening, somigliana
real(digits12) :: grav_ratio, sin2, termg, termr, grav, eccentricity

!  Parameters below from WGS-84 model software inside GPS receivers.
parameter(semi_major_axis = 6378.1370d3)    ! (m)
parameter(semi_minor_axis = 6356.7523142d3) ! (m)
parameter(grav_polar = 9.8321849378)        ! (m/s2)
parameter(grav_equator = 9.7803253359)      ! (m/s2)
parameter(earth_omega = 7.292115d-5)        ! (rad/s)
parameter(grav_constant = 3.986004418d14)   ! (m3/s2)
parameter(grav = 9.80665d0)                 ! (m/s2) WMO std g at 45 deg lat
parameter(eccentricity = 0.081819d0)        ! unitless
parameter(pi2 = 3.14159265358979d0/180.d0)

!  Derived geophysical constants
parameter(flattening = (semi_major_axis-semi_minor_axis) / semi_major_axis)

parameter(somigliana = (semi_minor_axis/semi_major_axis)*(grav_polar/grav_equator)-1.d0)

parameter(grav_ratio = (earth_omega*earth_omega * &
                semi_major_axis*semi_major_axis * semi_minor_axis)/grav_constant)

! To use geopotential height uncomment the following two lines:
!compute_geometric_height = geopot
!return

latr = lat * (pi2)        ! in radians
sin2  = sin(latr) * sin(latr)
termg = grav_equator * ( (1.d0+somigliana*sin2) / &
        sqrt(1.d0-eccentricity*eccentricity*sin2) )
termr = semi_major_axis / (1.d0 + flattening + grav_ratio - 2.d0*flattening*sin2)

compute_geometric_height = (termr*geopot) / ( (termg/grav) * termr - geopot )


end function compute_geometric_height

!------------------------------------------------------------------
! This is mass level
subroutine get_level_below_obs(which_vert, id, lon_lat_vert, ens_size, state_handle, &
                               ll, ul, lr, ur, dx, dy, dxm, dym, &
                               level_below, zloc, fail)

integer,  intent(in)  :: which_vert
integer,  intent(in)  :: id
real(r8), intent(in)  :: lon_lat_vert(3)
integer,  intent(in)  :: ens_size
type(ensemble_type), intent(in) :: state_handle
integer, dimension(2), intent(in) :: ll, ul, lr, ur ! (x,y) of each corner
real(r8), intent(in)  :: dx, dxm, dy, dym ! grid fractions to obs
integer,  intent(out) :: level_below(ens_size)
real(r8), intent(out) :: zloc(ens_size) ! vertical location of the obs for each ens member
logical,  intent(out) :: fail

integer :: e ! loop variable
real(r8) :: v_p(0:grid(id)%bt,ens_size)
real(r8) :: v_h(0:grid(id)%bt,ens_size)
logical :: lev0

fail = .false.

select case (which_vert)
   case(VERTISLEVEL)
      zloc(:) = lon_lat_vert(3); fail = .false.
   case(VERTISPRESSURE)
      call get_model_pressure_profile(id, ll, ul, lr, ur, dx, dy, dxm, dym, ens_size, state_handle, v_p)
      do e = 1, ens_size
         call pres_to_zk(lon_lat_vert(3), v_p(:,e), grid(id)%bt, zloc(e), level_below(e), lev0, fail)
         if (fail) return
         if (lev0) then
            print*, "pressure obs below lowest sigma"
            fail = .true.
            return
         endif
      enddo
   case(VERTISHEIGHT)
      call get_model_height_profile(ll, ul, lr, ur, dx, dy, dxm, dym, id, v_h, state_handle, ens_size)
      do e = 1, ens_size
         call height_to_zk(lon_lat_vert(3), v_h(:, e), grid(id)%bt, zloc(e), level_below(e), lev0, fail)
         if (fail) return
         if (lev0) then
            print*, "height obs below lowest sigma"
               fail = .true.
            return
         endif
      enddo
   case(VERTISSURFACE)
       zloc(:) = 1.0_r8
       ! HK todo call check to see if the station height is too far away from the model surface height
   case(VERTISUNDEF)
       zloc = 0.0_r8
   case default
       fail = .true.
end select

end subroutine get_level_below_obs

!------------------------------------------------------------------
! Calculate the model pressure profile on half (mass) levels,
! horizontally interpolated at the observation location.
subroutine get_model_pressure_profile(id, ll, ul, lr, ur, dx, dy, dxm, dym, &
                                      ens_size, state_handle, v_p)

integer,  intent(in)             :: id
integer,  intent(in)             :: ll(2), ul(2), lr(2), ur(2) ! (x,y) mass grid corners
real(r8), intent(in)             :: dx, dy, dxm, dym
integer,  intent(in)             :: ens_size
type(ensemble_type), intent(in)  :: state_handle
real(r8),            intent(out) :: v_p(0:grid(id)%bt, ens_size)

integer(i8)           :: ill, ilr, iul, iur
real(r8), dimension(ens_size) :: x_ill, x_ilr, x_iul, x_iur
real(r8), dimension(ens_size) :: pres1, pres2, pres3, pres4
real(r8), dimension(ens_size) :: lev2_pres1, lev2_pres2, lev2_pres3, lev2_pres4

integer :: var_id, levk
integer :: k(ens_size)

do levk=1, grid(id)%bt ! number of mass levels

   k(:) = levk
   pres1 = model_pressure_t(ll(1), ll(2), k, id, state_handle, ens_size)
   pres2 = model_pressure_t(lr(1), lr(2), k, id, state_handle, ens_size)
   pres3 = model_pressure_t(ul(1), ul(2), k, id, state_handle, ens_size)
   pres4 = model_pressure_t(ur(1), ur(2), k, id, state_handle, ens_size)

   v_p(levk, :) = interp_4pressure(pres1, pres2, pres3, pres4, dx, dxm, dy, dym, ens_size)

   if (levk == 2) then ! store result for extrapolation
     lev2_pres1(:) = pres1(:)
     lev2_pres2(:) = pres2(:)
     lev2_pres3(:) = pres3(:)
     lev2_pres4(:) = pres4(:)
   endif

enddo

var_id = get_varid_from_kind(wrf_dom(id), QTY_SURFACE_PRESSURE)

if (var_id > 0) then ! surface pressure in domain so get v_p(0,:) from surface pressure

   ill = get_dart_vector_index(ll(1), ll(2), 1, wrf_dom(id), var_id)
   ilr = get_dart_vector_index(lr(1), lr(2), 1, wrf_dom(id), var_id)
   iul = get_dart_vector_index(ul(1), ul(2), 1, wrf_dom(id), var_id)
   iur = get_dart_vector_index(ur(1), ur(2), 1, wrf_dom(id), var_id)

   x_ill = get_state(ill, state_handle)
   x_ilr = get_state(ilr, state_handle)
   x_iul = get_state(iul, state_handle)
   x_iur = get_state(iur, state_handle)

   v_p(0,:) = interp_4pressure(x_ill, x_ilr, x_iul, x_iur, dx, dxm, dy, dym, ens_size)

!!! Old code: has a check for 0.0 surface pressure
!!! https://github.com/NCAR/DART/blob/9729d784226295a197ca3bf00c917e4aaab5003b/models/wrf/model_mod.f90#L4600-L4606

else !extrapolate v_p(0:,) from pressure level 2 and v_p(1:,:)

  v_p(0,:) = extrap_4pressure(lev2_pres1(:), lev2_pres2(:), lev2_pres3(:), lev2_pres4(:), dx, dxm, dy, dym, ens_size, &
                              edgep=v_p(1,:))

endif

end subroutine get_model_pressure_profile

!------------------------------------------------------------------
! returns pressure at a point on the mass grid
function model_pressure_t(i,j,k,id,state_handle, ens_size)

integer,             intent(in) :: ens_size
integer,             intent(in) :: i,j,k(ens_size),id
type(ensemble_type), intent(in) :: state_handle
real(r8) :: model_pressure_t(ens_size)

integer(i8), dimension(ens_size) :: iqv, it
real(r8),    dimension(ens_size) :: qvf1, rho, x_iqv, x_it

integer :: var_idv, var_idt, e

! Adapted the code from WRF module_big_step_utilities_em.F ----
!         subroutine calc_p_rho_phi      Y.-R. Guo (10/20/2004)

! Simplification: alb*mub = (phb(i,j,k+1) - phb(i,j,k))/dnw(k)

var_idv = get_varid_from_kind(wrf_dom(id), QTY_VAPOR_MIXING_RATIO)
var_idt = get_varid_from_kind(wrf_dom(id), QTY_POTENTIAL_TEMPERATURE) ! HK original code type_t is this always QTY_POTENTIAL_TEMPERATURE
do e = 1, ens_size
  iqv = get_dart_vector_index(i,j,k(e), wrf_dom(id), var_idv)
  it  = get_dart_vector_index(i,j,k(e), wrf_dom(id), var_idt)
enddo

call get_state_array(x_iqv, iqv, state_handle)
call get_state_array(x_it, it, state_handle)

qvf1(:) = 1.0_r8 + x_iqv(:) / rd_over_rv

rho(:) = model_rho_t(i,j,k,id,state_handle, ens_size)

model_pressure_t(:) = ps0 * ( (gas_constant*(ts0+x_it)*qvf1) / &
     (ps0/rho(:)) )**cpovcv

end function model_pressure_t

!------------------------------------------------------------------
! returns surface pressure at a point on the mass grid
function model_pressure_s(i, j, id, state_handle, ens_size)

integer,             intent(in) :: i,j,id
type(ensemble_type), intent(in) :: state_handle
integer,             intent(in) :: ens_size

real(r8)              :: model_pressure_s(ens_size)

integer(i8) :: ips, imu
integer     :: var_id_psfc, var_id_mu
real(r8)    :: x_imu(1), x_ips(1)

var_id_psfc = get_varid_from_varname(wrf_dom(id), 'PSFC')
var_id_mu = get_varid_from_varname(wrf_dom(id), 'MU')

if ( var_id_PSFC > 0 ) then
   ips = get_dart_vector_index(i,j,1, wrf_dom(id), var_id_psfc)
   x_ips = get_state(ips, state_handle)
   model_pressure_s = x_ips

elseif (var_id_mu > 0) then
   imu = get_dart_vector_index(i,j,1, wrf_dom(id), var_id_mu)
   x_imu = get_state(imu, state_handle)
   model_pressure_s = stat_dat(id)%p_top + stat_dat(id)%mub(i,j) + x_imu

else
   call error_handler(E_ERR, 'model_pressure_s:', &
    'One of MU (QTY_PRESSURE) or PSFC (QTY_SURFACE_PRESSURE) must be in state vector to compute surface pressure', &
     source)
endif

end function model_pressure_s


!------------------------------------------------------------------
! Calculate the model level "zk" on half (mass) levels,
! corresponding to pressure "pres"
subroutine pres_to_zk(pres, mdl_v, n3, zk, level_below, lev0, fail)

real(r8), intent(in)  :: pres
real(r8), intent(in)  :: mdl_v(0:n3)
integer,  intent(in)  :: n3
real(r8), intent(out) :: zk
integer,  intent(out) :: level_below
logical,  intent(out) :: lev0
logical,  intent(out) :: fail

integer  :: k

lev0 = .false.
zk = MISSING_R8
fail = .false.

! if out of range completely, return missing_r8 and lev0 false
if (pres > mdl_v(0) .or. pres < mdl_v(n3)) then
   fail = .true.
   return
endif

! if above surface but below lowest sigma level, return the
! sigma value but set lev0 true
if(pres <= mdl_v(0) .and. pres > mdl_v(1)) then
  lev0 = .true.
  level_below = 1
  if (log_vert_interp) then
     zk = (log(mdl_v(0)) - log(pres))/(log(mdl_v(0)) - log(mdl_v(1)))
  else
     zk = (mdl_v(0) - pres)/(mdl_v(0) - mdl_v(1))
  endif
  return
 endif

! find the 2 sigma levels the value is between and return that
! as a real number, including the fraction between the levels.
do k = 1,n3-1
   if(pres <= mdl_v(k) .and. pres >= mdl_v(k+1)) then
      level_below = k
      if (log_vert_interp) then
         zk = real(k) + (log(mdl_v(k)) - log(pres))/(log(mdl_v(k)) - log(mdl_v(k+1)))
      else
         zk = real(k) + (mdl_v(k) - pres)/(mdl_v(k) - mdl_v(k+1))
      endif
      exit
   endif
enddo

end subroutine pres_to_zk

!------------------------------------------------------------------
! Calculate the model height profile on half (mass) levels,
! horizontally interpolated at the observation location.
subroutine get_model_height_profile(ll, ul, lr, ur, dx, dy, dxm, dym, id, v_h, state_handle, ens_size)

integer, dimension(2), intent(in)  :: ll, ul, lr, ur ! (x,y) mass grid corners
integer,  intent(in)  :: id
real(r8), intent(in)  :: dx,dy,dxm,dym
integer,  intent(in)  :: ens_size
real(r8), intent(out) :: v_h(0:grid(id)%bt, ens_size)
type(ensemble_type), intent(in)  :: state_handle
integer e !< for ensemble loop

real(r8)              :: fll(grid(id)%bts, ens_size), geop(ens_size), lat(ens_size)
integer(i8)           :: ill, iul, ilr, iur
integer               :: k, rc

real(r8), dimension(ens_size) :: x_ill, x_ilr, x_iul, x_iur
integer :: var_id

var_id = get_varid_from_kind(wrf_dom(id), QTY_GEOPOTENTIAL_HEIGHT)

do k = 1, grid(id)%bts  ! geopotential height (PH) is on bottom_top_stag

   ill = get_dart_vector_index(ll(1), ll(2), k, wrf_dom(id), var_id)
   iul = get_dart_vector_index(ul(1), ul(2), k, wrf_dom(id), var_id)
   ilr = get_dart_vector_index(lr(1), lr(2), k, wrf_dom(id), var_id)
   iur = get_dart_vector_index(ur(1), ur(2), k, wrf_dom(id), var_id)

   x_ill = get_state(ill, state_handle)
   x_ilr = get_state(ilr, state_handle)
   x_iul = get_state(iul, state_handle)
   x_iur = get_state(iur, state_handle)

   geop(:) = ( dym*( dxm*( stat_dat(id)%phb(ll(1),ll(2),k) + x_ill ) + &
                   dx*( stat_dat(id)%phb(lr(1),lr(2),k) + x_ilr ) ) + &
             dy*( dxm*( stat_dat(id)%phb(ul(1),ul(2),k) + x_iul ) + &
                   dx*( stat_dat(id)%phb(ur(1),ur(2),k) + x_iur ) ) )/gravity

   lat(:) = ( grid(id)%latitude(ll(1),ll(2)) + &
              grid(id)%latitude(lr(1),lr(2)) + &
              grid(id)%latitude(ul(1),ul(2)) + &
              grid(id)%latitude(ur(1),ur(2)) ) / 4.0_r8

   do e = 1, ens_size
      fll(k, e) = compute_geometric_height(geop(e), lat(e))
   enddo
end do

do k = 1, grid(id)%bt
   v_h(k, :) = 0.5_r8*(fll(k, :) + fll(k+1, :) )
end do

v_h(0, :) = dym*( dxm*stat_dat(id)%hgt(ll(1), ll(2)) + &
                dx*stat_dat(id)%hgt(lr(1), lr(2)) ) + &
          dy*( dxm*stat_dat(id)%hgt(ul(1), ul(2)) + &
                dx*stat_dat(id)%hgt(ur(1), ur(2)) )

 

end subroutine get_model_height_profile

!------------------------------------------------------------------
! Calculate the model level zk on half (mass) levels,
! corresponding to height obs_v.
subroutine height_to_zk(obs_v, mdl_v, n3, zk, level_below, lev0, fail)

real(r8), intent(in)  :: obs_v
integer,  intent(in)  :: n3
real(r8), intent(in)  :: mdl_v(0:n3)
real(r8), intent(out) :: zk
integer,  intent(out) :: level_below
logical,  intent(out) :: lev0
logical,  intent(out) :: fail

integer   :: k

zk = MISSING_R8
lev0 = .false.
fail = .false.

! if out of range completely, return missing_r8 and lev0 false
if (obs_v < mdl_v(0) .or. obs_v > mdl_v(n3)) then
   fail = .true.
   return
endif

! if above surface but below lowest 3-d height level, return the
! height value but set lev0 true
if(obs_v >= mdl_v(0) .and. obs_v < mdl_v(1)) then
  lev0 = .true.
  level_below = 1
  zk = (mdl_v(0) - obs_v)/(mdl_v(0) - mdl_v(1))
  return
endif

! find the 2 height levels the value is between and return that
! as a real number, including the fraction between the levels.
do k = 1,n3-1
   if(obs_v >= mdl_v(k) .and. obs_v <= mdl_v(k+1)) then
      level_below = k
      zk = real(k) + (mdl_v(k) - obs_v)/(mdl_v(k) - mdl_v(k+1))
      exit
   endif
enddo

end subroutine height_to_zk

!------------------------------------------------------------------
! Interpolate pressure inside quad
! Four corners of a quad:
!   p1 lower left
!   p2 lower right
!   p3 upper left
!   p4 upper right
! dx is the distance in x, dxm is 1.0-dx
! dy is distance in y, dym is 1.0-dy
function interp_4pressure(p1, p2, p3, p4, dx, dxm, dy, dym, ens_size)

integer, intent(in)                :: ens_size
real(r8), intent(in)               :: p1(ens_size), p2(ens_size), p3(ens_size), p4(ens_size)
real(r8), intent(in)               :: dx, dxm, dy, dym
real(r8)                           :: interp_4pressure(ens_size)

real(r8) :: l1(ens_size), l2(ens_size), l3(ens_size), l4(ens_size)

if (log_horz_interpQ) then
   l1 = log(p1)
   l2 = log(p2)
   l3 = log(p3)
   l4 = log(p4)
endif

! once we like the results, remove the log_horz_interpQ test.
if (log_horz_interpQ) then
   interp_4pressure = exp(dym*( dxm*l1 + dx*l2 ) + dy*( dxm*l3 + dx*l4 ))
else
   interp_4pressure = dym*( dxm*p1 + dx*p2 ) + dy*( dxm*p3 + dx*p4 )
endif

end function interp_4pressure

!------------------------------------------------------------------
! extrapolate quad where edgep is the edge pressure
function extrap_4pressure(p1, p2, p3, p4, dx, dxm, dy, dym, ens_size, edgep)

integer, intent(in)                :: ens_size
real(r8), intent(in)               :: p1(ens_size), p2(ens_size), p3(ens_size), p4(ens_size)
real(r8), intent(in)               :: dx, dxm, dy, dym
real(r8), intent(in)               :: edgep(ens_size)
real(r8)                           :: extrap_4pressure(ens_size)

real(r8) :: intermediate(ens_size)
real(r8) :: l1(ens_size), l2(ens_size), l3(ens_size), l4(ens_size)

if (log_horz_interpQ) then
   l1 = log(p1)
   l2 = log(p2)
   l3 = log(p3)
   l4 = log(p4)
endif

! once we like the results, remove the log_horz_interpQ test.
if (log_horz_interpQ) then
   intermediate = (3.0_r8*log(edgep) - &
              dym*( dxm*l1 + dx*l2 ) - dy*( dxm*l3 + dx*l4 ))/2.0_r8
   
   where (intermediate <= 0.0_r8)
      extrap_4pressure = edgep
   else where
       extrap_4pressure = exp(intermediate)
   end where
else
      extrap_4pressure = (3.0_r8*edgep - &
              dym*( dxm*p1 + dx*p2 ) - dy*( dxm*p3 + dx*p4 ))/2.0_r8
endif

end function extrap_4pressure

!------------------------------------------------------------------
! Calculate the total density on mass point (half (mass) levels, T-point).
function model_rho_t(i,j,k,id,state_handle, ens_size)

integer,             intent(in)  :: ens_size
integer,             intent(in)  :: i,j,k(ens_size),id
type(ensemble_type), intent(in)  :: state_handle
real(r8) :: model_rho_t(ens_size)

integer(i8), dimension(ens_size) :: imu,iph,iphp1
real(r8),    dimension(ens_size) :: ph_e, x_imu, x_iph, x_iphp1
integer :: var_id_mu, var_id_ph, e

! Adapted the code from WRF module_big_step_utilities_em.F ----
!         subroutine calc_p_rho_phi      Y.-R. Guo (10/20/2004)

! Simplification: alb*mub = (phb(i,j,k+1) - phb(i,j,k))/dnw(k)

var_id_mu = get_varid_from_varname(wrf_dom(id), 'MU')
var_id_ph = get_varid_from_kind(wrf_dom(id), QTY_GEOPOTENTIAL_HEIGHT)

do e = 1, ens_size
   imu   = get_dart_vector_index(i,j,1,      wrf_dom(id), var_id_mu)
   iph   = get_dart_vector_index(i,j,k(e),   wrf_dom(id), var_id_ph)
   iphp1 = get_dart_vector_index(i,j,k(e)+1, wrf_dom(id), var_id_ph)
enddo

call get_state_array(x_imu, imu, state_handle)
call get_state_array(x_iph, iph, state_handle)
call get_state_array(x_iphp1, iphp1, state_handle)

do e = 1, ens_size
   ph_e(e) = ( (x_iphp1(e) + stat_dat(id)%phb(i,j,k(e)+1)) &
             - (x_iph(e)   + stat_dat(id)%phb(i,j,k(e)  )) ) / stat_dat(id)%dnw(k(e))
enddo

! rho = - mu / dphi/deta
model_rho_t(:) = - (stat_dat(id)%mub(i,j)+x_imu) / ph_e

end function model_rho_t


!------------------------------------------------------------------
function density_interpolate(ens_size, state_handle, id, ll, ul, lr, ur, k, dxm, dx, dy, dym)

integer,             intent(in) :: ens_size
type(ensemble_type), intent(in) :: state_handle
integer,             intent(in) :: id
integer,             intent(in) :: ll(2), ul(2), lr(2), ur(2) ! (x,y) at  four corners
integer,             intent(in) :: k(ens_size) ! k may be different across the ensemble
real(r8),            intent(in) :: dxm, dx, dy, dym
real(r8) :: density_interpolate(ens_size)

real(r8), dimension(ens_size) :: rho1, rho2, rho3, rho4

rho1 = model_rho_t(ll(1), ll(2), k, id, state_handle, ens_size)
rho2 = model_rho_t(lr(1), lr(2), k, id, state_handle, ens_size)
rho3 = model_rho_t(ul(1), ul(2), k, id, state_handle, ens_size)
rho4 = model_rho_t(ur(1), ur(2), k, id, state_handle, ens_size)

density_interpolate = dym*( dxm*rho1(:) + dx*rho2(:) ) + dy*( dxm*rho3(:) + dx*rho4(:) )

end function density_interpolate

!------------------------------------------------------------------
! wrfinput land mask XLAND  1 = land, 2 = water
! obs_def_rttov_mod 0 = land, 1 = water, 2 = sea ice
function surface_type_interpolate(ens_size, id, ll, ul, lr, ur, dxm, dx, dy, dym)

integer,             intent(in) :: ens_size
integer,             intent(in) :: id
integer,             intent(in) :: ll(2), ul(2), lr(2), ur(2) ! (x,y) at  four corners
real(r8),            intent(in) :: dxm, dx, dy, dym
real(r8) :: surface_type_interpolate(ens_size) ! same across the ensemble

surface_type_interpolate(:) = -1 + dym*( dxm*stat_dat(id)%land(ll(1), ll(2))      + &
                                         dx*stat_dat(id)%land(lr(1), lr(2)) )     + &
                                         dy*( dxm*stat_dat(id)%land(ul(1), ul(2)) + &
                                         dx*stat_dat(id)%land(ur(1), ur(2)) )

end function surface_type_interpolate

!------------------------------------------------------------------
function surface_elevation_interpolate(ens_size, id, ll, ul, lr, ur, dxm, dx, dy, dym)

integer,             intent(in) :: ens_size
integer,             intent(in) :: id
integer,             intent(in) :: ll(2), ul(2), lr(2), ur(2) ! (x,y) at  four corners
real(r8),            intent(in) :: dxm, dx, dy, dym
real(r8) :: surface_elevation_interpolate(ens_size)

surface_elevation_interpolate(:) = dym*( dxm*stat_dat(id)%hgt(ll(1), ll(2))      + &
                                         dx*stat_dat(id)%hgt(lr(1), lr(2)) )     + &
                                         dy*( dxm*stat_dat(id)%hgt(ul(1), ul(2)) + &
                                         dx*stat_dat(id)%hgt(ur(1), ur(2)) )

end function surface_elevation_interpolate

!------------------------------------------------------------------
function geopotential_height_interpolate(ens_size, state_handle, qty, id, ll, ul, lr, ur, k, dxm, dx, dy, dym)

integer,             intent(in) :: ens_size
type(ensemble_type), intent(in) :: state_handle
integer,             intent(in) :: qty
integer,             intent(in) :: id
integer,             intent(in) :: ll(2), ul(2), lr(2), ur(2) ! (x,y) at  four corners
integer,             intent(in) :: k(ens_size) ! k may be different across the ensemble
real(r8),            intent(in) :: dxm, dx, dy, dym
real(r8) :: geopotential_height_interpolate(ens_size)

real(r8), dimension(ens_size) :: a1

a1 = simple_interpolation(ens_size, state_handle, qty, id, ll, ul, lr, ur, k, dxm, dx, dy, dym)

! phb is constant across the ensemble, so use k(1)
geopotential_height_interpolate = ( a1 + &
                                    dym* ( dxm*stat_dat(id)%phb(ll(1), ll(2), k(1) ) + &
                                    dx * stat_dat(id)%phb(lr(1), lr(2), k(1)) )      + &
                                    dy * ( dxm*stat_dat(id)%phb(ul(1), ul(2), k(1) ) + &
                                    dx * stat_dat(id)%phb(ur(1), ur(2), k(1)) )          )  / gravity

end function geopotential_height_interpolate

!------------------------------------------------------------------
function temperature_interpolate(ens_size, state_handle, id, ll, ul, lr, ur, k, dxm, dx, dy, dym)

integer,             intent(in) :: ens_size
type(ensemble_type), intent(in) :: state_handle
integer,             intent(in) :: id
integer,             intent(in) :: ll(2), ul(2), lr(2), ur(2) ! (x,y) at  four corners
integer,             intent(in) :: k(ens_size) ! k may be different across the ensemble
real(r8),            intent(in) :: dxm, dx, dy, dym
real(r8) :: temperature_interpolate(ens_size)

real(r8), dimension(ens_size) :: a1, pres, pres1, pres2, pres3, pres4

! In terms of perturbation potential temperature
a1 = simple_interpolation(ens_size, state_handle, QTY_POTENTIAL_TEMPERATURE, id, ll, ul, lr, ur, k, dxm, dx, dy, dym)

pres1 = model_pressure_t(ll(1), ll(2), k, id, state_handle, ens_size)
pres2 = model_pressure_t(lr(1), lr(2), k, id, state_handle, ens_size)
pres3 = model_pressure_t(ul(1), ul(2), k, id, state_handle, ens_size)
pres4 = model_pressure_t(ur(1), ur(2), k, id, state_handle, ens_size)

! Pressure at location
pres = dym*( dxm*pres1 + dx*pres2 ) + dy*( dxm*pres3 + dx*pres4 )

! Full sensible temperature field
temperature_interpolate = (ts0 + a1(:))*(pres(:)/ps0)**kappa

end function temperature_interpolate

!------------------------------------------------------------------
function pressure_interpolate(ens_size, state_handle, id, ll, ul, lr, ur, k, dxm, dx, dy, dym)

integer,             intent(in) :: ens_size
type(ensemble_type), intent(in) :: state_handle
integer,             intent(in) :: id
integer,             intent(in) :: ll(2), ul(2), lr(2), ur(2) ! (x,y) at  four corners
integer,             intent(in) :: k(ens_size) ! k may be different across the ensemble
real(r8),            intent(in) :: dxm, dx, dy, dym
real(r8) :: pressure_interpolate(ens_size)

real(r8), dimension(ens_size) :: pres1, pres2, pres3, pres4

pres1 = model_pressure_t(ll(1), ll(2), k, id, state_handle, ens_size)
pres2 = model_pressure_t(lr(1), lr(2), k, id, state_handle, ens_size)
pres3 = model_pressure_t(ul(1), ul(2), k, id, state_handle, ens_size)
pres4 = model_pressure_t(ur(1), ur(2), k, id, state_handle, ens_size)

! Pressure at location
pressure_interpolate = dym*( dxm*pres1 + dx*pres2 ) + dy*( dxm*pres3 + dx*pres4 )

end function pressure_interpolate

!------------------------------------------------------------------
function surface_pressure_interpolate(ens_size, state_handle, id, ll, ul, lr, ur, dxm, dx, dy, dym)

integer,             intent(in) :: ens_size
type(ensemble_type), intent(in) :: state_handle
integer,             intent(in) :: id
integer,             intent(in) :: ll(2), ul(2), lr(2), ur(2) ! (x,y) at  four corners
real(r8),            intent(in) :: dxm, dx, dy, dym
real(r8) :: surface_pressure_interpolate(ens_size)

real(r8), dimension(ens_size) :: pres1, pres2, pres3, pres4
integer :: e

pres1 = model_pressure_s(ll(1), ll(2), id, state_handle, ens_size)
pres2 = model_pressure_s(lr(1), lr(2), id, state_handle, ens_size)
pres3 = model_pressure_s(ul(1), ul(2), id, state_handle, ens_size)
pres4 = model_pressure_s(ur(1), ur(2), id, state_handle, ens_size)

! Pressure at location

do e = 1, ens_size
   ! HK todo original code comment:
   ! I'm not quite sure where this comes from, but I will trust them on it....
   if ( pres1 (e) /= 0.0_r8 .and. pres2(e) /= 0.0_r8 .and. pres3(e) /= 0.0_r8 .and. &
        pres4(e) /= 0.0_r8 ) then

    surface_pressure_interpolate(e) = dym*( dxm*pres1(e) + dx*pres2(e) ) + dy*( dxm*pres3(e) + dx*pres4(e) )
   endif  !HK todo initialize to missing_r8?
enddo


end function surface_pressure_interpolate

!------------------------------------------------------------------
function specific_humidity_interpolate(ens_size, state_handle, qty, id, ll, ul, lr, ur, k, dxm, dx, dy, dym)

integer,             intent(in) :: ens_size
type(ensemble_type), intent(in) :: state_handle
integer,             intent(in) :: qty
integer,             intent(in) :: id
integer,             intent(in) :: ll(2), ul(2), lr(2), ur(2) ! (x,y) at  four corners
integer,             intent(in) :: k(ens_size) ! k may be different across the ensemble
real(r8),            intent(in) :: dxm, dx, dy, dym
real(r8) :: specific_humidity_interpolate(ens_size)

real(r8), dimension(ens_size) :: a1

a1 = simple_interpolation(ens_size, state_handle, qty, id, ll, ul, lr, ur, k, dxm, dx, dy, dym)

specific_humidity_interpolate = a1(:) /(1.0_r8 + a1(:))

end function specific_humidity_interpolate

!------------------------------------------------------------------
function wind_interpolate(ens_size, state_handle, qty, id, k, xloc, yloc, i, j, dxm, dx, dy, dym, lon)

integer,             intent(in) :: ens_size
type(ensemble_type), intent(in) :: state_handle
integer,             intent(in) :: qty
integer,             intent(in) :: id
integer,             intent(in) :: k(ens_size) ! k may be different across the ensemble
real(r8),            intent(in) :: xloc, yloc ! location on mass grid
integer,             intent(in) :: i,j ! ll corners of mass grid
real(r8),            intent(in) :: dxm, dx, dy, dym
real(r8),            intent(in) :: lon ! Longitude of point in degrees
real(r8) :: wind_interpolate(ens_size)

real(r8), dimension(ens_size) :: u_wind_grid, v_wind_grid, u_wind, v_wind
real(r8) :: xloc_u, yloc_v  ! x ugrid, y vgrid
real(r8) :: dx_u, dxm_u, dy_v, dym_v
integer :: i_u, j_v
integer :: ll(2), ul(2), lr(2), ur(2) ! (x,y) at  four corners
integer :: e, rc

! HK TODO relationship between mass grid and u grid and v grid
! Original code adds 0.5 to xloc, yloc. But what if you are on the edge of a domain?
! https://github.com/NCAR/DART/blob/70e6af803a52d14b9f77f872c94b1fe11d5dc2d9/models/wrf/model_mod.f90#L1425-L1432

! xloc and yloc are indices on mass-grid.  If we are on a periodic longitude domain,
!   then xloc can range from [1 wes).  This means that simply adding 0.5 to xloc has
!   the potential to render xloc_u out of the valid mass-grid index bounds (>wes).
xloc_u = xloc + 0.5_r8
yloc_v = yloc + 0.5_r8

! HK TODO what about periodic_y?
if ( grid(id)%periodic_x .and. xloc_u > real(grid(id)%wes,r8) ) xloc_u = xloc_u - real(grid(id)%we,r8)

call toGrid(xloc_u,i_u,dx_u,dxm_u)
call toGrid(yloc_v,j_v,dy_v,dym_v)

call getCorners(i_u, j, id, QTY_U_WIND_COMPONENT, ll, ul, lr, ur, rc)
u_wind_grid = simple_interpolation(ens_size, state_handle, QTY_U_WIND_COMPONENT, id, ll, ul, lr, ur, k, dxm_u, dx_u, dy, dym)
call getCorners(i, j_v, id, QTY_V_WIND_COMPONENT, ll, ul, lr, ur, rc)
v_wind_grid = simple_interpolation(ens_size, state_handle, QTY_V_WIND_COMPONENT, id, ll, ul, lr, ur, k, dxm, dx, dy_v, dym_v)

do e = 1, ens_size
   call gridwind_to_truewind(lon, grid(id)%proj, u_wind_grid(e), v_wind_grid(e), u_wind(e), v_wind(e))
enddo

if ( qty == QTY_U_WIND_COMPONENT ) then
   wind_interpolate = u_wind
else
   wind_interpolate = v_wind
endif

end function wind_interpolate

!------------------------------------------------------------------
! If there are other domains in the state:
!      wrf domain id =/ state domain id
function get_wrf_domain(state_id)

integer, intent(in) :: state_id
integer :: get_wrf_domain

integer :: i

do i = 1, num_domains
   if (wrf_dom(i) == state_id) then
      get_wrf_domain = i
      return
   endif
enddo

end function get_wrf_domain

!------------------------------------------------------------------
subroutine convert_vertical_state(state_handle, num, locs, loc_qtys, loc_indx, &
                                  which_vert, istatus)

type(ensemble_type), intent(in)    :: state_handle
integer,             intent(in)    :: num
type(location_type), intent(inout) :: locs(num)     ! location for each state element
integer,             intent(in)    :: loc_qtys(num) ! qty for each state element
integer(i8),         intent(in)    :: loc_indx(num) ! index into the state vector
integer,             intent(in)    :: which_vert  ! vertical coordinate to be converted to
integer,             intent(out)   :: istatus

integer  :: ip, jp, kp, id ! x,y,z of state index
integer  :: var_id, state_id
real(r8) :: vert ! vertical after conversion
integer  :: i
integer  :: vert_coord_in
real(r8) :: lon_lat_vert(3)

istatus = 0 ! can not fail

do i = 1, num

   lon_lat_vert = get_location(locs(i))
   vert_coord_in = nint(query_location(locs(i)))

   if (vert_coord_in == which_vert) then ! conversion is already done
      cycle
   endif

   call get_model_variable_indices(loc_indx(i), ip, jp, kp, var_id=var_id, dom_id=state_id)
   id = get_wrf_domain(state_id)

   if (which_vert == VERTISLEVEL) then

      if (on_w_grid(state_id, var_id)) then
         vert = real(kp) - 0.5_r8
      else
          vert = real(kp)
      endif
   
   elseif (which_vert == VERTISPRESSURE) then
   
      vert = model_pressure(ip, jp, kp, id, var_id, state_id, state_handle)
   
   elseif (which_vert == VERTISHEIGHT) then

      vert = model_height(ip, jp, kp, id, loc_qtys(i), var_id, state_id, state_handle)

   elseif (which_vert == VERTISSCALEHEIGHT) then
   
      vert = -log(model_pressure(ip, jp, kp, id, var_id, state_id, state_handle) / &
               model_surface_pressure(ip, jp, id, var_id, state_id, state_handle))
   
   endif

   call set_vertical(locs(i), vert, which_vert)

enddo



end subroutine convert_vertical_state

!------------------------------------------------------------------
subroutine convert_vertical_obs(state_handle, num, locs, loc_qtys, loc_types, &
                                which_vert, istatus)

type(ensemble_type), intent(in)    :: state_handle
integer,             intent(in)    :: num
type(location_type), intent(inout) :: locs(num)
integer,             intent(in)    :: loc_qtys(num)
integer,             intent(in)    :: loc_types(num)
integer,             intent(in)    :: which_vert
integer,             intent(out)   :: istatus(num)

integer  :: vert_coord_in, id
real(r8) :: lon_lat_vert(3)
integer  :: ll(2), ul(2), lr(2), ur(2) !(x,y) of four corners
integer  :: i, j ! grid
real(r8) :: xloc, yloc ! WRF i,j in the grid
real(r8) :: dx, dxm, dy, dym ! grid fractions
real(r8) :: zout(1), zk(1), zk1(1), geop(1), zloc(1)
real(r8) :: pres1(1), pres2(1), pres3(1), pres4(1)
integer  :: ens_size, rc, k(1)
logical  :: fail
integer  :: ob ! loop variable
real(r8) :: dz, dzm ! dummys for toGrid call

integer, parameter :: FAILED_BOUNDS_CHECK = 144
integer, parameter :: CANNOT_INTERPOLATE_QTY = 155
integer, parameter :: NOT_IN_ANY_DOMAIN = 111
integer, parameter :: VERTICAL_LOCATION_FAIL = 166

ens_size = 1 ! working with the mean state

do ob = 1, num

   lon_lat_vert = get_location(locs(ob))
   vert_coord_in = nint(query_location(locs(ob)))
   
   if (vert_coord_in == which_vert) then ! conversion is already done
      istatus(ob) = 0
      cycle
   endif
   
   if (lon_lat_vert(3) == VERTISUNDEF) then ! no vertical, no conversion
      istatus(ob) = 0
      cycle
   endif
   
   if (lon_lat_vert(3) == MISSING_R8) then ! vertical is missing, no conversion
      call set_vertical(locs(ob), MISSING_R8, which_vert)
      istatus(ob) = 1  ! HK todo original code does not set success for this
      cycle
   endif
   
   ! convert to which_vert
   call get_domain_info(lon_lat_vert(1),lon_lat_vert(2),id,xloc,yloc)
   if (id == 0) then
      call set_vertical(locs(ob), MISSING_R8, which_vert) !HK original code - why set to missing?
      istatus(ob) = NOT_IN_ANY_DOMAIN
      cycle
   endif
   
   ! horizontal location mass point
   call toGrid(xloc,i,dx,dxm)
   call toGrid(yloc,j,dy,dym)
   
   if ( .not. within_bounds_horizontal(i, j, id, loc_qtys(ob)) ) then
      call set_vertical(locs(ob), MISSING_R8, which_vert)
      istatus(ob) = FAILED_BOUNDS_CHECK
      cycle
   endif
   
   call getCorners(i, j, id, QTY_POTENTIAL_TEMPERATURE, ll, ul, lr, ur, rc) !HK todo what qty to pick?
   
   ! vertical location mass level
   call get_level_below_obs(vert_coord_in, id, lon_lat_vert, ens_size, state_handle, ll, ul, lr, ur, dx, dy, dxm, dym, k, zloc, fail)
   if (fail) then
      istatus(ob) = VERTICAL_LOCATION_FAIL
      ! set vertical?
      cycle
   endif
   
   select case (which_vert)
   
      case (VERTISLEVEL)

           if (vert_coord_in == VERTISSURFACE) then
              zout = 1.0_r8
           else
              zout = zloc
           endif

      case (VERTISPRESSURE)

        if (vert_coord_in == VERTISSURFACE) then
           ! compute surface pressure at all neighboring mass points
           pres1 = model_pressure_s(ll(1), ll(2), id, state_handle, ens_size)
           pres2 = model_pressure_s(lr(1), lr(2), id, state_handle, ens_size)
           pres3 = model_pressure_s(ul(1), ul(2), id, state_handle, ens_size)
           pres4 = model_pressure_s(ur(1), ur(2), id, state_handle, ens_size)
           zout = dym*( dxm*pres1(1) + dx*pres2(1) ) + dy*( dxm*pres3(1) + dx*pres4(1) )

        else
           zk  = pressure_interpolate(ens_size, state_handle, id, ll, ul, lr, ur, k, dxm, dx, dy, dym)
           zk1 = pressure_interpolate(ens_size, state_handle, id, ll, ul, lr, ur, k+1, dxm, dx, dy, dym)
           zout = vertical_interpolation(ens_size, zloc, zk, zk1)
        endif

      case (VERTISHEIGHT)

         if (vert_coord_in == VERTISSURFACE) then
            ! a surface ob is assumed to have height as vertical coordinate.
            ! this code needs to be revised if this is not true
            ! (in that case uncomment lines below to get terrain height
            ! from model)
            zout = lon_lat_vert(3)
            !! or: directly interpolate terrain height at neighboring mass points
            !zout = dym*( dxm*stat_dat(id)%hgt(i,  j) + &
            !              dx*stat_dat(id)%hgt(i+1,j) ) + &
            !        dy*( dxm*stat_dat(id)%hgt(i,  j+1) + &
            !              dx*stat_dat(id)%hgt(i+1,j+1) )
         else
            ! adding 0.5 to get to the staggered vertical grid for height
            zloc = zloc + 0.5_r8 ! Adjust zloc for staggered
            call toGrid(zloc(1),k(1),dz,dzm)

            ! This method does not give bitwise with main
            !zk  = geopotential_height_interpolate(ens_size, state_handle, QTY_GEOPOTENTIAL_HEIGHT, id, ll, ul, lr, ur, k, dxm, dx, dy, dym)
            !zk1 = geopotential_height_interpolate(ens_size, state_handle, QTY_GEOPOTENTIAL_HEIGHT, id, ll, ul, lr, ur, k+1, dxm, dx, dy, dym)
            !geop = vertical_interpolation(ens_size, zloc, zk, zk1)
            !zout = compute_geometric_height(geop(1), lon_lat_vert(2))

            zk = interpolate_geometric_height(ens_size, state_handle, id, ll, ul, lr, ur, k, dxm, dx, dy, dym)
            zk1 = interpolate_geometric_height(ens_size, state_handle, id, ll, ul, lr, ur, k+1, dxm, dx, dy, dym)
            zout = vertical_interpolation(ens_size, zloc, zk, zk1)
        endif

      case (VERTISSCALEHEIGHT)
         if (vert_coord_in == VERTISSURFACE) then
            zout = -log(1.0_r8)
         elseif (vert_coord_in == VERTISPRESSURE) then
            ! surface pressure
            pres1 = model_pressure_s(ll(1), ll(2), id, state_handle, ens_size)
            pres2 = model_pressure_s(lr(1), lr(2), id, state_handle, ens_size)
            pres3 = model_pressure_s(ul(1), ul(2), id, state_handle, ens_size)
            pres4 = model_pressure_s(ur(1), ur(2), id, state_handle, ens_size)

            zout = -log(lon_lat_vert(3) / (dym*( dxm*pres1(1) + dx*pres2(1) ) + dy*( dxm*pres3(1) + dx*pres4(1) )))

         else ! vert_cood_in == VERTISHEIGHT

            zk  = pressure_interpolate(ens_size, state_handle, id, ll, ul, lr, ur, k, dxm, dx, dy, dym)
            zk1 = pressure_interpolate(ens_size, state_handle, id, ll, ul, lr, ur, k+1, dxm, dx, dy, dym)
            zout = vertical_interpolation(ens_size, zloc, zk, zk1)

            ! surface pressure
            pres1 = model_pressure_s(ll(1), ll(2), id, state_handle, ens_size)
            pres2 = model_pressure_s(lr(1), lr(2), id, state_handle, ens_size)
            pres3 = model_pressure_s(ul(1), ul(2), id, state_handle, ens_size)
            pres4 = model_pressure_s(ur(1), ur(2), id, state_handle, ens_size)

            zout = -log(zout / (dym*( dxm*pres1(1) + dx*pres2(1) ) + dy*( dxm*pres3(1) + dx*pres4(1) )))

         endif
        

   end select
 
   !HK original code uses set_location, you could use set_vertical here see #621
   locs(ob) = set_location(lon_lat_vert(1), lon_lat_vert(2), zout(1), which_vert)
   istatus(ob) = 0

enddo

end subroutine convert_vertical_obs


!------------------------------------------------------------------
function interpolate_geometric_height(ens_size, state_handle, id, ll, ul, lr, ur, k, dxm, dx, dy, dym)

integer,             intent(in) :: ens_size
type(ensemble_type), intent(in) :: state_handle
integer,             intent(in) :: id
integer,             intent(in) :: ll(2), ul(2), lr(2), ur(2) ! (x,y) at  four corners
integer,             intent(in) :: k(ens_size) ! k may be different across the ensemble
real(r8),            intent(in) :: dxm, dx, dy, dym
real(r8) :: interpolate_geometric_height(ens_size)

integer :: e ! loop variable
! lower left, upper left, lower right, upper right
integer(i8), dimension(ens_size)  :: ill, iul, ilr, iur  ! dart index at four corners
real(r8),    dimension(ens_size)  :: x_ill, x_iul, x_ilr, x_iur ! state value at four corners
real(r8),    dimension(ens_size)  :: geop_ll, geop_ul, geop_lr, geop_ur ! geopotential height at four corners
real(r8),    dimension(ens_size)  :: geomet_ll, geomet_ul, geomet_lr, geomet_ur ! geometric height at four corners
integer :: var_id

var_id = get_varid_from_kind(wrf_dom(id), QTY_GEOPOTENTIAL_HEIGHT)

do e = 1, ens_size
                              !   x,     y,     z,    domain,      variable
   ill(e) = get_dart_vector_index(ll(1), ll(2), k(e), wrf_dom(id), var_id)
   iul(e) = get_dart_vector_index(ul(1), ul(2), k(e), wrf_dom(id), var_id)
   ilr(e) = get_dart_vector_index(lr(1), lr(2), k(e), wrf_dom(id), var_id)
   iur(e) = get_dart_vector_index(ur(1), ur(2), k(e), wrf_dom(id), var_id)

enddo

call get_state_array(x_ill, ill, state_handle)
call get_state_array(x_iul, iul, state_handle)
call get_state_array(x_ilr, ilr, state_handle)
call get_state_array(x_iur, iur, state_handle)

do e = 1, ens_size
   geop_ll(e) = (stat_dat(id)%phb(ll(1), ll(2), k(e)) + x_ill(e)) / gravity
   geop_ul(e) = (stat_dat(id)%phb(ul(1), ul(2), k(e)) + x_iul(e)) / gravity
   geop_lr(e) = (stat_dat(id)%phb(lr(1), lr(2), k(e)) + x_ilr(e)) / gravity
   geop_ur(e) = (stat_dat(id)%phb(ur(1), ur(2), k(e)) + x_iur(e)) / gravity

   geomet_ll(e) = compute_geometric_height(geop_ll(e), grid(id)%latitude(ll(1), ll(2)))
   geomet_ul(e) = compute_geometric_height(geop_ul(e), grid(id)%latitude(ul(1), ul(2)))
   geomet_lr(e) = compute_geometric_height(geop_lr(e), grid(id)%latitude(lr(1), lr(2)))
   geomet_ur(e) = compute_geometric_height(geop_ur(e), grid(id)%latitude(ur(1), ur(2)))

enddo

interpolate_geometric_height = dym*( dxm*geomet_ll(:) + dx*geomet_lr(:) ) + dy*( dxm*geomet_ul(:) + dx*geomet_ur(:) )

end function interpolate_geometric_height

!------------------------------------------------------------------
! model height any grid u,v,w,t
! used in convert_vertical_state so ens_size = 1
function model_height(i, j, k, id, qty, var_id, state_id, state_handle)

type(ensemble_type), intent(in) :: state_handle
integer,             intent(in) :: i, j, k, id, qty
integer,             intent(in) :: var_id,state_id
real(r8)                        :: model_height

integer(i8)  :: i1, i2, i3, i4
integer      :: off
real(r8)     :: x_i1(1), x_i2(1), x_i3(1), x_i4(1)
real(r8)     :: geop, lat
integer      :: gz_id

gz_id = get_varid_from_kind(state_id, QTY_GEOPOTENTIAL_HEIGHT)

!HK todo for these special cases would it be better to check by variable name
! instead of QTY?
if( qty == QTY_PRESSURE .or. &             ! MU
        qty == QTY_SURFACE_PRESSURE .or. & ! PSFC SFC PRESSUR
        qty == QTY_SKIN_TEMPERATURE) then  ! TSK SURFACE SKIN TEMPERATURE

   model_height = stat_dat(id)%hgt(i,j)

elseif( qty == QTY_SOIL_TEMPERATURE  .or. &  ! TSLB SOIL TEMPERATURE
        qty == QTY_SOIL_MOISTURE .or. &  ! SMOIS SOIL MOISTURE
        qty == QTY_SOIL_LIQUID_WATER ) then  !  SH2O SOIL LIQUID WATER

   model_height = stat_dat(id)%hgt(i,j) - stat_dat(id)%zs(k)

elseif( qty == QTY_10M_U_WIND_COMPONENT .or. &
        qty == QTY_10M_V_WIND_COMPONENT  ) then

   model_height = stat_dat(id)%hgt(i,j) + 10.0_r8

elseif( qty == QTY_2M_TEMPERATURE  .or. &
        !qty == QTY_POTENTIAL_TEMPERATURE .or. & ! TH2 POT TEMP at 2 M
        qty == QTY_2M_SPECIFIC_HUMIDITY ) then  ! Q2 QV at 2 M

   model_height = stat_dat(id)%hgt(i,j) + 2.0_r8

! If W-grid (on ZNW levels), native to GZ
elseif( on_w_grid(state_id, var_id) ) then

   i1 = get_dart_vector_index(i,j,k, state_id, gz_id)
   x_i1 = get_state(i1, state_handle)

   geop = minval(stat_dat(id)%phb(i,j,k)+x_i1)/gravity
   model_height = compute_geometric_height(geop, grid(id)%latitude(i, j))

! If U-grid, then height is defined between U points, both in horizontal
!   and in vertical, so average -- averaging depends on longitude periodicity
elseif( on_u_grid(state_id, var_id) ) then

   if( i == grid(id)%wes ) then

      ! Check to see if periodic in longitude
      if ( grid(id)%periodic_x ) then

         ! We are at the seam in longitude, so take first and last mass points
         i1 = get_dart_vector_index(i-1,j,k  , state_id, gz_id)
         i2 = get_dart_vector_index(i-1,j,k+1, state_id, gz_id)
         i3 = get_dart_vector_index(1,  j,k  , state_id, gz_id)
         i4 = get_dart_vector_index(1,  j,k+1, state_id, gz_id)

         x_i1 = get_state(i1, state_handle)
         x_i2 = get_state(i2, state_handle)
         x_i3 = get_state(i3, state_handle)
         x_i4 = get_state(i4, state_handle)

! HK todo what is minval for? Is it just for converting an array to a scalar?
         geop = minval(( (stat_dat(id)%phb(i-1,j,k  ) + x_i1) &
                 +(stat_dat(id)%phb(i-1,j,k+1) + x_i2) &
                 +(stat_dat(id)%phb(1  ,j,k  ) + x_i3) &
                 +(stat_dat(id)%phb(1  ,j,k+1) + x_i4) )/(4.0_r8*gravity))
         
         lat = ( grid(id)%latitude(i-1,j)  &
                +grid(id)%latitude(i-1,j)  &
                +grid(id)%latitude(1  ,j)  &
                +grid(id)%latitude(1  ,j) ) / 4.0_r8

         model_height = compute_geometric_height(geop, lat)
         
      else

         ! If not periodic, then try extrapolating
         i1 = get_dart_vector_index(i-1,j,k  , state_id, gz_id)
         i2 = get_dart_vector_index(i-1,j,k+1, state_id, gz_id)

         x_i1 = get_state(i1, state_handle)
         x_i2 = get_state(i2, state_handle)
         x_i3 = get_state(i1 -1, state_handle)
         x_i4 = get_state(i2 -1, state_handle)


         geop = minval(( 3.0_r8*(stat_dat(id)%phb(i-1,j,k  )+x_i1) &
                 +3.0_r8*(stat_dat(id)%phb(i-1,j,k+1)+x_i2) &
                        -(stat_dat(id)%phb(i-2,j,k  )+x_i3) &
                        -(stat_dat(id)%phb(i-2,j,k+1)+x_i4) )/(4.0_r8*gravity))

         lat = ( 3.0_r8*grid(id)%latitude(i-1,j)  &
                +3.0_r8*grid(id)%latitude(i-1,j)  &
                       -grid(id)%latitude(i-2,j)  &
                       -grid(id)%latitude(i-2,j)) / 4.0_r8

         model_height = compute_geometric_height(geop, lat)

      endif

   elseif( i == 1 ) then

      ! Check to see if periodic in longitude
      if ( grid(id)%periodic_x ) then

         ! We are at the seam in longitude, so take first and last mass points
         off = grid(id)%we
         i1 = get_dart_vector_index(i  ,j,k  ,state_id, gz_id)
         i2 = get_dart_vector_index(i  ,j,k+1,state_id, gz_id)
         i3 = get_dart_vector_index(off,j,k  ,state_id, gz_id)
         i4 = get_dart_vector_index(off,j,k+1,state_id, gz_id)

         x_i1 = get_state(i1, state_handle)
         x_i2 = get_state(i2, state_handle)
         x_i3 = get_state(i3, state_handle)
         x_i4 = get_state(i4, state_handle)

         geop = minval(( (stat_dat(id)%phb(i  ,j,k  ) + x_i1) &
                 +(stat_dat(id)%phb(i  ,j,k+1) + x_i2) &
                 +(stat_dat(id)%phb(off,j,k  ) + x_i3) &
                 +(stat_dat(id)%phb(off,j,k+1) + x_i4) )/(4.0_r8*gravity))
         
         lat = ( grid(id)%latitude(i  ,j)  &
                +grid(id)%latitude(i  ,j)  &
                +grid(id)%latitude(off,j)  &
                +grid(id)%latitude(off,j)) / 4.0_r8

         model_height = compute_geometric_height(geop, lat)

      else

         ! If not periodic, then try extrapolating
         i1 = get_dart_vector_index(i,j,k  ,state_id, gz_id)
         i2 = get_dart_vector_index(i,j,k+1,state_id, gz_id)

         x_i1 = get_state(i1, state_handle)
         x_i2 = get_state(i2, state_handle)
         x_i3 = get_state(i1 +1, state_handle)
         x_i4 = get_state(i2 +1, state_handle)


         geop = minval(( 3.0_r8*(stat_dat(id)%phb(i  ,j,k  )+x_i1)   &
                 +3.0_r8*(stat_dat(id)%phb(i  ,j,k+1)+x_i2)   &
                        -(stat_dat(id)%phb(i+1,j,k  )+x_i3) &
                        -(stat_dat(id)%phb(i+1,j,k+1)+x_i4) )/(4.0_r8*gravity))

         lat = ( 3.0_r8*grid(id)%latitude(i  ,j)  &
                +3.0_r8*grid(id)%latitude(i  ,j)  &
                       -grid(id)%latitude(i+1,j)  &
                       -grid(id)%latitude(i+1,j)) / 4.0_r8

         model_height = compute_geometric_height(geop, lat)

      endif

   else

      i1 = get_dart_vector_index(i,j,k  ,state_id, gz_id)
      i2 = get_dart_vector_index(i,j,k+1,state_id, gz_id)

      x_i1 = get_state(i1, state_handle)
      x_i2 = get_state(i2, state_handle)
      x_i3 = get_state(i1 -1, state_handle)
      x_i4 = get_state(i2 -1, state_handle)


      geop = minval(( (stat_dat(id)%phb(i  ,j,k  )+x_i1)   &
              +(stat_dat(id)%phb(i  ,j,k+1)+x_i2)   &
              +(stat_dat(id)%phb(i-1,j,k  )+x_i3) &
              +(stat_dat(id)%phb(i-1,j,k+1)+x_i4) )/(4.0_r8*gravity))

      lat = (  grid(id)%latitude(i  ,j)  &
              +grid(id)%latitude(i  ,j)  &
              +grid(id)%latitude(i-1,j)  &
              +grid(id)%latitude(i-1,j)) / 4.0_r8

      model_height = compute_geometric_height(geop, lat)

   endif

! If V-grid, then pressure is defined between V points, both in horizontal
!   and in vertical, so average -- averaging depends on polar periodicity
elseif( on_v_grid(state_id, var_id) ) then

   if( j == grid(id)%sns ) then

      ! Check to see if periodic in latitude (polar)
      if ( grid(id)%polar ) then

         ! The upper corner is 180 degrees of longitude away
         off = i + grid(id)%we/2
         if ( off > grid(id)%we ) off = off - grid(id)%we

         i1 = get_dart_vector_index(off,j-1,k  , state_id, gz_id)
         i2 = get_dart_vector_index(off,j-1,k+1, state_id, gz_id)
         i3 = get_dart_vector_index(i  ,j-1,k  , state_id, gz_id)
         i4 = get_dart_vector_index(i  ,j-1,k+1, state_id, gz_id)

         x_i1 = get_state(i1, state_handle)
         x_i2 = get_state(i2, state_handle)
         x_i3 = get_state(i3, state_handle)
         x_i4 = get_state(i4, state_handle)

         geop = minval(( (stat_dat(id)%phb(off,j-1,k  )+x_i1) &
                 +(stat_dat(id)%phb(off,j-1,k+1)+x_i2) &
                 +(stat_dat(id)%phb(i  ,j-1,k  )+x_i3) &
                 +(stat_dat(id)%phb(i  ,j-1,k+1)+x_i4) )/(4.0_r8*gravity))
         
         lat = ( grid(id)%latitude(off,j-1)  &
                +grid(id)%latitude(off,j-1)  &
                +grid(id)%latitude(i  ,j-1)  &
                +grid(id)%latitude(i  ,j-1)) / 4.0_r8

         model_height = compute_geometric_height(geop, lat)

      else

         ! If not periodic, then try extrapolating
         i1 = get_dart_vector_index(i,j-1,k ,  state_id, gz_id)
         i2 = get_dart_vector_index(i,j-1,k+1, state_id, gz_id)
         i3 = get_dart_vector_index(i,j-2,k  , state_id, gz_id)
         i4 = get_dart_vector_index(i,j-2,k+1, state_id, gz_id)

         x_i1 = get_state(i1, state_handle)
         x_i2 = get_state(i2, state_handle)
         x_i3 = get_state(i3, state_handle)
         x_i4 = get_state(i4, state_handle)

         geop = minval(( 3.0_r8*(stat_dat(id)%phb(i,j-1,k  )+x_i1) &
                 +3.0_r8*(stat_dat(id)%phb(i,j-1,k+1)+x_i2) &
                        -(stat_dat(id)%phb(i,j-2,k  )+x_i3) &
                        -(stat_dat(id)%phb(i,j-2,k+1)+x_i4) )/(4.0_r8*gravity))

         lat = ( 3.0_r8*grid(id)%latitude(i,j-1)  &
                +3.0_r8*grid(id)%latitude(i,j-1)  &
                       -grid(id)%latitude(i,j-2)  &
                       -grid(id)%latitude(i,j-2)) / 4.0_r8

         model_height = compute_geometric_height(geop, lat)

      endif

   elseif( j == 1 ) then

      ! Check to see if periodic in latitude (polar)
      if ( grid(id)%polar ) then

         ! The lower corner is 180 degrees of longitude away
         off = i + grid(id)%we/2
         if ( off > grid(id)%we ) off = off - grid(id)%we

         i1 = get_dart_vector_index(off,j,k  , state_id, gz_id)
         i2 = get_dart_vector_index(off,j,k+1, state_id, gz_id)
         i3 = get_dart_vector_index(i  ,j,k  , state_id, gz_id)
         i4 = get_dart_vector_index(i  ,j,k+1, state_id, gz_id)

         x_i1 = get_state(i1, state_handle)
         x_i2 = get_state(i2, state_handle)
         x_i3 = get_state(i3, state_handle)
         x_i4 = get_state(i4, state_handle)

         geop = minval(( (stat_dat(id)%phb(off,j,k  )+x_i1) &
                 +(stat_dat(id)%phb(off,j,k+1)+x_i2) &
                 +(stat_dat(id)%phb(i  ,j,k  )+x_i3) &
                 +(stat_dat(id)%phb(i  ,j,k+1)+x_i4) )/(4.0_r8*gravity))
         
         lat = ( grid(id)%latitude(off,j)  &
                +grid(id)%latitude(off,j)  &
                +grid(id)%latitude(i  ,j)  &
                +grid(id)%latitude(i  ,j)) / 4.0_r8

         model_height = compute_geometric_height(geop, lat)

      else

         ! If not periodic, then try extrapolating
         i1 = get_dart_vector_index(i,j  ,k  , state_id, gz_id)
         i2 = get_dart_vector_index(i,j  ,k+1, state_id, gz_id)
         i3 = get_dart_vector_index(i,j+1,k  , state_id, gz_id)
         i4 = get_dart_vector_index(i,j+1,k+1, state_id, gz_id)

         x_i1 = get_state(i1, state_handle)
         x_i2 = get_state(i2, state_handle)
         x_i3 = get_state(i3, state_handle)
         x_i4 = get_state(i4, state_handle)

         geop = minval(( 3.0_r8*(stat_dat(id)%phb(i,j  ,k  )+x_i1) &
                 +3.0_r8*(stat_dat(id)%phb(i,j  ,k+1)+x_i2) &
                        -(stat_dat(id)%phb(i,j+1,k  )+x_i3) &
                        -(stat_dat(id)%phb(i,j+1,k+1)+x_i4) )/(4.0_r8*gravity))

         lat = ( 3.0_r8*grid(id)%latitude(i,j  )  &
                +3.0_r8*grid(id)%latitude(i,j  )  &
                       -grid(id)%latitude(i,j+1)  &
                       -grid(id)%latitude(i,j+1)) / 4.0_r8

         model_height = compute_geometric_height(geop, lat)

      endif

   else

      i1 = get_dart_vector_index(i,j  ,k  , state_id, gz_id)
      i2 = get_dart_vector_index(i,j  ,k+1, state_id, gz_id)
      i3 = get_dart_vector_index(i,j-1,k  , state_id, gz_id)
      i4 = get_dart_vector_index(i,j-1,k+1, state_id, gz_id)

      x_i1 = get_state(i1, state_handle)
      x_i2 = get_state(i2, state_handle)
      x_i3 = get_state(i3, state_handle)
      x_i4 = get_state(i4, state_handle)

      geop = minval(( (stat_dat(id)%phb(i,j  ,k  )+x_i1) &
              +(stat_dat(id)%phb(i,j  ,k+1)+x_i2) &
              +(stat_dat(id)%phb(i,j-1,k  )+x_i3) &
              +(stat_dat(id)%phb(i,j-1,k+1)+x_i4) )/(4.0_r8*gravity))

      lat = ( grid(id)%latitude(i,j  )  &
             +grid(id)%latitude(i,j  )  &
             +grid(id)%latitude(i,j-1)  &
             +grid(id)%latitude(i,j-1)) / 4.0_r8

      model_height = compute_geometric_height(geop, lat)

   endif

else

   i1 = get_dart_vector_index(i,j,k  , state_id, gz_id)
   i2 = get_dart_vector_index(i,j,k+1, state_id, gz_id)

   x_i1 = get_state(i1, state_handle)
   x_i2 = get_state(i2, state_handle)

   geop = minval(( (stat_dat(id)%phb(i,j,k  )+x_i1) &
           +(stat_dat(id)%phb(i,j,k+1)+x_i2) )/(2.0_r8*gravity))

   lat = grid(id)%latitude(i,j)

   model_height = compute_geometric_height(geop, lat)

endif

end function model_height

!------------------------------------------------------------------
pure function interp_pressure(p1, p2, use_log)

real(r8), intent(in) :: p1(1), p2(1)
logical,  intent(in) :: use_log
real(r8) :: interp_pressure

if (use_log) then
   interp_pressure = exp((log(p1(1)) + log(p2(1)))/2.0_r8)
else
   interp_pressure = (p1(1) + p2(1))/2.0_r8
endif

end function interp_pressure

!------------------------------------------------------------------
pure function extrap_pressure(p1, p2, use_log)

real(r8), intent(in) :: p1(1), p2(1)
logical,  intent(in) :: use_log

real(r8) :: extrap_pressure
real(r8) :: intermediate

if (use_log) then
   intermediate = (3.0_r8*log(p1(1)) - log(p2(1)))/2.0_r8
   if (intermediate <= 0.0_r8) then
      extrap_pressure = p1(1)
   else
      extrap_pressure = exp(intermediate)
   endif
else
   extrap_pressure = (3.0_r8*p1(1) - p2(1))/2.0_r8
endif

end function extrap_pressure

!------------------------------------------------------------------
! model pressure any grid (w,u,v,t)
! used in convert_vertical_state so ens_size = 1
function model_pressure(i, j, kp, id, var_id, state_id, state_handle)

type(ensemble_type), intent(in) :: state_handle
integer,             intent(in) :: i,j,kp,id
integer,             intent(in) :: var_id
integer,             intent(in) :: state_id
real(r8)                        :: model_pressure

real(r8) :: p(1), p_one(1), p_two(1) ! only using the mean, so calling model_pressure_t with ens_size=1
integer   :: k(1), off, n

k(1) = kp ! array version

! HK this is for soil variables
do n = 1, get_num_dims(state_id, var_id)
   if ( get_dim_name(state_id, var_id, n) == 'soil_layers_stag' ) then
      p = model_pressure_s(i, j, id, state_handle, 1)
      model_pressure = p(1)
      return
   endif
enddo

if (get_num_dims(state_id, var_id) == 2) then ! surface qty
   p = model_pressure_s(i, j, id, state_handle, 1)
   model_pressure = p(1)
   return
endif

if (on_w_grid(state_id, var_id)) then ! average in the vertical

   if (kp==1) then ! on boundary, extrapolate in the vertical

      p_one = model_pressure_t(i, j, k,   id, state_handle, 1)
      p_two = model_pressure_t(i, j, k+1, id, state_handle, 1)
      model_pressure = extrap_pressure(p_one, p_two, log_vert_interp)

   elseif ( kp==grid(id)%bts ) then ! on boundary, extrapolate in the vertical

      p_one = model_pressure_t(i, j, k-1, id, state_handle, 1)
      p_two = model_pressure_t(i, j, k-2, id, state_handle, 1)
      model_pressure = extrap_pressure(p_one, p_two, log_vert_interp)

   else  ! interpolate 

      p_one = model_pressure_t(i, j, k,  id, state_handle, 1)
      p_two = model_pressure_t(i, j, k-1,id, state_handle, 1)
      model_pressure = interp_pressure(p_one, p_two, log_vert_interp)

   endif

elseif (on_u_grid(state_id, var_id)) then ! average in the horizontal u direction

   if  (i==grid(id)%wes) then

      if ( grid(id)%periodic_x ) then

          ! We are at seam in longitude, take first and last M-grid points
          p_one = model_pressure_t(i-1,j, k, id, state_handle, 1)
          p_two = model_pressure_t(1,  j, k, id, state_handle, 1)
          model_pressure = interp_pressure(p_one, p_two, log_horz_interpM)

      else

          ! If not periodic, then try extrapolating
          p_one = model_pressure_t(i-1, j, k, id, state_handle, 1)
          p_two = model_pressure_t(i-2, j, k, id, state_handle, 1)
          model_pressure = extrap_pressure(p_one, p_two, log_horz_interpM)

      endif

   elseif (i==1) then

      if ( grid(id)%periodic_x ) then

         ! We are at seam in longitude, take first and last M-grid points
         p_one = model_pressure_t(i,           j, k, id, state_handle, 1)
         p_two = model_pressure_t(grid(id)%we, j, k, id, state_handle, 1)
         model_pressure = interp_pressure(p_one, p_two, log_horz_interpM)
      else

         ! If not periodic, then try extrapolating
         p_one = model_pressure_t(i,   j, k, id, state_handle, 1)
         p_two = model_pressure_t(i+1, j, k, id, state_handle, 1)
         model_pressure = extrap_pressure(p_one, p_two, log_horz_interpM)

      endif

   else
       p_one = model_pressure_t(i,   j, k, id, state_handle, 1)
       p_two = model_pressure_t(i-1, j, k, id, state_handle, 1)
       model_pressure = interp_pressure(p_one, p_two, log_horz_interpM)
   endif

elseif (on_v_grid(state_id, var_id)) then ! average in the horizontal v direction

   if (j==grid(id)%sns) then

      ! Check to see if periodic in latitude (polar)
      if ( grid(id)%polar ) then

         ! The upper corner is 180 degrees of longitude away
         off = i + grid(id)%we/2
         if ( off > grid(id)%we ) off = off - grid(id)%we

         p_one = model_pressure_t(off, j-1, k, id, state_handle, 1)
         p_two = model_pressure_t(i,   j-1, k, id, state_handle, 1)
         model_pressure = interp_pressure(p_one, p_two, log_horz_interpM)
   
         ! If not periodic, then try extrapolating
      else

         p_one = model_pressure_t(i,j-1,k,id, state_handle, 1)
         p_two = model_pressure_t(i,j-2,k,id, state_handle, 1)
         model_pressure = extrap_pressure(p_one, p_two, log_horz_interpM)

      endif

    elseif (j==1) then

       ! Check to see if periodic in latitude (polar)
       if ( grid(id)%polar ) then

          ! The lower corner is 180 degrees of longitude away
          off = i + grid(id)%we/2
          if ( off > grid(id)%we ) off = off - grid(id)%we

          p_one = model_pressure_t(off,j,k,id, state_handle, 1)
          p_two = model_pressure_t(i,  j,k,id, state_handle, 1)
          model_pressure = interp_pressure(p_one, p_two, log_horz_interpM)

          ! If not periodic, then try extrapolating
       else
          p_one = model_pressure_t(i, j, k, id, state_handle, 1)
          p_two = model_pressure_t(i, j+1, k, id, state_handle, 1)
          model_pressure = extrap_pressure(p_one, p_two, log_horz_interpM)
       endif

    else
      p_one = model_pressure_t(i, j   ,k, id, state_handle, 1) 
      p_two = model_pressure_t(i, j-1, k, id, state_handle, 1)
      model_pressure = interp_pressure(p_one, p_two, log_horz_interpM)

    endif

else ! on t_grid

   p = model_pressure_t(i, j, k, id, state_handle, 1)
   model_pressure = p(1) 

endif

end function model_pressure

!------------------------------------------------------------------
! model surface pressure any grid u,v,w,t
function model_surface_pressure(i, j, id, var_id, state_id, state_handle)

! Calculate the surface pressure at grid point (i,j), domain id.
! The grid is defined according to var_type.

type(ensemble_type), intent(in) :: state_handle
integer,             intent(in) :: i, j, id
integer,             intent(in) :: var_id
integer,             intent(in) :: state_id
real(r8)              :: model_surface_pressure

integer  :: off
integer  :: ens_size ! only using the mean, so calling model_pressure_s with ens_size=1
real(r8) :: pres1(1), pres2(1)

ens_size = 1

model_surface_pressure = missing_r8

! If U-grid, then pressure is defined between U points, so average --
!   averaging depends on longitude periodicity
if( on_u_grid(state_id, var_id) ) then

   if (i==grid(id)%wes) then

      ! Check to see if periodic in longitude
      if ( grid(id)%periodic_x ) then

         ! We are at seam in longitude, take first and last M-grid points
         pres1 = model_pressure_s(i-1,j,id, state_handle, ens_size)
         pres2 = model_pressure_s(1,  j,id, state_handle, ens_size)
         model_surface_pressure = interp_pressure(pres1, pres2, log_horz_interpM)
         
      else

         ! If not periodic, then try extrapolating
         pres1 = model_pressure_s(i-1,j,id, state_handle, ens_size)
         pres2 = model_pressure_s(i-2,j,id, state_handle, ens_size)
         model_surface_pressure = extrap_pressure(pres1, pres2, log_horz_interpM)

      endif

   elseif( i == 1 ) then

      ! Check to see if periodic in longitude
      if ( grid(id)%periodic_x ) then

         ! We are at seam in longitude, take first and last M-grid points
         pres1 = model_pressure_s(i,          j,id, state_handle, ens_size)
         pres2 = model_pressure_s(grid(id)%we,j,id, state_handle, ens_size)
         model_surface_pressure = interp_pressure(pres1, pres2, log_horz_interpM)
         
      else

         ! If not periodic, then try extrapolating
         pres1 = model_pressure_s(i,  j,id, state_handle, ens_size)
         pres2 = model_pressure_s(i+1,j,id, state_handle, ens_size)
         model_surface_pressure = extrap_pressure(pres1, pres2, log_horz_interpM)

      endif

   else

      pres1 = model_pressure_s(i,  j,id, state_handle, ens_size)
      pres2 = model_pressure_s(i-1,j,id, state_handle, ens_size)
      model_surface_pressure = interp_pressure(pres1, pres2, log_horz_interpM)

   endif

! If V-grid, then pressure is defined between V points, so average --
!   averaging depends on polar periodicity
elseif( on_v_grid(state_id, var_id) ) then

   if (j==grid(id)%sns) then

      ! Check to see if periodic in latitude (polar)
      if ( grid(id)%polar ) then

         ! The upper corner is 180 degrees of longitude away
         off = i + grid(id)%we/2
         if ( off > grid(id)%we ) off = off - grid(id)%we

         pres1 = model_pressure_s(off,j-1,id, state_handle, ens_size)
         pres2 = model_pressure_s(i  ,j-1,id, state_handle, ens_size)
         model_surface_pressure = interp_pressure(pres1, pres2, log_horz_interpM)

      ! If not periodic, then try extrapolating
      else

         pres1 = model_pressure_s(i,j-1,id, state_handle, ens_size)
         pres2 = model_pressure_s(i,j-2,id, state_handle, ens_size)
         model_surface_pressure = extrap_pressure(pres1, pres2, log_horz_interpM)

      endif

   elseif( j == 1 ) then

      ! Check to see if periodic in latitude (polar)
      if ( grid(id)%polar ) then

         ! The lower corner is 180 degrees of longitude away
         off = i + grid(id)%we/2
         if ( off > grid(id)%we ) off = off - grid(id)%we

         pres1 = model_pressure_s(off,j,id, state_handle, ens_size)
         pres2 = model_pressure_s(i,  j,id, state_handle, ens_size)
         model_surface_pressure = interp_pressure(pres1, pres2, log_horz_interpM)

      ! If not periodic, then try extrapolating
      else

         pres1 = model_pressure_s(i,j,  id, state_handle, ens_size)
         pres2 = model_pressure_s(i,j+1,id, state_handle, ens_size)
         model_surface_pressure = extrap_pressure(pres1, pres2, log_horz_interpM)

      endif

   else

      pres1 = model_pressure_s(i,j,  id, state_handle, ens_size)
      pres2 = model_pressure_s(i,j-1,id, state_handle, ens_size)
      model_surface_pressure = interp_pressure(pres1, pres2, log_horz_interpM)

   endif

else

   pres1 = model_pressure_s(i,j,id, state_handle, ens_size)
   model_surface_pressure = pres1(1)

endif

end function model_surface_pressure

!------------------------------------------------------------------
function convert_indices_to_lon_lat_lev(i, j, k, var_id, state_id)

integer, intent(in) :: i, j, k, var_id, state_id
type(location_type) :: convert_indices_to_lon_lat_lev

real(r8) :: long, lat, lev
integer :: dom_id

dom_id = get_wrf_domain(state_id)

if ( on_u_grid(state_id, var_id) ) then
   long = grid(dom_id)%longitude_u(i,j)
   lat = grid(dom_id)%latitude_u(i,j)
elseif ( on_v_grid(state_id, var_id) ) then
   long = grid(dom_id)%longitude_v(i,j)
   lat = grid(dom_id)%latitude_v(i,j)
else ! on mass grid
   long = grid(dom_id)%longitude(i,j)
   lat = grid(dom_id)%latitude(i,j)
endif

! dart expects longitude [0,360]
do while (long <   0.0_r8)
   long = long + 360.0_r8
end do
do while (long > 360.0_r8)
   long = long - 360.0_r8
end do


if ( on_w_grid(state_id, var_id) ) then
   lev = real(k) - 0.5_r8
else
   lev = real(k)
endif

convert_indices_to_lon_lat_lev = set_location(long,lat,lev, VERTISLEVEL)

end function convert_indices_to_lon_lat_lev

!------------------------------------------------------------------
! which grid a variable is on.
!   querying dimension here, could do by qty?
!------------------------------------------------------------------
function on_u_grid(state_id, ivar)
integer, intent(in) :: state_id, ivar
logical :: on_u_grid

on_u_grid = (get_dim_name(state_id, ivar, 1) == 'west_east_stag')

end function

!------------------------------------------------------------------
function on_v_grid(state_id, ivar)
integer, intent(in) :: state_id, ivar
logical :: on_v_grid

on_v_grid = (get_dim_name(state_id, ivar, 2) == 'south_north_stag')

end function

!------------------------------------------------------------------
function on_w_grid(state_id, ivar)
integer, intent(in) :: state_id, ivar
logical :: on_w_grid

if (get_num_dims(state_id, ivar) > 2) then
   on_w_grid = (get_dim_name(state_id, ivar, 3) == 'bottom_top_stag')
else
   on_w_grid = .false.
endif

end function on_w_grid

!------------------------------------------------------------------
function on_t_grid(state_id, ivar)
integer, intent(in) :: state_id, ivar
logical :: on_t_grid

on_t_grid = (get_dim_name(state_id, ivar, 1) == 'west_east') .and. &
            (get_dim_name(state_id, ivar, 2) == 'south_north')

end function on_t_grid


!------------------------------------------------------------------
!------------------------------------------------------------------

function within_bounds_horizontal(i, j, id, qty_in)

integer, intent(in) :: i, j
integer, intent(in) :: id
integer, intent(in) :: qty_in
logical :: within_bounds_horizontal

integer :: var_id, qty

! Force QTY_TEMPERATURE to QTY_POTENTIAL_TEMPERATURE
if (qty_in == QTY_TEMPERATURE) then
   qty = QTY_POTENTIAL_TEMPERATURE
else
   qty = qty_in
endif

var_id = get_varid_from_kind(wrf_dom(id), qty) 

within_bounds_horizontal = (bounds_check_lon(i, id, var_id) .and. bounds_check_lat(j, id, var_id))

end function within_bounds_horizontal

!------------------------------------------------------------------
function able_to_interpolate_qty(id, qty)

integer, intent(in) :: id
integer, intent(in) :: qty

logical :: able_to_interpolate_qty

select case (qty)
   case (QTY_U_WIND_COMPONENT, QTY_V_WIND_COMPONENT)
      able_to_interpolate_qty = qty_in_domain(id, QTY_U_WIND_COMPONENT) .and. &
                                qty_in_domain(id, QTY_V_WIND_COMPONENT)

   case (QTY_10M_U_WIND_COMPONENT, QTY_10M_V_WIND_COMPONENT)
      able_to_interpolate_qty = qty_in_domain(id, QTY_10M_U_WIND_COMPONENT) .and. &
                                qty_in_domain(id, QTY_10M_V_WIND_COMPONENT)

   case (QTY_DENSITY)
      able_to_interpolate_qty = qty_in_domain(id, QTY_GEOPOTENTIAL_HEIGHT) .and. &
                                qty_in_domain(id, QTY_PRESSURE)

   case (QTY_SURFACE_TYPE)
      able_to_interpolate_qty = .true.  ! land mask XLAND is static data

   case (QTY_LANDMASK)
      able_to_interpolate_qty = .true.  ! land mask XLAND is static data

   case (QTY_SURFACE_ELEVATION)
      able_to_interpolate_qty = .true.  ! terrain height HGT is static data

   case (QTY_TEMPERATURE, QTY_POTENTIAL_TEMPERATURE)
      able_to_interpolate_qty = qty_in_domain(id, QTY_TEMPERATURE) .or. &
                                qty_in_domain(id, QTY_POTENTIAL_TEMPERATURE)

   case default
     able_to_interpolate_qty = qty_in_domain(id, qty)

end select


end function able_to_interpolate_qty

!------------------------------------------------------------------
function qty_in_domain(id, qty)

integer, intent(in) :: id
integer, intent(in) :: qty

logical :: qty_in_domain

integer :: varid

varid = get_varid_from_kind(wrf_dom(id), qty)

if (varid > 0) then
   qty_in_domain = .true.
else
   qty_in_domain = .false.
endif

end function qty_in_domain

!------------------------------------------------------------------
!  returns closest domain id and horizontal mass point grid points (iloc,jloc)
subroutine get_domain_info(obslon,obslat,id,iloc,jloc,domain_id_start)

real(r8), intent(in)           :: obslon, obslat
integer,  intent(out)          :: id
real(r8), intent(out)          :: iloc, jloc
integer,  intent(in), optional :: domain_id_start ! HK this is used in wrf_dart_obs_preprocess.f90

integer :: n ! domain to start from

if (present(domain_id_start)) then
  n = domain_id_start
else
  n = num_domains
endif

do id = n, 1, -1

! From module_map_utils.f90
!       latlon_to_ij(proj, lat, lon, i, j)
!       ij_to_latlon(proj, i, j, lat, lon)
!
!       It is incumbent upon the calling routine to determine whether or
!       not the values returned are within your domain's bounds.  All values
!       of i, j, lat, and lon are REAL values.

   call latlon_to_ij(grid(id)%proj,min(max(obslat,-89.9999999_r8),89.9999999_r8),obslon,iloc,jloc)

   if (found_in_domain(id, iloc,jloc)) return

enddo

! domain not found
id=0

end subroutine get_domain_info

!------------------------------------------------------------------
function found_in_domain(id, i,j)
integer,  intent(in) :: id
real(r8), intent(in) :: i, j
logical :: found_in_domain

found_in_domain = .true.

if (id > 1) then

   found_in_domain = ( i >= 1.0_r8 .and. i <= real(grid(id)%we,r8) .and. &
                       j >= 1.0_r8 .and. j <= real(grid(id)%sn,r8) )

else ! have to check periodic

! Array bound checking depends on whether periodic or not -- these are
!   real-valued indices here, so we cannot use boundsCheck  :(

   if ( grid(id)%periodic_x .and. .not. grid(id)%periodic_y  ) then
      if ( grid(id)%polar ) then
         !   Periodic     X & M_grid ==> [1 we+1)
         !   Periodic     Y & M_grid ==> [0.5 sn+0.5]
         found_in_domain = ( i >= 1.0_r8 .and. i <  real(grid(id)%we,r8)+1.0_r8 .and. &
                             j >= 0.5_r8 .and. j <= real(grid(id)%sn,r8)+0.5_r8 )
      else
         !   Periodic     X & M_grid ==> [1 we+1)
         !   NOT Periodic Y & M_grid ==> [1 sn]
         found_in_domain = ( i >= 1.0_r8 .and. i <  real(grid(id)%we,r8)+1.0_r8 .and. &
                             j >= 1.0_r8 .and. j <= real(grid(id)%sn,r8) )
      endif
   
   elseif ( grid(id)%periodic_x .and. grid(id)%periodic_y ) then
         !   Periodic     X & M_grid ==> [1 we+1)
         !   Periodic     Y & M_grid ==> [1 sn+1]
         found_in_domain = ( i >= 1.0_r8 .and. i <  real(grid(id)%we,r8)+1.0_r8 .and. &
                             j >= 1.0_r8 .and. j <= real(grid(id)%sn,r8)+1.0_r8 )
   
   else
      if ( grid(id)%polar ) then
         !   NOT Periodic X & M_grid ==> [1 we]
         !   Periodic     Y & M_grid ==> [0.5 sn+0.5]
         found_in_domain = ( i >= 1.0_r8 .and. i <= real(grid(id)%we,r8) .and. &
                             j >= 0.5_r8 .and. j <= real(grid(id)%sn,r8)+0.5_r8 )
      else
         !   NOT Periodic X & M_grid ==> [1 we]
         !   NOT Periodic Y & M_grid ==> [1 sn]
         found_in_domain = ( i >= 1.0_r8 .and. i <= real(grid(id)%we,r8) .and. &
                             j >= 1.0_r8 .and. j <= real(grid(id)%sn,r8) )
      endif
   endif
endif


end function found_in_domain

!------------------------------------------------------------------
subroutine getCorners(i, j, id, qty, ll, ul, lr, ur, rc)

integer, intent(in)  :: i, j, id, qty
integer, dimension(2), intent(out) :: ll, ul, lr, ur ! (x,y) of each corner
integer, intent(out) :: rc

integer :: var_id

! set return code to 0, and change this if necessary
rc = 0

var_id = get_varid_from_kind(wrf_dom(id), qty)


!----------------
! LOWER LEFT ll
!----------------
! i and j are the lower left (ll) corner already
!
! NOTE :: once we allow for polar periodicity, the incoming j index could actually
!           be 0, which would imply a ll(2) value of 1, with a ll(1) value 180 degrees
!           of longitude away from the incoming i index!  But we have not included
!           this possibility yet.

! As of 22 Oct 2007, this option is not allowed!
!   Note that j = 0 can only happen if we are on the M (or U) wrt to latitude
if ( grid(id)%polar .and. j == 0 .and. .not. restrict_polar ) then

   ! j = 0 should be mapped to j = 1 (ll is on other side of globe)
   ll(2) = 1
   
   ! Need to map i index 180 degrees away
   ll(1) = i + grid(id)%we/2
   
   ! Check validity of bounds & adjust by periodicity if necessary
   if ( ll(1) > grid(id)%we ) ll(1) = ll(1) - grid(id)%we

   ! We shouldn't be able to get this return code if restrict_polar = .true.
    rc = 1
    print*, 'model_mod.f90 :: getCorners :: Tried to do polar bc -- rc = ', rc

else
   
   ll(1) = i
   ll(2) = j

endif


!----------------
! LOWER RIGHT lr
!----------------

! Most of the time, the lower right (lr) corner will simply be (i+1,j), but we need to check
! Summary of x-direction corners:
!   Periodic     & M_grid has ind = [1 wes)
!     ind = [1 we)    ==> ind_p_1 = ind + 1
!     ind = [we wes)  ==> ind_p_1 = 1
!   Periodic     & U_grid has ind = [1 wes)
!     ind = [1 we)    ==> ind_p_1 = ind + 1
!     ind = [we wes)  ==> ind_p_1 = wes       ( keep in mind that U(1) = U(wes) if periodic )
!   NOT Periodic & M_grid has ind = [1 we)
!     ind = [1 we-1)  ==> ind_p_1 = ind + 1
!     ind = [we-1 we) ==> ind_p_1 = we
!   NOT Periodic & U_grid has ind = [1 wes)
!     ind = [1 we)    ==> ind_p_1 = ind + 1
!     ind = [we wes)  ==> ind_p_1 = wes

if ( grid(id)%periodic_x ) then
  
   ! Check to see what grid we have, M vs. U
   if (on_u_grid(wrf_dom(id), var_id) ) then
      ! U-grid is always i+1 -- do this in reference to already adjusted ll points
      lr(1) = ll(1) + 1
      lr(2) = ll(2)
   else
      ! M-grid is i+1 except if we <= ind < wes, in which case it's 1
      if ( i < grid(id)%we ) then
         lr(1) = ll(1) + 1
      else
         lr(1) = 1
      endif
      lr(2) = ll(2)
   endif

else

  ! Regardless of grid, NOT Periodic always has i+1
   lr(1) = ll(1) + 1
   lr(2) = ll(2)

endif
      

!----------------
! UPPER LEFT ul
!----------------

!** NOTE: For now are disallowing observation locations that occur poleward of the
!           first and last M-grid gridpoints.  This need not be the case because
!           the information should be available for proper interpolation across the
!           poles, but it will require more clever thinking.  Hopefully this can
!           be added in later.

! Most of the time, the upper left (ul) corner will simply be (i,j+1), but we need to check
! Summary of y-direction corners:
!   Periodic     & M_grid has ind = [0 sns)  *though in practice, [1 sn)*
!     ind = [1 sn-1)  ==> ind_p_1 = ind + 1
!     ind = [sn-1 sn) ==> ind_p_1 = sn
!   Periodic     & V_grid has ind = [1 sns)
!     ind = [1 sn)    ==> ind_p_1 = ind + 1
!     ind = [sn sns)  ==> ind_p_1 = sns
!   NOT Periodic & M_grid has ind = [1 sn)
!     ind = [1 sn-1)  ==> ind_p_1 = ind + 1
!     ind = [sn-1 sn) ==> ind_p_1 = sn
!   NOT Periodic & V_grid has ind = [1 sns)
!     ind = [1 sn)    ==> ind_p_1 = ind + 1
!     ind = [sn sns)  ==> ind_p_1 = sns
!
! Hence, with our current polar obs restrictions, all four possible cases DO map into
!   ul = (i,j+1).  But this will not always be the case.

if ( grid(id)%polar ) then

   ! Check to see what grid we have, M vs. V
   if ( on_v_grid(wrf_dom(id), var_id) ) then
      ! V-grid is always j+1, even if we allow for full [1 sns) range
      ul(1) = ll(1)
      ul(2) = ll(2) + 1
   else
      ! M-grid changes depending on polar restriction
      if ( restrict_polar ) then
         ! If restricted, then we can simply add 1
         ul(1) = ll(1)
         ul(2) = ll(2) + 1
      else
         ! If not restricted, then we can potentially wrap over the north pole, which
         !   means that ul(2) is set to sn and ul(1) is shifted by 180 deg.

         if ( j == grid(id)%sn ) then
            ! j > sn should be mapped to j = sn (ul is on other side of globe)
            ul(2) = grid(id)%sn
   
            ! Need to map i index 180 degrees away
            ul(1) = ll(1) + grid(id)%we/2
   
            ! Check validity of bounds & adjust by periodicity if necessary
            if ( ul(1) > grid(id)%we ) ul(1) = ul(1) - grid(id)%we

            ! We shouldn't be able to get this return code if restrict_polar = .true.
             rc = 1
             print*, 'model_mod.f90 :: getCorners :: Tried to do polar bc -- rc = ', rc

         elseif ( j == 0 ) then
            ! In this case, we have place ll on the other side of the globe, so we
            !   cannot reference ul to ll
            ul(1) = i
            ul(2) = 1

         else
            ! We can confidently set to j+1
            ul(1) = ll(1)
            ul(2) = ll(2) + 1
         endif

      endif
   endif

elseif ( grid(id)%periodic_y ) then

   ! Check to see what grid we have, M vs. V
   if ( on_v_grid(wrf_dom(id), var_id) ) then
      ! V-grid is always j+1 -- do this in reference to already adjusted ll points
      ul(1) = ll(1)
      ul(2) = ll(2)+1
   else
      ! M-grid is j+1 except if we <= ind < wes, in which case it's 1
      if ( j < grid(id)%sn ) then
         ul(2) = ll(2) + 1
      else
         ul(2) = 1
      endif
      ul(1) = ll(1)
   endif

else

   ! Regardless of grid, NOT Periodic always has j+1
   ul(1) = ll(1)
   ul(2) = ll(2) + 1

endif
   

!----------------
! UPPER RIGHT ur
!----------------

!*** NOTE: For now are disallowing observation locations that occur poleward of the
!            first and last M-grid gridpoints.  This need not be the case because
!            the information should be available for proper interpolation across the
!            poles, but it will require more clever thinking.  Hopefully this can
!            be added in later.

! Most of the time, the upper right (ur) corner will simply be (i+1,j+1), but we need to check
!   In fact, we can largely get away with ur = (lr(1),ul(2)).  Where this will NOT work is
!   where we have had to re-map the i index to the other side of the globe (180 deg) due to
!   the polar boundary condition.  There are no situations where ur(2) will not be equal to
!   ul(2).

ur(2) = ul(2)

! Need to check if ur(1) .ne. lr(1)
if ( grid(id)%polar .and. .not. restrict_polar ) then

   ! Only if j == 0 or j == sn
   if ( j == 0 .or. j ==  grid(id)%sn) then
      ! j == 0 means that ll(1) = i + 180 deg, so we cannot use lr(1) -- hence, we will
      !   add 1 to ul(1), unless doing so spans the longitude seam point.
      ! j == sn means that ul(1) = i + 180 deg.  Here we cannot use lr(1) either because
      !   it will be half a domain away from ul(1)+1.  Be careful of longitude seam point.

      !   Here we need to check longitude periodicity and the type of grid
      if ( grid(id)%periodic_x ) then
  
         ! Check to see what grid we have, M vs. U
         if ( on_u_grid(wrf_dom(id), var_id) ) then
            ! U-grid is always i+1 -- do this in reference to already adjusted ll points
            ur(1) = ul(1) + 1
         else
            ! M-grid is i+1 except if we <= ind < wes, in which case it's 1
            if ( ul(1) < grid(id)%we ) then
               ur(1) = ul(1) + 1
            else
               ur(1) = 1
            endif
         endif

      else

         ! Regardless of grid, NOT Periodic always has i+1
         ur(1) = ul(1) + 1

      endif

   ! If not a special j value, then we are set for the ur(1) = lr(1)
   else

      ur(1) = lr(1)

   endif

! If not an unrestricted polar periodic domain, then we have nothing to worry about
else

   ur(1) = lr(1)

endif

end subroutine getCorners


!------------------------------------------------------------------
!  Transfer obs. x to grid j and calculate its
!  distance to grid j and j+1
subroutine toGrid (x, j, dx, dxm)

real(r8), intent(in)  :: x
real(r8), intent(out) :: dx, dxm
integer,  intent(out) :: j

j = int(x)
dx = x - real(j)  ! HK todo this should be real(j, KIND=r8)?
dxm= 1.0_r8 - dx

end subroutine toGrid

!------------------------------------------------------------------
subroutine setup_map_projection(id)

integer, intent(in)   :: id
logical, parameter    :: debug = .false.

real(r8) :: latinc,loninc,stdlon
real(r8) :: truelat1, truelat2

call map_init(grid(id)%proj)

! Populate the map projection structure

!nc -- added in case structures for CASSINI and CYL
!nc -- global wrfinput_d01 has truelat1 = 1.e20, so we need to change it where needed
!nc -- for PROJ_LATLON stdlon and truelat1 have different meanings --
!nc --   stdlon --> loninc  and  truelat1 --> latinc
!JPH --- this latinc/loninc calculations are only valid for global domains

!if ( wrf%dom(id)%scm ) then
!! JPH -- set to zero which should cause the map utils to return NaNs if called
!   latinc = 0.0_r8
!   loninc = 0.0_r8
!else
!   latinc = 180.0_r8/wrf%dom(id)%sn
!   loninc = 360.0_r8/wrf%dom(id)%we
!endif

latinc = 180.0_r8/size(grid(id)%longitude(:,1)) ! west_east
loninc = 360.0_r8/size(grid(id)%longitude(1,:)) ! north_south

if(grid(id)%map_proj == PROJ_LATLON) then !HK why are these different to the wrfinput file?
   truelat1 = latinc
   stdlon = loninc
else
   truelat1 = grid(id)%truelat1
   truelat2 = grid(id)%truelat2
   stdlon = grid(id)%stand_lon
endif

!nc -- specified inputs to hopefully handle ALL map projections -- hopefully map_set will
!        just ignore the inputs it doesn't need for its map projection of interest (?)
!
!   NOTE:: We are NOT yet supporting the Gaussian grid or the Rotated Lat/Lon, so we
!            are going to skip the entries:  nlon, nlat, ixdim, jydim, stagger, phi, lambda
!
!      + Gaussian grid uses nlat & nlon
!      + Rotated Lat/Lon uses ixdim, jydim, stagger, phi, & lambda
!
call map_set( proj_code=grid(id)%map_proj, &
              proj=grid(id)%proj, &
              lat1=grid(id)%latitude(1,1), &
              lon1=grid(id)%longitude(1,1), &
              lat0=90.0_r8, &
              lon0=0.0_r8, &
              knowni=1.0_r8, &
              knownj=1.0_r8, &
              dx=grid(id)%dx, &
              latinc=latinc, &
              loninc=loninc, &
              stdlon=stdlon, &
              truelat1=truelat1, &
              truelat2=truelat2  )

end subroutine setup_map_projection


!------------------------------------------------------------------
!------------------------------------------------------------------
! Bounds checking routines:
  ! Summary of Allowable REAL-VALUED Index Values ==> INTEGER Index Values 
  !
  ! In longitude (x) direction
  !   Periodic     & M_grid ==> [1 we+1)       ==> [1 wes)
  !   Periodic     & U_grid ==> [1 wes)        ==> [1 wes)
  !   NOT Periodic & M_grid ==> [1 we]         ==> [1 we)
  !   NOT Periodic & U_grid ==> [1.5 we+0.5]   ==> [1 wes)
  ! In latitude (y) direction
  !   Periodic     & M_grid ==> [0.5 sn+0.5]   ==> [0 sns) *though in practice, [1 sn)*
  !   Periodic     & V_grid ==> [1 sns]        ==> [1 sns) *though allowable range, [1.5 sn+.5]*
  !   NOT Periodic & M_grid ==> [1 sn]         ==> [1 sn)
  !   NOT Periodic & V_grid ==> [1.5 sn+0.5]   ==> [1 sns)
  ! In vertical (z) direction
  !                  M_grid ==> [1 bt]         ==> [1 bt)
  !                  W_grid ==> [1.5 bt+0.5]   ==> [1 bts)

!------------------------------------------------------------------

!------------------------------------------------------------------
! Determines whether real-valued location indices are
! within a sensible range based on the assumed (un)staggered grid and based on 
! whether the domain is assumed to be periodic in a given direction.
function bounds_check_lon(ind, id, var_id)  ! or variable?

integer,  intent(in)  :: ind ! index into the grid
integer,  intent(in)  :: id  ! domain id
integer,  intent(in)  :: var_id 

logical :: bounds_check_lon

! Consider cases in REAL-VALUED indexing:
!
! I. Longitude -- x-direction
!    A. PERIODIC (period_x = .true.)
!
!       Consider Mass-grid (& V-grid) longitude grid with 4 west-east gridpoints
!         Values  ::  [ -135 -45  45 135 ] .. {225}
!         Indices ::  [   1   2   3   4  ] .. {1,5}
!       Complementary U-grid
!         Values  ::  [ -180 -90  0  90  180 ]
!         Indices ::  [   1   2   3   4   5  ]
!
!       What are the allowable values for a real-valued index on each of these grids?
!       1. M-grid  --->  [1 5)       ---> [1 we+1)
!                  --->  [-135 225)  
!       2. U-grid  --->  [1 5)       ---> [1 wes)
!                  --->  [-180 180)
!       [Note that above "allowable values" reflect that one should be able to have
!        an observation anywhere on a given longitude circle -- the information 
!        exists in order to successfully interpolate to anywhere over [0 360).]
!
!       It is up to the routine calling "boundsCheck" to have handled the 0.5 offset
!         in indices between the M-grid & U-grid.  Hence, two examples: 
!          a. If there is an observation location at -165 longitude, then:
!             * An observation of TYPE_T (on the M-grid) would have ind = 4.667
!             * An observation of TYPE_U (on the U-grid) would have ind = 1.167
!          b. If there is an observation location at 0 longitude, then:
!             * An observation of TYPE_T (on the M-grid) would have ind = 2.5
!             * An observation of TYPE_U (on the U-grid) would have ind = 3.0
!
!    B. NOT periodic (period_x = .false.)
!
!       Consider Mass-grid (& V-grid) longitude grid with 4 west-east gridpoints
!         Values  ::  [  95  105 115 125 ] 
!         Indices ::  [   1   2   3   4  ] 
!       Complementary U-grid
!         Values  ::  [  90  100 110 120 130 ]
!         Indices ::  [   1   2   3   4   5  ]
!
!       What are the allowable values for a real-valued index on each of these grids?
!       1. M-grid  --->  [1 4]       ---> [1 we]
!                  --->  [95 125]  
!       2. U-grid  --->  [1.5 4.5]       ---> [1.5 we+0.5]
!                  --->  [95 125]
!       [Note that above "allowable values" reflect that one should only be able to
!        have an observation within the M-grid, since that is the only way to  
!        guarantee that the necessary information exists in order to successfully 
!        interpolate to a specified location.]
!
!       It is up to the routine calling "boundsCheck" to have handled the 0.5 offset
!         in indices between the M-grid & U-grid.  Hence, two examples: 
!          a. If there is an observation location at 96 longitude, then:
!             * An observation of TYPE_T (on the M-grid) would have ind = 1.1
!             * An observation of TYPE_U (on the U-grid) would have ind = 1.6
!          b. If there is an observation location at 124 longitude, then:
!             * An observation of TYPE_T (on the M-grid) would have ind = 3.9
!             * An observation of TYPE_U (on the U-grid) would have ind = 4.4
!

bounds_check_lon = .false.

! Next check periodicity
if ( grid(id)%periodic_x ) then
   
   ! If periodic in longitude, then no need to check staggering because both
   !   M and U grids allow integer indices from [1 wes)
   if ( ind >= 1 .and. ind < grid(id)%wes ) bounds_check_lon = .true.

else

   ! If NOT periodic in longitude, then we need to check staggering because
   !   M and U grids allow different index ranges

   ! Check staggering by comparing var_size(dim,type) to the staggered dimension 
   if ( on_u_grid(wrf_dom(id),var_id) ) then
      ! U-grid allows integer range of [1 wes)
      if ( ind >= 1 .and. ind < grid(id)%wes ) bounds_check_lon = .true.
   else  
      ! M & V-grid allow [1 we)
      if ( ind >= 1 .and. ind < grid(id)%we ) bounds_check_lon = .true.
   endif

endif

end function bounds_check_lon

!------------------------------------------------------------------
function bounds_check_lat(ind, id, var_id)

integer,  intent(in)  :: ind ! index into the grid
integer,  intent(in)  :: id  ! domain id
integer,  intent(in)  :: var_id  

logical :: bounds_check_lat
! II. Latitude -- y-direction
!    A. PERIODIC (polar = .true.)
!
!       Consider Mass-grid (& U-Grid) latitude grid with 4 south-north gridpoints
!         Values  :: [ -67.5 -22.5  22.5  67.5 ] 
!         Indices :: [   1     2     3     4   ] 
!       Complementary V-grid 
!         Values  :: [ -90   -45     0    45    90 ] 
!         Indices :: [   1     2     3     4     5 ] 
!
!       What are the allowable values for a real-valued index on each of these grids?
!       1. M-grid  --->  [0.5 4.5]   ---> [0.5 sn+0.5]
!                  --->  [-90 90]  
!       2. U-grid  --->  [1 5]       ---> [1 sns]
!                  --->  [-90 90]
!       [Note that above "allowable values" reflect that one should be able to have
!        an observation anywhere along a give latitude circle -- the information 
!        exists in order to successfully interpolate to anywhere over [-90 90]; 
!        however, in latitude this poses a special challenge since the seams join
!        two separate columns of data over the pole, as opposed to in longitude
!        where the seam wraps back on a single row of data.]  
!
!       It is up to the routine calling "boundsCheck" to have handled the 0.5 offset
!         in indices between the M-grid & V-grid.  Hence, two examples: 
!          a. If there is an observation location at -75 latitude, then:
!             * An observation of TYPE_T (on the M-grid) would have ind = 0.833
!             * An observation of TYPE_V (on the V-grid) would have ind = 1.333
!          b. If there is an observation location at 0 latitude, then:
!             * An observation of TYPE_T (on the M-grid) would have ind = 2.5
!             * An observation of TYPE_V (on the V-grid) would have ind = 3.0
!
!    B. NOT periodic (polar = .false.)
!
!       Consider Mass-grid (& U-Grid) latitude grid with 4 south-north gridpoints
!         Values  :: [ 10  20  30  40 ] 
!         Indices :: [  1   2   3   4 ] 
!       Complementary V-grid 
!         Values  :: [  5  15  25  35  45 ] 
!         Indices :: [  1   2   3   4   5 ] 
!
!       What are the allowable values for a real-valued index on each of these grids?
!       1. M-grid  --->  [1 4]   ---> [1 sn]
!                  --->  [10 40]  
!       2. U-grid  --->  [1.5 4.5]       ---> [1.5 sn+0.5]
!                  --->  [10 40]
!       [Note that above "allowable values" reflect that one should only be able to
!        have an observation within the M-grid, since that is the only way to  
!        guarantee that the necessary information exists in order to successfully 
!        interpolate to a specified location.]
!
!       It is up to the routine calling "boundsCheck" to have handled the 0.5 offset
!         in indices between the M-grid & V-grid.  Hence, two examples: 
!          a. If there is an observation location at 11 latitude, then:
!             * An observation of TYPE_T (on the M-grid) would have ind = 1.1
!             * An observation of TYPE_V (on the V-grid) would have ind = 1.6
!          b. If there is an observation location at 25 latitude, then:
!             * An observation of TYPE_T (on the M-grid) would have ind = 2.5
!             * An observation of TYPE_V (on the V-grid) would have ind = 3.0
!

if ( grid(id)%periodic_y ) then
        
     ! We need to check staggering because M and V grids allow different indices

!*** NOTE: For now are disallowing observation locations that occur poleward of the 
!            first and last M-grid gridpoints.  This means that this function will 
!            return false for polar observations.  This need not be the case because
!            the information should be available for proper interpolation across the
!            poles, but it will require more clever thinking.  Hopefully this can 
!            be added in later.  

    ! Check staggering by comparing var_size(dim,type) to the staggered dimension 
    if ( on_v_grid(wrf_dom(id), var_id) ) then
       ! V-grid allows integer range [1 sns)
       if ( ind >= 1 .and. ind < grid(id)%sns ) bounds_check_lat = .true.
    else  
       ! For now we will set a logical flag to more restrictively check the array
       !   bounds under our no-polar-obs assumptions
       if ( restrict_polar ) then
          ! M & U-grid allow integer range [1 sn) in practice (though properly, [0 sns) )
          if ( ind >= 1 .and. ind < grid(id)%sn ) bounds_check_lat = .true.
       else
          ! M & U-grid allow integer range [0 sns) in unrestricted circumstances
          if ( ind >= 0 .and. ind < grid(id)%sns ) bounds_check_lat = .true.
       endif
    endif
    
 else

    ! We need to check staggering because M and V grids allow different indices
    if ( on_v_grid(wrf_dom(id), var_id) ) then
       ! V-grid allows [1 sns)
       if ( ind >= 1 .and. ind < grid(id)%sns ) bounds_check_lat = .true.
    else 
       ! M & U-grid allow [1 sn)
       if ( ind >= 1 .and. ind < grid(id)%sn ) bounds_check_lat = .true.
    endif

 endif

end function bounds_check_lat


!------------------------------------------------------------------
function bounds_check_vertical(ind, id, var_id) 

integer,  intent(in)  :: ind ! index into the grid
integer,  intent(in)  :: id  ! domain id
integer,  intent(in)  :: var_id

logical :: bounds_check_vertical
! III. Vertical -- z-direction (periodicity not an issue)
!    
!    Consider Mass vertical grid with 4 bottom-top gridpoints
!      Values  :: [ 0.875 0.625 0.375 0.125 ]
!      Indices :: [   1     2     3     4   ]
!    Complementary W-grid
!      Values  :: [   1   0.75  0.50  0.25    0   ]
!      Indices :: [   1     2     3     4     5   ]
!
!    What are the allowable values for a real-valued index on each of these grids?
!    1. M-grid  --->  [1 4]           ---> [1 bt]
!               --->  [0.875 0.125]  
!    2. W-grid  --->  [1.5 4.5]       ---> [1.5 bt+0.5]
!               --->  [0.875 0.125]
!
!    [Note that above "allowable values" reflect that one should only be able to
!     have an observation within the M-grid, since that is the only way to  
!     guarantee that the necessary information exists in order to successfully 
!     interpolate to a specified location.]

bounds_check_vertical = .false.

if ( on_w_grid(wrf_dom(id), var_id) ) then
   ! W vertical grid allows [1 bts)
   if ( ind >= 1 .and. ind < grid(id)%bts ) bounds_check_vertical = .true.
else
   ! M vertical grid allows [1 bt)
   if ( ind >= 1 .and. ind < grid(id)%bt ) bounds_check_vertical = .true.
endif
  
end function bounds_check_vertical 

!------------------------------------------------------------------
! Summary of Allowable REAL-VALUED Index Values ==> INTEGER Index Values 
!
! In longitude (x) direction
!   Periodic     & M_grid ==> [1 we+1)       ==> [1 wes)
!   Periodic     & U_grid ==> [1 wes)        ==> [1 wes)
!   NOT Periodic & M_grid ==> [1 we]         ==> [1 we)
!   NOT Periodic & U_grid ==> [1.5 we+0.5]   ==> [1 wes)
! In latitude (y) direction
!   Periodic     & M_grid ==> [0.5 sn+0.5]   ==> [0 sns) *though in practice, [1 sn)*
!   Periodic     & V_grid ==> [1 sns]        ==> [1 sns) *though allowable range, [1.5 sn+.5]*
!   NOT Periodic & M_grid ==> [1 sn]         ==> [1 sn)
!   NOT Periodic & V_grid ==> [1.5 sn+0.5]   ==> [1 sns)
! In vertical (z) direction
!                  M_grid ==> [1 bt]         ==> [1 bt)
!                  W_grid ==> [1.5 bt+0.5]   ==> [1 bts)
!------------------------------------------------------------------

!------------------------------------------------------------------
! Vortex

subroutine vortex()

print*, 'Do vortex'

end subroutine vortex

end module model_mod
