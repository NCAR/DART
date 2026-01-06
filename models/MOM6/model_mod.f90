! DART software - Copyright UCAR. This open source software is provided
! by UCAR, "as is", without charge, subject to all terms of use at
! http://www.image.ucar.edu/DAReS/DART/DART_download
!

! @todo
!   * QUAD_LOCATED_CELL_CENTERS - what difference does this make?
!   * t_grid interp for thickess locate and evaluate
!   * vertical location for QTY_DRY_LAND

module model_mod

use        types_mod, only : r8, i8, MISSING_R8, vtablenamelength

use time_manager_mod, only : time_type, set_time, set_date, set_calendar_type

use     location_mod, only : location_type, get_close_type, &
                             loc_get_close_obs => get_close_obs, &
                             loc_get_close_state => get_close_state, &
                             set_location, set_location_missing, &
                             get_location, query_location, VERTISLEVEL, &
                             VERTISHEIGHT, set_vertical, get_dist

use    utilities_mod, only : error_handler, E_ERR, E_MSG, &
                             nmlfileunit, do_output, do_nml_file, do_nml_term,  &
                             find_namelist_in_file, check_namelist_read

use netcdf_utilities_mod, only : nc_add_global_attribute, nc_synchronize_file, &
                                 nc_add_global_creation_time, &
                                 nc_begin_define_mode, nc_end_define_mode, &
                                 nc_open_file_readonly, nc_close_file, &
                                 nc_get_variable, nc_get_variable_size, &
                                 NF90_MAX_NAME, nc_get_attribute_from_variable

use        quad_utils_mod,  only : quad_interp_handle, init_quad_interp, &
                                   set_quad_coords, quad_lon_lat_locate, &
                                   quad_lon_lat_evaluate, quad_interp_handle, &
                                   GRID_QUAD_FULLY_IRREGULAR, &
                                   QUAD_LOCATED_CELL_CENTERS

use state_structure_mod, only : add_domain, get_domain_size, &
                                get_model_variable_indices, &
                                get_kind_string, get_varid_from_kind, &
                                get_dart_vector_index

use distributed_state_mod, only : get_state, get_state_array

use obs_kind_mod, only : get_index_for_quantity, QTY_U_CURRENT_COMPONENT, &
                         QTY_V_CURRENT_COMPONENT, QTY_LAYER_THICKNESS, &
                         QTY_DRY_LAND, QTY_SALINITY, QTY_TEMPERATURE, &
                         QTY_POTENTIAL_TEMPERATURE

use ensemble_manager_mod, only : ensemble_type

! These routines are passed through from default_model_mod.
! To write model specific versions of these routines
! remove the routine from this use statement and add your code to
! this the file.
use default_model_mod, only : pert_model_copies, write_model_time, &
                              init_time => fail_init_time, &
                              init_conditions => fail_init_conditions, &
                              convert_vertical_obs, adv_1step, &
                              parse_variables, &
                              MAX_STATE_VARIABLE_FIELDS

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


character(len=256), parameter :: source   = "MOM6/model_mod.f90"
logical :: module_initialized = .false.
type(time_type) :: assimilation_time_step
integer :: dom_id ! used to access the state structure
integer(i8) :: model_size
integer :: nfields ! number of fields in the state vector
! Grid parameters, nz is number of layers
integer :: nx=-1, ny=-1, nz=-1     ! grid counts for each field
real(r8), allocatable :: geolon(:,:), geolat(:,:),     & ! T
                         geolon_u(:,:), geolat_u(:,:), & ! U
                         geolon_v(:,:), geolat_v(:,:)    ! V
logical, allocatable  :: mask(:,:), mask_u(:,:), mask_v(:,:) ! geolat/lon/u/v has missing values
type(quad_interp_handle) :: interp_t_grid,  &
                            interp_u_grid,  &
                            interp_v_grid

! Ocean vertical
real(r8), allocatable :: zstar(:) ! pseudo depth for each layer

! Ocean vs land
real(r8), allocatable :: wet(:,:), basin_depth(:,:)

! model_interpolate failure codes
integer, parameter :: NOT_IN_STATE = 12
integer, parameter :: THICKNESS_NOT_IN_STATE = 13
integer, parameter :: QUAD_LOCATE_FAILED = 14
integer, parameter :: THICKNESS_QUAD_EVALUATE_FAILED = 15
integer, parameter :: QUAD_EVALUATE_FAILED = 16
integer, parameter :: QUAD_ON_LAND = 17
integer, parameter :: QUAD_ON_BASIN_EDGE = 18
integer, parameter :: OBS_ABOVE_SURFACE = 20
integer, parameter :: OBS_TOO_DEEP = 22


! namelist
character(len=256) :: template_file = 'mom6.r.nc'
character(len=256) :: static_file = 'c.e22.GMOM.T62_g16.nuopc.001.mom6.static.nc'
character(len=256) :: ocean_geometry = 'ocean_geometry.nc'
integer  :: assimilation_period_days      = 1
integer  :: assimilation_period_seconds   = 0
character(len=vtablenamelength) :: model_state_variables(MAX_STATE_VARIABLE_FIELDS) = ' '
character(len=NF90_MAX_NAME) :: layer_name = 'Layer'
logical :: use_pseudo_depth = .false. ! use pseudo depth instead of sum(layer thickness) for vertical location

namelist /model_nml/ template_file, static_file, ocean_geometry, assimilation_period_days, &
                     assimilation_period_seconds, model_state_variables, layer_name, &
                     use_pseudo_depth


interface on_land
   module procedure on_land_point
   module procedure on_land_quad
end interface on_land

contains

!------------------------------------------------------------------
!
! Called to do one time initialization of the model. As examples,
! might define information about the model size or model timestep.
! In models that require pre-computed static data, for instance
! spherical harmonic weights, these would also be computed here.

subroutine static_init_model()

integer  :: iunit, io

module_initialized = .true.

call find_namelist_in_file("input.nml", "model_nml", iunit)
read(iunit, nml = model_nml, iostat = io)
call check_namelist_read(iunit, io, "model_nml")

! Record the namelist values used for the run
if (do_nml_file()) write(nmlfileunit, nml=model_nml)
if (do_nml_term()) write(     *     , nml=model_nml)

call set_calendar_type('gregorian')

! This time is both the minimum time you can ask the model to advance
! (for models that can be advanced by filter) and it sets the assimilation
! window.  All observations within +/- 1/2 this interval from the current
! model time will be assimilated. If this is not settable at runtime 
! feel free to hardcode it and remove from the namelist.
assimilation_time_step = set_time(assimilation_period_seconds, &
                                  assimilation_period_days)

! Define which variables are in the model state;
! parse_variables converts the character table that was read in from
! model_nml:model_state_variables to a state_var_type that can be passed
! to add_domain
dom_id = add_domain(template_file, parse_variables(model_state_variables))

model_size = get_domain_size(dom_id)

call read_horizontal_grid()
call setup_interpolation()

call read_num_layers()
call read_ocean_geometry() ! ocean vs. land and basin depth

end subroutine static_init_model

!------------------------------------------------------------------
! Returns the number of items in the state vector as an integer. 
function get_model_size()

integer(i8) :: get_model_size

if ( .not. module_initialized ) call static_init_model

get_model_size = get_domain_size(dom_id)

end function get_model_size


!------------------------------------------------------------------
! Given a state handle, a location, and a state quantity,
! interpolates the state variable fields to that location and returns
! the values in expected_obs. The istatus variables should be returned as
! 0 unless there is some problem in computing the interpolation in
! which case a positive istatus should be returned.

subroutine model_interpolate(state_handle, ens_size, location, qty_in, expected_obs, istatus)


type(ensemble_type), intent(in) :: state_handle
integer,             intent(in) :: ens_size
type(location_type), intent(in) :: location
integer,             intent(in) :: qty_in
real(r8),           intent(out) :: expected_obs(ens_size) !< array of interpolated values
integer,            intent(out) :: istatus(ens_size)

real(r8), parameter             :: CONCENTRATION_TO_PPT = 1000.0_r8

integer  :: qty ! local qty
integer  :: which_vert, four_ilons(4), four_ilats(4)
integer  :: lev(ens_size,2), levz(2) ! level below and above obs
integer  :: locate_status, quad_status
real(r8) :: lev_fract(ens_size), levz_fract ! fraction between bottom and top level
real(r8) :: lon_lat_vert(3)
real(r8) :: quad_vals(4, ens_size)
real(r8) :: expected(ens_size, 2)
real(r8) :: expected_pot_temp(ens_size), expected_salinity(ens_size), pressure_bars(ens_size)
type(quad_interp_handle) :: interp
integer :: varid, i, e, thick_id, corner
integer(i8) :: th_indx
real(r8) :: depth_at_x(ens_size), thick_at_x(ens_size) ! depth, layer thickness at obs lat lon
logical :: found(ens_size)

if ( .not. module_initialized ) call static_init_model

expected_obs(:) = MISSING_R8
istatus(:) = 1

if (qty_in == QTY_TEMPERATURE) then
   qty = QTY_POTENTIAL_TEMPERATURE  ! model has potential temperature
   if (get_varid_from_kind(dom_id, QTY_SALINITY) < 0) then ! Require salinity to convert to temperature
      istatus = NOT_IN_STATE
      return
   end if
else
   qty = qty_in
endif

varid = get_varid_from_kind(dom_id, qty)
if (varid < 0) then ! not in state
   istatus = NOT_IN_STATE
   return
endif

if (.not. use_pseudo_depth) then
   thick_id = get_varid_from_kind(dom_id, QTY_LAYER_THICKNESS)
   if (thick_id < 0) then ! thickness not in state
      istatus = THICKNESS_NOT_IN_STATE
      return
   endif
endif

! find which grid the qty is on
interp = get_interp_handle(qty)

! unpack the location type into lon, lat, vert, vert_type
lon_lat_vert = get_location(location)
which_vert   = nint(query_location(location))

! get the indices for the 4 corners of the quad in the horizontal
call quad_lon_lat_locate(interp, lon_lat_vert(1), lon_lat_vert(2), &
                         four_ilons, four_ilats, locate_status)
if (locate_status /= 0) then
  istatus(:) = QUAD_LOCATE_FAILED
  return
endif

! check if all four corners are in the ocean
if (on_land(four_ilons, four_ilats)) then
   istatus(:) = QUAD_ON_LAND
   return
endif

if (use_pseudo_depth) then 

   ! Get the bounding vertical levels and the fraction between bottom and top
   call find_level_bounds(lon_lat_vert(3), which_vert, levz, levz_fract, locate_status)
   if (locate_status /= 0) then
      istatus(:) = locate_status
      return
   endif   
   ! pseudo depth is the same for all ensemble members
   lev(:,1) = levz(1) ! layer_below
   lev(:,2) = levz(2) ! layer_above
   lev_fract(:) = levz_fract
   depth_at_x(:) = zstar(levz(1)) ! pseudo depth at obs lat lon
else

   ! find which layer the observation is in. Layer thickness is a state variable.
   ! HK @todo Do you need to use t_grid interp for thickness four_ilons, four_ilats?
   found(:) = .false.
   depth_at_x(:) = 0
   FIND_LAYER: do i = 2, nz

      do corner = 1, 4
         th_indx = get_dart_vector_index(four_ilons(corner), four_ilats(corner), i, dom_id, thick_id)
          quad_vals(corner, :) = get_state(th_indx, state_handle)
      enddo

      
      call quad_lon_lat_evaluate(interp, &
                                 lon_lat_vert(1), lon_lat_vert(2), & ! lon, lat of obs
                                 four_ilons, four_ilats, &
                                 ens_size, &
                                 quad_vals, & ! 4 corners x ens_size
                                 thick_at_x, &
                                 quad_status)
      if (quad_status /= 0) then
         istatus(:) = THICKNESS_QUAD_EVALUATE_FAILED
         return
      endif

      depth_at_x = depth_at_x + thick_at_x

      do e = 1, ens_size
         if (lon_lat_vert(3) < depth_at_x(e)) then
            lev(e,1) = i ! layer_below
            lev(e,2) = i-1 ! layer_above
            lev_fract(e) = (depth_at_x(e) - lon_lat_vert(3)) / thick_at_x(e)
            found(e) = .true.
            if (all(found)) exit FIND_LAYER
         endif
      enddo

   enddo FIND_LAYER

   if (any(found .eqv. .false.)) then
      istatus(:) = OBS_TOO_DEEP
      return
   endif

endif

if (on_basin_edge(four_ilons, four_ilats, ens_size, depth_at_x)) then
   istatus(:) = QUAD_ON_BASIN_EDGE
   return
endif


select case (qty_in)
   case (QTY_TEMPERATURE)
      ! convert from potential temperature to temperature

      call state_on_quad(four_ilons, four_ilats, lon_lat_vert, ens_size, lev, lev_fract, interp, state_handle, varid, expected_pot_temp, quad_status)
      if (quad_status /= 0) then
         istatus(:) = QUAD_EVALUATE_FAILED
         return
      endif
      call state_on_quad(four_ilons, four_ilats, lon_lat_vert, ens_size, lev, lev_fract, interp, state_handle, get_varid_from_kind(dom_id, QTY_SALINITY), expected_salinity, quad_status)
      if (quad_status /= 0) then
         istatus(:) = QUAD_EVALUATE_FAILED
         return
      endif

      pressure_bars =  0.059808_r8*(exp(-0.025_r8*lon_lat_vert(3)) - 1.0_r8)  &
                        + 0.100766_r8*lon_lat_vert(3) + 2.28405e-7_r8*lon_lat_vert(3)**2
      expected_obs = sensible_temp(expected_pot_temp, expected_salinity, pressure_bars*10.0_r8)

   case (QTY_SALINITY) ! convert from g of salt per kg of seawater (model) to kg of salt per kg of seawater (observation)
      call state_on_quad(four_ilons, four_ilats, lon_lat_vert, ens_size, lev, lev_fract, interp, state_handle, varid, expected_obs, quad_status)
      if (quad_status /= 0) then 
         istatus(:) = QUAD_EVALUATE_FAILED
         return
      endif
      expected_obs = expected_obs/CONCENTRATION_TO_PPT
 
   case default
      call state_on_quad(four_ilons, four_ilats, lon_lat_vert, ens_size, lev, lev_fract, interp, state_handle, varid, expected_obs, quad_status)
      if (quad_status /= 0) then
         istatus(:) = QUAD_EVALUATE_FAILED
         return
      endif
end select

istatus(:) = 0

end subroutine model_interpolate

!------------------------------------------------------------------
! Interpolate on the quad, between two levels
subroutine state_on_quad(four_ilons, four_ilats, lon_lat_vert, ens_size, lev, lev_fract, interp, state_handle, varid, expected_obs, quad_status)

integer,  intent(in) :: four_ilons(4), four_ilats(4) ! indices into lon, lat
real(r8), intent(in) :: lon_lat_vert(3) ! lon, lat, vert of obs
integer,  intent(in) :: ens_size
integer,  intent(in) :: lev(ens_size,2) ! levels below and above obs
real(r8), intent(in) :: lev_fract(ens_size) 
type(quad_interp_handle), intent(in) :: interp
type(ensemble_type), intent(in) :: state_handle
integer,  intent(in) :: varid ! which state variable
real(r8), intent(out) :: expected_obs(ens_size)
integer,  intent(out) :: quad_status

integer :: i, e
integer(i8) :: indx(ens_size)
real(r8) :: quad_vals(4, ens_size)
real(r8) :: expected(ens_size, 2) ! state value at level below and above obs

do i = 1, 2 
   ! corner1
   do e = 1, ens_size
      indx(e) = get_dart_vector_index(four_ilons(1), four_ilats(1), lev(e, i), dom_id, varid)
   enddo
   call get_state_array(quad_vals(1, :), indx, state_handle)

   ! corner2
   do e = 1, ens_size
      indx(e) = get_dart_vector_index(four_ilons(2), four_ilats(2), lev(e, i), dom_id, varid)
   enddo
   call get_state_array(quad_vals(2, :), indx, state_handle)

   ! corner3
   do e = 1, ens_size
      indx(e) = get_dart_vector_index(four_ilons(3), four_ilats(3), lev(e, i), dom_id, varid)
   enddo
   call get_state_array(quad_vals(3, :), indx, state_handle)

   ! corner4
   do e = 1, ens_size
      indx(e) = get_dart_vector_index(four_ilons(4), four_ilats(4), lev(e, i), dom_id, varid)
   enddo
   call get_state_array(quad_vals(4, :), indx, state_handle)

   call quad_lon_lat_evaluate(interp, &
                              lon_lat_vert(1), lon_lat_vert(2), & ! lon, lat of obs
                              four_ilons, four_ilats, &
                              ens_size, &
                              quad_vals, & ! 4 corners x ens_size
                              expected(:,i), &
                              quad_status)
   if (quad_status /= 0) return

enddo

! Interpolate between levels
! expected_obs = bot_val + lev_fract * (top_val - bot_val)
expected_obs = expected(:,1) + lev_fract(:) * (expected(:,2) - expected(:,1))

end subroutine state_on_quad

!------------------------------------------------------------------
! Returns the smallest increment in time that the model is capable 
! of advancing the state in a given implementation, or the shortest
! time you want the model to advance between assimilations.

function shortest_time_between_assimilations()

type(time_type) :: shortest_time_between_assimilations

if ( .not. module_initialized ) call static_init_model

shortest_time_between_assimilations = assimilation_time_step

end function shortest_time_between_assimilations



!------------------------------------------------------------------
! Given an integer index into the state vector, returns the
! associated location and optionally the physical quantity.

subroutine get_state_meta_data(index_in, location, qty)

integer(i8),         intent(in)  :: index_in
type(location_type), intent(out) :: location
integer,             intent(out), optional :: qty

real(r8) :: lat, lon
integer :: lon_index, lat_index, level, local_qty

if ( .not. module_initialized ) call static_init_model

call get_model_variable_indices(index_in, lon_index, lat_index, level, kind_index=local_qty)

call get_lon_lat(lon_index, lat_index, local_qty, lon, lat)

location = set_location(lon, lat, real(level,r8), VERTISLEVEL)

if (present(qty)) then
   qty = local_qty
   if (on_land(lon_index, lat_index)) qty = QTY_DRY_LAND
endif


end subroutine get_state_meta_data

!------------------------------------------------------------------
subroutine convert_vertical_state(state_handle, num, locs, loc_qtys, loc_indx, &
which_vert, istatus)

type(ensemble_type), intent(in)    :: state_handle
integer,             intent(in)    :: num
type(location_type), intent(inout) :: locs(num) !locations
integer,             intent(in)    :: loc_qtys(num) !qty at location
integer(i8),         intent(in)    :: loc_indx(num) !state index
integer,             intent(in)    :: which_vert
integer,             intent(out)   :: istatus

integer :: i,j,k
integer :: ii, layer ! loop variables
integer :: thick_id
integer(i8) :: indx
real(r8) :: depth(1)

! assert(which_vert == VERTISHEIGHT)

if (use_pseudo_depth) then
   do ii = 1, num
      if (loc_qtys(ii) == QTY_DRY_LAND) then
         call set_vertical(locs(ii), 0.0_r8, VERTISHEIGHT)
      else
         call get_model_variable_indices(loc_indx(ii), i, j, k)
         depth(1) = zstar(k) ! pseudo depth
         call set_vertical(locs(ii), depth(1), VERTISHEIGHT)
      endif
   enddo
   istatus = 0
   return
endif

! If not using pseudo depth, then sum the layer thicknesses to get vertical location
thick_id = get_varid_from_kind(dom_id, QTY_LAYER_THICKNESS)
if (thick_id < 0) then
   istatus = THICKNESS_NOT_IN_STATE
   return
endif

do ii = 1, num

   if (loc_qtys(ii) == QTY_DRY_LAND) then
      call set_vertical(locs(ii), 0.0_r8, VERTISHEIGHT)
   else

      call get_model_variable_indices(loc_indx(ii), i, j, k)

      depth = 0.0_r8
      do layer = 1, k
         indx = get_dart_vector_index(i, j, layer, dom_id, thick_id)
         depth = depth + get_state(indx, state_handle)
      enddo

      call set_vertical(locs(ii), depth(1), VERTISHEIGHT)

   endif

enddo

istatus = 0

end subroutine convert_vertical_state
!------------------------------------------------------------------
! Any model specific distance calcualtion can be done here
subroutine get_close_obs(gc, base_loc, base_type, locs, loc_qtys, loc_types, &
                         num_close, close_ind, dist, ens_handle)

type(get_close_type),          intent(in)    :: gc            ! handle to a get_close structure
integer,                       intent(in)    :: base_type     ! observation TYPE
type(location_type),           intent(inout) :: base_loc      ! location of interest
type(location_type),           intent(inout) :: locs(:)       ! obs locations
integer,                       intent(in)    :: loc_qtys(:)   ! QTYS for obs
integer,                       intent(in)    :: loc_types(:)  ! TYPES for obs
integer,                       intent(out)   :: num_close     ! how many are close
integer,                       intent(out)   :: close_ind(:)  ! incidies into the locs array
real(r8),            optional, intent(out)   :: dist(:)       ! distances in radians
type(ensemble_type), optional, intent(in)    :: ens_handle

character(len=*), parameter :: routine = 'get_close_obs'

call loc_get_close_obs(gc, base_loc, base_type, locs, loc_qtys, loc_types, &
                          num_close, close_ind, dist, ens_handle)

end subroutine get_close_obs


!------------------------------------------------------------------
! Any model specific distance calcualtion can be done here
subroutine get_close_state(gc, base_loc, base_type, locs, loc_qtys, loc_indx, &
                           num_close, close_ind, dist, ens_handle)

type(get_close_type),          intent(in)    :: gc           ! handle to a get_close structure
type(location_type),           intent(inout) :: base_loc     ! location of interest
integer,                       intent(in)    :: base_type    ! observation TYPE
type(location_type),           intent(inout) :: locs(:)      ! state locations
integer,                       intent(in)    :: loc_qtys(:)  ! QTYs for state
integer(i8),                   intent(in)    :: loc_indx(:)  ! indices into DART state vector
integer,                       intent(out)   :: num_close    ! how many are close
integer,                       intent(out)   :: close_ind(:) ! indices into the locs array
real(r8),            optional, intent(out)   :: dist(:)      ! distances in radians
type(ensemble_type), optional, intent(in)    :: ens_handle

character(len=*), parameter :: routine = 'get_close_state'

integer :: ii ! loop index
integer :: i, j, k, istatus
real(r8) :: lon_lat_vert(3)
integer(i8) :: ind


call loc_get_close_state(gc, base_loc, base_type, locs, loc_qtys, loc_indx, &
                            num_close, close_ind)

if (.not. present(dist)) return

! Put any land or sea floor points very far away
! so they are not updated by assimilation
do ii = 1, num_close

  ind = close_ind(ii)
  if(loc_qtys(ind) == QTY_DRY_LAND) then 
    dist(ii) = 1.0e9_r8
    cycle
  endif

  lon_lat_vert = get_location(locs(ind))
  if (query_location(locs(ind)) /= VERTISHEIGHT) then ! assuming VERTISHEIGHT
    call convert_vertical_state(ens_handle, 1, locs(ind:ind), loc_qtys(ind:ind), loc_indx(ind:ind), VERTISHEIGHT, istatus)
  endif
  call get_model_variable_indices(loc_indx(ind), i, j, k)
  if ( below_sea_floor(i,j,lon_lat_vert(3)) ) then 
    dist(ii) = 1.0e9_r8
    cycle
  endif
  dist(ii) = get_dist(base_loc, locs(ind), base_type, loc_qtys(ind))

enddo

end subroutine get_close_state


!------------------------------------------------------------------
! Does any shutdown and clean-up needed for model. Can be a NULL
! INTERFACE if the model has no need to clean up storage, etc.

subroutine end_model()


end subroutine end_model


!------------------------------------------------------------------
! write any additional attributes to the output and diagnostic files

subroutine nc_write_model_atts(ncid, domain_id)

integer, intent(in) :: ncid      ! netCDF file identifier
integer, intent(in) :: domain_id

if ( .not. module_initialized ) call static_init_model

! put file into define mode.

call nc_begin_define_mode(ncid)

call nc_add_global_creation_time(ncid)

call nc_add_global_attribute(ncid, "model_source", source )
call nc_add_global_attribute(ncid, "model", "MOM6")

call nc_end_define_mode(ncid)

! Flush the buffer and leave netCDF file open
call nc_synchronize_file(ncid)

end subroutine nc_write_model_atts

!------------------------------------------------------------
! Read lon, lat for T,U,V grids from mom6 static file
subroutine read_horizontal_grid()

integer :: ncid
integer :: nxy(2) ! (nx,ny)
real(r8) :: fillval

character(len=*), parameter :: routine = 'read_horizontal_grid'

ncid = nc_open_file_readonly(static_file)

call nc_get_variable_size(ncid, 'geolon', nxy)
nx = nxy(1)
ny = nxy(2)
allocate(geolon(nx,ny), geolat(nx,ny))      ! T grid
allocate(geolon_u(nx,ny), geolat_u(nx,ny))  ! U grid
allocate(geolon_v(nx,ny), geolat_v(nx,ny))  ! V grid
allocate(mask(nx,ny))  ! missing values
allocate(mask_u(nx,ny))  ! missing values
allocate(mask_v(nx,ny))  ! missing values

call nc_get_variable(ncid, 'geolon', geolon, routine)
call nc_get_variable(ncid, 'geolon_u', geolon_u, routine)
call nc_get_variable(ncid, 'geolon_v', geolon_v, routine)

call nc_get_variable(ncid, 'geolat', geolat, routine)
call nc_get_variable(ncid, 'geolat_u', geolat_u, routine)
call nc_get_variable(ncid, 'geolat_v', geolat_v, routine)

! mom6 has missing values in the grid
! and set missing value to a land point to prevent set_location erroring
mask(:,:) = .false.
mask_u(:,:) = .false.
mask_v(:,:) = .false.
call nc_get_attribute_from_variable(ncid, 'geolon', '_FillValue', fillval)
where (geolon == fillval) mask = .true.  
where (geolon == fillval) geolon = 72.51_r8
where (geolat == fillval) geolat = 42.56_r8

call nc_get_attribute_from_variable(ncid, 'geolon_u', '_FillValue', fillval)
where (geolon_u == fillval) mask_u = .true.  
where (geolon_u == fillval) geolon_u = 72.51_r8
where (geolat_u == fillval) geolat_u = 42.56_r8

call nc_get_attribute_from_variable(ncid, 'geolon_v', '_FillValue', fillval)
where (geolon_v == fillval) mask_v = .true.  
where (geolon_v == fillval) geolon_v = 72.51_r8
where (geolat_v == fillval) geolat_v = 42.56_r8

! mom6 example files have longitude > 360 and longitudes < 0
! DART uses [0,360]
geolon = mod(geolon, 360.0_r8)
geolon_u = mod(geolon_u, 360.0_r8)
geolon_v = mod(geolon_v, 360.0_r8)

where (geolon < 0.0) geolon = geolon + 360.0_r8
where (geolon_u < 0.0) geolon_u = geolon_u + 360.0_r8
where (geolon_v < 0.0) geolon_v = geolon_v + 360.0_r8

call nc_close_file(ncid)

end subroutine read_horizontal_grid

!------------------------------------------------------------
! Read number of vertical layers from mom6 template file
subroutine read_num_layers()

integer :: ncid

character(len=*), parameter :: routine = 'read_num_layers'

ncid = nc_open_file_readonly(template_file)

call nc_get_variable_size(ncid, layer_name, nz)
if (use_pseudo_depth) then
   allocate(zstar(nz))
   call nc_get_variable(ncid, layer_name, zstar, routine)
endif

call nc_close_file(ncid)

end subroutine read_num_layers
!------------------------------------------------------------
! ocean_geom are 2D state sized static data
! HK Do these arrays become too big in high res cases?
subroutine read_ocean_geometry()

integer :: ncid

character(len=*), parameter :: routine = 'read_ocean_geometry'

! Need nx, ny
if ( .not. module_initialized ) call static_init_model

ncid = nc_open_file_readonly(ocean_geometry)

allocate(wet(nx,ny), basin_depth(nx,ny))
call nc_get_variable(ncid, 'wet', wet, routine)
call nc_get_variable(ncid, 'D', basin_depth, routine)

call nc_close_file(ncid)

end subroutine read_ocean_geometry

!------------------------------------------------------------
! wet is a 2D array of ones and zeros
! 1 is ocean
! 0 is land
function on_land_quad(ilon, ilat)

integer :: ilon(4), ilat(4) ! these are indices into lon, lat
logical ::  on_land_quad

if ( wet(ilon(1), ilat(1)) + &
     wet(ilon(2), ilat(2)) + &
     wet(ilon(3), ilat(3)) + &
     wet(ilon(4), ilat(4))  < 4) then
   on_land_quad = .true.
else
   on_land_quad = .false.
endif

end function on_land_quad

!------------------------------------------------------------
function on_land_point(ilon, ilat)

integer :: ilon, ilat ! these are indices into lon, lat
logical :: on_land_point

if ( wet(ilon, ilat) == 0) then
   on_land_point = .true.
else
   on_land_point = .false.
endif


end function on_land_point

!------------------------------------------------------------
! basin_depth is a 2D array with the basin depth
function on_basin_edge(ilon, ilat, ens_size, depth)

! indices into lon, lat lev
integer, intent(in)  :: ilon(4), ilat(4)
integer, intent(in)  :: ens_size
real(r8), intent(in) :: depth(ens_size)
logical :: on_basin_edge

integer  :: i, e
real(r8) :: d(4) ! basin depth at each corner

d(1) = basin_depth(ilon(1), ilat(1))
d(2) = basin_depth(ilon(2), ilat(2))
d(3) = basin_depth(ilon(3), ilat(3))
d(4) = basin_depth(ilon(4), ilat(4))

do e = 1, ens_size
   do i = 1, 4
      if (d(i) < depth(e)) then
         on_basin_edge = .true.
         return
      endif
   enddo
enddo

! four points are in the ocean
on_basin_edge = .false.

end function on_basin_edge

!------------------------------------------------------------
function below_sea_floor(ilon, ilat, depth)

! indices into lon, lat lev
integer,  intent(in) :: ilon, ilat
real(r8), intent(in) :: depth
logical :: below_sea_floor

if (basin_depth(ilon, ilat) < depth) then
   below_sea_floor = .true.
else
   below_sea_floor = .false.
endif

end function below_sea_floor

!------------------------------------------------------------
! longitude and latitide values from indices
subroutine get_lon_lat(lon_indx, lat_indx, qty, lon, lat)

integer, intent(in) :: lon_indx, lat_indx, qty
real(r8) :: lon, lat

if (on_u_grid(qty)) then
   lon = geolon_u(lon_indx, lat_indx)
   lat = geolat_u(lon_indx, lat_indx)
elseif (on_v_grid(qty)) then
   lon = geolon_v(lon_indx, lat_indx)
   lat = geolat_v(lon_indx, lat_indx)
else ! T grid
   lon = geolon(lon_indx, lat_indx)
   lat = geolat(lon_indx, lat_indx)
endif

end subroutine get_lon_lat

!------------------------------------------------------------
function on_v_grid(qty)

integer, intent(in)  :: qty
logical :: on_v_grid

if (qty == QTY_V_CURRENT_COMPONENT) then
  on_v_grid = .true.
else
  on_v_grid = .false.
endif

end function on_v_grid

!----------------------------------------------------------
function on_u_grid(qty)

integer, intent(in)  :: qty
logical :: on_u_grid

if (qty == QTY_U_CURRENT_COMPONENT) then
  on_u_grid = .true.
else
  on_u_grid = .false.
endif

end function on_u_grid

!----------------------------------------------------------
function on_t_grid(qty)

integer, intent(in)  :: qty
logical :: on_t_grid

if (qty == QTY_U_CURRENT_COMPONENT .or. qty == QTY_V_CURRENT_COMPONENT) then
  on_t_grid = .false.
else
  on_t_grid = .true.
endif

end function on_t_grid

!------------------------------------------------------------
function on_layer(qty)

logical :: on_layer
integer :: qty

! Salt, Temp, u, v all on layer 
on_layer = .true.

end function on_layer

!------------------------------------------------------------
subroutine setup_interpolation()

! T
call init_quad_interp(GRID_QUAD_FULLY_IRREGULAR, nx, ny, &
                      QUAD_LOCATED_CELL_CENTERS, &
                      global=.true., spans_lon_zero=.true., pole_wrap=.true., &
                      interp_handle=interp_t_grid)

call set_quad_coords(interp_t_grid, geolon, geolat, mask)

! U
call init_quad_interp(GRID_QUAD_FULLY_IRREGULAR, nx, ny, &
                      QUAD_LOCATED_CELL_CENTERS, &
                      global=.true., spans_lon_zero=.true., pole_wrap=.true., &
                      interp_handle=interp_u_grid)

call set_quad_coords(interp_u_grid, geolon_u, geolat_u, mask_u)


! V
call init_quad_interp(GRID_QUAD_FULLY_IRREGULAR, nx, ny, &
                      QUAD_LOCATED_CELL_CENTERS, &
                      global=.true., spans_lon_zero=.true., pole_wrap=.true., &
                      interp_handle=interp_v_grid)

call set_quad_coords(interp_v_grid, geolon_v, geolat_v, mask_v)

end subroutine setup_interpolation

!------------------------------------------------------------
! return the appropriate quad_interp handle
function get_interp_handle(qty)

type(quad_interp_handle) :: get_interp_handle
integer, intent(in) :: qty

if (on_v_grid(qty)) then
  get_interp_handle = interp_v_grid
elseif (on_u_grid(qty)) then
  get_interp_handle = interp_u_grid
else
  get_interp_handle = interp_t_grid
endif

end function

!------------------------------------------------------------
! calculate sensible (in-situ) temperature from 
! local pressure, salinity, and potential temperature
elemental function sensible_temp(potemp, s, lpres)

real(r8), intent(in)  :: potemp ! potential temperature in C
real(r8), intent(in)  :: s      ! salinity Practical Salinity Scale 1978 (PSS-78)
real(r8), intent(in)  :: lpres  ! pressure in decibars
real(r8) :: sensible_temp ! in-situ (sensible) temperature (C)

integer  :: i,j,n
real(r8) :: dp,p,q,r1,r2,r3,r4,r5,s1,t,x

s1 = s - 35.0_r8
p  = 0.0_r8
t  = potemp

dp = lpres - p
n  = int (abs(dp)/1000.0_r8) + 1
dp = dp/n

do i=1,n
   do j=1,4

      r1 = ((-2.1687e-16_r8 * t + 1.8676e-14_r8) * t - 4.6206e-13_r8) * p
      r2 = (2.7759e-12_r8*t - 1.1351e-10_r8) * s1
      r3 = ((-5.4481e-14_r8 * t + 8.733e-12_r8) * t - 6.7795e-10_r8) * t
      r4 = (r1 + (r2 + r3 + 1.8741e-8_r8)) * p + (-4.2393e-8_r8 * t+1.8932e-6_r8) * s1
      r5 = r4 + ((6.6228e-10_r8 * t-6.836e-8_r8) * t + 8.5258e-6_r8) * t + 3.5803e-5_r8

      x  = dp*r5

      if (j == 1) then
         t = t + 0.5_r8 * x
         q = x
         p = p + 0.5_r8 * dp
      
      else if (j == 2) then
         t = t + 0.29298322_r8 * (x-q)
         q = 0.58578644_r8 * x + 0.121320344_r8 * q

      else if (j == 3) then
         t = t + 1.707106781_r8 * (x-q)
         q = 3.414213562_r8*x - 4.121320344_r8*q
         p = p + 0.5_r8*dp

      else ! j must == 4
         t = t + (x - 2.0_r8 * q) / 6.0_r8

      endif

   enddo ! j loop
enddo ! i loop

sensible_temp = t

end function sensible_temp

!--------------------------------------------------------------------
function read_model_time(filename)

character(len=*), intent(in) :: filename
type(time_type) :: read_model_time

integer :: ncid
character(len=*), parameter :: routine = 'read_model_time'
real(r8) :: days
type(time_type) :: mom6_time
integer :: dart_base_date_in_days, dart_days

dart_base_date_in_days = 584388 ! 1601 1 1 0 0
ncid = nc_open_file_readonly(filename, routine)

call nc_get_variable(ncid, 'Time', days, routine)

call nc_close_file(ncid, routine)

! MOM6 counts days from year 1
! DART counts days from 1601 
dart_days = int(days) - dart_base_date_in_days

read_model_time = set_time(0,dart_days)

end function read_model_time

!--------------------------------------------------------------------
!                  Surface
!                    --- 0 i=1  Interface pseudo-depth
! Layer pseudo-depth  .
!                    --- 1 i=2
!                     .
!                    --- 1 i=3
!

subroutine find_level_bounds(vert_loc, which_vert, lev, lev_fract, istatus)

real(r8), intent(in) :: vert_loc   ! observation location
integer,  intent(in) :: which_vert ! obs vertical coordinate
integer,  intent(out) :: lev(2)    ! bottom, top
real(r8), intent(out) :: lev_fract
integer,  intent(out) :: istatus

integer :: i

! HK assert(which_vert == height meters) ?

do i = 2, nz
   if(vert_loc < zstar(i)) then
      lev(1) = i   ! bottom
      lev(2) = i-1 ! top

      lev_fract = (zstar(lev(1)) - vert_loc) / (zstar(lev(2)) - zstar(lev(1)))
      istatus = 0
      return
   endif
enddo

istatus = OBS_TOO_DEEP

end subroutine find_level_bounds

!===================================================================
! End of model_mod
!===================================================================
end module model_mod

