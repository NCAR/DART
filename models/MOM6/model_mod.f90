! DART software - Copyright UCAR. This open source software is provided
! by UCAR, "as is", without charge, subject to all terms of use at
! http://www.image.ucar.edu/DAReS/DART/DART_download
!

! @todo
!   layer thickness - this is per ensemble member?
!    

module model_mod

use        types_mod, only : r8, i8, MISSING_R8, vtablenamelength

use time_manager_mod, only : time_type, set_time, set_date, set_calendar_type

use     location_mod, only : location_type, get_close_type, &
                             loc_get_close_obs => get_close_obs, &
                             loc_get_close_state => get_close_state, &
                             set_location, set_location_missing, &
                             get_location, query_location, VERTISHEIGHT

use    utilities_mod, only : error_handler, &
                             E_ERR, E_MSG, &
                             nmlfileunit, do_output, do_nml_file, do_nml_term,  &
                             find_namelist_in_file, check_namelist_read, &
                             to_upper

use netcdf_utilities_mod, only : nc_add_global_attribute, nc_synchronize_file, &
                                 nc_add_global_creation_time, &
                                 nc_begin_define_mode, nc_end_define_mode, &
                                 nc_open_file_readonly, nc_close_file, &
                                 nc_get_variable, nc_get_variable_size, &
                                 NF90_MAX_NAME

use        quad_utils_mod,  only : quad_interp_handle, init_quad_interp, &
                                   set_quad_coords, quad_lon_lat_locate, &
                                   quad_lon_lat_evaluate, quad_interp_handle, &
                                   GRID_QUAD_IRREG_SPACED_REGULAR, &
                                   QUAD_LOCATED_CELL_CENTERS

use state_structure_mod, only : add_domain, get_domain_size, &
                                get_model_variable_indices, &
                                get_kind_string, get_varid_from_kind, &
                                get_dart_vector_index

use distributed_state_mod, only : get_state

use obs_kind_mod, only : get_index_for_quantity, QTY_U_CURRENT_COMPONENT, &
                         QTY_V_CURRENT_COMPONENT

use ensemble_manager_mod, only : ensemble_type

! These routines are passed through from default_model_mod.
! To write model specific versions of these routines
! remove the routine from this use statement and add your code to
! this the file.
use default_model_mod, only : pert_model_copies, write_model_time, &
                              init_time => fail_init_time, &
                              init_conditions => fail_init_conditions, &
                              convert_vertical_obs, convert_vertical_state, adv_1step

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
! Grid parameters
! nx, ny and nz are the size of the dipole (or irregular) grids.
integer :: nx=-1, ny=-1, nz=-1     ! grid counts for each field
!netcdf bmom.e20.BMOM.f09_t061.long_run_mct.006.mom6.r.0006-01-01-00000 {
!dimensions:
!        lath = 458 ;
!        lonh = 540 ;
!        latq = 458 ;
!        lonq = 540 ;
!        Layer = 63 ;
!        Interface = 64 ;
!        Time = UNLIMITED ; // (1 currently)
real(r8), allocatable :: lath(:), lonh(:), latq(:), lonq(:), layer(:), interf(:)
real(r8), allocatable :: wet(:,:), basin_depth(:,:)
type(quad_interp_handle) :: interp_v_grid, interp_u_grid, interp_t_grid !HK do we need all three?

! DART state vector contents are specified in the input.nml:&model_nml namelist.
integer, parameter :: MAX_STATE_VARIABLES = 10
integer, parameter :: NUM_STATE_TABLE_COLUMNS = 3

! model_interpolate failure codes
integer, parameter :: NOT_IN_STATE = 12
integer, parameter :: QUAD_LOCATE_FAILED = 14
integer, parameter :: QUAD_ON_LAND = 16
integer, parameter :: QUAD_ON_BASIN_EDGE = 18
integer, parameter :: OBS_ABOVE_SURFACE = 20
integer, parameter :: OBS_TOO_DEEP = 22


! namelist
character(len=256) :: template_file = 'mom6.r.nc'
character(len=256) :: ocean_geometry = 'ocean_geometry.nc'
integer  :: assimilation_period_days      = -1
integer  :: assimilation_period_seconds   = -1
character(len=vtablenamelength) :: model_state_variables(MAX_STATE_VARIABLES * NUM_STATE_TABLE_COLUMNS ) = ' '

namelist /model_nml/ template_file, ocean_geometry, assimilation_period_days, &
                     assimilation_period_seconds, model_state_variables

contains

!------------------------------------------------------------------
!
! Called to do one time initialization of the model. As examples,
! might define information about the model size or model timestep.
! In models that require pre-computed static data, for instance
! spherical harmonic weights, these would also be computed here.

subroutine static_init_model()

integer  :: iunit, io
character(len=vtablenamelength) :: variable_table(MAX_STATE_VARIABLES, NUM_STATE_TABLE_COLUMNS)

integer :: state_qty_list(MAX_STATE_VARIABLES)
logical :: update_var_list(MAX_STATE_VARIABLES)

! identifiers for variable_table
integer, parameter :: VAR_NAME_INDEX = 1
integer, parameter :: VAR_QTY_INDEX = 2
integer, parameter :: VAR_UPDATE_INDEX = 3

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

! verify that the model_state_variables namelist was filled in correctly.
! returns variable_table which has variable names, kinds and update strings.
call verify_state_variables(model_state_variables, nfields, variable_table, state_qty_list, update_var_list)

! Define which variables are in the model state
dom_id = add_domain(template_file, nfields, &
                    var_names = variable_table(1:nfields, VAR_NAME_INDEX), &
                    kind_list = state_qty_list(1:nfields), &
                    update_list = update_var_list(1:nfields))

model_size = get_domain_size(dom_id)

call read_grid()

call setup_interpolation()

call read_ocean_geometry()

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

subroutine model_interpolate(state_handle, ens_size, location, qty, expected_obs, istatus)


type(ensemble_type), intent(in) :: state_handle
integer,             intent(in) :: ens_size
type(location_type), intent(in) :: location
integer,             intent(in) :: qty
real(r8),           intent(out) :: expected_obs(ens_size) !< array of interpolated values
integer,            intent(out) :: istatus(ens_size)

integer  :: which_vert, four_ilons(4), four_ilats(4)
integer  :: locate_status
real(r8) :: lon_fract, lat_fract, lev_fract
real(r8) :: lon_lat_vert(3)
real(r8) :: quad_vals(4, ens_size)
real(r8) :: expected(ens_size, 2) ! level below and above obs
type(quad_interp_handle) :: interp
integer :: varid, i
integer(i8) :: indx
integer :: lev(2) ! level bounds

if ( .not. module_initialized ) call static_init_model

expected_obs(:) = MISSING_R8
istatus(:) = 1

varid = get_varid_from_kind(dom_id, qty)
if (varid < 0) then ! not in state
   istatus = NOT_IN_STATE
   return
endif

! find which grid the qty is on
interp = get_interp_handle(qty)

! unpack the location type into lon, lat, vert, vert_type
lon_lat_vert = get_location(location)
which_vert   = nint(query_location(location))

! get the indices for the 4 corners of the quad in the horizontal, plus
! the fraction across the quad for the obs location
call quad_lon_lat_locate(interp, lon_lat_vert(1), lon_lat_vert(2), &
                         four_ilons, four_ilats, lon_fract, lat_fract, locate_status)

if (locate_status /= 0) then
  istatus(:) = QUAD_LOCATE_FAILED
  return
endif

! check if all four corners are in the ocean
if (on_land(four_ilons, four_ilats)) then
   istatus(:) = QUAD_ON_LAND
   return
endif

! find levels
! Get the bounding vertical levels and the fraction between bottom and top
call find_level_bounds(lon_lat_vert(3), which_vert, lev, lev_fract, locate_status)
if (locate_status /= 0) then
  istatus(:) = locate_status
  return
endif

if (on_basin_edge(four_ilons, four_ilats, lev)) then
   istatus(:) = QUAD_ON_BASIN_EDGE
   return
endif

! get values of state at four corners

do i = 1, 2
   !HK which corner of the quad is which?
   ! corner1
   indx = get_dart_vector_index(four_ilons(1), four_ilats(1), lev(i), dom_id, varid)
   quad_vals(1, :) = get_state(indx, state_handle)

   ! corner2
   indx = get_dart_vector_index(four_ilons(1), four_ilats(2), lev(i), dom_id, varid)
   quad_vals(2, :) = get_state(indx, state_handle)

   ! corner3
   indx = get_dart_vector_index(four_ilons(2), four_ilats(1), lev(i), dom_id, varid)
   quad_vals(3, :) = get_state(indx, state_handle)

   ! corner4
   indx = get_dart_vector_index(four_ilons(2), four_ilats(2), lev(i), dom_id, varid)
   quad_vals(4, :) = get_state(indx, state_handle)

   call quad_lon_lat_evaluate(interp, &
                              lon_fract, lat_fract, &
                              ens_size, &
                              quad_vals, & ! 4 corners x ens_size
                              expected(:,i), &
                              istatus)
enddo

! Interpolate between levels
! expected_obs = bot_val + lev_fract * (top_val - bot_val)
expected_obs = expected(:,1) + lev_fract * (expected(:,2) - expected(:,1))

end subroutine model_interpolate



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

real(r8) :: lat, lon, depth
integer :: lon_index, lat_index, depth_index, local_qty

if ( .not. module_initialized ) call static_init_model

call get_model_variable_indices(index_in, lon_index, lat_index, depth_index, kind_index=local_qty)

lon = get_lon(lon_index, local_qty)
lat = get_lat(lat_index, local_qty)
depth = get_depth(depth_index, local_qty)

location = set_location(lon, lat, depth, VERTISHEIGHT)

if (present(qty)) qty = local_qty

end subroutine get_state_meta_data


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


call loc_get_close_state(gc, base_loc, base_type, locs, loc_qtys, loc_indx, &
                            num_close, close_ind, dist, ens_handle)


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
! Read lon, lat, vertical from mom6 template file
subroutine read_grid()

integer :: ncid

character(len=*), parameter :: routine = 'read_grid'

ncid = nc_open_file_readonly(template_file)

call nc_get_variable_size(ncid, 'lonh', nx)
allocate(lonh(nx), lonq(nx))
call nc_get_variable(ncid, 'lonh', lonh, routine)
call nc_get_variable(ncid, 'lonq', lonq, routine)
! mom6 example file has -ve longitude
! DART uses [0,360]
where(lonh < 0 ) lonh = lonh + 360.0_r8
where(lonq < 0 ) lonq = lonq + 360.0_r8


call nc_get_variable_size(ncid, 'lath', ny)
allocate(lath(ny), latq(ny))
call nc_get_variable(ncid, 'lath', lath, routine)
call nc_get_variable(ncid, 'latq', latq, routine)

call nc_get_variable_size(ncid, 'Layer', nz)
allocate(layer(nz), interf(nz+1))
call nc_get_variable(ncid, 'Layer', layer, routine) ! Layer pseudo-depth
call nc_get_variable(ncid, 'Interface', interf, routine) ! Interface pseudo-depth

call nc_close_file(ncid)

end subroutine read_grid

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
function on_land(ilon, ilat)

integer :: ilon(4), ilat(4) ! these are indices into lon, lat
logical ::  on_land

if ( wet(ilon(1), ilat(1)) + &
     wet(ilon(1), ilat(2)) + &
     wet(ilon(2), ilat(1)) + &
     wet(ilon(2), ilat(2))  < 4) then
   on_land = .true.
else
   on_land = .false.
endif

end function on_land

!------------------------------------------------------------
! basin_depth is a 2D array with the basin depth
function on_basin_edge(ilon, ilat, lev)

! indices into lon, lat, lev
integer, intent(in) :: ilon(4), ilat(4), lev(2)
logical :: on_basin_edge

integer  :: i
real(r8) :: d(4) ! basin depth at each corner

d(1) = basin_depth(ilon(1), ilat(1))
d(2) = basin_depth(ilon(1), ilat(2))
d(3) = basin_depth(ilon(2), ilat(1))
d(4) = basin_depth(ilon(2), ilat(2))

do i = 1, 4
  if (d(i) < Layer(lev(1)) .or. d(i) < Layer(lev(2))) then
     on_basin_edge = .true.
     return
  endif
enddo

! four points are in the ocean
on_basin_edge = .false.

end function on_basin_edge

!------------------------------------------------------------
! longitude value from index
function get_lon(indx, qty)

integer, intent(in) :: indx, qty
real(r8) :: get_lon

if (on_u_grid(qty)) then
   get_lon = lonq(indx)  ! only u current on lonq
else
   get_lon = lonh(indx)
endif

end function get_lon

!------------------------------------------------------------
! latitude value from index
function get_lat(indx, qty)

integer, intent(in) :: indx, qty
real(r8) :: get_lat

if (on_v_grid(qty)) then ! only v current on latq
   get_lat = latq(indx)
else !on q grid
   get_lat = lath(indx)
endif

end function get_lat

!------------------------------------------------------------
! depth value from index
function get_depth(indx, qty)

integer, intent(in) :: indx, qty
real(r8) :: get_depth

if (on_layer(qty)) then
   get_depth = layer(indx)
else
   get_depth = interf(indx)
endif

end function get_depth

!------------------------------------------------------------
! Salt h
! Temp h
! u lath, lonq
! v latq, lonh
!----------------------------------------------------------
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

if (qty == QTY_U_CURRENT_COMPONENT .or. QTY_V_CURRENT_COMPONENT) then
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
! 1D arrays vs 2D arrays
! HK what grids do we have to deal with?
subroutine setup_interpolation()


call init_quad_interp(GRID_QUAD_IRREG_SPACED_REGULAR, nx, ny, &
                      QUAD_LOCATED_CELL_CENTERS, &
                      global=.true., spans_lon_zero=.true., pole_wrap=.true., &
                      interp_handle=interp_v_grid)

call set_quad_coords(interp_v_grid, lonq, lath)


call init_quad_interp(GRID_QUAD_IRREG_SPACED_REGULAR, nx, ny, &
                      QUAD_LOCATED_CELL_CENTERS, &
                      global=.true., spans_lon_zero=.true., pole_wrap=.true., &
                      interp_handle=interp_u_grid)

call set_quad_coords(interp_u_grid, lonh, latq)

call init_quad_interp(GRID_QUAD_IRREG_SPACED_REGULAR, nx, ny, &
                      QUAD_LOCATED_CELL_CENTERS, &
                      global=.true., spans_lon_zero=.true., pole_wrap=.true., &
                      interp_handle=interp_t_grid)

call set_quad_coords(interp_t_grid, lonh, lath)


end subroutine setup_interpolation

!------------------------------------------------------------
! return the appropriate quad_interp handle
function get_interp_handle(qty)

type(quad_interp_handle) :: get_interp_handle
integer, intent(in) :: qty

if (on_v_grid(qty)) then
  get_interp_handle = interp_v_grid
elseif (on_v_grid(qty)) then
  get_interp_handle = interp_u_grid
else
  get_interp_handle = interp_t_grid
endif

end function


!------------------------------------------------------------------
! Verify that the namelist was filled in correctly, and check
! that there are valid entries for the dart_kind.
! Returns a table with columns:
!
! netcdf_variable_name ; dart_qty_string ; update_string

subroutine verify_state_variables(state_variables, ngood, table, qty_list, update_var)

character(len=*),  intent(inout) :: state_variables(:)
integer,           intent(out) :: ngood
character(len=*),  intent(out) :: table(:,:)
integer,           intent(out) :: qty_list(:)   ! kind number
logical,           intent(out) :: update_var(:) ! logical update

integer :: nrows, i
character(len=NF90_MAX_NAME) :: varname, dartstr, update
character(len=256) :: string1, string2

if ( .not. module_initialized ) call static_init_model

nrows = size(table,1)

ngood = 0

MyLoop : do i = 1, nrows

   varname = trim(state_variables(3*i -2))
   dartstr = trim(state_variables(3*i -1))
   update  = trim(state_variables(3*i   ))
   
   call to_upper(update)

   table(i,1) = trim(varname)
   table(i,2) = trim(dartstr)
   table(i,3) = trim(update)

   if ( table(i,1) == ' ' .and. table(i,2) == ' ' .and. table(i,3) == ' ') exit MyLoop ! Found end of list.

   if ( table(i,1) == ' ' .or. table(i,2) == ' ' .or. table(i,3) == ' ' ) then
      string1 = 'model_nml:model_state_variables not fully specified'
      call error_handler(E_ERR,'verify_state_variables',string1)
   endif

   ! Make sure DART qty is valid

   qty_list(i) = get_index_for_quantity(dartstr)
   if( qty_list(i)  < 0 ) then
      write(string1,'(''there is no obs_kind <'',a,''> in obs_kind_mod.f90'')') trim(dartstr)
      call error_handler(E_ERR,'verify_state_variables',string1)
   endif
   
   ! Make sure the update variable has a valid name

   select case (update)
      case ('UPDATE')
         update_var(i) = .true.
      case ('NO_COPY_BACK')
         update_var(i) = .false.
      case default
         write(string1,'(A)')  'only UPDATE or NO_COPY_BACK supported in model_state_variable namelist'
         write(string2,'(6A)') 'you provided : ', trim(varname), ', ', trim(dartstr), ', ', trim(update)
         call error_handler(E_ERR,'verify_state_variables',string1, text2=string2)
   end select

   ngood = ngood + 1
enddo MyLoop


end subroutine verify_state_variables

!--------------------------------------------------------------------
function read_model_time(filename)

character(len=*), intent(in) :: filename
type(time_type) :: read_model_time

integer :: ncid
character(len=*), parameter :: routine = 'read_model_time'
real(r8) :: days

ncid = nc_open_file_readonly(filename, routine)

call nc_get_variable(ncid, 'Time', days, routine)

call nc_close_file(ncid, routine)

read_model_time = set_time(0,int(days))

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

if (vert_loc < interf(1)) then
  istatus = OBS_ABOVE_SURFACE
  return
endif

do i = 2, nz+1
   if(vert_loc < interf(i)) then
      lev(1) = i   ! bottom
      lev(2) = i-1 ! top

      lev_fract = (interf(lev(1)) - vert_loc) / (interf(lev(2)) - interf(lev(1)))
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

