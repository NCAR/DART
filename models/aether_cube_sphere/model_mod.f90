! DART software - Copyright UCAR. This open source software is provided
! by UCAR, "as is", without charge, subject to all terms of use at
! http://www.image.ucar.edu/DAReS/DART/DART_download
!

module model_mod

! This is a template showing the interfaces required for a model to be compliant
! with the DART data assimilation infrastructure. Do not change the arguments
! for the public routines.

use           netcdf

use        types_mod, only : r8, i8, MISSING_R8, vtablenamelength, DEG2RAD, RAD2DEG, radius => earth_radius, PI

use time_manager_mod, only : time_type, set_time

use     location_mod, only : location_type, get_close_type, &
                             loc_get_close_obs => get_close_obs, &
                             loc_get_close_state => get_close_state, &
                             set_location, set_location_missing, VERTISHEIGHT, query_location, &
                             get_location

use    utilities_mod, only : register_module, error_handler, &
                             E_ERR, E_MSG, &
                             nmlfileunit, do_output, do_nml_file, do_nml_term,  &
                             find_namelist_in_file, check_namelist_read, to_upper, &
                             find_enclosing_indices

use obs_kind_mod,     only : get_index_for_quantity

use netcdf_utilities_mod, only : nc_add_global_attribute, nc_synchronize_file, &
                                 nc_add_global_creation_time, &
                                 nc_begin_define_mode, nc_end_define_mode, nc_open_file_readonly, & 
                                 nc_close_file
use distributed_state_mod, only : get_state
use state_structure_mod, only : add_domain, get_dart_vector_index, get_domain_size, &
                                get_model_variable_indices, get_varid_from_kind

use transform_state_mod, only : file_type

use ensemble_manager_mod, only : ensemble_type

! These routines are passed through from default_model_mod.
! To write model specific versions of these routines
! remove the routine from this use statement and add your code to
! this the file.
use default_model_mod, only : pert_model_copies, read_model_time, write_model_time, &
                              init_time => fail_init_time, &
                              init_conditions => fail_init_conditions, &
                              convert_vertical_obs, convert_vertical_state, adv_1step

use quad_utils_mod,    only : in_quad, quad_bilinear_interp

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

! Routine for comprehensive test of interpolation
public :: test_grid_box

integer  :: iunit, io

character(len=256), parameter :: source   = "model_mod.f90"
logical :: module_initialized = .false.
integer :: dom_id ! used to access the state structure
type(time_type) :: assimilation_time_step

real(r8), parameter :: roundoff = 1.0e-12_r8

! Geometry variables that are used throughout the module
integer                          :: ncenter_columns, ncenter_altitudes ! The number of center_columns and altitudes are read from the geometry file
integer, parameter               :: nvertex_columns = 8
integer, parameter               :: nvertex_neighbors = 3
integer                          :: nquad_columns ! The number of quad_columns is read from the geometry file
integer, parameter               :: nquad_neighbors = 4

! Just like in cam-se, the aether cube_sphere filter input files are created to have a horizonal
! column dimension rather than being functions of latitude and longitude.
integer                          :: no_third_dimension = -99

integer :: inorth_pole_quad_column, isouth_pole_quad_column
real(r8), allocatable, dimension(:)       :: center_latitude, center_longitude, center_altitude
integer, allocatable, dimension(:, :)     :: quad_neighbor_indices, vertex_neighbor_indices

! Error codes
integer, parameter :: INVALID_VERT_COORD_ERROR_CODE = 15
integer, parameter :: INVALID_ALTITUDE_VAL_ERROR_CODE = 17
integer, parameter :: UNKNOWN_OBS_QTY_ERROR_CODE = 20

! Example Namelist
! Use the namelist for options to be set at runtime.
character(len=256) :: template_file = 'model_restart.nc'
integer  :: time_step_days      = 0
integer  :: time_step_seconds   = 3600

integer, parameter              :: MAX_STATE_VARIABLES     = 100
integer, parameter              :: NUM_STATE_TABLE_COLUMNS = 5
character(len=vtablenamelength) :: variables(NUM_STATE_TABLE_COLUMNS,MAX_STATE_VARIABLES) = ''

type :: var_type
    integer :: count
    character(len=64), allocatable :: names(:)
    integer,           allocatable :: qtys(:)
    real(r8),          allocatable :: clamp_values(:, :)
    logical,           allocatable :: updates(:)
end type var_type

namelist /model_nml/ template_file, time_step_days, time_step_seconds, variables

contains

!------------------------------------------------------------------
!
! Called to do one time initialization of the model. As examples,
! might define information about the model size or model timestep.
! In models that require pre-computed static data, for instance
! spherical harmonic weights, these would also be computed here.

subroutine static_init_model()

type(var_type) :: var

module_initialized = .true.

! Print module information to log file and stdout.
call register_module(source)

call find_namelist_in_file('input.nml', 'model_nml', iunit)
read(iunit, nml = model_nml, iostat = io)
call check_namelist_read(iunit, io, 'model_nml')

! Record the namelist values used for the run 
if (do_nml_file()) write(nmlfileunit, nml=model_nml)
if (do_nml_term()) write(     *     , nml=model_nml)

! This time is both the minimum time you can ask the model to advance
! (for models that can be advanced by filter) and it sets the assimilation
! window.  All observations within +/- 1/2 this interval from the current
! model time will be assimilated. If this is not settable at runtime 
! feel free to hardcode it and remove from the namelist.
assimilation_time_step = set_time(time_step_seconds, &
                                  time_step_days)

var = assign_var(variables, MAX_STATE_VARIABLES)

call read_geometry_file()

! Define which variables are in the model state
dom_id = add_domain(template_file, var%count, var%names, var%qtys, &
                    var%clamp_values, var%updates)

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
!

subroutine model_interpolate(state_handle, ens_size, location, qty, expected_obs, istatus)


type(ensemble_type), intent(in) :: state_handle
integer,             intent(in) :: ens_size
type(location_type), intent(in) :: location
integer,             intent(in) :: qty
real(r8),           intent(out) :: expected_obs(ens_size) !< array of interpolated values
integer,            intent(out) :: istatus(ens_size)

character(len=512) :: error_string_1

! Location values stored in a vector
real(r8), dimension(3)                           :: lon_lat_alt

! Logical needed for quad_bilinear_interp
logical                                          :: cyclic

! Vertical interpolation variables
integer(i8) :: state_index
integer     :: below_index, above_index, enclosing_status, which_vertical
real(r8)    :: fraction

real(r8), dimension(nvertex_neighbors, ens_size) :: vertex_temp_values
real(r8), dimension(nquad_neighbors, ens_size)   :: quad_temp_values
real(r8), dimension(ens_size)                    :: above_values, below_values

! Actual variables needed for this routine

integer :: icolumn, ineighbor, iens

logical :: inside

real(r8), dimension(nvertex_neighbors, 3)        :: vertex_xyz
real(r8), dimension(3)                           :: r
real(r8), dimension(3)                           :: weights
real(r8), dimension(nquad_neighbors)             :: x_neighbors_quad, y_neighbors_quad

! Begin local variables for model_interpolate
integer  :: varid

! Initialize module if not already done
if ( .not. module_initialized ) call static_init_model

! Set all obs to MISSING_R8 initially
expected_obs(:) = MISSING_R8
! JLA: Is this successful default
istatus(:) = 0

! End variables for horizontal interpolation
cyclic = .true.

! Determine the vertical location type
lon_lat_alt = get_location(location)
which_vertical = nint(query_location(location, 'WHICH_VERT'))

! Only height currently supported for observations
if (.not. which_vertical == VERTISHEIGHT ) then
   istatus(:) = INVALID_VERT_COORD_ERROR_CODE
   !!!write(error_string_1, *) 'unsupported vertical type: ', which_vertical
   !!!call error_handler(E_ERR, 'model_interpolate', error_string_1, source)
   return
endif

! See if the state contains the obs quantity
varid = get_varid_from_kind(dom_id, qty)
if(varid <= 0) then
   istatus = UNKNOWN_OBS_QTY_ERROR_CODE
   return
endif

! Find the bounding vertical levels and the fractional distance between
call find_enclosing_indices(ncenter_altitudes, center_altitude, lon_lat_alt(3), below_index, &
   above_index, fraction, enclosing_status)
if (enclosing_status /= 0) then
   istatus(:) = INVALID_ALTITUDE_VAL_ERROR_CODE
   return
endif

! If the vertical location is acceptable, then do the horizontal interpolation
! Find the enclosing triangle or quad
!!!call get_bounding_box(pt_lat, pt_lon, np, &
   !!!grid_face, grid_lat_ind, grid_lon_ind, grid_pt_lat, grid_pt_lon, num_bound_points)
! Map the grid_face, latitude index and longitude index to the one dimensional index used in the state vector
!JLA

inside = .false.
! This is trying to search every single column
do icolumn = 1, nvertex_columns
   do ineighbor = 1, nvertex_neighbors
      call latlon_to_xyz(center_latitude(vertex_neighbor_indices(icolumn, ineighbor)), &
         center_longitude(vertex_neighbor_indices(icolumn, ineighbor)), &
         vertex_xyz(ineighbor, 1), vertex_xyz(ineighbor, 2), vertex_xyz(ineighbor, 3))
   end do

   call latlon_to_xyz(lon_lat_alt(2), lon_lat_alt(1), r(1), r(2), r(3))

   call inside_triangle(vertex_xyz(1, :), vertex_xyz(2, :), vertex_xyz(3, :), r, lon_lat_alt(2), lon_lat_alt(1), inside, weights)

   if (inside) then

      ! Do above level
      do ineighbor = 1, nvertex_neighbors
         state_index = get_dart_vector_index(vertex_neighbor_indices(icolumn, ineighbor), above_index, no_third_dimension, dom_id, varid)
         vertex_temp_values(ineighbor, :) = get_state(state_index, state_handle)
      end do

      above_values = barycentric_average(ens_size, weights, vertex_temp_values)

      ! Do below level
      do ineighbor = 1, nvertex_neighbors
         state_index = get_dart_vector_index(vertex_neighbor_indices(icolumn, ineighbor), below_index, no_third_dimension, dom_id, varid)
         vertex_temp_values(ineighbor, :) = get_state(state_index, state_handle)
      end do

      below_values = barycentric_average(ens_size, weights, vertex_temp_values)

      call vert_interp(ens_size, below_values, above_values, fraction, expected_obs)

      exit
   end if
end do

! If the location is not inside any of the vertex triangles, do the quad loop

if (.not. inside) then

   do icolumn = 1, nquad_columns
      if (icolumn == inorth_pole_quad_column) then
         inside = .true.
      else if  (icolumn == isouth_pole_quad_column) then
         inside = .true.
      else
         do ineighbor = 1, nquad_neighbors
            x_neighbors_quad(ineighbor) = center_longitude(quad_neighbor_indices(icolumn, ineighbor))
            y_neighbors_quad(ineighbor) = center_latitude(quad_neighbor_indices(icolumn, ineighbor))
         end do
         inside = in_quad(lon_lat_alt(1), lon_lat_alt(2), x_neighbors_quad, y_neighbors_quad)
      end if

      if (inside) then
            
         ! Get quad temp_values for the above level
         do ineighbor = 1, nquad_neighbors

            state_index = get_dart_vector_index(quad_neighbor_indices(icolumn, ineighbor), above_index, no_third_dimension, dom_id, varid)

            quad_temp_values(ineighbor, :) = get_state(state_index, state_handle)
         end do

         do iens = 1, ens_size
            call quad_bilinear_interp(lon_lat_alt(1), lon_lat_alt(2), x_neighbors_quad, &
                                       y_neighbors_quad, cyclic, quad_temp_values(:,iens), &
                                       above_values(iens))
         enddo

         ! Get quad temp_values for the below level
         do ineighbor =1, nquad_neighbors
            state_index = get_dart_vector_index(quad_neighbor_indices(icolumn, ineighbor), below_index, no_third_dimension, dom_id, varid)
            quad_temp_values(ineighbor, :) = get_state(state_index, state_handle)
         end do

         do iens = 1, ens_size
            call quad_bilinear_interp(lon_lat_alt(1), lon_lat_alt(2), x_neighbors_quad, &
                                          y_neighbors_quad, cyclic, quad_temp_values(:,iens), &
                                          below_values(iens))
         enddo

         call vert_interp(ens_size, below_values, above_values, fraction, expected_obs)

         exit
      end if

   end do

end if

! All good.
istatus(:) = 0


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

! Local variables

integer :: lev_index, col_index
integer :: my_var_id, my_qty

if ( .not. module_initialized ) call static_init_model

call get_model_variable_indices(index_in, col_index, lev_index, no_third_dimension, &
                                var_id=my_var_id, kind_index=my_qty)

! should be set to the actual location using set_location()
location = set_location(center_longitude(col_index), center_latitude(col_index), center_altitude(lev_index), VERTISHEIGHT)

! should be set to the physical quantity, e.g. QTY_TEMPERATURE
if (present(qty)) qty = my_qty

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
call nc_add_global_attribute(ncid, "model", "template")

call nc_end_define_mode(ncid)

! Flush the buffer and leave netCDF file open
call nc_synchronize_file(ncid)

end subroutine nc_write_model_atts

!-----------------------------------------------------------------------
! Parse the table of variables characteristics into arrays for easier access.

function assign_var(variables, MAX_STATE_VARIABLES) result(var)

character(len=vtablenamelength), intent(in) :: variables(:, :)
integer, intent(in)                         :: MAX_STATE_VARIABLES

type(var_type) :: var
integer        :: ivar
character(len=vtablenamelength) :: table_entry

!-----------------------------------------------------------------------
! Codes for interpreting the NUM_STATE_TABLE_COLUMNS of the variables table
integer, parameter :: NAME_INDEX      = 1 ! ... variable name
integer, parameter :: QTY_INDEX       = 2 ! ... DART qty
integer, parameter :: MIN_VAL_INDEX   = 3 ! ... minimum value if any
integer, parameter :: MAX_VAL_INDEX   = 4 ! ... maximum value if any
integer, parameter :: UPDATE_INDEX    = 5 ! ... update (state) or not

! Loop through the variables array to get the actual count of the number of variables
do ivar = 1, MAX_STATE_VARIABLES
   ! If the element is an empty string, the loop has exceeded the extent of the variables
   if (variables(1, ivar) == '') then
      var%count = ivar-1
      exit
   endif 
enddo

! Allocate the arrays in the var derived type
allocate(var%names(var%count), var%qtys(var%count), var%clamp_values(var%count, 2), var%updates(var%count))

do ivar = 1, var%count
   
   var%names(ivar) = trim(variables(NAME_INDEX, ivar))

   table_entry = variables(QTY_INDEX, ivar)
   call to_upper(table_entry)

   var%qtys(ivar) = get_index_for_quantity(table_entry)

   if (variables(MIN_VAL_INDEX, ivar) /= 'NA') then
      read(variables(MIN_VAL_INDEX, ivar), '(d16.8)') var%clamp_values(ivar,1)
   else
      var%clamp_values(ivar,1) = MISSING_R8
   endif

   if (variables(MAX_VAL_INDEX, ivar) /= 'NA') then
      read(variables(MAX_VAL_INDEX, ivar), '(d16.8)') var%clamp_values(ivar,2)
   else
      var%clamp_values(ivar,2) = MISSING_R8
   endif

   table_entry = variables(UPDATE_INDEX, ivar)
   call to_upper(table_entry)

   if (table_entry == 'UPDATE') then
      var%updates(ivar) = .true.
   else
      var%updates(ivar) = .false.
   endif

enddo
   
end function assign_var

subroutine read_geometry_file()

   integer               :: dimid, varid
   character(len=256)      :: name

   type(file_type)       :: geometry_file
   character(len=256)    :: restart_directory, grid_directory, filter_directory
   namelist /directory_nml/ restart_directory, grid_directory, filter_directory

   call find_namelist_in_file('input.nml', 'directory_nml', iunit)
   read(iunit, nml = directory_nml, iostat = io)
   call check_namelist_read(iunit, io, 'directory_nml')

   geometry_file%file_path = trim(filter_directory) // 'geometry_file.nc'

   geometry_file%ncid = nc_open_file_readonly(geometry_file%file_path)

   ! attributes
   geometry_file%ncstatus = nf90_get_att(geometry_file%ncid, NF90_GLOBAL, 'index_of_north_pole_quad_column', inorth_pole_quad_column)
   geometry_file%ncstatus = nf90_get_att(geometry_file%ncid, NF90_GLOBAL, 'index_of_south_pole_quad_column', isouth_pole_quad_column)

   ! dimensions
   geometry_file%ncstatus = nf90_inq_dimid(geometry_file%ncid, 'center_altitudes', dimid)
   geometry_file%ncstatus = nf90_inquire_dimension(geometry_file%ncid, dimid, name, ncenter_altitudes)

   geometry_file%ncstatus = nf90_inq_dimid(geometry_file%ncid, 'quad_columns', dimid)
   geometry_file%ncstatus = nf90_inquire_dimension(geometry_file%ncid, dimid, name, nquad_columns)

   geometry_file%ncstatus = nf90_inq_dimid(geometry_file%ncid, 'center_columns', dimid)
   geometry_file%ncstatus = nf90_inquire_dimension(geometry_file%ncid, dimid, name, ncenter_columns)

   ! allocate arrays
   allocate(center_altitude(ncenter_altitudes))
   allocate(center_latitude(ncenter_columns))
   allocate(center_longitude(ncenter_columns))

   allocate(vertex_neighbor_indices(nvertex_columns, nvertex_neighbors))
   allocate(quad_neighbor_indices(nquad_columns, nquad_neighbors))

   ! variables
   geometry_file%ncstatus = nf90_inq_varid(geometry_file%ncid, 'center_altitude', varid)
   geometry_file%ncstatus = nf90_get_var(geometry_file%ncid, varid, center_altitude)

   geometry_file%ncstatus = nf90_inq_varid(geometry_file%ncid, 'center_longitude', varid)
   geometry_file%ncstatus = nf90_get_var(geometry_file%ncid, varid, center_longitude)

   geometry_file%ncstatus = nf90_inq_varid(geometry_file%ncid, 'center_latitude', varid)
   geometry_file%ncstatus = nf90_get_var(geometry_file%ncid, varid, center_latitude)

   geometry_file%ncstatus = nf90_inq_varid(geometry_file%ncid, 'vertex_neighbor_indices', varid)
   geometry_file%ncstatus = nf90_get_var(geometry_file%ncid, varid, vertex_neighbor_indices)

   geometry_file%ncstatus = nf90_inq_varid(geometry_file%ncid, 'quad_neighbor_indices', varid)
   geometry_file%ncstatus = nf90_get_var(geometry_file%ncid, varid, quad_neighbor_indices)

   call nc_close_file(geometry_file%ncid)

end subroutine read_geometry_file

! Barycentric procedures

subroutine inside_triangle(t1, t2, t3, r, lat, lon, inside, weights)

   ! given 3 corners of a triangle and an xyz point, compute whether
   ! the point is inside the triangle.  this assumes r is coplanar
   ! with the triangle - the caller must have done the lat/lon to
   ! xyz conversion with a constant radius and then this will be
   ! true (enough).  sets inside to true/false, and returns the
   ! weights if true.  weights are set to 0 if false.
   
   real(r8), intent(in)  :: t1(3), t2(3), t3(3)
   real(r8), intent(in)  :: r(3), lat, lon
   logical,  intent(out) :: inside
   real(r8), intent(out) :: weights(3)
   
   ! check for degenerate cases first - is the test point located
   ! directly on one of the vertices?  (this case may be common
   ! if we're computing on grid point locations.)
   if (all(abs(r - t1) < roundoff)) then
      inside = .true.
      weights = (/ 1.0_r8, 0.0_r8, 0.0_r8 /)
      return
   else if (all(abs(r - t2) < roundoff)) then
      inside = .true.
      weights = (/ 0.0_r8, 1.0_r8, 0.0_r8 /)
      return
   else if (all(abs(r - t3) < roundoff)) then
      inside = .true.
      weights = (/ 0.0_r8, 0.0_r8, 1.0_r8 /)
      return
   endif
   
   ! not a vertex. compute the weights.  if any are
   ! negative, the point is outside.  since these are
   ! real valued computations define a lower bound for
   ! numerical roundoff error and be sure that the
   ! weights are not just *slightly* negative.
   call get_3d_weights(r, t1, t2, t3, weights)
   
   if (any(weights < -roundoff)) then
      inside = .false.
      weights = 0.0_r8
      return
   endif
   
   ! truncate barely negative values to 0
   inside = .true.
   where (weights < 0.0_r8) weights = 0.0_r8
   return
   
end subroutine inside_triangle

subroutine get_3d_weights_old(p, v1, v2, v3, lat, lon, weights)

   ! Given a point p (x,y,z) inside a triangle, and the (x,y,z)
   ! coordinates of the triangle corner points (v1, v2, v3),
   ! find the weights for a barycentric interpolation.  this
   ! computation only needs two of the three coordinates, so figure
   ! out which quadrant of the sphere the triangle is in and pick
   ! the 2 axes which are the least planar:
   !  (x,y) near the poles,
   !  (y,z) near 0 and 180 longitudes near the equator,
   !  (x,z) near 90 and 270 longitude near the equator.
   ! (lat/lon are the coords of p. we could compute them here
   ! but since in all cases we already have them, pass them
   ! down for efficiency)
   
   real(r8), intent(in)  :: p(3)
   real(r8), intent(in)  :: v1(3), v2(3), v3(3)
   real(r8), intent(in)  :: lat, lon
   real(r8), intent(out) :: weights(3)
   
   real(r8) :: cxs(3), cys(3)
   
   ! above or below 45 in latitude, where -90 < lat < 90:
   if (lat >= 45.0_r8 .or. lat <= -45.0_r8) then
      cxs(1) = v1(1)
      cxs(2) = v2(1)
      cxs(3) = v3(1)
      cys(1) = v1(2)
      cys(2) = v2(2)
      cys(3) = v3(2)
      !call get_barycentric_weights(p(1), p(2), cxs, cys, weights)
      return
   endif
   
   ! nearest 0 or 180 in longitude, where 0 < lon < 360:
   if ( lon <= 45.0_r8 .or. lon >= 315.0_r8 .or. &
       (lon >= 135.0_r8 .and. lon <= 225.0_r8)) then
      cxs(1) = v1(2)
      cxs(2) = v2(2)
      cxs(3) = v3(2)
      cys(1) = v1(3)
      cys(2) = v2(3)
      cys(3) = v3(3)
      !call get_barycentric_weights(p(2), p(3), cxs, cys, weights)
      return
   endif
   
   ! last option, nearest 90 or 270 in lon:
   cxs(1) = v1(1)
   cxs(2) = v2(1)
   cxs(3) = v3(1)
   cys(1) = v1(3)
   cys(2) = v2(3)
   cys(3) = v3(3)
   !call get_barycentric_weights(p(1), p(3), cxs, cys, weights)
   
end subroutine get_3d_weights_old

subroutine get_3d_weights(p, v1, v2, v3, weights)

   ! MEG: 
   ! This replaces 'get_3d_weights_old' which to me it seems
   ! like it's trying to guess which side to draw flat based on the address.

   ! The code below looks at the triangle itself and picks the best 
   ! angle so it is flattest and least distorted. I also added a 'success' 
   ! option in 'get_barycentric_weights' to check for the collinearity issue. 

   ! Obviously, we can do better, but for now I think this keeps us going.

   real(r8), intent(in)  :: p(3), v1(3), v2(3), v3(3)
   real(r8), intent(out) :: weights(3)

   real(r8) :: e1(3), e2(3), n(3)
   real(r8) :: cxs(3), cys(3)
   logical  :: success

   ! Compute triangle edge vectors
   e1(1) = v2(1) - v1(1)
   e1(2) = v2(2) - v1(2)
   e1(3) = v2(3) - v1(3)

   e2(1) = v3(1) - v1(1)
   e2(2) = v3(2) - v1(2)
   e2(3) = v3(3) - v1(3)

   ! Cross product (normal vector)
   n(1) = e1(2)*e2(3) - e1(3)*e2(2)
   n(2) = e1(3)*e2(1) - e1(1)*e2(3)
   n(3) = e1(1)*e2(2) - e1(2)*e2(1)

   ! Try projection on XY plane
   if (abs(n(3)) >= abs(n(1)) .and. abs(n(3)) >= abs(n(2))) then
      cxs(1) = v1(1); cxs(2) = v2(1); cxs(3) = v3(1)
      cys(1) = v1(2); cys(2) = v2(2); cys(3) = v3(2)
      call get_barycentric_weights(p(1), p(2), cxs, cys, weights, success)
      if (success) return
   endif

   ! Try projection on YZ plane
   if (abs(n(1)) >= abs(n(2))) then
      cxs(1) = v1(2); cxs(2) = v2(2); cxs(3) = v3(2)
      cys(1) = v1(3); cys(2) = v2(3); cys(3) = v3(3)
      call get_barycentric_weights(p(2), p(3), cxs, cys, weights, success)
      if (success) return
   endif

   ! Fallback: try projection on XZ plane
   cxs(1) = v1(1); cxs(2) = v2(1); cxs(3) = v3(1)
   cys(1) = v1(3); cys(2) = v2(3); cys(3) = v3(3)
   call get_barycentric_weights(p(1), p(3), cxs, cys, weights, success)

   if (.not. success) then
      ! Tried all planes and it didn't work out, so ...
      print *, 'Carefull: get_3d_weights failed to compute weights.'
      weights(1) = -1.0_r8
      weights(2) = -1.0_r8
      weights(3) = -1.0_r8
   endif

end subroutine get_3d_weights


subroutine get_barycentric_weights(x, y, cxs, cys, weights, success)

   ! Computes the barycentric weights for a 2d interpolation point
   ! (x,y) in a 2d triangle with the given (cxs,cys) corners.
   
   real(r8), intent(in)  :: x, y, cxs(3), cys(3)
   real(r8), intent(out) :: weights(3)
   logical, intent(out)  :: success
   
   real(r8) :: denom
   
   ! Get denominator
   denom = (cys(2) - cys(3)) * (cxs(1) - cxs(3)) + &
      (cxs(3) - cxs(2)) * (cys(1) - cys(3))

   ! If the vertices are collinear then the triangle is degenerate
   if (abs(denom) < roundoff) then
      success = .false.
      weights = 0.0_r8
      return
   endif
 
   weights(1) = ((cys(2) - cys(3)) * (x - cxs(3)) + &
                (cxs(3) - cxs(2)) * (y - cys(3))) / denom
   
   weights(2) = ((cys(3) - cys(1)) * (x - cxs(3)) + &
                (cxs(1) - cxs(3)) * (y - cys(3))) / denom
   
   weights(3) = 1.0_r8 - weights(1) - weights(2)
   
   if (any(abs(weights) < roundoff)) then
      where (abs(weights) < roundoff) weights = 0.0_r8
      where (abs(1.0_r8 - abs(weights)) < roundoff) weights = 1.0_r8
   endif

   ! I'm here so it must be a good looking triangle
   success = .true.  
 
end subroutine get_barycentric_weights

!-----------------------------------------------------------------------

subroutine latlon_to_xyz(lat, lon, x, y, z)

   ! Given a lat, lon in degrees, return the cartesian x,y,z coordinate
   ! on the surface of a specified radius relative to the origin
   ! at the center of the earth.
   
   real(r8), intent(in)  :: lat, lon
   real(r8), intent(out) :: x, y, z
   
   real(r8) :: rlat, rlon
   
   rlat = lat * deg2rad
   rlon = lon * deg2rad
   
   x = radius * cos(rlon) * cos(rlat)
   y = radius * sin(rlon) * cos(rlat)
   z = radius * sin(rlat)
   
end subroutine latlon_to_xyz


!-----------------------------------------------------------------------
! Determines if the projection of a point p onto the plane of a triangle with vertices
! v1, v2 and v3 is inside the triangle or not. Computes the areas of each of the triangles
! between p and a pair of vertices. These should sum to the area of the triangle if p
! is inside and be larger than that if p is outside. 

function is_point_in_triangle(v1, v2, v3, p)

logical              :: is_point_in_triangle
real(r8), intent(in) :: v1(3), v2(3), v3(3), p(3)    

real(r8) :: a(3), b(3), perp(3), unit_perp(3), p_proj(3)
real(r8) :: offset, len_s1, len_s2, len_s3, len_p1, len_p2, len_p3, at, at1, at2, at3
real(r8) :: area_dif, dif_frac, threshold

! Get the projection of the point p onto the plane containing the triangle
! Start by getting perpendicular vector to plane by cross product
a = v1 - v2
b = v2 - v3
perp(1) = a(2) * b(3) - a(3) * b(2)
perp(2) = a(3) * b(1) - a(1) * b(3)
perp(3) = a(1) * b(2) - a(2) * b(1)
! Get unit vector in direction of perp
unit_perp = perp / sqrt(dot_product(perp, perp))
! Projection of vector from v1 to p on the unit perp vector is how much to move to get to plane
offset = dot_product((p-v1), unit_perp)
p_proj = p - offset * unit_perp

! Compute lengths of the sides
len_s1 = sqrt(dot_product(v1-v2, v1-v2))
len_s2 = sqrt(dot_product(v3-v2, v3-v2))
len_s3 = sqrt(dot_product(v1-v3, v1-v3))

! Compute the lengths from the point p
len_p1 = sqrt(dot_product(p_proj-v1, p_proj-v1))
len_p2 = sqrt(dot_product(p_proj-v2, p_proj-v2))
len_p3 = sqrt(dot_product(p_proj-v3, p_proj-v3))

! Area of triangle
at = heron(len_s1, len_s2, len_s3)

! Compute areas of sub triangles
at1 = heron(len_p1, len_p2, len_s1)
at2 = heron(len_p2, len_p3, len_s2)
at3 = heron(len_p3, len_p1, len_s3)

! Difference between sub triangles and the triangle area
area_dif = at1 + at2 + at3 - at

! Quadrilaterals on the interior of the cube sphere sides are really spherical quads,
! There sides are great circles. This routine assumes that the triangles composing the quads
! have straight sides in regular space. The algorithm finds points that are inside the 
! spherical quads. These quads actually 'bulge' out compared to the regular sides, so it is possible
! to have points that are inside the spherical quad but just barely outside of the regular
! quads. This threshold is tuned so that these points still show as inside. The tuning is for
! np = 18 (number of points along a grid face is 18). Fewer points might require a larger
! threshold while more points might be okay with a smaller one.
threshold = 0.002_r8

dif_frac = area_dif / at
is_point_in_triangle = abs(dif_frac) < threshold

end function is_point_in_triangle
!-----------------------------------------------------------------------

function lat_lon_to_xyz(lat, lon)

real(r8)             :: lat_lon_to_xyz(3)
real(r8), intent(in) :: lat, lon

lat_lon_to_xyz(1) = cos(lat) * cos(lon)
lat_lon_to_xyz(2) = cos(lat) * sin(lon)
lat_lon_to_xyz(3) = sin(lat)

end function lat_lon_to_xyz

!-----------------------------------------------------------------------

! Given the face (from 0 to 5) the number of lons/lats on a face (np) and
! the indices of the i and j grid point, returns the latitude and longitude of the point

subroutine grid_to_lat_lon(face, lat_ind, lon_ind, np, lat, lon)

integer,  intent(in)  :: face, lat_ind, lon_ind, np
real(r8), intent(out) :: lat, lon

real(r8) :: cube_side, del, half_del, x, y, blon, blat, rot_angle
real(r8) :: vect(3), RZ(3, 3), rot_vect(3)

! Cube side is divided into np-1 interior intervals of width 2sqrt(1/3) / np and
! two exterior intervals of half  width, sqrt(1/3) / np 
cube_side = 2.0_r8 * sqrt(1.0_r8 / 3.0_r8)
del = cube_side / np
half_del = del / 2.0_r8

x = sqrt(1.0_r8/3.0_r8) - (half_del + del * (lon_ind - 1))
if(face == 5) then 
   y =  sqrt(1.0_r8/3.0_r8) - (half_del + del * (lat_ind - 1))
else
   y = -sqrt(1.0_r8/3.0_r8) + (half_del + del * (lat_ind - 1))
endif

if(face < 4) then
   blon = atan2(sqrt(1.0_r8 / 3.0_r8), x)
   blat = atan2(y, sqrt(1.0_r8/3.0_r8 + x**2.0_r8))
   blon = blon - pi/4.0_r8

   ! Above is for face 0; add pi/2 for each additional face tangent to equator
   lon = blon + pi/2.0_r8 * face; 
   lat = blat;
elseif(face == 4 .or. face == 5) then
   ! Face 4 is tangent to south pole
   lon = atan2(y, x)
   lat = atan2(sqrt(1.0_r8/3.0_r8), sqrt(x**2.0_r8 + y**2.0_r8))

   if(face == 4) lat = -lat

   !  Get ready for rotation
   vect = lat_lon_to_xyz(lat, lon)

   ! Then rotate 45 degrees around Z
   rot_angle = -pi/4.0_r8;

   ! Create the rotation matrix
   RZ(1, 1:3) = [cos(rot_angle),  sin(rot_angle), 0.0_r8]
   RZ(2, 1:3) = [-sin(rot_angle), cos(rot_angle), 0.0_r8]
   RZ(3, 1:3) = [0.0_r8,          0.0_r8,         1.0_r8]
   rot_vect = matmul(RZ, vect)

   lat = asin(rot_vect(3))
   lon = atan2(rot_vect(2), rot_vect(1))
   ! Note that there are inconsistent treatments of the value near longitude
   ! 0 in the grid files for Aether. Some points have a value near or just less
   ! than 2pi, other points have values just greater than 0. This code 
   ! avoids values near to 2pi and have 0 instead.
   if(lon < 0.0_r8) lon =  lon + 2.0_r8*PI
   if(lon >= 2.0_r8*pi) lon = 0.0_r8

endif

end subroutine grid_to_lat_lon

!-----------------------------------------------------------------------

subroutine fix_face(face, lat_grid, lon_grid, np, f_face, f_lat_grid, f_lon_grid, edge, corner)

integer, intent(in)  :: face, lat_grid, lon_grid, np
integer, intent(out) :: f_face, f_lat_grid, f_lon_grid
logical, intent(out) :: edge, corner

integer :: left_neighbor(6),   right_neighbor(6)
integer :: bottom_neighbor(6), top_neighbor(6)
integer :: left_lon_grid(6),   right_lon_grid(6)
integer :: left_lat_grid(6),   right_lat_grid(6)
integer :: bottom_lon_grid(6), top_lon_grid(6)
integer :: bottom_lat_grid(6), top_lat_grid(6)
! For points past the edge of a face, finds the corresponding points on the adjacent face
! Need to do something more special for corner points

! Default is not a corner or an edge
corner = .false.
edge   = .false.

! Just return if no edge
if(lon_grid > 0 .and. lon_grid < np + 1 .and. lat_grid > 0 .and. lat_grid < np + 1) then
   f_face = face
   f_lon_grid = lon_grid
   f_lat_grid = lat_grid
   return
endif

if((lat_grid == 0 .or. lat_grid == np + 1) .and. (lon_grid == 0 .or. lon_grid == np + 1)) then
   corner = .true.
   f_face     = -99
   f_lon_grid = -99
   f_lat_grid = -99
   return
endif

! Otherwise, on an edge
edge = .true.
! Deal with each side of faces separately
if(lon_grid == 0) then
   ! On left edge not corner
   left_neighbor = [3, 0, 1, 2, 0, 0]
   f_face = left_neighbor(face + 1)
   left_lon_grid = [np, np, np, np, lat_grid, np+1-lat_grid]
   f_lon_grid = left_lon_grid(face + 1)
   left_lat_grid = [lat_grid, lat_grid, lat_grid, lat_grid, 1, np]
   f_lat_grid = left_lat_grid(face + 1)
elseif(lon_grid == np + 1) then
   ! On right edge not corner
   right_neighbor = [1, 2, 3, 0, 2, 2]
   f_face = right_neighbor(face + 1)
   right_lon_grid = [1, 1, 1, 1, np+1-lat_grid, lat_grid]
   f_lon_grid = right_lon_grid(face + 1)
   right_lat_grid = [lat_grid, lat_grid, lat_grid, lat_grid, 1, np]
   f_lat_grid = right_lat_grid(face + 1)
elseif(lat_grid == 0) then
   ! On bottom edge not corner
   bottom_neighbor = [4, 4, 4, 4, 3, 1]
   f_face = bottom_neighbor(face + 1)
   bottom_lon_grid = [1, lon_grid, np, np+1-lon_grid, np+1-lon_grid, lon_grid]
   f_lon_grid = bottom_lon_grid(face + 1)
   bottom_lat_grid = [lon_grid, np, np+1-lon_grid, 1, 1, np]
   f_lat_grid = bottom_lat_grid(face + 1)
elseif(lat_grid == np + 1) then
   ! On top edge not corner
   top_neighbor = [5, 5, 5, 5, 1, 3]
   f_face = top_neighbor(face + 1)
   top_lon_grid = [1, lon_grid, np, np+1-lon_grid, lon_grid, np+1-lon_grid]
   f_lon_grid = top_lon_grid(face + 1)
   top_lat_grid = [np+1-lon_grid, 1, lon_grid, np, 1, np]
   f_lat_grid = top_lat_grid(face + 1)
endif

end subroutine fix_face

!-----------------------------------------------------------------------

subroutine get_face(lat, lon_in, face, len)

real(r8), intent(in)  :: lat, lon_in
integer,  intent(out) :: face
real(r8), intent(out) :: len(2)

integer :: side, rside, rside2
real(r8) :: inv_sqrt_2, rlon, rlon2, gama, gamb, lon
real(r8) :: vec(3), rot_vec(3), rot_vec2(3), rot(3, 3), rot2(3, 3), lon_grid(2), lon_grid_m(2)
! Returns which face contains (lat, lon) and the length from the edge of the point
! along each of the great circle axes.

! Range adjustment
lon = lon_in
if(lon >= 2.0_r8*PI) lon = 0.0_r8

! Convert lat lon to x y z on unit sphere
vec = lat_lon_to_xyz(lat, lon)

! Get the longitudes for this point in the two rotated spaces

!====================================================================
! Following code shows individual rotations;
! Single matrix at the end multiplies them together off-line for single rotation
! Can collapse these to a single rotation vector for efficiency
! Rotation 90 degrees around y to put pole on equator
!RY = [cosd(90) 0 -sind(90); 0 1 0; sind(90) 0 cosd(90)];
! Rotate 45 degrees around x
!RX = [1 0 0; 0 cosd(45) sind(45); 0 -sind(45) cosd(45)];
! Then rotate 45 degrees around Z
!RZ = [cosd(45) sind(45) 0; -sind(45) cosd(45) 0; 0 0 1];
! Get longitude in the two rotated spaces
!rot_vec = RZ * RX * RY * vec;
!====================================================================
inv_sqrt_2 = 1.0_r8 / sqrt(2.0_r8);
rot(1, 1:3) = [0.5_r8,     0.5_r8,      -inv_sqrt_2]
rot(2, 1:3) = [0.5_r8,     0.5_r8,      inv_sqrt_2]
rot(3, 1:3) = [inv_sqrt_2, -inv_sqrt_2, 0.0_r8]
rot_vec = matmul(rot, vec)
! Compute the longitude in the rotated space
rlon = atan2(rot_vec(2), rot_vec(1))
if(rlon < 0.0_r8) rlon = rlon + 2.0_r8*PI

!====================================================================
! Can collapse these to a single rotation vector for efficiency
! Rotation 90 degrees around y to put pole on equator
!RY = [cosd(90) 0 -sind(90); 0 1 0; sind(90) 0 cosd(90)];
! Rotate -45 degrees around x
!RX2 = [1 0 0; 0 cosd(-45) sind(-45); 0 -sind(-45) cosd(-45)];
! Then rotate 45 degrees around Z
!RZ = [cosd(45) sind(45) 0; -sind(45) cosd(45) 0; 0 0 1];
! Get longitude in the two rotated spaces
!rot_vec2 = RZ * RX2 * RY * vec;
!====================================================================
rot2(1, 1:3) = [-0.5_r8,     0.5_r8,        -inv_sqrt_2]
rot2(2, 1:3) = [-0.5_r8,     0.5_r8,        inv_sqrt_2]
rot2(3, 1:3) = [inv_sqrt_2,  inv_sqrt_2,    0.0_r8]
rot_vec2 = matmul(rot2, vec)
! Compute the longitude in the rotated space
rlon2 = atan2(rot_vec2(2), rot_vec2(1))
if(rlon2 < 0.0_r8) rlon2 = rlon2 + 2.0_r8*PI

! Which non-polar side could we be on, 1 to 4
side = floor(lon / (PI/2.0_r8)) + 1.0_r8
! Which rotated 1 side are we on 
rside = floor(rlon / (PI/2.0_r8)) + 1.0_r8
! Which rotated 2 side
rside2 = floor(rlon2 / (PI/2.0_r8)) + 1.0_r8

! Figure out the face from here (0 to 5, 4 is south, 5 is north)
! These are consistent with the numbering on Aether grid files for the cubed sphere
if    ( side == 1 .and. rside  == 1) then
   face = 0; lon_grid(1) = lon;  lon_grid(2) = rlon
elseif( side == 2 .and. rside2 == 1) then
   face = 1; lon_grid(1) = lon;  lon_grid(2) = rlon2
elseif( side == 3 .and. rside  == 3) then 
   face = 2; lon_grid(1) = lon;  lon_grid(2) = rlon
elseif( side == 4 .and. rside2 == 3) then
   face = 3; lon_grid(1) = lon;  lon_grid(2) = rlon2
elseif(rside == 4 .and. rside2 == 4) then
   face = 4; lon_grid(1) = rlon; lon_grid(2) = rlon2
elseif(rside == 2 .and. rside2 == 2) then
   face = 5; lon_grid(1) = rlon; lon_grid(2) = rlon2
endif

! Can also use the fact that the projection is equidistant on the imbedded cube to get what fraction 
! across the imbedded rectangle we are
! Take the longitudes and turn them into a number between -sqrt(1/3) and sqrt(1/3)
lon_grid_m = mod(lon_grid, PI/2.0_r8)

! Use law of sines to go from lon back to position along edge of imbedded cube 
! The triangle of interest has a side of length 2sqrt(1/3) (1/2 of the planar diagonal of the imbedded cube)
! The angles adjacent to this side are the longitude and 45 degrees
! The angle opposite the side of length 2sqrt(1/3) is pi - (longitude + pi/4)
! The side opposite the longitude is how far along the side of the cube
! The cube side is 2sqrt(1/3), so the length along the side is between zero and this value
gama = PI - (PI/4.0_r8 + lon_grid_m(1))
len(1) = sqrt(2.0_r8/3.0_r8) * sin(lon_grid_m(1)) / sin(gama)

gamb = PI - (PI/4 + lon_grid_m(2))
len(2) = sqrt(2.0_r8/3.0_r8) * sin(lon_grid_m(2)) / sin(gamb)

! If we are on sides 2 or 3, the lengths need to be modified because the grid storage
! for Aether goes from smallest latitude to largest and the longitudes of the shifted
! poles are going the opposite way
if(face == 2 .or. face == 3) len(2) = 2.0_r8 * sqrt(1.0_r8/3.0_r8) - len(2)

! Same for face 4 (the bottom) but it's the other coordinate that's reversed
if(face == 4) len(1) = 2.0_r8*sqrt(1.0_r8/3.0_r8) - len(1)

end subroutine get_face

!-----------------------------------------------------------------------
subroutine get_corners(face, lat_grid, lon_grid, lat, lon, np,   &
   f_face, f_lat_grid, f_lon_grid, num_bound_points)
   
integer, intent(in) :: face, lat_grid, lon_grid, np
real(r8), intent(in) :: lat, lon
integer, intent(out) :: f_face(4), f_lat_grid(4), f_lon_grid(4), num_bound_points

integer :: corner, quad, i
integer :: quad_lon_grid(3, 4), quad_lat_grid(3, 4), quad_face(3, 4)
real(r8) :: pxyz(3), qxyz(4, 3), grid_pt_lat, grid_pt_lon

! Checks to see if the point under consideration is at a corner
! If it is, return the face, lat_index, and lon_index for each of the three bounding points

! Default is to find a triangle
num_bound_points = 3

if(face == 0) then
   if    (lat_grid == 0    .and. lon_grid == 0   ) then
      corner = 1
   elseif(lat_grid == 0    .and. lon_grid == np+1) then
      corner = 2
   elseif(lat_grid == np+1 .and. lon_grid == 0   ) then
      corner = 5
   else
      corner = 6;
   endif
elseif(face == 1) then
   if    (lat_grid == 0    .and. lon_grid == 0   ) then
      corner = 2
   elseif(lat_grid == 0    .and. lon_grid == np+1) then 
      corner = 3
   elseif(lat_grid == np+1 .and. lon_grid == 0   ) then
      corner = 6
   else
      corner = 7
   endif
elseif(face == 2) then
   if    (lat_grid == 0    .and. lon_grid == 0   ) then
      corner = 3
   elseif(lat_grid == 0    .and. lon_grid == np+1) then
      corner = 4
   elseif(lat_grid == np+1 .and. lon_grid == 0   ) then
      corner = 7
   else
      corner = 8
   endif
elseif(face == 3) then
   if    (lat_grid == 0    .and. lon_grid == 0   ) then
      corner = 4
   elseif(lat_grid == 0    .and. lon_grid == np+1) then
      corner = 1
   elseif(lat_grid == np+1 .and. lon_grid == 0   ) then
      corner = 8
   else
      corner = 5;
   endif
elseif(face == 4) then
   if    (lat_grid == 0    .and. lon_grid == 0   ) then
      corner = 1
   elseif(lat_grid == 0    .and. lon_grid == np+1) then
      corner = 4
   elseif(lat_grid == np+1 .and. lon_grid == 0   ) then
      corner = 2
   else
      corner = 3
   endif
elseif(face == 5) then
   if    (lat_grid == 0    .and. lon_grid == 0   ) then
      corner = 6
   elseif(lat_grid == 0    .and. lon_grid == np+1) then
      corner = 7
   elseif(lat_grid == np+1 .and. lon_grid == 0   ) then
      corner = 5
   else
      corner = 8;
   endif
endif

! Harvest the information on the grid points bounding the appropriate corner
! Arrays of info for adjacent quads for bulges (three of them, first index)
quad_lon_grid(1:3, 1:4)  = -99
quad_lat_grid(1:3, 1:4)  = -99
quad_face(1:3, 1:4) = -99

if(corner == 1) then
   f_face(1:3) =     [3,  0, 4]
   f_lon_grid(1:3) = [np, 1, 1]
   f_lat_grid(1:3) = [1,  1, 1]
   quad_face(1, 1:4) = [3,  0, 0, 3]
   quad_face(2, 1:4) = [0,  0, 4, 4]
   quad_face(3, 1:4) = [3,  3, 4, 4]
   quad_lat_grid(1, 1:4) = [1,  1, 2, 2]
   quad_lat_grid(2, 1:4) = [1,  1, 1, 2]
   quad_lat_grid(3, 1:4) = [1,  1, 1, 1]
   quad_lon_grid(1, 1:4) = [np,   1,  1, np]
   quad_lon_grid(2, 1:4) = [1,    2,  1, 1 ]
   quad_lon_grid(3, 1:4) = [np-1, np, 1, 2 ]
elseif(corner == 2) then
   f_face(1:3) =     [0,  1, 4 ]
   f_lon_grid(1:3) = [np, 1, 1 ]
   f_lat_grid(1:3) = [1,  1, np]
   quad_face(1, 1:4) = [0, 1, 1, 0]
   quad_face(2, 1:4) = [1, 1, 4, 4]
   quad_face(3, 1:4) = [0, 0, 4, 4]
   quad_lat_grid(1, 1:4) = [1, 1, 2,  2   ]
   quad_lat_grid(2, 1:4) = [1, 1, np, np  ]
   quad_lat_grid(3, 1:4) = [1, 1, np, np-1]
   quad_lon_grid(1, 1:4) = [np,   1,  1,  np]
   quad_lon_grid(2, 1:4) = [1,    2,  2,  1 ]
   quad_lon_grid(3, 1:4) = [np-1, np, 1,  1 ]
elseif(corner == 3) then
   f_face(1:3) =     [1,  2, 4 ]
   f_lon_grid(1:3) = [np, 1, np]
   f_lat_grid(1:3) = [1,  1, np]
   quad_face(1, 1:4) = [1, 2, 2, 1]
   quad_face(2, 1:4) = [2, 2, 4, 4]
   quad_face(3, 1:4) = [1, 1, 4, 4]
   quad_lat_grid(1, 1:4) = [1, 1, 2,    2 ]
   quad_lat_grid(2, 1:4) = [1, 1, np-1, np]
   quad_lat_grid(3, 1:4) = [1, 1, np,   np]
   quad_lon_grid(1, 1:4) = [np,   1,  1,  np  ]
   quad_lon_grid(2, 1:4) = [1,    2,  np, np  ]
   quad_lon_grid(3, 1:4) = [np-1, np, np, np-1]
elseif(corner == 4) then
   f_face(1:3) =     [2,  3, 4 ]
   f_lon_grid(1:3) = [np, 1, np]
   f_lat_grid(1:3) = [1,  1, 1 ]
   quad_face(1, 1:4) = [2, 3, 3,2]
   quad_face(2, 1:4) = [3, 3, 4, 4]
   quad_face(3, 1:4) = [2, 2, 4, 4]
   quad_lat_grid(1, 1:4) = [1, 1, 2, 2]
   quad_lat_grid(2, 1:4) = [1, 1, 1, 1]
   quad_lat_grid(3, 1:4) = [1, 1, 1, 2]
   quad_lon_grid(1, 1:4) = [np,   1,  1,    np]
   quad_lon_grid(2, 1:4) = [1,    2,  np-1, np]
   quad_lon_grid(3, 1:4) = [np-1, np, np,   np]
elseif(corner == 5) then
   f_face(1:3) =     [3,  0,  5 ]
   f_lon_grid(1:3) = [np, 1,  1 ]
   f_lat_grid(1:3) = [np, np, np]
   quad_face(1, 1:4) = [3, 0, 0, 3]
   quad_face(2, 1:4) = [0, 0, 5, 5]
   quad_face(3, 1:4) = [3, 3, 5, 5]
   quad_lat_grid(1, 1:4) = [np-1, np-1, np, np  ]
   quad_lat_grid(2, 1:4) = [np,   np,   np, np-1]
   quad_lat_grid(3, 1:4) = [np,   np,   np, np  ]
   quad_lon_grid(1, 1:4) = [np,   1,  1,  np]
   quad_lon_grid(2, 1:4) = [1,    2,  1,  1 ]
   quad_lon_grid(3, 1:4) = [np-1, np, 1,  2 ]
elseif(corner == 6) then
   f_face(1:3) =     [0,  1, 5 ]
   f_lon_grid(1:3) = [np, 1,  1]
   f_lat_grid(1:3) = [np, np, 1]
   quad_face(1, 1:4) = [0, 1, 1, 0]
   quad_face(2, 1:4) = [1, 1, 5, 5]
   quad_face(3, 1:4) = [0, 0, 5, 5]
   quad_lat_grid(1, 1:4) = [np-1, np-1, np, np]
   quad_lat_grid(2, 1:4) = [np,   np,   1,  1 ]
   quad_lat_grid(3, 1:4) = [np,   np,   1,  2 ]
   quad_lon_grid(1, 1:4) = [np,   1,  1,  np]
   quad_lon_grid(2, 1:4) = [1,    2,  2,  1 ]
   quad_lon_grid(3, 1:4) = [np-1, np, 1,  1 ]
elseif(corner == 7) then
   f_face(1:3) =     [1,  2,  5 ]
   f_lon_grid(1:3) = [np, 1,  np]
   f_lat_grid(1:3) = [np, np, 1 ]
   quad_face(1, 1:4) = [1, 2, 2, 1]
   quad_face(2, 1:4) = [2, 2, 5, 5]
   quad_face(3, 1:4) = [1, 1, 5, 5]
   quad_lat_grid(1, 1:4) = [np-1, np-1, np, np]
   quad_lat_grid(2, 1:4) = [np,   np,   1,   2]
   quad_lat_grid(3, 1:4) = [np,   np,   1,   1]
   quad_lon_grid(1, 1:4) = [np,   1,  1,  np  ]
   quad_lon_grid(2, 1:4) = [1,    2,  np, np  ]
   quad_lon_grid(3, 1:4) = [np-1, np, np, np-1]
elseif(corner == 8) then
   f_face(1:3) =     [2,  3,  5 ]
   f_lon_grid(1:3) = [np, 1,  np]
   f_lat_grid(1:3) = [np, np, np]
   quad_face(1, 1:4) = [2, 3, 3, 2]
   quad_face(2, 1:4) = [3, 3, 5, 5]
   quad_face(3, 1:4) = [2, 2, 5, 5];
   quad_lat_grid(1, 1:4) = [np-1, np-1, np, np  ]
   quad_lat_grid(2, 1:4) = [np,   np,   np, np  ]
   quad_lat_grid(3, 1:4) = [np,   np,   np, np-1]
   quad_lon_grid(1, 1:4) = [np,   1,  1,    np]
   quad_lon_grid(2, 1:4) = [1,    2,  np-1, np]
   quad_lon_grid(3, 1:4) = [np-1, np, np,   np]
endif

! Load up the array for the point
pxyz = lat_lon_to_xyz(lat, lon)

! Get lats and lons of the triangle vertices
do i = 1, 3
   call grid_to_lat_lon(f_face(i), f_lat_grid(i), f_lon_grid(i), np, grid_pt_lat, grid_pt_lon)
   ! Convert to x, y, z coords to check for whether points are in tris/quads
   qxyz(i, 1:3) = lat_lon_to_xyz(grid_pt_lat, grid_pt_lon)
enddo

! See if the point is in the triangle; if so, all is good
if(is_point_in_triangle(qxyz(1, :), qxyz(2, :), qxyz(3, :), pxyz)) return

! If it's not in the triangle, have to check the three adjacent quads at the corner
num_bound_points = 4

do quad = 1, 3
   ! Compute lat and lon for a quad
   do i = 1, 4
         call grid_to_lat_lon(quad_face(quad, i), quad_lat_grid(quad, i), quad_lon_grid(quad, i), &
            np, grid_pt_lat, grid_pt_lon)
      ! Convert to x, y, z coords to check for whether points are in tris/quads
      qxyz(i, 1:3) = lat_lon_to_xyz(grid_pt_lat, grid_pt_lon)
   enddo

   ! See if the point is inside this quad
   if(is_point_in_quad(qxyz, pxyz)) then
      f_face = quad_face(quad, 1:4)
      f_lat_grid = quad_lat_grid(quad, 1:4)
      f_lon_grid = quad_lon_grid(quad, 1:4)
      return
   endif
enddo

! Falling of the end is not happy
!!!fprintf('UNEXPECTED FAILURE IN GET_CORNERS.M\n');
!!!stop

end subroutine get_corners

!-----------------------------------------------------------------------

! Given the latitude and longitude of a point, returns the face, array indices, latitude
! and longitude of the bounding three or four grid points along the number of points
! (3 triangle; 4 quad). np is the number of grid points across each face of the cube sphere.

subroutine get_bounding_box(lat, lon, np, &
   grid_face, grid_lat_ind, grid_lon_ind, grid_pt_lat, grid_pt_lon, num_bound_points)

real(r8), intent(in)  :: lat, lon
integer,  intent(in)  :: np
real(r8), intent(out) :: grid_pt_lat(4), grid_pt_lon(4)
integer,  intent(out) :: grid_face(4), grid_lat_ind(4), grid_lon_ind(4), num_bound_points

real(r8) :: cube_side, del, half_del, len(2), qxyz(4, 3), pxyz(3)
integer  :: face, low_grid(2), hi_grid(2), i, my_pt, corner_index
integer  :: lat_grid(4), lon_grid(4), face1_pts(2), face2_pts(2), face1_count, face2_count
logical  :: on_edge, edge, corner

! Cube side is divided into np-1 interior intervals of width 2sqrt(1/3) / np and
! two exterior intervals with half the width, sqrt(1/3) / np 
cube_side = 2.0_r8 * sqrt(1.0_r8 / 3.0_r8) 
del = cube_side / np 
half_del = del / 2.0_r8

! Get the face and the length along the two imbedded cube faces for the point
call get_face(lat, lon, face, len);

! Figure out which interval this is in along each cube face; This gives 0 to np grid indices
low_grid(1) = floor((len(1) + half_del) / del)
low_grid(2) = floor((len(2) + half_del) / del)
hi_grid = low_grid + 1

! Get the indices for the lat and lon directions: Points go counterclockwise starting from lower left 
! For now assume this is a quad, but will correct below if it is a triangle
lat_grid(1) = low_grid(2); lat_grid(2) = hi_grid(2); lat_grid(3) = lat_grid(2); lat_grid(4) = lat_grid(1)
lon_grid(1) = low_grid(1); lon_grid(2) = lon_grid(1); lon_grid(3) = hi_grid(1); lon_grid(4) = lon_grid(3)

! If points are on the edge map to adjacent faces
on_edge = .false.
do i = 1, 4
   call fix_face(face, lat_grid(i), lon_grid(i), np, &
      grid_face(i), grid_lat_ind(i), grid_lon_ind(i), edge, corner)
   ! If any point is on an edge, on_edge is true
   if(edge) on_edge = .true.
   if(corner) then
      corner_index = i
      exit
   endif
enddo

! If it's at a corner, need to find the triangles in a different fashion
! It is possible that the point initially looks like it is in a corner due to the fact 
! that the edges of the grid are on great circles from the corresponding faces, but the
! grid point at the edge are not connected by these great circles.
if(corner) then
   call get_corners(face, lat_grid(corner_index), lon_grid(corner_index), lat, lon, np, &
      grid_face, grid_lat_ind, grid_lon_ind, num_bound_points)
   if(num_bound_points == 4) corner = .false.
else
   ! If not initially at a corner it's definitely in a quad
   num_bound_points = 4
endif

! Compute the lat and lon corresponding to these point
do i = 1, num_bound_points
   call grid_to_lat_lon(grid_face(i), grid_lat_ind(i), grid_lon_ind(i), np, &
      grid_pt_lat(i), grid_pt_lon(i))
enddo

! Make on_edge true only if we are on an edge but not at a corner
on_edge = (on_edge .and.  .not. corner)

if(on_edge) then
   ! If this is an edge, may need to revise box selection
   ! See if the point is in the box (approximately)
   ! Load up the arrays for the vertex points
   do i = 1, num_bound_points   
      ! Convert to x, y, z coords to check for whether points are in tris/quads
      qxyz(i, 1:3) = lat_lon_to_xyz(grid_pt_lat(i), grid_pt_lon(i))
   enddo

   ! Convert point to xyz
   pxyz = lat_lon_to_xyz(lat, lon);

   if(.not. is_point_in_quad(qxyz, pxyz)) then
      ! Not in this box, need to move 'equatorward'
      ! Find indices (from 1 to 4) of points on the same face
      face1_pts(1:2) = 0; face2_pts(1:2) = 0;
      face1_count = 0;    face2_count = 0;
      do i = 1, 4
         if(grid_face(i) == grid_face(1)) then
            face1_count = face1_count + 1;    face1_pts(face1_count) = i
         else
            face2_count = face2_count + 1;    face2_pts(face2_count) = i
         endif
      enddo

      ! First process points of the first face
      ! Are the latitudes or the longitudes on the edge
      if(grid_lon_ind(face1_pts(1)) == grid_lon_ind(face1_pts(2))) then
         ! Adjust the face1 latitudes
         do i = 1, 2
            my_pt = face1_pts(i)
            if(grid_lat_ind(my_pt) > np/2) then
               grid_lat_ind(my_pt) = grid_lat_ind(my_pt) - 1
            else
               grid_lat_ind(my_pt) = grid_lat_ind(my_pt) + 1
            endif
         enddo
      else
         ! Adjust the face1 longitudes
         do i = 1, 2
            my_pt = face1_pts(i)
            if(grid_lon_ind(my_pt) > np/2) then
               grid_lon_ind(my_pt) = grid_lon_ind(my_pt) - 1
            else
               grid_lon_ind(my_pt) = grid_lon_ind(my_pt) + 1
            endif
         enddo
      endif

      ! Do the same thing for face2 
      ! Are the latitudes or the longitudes on the edge
      if(grid_lon_ind(face2_pts(1)) == grid_lon_ind(face2_pts(2))) then
         ! Adjust the face1 latitudes
         do i = 1, 2
            my_pt = face2_pts(i);
            if(grid_lat_ind(my_pt) > np/2) then
               grid_lat_ind(my_pt) = grid_lat_ind(my_pt) - 1
            else
               grid_lat_ind(my_pt) = grid_lat_ind(my_pt) + 1
            endif
         enddo
      else
         ! Adjust the face1 longitudes
         do i = 1, 2
            my_pt = face2_pts(i);
            if(grid_lon_ind(my_pt) > np/2) then
               grid_lon_ind(my_pt) = grid_lon_ind(my_pt) - 1
            else
               grid_lon_ind(my_pt) = grid_lon_ind(my_pt) + 1
            endif
         enddo
      endif

      ! Compute the lat and lon corresponding to these point
      do i = 1, num_bound_points
         call grid_to_lat_lon(grid_face(i), grid_lat_ind(i), grid_lon_ind(i), np, &
            grid_pt_lat(i), grid_pt_lon(i))
      enddo
   endif

endif

end subroutine get_bounding_box

!-----------------------------------------------------------------------

subroutine test_grid_box

integer  :: np, i, j, num_bound_points
integer  :: grid_face(4), grid_lat_ind(4), grid_lon_ind(4)
real(r8) :: pt_lon_d, pt_lat_d, pt_lon, pt_lat 
real(r8) :: qxyz(4, 3), pxyz(3), grid_pt_lat(4), grid_pt_lon(4)
logical  :: inside

type(location_type) :: location
integer             :: qty, lon_count, lat_count, my_face, my_level, my_qty, my_lon_ind, my_lat_ind
integer(i8)         :: state_index, state_index2
real(r8)            :: lon_lat_hgt(3), my_lat, my_lon

! Parameter that has the number of model grid points along each dimension 
! This does not include halos; the points are offset from the boundaries
np = 18

! Open grid files for each face for comparison
!!!for face = 0:5
   !!!! Generate the grid file name
   !!!fname = strcat('grid_g000', num2str(face), '.nc');
   !!!xt = ncread(fname, 'Longitude');
   !!!glon(face + 1, :, :) = squeeze(xt(1, 3:end-2, 3:end-2));
   !!!yt = ncread(fname, 'Latitude');
   !!!glat(face + 1, :, :) = squeeze(yt(1, 3:end-2, 3:end-2));
!!!end

! Test points for the following:
! 1. Does the bounding box found contain the observed point?
! 2. Are the computed vertex latitudes and longitudes the same as those in the Aether grid files?

do lon_count = 0, 36000
   pt_lon_d = lon_count / 100.0_r8
   write(*, *) pt_lon_d
   do lat_count = -9000, 9000
   pt_lat_d = lat_count / 100.0_r8


! Convert to radians
pt_lon = DEG2RAD * pt_lon_d
pt_lat = DEG2RAD * pt_lat_d

! Get the x, y, z coords for this point
pxyz = lat_lon_to_xyz(pt_lat, pt_lon);

call get_bounding_box(pt_lat, pt_lon, np, &
   grid_face, grid_lat_ind, grid_lon_ind, grid_pt_lat, grid_pt_lon, num_bound_points)

do i = 1, num_bound_points
   ! Convert to x, y, z coords to check for whether points are in tris/quads
   qxyz(i, 1:3) = lat_lon_to_xyz(grid_pt_lat(i), grid_pt_lon(i));
enddo


! Get the latitude and longitude of the bounding points from get_state_meta_data as a confirmation test
do i = 1, num_bound_points
   state_index = get_state_index(grid_face(i), grid_lat_ind(i), grid_lon_ind(i), &
      1, 1, np, 44, 2)
   call get_state_meta_data(state_index, location, qty)
   lon_lat_hgt = get_location(location)       
   ! Deal with Aether file round off
   if(abs(lon_lat_hgt(1) - 360.0_r8) < 0.0001) lon_lat_hgt(1) = 0.0_r8 
   if(abs(RAD2DEG*grid_pt_lat(i) - lon_lat_hgt(2)) > 0.0001_r8 .or. &
      abs(RAD2DEG*grid_pt_lon(i) - lon_lat_hgt(1)) > 0.0001_r8) then
      write(*, *) grid_pt_lat(i), grid_pt_lon(i), lon_lat_hgt(2), lon_lat_hgt(1)
      stop
   endif
enddo

if(num_bound_points == 3) then
   ! See if the point is inside a local approximately tangent triangle
   inside = is_point_in_triangle(qxyz(1, :), qxyz(2, :), qxyz(3, :), pxyz);
else
   ! Or quadrilateral
   inside = is_point_in_quad(qxyz, pxyz);
endif

if(.not. inside) then
   write(*, *) 'ERROR: inside is false'
   write(*, *)pt_lat, pt_lon, num_bound_points
   stop
endif

enddo
enddo

!-------------------
! Block that loops through all state variables and confirms that the algorithms for mapping
! state vector (face/lon/lat) indices and get_state_meta_data correcty match up.
do my_qty = 1, 2
   do my_level = 1, 44
      do my_face = 0, 5
         do my_lat_ind = 1, np
            do my_lon_ind = 1, np
               state_index = get_state_index(my_face, my_lat_ind, my_lon_ind, &
                  my_level, my_qty, np, 44, 2)
               call get_state_meta_data(state_index, location, qty)
               lon_lat_hgt = get_location(location)       

               ! Want to compare the lat lon directly from code to that from get_state_meta_data
               call grid_to_lat_lon(my_face, my_lat_ind, my_lon_ind, np, my_lat, my_lon)
                
               ! ROUNDOFF FROM AETHER, NEED TO FIX somewhere.
               if(abs(lon_lat_hgt(1) - 360.0_r8) < 0.0001) lon_lat_hgt(1) = 0.0_r8 

               ! Check that things are consistent            
               if(abs(RAD2DEG*my_lat - lon_lat_hgt(2)) > 0.0001_r8 .or. &
                  abs(RAD2DEG*my_lon - lon_lat_hgt(1)) > 0.0001_r8) then
                  write(*, *) 'ERROR: grid points not appropriately mapping'
                  write(*, *) my_face, my_qty, my_level, my_lat_ind, my_lon_ind
                  stop
               endif
   
            enddo
         enddo
      enddo
   enddo
enddo
!-------------------------------

end subroutine test_grid_box

!-----------------------------------------------------------------------

function is_point_in_quad(v, p)

logical              :: is_point_in_quad
real(r8), intent(in) :: v(4, 3), p(3)

logical :: inside_t(4)

! See if the point is inside this quad; it's inside if it's in one or more contained triangles
inside_t(1) = is_point_in_triangle(v(1, :), v(2, :), v(3, :), p)
inside_t(2) = is_point_in_triangle(v(1, :), v(2, :), v(4, :), p)
inside_t(3) = is_point_in_triangle(v(1, :), v(3, :), v(4, :), p)
inside_t(4) = is_point_in_triangle(v(2, :), v(3, :), v(4, :), p)

is_point_in_quad = any(inside_t)

end function is_point_in_quad
   

!-----------------------------------------------------------------------

! Computes Herons formula to get area of triangle from lenghts of sides
! Super accuracy is not needed in the area calculation here

function heron(a, b, c)

real(r8)             :: heron
real(r8), intent(in) :: a, b, c

real(r8) :: s, arg

s = (a + b + c) /2
arg = (s * (s - a) * (s - b) * (s - c))

! Make sure we don't roundoff to a negative
if(arg <= 0.0_r8) then
   heron = 0.0_r8
else
   heron = sqrt(arg)
endif

end function heron

!-----------------------------------------------------------------------

function get_state_index(face, lat_ind, lon_ind, lev_ind, var_ind, np, n_lev, n_var)

integer             :: get_state_index
integer, intent(in) :: face, lat_ind, lon_ind, lev_ind, var_ind, np, n_lev, n_var

! Given the cube face, latitude (first) index, longitude (second) index on the face,
! the level index and the variable index, returns the state index for use
! by get_state_meta_data. Needs to know (are these in global storage) the number
! of lat and lon points across the face (np), the number of levels (n_lev) and
! the number of variables (n_var).

! This function makes the explicit assumption that the state is mapped to 
! the state index in the following fashion:
! From fastest to slowest varying index:
! 1. longitude index : 1 to np
! 2. latitude index  : 1 to np
! 3. face            : 0 to 5 consistent with Aether defs
! 4. level index     : 1 to n_lev
! 4. variable index  : 1 to number of variables

integer :: temp

get_state_index = lon_ind + np * ((lat_ind -1) + &
   np * (face + 6 * ((lev_ind - 1) + n_lev * (var_ind - 1))))

!!!temp = lon_ind + np *(lat_ind - 1) + np*np*face + &
   !!!np*np*6 * (lev_ind - 1) + np*np*6*n_lev * (var_ind - 1)
!!!if(get_state_index .ne. temp) then
   !!!write(*, *) 'get_state_index error ', get_state_index, temp
   !!!write(*, *) face, lat_ind, lon_ind, lev_ind, var_ind
   !!!stop
!!!endif

end function get_state_index

!-----------------------------------------------------------------------
! interpolate in the vertical between 2 arrays of items.

! vert_fracts: 0 is 100% of the first level and 
!              1 is 100% of the second level

subroutine vert_interp(nitems, levs1, levs2, vert_fract, out_vals)

   integer,  intent(in)  :: nitems
   real(r8), intent(in)  :: levs1(nitems)
   real(r8), intent(in)  :: levs2(nitems)
   real(r8), intent(in)  :: vert_fract
   real(r8), intent(out) :: out_vals(nitems)
   
   out_vals(:) = (levs1(:) * (1.0_r8 - vert_fract)) + &
                 (levs2(:) *           vert_fract )
   
end subroutine vert_interp


function barycentric_average(nitems, weights, vertex_temp_values) result (averaged_values)

   integer, intent(in)                                        :: nitems
   real(r8), dimension(nvertex_neighbors), intent(in)         :: weights
   real(r8), dimension(nvertex_neighbors, nitems), intent(in) :: vertex_temp_values

   real(r8), dimension(nitems)                :: averaged_values

   integer :: iweight, iitem

   averaged_values(:) = 0

   do iitem = 1, nitems
      do iweight = 1, nvertex_neighbors
         averaged_values(iitem) = averaged_values(iitem) + weights(iweight)*vertex_temp_values(iweight, iitem)
      end do
   end do

end function barycentric_average

!===================================================================
! End of model_mod
!===================================================================
end module model_mod
