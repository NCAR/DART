! DART software - Copyright UCAR. This open source software is provided
! by UCAR, "as is", without charge, subject to all terms of use at
! http://www.image.ucar.edu/DAReS/DART/DART_download
!

module model_mod

! This is a template showing the interfaces required for a model to be compliant
! with the DART data assimilation infrastructure. Do not change the arguments
! for the public routines.

use           netcdf

use        types_mod, only : r8, i8, MISSING_R8, vtablenamelength, DEG2RAD, radius => earth_radius

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
! For applications in which only perfect model experiments
! with identity observations (i.e. only the value of a particular
! state variable is observed), this can be a NULL INTERFACE.

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

! Begin variables for horiziontal interpolation

! End variables for horizontal interpolation
cyclic = .true.

lon_lat_alt = get_location(location)
which_vertical = nint(query_location(location))

! See if the state contains the obs quantity
varid = get_varid_from_kind(dom_id, qty)

if (varid > 0) then
    istatus = 0
else
    istatus = UNKNOWN_OBS_QTY_ERROR_CODE
endif

if ( .not. module_initialized ) call static_init_model

expected_obs(:) = MISSING_R8
istatus(:) = 1

! Find the vertical levels
if ( which_vertical == VERTISHEIGHT ) then
   call find_enclosing_indices(ncenter_altitudes, center_altitude, lon_lat_alt(3), below_index, &
                               above_index, fraction, enclosing_status)
   if (enclosing_status /= 0) then
      istatus(:) = INVALID_ALTITUDE_VAL_ERROR_CODE
   end if
else
   istatus(:) = INVALID_VERT_COORD_ERROR_CODE
   write(error_string_1, *) 'unsupported vertical type: ', which_vertical
   call error_handler(E_ERR, 'model_interpolate', error_string_1, source)
end if

! If the vertical location is acceptable, then do the horizontal interpolation
if (istatus(1) == 1) then

   ! Find the enclosing triangle or quad
   inside = .false.
   do icolumn = 1, nvertex_columns
      do ineighbor = 1, nvertex_neighbors
         call latlon_to_xyz(center_latitude(vertex_neighbor_indices(icolumn, ineighbor)), &
                           center_longitude(vertex_neighbor_indices(icolumn, ineighbor)), &
                           vertex_xyz(ineighbor, 1), &
                           vertex_xyz(ineighbor, 2), &
                           vertex_xyz(ineighbor, 3))
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

end if

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
   call get_3d_weights(r, t1, t2, t3, lat, lon, weights)
   
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

subroutine get_3d_weights(p, v1, v2, v3, lat, lon, weights)

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
      call get_barycentric_weights(p(1), p(2), cxs, cys, weights)
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
      call get_barycentric_weights(p(2), p(3), cxs, cys, weights)
      return
   endif
   
   ! last option, nearest 90 or 270 in lon:
   cxs(1) = v1(1)
   cxs(2) = v2(1)
   cxs(3) = v3(1)
   cys(1) = v1(3)
   cys(2) = v2(3)
   cys(3) = v3(3)
   call get_barycentric_weights(p(1), p(3), cxs, cys, weights)
   
end subroutine get_3d_weights

subroutine get_barycentric_weights(x, y, cxs, cys, weights)

   ! Computes the barycentric weights for a 2d interpolation point
   ! (x,y) in a 2d triangle with the given (cxs,cys) corners.
   
   real(r8), intent(in)  :: x, y, cxs(3), cys(3)
   real(r8), intent(out) :: weights(3)
   
   real(r8) :: denom
   
   ! Get denominator
   denom = (cys(2) - cys(3)) * (cxs(1) - cxs(3)) + &
      (cxs(3) - cxs(2)) * (cys(1) - cys(3))
   
   weights(1) = ((cys(2) - cys(3)) * (x - cxs(3)) + &
      (cxs(3) - cxs(2)) * (y - cys(3))) / denom
   
   weights(2) = ((cys(3) - cys(1)) * (x - cxs(3)) + &
      (cxs(1) - cxs(3)) * (y - cys(3))) / denom
   
   weights(3) = 1.0_r8 - weights(1) - weights(2)
   
   if (any(abs(weights) < roundoff)) then
      where (abs(weights) < roundoff) weights = 0.0_r8
      where (abs(1.0_r8 - abs(weights)) < roundoff) weights = 1.0_r8
   endif
   
end subroutine get_barycentric_weights

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
