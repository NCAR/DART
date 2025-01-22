! DART software - Copyright UCAR. This open source software is provided
! by UCAR, "as is", without charge, subject to all terms of use at
! http://www.image.ucar.edu/DAReS/DART/DART_download
!

module model_mod

use types_mod, only: r8, i8, deg2rad, missing_r8, missing_i, &
                     ps0, earth_radius, digits12, vtablenamelength, &
                     gas_constant, gas_constant_v, gravity, pi

use time_manager_mod, only: time_type, set_time, set_calendar_type, GREGORIAN, &
                            set_date, get_date

use location_mod, only: location_type, get_location, set_location, &
                        query_location, VERTISUNDEF, VERTISSURFACE, &
                        VERTISLEVEL, VERTISPRESSURE, VERTISHEIGHT, &
                        VERTISSCALEHEIGHT, vertical_localization_on, &
                        set_vertical_localization_coord, &
                        get_close_type, get_dist, is_vertical, &
                        loc_get_close_obs => get_close_obs, &
                        loc_get_close_state => get_close_state, &
                        get_maxdist, set_vertical
! VERTISUNDEF    = -2  ! has no specific vertical location (undefined) \ VERTISSURFACE     = -1  ! surface value  \ VERTISLEVEL       =  1  ! by level  \  VERTISPRESSURE    =  2  ! by pressure (in pascals) \ VERTISHEIGHT      =  3  ! by height (in meters) \ VERTISSCALEHEIGHT =  4  ! by scale height (unitless)
use utilities_mod, only: file_exist, open_file, close_file, &
                         register_module, error_handler, E_ERR, E_WARN, &
                         E_MSG, nmlfileunit, &
                         find_namelist_in_file, check_namelist_read, &
                         find_textfile_dims, file_to_text, &
                         do_nml_file, do_nml_term, scalar, &
                         find_enclosing_indices, to_upper, &
                         string_to_logical, string_to_real

use netcdf_utilities_mod, only: nc_get_variable, nc_get_variable_size, nc_create_file, &
                                nc_add_attribute_to_variable, &
                                nc_define_integer_variable, nc_define_double_variable, &
                                nc_define_real_variable, &
                                nc_define_real_scalar, &
                                nc_add_global_creation_time, &
                                nc_add_global_attribute, &
                                nc_define_dimension, nc_put_variable, &
                                nc_synchronize_file, nc_end_define_mode, &
                                nc_begin_define_mode, nc_open_file_readonly, &
                                nc_close_file, nc_variable_exists, nc_get_global_attribute, &
                                nc_get_dimension_size, nc_check

use mpi_utilities_mod, only: my_task_id, task_count, all_reduce_min_max

use random_seq_mod, only: random_seq_type, init_random_seq, random_gaussian

use obs_kind_mod, only: QTY_U_WIND_COMPONENT, QTY_V_WIND_COMPONENT, &
                        QTY_SURFACE_PRESSURE, QTY_TEMPERATURE, &
                        QTY_SPECIFIC_HUMIDITY, QTY_SURFACE_ELEVATION, &
                        QTY_PRESSURE, QTY_2M_TEMPERATURE, &
                        QTY_VAPOR_MIXING_RATIO, QTY_10M_U_WIND_COMPONENT, &
                        QTY_10M_V_WIND_COMPONENT, QTY_GEOPOTENTIAL_HEIGHT, &
                        QTY_LANDMASK, QTY_VERTLEVEL, QTY_GEOMETRIC_HEIGHT, &
                        QTY_SURFACE_TYPE, &
                        get_index_for_quantity, get_num_quantities, &
                        get_name_for_quantity, get_quantity_for_type_of_obs

use ensemble_manager_mod, only: ensemble_type, get_my_num_vars, get_my_vars

use distributed_state_mod, only: get_state

! These routines are passed through from default_model_mod.
! To write model specific versions of these routines
! remove the routine from this use statement and add your code to
! this the file.
! models/utilities/default_model_mod.f90
use default_model_mod, only: adv_1step, nc_write_model_vars, &
                             init_conditions => fail_init_conditions, &
                             init_time => fail_init_time

! assimilation_code/modules/io/state_structure_mod.f90
use state_structure_mod, only: add_domain, get_model_variable_indices, &
                               state_structure_info, get_domain_size, &
                               get_index_start, get_index_end, &
                               get_dart_vector_index, get_varid_from_kind, &
                               get_num_dims


use quad_utils_mod, only: quad_interp_handle, init_quad_interp, &
                          set_quad_coords, finalize_quad_interp, &
                          quad_lon_lat_locate, quad_lon_lat_evaluate, &
                          GRID_QUAD_IRREG_SPACED_REGULAR, &
                          QUAD_LOCATED_CELL_CENTERS

use netcdf

implicit none
private

!-----
! DART requires 18 specific public interfaces from model_mod.f90

! routines with code in this module
public ::  get_model_size, &  ! x same as in cam-se now, might need to rewrite for pangu, which doesn't call static_init_model          get_state_meta_data,    &
          get_state_meta_data, & ! x
          shortest_time_between_assimilations, &  ! x
          static_init_model, & ! x
          model_interpolate, & ! x
          nc_write_model_atts, & ! x
          pert_model_copies, & ! x
          get_close_obs, & ! x
          get_close_state, & ! x
          convert_vertical_obs, & ! x
          convert_vertical_state, & ! x
          read_model_time, & ! x
          write_model_time, & ! x
          end_model, & ! x
          domain_id
! pert_model_copies,             & !x let filter add gaussian noise to a single state vector to generate an ensemble

! required routines where the code is in another module
public ::  adv_1step, &
          init_time, &
          init_conditions, &
          nc_write_model_vars   !?  what is this?

!-----
! Here is the appropriate place for other users to make additional routines
!   contained within model_mod available for public use:
public ::   fill_default_state_table, &
          get_number_of_state_variables, &
          height_diff_check

! public parameters
public :: MAX_STATE_VARIABLES, &
          num_state_table_columns, &
          num_bounds_table_columns

! types
public :: model_static_data_for_dart

character(len=256), parameter :: source = "model_mod.f90"

! miscellaneous
integer, parameter :: MAX_STATE_VARIABLES = 100
integer, parameter :: num_state_table_columns = 5
integer, parameter :: num_bounds_table_columns = 4

! model_nml namelist variables and default values
integer :: domain_id = 1

character(len=256) :: input_template_filename = 'pginput_0001.nc'
character(len=32)  :: vertical_localization_coord = 'PRESSURE'
integer :: assimilation_period_seconds = 21600
integer :: debug_level = 0

! in converting to scale height for the vertical:
!  set this to .true. to compute the log of the pressure.
!  set this to .false. to additionally normalize by the surface
!    pressure at this location.  this option is backwards-compatible
!    with previous versions of this code.
logical :: no_normalization_of_scale_heights = .true.

! state_variables defines the contents of the state vector.
! each line of this input should have the form:
!
!    netcdf_variable_name, dart_quantity, clamp_min, clamp_max, update_variable
!
! all items must be strings (even if numerical values).
! for no clamping, use the string 'NA'
! to have the assimilation change the variable use 'UPDATE', else 'NO_UPDATE'
character(len=129) :: state_variables(MAX_STATE_VARIABLES* &
                                      num_state_table_columns) = ' '

! Just a global storage for output strings
character(len=512)  :: string1, string2, string3, errstring
logical, save      :: module_initialized = .false.

integer :: calendar_type = GREGORIAN
integer :: vertical_localization_type = VERTISPRESSURE
!integer :: vertical_localization_type = VERTISHEIGHT

! flag used to know if the vertical unit system has numbers
! that get larger as you move away from the earth's surface
! (e.g. height) or smaller (e.g. pressure)
logical :: higher_is_smaller

! things related to damping at the model top
logical  :: are_damping = .false.
real(r8) :: ramp_end         ! fixed top of ramp; the start (bottom) varies
logical  :: discarding_high_obs = .true.
real(kind=r8) :: no_assim_above_height = 20000_r8  ! approximate 50 hPa
real(kind=r8) :: no_assim_above_level = 13_r8
real(kind=r8) :: no_assim_above_pressure = 50_r8

! commonly used numbers that we'll set in static_init_model
real(r8) :: ref_model_top_pressure
real(r8) :: ref_surface_pressure
integer  :: ref_nlevels = 13

!> TODO: modify pert_model_copies subroutine
logical :: allow_perturbed_ics = .false.
type(time_type) :: assimilation_time_step

! Example Namelist
! Use the namelist for options to be set at runtime.
integer  :: time_step_days = 0
integer  :: time_step_seconds = 21600

namelist /model_nml/ input_template_filename, &
   vertical_localization_coord, &
   calendar_type, &
   state_variables, &
   assimilation_period_seconds

integer :: num_obs_kinds = 0

integer, parameter :: STAGGER_NONE = -1

type pangu_stagger
   integer, allocatable :: qty_stagger(:)
end type

type(pangu_stagger) :: grid_stagger

type(quad_interp_handle) :: interp_nonstaggered

real(kind=r8), PARAMETER    :: kappa = 2.0_r8/7.0_r8 ! gas_constant / cp
real(kind=r8), PARAMETER    :: ts0 = 300.0_r8        ! Base potential temperature for all levels.

TYPE model_static_data_for_dart

   integer  :: bt, sn, we
   real(r8) :: dt, dlon, dlat

   integer  :: domain_size
   integer  :: localization_coord
   real(r8), dimension(:, :), pointer :: hgt
   real(r8), dimension(:, :), pointer :: latitude
   real(r8), dimension(:, :), pointer :: longitude
   real(r8), dimension(:, :, :), pointer :: phb
   integer, dimension(:, :), pointer :: land

end type model_static_data_for_dart

type(model_static_data_for_dart) :: grid_data

contains

!------------------------------------------------------------------
!
! Called to do one time initialization of the model. As examples,
! might define information about the model size or model timestep.
! In models that require pre-computed static data, for instance
! spherical harmonic weights, these would also be computed here.

subroutine static_init_model()

   integer :: ncid
   integer  :: iunit, io

   logical, parameter    :: debug = .false.
   integer               :: nsize(5)
   integer               :: ref_nlevels

   character(len=*), parameter :: routine = 'static_init_model'
   
   module_initialized = .true.

   call find_namelist_in_file("input.nml", "model_nml", iunit)
   read (iunit, nml=model_nml, iostat=io)
   call check_namelist_read(iunit, io, "model_nml")
   ! Record the namelist values used for the run
   if (do_nml_file()) write (nmlfileunit, nml=model_nml)
   if (do_nml_term()) write (*, nml=model_nml)

   num_obs_kinds = get_num_quantities()

   ! set calendar type
   call set_calendar_type(calendar_type)
   ! read_grid_info
   if (file_exist(input_template_filename)) then  ! can be any .nc file converted from .npy
      ! contains necessary information
      call nc_check(nf90_open(input_template_filename, NF90_NOWRITE, ncid), &
                    'static_init_model', input_template_filename)
   else
      write (string1, *) 'Please put ', trim(input_template_filename), &
                         ' in the work directory.'
      call error_handler(E_ERR, routine, string1, source)
   end if
   if (debug) write (*, *) ' ncid is ', ncid

   !-------------------------------------------------------
   ! read dimensions
   ! returns numbers for bt, sn ,and we
   !-------------------------------------------------------
   call nc_get_variable_size(ncid, 'U', nsize(:))
   grid_data%we = nsize(1)
   grid_data%sn = nsize(2)
   grid_data%bt = nsize(3)

   !-------------------------------------------------------
   ! read model file attributes
   ! get grid_data%dlon, grid_data%dlat, grid_data%dt
   !-------------------------------------------------------
   call read_grid_file_attributes(ncid)

   !-------------------------------------------------------
   ! read static data
   ! get grid_data%longitude, grid_data%latitude, etc.
   !-------------------------------------------------------
   call read_model_static_data(ncid)

   ! init_globals equivalent
   !ref_nlevels = 13
   ref_model_top_pressure = 5000 ! Pa
   ref_surface_pressure = 100000 ! Pa

   !-------------------------------------------------------
   ! read the namelist &model_nml :: state_variables  and add_domain
   ! to set up what will be read into the model state vector
   call set_model_variable_info(input_template_filename, state_variables)

   call fill_model_stagger_info(grid_stagger)    ! doesn't seem necessary
   if (debug_level > 100) call state_structure_info(domain_id)

   ! convert from string in namelist to integer (e.g. VERTISxxx)
   ! and tell the dart code which vertical type we want to localize in.
   call set_vert_localization(vertical_localization_coord)

   call setup_interpolation() !grid is global

   !>TODO: set top limit where obs impacts are diminished to 0.
   !>TODO: discarding_high_obs

   ! set a flag based on the vertical localization coordinate selected
   call init_sign_of_vert_units()

   ! close data file, we have all we need
   call nc_check(nf90_close(ncid), 'static_init_model', 'close ', input_template_filename)

end subroutine static_init_model

!------------------------------------------------------------------
! Returns the number of items in the state vector as an integer.

function get_model_size()
   integer(i8) :: get_model_size
   integer     :: domain_id

   domain_id = 1
   if (.not. module_initialized) call static_init_model

   get_model_size = get_domain_size(domain_id)

end function get_model_size

!------------------------------------------------------------------
! Returns the smallest increment in time that the model is capable
! of advancing the state in a given implementation, or the shortest
! time you want the model to advance between assimilations.

function shortest_time_between_assimilations()

   type(time_type) :: shortest_time_between_assimilations
   integer :: model_dt, assim_dt

   ! We need to coordinate the desired assimilation window to be a
   ! multiple of the model time step (which has no precision past integer seconds).

   model_dt = nint(grid_data%dt)

   ! The integer arithmetic does its magic.
   assim_dt = (assimilation_period_seconds/model_dt)*model_dt

   shortest_time_between_assimilations = set_time(assim_dt)

end function shortest_time_between_assimilations

!--------------------------------------------------------------------
!> level 1 is at the lowest level, level N is the model top
!> our convention in this code is:  between levels a fraction of 0
!> is 100% level 1, and fraction of 1 is 100% level 2.

function check_good_levels(vert_value, valid_range, l1, l2, fract)
   real(r8), intent(in)  :: vert_value
   integer, intent(in)  :: valid_range
   integer, intent(out) :: l1
   integer, intent(out) :: l2
   real(r8), intent(out) :: fract
   logical               :: check_good_levels

   integer :: integer_level
   real(r8) :: fract_level

   ! be a pessimist, then you're never disappointed
   check_good_levels = .false.
   l1 = MISSING_I
   l2 = MISSING_I
   fract = MISSING_R8

   ! out of range checks
   if (vert_value < 1.0_r8 .or. vert_value > valid_range) return

   integer_level = floor(vert_value)
   fract_level = vert_value - integer_level

   ! cam levels start at the top so level 1 is
   ! the highest level and increases on the way down.

   !>might want to allow extrapolation - which means
   !>allowing out of range values here and handling
   !>them correctly in the calling and vert_interp() code.

   if (vert_value /= valid_range) then
      l1 = integer_level
      l2 = integer_level + 1
      fract = fract_level
   else
      ! equal to the largest level number
      l1 = integer_level - 1
      l2 = integer_level
      fract = 1.0_r8
   end if

   check_good_levels = .true.

end function check_good_levels

!-----------------------------------------------------------------------
!> Given an integer index into the state vector structure, returns the
!> associated location. A second intent(out) optional argument quantity
!> can be returned if the model has more than one type of field (for
!> instance temperature and zonal wind component).
!>
!> @param index_in the index into the DART state vector
!> @param location the location at that index
!> @param var_type_out the DART Quantity at that index
!>

subroutine get_state_meta_data(index_in, location, var_type_out)

   integer(i8), intent(in)  :: index_in
   type(location_type), intent(out) :: location
   integer, optional, intent(out) :: var_type_out

   integer     :: var_type, dart_type
   integer     :: iloc, jloc, vloc
   real(r8)    :: lon, lat, lev

   integer :: var_id, var_qty, nd

   real(r8) :: use_vert_val
   integer  :: use_vert_type
   character(len=*), parameter :: routine = 'get_state_meta_data'

   logical, parameter :: debug = .false.

   if (.not. module_initialized) call static_init_model

   ! from the dart index get the local variables indices
   call get_model_variable_indices(index_in, iloc, jloc, vloc, var_id=var_id, kind_index=var_qty)

   nd = get_num_dims(domain_id, var_id)

!                if(debug) write(*,*) ' ip, jp, kp for index ', iloc, jloc, vloc, index_in
!                if(debug) write(*,*) ' Var type: ',var_type

   ! first obtain lat/lon from (iloc,jloc)
   call get_model_horizontal_location(iloc, jloc, var_type, lon, lat)

   ! full 3d fields has lon/lat/level.
   ! 2d fields are either surface fields

   if (nd == 3) then
      use_vert_type = VERTISLEVEL
      use_vert_val = real(vloc, r8)
   else if (nd == 2) then
      if (is_surface_field(var_qty)) then
         use_vert_type = VERTISSURFACE
         use_vert_val = grid_data%hgt(iloc, jloc)
      else
         ! any 2d field not listed as a surface field (in is_surface_field() function)
         ! is assumed to be an integrated quantity with a vert type of VERTISUNDEF.
         use_vert_type = VERTISUNDEF
         use_vert_val = MISSING_R8
      end if
   else
      write (string1, *) 'state vector field not 2D or 3D and no code to handle other dimensionity'
      write (string2, *) 'dimensionality = ', nd, ' quantity type = ', trim(get_name_for_quantity(var_qty))
      call error_handler(E_ERR, routine, string1, source, text2=string2)
   end if

   if (debug) write (*, *) 'lon, lat, lev: ', lon, lat, use_vert_val
   location = set_location(lon, lat, use_vert_val, use_vert_type)

   ! return DART variable kind if requested
   if (present(var_type_out)) var_type_out = var_qty

end subroutine get_state_meta_data

!-----------------------------------------------------------------------
!>
!> Model interpolate will interpolate any DART state variable
!> to the given location.
!>
!> @param state_handle DART ensemble handle
!> @param ens_size DART ensemble size
!> @param location the location of interest
!> @param obs_qty the DART Quantity of interest
!> @param interp_vals the estimated value of the DART state at the location
!>          of interest (the interpolated value).
!> @param istatus interpolation status ... 0 == success, /=0 is a failure
!>
!> istatus = 2    asked to interpolate an unknown/unsupported quantity
!> istatus = 3    cannot locate horizontal quad
!> istatus = 4    cannot locate enclosing vertical levels
!> istatus = 5    cannot retrieve state vector values
!> istatus = 6    cannot get values at quad corners
!> istatus = 7    unused (error code available)
!> istatus = 8    cannot interpolate in the quad to get the values
!> istatus = 9    unused (error code available)
!> istatus = 10   cannot get vertical levels for an obs on pressure levels
!> istatus = 11   cannot get vertical levels for an obs on height levels
!> istatus = 12   cannot get values from obs quantity
!> istatus = 13   can not interpolate values of this quantity
!> istatus = 14   obs above user-defined assimilation top pressure
!> istatus = 15   can not get indices from given state vector index
!> istatus = 16   cannot do vertical interpolation for bottom layer
!> istatus = 17   cannot do vertical interpolation for top layer
!> istatus = 98   unknown error - shouldn't happen
!> istatus = 99   unknown error - shouldn't happen
!>

subroutine model_interpolate(state_handle, ens_size, location, obs_qty, interp_vals, istatus)

   type(ensemble_type), intent(in) :: state_handle
   integer, intent(in) :: ens_size
   type(location_type), intent(in) :: location
   integer, intent(in) :: obs_qty
   real(r8), intent(out) :: interp_vals(ens_size) !< array of interpolated values
   integer, intent(out) :: istatus(ens_size)

   character(len=*), parameter :: routine = 'model_interpolate:'

   ! local
   logical, parameter  :: debug = .false.
   integer             :: status1
   integer             :: status_array(ens_size)
   integer             :: four_lons(4), four_lats(4)
   integer             :: which_vert, varid
   real(r8)            :: lon_fract, lat_fract
   real(r8)            :: lon_lat_vert(3)
   real(r8)            :: quad_vals(ens_size, 4)
   type(quad_interp_handle) :: interp_handle

   integer  :: imem !< index varibles for loop
   if (.not. module_initialized) call static_init_model

   ! Initialize stuff
   istatus(:) = 99
   interp_vals(:) = missing_r8  !< array of obs_vals
   call ok_to_interpolate(obs_qty, varid, status1)

   if (status1 /= 0) then
      if (debug_level > 12) then
         write (string1, *) 'did not find observation quantity ', obs_qty, ' in the state vector'
         call error_handler(E_MSG, routine, string1, source)
      end if
      istatus(:) = status1   ! this quantity not in the state vector
      return
   end if

   ! Unravel location_type information
   lon_lat_vert = get_location(location)
   which_vert = nint(query_location(location))

   interp_handle = get_interp_handle()   ! all variables are not staggered
   call quad_lon_lat_locate(interp_handle, lon_lat_vert(1), lon_lat_vert(2), &
                            four_lons, four_lats, lon_fract, lat_fract, status1)

   if (status1 /= 0) then
      istatus(:) = 3  ! cannot locate enclosing horizontal quad
      return
   end if

!                         if (discarding_high_obs) then
!                            call obs_too_high(lon_lat_vert(3), which_vert, status1)
!                            if (status1 /= 0) then
!                               istatus(:) = status1
!                               return
!                            endif
!                         endif
   !*****************************************************************************
   call get_quad_vals(state_handle, ens_size, varid, obs_qty, four_lons, four_lats, &
                      lon_lat_vert, which_vert, quad_vals, status_array)
   if (any(status_array /= 0)) then
      istatus(:) = maxval(status_array)   ! cannot get the state values at the corners
      return
   end if

   do imem = 1, ens_size
      ! do the horizontal interpolation for each ensemble member.
      ! call one at a time to avoid creating temporary arrays
      call quad_lon_lat_evaluate(interp_handle, lon_fract, lat_fract, &
                                 quad_vals(imem, :), interp_vals(imem), status_array(imem))
   end do

   ! all interp values should be set by now.  set istatus
   istatus(:) = 0
end subroutine model_interpolate

subroutine get_quad_vals(state_handle, ens_size, varid, obs_qty, four_lons, four_lats, &
                         lon_lat_vert, which_vert, quad_vals, my_status)
   type(ensemble_type), intent(in) :: state_handle
   integer, intent(in) :: ens_size
   integer, intent(in) :: varid
   integer, intent(in) :: obs_qty
   integer, intent(in) :: four_lons(4)
   integer, intent(in) :: four_lats(4)
   real(r8), intent(in) :: lon_lat_vert(3)
   integer, intent(in) :: which_vert
   real(r8), intent(out) :: quad_vals(ens_size, 4) !< array of interpolated values
   integer, intent(out) :: my_status(ens_size)

   integer  :: icorner, numdims
   integer  :: level_one_array(ens_size)
   integer  :: four_levs1(ens_size, 4), four_levs2(ens_size, 4)
   real(r8) :: four_vert_fracts(ens_size, 4)

   character(len=*), parameter :: routine = 'get_quad_vals:'

   quad_vals(:, :) = MISSING_R8
   my_status(:) = 99

   ! need to consider the case for 2d vs 3d variables
   numdims = get_dims_from_qty(obs_qty, varid)
   ! now here potentially we have different results for different
   ! ensemble members.  the things that can vary are dimensioned by ens_size.

   if (numdims == 3) then

      ! build 4 columns to find vertical level numbers
      do icorner = 1, 4
         call find_vertical_levels(state_handle, ens_size, &
                                   four_lons(icorner), four_lats(icorner), lon_lat_vert(3), &
                                   which_vert, obs_qty, varid, &
                                   four_levs1(:, icorner), four_levs2(:, icorner), &
                                   four_vert_fracts(:, icorner), my_status)
         if (any(my_status /= 0)) return

      end do

      ! we have all the indices and fractions we could ever want.
      ! now get the data values at the bottom levels, the top levels,
      ! and do vertical interpolation to get the 4 values in the columns.
      ! the final horizontal interpolation will happen later.

      if (varid > 0) then

         call get_four_state_values(state_handle, ens_size, four_lons, four_lats, &
                                    four_levs1, four_levs2, four_vert_fracts, &
                                    varid, quad_vals, my_status)
      else ! get 3d special variables in another ways ( like QTY_PRESSURE )
         call get_four_nonstate_values(state_handle, ens_size, four_lons, four_lats, &
                                       four_levs1, four_levs2, four_vert_fracts, &
                                       obs_qty, quad_vals, my_status)
      end if

      if (any(my_status /= 0)) return

   else if (numdims == 2) then

      if (varid > 0) then
         level_one_array(:) = 1
         do icorner = 1, 4
            call get_values_from_varid(state_handle, ens_size, &
                                       four_lons(icorner), four_lats(icorner), &
                                       level_one_array, varid, quad_vals(:, icorner), my_status)
            if (any(my_status /= 0)) return

         end do

      else ! special 2d case
         do icorner = 1, 4
            call get_quad_values(ens_size, four_lons(icorner), four_lats(icorner), &
                                 obs_qty, obs_qty, quad_vals(:, icorner))
         end do
         ! apparently this can't fail
         my_status(:) = 0

      end if

   else
      write (string1, *) trim(get_name_for_quantity(obs_qty)), ' has dimension ', numdims
      call error_handler(E_ERR, routine, 'only supports 2D or 3D fields', &
                         source, text2=string1)
   end if

   ! when you get here, my_status() was set either by passing it to a
   ! subroutine, or setting it explicitly here.  if this routine returns
   ! the default value of 99 something went wrong in this logic.

end subroutine get_quad_vals

subroutine get_four_state_values(state_handle, ens_size, four_lons, four_lats, &
                                 four_levs1, four_levs2, four_vert_fracts, &
                                 varid, quad_vals, my_status)

   type(ensemble_type), intent(in) :: state_handle
   integer, intent(in) :: ens_size
   integer, intent(in) :: four_lons(4), four_lats(4)
   integer, intent(in) :: four_levs1(ens_size, 4), four_levs2(ens_size, 4)
   real(r8), intent(in) :: four_vert_fracts(ens_size, 4)
   integer, intent(in) :: varid
   real(r8), intent(out) :: quad_vals(ens_size, 4) !< array of interpolated values
   integer, intent(out) :: my_status(ens_size)

   integer  :: icorner
   real(r8) :: vals1(ens_size), vals2(ens_size)

   character(len=*), parameter :: routine = 'get_four_state_values:'

   do icorner = 1, 4
      call get_values_from_varid(state_handle, ens_size, &
                                 four_lons(icorner), four_lats(icorner), &
                                 four_levs1(:, icorner), varid, vals1, &
                                 my_status)

      if (any(my_status /= 0)) then
         my_status(:) = 16   ! cannot retrieve vals1 values
         return
      end if

      call get_values_from_varid(state_handle, ens_size, &
                                 four_lons(icorner), four_lats(icorner), &
                                 four_levs2(:, icorner), varid, vals2, my_status)
      if (any(my_status /= 0)) then
         my_status(:) = 17   ! cannot retrieve top values
         return
      end if

      call vert_interp(ens_size, vals1, vals2, four_vert_fracts(:, icorner), &
                       quad_vals(:, icorner))

   end do

end subroutine get_four_state_values

!-----------------------------------------------------------------------
!> interpolate in the vertical between 2 arrays of items.
!>
!> vert_fracts: 0 is 100% of the first level and
!>              1 is 100% of the second level

subroutine vert_interp(nitems, levs1, levs2, vert_fracts, out_vals)
   integer, intent(in)  :: nitems
   real(r8), intent(in)  :: levs1(nitems)
   real(r8), intent(in)  :: levs2(nitems)
   real(r8), intent(in)  :: vert_fracts(nitems)
   real(r8), intent(out) :: out_vals(nitems)

   out_vals(:) = (levs1(:)*(1.0_r8 - vert_fracts(:))) + &
                 (levs2(:)*vert_fracts(:))

end subroutine vert_interp

subroutine get_four_nonstate_values(state_handle, ens_size, four_lons, four_lats, &
                                    four_levs1, four_levs2, four_vert_fracts, &
                                    obs_qty, quad_vals, my_status)

   type(ensemble_type), intent(in) :: state_handle
   integer, intent(in) :: ens_size
   integer, intent(in) :: four_lons(4), four_lats(4)
   integer, intent(in) :: four_levs1(ens_size, 4), four_levs2(ens_size, 4)
   real(r8), intent(in) :: four_vert_fracts(ens_size, 4)
   integer, intent(in) :: obs_qty
   real(r8), intent(out) :: quad_vals(ens_size, 4) !< array of interpolated values
   integer, intent(out) :: my_status(ens_size)

   integer  :: icorner
   real(r8) :: vals1(ens_size), vals2(ens_size)

   character(len=*), parameter :: routine = 'get_four_nonstate_values:'

   do icorner = 1, 4
      call get_values_from_nonstate_fields(state_handle, ens_size, &
                                           four_lons(icorner), four_lats(icorner), &
                                           four_levs1(:, icorner), obs_qty, vals1, my_status)
      if (any(my_status /= 0)) then
         my_status(:) = 16   ! cannot retrieve vals1 values
         return
      end if

      call get_values_from_nonstate_fields(state_handle, ens_size, &
                                           four_lons(icorner), four_lats(icorner), &
                                           four_levs2(:, icorner), obs_qty, vals2, my_status)
      if (any(my_status /= 0)) then
         my_status(:) = 17   ! cannot retrieve top values
         return
      end if

      call vert_interp(ens_size, vals1, vals2, four_vert_fracts(:, icorner), &
                       quad_vals(:, icorner))

   end do

end subroutine get_four_nonstate_values

!-----------------------------------------------------------------------
!> figure out whether this is a 2d or 3d field based on the quantity.
!> if this field is in the state vector, use the state routines.
!> if it's not, there are cases for known other quantities we can
!> interpolate and return.  add any new non-state fields here.

function get_dims_from_qty(obs_quantity, var_id)
   integer, intent(in) :: obs_quantity
   integer, intent(in) :: var_id
   integer :: get_dims_from_qty
   character(len=*), parameter :: routine = 'get_dims_from_qty:'

   if (var_id > 0) then
      get_dims_from_qty = get_num_dims(domain_id, var_id)
   else
      select case (obs_quantity)
      case (QTY_SURFACE_ELEVATION)
         get_dims_from_qty = 2
         ! NCHEN provide the pressure levels for interpolations?
      case (QTY_PRESSURE, QTY_GEOMETRIC_HEIGHT)
         get_dims_from_qty = 3
      case default
         write (string1, *) 'we can not interpolate qty "', get_name_for_quantity(obs_quantity), &
            '" if the dimension is not known'
         call error_handler(E_ERR, routine, string1, source)
      end select
   end if

end function get_dims_from_qty

!-----------------------------------------------------------------------
!>
!>  This is for 2d special observations quantities not in the state
!>  All Pangu states are nonstaggered

subroutine get_quad_values(ens_size, lon_index, lat_index, obs_quantity, stagger_qty, vals)
   integer, intent(in) :: ens_size
   integer, intent(in) :: lon_index
   integer, intent(in) :: lat_index
   integer, intent(in) :: obs_quantity
   integer, optional, intent(in) :: stagger_qty
   real(r8), intent(out) :: vals(ens_size)
   character(len=*), parameter :: routine = 'get_quad_values'

   select case (obs_quantity)
   case (QTY_SURFACE_ELEVATION)

      vals = grid_data%hgt(lon_index, lat_index)

   case default
      write (string1, *) 'we can not interpolate qty', obs_quantity
      call error_handler(E_ERR, routine, string1, source)

   end select

end subroutine get_quad_values

!-----------------------------------------------------------------------
!> this routine converts the 3 index values and a quantity into a state vector
!> offset and gets the ensemble of state values for that offset.  this only
!> gets a single vertical location - if you need to get values which might
!> have different vertical locations in different ensemble members
!> see get_values_from_varid() below.

subroutine get_values_from_single_level(ens_handle, ens_size, qty, lon_index, lat_index, lev_index, &
                                        vals, my_status)
   type(ensemble_type), intent(in) :: ens_handle
   integer, intent(in) :: ens_size
   integer, intent(in) :: qty
   integer, intent(in) :: lon_index
   integer, intent(in) :: lat_index
   integer, intent(in) :: lev_index
   real(r8), intent(out) :: vals(ens_size)
   integer, intent(out) :: my_status
   character(len=*), parameter :: routine = 'get_values_from_single_level:'

   integer :: varid
   integer(i8) :: state_indx

   varid = get_varid_from_kind(domain_id, qty)
   if (varid < 0) then
      vals(:) = MISSING_R8
      my_status = 12
      return
   end if

   state_indx = get_dart_vector_index(lon_index, lat_index, lev_index, domain_id, varid)
   if (state_indx < 1 .or. state_indx > get_domain_size(domain_id)) then
      write (string1, *) 'state_index out of range: ', state_indx, ' not between ', 1, get_domain_size(domain_id)
      call error_handler(E_ERR, routine, string1, source, text2=string2, text3='should not happen')
   end if
   vals(:) = get_state(state_indx, ens_handle)

   my_status = 0

end subroutine get_values_from_single_level

!-----------------------------------------------------------------------
!> this routine takes care of getting the actual state values.  get_state()
!> communicates with other MPI tasks and can be expensive.
!>
!> all ensemble members have the same horizontal location, but different
!> ensemble members could have different vertical locations and
!> so be between different vertical layers.  this code tries to do the fewest
!> calls to get_state by only calling it for levels that are actually needed
!> and setting all members with those same levels in a single pass.
!>

subroutine get_values_from_varid(ens_handle, ens_size, lon_index, lat_index, lev_index, varid, &
                                 vals, my_status)
   type(ensemble_type), intent(in)  :: ens_handle
   integer, intent(in)  :: ens_size
   integer, intent(in)  :: lon_index
   integer, intent(in)  :: lat_index
   integer, intent(in)  :: lev_index(ens_size)
   integer, intent(in)  :: varid
   real(r8), intent(out) :: vals(ens_size)
   integer, intent(out) :: my_status(ens_size)

   integer(i8) :: state_indx
   integer  :: i, j
   real(r8) :: temp_vals(ens_size)
   logical  :: member_done(ens_size)
   character(len=*), parameter :: routine = 'get_values_from_varid:'

   ! as we get the values for each ensemble member, we set the 'done' flag
   ! and a good return code.
   my_status(:) = 12
   member_done(:) = .false.

   ! start with lev_index(1).  get the vals into a temp var.
   ! run through 2-N. any other member that has the same level
   ! set the outgoing values.  keep a separate flag for which
   ! member(s) have been done.  skip to the next undone member
   ! and get the state for that level.  repeat until all levels done.

   do i = 1, ens_size

      if (member_done(i)) cycle

      state_indx = get_dart_vector_index(lon_index, lat_index, lev_index(i), domain_id, varid)

      if (state_indx < 0) then
         write (string1, *) 'Should not happen: could not find dart state index from '
         write (string2, *) 'lon, lat, and lev index :', lon_index, lat_index, lev_index
         call error_handler(E_ERR, routine, string1, source, text2=string2)
         return
      end if

      temp_vals(:) = get_state(state_indx, ens_handle)    ! all the ensemble members for level (i)

      ! start at i, because my ensemble member is clearly at this level.
      ! then continue on to see if any other members are also at this level.
      do j = i, ens_size
         if (member_done(j)) cycle

         if (lev_index(j) == lev_index(i)) then
            vals(j) = temp_vals(j)
            member_done(j) = .true.
            my_status(j) = 0
         end if
      end do
   end do

end subroutine get_values_from_varid

!-----------------------------------------------------------------------
!> this is just for 3d fields

subroutine get_values_from_nonstate_fields(ens_handle, ens_size, lon_index, lat_index, &
                                           lev_index, obs_quantity, vals, my_status)
   type(ensemble_type), intent(in)  :: ens_handle
   integer, intent(in)  :: ens_size
   integer, intent(in)  :: lon_index
   integer, intent(in)  :: lat_index
   integer, intent(in)  :: lev_index(ens_size)
   integer, intent(in)  :: obs_quantity
   real(r8), intent(out) :: vals(ens_size)
   integer, intent(out) :: my_status(ens_size)

   integer  :: imember
   real(r8) :: vals_array(ref_nlevels, ens_size)
   character(len=*), parameter :: routine = 'get_values_from_nonstate_fields:'

   vals(:) = MISSING_R8
   my_status(:) = 99

   select case (obs_quantity)
   case (QTY_PRESSURE)

      call get_pressure_column(ens_size, vals_array)

      do imember = 1, ens_size
         vals(imember) = vals_array(lev_index(imember), imember)
      end do

   case (QTY_VERTLEVEL)
      vals(:) = lev_index(:)
      my_status(:) = 0

   case default
      write (string1, *) 'contact dart support. unexpected error for quantity ', obs_quantity
      call error_handler(E_MSG, routine, string1, source)

   end select

end subroutine get_values_from_nonstate_fields

!-----------------------------------------------------------------------
!> given a lon/lat index number, a quantity and a vertical value and type,
!> return which two levels these are between and the fraction across.
!>

subroutine find_vertical_levels(ens_handle, ens_size, lon_index, lat_index, vert_val, &
                                which_vert, obs_qty, var_id, levs1, levs2, vert_fracts, my_status)
   type(ensemble_type), intent(in)  :: ens_handle
   integer, intent(in)  :: ens_size
   integer, intent(in)  :: lon_index
   integer, intent(in)  :: lat_index
   real(r8), intent(in)  :: vert_val
   integer, intent(in)  :: which_vert
   integer, intent(in)  :: obs_qty
   integer, intent(in)  :: var_id
   integer, intent(out) :: levs1(ens_size)
   integer, intent(out) :: levs2(ens_size)
   real(r8), intent(out) :: vert_fracts(ens_size)
   integer, intent(out) :: my_status(ens_size)
   character(len=*), parameter :: routine = 'find_vertical_levels:'

   integer  :: l1, l2, imember, level_one, status1, k
   real(r8) :: fract1
   real(r8) :: surf_pressure(ens_size)
   real(r8) :: pressure_array(ref_nlevels, ens_size)
   real(r8) :: height_array(ref_nlevels, ens_size)

   ! assume the worst
   levs1(:) = MISSING_I
   levs2(:) = MISSING_I
   vert_fracts(:) = MISSING_R8
   my_status(:) = 98

   ! ref_nlevels is the number of vertical levels (midlayer points)

   level_one = 1

   select case (which_vert)

   case (VERTISPRESSURE)
      ! construct a pressure column here and find the model levels
      ! that enclose this value
      call get_pressure_column(ens_size, pressure_array)

      do imember = 1, ens_size
         call pressure_to_level(ref_nlevels, pressure_array(1:ref_nlevels, imember), vert_val, &
                                levs1(imember), levs2(imember), &
                                vert_fracts(imember), my_status(imember))

      end do

      if (debug_level > 100) then
         do k = 1, ens_size
            print *, 'ISPRESSURE levs1(k), levs2(k), vert_fracts(k), vert_val', &
               levs1(k), levs2(k), vert_fracts(k), vert_val
         end do
      end if

   case (VERTISHEIGHT)
      ! construct a height column here and find the model levels
      ! that enclose this value
      call pangu_height_levels(ens_handle, ens_size, lon_index, lat_index, ref_nlevels, obs_qty, &
                               height_array, my_status)

      !>@todo FIXME let successful members continue?
      if (any(my_status /= 0)) return

      if (debug_level > 400) then
         do k = 1, ref_nlevels
            print *, 'ISHEIGHT: ', k, height_array(k, 1)
         end do
      end if

      do imember = 1, ens_size
         call height_to_level(ref_nlevels, height_array(:, imember), vert_val, &
                              levs1(imember), levs2(imember), vert_fracts(imember), &
                              my_status(imember))
      end do

      !>@todo FIXME let successful members continue?
      if (any(my_status /= 0)) return

      if (debug_level > 100) then
         do k = 1, ens_size
            print *, 'ISHEIGHT ens#, levs1(#), levs2(#), vert_fracts(#), top/bot height(#)', &
               k, levs1(k), levs2(k), vert_fracts(k), height_array(levs2(k), k), height_array(levs1(k), k)
         end do
      end if

   case (VERTISLEVEL)
      ! this routine returns false if the level number is out of range.
      if (.not. check_good_levels(vert_val, ref_nlevels, l1, l2, fract1)) then
         my_status(:) = 8
         return
      end if

      ! because we're given a model level as input, all the ensemble
      ! members have the same outgoing values.
      levs1(:) = l1
      levs2(:) = l2
      vert_fracts(:) = fract1
      my_status(:) = 0

      if (debug_level > 100) then
         do k = 1, ens_size
            print *, 'ISLEVEL levs1(k), levs2(k), vert_fracts(k), vert_val', &
               levs1(k), levs2(k), vert_fracts(k), vert_val
         end do
      end if

      ! 2d fields
   case (VERTISUNDEF, VERTISSURFACE)
      if (get_dims_from_qty(obs_qty, var_id) == 2) then
         levs1(:) = ref_nlevels - 1
         levs2(:) = ref_nlevels
         vert_fracts(:) = 1.0_r8
         my_status(:) = 0
      else
         my_status(:) = 4 ! can not get vertical levels
      end if

   case default
      write (string1, *) 'unsupported vertical type: ', which_vert
      call error_handler(E_ERR, routine, string1, source)

   end select

   ! by this time someone has already set my_status(), good or bad.

end subroutine find_vertical_levels

!-----------------------------------------------------------------------
!> Compute the heights at pressure midpoints
!>
!> this version does all ensemble members at once.

subroutine pangu_height_levels(ens_handle, ens_size, lon_index, lat_index, nlevels, qty, height_array, my_status)
   type(ensemble_type), intent(in)  :: ens_handle
   integer, intent(in)  :: ens_size
   integer, intent(in)  :: lon_index
   integer, intent(in)  :: lat_index
   integer, intent(in)  :: nlevels
   integer, intent(in)  :: qty
   real(r8), intent(out) :: height_array(nlevels, ens_size)
   integer, intent(out) :: my_status(ens_size)

   integer  :: k, level_one, imember, status1
   real(r8) :: surface_elevation(1)
   real(r8) :: grav
   real(r8) :: surface_pressure(ens_size), mbar(nlevels, ens_size)
   real(r8) :: pressure(nlevels, ens_size)
   real(r8) :: perturb_ph(nlevels, ens_size)
   real(r8) :: tv(nlevels, ens_size)  ! Virtual temperature, top to bottom

   grav = 9.80665
   ! this is for surface obs
   level_one = 1

                      !! get the surface pressure from the ens_handle
   !call get_staggered_values_from_qty(ens_handle, ens_size, QTY_SURFACE_PRESSURE, &
   !                                   lon_index, lat_index, level_one, qty, surface_pressure, status1)
   !
   ! get the surface elevation from the phis, including stagger if needed
   call get_quad_values(1, lon_index, lat_index, QTY_SURFACE_ELEVATION, qty, surface_elevation)
                        !! Build the pressure columns for the entire ensemble
   !call build_cam_pressure_columns(ens_size, surface_pressure, nlevels, pressure)
   call get_values_from_single_level(ens_handle, ens_size, QTY_SURFACE_PRESSURE, &
                                     lon_index, lat_index, level_one, surface_pressure, status1)
   call get_pressure_column(ens_size, pressure)

   call compute_virtual_temperature(ens_handle, ens_size, lon_index, lat_index, nlevels, tv, status1)

   if (status1 /= 0) then
      my_status = status1
      return
   end if

   ! compute the height columns for each ensemble member - no mbar() argument here.
   ! (you cannot just pass 1.0 in for the mbar() array; it uses a different gas constant
   ! in the variable mean mass case.)
   do imember = 1, ens_size
      call build_heights(nlevels, surface_pressure(imember), surface_elevation(1), &
                         pressure(:, imember), tv(:, imember), height_array(:, imember))
   end do

!                if (debug_level > 100) then
!                 do imember = 1, ens_size
!                  print *, ''
!                  print *, 'geopotential, member: ', imember
!                  do k = 1, nlevels
!                    print*, 'tv(level)    ', k, tv(k, imember)
!                  enddo
!                  do k = 1, nlevels
!                    print*, 'height(level)', k, height_array(k, imember)
!                  enddo
!                 enddo
!                endif

   ! convert entire array to geometric height (from potential height)

   call gph2gmh(height_array, grid_data%latitude(lon_index, lat_index))
   print *, 'height(level)', 5, height_array(5, imember)

   do k = 1, nlevels
      call get_values_from_single_level(ens_handle, ens_size, QTY_GEOPOTENTIAL_HEIGHT, &
                                        lon_index, lat_index, k, perturb_ph(k, :), status1)

      ! convert perturbation geopotential height to full geopotential altitude
      perturb_ph(k, :) = (perturb_ph(k, :) + grid_data%phb(lon_index, lat_index, k))/grav
   end do
   height_array = perturb_ph
   call gph2gmh(height_array, grid_data%latitude(lon_index, lat_index))

   print *, 'height(level)', 5, height_array(5, imember)

   if (debug_level > 100) then
      do imember = 1, ens_size
         print *, ''
         print *, 'geometric, member: ', imember
         do k = 1, nlevels
            print *, 'height(level)', k, height_array(k, imember)
         end do
      end do
   end if

   my_status(:) = 0

end subroutine pangu_height_levels

!-----------------------------------------------------------------------
!> Compute the virtual temperature at the midpoints
!>
!> this version does all ensemble members at once.
!>

subroutine compute_virtual_temperature(ens_handle, ens_size, lon_index, lat_index, nlevels, tv, istatus)

   type(ensemble_type), intent(in)   :: ens_handle
   integer, intent(in)   :: ens_size
   integer, intent(in)   :: lon_index
   integer, intent(in)   :: lat_index
   integer, intent(in)   :: nlevels
   real(r8), intent(out)  :: tv(nlevels, ens_size)
   integer, intent(out)  :: istatus

   integer :: k
   real(r8) :: temperature(ens_size), specific_humidity(ens_size)

   !>@todo this should come from a model specific constant module.
   !> the forward operators and model_mod should use it.
   real(r8), parameter :: rd = 287.05_r8 ! dry air gas constant
   real(r8), parameter :: rv = 461.51_r8 ! wet air gas constant
   real(r8), parameter :: rr_factor = (rv/rd) - 1.0_r8

   ! construct a virtual temperature column, one for each ensemble member
   do k = 1, nlevels
      ! temperature
      call get_values_from_single_level(ens_handle, ens_size, QTY_TEMPERATURE, &
                                        lon_index, lat_index, k, temperature, istatus)
      if (istatus < 0) return

      ! specific humidity
      call get_values_from_single_level(ens_handle, ens_size, QTY_SPECIFIC_HUMIDITY, &
                                        lon_index, lat_index, k, specific_humidity, istatus)
      if (istatus < 0) return

      !>tv == virtual temperature.
      tv(k, :) = temperature(:)*(1.0_r8 + rr_factor*specific_humidity(:))
      !print*, 'tv(levels)', k,tv(k,1), temperature(1), specific_humidity(1)
   end do

end subroutine compute_virtual_temperature

subroutine build_heights(nlevels, p_surf, h_surf, pressure, virtual_temp, height_midpts, mbar)

   integer, intent(in)  :: nlevels                            ! Number of vertical levels
   real(r8), intent(in)  :: p_surf                             ! Surface pressure (pascals)
   real(r8), intent(in)  :: h_surf                             ! Surface height (m)
   real(r8), intent(in)  :: pressure(nlevels)                 ! Pressure
   real(r8), intent(in)  :: virtual_temp(nlevels)             ! Virtual Temperature
   real(r8), intent(out) :: height_midpts(nlevels)             ! Geopotential height at midpoints, top to bottom
   real(r8), intent(in), optional :: mbar(nlevels)            ! Factor to support for variable gas constant

   ! Local variables
   !>@todo FIXME can we use the types_mod values here?  or have a model constants module?
   real(r8), parameter :: const_r = 287.04_r8    ! Different than model_heights (dry air gas constant)
   real(r8), parameter :: universal_gas_constant = 8314.0_r8 ! [J/K/kmol]
   real(r8), parameter :: g0 = 9.80616_r8        ! Different than model_heights (gph2gmh:G) !

   integer  :: k, l

   ! an array now: real(r8), parameter :: rbyg=r/g0
   real(r8) :: pterm(nlevels)   ! vertical scratch space, to improve computational efficiency
   real(r8) :: r_g0_tv(nlevels) ! rbyg=r/g0 * tv
   real(r8) :: pm_ln(nlevels + 1) ! logs of midpoint pressures plus surface interface pressure

   ! cam uses a uniform gas constant value, but high top
   ! models like waccm change the gas constant with height.
   ! allow for the calling code to pass in an array of r.

   ! if mbar() array is present notice that the units are different
   ! for the gas constants, so an mbar() array of 1.0s will NOT give
   ! the same results as if it isn't present.

   if (present(mbar)) then
      r_g0_tv(:) = (universal_gas_constant/(mbar(:)*g0))*virtual_temp(:)
   else
      r_g0_tv(:) = (const_r/g0)*virtual_temp(:)
   end if

   ! calculate the log of the pressure column midpoints.
   ! items 1:nlevels are the midpoints, but NOTICE THAT
   ! the pressure at nlevels+1 is the pressure of the
   ! actual surface interface, not a midpoint!!

   ! The original routine that did this conversion allowed the bottom boundary of the lowest pressure
   ! level to be something other than the surface pressure and computed it with the following :
   ! p_surf * grid_data%hybi%vals(nlevels+1)   ! surface interface
   ! However, all modern SE models appear to have the lowest level boundary the same as the surface. This
   ! means that this can be replaced by just the surface pressure. If this is not true, careful thought is
   ! required, especially for the dry_mass_vertical_coordinate.

   ! Put the log of the surface pressure in the top entry of the log pressure column for the conversion
   pm_ln(1) = log(p_surf)

   ! Some weird vertical coord could have top pressure 0, so leave this check
   where (pressure > 0.0_r8)
      pm_ln(2:nlevels + 1) = log(pressure)
      else where (pressure <= 0.0_r8)
      pm_ln(1:nlevels) = 0
   end where

   !        height_midpts(1)=top  ->  height_midpts(nlevels)=bottom
   !
   ! level
   ! 1/2    ---------------------------------------------------------------
   ! 1      - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - top
   ! 3/2    ---------------------------------------------------------------
   !
   !                 ---------------------------------------------
   !        --------/                                             \--------
   !                - - - - - - - - - - - - - - - - - - - - - - - -
   ! NL     - - - /                                                 \ - - - bottom
   !              ---------------------------------------------------
   ! NL+1/2 -----/|||||||||||||||||||||||||||||||||||||||||||||||||||\-----

   ! now the finite differences.
   ! Eq 3.a.109 has 5 piecewise (PW) terms.  The numbers below refer to each of these terms
   ! in the order they're listed in the paper.

   !
   ! See 2nd PW term here: Eq 3.a.109  where l=K,k<K  h(k,l) = 1/2 * ln [  p(k+1) / p(k) ]
   !

   do k = 2, nlevels
      height_midpts(k) = h_surf + r_g0_tv(k)*0.5_r8*(pm_ln(k) - pm_ln(k + 1))
   end do
   height_midpts(1) = h_surf + r_g0_tv(1)*(pm_ln(1) - pm_ln(2))

   !
   ! See 4th PW term here: Eq 3.a.109  where l=K,k<K  h(k,l) = 1/2*ln[pi*pi/(p(k-1)*p(k))
   !

   do k = 2, nlevels
      height_midpts(k) = height_midpts(k) + r_g0_tv(1)* &
                         (pm_ln(1) - 0.5_r8*(pm_ln(2) + pm_ln(1)))
   end do

   !
   ! See 3rd PW term here:  Eqs 1.14 & 3.a.109 where l>K, k<K
   !                                h(k,l) = 1/2 * ln [ p(l+1)/p(l-1) ]

   ! don't recompute the same values multiple times;
   ! compute once and put into a temporary array.
   ! (see the double nested loops below with k and l)

   ! this is really a matrix multiply with a upper triangular
   ! matrix so it simplifies to a doubly nested loop.

   ! pterm(1) and (nlevels) are never used, but to prevent
   ! confusion when debugging set them to 0 so they don't
   ! look like uninitialized variables.
   pterm(nlevels) = 0.0_r8
   do k = 2, nlevels - 1
      pterm(k) = r_g0_tv(k)*0.5_r8*(pm_ln(k) - pm_ln(k + 2))
   end do
   pterm(1) = 0.0_r8

   do k = 3, nlevels
      do l = 2, k - 1
         height_midpts(k) = height_midpts(k) + pterm(l)
      end do
   end do

end subroutine build_heights

!-----------------------------------------------------------------------
!>  Convert a 2d array of geopotential altitudes to mean sea level altitudes.
!>  To avoid overflow with very high model tops, convert to km first, compute,
!>  then convert back.

subroutine gph2gmh(h, lat)
   real(r8), intent(inout) :: h(:, :)    ! geopotential altitude in m
   real(r8), intent(in)    :: lat       ! latitude in degrees.

   real(r8), parameter ::  be = 6356.7516_r8        ! min earth radius, km
   real(r8), parameter ::  ae = 6378.1363_r8        ! max earth radius, km
   real(r8), parameter ::  G = 0.00980665_r8        ! WMO reference g value, km/s**2, at 45.542N(S)

   real(r8) :: g0
   real(r8) :: r0
   real(r8) :: latr

   integer :: i, j

   latr = lat*DEG2RAD  ! convert to radians
   call compute_surface_gravity(latr, g0)

   ! compute local earth's radius using ellipse equation

   r0 = sqrt(ae**2*cos(latr)**2 + be**2*sin(latr)**2)

   ! Compute altitude above sea level
   do j = 1, size(h, 2)
      do i = 1, size(h, 1)
         h(i, j) = h(i, j)/1000.0_r8   ! m to km
         if (((g0*r0)/G) - h(i, j) > 0) &
            h(i, j) = (r0*h(i, j))/(((g0*r0)/G) - h(i, j))
         h(i, j) = h(i, j)*1000.0_r8   ! km to m
      end do
   end do

end subroutine gph2gmh

!-----------------------------------------------------------------------
!> This subroutine computes the Earth's gravity at any latitude.
!> The model assumes the Earth is an oblate spheriod rotating at
!> the Earth's spin rate.  The model was taken from
!> "Geophysical Geodesy, Kurt Lambeck, 1988".
!>
!>  input:    xlat, latitude in radians
!>  output:   galt, gravity at the given lat, km/sec**2

subroutine compute_surface_gravity(xlat, galt)
   real(r8), intent(in)  :: xlat
   real(r8), intent(out) :: galt

   real(r8), parameter :: xmu = 398600.4415_r8         ! km^3/s^2
   real(r8), parameter :: ae = 6378.1363_r8           ! km
   real(r8), parameter :: f = 1.0_r8/298.2564_r8
   real(r8), parameter :: xm = 0.003468_r8            !
   real(r8), parameter :: f2 = 5.3481622134089e-03_r8 ! f2 = -f + 5.0* 0.50*xm - 17.0/14.0*f*xm + 15.0/4.0*xm**2
   real(r8), parameter :: f4 = 2.3448248012911e-05_r8 ! f4 = -f**2* 0.50 + 5.0* 0.50*f*xm

   real(r8) :: g

   ! gravity at the equator, km/s2
   real(r8), parameter :: ge = xmu/ae**2/(1.0_r8 - f + 1.5_r8*xm - 15.0_r8/14.0_r8*xm*f)

   ! compute gravity at any latitude, km/s2
   g = ge*(1.0_r8 + f2*(sin(xlat))**2 - 1.0_r8/4.0_r8*f4*(sin(2.0_r8*xlat))**2)

   ! at a fixed altitude of 0.0, g and galt are the same
   galt = g

end subroutine compute_surface_gravity

!#######################################################################

subroutine get_model_horizontal_location(i, j, var_type, lon, lat)

   integer, intent(in)  :: i, j, var_type
   real(r8), intent(out) :: lon, lat
   integer               :: id

   ! given i, j indices into the horizontal grid return the lat/long.
   ! if u or v staggering use the staggered grids, otherwise use the mass
   ! grid.  this code has changed -- earlier versions only had the mass
   ! grid available and used it to compute cell midpoints and called them
   ! the staggered points.  now that all three grids are being read, look
   ! up the point locations directly from the appropriate array.

   id = 1

   lon = grid_data%longitude(i, j)
   lat = grid_data%latitude(i, j)

   do while (lon < 0.0_r8)
      lon = lon + 360.0_r8
   end do
   do while (lon > 360.0_r8)
      lon = lon - 360.0_r8
   end do

end subroutine get_model_horizontal_location

!#######################################################################

!------------------------------------------------------------------
! write any additional attributes to the output and diagnostic files

subroutine nc_write_model_atts(ncid, id)
   !-----------------------------------------------------------------
   ! Writes the model-specific attributes to a netCDF file

   integer, intent(in) :: ncid      ! netCDF file identifier

   character(len=*), parameter :: routine = 'nc_write_model_atts'

   ! currently unused, but if needed could be added back in.  these fields
   ! only appear to be supported in certain projections, so the code should
   ! test to be sure they exist before trying to read them from the netcdf file.

   integer :: id
   integer :: i, ret, tmp

   character(len=129) :: title

   character(len=129), allocatable, dimension(:) :: textblock
   integer :: ind
   character(len=NF90_MAX_NAME) :: attname, varname
   character(len=129) :: unitsval, descriptionval, coordinatesval, long_nameval, coordinate_char
   logical               :: debug = .false.
   character(len=256) :: filename

   id = 1

   ! use netcdf file id for identification
   write (filename, *) 'ncid', ncid

   !-------------------------------------------------------------------------------
   ! Put file into define mode and
   ! Write Global Attributes
   !-------------------------------------------------------------------------------
   call nc_begin_define_mode(ncid, routine)

   call nc_add_global_creation_time(ncid, routine)

   call nc_add_global_attribute(ncid, "model_source", source, routine)   
   call nc_add_global_attribute(ncid, "model", "pangu", routine)

   !-----------------------------------------------------------------
   ! Define the dimensions IDs
   !-----------------------------------------------------------------
   call nc_define_dimension(ncid, 'west_east', grid_data%we, routine)
   call nc_define_dimension(ncid, 'south_north', grid_data%sn, routine)
   call nc_define_dimension(ncid, 'bottom_top', grid_data%bt, routine)

   !----------------------------------------------------------------------------
   ! Create the Coordinate Variables and the Attributes
   ! The contents will be written in a later block of code.
   !----------------------------------------------------------------------------

   ! Mass Grid Longitudes
   call nc_define_real_variable(ncid, 'XLONG', (/'west_east  ', 'south_north'/), routine)
   call nc_add_attribute_to_variable(ncid, 'XLONG', 'long_name', 'longitude', routine)
   call nc_add_attribute_to_variable(ncid, 'XLONG', 'units', 'degrees_east', routine)

   ! Mass Grid Latitudes
   call nc_define_real_variable(ncid, 'XLAT', (/'west_east  ', 'south_north'/), routine)
   call nc_add_attribute_to_variable(ncid, 'XLAT', 'long_name', 'latitude', routine)
   call nc_add_attribute_to_variable(ncid, 'XLAT', 'units', 'degrees_north', routine)

   ! Vertical Grid
   ! define as
   call nc_define_real_variable(ncid, 'level', (/'bottom_top'/), routine)
   call nc_add_attribute_to_variable(ncid, 'level', 'long_name', 'level index', routine)
   call nc_add_attribute_to_variable(ncid, 'level', 'units', 'nondimension', routine)

   call nc_end_define_mode(ncid, routine)

   !----------------------------------------------------------------------------
   ! Fill the coordinate variables
   !----------------------------------------------------------------------------

   call nc_put_variable(ncid, 'XLONG', grid_data%longitude, routine)
   call nc_put_variable(ncid, 'XLAT', grid_data%latitude, routine)
   call nc_put_variable(ncid, 'level', (/(i, i=1, grid_data%bt)/), routine)

   ! flush any pending i/o to disk
   call nc_synchronize_file(ncid, routine)

end subroutine nc_write_model_atts

!-----------------------------------------------------------------------
!>
!> Does any shutdown and clean-up needed for model.
!>

subroutine end_model()

   ! deallocate arrays from grid and anything else

   deallocate (grid_data%longitude)
   deallocate (grid_data%latitude)
   deallocate (grid_data%hgt)
   deallocate (grid_data%land)
   deallocate (grid_data%phb)

   ! deallocate(phis)

   ! call free_std_atm_tables()

   call finalize_quad_interp(interp_nonstaggered)

end subroutine end_model

!#######################################################################

subroutine get_pressure_column(ens_size, pressure_array)
   !> replace get_model_pressure_profile_distrib

   integer, intent(in)  :: ens_size
   real(r8), intent(out) :: pressure_array(ref_nlevels, ens_size)

   integer :: j

   do j = 1, ens_size
      pressure_array(1:ref_nlevels, j) = (/100000, 92500, 85000, 70000, 60000, 50000, 40000, &
                                           30000, 25000, 20000, 15000, 10000, 5000/)
   end do

end subroutine get_pressure_column

!--------------------------------------------------------------------
!> for pressure and one flavor of scale height
!> smaller numbers are further away from the surface.
!> for height, levels,  and the other flavor of scale height
!> the opposite is true.  set this once at init time.

subroutine init_sign_of_vert_units()

   if (vertical_localization_type == VERTISHEIGHT .or. &
       vertical_localization_type == VERTISLEVEL) then
      higher_is_smaller = .false.

   else if (vertical_localization_type == VERTISSCALEHEIGHT) then
      ! FIXME: note from nick on scale height:
      !  If no_normalization_of_scale_heights is true, then SH=log(pressure),
      !  and scale height will decrease with increasing height.
      !  However, if it is false then SH= -1*log(pressure/surface_pressure)
      !  and it will increase with increasing height.
      ! cno_normalization_of_scale_heights == .true.
      higher_is_smaller = .true.

   else
      ! pressure
      higher_is_smaller = .true.

   end if

end subroutine init_sign_of_vert_units

!-----------------------------------------------------------------------
!> return the level indices and fraction across the level.
!> level 1 is model bottom, pressure is smallest at the top, l
!> level N is model top, pressure is largest at the bottom,
!> so the values *are* inverted
!> in the array.
!> fract = 0 means full lev1 value,
!> fract = 1 means full lev2 value.
!> return non-zero if value outside valid range.

subroutine pressure_to_level(nlevels, pressures, p_val, &
                             lev1, lev2, fract, my_status)

   integer, intent(in)  :: nlevels
   real(r8), intent(in)  :: pressures(:)
   real(r8), intent(in)  :: p_val
   integer, intent(out) :: lev1
   integer, intent(out) :: lev2
   real(r8), intent(out) :: fract
   integer, intent(out) :: my_status

   my_status = 0
   call find_enclosing_indices(nlevels, pressures(1:ref_nlevels), p_val, lev1, lev2, fract, my_status, &
                               inverted=.true., log_scale=.false.)

   if (my_status /= 0) my_status = 10

end subroutine pressure_to_level

!-----------------------------------------------------------------------
!> return the level indices and fraction across the level.
!> level 1 is model bottom, height is smallest at the bottom
!> level N is model top, height is largest at the top
!> so the values *are not* inverted in the array.
!> fract = 0 means full lev1 value,
!> fract = 1 means full lev2 value.
!> return non-zero if value outside valid range.

subroutine height_to_level(nlevels, heights, h_val, &
                           lev1, lev2, fract, my_status)

   integer, intent(in)  :: nlevels
   real(r8), intent(in)  :: heights(:)
   real(r8), intent(in)  :: h_val
   integer, intent(out) :: lev1
   integer, intent(out) :: lev2
   real(r8), intent(out) :: fract
   integer, intent(out) :: my_status

   character(len=*), parameter :: routine = 'height_to_level:'

   my_status = 0

   call find_enclosing_indices(nlevels, heights, h_val, lev1, lev2, fract, my_status, &
                               inverted=.false., log_scale=.false.)

   if (my_status /= 0) my_status = 11

end subroutine height_to_level
!--------------------------------------------------------------------

!> create an ensemble of states from a single state.  if this routine is
!> not specialized by the model_mod, the default is to add gaussian noise
!> to all parts of the state vector.  unless the model_mod wants to add
!> different amounts of noise to different parts of the state, the model_mod
!> can use this routine.

subroutine pert_model_copies(state_ens_handle, ens_size, pert_amp, interf_provided)

   type(ensemble_type), intent(inout) :: state_ens_handle
   integer, intent(in)    :: ens_size
   real(r8), intent(in)    :: pert_amp ! not used
   logical, intent(out)   :: interf_provided

   ! Perturbs model states for generating initial ensembles.
   ! Because this requires some care when using - see the comments in the
   ! code below - you must set a namelist variable to enable this functionality.

   ! Using this routine is not a substitute for a good set of real initial
   ! condition files.  Intended as a last resort, this routine should be used
   ! to start a period of free-running to allow the model to spin-up a set of
   ! internally consistent states with their own structure before assimilating
   ! a set of observations.  A good ensemble of boundary conditions is important
   ! to evolve the ensemble members differently, which is the goal.
   !

   real(r8)              :: pert_amount = 0.005   ! 0.5%

   real(r8)              :: pert_ampl, range
   real(r8)              :: minv, maxv, temp
   type(random_seq_type) :: random_seq
   integer               :: id, i, j, s, e
   logical, allocatable  :: within_range(:)
   integer               :: num_variables, number_of_state_variables
   real(r8), allocatable :: min_var(:), max_var(:)
   integer(i8)           :: start_ind, end_ind
   integer(i8), allocatable :: var_list(:)
   integer               :: count, copy

   ! generally you do not want to just perturb a single state to begin an
   ! experiment, especially for a regional weather model, because the
   ! resulting fields will have spread but they won't have organized features.
   ! we have had good luck with some global atmosphere models where there is
   ! a lot of model divergence; after a few days of running they evolve into
   ! plausible conditions that allow assimilation of real obs.
   !
   ! if you really need to start with a single state and proceed, the suggestion
   ! is to start with small magnitude perturbations and then get a good ensemble
   ! of boundary conditions and run the model for a while (many days) to let it
   ! evolve into plausible weather patterns.  then start assimilating real obs.
   !
   ! using this routine requires you to set the new namelist item
   ! 'allow_perturbed_ics' to true so you have to read the warnings here or
   ! in the html docs.
   !
   ! this code will add random noise field by field (T, U, V, etc), and new values
   ! will not exceed the original max or min values for each field.  this means
   ! it will not generate illegal values (e.g. negatives for percentages or
   ! number concentrations) but it also means that if all values in a field are
   ! identical (e.g. all 0.0) this routine will not change those values.  the code
   ! can easily be modified to set allowed min and max values here instead of
   ! using the incoming field min and max values; but you will have to modify
   ! the code below to enable that functionality.

   if (.not. allow_perturbed_ics) then
      call error_handler(E_ERR, 'pert_model_copies', &
                         'starting pangu model from a single vector requires additional steps', &
                         source, &
                         text2='see comments in pangu/model_mod.f90::pert_model_copies()')
   end if

   ! NOT REACHED unless allow_perturbed_ics is true in the namelist

   ! start of pert code
   interf_provided = .true.

   ! Make space for the state vector index numbers that are
   ! physically located on my task and get the global numbers.

   allocate (var_list(get_my_num_vars(state_ens_handle)))
   call get_my_vars(state_ens_handle, var_list)

   ! count up the total number of variables across all domains.
   num_variables = 0
   number_of_state_variables = get_number_of_state_variables(state_variables)
   num_variables = num_variables + number_of_state_variables

   ! get the global min/max on a variable by variable basis
   allocate (min_var(num_variables), max_var(num_variables))
   allocate (within_range(state_ens_handle%my_num_vars))

   count = 1
   do i = 1, number_of_state_variables

      start_ind = get_index_start(domain_id, i)
      end_ind = get_index_end(domain_id, i)

      ! at this point we only have 1 ensemble
      within_range = (var_list >= start_ind .and. var_list <= end_ind)
      min_var(count) = minval(state_ens_handle%copies(1, :), MASK=within_range)
      max_var(count) = maxval(state_ens_handle%copies(1, :), MASK=within_range)

      count = count + 1
   end do

   ! find the global min/max values across all tasks.
   call all_reduce_min_max(min_var, max_var, num_variables)

   deallocate (within_range)

   ! Now do the perturbing

   ! using task id as the seed for the random number generator is ok
   ! because pert_model_copies() is only called once on any single task.
   ! it perturbs all ensemble members for the items in the state vector
   ! that it owns. because the decomposition will be different with a
   ! different task count, you will NOT get the same result if you change
   ! the number of tasks.

   call init_random_seq(random_seq, my_task_id() + 1)

   count = 1 ! min and max are numbered 1 to n, where n is the total number of variables (all domains)
   do i = 1, number_of_state_variables

      start_ind = get_index_start(domain_id, i)
      end_ind = get_index_end(domain_id, i)

                        !! Option 1:
                        !! make the perturbation amplitude N% of the total
                        !! range of this variable.  values could vary a lot
                        !! over some of the types, like pressure
      range = max_var(count) - min_var(count)
      pert_ampl = pert_amount*range

      do j = 1, state_ens_handle%my_num_vars
         ! is this state variable index the current variable type we're perturbing?
         if (var_list(j) >= start_ind .and. var_list(j) <= end_ind) then
            do copy = 1, ens_size
                                                !! Option 2: perturb each value individually
                                                !! make the perturbation amplitude N% of this value
               !pert_ampl = pert_amount * state_ens_handle%copies(copy, j)
               state_ens_handle%copies(copy, j) = random_gaussian(random_seq, state_ens_handle%copies(copy, j), pert_ampl)
            end do

            ! keep variable from exceeding the original range
            state_ens_handle%copies(1:ens_size, j) = max(min_var(count), state_ens_handle%copies(1:ens_size, j))
            state_ens_handle%copies(1:ens_size, j) = min(max_var(count), state_ens_handle%copies(1:ens_size, j))

         end if
      end do

      count = count + 1

   end do

   deallocate (var_list, min_var, max_var)

end subroutine pert_model_copies

!#######################################################################

subroutine fill_model_stagger_info(stagger)
   type(pangu_stagger), intent(inout) :: stagger
   allocate (stagger%qty_stagger(0:get_num_quantities()))

   stagger%qty_stagger = STAGGER_NONE

end subroutine fill_model_stagger_info

subroutine setup_interpolation()

   ! mass points at cell centers
   call init_quad_interp(GRID_QUAD_IRREG_SPACED_REGULAR, grid_data%we, grid_data%sn, &
                         QUAD_LOCATED_CELL_CENTERS, &
                         global=.true., spans_lon_zero=.true., pole_wrap=.true., &
                         interp_handle=interp_nonstaggered)
   call set_quad_coords(interp_nonstaggered, grid_data%longitude(1:grid_data%we, 1), grid_data%latitude(1, 1:grid_data%sn))

   !call set_quad_coords(interp_nonstaggered, -180.0_r8, 0.25_r8, -90.0_r8, 0.25_r8)

end subroutine setup_interpolation

!-----------------------------------------------------------------------
!> Compute the heights at pressure midpoints
!>
!> this version does all ensemble members at once.

!    subroutine pangu_height_levels(ens_handle, ens_size, lon_index, lat_index, nlevels, qty, height_array, my_status)
!    type(ensemble_type), intent(in)  :: ens_handle
!    integer,             intent(in)  :: ens_size
!    integer,             intent(in)  :: lon_index
!    integer,             intent(in)  :: lat_index
!    integer,             intent(in)  :: nlevels
!    integer,             intent(in)  :: qty
!    real(r8),            intent(out) :: height_array(nlevels, ens_size)
!    integer,             intent(out) :: my_status(ens_size)
!
!    integer  :: k, level_one, imember, status1
!    real(r8) :: surface_elevation(1)
!    real(r8) :: surface_pressure(ens_size)
!    real(r8) :: pressure(nlevels, ens_size)
!    real(r8) :: tv(nlevels, ens_size)  ! Virtual temperature, top to bottom
!
!    ! this is for surface obs
!    level_one = 1
!
!    ! get the surface pressure from the ens_handle
!    call get_staggered_values_from_qty(ens_handle, ens_size, QTY_SURFACE_PRESSURE, &
!                                       lon_index, lat_index, level_one, qty, surface_pressure, status1)
!
!    ! get the surface elevation from the phis, including stagger if needed
!    call get_quad_values(1, lon_index, lat_index, QTY_SURFACE_ELEVATION, qty, surface_elevation)
!
!    call compute_virtual_temperature(ens_handle, ens_size, lon_index, lat_index, nlevels, tv, status1)
!
!    if (status1 /= 0) then
!       my_status = status1
!       return
!    endif
!
!    ! Build the pressure columns for the entire ensemble
!    call get_pressure_column(ens_size, pressure)
!
!
!   ! compute the height columns for each ensemble member - no mbar() argument here.
!   do imember = 1, ens_size
!      call build_heights(nlevels, surface_pressure(imember), surface_elevation(1), &
!                         pressure(:, imember), tv(:, imember), height_array(:, imember))
!   enddo
!
!    if (debug_level > 100) then
!     do imember = 1, ens_size
!      print *, ''
!      print *, 'geopotential, member: ', imember
!      do k = 1, nlevels
!        print*, 'tv(level)    ', k, tv(k, imember)
!      enddo
!      do k = 1, nlevels
!        print*, 'height(level)', k, height_array(k, imember)
!      enddo
!     enddo
!    endif
!
!    ! convert entire array to geometric height (from potential height)
!    call gph2gmh(height_array, grid_data%latitude(lon_index, lat_index))
!
!    if (debug_level > 100) then
!     do imember = 1, ens_size
!      print *, ''
!      print *, 'geometric, member: ', imember
!      do k = 1, nlevels
!        print*, 'height(level)', k, height_array(k, imember)
!      enddo
!     enddo
!    endif
!
!    my_status(:) = 0
!
!    end subroutine pangu_height_levels
!
function get_interp_handle()
   type(quad_interp_handle) :: get_interp_handle
   character(len=*), parameter :: routine = 'get_interp_handle:'

   get_interp_handle = interp_nonstaggered

end function get_interp_handle

subroutine interpolate_values(state_handle, ens_size, location, obs_qty, varid, &
                              interp_vals, istatus)

   type(ensemble_type), intent(in) :: state_handle
   integer, intent(in) :: ens_size
   type(location_type), intent(in) :: location
   integer, intent(in) :: obs_qty
   integer, intent(in) :: varid
   real(r8), intent(out) :: interp_vals(ens_size)
   integer, intent(out) :: istatus(ens_size)

   character(len=*), parameter :: routine = 'interpolate_values:'

   integer  :: which_vert, four_lons(4), four_lats(4)
   real(r8) :: lon_fract, lat_fract
   real(r8) :: lon_lat_vert(3), quad_vals(4, ens_size)
   type(quad_interp_handle) :: interp_handle

   interp_vals(:) = MISSING_R8
   istatus(:) = 99

   interp_handle = get_interp_handle()
   lon_lat_vert = get_location(location)
   which_vert = nint(query_location(location))

   call quad_lon_lat_locate(interp_handle, lon_lat_vert(1), lon_lat_vert(2), &
                            four_lons, four_lats, lon_fract, lat_fract, istatus(1))
   if (istatus(1) /= 0) then
      istatus(:) = 3  ! cannot locate enclosing horizontal quad
      return
   end if

   call get_quad_vals(state_handle, ens_size, varid, obs_qty, four_lons, four_lats, &
                      lon_lat_vert, which_vert, quad_vals, istatus)
   if (any(istatus /= 0)) return

   call quad_lon_lat_evaluate(interp_handle, lon_fract, lat_fract, ens_size, &
                              quad_vals, interp_vals, istatus)
   if (any(istatus /= 0)) then
      istatus(:) = 8   ! cannot evaluate in the quad
      return
   end if

end subroutine interpolate_values

subroutine obs_vertical_to_pressure(ens_handle, location, my_status)

   type(ensemble_type), intent(in)    :: ens_handle
   type(location_type), intent(inout) :: location
   integer, intent(out)   :: my_status

   integer  :: varid, ens_size, status(1), qty
   real(r8) :: pressure_array(ref_nlevels)

   character(len=*), parameter :: routine = 'obs_vertical_to_pressure'

   ens_size = 1

   qty = QTY_PRESSURE
   if (query_location(location) == VERTISSURFACE) then
      qty = QTY_SURFACE_PRESSURE
   end if

   call ok_to_interpolate(qty, varid, my_status)

   if (my_status /= 0) return

   call interpolate_values(ens_handle, ens_size, location, &
                           qty, varid, pressure_array(:), status(:))

   if (status(1) /= 0) then
      my_status = status(1)
      return
   end if

   call set_vertical(location, pressure_array(1), VERTISPRESSURE)

   my_status = 0

end subroutine obs_vertical_to_pressure

!--------------------------------------------------------------------

subroutine obs_vertical_to_height(ens_handle, location, my_status)
   type(ensemble_type), intent(in)    :: ens_handle
   type(location_type), intent(inout) :: location
   integer, intent(out)   :: my_status

   integer  :: varid, ens_size, status(1)
   real(r8) :: height_array(1)

   character(len=*), parameter :: routine = 'obs_vertical_to_height'

   ens_size = 1

   call ok_to_interpolate(QTY_GEOMETRIC_HEIGHT, varid, my_status)
   if (my_status /= 0) return

   call interpolate_values(ens_handle, ens_size, location, &
                           QTY_GEOMETRIC_HEIGHT, varid, height_array(:), status(:))
   if (status(1) /= 0) then
      my_status = status(1)
      return
   end if

   call set_vertical(location, height_array(1), VERTISHEIGHT)

   my_status = 0

end subroutine obs_vertical_to_height

!--------------------------------------------------------------------

subroutine obs_vertical_to_level(ens_handle, location, my_status)
   type(ensemble_type), intent(in)    :: ens_handle
   type(location_type), intent(inout) :: location
   integer, intent(out)   :: my_status

   integer  :: varid, ens_size, status(1)
   real(r8) :: level_array(1)

   ens_size = 1
   varid = -1

   call interpolate_values(ens_handle, ens_size, location, &
                           QTY_VERTLEVEL, varid, level_array(:), status(:))
   if (status(1) /= 0) then
      my_status = status(1)
      return
   end if

   call set_vertical(location, level_array(1), VERTISLEVEL)

   my_status = 0

end subroutine obs_vertical_to_level

!--------------------------------------------------------------------

subroutine obs_vertical_to_scaleheight(ens_handle, location, my_status)
   type(ensemble_type), intent(in)    :: ens_handle
   type(location_type), intent(inout) :: location
   integer, intent(out)   :: my_status

   integer  :: varid1, varid2, ens_size, status(1), ptype
   real(r8) :: pressure_array(1), surface_pressure_array(1)
   real(r8) :: scaleheight_val

   character(len=*), parameter :: routine = 'obs_vertical_to_scaleheight'

   ens_size = 1

   ! if this location is on the surface, use the surface pressure field
   ! in the computations below.  otherwise use the 3d pressure field.
   if (query_location(location) == VERTISSURFACE) then
      ptype = QTY_SURFACE_PRESSURE
   else
      ptype = QTY_PRESSURE
   end if

   ! there are 4 cases here.

   if (no_normalization_of_scale_heights) then

      ! take log of pressure, either surface pressure or regular pressure

      call ok_to_interpolate(ptype, varid1, my_status)
      if (my_status /= 0) return

      !>@todo FIXME IFF the obs location is already pressure, we can take it at
      !> face value here and not interpolate it.  however it won't fail if the
      !> pressure here is less than the ensemble mean pressure at this point.
      !> is that ok?

      if (ptype == QTY_PRESSURE .and. is_vertical(location, "PRESSURE")) then
         pressure_array(:) = query_location(location, "VLOC")
         my_status = 0
      else
         call interpolate_values(ens_handle, ens_size, location, ptype, varid1, &
                                 pressure_array(:), status(:))
         if (status(1) /= 0) then
            my_status = status(1)
            return
         end if
      end if

      scaleheight_val = log(pressure_array(1))

   else

      ! handle surface obs separately here.
      if (query_location(location) == VERTISSURFACE) then

         scaleheight_val = 0.0_r8  ! -log(1.0)

      else

         call ok_to_interpolate(QTY_PRESSURE, varid1, my_status)
         if (my_status /= 0) return

         !>@todo FIXME IFF the obs location is already pressure, we can take it at
         !> face value here and not interpolate it.  however, it can result in negative
         !> scale height values if the pressure is larger than the surface pressure at
         !> that location.  that's what the original cam model_mod did.  is that ok?

         if (ptype == QTY_PRESSURE .and. is_vertical(location, "PRESSURE")) then
            pressure_array(:) = query_location(location, "VLOC")
            my_status = 0
         else
            call interpolate_values(ens_handle, ens_size, location, QTY_PRESSURE, varid1, &
                                    pressure_array(:), status(:))
            if (status(1) /= 0) then
               my_status = status(1)
               return
            end if
         end if

         call ok_to_interpolate(QTY_SURFACE_PRESSURE, varid2, my_status)
         if (my_status /= 0) return

         call interpolate_values(ens_handle, ens_size, location, QTY_SURFACE_PRESSURE, varid2, &
                                 surface_pressure_array(:), status(:))
         if (status(1) /= 0) then
            my_status = status(1)
            return
         end if

         scaleheight_val = scale_height(pressure_array(1), surface_pressure_array(1), no_normalization_of_scale_heights)

      end if

   end if

   call set_vertical(location, scaleheight_val, VERTISSCALEHEIGHT)

   my_status = 0

end subroutine obs_vertical_to_scaleheight

subroutine convert_vertical_obs(ens_handle, num, locs, loc_qtys, loc_types, &
                                which_vert, my_status)

   type(ensemble_type), intent(in)    :: ens_handle
   integer, intent(in)    :: num
   type(location_type), intent(inout) :: locs(:)
   integer, intent(in)    :: loc_qtys(:)
   integer, intent(in)    :: loc_types(:)
   integer, intent(in)    :: which_vert
   integer, intent(out)   :: my_status(:)
   character(len=512) :: string1
   character(len=*), parameter :: routine = 'convert_vertical_obs'

   integer :: current_vert_type, i

   do i = 1, num
      current_vert_type = nint(query_location(locs(i)))

      if ((current_vert_type == which_vert) .or. &
          (current_vert_type == VERTISUNDEF)) then
         my_status(i) = 0
         cycle
      end if

      select case (which_vert)
      case (VERTISPRESSURE)
         call obs_vertical_to_pressure(ens_handle, locs(i), my_status(i))
      case (VERTISHEIGHT)
         call obs_vertical_to_height(ens_handle, locs(i), my_status(i))
      case (VERTISLEVEL)
         call obs_vertical_to_level(ens_handle, locs(i), my_status(i))
      case (VERTISSCALEHEIGHT)
         call obs_vertical_to_scaleheight(ens_handle, locs(i), my_status(i))
      case default
         write (string1, *) 'unable to convert vertical obs "', which_vert, '"'
         call error_handler(E_ERR, routine, string1, source)
      end select
   end do

end subroutine convert_vertical_obs

subroutine convert_vert_one_obs(ens_handle, loc, otype, vert_type, status1)
   type(ensemble_type), intent(in)    :: ens_handle
   type(location_type), intent(inout) :: loc
   integer, intent(in)    :: otype
   integer, intent(in)    :: vert_type
   integer, intent(out)   :: status1

   type(location_type) :: bl(1)
   integer :: bq(1), bt(1), status(1)

   ! these need to be arrays.  kinda a pain.
   bl(1) = loc
   bt(1) = otype
   bq(1) = get_quantity_for_type_of_obs(otype)

   call convert_vertical_obs(ens_handle, 1, bl, bq, bt, &
                             vert_type, status)

   status1 = status(1)

   if (status1 /= 0) return

   loc = bl(1)

end subroutine convert_vert_one_obs

subroutine convert_vertical_state(ens_handle, num, locs, loc_qtys, loc_indx, &
                                  which_vert, istatus)
   type(ensemble_type), intent(in)    :: ens_handle
   integer, intent(in)    :: num
   type(location_type), intent(inout) :: locs(:)
   integer, intent(in)    :: loc_qtys(:)
   integer(i8), intent(in)    :: loc_indx(:)
   integer, intent(in)    :: which_vert
   integer, intent(out)   :: istatus
   character(len=512) :: string1
   character(len=*), parameter :: routine = 'convert_vertical_state'

   integer :: current_vert_type, ens_size, i

   ens_size = 1
   !print *, "current vert type", nint(query_location(locs(10)))
   !print *, "which vert", which_vert
   do i = 1, num
      current_vert_type = nint(query_location(locs(i)))

      if (current_vert_type == which_vert) cycle
      if (current_vert_type == VERTISUNDEF) cycle

      select case (which_vert)
      case (VERTISPRESSURE)
         call state_vertical_to_pressure(ens_handle, ens_size, locs(i), loc_indx(i), loc_qtys(i))
      case (VERTISHEIGHT)
         call state_vertical_to_height(ens_handle, ens_size, locs(i), loc_indx(i), loc_qtys(i))
      case (VERTISLEVEL)
         call state_vertical_to_level(ens_size, locs(i), loc_indx(i), loc_qtys(i))
      case (VERTISSCALEHEIGHT)
         call state_vertical_to_scaleheight(ens_handle, ens_size, locs(i), loc_indx(i), loc_qtys(i))
      case default
         write (string1, *) 'unable to convert vertical state "', which_vert, '"'
         call error_handler(E_MSG, routine, string1, source)
      end select
   end do

   istatus = 0

end subroutine convert_vertical_state

!--------------------------------------------------------------------

subroutine state_vertical_to_pressure(ens_handle, ens_size, location, location_indx, qty)
   type(ensemble_type), intent(in)    :: ens_handle
   integer, intent(in)    :: ens_size
   type(location_type), intent(inout) :: location
   integer(i8), intent(in)    :: location_indx
   integer, intent(in)    :: qty

   integer  :: iloc, jloc, vloc, myqty, level_one, status1
   integer  :: my_status(ens_size)
   real(r8) :: pressure_array(ref_nlevels), surface_pressure(ens_size)
   real(r8) :: pres_arr(ref_nlevels, 1)

   call get_model_variable_indices(location_indx, iloc, jloc, vloc, kind_index=myqty)

   if (is_surface_field(myqty)) then

      level_one = 1
      call get_values_from_single_level(ens_handle, ens_size, QTY_SURFACE_PRESSURE, &
                                        iloc, jloc, level_one, surface_pressure, status1)

      if (status1 /= 0) then
         return
      end if
      call set_vertical(location, surface_pressure(1), VERTISPRESSURE)
   else
      call get_pressure_column(1, pres_arr)
      pressure_array = pres_arr(:, 1)

      call set_vertical(location, pressure_array(vloc), VERTISPRESSURE)
   end if

end subroutine state_vertical_to_pressure

!--------------------------------------------------------------------

subroutine state_vertical_to_height(ens_handle, ens_size, location, location_indx, qty)
   type(ensemble_type), intent(in)    :: ens_handle
   integer, intent(in)    :: ens_size
   type(location_type), intent(inout) :: location
   integer(i8), intent(in)    :: location_indx
   integer, intent(in)    :: qty

   integer  :: iloc, jloc, vloc, my_status(ens_size)
   real(r8) :: height_array(ref_nlevels, ens_size)

   ! build a height column and a pressure column and find the levels
   call get_model_variable_indices(location_indx, iloc, jloc, vloc)

   call pangu_height_levels(ens_handle, ens_size, iloc, jloc, ref_nlevels, &
                            qty, height_array, my_status)
   print *, "pangu state_vertical_to_height", height_array
   !>@todo FIXME this can only be used if ensemble size is 1
   call set_vertical(location, height_array(vloc, 1), VERTISHEIGHT)

end subroutine state_vertical_to_height

!--------------------------------------------------------------------

subroutine state_vertical_to_scaleheight(ens_handle, ens_size, location, location_indx, qty)
   type(ensemble_type), intent(in)    :: ens_handle
   integer, intent(in)    :: ens_size
   type(location_type), intent(inout) :: location
   integer(i8), intent(in)    :: location_indx
   integer, intent(in)    :: qty

   integer  :: iloc, jloc, vloc, level_one, status1, my_status(ens_size)
   real(r8) :: pressure_array(ref_nlevels)
   real(r8) :: pres_arr(ref_nlevels, 1)
   real(r8) :: surface_pressure(1), scaleheight_val

   !> this is currently only called with an ensemble size of 1 for
   !> vertical conversion.  since it is working only on state variables
   !> we don't expect it to ever fail.

   level_one = 1
   scaleheight_val = MISSING_R8

   if (no_normalization_of_scale_heights) then

      if (query_location(location) == VERTISSURFACE) then

         ! get the surface pressure from the ens_handle
         call get_model_variable_indices(location_indx, iloc, jloc, vloc)

         call get_values_from_single_level(ens_handle, ens_size, QTY_SURFACE_PRESSURE, &
                                           iloc, jloc, level_one, surface_pressure, status1)
         if (status1 /= 0) goto 200

         scaleheight_val = log(surface_pressure(1))

      else

         ! build a pressure column and and find the levels
         call get_model_variable_indices(location_indx, iloc, jloc, vloc)
         call get_pressure_column(1, pres_arr)
         pressure_array = pres_arr(:, 1)
         if (any(my_status /= 0)) goto 200

         scaleheight_val = log(pressure_array(vloc))

      end if

   else

      ! handle surface obs separately here.
      if (query_location(location) == VERTISSURFACE) then

         scaleheight_val = 0.0_r8   ! log(1.0)

      else

         ! build a pressure column and and find the levels
         call get_model_variable_indices(location_indx, iloc, jloc, vloc)

         call get_pressure_column(1, pres_arr)
         pressure_array = pres_arr(:, 1)

         if (any(my_status /= 0)) goto 200

         ! get the surface pressure from the ens_handle
         call get_values_from_single_level(ens_handle, ens_size, QTY_SURFACE_PRESSURE, &
                                           iloc, jloc, level_one, surface_pressure, status1)
         if (status1 /= 0) goto 200

         scaleheight_val = scale_height(pressure_array(vloc), surface_pressure(1), no_normalization_of_scale_heights)

      end if

   end if

200 continue   ! done

   call set_vertical(location, scaleheight_val, VERTISSCALEHEIGHT)

end subroutine state_vertical_to_scaleheight

!--------------------------------------------------------------------

subroutine state_vertical_to_level(ens_size, location, location_indx, qty)
   integer, intent(in)    :: ens_size
   type(location_type), intent(inout) :: location
   integer(i8), intent(in)    :: location_indx
   integer, intent(in)    :: qty

   integer  :: iloc, jloc, vloc

   !>@todo FIXME qty is currently unused.  if we need it, its here.
   !>if we really don't need it, we can remove it.  all the other
   !>corresponding routines like this use it.  (not clear what to
   !>return if field is W or something else with a vertical stagger.)

   call get_model_variable_indices(location_indx, iloc, jloc, vloc)

   call set_vertical(location, real(vloc, r8), VERTISLEVEL)

end subroutine state_vertical_to_level

subroutine get_close_obs(gc, base_loc, base_type, locs, loc_qtys, loc_types, &
                         num_close, close_ind, dist, ens_handle)

   ! The specific type of the base observation, plus the generic kinds list
   ! for either the state or obs lists are available if a more sophisticated
   ! distance computation is needed.

   type(get_close_type), intent(in)  :: gc
   type(location_type), intent(inout) :: base_loc, locs(:)
   integer, intent(in)  :: base_type, loc_qtys(:), loc_types(:)
   integer, intent(out) :: num_close, close_ind(:)
   real(r8), optional, intent(out) :: dist(:)
   type(ensemble_type), optional, intent(in)  :: ens_handle

   character(len=*), parameter :: routine = 'get_close_obs'

   integer :: i, status(1), this, vert_type
   real(r8) :: vert_value, extra_damping_dist
   real(r8), parameter :: LARGE_DIST = 999999.0  ! positive and large

   ! if absolute distances aren't needed, or vertical localization isn't on,
   ! the default version works fine since no conversion will be needed and
   ! there won't be any damping since there are no vert distances.
   if (.not. present(dist) .or. .not. vertical_localization_on()) then
      call loc_get_close_obs(gc, base_loc, base_type, locs, loc_qtys, loc_types, &
                             num_close, close_ind, dist, ens_handle)
      return
   end if

   if (.not. present(ens_handle)) then
      call error_handler(E_ERR, routine, &
                         'unexpected error: cannot convert distances without an ensemble handle', &
                         source)
   end if

   ! does the base obs need conversion first?
   vert_type = query_location(base_loc)

   if (vert_type /= vertical_localization_type) then
      call convert_vert_one_obs(ens_handle, base_loc, base_type, &
                                vertical_localization_type, status(1))
      if (status(1) /= 0) then
         num_close = 0
         return
      end if
   end if

   ! FIXME: is here where we need to compute start of ramp for this
   ! obs type?  should we cache these?  start with doing the computation
   ! each time, then make an array indexed by obs types with the
   ! start of the ramp and fill it in on demand.  have to call for
   ! maxdist(obs_type) and do the math, but just once per type.

   ! ok, distance is needed and we are localizing in the vertical.
   ! call default get close to get potentically close locations
   ! but call without distance so it doesn't do extra work.
   call loc_get_close_obs(gc, base_loc, base_type, locs, loc_qtys, loc_types, &
                          num_close, close_ind)

   ! compute distances, converting vertical first if need be.
   do i = 1, num_close
      this = close_ind(i)

      vert_type = query_location(locs(this))
      !print *, 'close_o, vval, vtype = ', i, query_location(locs(this), 'VLOC'), vert_type

      if (vert_type /= vertical_localization_type) then
         call convert_vertical_obs(ens_handle, 1, locs(this:this), &
                                   loc_qtys(this:this), loc_types(this:this), &
                                   vertical_localization_type, status)
         if (status(1) /= 0) then
            dist(i) = LARGE_DIST
            cycle
         end if

      end if

      dist(i) = get_dist(base_loc, locs(this))

      ! do not try to damp impacts when obs has "vert is undefined".
      ! the impact will go all the way to the model top.
      ! this is the agreed-on functionality.
      if (.not. are_damping .or. vert_type == VERTISUNDEF) cycle

      vert_value = query_location(locs(this), 'VLOC')
      if (above_ramp_start(vert_value, gc, base_type, ramp_end, dist(i), extra_damping_dist)) then
         dist(i) = dist(i) + extra_damping_dist
      end if
   end do

end subroutine get_close_obs

!----------------------------------------------------------------------------

subroutine get_close_state(gc, base_loc, base_type, locs, loc_qtys, loc_indx, &
                           num_close, close_ind, dist, ens_handle)

   ! The specific type of the base observation, plus the generic kinds list
   ! for either the state or obs lists are available if a more sophisticated
   ! distance computation is needed.

   type(get_close_type), intent(in)  :: gc
   type(location_type), intent(inout)  :: base_loc, locs(:)
   integer, intent(in)  :: base_type, loc_qtys(:)
   integer(i8), intent(in)  :: loc_indx(:)
   integer, intent(out) :: num_close, close_ind(:)
   real(r8), optional, intent(out) :: dist(:)
   type(ensemble_type), optional, intent(in)  :: ens_handle

   character(len=*), parameter :: routine = 'get_close_state'

   integer :: i, status, this, vert_type
   real(r8) :: vert_value, extra_damping_dist
   real(r8), parameter :: LARGE_DIST = 999999.0  ! positive and large

   ! if absolute distances aren't needed, or vertical localization isn't on,
   ! the default version works fine since no conversion will be needed and
   ! there won't be any damping since there are no vert distances.
   if (.not. present(dist) .or. .not. vertical_localization_on()) then
      call loc_get_close_state(gc, base_loc, base_type, locs, loc_qtys, loc_indx, &
                               num_close, close_ind, dist, ens_handle)
      return
   end if

   if (.not. present(ens_handle)) then
      call error_handler(E_ERR, routine, &
                         'unexpected error: cannot convert distances without an ensemble handle', &
                         source)
   end if

   ! does the base obs need conversion first?
   vert_type = query_location(base_loc)

   if (vert_type /= vertical_localization_type) then
      call convert_vert_one_obs(ens_handle, base_loc, base_type, &
                                vertical_localization_type, status)
      if (status /= 0) then
         num_close = 0
         return
      end if
   end if

   ! ok, distance is needed and we are localizing in the vertical.
   ! call default get close to get potentically close locations
   ! but call without distance so it doesn't do extra work.
   call loc_get_close_state(gc, base_loc, base_type, locs, loc_qtys, loc_indx, &
                            num_close, close_ind)

   ! compute distances, converting vertical first if need be.
   do i = 1, num_close
      this = close_ind(i)

      vert_type = query_location(locs(this))
      !print *, 'close_s, vval, vtype = ', i, query_location(locs(this), 'VLOC'), vert_type

      if (vert_type /= vertical_localization_type) then
         call convert_vertical_state(ens_handle, 1, locs(this:this), &
                                     loc_qtys(this:this), loc_indx(this:this), &
                                     vertical_localization_type, status)
         if (status /= 0) then
            dist(i) = LARGE_DIST
            cycle
         end if

      end if

      dist(i) = get_dist(base_loc, locs(this))

      ! do not try to damp impacts when obs has "vert is undefined".
      ! the impact will go all the way to the model top.
      ! this is the agreed-on functionality.
      if (.not. are_damping .or. vert_type == VERTISUNDEF) cycle

      vert_value = query_location(locs(this), 'VLOC')
      if (above_ramp_start(vert_value, gc, base_type, ramp_end, dist(i), extra_damping_dist)) then
         dist(i) = dist(i) + extra_damping_dist
      end if
   end do

end subroutine get_close_state

!--------------------------------------------------------------------
!--------------------------------------------------------------------

subroutine read_grid_file_attributes(ncid)

   ! ncid: input, file handl
   ! id:   input, domain id

   integer, intent(in)   :: ncid
   logical, parameter    :: debug = .false.

   ! get meta data and static data we need

   call nc_check(nf90_get_att(ncid, nf90_global, 'DLON', grid_data%dlon), &
                 'static_init_model', 'get_att DLON')
   call nc_check(nf90_get_att(ncid, nf90_global, 'DLAT', grid_data%dlat), &
                 'static_init_model', 'get_att DLAT')
   call nc_check(nf90_get_att(ncid, nf90_global, 'DT', grid_data%dt), &
                 'static_init_model', 'get_att DT')
   write (errstring, *) 'dt from wrfinput_d01 file is: ', grid_data%dt
   call error_handler(E_MSG, ' ', errstring)
   if (debug) write (*, *) ' dlon, dlat, dt are ', grid_data%dlon, &
      grid_data%dlat, grid_data%dt

   RETURN

end subroutine read_grid_file_attributes

subroutine read_model_static_data(ncid)

   ! ncid: input, file handle
   ! id:   input, domain id

   integer, intent(in)   :: ncid
   logical, parameter    :: debug = .false.
   integer               :: var_id

   allocate (grid_data%longitude(1:grid_data%we, 1:grid_data%sn))
   call nc_check(nf90_inq_varid(ncid, "XLONG", var_id), &
                 'read_model_static_data', 'inq_varid XLONG')
   call nc_check(nf90_get_var(ncid, var_id, grid_data%longitude), &
                 'read_model_static_data', 'get_var XLONG')

   allocate (grid_data%latitude(1:grid_data%we, 1:grid_data%sn))
   call nc_check(nf90_inq_varid(ncid, "XLAT", var_id), &
                 'read_model_static_data', 'inq_varid XLAT')
   call nc_check(nf90_get_var(ncid, var_id, grid_data%latitude), &
                 'read_model_static_data', 'get_var XLAT')

   allocate (grid_data%land(1:grid_data%we, 1:grid_data%sn))
   call nc_check(nf90_inq_varid(ncid, "XLAND", var_id), &
                 'read_model_static_data', 'inq_varid XLAND')
   call nc_check(nf90_get_var(ncid, var_id, grid_data%land), &
                 'read_model_static_data', 'get_var XLAND')

   if (debug) then
      write (*, *) ' corners of land '
      write (*, *) grid_data%land(1, 1), grid_data%land(grid_data%we, 1), &
         grid_data%land(1, grid_data%sn), grid_data%land(grid_data%we, &
                                                         grid_data%sn)
   end if

   if (debug) then
      write (*, *) ' corners of lat '
      write (*, *) grid_data%latitude(1, 1), grid_data%latitude(grid_data%we, 1), &
         grid_data%latitude(1, grid_data%sn), &
         grid_data%latitude(grid_data%we, grid_data%sn)
      write (*, *) ' corners of long '
      write (*, *) grid_data%longitude(1, 1), grid_data%longitude(grid_data%we, 1), &
         grid_data%longitude(1, grid_data%sn), &
         grid_data%longitude(grid_data%we, grid_data%sn)
   end if

   allocate (grid_data%hgt(1:grid_data%we, 1:grid_data%sn))
   call nc_check(nf90_inq_varid(ncid, "HGT", var_id), &
                 'read_model_static_data', 'inq_varid HGT')
   call nc_check(nf90_get_var(ncid, var_id, grid_data%hgt), &
                 'read_model_static_data', 'get_var HGT')

   ! get 3D base state geopotential

   allocate (grid_data%phb(1:grid_data%we, 1:grid_data%sn, 1:grid_data%bt))
   call nc_check(nf90_inq_varid(ncid, "PHB", var_id), &
                 'read_model_static_data', 'inq_varid PHB')
   call nc_check(nf90_get_var(ncid, var_id, grid_data%phb), &
                 'read_model_static_data', 'get_var PHB')
   if (debug) then
      write (*, *) ' corners of phb '
      write (*, *) grid_data%phb(1, 1, 1), grid_data%phb(grid_data%we, 1, 1), &
         grid_data%phb(1, grid_data%sn, 1), grid_data%phb(grid_data%we, &
                                                          grid_data%sn, 1)
      write (*, *) grid_data%phb(1, 1, grid_data%bt), &
         grid_data%phb(grid_data%we, 1, grid_data%bt), &
         grid_data%phb(1, grid_data%sn, grid_data%bt), &
         grid_data%phb(grid_data%we, grid_data%sn, grid_data%bt)
   end if

end subroutine read_model_static_data

!--------------------------------------------

!-----------------------------------------------------------------------
!>
!> Fill the array of requested variables, dart kinds, possible min/max
!> values and whether or not to update the field in the output file.
!> Then calls 'add_domain()' to tell the DART code which variables to
!> read into the state vector after this code returns.
!>
!>@param variable_array  the list of variables and kinds from model_mod_nml
!>@param nfields         the number of variable/Quantity pairs specified

subroutine set_model_variable_info(input_template_filename, variable_array)

   character(len=*), intent(in)  :: input_template_filename
   character(len=*), intent(in)  :: variable_array(:)

   character(len=*), parameter :: routine = 'set_model_variable_info:'

   integer :: i, nfields
   integer, parameter :: MAX_STRING_LEN = 128

   character(len=MAX_STRING_LEN) :: varname    ! column 1, NetCDF variable name
   character(len=MAX_STRING_LEN) :: dartstr    ! column 2, DART Quantity
   character(len=MAX_STRING_LEN) :: minvalstr  ! column 3, Clamp min val
   character(len=MAX_STRING_LEN) :: maxvalstr  ! column 4, Clamp max val
   character(len=MAX_STRING_LEN) :: updatestr  ! column 5, Update output or not

   character(len=vtablenamelength) :: var_names(MAX_STATE_VARIABLES) = ' '
   logical  :: update_list(MAX_STATE_VARIABLES) = .FALSE.
   integer  ::   kind_list(MAX_STATE_VARIABLES) = MISSING_I
   real(r8) ::  clamp_vals(MAX_STATE_VARIABLES, 2) = MISSING_R8

   nfields = 0
   ParseVariables: do i = 1, MAX_STATE_VARIABLES

      varname = variable_array(num_state_table_columns*i - 4)
      dartstr = variable_array(num_state_table_columns*i - 3)
      minvalstr = variable_array(num_state_table_columns*i - 2)
      maxvalstr = variable_array(num_state_table_columns*i - 1)
      updatestr = variable_array(num_state_table_columns*i)

      if (varname == ' ' .and. dartstr == ' ') exit ParseVariables ! Found end of list.

      if (varname == ' ' .or. dartstr == ' ') then
         string1 = 'model_nml:model "state_variables" not fully specified'
         call error_handler(E_ERR, routine, string1, source)
      end if

      ! Make sure DART kind is valid

      if (get_index_for_quantity(dartstr) < 0) then
         write (string1, '(3A)') 'there is no obs_kind "', trim(dartstr), '" in obs_kind_mod.f90'
         call error_handler(E_ERR, routine, string1, source)
      end if

      call to_upper(minvalstr)
      call to_upper(maxvalstr)
      call to_upper(updatestr)

      var_names(i) = varname
      kind_list(i) = get_index_for_quantity(dartstr)
      clamp_vals(i, 1) = string_to_real(minvalstr)
      clamp_vals(i, 2) = string_to_real(maxvalstr)
      update_list(i) = string_to_logical(updatestr, 'UPDATE')

      nfields = nfields + 1

   end do ParseVariables

   if (nfields == MAX_STATE_VARIABLES) then
      write (string1, '(2A)') 'WARNING: There is a possibility you need to increase ', &
         'MAX_STATE_VARIABLES in the global variables in model_mod.f90'

      write (string2, '(A,i4,A)') 'WARNING: you have specified at least ', nfields, &
         ' perhaps more'

      call error_handler(E_MSG, routine, string1, source, text2=string2)
   end if

   ! CAM only has a single domain (only a single grid, no nests or multiple grids)

   domain_id = add_domain(input_template_filename, nfields, var_names, kind_list, &
                          clamp_vals, update_list)

end subroutine set_model_variable_info

subroutine fill_default_state_table(default_table)

   character(len=129), intent(out) :: default_table(num_state_table_columns, MAX_STATE_VARIABLES)

   integer :: row

   default_table = 'NULL'
   row = 0

   ! fill default state variable table here.
   row = row + 1
   default_table(:, row) = (/'U                         ', &
                             'QTY_U_WIND_COMPONENT      ', &
                             'NA                        ', &
                             'NA                        ', &
                             'UPDATE                    '/)
   row = row + 1
   default_table(:, row) = (/'V                         ', &
                             'QTY_V_WIND_COMPONENT      ', &
                             'NA                        ', &
                             'NA                        ', &
                             'UPDATE                    '/)
   row = row + 1
   default_table(:, row) = (/'T                         ', &
                             'QTY_TEMPERATURE           ', &
                             'NA                        ', &
                             'NA                        ', &
                             'UPDATE                    '/)
   row = row + 1
   default_table(:, row) = (/'Q                          ', &
                             'QTY_SPECIFIC_HUMIDITY      ', &
                             'NA                         ', &
                             'NA                         ', &
                             'UPDATE                     '/)

   return

end subroutine fill_default_state_table

!--------------------------------------------
!--------------------------------------------

integer function get_number_of_state_variables(variable_array)

   character(len=*), intent(in) :: variable_array(num_state_table_columns &
                                                  *MAX_STATE_VARIABLES)

   integer :: num_vars, i
   logical :: debug = .false.

   integer, parameter :: MAX_STRING_LEN = 128

   character(len=MAX_STRING_LEN) :: varname    ! column 1, NetCDF variable name
   character(len=MAX_STRING_LEN) :: dartstr    ! column 2, DART Quantity

   num_vars = 0

   ParseVariables: do i = 1, MAX_STATE_VARIABLES
      varname = variable_array(num_state_table_columns*i - 4)
      dartstr = variable_array(num_state_table_columns*i - 3)
      if (varname == ' ' .and. dartstr == ' ') exit ParseVariables ! Found end of list.
      num_vars = num_vars + 1
   end do ParseVariables

   get_number_of_state_variables = num_vars

   return

end function get_number_of_state_variables

!--------------------------------------------------------------------
!--------------------------------------------------------------------

subroutine get_variable_metadata_from_file(ncid, wrf_var_name, description, units)
   !subroutine get_variable_metadata_from_file(ncid,wrf_var_name,description, coordinates,units)

   ! ncid: input, file handle
   ! id:   input, domain index

   integer, intent(in)               :: ncid
   character(len=*), intent(in)    :: wrf_var_name
   !character(len=129), intent(out)   :: description, coordinates, units
   character(len=129), intent(out)   :: description, units

   logical, parameter    :: debug = .false.
   integer               :: var_id

   call nc_check(nf90_inq_varid(ncid, trim(wrf_var_name), var_id), &
                 'get_variable_metadata_from_file', &
                 'inq_varid '//wrf_var_name)

   description = ''
   call nc_check(nf90_get_att(ncid, var_id, 'description', description), &
                 'get_variable_metadata_from_file', &
                 'get_att '//wrf_var_name//' '//description)

   !coordinates = ''
   !call nc_check( nf90_get_att(ncid, var_id, 'coordinates', coordinates), &
   !                        'get_variable_metadata_from_file', &
   !                        'get_att '//wrf_var_name//' '//coordinates)

   units = ''
   call nc_check(nf90_get_att(ncid, var_id, 'units', units), &
                 'get_variable_metadata_from_file', &
                 'get_att '//wrf_var_name//' '//units)

   return

end subroutine get_variable_metadata_from_file

        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!    height_diff_check - function that determines whether a pair of heights
!                        in meters are closer than the given tolerance.
!                        returns .TRUE. if closer or equal to limit
!
!    max_diff_meters   - maximum difference between 2 elevations (m)
!    height1           - first height (m)
!    height2           - second height (m)
!
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
function height_diff_check(max_diff_meters, height1, height2)
   real(r8), intent(in) :: max_diff_meters, height1, height2
   logical              :: height_diff_check

   height_diff_check = .true.

   if (abs(height1 - height2) > max_diff_meters) height_diff_check = .false.

end function height_diff_check

!-----------------------------------------------------------------------
!> return my_status /= 0 if obs is above a user-defined threshold.
!> intended to be quick (low-cost) and not exact.
!> This intentionally does NOT have a case for vert type of
!> SCALEHEIGHT - because this routine is only used to look at
!> observation locations.  we have not yet encountered obs
!> with that vertical type.

subroutine obs_too_high(vert_value, which_vert, my_status)
   real(r8), intent(in) :: vert_value
   integer, intent(in) :: which_vert
   integer, intent(out) :: my_status

   ! assume ok to begin with
   my_status = 0

   if (which_vert == VERTISPRESSURE) then
      ! lower pressures are higher; watch the less than/greater than tests
      if (vert_value < no_assim_above_pressure) my_status = 14
      return
   end if

   ! these are always ok
   if (which_vert == VERTISSURFACE .or. which_vert == VERTISUNDEF) return

   if (which_vert == VERTISHEIGHT) then
      if (vert_value > no_assim_above_height) my_status = 14
      return
   end if

   if (which_vert == VERTISLEVEL) then
      ! level 1 is bottom; watch less than/greater than in tests
      if (vert_value > no_assim_above_level) my_status = 14
      return
   end if

   ! for now we haven't run into observations where the vertical coordinate
   ! (of the OBS) is in scale height - but if we do it will fall into here.

   write (errstring, *) 'vertical type: ', which_vert
   call error_handler(E_ERR, 'obs_too_high', 'unrecognized vertical type', &
                      source, text2=errstring)

end subroutine obs_too_high

!-----------------------------------------------------------------------
!> return 0 (ok) if we know how to interpolate this quantity.
!> if it is a field in the state, return the variable id from
!> the state structure.  if not in the state, varid will return -1

subroutine ok_to_interpolate(obs_qty, varid, my_status)
   integer, intent(in)  :: obs_qty
   integer, intent(out) :: varid
   integer, intent(out) :: my_status

   ! See if the state contains the obs quantity
   varid = get_varid_from_kind(domain_id, obs_qty)

   ! in the state vector
   if (varid > 0) then
      my_status = 0
      return
   end if

   ! add any quantities that can be interpolated to this list if they
   ! are not in the state vector.
   select case (obs_qty)
   case (QTY_SURFACE_ELEVATION, &
         QTY_PRESSURE, &
         QTY_GEOMETRIC_HEIGHT, &
         QTY_VERTLEVEL)
      my_status = 0
   case default
      my_status = 2
   end select

end subroutine ok_to_interpolate

!-----------------------------------------------------------------------
!> convert from string to integer, and set in the dart code the
!> vertical type we are going to want to localize in.

subroutine set_vert_localization(typename)
   character(len=*), intent(in)  :: typename
   character(len=*), parameter :: routine = 'set_vert_localization'

   character(len=32) :: ucasename
   integer :: vcoord

   ucasename = typename
   call to_upper(ucasename)

   select case (ucasename)
   case ("PRESSURE")
      vcoord = VERTISPRESSURE
   case ("HEIGHT")
      vcoord = VERTISHEIGHT
   case ("SCALEHEIGHT", "SCALE_HEIGHT", "SCALE HEIGHT")
      vcoord = VERTISSCALEHEIGHT
   case ("LEVEL", "MODEL_LEVEL", "MODEL LEVEL")
      vcoord = VERTISLEVEL
   case default
      write (string1, *) 'unrecognized vertical localization coordinate type: '//trim(typename)
      write (string2, *) 'valid values are: PRESSURE, HEIGHT, SCALEHEIGHT, LEVEL'
      call error_handler(E_ERR, routine, string1, source, text2=string2)
   end select

   ! during assimilation, when get_close() is called to compute the separation distance
   ! between items, convert all state and obs to use this vertical type if vertical localization
   ! is enabled (usually true for cam).

   call set_vertical_localization_coord(vcoord)

   ! save in module global for later use.
   vertical_localization_type = vcoord

end subroutine set_vert_localization

function read_model_time(filename)

   character(len=*), intent(in) :: filename
   type(time_type)              :: read_model_time

   integer :: ncid
   integer :: litime, sub_date, sub_tod
   character(len=19) :: timestring
   integer :: iyear, imonth, iday, ihour, iminute, isecond, rem

   character(len=*), parameter :: routine = 'read_model_time'

   !SENote: Doesn't actually need model to be initialized
   !if ( .not. module_initialized ) call static_init_model

   if (.not. file_exist(filename)) then
      write (string1, *) trim(filename), ' does not exist.'
      call error_handler(E_ERR, routine, string1, source)
   end if

   ncid = nc_open_file_readonly(filename, routine)

   ! CAM initial files have two variables of length
   ! 'time' (the unlimited dimension): date, datesec

   call nc_get_variable(ncid, 'time', litime, routine)
   call nc_get_variable(ncid, 'date', sub_date, routine)

   write (timestring, '(I0)') litime

   call nc_get_variable(ncid, 'datesec', sub_tod, routine)

   !'date' is YYYYMMDD
   !'sub_tod' is seconds of current day
   iyear = sub_date/10000
   rem = sub_date - iyear*10000
   imonth = rem/100
   iday = rem - imonth*100

   ihour = sub_tod/3600
   rem = sub_tod - ihour*3600
   iminute = rem/60
   isecond = rem - iminute*60

   ! some cam files are from before the start of the gregorian calendar.
   ! since these are 'arbitrary' years, just change the offset.
   if (iyear < 1601) then
      write (string1, *) ' '
      write (string2, *) 'WARNING - ', trim(filename), ' changing year from ', &
         iyear, 'to', iyear + 1601

      call error_handler(E_MSG, routine, string1, source, &
                         text2=string2, text3='to make it a valid Gregorian date.')

      write (string1, *) ' '
      call error_handler(E_MSG, routine, string1, source)
      iyear = iyear + 1601
   end if

   read_model_time = set_date(iyear, imonth, iday, ihour, iminute, isecond)

   call nc_close_file(ncid, routine)

end function read_model_time

subroutine write_model_time(ncid, model_time)
   integer, intent(in) :: ncid
   type(time_type), intent(in) :: model_time

   integer :: iyear, imonth, iday, ihour, iminute, isecond
   integer :: sub_date(1), sub_tod(1), sub_time(1)

   character(len=*), parameter :: routine = 'write_model_time'

   call get_date(model_time, iyear, imonth, iday, ihour, iminute, isecond)

   sub_date = iyear*10000 + imonth*100 + iday
   sub_tod = ihour*3600 + iminute*60 + isecond
   sub_time = sub_date + sub_tod/86400

   if (.not. nc_variable_exists(ncid, "time")) then
      call nc_begin_define_mode(ncid, routine)
      call nc_define_real_variable(ncid, 'time', (/'time'/), routine)
      call nc_end_define_mode(ncid, routine)
      call nc_put_variable(ncid, 'time', sub_time, routine)
   end if

   ! if the file doesn't already have a "date" variable make one
   if (.not. nc_variable_exists(ncid, "date")) then
      call nc_begin_define_mode(ncid, routine)
      call nc_define_integer_variable(ncid, 'date', (/'time'/), routine)
      call nc_end_define_mode(ncid, routine)
      call nc_put_variable(ncid, 'date', sub_date, routine)
   end if

   ! if the file doesn't already have a "datesec" variable make one
   if (.not. nc_variable_exists(ncid, "datesec")) then
      call nc_begin_define_mode(ncid, routine)
      call nc_define_integer_variable(ncid, 'datesec', (/'time'/), routine)
      call nc_end_define_mode(ncid, routine)
      call nc_put_variable(ncid, 'datesec', sub_tod, routine)
   end if

end subroutine write_model_time

!--------------------------------------------------------------------
!> pressure gets smaller as you go up, everything else gets larger.
!> return true if this value is above the start of the ramp.
!> test_value and ramp_end need to already be in vert localization units

! FIXME: test this new code section carefully.
!
! right now the calling code is expecting extra_dist to be added
! to the original get_dist() value, so any scaling or modifications
! should happen in this routine.
!
! do we need the 2 locations here to compute the horizontal distance?
! or is having the total dist and the vertical separation enough?

function above_ramp_start(test_value, gc, obs_type, ramp_end, total_dist, extra_dist)
   real(r8), intent(in)  :: test_value
   type(get_close_type), intent(in)  :: gc
   integer, intent(in)  :: obs_type
   real(r8), intent(in)  :: ramp_end
   real(r8), intent(in)  :: total_dist
   real(r8), intent(out) :: extra_dist
   logical :: above_ramp_start

   real(r8) :: vert_localize_dist, ramp_start, norm, vert_norm, vert_only_dist
   real(r8) :: horiz_dist, ramp_dist, ramp_width
   type(location_type) :: this_loc, ramp_start_loc, loc1, loc2
   logical, save :: onetime = .true.

   ! do the easy cases first - either above the ramp end
   ! or below the ramp start.  leave the middle ground for
   ! last because we have to then compute a damping factor.

   ! FIXME: test this!!!
   ! is it above the ramp end? set damp dist to something
   ! large enough to turn off all impacts.  is vert_localize_dist enough?
   vert_localize_dist = get_maxdist(gc, obs_type)
   if (.false. .and. onetime) then
      print *, 'vert_localize_dist = ', vert_localize_dist
      onetime = .false.
   end if

   if (v_above(test_value, ramp_end)) then
      extra_dist = vert_localize_dist
      above_ramp_start = .true.
      return
   end if

   ! compute ramp start and see if we're lower than that.

   ! vert norm for this obs type
   loc1 = set_location(0.0_r8, 0.0_r8, 0.0_r8, vertical_localization_type)
   loc2 = set_location(0.0_r8, 0.0_r8, 1.0_r8, vertical_localization_type)
   norm = get_dist(loc1, loc2, obs_type)   ! units: rad/loc units
   vert_norm = 1.0_r8/norm               ! units now: loc units/rad

   ramp_start = v_down(ramp_end, vert_norm*vert_localize_dist)

   !print *, 'computing ramp start: ramp_end, vert_norm, vert_localize_dist', &
   !                    ramp_start, ramp_end, vert_norm, vert_localize_dist

   if (.not. v_above(test_value, ramp_start)) then
      extra_dist = 0.0_r8
      above_ramp_start = .false.
      return
   end if

   ! ok, we're somewhere inbetween.  compute horiz and vert distances
   ! and see what the ramping factor needs to be.

   !print *, 'test value within ramp range: ', ramp_start, test_value, ramp_end
   above_ramp_start = .true.

   ! see what the vertical separation is from obs to start of ramp
   this_loc = set_location(0.0_r8, 0.0_r8, test_value, vertical_localization_type)
   ramp_start_loc = set_location(0.0_r8, 0.0_r8, ramp_start, vertical_localization_type)

   ! do we need this?  i think so.   radians
   vert_only_dist = get_dist(ramp_start_loc, this_loc, obs_type)

   ! we need this to compute what?
   if (vert_only_dist > total_dist) then
      !print *, 'unexpected, vert larger than total:  ', vert_only_dist, total_dist
      !print *, 'obs_type, vert_norm = ', obs_type, vert_norm
      horiz_dist = 0.0_r8
   else
      horiz_dist = sqrt(total_dist**2 - vert_only_dist**2)
   end if

   ramp_dist = v_difference(test_value, ramp_start)
   ramp_width = v_difference(ramp_end, ramp_start)
   extra_dist = (ramp_dist/ramp_width)*vert_localize_dist

   ! DEBUG - disable for now
   if (.false. .and. above_ramp_start) then
      print *, 'ramp s/v/e: ', ramp_start, test_value, ramp_end
      print *, 'v, h:       ', vert_only_dist, horiz_dist
      print *, 'rampd, tot: ', ramp_dist, ramp_width
      print *, 'ed, return: ', extra_dist, above_ramp_start
   end if

end function above_ramp_start

!--------------------------------------------------------------------
! returns true if a is above b (higher in the atmosphere,
! further from the surface of the earth).

pure function v_above(a, b)
   real(r8), intent(in) :: a, b
   logical :: v_above

   if (higher_is_smaller) then
      v_above = (a < b)
   else
      v_above = (a > b)
   end if

end function v_above

!--------------------------------------------------------------------
! returns new value of moving b distance down in the atmosphere
! starting at a.  for height, this results in a smaller value
! (also one flavor of scale height), but for other vertical types
! this results in a larger value.

pure function v_down(a, b)
   real(r8), intent(in) :: a, b
   real(r8) :: v_down

   if (higher_is_smaller) then
      v_down = (a + b)
   else
      v_down = (a - b)
   end if

end function v_down

!--------------------------------------------------------------------
! returns difference of a and b
! (doesn't depend on the vertical_localization_type)

pure function v_difference(a, b)
   real(r8), intent(in) :: a, b
   real(r8) :: v_difference

   v_difference = abs(a - b)

end function v_difference

! add any 2d fields here that are surface quantities
function is_surface_field(qty)
   integer, intent(in) :: qty
   logical :: is_surface_field

   select case (qty)
   case (QTY_SURFACE_PRESSURE, QTY_SURFACE_ELEVATION)
      is_surface_field = .true.

   case default
      is_surface_field = .false.

   end select
end function is_surface_field
!--------------------------------------------------------------------
! Function to calculate scale height given a pressure and optionally
! a surface pressure.  (See the namelist item which controls whether to
! normalize the pressure value aloft with the surface pressure or not.
! We currently only use scale height for computing distances between
! two locations, so the surface pressure terms cancel out - exactly if
! the two locations are co-located horizontally, almost if they are not.
! Normalizing by the surface pressure means in areas of high orography
! the surface differences propagate all the way to the model top.
! To be backwards-compatible, do this normalization; the current thinking
! is we shouldn't do it both for scientific reasons and because it
! doubles the work if it's expensive to find the correct horizontal
! location, i.e. mpas irregular grids. In this model we always have
! the surface pressure at a location so it's not a performance issue.)
!
! Watch out for unusual cases that could crash the log() function
! We pass in the surface pressure here even if it isn't going to be
! used because in all the cases above we seem to have it (or the standard
! reference pressure) everywhere we are going to compute this value.
! The "skip_norm" parameter controls whether this code uses the
! surface pressure or not.

function scale_height(p_above, p_surface, skip_norm)
   real(r8), intent(in) :: p_above
   real(r8), intent(in) :: p_surface
   logical, intent(in) :: skip_norm
   real(r8)             :: scale_height

   real(r8), parameter :: tiny = epsilon(1.0_r8)
   real(r8) :: diff

   if (skip_norm) then
      scale_height = log(p_above)
      return
   end if

   diff = p_surface - p_above  ! should be positive

   if (abs(diff) < tiny) then
      ! surface obs will have (almost) identical values
      scale_height = 0.0_r8   ! -log(1.0_r8)

   else if (diff <= tiny .or. p_above <= 0.0_r8) then
      ! weed out bad cases
      scale_height = MISSING_R8

   else
      ! normal computation - should be safe now
      scale_height = -log(p_above/p_surface)

   end if

end function scale_height
!--------------------------------------------------------------------
!> if the namelist is set to not use this custom routine, the default
!> dart routine will add 'pert_amp' of noise to every field in the state
!> to generate an ensemble from a single member.  if it is set to true
!> this routine will be called.  the pert_amp will be ignored, and the
!> given list of quantities will be perturbed by the given amplitude
!> (which can be different for each field) to generate an ensemble.

!===================================================================
! End of model_mod
!===================================================================
end module model_mod

