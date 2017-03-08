! This code may (or may not) be part of the COAMPS distribution,
! So it is not protected by the DART copyright agreement.
!
! DART $Id$

module coamps_interp_mod

!------------------------------
! MODULE:       coamps_interp_mod
! AUTHOR:       T. R. Whitcomb
!               Naval Research Laboratory
! DART VERSION: Jamaica
!
! Module to handle interpolating the DART state vector to a point.
! This allows the use of observations that are not defined at 
! model grid points
!------------------------------ 
  
  use coamps_grid_mod,      only : check_ij_within_grid,           &
                                   coamps_grid,                    &
                                   dump_grid_info,                 &
                                   get_grid_delta_x,               &
                                   get_grid_delta_y,               &
                                   get_grid_dims,                  &
                                   get_grid_dsigmaw,               &
                                   get_grid_field_size,            & 
                                   get_grid_msigma,                &
                                   get_grid_num_levels,            &
                                   get_grid_wsigma,                &
                                   get_terrain_height_at_points,   &
                                   location_to_gridpt
  use coamps_intrinsic_mod, only : s2pint,                         &
                                   seaprs,                         &
                                   utom,                           &
                                   vtom,                           &
                                   z2zint
  use coamps_restart_mod,   only : get_num_vars,                   &
                                   get_restart_index_by_properties
  use coamps_util_mod,      only : check_alloc_status,             &
                                   check_dealloc_status,           &
                                   print_label_name,               &
                                   print_2d_real8_array

  use location_mod,       only : location_type,    &
                                 vert_is_height,   &
                                 vert_is_level,    &
                                 vert_is_pressure, &
                                 vert_is_surface
  use obs_kind_mod
  use types_mod,          only : MISSING_I,        &
                                 MISSING_R8,       &
                                 RAD2DEG,          &
                                 r8
  use utilities_mod,      only : do_output,        &
                                 E_ERR,            &
                                 E_MSG,            &
                                 error_handler,    &
                                 register_module 

  implicit none

  private

  !------------------------------
  ! BEGIN PUBLIC INTERFACE
  !------------------------------

  public :: interpolate

  !------------------------------
  ! END PUBLIC INTERFACE
  !------------------------------

  !------------------------------
  ! BEGIN EXTERNAL INTERFACES
  !------------------------------
  ! [none]
  !------------------------------
  ! END EXTERNAL INTERFACES
  !------------------------------

  !------------------------------
  ! BEGIN TYPES AND CONSTANTS 
  !------------------------------

  ! Atmospheric constants - specific heat, gas constant (dry air),
  ! and the initial pressure for calculating the Exner function
  real(kind=r8), parameter :: R   = real(287.0,  kind=r8)
  real(kind=r8), parameter :: Cp  = real(1004.0, kind=r8)
  real(kind=r8), parameter :: P00 = real(1000.0, kind=r8)

  ! For the height/sigma vertical interpolation, we have the option
  ! of either setting below-ground levels to the value of the lowest
  ! level or we can set the value as "missing" - pick the second one
  real(kind=r8), parameter :: MISSING_VALUE     = MISSING_R8
  integer,       parameter :: USE_MISSING_VALUE = 1

  ! Ordering of neighboring points
  integer, parameter :: NUM_NEIGHBORS        = 4
  integer, parameter :: NEIGHBOR_LOWER_LEFT  = 1
  integer, parameter :: NEIGHBOR_UPPER_LEFT  = 2
  integer, parameter :: NEIGHBOR_UPPER_RIGHT = 3
  integer, parameter :: NEIGHBOR_LOWER_RIGHT = 4

  ! The minimum number of levels needed for an interpolation
  integer, parameter :: MIN_LEVELS_NEEDED = 5

  ! The only real supported interpolation modes at this time are
  ! pressure, height, and sigma level - surface obs like SLP are
  ! handled separately
  integer, parameter :: INTERPOLATE_TO_PRESSURE = 1
  integer, parameter :: INTERPOLATE_TO_HEIGHT   = 2
  integer, parameter :: INTERPOLATE_TO_SIGMA    = 3
  integer, parameter :: INTERPOLATE_TO_SURFACE  = 4
  integer, parameter :: INTERPOLATE_TO_OTHER    = 5

  ! Pressure conversion constant
  real(kind=r8), parameter :: CONVERT_PA_TO_MB = real(0.01, kind=r8)

  ! Number of things we need to pull out of the state vector for the
  ! various types of interpolation.  
  !  Symbol Key:
  !   [entity]          entity is guaranteed to be present
  !   #entity#          entity is calculated during the course of
  !                     interpolation
  !   %entity%          entity is part of the model mean state
  ! Note that currently the variables in brackets are the *only*
  ! entries that are guaranteed to be present.
  !  Pressure:          target variable
  !                     #pressure#
  !                     perturbation Exner function
  !                     %mean Exner function%
  !  Height:            target variable
  !                     #geometric height#
  !                     [mass sigma levels]
  !                     [w sigma levels]
  !                     [terrain height]
  !  Sigma:             target variable
  !                     [mass sigma levels]
  !  SLP:               [mass sigma levels]
  !                     [w sigma levels]
  !                     potential temperature
  !                     %mean potential temperature%
  !                     perturbation Exner function
  !                     %mean Exner function (mass levels)%
  !                     %mean Exner function (w levels)%
  !                     [dsigmaw]
  integer, parameter :: NUM_VARS_NEEDED_FOR_P_INTERP = 4
  integer, parameter :: NUM_VARS_NEEDED_FOR_Z_INTERP = 2
  integer, parameter :: NUM_VARS_NEEDED_FOR_S_INTERP = 2
  integer, parameter :: NUM_VARS_NEEDED_FOR_SLP      = 7

  ! Alias the dimension types for the 2-D arrays so
  ! I can keep them straight: 
  !  "values" arrays index as (neighbor, level)
  !  availability array index as (level, variable)
  integer, parameter :: VALUES_DIM_NEIGHBOR       = 1
  integer, parameter :: VALUES_DIM_LEVEL          = 2
  integer, parameter :: AVAILABILITY_DIM_LEVEL    = 1
  integer, parameter :: AVAILABILITY_DIM_VARIABLE = 2
  !------------------------------
  ! END TYPES AND CONSTANTS 
  !------------------------------

  !------------------------------
  ! BEGIN MODULE VARIABLES
  !------------------------------

! version controlled file description for error handling, do not edit
character(len=256), parameter :: source   = &
   "$URL$"
character(len=32 ), parameter :: revision = "$Revision$"
character(len=128), parameter :: revdate  = "$Date$"

  ! Module-level aliases - allow us to access these items passed
  ! in to the driver function throughout the module
  real(kind=r8), dimension(:), pointer :: state_vector
  type(coamps_grid),           pointer :: interp_grid

  ! Store a column of values at each of the neighboring points
  ! These arrays for the interpolated variable values and the 
  ! corresponding model pressure are defined on sigma levels:
  !  (neighbor index, sigma level index)
  real(kind=r8), dimension(:,:), allocatable :: col_values
  real(kind=r8), dimension(:,:), allocatable :: vcoord_values
  real(kind=r8), dimension(:,:), allocatable :: vinterp_values
  real(kind=r8), dimension(:),   allocatable :: vinterp_levels
  real(kind=r8), dimension(:,:), allocatable :: mean_pres_m_values
  real(kind=r8), dimension(:,:), allocatable :: pert_pres_m_values
  real(kind=r8), dimension(:,:), allocatable :: theta_values
  real(kind=r8), dimension(:,:), allocatable :: thetabar_values
  real(kind=r8), dimension(:,:), allocatable :: mean_pres_w_values
  
  ! Store the column of values at each of the neighboring points
  ! but only at the levels that we have data available - this is
  ! a more general form since we only have two entries to worry
  ! about here
  real(kind=r8), dimension(:,:), allocatable :: available_target_vals
  real(kind=r8), dimension(:,:), allocatable :: available_vcoord_vals

  ! The coordinate variables can be either directly generated or 
  ! pulled directly out of the state vector.  For the target variable
  ! though, read it into a 2D field to make destaggering easier
  real(kind=r8), dimension(:,:), allocatable :: var_field
  
  ! Vertical coordinate information for both the COAMPS interp
  ! functions and conversion between vertical coordinate types
  real(kind=r8), dimension(:),   allocatable :: terrain_values
  real(kind=r8), dimension(:),   allocatable :: msigma_values
  real(kind=r8), dimension(:),   allocatable :: wsigma_values
  real(kind=r8), dimension(:),   allocatable :: dsigmaw_values
  
  ! Size of grid we'll be interpolating on - this includes both the
  ! x/y dimensions and the number of levels (both the total number
  ! of levels defined in the model and the number of levels the
  ! variable we're interpolating is defined on - these are not 
  ! guaranteed to be equal
  integer :: grid_i, grid_j
  integer :: field_size
  integer :: num_model_levels

  ! How many variables are in the state vector
  integer :: total_num_vars

  ! Data about the observation - note that the vertical direction is
  ! a "z", since the vertical coordinate is most likely pressure or
  ! height while the horizontal direction is i/j since that is 
  ! stored as non-integer grid indices, not actual x/y|lon/lat locs
  type(location_type) :: obs_loc
  real(kind=r8)       :: obs_loc_i
  real(kind=r8)       :: obs_loc_j
  real(kind=r8)       :: obs_loc_z

  ! Coordinates of the grid points near the observation
  integer, dimension(NUM_NEIGHBORS) :: neighbors_i
  integer, dimension(NUM_NEIGHBORS) :: neighbors_j

  ! What kind of level are we interpolating to?
  integer :: interpolate_to_level_type

  ! Availability matrix - this is an (# of levels) x (# of vars)
  ! array - if the value at (level, variable) is .true., then that
  ! variable is available at that level.  This is done to handle the
  ! case that variables may not be available at all levels and/or 
  ! different variables have different availabilities.  This also 
  ! allows for (accidental) duplication  of entries in the state
  ! vector definition - but beware that the later entry will over-
  ! write the value of that variable at that level.
  logical, dimension(:,:), allocatable :: vars_available
  logical, dimension(:), allocatable   :: levels_available
  integer                              :: num_levels_available

  ! Select out only the available levels
  logical, dimension(:,:), allocatable :: level_mask

  ! Since we don't care about keeping track of what's available 
  ! where once we've established the availability (once everything
  ! is read in we only care about which levels have all parts
  ! present and accounted for), use a counter to keep track of
  ! which column to put the availability of whatever we're reading
  ! into
  integer :: availability_index
  
  ! Strictly speaking, these are set at runtime.  However, since
  ! they are aliases for array indices, and don't "change", they
  ! should be treated like constants
  integer :: AVAILABLE_INDEX_TARGET
  integer :: AVAILABLE_INDEX_MEAN_EXNER
  integer :: AVAILABLE_INDEX_PERT_EXNER
  integer :: AVAILABLE_INDEX_PRESSURE
  integer :: AVAILABLE_INDEX_HEIGHT
  integer :: AVAILABLE_INDEX_SIGMA

  ! Weights for the interpolation
  real(kind=r8), dimension(NUM_NEIGHBORS) :: interp_weights

  logical, save :: module_initialized = .false.

  !------------------------------
  ! END MODULE VARIABLES
  !------------------------------

contains

  !------------------------------
  ! BEGIN PUBLIC ROUTINES
  !------------------------------

  ! interpolate
  ! -----------
  ! Driver for interpolation calculation
  !  PARAMETERS
  !   IN  state             big ol' DART state vector
  !   IN  grid              coamps_grid structure to interpolate on
  !   IN  obs_loc           DART location structure to interpolate to
  !   IN  obs_type          integer version of raw variable type
  !   OUT obs_value         result of interpolation
  !   OUT interp_worked     true if interpolation was successful
  subroutine interpolate(state, grid, obs_loc, obs_type, obs_value,&
                         interp_worked)
    real(kind=r8), dimension(:), target, intent(in)  :: state
    type(coamps_grid),           target, intent(in)  :: grid
    type(location_type),                 intent(in)  :: obs_loc
    integer,                             intent(in)  :: obs_type
    real(kind=r8),                       intent(out) :: obs_value
    logical,optional,                    intent(out) :: interp_worked

    ! Assume that the interpolation could go wrong in several places.
    ! Check them all - if they didn't work, throw up our hands and
    ! let the calling function deal with us not working properly
    logical :: obs_within_grid
    logical :: valid_level_type
    logical :: enough_levels_available
    logical :: level_in_range
    logical :: no_missing_values

    ! Module-level aliases
    interp_grid  => grid
    state_vector => state

    if (.not. module_initialized) then
      call initialize_module()
    end if

    call location_to_gridpt(grid, obs_loc, obs_loc_i, obs_loc_j, &
                            obs_loc_z)
    call check_ij_within_grid(grid, ceiling(obs_loc_i),          &
                              ceiling(obs_loc_j), obs_within_grid)
    if (.not. obs_within_grid) then
      interp_worked = .false.
      return
    end if

    ! Whether obs is on pressure, height, etc.
    call set_interpolation_level_type(obs_loc)
    call check_level_type(valid_level_type)
    if (.not. valid_level_type) then
      interp_worked = .false.
      return
    end if

    call allocate_availability_data()

    ! Reset since every call to interpolate() means we are interpolating
    ! new value of a potentially different type
    availability_index = 0

    call get_neighbors(obs_loc_i, obs_loc_j, neighbors_i, neighbors_j) 

    ! Terrain heights and sigma levels
    call get_vertical_data(interp_grid)

    ! The target variable is defined by the observation type and the
    ! coordinate variables are defined by the vertical level type
    ! in the observation location supplied to the routine. 
    call get_target_var(obs_type)
    call get_coordinate_vars()

    ! Now we can figure out how data points we have to work with
    call calculate_available_levels()
    call check_enough_levels_available(enough_levels_available)
    if (.not. enough_levels_available) then
      interp_worked = .false.
      return
    end if

    call allocate_available_values()
    call generate_level_mask()
    call get_available_values(col_values,    available_target_vals)
    call get_available_values(vcoord_values, available_vcoord_vals)

    ! Make sure that the vertical level we are trying to interpolate
    ! to is within the vertical levels that we actually have
    call check_level_in_available_range(level_in_range)
    if (.not. level_in_range) then
      interp_worked = .false.
      return
    end if

    ! Interpolate onto a single vertical level - we probably should
    ! not end up with any missing values since we already checked
    ! that the level was within range.
    call interpolate_to_level()
    call check_for_missing_values(no_missing_values)
    if (.not. no_missing_values) then
      interp_worked = .false.
      return
    end if

    ! Interpolate to observation location
    call calculate_interp_weights()
    call interpolate_to_point(obs_value) 
    interp_worked = .true.

    call print_interpolation_diagnostics(obs_type, obs_value)

    call deallocate_availability_data()
    call deallocate_available_values()

    nullify(state_vector)
    nullify(interp_grid)
  end subroutine interpolate
  
  !------------------------------
  ! END PUBLIC ROUTINES
  !------------------------------

  !------------------------------
  ! BEGIN PRIVATE ROUTINES
  !------------------------------

  ! get_neighbors
  ! -------------
  ! Given a i/j location, returns the i/j information of the four
  ! neighbor points.  This assumes that the i/j coordinate supplied
  ! is within the bounds of the grid.
  !  PARAMETERS
  !   IN  loc_ii            target x-coordinate index
  !   IN  loc_jj            target y-coordinate index
  !   OUT adj_ii            neighbors's x-coordinate indices
  !   OUT adj_jj            neighbors's y-coordinate indices
  subroutine get_neighbors(loc_ii, loc_jj, adj_ii, adj_jj)
    real(kind=r8), intent(in)                      :: loc_ii
    real(kind=r8), intent(in)                      :: loc_jj
    integer, dimension(NUM_NEIGHBORS), intent(out) :: adj_ii
    integer, dimension(NUM_NEIGHBORS), intent(out) :: adj_jj

    integer :: jj_lower, jj_upper
    integer :: ii_left, ii_right

    ii_left  = floor(loc_ii)
    ii_right = ceiling(loc_ii)
    jj_lower = floor(loc_jj)
    jj_upper = ceiling(loc_jj)

    ! Ordering goes clockwise starting at lower left 
    adj_ii(NEIGHBOR_LOWER_LEFT) = ii_left
    adj_jj(NEIGHBOR_LOWER_LEFT) = jj_lower
    adj_ii(NEIGHBOR_UPPER_LEFT) = ii_left
    adj_jj(NEIGHBOR_UPPER_LEFT) = jj_upper
    adj_ii(NEIGHBOR_UPPER_RIGHT) = ii_right
    adj_jj(NEIGHBOR_UPPER_RIGHT) = jj_upper
    adj_ii(NEIGHBOR_LOWER_RIGHT) = ii_right
    adj_jj(NEIGHBOR_LOWER_RIGHT) = jj_lower
  end subroutine get_neighbors

  ! calculate_interp_weights
  ! ------------------------
  ! Uses the module-level observation location to calculate the 
  ! weights for the interpolation scheme.
  ! Currently use bilinear interpolation using the 4 nearest points
  !  PARAMETERS
  !   [none]
  subroutine calculate_interp_weights()

    integer       :: whole_x, whole_y
    real(kind=r8) :: frac_x, frac_y

    real(kind=r8) :: grid_delta_x, grid_delta_y

    whole_x = int(obs_loc_i)
    whole_y = int(obs_loc_j)

    frac_x = obs_loc_i - real(whole_x, kind=r8)
    frac_y = obs_loc_j - real(whole_y, kind=r8)

    ! Bilinear interpolation weights is based on a position within
    ! the unit square: need to account for possible anisotropic spacing
    call get_grid_delta_x(interp_grid, grid_delta_x)
    call get_grid_delta_y(interp_grid, grid_delta_y)
    frac_x = frac_x / grid_delta_x
    frac_y = frac_y / grid_delta_y

    interp_weights(NEIGHBOR_LOWER_LEFT)  = (1 - frac_x) * (1 - frac_y)
    interp_weights(NEIGHBOR_UPPER_LEFT)  = (1 - frac_x) * (frac_y)
    interp_weights(NEIGHBOR_UPPER_RIGHT) = (frac_x)     * (frac_y)
    interp_weights(NEIGHBOR_LOWER_RIGHT) = (frac_x)     * (1 - frac_y)
  end subroutine calculate_interp_weights

  ! initialize_module
  ! -----------------
  ! One-time initialization - this assumes that the "interp_grid"
  ! pointer has already been set to point to the grid used for this
  ! interpolation.  The grid used by the model does not change from 
  ! timestep to timestep (since we do not support nested or
  ! moving grids) so we only need to get the domain parameters once
  ! and don't need to reallocate space every time we want to 
  ! interpolate something.
  !  PARAMETERS
  !   [none]
  subroutine initialize_module()

    call register_module(source, revision, revdate)

    ! Collect information on the grid and the state vector
    call get_grid_dims(interp_grid, grid_i, grid_j)
    call get_grid_field_size(interp_grid, field_size)
    call get_grid_num_levels(interp_grid, num_model_levels)
    call get_num_vars(total_num_vars)

    ! Now that we know how large the grid is we can set up the 
    ! storage arrays.  There is no corresponding deallocate()
    ! function for these items as we assume that the only time
    ! they need to be destroyed is on program termination.
    call allocate_fields()
    call allocate_columns()
    call allocate_coordinates()

    module_initialized = .true.
  end subroutine initialize_module

  ! allocate_fields
  ! ---------------
  ! Set up dynamic memory for the 2-D fields - only need grid
  ! size information to do this
  !  PARAMETERS
  !   [none]
  subroutine allocate_fields()
    
    character(len=15), parameter :: routine = 'allocate_fields'
    integer :: alloc_status

    allocate(var_field(grid_i,grid_j), stat=alloc_status)
    call check_alloc_status(alloc_status, routine, source, revision,&
                            revdate)
  end subroutine allocate_fields

  ! allocate_columns
  ! ----------------
  ! Set up dynamic memory for the vertical columns of values at the
  ! neighboring points
  !  PARAMETERS
  !   [none]
  subroutine allocate_columns()

    character(len=*), parameter :: routine = 'allocate_columns'
    integer :: alloc_status

    allocate(col_values(NUM_NEIGHBORS, num_model_levels),           &
             stat=alloc_status)
    call check_alloc_status(alloc_status, routine, source, revision,&
                            revdate, 'col_values')
    allocate(mean_pres_m_values(NUM_NEIGHBORS, num_model_levels),   &
             stat=alloc_status)
    call check_alloc_status(alloc_status, routine, source, revision,&
                            revdate, 'mean_pres_m_values')
    allocate(pert_pres_m_values(NUM_NEIGHBORS, num_model_levels),   &
             stat=alloc_status)
    call check_alloc_status(alloc_status, routine, source, revision,&
                            revdate, 'pert_pres_m_values')
    allocate(vcoord_values(NUM_NEIGHBORS, num_model_levels),        &
             stat=alloc_status)
    call check_alloc_status(alloc_status, routine, source, revision,&
                            revdate, 'vcoord_values')

    ! After the vertical interpolation to a single level we will have one
    ! value for each neighboring point
    allocate(vinterp_values(NUM_NEIGHBORS, 1), stat=alloc_status)
    call check_alloc_status(alloc_status, routine, source, revision,&
                            revdate, 'vinterp_values')
    allocate(vinterp_levels(1), stat=alloc_status)
    call check_alloc_status(alloc_status, routine, source, revision,&
                            revdate, 'vinterp_levels')
    
    ! These only get used for MSLP computation
    allocate(theta_values(NUM_NEIGHBORS, num_model_levels),         &
             stat=alloc_status)
    call check_alloc_status(alloc_status, routine, source, revision,&
                            revdate, 'theta_values')
    allocate(thetabar_values(NUM_NEIGHBORS, num_model_levels),      &
             stat=alloc_status)
    call check_alloc_status(alloc_status, routine, source, revision,&
                            revdate, 'thetabar_values')
    allocate(mean_pres_w_values(NUM_NEIGHBORS, num_model_levels+1), &
             stat=alloc_status)
    call check_alloc_status(alloc_status, routine, source, revision,&
                            revdate, 'mean_pres_w_values')
  end subroutine allocate_columns

  ! allocate_coordinates
  ! --------------------
  ! Set up dynamic storage for the coordinate information (to help
  ! move around between pressure levels, height levels, sigma levels
  ! etc.)
  !  PARAMETERS
  !   [none]
  subroutine allocate_coordinates()

    character(len=*), parameter :: routine = 'allocate_coordinates'
    integer :: alloc_status

    allocate(terrain_values(NUM_NEIGHBORS),  stat=alloc_status)
    call check_alloc_status(alloc_status, routine, source, revision, &
                            revdate, 'terrain_values')
    allocate(msigma_values(num_model_levels),      stat=alloc_status)
    call check_alloc_status(alloc_status, routine, source, revision, &
                            revdate, 'msigma_values')
    allocate(wsigma_values(num_model_levels + 1),  stat=alloc_status)
    call check_alloc_status(alloc_status, routine, source, revision, &
                            revdate, 'wsigma_values')
    allocate(dsigmaw_values(num_model_levels + 1), stat=alloc_status)
    call check_alloc_status(alloc_status, routine, source, revision, &
                            revdate, 'dsigmaw_values')
  end subroutine allocate_coordinates

  ! allocate_available_values
  ! -------------------------
  ! Set up dynamic storage for the target variable values and the
  ! vertical coordinate values *only* at the levels where they are
  ! available
  !  PARAMETERS
  !   [none]
  subroutine allocate_available_values()

    character(len=*), parameter :: routine = 'allocate_available_values'
    integer :: alloc_status

    allocate(available_target_vals(NUM_NEIGHBORS, num_levels_available), &
             stat=alloc_status)
    call check_alloc_status(alloc_status, routine, source, revision,     &
                            revdate, 'available_target_vals')
    allocate(available_vcoord_vals(NUM_NEIGHBORS, num_levels_available), &
             stat=alloc_status)
    call check_alloc_status(alloc_status, routine, source, revision,     &
                            revdate, 'available_vcoord_vals')
  end subroutine allocate_available_values

  ! deallocate_available_values
  ! ---------------------------
  ! Clean up the dynamic storage for the target variable and vertical
  ! coordinate variable - need to do this after each call to the
  ! interpolation driver since these will potentially change every
  ! time
  !  PARAMETERS
  !   [none]
  subroutine deallocate_available_values()

    character(len=*), parameter :: routine = 'deallocate_available_values'
    integer :: dealloc_status

    deallocate(available_target_vals, stat=dealloc_status)
    call check_dealloc_status(dealloc_status, routine, source,          &
                              revision, revdate, 'avaiable_target_vals')
    deallocate(available_vcoord_vals, stat=dealloc_status)
    call check_dealloc_status(dealloc_status, routine, source,          &
                              revision, revdate, 'available_vcoord_vals')
  end subroutine deallocate_available_values
  
  ! allocate_availability_data
  ! --------------------------
  ! Set up dynamic storage for the main availability array, the
  ! array of which levels contain available data, and the mask that
  ! is used to pull out only the available levels from the values
  ! we read from the restart file
  !  PARAMETERS
  !   [none]
  subroutine allocate_availability_data()

    character(len=*), parameter :: routine = 'allocate_availability_data'
    integer :: alloc_status
    integer :: num_vars_needed

    select case (interpolate_to_level_type)
    case (INTERPOLATE_TO_PRESSURE)
      num_vars_needed = NUM_VARS_NEEDED_FOR_P_INTERP
    case (INTERPOLATE_TO_HEIGHT)
      num_vars_needed = NUM_VARS_NEEDED_FOR_Z_INTERP
    case (INTERPOLATE_TO_SIGMA)
      num_vars_needed = NUM_VARS_NEEDED_FOR_S_INTERP
    case (INTERPOLATE_TO_SURFACE)
      num_vars_needed = NUM_VARS_NEEDED_FOR_SLP
    !case (INTERPOLATE_TO_OTHER)
    !case default
    end select

    allocate(vars_available(num_model_levels, num_vars_needed),     &
             stat=alloc_status)
    call check_alloc_status(alloc_status, routine, source, revision,&
                            revdate, 'vars_available')
    allocate(levels_available(num_model_levels), stat=alloc_status)
    call check_alloc_status(alloc_status, routine, source, revision,&
                            revdate, 'levels_available')
    allocate(level_mask(NUM_NEIGHBORS, num_model_levels),           &
             stat=alloc_status)
    call check_alloc_status(alloc_status, routine, source, revision,&
                            revdate, 'level_mask')
  end subroutine allocate_availability_data

  ! deallocate_availability_data
  ! ----------------------------
  ! Clean up the dynamic storage for the availability information
  ! This needs to be done every time interpolate() is called since it
  ! is potentially different each time
  !  PARAMETERS
  !   [none]
  subroutine deallocate_availability_data()

    character(len=*), parameter :: routine = 'deallocate_availability_data'
    integer :: dealloc_status

    deallocate(vars_available, stat=dealloc_status)
    call check_dealloc_status(dealloc_status, routine, source,     &
                              revision, revdate, 'vars_available')
    deallocate(levels_available, stat=dealloc_status)
    call check_dealloc_status(dealloc_status, routine, source,     &
                              revision, revdate, 'levels_available')
    deallocate(level_mask, stat = dealloc_status)
    call check_dealloc_status(dealloc_status, routine, source,     &
                              revision, revdate, 'level_mask')
  end subroutine deallocate_availability_data

  ! set_interpolation_level_type
  ! ----------------------------
  ! Sets the module-level vertical coordinate type based on the 
  ! observation location
  !  PARAMETERS
  !   IN  obs_loc           location_type of where the ob is
  subroutine set_interpolation_level_type(obs_loc)
    type(location_type), intent(in) :: obs_loc

    ! Note that the DART vertical location type of "model level" is
    ! a sigma level in COAMPS
    if      ( vert_is_pressure(obs_loc) ) then
      interpolate_to_level_type = INTERPOLATE_TO_PRESSURE
    else if ( vert_is_height(obs_loc)   ) then
      interpolate_to_level_type = INTERPOLATE_TO_HEIGHT
    else if ( vert_is_level(obs_loc)    ) then
      interpolate_to_level_type = INTERPOLATE_TO_SIGMA
    else if ( vert_is_surface(obs_loc)  ) then
      interpolate_to_level_type = INTERPOLATE_TO_SURFACE
    else
      interpolate_to_level_type = INTERPOLATE_TO_OTHER
    end if
  end subroutine set_interpolation_level_type

  ! get_vertical_data
  ! -----------------
  ! Pulls out from the grid the model vertical coordinate information 
  ! that we need for interpolation
  !  PARAMETERS
  !   IN  grid              coamps_grid structure that we're using
  subroutine get_vertical_data(grid)
    type(coamps_grid), intent(in) :: grid

    ! Most of these get passed directly into COAMPS interpolation
    ! routines, but we'll use the terrain height, sigma level, and
    ! the total depth of the atmosphere (wsigma(1)) to back out
    ! actual height
    call get_terrain_height_at_points(grid, neighbors_i,         &
                                      neighbors_j, terrain_values) 
    call get_grid_msigma(grid, msigma_values)
    call get_grid_wsigma(grid, wsigma_values)
    call get_grid_dsigmaw(grid, dsigmaw_values)
  end subroutine get_vertical_data

  ! get_target_var
  ! --------------
  ! Generates the data for the target variable by either reading it
  ! from the state vector (for most items) or calculating it
  ! (geopotential height)
  !  PARAMETERS
  !   IN  obs_type          integer representation of raw obs type
  subroutine get_target_var(obs_type)
    integer, intent(in) :: obs_type

    ! If we're searching for something in the state vector, assume
    ! that it is on mass levels and is *not* a mean state variable
    logical, parameter :: ISNT_MEAN  = .false.
    logical, parameter :: IS_M_LEVEL = .true.

    ! Keep track of which sigma level and neighbor point we're at
    integer :: cur_level
    integer :: cur_neighbor

    ! Index of the current field do get from the restart file in the
    ! array of results from find_matching_levels
    !integer :: cur_match
    integer :: cur_field_index

    call get_next_availability_index(AVAILABLE_INDEX_TARGET)

    if (obs_type .eq. QTY_GEOPOTENTIAL_HEIGHT) then
      ! If we're observing geopotential height, there's no need to
      ! go looking around in the state vector - we have all the 
      ! data we need contained in the grid, so just calculate these
      ! right here
      call calculate_heights(terrain_values, col_values)

      ! The height is calculated entirely from the COAMPS
      ! grid variables, so it is guaranteed available at
      ! every sigma level
      vars_available(:,AVAILABLE_INDEX_TARGET) = .true.
    else
      ! If we got here, we're observing something contained in
      ! the state vector directly.
      do cur_level = 1, num_model_levels
        call get_restart_index_by_properties(obs_type, ISNT_MEAN,  &
                                             IS_M_LEVEL, cur_level,&
                                             cur_field_index       )
                                             
        ! Negative index -> not found
        if (cur_field_index < 0) then
          call mark_unavailable(cur_level, AVAILABLE_INDEX_TARGET)
          cycle
        end if

        ! Pull out the entire 2D field - this makes the destaggering
        ! *much* easier.  
        call extract_field_from_state_vector(cur_field_index, &
                                             var_field)

        call destagger_field(var_field, obs_type)

        call extract_neighbors_from_field(var_field,             &
                                          col_values(:,cur_level)) 

        call mark_available(cur_level, AVAILABLE_INDEX_TARGET)
      end do
    end if
  end subroutine get_target_var

  ! get_coordinate_vars
  ! -------------------
  ! Takes the module-level interpolation level type and prepares
  ! the appropriate variables for interpolation
  !  PARAMETERS
  !   [none]
  subroutine get_coordinate_vars()

    integer :: cur_neighbor

    select case (interpolate_to_level_type)
    case (INTERPOLATE_TO_PRESSURE)
      call get_next_availability_index(AVAILABLE_INDEX_MEAN_EXNER)
      call get_next_availability_index(AVAILABLE_INDEX_PERT_EXNER)
      call get_next_availability_index(AVAILABLE_INDEX_PRESSURE)

      ! Do these as function calls for greater flexibility
      ! when we pull the mean state out of the state vector
      ! Note that these functions handle writing the availability 
      ! themselves based on the module-level variables we just got
      ! values for
      call get_mean_exner_values()
      call get_pert_exner_values()
      call convert_exner_to_pressure()
    case (INTERPOLATE_TO_HEIGHT)
      call get_next_availability_index(AVAILABLE_INDEX_HEIGHT)

      ! Data used to calculate heights is read out of the datahd
      ! file and is therefore *guaranteed* available
      call calculate_heights(terrain_values, vcoord_values)
      vars_available(:,AVAILABLE_INDEX_HEIGHT) = .true.
    case (INTERPOLATE_TO_SIGMA)
      call get_next_availability_index(AVAILABLE_INDEX_SIGMA)

      ! Data pulled from the datahd file is always available
      vcoord_values = spread(msigma_values, VALUES_DIM_NEIGHBOR, &
                             NUM_NEIGHBORS)
      vars_available(:,AVAILABLE_INDEX_SIGMA) = .true.
    case (INTERPOLATE_TO_SURFACE)
      ! Skip handling this for now
    case (INTERPOLATE_TO_OTHER)
      ! Skip handling this for now
    case default
      call error_handler(E_ERR, 'get_coordinate_vars', 'Interpolate&
        & level type not found!', source, revision, revdate)
    end select
  end subroutine get_coordinate_vars

  ! calculate_heights
  ! -----------------
  ! Given an 1D array of terrain heights of length N, uses the 
  ! module-level msigma values (of length K) and wsigma values to
  ! return a matrix that is N rows by K columns of the height of
  ! that sigma level given the terrain height at that point.
  ! This can be used to calculate height as either a coordinate
  ! or a target variable
  !  PARAMETERS
  !   IN  terrain_heights   1-D vector of the surface height
  !   OUT height_array      Matrix giving the geometric height
  !                         corresponding to a particular point
  !                         and sigma level
  subroutine calculate_heights(terrain_heights, height_array)
    real(kind=r8), dimension(:), intent(in)    :: terrain_heights
    real(kind=r8), dimension(:,:), intent(out) :: height_array

    ! How many sigma levels are in the model and how many points
    ! we need to calculate the height at
    integer :: num_levels
    integer :: num_points

    integer :: cur_level, cur_point

    ! These are just aliases to help readability
    integer, parameter :: MODEL_TOP = 1
    real(kind=r8) ::  Z, Zs, H, sigma
   

    num_levels = size(msigma_values)
    num_points = size(terrain_heights)

    ! total depth of the model atmosphere
    H = wsigma_values(MODEL_TOP)

    do cur_level = 1, num_levels
      ! No need to check this at every point at this level
      sigma = msigma_values(cur_level)

      do cur_point = 1, num_points
        Zs = terrain_heights(cur_point)
        
        ! sigma = H * (Z-Zs) / (H-Zs)
        Z = sigma * ( (H - Zs) / H ) + Zs
        height_array(cur_point, cur_level) = Z
      end do
    end do
  end subroutine calculate_heights

  ! get_mean_exner_value
  ! --------------------
  ! Get the value of the mean exner function on mass levels. This
  ! currently just reads it out of the state vector but will 
  ! eventually pull it out of the model mean state structure when
  ! I get around to adding that in.
  !  PARAMETERS
  !   [none]
  subroutine get_mean_exner_values()
    
    ! Search parameters for this variable
    logical, parameter :: IS_MEAN    = .true.
    logical, parameter :: IS_M_LEVEL = .true.

    ! Will need to loop over sigma levels
    integer :: cur_level
    
    ! Where the mean exner field is in the long vector
    integer :: mean_exner_index

    ! Don't need to do the field reading here since we don't need
    ! to destagger anything - can just read the points straight out
    ! of the long state vector
    do cur_level = 1, num_model_levels 
      call get_restart_index_by_properties(QTY_EXNER_FUNCTION, &
                                           IS_MEAN, IS_M_LEVEL, &
                                           cur_level,           &
                                           mean_exner_index     )
      
      ! See if there is an entry at this level
      if (mean_exner_index < 0) then
        call mark_unavailable(cur_level, AVAILABLE_INDEX_MEAN_EXNER)
        cycle
      end if

      ! Recall that the values arrays are indexed (neighbor, level)
      call extract_points_from_state_vector(mean_exner_index,   &
                                      neighbors_i, neighbors_j, &
                                      mean_pres_m_values(:,cur_level))
      call mark_available(cur_level, AVAILABLE_INDEX_MEAN_EXNER) 
    end do
  end subroutine get_mean_exner_values

  ! get_pert_exner_value
  ! --------------------
  ! Get the value of the pert exner function on mass levels.
  !  PARAMETERS
  !   [none]
  subroutine get_pert_exner_values()
    
    ! Search parameters for this variable
    logical, parameter :: IS_MEAN    = .false.
    logical, parameter :: IS_M_LEVEL = .true.

    ! Will need to loop over sigma levels
    integer :: cur_level
    
    ! Where the pert exner field is in the long vector
    integer :: pert_exner_index

    ! Don't need to do the field reading here since we don't need
    ! to destagger anything - can just read the points straight out
    ! of the long state vector
    do cur_level = 1, num_model_levels 
      call get_restart_index_by_properties(QTY_EXNER_FUNCTION, &
                                           IS_MEAN, IS_M_LEVEL, &
                                           cur_level,           &
                                           pert_exner_index     )
      
      ! See if there is an entry at this level
      if (pert_exner_index < 0) then
        call mark_unavailable(cur_level, AVAILABLE_INDEX_PERT_EXNER)
        cycle
      end if

      ! Recall that the values arrays are indexed (neighbor, level)
      call extract_points_from_state_vector(pert_exner_index,   &
                                      neighbors_i, neighbors_j, &
                                      pert_pres_m_values(:,cur_level))
      call mark_available(cur_level, AVAILABLE_INDEX_PERT_EXNER) 
    end do
  end subroutine get_pert_exner_values

  ! convert_exner_to_pressure
  ! -------------------------
  ! Uses the module-level exner values to convert to pressure
  !  exner = (p/p00)^(R/Cp)
  !  PARAMETERS
  !   [none]
  subroutine convert_exner_to_pressure()
  
    ! No need for loops here - can just do this in one fell swoop
    vcoord_values = (mean_pres_m_values + pert_pres_m_values)
    vcoord_values = vcoord_values ** (Cp / R)
    vcoord_values = vcoord_values * P00

    ! That will calculate even if we have junk data - make sure
    ! that the availability matrix reflects that we have both the
    ! mean and perturbation exner functions and hence, the pressure
    vars_available(:, AVAILABLE_INDEX_PRESSURE) =           &
              vars_available(:, AVAILABLE_INDEX_MEAN_EXNER) &
              .and.                                         &
              vars_available(:, AVAILABLE_INDEX_PERT_EXNER)
  end subroutine convert_exner_to_pressure

  ! extract_field_from_state_vector
  ! -------------------------------
  ! Given the index of a field within the large state vector, turns
  ! it into an offset and reads in the entire field from the state
  ! vector into a supplied array.
  !  PARAMETERS
  !   IN  field_index       position of this field in the long
  !                         state vector
  !   OUT field_data        array populated from the state vector
  subroutine extract_field_from_state_vector(field_index, field_data)
    integer, intent(in)                        :: field_index
    real(kind=r8), dimension(:,:), intent(out) :: field_data

    ! Where the field starts in terms of r8 blocks
    integer :: field_offset

    ! Indices for the where we are in the 2D field and 1D state vec 
    integer :: ii, jj, nn
    
    ! Assert that the index that is supplied is actually within the
    ! bounds of the state vector
    call check_field_index(field_index)

    call restart_index_to_vector_offset(field_index, field_offset)
    do jj = 1, grid_j
      do ii = 1, grid_i
        nn = field_offset + ( (jj - 1) * grid_i + ii )
        field_data(ii, jj) = state_vector(nn)
      end do
    end do

  end subroutine extract_field_from_state_vector

  ! extract_points_from_state_vector
  ! --------------------------------
  ! The same idea as extract_field_from_state_vector, but instead of
  ! taking the entire 2D field, just pulls out selected points given
  ! an index of a field within the large state vector
  subroutine extract_points_from_state_vector(field_index, i_coords,&
                                              j_coords, field_data  )
    integer,                     intent(in)  :: field_index
    integer,       dimension(:), intent(in)  :: i_coords
    integer,       dimension(:), intent(in)  :: j_coords
    real(kind=r8), dimension(:), intent(out) :: field_data

    ! where the field starts in terms of r8 blocks
    integer :: field_offset

    ! How many coordinates we are working with - these should be
    ! the same
    integer :: num_i_coords, num_j_coords

    ! Keep track of the current coordinate we're processing,
    ! its corresponding i/j indices, as well as where that is
    ! in the long state vector
    integer :: cur_coord_index
    integer :: ii, jj, nn

    ! Make sure that we have enough coordinate data and a place to 
    ! put the results
    num_i_coords = size(i_coords)
    num_j_coords = size(j_coords)
    if ( num_i_coords .ne. num_j_coords) then
      call error_handler(E_ERR, 'extract_points_from_state_vector',&
                         'Number of i/j coordinates do not match!',&
                         source, revision, revdate                 )
    end if
    if (size(field_data) < num_i_coords) then
      call error_handler(E_ERR, 'extract_points_from_state_vector',&
                         'Not enough storage space provided!',     &
                         source, revision, revdate                 )
    end if

    ! Assert that the index supplied is actually in the state vector
    call check_field_index(field_index)

    call restart_index_to_vector_offset(field_index, field_offset)

    ! note that num_i_coords and num_j_coords are the same if we
    ! have gotten this far, so only need to use one of the two
    do cur_coord_index = 1, num_i_coords
      ii = i_coords(cur_coord_index)
      jj = j_coords(cur_coord_index)
      nn = field_offset + ( (jj - 1) * grid_i + ii )
      
      field_data(cur_coord_index) = state_vector(nn)
    end do
  end subroutine extract_points_from_state_vector


  ! check_field_index
  ! -----------------
  ! Checks a supplied field index against the known size of the
  ! state vector to make sure it's within the bounds of that vector
  ! Since this code is called internally from here, it should really
  ! never fail.  If it does, though, complain loudly.
  !  PARAMETERS
  !   IN  field_index       index of the field within the state vec
  subroutine check_field_index(field_index)
    integer, intent(in) :: field_index

    if ((field_index < 1) .or. (field_index > total_num_vars)) then
      call error_handler(E_ERR, 'check_field_index',               &
                         'Field index out of range', source,       &
                         revision, revdate)
    end if
  end subroutine check_field_index


  ! restart_index_to_vector_offset
  ! ------------------------------
  ! Calculates where a particular entry in the restart file begins
  ! given the size of the grid.
  !  PARAMETERS
  !   IN  var_index         where the variable is in the restart
  !                         array
  !   OUT var_offset        beginning index for this field in the
  !                         state vector
  subroutine restart_index_to_vector_offset(var_index, var_offset)
    integer, intent(in)  :: var_index
    integer, intent(out) :: var_offset

    var_offset = (var_index - 1) * field_size 
  end subroutine restart_index_to_vector_offset

  ! destagger_field
  ! ---------------
  ! Given a 2D field and an observation type, destaggers the field
  ! if necessary using the COAMPS destaggering functions.  This
  ! subroutine assumes that the field dimensions is given by the
  ! module levels grid_i and grid_j (i.e. that we're destaggering
  ! fields from *this* grid)
  !  PARAMETERS
  ! INOUT field             the 2-D field to destagger
  !   IN  field_type        index representation of the raw type
  !                         contained in the field
  subroutine destagger_field(field, field_type)
    real(kind=r8), dimension(:,:), intent(inout) :: field
    integer,                       intent(in)    :: field_type

    ! Since we're processing fields, there is only 1 vertical level
    integer, parameter :: NUM_VERT_LEVELS = 1

    ! Only need to worry about u/v winds - everything else is
    ! returned as is
    if (field_type .eq. QTY_U_WIND_COMPONENT) then
      call utom(field, grid_i, grid_j, NUM_VERT_LEVELS)
    else if (field_type .eq. QTY_V_WIND_COMPONENT) then
      call vtom(field, grid_i, grid_j, NUM_VERT_LEVELS)
    end if
  end subroutine destagger_field

  ! extract_neighbors_from_field
  ! ----------------------------
  ! Given a 2D field and which level it represents, extracts the
  ! neighboring points that will be used for the interpolation
  ! Pass array slices into this array to avoid needing to carry 
  ! around what the level is
  !  PARAMETERS
  !   IN  field             the entire 2d field to pull data from
  ! INOUT neighbors         the array containing the value of field
  !                         at the neighboring points
  subroutine extract_neighbors_from_field(field, neighbors)
    real(kind=r8), dimension(:,:), intent(in)    :: field
    real(kind=r8), dimension(:),   intent(inout) :: neighbors

    ! which neighbor to process and what it's grid coordinates are
    integer :: cur_neighbor_num
    integer :: cur_neighbor_ii
    integer :: cur_neighbor_jj

    do cur_neighbor_num = 1, NUM_NEIGHBORS
      cur_neighbor_ii = neighbors_i(cur_neighbor_num)
      cur_neighbor_jj = neighbors_j(cur_neighbor_num)
      neighbors(cur_neighbor_num) = field( cur_neighbor_ii, &
                                           cur_neighbor_jj  )
    end do                                           
  end subroutine extract_neighbors_from_field

  ! mark_available
  ! --------------
  ! Given a level and variable index, marks the variable available
  ! at that particular level
  !  PARAMETERS
  !   IN  level             sigma level of the variable
  !   IN  var_index         variable's index for availability matrix
  subroutine mark_available(level, var_index)
    integer, intent(in) :: level
    integer, intent(in) :: var_index

    ! Each variable has its own column in the availability matrix
    vars_available(level, var_index) = .true.
  end subroutine mark_available

  ! mark_unavailable
  ! ----------------
  ! Given a level and variable index, marks the variable unavailable
  ! at that particular level
  !  PARAMETERS
  !   IN  level             sigma level of the variable
  !   IN  var_index         variable's index for availability matrix
  subroutine mark_unavailable(level, var_index)
    integer, intent(in) :: level
    integer, intent(in) :: var_index

    ! Each variable has its own column in the availability matrix
    vars_available(level, var_index) = .false.
  end subroutine mark_unavailable

  ! calculate_available_levels
  ! --------------------------
  ! Once all the availability statistics have been populated, this
  ! routine will calculate the levels that have all the data present
  !  PARAMETERS
  !   [none]
  subroutine calculate_available_levels()

    integer :: cur_level

    ! A level is defined as "available" if every variable needed
    ! for the interpolation is present at that level
    do cur_level = 1, size(vars_available,dim=AVAILABILITY_DIM_LEVEL)
      levels_available(cur_level) = all(vars_available(cur_level,:))
    end do

    num_levels_available = count(levels_available)
  end subroutine calculate_available_levels

  ! generate_level_mask
  ! -------------------
  ! Once the available levels have been calculated, generate a mask
  ! that will be used to pull out only the available levels from the
  ! data arrays.
  !  PARAMETERS
  !   [none]
  subroutine generate_level_mask()

    level_mask = spread(levels_available, VALUES_DIM_NEIGHBOR, &
                        NUM_NEIGHBORS)
  end subroutine generate_level_mask

  ! interpolate_to_level
  ! --------------------
  ! Interpolates an array to a single level using the interpolation
  ! functions from the COAMPS utility package.
  subroutine interpolate_to_level()
    
    ! Array size information for the interpolators
    integer :: num_levels_in
    integer :: num_levels_out
 
    num_levels_in  = num_levels_available
    num_levels_out = 1

    select case (interpolate_to_level_type)
    case (INTERPOLATE_TO_PRESSURE)
      ! DART location module stores pressure in Pa and the COAMPS
      ! pressure calculated from the Exner function is in mb
      vinterp_levels(1) = obs_loc_z * CONVERT_PA_TO_MB

      call s2pint(available_target_vals, vinterp_values,        &
                  available_vcoord_vals, vinterp_levels,        &
                  num_levels_in, num_levels_out, NUM_NEIGHBORS, &
                  USE_MISSING_VALUE, MISSING_VALUE)
    case (INTERPOLATE_TO_HEIGHT, INTERPOLATE_TO_SIGMA)
      ! DART location module stores height in meters so no conversion
      ! is necessary
      vinterp_levels(1) = obs_loc_z

      ! Same function call will work on either sigma interpolation
      ! or height interpolation - only difference is the values in
      ! the vertical coordinate variable
      call z2zint(available_target_vals, vinterp_values,         &
                  available_vcoord_vals, vinterp_levels,         &
                  num_levels_in, num_levels_out, NUM_NEIGHBORS,  &
                  USE_MISSING_VALUE, MISSING_VALUE)
    end select
  end subroutine interpolate_to_level

  ! interpolate_to_point
  ! --------------------
  ! Interpolates an array to a single point.  Uses module-level
  ! variables for arrays holding the data and the weights
  !  PARAMETERS
  !   OUT obs_value         result of interpolation
  subroutine interpolate_to_point(obs_value)
    real(kind=r8), intent(out) :: obs_value

    ! dot_product expects a 1-D vector - vinterp_values is a 2-D
    ! matrix with a singleton dimension.  Could also do this with 
    ! reshape(vinterp_values, (/ NUM_NEIGHBORS /))
    ! instead of pack()
    obs_value = dot_product(interp_weights,              &
                            pack(vinterp_values, .true.) )
  end subroutine interpolate_to_point

  ! get_next_availability_index
  ! ---------------------------
  ! Return where to write a variable availability within the 
  ! availability array and change the master index to reflect
  ! the request.
  !  PARAMETERS
  !   OUT next_index        The next available index in the
  !                         availability array
  subroutine get_next_availability_index(next_index)
    integer, intent(out) :: next_index

    availability_index = availability_index + 1
    next_index         = availability_index
  end subroutine get_next_availability_index

  ! get_available_values
  ! --------------------
  ! Uses the module-level mask for which levels are available to
  ! pull data from the supplied array,
  !  PARAMETERS
  !   [none]
  subroutine get_available_values(all_values, available_values)
    real(kind=r8), dimension(:,:), intent(in)  :: all_values
    real(kind=r8), dimension(:,:), intent(out) :: available_values

    integer, dimension(2) :: new_shape

    ! The new array will still be defined as (levels, neighbors) but
    ! now the extent of levels is defined as how many are available
    new_shape = (/ NUM_NEIGHBORS, num_levels_available /)

    ! This works since the first dimension of the array does not 
    ! change size.  When we re-order the vector created by pack, 
    ! we get all the neighbor values at only the levels defined by 
    ! .true. values in the level mask
    available_values = reshape(pack(all_values,level_mask),new_shape)
  end subroutine get_available_values

  ! check_level_type
  ! ----------------
  ! Ensures that the level type defined by the observation location
  ! is one that we actually will interpolate to
  !  PARAMETERS
  !   OUT valid_level_type  true if we can use this level type
  subroutine check_level_type(valid_level_type)
    logical, intent(out) :: valid_level_type

    character(len=128) :: message

    if (interpolate_to_level_type .eq. INTERPOLATE_TO_PRESSURE .or.&
        interpolate_to_level_type .eq. INTERPOLATE_TO_HEIGHT   .or.&
        interpolate_to_level_type .eq. INTERPOLATE_TO_SIGMA ) then
      valid_level_type = .true.
    else
      write (message,*) 'Level type ', interpolate_to_level_type, &
                        'is not supported.'
      call error_handler(E_MSG, 'check_level_type', trim(message),&
                         source, revision, revdate)
      valid_level_type = .false.
    end if
  end subroutine check_level_type

  ! check_enough_levels_available
  ! -----------------------------
  ! Ensures that we have enough vertical levels to do a meaningful
  ! interpolation.  Since we're currently not supporting surface
  ! level types, use the module-level parameter for the threshold
  !  PARAMETERS
  !   OUT enough_levels_available   true if there are enough levels
  !                                 to do the interpolation
  subroutine check_enough_levels_available(enough_levels_available)
    logical, intent(out) :: enough_levels_available

    character(len=128) :: message

    if (num_levels_available < MIN_LEVELS_NEEDED) then
      write (message,*) 'There are ', num_levels_available,        &
                        'levels available, which is less than the',&
                        MIN_LEVELS_NEEDED, 'needed to interpolate.'
      call error_handler(E_MSG, 'check_enough_levels_available', &
                         trim(message), source, revision, revdate)
      enough_levels_available = .false.
    else
      enough_levels_available = .true.
    end if
  end subroutine check_enough_levels_available

  ! check_level_in_available_range
  ! ------------------------------
  ! Ensures that the vertical level we are trying to interpolate
  ! to lies within the available levels that we have for the
  ! appropriate vertical coordinate type.
  !  PARAMETERS
  !   OUT level_in_range    true if the interpolating level is
  !                         within the range of available levels
  subroutine check_level_in_available_range(level_in_range)
    logical, intent(out) :: level_in_range

    character(len=256) :: message

    real(kind=r8) :: scaled_level

    ! If we're interpolating onto pressure levels, DART stores
    ! the pressure in Pa and we calculate in hPa
    if (interpolate_to_level_type == INTERPOLATE_TO_PRESSURE) then
      scaled_level = obs_loc_z * CONVERT_PA_TO_MB
    else
      scaled_level = obs_loc_z
    end if

    ! We're OK if the level is below the highest available and
    ! higher than the lowest available
    if (scaled_level .le. maxval(available_vcoord_vals) .and.      &
        scaled_level .ge. minval(available_vcoord_vals)     ) then
      level_in_range = .true.
    else
      write (message,*) 'Trying to interpolate to vertical level', &
                        scaled_level, 'but the available range is',&
                        minval(available_vcoord_vals),'to',        &
                        maxval(available_vcoord_vals)
      call error_handler(E_MSG, 'check_level_in_available_range',  &
                         trim(message), source, revision, revdate)
      level_in_range = .false.
    end if
  end subroutine check_level_in_available_range

  ! check_for_missing_values
  ! ------------------------
  ! Looks at the results of the level interpolation and ensures 
  ! that there are no values labeled as "missing"
  !  PARAMETER
  !   OUT no_missing_values     true if all the values are present
  subroutine check_for_missing_values(no_missing_values)
    logical, intent(out) :: no_missing_values
    
    character(len=128) :: message

    integer :: num_missing_values

    num_missing_values = count(int(vinterp_values) == MISSING_VALUE)
    if (num_missing_values .eq. 0) then
      no_missing_values = .true.
    else
      write (message,*) 'Missing values were found after level ',&
                        'interpolation'
      call error_handler(E_MSG, 'check_for_missing_values',      &
                         trim(message), source, revision, revdate)
    end if
  end subroutine check_for_missing_values

  ! print_interpolation_diagnostics
  ! -------------------------------
  ! Print some data about the interpolation
  !  PARAMETERS
  !   IN  obs_type          Type of the observation (integer)
  !   IN  obs_value         Result of the interpolation
  subroutine print_interpolation_diagnostics(obs_type, obs_value)
    integer,       intent(in) :: obs_type
    real(kind=r8), intent(in) :: obs_value
    
    if (do_output()) then
      call print_label_name('Interpolation Results')
      write (*,'(A15T25F15.6)') 'I Coordinate  :', obs_loc_i
      write (*,'(A15T25F15.6)') 'J Coordinate  :', obs_loc_j
      write (*,'(A15T25F10.2)') 'Vertical      :', obs_loc_z
      write (*,'(A15T25I10)'  ) 'Vertical Type :', interpolate_to_level_type
      write (*,'(A15T25I10)'  ) 'Variable Type :', obs_type
      write (*,'(A15T25F15.6)') 'Value         :', obs_value
      write (*,*)
    end if
  end subroutine print_interpolation_diagnostics

  !------------------------------
  ! END PRIVATE ROUTINES
  !------------------------------

end module coamps_interp_mod

! <next few lines under version control, do not edit>
! $URL$
! $Id$
! $Revision$
! $Date$
