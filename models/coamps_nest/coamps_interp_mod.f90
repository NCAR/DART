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
  
    use coamps_domain_mod,    only : coamps_domain, location_to_nest_point, &
                                     get_domain_msigma, get_domain_wsigma,  &
                                     get_domain_dsigmaw,                    &
                                     get_domain_num_levels,                 &
                                     nest_point_to_latlon

    use coamps_nest_mod,      only : coamps_nest, get_terrain_height_at_points, &
                                     get_i_coord, get_j_coord, get_nest,        &
                                     nest_point, nest_index_2d_to_1d,           &
                                     get_nest_i_width, get_nest_j_width,        &
                                     get_nest_size, get_nest_level,             &
                                     in_this_nest

    use coamps_statevar_mod,  only : state_variable, get_var_substate, is_sigma_level

    use coamps_statevec_mod,  only : state_vector, find_state_variable

    use coamps_intrinsic_mod, only : s2pint,                         &
                                     sfcp,                           &
                                     utom,                           &
                                     vtom,                           &
                                     vor,                            &
                                     z2zint

    use coamps_util_mod,      only : check_alloc_status,             &
                                     check_dealloc_status,           &
                                     trace_message,           &
                                     print_label_name
  
    use location_mod,         only : location_type,    &
                                     query_location,   &
                                     vert_is_height,   &
                                     vert_is_undef,    &
                                     vert_is_level,    &
                                     vert_is_pressure, &
                                     vert_is_surface
    use obs_kind_mod
    use types_mod,            only : MISSING_I,        &
                                     MISSING_R8,       &
                                     r8
    use utilities_mod,        only : do_output,        &
                                     E_ERR,            &
                                     E_MSG,            &
                                     E_WARN,           &
                                     error_handler,    &
                                     register_module 
  
    implicit none
  
    private
  
    !------------------------------
    ! BEGIN PUBLIC INTERFACE
    !------------------------------
  
    public :: interpolate
    public :: set_interp_diag
  
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

    ! Set below-ground levels as "missing" in the vertical interpolation
    integer, parameter :: USE_MISSING_VALUE = 1
  
    ! Number of neighboring points to use for horizontal interpolation
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
    integer, parameter :: INTERPOLATE_TO_UNDEF    = 5
    integer, parameter :: INTERPOLATE_TO_OTHER    = 6

    ! Pressure conversion constant
    real(kind=r8), parameter :: CONVERT_PA_TO_MB = 0.01_r8
    real(kind=r8), parameter :: CONVERT_MB_TO_PA = 100.0_r8

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
  
    integer, parameter :: SINGLE_LEVEL = 1
    integer, parameter :: SINGLE_POINT = 1
    real(kind=r8), parameter :: ONE_R8 = 1.0_r8

    ! Interpolation data type
    type :: coamps_interpolator
        private

        ! Shortcut references
        real(kind=r8), dimension(:), pointer :: model_state 
        type(coamps_domain),         pointer :: model_domain
        type(coamps_nest),           pointer :: interp_nest
        type(state_vector),          pointer :: state_definition

        ! Location (horizontal and vertical) of where we will interpolate
        type(nest_point)  :: interp_point
        real(kind=r8)     :: interp_level
        integer           :: interp_level_type

        ! Neighboring horizontal points
        integer, dimension(NUM_NEIGHBORS) :: neighbors_i
        integer, dimension(NUM_NEIGHBORS) :: neighbors_j

        ! Raw values - indexed by (neighbor, sigma level)
        real(kind=r8), dimension(:,:), pointer :: target_values
        real(kind=r8), dimension(:,:), pointer :: vcoord_values

        ! Variable availability is stored as an (# of levels) x (# of vars)
        ! array - if the value at (level, variable) is .true., then that
        ! variable is available at that level.  Level availability is the 
        ! 1-D analog, .true. when AND(vars_available(level, :)) is true
        logical, dimension(:,:), pointer :: vars_available
        logical, dimension(:),   pointer :: levels_available
        integer                          :: num_levels_available

        ! Filtered raw values only on levels that have all data available
        real(kind=r8), dimension(:,:), pointer :: available_target_values
        real(kind=r8), dimension(:,:), pointer :: available_vcoord_values

        ! Results of vertical interpolation
        real(kind=r8), dimension(:,:),          pointer :: vinterp_values
        real(kind=r8), dimension(SINGLE_LEVEL)          :: vinterp_level

        ! Weights for the horizontal interpolation
        real(kind=r8), dimension(NUM_NEIGHBORS) :: interp_weights

        ! Surface elevation values at neighboring points
        real(kind=r8), dimension(NUM_NEIGHBORS) :: zsfc_values
    end type

    ! Flag to output interpolation diagnostics
    logical :: output_interpolation

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
  
    integer, save :: NUM_MODEL_LEVELS
  
    ! Since the number of variables required for each type of interpolation
    ! is different, just increment an index as needed, storing the value in
    ! differently named fields
    integer :: cur_availability_index
  
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
    !   IN  domain            COAMPS domain to interpolate on
    !   IN  state_def         COAMPS state vector definition
    !   IN  obs_loc           DART location structure to interpolate to
    !   IN  obs_kind          integer version of raw variable type
    !   OUT obs_value         result of interpolation
    !   OUT interp_worked     true if interpolation was successful
    subroutine interpolate(state, domain, state_def, obs_loc, obs_kind, &
                           obs_value, interp_worked)
        real(kind=r8), dimension(:), intent(in)  :: state
        type(coamps_domain),         intent(in)  :: domain
        type(state_vector),          intent(in)  :: state_def
        type(location_type),         intent(in)  :: obs_loc
        integer,                     intent(in)  :: obs_kind
        real(kind=r8),               intent(out) :: obs_value
        logical, optional,           intent(out) :: interp_worked

        type(coamps_interpolator) :: interpolator
        integer :: ii,k
        logical :: is_success

        if (.not. module_initialized) call initialize_module(domain)

        call initialize_interpolator(interpolator, state, domain, state_def)

        call set_interpolation_location(interpolator, obs_loc)

        if (.not. has_valid_location(interpolator)) then
            interp_worked = .false.
            return
        end if

        ! Figure out which points out of an entire field need to be read 
        call get_nearest_neighbors(interpolator) 

        ! The target variable is defined by the observation type and the
        ! coordinate variables are defined by the observation location's 
        ! vertical level type
        call reset_availability_index()
        call initialize_availability_data(interpolator)

        ! Weights for bilinear horizontal interpolation
        call calculate_interp_weights(interpolator)

        select case (obs_kind)
        case (QTY_SURFACE_PRESSURE)
          call calculate_surface_pressure(interpolator)

        case (QTY_SURFACE_ELEVATION)
          call calculate_surface_heights(interpolator)

        case default
           
          ! Try to find if the state variable is defined on the same level 
          ! as the observation.  If it is not, interpolate in the vertical. 
          if( vert_is_height(obs_loc) ) then
             call calculate_height_level_var(interpolator, obs_kind, &
                     query_location(obs_loc, 'VLOC'), is_success)

          !elseif( vert_is_pressure(obs_loc) ) then
          !elseif( vert_is_surface(obs_loc) ) then
          elseif( vert_is_undef(obs_loc) ) then
             call calculate_undef_level_var(interpolator, obs_kind, &
                     query_location(obs_loc, 'VLOC'), is_success)
          end if 

          if(.not. is_success) then

          call get_target_var(interpolator, obs_kind)
          call get_coordinate_vars(interpolator)

          ! Interpolation is spotty if there aren't enough vertical levels,
          ! so declare failure rather than returning a (probably bad) result
          call calculate_available_levels(interpolator)
          if (.not. enough_levels_available(interpolator)) then
            interp_worked = .false.
            return
          end if

          call initialize_available_values(interpolator)
          call collect_available_values(interpolator)

          ! Avoid extrapolation and only interpolate values at vertical levels
          ! bounded by the available levels
          if (.not. interp_level_in_available_range(interpolator)) then
            interp_worked = .false.
            return
          end if

          ! Interpolate everything to a single level...
          call interpolate_to_level(interpolator)
          if (.not. no_missing_values(interpolator)) then
            interp_worked = .false.
            return
          end if
          end if

        end select

        ! ... then interpolate neighbors to a single point
        obs_value = interpolate_to_point(interpolator) 

        interp_worked = .true.

        if(output_interpolation) &
        call print_interpolation_diagnostics(interpolator, obs_kind, obs_value)

    end subroutine interpolate

  ! set_interp_diag
  ! ---------------
  ! Sets flag to output interpolation diagnostics
  ! IN  flag to output interpolation diagnostics
  subroutine set_interp_diag(interp_diag)
    logical, intent(in) :: interp_diag
    output_interpolation=interp_diag
  end subroutine set_interp_diag

    !------------------------------
    ! END PUBLIC ROUTINES
    !------------------------------

    !------------------------------
    ! BEGIN PRIVATE ROUTINES
    !------------------------------

    ! set_interpolation_location
    ! --------------------------
    ! Populate the location component of an interpolation structure given
    ! a DART observation.  
    subroutine set_interpolation_location(interpolator, obs_location)
        type(coamps_interpolator),  intent(inout) :: interpolator
        type(location_type),        intent(in)    :: obs_location

        ! Take the finest nest available for this point
        call location_to_nest_point(interpolator%model_domain, &
                                    obs_location,              &
                                    interpolator%interp_point, &
                                    interpolator%interp_level)
        call set_interpolation_level_type(interpolator, obs_location)

        ! Save an alias - this will eliminate a lot of calls to "get_nest"
        interpolator%interp_nest => get_nest(interpolator%interp_point)
    end subroutine set_interpolation_location

    ! has_valid_location
    ! ------------------
    ! Check the correctness of a proposed interpolation location 
    function has_valid_location(interpolator)
        type(coamps_interpolator), intent(in)  :: interpolator
        logical                                :: has_valid_location

        ! Need to check the horizontal and vertical components separately
        if (is_valid_level_type(interpolator%interp_level_type) .and. &
            in_this_nest(interpolator%interp_point)) then
            has_valid_location = .true.
        else
            has_valid_location = .false.
        end if

    end function has_valid_location

    ! get_nearest_neighbors
    ! ---------------------
    ! Calculates the four nearest neighbors of a given point.  This
    ! assumes that the i/j coordinates supplied is within the bounds
    ! of the grid, so that the neighboring points will be as well.
    subroutine get_nearest_neighbors(interpolator)
        type(coamps_interpolator), intent(inout) :: interpolator

        integer :: jj_lower, jj_upper
        integer :: ii_left, ii_right

        ii_left  = floor(get_i_coord(interpolator%interp_point))
        ii_right = ceiling(get_i_coord(interpolator%interp_point))
        jj_lower = floor(get_j_coord(interpolator%interp_point))
        jj_upper = ceiling(get_j_coord(interpolator%interp_point))

        ! Ordering goes clockwise starting at lower left 
        interpolator%neighbors_i(NEIGHBOR_LOWER_LEFT)  = ii_left
        interpolator%neighbors_j(NEIGHBOR_LOWER_LEFT)  = jj_lower

        interpolator%neighbors_i(NEIGHBOR_UPPER_LEFT)  = ii_left
        interpolator%neighbors_j(NEIGHBOR_UPPER_LEFT)  = jj_upper

        interpolator%neighbors_i(NEIGHBOR_UPPER_RIGHT) = ii_right
        interpolator%neighbors_j(NEIGHBOR_UPPER_RIGHT) = jj_upper

        interpolator%neighbors_i(NEIGHBOR_LOWER_RIGHT) = ii_right
        interpolator%neighbors_j(NEIGHBOR_LOWER_RIGHT) = jj_lower
    end subroutine get_nearest_neighbors

    ! calculate_interp_weights
    ! ------------------------
    ! Calculate how each neighbor should be weighted for horizontal
    ! interpolation.
    subroutine calculate_interp_weights(interpolator)
        type(coamps_interpolator), intent(inout) :: interpolator

        real(kind=r8) :: obs_x, obs_y
        real(kind=r8) :: frac_x, frac_y

        obs_x = get_i_coord(interpolator%interp_point)
        obs_y = get_j_coord(interpolator%interp_point)

        ! Bilinear interpolation weights are based on a position within
        ! the unit square.  Map the square [x,y]->[x+dx,y+dy] to 
        ! [0,0]->[1,1] by shifting and scaling to account for any anisotropy
        frac_x = obs_x - real(int(obs_x), kind=r8)
        frac_y = obs_y - real(int(obs_y), kind=r8)

        interpolator%interp_weights(NEIGHBOR_LOWER_LEFT)  = (1 - frac_x) * &
                                                            (1 - frac_y)

        interpolator%interp_weights(NEIGHBOR_UPPER_LEFT)  = (1 - frac_x) * &
                                                            (    frac_y)

        interpolator%interp_weights(NEIGHBOR_UPPER_RIGHT) = (    frac_x) * &
                                                            (    frac_y)

        interpolator%interp_weights(NEIGHBOR_LOWER_RIGHT) = (    frac_x) * &
                                                            (1 - frac_y)
    end subroutine calculate_interp_weights

    ! initialize_module
    ! -----------------
    ! Handle one-time initialization.   The number and height of the model 
    ! levels does not vary between nests, so we can get that information
    ! once. 
    subroutine initialize_module(domain)
        type(coamps_domain), intent(in) :: domain

        output_interpolation = .false.

        call register_module(source, revision, revdate)

        NUM_MODEL_LEVELS = get_domain_num_levels(domain)

        module_initialized = .true.
    end subroutine initialize_module

    ! initialize_interpolator
    ! -----------------------
    ! Set up the interpolator object
    subroutine initialize_interpolator(interpolator, state, domain, state_def)
        type(coamps_interpolator),           intent(inout) :: interpolator
        real(kind=r8), dimension(:), target, intent(in)    :: state
        type(coamps_domain),         target, intent(in)    :: domain
        type(state_vector),          target, intent(in)    :: state_def

        interpolator%model_state      => state
        interpolator%model_domain     => domain
        interpolator%state_definition => state_def

        call allocate_raw_values(interpolator)
    end subroutine initialize_interpolator


    ! finalize_interpolator
    ! ---------------------
    ! Clean up the interpolator object
    subroutine finalize_interpolator(interpolator)
        type(coamps_interpolator), intent(inout) :: interpolator

        nullify(interpolator%model_state)
        nullify(interpolator%model_domain)
        nullify(interpolator%state_definition)

        call deallocate_raw_values(interpolator)
        call deallocate_available_values(interpolator)
        call deallocate_availability_data(interpolator)
    end subroutine finalize_interpolator

    ! allocate_raw_values
    ! -------------------
    ! Set up dynamic memory for the vertical raw_values of values at the
    ! neighboring points
    subroutine allocate_raw_values(interpolator)
        type(coamps_interpolator), intent(inout) :: interpolator

        character(len=*), parameter :: routine = 'allocate_raw_values'
        integer                     :: alloc_status

        allocate(interpolator%target_values(NUM_NEIGHBORS, num_model_levels),&
                 stat=alloc_status)
        call check_alloc_status(alloc_status, routine, source, revision,     &
                                revdate, 'target_values')
        allocate(interpolator%vcoord_values(NUM_NEIGHBORS, num_model_levels),&
                 stat=alloc_status)
        call check_alloc_status(alloc_status, routine, source, revision,     &
                                revdate, 'vcoord_values')

        ! After the vertical interpolation to a single level we will have one
        ! value for each neighboring point
        allocate(interpolator%vinterp_values(NUM_NEIGHBORS, SINGLE_LEVEL),   &
                 stat=alloc_status)
        call check_alloc_status(alloc_status, routine, source, revision,     &
                                revdate, 'vinterp_values')

    end subroutine allocate_raw_values

    ! deallocate_raw_values
    ! ------------------
    ! Clean up dynamic memory for the vertical raw_values of values at the
    ! neighboring points
    subroutine deallocate_raw_values(interpolator)
        type(coamps_interpolator), intent(inout) :: interpolator

        character(len=*), parameter :: routine = 'deallocate_raw_values'
        integer                     :: dealloc_status

        deallocate(interpolator%target_values, stat=dealloc_status)
        call check_dealloc_status(dealloc_status, routine, source, revision, &
                                revdate, 'target_values')
        deallocate(interpolator%vcoord_values, stat=dealloc_status)
        call check_dealloc_status(dealloc_status, routine, source, revision, &
                                revdate, 'vcoord_values')

        ! After the vertical interpolation to a single level we will have one
        ! value for each neighboring point
        deallocate(interpolator%vinterp_values, stat=dealloc_status)
        call check_dealloc_status(dealloc_status, routine, source, revision, &
                                revdate, 'vinterp_values')

    end subroutine deallocate_raw_values

    ! deallocate_available_values
    ! ---------------------------
    ! Clean up the dynamic storage for the target variable and vertical
    ! coordinate variable - need to do this after each call to the
    ! interpolation driver since these will potentially change every
    ! time
    subroutine deallocate_available_values(interpolator)
        type(coamps_interpolator), intent(inout) :: interpolator

        character(len=*), parameter :: routine = 'deallocate_available_values'
        integer :: dealloc_status

        deallocate(interpolator%available_target_values, stat=dealloc_status)
        call check_dealloc_status(dealloc_status, routine, source, revision, &
                                  revdate, 'available_target_values')
        deallocate(interpolator%available_vcoord_values, stat=dealloc_status)
        call check_dealloc_status(dealloc_status, routine, source, revision, &
                                  revdate, 'available_vcoord_values')
    end subroutine deallocate_available_values
  
    ! initialize_availability_data
    ! --------------------------
    ! Set up dynamic storage for the main availability array, the
    ! array of which levels contain available data, and the mask that
    ! is used to pull out only the available levels from the values
    ! we read from the restart file
    !  PARAMETERS
    !   [none]
    subroutine initialize_availability_data(interpolator)
        type(coamps_interpolator), intent(inout) :: interpolator

        integer :: num_vars_needed

        character(len=*), parameter :: routine = 'initialize_availability_data'
        integer                     :: alloc_status

        select case (interpolator%interp_level_type)
        case (INTERPOLATE_TO_PRESSURE)
            num_vars_needed = NUM_VARS_NEEDED_FOR_P_INTERP
        case (INTERPOLATE_TO_HEIGHT)
            num_vars_needed = NUM_VARS_NEEDED_FOR_Z_INTERP
        case (INTERPOLATE_TO_SIGMA)
            num_vars_needed = NUM_VARS_NEEDED_FOR_S_INTERP
        case (INTERPOLATE_TO_SURFACE)
            num_vars_needed = NUM_VARS_NEEDED_FOR_SLP
        case (INTERPOLATE_TO_UNDEF)
            num_vars_needed = SINGLE_LEVEL
        !case (INTERPOLATE_TO_OTHER)
        !case default
        end select

        allocate(interpolator%vars_available(num_model_levels, num_vars_needed),   &
                 stat=alloc_status)
        call check_alloc_status(alloc_status, routine, source, revision,           &
                                revdate, 'vars_available')
        allocate(interpolator%levels_available(num_model_levels), stat=alloc_status)
        call check_alloc_status(alloc_status, routine, source, revision, &
                                revdate, 'levels_available')
    end subroutine initialize_availability_data

    ! deallocate_availability_data
    ! ----------------------------
    ! Clean up the dynamic storage for the availability information
    ! This needs to be done every time interpolate() is called since it
    ! is potentially different each time
    subroutine deallocate_availability_data(interpolator)
        type(coamps_interpolator), intent(inout) :: interpolator

        character(len=*), parameter :: routine = 'deallocate_availability_data'
        integer :: dealloc_status

        deallocate(interpolator%vars_available, stat=dealloc_status)
        call check_dealloc_status(dealloc_status, routine, source,     &
                                  revision, revdate, 'vars_available')
        deallocate(interpolator%levels_available, stat=dealloc_status)
        call check_dealloc_status(dealloc_status, routine, source,     &
                                  revision, revdate, 'levels_available')
    end subroutine deallocate_availability_data

    ! set_interpolation_level_type
    ! ----------------------------
    ! Sets the module-level vertical coordinate type based on the 
    ! observation location
    subroutine set_interpolation_level_type(interpolator, obs_loc)
        type(coamps_interpolator), intent(inout) :: interpolator
        type(location_type),       intent(in)    :: obs_loc

        ! Note that the DART vertical location type of "model level" is
        ! a sigma level in COAMPS
        if      ( vert_is_pressure(obs_loc) ) then
            interpolator%interp_level_type = INTERPOLATE_TO_PRESSURE
        else if ( vert_is_height(obs_loc)   ) then
            interpolator%interp_level_type = INTERPOLATE_TO_HEIGHT
        else if ( vert_is_level(obs_loc)    ) then
            interpolator%interp_level_type = INTERPOLATE_TO_SIGMA
        else if ( vert_is_surface(obs_loc)  ) then
            interpolator%interp_level_type = INTERPOLATE_TO_SURFACE
        else if ( vert_is_undef(obs_loc)  ) then
            interpolator%interp_level_type = INTERPOLATE_TO_UNDEF
        else
            interpolator%interp_level_type = INTERPOLATE_TO_OTHER
        end if
    end subroutine set_interpolation_level_type

    ! get_target_var
    ! --------------
    ! Generates the data for the target variable by either reading it
    ! from the state vector (for most items) or calculating it
    ! (geopotential height)
    subroutine get_target_var(interpolator, obs_kind)
        type(coamps_interpolator), intent(inout) :: interpolator
        integer,                   intent(in)    :: obs_kind

        integer :: AVAILABLE_INDEX_TARGET

        ! If we're searching for something in the state vector, assume
        ! that it is on mass levels and is *not* a mean state variable
        logical, parameter :: IS_NOT_MEAN  = .false.
        logical, parameter :: IS_M_LEVEL   = .true.

        AVAILABLE_INDEX_TARGET = get_next_availability_index()

        ! Not all the observed variables must come from the state vector
        select case (obs_kind)
          case (QTY_GEOPOTENTIAL_HEIGHT)
              call calculate_heights(interpolator, AVAILABLE_INDEX_TARGET,      &
                                     interpolator%target_values)
          case (QTY_VERTLEVEL)
             interpolator%target_values  =                             &
                spread(get_domain_msigma(interpolator%model_domain),   &
                       VALUES_DIM_NEIGHBOR, NUM_NEIGHBORS)

             interpolator%vars_available(:, AVAILABLE_INDEX_TARGET) = .true.

          case default

              call get_matching_var_values(interpolator, obs_kind, IS_NOT_MEAN, &
                                           IS_M_LEVEL, AVAILABLE_INDEX_TARGET,  &
                                           interpolator%target_values)

        end select
    end subroutine get_target_var

    ! get_coordinate_vars
    ! -------------------
    ! Calculate the vertical coordinate conditioned on the type of level
    ! (e.g. pressure, height, sigma level) we are interpolating to
    subroutine get_coordinate_vars(interpolator)
        type(coamps_interpolator), intent(inout) :: interpolator
    
        integer :: AVAILABLE_INDEX_PRESSURE
        integer :: AVAILABLE_INDEX_HEIGHT
        integer :: AVAILABLE_INDEX_SIGMA

        select case (interpolator%interp_level_type)
        case (INTERPOLATE_TO_PRESSURE)
            AVAILABLE_INDEX_PRESSURE   = get_next_availability_index()

            call calculate_pressure(interpolator, AVAILABLE_INDEX_PRESSURE,  &
                                    interpolator%vcoord_values )

        case (INTERPOLATE_TO_HEIGHT)
            AVAILABLE_INDEX_HEIGHT = get_next_availability_index()

            call calculate_heights(interpolator, AVAILABLE_INDEX_HEIGHT,     &
                                   interpolator%vcoord_values)

        case (INTERPOLATE_TO_SIGMA)
            AVAILABLE_INDEX_SIGMA = get_next_availability_index()

            call calculate_sigma(interpolator, AVAILABLE_INDEX_SIGMA,        &
                                 interpolator%vcoord_values)
        case (INTERPOLATE_TO_SURFACE)
            ! Skip handling this for now
        case (INTERPOLATE_TO_UNDEF)
            ! Skip handling this for now
        case (INTERPOLATE_TO_OTHER)
            ! Skip handling this for now
        case default
            call error_handler(E_ERR, 'get_coordinate_vars', 'Interpolate&
            & level type not found!', source, revision, revdate)
        end select
    end subroutine get_coordinate_vars

    ! calculate_surface_heights
    ! ------------------
    ! Calculates the terrain neihboring terrain heights of a point
    subroutine calculate_surface_heights(interpolator)
        type(coamps_interpolator),     intent(inout)  :: interpolator

        character(len=*), parameter :: routine = 'calculate_surface_heights'

        call get_terrain_height_at_points(get_nest(interpolator%interp_point),                      &
                                          interpolator%neighbors_i, interpolator%neighbors_j,       &
                                          interpolator%vinterp_values(:, SINGLE_LEVEL) )
    end subroutine calculate_surface_heights

    ! calculate_surface_pressure
    ! ------------------
    ! Calculate the surface pressure, keeping track of which
    ! levels have pressure information available
    subroutine calculate_surface_pressure(interpolator)
        type(coamps_interpolator),     intent(inout)  :: interpolator

        real(kind=r8), dimension(NUM_NEIGHBORS) :: zsfc
        real(kind=r8), dimension(SINGLE_POINT, SINGLE_POINT) :: sfc_pres 

        real(kind=r8), dimension(:,:), pointer :: mean_exner_values
        real(kind=r8), dimension(:,:), pointer :: mean_theta_values
        real(kind=r8), dimension(:,:), pointer :: theta_values
        real(kind=r8), dimension(:,:), pointer :: exner_values

        logical, parameter :: IS_NOT_MEAN  = .false.
        logical, parameter :: IS_MEAN      = .true.
        logical, parameter :: IS_M_LEVEL   = .true.
        logical, parameter :: IS_W_LEVEL   = .false.

        integer :: n, k

        character(len=*), parameter :: routine = 'calculate_surface_pressure'
        integer                     :: alloc_status
        integer                     :: dealloc_status

        allocate(exner_values(NUM_NEIGHBORS, num_model_levels), stat=alloc_status)
        call check_alloc_status(alloc_status, routine, source, revision, revdate, 'exner_values')

        allocate(theta_values(NUM_NEIGHBORS, num_model_levels), stat=alloc_status)
        call check_alloc_status(alloc_status, routine, source, revision, revdate, 'theta_values')

        allocate(mean_exner_values(NUM_NEIGHBORS, num_model_levels + 1), stat=alloc_status)
        call check_alloc_status(alloc_status, routine, source, revision, revdate, 'mean_exner_values')

        allocate(mean_theta_values(NUM_NEIGHBORS, num_model_levels), stat=alloc_status)
        call check_alloc_status(alloc_status, routine, source, revision, revdate, 'mean_theta_values')

        call get_matching_var_values(interpolator, QTY_EXNER_FUNCTION, IS_NOT_MEAN, &
                                     IS_M_LEVEL, matching_values = exner_values)

        call get_matching_var_values(interpolator, QTY_POTENTIAL_TEMPERATURE, IS_NOT_MEAN, &
                                     IS_M_LEVEL, matching_values = theta_values)

        call get_matching_var_values(interpolator, QTY_EXNER_FUNCTION, IS_MEAN, &
                                     IS_W_LEVEL, matching_values = mean_exner_values)

        call get_matching_var_values(interpolator, QTY_POTENTIAL_TEMPERATURE, IS_MEAN, &
                                     IS_M_LEVEL, matching_values = mean_theta_values)

        call get_terrain_height_at_points(get_nest(interpolator%interp_point),                      &
                                          interpolator%neighbors_i, interpolator%neighbors_j, zsfc) 

        do n=1,NUM_NEIGHBORS

          call sfcp(theta_values(n,:), exner_values(n,:),              &
                    mean_theta_values(n,:), mean_exner_values(n,:),    &
                    get_domain_dsigmaw(interpolator%model_domain),     &
                    get_domain_wsigma(interpolator%model_domain),      &
                    zsfc(n), SINGLE_POINT, SINGLE_POINT,               &
                    num_model_levels, sfc_pres)

                    interpolator%vinterp_values(n, SINGLE_LEVEL) =     &
                       sfc_pres(SINGLE_POINT, SINGLE_POINT) * CONVERT_MB_TO_PA

        end do

        deallocate(exner_values, stat=dealloc_status)
        call check_dealloc_status(dealloc_status, routine, source, revision, revdate, 'exner_values')

        deallocate(theta_values, stat=dealloc_status)
        call check_dealloc_status(dealloc_status, routine, source, revision, revdate, 'theta_values')

        deallocate(mean_exner_values, stat=dealloc_status)
        call check_dealloc_status(dealloc_status, routine, source, revision, revdate, 'mean_exner_values')

        deallocate(mean_theta_values, stat=dealloc_status)
        call check_dealloc_status(dealloc_status, routine, source, revision, revdate, 'mean_theta_values')

    end subroutine calculate_surface_pressure

    ! calculate_undef_level_var
    ! ------------------
    ! Calculate the surface pressure, keeping track of which
    ! levels have pressure information available
    subroutine calculate_undef_level_var(interpolator, var_kind, vert_value, is_succes)
        type(coamps_interpolator), intent(inout)  :: interpolator
        integer,                   intent(in)     :: var_kind
        real(kind=r8),             intent(in)     :: vert_value
        logical,                   intent(out)    :: is_succes

        real(kind=r8), dimension(NUM_NEIGHBORS)   :: matching_values
        logical,       dimension(NUM_NEIGHBORS)   :: is_available_var

        character(len=*), parameter :: LEVEL_TYPE = 'U'
        integer :: n

        character(len=*), parameter :: routine = 'calculate_undef_level_var'
        integer                     :: alloc_status
        integer                     :: dealloc_status

        call get_single_level_var_values(interpolator, var_kind, LEVEL_TYPE,   &
                                         vert_value, matching_values, is_succes)
        if(.not. is_succes) return

        do n=1,NUM_NEIGHBORS
           interpolator%vinterp_values(n, SINGLE_LEVEL) =  matching_values(n)
        end do

    end subroutine calculate_undef_level_var

    ! calculate_height_level_var
    ! ------------------
    ! Calculate the surface pressure, keeping track of which
    ! levels have pressure information available
    subroutine calculate_height_level_var(interpolator, var_kind, vert_value, is_succes)
        type(coamps_interpolator), intent(inout)  :: interpolator
        integer,                   intent(in)     :: var_kind
        real(kind=r8),             intent(in)     :: vert_value
        logical,                   intent(out)    :: is_succes

        real(kind=r8), dimension(NUM_NEIGHBORS)   :: matching_values
        logical,       dimension(NUM_NEIGHBORS)   :: is_available_var

        character(len=*), parameter :: LEVEL_TYPE = 'Z'
        integer :: n

        character(len=*), parameter :: routine = 'calculate_height_level_var'
        integer                     :: alloc_status
        integer                     :: dealloc_status

        call get_single_level_var_values(interpolator, var_kind, LEVEL_TYPE,   &
                                         vert_value, matching_values, is_succes)
        if(.not. is_succes) return

        do n=1,NUM_NEIGHBORS
           interpolator%vinterp_values(n, SINGLE_LEVEL) =  matching_values(n)
        end do

    end subroutine calculate_height_level_var

    ! calculate_pressure
    ! ------------------
    ! Calculate the pressure on mass sigma levels, keeping track of which
    ! levels have pressure information available
    subroutine calculate_pressure(interpolator, availability_index, pressure)
        type(coamps_interpolator),     intent(inout)  :: interpolator
        integer,                       intent(in)     :: availability_index
        real(kind=r8), dimension(:,:), intent(out)    :: pressure        

        real(kind=r8), dimension(:,:), pointer :: mean_exner_values
        real(kind=r8), dimension(:,:), pointer :: pert_exner_values

        integer :: mean_availability_index
        integer :: pert_availability_index

        character(len=*), parameter :: routine = 'calculate_pressure'
        integer                     :: alloc_status
        integer                     :: dealloc_status

        allocate(mean_exner_values(NUM_NEIGHBORS, num_model_levels), stat=alloc_status)
        call check_alloc_status(alloc_status, routine, source, revision, revdate, 'mean_exner_values')

        allocate(pert_exner_values(NUM_NEIGHBORS, num_model_levels), stat=alloc_status)
        call check_alloc_status(alloc_status, routine, source, revision, revdate, 'pert_exner_values')

        mean_availability_index = get_next_availability_index()
        call get_mean_exner_values(interpolator, mean_availability_index, mean_exner_values)

        pert_availability_index = get_next_availability_index()
        call get_pert_exner_values(interpolator, pert_availability_index, pert_exner_values )

        call convert_exner_to_pressure(interpolator, availability_index, &
                                       mean_exner_values,                &
                                       mean_availability_index,          &
                                       pert_exner_values,                &
                                       pert_availability_index, pressure)
                                       
        deallocate(mean_exner_values, stat=dealloc_status)
        call check_dealloc_status(dealloc_status, routine, source, revision, revdate, 'pert_exner_values')

        deallocate(pert_exner_values, stat=dealloc_status)
        call check_dealloc_status(dealloc_status, routine, source, revision, revdate, 'pert_exner_values')

    end subroutine calculate_pressure

    ! calculate_heights
    ! -----------------
    ! Use a nest's terrain height and the msigma and wsigma values from the 
    ! COAMPS vertical distribution to calculate the geometric height of each 
    ! (mass) sigma level.
    subroutine calculate_heights(interpolator, availability_index, heights)
        type(coamps_interpolator),     intent(inout) :: interpolator
        integer,                       intent(in)    :: availability_index
        real(kind=r8), dimension(:,:), intent(out)   :: heights

        ! Surface height
        real(kind=r8), dimension(NUM_NEIGHBORS) :: Zs

        ! Total depth of the model atmosphere
        real(kind=r8) :: H

        ! COAMPS indexes vertical levels from the top down
        integer, parameter :: MODEL_TOP = 1

        H = get_domain_wsigma(interpolator%model_domain, MODEL_TOP)

        call get_terrain_height_at_points(get_nest(interpolator%interp_point),&
                                          interpolator%neighbors_i,           &
                                          interpolator%neighbors_j,           &
                                          Zs) 

        ! sigma = H * (Z-Zs) / (H-Zs), so
        ! Z = (sigma * (1 - Zs/H)) + Zs
        heights =                                                          &
            (                                                              &
                spread(get_domain_msigma(interpolator%model_domain),       &
                       VALUES_DIM_NEIGHBOR, NUM_NEIGHBORS)                 &
              * (1 - (spread(Zs, VALUES_DIM_LEVEL, NUM_MODEL_LEVELS) / H)) &
            )                                                              &
           +                                                               &
            (                                                              &
                spread(Zs, VALUES_DIM_LEVEL, NUM_MODEL_LEVELS)             &
            )

        ! Since this data comes from the model domain, it's guaranteed to be
        ! available at all levels
        interpolator%vars_available(:, availability_index) = .true.

    end subroutine calculate_heights

    ! calculate_sigma
    ! ---------------
    ! Return the sigma values at each mass sigma level
    subroutine calculate_sigma(interpolator, availability_index, sigma)
        type(coamps_interpolator),     intent(inout) :: interpolator
        integer,                       intent(in)    :: availability_index
        real(kind=r8), dimension(:,:), intent(out)   :: sigma

        sigma = spread(get_domain_msigma(interpolator%model_domain),&
                       VALUES_DIM_NEIGHBOR, NUM_NEIGHBORS)

        ! Since this data comes from the model domain, it's guaranteed to be
        ! available at all levels
        interpolator%vars_available(:, availability_index) = .true.

    end subroutine calculate_sigma

    ! get_mean_theta_values
    ! ---------------------
    ! Get the value of the mean potential on all available mass levels. 
    subroutine get_mean_theta_values(interpolator, availability_index, mean_exner)
        type(coamps_interpolator),     intent(inout)  :: interpolator
        integer,                       intent(in)     :: availability_index
        real(kind=r8), dimension(:,:), intent(out)    :: mean_exner
        
        logical, parameter :: IS_MEAN_VARIABLE = .true.
        logical, parameter :: IS_ON_MASS_LEVEL = .true.

        call get_matching_var_values(interpolator, QTY_POTENTIAL_TEMPERATURE,       &
                                     IS_MEAN_VARIABLE, IS_ON_MASS_LEVEL,      &
                                     availability_index, mean_exner)
    end subroutine get_mean_theta_values

    ! get_mean_exner_values
    ! ---------------------
    ! Get the value of the mean exner function on all available mass levels. 
    subroutine get_mean_exner_values(interpolator, availability_index, mean_exner)
        type(coamps_interpolator),     intent(inout)  :: interpolator
        integer,                       intent(in)     :: availability_index
        real(kind=r8), dimension(:,:), intent(out)    :: mean_exner
        
        logical, parameter :: IS_MEAN_VARIABLE = .true.
        logical, parameter :: IS_ON_MASS_LEVEL = .true.

        call get_matching_var_values(interpolator, QTY_EXNER_FUNCTION,       &
                                     IS_MEAN_VARIABLE, IS_ON_MASS_LEVEL,      &
                                     availability_index, mean_exner)
    end subroutine get_mean_exner_values

    ! get_pert_exner_values
    ! ---------------------
    ! Get the value of the perturbation exner function on all available 
    ! mass levels.
    subroutine get_pert_exner_values(interpolator, availability_index, pert_exner)
        type(coamps_interpolator),     intent(inout)  :: interpolator
        integer,                       intent(in)     :: availability_index
        real(kind=r8), dimension(:,:), intent(out)    :: pert_exner

        ! Search parameters for this variable
        logical, parameter :: IS_NOT_MEAN_VARIABLE = .false.
        logical, parameter :: IS_ON_MASS_LEVEL     = .true.

        call get_matching_var_values(interpolator, QTY_EXNER_FUNCTION,       &
                                     IS_NOT_MEAN_VARIABLE, IS_ON_MASS_LEVEL,  &
                                     availability_index, pert_exner)
    end subroutine get_pert_exner_values

    ! convert_exner_to_pressure
    ! -------------------------
    ! Uses the module-level exner values to convert to pressure
    !  exner = (p/p00)^(R/Cp)
    subroutine convert_exner_to_pressure(interpolator,                &
                                       pres_availability_index,       &
                                       mean_exner_values,             &
                                       mean_availability_index,       &
                                       pert_exner_values,             &
                                       pert_availability_index,       &
                                       pressure)
        type(coamps_interpolator),     intent(inout) :: interpolator
        integer,                       intent(in)    :: pres_availability_index
        real(kind=r8), dimension(:,:), intent(in)    :: mean_exner_values
        integer,                       intent(in)    :: mean_availability_index
        real(kind=r8), dimension(:,:), intent(in)    :: pert_exner_values
        integer,                       intent(in)    :: pert_availability_index
        real(kind=r8), dimension(:,:), intent(out)   :: pressure

        ! Atmospheric constants - specific heat, gas constant (dry air),
        ! and the initial pressure for calculating the Exner function
        real(kind=r8), parameter :: R   = real(287.0,  kind=r8)
        real(kind=r8), parameter :: Cp  = real(1004.0, kind=r8)
        real(kind=r8), parameter :: P00 = real(1000.0, kind=r8)

        pressure = ((mean_exner_values + pert_exner_values) ** (Cp/R)) * P00

        ! That will calculate even if we have junk data - make sure
        ! that the availability matrix reflects that we have both the
        ! mean and perturbation exner functions and hence, the pressure
        interpolator%vars_available(:, pres_availability_index) =   &
            interpolator%vars_available(:, mean_availability_index) &
            .and.                                                   &
            interpolator%vars_available(:, pert_availability_index)
    end subroutine convert_exner_to_pressure

    ! get_matching_var_values
    ! -----------------------
    ! Reads all available levels of variables matching the specified
    ! kind, level type, and mean type at the neighboring points.
    subroutine get_matching_var_values(interpolator, var_kind, is_mean,       &
                                       is_mass_level, var_availability_index, &
                                       matching_values)
        type(coamps_interpolator),      intent(inout)  :: interpolator
        integer,                        intent(in)     :: var_kind
        logical,                        intent(in)     :: is_mean
        logical,                        intent(in)     :: is_mass_level
        integer, optional,              intent(in)     :: var_availability_index
        real(kind=r8), dimension(:,:),  intent(out)    :: matching_values

        integer                       :: num_levels, cur_level_num
        type(state_variable), pointer :: matching_var
        real(kind=r8), dimension(NUM_NEIGHBORS) :: neighbors
        real(kind=r8), allocatable, dimension(:)       :: var_field

        character(len=*), parameter :: MASS_LEVEL = 'M'
        character(len=*), parameter :: W_LEVEL    = 'W'
        character(len=1)            :: level_type

        character(len=*), parameter :: routine = 'get_matching_var_values'
        integer                     :: alloc_status
        integer                     :: dealloc_status

        nullify(matching_var)

        if(is_mass_level) then
          num_levels = num_model_levels
          level_type = MASS_LEVEL
        else
          num_levels = num_model_levels + 1
          level_type = W_LEVEL
        end if

        do cur_level_num = 1, num_levels
            matching_var => find_state_variable(interpolator%state_definition, &
                                               interpolator%interp_nest,       &
                                               var_kind,  is_mean,             &
                                                 level_type, cur_level_num)

            if(present(var_availability_index)) then
              if (.not. associated(matching_var)) then
                call mark_unavailable(interpolator, cur_level_num, var_availability_index)
                cycle
              end if
            end if

            call extract_neighbors(interpolator, matching_var, var_kind, neighbors)

            matching_values(:, cur_level_num) = neighbors(:)

            if(present(var_availability_index)) &
              call mark_available(interpolator, cur_level_num, var_availability_index)
        end do

        nullify(matching_var)

    end subroutine get_matching_var_values

    ! get_single_level_var_values
    ! -----------------------
    ! Reads all available levels of variables matching the specified
    ! kind, level type, and mean type at the neighboring points.
    subroutine get_single_level_var_values(interpolator, var_kind, level_type,   &
                                           vert_value, matching_values, is_available)
        type(coamps_interpolator),      intent(inout)  :: interpolator
        integer,                        intent(in)     :: var_kind
        character(len=*),               intent(in)     :: level_type
        real(kind=r8),                  intent(in)     :: vert_value
        real(kind=r8), dimension(:),    intent(out)    :: matching_values
        logical,                        intent(out)    :: is_available

        type(state_variable), pointer :: matching_var
        logical, parameter :: is_mean = .false.
        integer, parameter :: SINGLE_LEVEL = 1

        nullify(matching_var)
        matching_var => find_state_variable(interpolator%state_definition,  &
                                            interpolator%interp_nest,       &
                                            var_kind,  is_mean,             &
                                            level_type, SINGLE_LEVEL,       &
                                            vert_value = vert_value)
        if(.not.associated(matching_var)) then
          is_available = .false.
          return
        end if

        call extract_neighbors(interpolator, matching_var, var_kind, matching_values)
        is_available = .true.

    end subroutine get_single_level_var_values

    ! destagger_variable
    ! ---------------
    ! Given a 1D field and an observation type, destaggers the field
    ! if necessary using the COAMPS destaggering functions.
    subroutine destagger_variable(interpolator, var_values, var_kind, var_field)
        type(coamps_interpolator),           intent(in)  :: interpolator
        integer,                             intent(in)  :: var_kind
        real(kind=r8), dimension(:), intent(in)    :: var_values
        real(kind=r8), dimension(:), intent(inout) :: var_field

        ! Don't strictly need this but it makes typing much easier
        type(coamps_nest)                      :: nest

        ! Since we're processing fields, there is only 1 vertical level
        integer, parameter :: NUM_VERT_LEVELS = 1

        character(len=*), parameter :: routine = 'destagger_variable'
        integer                     :: alloc_status
        integer                     :: dealloc_status

        nest = get_nest(interpolator%interp_point)

        if ( (var_kind .eq. QTY_U_WIND_COMPONENT)  .or. &
             (var_kind .eq. QTY_V_WIND_COMPONENT)) then
                var_field = var_values
        end if

        select case (var_kind)
        case (QTY_U_WIND_COMPONENT)
            call utom(var_field, get_nest_i_width(nest),     &
                      get_nest_j_width(nest), NUM_VERT_LEVELS, .true.)
        case (QTY_V_WIND_COMPONENT)
            call vtom(var_field, get_nest_i_width(nest),     &
                      get_nest_j_width(nest), NUM_VERT_LEVELS, .true.)
        case default
            var_field(:) = var_values(:)
        end select

    end subroutine destagger_variable

    ! extract_neighbors_ustag
    ! -----------------
    ! Given an entire field, extracts the neighboring points specified by
    ! their i-j coordinate for variables staggered on the u grid
    subroutine extract_neighbors_ustag(interpolator, var_values, neighbors)
        type(coamps_interpolator),      intent(in)  :: interpolator
        real(kind=r8), dimension(:),    intent(in)  :: var_values
        real(kind=r8), dimension(:),    intent(out) :: neighbors

        integer :: cur_neighbor
        integer :: ii_west, ii_east, indx_west, indx_east

        do cur_neighbor = 1, NUM_NEIGHBORS
          ii_west = interpolator%neighbors_i(cur_neighbor) - 1
          ii_east = interpolator%neighbors_i(cur_neighbor)

          indx_west = nest_index_2d_to_1d(interpolator%interp_nest,    &
                      ii_west, interpolator%neighbors_j(cur_neighbor))

          indx_east = nest_index_2d_to_1d(interpolator%interp_nest,    &
                      ii_east, interpolator%neighbors_j(cur_neighbor))

          neighbors(cur_neighbor) = 0.5_r8 *   &
               (var_values(indx_west) + var_values(indx_east))   
        end do                                           
    end subroutine extract_neighbors_ustag

    ! extract_neighbors_vstag
    ! -----------------
    ! Given an entire field, extracts the neighboring points specified by
    ! their i-j coordinate for variables staggered on the u grid
    subroutine extract_neighbors_vstag(interpolator, var_values, neighbors)
        type(coamps_interpolator),      intent(in)  :: interpolator
        real(kind=r8), dimension(:),    intent(in)  :: var_values
        real(kind=r8), dimension(:),    intent(out) :: neighbors

        integer :: cur_neighbor
        integer :: jj_south, jj_north, indx_north, indx_south

        do cur_neighbor = 1, NUM_NEIGHBORS
          jj_south = interpolator%neighbors_j(cur_neighbor) - 1
          jj_north = interpolator%neighbors_j(cur_neighbor)

          indx_south = nest_index_2d_to_1d(interpolator%interp_nest,    &
                      interpolator%neighbors_i(cur_neighbor), jj_south)

          indx_north = nest_index_2d_to_1d(interpolator%interp_nest,    &
                      interpolator%neighbors_i(cur_neighbor), jj_north)

          neighbors(cur_neighbor) = 0.5_r8 *   &
               (var_values(indx_south) + var_values(indx_north))   
        end do                                           

    end subroutine extract_neighbors_vstag

    ! extract_neighbors_tstag
    ! -----------------
    ! Given an entire field, extracts the neighboring points specified by
    ! their i-j coordinate
    subroutine extract_neighbors_tstag(interpolator, var_values, neighbors)
        type(coamps_interpolator),      intent(in)  :: interpolator
        real(kind=r8), dimension(:),    intent(in)  :: var_values
        real(kind=r8), dimension(:),    intent(out) :: neighbors

        integer :: cur_neighbor

        do cur_neighbor = 1, NUM_NEIGHBORS
            neighbors(cur_neighbor) =                                    &
                var_values(nest_index_2d_to_1d(interpolator%interp_nest, &
                                interpolator%neighbors_i(cur_neighbor),  &
                                interpolator%neighbors_j(cur_neighbor)))
        end do                                           

    end subroutine extract_neighbors_tstag

    ! extract_neighbors
    ! -----------------
    ! Given an entire field, extracts the neighboring points specified by
    ! their i-j-k coordinate
    subroutine extract_neighbors(interpolator, matching_var, var_kind, neighbors)
        type(coamps_interpolator),     intent(in)    :: interpolator
        type(state_variable),          intent(in)    :: matching_var
        integer,                       intent(in)    :: var_kind
        real(kind=r8), dimension(:),   intent(out)   :: neighbors

            ! Take care for U and V wind components.  If the variable is on a sigma level
            ! it is staggered.  If it is on any other vert_type it is unstaggered.
            select case (var_kind)
              case (QTY_U_WIND_COMPONENT)

                if(is_sigma_level(matching_var)) then
                call extract_neighbors_ustag(interpolator,                              &
                     get_var_substate(matching_var, interpolator%model_state), neighbors)
                else
                  call extract_neighbors_tstag(interpolator,                              &
                       get_var_substate(matching_var, interpolator%model_state), neighbors)
                end if

              case (QTY_V_WIND_COMPONENT)

                if(is_sigma_level(matching_var)) then
                call extract_neighbors_vstag(interpolator,                              &
                     get_var_substate(matching_var, interpolator%model_state), neighbors)
                else
                  call extract_neighbors_tstag(interpolator,                              &
                       get_var_substate(matching_var, interpolator%model_state), neighbors)
                end if

              case default

                call extract_neighbors_tstag(interpolator,                              &
                     get_var_substate(matching_var, interpolator%model_state), neighbors)

            end select

    end subroutine extract_neighbors

    ! mark_available
    ! --------------
    ! Given a level and variable index, marks the variable available
    ! at that particular level
    subroutine mark_available(interpolator, level, var_index)
        type(coamps_interpolator), intent(inout) :: interpolator
        integer,                   intent(in)    :: level
        integer,                   intent(in)    :: var_index

        interpolator%vars_available(level, var_index) = .true.
    end subroutine mark_available

    ! mark_unavailable
    ! ----------------
    ! Given a level and variable index, marks the variable unavailable
    ! at that particular level
    subroutine mark_unavailable(interpolator, level, var_index)
        type(coamps_interpolator), intent(inout) :: interpolator
        integer,                   intent(in)    :: level
        integer,                   intent(in)    :: var_index

        interpolator%vars_available(level, var_index) = .false.
    end subroutine mark_unavailable

    ! calculate_available_levels
    ! --------------------------
    ! Once all the availability statistics have been populated, this routine 
    ! will calculate the levels that have all necessary data present
    subroutine calculate_available_levels(interpolator)
        type(coamps_interpolator), intent(inout) :: interpolator

        ! A level is defined as available if every variable required
        ! for the interpolation is available at that level
        interpolator%levels_available(:) = all(interpolator%vars_available, &
                                               AVAILABILITY_DIM_VARIABLE)

        interpolator%num_levels_available = &
            count(interpolator%levels_available)

    end subroutine calculate_available_levels

    ! initialize_available_values
    ! ---------------------------
    ! Runtime initialization of available values - this will depend on
    ! the variable being interpolated and the type of level it's being
    ! interpolated to
    subroutine initialize_available_values(interpolator)
        type(coamps_interpolator), intent(inout) :: interpolator

        integer :: num_levels

        character(len=*), parameter :: routine = 'initialize_available_values'
        integer                     :: alloc_status

        num_levels = interpolator%num_levels_available

        allocate(interpolator%available_target_values(NUM_NEIGHBORS, &
                                                      num_levels),   &
                 stat=alloc_status)
        call check_alloc_status(alloc_status, routine, source, revision, &
                                revdate, 'available_target_values')

        allocate(interpolator%available_vcoord_values(NUM_NEIGHBORS, &
                                                      num_levels),   &
                 stat=alloc_status)
        call check_alloc_status(alloc_status, routine, source, revision, &
                                revdate, 'available_vcoord_vals')

    end subroutine initialize_available_values


    ! collect_available_values
    ! ------------------------
    subroutine collect_available_values(interpolator)
        type(coamps_interpolator), intent(inout) :: interpolator

        call get_available_values(interpolator, interpolator%target_values,  &
                                  interpolator%available_target_values)

        call get_available_values(interpolator, interpolator%vcoord_values,  &
                                  interpolator%available_vcoord_values )

    end subroutine collect_available_values


    ! interpolate_to_level
    ! --------------------
    ! Interpolates an array to a single level using the interpolation
    ! functions from the COAMPS utility package.
    subroutine interpolate_to_level(interpolator)
        type(coamps_interpolator), intent(inout) :: interpolator

        ! Array size information for the interpolators
        integer :: num_levels_in
        integer :: num_levels_out

        num_levels_in  = interpolator%num_levels_available
        num_levels_out = 1

        select case (interpolator%interp_level_type)
        case (INTERPOLATE_TO_PRESSURE)
            ! DART location module stores pressure in Pa and the COAMPS
            ! pressure calculated from the Exner function is in mb
            interpolator%vinterp_level =  interpolator%interp_level &
                                        * CONVERT_PA_TO_MB

            call s2pint(interpolator%available_target_values,         & 
                        interpolator%vinterp_values,                  &
                        interpolator%available_vcoord_values,         &
                        interpolator%vinterp_level,                   &
                        num_levels_in, num_levels_out, NUM_NEIGHBORS, &
                        USE_MISSING_VALUE, MISSING_R8)
        case (INTERPOLATE_TO_HEIGHT, INTERPOLATE_TO_SIGMA)
            ! DART location module stores height in meters so no conversion
            ! is necessary
            interpolator%vinterp_level = interpolator%interp_level

            ! Same function call will work on either sigma interpolation
            ! or height interpolation - only difference is the values in
            ! the vertical coordinate variable
            call z2zint(interpolator%available_target_values,         & 
                        interpolator%vinterp_values,                  &
                        interpolator%available_vcoord_values,         &
                        interpolator%vinterp_level,                   &
                        num_levels_in, num_levels_out, NUM_NEIGHBORS, &
                        USE_MISSING_VALUE, MISSING_R8)
        end select
    end subroutine interpolate_to_level

    ! interpolate_to_point
    ! --------------------
    ! Interpolates an array to a single point.  Uses module-level
    ! variables for arrays holding the data and the weights
    !  PARAMETERS
    !   OUT obs_value         result of interpolation
    function interpolate_to_point(interpolator)
        type(coamps_interpolator), intent(in)  :: interpolator
        real(kind=r8)                          :: interpolate_to_point

        ! Need to remove the singleton dimension of vinterp_values since 
        ! dot_product wants a 1-d vector for both arguments
        interpolate_to_point = dot_product(interpolator%interp_weights,      &
                                           pack(interpolator%vinterp_values, &
                                                .true.))
    end function interpolate_to_point

    ! get_next_availability_index
    ! ---------------------------
    ! Return where to write a variable's availability within the availability 
    ! array and change the master index to reflect the request.
    function get_next_availability_index()
        integer :: get_next_availability_index

        cur_availability_index      = cur_availability_index + 1
        get_next_availability_index = cur_availability_index
    end function get_next_availability_index

    ! get_available_values
    ! --------------------
    ! Condense the raw values on all levels to only those levels that are 
    ! available for *all* variables
    subroutine get_available_values(interpolator, raw_values, available_values)
        type(coamps_interpolator),     intent(in)    :: interpolator
        real(kind=r8), dimension(:,:), intent(in)    :: raw_values
        real(kind=r8), dimension(:,:), intent(inout) :: available_values

        ! This works since the size of the first array dimension (i.e. the
        ! number of neighbors) does not change size, so the result is only
        ! the values selected by the levels_available mask after it is tiled
        ! across all neighbors)  
        available_values(:,:) =                                            &
            reshape(pack(raw_values, spread(interpolator%levels_available, &
                                            VALUES_DIM_NEIGHBOR,           &
                                            NUM_NEIGHBORS)),               &
                    (/ NUM_NEIGHBORS, interpolator%num_levels_available /))
    end subroutine get_available_values

    ! is_valid_level_type
    ! ----------------
    ! Ensures that the level type defined by the observation location
    ! is one that we actually will interpolate to
    function is_valid_level_type(interp_level_type)
        integer, intent(in)  :: interp_level_type
        logical              :: is_valid_level_type

        character(len=128) :: message

        if (interp_level_type .eq. INTERPOLATE_TO_PRESSURE .or.&
            interp_level_type .eq. INTERPOLATE_TO_HEIGHT   .or.&
            interp_level_type .eq. INTERPOLATE_TO_SURFACE  .or.&
            interp_level_type .eq. INTERPOLATE_TO_UNDEF    .or.&
            interp_level_type .eq. INTERPOLATE_TO_SIGMA ) then
            is_valid_level_type = .true.
        else
            write (message, '(A,I2,1x,A)') 'Level type ', interp_level_type, &
                               'is not supported.'
            call error_handler(E_MSG, 'is_valid_level_type', trim(message),&
                               source, revision, revdate)
            is_valid_level_type = .false.
        end if
    end function is_valid_level_type

    ! enough_levels_available
    ! -----------------------
    ! Ensures that we have enough vertical levels to do a meaningful
    ! interpolation.  Since we're currently not supporting surface
    ! level types, use the module-level parameter for the threshold
    function enough_levels_available(interpolator)
        type(coamps_interpolator), intent(in)  :: interpolator
        logical                                :: enough_levels_available

        character(len=128) :: message

        if (interpolator%num_levels_available < MIN_LEVELS_NEEDED) then
            write (message,'(3(A,1x,I2,1x))') 'There are only',            &
                              interpolator%num_levels_available,          &
                              'levels available, which is less than the', &
                              MIN_LEVELS_NEEDED, 'needed to interpolate.'
            call error_handler(E_MSG, 'enough_levels_available', &
                               trim(message), source, revision, revdate)
            enough_levels_available = .false.
        else
            enough_levels_available = .true.
        end if
    end function enough_levels_available

    ! interp_level_in_available_range
    ! -------------------------------
    ! Ensures that the vertical level we are trying to interpolate
    ! to lies within the available levels that we have for the
    ! appropriate vertical coordinate type.
    function interp_level_in_available_range(interpolator)
        type(coamps_interpolator), intent(in) :: interpolator
        logical                               :: interp_level_in_available_range

        character(len=256) :: message

        real(kind=r8) :: min_maxlevel_available, max_minlevel_available
        real(kind=r8) :: scaled_level

        min_maxlevel_available = minval(maxval(interpolator%available_vcoord_values,2))
        max_minlevel_available = maxval(minval(interpolator%available_vcoord_values,2))

        ! If we're interpolating onto pressure levels, DART stores
        ! the pressure in Pa and we calculate in hPa
        if (interpolator%interp_level_type == INTERPOLATE_TO_PRESSURE) then
            scaled_level = interpolator%interp_level * CONVERT_PA_TO_MB
        else
            scaled_level = interpolator%interp_level
        end if

        ! We're OK if the level is below the highest available and
        ! higher than the lowest available
        if (       (scaled_level <= min_maxlevel_available)       &
             .and. (scaled_level >= max_minlevel_available)) then 
            interp_level_in_available_range = .true.
        else
            !write (message,*) 'Trying to interpolate to vertical level',  &
            !                  scaled_level, 'but the available range is', &
            !                  max_minlevel_available, 'to', min_maxlevel_available
            !call error_handler(E_MSG, 'level_is_in_available_range',      &
            !                   trim(message), source, revision, revdate)
            interp_level_in_available_range = .false.
        end if

    end function interp_level_in_available_range

    ! no_missing_values
    ! -----------------
    ! Looks at the results of the level interpolation and returns true if 
    ! there are no values labeled as "missing"
    function no_missing_values(interpolator)
        type(coamps_interpolator), intent(in)  :: interpolator
        logical                                :: no_missing_values

        character(len=128) :: message

        integer :: num_missing_values

        ! Yes, I know I'm using an equality operator with floating points...
        num_missing_values = count(interpolator%vinterp_values == MISSING_R8)

        if (num_missing_values == 0) then
            no_missing_values = .true.
        else
            write (message,*) 'Missing values were found after level ', &
                              'interpolation'
            call error_handler(E_WARN, 'no_missing_values',             &
                               trim(message), source, revision, revdate)
            no_missing_values = .false.
        end if
    end function no_missing_values

    ! print_interpolation_diagnostics
    ! -------------------------------
    ! Print some data about the interpolation
    !  PARAMETERS
    !   IN  obs_kind          Type of the observation (integer)
    !   IN  obs_value         Result of the interpolation
    subroutine print_interpolation_diagnostics(interpolator, obs_kind, obs_value)
        type(coamps_interpolator), intent(in) :: interpolator
        integer,                   intent(in) :: obs_kind
        real(kind=r8),             intent(in) :: obs_value

        real(kind=r8)       :: lat, lon

        call nest_point_to_latlon(interpolator%model_domain, interpolator%interp_point, lat, lon)

        if (do_output()) then
            call print_label_name('Interpolation Results')
            write (*,'(A15,T25,2(F15.6,1x))') 'I/J Coordinate  :', &
                get_i_coord(interpolator%interp_point)           , &
                get_j_coord(interpolator%interp_point)
            write (*,'(A15,T25,2(F15.6,1x))') 'Lat/Lon Coord   :', &
                lat, lon
            write (*,'(A15,T25,F10.2)') 'Vertical      :', &
                interpolator%interp_level
            write (*,'(A15,T25,I10)'  ) 'Vertical Type :', &
                interpolator%interp_level_type
            write (*,'(A15,T25,I10)'  ) 'On Nest Level :', &
                get_nest_level(interpolator%interp_point)
            write (*,'(A15,T25,A30)'  ) 'Variable Type :', trim(get_name_for_quantity(obs_kind))
            write (*,'(A15,T25,F15.6)') 'Value         :', obs_value
            write (*,*)
        end if
    end subroutine print_interpolation_diagnostics

    ! reset_availability_index
    ! ------------------------
    ! Resets the module-wide availability index for a new interpolation
    subroutine reset_availability_index()
        cur_availability_index = 0
    end subroutine reset_availability_index

    !------------------------------
    ! END PRIVATE ROUTINES
    !------------------------------

end module coamps_interp_mod

! <next few lines under version control, do not edit>
! $URL$
! $Id$
! $Revision$
! $Date$
