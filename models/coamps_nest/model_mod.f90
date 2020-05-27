! This code may (or may not) be part of the COAMPS distribution,
! So it is not protected by the DART copyright agreement.
!
! DART $Id$

module model_mod

!------------------------------
! MODULE:       model_mod
! AUTHOR:       T. R. Whitcomb
!               Naval Research Laboratory
! DART VERSION: Jamaica
!
! Module for the interface between DART and the U. S. Navy's COAMPS
! mesoscale model.  COAMPS was developed by the Naval Research Laboratory,
! Monterey, California.  COAMPS is a registered trademark of the Naval
! Research Laboratory.
!------------------------------ 

    use coamps_domain_mod,   only : coamps_domain,            &
                                    initialize_domain,        &
                                    dump_domain_info,         &
                                    get_domain_nest,          &
                                    get_nest_count,           &
                                    get_domain_num_levels,    &
                                    get_domain_wsigma,        &
                                    get_domain_msigma,        &
                                    get_domain_dsigmaw,       &
                                    latlon_to_nest_point,     &
                                    nest_point_to_latlon


    use coamps_statevec_mod, only : state_vector,             &
                                    state_iterator,           &
                                    get_iterator,             &
                                    has_next,                 &
                                    get_next,                 &
                                    initialize_state_vector,  &
                                    find_state_variable,      &
                                    get_num_fields,           &
                                    get_var_by_index,         &
                                    dump_state_vector,         &
                                    get_total_size


    use coamps_statevar_mod, only : state_variable,           &
                                    get_var_substate,             &
                                    get_nest_number,          &
                                    get_vert_type,            &
                                    get_vert_loc,             &
                                    get_state_begin,          &
                                    get_state_end,            &
                                    is_var_at_index,          &
                                    dump_state_variable,          &
                                    get_var_kind

    use coamps_nest_mod,     only : coamps_nest,              &
                                    get_terrain,              &
                                    get_terrain_height_at_points, &
                                    get_nest_i_width,         &
                                    get_nest_j_width,         &
                                    get_nest_delta_x,         &
                                    get_nest_delta_y,         &
                                    get_nest,                 &
                                    make_nest_point,          &
                                    nest_point,               &
                                    get_nest_latlon,          &
                                    dump_nest_info,          &
                                    nest_index_1d_to_3d

    use coamps_intrinsic_mod, only : vor,              &
                                     z2zint

    use coamps_interp_mod,   only : interpolate,              &
                                    set_interp_diag

    use coamps_util_mod,     only : check_alloc_status,       &
                                    set_debug_level,          &
                                    timestamp_message,        &
                                    dump_data_file,           &
                                    check_dealloc_status

    use coamps_netcdf_mod,   only : nc_write_statearray_atts, &
                                    nc_write_prognostic_atts, &
                                    nc_write_statearray_data, &
                                    nc_write_prognostic_data

    use coamps_pert_mod,     only : perturb_state

    use location_mod,        only : get_close_type,                &
                                    get_dist,                      &
                                    get_location,                  &
                                    location_type,                 &
                                    loc_get_close_maxdist_init  => &
                                        get_close_maxdist_init,    &
                                    loc_get_close_obs           => &
                                        get_close_obs,             &
                                    loc_get_close_obs_init      => &
                                        get_close_obs_init,        &
                                    set_location,                  &
                                    vert_is_pressure,              &
                                  horiz_dist_only,                 &
                                  query_location,                  &
                                  VERTISLEVEL,                     &
                                  VERTISPRESSURE,                  &
                                  VERTISHEIGHT,                    &
                                  VERTISSURFACE,                   &
                                  VERTISUNDEF
    use obs_kind_mod

    use time_manager_mod,    only : set_time,                      &
                                    set_time_missing,              &
                                    time_type

    use types_mod,           only : MISSING_R8,                    &
                                    DEG2RAD,                       &
                                    r8

    use utilities_mod,       only : check_namelist_read,           &
                                    do_output,                     &
                                    E_ERR,                         &
                                    E_MSG,                         &
                                    error_handler,                 &
                                    find_namelist_in_file,         &
                                    get_unit,                      &
                                    register_module

    implicit none

    private

    !------------------------------
    ! BEGIN PUBLIC INTERFACE
    !------------------------------

    ! Initialization/finalization
    public :: static_init_model 
    public :: end_model      

    ! NetCDF
    public :: nc_write_model_atts 
    public :: nc_write_model_vars 

    ! Ensemble generation
    public :: pert_model_state 

    ! Forward operator
    public :: model_interpolate
    public :: ens_mean_for_model

    ! Localization
    public :: get_close_maxdist_init
    public :: get_close_obs_init
    public :: get_close_obs

    ! Information about model setup
    public :: get_model_size 
    public :: get_state_meta_data 
    public :: get_model_time_step 
    public :: get_coamps_domain 

    ! Null interfaces
    public :: init_conditions
    public :: init_time      
    public :: adv_1step 

    !------------------------------
    ! END PUBLIC INTERFACE
    !------------------------------

    !------------------------------
    ! BEGIN EXTERNAL INTERFACE
    !------------------------------
  !  [none]
    !------------------------------
    ! END EXTERNAL INTERFACE
    !------------------------------

    !------------------------------
    ! BEGIN TYPES AND CONSTANTS
    !------------------------------
    !  [none]
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

    ! Main model_mod namelist - not too much here as we read most of
    ! the data we need in from the COAMPS files themselves
    character(len=10) :: cdtg                 = '1999083100' ! Date-time group
    integer           :: y_bound_skip         = 3            ! How many x and y boundary
    integer           :: x_bound_skip         = 3            ! points to skip when
                                                             ! perturbing the model
                                                             ! state
    logical           :: need_mean            = .true.       ! Do we need the ensemble
                                                             ! mean for for forward
                                                             ! operator computation?
    character(len=80) :: dsnrff               = './'         ! Path to data files
    logical           :: output_interpolation = .false. 
    logical           :: output_state_vector  = .false.
    integer           :: debug                = 0            ! increase for debug messages

    namelist /model_nml/ cdtg, y_bound_skip, x_bound_skip, need_mean, dsnrff, &
                         output_interpolation, output_state_vector, debug

    ! Locations of state variables
    integer, dimension(:), allocatable :: all_vars

    ! Grid information structure
    type(coamps_domain) :: domain
    type(state_vector)  :: state_definition
    type(state_vector)  :: state_layout_3D

    ! Ensemble mean
    real(kind=r8), dimension(:), allocatable :: ensemble_mean

    !------------------------------
    ! END MODULE VARIABLES
    !------------------------------

contains

    !------------------------------
    ! BEGIN PUBLIC ROUTINES
    !------------------------------

    ! static_init_model
    ! -----------------
    ! One-time initialization of the model.  For COAMPS, this:
    !  1. Reads in the model_mod namelist
    !  2. Initializes the pressure levels for the state vector
    !  3. Generate the location data for each member of the state
    !  PARAMETERS
    !   [none]
    subroutine static_init_model()

        character(len=*), parameter :: STATE_VEC_DEF_FILE = 'state.vars'

        call register_module(source, revision, revdate)

        call read_model_namelist()
        call set_debug_level(debug)

        call initialize_domain(cdtg, domain)
        call set_interp_diag(output_interpolation)

        call initialize_state_vector(state_definition, STATE_VEC_DEF_FILE, domain)
        call initialize_state_vector(state_layout_3D, STATE_VEC_DEF_FILE, domain, .true.)

        if (do_output()) call dump_state_vector(state_layout_3D)

        if (do_output()) call dump_domain_info(domain)

        call allocate_metadata_arrays()

        call populate_metadata_arrays()

    end subroutine static_init_model

    ! end_model
    ! ---------
    ! Clean up the workspace once the program is ending
    !  PARAMETERS
    !   [none]
    subroutine end_model()
        real(kind=r8) x

        x = MISSING_R8

        call deallocate_metadata_arrays()
    end subroutine end_model

    ! nc_write_model_atts
    ! -------------------
    ! Write model-specific global attributes to a NetCDF file.
    !  PARAMETERS
    !   IN  ncFileID          numeric ID of an *open* NetCDF file
    !   OUT ierr              0 if writing was successful
    function nc_write_model_atts( ncFileID ) result (ierr)
        integer, intent(in)         :: ncFileID
        integer                     :: ierr          

        ! Error handling
        character(len=*), parameter :: routine = 'nc_write_model_atts'

        ! Assume success - the function will abort if there is an error
        ierr = 0

        !FIXME Add support for state vector output
        !if ( output_state_vector ) then
        !  call nc_write_statearray_atts( ncFileID, all_locs, all_types)
        !else
          call nc_write_prognostic_atts( ncFileID, state_layout_3D )
        !endif

    end function nc_write_model_atts

    ! nc_write_model_vars
    ! -------------------
    ! Writes the model variables (for now just the whole state vector)
    ! to the NetCDF file
    !  PARAMETERS
    !   IN  ncFileID          numeric *open* NetCDF file ID
    !   IN  statevec          large DART state vector
    !   IN  copyindex         which 'copy' to write - this is used for
    !                         indexing ensemble members
    !   IN  timeindex         which time to write (as an index)
    !   OUT ierr              0 if successful
    function nc_write_model_vars(ncFileID, statevec, copyindex, timeindex ) result (ierr)          
        integer,                intent(in) :: ncFileID
        real(r8), dimension(:), intent(in) :: statevec
        integer,                intent(in) :: copyindex
        integer,                intent(in) :: timeindex
        integer                            :: ierr   

        ! Error handling
        character(len=*), parameter :: routine = 'nc_write_model_vars'

        ! Assume success - the function will abort if there is an error
        ierr = 0

        if ( output_state_vector ) then
          call nc_write_statearray_data( ncFileID, statevec, copyindex, timeindex)
        else
          call nc_write_prognostic_data( ncFileID, state_layout_3D, statevec, copyindex, timeindex)
        endif

    end function nc_write_model_vars

    ! pert_model_state
    ! ----------------
    ! Perturb the model state, field by field.  This can be done 3
    ! different ways:
    !  1. No perturbation
    !  2. Uniform perturbation - each element of the field has the
    !                            same additive perturbation
    !  3. Individual perturbation - each element of the field has a 
    !                               different additive perturbation
    ! The perturbation magnitude and option are supplied out of the
    ! dynamic restart vector definition - this allows us to supply a
    ! variance appropriate for each type of variable at each level.
    !  PARAMETERS
    !   IN  state             DART state vector
    !   OUT pert_state        state vector after perturbations
    !   OUT interf_provided   true if this routine did the perturbation
    subroutine pert_model_state(state, pert_state, interf_provided)
        real(kind=r8), dimension(:), intent(in)  :: state
        real(kind=r8), dimension(:), intent(out) :: pert_state
        logical,                     intent(out) :: interf_provided

        call error_handler(E_MSG, 'pert_model_state', 'Perturbing '// &
                           'model state', source, revision, revdate) 
        
        pert_state = perturb_state(state, state_definition, x_bound_skip, y_bound_skip)
        interf_provided = .true.
    end subroutine pert_model_state
    
    ! model_interpolate
    ! -----------------
    ! Given the DART state vector, a location, and a raw variable type
    ! interpolates the state to that location
    ! This implementation currently only supports pressure coordinates.
    !  PARAMETERS
    !   IN  x                 full state vector
    !   IN  location          DART location structure for interpolation
    !                         target location
    !   IN  obs_kind          raw variable kind
    !   OUT obs_val           interpolated value
    !   OUT interp_status     status of interpolation (0 is success)
    subroutine model_interpolate(x, location, obs_kind, obs_val, interp_status)
        real(r8), dimension(:), intent(in)  :: x
        type(location_type),    intent(in)  :: location
        integer,                intent(in)  :: obs_kind
        real(r8),               intent(out) :: obs_val
        integer,                intent(out) :: interp_status

        integer :: which_vert
        logical :: interp_worked
        logical :: in_domain

        type(coamps_nest)      :: nest
        type(location_type)    :: cur_loc 
        type(nest_point)       :: nest_pt
        type(state_variable)   :: u_var, v_var

        integer                :: loc_which
        integer                :: i, j, nx, ny, nz
        integer,  dimension(2) :: ij
        real(r8), dimension(3) :: loc_array
        real(r8)               :: delx, dely, ztop

        real(r8), dimension(1) :: terrain

        !real(r8), allocatable  :: heights(:, :)
        real(r8), allocatable  :: heights(:)
        real(r8), allocatable  :: fm(:, :)
        real(r8), allocatable  :: vort_z(:, :, :)
        real(r8), allocatable  :: vort_z_zlev(:, :)
        real(r8), allocatable  :: vort_z_smooth(:, :)
        logical,  allocatable  :: mask(:,:)

        real(r8), parameter    :: vort_srch_radius = 5.0_r8
        real(r8), parameter    :: vrtx_scale = 200000.0_r8
        real(r8), parameter    :: vort_max_lev = 100.0_r8

        integer,  parameter    :: USE_MISSING_VALUE = 1
        integer,  parameter    :: DART_LOC_LON  = 1
        integer,  parameter    :: DART_LOC_LAT  = 2
        integer,  parameter    :: DART_LOC_VERT = 3

        character(len=*), parameter :: routine = 'model_interpolate'
        integer                     :: alloc_status, dealloc_status

        select case (obs_kind)    
        case (QTY_VORTEX_LAT, QTY_VORTEX_LON)

          ! obs_loc
          loc_array = get_location(location)
          loc_which = nint(query_location(location))

          ztop = get_domain_wsigma(domain, 1)

          call latlon_to_nest_point(domain, loc_array(DART_LOC_LAT), loc_array(DART_LOC_LON),&
                                    nest_pt, in_domain)
          if(.not. in_domain) then
            obs_val = MISSING_R8
            interp_status = 1
            return
          end if

          nest = get_nest(nest_pt)

          nx = get_nest_i_width(nest)
          ny = get_nest_j_width(nest)
          nz = get_domain_num_levels(domain)

          delx = get_nest_delta_x(nest)
          dely = get_nest_delta_y(nest)

          ! allocate arrays
          allocate(vort_z(nx, ny, nz), stat=alloc_status )
          call check_alloc_status(alloc_status, routine, source, revision, &
                                  revdate, 'model_interpolate: vort_z')

          allocate(vort_z_zlev(nx, ny), stat=alloc_status )
          call check_alloc_status(alloc_status, routine, source, revision, &
                                  revdate, 'model_interpolate: vort_z_zlev')

          allocate(vort_z_smooth(nx, ny), stat=alloc_status )
          call check_alloc_status(alloc_status, routine, source, revision, &
                                  revdate, 'model_interpolate: vort_z_smooth')

          !allocate(heights(nx*ny, nz), stat=alloc_status )
          allocate(heights(nz), stat=alloc_status )
          call check_alloc_status(alloc_status, routine, source, revision, &
                                  revdate, 'model_interpolate: heights')


          ! (1) Get u and v location in the state vector 
          u_var = find_state_variable(state_layout_3D, nest, &
                   QTY_U_WIND_COMPONENT, .false., 'M', nz)

          v_var = find_state_variable(state_layout_3D, nest, &
                   QTY_V_WIND_COMPONENT, .false., 'M', nz)

          ! (2) Calculate vorticity
          call timestamp_message('VRTX: Before calculate vorticity')
          allocate(fm(nx, ny), stat=alloc_status )
          call check_alloc_status(alloc_status, routine, source, revision, &
                                  revdate, 'model_interpolate: fm')

          fm(:,:) = 1.0_r8

          call vor(reshape(get_var_substate(u_var, x),(/nx, ny, nz/)),                         &
                   reshape(get_var_substate(v_var, x),(/nx, ny, nz/)),                         &
                   nx, ny, nz, delx, dely, fm,                              &
                   get_domain_msigma(domain), get_domain_wsigma(domain),    &
                   get_domain_dsigmaw(domain), get_terrain(nest), 1, vort_z )

          deallocate(fm, stat=dealloc_status)
          call check_dealloc_status(dealloc_status, routine, source, revision, &
                                    revdate, 'model_interpolate: fm')

          call timestamp_message('VRTX: After calculate vorticity')

          ! (3) Compute sigma level heights.  Since we want to interpolate to a 
          ! a height above ground level we will not add the surface height.
          !call timestamp_message('VRTX: Before calculate heights')
          !heights = spread(get_domain_msigma(domain), 1, nx*ny) * &
          !          ( 1 - spread(reshape(get_terrain(nest)/ztop, (/nx*ny/)), 2, nz) )
          !call timestamp_message('VRTX: After calculate heights')

          ! (4) Interpolate vorticity to the deisred level
          call timestamp_message('VRTX: Before interpolate vorticity')
          do i=1,nx ; do j=1,ny

            call get_terrain_height_at_points(nest, (/i/), (/j/), terrain) 

            heights(:) = get_domain_msigma(domain) * (1 - terrain(1)/ztop)

            call z2zint(vort_z(i, j, :), vort_z_zlev(i, j), heights(:),   &
                        (/vort_max_lev/), nz, 1, 1, USE_MISSING_VALUE, MISSING_R8 )

          end do ; end do

!          call z2zint(reshape(vort_z, (/nx*ny, nz/)), reshape(vort_z_zlev, (/nx*ny, nz/)), &
!                      reshape(heights,(/nx*ny, nz/)), (/vort_max_lev/), nz, 1, nx*ny,      &
!                      USE_MISSING_VALUE, MISSING_R8)
          call timestamp_message('VRTX: After interpolate vorticity')

          ! (5) Smooth vorticity
          call timestamp_message('VRTX: Before smooth vorticity')
          call kernal_smoother(vort_z_zlev, vort_z_smooth, &
                               ceiling(vrtx_scale/delx), nx, ny)

          call timestamp_message('VRTX: After smooth vorticity')
         
          ! (6) Search for vorticity maximum within a given radius of obs.
          call timestamp_message('VRTX: Before maximum vorticity search')

          allocate(mask(nx,ny), stat=alloc_status )
          call check_alloc_status(alloc_status, routine, source, revision, &
                                  revdate, 'model_interpolate: mask')

          do i=1,nx ; do j=1,ny
            call nest_point_to_latlon(domain, make_nest_point(nest, i, j),  &
                                      loc_array(DART_LOC_LAT), loc_array(DART_LOC_LON))

            cur_loc = set_location(loc_array(DART_LOC_LON), loc_array(DART_LOC_LAT), &
                                   loc_array(DART_LOC_VERT), loc_which)

            mask(i,j) = get_dist(location, cur_loc, no_vert=.true.) &
                                 <= (vort_srch_radius*DEG2RAD)
          end do ; end do
          ij = maxloc(vort_z_smooth, mask)

          call timestamp_message('VRTX: After maximum vorticity search')

          ! (7) Set desired return value
          call nest_point_to_latlon(domain, make_nest_point(nest, ij(1), ij(2)),  &
                                    loc_array(DART_LOC_LAT), loc_array(DART_LOC_LON))
          if(obs_kind == QTY_VORTEX_LAT) then
            obs_val = loc_array(DART_LOC_LAT)
          else
            obs_val = loc_array(DART_LOC_LON)
          endif
          interp_worked = .true.

          ! deallocate arrays
          deallocate(vort_z_smooth, stat=dealloc_status)
          call check_dealloc_status(dealloc_status, routine, source, revision, &
                                    revdate, 'model_interpolate: vort_z_smooth')

          deallocate(heights, stat=dealloc_status)
          call check_dealloc_status(dealloc_status, routine, source, revision, &
                                    revdate, 'model_interpolate: heights')

          deallocate(vort_z, stat=dealloc_status)
          call check_dealloc_status(dealloc_status, routine, source, revision, &
                                    revdate, 'model_interpolate: vort_z')

          deallocate(vort_z_zlev, stat=dealloc_status)
          call check_dealloc_status(dealloc_status, routine, source, revision, &
                                    revdate, 'model_interpolate: vort_z_zlev')

          deallocate(mask, stat=dealloc_status)
          call check_dealloc_status(dealloc_status, routine, source, revision, &
                                    revdate, 'model_interpolate: mask')

        case default 

          call interpolate(x, domain, state_definition, location, &
                           obs_kind, obs_val, interp_worked)

        end select

        if (interp_worked) then
            interp_status = 0
        else
            obs_val = MISSING_R8
            interp_status = 1
        end if

        return
    end subroutine model_interpolate

    ! ens_mean_for_model
    ! ------------------
    ! Allow the ensemble mean to be passed in and stored if we need it
    ! (can be handy for forward operators)
    !  PARAMETERS
    !   IN  ens_mean          ensemble mean state vector
    subroutine ens_mean_for_model(ens_mean)
        real(kind=r8), intent(in), dimension(:) :: ens_mean

        if (need_mean) then
            ensemble_mean = ens_mean
        end if
    end subroutine ens_mean_for_model

    ! get_close_maxdist_init
    ! ----------------------
    ! Set the maximum distance for the processor for finding nearby 
    ! points. Wrapper for location module's get_close_maxdist_init 
    ! subroutine.
    !  PARAMETERS
    ! INOUT gc                get_close_type structure to initialize
    !   IN  maxdist           the maximum distance to process  
    subroutine get_close_maxdist_init (gc, maxdist, maxdist_array)
        type(get_close_type), intent(inout) :: gc
        real(r8), intent(in)                :: maxdist
        real(r8), intent(in), optional      :: maxdist_array(:)

        call loc_get_close_maxdist_init(gc, maxdist, maxdist_array)
    end subroutine get_close_maxdist_init

    ! get_close_obs_init
    ! ------------------
    ! Initializes part of get_close accelerator that depends on the
    ! particular observation(s).  Wrapper for location module's
    ! get_close_obs_init subroutine.
    !  PARAMETERS
    ! INOUT  gc               get_close_type accelerator to initialize
    !   IN   num              number of observations in the set
    !   IN   obs              set of observation locations
    subroutine get_close_obs_init(gc, num, obs)
        type(get_close_type), intent(inout) :: gc
        integer, intent(in)                 :: num
        type(location_type), intent(in)     :: obs(num)

        call loc_get_close_obs_init(gc, num, obs)
    end subroutine get_close_obs_init

    ! get_close_obs
    ! -------------
    ! Gets the number of close observations.  Wrapper for location
    ! module's get_close_obs subroutine.
    !  PARAMETERS
    !   IN  gc                get_close_type accelerator
    !   IN  base_obs_loc      location of the base observation
    !   IN  obs_loc           location of all the observations
    !   IN  base_obs_kind     raw type of the base observation
    !   IN  obs_kind          raw type of all the observations
    !   OUT num_close         how many observations are close to the
    !                         base observation
    !   OUT close_ind         which of the observations are close to
    !                         the base observation
    !   OUT dist              OPTIONAL distance from the observations
    !                         to the base observation
    subroutine get_close_obs(gc, base_obs_loc, base_obs_kind, obs_loc, &
                             obs_kind, num_close, close_ind, dist)

    type(get_close_type), intent(in)    :: gc
    type(location_type), intent(inout)  :: base_obs_loc, obs_loc(:)
    integer, intent(in)                 :: base_obs_kind, obs_kind(:)
    integer, intent(out)                :: num_close, close_ind(:)
    real(r8), optional, intent(out)     :: dist(:)

    integer                                 :: t_ind, istatus1, istatus2, k
    integer                                 :: base_which, local_obs_which
	real(r8), dimension(3)                  :: base_array, local_obs_array
	real(r8)                                :: obs_val
    type(location_type)                     :: local_obs_loc

    ! Initialize variables to missing status
    num_close = 0 ; close_ind = -99 ; dist = 1.0e9
    istatus1  = 0 ; istatus2  = 0

    base_which =nint(query_location(base_obs_loc)) 
    base_array = get_location(base_obs_loc)

    ! Consider horizontal localization first.  Also consider undefined vertical coordinate.
    if(horiz_dist_only .or. base_which == VERTISUNDEF) then
      if(present(dist)) then
        call loc_get_close_obs(gc, base_obs_loc, base_obs_kind, obs_loc, obs_kind, &
                               num_close, close_ind, dist)
      else    
         call loc_get_close_obs(gc, base_obs_loc, base_obs_kind, obs_loc, obs_kind, &
                                num_close, close_ind)
      end if
    else
      ! Convert base_obs vertical coordinate to requested vertical coordinate if necessary

      if (base_which /= VERTISLEVEL) then
        if(base_which == VERTISSURFACE) then
          base_array(3)=0.0_r8 ; base_which=VERTISLEVEL
        else
          call model_interpolate(ensemble_mean, base_obs_loc, QTY_VERTLEVEL, obs_val, istatus1)
          base_array(3)=obs_val ; base_which=VERTISLEVEL
        endif
        base_obs_loc = set_location(base_array(1),  base_array(2), base_array(3), base_which)
      elseif (base_array(3) == MISSING_R8) then
        istatus1 = 1
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
          t_ind           = close_ind(k)
          local_obs_loc   = obs_loc(t_ind)
          local_obs_which = nint(query_location(local_obs_loc))
          local_obs_array = get_location(local_obs_loc)

          ! Convert local_obs vertical coordinate to requested vertical coordinate if necessary.
          ! This should only be necessary for obs priors, as state location information already
          ! contains the correct vertical coordinate (filter_assim's call to get_state_meta_data).
          if (local_obs_which /= VERTISLEVEL) then
            if(local_obs_which == VERTISSURFACE) then
              local_obs_array(3)=obs_val ; local_obs_which=VERTISLEVEL
            elseif(local_obs_which == VERTISUNDEF) then
              local_obs_which=VERTISUNDEF
            else
              call model_interpolate(ensemble_mean, obs_loc(t_ind), QTY_VERTLEVEL, obs_val, istatus2)
              local_obs_array(3)=obs_val ; local_obs_which=VERTISLEVEL
            end if

            ! Store the "new" location into the original full local array
            local_obs_loc = set_location(local_obs_array(1),local_obs_array(2), &
                                         local_obs_array(3),local_obs_which)
            obs_loc(t_ind) = local_obs_loc
          endif

          ! Compute distance - set distance to a very large value if vert coordinate is missing
          ! or vert_interpolate returned error (istatus2=1)
          if ((local_obs_array(3) == missing_r8) .or. (istatus2 == 1)) then
            dist(k) = 1.0e9
          else
            dist(k) = get_dist(base_obs_loc, local_obs_loc, base_obs_kind, obs_kind(t_ind))
          endif

        end do
      end if
    end if

    end subroutine get_close_obs

    ! get_model_size
    ! --------------
    ! Returns the size of the DART state vector
    !  PARAMETERS
    !   OUT get_model_size    length of the DART state vector
    function get_model_size()
        integer :: get_model_size

        get_model_size = get_total_size(state_definition)
    end function get_model_size

    ! get_state_meta_data
    ! -------------------
    ! Get the location and variable type for a particular index in
    ! the state vector
    !  PARAMETERS
    !   IN  index_in          position in state vector to query
    !   OUT location          DART location_type for that index
    !   OUT var_type          OPTIONAL numeric variable type
    subroutine get_state_meta_data(index_in, location, var_type)
        integer,                       intent(in)  :: index_in
        type(location_type), optional, intent(out) :: location
        integer,             optional, intent(out) :: var_type

        type(state_variable)           :: cur_var
        type(coamps_nest)              :: cur_nest
        type(nest_point)               :: cur_point
        integer, dimension(3)          :: ijk
        integer                        :: index_loc
        real(kind=r8)                  :: lat, lon

        cur_var = get_var_by_index(state_layout_3D, all_vars(index_in))

        if (present(location)) then 
            cur_nest = get_domain_nest(domain, get_nest_number(cur_var))

            index_loc = index_in - get_state_begin(cur_var) + 1
            ijk = nest_index_1d_to_3d(cur_nest, index_loc)

            call get_nest_latlon(cur_nest, ijk(1), ijk(2), lat, lon)        

            location = set_location(lon, lat, &
                                    get_vert_loc(cur_var, domain, ijk(3)), &
                                    get_vert_type(cur_var))
        end if

        if (present(var_type)) var_type = get_var_kind(cur_var)

    end subroutine get_state_meta_data

    ! get_model_time_step
    ! -------------------
    ! Returns the smallest increment in time that the model is capable
    ! of advancing the state - just call it a minute for now
    !  PARAMETERS
    !   OUT get_model_time_step  model time step as a DART time_type
    function get_model_time_step()
        type(time_type) :: get_model_time_step

        get_model_time_step = set_time(60,0)
    end function get_model_time_step

    ! init_conditions
    ! ---------------
    ! NULL INTERFACE
    ! (sets up initial conditions for the model, but we're using
    ! already existing restart file data)
    !  PARAMETERS
    !   OUT x                 state vector initial condition 
    subroutine init_conditions(x)
        real(r8), intent(out) :: x(:)

        call error_handler(E_ERR, 'init_conditions', 'inoperable', &
                           source, revision, revdate) 
        x = MISSING_R8

    end subroutine init_conditions

    ! init_time
    ! ---------
    ! NULL INTERFACE
    ! (sets up initial time information for the model, but we're using
    ! already existing restart file data)
    !  PARAMETERS
    !   OUT time              time initial condition
    subroutine init_time(time)
        type(time_type), intent(out) :: time

      time = set_time_missing()

    end subroutine init_time

    ! adv_1step
    ! ---------
    ! NULL INTERFACE
    ! (advances the model with function calls, but we're doing it
    ! asynchronously)
    !  PARAMETERS
    ! INOUT x                 (in)  state vector analysis
    !                         (out) state vector forecast 
    ! INOUT time              (in)  analysis time
    !                         (out) forecast time 
    subroutine adv_1step(x, time)
        real(r8), intent(inout) :: x(:)
        type(time_type), intent(in) :: time
        x = MISSING_R8  ! Just to satisfy compiler/complainer

    end subroutine adv_1step

    ! get_coamps_domain
    ! ---------
    ! Returns the module defined coamps_domain
    function get_coamps_domain()
      type(coamps_domain) :: get_coamps_domain
      get_coamps_domain = domain
    end function get_coamps_domain

    !------------------------------
    ! END PUBLIC ROUTINES
    !------------------------------

    !------------------------------
    ! BEGIN PRIVATE ROUTINES
    !------------------------------

    ! read_model_namelist
    ! -------------------
    ! Read in parameters from the model_nml namelist
    subroutine read_model_namelist

        character(len=*), parameter :: NAMELIST_FILE  = 'input.nml'
        character(len=*), parameter :: MODEL_NAMELIST = 'model_nml'

        integer :: io_status
        integer :: nml_unit

        call find_namelist_in_file(NAMELIST_FILE, MODEL_NAMELIST, nml_unit)
        read (nml_unit,nml=model_nml,iostat=io_status)
        call check_namelist_read(nml_unit, io_status, MODEL_NAMELIST)
        if (do_output()) write (*, model_nml)
    end subroutine read_model_namelist

    ! allocate_metadata_arrays
    ! ------------------------
    ! Initialize storage for state vector locations, types, and possibly the
    ! ensemble mean
    subroutine allocate_metadata_arrays

        character(len=*), parameter :: routine = 'allocate_metadata_arrays'
        integer                     :: alloc_status

        integer                     :: nvars

        allocate( all_vars(get_model_size()), stat=alloc_status )
        call check_alloc_status(alloc_status, routine, source, revision, &
                                revdate, 'all_vars')

        if (need_mean) then
            allocate( ensemble_mean(get_model_size()), stat=alloc_status )
            call check_alloc_status(alloc_status, routine, source, revision, &
            revdate, 'ensemble_mean')
        end if

    end subroutine allocate_metadata_arrays

    ! deallocate_metadata_arrays
    ! --------------------------
    ! Finalize storage for state vector lcoations, types and possibly the
    ! ensemble mean
    subroutine deallocate_metadata_arrays

        character(len=*), parameter :: routine = 'deallocate_metadata_arrays'
        integer                     :: dealloc_status

        deallocate(all_vars, stat=dealloc_status)
        call check_dealloc_status(dealloc_status, routine, source, revision, &
                                  revdate, 'all_vars')

        if (need_mean) then
            deallocate(ensemble_mean, stat=dealloc_status)
            call check_dealloc_status(dealloc_status, routine, source,       &
                                      revision, revdate, 'all_locs')
        end if
    end subroutine deallocate_metadata_arrays

    ! populate_metadata_arrays
    ! ------------------------
    ! Populates the various arrays of metadata for the model state
    subroutine populate_metadata_arrays()

        type(state_variable), pointer :: cur_var
        type(state_iterator)          :: iterator
        integer :: state_ii

        state_ii = 0
        iterator  = get_iterator(state_layout_3D)
        do while (has_next(iterator))
            state_ii = state_ii + 1
            cur_var   =>  get_next(iterator)
            all_vars(get_state_begin(cur_var):get_state_end(cur_var)) = state_ii
        end do

    end subroutine populate_metadata_arrays

    ! index_to_var
    ! ------------------------
    ! returns the state_variable at a given index point
    function index_to_var(index_in)
      integer,            intent(in) :: index_in
      type(state_variable)           :: index_to_var

      type(state_iterator) :: var_iterator
      integer              :: ii

      var_iterator  = get_iterator(state_layout_3D)
      do while (has_next(var_iterator))
         index_to_var = get_next(var_iterator)
         if(is_var_at_index(index_to_var, index_in)) exit 
      end do

    end function index_to_var

    subroutine kernal_smoother(fin, fout, half_width, nx, ny)
      real(kind=r8), intent(in)  :: fin(:, :)
      real(kind=r8), intent(out) :: fout(:, :)
      integer,       intent(in)  :: half_width
      integer,       intent(in)  :: nx
      integer,       intent(in)  :: ny

      real(kind=r8), allocatable :: kernal(:,:)
      real(kind=r8) :: LL, kernal_sum
      integer :: i, j, ii, jj, n_width

        character(len=*), parameter :: routine = 'kernal_smoother'
        integer                     :: alloc_status, dealloc_status

      n_width = 2*half_width + 1

      allocate(kernal(n_width, n_width), stat=alloc_status )
      call check_alloc_status(alloc_status, routine, source, revision, revdate, 'kernal')

      ! Setup kernal
      kernal_sum = 0.0_r8 
      do i=1,n_width ; do j=1,n_width
        LL = sqrt(real((half_width + 1 - j)**2 + (half_width + 1 - i)**2, kind=r8))  
        if(LL <= real(half_width,kind=r8)) then
          kernal(i, j) = 1.0_r8 - LL/real(half_width,kind=r8)
        else
          kernal(i, j) = 0.0_r8
        end if
        kernal_sum = kernal(i, j) + kernal_sum
      end do ; end do

      kernal(:,:) = kernal(:,:) / kernal_sum

      fout(:, :) = 0.0_r8
      do i=1+half_width,nx-half_width ; do j=1+half_width,ny-half_width
        do ii=-half_width,half_width ; do jj=-half_width,half_width
          fout(i, j) = fout(i, j) + fin(i+ii, j+jj)*kernal(ii+half_width+1, jj+half_width+1)
        end do ; end do
      end do ; end do

      deallocate(kernal, stat=dealloc_status )
      call check_dealloc_status(dealloc_status, routine, source, revision, revdate, 'kernal')

    end subroutine kernal_smoother
     
    !------------------------------
    ! END PRIVATE ROUTINES
    !------------------------------

end module model_mod

! <next few lines under version control, do not edit>
! $URL$
! $Id$
! $Revision$
! $Date$
