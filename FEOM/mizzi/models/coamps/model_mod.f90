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

  use coamps_grid_mod,     only : coamps_grid,                  &
                                  dump_grid_info,               &
                                  get_grid_dims,                &
                                  get_grid_dsigmaw,             &
                                  get_grid_field_size,          &
                                  get_grid_msigma,              &
                                  get_grid_num_levels,          &
                                  get_grid_wsigma,              &
                                  get_terrain_height_at_points, &
                                  gridpt_to_latlon,             &
                                  location_to_gridpt
  use coamps_interp_mod,   only : interpolate
  use coamps_restart_mod,  only : dump_restart_vars,            &
                                  get_num_vars,                 &
                                  get_pert_magnitude_by_index,  &
                                  get_pert_type_by_index,       &
                                  get_restart_grid,             &
                                  get_var_type_by_index,        &
                                  get_vert_coord_by_index,      &
                                  initialize_restart_info,      &
                                  PERT_TYPE_INDIVID,            &  
                                  PERT_TYPE_NOPERTS,            &
                                  PERT_TYPE_UNIFORM
  use coamps_util_mod,     only : check_alloc_status,           &
                                  check_dealloc_status

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
                                  vert_is_pressure
  use netcdf
  use obs_kind_mod
  use random_seq_mod,      only : init_random_seq,               &
                                  random_gaussian,               &
                                  random_seq_type
  use time_manager_mod,    only : set_time,                      &
                                  time_type
  use typeSizes
  use types_mod,           only : MISSING_R8,                    &
                                  r8
  use utilities_mod,       only : check_namelist_read,           &
                                  do_output,                     &
                                  E_ERR,                         &
                                  E_MSG,                         &
                                  error_handler,                 &
                                  find_namelist_in_file,         &
                                  get_unit,                      &
                                  nc_check,                      &
                                  register_module,               &
                                  do_nml_file,                   &
                                  do_nml_term,                   &
                                  nmlfileunit

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
  character(len=10) :: cdtg = '1999083100' ! Date-time group
  integer           :: y_bound_skip = 3    ! How many x and y boundary
  integer           :: x_bound_skip = 3    ! points to skip when
                                           ! perturbing the model
                                           ! state
  logical           :: need_mean = .false. ! Do we need the ensemble
                                           ! mean for for forward
                                           ! operator computation?
  namelist /model_nml/ cdtg, y_bound_skip, x_bound_skip, need_mean

  ! Locations of state variables
  type(location_type), dimension(:), allocatable :: all_locs
  integer, dimension(:), allocatable             :: all_types

  ! Grid information structure
  type(coamps_grid) :: restart_grid

  ! Ensemble mean
  real(r8), dimension(:), allocatable :: ensemble_mean

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

    integer :: nml_unit
    
    character(len=*), parameter :: routine = 'static_init_model'
    integer :: io_status, alloc_status

    integer :: num_vars, max_i, max_j

    integer :: state_ii, ii,jj, var_index

    real(kind=r8) :: lon, lat
    real(kind=r8) :: vert_loc
    integer       :: vert_type


    call register_module(source, revision, revdate)

    ! Read in the namelist information
    call find_namelist_in_file('input.nml', 'model_nml', nml_unit)
    read (nml_unit,nml=model_nml,iostat=io_status)
    call check_namelist_read(nml_unit, io_status, 'model_nml')
    if (do_nml_file()) write (nmlfileunit,model_nml)
    if (do_nml_term()) write (     *     ,model_nml)

    call initialize_restart_info(cdtg, 'restart.vars')

    call get_restart_grid(restart_grid)
    if (do_output()) call dump_grid_info(restart_grid)
    
    allocate( all_locs(get_model_size()), stat=alloc_status )
    call check_alloc_status(alloc_status, routine, source, revision, &
                            revdate, 'all_locs')
    allocate( all_types(get_model_size()), stat=alloc_status )
    call check_alloc_status(alloc_status, routine, source, revision, &
                            revdate, 'all_types')
    if (need_mean) then
       allocate( ensemble_mean(get_model_size()), stat=alloc_status )
       call check_alloc_status(alloc_status, routine, source, revision, &
                               revdate, 'ensemble_mean')
    end if

    ! Calculate loop limits
    call get_num_vars(num_vars)
    call get_grid_dims(restart_grid, max_i, max_j)

    ! Location is easy provided the i/j coordinate and
    ! the sigma level, and the type is given by the restart file
    ! information structure
    state_ii = 1
    do var_index = 1, num_vars
       do jj = 1, max_j
          do ii =1, max_i
             call gridpt_to_latlon(grid=restart_grid,    &
                                   ii=real(ii, kind=r8), &
                                   jj=real(jj, kind=r8), &
                                   lat=lat, lon=lon)

             ! Store the variable type for easy access from
             ! get_state_meta_data.
             call get_var_type_by_index(var_index,                 &
                                        all_types(state_ii))

             ! Vertical coordinate information is based on the 
             ! state vector definition since it's a field-wide
             ! constant
             call get_vert_coord_by_index(var_index, vert_type,    &
                                          vert_loc)

             all_locs(state_ii) = set_location(lon, lat, vert_loc, &
                                               vert_type)

             state_ii = state_ii + 1
          enddo
       enddo
    enddo
  end subroutine static_init_model

  ! end_model
  ! ---------
  ! Clean up the workspace once the program is ending
  !  PARAMETERS
  !   [none]
  subroutine end_model()

    character(len=*), parameter :: routine = 'end_model'
    integer                     :: dealloc_status

    deallocate(all_locs, stat=dealloc_status)
    call check_dealloc_status(dealloc_status, routine, source, revision, &
                              revdate, 'all_locs')
    deallocate(all_types, stat=dealloc_status)
    call check_dealloc_status(dealloc_status, routine, source, revision, &
                              revdate, 'all_locs')
    if (need_mean) then
       deallocate(ensemble_mean, stat=dealloc_status)
        call check_dealloc_status(dealloc_status, routine, source,       &
                                  revision, revdate, 'all_locs')
    end if
  end subroutine end_model

  ! nc_write_model_atts
  ! -------------------
  ! Write model-specific global attributes to a NetCDF file.
  !  PARAMETERS
  !   IN  ncFileID          numeric ID of an *open* NetCDF file
  !   OUT ierr              0 if writing was successful
  function nc_write_model_atts( ncFileID ) result (ierr)
    integer, intent(in)         :: ncFileID      ! netCDF file 
    integer                     :: ierr          

    integer                     :: model_size 
    integer                     :: grid_i_size
    integer                     :: grid_j_size
    type(location_type)         :: loc        
    real(kind=r8), dimension(3) :: loc3d      ! lon/lat/hgt
    integer                     :: vartype    ! state variable type
    integer                     :: ii         
 
    ! NetCDF variables                                               
    integer :: n_dims, n_vars, n_atts  ! NetCDF counts               
    integer :: state_varid             ! State vars     
    integer :: state_coord_varid       ! State vector coordinate var 
    integer :: type_varid              ! State var type variable
    integer :: time_dimid, state_dimid ! Time and state dimensions   
    integer :: mem_dimid               ! Ensemble member dimension   
    integer :: ulim_dimid              ! unlimited dimension         
    integer :: lat_dimid,lat_varid     ! Latitude coordinate
    integer :: lon_dimid,lon_varid     ! Longitude coordinate
    integer :: vert_dimid,vert_varid   ! Vertical coordinate

    ! Date and time
    integer, dimension(8)        :: dt_values 
    character(len=NF90_MAX_NAME) :: dt_string 

    ! Error handling
    character(len=*), parameter :: routine = 'nc_write_model_atts'

    call get_grid_dims(restart_grid, grid_i_size, grid_j_size)

    ! Default to no errors                                           
    ierr = 0 

    ! Ensure that we're dealing with an open & current  NetCDF file  
    call nc_check(nf90_inquire(ncFileID, n_dims, n_vars, n_atts, &
                               ulim_dimid),                      &
                  routine)  
    call nc_check(nf90_sync(ncFileID), routine) 
 
    ! Go into define mode                                            
    call nc_check(nf90_redef(ncFileID), routine) 
 
    ! Find the ensemble member and time dimensions                    
    call nc_check(nf90_inq_dimid(ncFileID,'copy',dimid=mem_dimid), &
                  routine) 
    call nc_check(nf90_inq_dimid(ncFileID,'time',dimid=time_dimid),&
                  routine) 
 
    ! Get the 'creation date' using the intrinsic F90 function       
    call DATE_AND_TIME(values=dt_values) 
    write (dt_string, '(''YYYY MM DD HH MM SS = '',i4,5(1x,i2.2))')& 
         dt_values(1), dt_values(2), dt_values(3), dt_values(5),   &
         dt_values(6), dt_values(7) 
 
    ! Global attributes                                   
    call nc_check(nf90_put_att(ncFileID, NF90_GLOBAL,           &
                               'creation_date', dt_string),     &
                  routine) 
    call nc_check(nf90_put_att(ncFileID, NF90_GLOBAL,           &
                               'model_source', source),         &
                  routine) 
    call nc_check(nf90_put_att(ncFileID, NF90_GLOBAL,           &
                               'model_revision', revision),     & 
                  routine)  
    call nc_check(nf90_put_att(ncFileID, NF90_GLOBAL,           &
                               'model_revdate', revdate),       &
                  routine) 
    call nc_check(nf90_put_att(ncFileID, NF90_GLOBAL, 'model',  & 
                               'COAMPS'),                       &
                  routine)
    call nc_check(nf90_put_att(ncFileID, NF90_GLOBAL, 'grid_x', &
                               grid_i_size),                    &
                  routine)
    call nc_check(nf90_put_att(ncFileID, NF90_GLOBAL, 'grid_y', &
                               grid_j_size),                    &
                  routine)

    ! State dimension                                     
    model_size = get_model_size() 
    call nc_check(nf90_def_dim(ncid=ncFileID, name='StateVariable',&
                               len=model_size, dimid=state_dimid), &
                  routine)  
 
    ! Location dimensions
    call nc_check(nf90_def_dim(ncid=ncFileID, name='lat',         &
                               len=model_size, dimid=lat_dimid),  &
                  routine)
    call nc_check(nf90_def_dim(ncid=ncFileID, name='lon',         &
                               len=model_size, dimid=lon_dimid),  &
                  routine)
    call nc_check(nf90_def_dim(ncid=ncFileID, name='height',      &
                               len=model_size, dimid=vert_dimid), &
                  routine)

    ! Location coordinates - latitude
    call nc_check(nf90_def_var(ncid=ncFileID, name='lat',           &
                               xtype=nf90_double, dimids=lat_dimid, &
                               varid=lat_varid),                    &
                  routine) 
    call nc_check(nf90_put_att(ncFileID, lat_varid, 'long_name',    &
                               'Latitude'),                         &
                  routine)
    call nc_check(nf90_put_att(ncFileID, lat_varid, 'units',        &
                               'degrees_north'),                    &
                  routine)
    call nc_check(nf90_put_att(ncFileID, lat_varid, 'valid_range',  &
                               (/ -90.0_r8, 90.0_r8 /)),            &
                  routine)

    ! Location coordinates - longitude
    call nc_check(nf90_def_var(ncid=ncFileID, name='lon',           &
                               xtype=nf90_double, dimids=lon_dimid, &
                               varid=lon_varid),                    &
                  routine) 
    call nc_check(nf90_put_att(ncFileID, lon_varid, 'long_name',    &
                               'Longitude'),                        &
                  routine)
    call nc_check(nf90_put_att(ncFileID, lon_varid, 'units',        &
                               'degrees_east'),                     &
                  routine) 
    call nc_check( nf90_put_att(ncFileID, lon_varid, 'valid_range', &
                                (/ 0.0_r8, 360.0_r8 /)),            &
                   routine)

    ! Location coordinates - height
    call nc_check(nf90_def_var(ncid=ncFileID, name='height',        &
                               xtype=nf90_double, dimids=vert_dimid,&
                               varid=vert_varid),                   &
                  routine) 
    call nc_check(nf90_put_att(ncFileID, vert_varid, 'long_name',   &
                               'Sigma Height'),                     &
                  routine)
    call nc_check(nf90_put_att(ncFileID, vert_varid, 'units',       &
                               'meters'),                           &
                  routine)
 
    ! Location coordinates - state    
    call nc_check(nf90_def_var(ncid=ncFileID, name='StateVariable', & 
                               xtype=nf90_int, dimids=state_dimid,  &
                                varid=state_coord_varid),           &
                  routine) 
    call nc_check(nf90_put_att(ncFileID, state_coord_varid,         &
                               'long_name',                         &
                               ' State Variable Coordinate'),       &
                  routine) 
    call nc_check(nf90_put_att(ncFileID, state_coord_varid, 'units',& 
                               'index'),                            &
                  routine) 
    call nc_check(nf90_put_att(ncFileID, state_coord_varid,         & 
                               'valid_range', (/1, model_size /)),  &
                  routine) 
    
    ! State variable              
    call nc_check(nf90_def_var(ncid=ncFileID, name='state',         &
                               xtype=nf90_double,                   &
                               dimids=(/ state_dimid, mem_dimid,    &
                               time_dimid /), varid=state_varid),   &
                  routine)  
    call nc_check(nf90_put_att(ncFileID, state_varid, 'long_name',  & 
                               'model state'),                      &
                  routine) 

    ! Variable types
    call nc_check(nf90_def_var(ncFileID, name='type',               &
                               xtype=nf90_int, dimids=state_dimid,  &
                               varid=type_varid),                   &
                  routine)
    call nc_check(nf90_put_att(ncFileID, type_varid, 'long_name',   &
                               'Variable Type'),                    &
                  routine)

    ! Need to get out of define mode to fill the variables 
    call nc_check(nf90_enddef(ncFileID), routine) 
 
    ! Fill the state variable coordinate                             
    call nc_check(nf90_put_var(ncFileID, state_coord_varid,          &
                               (/ (ii,ii=1, model_size) /)),         &
                  routine)  
  
    ! Fill the location variables
    do ii=1, model_size
       call get_state_meta_data(ii,loc,vartype)
       loc3d = get_location(loc)
       call nc_check(nf90_put_var(ncFileID, lon_varid, loc3d(1),  &
                                  (/ ii /)),                      &
                     routine)
       call nc_check(nf90_put_var(ncFileID, lat_varid, loc3d(2),  &
                                  (/ ii /)),                      &
                     routine)
       call nc_check(nf90_put_var(ncFileID, vert_varid, loc3d(3), &
                                  (/ ii /)),                      &
                     routine)
       call nc_check(nf90_put_var(ncFileID, type_varid, vartype,  &
                                  (/ ii /)),                      &
                     routine)
    end do 
 
    ! Sync up the buffer with the NetCDF file on disk, but don't     
    ! close it                                                       
    call nc_check(nf90_sync(ncFileID), routine) 
    if (do_output()) then
       write (*,*) 'Model attributes written and file synced...' 
    end if
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
  function nc_write_model_vars(ncFileID, statevec, copyindex, &
                               timeindex ) result (ierr)          
    integer,                intent(in) :: ncFileID
    real(r8), dimension(:), intent(in) :: statevec
    integer,                intent(in) :: copyindex
    integer,                intent(in) :: timeindex
    integer                            :: ierr   

    ! NetCDF variables
    integer :: n_dims, n_vars, n_atts ! NetCDF counts
    integer :: state_varid            ! state variable ID
    integer :: ulim_dimid             ! unlimited dimension (time)

    ! Error handling
    character(len=*), parameter :: routine = 'nc_write_model_atts'

    ! Default to no errors - since any errors encountered here will
    ! just throw us out of the program
    ierr = 0

    ! Ensure that we're dealing with a valid NetCDF file
    call nc_check(nf90_inquire(ncFileID, n_dims, n_vars, n_atts,  &
                               ulim_dimid),                       & 
                  routine)

    ! Get the numerical index of the state variable, then dump the
    ! state vector to the corresponding ensemble member and time
    call nc_check(nf90_inq_varid(ncFileID, 'state', state_varid), &
                  routine) 
    call nc_check(nf90_put_var(ncFileID, state_varid, statevec,   &
                               start=(/ 1,copyindex,timeindex /)),&
                  routine)
 
    ! Flush the buffer to disk                                       
    call nc_check( nf90_sync(ncFileId), routine) 
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
    ! Avoid having to link *everything* with the MPI utilities
    use mpi_utilities_mod, only : my_task_id

    real(r8), intent(in)  :: state(:)
    real(r8), intent(out) :: pert_state(:)
    logical,  intent(out) :: interf_provided

    type(random_seq_type), save :: random_sequence
    logical,               save :: initialized_random_sequence = .false. 
    integer,               save :: rand_seq_seed

    integer, parameter          :: DT_MILLISECOND = 8
    integer, dimension(8), save :: datetime

    real(r8)              :: pert, scale
    integer               :: ii,jj,kk
    integer               :: ijarea,gridii,gridjj
    integer               :: fullii

    integer :: num_vars

    character(len=*), parameter :: routine = 'pert_model_state'
    integer                     :: alloc_status, dealloc_status
    
    ! Base the perturbation on the restart.vars file
    real(kind=r8) :: pertsize
    integer       :: perttype

    ! Do the perturbations based on a 2-D grid: this'll make it
    ! easier to skip boundaries
    real(r8), dimension(:,:), allocatable :: temp_state
    real(r8), dimension(:,:), allocatable :: temp_pert

    if (.not. initialized_random_sequence) then

        ! Seed the random number generator based on the number
        ! SSmmmP, where:
        !  SS  = current time's seconds value
        !  mmm = current time's milliseconds value
        !  P   = current process ID
        ! This should give us unique perturbations between different
        ! processors and different runs of the software.
        call DATE_AND_TIME(values=datetime)
        rand_seq_seed = 1E1 * datetime(DT_MILLISECOND) + & 
                        1E0 * my_task_id()
        call init_random_seq(r=random_sequence, seed=rand_seq_seed)
    end if

    call error_handler(E_MSG, routine, 'Perturbing model state', &
                       source, revision, revdate) 

    call get_grid_field_size(restart_grid, ijarea)
    call get_grid_dims(restart_grid, gridii, gridjj)
    call get_num_vars(num_vars)

    allocate( temp_state(gridii,gridjj), stat=alloc_status )
    call check_alloc_status(alloc_status, routine, source, revision, &
                            revdate, 'temp_state')
    allocate( temp_pert(gridii,gridjj) , stat=alloc_status)
    call check_alloc_status(alloc_status, routine, source, revision, &
                            revdate, 'temp_pert')

    ! Each field in the restart vector may be treated differently
    do kk=1,num_vars
       call get_pert_type_by_index(kk,perttype)
       call get_pert_magnitude_by_index(kk, pertsize)
 
       ! Pull out just this field from the state vector      
       fullii = 1 + (kk - 1) * ijarea
       temp_state = reshape(state(fullii:(kk * ijarea)), &
                           (/ gridii, gridjj /))

       ! Apply the perturbation - may not perturb the boundaries
       if (perttype .eq. PERT_TYPE_NOPERTS) then

          ! Perturbed state is the same as the state
          do jj=1,gridjj
             do ii=1,gridii
                temp_pert(ii,jj) = temp_state(ii,jj)
             end do
          end do

       elseif (perttype .eq. PERT_TYPE_UNIFORM) then

          ! Single perturbation value for this field
          pert = random_gaussian(random_sequence, real(0.0,kind=r8), &
                                 pertsize)

          do jj=(1 + y_bound_skip), (gridjj - y_bound_skip)
             do ii=(1 + x_bound_skip),(gridii - x_bound_skip)
                temp_pert(ii,jj) = temp_state(ii,jj) + pert
             enddo
          end do

       elseif (perttype .eq. PERT_TYPE_INDIVID) then

          ! Each point in the field gets its own perturbation
          do jj=(1 + y_bound_skip), (gridjj - y_bound_skip)
             do ii=(1 + x_bound_skip), (gridii - x_bound_skip)
                pert = random_gaussian(random_sequence, real(0.0,kind=r8),&
                                       pertsize)
                temp_pert(ii,jj) = temp_state(ii,jj) + pert
             enddo
          end do

       end if

       pert_state(fullii:(kk*ijarea)) = reshape(temp_pert,  &
                                                (/ ijarea /))
     enddo

    interf_provided = .true.

    deallocate(temp_state, stat=dealloc_status)
    call check_dealloc_status(dealloc_status, routine, source, revision, &
                              revdate, 'temp_state')
    deallocate(temp_pert,  stat=dealloc_status)
    call check_dealloc_status(dealloc_status, routine, source, revision, &
                              revdate, 'temp_pert')
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
  !   IN  obs_type          raw variable type
  !   OUT obs_val           interpolated value
  !   OUT interp_status     status of interpolation (0 is success)
  subroutine model_interpolate(x, location, obs_type, obs_val, &
                               interp_status)
    real(r8), dimension(:), intent(in) :: x
    type(location_type),    intent(in) :: location
    integer,                intent(in) :: obs_type
    real(r8),              intent(out) :: obs_val
    integer,               intent(out) :: interp_status

    logical :: successful_interpolation

    ! Just hand this off to the COAMPS interpolation module - it will
    ! handle getting all the proper variables out - note that we are
    ! currently not supporting sea level pressure
    call interpolate(x, restart_grid, location, obs_type, obs_val, &
                     successful_interpolation)
    if (successful_interpolation) then
      interp_status = 0
    else
      obs_val = MISSING_R8
      interp_status = 1
    end if
  end subroutine model_interpolate

  ! ens_mean_for_model
  ! ------------------
  ! Allow the ensemble mean to be passed in and stored if we need it
  ! (can be handy for forward operators)
  !  PARAMETERS
  !   IN  ens_mean          ensemble mean state vector
  subroutine ens_mean_for_model(ens_mean)
    real(r8), intent(in), dimension(:) :: ens_mean

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
  subroutine get_close_maxdist_init (gc, maxdist)
    type(get_close_type), intent(inout) :: gc
    real(r8), intent(in)                :: maxdist

    call loc_get_close_maxdist_init(gc, maxdist)
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
  !   IN  obs               location of all the observations
  !   IN  base_obs_kind     raw type of the base observation
  !   IN  obs_kind          raw type of all the observations
  !   OUT num_close         how many observations are close to the
  !                         base observation
  !   OUT close_ind         which of the observations are close to
  !                         the base observation
  !   OUT dist              OPTIONAL distance from the observations
  !                         to the base observation
  subroutine get_close_obs(gc, base_obs_loc, base_obs_kind, obs,&
       & obs_kind, num_close, close_ind, dist)
    type(get_close_type), intent(in) :: gc
    type(location_type), intent(in)  :: base_obs_loc, obs(:)
    integer, intent(in)              :: base_obs_kind, obs_kind(:)
    integer, intent(out)             :: num_close, close_ind(:)
    real(r8), optional, intent(out)  :: dist(:)

    if (present(dist)) then
       call loc_get_close_obs(gc, base_obs_loc, base_obs_kind, obs,&
            & obs_kind, num_close, close_ind, dist) 
    else
       call loc_get_close_obs(gc, base_obs_loc, base_obs_kind, obs,&
            & obs_kind, num_close, close_ind)
    end if
  end subroutine get_close_obs

  ! get_model_size
  ! --------------
  ! Returns the size of the DART state vector
  !  PARAMETERS
  !   OUT get_model_size    length of the DART state vector
  function get_model_size()
    integer :: get_model_size

    integer :: level_gridpoints, total_levels

    call get_grid_field_size(restart_grid, level_gridpoints)
    call get_num_vars(total_levels)
    get_model_size = level_gridpoints * (total_levels)
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
    integer,             intent(in)            :: index_in
    type(location_type), intent(out)           :: location
    integer,             intent(out), optional :: var_type

    ! All this has been pre-calculated
    location = all_locs(index_in)
    if (present(var_type)) then
       var_type = all_types(index_in)
    end if
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
  end subroutine adv_1step

  !------------------------------
  ! END PUBLIC ROUTINES
  !------------------------------

  !------------------------------
  ! BEGIN PRIVATE ROUTINES
  !------------------------------
  ! [none]
  !------------------------------
  ! END PRIVATE ROUTINES
  !------------------------------

end module model_mod

! <next few lines under version control, do not edit>
! $URL$
! $Id$
! $Revision$
! $Date$
