module model_mod

! This module provides routines to work with COSMO data
! files in the DART framework
!
! Author: Jan D. Keller
!         Meteorological Institute, University of Bonn, Germany
!         2011-09-15
!
! <next few lines under version control, do not edit>
! $URL$
! $Id$
! $Revision$
! $Date$

  use        types_mod, only : r4, r8, digits12, SECPERDAY, MISSING_R8,          &
                               rad2deg, deg2rad, PI
  
  use time_manager_mod, only : time_type, set_time, set_date, get_date, get_time,&
                               print_time, print_date, set_calendar_type,        &
                               operator(*),  operator(+), operator(-),           &
                               operator(>),  operator(<), operator(/),           &
                               operator(/=), operator(<=)

  use   cosmo_data_mod, only : cosmo_meta,cosmo_hcoord,cosmo_non_state_data,     &
                               get_cosmo_info,get_data_from_binary,              &
                               set_vertical_coords,grib_header_type,             &
                               model_dims,record_length

  use     location_mod, only : location_type, get_dist, query_location,          &
                               get_close_maxdist_init, get_close_type,           &
                               set_location, get_location, horiz_dist_only,      & 
                               vert_is_undef,        VERTISUNDEF,                &
                               vert_is_surface,      VERTISSURFACE,              &
                               vert_is_level,        VERTISLEVEL,                &
                               vert_is_pressure,     VERTISPRESSURE,             &
                               vert_is_height,       VERTISHEIGHT,               &
                               vert_is_scale_height, VERTISSCALEHEIGHT,          &
                               get_close_obs_init, get_close_obs

  use    utilities_mod, only : register_module, error_handler,                   &
                               E_ERR, E_WARN, E_MSG, logfileunit, get_unit,      &
                               nc_check, do_output, to_upper,                    &
                               find_namelist_in_file, check_namelist_read,       &
                               open_file, close_file, file_exist,                &
                               find_textfile_dims, file_to_text

  use     obs_kind_mod, only : KIND_U_WIND_COMPONENT,                            &
                               KIND_V_WIND_COMPONENT,                            &
                               KIND_VERTICAL_VELOCITY,                           &
                               KIND_TEMPERATURE,                                 &
                               KIND_PRESSURE,                                    &
                               KIND_PRESSURE_PERTURBATION,                       &
                               KIND_SPECIFIC_HUMIDITY,                           &
                               KIND_CLOUD_LIQUID_WATER,                          &
                               KIND_CLOUD_ICE,                                   &
                               KIND_SURFACE_ELEVATION,                           &
                               KIND_SURFACE_GEOPOTENTIAL

  use    random_seq_mod, only: random_seq_type, init_random_seq, random_gaussian

  use byte_mod, only: to_float1,from_float1,word_to_byte,byte_to_word_signed,concat_bytes1

  use netcdf 

  implicit none

! version controlled file description for error handling, do not edit
character(len=128), parameter :: &
   source   = "$URL$", &
   revision = "$Revision$", &
   revdate  = "$Date$"

  public  :: get_model_size
  public  :: static_init_model
  public  :: get_state_meta_data
  public  :: get_model_time_step
  public  :: model_interpolate
  public  :: init_conditions
  public  :: init_time
  public  :: adv_1step
  public  :: end_model
  public  :: nc_write_model_atts
  public  :: nc_write_model_vars
  public  :: pert_model_state
  public  :: get_close_maxdist_init
  public  :: get_close_obs_init
  public  :: get_close_obs
  public  :: ens_mean_for_model

!  public  :: grib_to_sv
!  public  :: sv_to_grib

  private :: set_allowed_state_vector_vars
  private :: ll_to_xyz_vector
  private :: ll_to_xyz_single
  private :: get_enclosing_grid_box
  private :: bilinear_interpolation
  private :: linear_interpolation
!  public :: linear_interpolation
!  private :: data_to_state_vector
  private :: get_vertical_boundaries

  public  :: get_state_time
  public  :: get_state_vector
  public  :: write_grib_file
  public  :: get_cosmo_filename
  public  :: write_state_times

  INTERFACE sv_to_field
    MODULE PROCEDURE sv_to_field_2d
    MODULE PROCEDURE sv_to_field_3d
  END INTERFACE

  type dart_variable_info
    character(len=16)    :: varname_short
    character(len=256)   :: varname_long
    character(len=32)    :: units
    logical              :: is_present
    integer              :: nx
    integer              :: ny
    integer              :: nz
    real(r8),allocatable :: vertical_level(:)
    integer              :: vertical_coordinate
    integer              :: horizontal_coordinate
    integer,allocatable  :: state_vector_sindex(:) ! starting index in state vector for every vertical level
    integer,allocatable  :: cosmo_state_index(:)   ! index in cosmo state of every vertical level
  end type dart_variable_info

  integer,parameter              :: n_state_vector_vars=8
  integer,parameter              :: n_non_state_vars=1

  character(len=256)             :: string
  logical, save                  :: module_initialized = .FALSE.

  type(cosmo_meta),allocatable   :: cosmo_slabs(:)
  type(cosmo_hcoord)             :: cosmo_lonlat(3) ! 3 is for the stagger
  integer                        :: nslabs

  character(len=256)             :: cosmo_filename               = "test.grb"
  integer                        :: model_dt                     = 40
  logical                        :: output_state_vector          = .FALSE.
  real(r8)                       :: model_perturbation_amplitude = 0.1
  logical                        :: verbose                      = .false.

  namelist /model_nml/  &
    cosmo_filename, model_dt, model_perturbation_amplitude, output_state_vector,&
    model_dims,record_length,verbose

  integer                        :: model_size
  type(time_type)                :: model_timestep ! smallest time to adv model

  integer, parameter             :: n_max_kinds=200

  integer                        :: allowed_state_vector_vars(n_state_vector_vars)
  integer                        :: allowed_non_state_vars(1:n_max_kinds)
  logical                        :: is_allowed_state_vector_var(n_max_kinds)
  logical                        :: is_allowed_non_state_var(n_max_kinds)
  type(dart_variable_info)       :: state_vector_vars(1:n_max_kinds)
  type(cosmo_non_state_data)     :: non_state_data

  real(r8),allocatable           :: state_vector(:)

  real(r8),allocatable           :: ens_mean(:)  

  type(random_seq_type) :: random_seq

  character(len=256)    :: error_string,error_string2

  type(time_type)                :: cosmo_fc_time
  type(time_type)                :: cosmo_an_time

  type(grib_header_type),allocatable  :: grib_header(:)

contains

  function get_model_size()

    integer :: get_model_size
    
    if ( .not. module_initialized ) call static_init_model
    
    get_model_size = model_size

  end function get_model_size




  subroutine static_init_model()

    integer                       :: iunit,io,islab,ikind,sv_length
    real(r8),allocatable          :: data(:,:)

    real(r8),parameter            :: g = 9.80665_r8

    if ( module_initialized ) return ! only need to do this once.
    
    module_initialized=.TRUE.

    call set_allowed_state_vector_vars()

    ! read the DART namelist for this model
    call find_namelist_in_file('input.nml', 'model_nml', iunit)
    read(iunit, nml = model_nml, iostat = io)
    call check_namelist_read(iunit, io, 'model_nml')

    call set_calendar_type('Gregorian')

    call get_cosmo_info(cosmo_filename,cosmo_slabs,cosmo_lonlat,grib_header,&
                        is_allowed_state_vector_var,cosmo_fc_time)

    state_vector_vars(:)%is_present=.false.
    non_state_data%orography_present=.false.
    non_state_data%pressure_perturbation_present=.false.

    model_size = maxval(cosmo_slabs(:)%dart_eindex)
    nslabs     = size(cosmo_slabs,1)

    sv_length = 0
    do islab = 1,nslabs
      ikind = cosmo_slabs(islab)%dart_kind
      if (ikind>0) then 
        if (is_allowed_state_vector_var(ikind).OR.(ikind==KIND_PRESSURE_PERTURBATION)) then
          sv_length = sv_length+cosmo_slabs(islab)%dims(1)*cosmo_slabs(islab)%dims(2)
        end if
      end if
    end do

    allocate(state_vector(1:sv_length))

    ! cycle through all GRIB records
    ! one record corresponds to one horizontal field
    do islab=1,nslabs
      ikind=cosmo_slabs(islab)%dart_kind

      ! check if variable is a possible state vector variable
      if (ikind>0) then 
        if (is_allowed_state_vector_var(ikind)) then
          
          ! check if state vector variable information has not already been read
          ! e.g. for another vertical level
          if (.not. state_vector_vars(ikind)%is_present) then
            ! assign the variable information
            state_vector_vars(ikind)%is_present   = .true.
            state_vector_vars(ikind)%varname_short= cosmo_slabs(islab)%varname_short
            state_vector_vars(ikind)%varname_long = cosmo_slabs(islab)%varname_long
            state_vector_vars(ikind)%units        = cosmo_slabs(islab)%units
            state_vector_vars(ikind)%nx           = cosmo_slabs(islab)%dims(1)
            state_vector_vars(ikind)%ny           = cosmo_slabs(islab)%dims(2)
            state_vector_vars(ikind)%nz           = cosmo_slabs(islab)%dims(3)
            state_vector_vars(ikind)%horizontal_coordinate=cosmo_slabs(islab)%hcoord_type
            if (state_vector_vars(ikind)%nz>1) then
              state_vector_vars(ikind)%vertical_coordinate=VERTISLEVEL
            else
              state_vector_vars(ikind)%vertical_coordinate=VERTISSURFACE
            end if
            
            allocate(state_vector_vars(ikind)%vertical_level(     1:state_vector_vars(ikind)%nz))
            allocate(state_vector_vars(ikind)%state_vector_sindex(1:state_vector_vars(ikind)%nz))
            allocate(state_vector_vars(ikind)%cosmo_state_index(  1:state_vector_vars(ikind)%nz))
            
          end if
          
          ! set vertical information for this vertical level (record/slab)
          state_vector_vars(ikind)%vertical_level(     cosmo_slabs(islab)%ilevel)=cosmo_slabs(islab)%dart_level
          state_vector_vars(ikind)%state_vector_sindex(cosmo_slabs(islab)%ilevel)=cosmo_slabs(islab)%dart_sindex
          state_vector_vars(ikind)%cosmo_state_index(  cosmo_slabs(islab)%ilevel)=islab
          
        end if
        
      ! check for non state vector data (e.g. surface elevation) needed to run DART
        if (is_allowed_non_state_var(ikind)) then
          if (ikind==KIND_SURFACE_ELEVATION) then
            allocate(data(1:cosmo_slabs(islab)%dims(1),1:cosmo_slabs(islab)%dims(2)))
            data=get_data_from_binary(cosmo_filename,grib_header(islab),cosmo_slabs(islab)%dims(1),cosmo_slabs(islab)%dims(2))
            if (.not. allocated(non_state_data%surface_orography)) then
              allocate(non_state_data%surface_orography(1:cosmo_slabs(islab)%dims(1),1:cosmo_slabs(islab)%dims(2)))
            end if
            non_state_data%surface_orography(:,:)=data(:,:)
            deallocate(data)
            non_state_data%orography_present=.true.
          end if
          if ((ikind==KIND_SURFACE_GEOPOTENTIAL).and.(.not. allocated(non_state_data%surface_orography))) then
            allocate(data(1:cosmo_slabs(islab)%dims(1),1:cosmo_slabs(islab)%dims(2)))
            data=get_data_from_binary(cosmo_filename,grib_header(islab),cosmo_slabs(islab)%dims(1),cosmo_slabs(islab)%dims(2))
            allocate(non_state_data%surface_orography(1:cosmo_slabs(islab)%dims(1),1:cosmo_slabs(islab)%dims(2)))
            non_state_data%surface_orography(:,:)=data(:,:)/g
            deallocate(data)
            non_state_data%orography_present=.true.
          end if
          if (ikind==KIND_PRESSURE_PERTURBATION) then
            allocate(data(1:cosmo_slabs(islab)%dims(1),1:cosmo_slabs(islab)%dims(2)))
            data=get_data_from_binary(cosmo_filename,grib_header(islab),cosmo_slabs(islab)%dims(1),cosmo_slabs(islab)%dims(2))
            if (.not. allocated(non_state_data%pressure_perturbation)) then
              allocate(non_state_data%pressure_perturbation(1:cosmo_slabs(islab)%dims(1),1:cosmo_slabs(islab)%dims(2),1:cosmo_slabs(islab)%dims(3)))
            end if
            non_state_data%pressure_perturbation(:,:,cosmo_slabs(islab)%ilevel)=data(:,:)
            deallocate(data)
            
            if (.not. state_vector_vars(KIND_PRESSURE)%is_present) then
              ! assign the pressure variable information
              state_vector_vars(KIND_PRESSURE)%is_present   = .true.
              state_vector_vars(KIND_PRESSURE)%varname_short= cosmo_slabs(islab)%varname_short
              state_vector_vars(KIND_PRESSURE)%varname_long = cosmo_slabs(islab)%varname_long
              state_vector_vars(KIND_PRESSURE)%units        = cosmo_slabs(islab)%units
              state_vector_vars(KIND_PRESSURE)%nx           = cosmo_slabs(islab)%dims(1)
              state_vector_vars(KIND_PRESSURE)%ny           = cosmo_slabs(islab)%dims(2)
              state_vector_vars(KIND_PRESSURE)%nz           = cosmo_slabs(islab)%dims(3)
              state_vector_vars(KIND_PRESSURE)%horizontal_coordinate=cosmo_slabs(islab)%hcoord_type
              state_vector_vars(KIND_PRESSURE)%vertical_coordinate=VERTISLEVEL
              allocate(state_vector_vars(KIND_PRESSURE)%vertical_level(     1:state_vector_vars(KIND_PRESSURE)%nz))
              allocate(state_vector_vars(KIND_PRESSURE)%state_vector_sindex(1:state_vector_vars(KIND_PRESSURE)%nz))
              allocate(state_vector_vars(KIND_PRESSURE)%cosmo_state_index(  1:state_vector_vars(KIND_PRESSURE)%nz))
            end if
            
            ! set vertical information for this vertical level (record/slab)
            state_vector_vars(KIND_PRESSURE)%vertical_level(     cosmo_slabs(islab)%ilevel)=cosmo_slabs(islab)%dart_level
            state_vector_vars(KIND_PRESSURE)%state_vector_sindex(cosmo_slabs(islab)%ilevel)=cosmo_slabs(islab)%dart_sindex
            state_vector_vars(KIND_PRESSURE)%cosmo_state_index(  cosmo_slabs(islab)%ilevel)=islab
            
            non_state_data%pressure_perturbation_present=.true.
          end if
          
        end if
      end if
    end do

    ! set up the vertical coordinate system information
    !   search for one 3D variable (U-wind component should be contained in every analysis file)
    setlevel : do islab=1,nslabs
      if (cosmo_slabs(islab)%dart_kind==KIND_U_WIND_COMPONENT) then

        ! calculate the vertical coordinates for every grid point
        call set_vertical_coords(grib_header(islab),non_state_data,state_vector_vars(KIND_PRESSURE)%state_vector_sindex(:),state_vector)

        exit setlevel
      end if
    end do setlevel

    return
  end subroutine static_init_model



  subroutine get_state_meta_data(index_in,location,var_type)

    integer, intent(in)            :: index_in
    type(location_type)            :: location
    integer, optional, intent(out) :: var_type

    integer                        :: islab,var,hindex,dims(3)
    real(r8)                       :: lon,lat,vloc

    if (.NOT. module_initialized) CALL static_init_model()
    
    var=-1

    findindex : DO islab=1,nslabs
      IF ((index_in >= cosmo_slabs(islab)%dart_sindex) .AND. (index_in <= cosmo_slabs(islab)%dart_eindex)) THEN
        var      = islab
        hindex   = index_in-cosmo_slabs(islab)%dart_sindex+1
        var_type = cosmo_slabs(islab)%dart_kind
        dims     = cosmo_slabs(islab)%dims
        vloc     = cosmo_slabs(islab)%dart_level
        lon      = cosmo_lonlat(cosmo_slabs(islab)%hcoord_type)%lon(hindex)
        lat      = cosmo_lonlat(cosmo_slabs(islab)%hcoord_type)%lat(hindex)
        location = set_location(lon,lat,vloc,VERTISLEVEL)
        EXIT findindex
      END IF
    END DO findindex

    IF( var == -1 ) THEN
      write(string,*) 'Problem, cannot find base_offset, index_in is: ', index_in
      call error_handler(E_ERR,'get_state_meta_data',string,source,revision,revdate)
    ENDIF

  end subroutine get_state_meta_data



  function get_model_time_step()
  ! Returns the smallest increment of time that we want to advance the model.
  ! This defines the minimum assimilation interval.
  ! It is NOT the dynamical timestep of the model.

    type(time_type) :: get_model_time_step
    
    if ( .not. module_initialized ) call static_init_model
    
    model_timestep      = set_time(model_dt)
    get_model_time_step = model_timestep
    return

  end function get_model_time_step



  subroutine model_interpolate(x, location, obs_type, interp_val, istatus)

  ! Error codes:
  ! istatus = 99 : unknown error
  ! istatus = 10 : observation type is not in state vector
  ! istatus = 15 : observation lies outside the model domain (horizontal)
  ! istatus = 16 : observation lies outside the model domain (vertical)
  ! istatus = 19 : observation vertical coordinate is not supported

  ! Passed variables
  
  real(r8),            intent(in)  :: x(:)
  type(location_type), intent(in)  :: location
  integer,             intent(in)  :: obs_type
  real(r8),            intent(out) :: interp_val
  integer,             intent(out) :: istatus
  
  ! Local storage
  
  real(r8)             :: point_coords(1:3)
  
  integer              :: i,j,hbox(2,2),n,vbound(2),sindex
  real(r8)             :: hbox_weight(2,2),hbox_val(2,2),hbox_lon(2,2),hbox_lat(2,2)
  real(r8)             :: vbound_weight(2),val1,val2

  IF ( .not. module_initialized ) call static_init_model
  
  interp_val = MISSING_R8     ! the DART bad value flag
  istatus = 99                ! unknown error

  ! FIXME ... want some sort of error message here?
  if ( .not. state_vector_vars(obs_type)%is_present) then
     istatus=10
     return
  end if

  ! horizontal interpolation

  n = size(cosmo_lonlat(state_vector_vars(obs_type)%horizontal_coordinate)%lon,1)

  point_coords(1:3) = get_location(location)

  ! Find grid indices of box enclosing the observation location
  call get_enclosing_grid_box_lonlat(cosmo_lonlat(state_vector_vars(obs_type)%horizontal_coordinate)%lon,&
                                     cosmo_lonlat(state_vector_vars(obs_type)%horizontal_coordinate)%lat,&
                                     point_coords(1:2),n,state_vector_vars(obs_type)%nx,                 &
                                     state_vector_vars(obs_type)%ny, hbox, hbox_weight)
   
  if (hbox(1,1)==-1) then
     istatus=15
     return
  end if

  ! determine vertical level above and below obsevation
  call get_vertical_boundaries(hbox, hbox_weight, obs_type, query_location(location,'which_vert'),&
                               point_coords(3), vbound, vbound_weight, istatus)

  ! check if observation is in vertical domain and vertical coordinate system is supported
  ! FIXME istatus value?
  if (vbound(1)==-1) then
     return
  end if
   
  ! Perform a bilinear interpolation from the grid box to the desired location
  ! for the level above and below the observation

  sindex=state_vector_vars(obs_type)%state_vector_sindex(vbound(1))

  do i=1,2
  do j=1,2
     hbox_val(i,j)=x(sindex+hbox(i,j)-1)
     hbox_lon(i,j)=cosmo_lonlat(state_vector_vars(obs_type)%horizontal_coordinate)%lon(hbox(i,j))
     hbox_lat(i,j)=cosmo_lonlat(state_vector_vars(obs_type)%horizontal_coordinate)%lat(hbox(i,j))
  enddo
  enddo

  call bilinear_interpolation(hbox_val,hbox_lon,hbox_lat,point_coords,val1)

  sindex=state_vector_vars(obs_type)%state_vector_sindex(vbound(2))
  do i=1,2
  do j=1,2
     hbox_val(i,j)=x(sindex+hbox(i,j)-1)
  end do
  end do

  call bilinear_interpolation(hbox_val,hbox_lon,hbox_lat,point_coords,val2)

  ! vertical interpolation of horizontally interpolated values

  interp_val=val1*vbound_weight(1)+val2*vbound_weight(2)
  istatus=0
   
  end subroutine model_interpolate



  subroutine init_conditions(x)
  !------------------------------------------------------------------
  ! Returns a model state vector, x, that is some sort of appropriate
  ! initial condition for starting up a long integration of the model.
  ! At present, this is only used if the namelist parameter 
  ! start_from_restart is set to .false. in the program perfect_model_obs.

  real(r8), intent(out) :: x(:)

  if ( .not. module_initialized ) call static_init_model

  x = 0.0_r8  ! suppress compiler warnings about unused variables

  write(string,*) 'Cannot initialize COSMO state via subroutine call; start_from_restart cannot be F'
  call error_handler(E_ERR,'init_conditions',string,source,revision,revdate)

  end subroutine init_conditions



  subroutine init_time(time)
  !------------------------------------------------------------------
  ! Companion interface to init_conditions. Returns a time that is somehow 
  ! appropriate for starting up a long integration of the model.
  ! At present, this is only used if the namelist parameter 
  ! start_from_restart is set to .false. in the program perfect_model_obs.

  type(time_type), intent(out) :: time

  time = set_time(0,0) ! suppress compiler warnings about unused variables

  write(string,*) 'Cannot initialize COSMO time via subroutine call; start_from_restart cannot be F'
  call error_handler(E_ERR,'init_time',string,source,revision,revdate)
  
  end subroutine init_time


  
  subroutine adv_1step(x, time)
  !------------------------------------------------------------------
  ! As COSMO can only be advanced as a separate executable,
  ! this is a NULL INTERFACE.
  !------------------------------------------------------------------
    
  real(r8),        intent(inout) :: x(:)
  type(time_type), intent(in)    :: time
    
  if ( .not. module_initialized ) call static_init_model
    
  if (do_output()) then
      call print_time(time,'NULL interface adv_1step (no advance) DART time is')
      call print_time(time,'NULL interface adv_1step (no advance) DART time is',logfileunit)
  endif

  write(string,*) 'Cannot advance COSMO with a subroutine call; async cannot equal 0'
  call error_handler(E_ERR,'adv_1step',string,source,revision,revdate)

  end subroutine adv_1step


  
  subroutine end_model()

  deallocate(cosmo_slabs)
  deallocate(state_vector)

  end subroutine end_model



  function nc_write_model_atts( ncFileID ) result (ierr)
    
    integer, intent(in)  :: ncFileID      ! netCDF file identifier
    integer              :: ierr          ! return value of function

    integer              :: nDimensions, nVariables, nAttributes, unlimitedDimID, TimeDimID
    integer              :: StateVarDimID   ! netCDF pointer to state variable dimension (model size)
    integer              :: MemberDimID     ! netCDF pointer to dimension of ensemble    (ens_size)
    integer              :: LineLenDimID
    integer              :: StateVarVarID,StateVarID,VarID
    integer              :: ikind,ndims,idim,dims(100),nx,ny,nz,i
    character(len=6)     :: ckind

    integer              :: lonDimID, latDimID, levDimID, wlevDimID
    integer              :: lonVarID, latVarID, ulonVarID, ulatVarID, vlonVarID, vlatVarID
    integer              :: levVarID, wlevVarID

    character(len=128)   :: filename
    real(r8)             :: levs(1:500),wlevs(1:501)
    real(r8),allocatable :: data2d(:,:)

    character(len=8)      :: crdate      ! needed by F90 DATE_AND_TIME intrinsic
    character(len=10)     :: crtime      ! needed by F90 DATE_AND_TIME intrinsic
    character(len=5)      :: crzone      ! needed by F90 DATE_AND_TIME intrinsic
    integer, dimension(8) :: values      ! needed by F90 DATE_AND_TIME intrinsic

    logical :: has_std_latlon, has_ustag_latlon, has_vstag_latlon
   
    if ( .not. module_initialized ) call static_init_model

    ierr = -1 ! assume things go poorly

    has_std_latlon   = .FALSE.
    has_ustag_latlon = .FALSE.
    has_vstag_latlon = .FALSE.

    if (allocated(cosmo_lonlat(1)%lon) .and. allocated(cosmo_lonlat(1)%lat)) has_std_latlon   = .TRUE.
    if (allocated(cosmo_lonlat(2)%lon) .and. allocated(cosmo_lonlat(2)%lat)) has_ustag_latlon = .TRUE.
    if (allocated(cosmo_lonlat(3)%lon) .and. allocated(cosmo_lonlat(3)%lat)) has_vstag_latlon = .TRUE.

    write(filename,*) 'ncFileID', ncFileID

    !-------------------------------------------------------------------------------
    ! make sure ncFileID refers to an open netCDF file, 
    ! and then put into define mode.
    !-------------------------------------------------------------------------------
    
    call nc_check(nf90_Inquire(ncFileID,nDimensions,nVariables,nAttributes,unlimitedDimID),&
                                       'nc_write_model_atts', 'inquire '//trim(filename))
    call nc_check(nf90_Redef(ncFileID),'nc_write_model_atts',   'redef '//trim(filename))
    
    !-------------------------------------------------------------------------------
    ! We need the dimension ID for the number of copies/ensemble members, and
    ! we might as well check to make sure that Time is the Unlimited dimension. 
    ! Our job is create the 'model size' dimension.
    !-------------------------------------------------------------------------------
    
    call nc_check(nf90_inq_dimid(ncid=ncFileID, name='NMLlinelen', dimid=LineLenDimID), &
     'nc_write_model_atts','inq_dimid NMLlinelen')
    call nc_check(nf90_inq_dimid(ncid=ncFileID, name='copy', dimid=MemberDimID), &
     'nc_write_model_atts', 'copy dimid '//trim(filename))
    call nc_check(nf90_inq_dimid(ncid=ncFileID, name='time', dimid=  TimeDimID), &
     'nc_write_model_atts', 'time dimid '//trim(filename))
    
    if ( TimeDimID /= unlimitedDimId ) then
      write(error_string,*)'Time Dimension ID ',TimeDimID, &
       ' should equal Unlimited Dimension ID',unlimitedDimID
!      call error_handler(E_ERR,'nc_write_model_atts', error_string, source, revision, revdate)
    endif

    !-------------------------------------------------------------------------------
    ! Define the model size / state variable dimension / whatever ...
    !-------------------------------------------------------------------------------
    call nc_check(nf90_def_dim(ncid=ncFileID, name='StateVariable', len=model_size, &
     dimid = StateVarDimID),'nc_write_model_atts', 'state def_dim '//trim(filename))
    
    !-------------------------------------------------------------------------------
    ! Write Global Attributes 
    !-------------------------------------------------------------------------------

     call DATE_AND_TIME(crdate,crtime,crzone,values)
     write(string,'(''YYYY MM DD HH MM SS = '',i4,5(1x,i2.2))') &
      values(1), values(2), values(3), values(5), values(6), values(7)
    
     call nc_check(nf90_put_att(ncFileID, NF90_GLOBAL, 'creation_date' ,string    ), &
                   'nc_write_model_atts', 'creation put '//trim(filename))
     call nc_check(nf90_put_att(ncFileID, NF90_GLOBAL, 'model_source'  ,source  ), &
                   'nc_write_model_atts', 'source put '//trim(filename))
     call nc_check(nf90_put_att(ncFileID, NF90_GLOBAL, 'model_revision',revision), &
                   'nc_write_model_atts', 'revision put '//trim(filename))
     call nc_check(nf90_put_att(ncFileID, NF90_GLOBAL, 'model_revdate' ,revdate ), &
                   'nc_write_model_atts', 'revdate put '//trim(filename))
     call nc_check(nf90_put_att(ncFileID, NF90_GLOBAL, 'model',  'cosmo' ), &
                   'nc_write_model_atts', 'model put '//trim(filename))
    
    !-------------------------------------------------------------------------------
    ! Here is the extensible part. The simplest scenario is to output the state vector,
    ! parsing the state vector into model-specific parts is complicated, and you need
    ! to know the geometry, the output variables (PS,U,V,T,Q,...) etc. We're skipping
    ! complicated part.
    !-------------------------------------------------------------------------------
    
    if ( output_state_vector ) then
      
      !----------------------------------------------------------------------------
      ! Create a variable for the state vector
      !----------------------------------------------------------------------------
      
      ! Define the state vector coordinate variable and some attributes.
      call nc_check(nf90_def_var(ncid=ncFileID,name='StateVariable', xtype=nf90_int, &
                    dimids=StateVarDimID, varid=StateVarVarID), 'nc_write_model_atts', &
                    'statevariable def_var '//trim(filename))
      call nc_check(nf90_put_att(ncFileID,StateVarVarID,'long_name','State Variable ID'),&
                    'nc_write_model_atts','statevariable long_name '//trim(filename))
      call nc_check(nf90_put_att(ncFileID, StateVarVarID, 'units','indexical'), &
                    'nc_write_model_atts', 'statevariable units '//trim(filename))
      call nc_check(nf90_put_att(ncFileID,StateVarVarID,'valid_range',(/ 1,model_size /)),&
                    'nc_write_model_atts', 'statevariable valid_range '//trim(filename))
      
      ! Define the actual (3D) state vector, which gets filled as time goes on ... 
      call nc_check(nf90_def_var(ncid=ncFileID, name='state', xtype=nf90_real, &
                    dimids=(/StateVarDimID,MemberDimID,unlimitedDimID/),varid=StateVarID),&
                    'nc_write_model_atts','state def_var '//trim(filename))
      call nc_check(nf90_put_att(ncFileID,StateVarID,'long_name','model state or fcopy'),&
                    'nc_write_model_atts', 'state long_name '//trim(filename))
                   
      ! Leave define mode so we can fill the coordinate variable.
      call nc_check(nf90_enddef(ncfileID),'nc_write_model_atts','state enddef '//trim(filename))
      
      ! Fill the state variable coordinate variable
      call nc_check(nf90_put_var(ncFileID, StateVarVarID, (/ (i,i=1,model_size) /) ), &
                    'nc_write_model_atts', 'state put_var '//trim(filename))
      
    else
      
      !----------------------------------------------------------------------------
      ! We need to output the prognostic variables.
      !----------------------------------------------------------------------------
      ! Define the new dimensions IDs
      !----------------------------------------------------------------------------

      findnxny : do ikind=1,n_max_kinds
        if (state_vector_vars(ikind)%is_present) then
          nx=state_vector_vars(ikind)%nx
          ny=state_vector_vars(ikind)%ny
          exit findnxny
        end if
      end do findnxny

      findnz : do ikind=1,n_max_kinds
        if (state_vector_vars(ikind)%is_present) then
          if ((state_vector_vars(ikind)%nz>1) .and. (ikind .ne. KIND_VERTICAL_VELOCITY)) then
            nz=state_vector_vars(ikind)%nz
            exit findnz
          end if
        end if
      end do findnz

      call nc_check(nf90_def_dim(ncid=ncFileID, name='lon', len=nx, dimid = lonDimID), &
                     'nc_write_model_atts', 'lon def_dim '//trim(filename))
      call nc_check(nf90_def_dim(ncid=ncFileID, name='lat', len=ny, dimid = latDimID), &
                     'nc_write_model_atts', 'lat def_dim '//trim(filename))
      call nc_check(nf90_def_dim(ncid=ncFileID, name='lev', len=nz, dimid = levDimID), &
                     'nc_write_model_atts', 'lev def_dim '//trim(filename))
      call nc_check(nf90_def_dim(ncid=ncFileID, name='wlev', len=nz+1, dimid = wlevDimID), &
                     'nc_write_model_atts', 'lev def_dim '//trim(filename))

      if ( has_std_latlon ) then 
      ! Standard Grid Longitudes
      call nc_check(nf90_def_var(ncFileID,name='LON', xtype=nf90_real, &
                    dimids=(/ lonDimID, latDimID /), varid=lonVarID),&
                    'nc_write_model_atts', 'LON def_var '//trim(filename))
      call nc_check(nf90_put_att(ncFileID,  lonVarID, 'long_name', 'longitudes of grid'), &
                    'nc_write_model_atts', 'LON long_name '//trim(filename))
      call nc_check(nf90_put_att(ncFileID,  lonVarID, 'cartesian_axis', 'X'),  &
                    'nc_write_model_atts', 'LON cartesian_axis '//trim(filename))
      call nc_check(nf90_put_att(ncFileID,  lonVarID, 'units', 'degrees_east'), &
                    'nc_write_model_atts', 'LON units '//trim(filename))
      call nc_check(nf90_put_att(ncFileID,  lonVarID, 'valid_range', (/ -180.0_r8, 360.0_r8 /)), &
                    'nc_write_model_atts', 'LON valid_range '//trim(filename))
      ! Standard Grid Latitudes
      call nc_check(nf90_def_var(ncFileID,name='LAT', xtype=nf90_real, &
                    dimids=(/ lonDimID, latDimID /), varid=latVarID),&
                    'nc_write_model_atts', 'LAT def_var '//trim(filename))
      call nc_check(nf90_put_att(ncFileID,  latVarID, 'long_name', 'latitudes of grid'), &
                    'nc_write_model_atts', 'LAT long_name '//trim(filename))
      call nc_check(nf90_put_att(ncFileID,  latVarID, 'cartesian_axis', 'Y'),  &
                    'nc_write_model_atts', 'LAT cartesian_axis '//trim(filename))
      call nc_check(nf90_put_att(ncFileID,  latVarID, 'units', 'degrees_east'), &
                    'nc_write_model_atts', 'LAT units '//trim(filename))
      call nc_check(nf90_put_att(ncFileID,  latVarID, 'valid_range', (/ -180.0_r8, 360.0_r8 /)), &
                    'nc_write_model_atts', 'LAT valid_range '//trim(filename))
      endif


      if ( has_ustag_latlon ) then 
      ! U Grid Longitudes
      call nc_check(nf90_def_var(ncFileID,name='ULON', xtype=nf90_real, &
                    dimids=(/ lonDimID, latDimID /), varid=ulonVarID),&
                    'nc_write_model_atts', 'ULON def_var '//trim(filename))
      call nc_check(nf90_put_att(ncFileID,  ulonVarID, 'long_name', 'longitudes for U-wind'), &
                    'nc_write_model_atts', 'ULON long_name '//trim(filename))
      call nc_check(nf90_put_att(ncFileID,  ulonVarID, 'cartesian_axis', 'X'),  &
                    'nc_write_model_atts', 'ULON cartesian_axis '//trim(filename))
      call nc_check(nf90_put_att(ncFileID,  ulonVarID, 'units', 'degrees_east'), &
                    'nc_write_model_atts', 'ULON units '//trim(filename))
      call nc_check(nf90_put_att(ncFileID,  ulonVarID, 'valid_range', (/ -180.0_r8, 360.0_r8 /)), &
                    'nc_write_model_atts', 'ULON valid_range '//trim(filename))
      ! U Grid Latitudes
      call nc_check(nf90_def_var(ncFileID,name='ULAT', xtype=nf90_real, &
                    dimids=(/ lonDimID, latDimID /), varid=ulatVarID),&
                    'nc_write_model_atts', 'ULAT def_var '//trim(filename))
      call nc_check(nf90_put_att(ncFileID,  ulatVarID, 'long_name', 'latitudes for U-wind'), &
                    'nc_write_model_atts', 'ULAT long_name '//trim(filename))
      call nc_check(nf90_put_att(ncFileID,  ulatVarID, 'cartesian_axis', 'Y'),  &
                    'nc_write_model_atts', 'ULAT cartesian_axis '//trim(filename))
      call nc_check(nf90_put_att(ncFileID,  ulatVarID, 'units', 'degrees_east'), &
                    'nc_write_model_atts', 'ULAT units '//trim(filename))
      call nc_check(nf90_put_att(ncFileID,  ulatVarID, 'valid_range', (/ -180.0_r8, 360.0_r8 /)), &
                    'nc_write_model_atts', 'ULAT valid_range '//trim(filename))
      endif


      if ( has_vstag_latlon ) then 
      ! V Grid Longitudes
      call nc_check(nf90_def_var(ncFileID,name='VLON', xtype=nf90_real, &
                    dimids=(/ lonDimID, latDimID /), varid=vlonVarID),&
                    'nc_write_model_atts', 'VLON def_var '//trim(filename))
      call nc_check(nf90_put_att(ncFileID,  vlonVarID, 'long_name', 'longitudes for V-wind'), &
                    'nc_write_model_atts', 'VLON long_name '//trim(filename))
      call nc_check(nf90_put_att(ncFileID,  vlonVarID, 'cartesian_axis', 'X'),  &
                    'nc_write_model_atts', 'VLON cartesian_axis '//trim(filename))
      call nc_check(nf90_put_att(ncFileID,  vlonVarID, 'units', 'degrees_east'), &
                    'nc_write_model_atts', 'VLON units '//trim(filename))
      call nc_check(nf90_put_att(ncFileID,  vlonVarID, 'valid_range', (/ -180.0_r8, 360.0_r8 /)), &
                    'nc_write_model_atts', 'VLON valid_range '//trim(filename))
      ! V Grid Latitudes
      call nc_check(nf90_def_var(ncFileID,name='VLAT', xtype=nf90_real, &
                    dimids=(/ lonDimID, latDimID /), varid=vlatVarID),&
                    'nc_write_model_atts', 'VLAT def_var '//trim(filename))
      call nc_check(nf90_put_att(ncFileID,  vlatVarID, 'long_name', 'latitudes for V-wind'), &
                    'nc_write_model_atts', 'VLAT long_name '//trim(filename))
      call nc_check(nf90_put_att(ncFileID,  vlatVarID, 'cartesian_axis', 'Y'),  &
                    'nc_write_model_atts', 'VLAT cartesian_axis '//trim(filename))
      call nc_check(nf90_put_att(ncFileID,  vlatVarID, 'units', 'degrees_east'), &
                    'nc_write_model_atts', 'VLAT units '//trim(filename))
      call nc_check(nf90_put_att(ncFileID,  vlatVarID, 'valid_range', (/ -180.0_r8, 360.0_r8 /)), &
                    'nc_write_model_atts', 'VLAT valid_range '//trim(filename))
      endif

      ! Standard Z Levels
      call nc_check(nf90_def_var(ncFileID,name='LEV', xtype=nf90_real, &
                    dimids=(/ levDimID /), varid=levVarID),&
                    'nc_write_model_atts', 'LEV def_var '//trim(filename))
      call nc_check(nf90_put_att(ncFileID,  levVarID, 'long_name', 'standard hybrid model levels'), &
                    'nc_write_model_atts', 'LEV long_name '//trim(filename))
      call nc_check(nf90_put_att(ncFileID,  levVarID, 'cartesian_axis', 'Z'),  &
                    'nc_write_model_atts', 'LEV cartesian_axis '//trim(filename))
      call nc_check(nf90_put_att(ncFileID,  levVarID, 'units', 'model level'), &
                    'nc_write_model_atts', 'LEV units '//trim(filename))
      call nc_check(nf90_put_att(ncFileID,  levVarID, 'valid_range', (/ 1._r8,float(nz)+1._r8 /)), &
                    'nc_write_model_atts', 'LEV valid_range '//trim(filename))

      ! W-wind Z Levels
      call nc_check(nf90_def_var(ncFileID,name='WLEV', xtype=nf90_real, &
                    dimids=(/ wlevDimID /), varid=wlevVarID),&
                    'nc_write_model_atts', 'WLEV def_var '//trim(filename))
      call nc_check(nf90_put_att(ncFileID,  wlevVarID, 'long_name', 'standard model levels for W-wind'), &
                    'nc_write_model_atts', 'WLEV long_name '//trim(filename))
      call nc_check(nf90_put_att(ncFileID,  wlevVarID, 'cartesian_axis', 'Z'),  &
                    'nc_write_model_atts', 'WLEV cartesian_axis '//trim(filename))
      call nc_check(nf90_put_att(ncFileID,  wlevVarID, 'units', 'model level'), &
                    'nc_write_model_atts', 'WLEV units '//trim(filename))
      call nc_check(nf90_put_att(ncFileID,  wlevVarID, 'valid_range', (/ 1._r8,float(nz)+1._r8 /)), &
                    'nc_write_model_atts', 'WLEV valid_range '//trim(filename))

      ! DEBUG block to check shape of netCDF variables
      if ( 1 == 0 ) then
         write(*,*)'lon   dimid is ',lonDimID
         write(*,*)'lat   dimid is ',latDimID
         write(*,*)'lev   dimid is ',levDimID
         write(*,*)'wlev  dimid is ',wlevDimID
         write(*,*)'unlim dimid is ',unlimitedDimID
         write(*,*)'copy  dimid is ',MemberDimID
      endif

      do ikind=1,n_max_kinds
        if (state_vector_vars(ikind)%is_present) then

          error_string = trim(filename)//' '//trim(state_vector_vars(ikind)%varname_short)

          dims(1)=lonDimID
          dims(2)=latDimID

          idim=3
          if (state_vector_vars(ikind)%nz>1) then
            dims(idim)=levDimID
            if (ikind==KIND_VERTICAL_VELOCITY) then
              wlevs(1:nz+1)=state_vector_vars(ikind)%vertical_level(1:nz+1)
              dims(idim)=wlevDimID
            else
              levs(1:nz)=state_vector_vars(ikind)%vertical_level(1:nz)
            end if
            idim=idim+1
          end if

          ! Create a dimension for the ensemble 
          dims(idim) = memberDimID 
          idim=idim+1

          ! Put ensemble member dimension here
          dims(idim) = unlimitedDimID
          ndims=idim

          ! check shape of netCDF variables
          if ( verbose ) &
          write(*,*)trim(state_vector_vars(ikind)%varname_short),' has netCDF dimIDs ',dims(1:ndims)

          call nc_check(nf90_def_var(ncid=ncFileID, name=trim(state_vector_vars(ikind)%varname_short), xtype=nf90_real, &
                        dimids = dims(1:ndims), varid=VarID),&
                        'nc_write_model_atts', trim(error_string)//' def_var' )

          call nc_check(nf90_put_att(ncFileID, VarID, 'long_name', trim(state_vector_vars(ikind)%varname_long)), &
                        'nc_write_model_atts', trim(error_string)//' put_att long_name' )

          write(ckind,'(I6)') ikind
          call nc_check(nf90_put_att(ncFileID, VarID, 'DART_kind', trim(ckind)), &
                        'nc_write_model_atts', trim(error_string)//' put_att dart_kind' )

          call nc_check(nf90_put_att(ncFileID, VarID, 'units', trim(state_vector_vars(ikind)%units)), &
                        'nc_write_model_atts', trim(error_string)//' put_att units' )

        end if
      end do

      ! Leave define mode so we can fill the coordinate variable.
      call nc_check(nf90_enddef(ncfileID),'nc_write_model_atts','prognostic enddef '//trim(filename))

      !----------------------------------------------------------------------------
      ! Fill the coordinate variables - reshape the 1D arrays to the 2D shape
      !----------------------------------------------------------------------------
      allocate(data2d(nx,ny))

      if (has_std_latlon) then
      data2d = reshape(cosmo_lonlat(1)%lon, (/ nx, ny /) )
      call nc_check(nf90_put_var(ncFileID, lonVarID, data2d), &
                    'nc_write_model_atts', 'LON put_var '//trim(filename))

      data2d = reshape(cosmo_lonlat(1)%lat, (/ nx, ny /) )
      call nc_check(nf90_put_var(ncFileID, latVarID, data2d ), &
                    'nc_write_model_atts', 'LAT put_var '//trim(filename))
      endif
      
      if (has_ustag_latlon) then
      data2d = reshape(cosmo_lonlat(2)%lon, (/ nx, ny /) )
      call nc_check(nf90_put_var(ncFileID, ulonVarID, data2d ), &
                    'nc_write_model_atts', 'ULON put_var '//trim(filename))

      data2d = reshape(cosmo_lonlat(2)%lat, (/ nx, ny /) )
      call nc_check(nf90_put_var(ncFileID, ulatVarID, data2d ), &
                    'nc_write_model_atts', 'ULAT put_var '//trim(filename))
      endif
      
      if (has_vstag_latlon) then
      data2d = reshape(cosmo_lonlat(3)%lon, (/ nx, ny /) )
      call nc_check(nf90_put_var(ncFileID, vlonVarID, data2d ), &
                    'nc_write_model_atts', 'VLON put_var '//trim(filename))

      data2d = reshape(cosmo_lonlat(3)%lat, (/ nx, ny /) )
      call nc_check(nf90_put_var(ncFileID, vlatVarID, data2d ), &
                    'nc_write_model_atts', 'VLAT put_var '//trim(filename))
      endif

      deallocate(data2d)

      call nc_check(nf90_put_var(ncFileID, levVarID, levs(1:nz) ), &
                    'nc_write_model_atts', 'LEV put_var '//trim(filename))

      call nc_check(nf90_put_var(ncFileID, wlevVarID, wlevs(1:nz+1) ), &
                    'nc_write_model_atts', 'WLEV put_var '//trim(filename))
      
    end if

    !-------------------------------------------------------------------------------
    ! Flush the buffer and leave netCDF file open
    !-------------------------------------------------------------------------------
    call nc_check(nf90_sync(ncFileID), 'nc_write_model_atts', 'atts sync')
    
    ierr = 0 ! If we got here, things went well.
    
  end function nc_write_model_atts
 

 
  function nc_write_model_vars( ncFileID, state_vec, copyindex, timeindex ) result (ierr)         
    !------------------------------------------------------------------
    ! TJH 24 Oct 2006 -- Writes the model variables to a netCDF file.
    !
    ! TJH 29 Jul 2003 -- for the moment, all errors are fatal, so the
    ! return code is always '0 == normal', since the fatal errors stop execution.
    !
    ! For the lorenz_96 model, each state variable is at a separate location.
    ! that's all the model-specific attributes I can think of ...
    !
    ! assim_model_mod:init_diag_output uses information from the location_mod
    !     to define the location dimension and variable ID. All we need to do
    !     is query, verify, and fill ...
    !
    ! Typical sequence for adding new dimensions,variables,attributes:
    ! NF90_OPEN             ! open existing netCDF dataset
    !    NF90_redef         ! put into define mode
    !    NF90_def_dim       ! define additional dimensions (if any)
    !    NF90_def_var       ! define variables: from name, type, and dims
    !    NF90_put_att       ! assign attribute values
    ! NF90_ENDDEF           ! end definitions: leave define mode
    !    NF90_put_var       ! provide values for variable
    ! NF90_CLOSE            ! close: save updated netCDF dataset
    
    integer,                intent(in) :: ncFileID      ! netCDF file identifier
    real(r8), dimension(:), intent(in) :: state_vec
    integer,                intent(in) :: copyindex
    integer,                intent(in) :: timeindex
    integer                            :: ierr          ! return value of function
    
    integer, dimension(NF90_MAX_VAR_DIMS) :: dimIDs, mystart, mycount
    character(len=NF90_MAX_NAME)          :: varname 
    integer :: i,ikind, VarID, ncNdims, dimlen,ndims,vardims(3)
    integer :: TimeDimID, CopyDimID
    
    real(r8), allocatable, dimension(:,:)     :: data_2d_array
    real(r8), allocatable, dimension(:,:,:)   :: data_3d_array
    
    character(len=128) :: filename
    
    if ( .not. module_initialized ) call static_init_model
    
    ierr = -1 ! assume things go poorly
    
    !--------------------------------------------------------------------
    ! we only have a netcdf handle here so we do not know the filename
    ! or the fortran unit number.  but construct a string with at least
    ! the netcdf handle, so in case of error we can trace back to see
    ! which netcdf file is involved.
    !--------------------------------------------------------------------

    write(filename,*) 'ncFileID', ncFileID
    
    !-------------------------------------------------------------------------------
    ! make sure ncFileID refers to an open netCDF file, 
    !-------------------------------------------------------------------------------
    
    call nc_check(nf90_inq_dimid(ncFileID, 'copy', dimid=CopyDimID), &
     'nc_write_model_vars', 'inq_dimid copy '//trim(filename))
    
    call nc_check(nf90_inq_dimid(ncFileID, 'time', dimid=TimeDimID), &
     'nc_write_model_vars', 'inq_dimid time '//trim(filename))
    
    if ( output_state_vector ) then
      
      call nc_check(NF90_inq_varid(ncFileID, 'state', VarID), &
       'nc_write_model_vars', 'state inq_varid '//trim(filename))
      call nc_check(NF90_put_var(ncFileID,VarID,state_vec,start=(/1,copyindex,timeindex/)),&
       'nc_write_model_vars', 'state put_var '//trim(filename))
      
    else

      !----------------------------------------------------------------------------
      ! We need to process the prognostic variables.
      !----------------------------------------------------------------------------
      
      do ikind=1,n_max_kinds
        if (state_vector_vars(ikind)%is_present) then
          
          varname      = trim(state_vector_vars(ikind)%varname_short)
          error_string = trim(filename)//' '//trim(varname)

          ! Ensure netCDF variable is conformable with progvar quantity.
          ! The TIME and Copy dimensions are intentionally not queried
          ! by looping over the dimensions stored in the progvar type.
          
          call nc_check(nf90_inq_varid(ncFileID, varname, VarID), &
                        'nc_write_model_vars', 'inq_varid '//trim(error_string))
          
          call nc_check(nf90_inquire_variable(ncFileID,VarId,dimids=dimIDs,ndims=ncNdims), &
                        'nc_write_model_vars', 'inquire '//trim(error_string))

          mystart(:)=1
          mycount(:)=1
          
          if (state_vector_vars(ikind)%nz==1) then
            ndims=2
          else
            ndims=3
          end if

          vardims(1)=state_vector_vars(ikind)%nx
          vardims(2)=state_vector_vars(ikind)%ny
          vardims(3)=state_vector_vars(ikind)%nz

          DimCheck : do i = 1,ndims
            
            write(error_string,'(a,i2,A)') 'inquire dimension ',i,trim(varname)
            call nc_check(nf90_inquire_dimension(ncFileID, dimIDs(i), len=dimlen), &
             'nc_write_model_vars', trim(error_string))
            
            if ( dimlen /= vardims(i) ) then
              write(error_string,*) trim(varname),' dim/dimlen ',i,dimlen,' not ',vardims(i)
              write(error_string2,*)' but it should be.'
              call error_handler(E_ERR, 'nc_write_model_vars', trim(error_string), &
                              source, revision, revdate, text2=trim(error_string2))
            endif
            
            mycount(i) = dimlen
            
          end do DimCheck
        
          where(dimIDs == CopyDimID) mystart = copyindex
          where(dimIDs == CopyDimID) mycount = 1
          where(dimIDs == TimeDimID) mystart = timeindex
          where(dimIDs == TimeDimID) mycount = 1

          if (ndims==2) then
            allocate(data_2d_array(vardims(1),vardims(2)))
            call sv_to_field(data_2d_array,state_vec,state_vector_vars(ikind))
            call nc_check(nf90_put_var(ncFileID, VarID, data_2d_array, &
                          start = mystart(1:ncNdims), count=mycount(1:ncNdims)), &
                          'nc_write_model_vars', 'put_var '//trim(error_string2))
            deallocate(data_2d_array)

          elseif (ndims==3) then
            allocate(data_3d_array(vardims(1),vardims(2),vardims(3)))
            call sv_to_field(data_3d_array,state_vec,state_vector_vars(ikind))
            call nc_check(nf90_put_var(ncFileID, VarID, data_3d_array, &
                          start = mystart(1:ncNdims), count=mycount(1:ncNdims)), &
                          'nc_write_model_vars', 'put_var '//trim(error_string2))
            deallocate(data_3d_array)

          else
             write(error_string, *) 'no support for data array of dimension ', ncNdims
             call error_handler(E_ERR,'nc_write_model_vars', error_string, &
                           source,revision,revdate)
          end if

        end if

      end do

    end if

    return

  end function nc_write_model_vars



  subroutine pert_model_state(state, pert_state, interf_provided)
    !------------------------------------------------------------------
    ! Perturbs a model state for generating initial ensembles.
    ! The perturbed state is returned in pert_state.
    ! A model may choose to provide a NULL INTERFACE by returning
    ! .false. for the interf_provided argument. This indicates to
    ! the filter that if it needs to generate perturbed states, it
    ! may do so by adding a perturbation to each model state 
    ! variable independently. The interf_provided argument
    ! should be returned as .true. if the model wants to do its own
    ! perturbing of states.
    !------------------------------------------------------------------
    ! Currently only implemented as rondom perturbations
    !------------------------------------------------------------------    

    real(r8), intent(in)  :: state(:)
    real(r8), intent(out) :: pert_state(:)
    logical,  intent(out) :: interf_provided

    real(r8)              :: stddev,mean

    integer               :: ikind,ilevel,i,istart,iend
    logical, save         :: random_seq_init = .false.
    
    if ( .not. module_initialized ) call static_init_model
    
    interf_provided = .true.

    ! Initialize my random number sequence (no seed is submitted here!)
    if(.not. random_seq_init) then
      call init_random_seq(random_seq)
      random_seq_init = .true.
    endif
    
    ! add some uncertainty to every state vector element
    do ikind=1,size(state_vector_vars)
      if (state_vector_vars(ikind)%is_present) then
        do ilevel=1,state_vector_vars(ikind)%nz
          istart=state_vector_vars(ikind)%state_vector_sindex(ilevel)
          iend=istart+(state_vector_vars(ikind)%nx*state_vector_vars(ikind)%ny)-1

          mean=sum(abs(state(istart:iend)))/float(iend-istart+1)
          stddev=sqrt(sum((state(istart:iend)-mean)**2))/float(iend-istart+1)

          do i=istart,iend
            pert_state(i) = random_gaussian(random_seq, state(i),model_perturbation_amplitude*stddev)
          end do
          if ((ikind==KIND_SPECIFIC_HUMIDITY) .or. &
              (ikind==KIND_CLOUD_LIQUID_WATER) .or. &
              (ikind==KIND_CLOUD_ICE)) then
            where (pert_state(istart:iend)<0.)
              pert_state(istart:iend)=0.
            end where
          end if
        end do
      end if
    enddo

    return
    
  end subroutine pert_model_state



  subroutine ens_mean_for_model(filter_ens_mean)

    real(r8), dimension(:), intent(in) :: filter_ens_mean
    
    if ( .not. module_initialized ) call static_init_model
    
    allocate(ens_mean(1:model_size))
    ens_mean(:) = filter_ens_mean(:)

!  write(string,*) 'COSMO has no ensemble mean in storage.'
!  call error_handler(E_ERR,'ens_mean_for_model',string,source,revision,revdate)

  end subroutine ens_mean_for_model



  subroutine set_allowed_state_vector_vars()
    ! set the information on which variables should go into the state vector

    is_allowed_state_vector_var(:)=.FALSE.
    is_allowed_non_state_var(:)=.FALSE.

    allowed_state_vector_vars(1)=KIND_U_WIND_COMPONENT
     is_allowed_state_vector_var(KIND_U_WIND_COMPONENT)=.TRUE.
    allowed_state_vector_vars(1)=KIND_U_WIND_COMPONENT
     is_allowed_state_vector_var(KIND_U_WIND_COMPONENT)=.TRUE.
    allowed_state_vector_vars(2)=KIND_V_WIND_COMPONENT
     is_allowed_state_vector_var(KIND_V_WIND_COMPONENT)=.TRUE.
    allowed_state_vector_vars(3)=KIND_VERTICAL_VELOCITY
     is_allowed_state_vector_var(KIND_VERTICAL_VELOCITY)=.TRUE.
    allowed_state_vector_vars(4)=KIND_TEMPERATURE
     is_allowed_state_vector_var(KIND_TEMPERATURE)=.TRUE.
    allowed_state_vector_vars(5)=KIND_PRESSURE
     is_allowed_state_vector_var(KIND_PRESSURE)=.TRUE.
    allowed_state_vector_vars(6)=KIND_SPECIFIC_HUMIDITY
     is_allowed_state_vector_var(KIND_SPECIFIC_HUMIDITY)=.TRUE.
    allowed_state_vector_vars(7)=KIND_CLOUD_LIQUID_WATER
     is_allowed_state_vector_var(KIND_CLOUD_LIQUID_WATER)=.TRUE.
    allowed_state_vector_vars(8)=KIND_CLOUD_ICE
     is_allowed_state_vector_var(KIND_CLOUD_ICE)=.TRUE.

    ! set the information which variables are needed but will not go into the state vector
    allowed_non_state_vars(1)=KIND_SURFACE_ELEVATION
     is_allowed_non_state_var(KIND_SURFACE_ELEVATION)=.TRUE.
    allowed_non_state_vars(2)=KIND_SURFACE_GEOPOTENTIAL
     is_allowed_non_state_var(KIND_SURFACE_GEOPOTENTIAL)=.TRUE.
    allowed_non_state_vars(3)=KIND_PRESSURE_PERTURBATION
     is_allowed_non_state_var(KIND_PRESSURE_PERTURBATION)=.TRUE.

    return

  end subroutine set_allowed_state_vector_vars



  function ll_to_xyz_vector(lon,lat) RESULT (xyz)
    
    ! Passed variables
    
    real(r8),allocatable :: xyz(:,:)      ! result: x,z,y-coordinates
    real(r8),intent(in)  :: lat(:),lon(:) ! input:  lat/lon coordinates in degrees

    real(r8)             :: radius
    integer              :: n

    ! define output vector size to be the same as the input vector size
    ! second dimension (3) is x,y,z

    n=SIZE(lat,1)
    ALLOCATE(xyz(1:n,1:3))

    ! as we are interested in relative distances we set the radius to 1 - may be changed later

    radius=1.0_r8

    ! caclulate the x,y,z-coordinates

    xyz(1:n,1)=radius*sin(lat(1:n)*deg2rad)*cos(lon(1:n)*deg2rad)
    xyz(1:n,2)=radius*sin(lat(1:n)*deg2rad)*sin(lon(1:n)*deg2rad)
    xyz(1:n,3)=radius*cos(lat(1:n)*deg2rad)
    
    return
  end function ll_to_xyz_vector



  function ll_to_xyz_single(lon,lat) result (xyz)
    
    ! Passed variables
    
    real(r8)             :: xyz(1:3) ! result: x,z,y-coordinates
    real(r8),intent(in)  :: lat,lon  ! input:  lat/lon coordinates in degrees

    real(r8)             :: radius

    ! as we are interested in relative distances we set the radius to 1 - may be changed later

    radius=1.0_r8

    ! caclulate the x,y,z-coordinates

    xyz(1)=radius*sin(lat*deg2rad)*cos(lon*deg2rad)
    xyz(2)=radius*sin(lat*deg2rad)*sin(lon*deg2rad)
    xyz(3)=radius*cos(lat*deg2rad)
    
    return
  end function ll_to_xyz_single



  subroutine get_enclosing_grid_box(p,g,n,nx,ny,b,bw)
    
    integer,intent(in)   :: n,nx,ny
    real(r8),intent(in)  :: p(1:3),g(1:n,1:3)
    integer,intent(out)  :: b(1:2,1:2)
    real(r8),intent(out) :: bw(1:2,1:2)

!    real(r8)             :: work(1:nx,1:ny,1:3),dist(1:nx,1:ny),boxdist(1:2,1:2)
    real(r8)             :: work(1:nx+2,1:ny+2,1:3),dist(1:nx+2,1:ny+2),boxdist(1:2,1:2)
    integer              :: i,j,minidx(2),boxidx(2),xb,yb

    real(r8) :: sqrt2

    sqrt2 = sqrt(2.0_r8)

    work(2:nx+1,2:ny+1,1:3)=RESHAPE( g, (/ nx,ny,3 /))

    do i=2,nx+1
      work(i,   1,1:3)=work(i,   2,1:3)-(work(i, 3,1:3)-work(i,   2,1:3))
      work(i,ny+2,1:3)=work(i,ny+1,1:3)-(work(i,ny,1:3)-work(i,ny+1,1:3))
    end do

    do j=2,ny+1
      work(   1,j,1:3)=work(   2,j,1:3)-(work( 3,j,1:3)-work(   2,j,1:3))
      work(nx+2,j,1:3)=work(nx+1,j,1:3)-(work(nx,j,1:3)-work(nx+1,j,1:3))
    end do

    work(   1,   1,1:3) = work(   2,   2,1:3) - 0.5_r8*(sqrt2*(work(   2,   2,1:3)-work(   1,   2,1:3)) + sqrt2*(work(   2,   2,1:3)-work(   2,   1,1:3)))
    work(   1,ny+2,1:3) = work(   2,ny+1,1:3) - 0.5_r8*(sqrt2*(work(   2,ny+1,1:3)-work(   1,ny+1,1:3)) + sqrt2*(work(   2,ny+1,1:3)-work(   2,ny+2,1:3)))
    work(nx+2,   1,1:3) = work(nx+1,   2,1:3) - 0.5_r8*(sqrt2*(work(nx+1,   2,1:3)-work(nx+2,   2,1:3)) + sqrt2*(work(nx+1,   2,1:3)-work(nx+1,   1,1:3)))
    work(nx+2,ny+2,1:3) = work(nx+1,ny+1,1:3) - 0.5_r8*(sqrt2*(work(nx+1,ny+1,1:3)-work(nx+2,ny+1,1:3)) + sqrt2*(work(nx+1,ny+1,1:3)-work(nx+1,ny+2,1:3)))

    do i=1,nx+2
    do j=1,ny+2
        dist(i,j)=sqrt(sum((work(i,j,:)-p(:))**2))
    end do
    end do

    minidx(:)=minloc(dist)

    ! watch for out of area values

    if (minidx(1)==1 .or. minidx(1)==(nx+2) .or. minidx(2)==1 .or. minidx(2)==(ny+2)) then
      b(:,:)=-1
      return
    end if


    do i=0,1
    do j=0,1
        boxdist(i+1,j+1)=sum(dist(minidx(1)+i-1:minidx(1)+i,minidx(2)+j-1:minidx(2)+j))
    end do
    end do 

    boxidx=minloc(boxdist)-1

    xb=minidx(1)+(2*(boxidx(1)-0.5))
    yb=minidx(2)+(2*(boxidx(2)-0.5))

    if (xb==1 .or. xb==(nx+2) .or. yb==1 .or. yb==(ny+2)) then
      b(:,:)=-1
      return
    else
      do i=1,2
      do j=1,2
          b(i,j)=((minidx(2)+(j-1)*(boxidx(2)-0.5)*2)*ny)+(minidx(1)+(i-1)*(2*(boxidx(1)-0.5)))
      end do
      end do

      do i=1,2
      do j=1,2
          boxdist(i,j)=dist(mod(b(i,j),ny),b(i,j)/ny)
      end do
      end do

      bw(:,:)=1./boxdist(:,:)
!      bw(:,:)=(((1.-boxdist(:,:))/(1.1*maxval(boxdist)))**2)/((boxdist(:,:)/(1.1*maxval(boxdist)))**2)
      bw=bw/sum(bw)
      b(:,:)=b(:,:)-1
    end if

  end subroutine get_enclosing_grid_box



  subroutine get_enclosing_grid_box_lonlat(lon,lat,p,n,nx,ny,b,bw)
    
    integer, intent(in)  :: n,nx,ny
    real(r8),intent(in)  :: p(1:2),lon(1:n),lat(1:n)
    integer, intent(out) :: b(1:2,1:2)
    real(r8),intent(out) :: bw(1:2,1:2)

!    real(r8)            :: work(1:nx,1:ny,1:3),dist(1:nx,1:ny),boxdist(1:2,1:2)
    real(r8)             :: work(1:nx+2,1:ny+2,1:2),dist(1:nx+2,1:ny+2),boxdist(1:2,1:2),pw(2)

    integer  :: i,j,minidx(2),boxidx(2),xb,yb,bx(2,2),by(2,2)
    real(r8) :: sqrt2

    sqrt2 = sqrt(2.0_r8)

    work(2:nx+1,2:ny+1,1)=reshape(lon,(/ nx,ny /))*deg2rad
    work(2:nx+1,2:ny+1,2)=reshape(lat,(/ nx,ny /))*deg2rad
    pw=p*deg2rad

    do i=2,nx+1
      work(i,   1,1:2)=work(i,   2,1:2)-(work(i, 3,1:2)-work(i,   2,1:2))
      work(i,ny+2,1:2)=work(i,ny+1,1:2)-(work(i,ny,1:2)-work(i,ny+1,1:2))
    end do

    do j=2,ny+1
      work(   1,j,1:2)=work(   2,j,1:2)-(work( 3,j,1:2)-work(   2,j,1:2))
      work(nx+2,j,1:2)=work(nx+1,j,1:2)-(work(nx,j,1:2)-work(nx+1,j,1:2))
    end do

    work(   1,   1,1:2) = work(   2,   2,1:2) - 0.5_r8*(sqrt2*(work(   2,   2,1:2)-work(   1,   2,1:2))+sqrt2*(work(   2,   2,1:2)-work(   2,   1,1:2)))
    work(   1,ny+2,1:2) = work(   2,ny+1,1:2) - 0.5_r8*(sqrt2*(work(   2,ny+1,1:2)-work(   1,ny+1,1:2))+sqrt2*(work(   2,ny+1,1:2)-work(   2,ny+2,1:2)))
    work(nx+2,   1,1:2) = work(nx+1,   2,1:2) - 0.5_r8*(sqrt2*(work(nx+1,   2,1:2)-work(nx+2,   2,1:2))+sqrt2*(work(nx+1,   2,1:2)-work(nx+1,   1,1:2)))
    work(nx+2,ny+2,1:2) = work(nx+1,ny+1,1:2) - 0.5_r8*(sqrt2*(work(nx+1,ny+1,1:2)-work(nx+2,ny+1,1:2))+sqrt2*(work(nx+1,ny+1,1:2)-work(nx+1,ny+2,1:2)))

    do i=1,nx+2
    do j=1,ny+2
!      dist(i,j)=sqrt(sum((work(i,j,:)-p(:))**2))
       dist(i,j) = 6173.0_r8*acos(cos(work(i,j,2)-pw(2))-cos(work(i,j,2))*cos(pw(2))*(1-cos(work(i,j,1)-pw(1))))
    end do
    end do

    minidx(:)=minloc(dist)

    ! watch for out of area values

    if (minidx(1)==1 .or. minidx(1)==(nx+2) .or. minidx(2)==1 .or. minidx(2)==(ny+2)) then
      b(:,:)=-1
      return
    end if

!   open(21,file='/daten02/jkeller/testbox.bin',form='unformatted')
!   iunit = open_file('testbox.bin',form='unformatted',action='write') 
!   write(iunit) nx
!   write(iunit) ny

    do i=0,1
    do j=0,1
        boxdist(i+1,j+1)=sum(dist(minidx(1)+i-1:minidx(1)+i,minidx(2)+j-1:minidx(2)+j))/4.0_r8
!       write(*,'(4(I5))') minidx(1)+i-1,minidx(1)+i,minidx(2)+j-1,minidx(2)+j
!       write(iunit) (minidx(2)+j-1),minidx(1)+i-1,&
!                 (minidx(2)+j-1),minidx(1)+i,&
!                 (minidx(2)+j),minidx(1)+i-1,&
!                 (minidx(2)+j),minidx(1)+i
!       write(iunit) boxdist(i+1,j+1)
    end do
    end do 

    boxidx=minloc(boxdist)-1

    xb=minidx(1)+(2*(boxidx(1)-0.5_r8))
    yb=minidx(2)+(2*(boxidx(2)-0.5_r8))

    if (xb==1 .or. xb==(nx+2) .or. yb==1 .or. yb==(ny+2)) then
      b(:,:)=-1
      return
    else
      do i=1,2
      do j=1,2
          bx(i,j)=minidx(1)+(i-1)*(2*(boxidx(1)-0.5_r8))
          by(i,j)=minidx(2)+(j-1)*(2*(boxidx(2)-0.5_r8))
      end do
      end do

      do i=1,2
      do j=1,2
          boxdist(i,j)=dist(bx(i,j),by(i,j))
      end do
      end do

      bw(:,:)=1.0_r8/boxdist(:,:)
      bw=bw/sum(bw)
      bx=bx-1
      by=by-1
      b(:,:)=(by-1)*nx+bx
    end if

    return

  end subroutine get_enclosing_grid_box_lonlat


  subroutine bilinear_interpolation(bv,blo,bla,p,v)

    ! Passed variables
    
    real(r8),intent(in)  :: bv(2,2),blo(2,2),bla(2,2)
    real(r8),intent(in)  :: p(3)
    real(r8),intent(out) :: v
    
    ! Local storage
    
    real(r8)             :: x1,lo1,la1
    real(r8)             :: x2,lo2,la2
    real(r8)             :: d1,d2,d
    
!    write(*,'(3(F8.5,1X))') bv(1,1),blo(1,1),bla(1,1)
!    write(*,'(3(F8.5,1X))') bv(2,1),blo(2,1),bla(2,1)

    call linear_interpolation(p(1),p(2),bv(1,1),blo(1,1),bla(1,1),&
                                        bv(2,1),blo(2,1),bla(2,1),&
                                        x1,lo1,la1)

!    write(*,'(3(F8.5,1X))') x1,lo1,la1

    call linear_interpolation(p(1),p(2),bv(1,2),blo(1,2),bla(1,2),&
                                        bv(2,2),blo(2,2),bla(2,2),&
                                        x2,lo2,la2)

!    write(*,'(3(F8.5,1X))') x2,lo2,la2

    d1=sqrt((lo1-p(1))**2+(la1-p(2))**2)
    d2=sqrt((lo2-p(1))**2+(la2-p(2))**2)
    d =sqrt((lo1-lo2 )**2+(la1-la2 )**2)

    v=(1.0_r8-(d1/d))*x1+(1.0_r8-(d2/d))*x2

    return

  end subroutine bilinear_interpolation



  subroutine linear_interpolation(lop,lap,x1,lo1,la1,x2,lo2,la2,x,lo,la)

    real(r8),intent(in)  :: lo1,lo2,la1,la2,x1,x2,lop,lap
    real(r8),intent(out) :: lo,la,x

    real(r8)             :: m1,m2,n1,n2,d1,d2,d,mylo1,mylo2,mylop,w1,w2

    mylo1=lo1
    mylo2=lo2
    mylop=lop

    if (lo1>180.0_r8) mylo1=lo1-360.0_r8
    if (lo2>180.0_r8) mylo2=lo2-360.0_r8
    if (lop>180.0_r8) mylop=lop-360.0_r8

    m1=(la2-la1)/(mylo2-mylo1)
    if (m1 .ne. 0.0_r8) then
      n1=la1-mylo1*m1
      m2=-1.0_r8/m1
      n2=lap-mylop*m2
      lo=(n2-n1)/(m1-m2)
      la=lo*m1+n1
      d1=sqrt((mylo1-lo)**2+(la1-la)**2)
      d2=sqrt((mylo2-lo)**2+(la2-la)**2)
      d =sqrt((mylo1-mylo2)**2+(la1-la2)**2)
    else
      la=la1
      lo=mylop

      d1=sqrt((mylo1-lo)**2+(la1-la)**2)
      d2=sqrt((mylo2-lo)**2+(la2-la)**2)
      d =sqrt((mylo1-mylo2)**2+(la1-la2)**2)
    end if

    if (lo < 0.0_r8) lo=lo+360.0_r8

    w1=abs(1.0_r8-(d1/d))
    w2=abs(1.0_r8-(d2/d))
    x=w1*x1+w2*x2

    return

  end subroutine linear_interpolation



  subroutine get_vertical_boundaries(hb,hw,otype,vcs,p,b,w,istatus)

    real(r8), intent(in)  :: hw(2,2),p,vcs
    integer,  intent(in)  :: hb(2,2),otype
    integer,  intent(out) :: b(2),istatus
    real(r8), intent(out) :: w(2)

    integer               :: k,nlevel,x1,x2,x3,x4,y1,y2,y3,y4
    real(r8)              :: u,l
    real(r8),allocatable  :: klevel(:),hlevel(:),plevel(:)

    b(:)=-1

    ! coordinate system not implemented
    if ( (nint(vcs) == VERTISUNDEF)        .or. &
         (nint(vcs) == VERTISSURFACE)      .or. &
         (nint(vcs) == VERTISSCALEHEIGHT) ) then
      istatus=19
      return
    end if

! TJH    write(*,*)'non_state_data%pfl min max ',minval(non_state_data%pfl),maxval(non_state_data%pfl)
! TJH    write(*,*)' mean is ',sum(non_state_data%pfl)/(665.0_r8*657.0_r8*40.0_r8)

    x1 = mod(hb(1,1),size(non_state_data%pfl,1))
    x2 = mod(hb(2,1),size(non_state_data%pfl,1))
    x3 = mod(hb(1,2),size(non_state_data%pfl,1))
    x4 = mod(hb(2,2),size(non_state_data%pfl,1))
    y1 =     hb(1,1)/size(non_state_data%pfl,1)
    y2 =     hb(2,1)/size(non_state_data%pfl,1)
    y3 =     hb(1,2)/size(non_state_data%pfl,1)
    y4 =     hb(2,2)/size(non_state_data%pfl,1)

! TJH    write(*,*)'hb is ',hb
! TJH    write(*,*)'x  is ',x1,x2,x3,x4
! TJH    write(*,*)'y  is ',y1,y2,y3,y4
! TJH    write(*,*)'hw is ',hw
    
    if (otype .ne. KIND_VERTICAL_VELOCITY) then
      ! The variable exists on the 'full' levels
      nlevel=non_state_data%nfl
      allocate(klevel(1:nlevel))
      allocate(hlevel(1:nlevel))
      allocate(plevel(1:nlevel))
      klevel=state_vector_vars(otype)%vertical_level(:)
      hlevel=hw(1,1)*non_state_data%hfl(x1,y1,:)+&
             hw(2,1)*non_state_data%hfl(x2,y2,:)+&
             hw(1,2)*non_state_data%hfl(x3,y3,:)+&
             hw(2,2)*non_state_data%hfl(x4,y4,:)
      plevel=hw(1,1)*non_state_data%pfl(x1,y1,:)+&
             hw(2,1)*non_state_data%pfl(x2,y2,:)+&
             hw(1,2)*non_state_data%pfl(x3,y3,:)+&
             hw(2,2)*non_state_data%pfl(x4,y4,:)

! TJH write(*,*)non_state_data%pfl(x1,y1,:)
! TJH write(*,*)non_state_data%pfl(x2,y2,:)
! TJH write(*,*)non_state_data%pfl(x3,y3,:)
! TJH write(*,*)non_state_data%pfl(x4,y4,:)

    else
      ! The variable exists on the 'half' levels
      nlevel=non_state_data%nhl
      allocate(klevel(1:nlevel))
      allocate(hlevel(1:nlevel))
      allocate(plevel(1:nlevel))
      klevel=state_vector_vars(otype)%vertical_level(:)
      hlevel=hw(1,1)*non_state_data%hhl(x1,y1,:)+&
             hw(2,1)*non_state_data%hhl(x2,y2,:)+&
             hw(1,2)*non_state_data%hhl(x3,y3,:)+&
             hw(2,2)*non_state_data%hhl(x4,y4,:)
      plevel=hw(1,1)*non_state_data%phl(x1,y1,:)+&
             hw(2,1)*non_state_data%phl(x2,y2,:)+&
             hw(1,2)*non_state_data%phl(x3,y3,:)+&
             hw(2,2)*non_state_data%phl(x4,y4,:)
    end if

    u = -1.0_r8
    l = -1.0_r8

    do k=1,nlevel-1

      ! Find the bounding levels for the respective coordinate system 
      if (nint(vcs) == VERTISLEVEL) then
        u=klevel(k+1)
        l=klevel(k)
      end if
      if (nint(vcs) == VERTISPRESSURE) then
      ! write(*,*)' vert is pressure '
      ! write(*,*)'plevel is ',plevel
        u=plevel(k+1)
        l=plevel(k)
      end if
      if (nint(vcs) == VERTISHEIGHT) then
      ! write(*,*)' vert is height '
      ! write(*,*)'hlevel is ',hlevel
        u=hlevel(k+1)
        l=hlevel(k)
      end if

! TJH write(*,*)'u p l',u,p,l

      if (u>=p .and. l<=p) then
        b(1)=k
        b(2)=k+1
        w(1)=1.0_r8-(p-l)/(u-l)
        w(2)=1.0_r8-(u-p)/(u-l)
        return
      end if

    end do

    istatus=16 ! out of domain
    return
  end subroutine get_vertical_boundaries



  subroutine sv_to_field_2d(f,x,v)

    real(r8),intent(out)                :: f(:,:)
    real(r8),intent(in)                 :: x(:)
    type(dart_variable_info),intent(in) :: v

    integer                             :: is,ie

    is=v%state_vector_sindex(1)
    ie=is+v%nx*v%ny-1
    f(:,:) = reshape(x(is:ie),(/ v%nx,v%ny /))

    return

  end subroutine sv_to_field_2d



  subroutine sv_to_field_3d(f,x,v)

    real(r8),intent(out)                :: f(:,:,:)
    real(r8),intent(in)                 :: x(:)
    type(dart_variable_info),intent(in) :: v

    integer                             :: is,ie,iz

    do iz=1,v%nz
      is=v%state_vector_sindex(iz)
      ie=is+v%nx*v%ny-1
      f(:,:,iz) = reshape(x(is:ie),(/ v%nx,v%ny /))
    end do
    
    return

  end subroutine sv_to_field_3d



  function get_state_time() result (time)
    type(time_type) :: time
    
    if ( .not. module_initialized ) call static_init_model
    time=cosmo_fc_time

    return

  end function get_state_time



  function get_state_vector() result (sv)

    real(r8)             :: sv(1:model_size)

    integer              :: islab,ikind,nx,ny,sidx,eidx
    real(r8),allocatable :: mydata(:,:)

    if ( .not. module_initialized ) call static_init_model

    call set_allowed_state_vector_vars()

    do islab=1,nslabs
      ikind=cosmo_slabs(islab)%dart_kind
      if (ikind>0) then
        if (is_allowed_state_vector_var(ikind)) then
          nx=state_vector_vars(ikind)%nx
          ny=state_vector_vars(ikind)%ny
          allocate(mydata(1:nx,1:ny))
          mydata=get_data_from_binary(cosmo_filename,grib_header(islab),nx,ny)
          sidx=cosmo_slabs(islab)%dart_sindex
          eidx=cosmo_slabs(islab)%dart_eindex
          state_vector(sidx:eidx)=reshape(mydata,(/ (nx*ny) /))
          deallocate(mydata)
        end if
      end if
    end do

    sv(:)=state_vector(:)

    return

  end function get_state_vector



  subroutine write_grib_file(sv,nfile)

    real(r8),intent(in)           :: sv(:)
    character(len=128),intent(in) :: nfile

    integer                       :: istat = 0
    integer                       :: islab,ipos,lpos,bpos,griblen
    integer                       :: mylen,hlen,ix,iy,nx,ny,idx
    integer                       :: dval,ibsf,idsf
    real(r8)                      :: bsf,dsf
    integer(kind=1)               :: bin4(4)
    integer(kind=1),allocatable   :: bytearr(:),tmparr(:)
    real(r8),allocatable          :: mydata(:,:)
    real(r8)                      :: ref_value
    integer                       :: gribunitin, gribunitout,irec
    integer                       :: recpos(nslabs+1),bytpos(nslabs+1),myrlen,myblen

    logical                       :: write_from_sv = .false.

    if ( .not. module_initialized ) call static_init_model

    mylen=0

    ! read information on record and byte positions in GRIB file
    myrlen=(grib_header(2)%start_record-grib_header(1)%start_record)
    myblen=(grib_header(2)%start_record-grib_header(1)%start_record)*4+grib_header(2)%data_offset
    DO islab=1,nslabs
       recpos(islab)=grib_header(islab)%start_record
       bytpos(islab)=(grib_header(islab)%start_record-1)*4+1+grib_header(islab)%data_offset+4
    END DO
    recpos(nslabs+1)=recpos(nslabs)+myrlen
    bytpos(nslabs+1)=bytpos(nslabs)+myblen

    ! generate byte_array to write
    ALLOCATE(bytearr(1:bytpos(nslabs+1)+100))
    bytearr(:)=0

    ! read all data from the input GRIB file
    gribunitin  = get_unit()
    OPEN(gribunitin,FILE=TRIM(cosmo_filename),FORM='UNFORMATTED',ACCESS='DIRECT',RECL=record_length)
    irec=1
    lpos=1

    do while (istat==0)
      read(gribunitin,rec=irec,iostat=istat) bin4
      irec=irec+1
      lpos=lpos+4
      bytearr(lpos:lpos+3)=bin4
    end do
    griblen=lpos-1

    call close_file(gribunitin)

    ipos=bytpos(1)

    if (verbose) write(*,*)'number of slabs is ',nslabs

    do islab=1,nslabs

      if (verbose) write(*,'(A8,A,I4,A,I4,A,I12)')cosmo_slabs(islab)%varname_short," is GRIB record ",islab," of ",nslabs,", byte position is ",ipos

      nx=cosmo_slabs(islab)%dims(1)
      ny=cosmo_slabs(islab)%dims(2)

      idx=cosmo_slabs(islab)%dart_sindex

      ! check if variable is in state vector
      
      write_from_sv = .false.
      if (idx >= 0) write_from_sv = .true.

      if (write_from_sv) then

        ! if variable is in state vector

        if (verbose) write(*,*)'         ... data is written to GRIB file from state vector'

        if (verbose) write(*,'(A8,A,I4,A,I4,A,2(1x,I12))')cosmo_slabs(islab)%varname_short," is GRIB record ",islab," of ",nslabs,", i1/i2 are ",idx,(idx+nx*ny-1)

        allocate(mydata(1:nx,1:ny))
        mydata(:,:)=reshape(sv(idx:(idx+nx*ny-1)),(/ nx,ny /))

        if (cosmo_slabs(islab)%dart_kind==KIND_PRESSURE_PERTURBATION) then
          mydata(:,:)=mydata(:,:)-non_state_data%p0fl(:,:,cosmo_slabs(islab)%ilevel)
        end if

        ref_value=minval(mydata)
        bin4(1:4)=from_float1(ref_value,cosmo_slabs(islab)%ref_value_char)

        hlen=size(grib_header(islab)%pds)+size(grib_header(islab)%gds)

        ! offset to binary data section in GRIB data, 8 because of indicator section
        bpos=ipos+hlen+8

        ! set new reference value in GRIB data
        bytearr(bpos+6:bpos+9)=bin4

        ! get the binary scale factor
        CALL byte_to_word_signed(bytearr(bpos+4:bpos+5),ibsf,2)
        bsf=FLOAT(ibsf)
        
        ! get the decimal scale factor
        CALL byte_to_word_signed(grib_header(islab)%pds(27:28),idsf,2)
        dsf=FLOAT(idsf)
        
        ! allocate a temporal array to save the binary data values
        allocate(tmparr((nx*ny*2)))
        lpos=1
        DO iy=1,ny
        DO ix=1,nx
          dval=int((mydata(ix,iy)-ref_value)*((10.**dsf)/(2.**bsf)))
          tmparr(lpos:lpos+1)=word_to_byte(dval)
          lpos=lpos+2
        END DO
        END DO

        deallocate(mydata)

        ! overwrite the old with new data
        bytearr(bpos+11:(bpos+nx*ny*2))=tmparr(1:(nx*ny*2))
         
        deallocate(tmparr)
       
        ipos=bytpos(islab+1)

      else

        if (verbose) write(*,*)'         ... data is copied from old grib file'
        
        ! if variable is not in state vector then skip this slab
        
        ipos=bytpos(islab+1)
                
      end if

    END DO
    
    ! write the new GRIB file
    gribunitout = get_unit()
    OPEN(gribunitout,FILE=TRIM(nfile),FORM='UNFORMATTED',ACCESS='stream')
    WRITE(gribunitout) bytearr(5:griblen)

    call close_file(gribunitout)

    return

  end subroutine write_grib_file


  function get_cosmo_filename()
    character(len=256) :: get_cosmo_filename
    character(len=256) :: lj_filename
    
    lj_filename        = adjustl(cosmo_filename)
    get_cosmo_filename = trim(lj_filename)
    
  end function get_cosmo_filename
  
  
  subroutine write_state_times(iunit, statetime, advancetime)
    integer,         intent(in) :: iunit
    type(time_type), intent(in) :: statetime, advancetime
    
    character(len=32) :: timestring 
    integer           :: iyear, imonth, iday, ihour, imin, isec
    integer           :: ndays, nhours, nmins, nsecs
    type(time_type)   :: interval
    
    call get_date(statetime, iyear, imonth, iday, ihour, imin, isec)
    write(timestring, "(I4,5(1X,I2))") iyear, imonth, iday, ihour, imin, isec
    write(iunit, "(A)") trim(timestring)
    
    call get_date(advancetime, iyear, imonth, iday, ihour, imin, isec)
    write(timestring, "(I4,5(1X,I2))") iyear, imonth, iday, ihour, imin, isec
    write(iunit, "(A)") trim(timestring)
    
    interval = advancetime - statetime
    call get_time(interval, nsecs, ndays)
    nhours = nsecs / (60*60)
    nsecs  = nsecs - (nhours * 60*60)
    nmins  = nsecs / 60
    nsecs  = nsecs - (nmins * 60)
    
    write(timestring, "(I4,3(1X,I2))") ndays, nhours, nmins, nsecs
    write(iunit, "(A)") trim(timestring)
    
  end subroutine write_state_times
  
  
  
end module model_mod
