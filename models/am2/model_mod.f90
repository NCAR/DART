! DART software - Copyright UCAR. This open source software is provided
! by UCAR, "as is", without charge, subject to all terms of use at
! http://www.image.ucar.edu/DAReS/DART/DART_download
!
! $Id$

module model_mod

!----------------------------------------------------------------------
! purpose: interface between AM2 and DART
!              Translate to/from state_vector and restart file(s)
!              Initialize model
!              Generate expected obs from model state (model_interpolate)
!              Find state variables (or obs) that are close to a given base observation.
!                  (get_close_obs)
! author: Robert Pincus, CIRES/NOAA ESRL PSD1
!==============================================================================================

  use netcdf
  use types_mod,         only : r8, MISSING_R8
  use time_manager_mod,  only : time_type, set_time, print_time, set_calendar_type, GREGORIAN
  use utilities_mod,     only : open_file, close_file, find_namelist_in_file, check_namelist_read, &
                                register_module, error_handler, file_exist, E_ERR, E_WARN, E_MSG,  &
                                nmlfileunit, do_output, do_nml_file, do_nml_term, logfileunit
  use mpi_utilities_mod, only : my_task_id, task_count
  use location_mod,      only : location_type,      get_close_maxdist_init,      &
                                get_close_obs_init, get_close_obs, set_location, &
                                get_location, vert_is_level, vert_is_pressure, vert_is_height
  use obs_kind_mod,      only : QTY_U_WIND_COMPONENT, QTY_V_WIND_COMPONENT,  &
                                QTY_SURFACE_PRESSURE, QTY_TEMPERATURE,       &
                                QTY_SPECIFIC_HUMIDITY, QTY_PRESSURE,         &
                                QTY_CLOUD_LIQUID_WATER, QTY_CLOUD_ICE,       &
                                get_index_for_type_of_obs, get_quantity_for_type_of_obs
                                ! We'll need to add a kind_cloud_fraction to correspond to AM2 prog var
  use location_mod,      only:  VERTISSURFACE, VERTISLEVEL
  use    netcdf_utilities_mod, only : nc_check
  implicit none
  private

  !
  ! Standard DART interface
  !
  public :: model_type,             &
            get_model_size,         &
            adv_1step,              &
            get_state_meta_data,    &
            model_interpolate,      &
            get_model_time_step,    &
            end_model,              &
            static_init_model,      &
            init_time,              &
            init_conditions,        &
            nc_write_model_atts,    &
            nc_write_model_vars,    &
            pert_model_state,       &
            get_close_maxdist_init, &
            get_close_obs_init,     &
            get_close_obs,          &
            ens_mean_for_model
  !
  ! Procedures from CAM model_mod used in CAM utilities
  !
  public :: prog_var_to_vector,     &
            vector_to_prog_var,     &
            init_model_instance,    &
            end_model_instance,     &
            read_model_init,        &
            write_model_init

! version controlled file description for error handling, do not edit
character(len=256), parameter :: source   = &
   "$URL$"
character(len=32 ), parameter :: revision = "$Revision$"
character(len=128), parameter :: revdate  = "$Date$"

  !==============================================================================================
  !
  ! Global declarations
  !

  type model_type
    !
    ! "Prognostic variables"
    !   Surface pressure may be computed from pressure thickness (delp) on the fly
    !
    real, pointer, dimension(:, :)       :: ps => null()
    real, pointer, dimension(:, :, :)    :: u => null(), v => null(), T => null(), delp => null()
    ! Tracers; could generalize to 2D and 3D collections
    real, pointer, dimension(:, :, :, :) :: tracers => null()
    real, pointer, dimension(:)          :: tracerTypes  => null() ! Dart obs type associated with each tracer
  end type model_type

  !----------------------------------------------------------------------
  ! Namelist variables with default values follow
  !----------------------------------------------------------------------

  integer, parameter :: maxnum_tracers = 10
  character(len = nf90_max_name), dimension(maxnum_tracers) :: &
    tracer_names = '', tracer_files = '', tracer_obs_kind_names = "", tracer_config_files = ''
  integer, dimension(maxnum_tracers) :: tracer_obs_kinds = -1

  ! output_state_vector = .true.     results in a "state-vector" netCDF file
  ! output_state_vector = .false.    results in a "prognostic-var" netCDF file
  logical :: output_state_vector = .false.

  ! File where basic info about model configuration can be found
  character(len = nf90_max_name) :: model_config_file  = 'fv_rst.res.nc', &
                                    model_restart_file = 'fv_rst.res.nc'

  ! Define location restrictions on which observations are assimilated
  ! (values are calculated anyway, but istatus is set to 2)
  real(r8) :: max_obs_lat_degree        = 90.0_r8
  real(r8) :: highest_obs_pressure_mb   = 150.0_r8
  real(r8) :: highest_state_pressure_mb = 150.0_r8

  ! Specify shortest time step that the model will support
  ! This may be limited below by the model itself
  integer :: Time_step_seconds = 21600, Time_step_days = 0

  namelist /model_nml/                                                      &
    tracer_names, tracer_files, tracer_config_files,                        &
    tracer_obs_kind_names, output_state_vector,                             &
    model_config_file, model_restart_file,                                  &
    highest_obs_pressure_mb, highest_state_pressure_mb, max_obs_lat_degree, &
    Time_step_seconds, Time_step_days

  !----------------------------------------------------------------------
  !
  ! Useful global storage
  !

  character(len=256) :: string1
  logical, save :: module_initialized = .false.
  type(time_type) :: Time_step_atmos

  !
  ! Model top pressure is used for converting pressure thickness delp to surface pressure
  !   It should ideally be availible in phalf but the restart files I have don't have
  !   valid values in that dimension variable.
  !   The dimension variables have units of mb but delp has units of Pa.
  !
  real(kind = r8), parameter :: mb_to_pa = 100., model_top_pressure = 1. * mb_to_pa
  real(kind = r8), parameter :: gravity = 9.81, Rdgas = 287. ! MKS units

  integer :: num_tracers, num_lats, num_lons, num_levels, model_size

  !
  ! We assume here that all 2D prognostic variables are functions of lat, lon
  !
  integer, parameter :: num_3d_prog_vars = 3, num_2d_prog_vars = 1
  character(len = 1), dimension(num_3d_prog_vars), parameter :: &
         names_3d_prog_vars = (/ "U", "V", "T" /)
  character(len = 2), dimension(num_2d_prog_vars), parameter :: &
         names_2d_prog_vars = (/ "PS" /)
  integer, dimension(num_3d_prog_vars),            parameter :: &
         kinds_3d_prog_vars = (/ QTY_U_WIND_COMPONENT, QTY_V_WIND_COMPONENT, QTY_TEMPERATURE /)
  integer, dimension(num_2d_prog_vars),            parameter ::  &
         kinds_2d_prog_vars = (/ QTY_SURFACE_PRESSURE /)

  !
  ! Information related to dimensions (coordinates)
  !
  integer, parameter :: num_dims = 6
  character(len = 5), dimension(num_dims), parameter :: &
     dim_names = (/ 'lat  ', 'latu ', 'lon  ', 'lonv ', 'pfull', 'phalf' /)
  !
  ! Dimension ids and dimension variable ids for each of the dimensions in the model_config_file
  !
  integer, dimension(num_dims) :: dim_lens, dim_var_ids
  real, allocatable, dimension(:) :: lat, latu, lon, lonv, pfull, phalf

  !
  ! Surface geopotential and ak, bk terms, read from model_config_file,
  !   used for pressure computations
  !
  real, dimension(:, :), allocatable :: surface_geopotential
  real, dimension(:),    allocatable :: ak, bk, akmid, bkmid

  !-----------------------------------------------------------------------

contains

  ! ----------------------------------------------------------------------------
  !
  !  Public procedures
  !
  ! ----------------------------------------------------------------------------

  subroutine static_init_model()
    !
    ! Initializes class data for  model
    !

    ! Local variables
    !
    integer :: iunit, io, ens_member, num_tasks, my_task
    integer :: ncfileid, ncvarid, i, index
    logical :: do_out
    ! --------------------------------------------------------------------------

    if ( module_initialized ) return ! only need to do this once.

    ! Since this routine calls other routines that could call this routine
    ! we'll say we've been initialized pretty dang early.
    module_initialized = .true.

    call register_module(source, revision, revdate)

    ! Calendar information is not passed to the model; it must be set in the model namelist
    call set_calendar_type(GREGORIAN)

    call find_namelist_in_file("input.nml", "model_nml", iunit)
    read(iunit, nml = model_nml, iostat = io)
    call check_namelist_read(iunit, io, "model_nml")

    ! set the printed output logical variable to reduce printed output;
    ! depends on whether this is being called by trans_... (read ens member # from file 'element' )
    ! or by filter (multiple processes, printout controlled by do_output())
    if (file_exist('element')) then
      ! debug; fix this ugliness
      !
      ! Inherited from CAM model_mod; "element" exists for all processes save ???
      !
       open(unit = 99, file='element', form = 'formatted')
       read(99,*) ens_member
       close(99)
       do_out = .false.
       if (ens_member == 1) do_out = .true.
    else
       do_out = do_output()
       ! static_init_model is called once(?) for each task(?).
       ! There may be more or fewer ensemble members than tasks.
       ! No problem if there are fewer.
       ! In pert_model_state generate a unique ens_member from my_task and globally stored info
       !    about previous calls to pert_model_state.
       num_tasks = task_count()
       my_task = my_task_id()
    end if

    ! Record the namelist values
    if (do_nml_file()) write(nmlfileunit, nml = model_nml)
    if (do_nml_term()) write(    *      , nml = model_nml)

    ! Set the model minimum time step from the namelist seconds and days input
    Time_step_atmos = set_time(Time_step_seconds, Time_step_days)
    if (do_out) call print_time(Time_step_atmos)

    !
    ! Open config file; read number of and values of our six dimensions
    !   Checking the tracers would require having a "config" file for each of the tracer fields
    !
    call nc_check(nf90_open(path = trim(model_config_file), mode = nf90_nowrite, ncid = ncfileid), &
                 'static_init_model', 'opening '// trim(model_config_file))
    call read_dimension_info(ncFileId)

    !
    ! Read in ak, bk, surface geopotential
    !   The latter is a function of time but there's only one value of time
    !   We could check the dimension order for an extra warm-and-fuzzy feeling
    !
    allocate(surface_geopotential(num_lons, num_lats), ak(num_levels + 1), bk(num_levels + 1), &
             akmid(num_levels), bkmid(num_levels))

    call nc_check(nf90_inq_varid(ncfileid, 'ak', ncvarid), 'static_init_model', 'looking for varid of ak')
    call nc_check(nf90_get_var(ncfileid, ncvarid, ak),     'static_init_model', 'reading ak')

    call nc_check(nf90_inq_varid(ncfileid, 'bk', ncvarid), 'static_init_model', 'looking for varid of bk')
    call nc_check(nf90_get_var(ncfileid, ncvarid, bk),     'static_init_model', 'reading bk')

    call nc_check(nf90_inq_varid(ncfileid,  'Surface_geopotential', ncvarid), &
                 'static_init_model', 'looking for varid of Surface_geopotential')
    call nc_check(nf90_get_var(ncfileid, ncvarid, surface_geopotential), &
                 'static_init_model', 'reading Surface_geopotential')

    !
    ! Compute akmid, bkmid for computing pressure at layers
    !
    akmid(:num_levels) = ak(:num_levels) + (ak(2:) - ak(:num_levels))/2.
    bkmid(:num_levels) = bk(:num_levels) + (bk(2:) - bk(:num_levels))/2.

    !
    ! Count the number of tracer names provided; check to be sure each tracer name
    !   has a corresponding file to read the tracer from
    !
    num_tracers = count(len_trim(tracer_names) > 0)
    if(count(len_trim(tracer_files) > 0) /= num_tracers) &
      call error_handler(E_ERR, 'static_model_init', &
           "Different number of model tracer names, tracer files provided", source, revision, revdate)
    if(count(len_trim(tracer_obs_kind_names) > 0) /= num_tracers) &
      call error_handler(E_ERR, 'static_model_init', &
           "Different number of model tracer names, tracer kinds provided", source, revision, revdate)
    do i = 1, num_tracers
      index = get_index_for_type_of_obs(tracer_obs_kind_names(i))
      if(index > 0) then
        tracer_obs_kinds(i) = get_quantity_for_type_of_obs(index)
      else
        call error_handler(E_ERR, 'static_model_init', &
           "Tracer type " // trim(tracer_obs_kind_names(i)) // " unknown" , source, revision, revdate)
      end if
    end do

    model_size = num_2d_prog_vars * (num_lats * num_lons)              + &
                 num_3d_prog_vars * (num_lats * num_lons * num_levels) + &
                 num_tracers      * (num_lats * num_lons * num_levels)

    call nc_check(nf90_close(ncfileid), &
                 'static_init_model', 'closing '// trim(model_config_file))
  end subroutine static_init_model

  ! ----------------------------------------------------------------------------

  integer function get_model_size()

     if ( .not. module_initialized ) call static_init_model

     get_model_size = model_size

  end function get_model_size

  ! ----------------------------------------------------------------------------

  subroutine get_state_meta_data(index_in, location, var_type)
    integer,             intent(in)  :: index_in
    type(location_type), intent(out) :: location
    integer,             intent(out), optional :: var_type
    !------------------------------------------------------------------
    !
    ! Given an integer index into the state vector structure, returns the
    ! associated location and DART observation kind

    integer  :: i, start, finish, field_number, local_type, &
                local_index, lat_index, lon_index, level_index, &
                which_vert
    real(r8) :: local_lat, local_lon, vert_loc

    if ( .not. module_initialized ) call static_init_model

    start        = 1
    field_number = 0
    local_index  = 0

    ! Figure out which field the index points to and the index within that field
    !   Walk through each field, compute its start and end point, and see if index_in
    !   lies between them.
    !   Assign the observation type while we're at it.

    do i = 1, num_2d_prog_vars
      finish = start + (num_lons * num_lats) - 1
      if(index_in >= start .and. index_in <= finish) then
         field_number = i
         local_type = kinds_2d_prog_vars(i)
         local_index = index_in - start + 1
         exit
      end if
      start = finish + 1
    end do

    if(field_number <= 0) then
      do i = 1, num_3d_prog_vars
        finish = start + (num_levels * num_lons * num_lats) - 1
        if(index_in >= start .and. index_in <= finish) then
          field_number = num_2d_prog_vars + i
          local_type = kinds_3d_prog_vars(i)
          local_index = index_in - start + 1
          exit
        end if
        start = finish + 1
      end do
    end if

    if(field_number <= 0) then
      do i = 1, num_tracers
        finish = start + (num_levels * num_lons * num_lats) - 1
        if(index_in >= start .and. index_in <= finish) then
          field_number = num_2d_prog_vars + num_3d_prog_vars + i
          local_type = tracer_obs_kinds(i)
          local_index = index_in - start + 1
          exit
        end if
        start = finish + 1
      end do
    end if

    if(field_number <=0) &
      call error_handler(E_ERR,'get_state_meta_data', "State vector index out of bounds", &
                         source, revision, revdate)

    !
    ! Compute the x, y, z index. The state vector is stored in order level, lon, lat
    !   (or lon, lat for 2D fields)
    !
    if(field_number <= num_2d_prog_vars) then
      lat_index = (local_index-1)/num_lons + 1
      lon_index = mod(local_index-1, num_lons) + 1
      level_index = 0
      which_vert = VERTISSURFACE
    else
      lat_index = (local_index-1)/(num_lons * num_levels) + 1
      local_index = mod(local_index-1, num_lons * num_levels) + 1
      lon_index = (local_index-1)/num_levels + 1
      level_index = mod(local_index-1, num_levels) + 1
      which_vert = VERTISLEVEL
    end if
    if(lat_index < 0 .or. lat_index > num_lats) &
      call error_handler(E_ERR,'get_state_meta_data', "Lat calculation out of bounds", source, revision, revdate)
    if(lon_index < 0 .or. lon_index > num_lons) &
      call error_handler(E_ERR,'get_state_meta_data', "Lon calculation out of bounds", source, revision, revdate)
    !
    ! Convert indexes to lat, lon, level values
    !
    local_lat = lat(lat_index); if(local_type == QTY_U_WIND_COMPONENT) local_lat = latu(lat_index)
    local_lon = lon(lon_index); if(local_type == QTY_V_WIND_COMPONENT) local_lon = lonv(lon_index)

    location = set_location(local_lon, local_lat, real(level_index, r8),  which_vert)
    if(present(var_type)) var_type = local_type
  end subroutine get_state_meta_data

  ! ----------------------------------------------------------------------------

  function nc_write_model_atts( ncFileID ) result (ierr)
    ! Writes model-specific attributes to a netCDF file.
    ! TJH Fri Aug 29 MDT 2003
    ! Modified for AM2 by Robert Pincus, Dec. 2007
    !
    integer, intent(in)  :: ncFileID      ! netCDF file identifier
    integer              :: ierr          ! return value of function
    !-----------------------------------------------------------------------------------------
    integer :: unlimitedDimID
    integer :: MemberDimID, StateVarDimID, TimeDimID, & ! ScalarDimID, &
               LatDimId, LatuDimId, LonDimId, LonvDimId, pfullDimId, phalfDimId
    integer :: StateVarID, StateVarVarID
    integer :: configFileId, configVarId, diagVarId, dim
    integer :: i, tracer, trcrFileId
    character(len=129)    :: errstring
    character(len=8)      :: crdate      ! needed by F90 DATE_AND_TIME intrinsic
    character(len=10)     :: crtime      ! needed by F90 DATE_AND_TIME intrinsic
    character(len=5)      :: crzone      ! needed by F90 DATE_AND_TIME intrinsic
    integer, dimension(8) :: values      ! needed by F90 DATE_AND_TIME intrinsic
    character(len=NF90_MAX_NAME) :: str1
    !-------------------------------------------------------------------------------

    if ( .not. module_initialized ) call static_init_model

    ierr = -1     ! assume it's not going to work

    !
    ! Make sure ncFileID refers to an open netCDF file, and then put into define mode.
    !    More dimensions, variables and attributes will be added in this routine.
    !
    write(errstring,*) 'ncFileID is', ncFileID
    call nc_check(nf90_Inquire(ncFileID, unlimitedDimID = unlimitedDimID), &
                  'nc_write_model_atts', 'Inquire '// trim(errstring))
    call nc_check(nf90_Redef(ncFileID), 'nc_write_model_atts', 'Redef '// trim(errstring))

    !
    ! Write DART rev info and current date & time as global attributes
    !
    call DATE_AND_TIME(crdate,crtime,crzone,values)
    write(str1,'("YYYY MM DD HH MM SS = ",i4, 5(1x,i2.2))') &
                      values(1), values(2), values(3), values(5), values(6), values(7)

    call nc_check(nf90_put_att(ncFileID, NF90_GLOBAL, "creation_date",str1),        &
                  'nc_write_model_atts', 'put_att creation_date'//trim(str1))
    call nc_check(nf90_put_att(ncFileID, NF90_GLOBAL, "model_revision",revision),   &
                  'nc_write_model_atts', 'put_att model_revision'//trim(revision))
    call nc_check(nf90_put_att(ncFileID, NF90_GLOBAL, "model_revdate",revdate),     &
                  'nc_write_model_atts', 'put_att model_revdate'//trim(revdate))
    call nc_check(nf90_put_att(ncFileID, NF90_GLOBAL, "model","AM2"),               &
                  'nc_write_model_atts','put_att model AM2')

    !
    ! Dimension ids for existing dimensions so we can define variables properly
    !
    call nc_check(nf90_inq_dimid(ncid = ncFileID, name = "copy", dimid = MemberDimID), &
                  'nc_write_model_atts', 'inq_dimid copy')
    call nc_check(nf90_inq_dimid(ncid = ncFileID, name = "time", dimid =   TimeDimID), &
                  'nc_write_model_atts', 'inq_dimid time')
    if ( TimeDimID /= unlimitedDimId ) then
      write(errstring, *)'Time dimension ID ', TimeDimID, 'must match Unlimited Dimension ID ', unlimitedDimId
      call error_handler(E_ERR,'nc_write_model_atts', errstring, source, revision, revdate)
    end if

    !
    ! Define the new dimensions IDs
    !

    if ( output_state_vector ) then
      !
      ! Create a dimension that corresponds to the length of the state vector
      !
      call nc_check(nf90_def_dim(ncid=ncFileID, name="StateVariable",  &
                              len=model_size, dimid = StateVarDimID),  &
                    'nc_write_model_atts', 'def_dim StateVariable')
      !
      ! Define the state vector coordinate variable and its attributes
      !
      call nc_check(nf90_def_var(ncid=ncFileID, name="StateVariable", xtype=nf90_int,          &
                 dimids=StateVarDimID, varid=StateVarVarID),                                   &
                    'nc_write_model_atts','def_var  state vector')
      call nc_check(nf90_put_att(ncFileID, StateVarVarID, "long_name", "State Variable ID"),   &
                    'nc_write_model_atts','put_att long_name state vector ')
      call nc_check(nf90_put_att(ncFileID, StateVarVarID, "units",     "indexical"),           &
                    'nc_write_model_atts','put_att units state vector ' )
      call nc_check(nf90_put_att(ncFileID, StateVarVarID, "valid_range", (/ 1, model_size /)), &
                    'nc_write_model_atts','put_att valid range state vector ')
    else
      !
      ! We'll copy the dimension definitions (including all the attributes) for the num_dims
      !   dimensions we're using from the config file into the diagnostics file.
      ! Variable dim_names holds the dimension names
      ! We could also do this by looping over all the dimensions in the config file, skipping time
      !
      call nc_check(nf90_open(path = trim(model_config_file), mode = nf90_nowrite, ncid = configFileId), &
                   'nc_write_model_atts', 'opening '// trim(model_config_file) )
      do dim = 1, num_dims
        call copy_dim_var_pair(ncFileId, trim(dim_names(dim)), dim_lens(dim), configFileId, 'nc_write_model_atts')
      end do
    end if

    !
    ! Create the variables and their attributes
    !
    if ( output_state_vector ) then
      !
      ! Define the state vector and its attributes
      !
      call nc_check(nf90_def_var(ncid=ncFileID, name="state", xtype=nf90_real,                 &
                 dimids = (/ StateVarDimID, MemberDimID, unlimitedDimID /), varid=StateVarID), &
                    'nc_write_model_atts','def_var state vector')
      call nc_check(nf90_put_att(ncFileID, StateVarID, "long_name", "model state or fcopy"),   &
                    'nc_write_model_atts','put_att long_name model state or fcopy ')
      call nc_check(nf90_put_att(ncFileID, StateVarId, "vector_to_prog_var","AM2"),            &
                    'nc_write_model_atts','put_att vector_to_prog_var AM2 ')
    else
      !
      ! We'll need all the dimension ids except phalf
      !
      call nc_check(nf90_inq_dimid(ncFileId, 'lon',   LonDimId),   'nc_write_model_atts', 'finding dimid for lon')
      call nc_check(nf90_inq_dimid(ncFileId, 'lonv',  LonvDimId),  'nc_write_model_atts', 'finding dimid for lonv')
      call nc_check(nf90_inq_dimid(ncFileId, 'lat',   LatDimId),   'nc_write_model_atts', 'finding dimid for lat')
      call nc_check(nf90_inq_dimid(ncFileId, 'latu',  LatuDimId),  'nc_write_model_atts', 'finding dimid for latu')
      call nc_check(nf90_inq_dimid(ncFileId, 'pfull', pfullDimId), 'nc_write_model_atts', 'finding dimid for pfull')

      !
      ! We'l write out surface pressure regardless of whether the model is using delp or ps internally
      !    Then there are three more prognostic variables (u, v, t) and tracers
      !    Coordinate order in the diagnostic file is [lev], lon, lat, copy, time
      !
      ! -------------------
      !
      ! Surface pressure - copy units from delp field, write long name, and ignore any other attributes
      !
      call nc_check(nf90_def_var(ncFileID, "PS", nf90_real,              &
                    (/ lonDimId, latDimId, MemberDimID, TimeDimId /),  diagVarId), &
                    'nc_write_model_atts','def_var PS')
      call nc_check(nf90_inq_varid(configFileId, "DELP", configVarId),  &
                    'nc_write_model_atts', 'getting varid from config file for DELP' )
      call nc_check(nf90_copy_att(configFileId, configVarId, "units", ncFileId, diagVarId),  &
                    'nc_write_model_atts', 'copying units attribute for DELP')
      call nc_check(nf90_put_att(ncFileId, diagVarId, "long_name", "Surface pressure"),  &
                    'nc_write_model_atts', 'writing long_name attribute for DELP' )
      ! -------------------
      !
      ! The other prognostic variables (T, U, V)
      !
      call define_3d_real_var(ncFileId, configFileId, &
                              "U", (/ pfullDimId, lonDimId, latuDimId, MemberDimID, TimeDimId /), &
                              'nc_write_model_atts')
      call define_3d_real_var(ncFileId, configFileId, &
                              "V", (/ pfullDimId, lonvDimId, latDimId, MemberDimID, TimeDimId /), &
                              'nc_write_model_atts')
      call define_3d_real_var(ncFileId, configFileId, &
                              "T", (/ pfullDimId, lonDimId,  latDimId, MemberDimID, TimeDimId /), &
                              'nc_write_model_atts')

          call nc_check(nf90_close(configFileId), &
                   'nc_write_model_atts', 'closing '// trim(model_config_file))
      ! -------------------
      !
      ! Tracers, which we assume are all defined on levels (not interfaces) and
      !   on the regular (not staggered) grid
      !
      do tracer = 1, num_tracers
             call nc_check(nf90_open(path = trim(tracer_config_files(tracer)), mode = nf90_nowrite, &
                  ncid = trcrFileId),'nc_write_model_atts', 'opening '// trim(tracer_config_files(tracer)))

             call define_3d_real_var(ncFileId, trcrFileId,tracer_names(tracer), &
                  (/ pfullDimId, lonDimId, latDimId, MemberDimID, TimeDimId /), 'nc_write_model_atts')

             call nc_check(nf90_close(trcrFileId),'nc_write_model_atts', &
                  'closing' // trim(tracer_config_files(tracer)))
      end do
      ! -------------------

    end if

    call nc_check(nf90_enddef(ncfileID), 'nc_write_model_atts','enddef ')

    !
    ! Fill the coordinate variables
    !
    if ( output_state_vector ) then
      ! Fill the state variable coordinate variable
      call nc_check(nf90_put_var(ncFileID, StateVarVarID, (/ (i,i=1,model_size) /) ), &
                    'nc_write_model_atts','put_var StateVar ')
    else
      call write_1D_values(ncFileId, "pfull", pfull, 'nc_write_model_atts')
      call write_1D_values(ncFileId, "phalf", phalf, 'nc_write_model_atts')
      call write_1D_values(ncFileId, "lat",   lat,   'nc_write_model_atts')
      call write_1D_values(ncFileId, "latu",  latu,  'nc_write_model_atts')
      call write_1D_values(ncFileId, "lon",   lon,   'nc_write_model_atts')
      call write_1D_values(ncFileId, "lonv",  lonv,  'nc_write_model_atts')
    end if

    !
    ! Flush the buffer and leave netCDF file open
    !
    call nc_check(nf90_sync(ncFileID),'nc_write_model_atts', 'sync ')
    write (*, *) 'nc_write_model_atts: netCDF file ',ncFileID,' is synched ...'

    ierr = 0

  end function nc_write_model_atts

  ! ----------------------------------------------------------------------------

  integer function nc_write_model_vars(ncFileID, statevec, copyindex, timeindex ) result (ierr)
    !
    ! Writes the state vector to a diagnostics netcdf file.
    !
    integer,                intent(in) :: ncFileID      ! netCDF file identifier
    real(r8), dimension(:), intent(in) :: statevec
    integer,                intent(in) :: copyindex, timeindex

    integer :: varId, start, finish, i
    !
    ! We'll unpack the state vector into this temporary array before we write it
    !   We rely on the number of lats and lons being the same on the regular and staggered grids
    !   Dimension order is consistent with storage order in the state vector and
    !   the order desired in the diagnostics file (it's not the same as the order in the restart file)
    !
    real,    dimension(num_levels, num_lons, num_lats) :: tempField
    logical, dimension(num_levels, num_lons, num_lats) :: allTrue
    ! --------------------

     if ( .not. module_initialized ) call static_init_model

    allTrue(:, :, :) = .true.
    ierr = -1 ! Assume the worst
    if ( output_state_vector ) then
       call nc_check(NF90_inq_varid(ncFileID, "state", varId), "nc_write_model_vars", "getting varid for state vector"  )
       call nc_check(NF90_put_var(ncFileID, varId, statevec,  start=(/ 1, copyindex, timeindex /)), &
                  "nc_write_model_vars", "writing state vector")
    else
      !
      ! The unpacking code is copied from vector_to_prog_var
      !
      ! 2D field - surface pressure
      !
      start = 1; finish = start + (num_lons * num_lats) - 1
      tempField(1, :, :) = unpack(statevec(start:finish), allTrue(1, :, :), field = 0._r8)
      call nc_check(nf90_inq_varid(ncFileID, "PS", varID), "nc_write_model_vars", "Getting varid for PS")
      call nc_check(nf90_put_var(ncFileID, varID, tempField(1, :, :),start=(/ 1, 1, copyindex, &
timeindex /)), "nc_write_model_vars", "Writing PS")

      !
      ! 3D fields - u, v, and T
      !   These were stored in dim order level, lon, lat in the state vector and need
      !   to be mapped back to lon, lat, level
      !
      do i = 1, num_3d_prog_vars
        start = finish + 1; finish = start + (num_levels * num_lons * num_lats) - 1
        tempField(:, :, :) = unpack(statevec(start:finish), allTrue, field = 0._r8)
        call nc_check(nf90_inq_varid(ncFileID, names_3d_prog_vars(i), varID), &
            "nc_write_model_vars", "Getting varid for " // names_3d_prog_vars(i))
        call nc_check(nf90_put_var(ncFileID, varID, tempField,start=(/1,1,1,copyindex,timeindex/)), &
            "nc_write_model_vars", "Writing " // names_3d_prog_vars(i))
      end do

      !
      ! Tracers
      !   These were stored in dim order level, lon, lat, tracer_num in the state vector and need
      !   to be mapped back to lon, lat, level, tracer_num
      !
      do i = 1, num_tracers
        start = finish + 1; finish =  start + (num_levels * num_lons * num_lats ) - 1
        tempField(:, :, :) = unpack(statevec(start:finish), allTrue, field = 0._r8)
        call nc_check(nf90_inq_varid(ncFileID, tracer_names(i), varID), "nc_write_model_vars", &
             "Getting varid for " // tracer_names(i))
        call nc_check(nf90_put_var(ncFileID, varID, tempField,start=(/1,1,1,copyindex,timeindex/)), &
             "nc_write_model_vars", "Writing " // tracer_names(i))
      end do
    end if
    ierr = 0
  end function nc_write_model_vars

  ! ----------------------------------------------------------------------------
  !
  ! Stubs for unimplemented public procedures
  !
  ! ----------------------------------------------------------------------------

  subroutine adv_1step(x, time)
    real(r8),        intent(inout) :: x(:)
    type(time_type), intent(in)    :: time
    !
    ! Does a single timestep advance of the model. The input value of
    ! the vector x is the starting condition and x is updated to reflect
    ! the changed state after a timestep. The time argument is intent
    ! in and is used for models that need to know the date/time to
    ! compute a timestep, for instance for radiation computations.
    ! This interface is only called if the namelist parameter
    ! async is set to 0 in perfect_model_obs of filter or if the
    ! program integrate_model is to be used to advance the model
    ! state as a separate executable. If one of these options
    ! is not going to be used (the model will only be advanced as
    ! a separate model-specific executable), this can be a
    ! NULL INTERFACE.

     if ( .not. module_initialized ) call static_init_model

     if (do_output()) then
        call print_time(time,'NULL interface adv_1step (no advance) DART time is')
        call print_time(time,'NULL interface adv_1step (no advance) DART time is',logfileunit)
     endif

     write(string1,*)'DART should not be trying to advance AM2'
     call error_handler(E_ERR,'adv_1step',string1,source,revision,revdate)

     x(:) = MISSING_R8   ! just to satisfy compiler

  end subroutine adv_1step

  ! ----------------------------------------------------------------------------

  subroutine init_time(time)
     type(time_type), intent(out) :: time

     ! Companion interface to init_conditions. Returns a time that is somehow
     ! appropriate for starting up a long integration of the model.
     ! At present, this is only used if the namelist parameter
     ! start_from_restart is set to .false. in the program perfect_model_obs.

     if ( .not. module_initialized ) call static_init_model

     ! for now, just set to 0
     time = set_time(0, 0)

  end subroutine init_time

  ! ----------------------------------------------------------------------------

  subroutine init_conditions(x)
     real(r8), intent(out) :: x(:)

    ! Returns a model state vector, x, that is some sort of appropriate
    ! initial condition for starting up a long integration of the model.
    ! At present, this is only used if the namelist parameter
    ! start_from_restart is set to .false. in the program perfect_model_obs.
    if ( .not. module_initialized ) call static_init_model

    write(string1,*)'DART cannot initialize AM2'
    call error_handler(E_ERR,'init_conditions',string1,source,revision,revdate)

     x(:) = MISSING_R8   ! just to satisfy compiler

  end subroutine init_conditions

  ! ----------------------------------------------------------------------------

  subroutine model_interpolate(x, location, itype, obs_val, istatus)
    real(r8),            intent(in) :: x(:)
    type(location_type), intent(in) :: location
    integer,             intent(in) :: itype
    real(r8),           intent(out) :: obs_val
    integer,            intent(out) :: istatus
    !
    ! Given a state vector, a location, and a model state variable type,
    ! interpolates the state variable field to that location and returns
    ! the value in obs_val. The istatus variable should be returned as
    ! 0 unless there is some problem in computing the interpolation in
    ! which case an alternate value should be returned. The itype variable
    ! is a model specific integer that specifies the type of field (for
    ! instance temperature, zonal wind component, etc.).
    !
    ! As per CAM model_mod:
    !
    ! istatus   meaning                  return expected obs?   assimilate?
    ! 0         obs and model are fine;  yes                    yes
    ! 1         fatal problem;           no                     no
    ! 2         exclude valid obs        yes                    no
    ! 3         unfamiliar obs type      no                     no

    integer :: variableStart, lat_index, lon_index_below, lon_index_above, height_index
    integer :: i, j, k
    real(kind = r8)                  :: lon_weight, lat_weight, this_lon
    real(kind = r8), dimension(3)    :: lon_lat_height
    real(kind = r8), dimension(2, 2) :: cornerValues, cornerLogPressures, psurf, weights ! Dimension order lon, lat
    real(kind = r8), dimension(2, 2, num_levels) &
                                     :: cornerTemperatures, cornerHumidities     ! Dimension order lon, lat, level
    real(kind = r8), dimension(2, 2, num_levels + 1) &
                                     :: cornerHeights, cornerLogPressureProfiles ! Dimension order lon, lat, level
    integer        , dimension(2, 2) :: cornerStarts, cornerTempStarts
    real(kind = r8), dimension(:), &
                      allocatable    :: theseLons, theseLats

    ! Notes after talking to Jeff A
    !   We interpolate the KINDS we have: U, V, T, tracers, and PS; we also allow for surface height
    !   **adding in QTY_PRESSURE**
    !   Vertical coordinates can be expressed as level, pressure (most usual), and height
    !   We'll interpolate linearly in lat, lon, and vertical coordinate since the
    !   errors this introduces are small and can be lumped in to "representativeness"
    !   When the vertical coordinate is pressure the vertical indicies bounding the
    !   observation location might be different at te different corners - I think the
    !   easiest way to get around that is to interpolate in the vertical at each of the four
    !   corners first, then interpolate in the horizontal

     if ( .not. module_initialized ) call static_init_model

    ! -------------------------------------------------------------
    !
    ! Which variable are we looking for?
    !   State vector is ordered ps (2D), u, v, t, tracers
    !
    variableStart = -1
    if(itype == QTY_SURFACE_PRESSURE .or. itype == QTY_PRESSURE) then
      variableStart = 1
    else if (itype == QTY_U_WIND_COMPONENT) then
      variableStart = 1 + (num_lons * num_lats)
    else if (itype == QTY_V_WIND_COMPONENT) then
      variableStart = 1 + (num_lons * num_lats) + 1 * (num_lons * num_lats * num_levels)
    else if (itype == QTY_TEMPERATURE) then
      variableStart = 1 + (num_lons * num_lats) + 2 * (num_lons * num_lats * num_levels)
    else
      do i = 1, num_tracers
        if(itype == tracer_obs_kinds(i)) &
          variableStart = 1 + (num_lons) * num_lats + (i - 1 + 3) * (num_lons * num_lats * num_levels)
      end do
    end if

    if(variableStart <= 0) then
      !
      ! We don't know how to interpolate this kind of observation
      !
      istatus = 3
      obs_val = MISSING_R8
    else
      !
      ! Horizontal interpolation
      !
      allocate(theseLons(num_lons), theseLats(num_lats))
      !
      ! U and V are on staggered grids
      !
      if(itype == QTY_V_WIND_COMPONENT) then
        theseLons(:) = lonv(:)
      else
        theseLons(:) = lon(:)
      end if
      if(itype == QTY_U_WIND_COMPONENT) then
        theseLats(:) = (/ latu(:), real(MISSING_R8) /)
      else
        theseLats(:) = lat(:)
      end if

      lon_lat_height = get_location(location)

      !
      ! Choose lat_index such that (lat_index) <= lat
      !   There's no doubt a smart way to interpolate at the poles, but I'm going to punt
      !
      if(lon_lat_height(2) < theseLats(1) .or. lon_lat_height(2) >= theseLats(num_lats)) then
        !
        ! Latitude is out of bounds
        !
        lat_index = -1
        istatus = 3
        obs_val = MISSING_R8
      else
        lat_index = findIndex(lon_lat_height(2), theseLats(:))
        lat_weight = 1._r8 - (lon_lat_height(2) - theselats(lat_index)) / (theselats(lat_index + 1) - theselats(lat_index))
      end if

      !
      ! Choose upper and lower lon indicies. This is more complicated than lat because we have to
      !   enforce periodicity
      !
      this_lon = lon_lat_height(1)
      if(this_lon < theseLons(1) .or. this_lon >= theseLons(num_lons)) then
        !
        ! Longitude is between last and first lon (or out of bounds)
        !
        if(this_lon < 0 .or. this_lon > 360) then
          lon_index_below = -1; lon_index_above = -1
          istatus = 3
          obs_val = MISSING_R8
        else
          lon_index_below = num_lons
          lon_index_above = 1
          if(this_lon <= theseLons(1)) this_lon = this_lon + 360.
          lon_weight = 1._r8 - (this_lon - theseLons(lon_index_below)) / (theseLons(lon_index_above) + 360. - theseLons(lon_index_below))
        end if
      else
        lon_index_below = findIndex(lon_lat_height(1),  theseLons(:))
        lon_index_above = lon_index_below + 1
        lon_weight = 1._r8 - (this_lon - theseLons(lon_index_below)) / (theseLons(lon_index_above) - theseLons(lon_index_below))
      end if

      deallocate(theseLons, theseLats)

      if(lat_index > 0 .and. lon_index_below > 0) then
        weights(1, 1) =          lon_weight  *          lat_weight
        weights(2, 1) = (1._r8 - lon_weight) *          lat_weight
        weights(1, 2) =          lon_weight  * (1._r8 - lat_weight)
        weights(2, 2) = (1._r8 - lon_weight) * (1._r8 - lat_weight)

        !
        ! We don't always need psurf for the interpolation but it's cheap enough to compute
        !
        psurf(1, 1) = x(1 + (lat_index - 1) * num_lons + lon_index_below - 1)
        psurf(2, 1) = x(1 + (lat_index    ) * num_lons + lon_index_below - 1)
        psurf(1, 2) = x(1 + (lat_index - 1) * num_lons + lon_index_above - 1)
        psurf(2, 2) = x(1 + (lat_index    ) * num_lons + lon_index_above - 1)

        !
        ! Find the four values that bracket the observation location in the horizontal
        !   Storage order is level, lon, lat
        !
        if(itype == QTY_SURFACE_PRESSURE) then
          cornerValues(:, :) = psurf(:, :)   
        else
          !
          ! Find the four starting locations of the 4 columns that bracket the observation location
          !   and interpolate within each column
          !
          cornerStarts(1, 1) = (variableStart + (lat_index - 1) * (num_lons * num_levels) + (lon_index_below - 1) * num_levels)
          cornerStarts(2, 1) = (variableStart + (lat_index - 1) * (num_lons * num_levels) + (lon_index_above - 1) * num_levels)
          cornerStarts(1, 2) = (variableStart + (lat_index    ) * (num_lons * num_levels) + (lon_index_below - 1) * num_levels)
          cornerStarts(2, 2) = (variableStart + (lat_index    ) * (num_lons * num_levels) + (lon_index_above - 1) * num_levels)

          !
          ! We know how to interpolate vertically in level, pressure, and height - are we missing any possibilities?
          !
          if(vert_is_level(location)) then
            if(itype == QTY_PRESSURE) then
               forall(i = 1:2, j = 1:2)
                  cornerValues(i, j) = interpolate1D(desiredLocation = lon_lat_height(3),                 &
                                                   values = akmid(:) + bkmid(:) * psurf(i,j),             &
                                                   locations = (/ (real(k, kind = r8), k = 1, num_levels) /) )
               end forall
            else 
               forall(i = 1:2, j = 1:2)
                  cornerValues(i, j) = interpolate1D(desiredLocation = lon_lat_height(3),                                &
                                                   values = x(cornerStarts(i, j):cornerStarts(i, j) + num_levels - 1), &
                                                   locations = (/ (real(k, kind = r8), k = 1, num_levels) /) )
               end forall
            end if
          else if (vert_is_pressure(location)) then
            !
            ! Better be sure this is in the right units
            !
            forall(i = 1:2, j = 1:2)
               cornerValues(i, j) = interpolate1D(desiredLocation = lon_lat_height(3),                                &
                                                  values = x(cornerStarts(i, j):cornerStarts(i, j) + num_levels - 1), &
                                                  ! Compute pressure at layer midpoints in this column on the fly
                                                  locations = akmid(:) + bkmid(:) * psurf(i, j) )

            end forall
          else if (vert_is_height(location)) then
            !
            ! First compute the height at each interface pressure according to the hypsometric equation
            !   For the moment we use dry air temp but need to use virtual temperature
            !
            !
            ! Temperature profile in each column surrounding our location
            !
            variableStart = 1 + (num_lons * num_lats) + 2 * (num_lons * num_lats * num_levels)
            cornerTempStarts(1, 1) = (variableStart + (lat_index - 1) * (num_lons * num_levels) + (lon_index_below - 1) * num_levels)
            cornerTempStarts(2, 1) = (variableStart + (lat_index - 1) * (num_lons * num_levels) + (lon_index_above - 1) * num_levels)
            cornerTempStarts(1, 2) = (variableStart + (lat_index    ) * (num_lons * num_levels) + (lon_index_below - 1) * num_levels)
            cornerTempStarts(2, 2) = (variableStart + (lat_index    ) * (num_lons * num_levels) + (lon_index_above - 1) * num_levels)
            forall(i = 1:2, j = 1:2) &
              cornerTemperatures(i, j, :) = x(cornerTempStarts(i, j):cornerTempStarts(i, j) + num_levels)

            !
            ! Correct temperature to virtual temperature
            !
            variableStart = -1
            do i = 1, num_tracers
              if(tracer_obs_kinds(i) == QTY_SPECIFIC_HUMIDITY) &
                variableStart = 1 + (num_lons) * num_lats + (i - 1 + 3) * (num_lons * num_lats * num_levels)
            end do
            if(variableStart > 0) then
              cornerTempStarts(1, 1) = (variableStart + (lat_index - 1) * (num_lons * num_levels) + (lon_index_below - 1) * num_levels)
              cornerTempStarts(2, 1) = (variableStart + (lat_index - 1) * (num_lons * num_levels) + (lon_index_above - 1) * num_levels)
              cornerTempStarts(1, 2) = (variableStart + (lat_index    ) * (num_lons * num_levels) + (lon_index_below - 1) * num_levels)
              cornerTempStarts(2, 2) = (variableStart + (lat_index    ) * (num_lons * num_levels) + (lon_index_above - 1) * num_levels)
              forall(i = 1:2, j = 1:2) &
                cornerHumidities(i, j, :) = x(cornerTempStarts(i, j):cornerTempStarts(i, j) + num_levels)
              !
              ! We ignore the distinction between specific humidity and mixing ratio -
              !   this should really be q/(1 - q) instead of q
              !
              cornerTemperatures(:, :, :) = cornerTemperatures(:, :, :) * (1. + 0.61 * (cornerHumidities(:, :, :)/(1.-cornerHumidities(:,:,:))))
            end if
            !
            ! Compute geometric height at each  level, remembering that the layers are ordered top to bottom
            !
            forall(k = 1:num_levels + 1) &
              cornerLogPressureProfiles(:, :, k) = log(ak(k) + bk(k) * psurf(:, :))
            cornerHeights(1, 1, num_levels + 1) = surface_geopotential(lon_index_below, lat_index)
            cornerHeights(2, 1, num_levels + 1) = surface_geopotential(lon_index_above, lat_index)
            cornerHeights(1, 2, num_levels + 1) = surface_geopotential(lon_index_below, lat_index + 1)
            cornerHeights(2, 2, num_levels + 1) = surface_geopotential(lon_index_above, lat_index + 1)
            cornerHeights(:, :, num_levels + 1) = cornerHeights(:, :, num_levels + 1)/gravity
            do k = num_levels, 2, -1
              cornerHeights(:, :, k) = cornerHeights(:, :, k + 1) +                  &
                                       Rdgas/gravity * cornerTemperatures(:, :, k) * &
                                       (cornerLogPressureProfiles(:, :, k + 1) - cornerLogPressureProfiles(:, :, k))
            end do
            
            cornerHeights(:, :, 1) = huge(cornerHeights)
            !
            ! Compute the pressure at the desired height, interpolating linearly in ln(p).
            !
            forall(i = 1:2, j = 1:2)
              cornerLogPressures(i, j) = interpolate1D(desiredLocation = lon_lat_height(3),          &
                                                       values = cornerLogPressureProfiles(i, j,num_levels+1:1:-1),  & 
                                                       locations =  cornerHeights(i, j,num_levels+1:1:-1))
            end forall
            if(itype == QTY_PRESSURE) then
               forall(i = 1:2, j = 1:2)
                  cornerValues(i,j) = exp(cornerLogPressures(i,j))
               end forall
               !
            else
               !
               ! Now interpolate the values as if pressure had been supplied
               !   It might make sense to interpolate in ln(p) but this would be inconsistent with how
               !   we do the interpolation when p is supplied.
               !
               forall(i = 1:2, j = 1:2)
                  cornerValues(i, j) = interpolate1D(desiredLocation = exp(cornerLogPressures(i, j)),                    &
                                                 values = x(cornerStarts(i, j):cornerStarts(i, j) + num_levels - 1), &
                                                 ! Compute pressure at layer midpoints in this column on the fly
                                                 locations = akmid(:) + bkmid(:) * psurf(i, j) )
               end forall
            end if !QTY_PRESSURE loop
            !
          else
            !
            ! The requested vertical coordinate isn't a pressure, height, or level
            !
            cornerValues(:, :) = -huge(cornerValues)
            istatus = 1
          end if ! Interpolate in level, pressure, or height
        end if ! Surface pressure or 3D variable

        if(any(cornerValues(:, :) <= -huge(cornerValues))) then
          !
          ! One or more of our horizontal values isn't valid
          !   Likely the vertical location isn't within the range of pressures
          !
          obs_val = MISSING_R8
          istatus = 1
        else
          obs_val = sum(weights(:, :) * cornerValues(:, :))
          !
          ! Set istatus to ensure we want to assimilate this obs
          !
          if(vert_is_pressure(location) .and.                             &
             (lon_lat_height(3) < highest_obs_pressure_mb * mb_to_pa .or. &
              abs(lon_lat_height(2)) > max_obs_lat_degree) )  then
            istatus = 2
          else
            istatus = 0
          end if
        end if
      end if ! Check for valid latitude
    end if ! Check for valid variable


  end subroutine model_interpolate


  ! ----------------------------------------------------------------------------

  function get_model_time_step()
     type(time_type) :: get_model_time_step

     ! Returns the the time step of the model; the smallest increment
     ! in time that the model is capable of advancing the state in a given
     ! implementation. This interface is required for all applications.

     if ( .not. module_initialized ) call static_init_model

     ! Time_step_atmos is global static storage
     get_model_time_step =  Time_step_atmos

  end function get_model_time_step

  ! ----------------------------------------------------------------------------

  subroutine end_model()
    !------------------------------------------------------------------
    !
    ! Does any shutdown and clean-up needed for model. Can be a NULL
    ! INTERFACE if the model has no need to clean up storage, etc.

     if ( .not. module_initialized ) call static_init_model

    ! good style ... perhaps you could deallocate stuff (from static_init_model?).
    ! deallocate(state_loc)

  end subroutine end_model

  ! ----------------------------------------------------------------------------

  subroutine pert_model_state(state, pert_state, interf_provided)
    real(r8), intent(in)  :: state(:)
    real(r8), intent(out) :: pert_state(:)
    logical,  intent(out) :: interf_provided
    !------------------------------------------------------------------
    !
    ! Perturbs a model state for generating initial ensembles.
    ! The perturbed state is returned in pert_state.
    ! A model may choose to provide a NULL INTERFACE by returning
    ! .false. for the interf_provided argument. This indicates to
    ! the filter that if it needs to generate perturbed states, it
    ! may do so by adding an O(0.1) magnitude perturbation to each
    ! model state variable independently. The interf_provided argument
    ! should be returned as .true. if the model wants to do its own
    ! perturbing of states.

     if ( .not. module_initialized ) call static_init_model

     interf_provided = .false.

  end subroutine pert_model_state

  ! ----------------------------------------------------------------------------

  subroutine ens_mean_for_model(ens_mean)
     real(r8), intent(in) :: ens_mean(:)

     if ( .not. module_initialized ) call static_init_model

  end subroutine ens_mean_for_model

  ! ----------------------------------------------------------------------------
  !
  !  Public procedures that aren't part of the standard DART interface
  !
  ! ----------------------------------------------------------------------------

  subroutine init_model_instance(var)
    type(model_type), intent(out) :: var

    ! Initializes an instance of a model state variable
    !   In our case this means storage allocation

    if ( .not. module_initialized ) call static_init_model

    call end_model_instance(var)
    allocate(var%u   (num_lons, num_lats, num_levels), &
             var%v   (num_lons, num_lats, num_levels), &
             var%T   (num_lons, num_lats, num_levels), &
             var%delp(num_lons, num_lats, num_levels))
    if(num_tracers > 0) &
      allocate(var%tracers(num_lons, num_lats, num_levels, num_tracers), &
               var%tracerTypes(                            num_tracers))

  end subroutine init_model_instance

  ! ----------------------------------------------------------------------------

   subroutine end_model_instance(var)
     type(model_type), intent(inout) :: var

     ! Ends an instance of a model state variable
     !   i.e. frees allocated storage

     if ( .not. module_initialized ) call static_init_model

     if(associated(var%ps  )) deallocate(var%ps)
     if(associated(var%u   )) deallocate(var%u)
     if(associated(var%v   )) deallocate(var%v)
     if(associated(var%T   )) deallocate(var%T)
     if(associated(var%delp)) deallocate(var%delp)
     if(associated(var%tracers))     deallocate(var%tracers)
     if(associated(var%tracerTypes)) deallocate(var%tracerTypes)

  end subroutine end_model_instance

  ! ----------------------------------------------------------------------------

  subroutine read_model_init(rst_file_name, trc_file_name, var)
    character(len = *), intent(in)  :: rst_file_name, trc_file_name
    type(model_type),   intent(inout) :: var

    integer :: ncfileid, delpvarid, uvarid, vvarid, tvarid
    integer :: ncfileid_t, varID
    integer :: i

    if ( .not. module_initialized ) call static_init_model

    ! Do restart file first
    call nc_check(nf90_open(path = trim(rst_file_name), mode = nf90_nowrite, ncid = ncfileid), &
             'read_model_init', 'opening '// trim(rst_file_name))

    call nc_check(nf90_inq_varid(ncfileid, "DELP", delpvarid),'read_model_init','inquiring delp varid')
    call nc_check(nf90_get_var(ncfileid, delpvarid, var%delp),'read_model_init','getting delp var')

    call nc_check(nf90_inq_varid(ncfileid, "U", uvarid),'read_model_init','inquiring u varid')
    call nc_check(nf90_get_var(ncfileid, uvarid, var%u),'read_model_init','getting u var')

    call nc_check(nf90_inq_varid(ncfileid, "V", vvarid),'read_model_init','inquiring v varid')
    call nc_check(nf90_get_var(ncfileid, vvarid, var%v),'read_model_init','getting v var')

    call nc_check(nf90_inq_varid(ncfileid, "T", tvarid),'read_model_init','inquiring T varid')
    call nc_check(nf90_get_var(ncfileid, tvarid, var%T),'read_model_init','getting T var')

    call nc_check(nf90_close(ncfileid),'read_model_init','closing restart file')

    if(num_tracers > 0) then
      call nc_check(nf90_open(path = trim(trc_file_name), mode = nf90_nowrite, ncid = ncfileid_t), &
               'read_model_init', 'opening '// trim(tracer_files(1)))

      do i = 1, num_tracers
         call nc_check(nf90_inq_varid(ncfileid_t, tracer_names(i), varID), "read_model_init", &
              "inquiring varid for " // tracer_names(i))
         call nc_check(nf90_get_var(ncfileid_t, varID, var%tracers(:,:,:,i)), "read_model_init", &
              "getting " // tracer_names(i))
      end do
      call nc_check(nf90_close(ncfileid_t),'read_model_init','closing tracer file')
    end if

  end subroutine read_model_init

  ! ----------------------------------------------------------------------------

  subroutine write_model_init(rst_file_name, trc_file_name, var)
   character(len = *), intent(in) :: rst_file_name, trc_file_name
   type(model_type),   intent(in) :: var

   integer :: ncfileid, ncfileid_t, ncvarID
   integer :: i, j, k, t

   real, dimension(:, :, :), allocatable :: tempVar

   if ( .not. module_initialized ) call static_init_model

   allocate(tempVar(num_lons, num_lats, num_levels))

   call nc_check(nf90_open(path = trim(rst_file_name), mode = nf90_write, ncid = ncfileid), &
                 'write_model_init', 'opening ' // trim(rst_file_name))

   !
   ! We don't write the first row of u, which should be all 0s
   !
   call nc_check(nf90_inq_varid(ncfileid, "U", ncVarId),'write_model_init','inquiring u varid')
   call nc_check(nf90_put_var(ncfileid, ncVarId, var%u(:,2:num_lats,:), start = (/ 1, 2, 1/)),& 
                                                        'write_model_init','putting u var')
   call nc_check(nf90_put_var(ncfileid, ncVarId, RESHAPE( (/(0,i=1,num_lons*num_levels)/), &
                                                          (/num_lons,num_levels/)), &
                              start = (/ 1, 1, 1/), count = (/num_lons, 1, num_levels/)), &
                                                        'write_model_init','putting u var')

   !
   ! Quoth Steve Klein, quoting S.J. Lin - v in first and last row is ignored. Fine - we'll write
   !   what we've got anyway
   !
   call nc_check(nf90_inq_varid(ncfileid, "V", ncVarId),'write_model_init','inquiring v varid')
   call nc_check(nf90_put_var(ncfileid, ncVarId, var%v),'write_model_init','putting v var')

   !
   ! Quoth Steve Klein, quoting S.J. Lin - T, delP, and all tracers must have identical values
   !   for all lons in first and last lats
   ! We'll use the average value for lack of a better idea
   !
   tempVar(:, :, :) = var%delp(:, :, :)
   forall(k = 1:num_Levels)
     tempVar(:, 1,        k)  = sum(tempVar(:, 1,       k)) / real(num_lons)
     tempVar(:, num_lats, k) = sum(tempVar(:, num_lats, k)) / real(num_lons)
   end forall
   call nc_check(nf90_inq_varid(ncfileid, "DELP", ncVarId),'write_model_init','inquiring delp varid')
   call nc_check(nf90_put_var(ncfileid, ncVarId, tempVar),'write_model_init','putting delp var')

   tempVar(:, :, :) = var%t(:, :, :)
   forall(k = 1:num_Levels)
     tempVar(:, 1,        k)  = sum(tempVar(:, 1,       k)) / real(num_lons)
     tempVar(:, num_lats, k) = sum(tempVar(:, num_lats, k)) / real(num_lons)
   end forall
   call nc_check(nf90_inq_varid(ncfileid, "T", ncVarId),'write_model_init','inquiring T varid')
   call nc_check(nf90_put_var(ncfileid, ncVarId, tempVar),'write_model_init','putting T var')

   call nc_check(nf90_close(ncfileid),'write_model_init','closing restart file')

   if(num_tracers > 0) then
     call nc_check(nf90_open(path = trim(trc_file_name), mode = nf90_write, ncid = ncfileid_t), &
              'write_model_init', 'opening '// trim(trc_file_name))

     do i = 1, num_tracers
       tempVar = var%tracers(:, :, :, i)
       forall(k = 1:num_Levels)
          tempVar(:, 1,        k) = sum(tempVar(:, 1,       k)) / real(num_lons)
          tempVar(:, num_lats, k) = sum(tempVar(:, num_lats, k)) / real(num_lons)
        end forall
        call nc_check(nf90_inq_varid(ncfileid_t, tracer_names(i), ncVarId), &
                      "write_model_init", "inquiring varid for " // tracer_names(i))
        call nc_check(nf90_put_var(ncfileid_t, ncVarId, tempVar), &
                      "write_model_init", "putting " // tracer_names(i))
     end do

     call nc_check(nf90_close(ncfileid_t),'write_model_init','closing tracer file')
   end if

  end subroutine write_model_init

  ! ----------------------------------------------------------------------------

  subroutine prog_var_to_vector(model_var, state_vector)
    type(model_type),              intent(in ) :: model_var
    real(kind = r8), dimension(:), intent(out) :: state_vector

    ! -----------------------------------------------
    real, dimension(:, :),    allocatable :: psurf
    real, dimension(:,:,:),   allocatable :: u_var, v_var, t_var
    real, dimension(:,:,:,:), allocatable :: tracers
    integer :: i, j, k, t

    if ( .not. module_initialized ) call static_init_model

    if(size(state_vector) /= model_size)       &
      call error_handler(E_ERR, 'prog_var_to_vector', "State vector is incorrect size for model_type", &
                         source, revision, revdate)

    allocate(psurf(num_lons, num_lats),              &
             u_var(num_levels, num_lons, num_lats), &
             v_var(num_levels, num_lons, num_lats ), &
             t_var(num_levels, num_lons, num_lats ), &
             tracers(num_levels, num_lons, num_lats, num_tracers))

    if (associated(model_var%ps)) then
      psurf(:, :) = model_var%ps
    else if (associated(model_var%delp)) then
      psurf(:, :) = sum(model_var%delp, dim = 3) + model_top_pressure
    else
      call error_handler(E_ERR, 'prog_var_to_vector', "Neither delp nor ps is present in model_var", &
                         source, revision, revdate)
    end if

    forall(i = 1:num_levels, j=1:num_lons, k=1:num_lats )
          u_var(i, j, k) = model_var%u(j, k, i)
          v_var(i, j, k) = model_var%v(j, k, i)
          t_var(i, j, k) = model_var%t(j, k, i)
    end forall
    if(num_tracers > 0)                                                         &
      forall(i = 1:num_levels, j=1:num_lons, k=1:num_lats, t = 1:num_tracers ) &
        tracers(i, j, k, t) = model_var%tracers(j, k, i, t)
    !
    ! 2D field ps is ordered lon, lat; 3D fields are ordered lon, lat, level but are
    !   reordered to level, lon, lat before packing into vectors
    ! Tracers are ordered lon, lat, level, tracer_num, pack as level, lon, lat, tracer_num.
    if(num_tracers > 0) then
      state_vector(:) = (/ pack(psurf, .true.), &
                           pack(u_var,  .true.), &
                           pack(v_var,  .true.), &
                           pack(t_var,  .true.), &
                           pack(tracers,  .true.)  /)
    else
      state_vector(:) = (/ pack(psurf, .true.), &
                           pack(u_var,  .true.), &
                           pack(v_var,  .true.), &
                           pack(t_var,  .true.)  /)
    end if
    deallocate(psurf, u_var, v_var, t_var, tracers)
  end subroutine prog_var_to_vector

   ! ----------------------------------------------------------------------------

  subroutine vector_to_prog_var(state_vector, model_var)
    real(kind = r8), dimension(:), intent(in ) :: state_vector
    type(model_type),              intent(inout) :: model_var

    ! -----------------------------------------------
    !
    ! Local variables
    !
    real,    dimension(:, :),    allocatable :: psurf
    logical, dimension(:, :, :), allocatable :: allTrue
    integer :: start, finish

    real, dimension(:, :, :),    allocatable :: u_var, v_var, t_var
    real, dimension(:, :, :, :), allocatable :: tracers
    integer :: i, j, k, t

    ! -----------------------------------------------

    if ( .not. module_initialized ) call static_init_model

    if(size(state_vector) /= model_size)       &
      call error_handler(E_ERR, 'prog_var_to_vector', "State vector is incorrect size for model_type", &
                         source, revision, revdate)

    allocate(allTrue(num_levels, num_lons, num_lats))
    allocate(psurf(num_lons, num_lats))

    allTrue(:, :, :) = .true.
    !
    ! 2D field - surface pressure (may be mapped to delp)
    !
    start = 1; finish = start + (num_lons * num_lats) - 1
    psurf = unpack(state_vector(start:finish), allTrue(1, :, :), field = 0._r8)
    if (associated(model_var%ps)) then
      model_var%ps = psurf(:, :)
    else if (associated(model_var%delp)) then
      !
      ! Compute delp from fixed ak and bk terms and surface pressure
      !
      do k = 1, num_levels
        model_var%delp(:, :, k) =  (ak(k+1) - ak(k)) + (bk(k+1) - bk(k)) * psurf(:, :)
      end do
    else
      call error_handler(E_ERR, 'vector_to_prog_var', "Neither delp nor ps is present in model_var", &
                         source, revision, revdate)
    end if

    !
    ! 3D fields - u, v, and T
    !   These were stored in dim order level, lon, lat in the state vector and need
    !   to be mapped back to lon, lat, level
    !

    allocate(u_var(num_levels, num_lons, num_lats), &
             v_var(num_levels, num_lons, num_lats),  &
             t_var(num_levels, num_lons, num_lats) )

    start = finish + 1; finish = start + (num_levels * num_lons * num_lats) - 1
    u_var = unpack(state_vector(start:finish), allTrue(:, :, 1:num_lats), field = 0._r8)

    start = finish + 1; finish = start + (num_levels * num_lons * num_lats) - 1
    v_var = unpack(state_vector(start:finish), allTrue, field = 0._r8)

    start = finish + 1; finish = start + (num_levels * num_lons * num_lats) - 1
    t_var = unpack(state_vector(start:finish), allTrue, field = 0._r8)

    !Fill model_var components
    forall(i = 1:num_lons, j=1:num_lats, k=1:num_levels)
          model_var%u(i, j, k) = u_var(k, i, j)
          model_var%v(i, j, k) = v_var(k, i, j)
          model_var%t(i, j, k) = t_var(k, i, j)
    end forall

    !
    ! Tracers
    !   These were stored in dim order level, lon, lat, tracer_num in the state vector and need
    !   to be mapped back to lon, lat, level, tracer_num
    !
    allocate(tracers(num_levels, num_lons, num_lats, num_tracers))
    if(num_tracers > 0) then
      start = finish + 1; finish =  start + (num_levels * num_lons * num_lats * num_tracers) - 1
      if(finish /= model_size) &
        call error_handler(E_ERR, 'vector_to_prog_var', "Mismatch between model size and state vector", &
                           source, revision, revdate)
      tracers = unpack(state_vector(start:finish), &
                       spread(allTrue, dim = 4, nCopies = num_tracers), field = 0._r8)
      !Fill model_var%tracers
      forall(i = 1:num_lons, j=1:num_lats, k=1:num_levels, t = 1:num_tracers) &
        model_var%tracers(i, j, k, t) = tracers(k, i, j, t)
    end if


    deallocate(allTrue)
    deallocate(psurf, u_var, v_var, t_var, tracers)
  end subroutine vector_to_prog_var

  ! ----------------------------------------------------------------------------
  !
  !  Private procedures
  !
  ! ----------------------------------------------------------------------------

  subroutine read_dimension_info(ncFileId)
    integer, intent(in) :: ncFileId
    !
    ! Called by static_model_init; fills in global variables related to the coordinates
    ! DimStrings must be in this order 
    integer :: i, temp_size
    integer, dimension(num_dims) :: dim_ids
    character(len=5), dimension(6) :: DimStrings = &
       (/ 'lat  ', 'latu ', 'lon  ', 'lonv ', 'pfull', 'phalf' /)
    logical :: inOrder = .true.

    do i = 1, num_dims
      call nc_check(nf90_inq_dimid(ncfileid, trim(dim_names(i)), dim_ids(i)), &
                   'static_init_model', 'looking for dimension id for '// trim(dim_names(i)))
      call nc_check(nf90_Inquire_Dimension(ncfileid,  dim_ids(i), len = dim_lens(i)), &
                   'static_init_model', 'looking for size of '// trim(dim_names(i)))
      call nc_check(nf90_inq_varid(ncfileid, trim(dim_names(i)), dim_var_ids(i)), &
                   'static_init_model', 'looking for variable id for '// trim(dim_names(i)))
    end do

    do i = 1, size(dim_names)
       if ( trim(dim_names(i)) /= trim(DimStrings(i)) ) inOrder = .false.
    end do

    if ( .not. inOrder ) then
         call error_handler(E_ERR, 'read_dimension_info', &
            "Mapping between dim names and variables is messed up", source, revision, revdate)
    endif

    if(dim_lens(1) /= dim_lens(2) .or. dim_lens(3) /= dim_lens(4) .or. &
       dim_lens(5) /= dim_lens(6) - 1) &
       call error_handler(E_ERR, 'read_dimension_info', &
           "Dimension sizes aren't what we expect", source, revision, revdate)

    num_lats = dim_lens(1)
    allocate(lat(num_lats), latu(num_lats))
    call nc_check(nf90_get_var(ncfileid, dim_var_ids(1), lat), &
                  'static_init_model', 'reading values of '// trim(dim_names(1)))
    call nc_check(nf90_get_var(ncfileid, dim_var_ids(2), latu), &
                  'static_init_model', 'reading values of '// trim(dim_names(2)))

    num_lons = dim_lens(3)
    allocate(lon(num_lons), lonv(num_lons))
    call nc_check(nf90_get_var(ncfileid, dim_var_ids(3), lon), &
                  'static_init_model', 'reading values of '// trim(dim_names(3)))
    call nc_check(nf90_get_var(ncfileid, dim_var_ids(4), lonv), &
                  'static_init_model', 'reading values of '// trim(dim_names(4)))

    num_levels = dim_lens(5)
    allocate(pfull(num_levels), phalf(num_levels + 1))
    call nc_check(nf90_get_var(ncfileid, dim_var_ids(5), pfull), &
                  'static_init_model', 'reading values of '// trim(dim_names(5)))
    call nc_check(nf90_get_var(ncfileid, dim_var_ids(6), phalf), &
                  'static_init_model', 'reading values of '// trim(dim_names(6)))

  end subroutine read_dimension_info

  ! ----------------------------------------------------------------------------

  subroutine copy_dim_var_pair(newFileId, dimName, dimLength, oldFileId, routineName)
    integer,            intent(in) :: newFileId, dimLength, oldFileId
    character(len = *), intent(in) :: dimName, routineName
    !
    ! Create a dimension and corresponding 1d double variable in newFileId
    !   using the definitions in oldFileId
    !

    integer :: att, num_atts, newDimId, oldVarId, newVarId
    character(len = NF90_MAX_NAME) :: attName

    call nc_check(nf90_def_dim(newFileId, trim(dimName), dimLength, newDimId),  &
                  trim(routineName), 'def_dim ' // trim(dimName) )
    call nc_check(nf90_def_var(newFileId, trim(dimName), nf90_double, newDimId, newVarId),  &
                  trim(routineName), 'def_var ' // trim(dimName) )
    !
    ! Copy over all the attributes from the old file to the newnostics file
    !
    call nc_check(nf90_inq_varid(oldFileId, trim(dimName), oldVarId),  &
                  trim(routineName), 'getting varid from old file for ' // trim(dimName) )
    call nc_check(nf90_inquire_variable(oldFileId, oldVarId, nAtts = num_atts),  &
                  trim(routineName), 'reading number of attributes from ' // trim(dimName) )
    do att = 1, num_atts
      call nc_check(nf90_inq_attname(oldFileId, oldVarId, att, attName),  &
                    trim(routineName), 'reading number of attributes from ' // trim(dimName) )
      call nc_check(nf90_copy_att(oldFileId, oldVarId, trim(attName), newFileId, newVarId),  &
                    trim(routineName), 'copying attribute ' // trim(attName) // ' for variable ' // trim(dimName) )
    end do
  end subroutine copy_dim_var_pair

  ! ----------------------------------------------------------------------------

  subroutine define_3d_real_var(newFileId, oldFileId, varName, dimIds, routineName)
    integer,               intent(in) :: newFileId, oldFileId
    character(len = *),    intent(in) :: varName, routineName
    integer, dimension(:), intent(in) :: dimIds
    !
    ! Define a new variable in newFileId, then copy all the attributes from the same variable
    !   in oldFileId
    !

    integer :: att, num_atts, oldVarId, newVarId
    character(len = NF90_MAX_NAME) :: attName

    call nc_check(nf90_def_var(newFileId, trim(varName), nf90_real, dimIds, newVarId), &
                  trim(routineName),'def_var ' // trim(varName))

    call nc_check(nf90_inq_varid(oldFileId, trim(varName), oldVarId),  &
                  trim(routineName), 'getting varid from old file for ' // trim(varName) )
    call nc_check(nf90_inquire_variable(oldFileId, oldVarId, nAtts = num_atts),  &
                  trim(routineName), 'reading number of attributes from ' // trim(varName) )
    do att = 1, num_atts
      call nc_check(nf90_inq_attname(oldFileId, oldVarId, att, attName),  &
                    trim(routineName), 'reading attributes from ' // trim(varName) )
      call nc_check(nf90_copy_att(oldFileId, oldVarId, trim(attName), newFileId, newVarId),  &
                    trim(routineName), 'copying attribute ' // trim(attName) // ' for variable ' // trim(varName) )
    end do

  end subroutine define_3d_real_var

  ! ----------------------------------------------------------------------------

  subroutine write_1D_values(ncFileId, varName, values, routineName)
    integer,            intent(in) :: ncFileId
    character(len = *), intent(in) :: varName, routineName
    real, dimension(:), intent(in) :: values

    integer :: varId

    call nc_check(nf90_inq_varid(ncFileId, trim(varName), varId), trim(routineName), 'getting varid for ' // trim(varName) )
    call nc_check(nf90_put_var(ncFileId, varId, values), trim(routineName), 'writing ' // trim(varName) )

  end subroutine write_1D_values

  !------------------------------------------------------------------------------------------
  pure function interpolate1D(desiredLocation, values, locations)
    real(kind = r8),               intent(in) :: desiredLocation
    real(kind = r8), dimension(:), intent(in) :: values, locations
    real(kind = r8)                           :: interpolate1D
    !
    ! Given a set of values at a set of locations (ordered from lowest to highest coordinate value)
    !   return the value at an arbirtrary point using linear interoplation
    !
    real(kind =r8) :: fraction
    integer        :: n, vIndex ! Such that locations(vIndex) <= desiredLocation < location(vIndex+1)
    ! ----------------------

    n = size(values) ! = size(locations)
    if(desiredLocation < locations(1) .or. desiredLocation > locations(n)) then
      !
      ! Desired location is out of bounds
      !
      interpolate1D = -huge(interpolate1D)
    else if (abs(desiredLocation - locations(n)) <= spacing(desiredLocation)) then
      !
      ! Desired location is exactly at the last location
      !
      interpolate1D = values(n)
    else
      vIndex = findIndex(desiredLocation, locations)
      fraction = 1._r8 - (desiredLocation         - locations(vIndex)) / &
                         (locations(vIndex + 1)   - locations(vIndex))
      interpolate1D = fraction * values(vIndex) + (1._r8 - fraction) * values(vIndex + 1)
    end if
  end function interpolate1D
  !------------------------------------------------------------------------------------------

  pure function findIndex(value, table, firstGuess)
    real(kind=r8),               intent( in) :: value
    real(kind=r8), dimension(:), intent( in) :: table
    integer,          optional,  intent( in) :: firstGuess
    integer                                  :: findIndex
    !
    ! Find the index i into the table such that table(i) <= value < table(i+i)
    !   This is modeled after routine "hunt" from Numerical Recipes, 2nd ed.,
    !   pg 112. Here we know that the values in the table are always increasing,
    !   that every value should be spanned by the table entries, and the firstGuess
    !   always makes sense.

    ! Local variables
    integer :: lowerBound, upperBound, midPoint
    integer :: increment

    ! Hunting; only done if a first guess is supplied
    !  Move upper and lower bounds around until the value is spanned by
    !   table(lowerBound) and table(upperBound). Make the interval twice as
    !   big at each step
    if(present(firstGuess)) then
      lowerBound = firstGuess
      increment = 1
      huntingLoop: do
        upperBound = min(lowerBound + increment, size(table))
        if(lowerBound == size(table) .or. &
           (table(lowerBound) <= value .and. table(upperBound) > value)) exit huntingLoop
        if(table(lowerBound) > value) then
          upperBound = lowerBound
          lowerBound = max(upperBound - increment, 1)
        else
          ! Both table(lowerBound) and table(upperBound) are <= value
          lowerBound = upperBound
        end if
        increment = increment * 2
      end do huntingLoop
    else
      lowerBound = 0; upperBound = size(table)
    end if

    ! Bisection: figure out which half of the remaining interval holds the
    !   desired value, discard the other half, and repeat
    bisectionLoop: do
      if(lowerBound == size(table) .or. upperBound <= lowerBound + 1) exit bisectionLoop
      midPoint = (lowerBound + upperBound)/2
      if(value >= table(midPoint)) then
        lowerBound = midPoint
      else
        upperBound = midPoint
      end if
    end do bisectionLoop

    findIndex = lowerBound
  end function findIndex

  !------------------------------------------------------------------------------------------

end module model_mod

! <next few lines under version control, do not edit>
! $URL$
! $Id$
! $Revision$
! $Date$
