! DART software - Copyright UCAR. This open source software is provided
! by UCAR, "as is", without charge, subject to all terms of use at
! http://www.image.ucar.edu/DAReS/DART/DART_download

module model_mod

use        types_mod, only : i4, r8, i8, MISSING_R8, metadatalength, vtablenamelength
use time_manager_mod, only : time_type, set_calendar_type, get_time, set_date, get_date
use     location_mod, only : location_type, get_close_type, get_close_obs, get_dist, &
                             convert_vertical_obs, convert_vertical_state, &
                             set_location, set_location_missing, VERTISLEVEL, &
                             get_location, loc_get_close_state => get_close_state
use    utilities_mod, only : register_module, error_handler, E_ERR, E_MSG, logfileunit, &
                             nmlfileunit, do_output, do_nml_file, do_nml_term, &
                             find_namelist_in_file, check_namelist_read,to_upper, &
                             file_exist
use netcdf_utilities_mod, only : nc_add_global_attribute, nc_synchronize_file, &
                                 nc_add_global_creation_time, nc_begin_define_mode, &
                                 nc_end_define_mode, nc_check
use state_structure_mod, only : add_domain, get_domain_size
use ensemble_manager_mod, only : ensemble_type
use distributed_state_mod, only : get_state
use default_model_mod, only : pert_model_copies, nc_write_model_vars, init_conditions, &
                              init_time, adv_1step, shortest_time_between_assimilations
use         dart_cice_mod, only : get_horiz_grid_dims, get_ncat_dim, read_horiz_grid
use   state_structure_mod, only : state_structure_info,get_index_start, get_num_variables, &
                                  get_dart_vector_index, get_model_variable_indices
use          obs_kind_mod, only : QTY_SEAICE_AGREG_CONCENTR  , &
                                  QTY_SEAICE_AGREG_VOLUME    , &
                                  QTY_SEAICE_AGREG_SNOWVOLUME, &
                                  QTY_SEAICE_AGREG_THICKNESS , &
                                  QTY_SEAICE_AGREG_SNOWDEPTH , &
                                  QTY_SEAICE_CATEGORY        , &
                                  QTY_U_SEAICE_COMPONENT     , &
                                  QTY_V_SEAICE_COMPONENT     , &
                                  QTY_SEAICE_ALBEDODIRVIZ    , &
                                  QTY_SEAICE_ALBEDODIRNIR    , &
                                  QTY_SEAICE_ALBEDOINDVIZ    , &
                                  QTY_SEAICE_ALBEDOINDNIR    , &
                                  QTY_SEAICE_CONCENTR        , &
                                  QTY_SEAICE_VOLUME          , &
                                  QTY_SEAICE_SNOWVOLUME      , &
                                  QTY_SEAICE_SURFACETEMP     , &
                                  QTY_SEAICE_FIRSTYEARAREA   , &
                                  QTY_SEAICE_ICEAGE          , &
                                  QTY_SEAICE_LEVELAREA       , &
                                  QTY_SEAICE_LEVELVOLUME     , &
                                  QTY_SEAICE_MELTPONDAREA    , &
                                  QTY_SEAICE_MELTPONDDEPTH   , &
                                  QTY_SEAICE_MELTPONDLID     , &
                                  QTY_SEAICE_MELTPONDSNOW    , &
                                  QTY_SEAICE_SALINITY001     , &
                                  QTY_SEAICE_SALINITY002     , &
                                  QTY_SEAICE_SALINITY003     , &
                                  QTY_SEAICE_SALINITY004     , &
                                  QTY_SEAICE_SALINITY005     , &
                                  QTY_SEAICE_SALINITY006     , &
                                  QTY_SEAICE_SALINITY007     , &
                                  QTY_SEAICE_SALINITY008     , &
                                  QTY_SEAICE_ICEENTHALPY001  , &
                                  QTY_SEAICE_ICEENTHALPY002  , &
                                  QTY_SEAICE_ICEENTHALPY003  , &
                                  QTY_SEAICE_ICEENTHALPY004  , &
                                  QTY_SEAICE_ICEENTHALPY005  , &
                                  QTY_SEAICE_ICEENTHALPY006  , &
                                  QTY_SEAICE_ICEENTHALPY007  , &
                                  QTY_SEAICE_ICEENTHALPY008  , &
                                  QTY_SEAICE_SNOWENTHALPY001 , &
                                  QTY_SEAICE_SNOWENTHALPY002 , &
                                  QTY_SEAICE_SNOWENTHALPY003 , &
                                  QTY_DRY_LAND               , &
                                  QTY_SOM_TEMPERATURE        , &
                                  QTY_SEAICE_FY              , &
                                  QTY_SEAICE_AGREG_FY        , &
                                  QTY_SEAICE_AGREG_SURFACETEMP,&
                                  get_index_for_quantity     , &
                                  get_name_for_quantity

use netcdf

implicit none
private

! required routines by DART code - will be called from filter and other
! DART executables. interfaces to these routines are fixed and cannot
! be changed in any way.
public :: get_model_size,         &
          adv_1step,              &
          get_state_meta_data,    &
          model_interpolate,      &
          shortest_time_between_assimilations, &
          end_model,              &
          static_init_model,      &
          nc_write_model_atts, &
          init_time, &
          init_conditions, &
          check_sfctemp_var

! required routines where code is in other modules
public :: nc_write_model_vars,    &
          pert_model_copies,      &
          get_close_obs,          &
          get_close_state,        &
          convert_vertical_obs,   &
          convert_vertical_state, &
          read_model_time, &
          write_model_time

character(len=256), parameter :: source   = 'icepack/model_mod.f90'

logical, save :: module_initialized = .false.

! message strings
character(len=512) :: string1
character(len=512) :: string2
character(len=512) :: string3

! DART state vector contents are specified in the input.nml:&model_nml namelist.
integer, parameter :: max_state_variables = 10
integer, parameter :: num_state_table_columns = 3
character(len=vtablenamelength) :: variable_table( max_state_variables, num_state_table_columns )
integer :: state_kinds_list( max_state_variables )
logical :: update_var_list( max_state_variables )

integer, parameter :: VAR_NAME_INDEX = 1
integer, parameter :: VAR_QTY_INDEX = 2
integer, parameter :: VAR_UPDATE_INDEX = 3

integer  :: model_size

real(r8), allocatable :: TLAT(:), TLON(:)

type(time_type) :: model_time, model_timestep

integer :: Nx=-1
integer :: Ncat=-1
integer :: domain_id,nfields

! Items in the model_nml
real(r8) :: model_perturbation_amplitude = 0.01
character(len=metadatalength) :: model_state_variables(max_state_variables * num_state_table_columns ) = ' '
integer  :: debug = 1
integer  :: grid_oi = 3

namelist /model_nml/  &
   model_perturbation_amplitude, &
   model_state_variables,        &
   debug,                        &
   grid_oi

contains

!------------------------------------------------------------------
! Called to do one time initialization of the model. Reads the
! namelist, defines information about the model size and model
! timestep, initializes module variables, and calls add_domain()
! to set what data should be read into the state

subroutine static_init_model()

integer  :: iunit, io, ss, dd

if ( module_initialized ) return ! only need to do this once
 
module_initialized = .true.

call find_namelist_in_file("input.nml", "model_nml", iunit)
read(iunit, nml = model_nml, iostat = io)
call check_namelist_read(iunit, io, "model_nml")

call error_handler(E_MSG,'static_init_model','model_nml values are',' ',' ',' ')
if (do_nml_file()) write(nmlfileunit, nml=model_nml)
if (do_nml_term()) write(*, nml=model_nml)

call set_calendar_type('Gregorian')

model_timestep = shortest_time_between_assimilations()

call get_time(model_timestep,ss,dd)

write(string1, *) 'assimilation period is ', dd,' days ', ss,' seconds'
call error_handler(E_MSG, 'static_init_model', string1, source)

call get_horiz_grid_dims(Nx)
call get_ncat_dim(Ncat)

call verify_state_variables(model_state_variables, nfields, variable_table, &
                            state_kinds_list, update_var_list)

allocate(TLAT(Nx), TLON(Nx))

call read_horiz_grid(Nx, TLAT, TLON)
 
if (do_output()) write(logfileunit, *) 'Using grid : Nx, Ncat = ', Nx, Ncat
if (do_output()) write(*, *) 'Using grid : Nx, Ncat = ', Nx, Ncat

domain_id = add_domain('cice.r.nc', nfields, &
                       var_names = variable_table(1:nfields, VAR_NAME_INDEX), &
                       kind_list = state_kinds_list(1:nfields), &
                       update_list = update_var_list(1:nfields))

if (debug > 2) call state_structure_info(domain_id)

model_size = get_domain_size(domain_id)
if (do_output()) write(*, *) 'model_size = ', model_size


end subroutine static_init_model

!------------------------------------------------------------------
! Returns the number of items in the state vector as an integer. 

function get_model_size()

integer(i8) :: get_model_size

get_model_size = model_size

end function get_model_size

!------------------------------------------------------------------
! Given a state handle, a location, and a model state variable type,
! interpolates the state variable fields to that location and returns
! the values in expected_obs. The istatus variables should be
! returned as 0 unless there is some problem in computing the 
! interpolation, in which case an alternate value should be returned.

subroutine model_interpolate(state_handle, ens_size, location, obs_type, expected_obs, istatus, thick_flag)

type(ensemble_type), intent(in) :: state_handle
integer,             intent(in) :: ens_size
type(location_type), intent(in) :: location
integer,             intent(in) :: obs_type
real(r8),           intent(out) :: expected_obs(ens_size) ! array of interpolated values
integer,            intent(out) :: istatus(ens_size)
logical,optional,    intent(inout) :: thick_flag

!local vars
real(r8)       :: loc_array(3), llon, llat
integer(i8)    :: base_offset
integer        :: cat_index, cat_signal, icat, cat_signal_interm
real(r8)       :: expected_aggr_conc(ens_size)
integer        :: set_obstype

!Fei---need aicen*fyn to calculate the aggregate FY concentration
real(r8)       :: expected_conc(ens_size)
real(r8)       :: expected_fy(ens_size)
real(r8)       :: expected_tsfc(ens_size)
real(r8)       :: temp(ens_size)
real(r8)       :: temp1(ens_size)

if ( .not. module_initialized ) call static_init_model

expected_obs(:) = MISSING_R8 ! represents a bad value in DART
istatus(:) = 99

loc_array = get_location(location)
llon    = loc_array(1)
llat    = loc_array(2)
cat_index = int(loc_array(3))

if (obs_type == QTY_SEAICE_CATEGORY) then
   if (cat_index <= Ncat) then
      istatus      = 0
      expected_obs = cat_index
      RETURN
   endif
endif
if (debug > 1) then
   print *, 'requesting interpolation of ', obs_type, ' at ', llon, llat, cat_index
endif

SELECT CASE (obs_type)
   CASE (QTY_SEAICE_AGREG_THICKNESS )  ! these kinds require aggregating 3D vars to make a 2D var
      if (any(variable_table(:,1)=='hi')) then
        cat_signal = 1 ! for extra special procedure to aggregate
        !base_offset = get_index_start(domain_id, get_varid_from_kind(QTY_SEAICE_AGREG_THICKNESS))
        thick_flag = .true.
        base_offset = cat_index
        set_obstype = obs_type
      else
        set_obstype = QTY_SEAICE_VOLUME
        cat_signal = 1 ! for extra special procedure to aggregate
        !base_offset = get_index_start(domain_id, get_varid_from_kind(QTY_SEAICE_VOLUME))
        base_offset = cat_index
      endif
   CASE (QTY_SEAICE_AGREG_SNOWDEPTH )  ! these kinds require aggregating 3D vars to make a 2D var
      if (any(variable_table(:,1)=='hs')) then
        cat_signal = 1 ! for extra special procedure to aggregate
        !base_offset = get_index_start(domain_id, get_varid_from_kind(QTY_SEAICE_AGREG_SNOWDEPTH))
        base_offset = cat_index
        thick_flag = .true.
        set_obstype = obs_type
      else
        set_obstype = QTY_SEAICE_SNOWVOLUME
        cat_signal = 1 ! for extra special procedure to aggregate
        !base_offset = get_index_start(domain_id, get_varid_from_kind(QTY_SEAICE_SNOWVOLUME))
        base_offset = cat_index
      endif
   CASE (QTY_SEAICE_AGREG_CONCENTR )   ! these kinds require aggregating a 3D var to make a 2D var
      cat_signal = 0 ! for aggregate variable, send signal to lon_lat_interp
      set_obstype = obs_type
      base_offset = get_index_start(domain_id, get_varid_from_kind(QTY_SEAICE_CONCENTR))
   CASE (QTY_SEAICE_AGREG_VOLUME   )   ! these kinds require aggregating a 3D var to make a 2D var
      cat_signal = 0 ! for aggregate variable, send signal to lon_lat_interp
      set_obstype = obs_type
      base_offset = get_index_start(domain_id, get_varid_from_kind(QTY_SEAICE_VOLUME))
   CASE (QTY_SEAICE_AGREG_SNOWVOLUME ) ! these kinds require aggregating a 3D var to make a 2D var
      cat_signal = 0 ! for aggregate variable, send signal to lon_lat_interp
      set_obstype = obs_type
      base_offset = get_index_start(domain_id, get_varid_from_kind(QTY_SEAICE_SNOWVOLUME))
   CASE (QTY_SEAICE_AGREG_SURFACETEMP) ! FEI need aicen to average the temp, have not considered open water temp yet
      if (any(variable_table(:,1)=='Tsfc')) then
        cat_signal = 1 ! for extra special procedure to aggregate
        base_offset = get_index_start(domain_id, get_varid_from_kind(QTY_SEAICE_AGREG_SURFACETEMP))
        thick_flag = .true.
        set_obstype = obs_type
      else
        cat_signal = -3
        set_obstype = QTY_SEAICE_SURFACETEMP
        base_offset = get_index_start(domain_id, get_varid_from_kind(QTY_SEAICE_SURFACETEMP))
      endif
   CASE (QTY_SOM_TEMPERATURE) ! these kinds are 1d variables
      cat_signal = 1 ! for extra special procedure to aggregate
      set_obstype = obs_type
      !base_offset = get_index_start(domain_id,get_varid_from_kind(QTY_SOM_TEMPERATURE))
      base_offset = cat_index
   CASE (QTY_SEAICE_CONCENTR       , &  ! these kinds have an additional dim for category
         QTY_SEAICE_FY       , &
         QTY_SEAICE_VOLUME         , &
         QTY_SEAICE_SNOWVOLUME     , &
         QTY_SEAICE_SURFACETEMP    , &
         QTY_SEAICE_FIRSTYEARAREA  , &
         QTY_SEAICE_ICEAGE         , &
         QTY_SEAICE_LEVELAREA      , &
         QTY_SEAICE_LEVELVOLUME    , &
         QTY_SEAICE_MELTPONDAREA   , &
         QTY_SEAICE_MELTPONDDEPTH  , &
         QTY_SEAICE_MELTPONDLID    , &
         QTY_SEAICE_MELTPONDSNOW   , &
         QTY_SEAICE_SALINITY001    , &
         QTY_SEAICE_SALINITY002    , &
         QTY_SEAICE_SALINITY003    , &
         QTY_SEAICE_SALINITY004    , &
         QTY_SEAICE_SALINITY005    , &
         QTY_SEAICE_SALINITY006    , &
         QTY_SEAICE_SALINITY007    , &
         QTY_SEAICE_SALINITY008    , &
         QTY_SEAICE_ICEENTHALPY001 , &
         QTY_SEAICE_ICEENTHALPY002 , &
         QTY_SEAICE_ICEENTHALPY003 , &
         QTY_SEAICE_ICEENTHALPY004 , &
         QTY_SEAICE_ICEENTHALPY005 , &
         QTY_SEAICE_ICEENTHALPY006 , &
         QTY_SEAICE_ICEENTHALPY007 , &
         QTY_SEAICE_ICEENTHALPY008 , &
         QTY_SEAICE_SNOWENTHALPY001, &
         QTY_SEAICE_SNOWENTHALPY002, &
         QTY_SEAICE_SNOWENTHALPY003  )
      ! move pointer to the particular category
      ! then treat as 2d field in lon_lat_interp
      
      base_offset = get_index_start(domain_id, get_varid_from_kind(obs_type))
      base_offset = base_offset + (cat_index-1)
      base_offset = cat_index
      set_obstype = obs_type
      cat_signal = 1 ! now same as boring 2d field
  CASE DEFAULT
      ! Not a legal type for interpolation, return istatus error
      istatus = 15
      return
END SELECT

if (cat_signal == -2) then
   temp = 0.0_r8
   temp1= 0.0_r8
   do icat = 1,Ncat
      !reads in aicen
      cat_signal_interm = 1
      base_offset = get_index_start(domain_id,get_varid_from_kind(QTY_SEAICE_CONCENTR))
      base_offset = base_offset + (icat-1) * Nx
      call lon_lat_interpolate(state_handle, ens_size, base_offset, llon, llat, set_obstype, cat_signal_interm, expected_conc, istatus)
      !reads in fyn
      cat_signal_interm = 1
      base_offset = get_index_start(domain_id,get_varid_from_kind(QTY_SEAICE_FY))
      base_offset = base_offset + (icat-1) * Nx
      call lon_lat_interpolate(state_handle, ens_size, base_offset, llon, llat, set_obstype, cat_signal_interm, expected_fy, istatus)
   temp = temp + expected_conc * expected_fy  !sum(aicen*fyn) = FY % over ice
   temp1= temp1+ expected_conc                        !sum(aicen) = aice

      if ((any(expected_conc<0.0) .or. any(expected_conc>1.0)) .and. (debug > 1)) then
      print*,'obstype FY expected sicn:',expected_conc
      print*,'FY sicn lat lon:',llat,llon
      endif
      if ((any(expected_fy>1.0) .or. any(expected_fy<0.0)) .and. (debug > 1)) then
      print*,'obstype FY expected fyn:',expected_fy,llat,llon
      print*,'FY fyn lat lon:',llat,llon
      endif

   end do
   expected_obs = temp/max(temp1,1.0e-8)  !sum(aicen*fyn)/aice = FY % in the gridcell
else if (cat_signal == -3 ) then
   temp = 0.0_r8
   temp1= 0.0_r8
   do icat = 1,Ncat
      !reads in aicen
      cat_signal_interm = 1
      base_offset = get_index_start(domain_id,get_varid_from_kind(QTY_SEAICE_CONCENTR))
      base_offset = base_offset + (icat-1) * Nx
      call lon_lat_interpolate(state_handle, ens_size, base_offset, llon, llat, set_obstype, cat_signal_interm, expected_conc, istatus)
      !reads in Tsfcn
      cat_signal_interm = 1
      base_offset = get_index_start(domain_id,get_varid_from_kind(QTY_SEAICE_SURFACETEMP))
      base_offset = base_offset + (icat-1) * Nx
      call lon_lat_interpolate(state_handle, ens_size, base_offset, llon, llat, set_obstype, cat_signal_interm, expected_tsfc, istatus)
      if ((any(expected_conc<0.0) .or. any(expected_conc>1.0)) .and. (debug > 1)) then
      print*,'obstype TSFC expected sicn:',expected_conc
      print*,'TSFC sicn lat lon:',llat,llon
      endif
      if ((any(expected_tsfc>50.0) .or. any(expected_tsfc<-100.0)) .and. (debug > 1)) then
      print*,'obstype TSFC expected tsfcn:',expected_tsfc
      print*,'TSFC tsfcn lat lon:',llat,llon
      endif
      temp = temp + expected_conc * expected_tsfc  !sum(aicen*Tsfcn)
      temp1= temp1+ expected_conc                  !sum(aicen) = aice
   end do
   expected_obs = temp/max(temp1,1.0e-8)  !sum(aicen*Tsfcn)/aice = Tsfc ;averaged temperature over sea-ice covered portion
   if ((any(expected_obs>50.0) .or. any(expected_obs<-100.0)) .and. (debug > 1)) then
      print*,'obstype TSFC expected obs:',expected_obs
      print*,'TSFC tsfc lat lon:' ,llat,llon
      print*,'temp:',temp
      print*,'temp1:',temp1
   endif
else
    call lon_lat_interpolate(state_handle, ens_size, base_offset, llon, llat, set_obstype, cat_signal, expected_obs, istatus)
    
      if (any(expected_obs<0.0) .and. (debug > 1)) then
      print*,'obstype SIC expected concs:',expected_obs
      print*,'SIC sic negative lat lon:',llat,llon
      endif
      if (any(expected_obs>1.0) .and. (debug > 1)) then
      print*,'obstype SIC expected concs:',expected_obs
      print*,'SIC sic positive lat lon:',llat,llon
      endif
endif

if (cat_signal == -1) then
      ! we need to know the aggregate sea ice concentration for these special cases
      base_offset = get_index_start(domain_id, get_varid_from_kind(QTY_SEAICE_CONCENTR))
      base_offset = base_offset + (cat_index-1)
      call lon_lat_interpolate(state_handle, ens_size, base_offset, llon, llat, set_obstype, cat_signal, expected_aggr_conc, istatus)
      expected_obs = expected_obs/max(expected_aggr_conc,1.0e-8)  ! hope this is allowed so we never divide by zero

      if ((any(expected_aggr_conc<0.0) .or. any(expected_aggr_conc>1.0)) .and. (debug > 1)) then
      print*,'obstype SIT expected conc:',expected_aggr_conc
      print*,'SIT sic lat lon:',llat,llon
      endif

endif

if (debug > 1) print *, 'interp val, istatus = ', expected_obs, istatus, size(expected_obs)

! This should be the result of the interpolation of a
! given obs_type of variable at the given location.

! The return code for successful return should be 0. 
! Any positive number is an error.
! Negative values are reserved for use by the DART framework.

end subroutine model_interpolate

!------------------------------------------------------------------

subroutine lon_lat_interpolate(state_handle, ens_size, offset, lon, lat, var_type, cat_signal, expected_obs, istatus)
type(ensemble_type), intent(in)  :: state_handle
integer,             intent(in)  :: ens_size
integer(i8),         intent(in)  :: offset
real(r8),            intent(in)  :: lon, lat
integer,             intent(in)  :: var_type
integer,             intent(in)  :: cat_signal
real(r8),            intent(out) :: expected_obs(ens_size)
integer,             intent(out) :: istatus(ens_size)

real(r8) :: work_expected_obs(ens_size)
integer  :: e, iterations, Niterations
integer(i8) :: state_index

if ( .not. module_initialized ) call static_init_model

istatus = 0
if (var_type == 14) then
  e = 1
else if (var_type == 15) then
  e = 2
else if (var_type == 16) then
  e = 3
endif
if ( cat_signal < 1 )  then
   Niterations = Ncat ! only iterate if aggregating over all types
else
   Niterations = 1 ! no need to iterate
endif
work_expected_obs = 0.0_r8
expected_obs = 0.0_r8
do iterations = 1, Niterations

   ! FIXME: this should use the state structure routine 'get_dart_vector_index'
   ! to get the start of the next category layer.  this code assumes it knows
   ! exactly how the state vector is laid out (reasonable, but might not be true
   ! in future versions of dart.)
   !next_offset = offset + (iterations-1)*Nx 
   !print*,'offset',offset
   state_index = get_dart_vector_index(grid_oi,int(offset,i4),1, domain_id, e)
   work_expected_obs = get_state(state_index,state_handle)
   expected_obs = expected_obs+work_expected_obs
enddo
end subroutine lon_lat_interpolate

!------------------------------------------------------------------
! Given an integer index into the state vector structure, returns 
! the associated location.

subroutine get_state_meta_data(index_in, location, var_type)

integer(i8),         intent(in)  :: index_in
type(location_type), intent(out) :: location
integer,             intent(out), optional :: var_type

real(r8) :: lat, lon, rcat
integer  :: ni_index, hold_index, cat_index, local_var, var_id
! these should be set to the actual location and state quantity

if ( .not. module_initialized ) call static_init_model

call get_model_variable_indices(index_in, ni_index, cat_index, hold_index, var_id=var_id)
call get_state_kind(var_id, local_var)

lon = TLON(ni_index)
lat = TLAT(ni_index)

if (debug > 5) print *, 'lon, lat, cat_index = ', lon, lat, cat_index
rcat     = cat_index*1.0_r8
location = set_location(lon, lat, rcat, VERTISLEVEL)

if (present(var_type)) then
   var_type = local_var
endif

end subroutine get_state_meta_data

!------------------------------------------------------------------
! Given an integer index into the state vector structure, returns the kind,
! and both the starting offset for this kind, as well as the offset into
! the block of this kind.

subroutine get_state_kind(var_ind, var_type)
 integer, intent(in)  :: var_ind
 integer, intent(out) :: var_type

if ( .not. module_initialized ) call static_init_model

var_type = state_kinds_list(var_ind)

end subroutine get_state_kind

!------------------------------------------------------------------
! Does any shutdown and clean-up needed for model

subroutine end_model()

deallocate(TLAT,TLON)

end subroutine end_model

!------------------------------------------------------------------
! write any additional attributes to the output and diagnostic files

subroutine nc_write_model_atts(ncid, domain_id)

integer, intent(in) :: ncid      ! netCDF file identifier
integer, intent(in) :: domain_id

integer :: NGridDimID
integer :: tlonVarID, tlatVarID

character(len=256) :: filename

if ( .not. module_initialized ) call static_init_model

! put file into define mode.

write(filename,*) 'ncid', ncid

call nc_begin_define_mode(ncid)

call nc_add_global_creation_time(ncid)

call nc_add_global_creation_time(ncid)

call nc_add_global_attribute(ncid, "model_source", source )
call nc_add_global_attribute(ncid, "model", "CICE-SCM")

call nc_check(nf90_def_dim(ncid, name='ni', &
       len = Nx, dimid = NGridDimID),'nc_write_model_atts', 'ni def_dim '//trim(filename))

call nc_check(nf90_def_var(ncid,name='TLON', xtype=nf90_real, &
              dimids=(/ NGridDimID /), varid=tlonVarID),&
              'nc_write_model_atts', 'TLON def_var '//trim(filename))
call nc_check(nf90_def_var(ncid,name='TLAT', xtype=nf90_real, &
              dimids=(/ NGridDimID /), varid=tlatVarID),&
              'nc_write_model_atts', 'TLAT def_var '//trim(filename))

call nc_end_define_mode(ncid)

call nc_check(nf90_put_var(ncid, tlonVarID, TLON ), &
             'nc_write_model_atts', 'TLON put_var '//trim(filename))
call nc_check(nf90_put_var(ncid, tlatVarID, TLAT ), &
             'nc_write_model_atts', 'TLAT put_var '//trim(filename))

! Flush the buffer and leave netCDF file open
call nc_synchronize_file(ncid)

end subroutine nc_write_model_atts

!------------------------------------------------------------------
! Given a kind, return what variable number it is

function get_varid_from_kind(dart_kind)

integer, intent(in) :: dart_kind
integer             :: get_varid_from_kind

integer :: i

do i = 1, get_num_variables(domain_id)
   if (dart_kind == state_kinds_list(i)) then
      get_varid_from_kind = i
      return
   endif
end do

if (debug > 1) then
   write(string1, *) 'Kind ', dart_kind, ' not found in state vector'
   write(string2, *) 'AKA ', get_name_for_quantity(dart_kind), ' not found in state vector'
   call error_handler(E_MSG, 'get_varid_from_kind', string1, source, text2=string2)
endif

get_varid_from_kind = -1

end function get_varid_from_kind

!------------------------------------------------------------------

subroutine verify_state_variables(state_variables, ngood, table, kind_list, update_var)

character(len=*),  intent(inout) :: state_variables(:)
integer,           intent(out) :: ngood
character(len=*),  intent(out) :: table(:,:)
integer,           intent(out) :: kind_list(:)   ! kind number
logical, optional, intent(out) :: update_var(:) ! logical update

integer :: nrows, i
character(len=NF90_MAX_NAME) :: varname, dartstr, update

if ( .not. module_initialized ) call static_init_model

nrows = size(table,1)

ngood = 0

! Including this check as well because the conditionals below do not pass and cause the
! code to error out when the model_state_variables namelist entry is completely empty
! (model_state_variables = '')
if ( state_variables(1) == ' ' ) then ! no model_state_variables namelist entry provided
   string1 = 'model_nml:model_state_variables not specified'
   call error_handler(E_ERR, 'verify_state_variables', string1, source)
endif

MyLoop : do i = 1, nrows

   varname = trim(state_variables(3*i -2))
   dartstr = trim(state_variables(3*i -1))
   update  = trim(state_variables(3*i   ))

   call to_upper(update)

   table(i,1) = trim(varname)
   table(i,2) = trim(dartstr)
   table(i,3) = trim(update)

   if ( table(i,1) == ' ' .and. table(i,2) == ' ' .and. table(i,3) == ' ') exit MyLoop

   if ( table(i,1) == ' ' .or. table(i,2) == ' ' .or. table(i,3) == ' ' ) then
      string1 = 'model_nml:model_state_variables not fully specified'
      call error_handler(E_ERR, 'verify_state_variables', string1, source)
   endif

   ! Make sure DART kind is valid
   kind_list(i) = get_index_for_quantity(dartstr)
   if( kind_list(i)  < 0 ) then
      write(string1,'(''there is no obs_kind <'',a,''> in obs_kind_mod.f90'')') trim(dartstr)
      call error_handler(E_ERR, 'verify_state_variables', string1, source)
   endif

   ! Make sure the update variable has a valid name
   if ( present(update_var) )then
      SELECT CASE (update)
         CASE ('UPDATE')
            update_var(i) = .true.
         CASE ('NO_COPY_BACK')
            update_var(i) = .false.
         CASE DEFAULT
            write(string1,'(A)')  'only UPDATE or NO_COPY_BACK supported in model_state_variable namelist'
            write(string2,'(6A)') 'you provided : ', trim(varname), ', ', trim(dartstr), ', ', trim(update)
            call error_handler(E_ERR, 'verify_state_variables', string1, source, text2=string2)
      END SELECT
   endif

   ! Record the contents of the DART state vector
   if (do_output()) then
      write(string1,'(A,I2,6A)') 'variable ',i,' is ',trim(varname), ', ', trim(dartstr), ', ', trim(update)
      call error_handler(E_MSG, 'verify_state_variables', string1, source)
   endif

   ngood = ngood + 1

enddo MyLoop

end subroutine verify_state_variables

!------------------------------------------------------------------
! Given a DART location (referred to as "base") and a set of
! candidate locations & kinds (locs, loc_qtys/indx), returns the
! subset close to the base, their indices, and their distances to
! the base

subroutine get_close_state(filt_gc, base_loc, base_type, locs, loc_qtys, loc_indx, &
                         num_close, close_indices, distances, state_handle)

type(get_close_type), intent(in)    :: filt_gc
type(location_type),  intent(inout) :: base_loc
integer,              intent(in)    :: base_type
type(location_type),  intent(inout) :: locs(:)
integer,              intent(in)    :: loc_qtys(:)
integer(i8),          intent(in)    :: loc_indx(:)
integer,              intent(out)   :: num_close
integer,              intent(out)   :: close_indices(:)
real(r8),             intent(out), optional :: distances(:)
type(ensemble_type),  intent(in),  optional :: state_handle

integer :: t_ind, k

! Initialize variables to missing status
num_close = 0
close_indices = -99
if (present(distances)) distances(:) = 1.0e9   !something big and positive (far away)

! Get all the potentially close obs but no dist (optional argument dist(:)
! is not present) This way, we are decreasing the number of distance
! computations that will follow.  This is a horizontal-distance operation and
! we don't need to have the relevant vertical coordinate information yet
! (for obs).
call loc_get_close_state(filt_gc, base_loc, base_type, locs, loc_qtys, loc_indx, &
                         num_close, close_indices)

! Loop over potentially close subset of obs priors or state variables
if (present(distances)) then
   do k = 1, num_close

      t_ind = close_indices(k)

      ! if dry land, leave original 1e9 value.  otherwise, compute real dist.
      distances(k) = get_dist(base_loc,      locs(t_ind), &
                                 base_type, loc_qtys(t_ind))
   enddo
endif

end subroutine get_close_state

!------------------------------------------------------------------

function read_model_time(filename)

character(len=256) :: filename
type(time_type) :: read_model_time

integer :: ncid         ! netcdf file id
integer :: nyr      , & ! year number, in cice restart
           month    , & ! month number, 1 to 12, in cice restart
           mday     , & ! day of the month, in cice restart
           sec          ! elapsed seconds into date, in cice restart
integer :: hour     , & ! hour of the day, needed for dart set_date
           minute   , & ! minute of the hour, needed for dart set_date
           secthismin

if ( .not. module_initialized ) call static_init_model

if ( .not. file_exist(filename) ) then
   write(string1,*) 'cannot open file ', trim(filename),' for reading.'
   call error_handler(E_ERR, 'read_model_time', string1, source)
endif

call nc_check( nf90_open(trim(filename), NF90_NOWRITE, ncid), &
                  'read_model_time', 'open '//trim(filename))
call nc_check( nf90_get_att(ncid, NF90_GLOBAL, 'nyr'  , nyr), &
                  'read_model_time', 'get_att nyr')
call nc_check( nf90_get_att(ncid, NF90_GLOBAL, 'month' , month), &
                  'read_model_time', 'get_att month')
call nc_check( nf90_get_att(ncid, NF90_GLOBAL, 'mday'   , mday), &
                  'read_model_time', 'get_att mday')
call nc_check( nf90_get_att(ncid, NF90_GLOBAL, 'sec', sec), &
                  'read_model_time', 'get_att sec')

if (nyr == 0) then
  call error_handler(E_ERR, 'read_model_time', &
                     'A model time with year 0 is not supported', & 
                     source)
endif

hour       = int(sec/3600)
minute     = int((sec-hour*3600)/60)
secthismin = int(sec-hour*3600-minute*60)

read_model_time = set_date(nyr, month, mday, hour, minute, secthismin)

end function read_model_time

!------------------------------------------------------------------

subroutine write_model_time(ncid, model_time, adv_to_time)

integer,         intent(in)           :: ncid
type(time_type), intent(in)           :: model_time
type(time_type), intent(in), optional :: adv_to_time

character(len=16), parameter :: routine = 'write_model_time'

integer :: iyear, imonth, iday, ihour, imin, isec
integer :: seconds

if ( .not. module_initialized ) call static_init_model

if (present(adv_to_time)) then
   call get_date(adv_to_time, iyear, imonth, iday, ihour, imin, isec)
   write(string1,*)'CICE/DART not configured to advance CICE.'
   write(string2,*)'called with optional advance_to_time of'
   write(string3,'(i4.4,5(1x,i2.2))')iyear,imonth,iday,ihour,imin, isec
   call error_handler(E_ERR, routine, string1, source, text2=string2, &
                      text3=string3)
endif

call get_date(model_time, iyear, imonth, iday, ihour, imin, isec)

seconds = (ihour*60 + imin)*60 + isec

call nc_begin_define_mode(ncid)
call nc_add_global_attribute(ncid, 'nyr'   , iyear)
call nc_add_global_attribute(ncid, 'month' , imonth)
call nc_add_global_attribute(ncid, 'mday'  , iday)
call nc_add_global_attribute(ncid, 'sec'   , seconds)
call nc_end_define_mode(ncid)

end subroutine write_model_time

!-----------------------------------------------------------------
! Check which surface temperature state variable is in restart

subroutine check_sfctemp_var(flag)

logical, intent(inout) :: flag

if (any(variable_table(:,1)=='Tsfc')) then
  flag = .true.
else
  flag = .false.
endif

end subroutine check_sfctemp_var

!===================================================================
! End of model_mod
!===================================================================

end module model_mod