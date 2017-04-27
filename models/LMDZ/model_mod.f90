! DART software - Copyright UCAR. This open source software is provided
! by UCAR, "as is", without charge, subject to all terms of use at
! http://www.image.ucar.edu/DAReS/DART/DART_download

module model_mod

!--------------------------------------------------------------------------------------------
!                Assimilation interface for LMDZ  model
!-------------------------------------------------------------------------------------------
!---------------- m o d u l e   i n f o r m a t i o n --------------------------------------
!-------------------------------------------------------------------------------------------


use           netcdf
use        types_mod, only : r8, MISSING_I,MISSING_R8,  gravity_const => gravity
use time_manager_mod, only : time_type, set_time,set_date,get_date, set_calendar_type, operator(+)
use    utilities_mod, only : register_module, error_handler, nc_check,                           &
                             E_ERR, E_MSG, nmlfileunit, do_output, do_nml_file, do_nml_term,     &
                             find_namelist_in_file, check_namelist_read,logfileunit,file_exist,  &
                             get_unit
use mpi_utilities_mod, only : my_task_id, task_count
!-------------------------------------------------------------------------
use     location_mod,  only : location_type, get_location, set_location, query_location,         &
                              LocationDims, LocationName, LocationLName,horiz_dist_only,         &
                              vert_is_level, vert_is_pressure, vert_is_height,vert_is_surface,   &
                              VERTISUNDEF, VERTISSURFACE, VERTISLEVEL,                           &
                              VERTISPRESSURE, VERTISHEIGHT,                                      &
                              get_close_type, get_close_maxdist_init, get_close_obs_init,        &
                              get_dist,loc_get_close_obs => get_close_obs

! BEGIN DART PREPROCESS USED KINDS
use     obs_kind_mod, only : QTY_U_WIND_COMPONENT, QTY_V_WIND_COMPONENT,QTY_PRESSURE,         &
                             QTY_SURFACE_PRESSURE, QTY_TEMPERATURE,QTY_SPECIFIC_HUMIDITY,     &
                             QTY_CLOUD_LIQUID_WATER, QTY_SURFACE_ELEVATION
use   random_seq_mod, only : random_seq_type, init_random_seq, random_gaussian

! end of use statements
!==============================================================================================
!
! LMDZ global/module declarations

implicit none
private

! The first block are the 16 required interfaces.  The following block
! are additional useful interfaces that utility programs can call.
public ::                                                            &
   static_init_model, get_model_size, get_model_time_step,           &
   pert_model_state, get_state_meta_data, model_interpolate,         &
   nc_write_model_atts, nc_write_model_vars,                         &
   init_conditions, init_time, adv_1step, end_model,                 &
   get_close_maxdist_init, get_close_obs_init, get_close_obs,        &
   ens_mean_for_model

public ::                                                            &
   data_2d_type,data_3d_type,PS,T,U,V,Q,CLDLIQ, prog_var_to_vector,  &
    vector_to_prog_var,  read_lmdz_init, read_lmdz_init_size,        &
   init_model_instance, end_model_instance, write_lmdz_init, coord_index
   


!----------------------------------------------------------------------------
!----------------------------------------------------------------------------

! version controlled file description for error handling, do not edit
character(len=128), parameter :: &
   source   = "$URL$", &
   revision = "$Revision$", &
   revdate  = "$Date$"

! Ensemble mean is used so that the same "state" will be used for the heigh calculations
! on all processors, for all ensemble members.
! This is allocated in static_init_model().
real(r8), allocatable :: ens_mean(:)

!----------------------------------------------------------------------
! Global storage for describing LMDZ  model class
!----------------------------------------------------------------------

!-----------------------------------------------------------------------
! Definition of variable types
! Values will be defined in order_state_fields
! All fields will have entries in the TYPE_xD corresponding to their orders
! in state_names_Xd.  The explicitly named TYPE_s are for convenience
integer :: TYPE_PS   = MISSING_I,         &
           TYPE_T    = MISSING_I,         &
           TYPE_U    = MISSING_I,         &
           TYPE_V    = MISSING_I,         &
           TYPE_Q    = MISSING_I,         &
           TYPE_CLDLIQ = MISSING_I
!-----------------------------------------------------------------------

type(time_type)                  :: time_step
type(location_type), allocatable :: state_loc(:)

! temporary output
integer :: num_calced = 0, num_searched = 0

!---------------Variables and grid structure of start.nc ----------------------
!               
!
!          (PS,T,Q,CLDLIQ)               (PS,T,Q,CLDLIQ)
!               *--------------(U)--------------*
!               !                               !
!               !                               !
!               !                               !
!               !                               !
!  rlatv------>(V)              *              (V)       
! (~slat)       !          (rlonu,rlatv)        !       
!               !                               !
!               !                               !
!               !                               !
!  rlatu-__---->*-------------(U)---------------*
!  (~lat)  (PS,T,Q,CLDLIQ)     ^          (PS,T,Q,CLDLIQ)
!               ^              !
!               !              !
!               !              !
!           rlonv (~lon)    rlonu (~slon)
!---------------------------------------------------------------------------------            
! (PS,PHIS)(lon,lat)  ! (T,Q,CLDLIQ)(lon,lat,sigs),! U(slon,lat,sigs),! V(lon,slat,sigs)
!--------------------------------------------------------------------------------
! PHIS  : Surface Geopotential 
! PS    : Surface Pressure (Pa)
! T     : Air Temparature  (K)
! U     : Zonal Wind Comp  (m/s)
! V     : Meridional Wind  (m/s)
! Q     : Specific Humidity (kg/kg)   ??
! CLDLIQ: Cloud liq. water
!------------------------------------------------------------------
!------------------------------------------------------------------
! In start.nc lon~(-180,180), lat~(90,-90), level(bottom,top) but
! for  DART  it has been converted to lon~(0,360), lat~(-90, 90), level(top, bottom)
!----------------------------------------------------------------
!----------------------------------------------------------------

type grid_1D_type
    !private
    integer                      :: dim_id
    integer                      :: length
    real(r8)                     :: resolution
    real(r8), pointer            :: vals(:)
    character (len=128)          :: var_name
    integer                      :: num_atts
    character (len=128), pointer :: atts_names(:)
    character (len=128), pointer :: atts_vals(:)
end type grid_1D_type

!-------
type data_2d_type
   !private
   integer               :: var_id
   integer               :: length
   real(r8), pointer     :: vals(:,:)
   character (len=128)   :: atts_names  ! start.nc have only one attribute
   character (len=128)   :: atts_vals  
end type data_2d_type
!-------
type data_3d_type
   !private
   integer               :: var_id
   integer               :: length
   real(r8), pointer     :: vals(:,:,:)
   character (len=128)   :: atts_names  ! start.nc have only one attribute
   character (len=128)   :: atts_vals   ! start.nc have only one attribute  
end type data_3d_type

! renaming of coordinate variables in  start.nc: lon=rlonv, lat=rlatu, slon=rlonu, slat=rlatv 
! here slat & slon represents staggered coordinates.
 type(grid_1D_type)    ::  slon, lat, lon, slat, sig, sigs, ap,apm,bp, bpm, presnivs,temps 
 type(data_3d_type)    ::  T,U,V,Q,CLDLIQ
 type(data_2d_type)    ::  PHIS,PS

! end of model section
!----------------------------------------------------------------------
! Derived parameters

! make sure static init code only called once
logical :: module_initialized = .false.
!-------

type(time_type) :: Time_step_atmos
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - --
! Random sequence and init for pert_model_state
type(random_seq_type)   :: random_seq
! Variable for keeping track of which ensemble member is to be perturbed
! by pert_model_state, which is called by filter for each ensemble member
! for a cold start.
integer                 :: num_tasks
integer                 :: my_task
integer                 :: ens_member = 0
logical                 :: do_out

! common message string used by many subroutines
character(len=150)               :: msgstring, string2
!----
character (len=8), allocatable   :: cflds(:)
integer                          :: nflds         ! # fields to read
integer                          :: num_dims
integer                          :: model_size

!----- used in data conversion from [-180 180] to [0 360] Formate
integer :: botm_positive_lon_index
integer :: botm_positive_slon_index
integer :: botm_negative_lon_index
integer :: botm_negative_slon_index

! Surface pressures, used by vertical interpolation routines.
logical               :: alloc_ps=.true.              ! Flag whether to alloc space for ps[_stagr] 
real(r8), allocatable :: ps_ens_mean(:, :)            ! ps ensemble mean
real(r8), allocatable :: ps_ens_mean_stagr_lonu(:, :) ! ps used to calc P profiles & heights on grid staggered
                                                      !    East-West (i.e. at location  of  UCOV) relative to ps
real(r8), allocatable :: ps_ens_mean_stagr_latv(:, :) ! ps used to calc P profiles & heights on grid staggered
                                                      !    North-South (i.e. at location  of  VCOV)

! height
! Surface potential; used for calculation of geometric heights.
logical               :: alloc_phis=.true.     ! Flag whether to allocate space for phis[_stagr] 
real(r8), allocatable :: phis_stagr_lonu(:, :) ! surface geopotential staggered as for ps
real(r8), allocatable :: phis_stagr_latv(:, :) ! surface geopotential staggered as for ps


integer :: ii
! columns of pressure and model level heights for use in convert_vert
real(r8), allocatable :: p_col(:), model_h(:)


! array for the linking of obs_kinds (QTY_) to model field TYPE_s
! It's filled in map_kinds.The max size of QTY_ should come from obs_kind_mod
! These should be dimensioned the same size as the total of state_names_Nd.
integer, dimension(100) :: dart_to_lmdz_kinds = (/(MISSING_I,ii=1,100)/)
integer, dimension(100) :: lmdz_to_dart_kinds = (/(MISSING_I,ii=1,100)/)


!----------------------------------------------------------------------
! Namelist variables with default values follow
!-----------------------------------------------------------------------
!define names list for  model_nml in input.nml file
! Special for an experiment.  Specify one string kind e.g QTY_CLOUD_LIQUID and 
! observations of that kind will only impact other obs and state vars of that
! same kind.  All other kinds of obs and state vars will not be impacted
! by obs of this kind.  A null string means behave as normal.  Kind strings
! are limited by compilers to be 32 chars, since they are declared as params.
character(len=32) :: impact_only_same_kind = ' '
integer           :: impact_kind_index = -1

logical :: print_details = .false.
logical :: write_grads   = .false.


! output_state_vector = .true.     results in a "state-vector" netCDF file
! output_state_vector = .false.    results in a "prognostic-var" netCDF file
logical  :: output_state_vector = .false.
character(len = 128) :: model_config_file = 'start.nc'

integer, parameter :: MAX_STATE_NAMES = 100

! list of fields which this code needs to perturb because they're
! constant valued model parameters and show no spread when start_from_restart = .true.
character (len=8),dimension(MAX_STATE_NAMES) :: pert_names = (/('',ii=1,MAX_STATE_NAMES)/)
real(r8)         ,dimension(MAX_STATE_NAMES) :: pert_sd    = (/(-888888.0d0,ii=1,MAX_STATE_NAMES)/)
real(r8)         ,dimension(MAX_STATE_NAMES) :: pert_base_vals = (/(-888888.0d0,ii=1,MAX_STATE_NAMES)/)

! Specify shortest time step that the model will support
! This is limited below by CAMs fixed time step but is also impacted
! by numerical stability concerns for repeated restarting in leapfrog.
integer  :: time_step_days      = 1
integer  :: time_step_seconds   = 0


! Define location restrictions on which observations are assimilated
! (values are calculated anyway, but istatus is set to 2)
real(r8) :: max_obs_lat_degree        = 90.0_r8
real(r8) :: highest_obs_pressure_mb   = 150.0_r8
real(r8) :: highest_state_pressure_mb = 150.0_r8

! exclude some upper levs fron initial perturbations let model balance in free run
integer  :: exclude_pert_upper_levs = 1
! These are not namelist variables, but are related, and calculated from
! highest_obs_pressure_mb
real(r8) :: highest_obs_level         = MISSING_R8
real(r8) :: highest_obs_height_m      = MISSING_R8
!---- end of namelist (found in file input.nml) ----

namelist /model_nml/ model_config_file,time_step_seconds, time_step_days, write_grads, &
                     impact_only_same_kind, print_details, max_obs_lat_degree, highest_obs_pressure_mb,     &
                      highest_state_pressure_mb, pert_names,pert_sd,pert_base_vals, exclude_pert_upper_levs    


contains

!==================================================================

subroutine static_init_model()
!------------------------------------------------------------------
!
! Called to do one time initialization of the model. As examples,
! might define information about the model size or model timestep.
! In models that require pre-computed static data, for instance
! spherical harmonic weights, these would also be computed here.
! Can be a NULL INTERFACE for the simplest models.

 real(r8)        :: x_loc
 integer         :: i, j, iunit, io, ncfileid
 integer         :: max_levs, topog_lons, topog_lats
 type(time_type) :: model_time 

!------------------------------------------------------------------
! only execute this code once
if (module_initialized) return
! Print module information to log file and stdout.
call register_module(source, revision, revdate)

! setting calendar type
! calendar types listed in time_manager_mod.f90
! this information is NOT passed to LMDZ; it must be set in the LMDZ namelist

call set_calendar_type('GREGORIAN')

! Read a the namelist entry 
call find_namelist_in_file("input.nml", "model_nml", iunit)
read(iunit, nml = model_nml, iostat = io)
call check_namelist_read(iunit, io, "model_nml")

! set the printed output logical variable to reduce printed output;
! depends on whether this is being called by dart_to_cam (read ens member # from
! file 'element' )  or by filter (multiple processes, printout controlled by do_output())

if (file_exist('element')) then
   iunit = get_unit()
   open(unit = iunit, file='element', form = 'formatted')
   read(iunit,*) ens_member
   close(iunit)
   do_out = .false.
   if (ens_member == 1) do_out = .true.
else
   do_out = do_output()
   !write(*,*) 'do_out = ',do_out
   ! static_init_model is called once(?) for each task(?).
   ! There may be more or fewer ensemble members than tasks.
   ! No problem if there are fewer.
   ! In pert_model_state generate a unique ens_member from my_task and globally
   ! stored info  about previous calls to pert_model_state.
   num_tasks = task_count()
   my_task = my_task_id()
end if


! Record (write in a file) the namelist values used for the run ...
if (do_nml_file()) write(nmlfileunit, nml=model_nml)
if (do_nml_term()) write(     *     , nml=model_nml)

! Set the model minimum time step from the namelist seconds and days input
Time_step_atmos = set_time(time_step_seconds, time_step_days)

! read LMDZ 'initial' file domain info

call nc_check(nf90_open(path = trim(model_config_file), mode =nf90_nowrite , ncid = ncfileid), &
             'static_init_model', 'opening '//trim(model_config_file))

call nc_check(nf90_inquire(ncfileid, num_dims), 'read_lmdz_init_size', 'inquire num_dims')

write(*,*)'Num of coordinate var. in start.nc=',num_dims

call read_lmdz_init_size(ncfileid)

  PS%length     = lon%length*lat%length
  U%length      = slon%length*lat%length*sigs%length
  V%length      = lon%length*slat%length*sigs%length
  T%length      = lon%length*lat%length*sigs%length
  Q%length      = lon%length*lat%length*sigs%length
  CLDLIQ%length = lon%length*lat%length*sigs%length

! Lenght of state vector
  model_size =  PS%length + U%length + V%length + T%length + Q%length + CLDLIQ%length
print*, 'Model Size =',model_size

call read_lmdz_coord(ncfileid, lon,      'rlonv     ')
call read_lmdz_coord(ncfileid, lat,      'rlatu     ')
call read_lmdz_coord(ncfileid, slon,     'rlonu     ')
call read_lmdz_coord(ncfileid, slat,     'rlatv     ')
call read_lmdz_coord(ncfileid, presnivs, 'presnivs  ')
call read_lmdz_coord(ncfileid, ap   ,    'ap        ')
call read_lmdz_coord(ncfileid, bp   ,    'bp        ')
call read_lmdz_coord(ncfileid, sigs ,    'nivsigs   ')
call read_lmdz_coord(ncfileid, temps   , 'temps     ')


call change_lon_lat_lev_to_dart()
!-----------------------------------------------------------------------
!find hybrid layer coefficient at mid point of layers
allocate(apm%vals(sigs%length))
allocate(bpm%vals(sigs%length))
call hybrid_coefi_mid_layer(apm, bpm, sig%length)
!-----------------------------------------------------------------------

!----
max_levs = sig%length  !(sig%length = sigs%length + 1)

nflds = 6  ! variables in state vector
allocate (cflds(nflds))
! Order of state variables in statevector 
cflds=(/'PS      ', 'T       ', 'U       ', 'V       ','Q       ','CLDLIQ  '/)
call order_state_fields()
!----------

! start.nc have only one attribute 'title'
call nc_read_model_atts_2d( ncfileid ,'title', PS ,      'ps      ')
call nc_read_model_atts_3d( ncfileid ,'title', U  ,      'ucov    ')
call nc_read_model_atts_3d( ncfileid ,'title', V  ,      'vcov    ')
call nc_read_model_atts_3d( ncfileid ,'title', T  ,      'teta    ')
call nc_read_model_atts_3d( ncfileid ,'title', Q  ,      'H2Ov    ')
call nc_read_model_atts_3d( ncfileid ,'title', CLDLIQ  , 'H2Ol    ')
!call nc_read_model_atts( ncfileid ,'title',PHIS ,'phisinit')

! Read surface geopotential for use in vertical interpolation in height

topog_lats= lat%length
topog_lons= lon%length
allocate(PHIS%vals (topog_lons, topog_lats))


allocate (p_col(max_levs), model_h(max_levs))
allocate(ens_mean(model_size))

!--get PHIS data
call read_lmdz_horiz (ncfileid, phis , topog_lons, topog_lats, 'phisinit')

!-------
call nc_check(nf90_close(ncfileid), &
              'static_init_model', 'closing '//trim(model_config_file))

!------------------------------------------------------------------------
! arrays for the linking of obs_kinds (QTY_) to model field TYPE_s; 
! Makes an array of 'locations within the state vector'
! of  all the available obs kinds that come from obs_kind_mod.

call map_kinds()

!if (len_trim(impact_only_same_kind) > 0) then
!   impact_kind_index = get_index_for_quantity(impact_only_same_kind)
!endif

! make sure we only come through here once

module_initialized = .true.

!=======================================================================
end subroutine 



subroutine read_lmdz_init_size(ncfileid)

!=======================================================================
integer,  intent(in)  :: ncfileid

!   dim_ids(i) = i   ! Dimension ids are sequential integers on the NetCDF file.

   call nc_check(nf90_inquire_dimension(ncfileid, 2, slon%var_name , slon%length), &
                 'read_lmdz_init_size', 'inquire for '//trim(slon%var_name))
   call nc_check(nf90_inquire_dimension(ncfileid, 3, lat%var_name , lat%length), &
                 'read_lmdz_init_size', 'inquire for '//trim(lat%var_name))
   call nc_check(nf90_inquire_dimension(ncfileid, 4, lon%var_name , lon%length), &
                 'read_lmdz_init_size', 'inquire for '//trim(lon%var_name))
   call nc_check(nf90_inquire_dimension(ncfileid, 5, slat%var_name , slat%length), &
                 'read_lmdz_init_size', 'inquire for '//trim(slat%var_name))
   call nc_check(nf90_inquire_dimension(ncfileid, 6, sigs%var_name , sigs%length), &
                 'read_lmdz_init_size', 'inquire for '//trim(sigs%var_name))
   call nc_check(nf90_inquire_dimension(ncfileid, 7, sig%var_name , sig%length), &
                 'read_lmdz_init_size', 'inquire for '//trim(sig%var_name))
   !call nc_check(nf90_inquire_dimension(ncfileid, 15, temps%var_name , temps%length), &
   !              'read_lmdz_init_size', 'inquire for '//trim(temps%var_name))
    !if (print_details .and. do_out) write(*,*) 'Dims info = ',i, trim(dim_names(i)), dim_sizes(i)

    ap%length       = sig%length
    bp%length       = sig%length
    presnivs%length = sigs%length
    temps%length    = 1

    write(*,"(A10,I3)") slon%var_name , slon%length
    write(*,"(A10,I3)") lat%var_name , lat%length
    write(*,"(A10,I3)") lon%var_name , lon%length
    write(*,"(A10,I3)") slat%var_name , slat%length
    write(*,"(A10,I3)") sigs%var_name  , sigs%length
    write(*,"(A10,I3)") sig%var_name   , sig%length

end subroutine


subroutine read_lmdz_coord(ncfileid, var, cfield)
!=======================================================================
integer,            intent(in)  :: ncfileid   ! file and field IDs
character (len=8),  intent(in)  :: cfield
type(grid_1d_type), intent(out) :: var
!-----------
!----------------------------------------------------------------------
! Local workspace
integer :: i, coord_size             ! grid/array indices
integer :: ncfldid    ! file and field IDs
integer :: ncerr      ! other nc errors; don't abort
integer, dimension(nf90_max_var_dims) :: coord_dims
integer :: num_atts, keep_atts, alen
integer                                     :: att_type
!character (len=nf90_max_name)               :: att_name
!character (len=nf90_max_name), pointer      :: att_vals(:)
character (len=nf90_max_name) , pointer     :: att_names(:)
character (len=nf90_max_name), pointer      :: att_vals(:)
real(r8)                                    :: resol, resol_1, resol_n

  call nc_check(nf90_inq_varid(ncfileid,trim(cfield), ncfldid), &
                 'read_lmdz_coord', 'inq_varid '//trim(cfield))
  call nc_check( nf90_inquire_variable(ncfileid, ncfldid, dimids = coord_dims, nAtts =num_atts), &
               'read_lmdz_coord', 'inqure_variable '//trim('cfield'))

  allocate(att_names(num_atts), att_vals(num_atts))
!-------
 do i=1,num_atts

  call nc_check(nf90_inq_attname(ncfileid, ncfldid, i, att_names(i)), &
                 'read_lmdz_coord', 'inq_attname '//trim('att_name'))
  !var%atts_names(i) = att_name

  call nc_check(nf90_get_att(ncfileid, ncfldid, att_names(i),att_vals(i)), &
                    'read_lmdz_coord', 'get_att '//trim('att_vals' ))

 end do

 call init_grid_1d_instance(var, var%length, num_atts)


 var%var_name = cfield
 var%dim_id = coord_dims(1)
 var%atts_vals = att_vals 
 var%atts_names = att_names

 call nc_check(nf90_get_var(ncfileid, ncfldid, var%vals, start=(/1/) &
    ,count=(/var%length/) ), 'read_lmdz_coord' ,'get_var '//cfield)

   !PRINT*,'reading ',cfield,' using id ',ncfldid,' size ',var%length
   !WRITE(*,*) 'first, last val: ', var%vals(1),var%vals(var%length)
   !PRINT*,var%length
! radian to degree
if (cfield(1:2) == 'rl') then
   call rad_to_degree(var)
end if


if (cfield(1:2) == 'ap'.OR.cfield(1:2) == 'bp') then
   var%resolution = MISSING_R8
else
   resol_1 = var%vals(2) - var%vals(1)
   resol = 1.0_r8/resol_1
   var%resolution = resol_1

   ! Test all other coordinate spacings.  If any of them differ from the first
   ! by more than epsilon (smallest meaningful number relative to the coordinate
   ! spacings)
   ! then spacing is irregular.
   Res: do i = 3, var%length 
      resol_n = var%vals(i) - var%vals(i-1)
      if (((resol_n - resol_1) *resol) > epsilon(resol_n)) then
         var%resolution = -1.0_r8
      !print*,var%var_name , resol_n , var%resolution ;pause
         exit Res
      endif
   end do Res
endif

 deallocate(att_names, att_vals)

end subroutine


   subroutine init_grid_1d_instance(var, length, num_atts)
!=======================================================================
! subroutine init_grid_1d_instance(var)
!
! Initializes an instance of a lmdz grid variable

type(grid_1d_type), intent(out) :: var
integer,       intent(in ) :: length, num_atts

! Initialize the storage space and return
allocate( var%vals      (length))
var%vals = 0.0_r8

allocate( var%atts_names(num_atts))
allocate( var%atts_vals (num_atts))
var%num_atts = num_atts

!var%length = length


end subroutine init_grid_1d_instance


subroutine end_grid_1d_instance(var)
!=======================================================================
! subroutine end_grid_1d_instance(var)
!
! Ends an instance of a lmdz  grid_1d variable

 type(grid_1d_type), intent(inout) :: var

 deallocate(var%vals, var%atts_names, var%atts_vals)

end subroutine end_grid_1d_instance



   subroutine nc_read_model_atts_3d(  ncfileid , att, var, cfield)
!=======================================================================
! Local workspace
integer :: i, nchars, ierr
integer :: ncfileid, ncfldid, ncattid, att_type
!----------------------------------------------------------------------

character (len=8),  intent(in)  :: cfield
character (len=*),  intent(in)  :: att
type(data_3d_type)              :: var
character (len=128)             :: att_vals
!----------------------------------------------------------------------

    call nc_check(nf90_inq_varid(ncfileid, cfield, ncfldid),'nc_read_model_atts', &
                 'inq_varid '//trim(cfield))
 
   ierr = nf90_inquire_attribute(ncfileid, ncfldid, trim(att), att_type, nchars, ncattid)

   if (ierr == nf90_noerr) then
      call nc_check(nf90_get_att(ncfileid, ncfldid, trim(att) ,att_vals ), &
                    'nc_read_model_atts', 'get_att '//trim(att))
      att_vals(nchars+1:128) = ' '
   else
      WRITE(*,*) ncfldid, 'NOT AVAILABLE!'
   end if

    var%atts_names = att 
    var%atts_vals  = att_vals

end subroutine nc_read_model_atts_3d


   subroutine nc_read_model_atts_2d(  ncfileid , att, var, cfield)
!=======================================================================
! Local workspace
integer :: i, nchars, ierr
integer :: ncfileid, ncfldid, ncattid, att_type
!----------------------------------------------------------------------

character (len=8),  intent(in)  :: cfield
character (len=*),  intent(in)  :: att
type(data_2d_type)              :: var
character (len=128)             :: att_vals
!----------------------------------------------------------------------

    call nc_check(nf90_inq_varid(ncfileid, cfield,ncfldid),'nc_read_model_atts', &
                 'inq_varid '//trim(cfield))

   ierr = nf90_inquire_attribute(ncfileid, ncfldid, trim(att), att_type, nchars,ncattid)

   if (ierr == nf90_noerr) then
      call nc_check(nf90_get_att(ncfileid, ncfldid, trim(att) ,att_vals ), &
                    'nc_read_model_atts', 'get_att '//trim(att))
      att_vals(nchars+1:128) = ' '
   else
      WRITE(*,*) ncfldid, 'NOT AVAILABLE!'
   end if

    var%atts_names = att
    var%atts_vals  = att_vals

end subroutine nc_read_model_atts_2d




   subroutine read_lmdz_horiz(ncfileid, var, dim1, dim2, cfield)
!======================================================
! should be called with cfield = a 2D record variable  (time,lat,lon):

implicit none
!------------------------------------------------------
integer,                         intent(in)  :: ncfileid, dim1, dim2
type(data_2d_type)                           :: var
character (len=8),               intent(in)  :: cfield

!------------------------------------------------------
integer :: ncfldid
integer :: i, j
!if (print_details .and. do_out) PRINT*,'read_lmdz_horiz; reading ',cfield


call nc_check(nf90_inq_varid(ncfileid, trim(cfield), ncfldid), &
              'read_lmdz_horiz', 'inq_varid '//trim(cfield))
call nc_check(nf90_get_var(ncfileid, ncfldid, var%vals, start=(/1,1,1/), &
           count=(/dim1, dim2, 1/)), 'read_lmdz_horiz', trim(cfield))

call convert_grid_2d_data_to_dart(dim1,dim2,botm_positive_lon_index,var%vals)

if (cfield == 'phisinit') then

 if (alloc_phis) then
     allocate ( phis_stagr_lonu( slon%length,  lat%length) )   ! at UCOV on grid C
     allocate ( phis_stagr_latv( lon%length , slat%length) )   ! at VCOV on grid C
     alloc_phis = .false. 


     do i = 1, slon%length
       do j = 1, lat%length
        phis_stagr_lonu(i, j) = 0.5 * (PHIS%vals(i, j) + PHIS%vals(i+1, j))       
       end do
     end do 

     do i = 1, lon%length
       do j = 1, slat%length
        phis_stagr_latv(i, j) = 0.5 * (PHIS%vals(i, j) + PHIS%vals(i, j+1))       
       end do
     end do 

  end if ! alloc_phis 
end if ! cfield 
end subroutine read_lmdz_horiz

 



 subroutine map_kinds()
!=======================================================================
! subroutine map_kinds()
! Borrowed from CAM

! ? Should this be a function instead; removes need to dimension obs_loc_in
! arbitrarily
!   and wastefully.  But then it's called millions of times, instead of
!   accessing an
!   array that's defined once.

! Makes an array of 'locations within the state vector'
! of  all the available obs kinds that come from obs_kind_mod. 
! The obs kind that's needed will be the index into this array,
! the corresponding value will be the position of that field (not individual
! variable) 
! within the state vector according to state_name_Xd.  
! This subroutine will be called from static_init_model, so it will not have to
! be 
! recomputed for every obs.
! Also maps the local model_mod TYPE_s onto the DART QTY_s by the same
! mechanism.

! other QTY_ possibilities are listed after the 'use obs_kind_mod' statement

!integer :: i


dart_to_lmdz_kinds(QTY_SURFACE_PRESSURE)   = TYPE_PS
dart_to_lmdz_kinds(QTY_TEMPERATURE)        = TYPE_T
dart_to_lmdz_kinds(QTY_U_WIND_COMPONENT)   = TYPE_U
dart_to_lmdz_kinds(QTY_V_WIND_COMPONENT)   = TYPE_V
dart_to_lmdz_kinds(QTY_SPECIFIC_HUMIDITY)  = TYPE_Q
dart_to_lmdz_kinds(QTY_CLOUD_LIQUID_WATER) = TYPE_CLDLIQ

if (TYPE_PS /= MISSING_I)      lmdz_to_dart_kinds(TYPE_PS)     = QTY_SURFACE_PRESSURE
if (TYPE_T /= MISSING_I)       lmdz_to_dart_kinds(TYPE_T)      = QTY_TEMPERATURE
if (TYPE_U /= MISSING_I)       lmdz_to_dart_kinds(TYPE_U)      = QTY_U_WIND_COMPONENT
if (TYPE_V /= MISSING_I)       lmdz_to_dart_kinds(TYPE_V)      = QTY_V_WIND_COMPONENT
if (TYPE_Q /= MISSING_I)       lmdz_to_dart_kinds(TYPE_Q)      = QTY_SPECIFIC_HUMIDITY
if (TYPE_CLDLIQ /= MISSING_I)  lmdz_to_dart_kinds(TYPE_CLDLIQ) = QTY_CLOUD_LIQUID_WATER 


!if (print_details .and. do_out) then
!   write(*,*) 'OBS_KIND   FIELD_TYPE'
!   do i=1,100
!      if (dart_to_lmdz_kinds(i) /= MISSING_I) write(*,'(2I8)') i,
!dart_to_lmdz_kinds(i)
!   end do
!end if

! In the future, if fields are not ordered nicely, or if users are specifying
! correspondence of obs fields with state fields, I may want code like:
! The max size of QTY_ should come from obs_kind_mod
! do i=1,state_num_3d
!    if (state_names_3d(i)(1:1) == 'T' .and. &
!        QTY_TEMPERATURE <= 100) ) dart_to_lmdz_kinds(QTY_TEMPERATURE) =
!        TYPE_3D(i)
! end do 

return

end subroutine map_kinds


   function get_model_size()
!=======================================================================

integer :: get_model_size

get_model_size = model_size

end function get_model_size

   subroutine init_model_instance(PS_local, T_local, U_local, V_local, Q_local, CLDLIQ_local)
!=======================================================================
! subroutine init_model_instance(var)
!
! Initializes an instance of a lmdz model state variable

! Initialize the storage space and return
! keep some others name of variabls
    type(data_2D_type), intent(in) :: PS_local
    type(data_3D_type), intent(in) :: U_local,V_local,T_local,Q_local,CLDLIQ_local

    allocate(PS_local%vals(lon%length,lat%length))
    allocate(T_local%vals(lon%length,lat%length,sigs%length))
    allocate(U_local%vals(slon%length,lat%length,sigs%length))
    allocate(V_local%vals(lon%length,slat%length,sigs%length))
    allocate(Q_local%vals(lon%length,lat%length,sigs%length))
    allocate(CLDLIQ_local%vals(lon%length,lat%length,sigs%length))

end subroutine init_model_instance


   subroutine end_model_instance(PS_local, T_local, U_local, V_local, Q_local, CLDLIQ_local)
!=======================================================================
! Initializes an instance of a lmdz model state variable

! Initialize the storage space and return
! keep some others name of variables
    type(data_2D_type), intent(in) :: PS_local
    type(data_3D_type), intent(in) ::  U_local,V_local,T_local,Q_local,CLDLIQ_local

    deallocate(PS_local%vals)
    deallocate(U_local%vals)
    deallocate(V_local%vals)
    deallocate(T_local%vals)
    deallocate(Q_local%vals)
    deallocate(CLDLIQ_local%vals)

end subroutine end_model_instance



   subroutine read_lmdz_init(file_name, model_time)
!=======================================================================

character(len = *),        intent(in)    :: file_name
type(time_type), optional, intent(out)   :: model_time

! Local workspace
integer :: ncfileid, ncfldid, dimid, varid, dimlen

integer(kind=4),save :: iyear, ayear, imonth, iday, ihour, imin, isec, rem
integer, allocatable, dimension(:) :: datetmp, datesec

integer :: i, j
real(r8) , allocatable ::  tmp_2d(:,:), tmp_3d(:,:,:)

!----------------------------------------------------------------------
call nc_check(nf90_open(path = trim(file_name), mode = nf90_nowrite, ncid = ncfileid), &
      'read_lmdz_init', 'opening '//trim(file_name))

!----PS--
call nc_check(nf90_inq_varid(ncfileid,'ps', ncfldid), &
                'read_lmdz_init', 'inq_varid '//trim('ps'))

call nc_check(nf90_get_var(ncfileid, ncfldid, PS%vals ,start=(/1,1,1/)  &
                           ,count=(/lon%length,lat%length, 1/) ), &
                           'read_lmdz_init', 'get_var '//trim('ps'))
call convert_grid_2d_data_to_dart(lon%length,lat%length,botm_positive_lon_index,PS%vals)
!----T--

call nc_check(nf90_inq_varid(ncfileid,'teta', ncfldid), &
                'read_lmdz_init', 'inq_varid '//trim('teta'))
call nc_check(nf90_get_var(ncfileid, ncfldid, T%vals ,start=(/1,1,1,1/)  &
                           ,count=(/lon%length,lat%length, sigs%length,1/) ), &
                           'read_lmdz_init', 'get_var '//trim('teta'))

call convert_grid_3d_data_to_dart(lon%length,lat%length,sigs%length,botm_positive_lon_index,T%vals)

!-----Q--
call nc_check(nf90_inq_varid(ncfileid,'H2Ov', ncfldid), &
                'read_lmdz_init', 'inq_varid '//trim('H2Ov'))
call nc_check(nf90_get_var(ncfileid, ncfldid, Q%vals ,start=(/1,1,1,1/)  &
                           ,count=(/lon%length,lat%length, sigs%length,1/) ), &
                           'read_lmdz_init', 'get_var '//trim('H2Ov'))
call convert_grid_3d_data_to_dart(lon%length,lat%length,sigs%length,botm_positive_lon_index,Q%vals)

!-----CLDLIQ--
call nc_check(nf90_inq_varid(ncfileid,'H2Ol', ncfldid), &
                'read_lmdz_init', 'inq_varid '//trim('H2Ol'))
call nc_check(nf90_get_var(ncfileid, ncfldid, CLDLIQ%vals ,start=(/1,1,1,1/)  &
                           ,count=(/lon%length,lat%length, sigs%length,1/) ), &
                           'read_lmdz_init', 'get_var '//trim('H2Ol'))
call convert_grid_3d_data_to_dart(lon%length,lat%length,sigs%length,botm_positive_lon_index,CLDLIQ%vals)
!----U--

call nc_check(nf90_inq_varid(ncfileid,'ucov', ncfldid), &
                'read_lmdz_init', 'inq_varid '//trim('ucov'))
call nc_check(nf90_get_var(ncfileid, ncfldid, U%vals ,start=(/1,1,1,1/)  &
                           ,count=(/slon%length,lat%length, sigs%length,1/) ), &
                           'read_lmdz_init', 'get_var '//trim('ucov'))
call convert_grid_3d_data_to_dart(slon%length,lat%length,sigs%length,botm_positive_slon_index,U%vals)

!----V--

call nc_check(nf90_inq_varid(ncfileid,'vcov', ncfldid), &
                'read_lmdz_init', 'inq_varid '//trim('vcov'))
call nc_check(nf90_get_var(ncfileid, ncfldid, V%vals ,start=(/1,1,1,1/)  &
                           ,count=(/lon%length,slat%length, sigs%length,1/) ), &
                           'read_lmdz_init', 'get_var '//trim('vcov'))
call convert_grid_3d_data_to_dart(lon%length,slat%length,sigs%length,botm_positive_lon_index,V%vals)

!---
call nc_check(nf90_close(ncfileid), 'read_lmdz_init', 'closing '//trim(file_name))


! Read the time of the current state.
! extarct date information from  time unit atrribute 
!e.g temps:unit='days since 2009-05-13 00:00:00' 

if (present( model_time)) then
 !conversion form Charector to inter 
 read(temps%atts_vals(2)(12:15),'(I4)') iyear
 read(temps%atts_vals(2)(17:18),'(I2)') imonth
 read(temps%atts_vals(2)(20:21),'(I2)') iday
 read(temps%atts_vals(2)(23:24),'(I2)') ihour
 read(temps%atts_vals(2)(26:27),'(I2)') imin
 read(temps%atts_vals(2)(29:30),'(I2)') isec

 !temps%vals is fraction of day for 6 hourly run its value is {0,.25,0.5,.75,1}

 model_time = set_date(iyear,imonth,iday,ihour,imin,isec)
 
 PRINT*,"*********************************************************"
 PRINT*,'MODEL RESTART BaseTime',iyear,imonth,iday,ihour,imin,isec
 PRINT*,"*********************************************************"

 ! the time is in days - multiply by 86400 to get total seconds
 model_time = model_time + set_time(int(24 * 60 * 60 * temps%vals(1)))
 call get_date(model_time, iyear,imonth,iday,ihour,imin,isec)

 PRINT*,"*********************************************************"
 PRINT*,'MODEL RESTART Time    ',iyear,imonth,iday,ihour,imin,isec
 PRINT*,"*********************************************************"

end if

end subroutine read_lmdz_init




  subroutine  prog_var_to_vector(x,PS_local,T_local,U_local,V_local,Q_local,CLDLIQ_local)
!=======================================================================
! subroutine prog_var_to_vector(var, x)
real(r8),          intent(out) :: x(:)
type(data_2D_type), intent(in) :: PS_local
type(data_3D_type), intent(in) :: U_local,V_local,T_local,Q_local,CLDLIQ_local

integer :: i, j, k,indx
!---

indx=0

! Don't change the order of do loops 
!Rivision num: 2


!---PS---
do j=1,lat%length
  do i=1,lon%length
      indx = indx + 1
      x(indx)=PS_local%vals(i,j)
!      write(101,*)indx,x(indx),i,j
  end do     
end do


!---T---
do j = 1, lat%length
  do i = 1, lon%length
    do k = 1, sigs%length
      indx = indx + 1
      x(indx) = T_local%vals(i,j,k)
!      write(102,*)indx,x(indx),i,j,k
    end do     
  end do     
end do

!---U---
do j = 1, lat%length
  do i = 1, slon%length
    do k = 1, sigs%length
      indx = indx + 1
      x(indx) = U_local%vals(i,j,k)
      !write(103,*)indx, x(indx), i, j, k
    end do     
  end do     
end do

!---V---
do j = 1, slat%length
  do i = 1, lon%length
    do k = 1, sigs%length
      indx = indx + 1
      x(indx) = V_local%vals(i,j,k)
      !write(104,*)indx,x(indx),i,j,k
    end do     
  end do     
end do

!-----Q
do j = 1, lat%length
  do i = 1, lon%length
    do k = 1, sigs%length
      indx = indx + 1
      x(indx) = Q_local%vals(i,j,k)
      !write(105,*)indx,x(indx),i,j,k
    end do     
  end do     
end do

!-----CLDLIQ
do j = 1, lat%length
  do i = 1, lon%length
    do k = 1, sigs%length
      indx = indx + 1
      x(indx) = CLDLIQ_local%vals(i,j,k)
      !write(106,*)indx,x(indx),i,j,k
    end do     
  end do     
end do

!do i=1,indx
!  write(111,*)i,x(i)
!end do

end subroutine prog_var_to_vector 


 subroutine vector_to_prog_var(x,PS_local,T_local,U_local,V_local,Q_local,CLDLIQ_local)
!=======================================================================
real(r8),           intent(in)  :: x(model_size)
type(data_2D_type), intent(out) :: PS_local
type(data_3D_type), intent(out) :: U_local,V_local,T_local,Q_local,CLDLIQ_local

integer :: i, j, k,indx
!---

indx=0

! Don't change the order of do loops 
!Rivision num: 2


!---PS---
do j=1,lat%length
  do i=1,lon%length
      indx = indx + 1
      PS_local%vals(i,j) = x(indx)
      !write(101,*)indx,x(indx),i,j
  end do
end do


!---T---
do j = 1, lat%length
  do i = 1, lon%length
    do k = 1, sigs%length
      indx = indx + 1
      T_local%vals(i,j,k) = x(indx)
      !write(102,*)indx,x(indx),i,j,k,lat%vals(j),lon%vals(i)
    end do
  end do
end do

!---U---
do j = 1, lat%length
  do i = 1, slon%length
    do k = 1, sigs%length
      indx = indx + 1
      U_local%vals(i,j,k) = x(indx)
      !write(103,*)indx, x(indx), i, j, k
    end do
  end do
end do

!---V---
do j = 1, slat%length
  do i = 1, lon%length
    do k = 1, sigs%length
      indx = indx + 1
      V_local%vals(i,j,k) =  x(indx)
      !write(104,*)indx,x(indx),i,j,k
    end do
  end do
end do

!-----Q
do j = 1, lat%length
  do i = 1, lon%length
    do k = 1, sigs%length
      indx = indx + 1
      Q_local%vals(i,j,k) = x(indx)
      !write(105,*)indx,x(indx),i,j,k
    end do
  end do
end do

!-----CLDLIQ
do j = 1, lat%length
  do i = 1, lon%length
    do k = 1, sigs%length
      indx = indx + 1
      CLDLIQ_local%vals(i,j,k) = x(indx)
      !write(106,*)indx,x(indx),i,j,k
    end do
  end do
end do

!do i=1,indx
!  write(111,*)i,x(i)
!end do

end subroutine vector_to_prog_var



   subroutine adv_1step(x, Time)
!=======================================================================
! subroutine adv_1step(x, Time)
!

real(r8), intent(inout) :: x(:)

! Time is needed for more general models like this; need to add in to
! low-order models
type(time_type), intent(in) :: Time

! make it an error by default; comment these calls out to actually
! test assimilations with null advance.

call error_handler(E_ERR,'adv_1step', &
                  'LMDZ model cannot be called as a subroutine; async cannot =0', &
                  source, revision, revdate)

end subroutine adv_1step


   function get_model_time_step()
!=======================================================================
! function get_model_time_step()
!
! Returns the the time step of the model. In the long run should be repalced
! by a more general routine that returns details of a general time-stepping
! capability.

type(time_type) :: get_model_time_step

! Time_step_atmos is global static storage
get_model_time_step =  Time_step_atmos

end function get_model_time_step


  subroutine end_model()
!=======================================================================
! subroutine end_model()
deallocate(ens_mean)

! Deallocate other variables?

end subroutine end_model


   subroutine init_conditions(x)
!=======================================================================
! subroutine init_conditions(x)
!
! Reads in restart initial conditions  -- noop for LMDZ

real(r8), intent(inout) :: x(:)

call error_handler(E_ERR,"init_conditions", &
                  "WARNING!!  LMDZ model has no built-in default state", &
                  source, revision, revdate, &
                  text2="cannot run with 'start_from_restart = .false.'", &
                  text3="use 'lmdz_to_dart' to create a LMDZ state vector file")

end subroutine init_conditions

   subroutine init_time(time)
!=======================================================================
! subroutine init_time(time)
!
! For now returns value of Time_init which is set in initialization routines.

type(time_type), intent(out) :: time
! WARNING: CURRENTLY SET TO 0
time = set_time(0, 0)
end subroutine init_time


 subroutine hybrid_coefi_mid_layer(hyam, hybm, max_levs)
!=======================================================================
!compute hybrid levels coefficient at mid layer from 'ap' and 'bp'
 integer              :: i,j
 type(grid_1D_type)   :: hyam,hybm
 integer , intent(in) :: max_levs
 
  do i=1,max_levs-1

    hyam%vals(i) = 0.5 * (ap%vals(i) + ap%vals(i+1)) 
    hybm%vals(i) = 0.5 * (bp%vals(i) + bp%vals(i+1)) 

  end do
 
end subroutine hybrid_coefi_mid_layer


 subroutine plevs_lmdz (pres_surf, num_levs, pres_mid )
!==============lmdz=========================================================
! Find Pressures at model-layers for given surface pressure at one grid point
! and apm,bpm hybrid coefficient

!-----------------------------------------------------------------------
real(r8), intent(in)  :: pres_surf        ! Surface pressure (pascals)
integer,  intent(in)  :: num_levs
real(r8), intent(out) :: pres_mid(:)   ! Pressure at model levels
!-----------------------------------------------------------------------

!---------------------------Local workspace-----------------------------
  integer k       
!-----------------------------------------------------------------------
!
! Set midpoint pressures and layer thicknesses
!
! coef
do k=1,num_levs
   pres_mid(k) = apm%vals(k) + bpm%vals(k)*pres_surf
end do

return

end subroutine plevs_lmdz


   function index_from_grid(lev_ind, lon_ind, lat_ind, ifld)
!=======================================================================
!Function to generate the state vector index corresponding to the grid location
!and state variable
!------------------

integer, intent(in) :: lev_ind, lon_ind, lat_ind, ifld

integer :: index_from_grid
!-----

index_from_grid = 0

select case (ifld)

 case(1)   !PS
  index_from_grid = (lat_ind-1)*lon%length + lon_ind  
 case(2)   !T
  index_from_grid = PS%length + (lat_ind -1 )*lon%length*sigs%length +  & 
                                 (lon_ind-1)*sigs%length + lev_ind 
 case(3)   !U
  index_from_grid = PS%length + T%length + (lat_ind -1 )*slon%length*sigs%length + & 
                                 (lon_ind-1)*sigs%length + lev_ind
 case(4)   !V
  index_from_grid = PS%length + T%length + U%length + (lat_ind -1 )*lon%length*sigs%length + &
                                 (lon_ind-1)*sigs%length + lev_ind
 case(5)   !Q 
  index_from_grid = PS%length + T%length + U%length + V%length + (lat_ind -1 )*lon%length*sigs%length + &
                                 (lon_ind-1)*sigs%length + lev_ind
 case(6)   !CLDLIQ
  index_from_grid = PS%length + T%length + U%length + V%length + Q%length + & 
                                 (lat_ind -1 )*lon%length*sigs%length + (lon_ind-1)*sigs%length + lev_ind
end select 

if(ifld>6.OR.ifld<1) stop 'stoping in index_from_grid : variable not found for ifld'
!print*,'PS',PS%length
!print*,'T',T%length
!print*,'U',U%length
!print*,'V',V%length
!print*,'Q',Q%length
end function index_from_grid


   subroutine coord_val(dim_name, indx, lon_val, lat_val, lev_val)
!==========================================================================
! given the name of the coordinate to be searched and the index into that array, 
! returns the coordinate value  in either lon_val, lat_val, or lev_val.
! All 3 _val arguments are present so that this routine can return the value
! in the coordinate that the calling routine wants it to be, and
! search/placement doesn't
! have to be done there.
! Used by get_state_meta_data and model_interpolate.

character (len=8), intent(in)    :: dim_name
integer,           intent(in)    :: indx
real(r8),          intent(inout) :: lon_val, lat_val, lev_val

if (dim_name == 'slon    ') lon_val = slon%vals(indx)
if (dim_name == 'lat     ') lat_val = lat%vals(indx)
if (dim_name == 'lon     ') lon_val = lon%vals(indx)
if (dim_name == 'slat    ') lat_val = slat%vals(indx)
if (dim_name == 'presnivs') lev_val = presnivs%vals(indx)

return

end subroutine coord_val


   subroutine gravity(xlat,alt,galt)
!=============================================================
! This subroutine computes the Earth's gravity at any altitude
! and latitude.  The model assumes the Earth is an oblate 
! spheriod rotating at a the Earth's spin rate.  The model
! was taken from "Geophysical Geodesy, Kurt Lambeck, 1988".
!
!  input:    xlat, latitude in radians
!            alt,  altitude above the reference ellipsiod, km
!  output:   galt, gravity at the given lat and alt, cm/sec
!
! Compute acceleration due to the Earth's gravity at any latitude/altitude
! author     Bill Schreiner   5/95
! ------------------------------------------------------------------------
! Borrowed from CAM

  implicit none
  real (r8) :: xmu, ae, f, w, xm, f2, f4, ge, g, galt, xlat,alt
!
      xmu = 398600.4415_r8       ! km^3/s^2
      ae = 6378.1363_r8          ! km
      f = 1.0_r8/298.2564_r8
      w = 7.292115d-05          ! rad/s
      xm = 0.003468_r8           !
!     f2 = -f + 5.0/2.0*xm - 17.0/14.0*f*xm + 15.0/4.0*xm**2
!     f4 = -f**2/2.0 + 5.0/2.0*f*xm
      f2 = 5.3481622134089D-03
      f4 = 2.3448248012911D-05
!
! compute gravity at the equator, km/s2
!
      ge = xmu/ae**2/(1.0_r8 - f + 1.5_r8*xm - 15.0/14.0*xm*f)
!
! compute gravity at any latitude, km/s2
!
      g = ge*(1.0_r8 + f2*(dsin(xlat))**2 - 1.0/4.0*f4*(dsin(2.0*xlat))**2)
!
! compute gravity at any latitude and at any height, km/s2
!
      galt = g - 2.0*ge*alt/ae*(1.0 + f + xm + (-3.0*f + 5.0/2.0*xm)*  &
                             (dsin(xlat))**2) + 3.0*ge*alt**2/ae**2

!liu     galt = galt*1.0d5              ! convert from km/s2 to cm/s2
!
end subroutine gravity

 
  function gph2gmh (h, lat)
 !======================================================
! Barrowed from CAM

      implicit none

      real(r8) ::  h, lat, gph2gmh
      real(r8) ::  be, ae, pi, G, g0, r0, latr

      be = 6356.7516_r8             ! min earth radius, km
      ae = 6378.1363_r8             ! max earth radius, km

      pi = 3.14159265358979_r8
      latr = lat * (pi/180.0_r8)           ! in radians

  ! These are the equations for g0 and r0 in Chris Rocken's paper.
  ! I am not using them because they imply a standard ellipsoid which
  ! is quite different from our standard ellipsoid. D. Hunt 10/28/96
  !G = 0.0098
  !g0 = 0.001 * 9.780356 * (1+0.0052885 * (sin(lat))**2 - 5.9e-6 *
  !(sin(2*lat))**2)
  !r0 = (2*g0)/(3.085462e-6 + 2.27e-9 * cos(2*lat) - 2e-12*cos(4*lat))

      G = 0.00980665_r8          ! WMO reference g value, km/s**2, at 45.542N(S)

      g0 = 0.0_r8
      call gravity (latr, 0.0_r8, g0)
! liu    g0 = g0 * 0.00001_r8             ! convert to km/s**2

! compute local earth's radius using ellipse equation
!
      r0 = dsqrt ( ae**2 * dcos(latr)**2 + be**2 * dsin(latr)**2)

!     if (h.eq.-999.0_r8) then
!        z = -999.0_r8
!     else 
! Compute altitude above sea level

         gph2gmh = (r0 * h) / (((g0*r0)/G) - h)
!     end if

end function gph2gmh




   subroutine get_val(val, x, lon_index, lat_index, level, obs_kind, istatus)
!=======================================================================
! Extract the value of field at a specified location from the DART state vector

real(r8), intent(out) :: val
real(r8), intent(in) :: x(:)
integer, intent(in) :: lon_index, lat_index, level, obs_kind
integer, intent(out) :: istatus

integer :: indx, field_type

! No errors to start with
istatus = 0

field_type = dart_to_lmdz_kinds(obs_kind)

if (field_type <= 0 .or. field_type > 6) then
   istatus = 1
   val = 0.0_r8
   return
end if
indx = index_from_grid(level, lon_index, lat_index, field_type)
val  = x(indx)
return

end subroutine get_val


   subroutine set_ps_ens_mean_arrays(vec)
 !======================================================
real(r8), intent(in)   :: vec(:)
integer :: ind, i, j, fld_index , ps_ens_mean_length 

!$Rivision num : 1

if (alloc_ps) then

  allocate ( ps_ens_mean ( lon%length, lat%length ) )
  allocate ( ps_ens_mean_stagr_lonu ( slon%length, lat%length ))
  allocate ( ps_ens_mean_stagr_latv ( lon%length, slat%length ))

  alloc_ps = .false.

! fill ps_ens_mean arrays from state vector 'vec' 
  ! for PS fld_index=1
  fld_index=1
  ind = index_from_grid(1,1,1,fld_index) -1

    do j = 1, lat%length
      do i = 1, lon%length
         ind = ind + 1
         ps_ens_mean(i, j) = vec(ind)
      end do
    end do


! ps_ens_mean at stagered lat & lon grid,  check for correctnesss !!!
    do i = 1, slon%length
      do j = 1, lat%length
        ps_ens_mean_stagr_lonu(i, j) = 0.5 * ( ps_ens_mean (i, j) + ps_ens_mean (i + 1, j))
      end do
    end do

    do i = 1, lon%length
      do j = 1, slat%length
        ps_ens_mean_stagr_latv(i, j) = 0.5 * ( ps_ens_mean (i, j) + ps_ens_mean (i, j + 1))
      end do
    end do

end if  ! (alloc_ps)

end subroutine set_ps_ens_mean_arrays

    subroutine coord_index(dim_name, val, indx, other_indx)
!==========================================================================
! subroutine coord_index(dim_name, indx, val, indx, other_indx)

! Given the name of the coordinate to be searched and the value, 
! Returns the index of the closest coordinate value.  
! Optionally returns the next closest index too, which may be < or > the
! closest.
! Used by get_state_meta_data.

character (len=8), intent(in)  :: dim_name
real(r8),          intent(in)  :: val
integer,           intent(out) :: indx
integer, optional, intent(out) :: other_indx

real(r8), pointer :: coord(:)
real(r8)          :: diff_upper, diff_lower, val_local, resol
integer           :: coord_len, i

val_local = val

if (dim_name == 'lon     ') then
   coord     => lon%vals
   coord_len =  lon%length
   resol     =  lon%resolution
elseif (dim_name == 'lat     ') then
   coord     => lat%vals
   coord_len =  lat%length
   resol     =  lat%resolution
elseif (dim_name == 'slon    ') then
   coord     => slon%vals
   coord_len =  slon%length
   resol     =  slon%resolution
elseif (dim_name == 'slat    ') then
   coord     => slat%vals
   coord_len =  slat%length
   resol     =  slat%resolution
elseif (dim_name == 'presnivs') then
   coord     => presnivs%vals
   coord_len =  presnivs%length
   resol     =  presnivs%resolution
else
   ! should not happen; fatal error.
   write(msgstring, *) 'unexpected dim_name, ', trim(dim_name)
   call error_handler(E_ERR, 'coord_index', msgstring,source,revision,revdate)
end if

! further check?  for blunders check that coord(1) - val is smaller than
! coord(2) - coord(1), etc.
! Assumes that coordinates are monotonic; not true for hyam, hyai.  But we don't
! reference them.
! The first 2 if blocks work for latitudes and levels.  Longitudes must be
! handled in the calling
!    routine.
if (val_local <= coord(1)) then
   indx = 1
   if (present(other_indx)) other_indx = 1
   return
elseif (val_local >= coord(coord_len)) then
   indx = coord_len
   if (present(other_indx)) other_indx = coord_len
   return
else
   if (resol > 0.0_r8) then
      ! temp output
      num_calced = num_calced + 1
      ! regularly spaced; calculate the index
      indx = NINT((val_local - coord(1))/resol) + 1
      if (present(other_indx)) then
         if (val_local > coord(indx)) then
            other_indx = indx + 1
         else
            other_indx = indx - 1
         endif
      endif
   else
      ! temp output
      num_searched = num_searched + 1
      ! IRregularly spaced; search for the index
      ! Replace with a binary search?
      do i=1, coord_len - 1
         diff_upper = coord(i+1) - val_local
         if (diff_upper >= 0.0_r8) then
            diff_lower = val_local - coord(i)
            ! Alway return the closer coord point in the first (non-optional)
            ! argument
            if (diff_upper > diff_lower) then
               indx = i
               if (present(other_indx)) other_indx = i + 1
            else
               indx = i + 1
               if (present(other_indx)) other_indx = i
            end if
            return
         end if
      end do
   endif
end if

end subroutine coord_index


   subroutine ens_mean_for_model(filter_ens_mean)
!======================================================

real(r8), intent(in) :: filter_ens_mean(:)

 ens_mean = filter_ens_mean
 ! Fill ps_ens_mean, ps_ens_mean_stagr_lxx if not filled yet.
 ! WATCH OUT that it's not still filled with something other than ens_mean
 call set_ps_ens_mean_arrays(ens_mean)


end subroutine ens_mean_for_model


  subroutine get_state_meta_data(index_in, location, var_kind)
  !subroutine get_state_meta_data(index_in, index_1, index_2, index_3,var_kind)
!=======================================================================
! subroutine get_state_meta_data(index_in, location, var_kind, set_loc)
!
! Given an integer index into the state vector structure, returns the
! associated location. 
! The location may have components that are MISSING_R8 values, since some fields
! don't have locations in all three dimensions, i.e. PS has no vertical level,
!The which_vert should take care of the vertical coordinate (it will be ignored),
!  but the others will require more interesting  fixes. 

integer,             intent(in)  :: index_in
type(location_type), intent(out) :: location
integer, optional,   intent(out) :: var_kind

integer  :: which_vert
integer  :: i, indx, index_1, index_2, index_3, nfld
integer  :: box, slice
real(r8) :: lon_val, lat_val, lev_val

! In order to find what variable this is, and its location, I must subtract the
! individual 
! components of the state vector, since they may have varying sizes.
! Save the original index.   
! index_in will be < 0 if it's an identity obs (called from convert_vert)
indx = abs(index_in)
which_vert = MISSING_I
index_1 = 0; index_2 = 0; index_3 = 0; nfld = 0
lon_val = MISSING_R8; lat_val = MISSING_R8; lev_val = MISSING_R8

!-----look for PS (2D)---
nfld = nfld + 1                              ! lon_index = index_2
                                             ! lat_index = index_3 
   slice = PS%length                         ! lev_index = index_1
   if (indx > slice ) then
      indx = indx - slice
   else
      ! We've found the desired field. 
      ! Now find lat ,lon indx if called by
      ! assim_tools_mod:filter_assim

         ! # second dimension rows to subtract off; temporary value for index_2
         index_3 = (indx -1) / lon%length
         index_2 = indx - (index_3 * lon%length )
         call coord_val('lon     ', index_2, lon_val, lat_val, lev_val)
          
         ! index_2 of the variable in question is 1 more than the # subtracted
         ! of to get index_1
         index_3 = index_3 + 1

         call coord_val('lat     ', index_3, lon_val, lat_val, lev_val)
         

         which_vert = VERTISSURFACE  ! or =-1 ! look in location_mod  for this variable
 
      goto 10
   end if

!------T (3D)---------
   nfld = nfld + 1
   box = T%length 
   if (indx > box ) then
      indx = indx - box
   else
      ! We've found the desired field. 
      ! Now find lat, lon, lev of indx if called by
      ! assim_tools_mod:filter_assim

         ! # of (first x second dimension) slices to subtract off in order to
         ! find the current slice
         ! debug; try printing out all this index info when there's just 1 obs
         ! to be processed.
         slice = lon%length * sigs%length
         index_3 = (indx -1) / slice         ! temporary value used to find index_2 and index_1
         index_2 = (indx -1 - (index_3 * slice)) / sigs%length    ! same for  index_2 to find index_1
         index_1 = (indx - (index_3 * slice) - (index_2 *  sigs%length))

   !     Don't return the value of lev(index_1), which = (ap+1000*bp)
   !     Return index_1 as the vertical location
         lev_val = real(index_1)

         ! index_2 of the variable in question is one more than the number
         ! subtracted off to get index_1
         index_2 = index_2 + 1
         call coord_val('lon     ', index_2, lon_val, lat_val, lev_val)

         ! index_3 of the variable in question is one more than the number
         ! subtracted off to get index_1
         index_3 = index_3 + 1
         call coord_val('lat     ', index_3, lon_val, lat_val, lev_val)

         which_vert = VERTISLEVEL  !or, =1 

      goto 10
   end if
!------U (3D)---------
   nfld = nfld + 1
   box = U%length
   if (indx > box ) then
      indx = indx - box
   else
         slice = slon%length * sigs%length
         index_3 = (indx -1) / slice        
         index_2 = (indx -1 - (index_3 * slice)) / sigs%length   
         index_1 = (indx - (index_3 * slice) - (index_2 *  sigs%length))

         lev_val = real(index_1)

         index_2 = index_2 + 1
         call coord_val('slon    ', index_2, lon_val, lat_val, lev_val)

         index_3 = index_3 + 1
         call coord_val('lat     ', index_3, lon_val, lat_val, lev_val)

         which_vert = VERTISLEVEL  !or, =1 

      goto 10
   end if

!------V (3D)---------
   nfld = nfld + 1
   box = V%length
   if (indx > box ) then
      indx = indx - box
   else
         slice = lon%length * sigs%length
         index_3 = (indx -1) / slice         
         index_2 = (indx -1 - (index_3 * slice)) / sigs%length    
         index_1 = (indx - (index_3 * slice) - (index_2 *  sigs%length))

         lev_val = real(index_1)

         index_2 = index_2 + 1
         call coord_val('lon     ', index_2, lon_val, lat_val, lev_val)

         index_3 = index_3 + 1
         call coord_val('slat    ', index_3, lon_val, lat_val, lev_val)

         which_vert = VERTISLEVEL  !or, =1 

      goto 10
   end if

!------Q (3D)---------
   nfld = nfld + 1
   box = Q%length
   if (indx > box ) then
      indx = indx - box
   else
         slice = lon%length * sigs%length
         index_3 = (indx -1) / slice    
         index_2 = (indx -1 - (index_3 * slice)) / sigs%length    
         index_1 = (indx - (index_3 * slice) - (index_2 *  sigs%length))
         lev_val = real(index_1)
         index_2 = index_2 + 1
         call coord_val('lon     ', index_2, lon_val, lat_val, lev_val)
         index_3 = index_3 + 1
         call coord_val('lat     ', index_3, lon_val, lat_val, lev_val)

         which_vert = VERTISLEVEL  !or, =1 

      goto 10
   end if

!------CLDLIQ (3D)---------
   nfld = nfld + 1
   box = CLDLIQ%length
   if (indx > box ) then
      indx = indx - box
   else
         slice = lon%length * sigs%length
         index_3 = (indx -1) / slice         
         index_2 = (indx -1 - (index_3 * slice)) / sigs%length    
         index_1 = (indx - (index_3 * slice) - (index_2 *  sigs%length))
         lev_val = real(index_1)
         index_2 = index_2 + 1
         call coord_val('lon     ', index_2, lon_val, lat_val, lev_val)
         index_3 = index_3 + 1
         call coord_val('lat     ', index_3, lon_val, lat_val, lev_val)

         which_vert = VERTISLEVEL  !or, =1 

      goto 10
   end if

!-------
10 continue

! This will malfunction for fields that are filled with MISSING_r8 for lat_val
! or lon_val.
if (lon_val == MISSING_r8 .or. lat_val == MISSING_r8 ) then
   write(msgstring, *) 'Field ',cflds(nfld),' has no lon or lat dimension.  ', &
         'What should be specified for it in the call to location?'
   call error_handler(E_ERR, 'get_state_meta_data', msgstring, source, revision, revdate)
else
 !write(*,*)'index_in,indx, nfld =',index_in,indx, nfld
 !write(*,* )'set_meta_data',index_2,index_3,index_1,lon_val,lat_val,lev_val
   location = set_location(lon_val, lat_val, lev_val, which_vert)
endif

! If the type is wanted, return it
if (present(var_kind)) then
   if (index_in < 0) then ! (Identity Obs)
      ! used by convert_vert which wants the LMDZ field index, not the DART QTY_ 
      var_kind = nfld
   else if (index_in > 0) then
      ! used by call from assim_tools_mod:filter_assim, which wants the DART
      ! QTY_
      var_kind = lmdz_to_dart_kinds(nfld)
   end if
end if



end subroutine get_state_meta_data


   subroutine get_val_level(val, x, lon_index, lat_index, level, obs_kind, istatus)
!=======================================================================
!   subroutine get_val_level(val, x, lon_index, lat_index, level, obs_kind,
!   istatus)
!
! Gets the value on level for variable obs_kind
! at lon_index, lat_index horizontal grid point
! This version excludes observations below lowest level and above
! highest level.

! $Rivision num : 2

real(r8), intent(out) :: val
real(r8), intent(in)  :: x(:)
integer,  intent(in)  :: lon_index, lat_index, level, obs_kind
integer,  intent(out) :: istatus

integer               :: vstatus, num_levs
real(r8)              :: p_surf, threshold

istatus = 0
vstatus = 0
! This assumes that all variables are defined on model levels, not on interface
! levels.
num_levs = sigs%length 
if (level > num_levs .or. level < 1) then
   ! Exclude obs below the model's lowest level and above the highest level
   istatus = 1
   val = MISSING_R8
else
   if (highest_obs_level == MISSING_R8) then  ! Enter only in first call
      ! To do this completely right; p_surf would depend on whether obs_kind was
      ! on a staggered grid
      ! and highest_obs_level would be recalculated for each location passed in.
      ! Since the hybrid coord system is pure pressure at the top levels, I'll
      ! ignore these for now.
      ! Give highest_obs_pressure_mb in mb determine highest_obs_level in model level
      ! here lev_index=-1 is ignored to get value of p_surf because it don't have levels.
      call get_val(p_surf, x, lon_index, lat_index, -1, QTY_SURFACE_PRESSURE,vstatus) 
      call plevs_lmdz (p_surf, num_levs, p_col)
      highest_obs_level = 1.0_r8
      threshold = highest_obs_pressure_mb*100.0_r8  ! conversion in Pa
      do while ((p_col(nint(highest_obs_level))) < threshold )  ! level: from top to bottom
         highest_obs_level = highest_obs_level + 1.0_r8
      end do
   end if

   ! exclude level above highest_obs_level (level=1 means top)
   if (level < nint(highest_obs_level)) then   ! level_index: top=1,bottom=last_index
      ! Exclude from assimilation the obs above a user specified level
      ! but still calculate the expected obs.
      istatus = 2
   else
      istatus = 0
   end if

   if (obs_kind == QTY_PRESSURE) then
      ! Can't get the value from get_val because 3d pressure is not a model
      ! variable.
      ! Can calculate it from ps.
      ! ps is on A-grid, so no need to check for staggered grids
      call get_val(p_surf, x, lon_index, lat_index, -1, QTY_SURFACE_PRESSURE, vstatus)
      if (vstatus > 0) then
         val = MISSING_R8
         istatus = 1
         return
      end if
      ! Next, get the values on the levels for this ps
      call plevs_lmdz (p_surf, num_levs, p_col)

      val = p_col(level)
   else
       call get_val(val, x, lon_index, lat_index, level, obs_kind, vstatus)
   end if

   if (vstatus /= 0) then
      istatus = 1
      val = MISSING_R8
   end if
end if

end subroutine get_val_level



   subroutine get_val_pressure(val, x, lon_index, lat_index, pressure, obs_kind, istatus)
!=======================================================================
!
! Gets the vertically interpolated value on pressure for variable obs_kind
! at lon_index, lat_index horizontal grid point
!
! This version excludes observations below lowest level pressure and above
! highest level pressure.


real(r8), intent(out) :: val
real(r8), intent(in)  :: x(:), pressure
integer,  intent(in)  :: lon_index, lat_index, obs_kind
integer,  intent(out) :: istatus

real(r8)              :: bot_val, top_val, p_surf, frac
integer               :: top_lev, bot_lev, i, vstatus, num_levs

! No errors to start with
istatus = 0
vstatus = 0

! Need to get the surface pressure at this point. 
! Find out whether the observed field is a staggered in LMDZ.
! This version excludes observations below lowest level pressure and above
! highest level pressure.

if  (obs_kind == QTY_U_WIND_COMPONENT ) then
      p_surf = ps_ens_mean_stagr_lonu(lon_index, lat_index)

elseif (obs_kind == QTY_V_WIND_COMPONENT ) then
     p_surf = ps_ens_mean_stagr_latv(lon_index, lat_index)
else
   ! A-grid ps can be retrieved from state vector, which was used to define ps
   ! on entry to   model_interpolate.
   p_surf = ps_ens_mean(lon_index, lat_index)
end if
 

num_levs = sigs%length
call plevs_lmdz(p_surf, num_levs, p_col)

if (pressure <= p_col(1) .or. pressure >= p_col(num_levs)) then
   ! Exclude obs below the model's lowest level and above the highest level
   ! We *could* possibly use ps and p(num_levs) to interpolate for points below
   ! the lowest level.
   istatus = 1
   val = MISSING_R8
else
   ! Interpolate in vertical to get two bounding levels
   if (pressure < highest_obs_pressure_mb * 100.0_r8) then
      ! Exclude from assimilation the obs above a user specified level
      ! but still calculate the expected obs.
      istatus = 2
   else
      istatus = 0
   end if
   ! Search down through pressures
   do i = 2, num_levs
      if (pressure < p_col(i)) then
         top_lev = i -1
         bot_lev = i
         frac = (p_col(i) - pressure) / (p_col(i) - p_col(i - 1))
         goto 21
      end if
   end do

   21 continue

   ! Pobs
   if (obs_kind == QTY_PRESSURE) then
      ! can't get pressure on levels from state vector; get it from p_col  instead 
      bot_val = p_col(bot_lev)
      top_val = p_col(top_lev)
      val = (1.0_r8 - frac) * bot_val + frac * top_val

   else
   ! Pobs end

                        call get_val(bot_val, x, lon_index, lat_index, bot_lev, obs_kind, vstatus)
      if (vstatus == 0) call get_val(top_val, x, lon_index, lat_index, top_lev, obs_kind, vstatus)
      if (vstatus == 0) then
         val = (1.0_r8 - frac) * bot_val + frac * top_val
      else
         istatus = 1
         val = MISSING_R8
      end if
   end if

end if

end subroutine 

   subroutine get_interp_prof (prof, vec, num_levs, lon_index, lat_index, stagr_lon, stagr_lat, &
                            kind_lmdz, vstatus)
!=====================================================================

real(r8), intent(out) :: prof(num_levs)
integer,  intent(out) :: vstatus
! type(model_type), intent(in) :: state
real(r8), intent(in)  :: vec(:)
integer,  intent(in)  :: kind_lmdz, num_levs, lon_index, lat_index
logical,  intent(in)  :: stagr_lon, stagr_lat

real(r8)  :: var(num_levs,0:1,0:1), weight
integer   :: k, lon_index_local, vec_ind

! So far only called from model_heights to get T and Q profiles for Tv
! calculation.
! If the point at which we need T and Q is a staggered US or VS point, 
! then T and Q must be interpolated to that point.
! lon_index and lat_index are the indices of the point at which we need the
! profiles,
!   not necessarily one of the points where the profiles exist.
! stagr_xx tells where to look for profiles and what interpolating to do.

! 
!  Grid index for var(num_levs,0:1,0:1)
!
!   (0,1)-------------(1,1)
!    !                  !
!    !                  !                (180)
!   (slat)     .        !                  ^
!    !                  !                  !
!    !                  !                  !
!   (0,0)-----(slon)---(1,0)   (0)--->-- index -----> (360)
!                                          !
! $Rivision num : 1                        !
!                                        (-180)
! 

vstatus = 0
var = 0._r8

vec_ind = index_from_grid(1,lon_index   ,lat_index   , kind_lmdz)
var(1:num_levs,0,0) = vec(vec_ind:vec_ind +num_levs -1)
weight = 1.0_r8

if (stagr_lat) then   ! interpolation at VCOV point
   ! Find the profile to the north
   vec_ind = index_from_grid(1,lon_index   , lat_index +1, kind_lmdz)
   var(1:num_levs,0,1) = vec(vec_ind:vec_ind +num_levs -1)
   ! This weighting is not strictly correct for Gaussian latitudes, 
   ! but is correct for equally spaced latitudes.
   weight = weight / 2.0_r8
end if

if (stagr_lon) then  ! at UCOV
   ! Find the profile to the east
   vec_ind = index_from_grid(1,lon_index + 1 ,lat_index, kind_lmdz)
   var(1:num_levs,1,0) = vec(vec_ind : vec_ind  + num_levs -1)
   weight = weight / 2.0_r8
end if

if (stagr_lon .and. stagr_lat) then  ! at  Middle of grid
   ! Find the profile to the northeast (corner of grid)
   vec_ind = index_from_grid(1,lon_index +1 ,lat_index +1, kind_lmdz)
   var(1:num_levs,1,1) = vec(vec_ind:vec_ind +num_levs -1)
end if

do k = 1,num_levs
   prof(k) = weight * (var(k,0,0) + var(k,1,0) + var(k,0,1) + var(k,1,1))
end do

end subroutine get_interp_prof


!===============================================================================
   subroutine model_heights(vec, p_surf, lon_index, lat_index, num_levs, &
                            stagr_lon, stagr_lat, model_h, istatus)
!===============================================================================
! This routine calculates geometrical height (m) at mid-layers of the LMDZ model
!
! Kevin Raeder converted to single column version 4/28/2006
!              removed longitude dimension entirely and extra arrays 10/2006

implicit none

real(r8),         intent(in)  :: vec(:)
real(r8),         intent(in)  :: p_surf
integer,          intent(in)  :: lon_index, lat_index, num_levs
logical,          intent(in)  :: stagr_lon, stagr_lat
integer,          intent(out) :: istatus

! OUTPUT: geometrical height at midlayer (m)  hui liu /03/29/2004 added.
real(r8),      intent(out) ::  model_h(num_levs)

! local variables; ps must be dimensioned as an array because dcz2 has it that
! way
real (r8):: phi(num_levs), tv(num_levs), q_prof(num_levs), t_prof(num_levs)
real (r8):: pmln(num_levs+1), hyba(num_levs+1,2), hybb(num_levs+1,2), pterm(num_levs)
real (r8):: phi_surf, ht_tmp, rd, rv, rr_factor, local_lat
integer :: k, i, vstatus

istatus = 0
vstatus = 0
! Scratch arrays

rd = 287.05_r8
rv = 461.51_r8
rr_factor = (rv/rd) - 1.0_r8

DO k = 1,num_levs
   model_h(k) = 0.0_r8
   phi(k)     = 0.0_r8
   PTERM(k)   = 0.0_r8
END DO

! copy to temporary arrays

!    All arrays except hyba, hybb are oriented top to bottom
!    modified to be consistent with CAM3.0 where hyai, hyam are top to bottom
!    H Liu, 04/05/2004

! interface=1 extra level at model bottom
k = num_levs +1
hyba(1,1) = ap%vals(k)
hybb(1,1) = bp%vals(k)
!   hyam(26) = 0 -> hyba(2,2) = 0, so it would be safe to set hyba(1,2) = 0.
!   This element is referenced below, but not ultimately used.
hyba(1,2) = 0.0_r8
hybb(1,2) = 1.0_r8
! hybX go from bottom to top; b coeffs multiply sigma, and coord is pure sigma
!      at the bottom, so hybb = 1.0 there.

! mid-points=2; note that hyXm(num_levs + 1) is not defined (= MISSING_R8)
do k = 2,num_levs +1
   i = num_levs +2 - k
   hyba(k,1) = ap%vals(i)
   hybb(k,1) = bp%vals(i)
   hyba(k,2) = apm%vals(i)
   hybb(k,2) = bpm%vals(i)
end do


if(stagr_lat) then
   phi_surf = phis_stagr_latv(lon_index, lat_index)
   local_lat = slat%vals(lat_index)

elseif (stagr_lon) then
   phi_surf = phis_stagr_lonu(lon_index, lat_index)
   local_lat = lat%vals(lat_index)
else
   phi_surf = PHIS%vals(lon_index, lat_index)
   local_lat = lat%vals(lat_index)
end if

call get_interp_prof (q_prof,vec,num_levs, lon_index, lat_index, stagr_lon, stagr_lat, &
                      TYPE_Q, vstatus)
call get_interp_prof (t_prof,vec,num_levs, lon_index, lat_index, stagr_lon, stagr_lat, &
                      TYPE_T, vstatus)
! Calculate tv for this column, for use by dcz2
if (vstatus == 0) then
   do k = 1, num_levs
      tv(k) = t_prof(k)*(1.0_r8 + rr_factor*q_prof(k))
   end do
elseif (vstatus > 0) then
   istatus = 1
   return
end if

call dcz2(p_surf, phi_surf, tv ,hyba, hybb, num_levs, pmln, pterm, phi)

do k = 1,num_levs
   ht_tmp = phi(k) * 0.001_r8        ! convert to km for following call only
   model_h(k) = gph2gmh (ht_tmp, local_lat) * 1000.0_r8  ! convert back to meter
end do


  end subroutine model_heights




! height
!=====================================================================
   subroutine dcz2(p_surf,phis0,tv,hyba,hybb,kmax,pmln, pterm,z2)
!=====================================================================
!       Purpose:
!         To compute geopotential height for a CCM2 hybrid coordinate
!         vertical slice.  Since the vertical integration matrix is a
!         function of latitude and longitude, it is not explicitly
!         computed as for sigma coordinates.  The integration algorithm
!         is derived from Boville's mods in the ibm file hybrid 1mods
!         (6/17/88).  All vertical slice arrays are oriented top to
!         bottom as in CCM2.  This field is on full model levels (aka
!         "midpoints") not half levels.
!
!       Equation references are to "Hybrid Coordinates for CCM1"
!
!----------------------------Code History-------------------------------
!       Original: Jun 25 1991  L Buja
!       Upgrades: Apr 21 1993  L Buja
!                     Truesdale reports difference from CCM2 Z2.
!                     - Make pterm a scratch array (was automatic)
!                     - Make g0=9.80616 (was 9.81).  This appears to
!                        affect only the lowest layer.
!       Change  : Feb    1999: D Shea
!                     NCL changes + error
!                     - The "Invert vertical loop" has to have
!                        the argument checked for >0.0
!-----------------------------------------------------------------------
      implicit none
!-----------------------------Parameters--------------------------------
      real(r8) :: r,g0,rbyg
      parameter (r=287.04_r8,g0=9.80616_r8,rbyg=r/g0)
!-----------------------------Arguments---------------------------------
! Input
!
! Number of vertical levels
      integer kmax

! Surface pressure           (pascals)
      real(r8) :: p_surf

! Surface geoptential
      real(r8) :: phis0

! Virtual temperature, top to bottom
      real(r8) TV(KMAX)

! Hybrid base pressure       (pascals)
      real(r8) :: HPRB

! Hybrid coord coeffs for base pressure
      real(r8) :: HYBA(KMAX+1,2)

!  All arrays except hyba, hybb are oriented top to bottom
!  ground to top, first subscript:

!  = 1 for layer interfaces 
!  = 2 for layer midpoints 
!  Lowest level is ground for both layer locations

! Hybrid coord coeffs for surf pressure (in same format as hyba)
      real(r8) ::  hybb(kmax+1,2)

! vertical slice scratch space used to
!   hold logs of midpoint pressures
      real(r8) ::  pmln(kmax+1)

! Note: These scratch vertical slices are used to improve computaional
! efficiency

! Vertical scratch space.
      real(r8)::   pterm(kmax)
! temporary
      real(r8)::   arg
!
! Output ---------------------------------
!
! Geopotential height, top to bottom
      real(r8)::   Z2(KMAX)
!
!--------------------------Local variables------------------------------
! indexes
      integer i,k,l,num
!-----------------------------------------------------------------------
!
      DATA NUM/0/

      NUM = NUM + 1
!

!       Compute intermediate quantities using scratch space

!       Invert vertical loop
!       Compute top only if top interface pressure is nonzero.
!       Implemented by setting loop limit klim
!
!       hyba, hybb are bottom to top, starting at ground.
!       pmln(i,k) is the mid-point pressure of layer k.
!       SHEA MODIFICATION

      DO K = KMAX + 1,1,-1
         i = KMAX-K+2
         ARG = HYBA(i,2) + P_surf *HYBB(i,2)
         IF (ARG.GT.0.0_r8) THEN
             PMLN(K) = DLOG(ARG)
         ELSE
             PMLN(K) = 0.0_r8
         END IF
      END DO

!       Initialize Z2 to sum of ground height and thickness of
!        top half-layer  (i.e. (phi)sfc in equation 1.14)
!       (Z2(i,1)=top  ->  Z2(i,kmax)=bottom
!       Eq 3.a.109.2  where l=K,k<K  h(k,l) = 1/2 * ln [  p(k+1) / p(k) ]

      DO K = 2,KMAX - 1
         pterm(k) = rbyg*tv(k)*0.5_r8* (pmln(k+1)-pmln(k-1))
      END DO


! 
      DO K = 1,KMAX - 1
         z2(k) = phis0/g0 + rbyg*tv(k)*0.5_r8* (PMLN(K+1)-PMLN(K))
      END DO



!       Eq 3.a.109.5  where l=K,k=K  h(k,l) = ln [ pi / (p(k)) ]

      K = KMAX
!
      z2(K) = phis0/g0 + rbyg*tv(k)* (dlog(p_surf*hybb(1,1))-pmln(k))

!       Eq 3.a.109.4  where l=K,k<K  h(k,l) = 1/2*ln[pi*pi/(p(k-1)*p(k))


! 
      do k = 1,kmax - 1
          l = kmax
          z2(k) = z2(k) + rbyg*tv(l)* (dlog(p_surf*hybb(1,1))-0.5_r8* &
                                       (pmln(l-1)+pmln(l)))
      end do

!       Add thickness of the remaining full layers
!        (i.e., integrate from ground to highest layer interface)

!       Eqs 1.14 & 3.a.109.3 where l>K, k<K
!                                h(k,l) = 1/2 * ln [ p(l+1)/p(l-1) ]
! 
      DO K = 1,KMAX - 2
          DO L = K + 1,KMAX - 1
             Z2(K) = Z2(K) + PTERM(L)
          END DO
      END DO

      RETURN

end subroutine dcz2 

   subroutine get_val_height(val, vec, lon_index, lat_index, height, obs_kind, istatus)
!=======================================================================
!
! Gets the vertically interpolated value on height for variable obs_kind
! at lon_index, lat_index horizontal grid point
!
! written by Kevin Raeder, based on code from Hui Liu 4/28/2006 and
! get_val_pressure from Jeff Anderson
!
! This version excludes observations below lowest level height and above
! highest level height.
! 

real(r8), intent(out) :: val
real(r8), intent(in)  :: vec(:), height
integer,  intent(in)  :: lon_index, lat_index, obs_kind
integer,  intent(out) :: istatus

integer  :: top_lev, bot_lev, i, vstatus, num_levs
real(r8) :: bot_val, top_val, frac
real(r8) :: p_surf
logical  :: stagr_lon, stagr_lat

! No errors to start with
istatus = 0
vstatus = 0
stagr_lon = .false.
stagr_lat = .false.

num_levs = sigs%length

if  (obs_kind == QTY_U_WIND_COMPONENT ) then
      p_surf = ps_ens_mean_stagr_lonu(lon_index, lat_index)
     stagr_lon = .true.

elseif (obs_kind == QTY_V_WIND_COMPONENT ) then
     p_surf = ps_ens_mean_stagr_latv(lon_index, lat_index)
     stagr_lat = .true.
else
   ! A-grid ps can be retrieved from state vector, which was used to define ps
   ! on entry to   model_interpolate.
   p_surf = ps_ens_mean(lon_index, lat_index)
end if

! Next, get the heights on the levels for this ps

! We want to use the new vec for each new ob on height because the state was  updated 
! for all previous obs, and we want to use the most up to date state to get the  best location.
call model_heights(vec, p_surf, lon_index, lat_index, num_levs, stagr_lon, stagr_lat, &
                   model_h, vstatus)

!print*,model_h
! debug
 write(logfileunit,'(A,6F7.0,/(10F7.0))') 'heights = ',(model_h(i),i=1,num_levs)

! Interpolate in vertical to get two bounding levels
if (height >= model_h(1) .or. height <= model_h(num_levs)) then
   ! Exclude obs below the model's lowest level and above the highest level
   istatus = 1
   val = MISSING_R8

else

      call plevs_lmdz (p_surf, num_levs, p_col)
      do i=1,num_levs
         if (p_col(i) > highest_obs_pressure_mb*100.0_r8) then
            highest_obs_height_m = model_h(i)
            go to 10
         end if
      end do

10 if (height > highest_obs_height_m ) then
      ! Exclude from assimilation the obs above a user specified level
      ! but still calculate the expected obs.
      istatus = 2
   else
      istatus = 0
   end if

   ! Search down through heights
   do i = 2, num_levs
      if (height > model_h(i)) then
         top_lev = i -1
         bot_lev = i
         frac = (model_h(i) - height      ) / &
                (model_h(i) - model_h(i-1))
! check
         goto 21
      end if
   end do

   21 continue

   ! Pobs
   if (obs_kind == QTY_PRESSURE) then
      ! Observing a pressure on a height surface sounds silly.  But for
      ! completeness:
      ! get_val_height is called for 4 different columns, which will have
      ! different p_cols for each.
      ! It's also requested by obs_def_gps_mod.

      ! Next, get the values on the levels for this ps
      ! ps is on A-grid, so no need to check for staggered grids
      call plevs_lmdz (p_surf, num_levs, p_col)

      bot_val = p_col(bot_lev)
      top_val = p_col(top_lev)
   else

   ! Pobs end
                        call get_val(bot_val, vec, lon_index, lat_index, bot_lev, obs_kind, vstatus)
      if (vstatus == 0) call get_val(top_val, vec, lon_index, lat_index, top_lev, obs_kind, vstatus)
   ! Pobs
   end if

   if (vstatus == 0) then
      val = (1.0_r8 - frac) * bot_val + frac * top_val
   else
      istatus = 1
      val = MISSING_R8
   end if
end if

end subroutine get_val_height


!=======================================================================
   subroutine convert_vert (old_array, old_which, new_array, new_which, dart_kind)
!=======================================================================
! subroutine convert_vert(old_loc, new_loc, dart_kind)
!
! Uses model information and subroutines to convert the vertical location of an  ob 
! (prior, model state variable, or actual ob) into the standard vertical
! coordinate (pressure).
! Called by model_mod:get_close_obs.
!-----Code history
! Kevin Raeder 10/26/2006 for CAM
! Modified by T. Singh on  20/01/2014 for LMDZ5
!---
integer,                intent(in)    :: dart_kind, old_which
integer,                intent(out)   :: new_which
real(r8), dimension(3), intent(in)    :: old_array
real(r8), dimension(3), intent(inout) :: new_array

integer   :: i, num_levs, top_lev, bot_lev
integer   :: lon_index, lat_index
integer   :: rank_kind, lmdz_kind, istatus
real(r8)  :: p_surf,   frac
logical   :: stagr_lon, stagr_lat
type(location_type)   :: dum_loc


! set good initial values, only differences will be changed.
stagr_lon = .false.
stagr_lat = .false.

! these should be set by the code below; it's an error if not.
lon_index       = MISSING_I
lat_index       = MISSING_I
new_array       = MISSING_R8
new_which       = MISSING_I

if (old_which == VERTISPRESSURE .or. old_which == VERTISHEIGHT  .or. &
    old_which == VERTISLEVEL    .or. old_which == VERTISSURFACE .or. &
    old_which == VERTISUNDEF   ) then
   !  proceed
else
   ! make this a fatal error - there should be no other options for vert.
   write(msgstring,'(''obs at '',3(F9.5,1x),I2,'' has bad vertical type'')') &
                   old_array, old_which
   call error_handler(E_ERR, 'convert_vert', msgstring,source,revision,revdate)
end if

! Find the nfld of this dart-KIND

if (dart_kind > 0) then
   ! non-identity obs
   lmdz_kind = dart_to_lmdz_kinds(dart_kind)
else if (dart_kind < 0) then
   ! identity obs; dart_kind = -1*state_vector_index
   ! Value returned in cam_kind will be the nfld value of this field, not the
   ! usual dart_kind.
   call get_state_meta_data(dart_kind, dum_loc, lmdz_kind)
end if

! Find the index of this kind within its group of same-rank fields
rank_kind = lmdz_kind

!Find closest lat_index & lat_index form given lat lat value  

if (rank_kind == TYPE_U) then  
         call coord_index('slon    ', old_array(1), lon_index)
         call coord_index('lat     ', old_array(2), lat_index)
         p_surf = ps_ens_mean_stagr_lonu(lon_index, lat_index)
         stagr_lon=.true. 
elseif (rank_kind == TYPE_V) then 
         call coord_index('lon     ', old_array(1), lon_index)
         call coord_index('slat    ', old_array(2), lat_index)
         p_surf = ps_ens_mean_stagr_latv(lon_index, lat_index)
         stagr_lat=.true.
else   ! for all others e.g. TYPE_T, TYPE_Q, TYPE_PS, TYPE_CLDLIQ
         call coord_index('lon     ', old_array(1), lon_index)
         call coord_index('lat     ', old_array(2), lat_index)
         p_surf = ps_ens_mean(lon_index, lat_index)
         stagr_lat=.true.
end if
!----------
! Need the vertical pressure structure for this column
! This routine will be called for :
!   model grid points (from get_close_obs) (just one column of the state vector
!   is the correct one),
!   expected obs (4 times from model_interpolate) (the correct column is an
!   interpolation of the 
!      surrounding 4 columns).  

! Convert vertical coordinate from one of the following to pressure.
! integer, parameter :: VERTISUNDEF    = -2 ! has no vertical location
! (undefined)
! integer, parameter :: VERTISSURFACE  = -1 ! surface value
! integer, parameter :: VERTISLEVEL    =  1 ! by level
! integer, parameter :: VERTISPRESSURE =  2 ! by pressure
! integer, parameter :: VERTISHEIGHT   =  3 ! by height

if (old_which == VERTISUNDEF) then
   new_array(3) = MISSING_R8
   new_which    = old_which

elseif (old_which == VERTISSURFACE ) then
   ! surface field; change which_vert for the distance calculation
   new_array(3) =  p_surf
   new_which    = 2

elseif (old_which == VERTISLEVEL ) then
   num_levs = sigs%length
   call plevs_lmdz(p_surf, num_levs, p_col)
   new_array(3) = p_col(nint(old_array(3)))
   new_which    = 2

elseif (old_which == VERTISHEIGHT) then

   num_levs = sigs%length
   call plevs_lmdz (p_surf, num_levs, p_col)
   call model_heights(ens_mean, p_surf, lon_index, lat_index, num_levs, stagr_lon, stagr_lat,  &
                      model_h, istatus)
 ! Search down through heights
   ! This assumes linear relationship of pressure to height over each model
   ! layer, 

   bot_lev = 2
   do while (old_array(3) <= model_h(bot_lev) .and. bot_lev <= num_levs)
      bot_lev = bot_lev + 1
   end do

   top_lev = bot_lev - 1
  ! write warning message if not found within model level heights.
   ! maybe this should return failure somehow?
   if (top_lev == 1 .and. old_array(3) > model_h(1)) then
      ! above top of model
      frac = 1.0_r8
      write(msgstring, *) 'ob height ',old_array(2),' above LMDZ levels at '  &
                          ,old_array(1) ,old_array(2) ,' for ob type',dart_kind
      call error_handler(E_MSG, 'convert_vert', msgstring,source,revision,revdate)


   else if (bot_lev <= num_levs) then ! within model levels
      frac = (old_array(3) - model_h(bot_lev)) / (model_h(top_lev) - model_h(bot_lev))

   else
      ! below bottom of model
      frac = 0.0_r8
      write(msgstring, *) 'ob height ',old_array(3),' below LMDZ levels at ' &
                          ,old_array(1) ,old_array(2) ,' for ob type',dart_kind
      call error_handler(E_MSG, 'convert_vert', msgstring,source,revision,revdate)
   endif


   new_array(3) = (1.0_r8 - frac) * p_col(bot_lev) + frac * p_col(top_lev)
   new_which    = 2

else
   write(msgstring, *) 'model which_vert = ',old_which,' not handled in convert_vert '
   call error_handler(E_ERR, 'convert_vert', msgstring,source,revision,revdate)
end if

return

  end subroutine convert_vert

   subroutine get_close_obs(gc, base_obs_loc, base_obs_kind, obs_loc, obs_kind, &
                            num_close, close_ind, dist)
!----------------------------------------------------------------------------

! get_close_obs will be getting an ob, with its location, and its horizontal
! distances  to an array of other locations (and the locations).
!  These locations were picked out based on the efficient search/box algorithm.
!    It converts vertical coordinates as needed, 
!    It tests for being above the highest_obs_pressure_mb threshold, and
!    increases the vertical distance based on height above highest_.
!    It calls location_mod/threed_sphere:get_close_obs, 
!       to which it sends this (converted) array of locations.
!    It gets back the new total distances and arrays of those locations that are
!    "close" to the base ob.
!    get_close_obs will use the ensemble average passed through new interface;
!    ens_mean_for_model  Convert the obs and/or state vertical location(s) to a standard (pressure)
!    vertical location   3D model_h would be useful here; calc once and use over and over.
!       reinstall height/lon slice code for model_heights to facilitate that
!    3D pressures also useful here; 
!       reinstall height/lon slice code for plevs_lmdz to facilitate that
!    throw away ens_mean after it's been used (or don't worry about it for now).
! 
! The kinds are available to do more sophisticated distance computations if
! needed

implicit none

type(get_close_type), intent(in)  :: gc
type(location_type),  intent(in)  :: base_obs_loc, obs_loc(:)
integer,              intent(in)  :: base_obs_kind, obs_kind(:)
integer,              intent(out) :: num_close, close_ind(:)
real(r8),             intent(out) :: dist(:)
integer                :: k, t_ind
integer                :: base_which, local_base_which, obs_which,local_obs_which
real(r8), dimension(3) :: base_array, local_base_array, obs_array, local_obs_array
real(r8)               :: increment, threshold, thresh_wght
type(location_type)    :: local_base_obs_loc, local_obs_loc

!If base_obs vert type is not pressure; convert it to pressure

   base_array = get_location(base_obs_loc)


base_which = nint(query_location(base_obs_loc))


if (base_which == VERTISPRESSURE) then
   local_base_obs_loc = base_obs_loc
   local_base_array   = get_location(base_obs_loc)  ! needed in num_close loop
   local_base_which   = base_which
else
   base_array = get_location(base_obs_loc)

   call convert_vert(base_array, base_which, local_base_array, local_base_which, base_obs_kind)
   local_base_obs_loc = set_location(base_array(1), base_array(2), local_base_array(3), &
                                     local_base_which)
end if


!print*,'local_base_array...TKKKKK',local_base_array ;stop

!! DEBUG: comment this in if you want to bypass the top level damping code
!below.
!call loc_get_close_obs(gc, base_obs_loc, base_obs_kind, obs_loc, obs_kind, &
!                       num_close, close_ind, dist)
!return

! Get all the potentially close obs but no dist (optional argument dist(:) is
! not present)

call loc_get_close_obs(gc, local_base_obs_loc, base_obs_kind, obs_loc, obs_kind, &
                       num_close, close_ind)

threshold = highest_state_pressure_mb *100._r8
if (threshold > 0.0_r8) thresh_wght = 1._r8/(threshold * threshold)


do k = 1, num_close

   t_ind = close_ind(k)
   obs_array = get_location(obs_loc(t_ind))
   obs_which = nint(query_location(obs_loc(t_ind)))
   
   if (obs_which == VERTISPRESSURE ) then
      ! put the vertical (pressure) of the state/ob in local storage
      local_obs_array(3) = obs_array(3)
      local_obs_which    = obs_which

   else
      ! Convert vertical coordinate of obs_loc to pressure.
      ! If horiz_dist_only is true, the vertical location and which vert aren't
      ! used by get_dist, 
      ! but need to be defined for set_loc and are used in the damping section
      ! below no matter what.
      call convert_vert(obs_array, obs_which, local_obs_array, local_obs_which, obs_kind(t_ind))

      ! obs_which = -2 (VERTISUNDEF) mean this ob is vertically close to
      ! base_obs, no matter what.
      if (local_obs_array(3) == MISSING_R8) then
         local_obs_array(3) = local_base_array(3)
         local_obs_which = local_base_which
      end if
   end if

   local_obs_loc = set_location(obs_array(1), obs_array(2), local_obs_array(3), &
                                   local_obs_which)
!  nsc fri, 13mar09
!  allow a namelist specified kind string to restrict the impact of those
!  obs kinds to only other obs and state vars of the same kind.
   if ((impact_kind_index >= 0)                .and. &
       (impact_kind_index == base_obs_kind)    .and. &
       (impact_kind_index /= obs_kind(t_ind))) then
      dist(k) = 999999._r8     ! arbitrary very large distance
   else if (local_base_which == VERTISUNDEF) then
      ! The last argument, no_vert = .true., makes get_dist calculate horzontal
      ! distance only.
      dist(k) = get_dist(local_base_obs_loc, local_obs_loc, base_obs_kind, obs_kind(t_ind),.true.)
      ! Then no damping can be done since vertical distance is undefined.
      ! ? Is this routine called *both* to get model points close to a real obs,
      !   AND ob close to a model point?  I want damping in the latter case,
      !   even if ob has which_vert = VERTISUNDEF.
      !   I think that testing on local_base_which will do that.
   else

      dist(k) = get_dist(local_base_obs_loc, local_obs_loc, base_obs_kind, obs_kind(t_ind))
      ! Damp the influence of obs (below the namelist variable
      ! highest_obs_pressure_mb) 
      ! on variables above highest_state_pressure_mb.  
      ! This section could also change the distance based on the QTY_s of the
      ! base_obs and obs.

      ! dist = 0 for some for synthetic obs.
      ! Additive increase, based on height above threshold, works better than
      ! multiplicative

      ! See model_mod circa 1/1/2007 for other damping algorithms.

      increment = threshold - local_obs_array(3)
      ! This if-test handles the case where no damping is performed, i.e. 
      ! highest_state_pressure_mb = 0 and threshold = 0.
      if (increment > 0) then
          dist(k) = dist(k) + increment * increment * thresh_wght
      ! too sharp      dist(k) = dist(k) + increment / threshold
      end if
   endif

end do

end subroutine get_close_obs

!--------------------------------------------------------------------------------------

   subroutine model_interpolate(x, location, obs_type, interp_val, istatus)
!=======================================================================
!

real(r8),            intent(in) :: x(:)
type(location_type), intent(in) :: location
integer,             intent(in) :: obs_type
integer,            intent(out) :: istatus
real(r8),           intent(out) :: interp_val

integer  :: i, vstatus
real(r8) :: bot_lon, top_lon, delta_lon,lon_above,temp_lon,             &
            lon_below, lat_below, lat_above, lev_below,                 &
            lon_fract, lat_fract, vals(2, 2),  a(2),                    &
            lon_lat_lev(3), level, pressure, height

integer  :: lmdz_type,  &
            lon_ind_below, lon_ind_above, lat_ind_below, lat_ind_above, &
            num_lons
character (len=8)   :: lon_name, lat_name, lev_name


! istatus   meaning                  return expected obs?   assimilate?
! 0         obs and model are fine;  yes                    yes
! 1         fatal problem;           no                     no
! 2         exclude valid obs        yes                    no
! 3         unfamiliar obs type      no                     no

! These are fields which were observed, and will have 3d locations, but the 
! corresponding state-vector component could, conceivably, be missing one of the
! dimensions.
! The only use for such fields I have thought of is parameterization tuning.
! Such fields would not have observations associated with them.
! for now I will assume that observed fields are not missing any dimensions.
! PS is missing lev, although we have never assimilated those obs.

! model_interpolate will continue to use state passed to it;
!    recalc p_col and model_h columns as needed.
!    no need to convert to a standard vert coord; no distance calc involved.

! Start with no errors in 
istatus = 0
vstatus = 0
vals = MISSING_R8

! Always fill the ps arrays with the state vector here, since most obs and
! vertical locations  will need that info.  "Always" allows ens_mean_for_model to set ps arrays
! once for all   of the get_close_obs calls, without having to remove the ps arrays contents
! at the end  of the get_close_obs calls, which is hard to identify.

call set_ps_ens_mean_arrays(x)

! Get the observation (horizontal) position, in degrees 

lon_lat_lev = get_location(location)
!PRINT*,'++++++++++++++++++++++++++++++++++++++++++++'
!PRINT*,'Given location for interpolation',lon_lat_lev

! check whether model_mod can interpolate the requested variable
! Pressure (3d) can't be specified as a state vector field (so lmdz_type will =  MISSING_I), 
! but can be calculated for LMDZ, so obs_type = QTY_PRESSURE is acceptable.
lmdz_type = dart_to_lmdz_kinds(obs_type)

if (lmdz_type == MISSING_I .and.(obs_type .ne. QTY_PRESSURE) ) then
   istatus = 3
! should be MISSING_R8 ?
   interp_val = MISSING_R8
! check
   write(*,*) 'Wrong type of obs = ', obs_type
   return
end if

! Get lon and lat grid specs

! Set [lon,lat,lev] names to a default, which will be overwritten for variables
! in the state vector, but not for other acceptable variables (3D pressure,
! surface  elevation, ...?)
lon_name = 'lon     '
lat_name = 'lat     '
! ? How to separate the 3D from 2D 'other' variables?
!   Can't do it automatically/generically because they're not part of state
!   vector and that info isn't coming from DART.

if (obs_type .eq. QTY_PRESSURE) then
   lev_name = 'lev     '
endif

if (lmdz_type == MISSING_I .and. (obs_type .eq. QTY_PRESSURE) ) then
   ! use defaults lon_name and lat_name set above

elseif(lmdz_type==TYPE_PS)then
   lon_name = 'lon    '   
   lat_name = 'lat     '
   lev_name = 'none    '

elseif(lmdz_type == TYPE_U) then
   lon_name = 'slon    ' 
   lat_name = 'lat     '
   lev_name = 'lev     '

elseif(lmdz_type == TYPE_V) then
   lon_name = 'lon     ' 
   lat_name = 'slat    '
   lev_name = 'lev     '

elseif(lmdz_type == TYPE_T .or. lmdz_type == TYPE_Q .or. lmdz_type == TYPE_CLDLIQ ) then
   lon_name = 'lon     ' 
   lat_name = 'lat     '
   lev_name = 'lev     '

else
   istatus = 3
   interp_val = MISSING_R8
   write(*,*) 'Unexpected state type value, lmdz_type = ', lmdz_type
   return

endif


! Find the location of above& below lon index for Obs lon (lon_lat_lev(1))
! Compute bracketing lon indices
if (lon_name == 'lon     ') then
   num_lons  = lon%length
   bot_lon   = lon%vals(1)
   top_lon   = lon%vals(num_lons)
elseif (lon_name == 'slon    ') then
   num_lons  = slon%length
   bot_lon   = slon%vals(1)
   top_lon   = slon%vals(num_lons) 
end if


! LMDZ (s)lon  starts from -180 to 180 so it has been rearrange from 0 to 360
! formate but it is not necessary that (s)lon will start from 0 or end at 360.
! suppose  bot_lon = 1.5 and and top_lon = 359.5  and obs_lon= 0.5 or 359.5 then
! to interpolated it properly following two 'if blocks' needed.   
if( lon_lat_lev(1) <= bot_lon ) then
   lon_fract = ( (lon_lat_lev(1) + 360.0_r8) - top_lon ) / ( (bot_lon + 360.0_r8) - top_lon )
   lon_ind_below = num_lons 
   lon_ind_above = 1 
   lon_below = top_lon 
   lon_above = bot_lon
   !print*,'TKKK  fract 1',lon_fract
elseif(  lon_lat_lev(1) >= top_lon ) then 
   lon_fract = (lon_lat_lev(1) - top_lon ) / ( (bot_lon + 360.0_r8) - top_lon )
   lon_ind_below = num_lons 
   lon_ind_above = 1 
   lon_below = top_lon 
   lon_above = bot_lon

   !print*,'TKKK  fract 2',lon_fract
else
  call coord_index(lon_name, lon_lat_lev(1), lon_ind_above, lon_ind_below)
!
  if (lon_ind_above == lon_ind_below) then
    if (lon_ind_above == 1) then
      lon_fract = 0.0_r8
    else                     !both must be equal to the max (s)lon index
      lon_fract = 1.0_r8
    end if
  else
   if (lon_ind_above < lon_ind_below) then
      ! switch order
      i = lon_ind_above
      lon_ind_above = lon_ind_below
      lon_ind_below = i
   end if
   ! only lon_xxx is changed by these calls
   call coord_val(lon_name, lon_ind_below, lon_below, lat_below, lev_below)
   !return only lon_below
   call coord_val(lon_name, lon_ind_above, lon_above, lat_below, lev_below)

   lon_fract = (lon_lat_lev(1) - lon_below) / (lon_above - lon_below)
   !print*,'TKKK',lon_lat_lev(1),lon_ind_above,lon_ind_below,lon_below, lon_above,lon_fract
 end if
end if
! Next, compute neighboring lat rows

call coord_index(lat_name, lon_lat_lev(2), lat_ind_above, lat_ind_below)

if (lat_ind_above == lat_ind_below) then
   if (lat_ind_above == 1) then
      lat_fract = 0.0_r8
   else                     !both must be equal to the max (s)lat index
      lat_fract = 1.0_r8
   end if
else
   if (lat_ind_above < lat_ind_below) then
      ! switch order
      i = lat_ind_above
      lat_ind_above = lat_ind_below
      lat_ind_below = i
   end if
   ! only lat_xxx is changed by these calls
   call coord_val(lat_name, lat_ind_below, lon_below, lat_below, lev_below)  ! return only lat_below
   call coord_val(lat_name, lat_ind_above, lon_below, lat_above, lev_below)
   lat_fract = (lon_lat_lev(2) - lat_below) / (lat_above - lat_below)
end if
! Now, need to find the values for the four corners
! determine the vertical coordinate: model level, pressure, or height

if (vert_is_level(location)) then
 ! Case 1: model level specified in vertical
 ! Pobs
 level = lon_lat_lev(3)
      call get_val_level            &
      (vals(1, 1), x, lon_ind_below, lat_ind_below, nint(level), obs_type, vstatus)

   if (vstatus /= 1) call get_val_level   &
      (vals(1, 2), x, lon_ind_below, lat_ind_above, nint(level), obs_type, vstatus)

   if (vstatus /= 1) call get_val_level   &
      (vals(2, 1), x, lon_ind_above, lat_ind_below, nint(level), obs_type, vstatus)

   if (vstatus /= 1) call get_val_level   &
      (vals(2, 2), x, lon_ind_above, lat_ind_above, nint(level), obs_type, vstatus)
   ! Pobs end

elseif (vert_is_pressure(location)) then
   ! which_vert is pressure for this obs
   pressure = lon_lat_lev(3)
      call get_val_pressure                  &
      (vals(1,1),x,lon_ind_below,lat_ind_below,pressure,obs_type,vstatus)

   if (vstatus /= 1) call get_val_pressure   &
      (vals(1,2),x,lon_ind_below,lat_ind_above,pressure,obs_type,vstatus)

   if (vstatus /= 1) call get_val_pressure   &
      (vals(2,1),x,lon_ind_above,lat_ind_below,pressure,obs_type,vstatus)

   if (vstatus /= 1) call get_val_pressure   &
      (vals(2,2),x,lon_ind_above,lat_ind_above,pressure,obs_type,vstatus)

elseif (vert_is_height(location)) then
   ! which_vert is height for this obs
   height = lon_lat_lev(3)
      call get_val_height                  &
      (vals(1, 1), x, lon_ind_below, lat_ind_below, height, obs_type, vstatus)
   if (vstatus /= 1) call get_val_height   &
      (vals(1, 2), x, lon_ind_below, lat_ind_above, height, obs_type, vstatus)
   if (vstatus /= 1) call get_val_height   &
      (vals(2, 1), x, lon_ind_above, lat_ind_below, height, obs_type, vstatus)
   if (vstatus /= 1) call get_val_height   &
      (vals(2, 2), x, lon_ind_above, lat_ind_above, height, obs_type, vstatus)

elseif (vert_is_surface(location)) then
   ! location_mod:interactive_location asks for surface obs to have vertical
   ! coord = ps(hPa)
   ! The 'lev' argument is set to 1 because there is no level for these types,
   ! and 'lev' will be ignored.

                     call get_val(vals(1,1),x, lon_ind_below, lat_ind_below, 1, obs_type, vstatus)
   if (vstatus /= 1) call get_val(vals(1,2),x, lon_ind_below, lat_ind_above, 1, obs_type, vstatus)
   if (vstatus /= 1) call get_val(vals(2,1),x, lon_ind_above, lat_ind_below, 1, obs_type, vstatus)
   if (vstatus /= 1) call get_val(vals(2,2),x, lon_ind_above, lat_ind_above, 1, obs_type, vstatus)

else
   write(*,*) '   No vert option chosen!'

end if
   
! lat is already converted to degrees by get_location
if (abs(lon_lat_lev(2)) > max_obs_lat_degree .and. vstatus /= 1) then
   istatus = 4
else
   istatus = vstatus
end if

! indices of vals are (longitude, latitude)
 if (istatus /= 1) then
   do i = 1, 2
      a(i) = lon_fract * vals(2, i) + (1.0_r8 - lon_fract) * vals(1, i)
   end do
   interp_val = lat_fract * a(2) + (1.0_r8 - lat_fract) * a(1)
else
   interp_val = MISSING_R8
end if

!!PRINT*,'*****************************************************************'
!PRINT*,'LEV,obs_type,lmdz_type',obs_type,lmdz_type
!PRINT*,'BOX around Given location lon1,lon2',lon_ind_below ,lon_ind_above 
!PRINT*,'BOX around Given location lon1,lon2',lon_below ,lon_above 
!PRINT*,'BOX around Given location lat1,lat2',lat_ind_below ,lat_ind_above 
!PRINT*,'BOX around Given location lat1,lat2',lat_below ,lat_above 
!PRINT*,'BOX corner values',vals(1,1),vals(1,2),vals(2,1),vals(2,2)
!PRINT*,'Interpolated Value, lon_fract ',interp_val, lon_fract

end subroutine model_interpolate



  subroutine pert_model_state(state, pert_state, interf_provided)
!=======================================================================
! subroutine pert_model_state(state, pert_state, interf_provided)
!
! Perturbs a model state for generating initial ensembles
! Returning interf_provided means go ahead and do this with
! small independent perturbations.
!
! If module storage variable 'pert_sd' is positive, then we will randomly
! perturb the fields (based on the values of pert_sd) for each of
! the variable names listed in pert_names.
!
! If 'pert_sd' is negative (which includes MISSING) then the field (only one)
! listed in pert_names is set to a different constant value for each 
! ensemble member.  Those values come from 'pert_base_vals'.

! added to give each ens member a different sequence when perturbing model
! parameter fields

real(r8), intent(in)    :: state(:)
real(r8), intent(out)   :: pert_state(:)
logical,  intent(out)   :: interf_provided

type(random_seq_type)   :: random_seq
type(data_3d_type)      :: T_temp,U_temp,V_temp,Q_temp,CLDLIQ_temp
type(data_2d_type)      :: PS_temp
integer                 :: i, j, k, m, pert_fld, mode, field_num
integer                 :: dim1, dim2, dim3, member
real(r8)                :: pert_val
integer, save           :: seed = -1

! perturb model parameters for the filter_ics.
! Use the (single) state value as the "ens_mean" here.

interf_provided = .true.

call init_model_instance(PS_temp,T_temp,U_temp,V_temp,Q_temp,CLDLIQ_temp)
call vector_to_prog_var(state,PS_temp,T_temp,U_temp,V_temp,Q_temp,CLDLIQ_temp)

! seed the random number generator differently on each task, and if one task
! has more than one ensemble member, make sure the seed next time through is
! a different value.  note that seed has the 'save' attribute so it maintains
! the value from the last time this routine was called.
call init_random_seq(random_seq, seed * (my_task_id()+1))
seed = seed - 1000

print*,nflds,cflds

pert_fld = 1
do while (pert_names(pert_fld) /= '        ')

 !  WRITE(*,*) 'Perturbing ',pert_names(pert_fld),' for ens_member ',ens_member

  if (pert_names(pert_fld) == cflds(1)) then  !PS
      WRITE(*,*) '   Found match  ',cflds(1)
      dim1 = lon%length 
      dim2 = lat%length


      ! Choose mode of perturbations/resets; 
      if (pert_sd(pert_fld) <= 0.0_r8 ) then
         ! Set each ensemble member to it's own constant value, 
         ! as found in pert_base_vals(ens_member).
         ! This only works when setting a single field = to a different
         ! constant value for each ensemble member.
         ! Could add more fields by overloading pert_base_vals and 
         ! adding code to find those values.
           mode = ens_member
      else
         ! Set each *field* to it's own pert_base_val +/- pert_sd
          mode = pert_fld
      endif

      if ( pert_base_vals(mode) /= -888888.0d0 ) then
          WRITE(*,*) '   around new base value ',pert_base_vals(mode)
          ! reset base values to value provided in namelist.
          PS_temp%vals(1:dim1,1:dim2) = pert_base_vals(mode)
      endif

      WRITE(*,'(A,A8,A3,1p,2E12.4)') '    first and last state for ',cflds(1),' = ', &
                   PS_temp%vals(1,1),PS_temp%vals(dim1,dim2)

      if (pert_sd(pert_fld) > 0.0_r8 ) then
         ! randomly perturb each point around its base value.
         do j = 1, dim2
           do i = 1, dim1
            pert_val = random_gaussian(random_seq, PS_temp%vals(i,j), pert_sd(mode))
            PS_temp%vals(i,j) = pert_val
           end do
         end do
      endif

      WRITE(*,'(A,1p,2E12.4)') ' new first and last state = ', &
             PS_temp%vals(1,1),PS_temp%vals(dim1,dim2)

  elseif(pert_names(pert_fld) == cflds(2))then    ! for T

        dim1 = lon%length
        dim2 = lat%length
        dim3 = sigs%length 

        ! Choose mode of perturbations/resets; 
        if (pert_sd(pert_fld) <= 0.0_r8 ) then
             mode = ens_member
        else
             mode = pert_fld
        endif

        WRITE(*,'(A,A8,A,I4,A3,1p,2E12.4)') '    first and last state for ',cflds(2), &
              ' member ',ens_member,' = ', T_temp%vals(1,1,1),T_temp%vals(dim1,dim2,dim3)
 
        ! reset base values to value provided in namelist.
        if ( pert_base_vals(mode) /= -888888.0d0 ) then
          WRITE(*,*) '  perturbed around new base value ',pert_base_vals(mode) 
          T_temp%vals(1:dim1,1:dim2,1:dim3) = pert_base_vals(mode)
        endif

        ! randomly perturb each point around its base value.
        if (pert_sd(pert_fld) > 0.0_r8 ) then
          WRITE(*,*) 'Perturbing ',cflds(2),' of member ',ens_member, &
                          ' by st dev ',pert_sd(mode)
          do j = 1, dim2 
            do i = 1, dim1 
             do k = exclude_pert_upper_levs , dim3
               ! pert_val = rand#(O(0-1)) * standard dev  + mean
               pert_val = random_gaussian(random_seq, T_temp%vals(i,j,k), pert_sd(mode))
               T_temp%vals(i,j,k) = pert_val
             end do
            end do
          end do                                    

          WRITE(*,*) 'Last pert_val for member ',ens_member,' is ',pert_val
 
       endif  !pert_sd(pert_fld) > 0.0_r8

       WRITE(*,'(A,I4,A3,1p,2E12.4)') ' new first and last state for member ',ens_member, &
            ' = ', T_temp%vals(1,1,1), T_temp%vals(dim1,dim2,dim3)

  elseif(pert_names(pert_fld) == cflds(3))then    ! for U

        dim1 = slon%length
        dim2 = lat%length
        dim3 = sigs%length

        ! Choose mode of perturbations/resets; 
        if (pert_sd(pert_fld) <= 0.0_r8 ) then
             mode = ens_member
        else
          ! Set each *field* to it's own pert_base_val +/- pert_sd
             mode = pert_fld
        endif

        WRITE(*,'(A,A8,A,I4,A3,1p,2E12.4)') '    first and last state for ',cflds(3), &
              ' member ',ens_member,' = ', U_temp%vals(1,1,1),U_temp%vals(dim1,dim2,dim3)

        ! reset base values to value provided in namelist.
        if ( pert_base_vals(mode) /= -888888.0d0 ) then
          WRITE(*,*) '  perturbed around new base value ',pert_base_vals(mode)
          U_temp%vals(1:dim1,1:dim2,1:dim3) = pert_base_vals(mode)
        endif

        ! randomly perturb each point around its base value.
        if (pert_sd(pert_fld) > 0.0_r8 ) then
          WRITE(*,*) 'Perturbing ',cflds(3),' of member ',ens_member, &
                          ' by st dev ',pert_sd(mode)
          do j = 1, dim2
            do i = 1, dim1
             do k = exclude_pert_upper_levs , dim3
               ! pert_val = rand#(O(0-1)) * standard dev  + mean
               pert_val = random_gaussian(random_seq, U_temp%vals(i,j,k), pert_sd(mode))
               U_temp%vals(i,j,k) = pert_val
             end do
            end do
          end do

          WRITE(*,*) 'Last pert_val for member ',ens_member,' is ',pert_val

       endif  !pert_sd(pert_fld) > 0.0_r8

       WRITE(*,'(A,I4,A3,1p,2E12.4)') ' new first and last state for member ',ens_member, &
            ' = ', U_temp%vals(1,1,1), U_temp%vals(dim1,dim2,dim3)

  elseif(pert_names(pert_fld) == cflds(4))then    ! for V

        dim1 = lon%length
        dim2 = slat%length
        dim3 = sigs%length

        ! Choose mode of perturbations/resets; 
        if (pert_sd(pert_fld) <= 0.0_r8 ) then
             mode = ens_member
        else
          ! Set each *field* to it's own pert_base_val +/- pert_sd
             mode = pert_fld
        endif

        WRITE(*,'(A,A8,A,I4,A3,1p,2E12.4)') '    first and last state for ',cflds(4), &
              ' member ',ens_member,' = ', V_temp%vals(1,1,1),V_temp%vals(dim1,dim2,dim3)

        ! reset base values to value provided in namelist.
        if ( pert_base_vals(mode) /= -888888.0d0 ) then
          WRITE(*,*) '  perturbed around new base value ',pert_base_vals(mode)
          V_temp%vals(1:dim1,1:dim2,1:dim3) = pert_base_vals(mode)
        endif

        ! randomly perturb each point around its base value.
        if (pert_sd(pert_fld) > 0.0_r8 ) then
          WRITE(*,*) 'Perturbing ',cflds(4),' of member ',ens_member, &
                          ' by st dev ',pert_sd(mode)
          do j = 1, dim2
            do i = 1, dim1
             do k = exclude_pert_upper_levs, dim3
               ! pert_val = rand#(O(0-1)) * standard dev  + mean
               pert_val = random_gaussian(random_seq, V_temp%vals(i,j,k),pert_sd(mode))
               V_temp%vals(i,j,k) = pert_val
             end do
            end do
          end do

          WRITE(*,*) 'Last pert_val for member ',ens_member,' is ',pert_val

       endif  !pert_sd(pert_fld) > 0.0_r8

       WRITE(*,'(A,I4,A3,1p,2E12.4)') ' new first and last state for member',ens_member, &
            ' = ', V_temp%vals(1,1,1), V_temp%vals(dim1,dim2,dim3)

  elseif(pert_names(pert_fld) == cflds(5))then    ! for Q

        dim1 = lon%length
        dim2 = lat%length
        dim3 = sigs%length

        ! Choose mode of perturbations/resets; 
        if (pert_sd(pert_fld) <= 0.0_r8 ) then
             mode = ens_member
        else
          ! Set each *field* to it's own pert_base_val +/- pert_sd
             mode = pert_fld
        endif

        WRITE(*,'(A,A8,A,I4,A3,1p,2E12.4)') '    first and last state for ',cflds(5), &
              ' member ',ens_member,' = ', Q_temp%vals(1,1,1),Q_temp%vals(dim1,dim2,dim3)

        ! reset base values to value provided in namelist.
        if ( pert_base_vals(mode) /= -888888.0d0 ) then
          WRITE(*,*) '  perturbed around new base value ',pert_base_vals(mode)
          Q_temp%vals(1:dim1,1:dim2,1:dim3) = pert_base_vals(mode)
        endif

        ! randomly perturb each point around its base value.
        if (pert_sd(pert_fld) > 0.0_r8 ) then
          WRITE(*,*) 'Perturbing ',cflds(5),' of member ',ens_member, &
                          ' by st dev ',pert_sd(mode)
          do j = 1, dim2
            do i = 1, dim1
             do k = exclude_pert_upper_levs, dim3
               ! pert_val = rand#(O(0-1)) * standard dev  + mean
               pert_val = random_gaussian(random_seq, Q_temp%vals(i,j,k),pert_sd(mode))
               Q_temp%vals(i,j,k) = pert_val
             end do
            end do
          end do

          WRITE(*,*) 'Last pert_val for member ',ens_member,' is ',pert_val

       endif  !pert_sd(pert_fld) > 0.0_r8

       WRITE(*,'(A,I4,A3,1p,2E12.4)') ' new first and last state for member',ens_member, &
            ' = ', Q_temp%vals(1,1,1), Q_temp%vals(dim1,dim2,dim3)

 elseif(pert_names(pert_fld) == cflds(6))then    ! for CLDLIQ

        dim1 = lon%length
        dim2 = lat%length
        dim3 = sigs%length

        ! Choose mode of perturbations/resets; 
        if (pert_sd(pert_fld) <= 0.0_r8 ) then
             mode = ens_member
        else
          ! Set each *field* to it's own pert_base_val +/- pert_sd
             mode = pert_fld
        endif

        WRITE(*,'(A,A8,A,I4,A3,1p,2E12.4)') '    first and last state for ',cflds(6), &
              ' member ',ens_member,' = ', CLDLIQ_temp%vals(1,1,1), CLDLIQ_temp%vals(dim1,dim2,dim3)

        ! reset base values to value provided in namelist.
        if ( pert_base_vals(mode) /= -888888.0d0 ) then
          WRITE(*,*) '  perturbed around new base value ',pert_base_vals(mode)
          CLDLIQ_temp%vals(1:dim1,1:dim2,1:dim3) = pert_base_vals(mode)
        endif

        ! randomly perturb each point around its base value.
        if (pert_sd(pert_fld) > 0.0_r8 ) then
          WRITE(*,*) 'Perturbing ',cflds(6),' of member ',ens_member,' by st dev ',pert_sd(mode)
          do j = 1, dim2
            do i = 1, dim1
             do k = exclude_pert_upper_levs , dim3
               ! pert_val = rand#(O(0-1)) * standard dev  + mean
               pert_val = random_gaussian(random_seq, CLDLIQ_temp%vals(i,j,k),pert_sd(mode))
               CLDLIQ_temp%vals(i,j,k) = pert_val
             end do
            end do
          end do

          WRITE(*,*) 'Last pert_val for member ',ens_member,' is ',pert_val

       endif  !pert_sd(pert_fld) > 0.0_r8

       WRITE(*,'(A,I4,A3,1p,2E12.4)') ' new first and last state for member',ens_member, &
            ' = ', CLDLIQ_temp%vals(1,1,1), CLDLIQ_temp%vals(dim1,dim2,dim3)

   end if   !(pert_names(pert_fld)

pert_fld = pert_fld +1

end do  ! do while


call prog_var_to_vector(pert_state,PS_temp,T_temp,U_temp,V_temp,Q_temp,CLDLIQ_temp)
call end_model_instance(PS_temp,T_temp,U_temp,V_temp,Q_temp,CLDLIQ_temp)

end subroutine pert_model_state


   subroutine write_lmdz_coord_def(ncFileID, c_name, coord, dim_id, c_id)
!=======================================================================================

character (len=8),  intent(in)  :: c_name
integer,            intent(in)  :: ncFileID, dim_id
type(grid_1d_type), intent(in)  :: coord
integer,            intent(out) :: c_id

integer  :: i
!integer  :: nch

call nc_check(nf90_def_var(ncFileID, name=c_name, xtype=nf90_double, dimids=dim_id, &
                        varid=c_id), 'write_lmdz_coord_def', 'def_var'//trim(c_name))

!if (print_details .and. do_out) write(*,'(/A,A)') 'write_lmdz_coord_def;  ', trim(c_name)

do i=1,coord%num_atts
!   if (print_details .and. do_out) then
!!      nch = len_trim(coord%atts_vals(i))
!!                 i,trim(coord%atts_names(i)),' ', coord%atts_vals(i)(1:nch)
!      write(*,*) '   i, att_name, att_val', &
!                 i,trim(coord%atts_names(i)),' ', trim(coord%atts_vals(i))
!   endif
   call nc_check(nf90_put_att(ncFileID, c_id, coord%atts_names(i), coord%atts_vals(i)), &
                 'write_lmdz_coord_def', 'put_att '//trim(coord%atts_names(i)))
end do

return

end subroutine write_lmdz_coord_def


   subroutine write_lmdz_init(file_name, PS_local,T_local,U_local,V_local,Q_local,  &
                                                          CLDLIQ_local, model_time)
!=======================================================================

! write LMDZ 'initial' file fields that have been updated

character (len = *), intent(in)           :: file_name
type(data_2d_type),  intent(inout)        :: PS_local 
type(data_3D_type),  intent(inout)        :: T_local,U_local,V_local,Q_local,CLDLIQ_local
type(time_type),     intent(in), optional :: model_time

integer               :: i, k, n, m, ifld, ncfileid, ncfldid, dim1, dim2, dim3
integer               :: iyear, imonth, iday, ihour, imin, isec
integer               :: dimid, dimlen, varid
integer, allocatable, dimension(:) :: datetmp, datesec
real(r8), allocatable :: temp_3d(:,:,:), temp_2d(:,:)
!character*30 unites
character(len=30) unites
integer status

integer :: xtype, len, attnum

! Read LMDZ 'initial' file domain info
call nc_check(nf90_open(path = trim(file_name), mode = nf90_write, ncid = ncfileid), &
           'write_lmdz_init', 'opening '//trim(file_name))

! The temporary arrays into which fields are read are dimensioned by the largest
! values of  the sizes of the dimensions listed in coord_RANKd
! debug
! ! ! This may not work for writing fields that have smaller sizes!
!     Unless array section is explicitly given in array argument

!------PS----
   ifld = 1
   dim1 = lon%length
   dim2 = lat%length

call convert_grid_2d_data_to_lmdz(dim1,dim2,botm_negative_lon_index,PS_local%vals)

PS_local%vals(1,:) = PS_local%vals(dim1,:) !periodicity

 ! special code:  check and error out if the PS field has gone negative
 if (minval(PS_local%vals(:,:)) < 0._r8) then
      write(msgstring, *)'PS has negative values; should not happen'
      call error_handler(E_ERR, 'write_lmdz_init', msgstring, source, revision, revdate)
 endif


 call nc_check(nf90_inq_varid(ncfileid,'ps' , ncfldid),           &
                 'write_lmdz_init','inq_varid '//trim(cflds(ifld)))

 call nc_check(nf90_put_var(ncfileid, ncfldid, PS_local%vals(1:dim1,1:dim2),  &
                           start=(/1, 1, 1/), count = (/dim1, dim2, 1/)), &
                 'write_lmdz_init','put_var '//trim(cflds(ifld)))

!--------T----
   ifld = 2
   dim1 = lon%length
   dim2 = lat%length
   dim3 = sigs%length
 call convert_grid_3d_data_to_lmdz(dim1,dim2,dim3,botm_negative_lon_index,T_local%vals)

  T_local%vals(1,:,:) = T_local%vals(dim1,:,:) !periodicity

 if (minval(T_local%vals(:,:,:)) < 0._r8) then
     write(msgstring, *)'T has negative values; should not happen'
     call error_handler(E_ERR, 'write_lmdz_init', msgstring, source,revision, revdate)
 endif

   call nc_check(nf90_inq_varid(ncfileid,'teta', ncfldid),   &
                 'write_lmdz_init', 'inq_varid '//trim(cflds(ifld)))
   call nc_check(nf90_put_var(ncfileid, ncfldid                          &
              ,T_local%vals(1:dim1, 1:dim2, 1:dim3) ,start=(/1,1,1,1/) , &
    count=(/dim1, dim2, dim3, 1/) ), 'write_lmdz_init', 'put_var '//trim(cflds(ifld)))

!--------U----
   ifld = 3
   dim1 = slon%length
   dim2 = lat%length
   dim3 = sigs%length
   
   call convert_grid_3d_data_to_lmdz(dim1,dim2,dim3,botm_negative_slon_index,U_local%vals)

   U_local%vals(1,:,:) = U_local%vals(dim1,:,:) !periodicity

   call nc_check(nf90_inq_varid(ncfileid, 'ucov', ncfldid),   &
                 'write_lmdz_init', 'inq_varid '//trim(cflds(ifld)))
   call nc_check(nf90_put_var(ncfileid, ncfldid                          &
              ,U_local%vals(1:dim1, 1:dim2, 1:dim3) ,start=(/1,1,1,1/) , &
    count=(/dim1, dim2, dim3, 1/) ), 'write_lmdz_init', 'put_var '//trim(cflds(ifld)))

!--------V----
   ifld = 4
   dim1 = lon%length
   dim2 = slat%length
   dim3 = sigs%length

   call convert_grid_3d_data_to_lmdz(dim1,dim2,dim3,botm_negative_lon_index,V_local%vals)

   V_local%vals(1,:,:) = V_local%vals(dim1,:,:) !periodicity

   call nc_check(nf90_inq_varid(ncfileid, 'vcov', ncfldid),   &
                 'write_lmdz_init', 'inq_varid '//trim(cflds(ifld)))
   call nc_check(nf90_put_var(ncfileid, ncfldid                          &
              ,V_local%vals(1:dim1, 1:dim2, 1:dim3) ,start=(/1,1,1,1/) , &
    count=(/dim1, dim2, dim3, 1/) ), 'write_lmdz_init', 'put_var '//trim(cflds(ifld)))
!--------Q----
   ifld = 5
   dim1 = lon%length
   dim2 = lat%length
   dim3 = sigs%length

   call convert_grid_3d_data_to_lmdz(lon%length,lat%length,sigs%length,botm_negative_lon_index,Q_local%vals)

    Q_local%vals(1,:,:) = Q_local%vals(dim1,:,:) !periodicity

   call nc_check(nf90_inq_varid(ncfileid,'H2Ov', ncfldid),   &
                 'write_lmdz_init', 'inq_varid '//trim(cflds(ifld)))
   call nc_check(nf90_put_var(ncfileid, ncfldid                          &
              ,Q_local%vals(1:dim1, 1:dim2, 1:dim3) ,start=(/1,1,1,1/) , &
    count=(/dim1, dim2, dim3, 1/) ), 'write_lmdz_init', 'put_var '//trim(cflds(ifld)))

!--------CLDLIQ----
   ifld = 6
   dim1 = lon%length
   dim2 = lat%length
   dim3 = sigs%length
   
   call convert_grid_3d_data_to_lmdz(lon%length,lat%length,sigs%length,botm_negative_lon_index,CLDLIQ_local%vals)

   CLDLIQ_local%vals(1,:,:) = CLDLIQ_local%vals(dim1,:,:) !periodicity

   call nc_check(nf90_inq_varid(ncfileid,'H2Ol' , ncfldid),        &
                 'write_lmdz_init', 'inq_varid '//trim(cflds(ifld)))
   call nc_check(nf90_put_var(ncfileid, ncfldid                              &
              ,CLDLIQ_local%vals(1:dim1, 1:dim2,1:dim3) ,start=(/1,1,1,1/) , &
    count=(/dim1, dim2, dim3, 1/) ), 'write_lmdz_init', 'put_var'//trim(cflds(ifld)))


!---------
  if (present( model_time)) then
    call get_date(model_time, iyear, imonth, iday, ihour, imin, isec)
!  
    call nc_check(nf90_inq_varid(ncfileid,'temps' , ncfldid),        &
                 'write_lmdz_init', 'inq_varid '//trim('temps'))
!
      write(unites,200)iyear,imonth,iday
200   format('days since ',i4,'-',i2.2,'-',i2.2,' 00:00:00')
!
      call nc_check(nf90_put_att(ncFileID, ncfldid, "units",unites), &
     'write_lmdz_init', 'nf90_put_att '//trim('temps'))
!
  end if 
!---------------
  call nc_check(nf90_close(ncfileid), &
              'static_init_model', 'closing '//trim(model_config_file))
 
end subroutine write_lmdz_init



   function nc_write_model_atts( ncFileID ) result (ierr)
!=======================================================================
! function nc_write_model_atts( ncFileID ) result (ierr)
!
! Writes the model-specific attributes to a netCDF file.
! TJH Fri Aug 29 MDT 2003
!

integer, intent(in)  :: ncFileID      ! netCDF file identifier
integer              :: ierr          ! return value of function

!-----------------------------------------------------------------------------------------

integer :: nDimensions, nVariables, nAttributes, unlimitedDimID
! merge
integer :: MemberDimID, StateVarDimID, TimeDimID,ScalarDimID
integer :: xVarID,StateVarID, StateVarVarID
integer :: P_id(num_dims)
integer :: i, ifld, dim_id, g_id
integer :: grid_id(10)   ! manual choice for 10.! num of coord var 10 is sufficient
character(len=8)      :: crdate      ! needed by F90 DATE_AND_TIME intrinsic
character(len=10)     :: crtime      ! needed by F90 DATE_AND_TIME intrinsic
character(len=5)      :: crzone      ! needed by F90 DATE_AND_TIME intrinsic
integer, dimension(8) :: values      ! needed by F90 DATE_AND_TIME intrinsic
character(len=NF90_MAX_NAME) :: str1


ierr = 0     ! assume normal termination

!-------------------------------------------------------------------------------
! Make sure ncFileID refers to an open netCDF file, and then put into define mode.
! nf90_Inquire  returns all but the ncFileID; these were defined in the calling
! routine. More dimensions, variables and attributes will be added in this routine.
!-------------------------------------------------------------------------------

write(msgstring,*) 'ncFileID', ncFileID

call nc_check(nf90_Inquire(ncFileID, nDimensions, nVariables, nAttributes, unlimitedDimID), &
              'nc_write_model_atts', 'Inquire '//trim(msgstring))
call nc_check(nf90_Redef(ncFileID), 'nc_write_model_atts', 'Redef'//trim(msgstring))

!-------------------------------------------------------------------------------
! We need the dimension ID for the number of copies
!--------------------------------------------------------------------------

call nc_check(nf90_inq_dimid(ncid=ncFileID, name="copy", dimid=MemberDimID), &
              'nc_write_model_atts', 'inq_dimid copy')
call nc_check(nf90_inq_dimid(ncid=ncFileID, name="time", dimid=  TimeDimID), &
              'nc_write_model_atts', 'inq_dimid time')

if ( TimeDimID /= unlimitedDimId ) then
  write(msgstring,*)'Time dimension ID ',TimeDimID,'must match Unlimited Dimension ID ',unlimitedDimId
  call error_handler(E_ERR,'nc_write_model_atts', msgstring, source, revision, revdate)
end if

!-------------------------------------------------------------------------------
! Define the model size, state variable dimension ... whatever ...
!-------------------------------------------------------------------------------
call nc_check(nf90_def_dim(ncid=ncFileID, name="StateVariable",  &
                        len=model_size, dimid = StateVarDimID),  &
              'nc_write_model_atts', 'def_dim StateVariable')
!-------------------------------------------------------------------------------
! Write Global Attributes
!-------------------------------------------------------------------------------

call DATE_AND_TIME(crdate,crtime,crzone,values)
! write(str1,'(''YYYY MM DD HH MM SS = '',i4,5(1x,i2.2))') &
write(str1,'("YYYY MM DD HH MM SS = ",i4,5(1x,i2.2))') &
                  values(1), values(2), values(3), values(5), values(6), values(7)

call nc_check(nf90_put_att(ncFileID, NF90_GLOBAL, "creation_date",str1), &
              'nc_write_model_atts', 'put_att creation_date'//trim(str1))
call nc_check(nf90_put_att(ncFileID, NF90_GLOBAL, "model_revision",revision), &
              'nc_write_model_atts', 'put_att model_revision'//trim(revision))
call nc_check(nf90_put_att(ncFileID, NF90_GLOBAL, "model_revdate",revdate),&
              'nc_write_model_atts', 'put_att model_revdate'//trim(revdate))
call nc_check(nf90_put_att(ncFileID, NF90_GLOBAL, "model","LMDZ"), &
              'nc_write_model_atts','put_att model LMDZ')

!-------------------------------------------------------------------------------
! Define the new dimensions IDs
!-------------------------------------------------------------------------------

! All the dim_ids on start.nc were read in, and will be written out here,
! and some will be used to define the coordinate variables below.
! They have different dimids for this file than they had for start.nc
! P_id serves as a map between the 2 sets.

if (print_details .and. do_out) write(*,*) ' dimens,       name,  size, lmdz dim_id, P[oste]rior id'

call nc_check(nf90_def_dim (ncid=ncFileID, name='lon     ', len=lon%length,           &
                    dimid=P_id(1)), 'nc_write_model_atts','def_dim'//trim('lon'))

call nc_check(nf90_def_dim (ncid=ncFileID, name='lat     ', len=lat%length,           &
                    dimid=P_id(2)), 'nc_write_model_atts','def_dim'//trim('lat'))

call nc_check(nf90_def_dim (ncid=ncFileID, name='slon    ', len=slon%length,           &
                    dimid=P_id(3)), 'nc_write_model_atts','def_dim'//trim('slon'))

call nc_check(nf90_def_dim (ncid=ncFileID, name='slat    ', len=slat%length,           &
                    dimid=P_id(4)), 'nc_write_model_atts','def_dim'//trim('slat'))

call nc_check(nf90_def_dim (ncid=ncFileID, name='sigs    ', len=sigs%length,           &
                    dimid=P_id(5)), 'nc_write_model_atts','def_dim'//trim('sigs'))

call nc_check(nf90_def_dim(ncid=ncFileID, name="scalar",   len = 1,   dimid = ScalarDimID) &
             ,'nc_write_model_atts', 'def_dim scalar')

!-------------------------------------------------------------------------------
! Create the (empty) Coordinate Variables and their attributes
!-------------------------------------------------------------------------------
! grid longitudes, latitudes, levels, and other coordinates.
! grid_id() is filled here; it's the dimid of the desired coordinate *on this
! P_Diag.nc file*.  



 call write_lmdz_coord_def(ncFileID,'lon     ',lon , P_id(1) , grid_id(1))

 call write_lmdz_coord_def(ncFileID,'lat     ',lat , P_id(2) , grid_id(2))
 call write_lmdz_coord_def(ncFileID,'slon    ',slon ,P_id(3) , grid_id(3))
 call write_lmdz_coord_def(ncFileID,'slat    ',slat ,P_id(4) , grid_id(4))
 call write_lmdz_coord_def(ncFileID,'sigs    ',sigs ,P_id(5) , grid_id(5))
! apm & bpm have dim_id of sigs i.e.P_id(5)
 call write_lmdz_coord_def(ncFileID,'apm     ',apm  ,P_id(5) , grid_id(6))
 call write_lmdz_coord_def(ncFileID,'bpm     ',bpm  ,P_id(5) , grid_id(7))

if ( output_state_vector ) then

   !----------------------------------------------------------------------------
   ! Create attributes for the state vector
   !----------------------------------------------------------------------------

   ! Define the state vector coordinate variable
   call nc_check(nf90_def_var(ncid=ncFileID,name="StateVariable", xtype=nf90_int,           &
              dimids=StateVarDimID, varid=StateVarVarID),&
                 'nc_write_model_atts','def_var  state vector')
   call nc_check(nf90_put_att(ncFileID, StateVarVarID, "long_name", "StateVariable ID"),   &
                 'nc_write_model_atts','put_att long_name state vector ')
   call nc_check(nf90_put_att(ncFileID, StateVarVarID, "units","indexical"),           &
                 'nc_write_model_atts','put_att units state vector ' )
   call nc_check(nf90_put_att(ncFileID, StateVarVarID, "valid_range", (/ 1,model_size /)), &
                 'nc_write_model_atts','put_att valid range state vector ')
   ! Define the actual state vector
   call nc_check(nf90_def_var(ncid=ncFileID, name="state", xtype=nf90_real,&
              dimids = (/ StateVarDimID, MemberDimID, unlimitedDimID /),varid=StateVarID), &
                 'nc_write_model_atts','def_var state vector')
   call nc_check(nf90_put_att(ncFileID, StateVarID, "long_name", "model state or fcopy"),   &
                 'nc_write_model_atts','put_att long_name model state or fcopy')
   call nc_check(nf90_put_att(ncFileID, StateVarId, "vector_to_prog_var","LMDZ"),&
                 'nc_write_model_atts','put_att vector_to_prog_var LMDZ ')

   ! Leave define mode so we can fill 
   call nc_check(nf90_enddef(ncfileID), 'nc_write_model_atts','enddef ')

   ! Fill the state variable coordinate variable
   call nc_check(nf90_put_var(ncFileID, StateVarVarID, (/ (i,i=1,model_size) /)),         &
                 'nc_write_model_atts','put_var StateVar ')

else
!--------PS-------
     ifld =  1
      call nc_check(nf90_def_var(ncid=ncFileID, name=trim(cflds(ifld)), xtype=nf90_real, &
                 dimids = (/ P_id(1), P_id(2), &
                             MemberDimID, unlimitedDimID /), varid  = xVarID),&
                 'nc_write_model_atts','def_var 2d '//trim(cflds(ifld)))

      call nc_check(nf90_put_att(ncFileID, xVarID, "long_name", 'Surface Pressure'), &
                 'nc_write_model_atts','put_att long_name ')
      call nc_check(nf90_put_att(ncFileID, xVarID, "units", 'Pa'),&
                 'nc_write_model_atts','put_att units ')
!-------T----------
      ifld = 2
      call nc_check(nf90_def_var &
           (ncid=ncFileID, name=trim(cflds(ifld)), xtype=nf90_real, &
            dimids = (/ P_id(1), P_id(2), P_id(5),  &
                        MemberDimID, unlimitedDimID /), &
            varid  = xVarID),&
                 'nc_write_model_atts','def_var 3d'//trim(cflds(ifld)))
      call nc_check(nf90_put_att(ncFileID, xVarID, "long_name", 'Air Temperature'),      &
                 'nc_write_model_atts','put_att long_name' )
      call nc_check(nf90_put_att(ncFileID, xVarID, "units", 'deg K'), &
                 'nc_write_model_atts','put_att units ')

!------U--------
      ifld = 3
      call nc_check(nf90_def_var & 
           (ncid=ncFileID, name=trim(cflds(ifld)), xtype=nf90_real, &
            dimids = (/ P_id(3), P_id(2), P_id(5),  &
                        MemberDimID, unlimitedDimID /), &
            varid  = xVarID),&
                 'nc_write_model_atts','def_var 3d'//trim(cflds(ifld)))
      call nc_check(nf90_put_att(ncFileID, xVarID, "long_name", 'Zonal WinD'),      &
                 'nc_write_model_atts','put_att long_name ')       
      call nc_check(nf90_put_att(ncFileID, xVarID, "units", 'm/s'), &
                 'nc_write_model_atts','put_att units ')

!-------V-------    

ifld = 4
      call nc_check(nf90_def_var &
           (ncid=ncFileID, name=trim(cflds(ifld)), xtype=nf90_real, &
            dimids = (/ P_id(1), P_id(4), P_id(5),  &
                        MemberDimID, unlimitedDimID /), &
            varid  = xVarID),&
                 'nc_write_model_atts','def_var 3d'//trim(cflds(ifld)))
      call nc_check(nf90_put_att(ncFileID, xVarID, "long_name", 'Meridional WinD'), &
                 'nc_write_model_atts','put_att long_name ')
      call nc_check(nf90_put_att(ncFileID, xVarID, "units", 'm/s'), &
                 'nc_write_model_atts','put_att units ')

!-------Q--------
ifld = 5
      call nc_check(nf90_def_var &
           (ncid=ncFileID, name=trim(cflds(ifld)), xtype=nf90_real, &
            dimids = (/ P_id(1), P_id(2), P_id(5),  &
                        MemberDimID, unlimitedDimID /), &
            varid  = xVarID),&
                 'nc_write_model_atts','def_var 3d'//trim(cflds(ifld)))
      call nc_check(nf90_put_att(ncFileID, xVarID, "long_name", 'Specific Humidity'), &
                 'nc_write_model_atts','put_att long_name ')
      call nc_check(nf90_put_att(ncFileID, xVarID, "units", ' '), &
                 'nc_write_model_atts','put_att units ')


!-----CLDLIQ--
ifld = 6
      call nc_check(nf90_def_var &
           (ncid=ncFileID, name=trim(cflds(ifld)), xtype=nf90_real, &
            dimids = (/ P_id(1), P_id(2), P_id(5),  &
                        MemberDimID, unlimitedDimID /), &
            varid  = xVarID),&
                 'nc_write_model_atts','def_var 3d'//trim(cflds(ifld)))
      call nc_check(nf90_put_att(ncFileID, xVarID, "long_name", 'CLDLIQ'), &
                 'nc_write_model_atts','put_att long_name ')
      call nc_check(nf90_put_att(ncFileID, xVarID, "units", ' '), &
                 'nc_write_model_atts','put_att units ')

  ! Leave define mode so we can fill variables
   call nc_check(nf90_enddef(ncfileID),                         &
                 'nc_write_model_atts','enddef ')


!--------
!-------------------------------------------------------------------------------
! Fill the coordinate variables
! Each 'vals' vector has been dimensioned to the right size for its coordinate.  
! The default values of 'start' and 'count'  write out the whole thing.
!-------------------------------------------------------------------------------
if (print_details .and. do_out) write(*,*) 'nc_write_model_atts; filling coords'



 call nc_check(nf90_put_var(ncFileID, grid_id(1),  lon%vals) &
                 ,'nc_write_model_atts', 'put_var lon')

 call nc_check(nf90_put_var(ncFileID, grid_id(2),  lat%vals) &
                 ,'nc_write_model_atts', 'put_var lat')

 call nc_check(nf90_put_var(ncFileID, grid_id(3),  slon%vals) &
                 ,'nc_write_model_atts', 'put_var slon')

 call nc_check(nf90_put_var(ncFileID, grid_id(4),  slat%vals) &
                 ,'nc_write_model_atts', 'put_var slat')

 call nc_check(nf90_put_var(ncFileID, grid_id(5),  sigs%vals) &
                 ,'nc_write_model_atts', 'put_var sigs')

 call nc_check(nf90_put_var(ncFileID, grid_id(6),  apm%vals) &
                 ,'nc_write_model_atts', 'put_var apm')

 call nc_check(nf90_put_var(ncFileID, grid_id(7),  bpm%vals) &
                 ,'nc_write_model_atts', 'put_var bpm')
end if


!-------------------------------------------------------------------------------
! Flush the buffer and leave netCDF file open
!-------------------------------------------------------------------------------
call nc_check(nf90_sync(ncFileID),'nc_write_model_atts', 'sync ')

write (*,*)'nc_write_model_atts: netCDF file ',ncFileID,' is synched ...'


end function  nc_write_model_atts

 


  function nc_write_model_vars( ncFileID, statevec, copyindex, timeindex ) result (ierr)
!=======================================================================
! function nc_write_model_vars( ncFileID, statevec, copyindex, timeindex )
! result (ierr)
!
! Writes the model-specific variables to a netCDF file
! TJH 25 June 2003
!

integer,                intent(in) :: ncFileID      ! netCDF file identifier
real(r8), dimension(:), intent(in) :: statevec
integer,                intent(in) :: copyindex
integer,                intent(in) :: timeindex
integer                            :: ierr          ! return value of function

!-----------------------------------------------------------------------------------------
!type(model_type) :: Var

integer :: nDimensions, nVariables, nAttributes, unlimitedDimID
integer :: StateVarID, ncfldid
integer :: ifld, ii,dim1,dim2,dim3

character (len=8) :: cfield

type(data_2d_type) :: PS_local
type(data_3d_type) :: T_local,U_local,V_local,Q_local,CLDLIQ_local

ierr = 0     ! assume normal termination

!-------------------------------------------------------------------------------
! make sure ncFileID refers to an open netCDF file, 
! then get all the Variable ID's we need.
!-------------------------------------------------------------------------------

call nc_check(nf90_Inquire(ncFileID, nDimensions, nVariables, nAttributes,unlimitedDimID), &
              'nc_write_model_vars','Inquire ')

if ( output_state_vector ) then

   call nc_check(nf90_inq_varid(ncFileID, "state", StateVarID),'nc_write_model_vars ','inq_varid state' )
   call nc_check(nf90_put_var(ncFileID, StateVarID, statevec,  &
                start=(/ 1, copyindex, timeindex /)),'nc_write_model_vars','put_var state')              

else

   !----------------------------------------------------------------------------
   ! Fill the variables
   !----------------------------------------------------------------------------
                                                                                       
call init_model_instance(PS_local,T_local,U_local,V_local,Q_local,CLDLIQ_local)     ! Explicity released at end of routine. 

call vector_to_prog_var(statevec,PS_local,T_local,U_local,V_local,Q_local,CLDLIQ_local)
  
!------------PS--------
  ifld =  1
  dim1=lon%length
  dim2=lat%length
  cfield = trim(cflds(ifld))

   call nc_check(nf90_inq_varid(ncFileID, cfield, ncfldid), &
                    'nc_write_model_vars ','inq_varid 2d '//cfield)
   call nc_check(nf90_put_var(ncFileID, ncfldid, &
                   PS_local%vals(1:dim1,1:dim2), &
                    start=   (/ 1              ,1              , copyindex, timeindex /),         &
                    count=   (/dim1, dim2, 1  , 1/) ),                                       &
                    'nc_write_model_vars ','put_var 2d '//cfield)
!-----------T---------
  ifld = 2
  dim1 = lon%length
  dim2 = lat%length
  dim3 = sigs%length
  cfield = trim(cflds(ifld))

      call nc_check(nf90_inq_varid(ncFileID, cfield, ncfldid), &
                    'nc_write_model_vars ','inq_varid 3d '//cfield)
      call nc_check(nf90_put_var(ncFileID, ncfldid,    &
                 T_local%vals(1:dim1, 1:dim2,  1:dim3) &
                 ,start=   (/1              ,1              ,1  ,copyindex,timeindex/)&
                 ,count=   (/  dim1,  dim2,  dim3, 1 , 1/) ),     &
                    'nc_write_model_vars ','put_var 3d '//cfield)

!-----------U---------
  ifld = 3
  dim1 = slon%length
  dim2 = lat%length
  dim3 = sigs%length
  cfield = trim(cflds(ifld))

      call nc_check(nf90_inq_varid(ncFileID, cfield, ncfldid), &
                    'nc_write_model_vars ','inq_varid 3d '//cfield)
      call nc_check(nf90_put_var(ncFileID, ncfldid,    &
                 U_local%vals(1:dim1, 1:dim2,  1:dim3) &
                 ,start=   (/1              ,1              ,1  ,copyindex,timeindex/)&
                 ,count=   (/  dim1,  dim2,  dim3, 1 , 1/) ),     &
                    'nc_write_model_vars ','put_var 3d '//cfield)
!-----------V---------
  ifld = 4
  dim1 = lon%length
  dim2 = slat%length
  dim3 = sigs%length
  cfield = trim(cflds(ifld))

      call nc_check(nf90_inq_varid(ncFileID, cfield, ncfldid), &
                    'nc_write_model_vars ','inq_varid 3d '//cfield)
      call nc_check(nf90_put_var(ncFileID, ncfldid,    &
                 V_local%vals(1:dim1, 1:dim2,  1:dim3) &
                 ,start=   (/1              ,1              ,1  ,copyindex,timeindex/)&
                 ,count=   (/  dim1,  dim2,  dim3, 1 , 1/) ),     &
                    'nc_write_model_vars ','put_var 3d '//cfield)
!-----------Q---------
  ifld = 5
  dim1 = lon%length
  dim2 = lat%length
  dim3 = sigs%length
  cfield = trim(cflds(ifld))

      call nc_check(nf90_inq_varid(ncFileID, cfield, ncfldid), &
                    'nc_write_model_vars ','inq_varid 3d '//cfield)
      call nc_check(nf90_put_var(ncFileID, ncfldid,    &
                 Q_local%vals(1:dim1, 1:dim2,  1:dim3) &
                 ,start=   (/1              ,1              ,1  ,copyindex,timeindex/)&
                 ,count=   (/  dim1,  dim2,  dim3, 1 , 1/) ),     &
                    'nc_write_model_vars ','put_var 3d '//cfield)

!-----------CLDLIQ---------
  ifld = 6
  dim1 = lon%length
  dim2 = lat%length
  dim3 = sigs%length
  cfield = trim(cflds(ifld))

      call nc_check(nf90_inq_varid(ncFileID, cfield, ncfldid), &
                    'nc_write_model_vars ','inq_varid 3d '//cfield)
      call nc_check(nf90_put_var(ncFileID, ncfldid,    &
                 CLDLIQ_local%vals(1:dim1, 1:dim2,  1:dim3) &
                 ,start=   (/1              ,1              ,1  ,copyindex,timeindex/)&
                 ,count=   (/  dim1,  dim2,  dim3, 1 , 1/) ),     &
                    'nc_write_model_vars ','put_var 3d '//cfield)

end if

!-------------------------------------------------------------------------------
! Flush the buffer and leave netCDF file open
!-------------------------------------------------------------------------------

!write (*,*)'Finished filling variables ...'
call nc_check(nf90_sync(ncFileID),'nc_write_model_vars ','sync ')
!write (*,*)'netCDF file is synched ...'

! temporary output
!print*,'num_calced, num_searched = ',num_calced, num_searched

call end_model_instance(PS_local,T_local,U_local,V_local,Q_local,CLDLIQ_local)   ! should avoid any memory leaking


end function nc_write_model_vars







  subroutine order_state_fields()
!=======================================================================
TYPE_PS     = 1
TYPE_T      = 2
TYPE_U      = 3
TYPE_V      = 4
TYPE_Q      = 5
TYPE_CLDLIQ = 6
end subroutine order_state_fields


  subroutine change_lon_lat_lev_to_dart()
!=======================================================================

integer :: i, j
real(r8), allocatable :: pres_tmp(:),ap_tmp(:),bp_tmp(:),tmp_lon(:), &
                         tmp_slon(:),tmp_lat(:),tmp_slat(:)


!--------------------lon ..convert from (-180,180) to (0,360)formate ----------
allocate(tmp_lon(lon%length))
tmp_lon(:) = MISSING_R8

botm_positive_lon_index = -1  ! initilise  
 do i=1, lon%length 
   if( lon%vals(i) >= 0.0_r8) then 
     botm_positive_lon_index = i
     go to 101
   endif
 end do

101 j = 1

   do i = botm_positive_lon_index, lon%length
     !print*,lon%vals(i)
     tmp_lon(j) = lon%vals(i) 
     j = j + 1
   end do

    botm_negative_lon_index = j  ! save for reverse conversion 
   
   do i = 1 , botm_positive_lon_index - 1
     tmp_lon(j) = lon%vals(i) + 360.0_r8
    !print*,tmp_lon(j)
     j = j + 1
   end do
   lon%vals =tmp_lon
   deallocate(tmp_lon)

!--------------------slon ..convert from (-180,180) to (0,360)formate ----------
allocate(tmp_slon(slon%length))
tmp_slon(:) = MISSING_R8

botm_positive_slon_index = -1  ! initilise  
 do i=1, slon%length
   if( slon%vals(i) >= 0.0_r8) then  
     botm_positive_slon_index = i
     go to 102
   endif
 end do

102 j=1
   do i =botm_positive_slon_index, slon%length
     !print*,slon%vals(i)
     tmp_slon(j) = slon%vals(i)
     j = j + 1
   end do

   botm_negative_slon_index = j  ! save for reverse conversion 

   do i = 1 , botm_positive_slon_index - 1
     tmp_slon(j) = slon%vals(i) + 360.0_r8
     !print*,tmp_slon(j)
     j = j + 1
   end do
   slon%vals =tmp_slon
   deallocate(tmp_slon)

!---- change lat = (90,-90) to lat = (-90,90)
allocate(tmp_lat(lat%length))
 j = 1
 do i = lat%length, 1, -1
    tmp_lat(j)=lat%vals(i)
    j = j + 1
 end do
 lat%vals =tmp_lat
 deallocate(tmp_lat)

!print*,lat%vals
!---- change slat = (90,-90) to slat = (-90,90)
allocate(tmp_slat(slat%length))
 j = 1
 do i = slat%length, 1, -1
    tmp_slat(j)=slat%vals(i)
    j = j + 1
 end do

slat%vals =tmp_slat
deallocate(tmp_slat)
!print*,slat%vals
!__________________________________________________________________
allocate(pres_tmp(presnivs%length))
allocate(ap_tmp(ap%length))
allocate(bp_tmp(bp%length))

!change vertical coordinate index from top to bottom 
j=1
do i=presnivs%length,1,-1
    pres_tmp(j)=presnivs%vals(i)
    j=j+1
end do
presnivs%vals=pres_tmp

j=1
do i=ap%length,1,-1
   ap_tmp(j)=ap%vals(i)
   bp_tmp(j)=bp%vals(i)
   j = j + 1
end do

ap%vals = ap_tmp
bp%vals = bp_tmp

deallocate(pres_tmp, ap_tmp, bp_tmp)

end subroutine change_lon_lat_lev_to_dart



subroutine convert_grid_3d_data_to_dart(lon_len,lat_len,lev_len,indx,var)
!****************************************************************
 implicit none
 integer, intent(in)    :: indx
 integer, intent(in)    :: lon_len,lat_len,lev_len
 real(r8),intent(inout) :: var(:,:,:)
 real(r8),allocatable   :: var_tmp(:,:,:)
 integer                :: i,j

 allocate(var_tmp(lon_len,lat_len,lev_len))
 var_tmp(:,:,:) = MISSING_R8
!-------change lon=(-180,180) to lon=(0,360)
 j = 1
 do i =indx, lon_len
   var_tmp(j, :, :) = var(i, :, :)
    j = j + 1
 end do
 do i = 1 , indx - 1
    var_tmp(j, :, :) = var(i, :, :)
   j = j + 1
 end do
 var=var_tmp
 var_tmp(:,:,:) = MISSING_R8
!---- change lat = (90,-90) to lat = (-90,90)

 j = 1
 do i = lat_len, 1, -1
    var_tmp(:,j,:)=var(:,i,:)
    j = j + 1
 end do

 var=var_tmp
 var_tmp(:,:,:) = MISSING_R8
!-------change model level form (bottom,top) to (top,bottom) index
   j = 1
 do i = lev_len, 1, -1
    var_tmp(:,:,j)=var(:,:,i)
    j = j + 1
 end do

  var=var_tmp

 deallocate(var_tmp)
  end subroutine

!****************************************************************


subroutine convert_grid_2d_data_to_dart(lon_len,lat_len,indx,var)
!****************************************************************
 implicit none
 integer, intent(in)    :: indx
 integer, intent(in)    :: lon_len,lat_len
 real(r8),intent(inout) :: var(:,:)
 real(r8),allocatable   :: var_tmp(:,:)
 integer                :: i,j

 allocate(var_tmp(lon_len,lat_len))
 var_tmp(:,:) = 0._r8

 j = 1
 do i =indx, lon_len
   var_tmp(j, :) = var(i, :)
    j = j + 1
 end do
 do i = 1 , indx - 1
    var_tmp(j, :) = var(i, :)
    j = j + 1
 end do

 var=var_tmp
!---- change lat = (90,-90) to lat = (-90,90)
 j = 1
 do i = lat_len, 1, -1
    var_tmp(:,j)=var(:,i)
    j = j + 1
 end do

  var=var_tmp

 deallocate(var_tmp)
 end subroutine

 subroutine convert_grid_3d_data_to_lmdz(lon_len,lat_len,lev_len,indx,var)
!****************************************************************
 implicit none
 integer, intent(in)    :: indx
 integer, intent(in)    :: lon_len,lat_len,lev_len
 real(r8),intent(inout) :: var(:,:,:)
 real(r8),allocatable   :: var_tmp(:,:,:)
 integer                :: i,j

 allocate(var_tmp(lon_len,lat_len,lev_len))
 var_tmp(:,:,:) = MISSING_R8
!-------change lon=(-180,180) to lon=(0,360)
 j = 1
 do i = indx, lon_len
   var_tmp(j, :, :) = var(i, :, :)
    j = j + 1
 end do
 do i = 1 ,  indx  - 1
    var_tmp(j, :, :) = var(i, :, :)
   j = j + 1
 end do
 var=var_tmp
 var_tmp(:,:,:) = MISSING_R8
!---- change lat = (90,-90) to lat = (-90,90)

 j = 1
 do i = lat_len, 1, -1
    var_tmp(:,j,:)=var(:,i,:)
    j = j + 1
 end do

 var=var_tmp
 var_tmp(:,:,:) = MISSING_R8
!-------change model level form (bottom,top) to (top,bottom) index

   j = 1
 do i = lev_len, 1, -1
    var_tmp(:,:,j)=var(:,:,i)
    j = j + 1
 end do

  var=var_tmp
  deallocate(var_tmp)
end subroutine


 subroutine convert_grid_2d_data_to_lmdz(lon_len,lat_len,indx,var)
!****************************************************************
 implicit none
 integer, intent(in)    :: lon_len,lat_len,indx
 real(r8),intent(inout) :: var(:,:)
 real(r8),allocatable   :: var_tmp(:,:)
 integer                :: i,j

 allocate(var_tmp(lon_len,lat_len))
 var_tmp(:,:) = 0._r8

!-------change data  lon=(0,360) to lon=(-180,180)

 j = 1
do i = indx, lon_len
   var_tmp(j, :) = var(i, :)
    j = j + 1
 end do
 do i = 1 , indx - 1
    var_tmp(j, :) = var(i, :)
    j = j + 1
 end do

 var=var_tmp
!---- change lat = (-90,90) to lat = (90,-90)
 j = 1
 do i = lat_len, 1, -1
    var_tmp(:,j)=var(:,i)
    j = j + 1
 end do

  var=var_tmp

 deallocate(var_tmp)
 end subroutine


 
  subroutine scopy(n,sx,incx,sy,incy)
!======================================================
  IMPLICIT NONE

  integer  n,incx,incy,ix,iy,i
  real sx((n-1)*incx+1),sy((n-1)*incy+1)

   iy=1
   ix=1
     do  i=1,n
         sy(iy)=sx(ix)
         ix=ix+incx
         iy=iy+incy
     end do

end subroutine scopy


  subroutine rad_to_degree(coord)
!======================================================
!convesion form radian to degree
  type(grid_1D_type) :: coord 
  real(r8)           :: pi
  integer            :: i

  pi = 4._r8 * atan (1._r8)

  do i= 1, size(coord%vals)
    coord%vals(i) = coord%vals(i) * (180._r8 / pi) 
  end do
 
end subroutine rad_to_degree

  subroutine degree_to_rad(coord)
!======================================================
!convesion form radian to degree
  type(grid_1D_type) :: coord 
  real(r8)           :: pi
  integer            :: i

  pi = 4._r8 * atan (1._r8)

  do i= 1, size(coord%vals)
    coord%vals(i) = coord%vals(i) * ( pi / 180._r8 ) 
  end do
 
end subroutine degree_to_rad



!======================================================
end module

! <next few lines under version control, do not edit>
! $URL$
! $Id$
! $Revision$
! $Date$
