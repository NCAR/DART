! DART software - Copyright UCAR. This open source software is provided
! by UCAR, "as is", without charge, subject to all terms of use at
! http://www.image.ucar.edu/DAReS/DART/DART_download

module model_mod

use             types_mod, only : r4, r8, i8, MISSING_R8, obstypelength

use      time_manager_mod, only : time_type, set_time, set_date, set_calendar_type, &
                                  get_time, get_date, print_time, print_date, &
                                  operator(*),  operator(+), operator(-), &
                                  operator(>),  operator(<), operator(/), &
                                  operator(/=), operator(<=)

use          location_mod, only : location_type, get_close_type, get_dist, &
                                  get_close_obs, get_close_state, &
                                  convert_vertical_obs, convert_vertical_state, &
                                  set_location, set_location_missing, &
                                  query_location, write_location, &
                                  get_location, is_vertical,  &
                                  VERTISSURFACE, VERTISHEIGHT

use         utilities_mod, only : error_handler, do_output, &
                                  E_ERR, E_MSG, file_exist, get_unit, &
                                  logfileunit, nmlfileunit, to_upper, &
                                  do_nml_file, do_nml_term, &
                                  find_namelist_in_file, check_namelist_read, &
                                  find_textfile_dims, file_to_text

use  netcdf_utilities_mod, only : nc_check, nc_synchronize_file, &
                                  nc_add_global_attribute, &
                                  nc_add_global_creation_time, &
                                  nc_begin_define_mode, &
                                  nc_end_define_mode, &
                                  nc_get_global_attribute, &
                                  nc_open_file_readonly, nc_close_file

use obs_def_utilities_mod, only : track_status

use     obs_kind_mod,      only : get_index_for_quantity, &
                                  get_name_for_quantity, &
                                  QTY_SNOWCOVER_FRAC, &
                                  QTY_GEOPOTENTIAL_HEIGHT

use  ensemble_manager_mod, only : ensemble_type

use distributed_state_mod, only : get_state

use     default_model_mod, only : adv_1step, nc_write_model_vars

use        noah_hydro_mod, only : configure_lsm, get_noah_timestepping, &
                                  num_soil_layers, lsm_namelist_filename, &
                                  soil_layer_thickness, get_lsm_domain_info, &
                                  wrf_static_data, get_lsm_domain_filename, &
                                  read_noah_global_atts, write_noah_global_atts, &
                                  num_soil_nitrogen_layers

use   state_structure_mod, only : add_domain,      get_domain_size,   &
                                  get_index_start, get_index_end,     &
                                  get_num_domains, get_num_variables, &
                                  get_num_dims,    get_dim_name,      &
                                  get_dim_length,  get_variable_name, &
                                  get_model_variable_indices,         &
                                  get_varid_from_kind,                &
                                  get_dart_vector_index,              &
                                  get_variable_size,                  &
                                  state_structure_info

use     mpi_utilities_mod, only : my_task_id

use           options_mod, only : set_missing_ok_status

use        random_seq_mod, only : random_seq_type, init_random_seq, random_gaussian

use             map_utils, only : proj_info, map_init, map_set, &
                                  latlon_to_ij, ij_to_latlon

use misc_definitions_module, only : map_latlon       => PROJ_LATLON, &
                                    map_lambert      => PROJ_LC, &
                                    map_polar_stereo => PROJ_PS, &
                                    map_mercator     => PROJ_MERC, &
                                    map_cylindrical  => PROJ_CYL, &
                                    map_cassini      => PROJ_CASSINI

use typesizes
use netcdf

implicit none
private

! required by DART code - will be called from filter and other
! DART executables.  interfaces to these routines are fixed and
! cannot be changed in any way.
public :: static_init_model, &
          init_time, &
          init_conditions, &
          pert_model_copies, &
          get_model_size, &
          read_model_time, &
          write_model_time, &
          shortest_time_between_assimilations, &
          get_state_meta_data, &
          model_interpolate, &
          get_close_obs, &
          get_close_state, &
          nc_write_model_atts, &
          end_model

! Also required, but the default routines are sufficient.

public :: nc_write_model_vars, &
          adv_1step, &
          convert_vertical_obs, &
          convert_vertical_state

character(len=*), parameter :: source   = 'noah_model_mod.f90'

!------------------------------------------------------------------
! The NSOLDX parameter comes from the NOAH source code. We need it
! because we have to read the NOAH namelist for timestep information.
!------------------------------------------------------------------

integer, parameter :: NSOLDX = 100
character(len=512) :: string1, string2, string3

! Codes for interpreting the columns of the variable_table
integer, parameter :: VT_VARNAMEINDX  = 1 ! ... variable name
integer, parameter :: VT_KINDINDX     = 2 ! ... DART kind
integer, parameter :: VT_MINVALINDX   = 3 ! ... minimum value if any
integer, parameter :: VT_MAXVALINDX   = 4 ! ... maximum value if any
integer, parameter :: VT_STATEINDX    = 5 ! ... update (state) or not
integer, parameter :: MAX_STATE_VARIABLES     = 40
integer, parameter :: NUM_STATE_TABLE_COLUMNS = 5

integer :: domain_count
integer :: idom, idom_lsm = -1

!------------------------------------------------------------------
! things which can/should be in the DART model_nml
! The variables in the noah restart file that are used to create the
! DART state vector are specified in the input.nml:model_nml namelist.
! For example:
!
! lsm_variables = 'SOIL_T',   'QTY_SOIL_TEMPERATURE',   '0.0',  'NA', 'UPDATE',
!                 'SMC',      'QTY_SOIL_MOISTURE',      '0.0',  'NA', 'UPDATE',
!                 'SOIL_W',   'QTY_SOIL_LIQUID_WATER',  '0.0',  'NA', 'UPDATE',
!                 'SKINTEMP', 'QTY_SKIN_TEMPERATURE',   '0.0',  'NA', 'UPDATE',
!                 'SNODEP',   'QTY_SNOW_THICKNESS',     '0.0',  'NA', 'UPDATE',
!                 'WEASD',    'QTY_SNOW_WATER',         '0.0',  'NA', 'UPDATE',
!                 'CANWAT',   'QTY_CANOPY_WATER',       '0.0',  'NA', 'UPDATE',
!                 'QFX',      'QTY_LATENT_HEAT_FLUX',   '0.0',  'NA', 'UPDATE',
!                 'HFX',      'QTY_SENSIBLE_HEAT_FLUX', '0.0',  'NA', 'UPDATE',
!                 'GRDFLX',   'QTY_GROUND_HEAT_FLUX'    '0.0',  'NA', 'UPDATE',
!
!------------------------------------------------------------------

character(len=256)    :: domain_shapefiles(3)         = ''
character(len=256)    :: lsm_model_choice             = ''
integer               :: assimilation_period_days     = 0
integer               :: assimilation_period_seconds  = 60
real(r8)              :: model_perturbation_amplitude = 0.002
character(len=256)    :: perturb_distribution         = 'lognormal'
integer               :: debug    = 0  ! turn up for more and more debug messages
character(len=obstypelength) :: lsm_variables(NUM_STATE_TABLE_COLUMNS,MAX_STATE_VARIABLES) = ' '

!nc -- we are adding these to the model.nml until they appear in the NetCDF files
logical :: polar      = .false.    ! wrap over the poles
logical :: periodic_x = .false.    ! wrap in longitude or x
logical :: periodic_y = .false.    ! used for single column model, wrap in y

namelist /model_nml/ domain_shapefiles, &
                     lsm_model_choice, &
                     assimilation_period_days, &
                     assimilation_period_seconds, &
                     model_perturbation_amplitude, &
                     perturb_distribution, &
                     debug, &
                     lsm_variables, &
                     polar, periodic_x, periodic_y

!------------------------------------------------------------------

type domain_locations
   private
   type(location_type), allocatable :: location(:,:)
end type domain_locations

type(domain_locations), allocatable :: domain_info(:)

! Thinking about supporting multiple (nested) domains for noah
! This follows the WRF mechanism. The nests must be numbered
! from the outermost (as domain 1) to innermost (as domain N).

type lsm_dom
   private
   type(wrf_static_data), pointer :: dom(:)
end type lsm_dom

type(lsm_dom) :: lsm

!------------------------------------------------------------------
! module storage
!------------------------------------------------------------------

integer            :: model_size       ! the state vector length
logical, save      :: module_initialized = .false.
character(len=32)  :: calendar = 'Gregorian'

real(r8), allocatable :: soil_depth(:)
real(r8), allocatable :: soil_nitrogen_depth(:)
real(r8), allocatable :: xlong(:,:), xlat(:,:)
integer :: south_north, west_east

!==================================================================
CONTAINS
!==================================================================

!------------------------------------------------------------------
!> one-time initialization of the model

subroutine static_init_model()

character(len=*), parameter :: routine = 'static_init_model'

integer  :: iunit, io, domainID
integer  :: n_lsm_fields
integer  :: i

character(len=obstypelength) :: var_names(MAX_STATE_VARIABLES)
real(r8) :: var_ranges(MAX_STATE_VARIABLES,2)
logical  :: var_update(MAX_STATE_VARIABLES)
integer  :: var_qtys(  MAX_STATE_VARIABLES)

character(len=256) :: filename

if ( module_initialized ) return ! only need to do this once.

! Since this routine calls other routines that could call this routine
! we'll say we've been initialized pretty dang early.
module_initialized = .true.

! Read the DART namelist
call find_namelist_in_file('input.nml', 'model_nml', iunit)
read(iunit, nml = model_nml, iostat = io)
call check_namelist_read(iunit, io, 'model_nml')

! Record the DART namelist values used for the run ...
if (do_nml_file()) write(nmlfileunit, nml=model_nml)
if (do_nml_term()) write(     *     , nml=model_nml)

! It is OK for this model to have MISSING_R8 in the DART state.
call set_missing_ok_status(.true.)

call set_calendar_type( calendar )

! Determine the composition of the DART state vector.

if ( .not. file_exist(domain_shapefiles(1)) ) then
   write(string1,*) 'Domain 1 shapefile "', trim(domain_shapefiles(1)),'" does not exist.'
   write(string2,*) 'There must be a domain 1.'
   call error_handler(E_ERR,routine,string1,source,text2=string2)
endif

DOMAINS: do domainID = 1,size(domain_shapefiles)

   if (len_trim(domain_shapefiles(domainID)) == 0) exit DOMAINS

   filename = domain_shapefiles(domainID)

   if ( .not. file_exist(domain_shapefiles(domainID)) ) then
      write(string1,'("Domain shapefile ",i3, A," does not exist.")')  &
                     domainID, '"'//trim(domain_shapefiles(domainID))//'"'
      call error_handler(E_ERR,routine,string1,source)
   else
      call configure_lsm(lsm_model_choice,domain_shapefiles(domainID))
      call read_noah_global_atts(domain_shapefiles(domainID))
      call verify_variables(lsm_variables, domain_shapefiles(domainID), n_lsm_fields, &
                       var_names, var_qtys, var_ranges, var_update)
      idom_lsm = add_domain(domain_shapefiles(domainID), n_lsm_fields, var_names, &
                         kind_list=var_qtys, &
                        clamp_vals=var_ranges(1:n_lsm_fields,:), &
                       update_list=var_update)
      if (debug > 99) call state_structure_info(idom_lsm)

      call check_vertical_dimension(idom_lsm)

   endif

enddo DOMAINS

! Each domain has its own array of locations.
! step 1: finish reading attributes from wrfinput.nc
! step 2: finish setup_map_projection, which contains map_init and map_set

domain_count = get_num_domains()

allocate(    lsm%dom(domain_count))
allocate(domain_info(domain_count))

call configure_domains()

!soil_layer_thickness is the thickness of each soil layer,
! soil_depth is the midpoint of each soil layer.

allocate(soil_depth(num_soil_layers))
soil_depth(1) = soil_layer_thickness(1)/2.0_r8
do i = 2,num_soil_layers
   soil_depth(i) = soil_depth(i-1) + &
                   (soil_layer_thickness(i)+soil_layer_thickness(i-1))/2.0_r8
enddo

! soil nitrogen support implemented for UT Austin research code.

if (num_soil_nitrogen_layers > 0) then
   allocate(soil_nitrogen_depth(num_soil_nitrogen_layers))
   
   soil_nitrogen_depth(1) = 0.005_r8 
   soil_nitrogen_depth(2) = soil_layer_thickness(1)/2.0_r8+soil_nitrogen_depth(1)
   
   do i = 3,num_soil_nitrogen_layers
      soil_nitrogen_depth(i) = soil_nitrogen_depth(i-1) + &
                      (soil_layer_thickness(i-1)+soil_layer_thickness(i-2))/2.0_r8
   enddo
endif

! Now that we know the composition of the DART state, determine the model size.

model_size = 0
do idom = 1,domain_count
   model_size = model_size + get_domain_size(idom)
enddo

end subroutine static_init_model


!------------------------------------------------------------------
!> Returns a time that is somehow appropriate for starting up a long
!> integration of the model.  At present, this is only used if the namelist
!> parameter start_from_restart is set to .false. in the program perfect_model_obs.

subroutine init_time(time)

type(time_type), intent(out) :: time

character(len=*), parameter :: routine = 'init_time'

if ( .not. module_initialized ) call static_init_model

! for now, just set to 0
time = set_time(0,0)

write(string1,*) 'no good way to specify initial time'
call error_handler(E_ERR,routine,string1,source)

end subroutine init_time



!------------------------------------------------------------------
!> Returns a model state vector, x, that is some sort of appropriate
!> initial condition for starting up a long integration of the model.
!> At present, this is only used if the namelist parameter
!> start_from_restart is set to .false. in the program perfect_model_obs.

subroutine init_conditions(x)

real(r8), intent(out) :: x(:)

character(len=*), parameter :: routine = 'init_conditions'

if ( .not. module_initialized ) call static_init_model

write(string1,*) 'PROBLEM: no known way to set arbitrary initial conditions.'
write(string2,*) 'start_from_restart must be .true. in this model.'
call error_handler(E_ERR, routine, string1, source, text2=string2)

x = MISSING_R8

end subroutine init_conditions


!------------------------------------------------------------------
!> Returns the size of the model as a long integer

function get_model_size()

integer(i8) :: get_model_size

if ( .not. module_initialized ) call static_init_model

get_model_size = model_size

end function get_model_size


!-----------------------------------------------------------------------
!> The LSM restart files have "time".
!> We are always using the 'most recent' which is, by defn, the last one.
!>
!> LSM restart filename is RESTART.2004010102_DOMAIN1 has
!>     Times = '2004-01-01_02:00:00' ;
!> The data is valid @ 2004-01-01_02:00:00
!>
!> This routine is identical to the wrf_hydro:model_mod:read_model_time()

function read_model_time(filename)

type(time_type) :: read_model_time
character(len=*),  intent(in)  :: filename    !! the file name

character(len=*), parameter :: routine = 'read_model_time'

integer, parameter :: STRINGLENGTH = 19
character(len=STRINGLENGTH), allocatable, dimension(:) :: datestring
character(len=STRINGLENGTH)                            :: datestring_scalar
integer               :: year, month, day, hour, minute, second
integer               :: DimID, VarID, strlen, ntimes
logical               :: isLsmFile
integer :: ncid, io

write(string1,*)trim(routine),' "'//trim(filename)//'"'
ncid = nc_open_file_readonly(filename, string1)

! Test if "Time" is a dimension in the file.
isLsmFile = nf90_inq_dimid(ncid, 'Time', DimID) == NF90_NOERR

if(isLsmFile) then ! Get the time from the LSM restart file

   ! TJH ... my preference is to read the dimension IDs for the Times variable
   ! and cycle over the dimension IDs to get the lengths and check compatibility

   ! TJH This has the assumption that Times(DateStrLen,Time)
   ! Get the dimensions for the strings of times
   io = nf90_inq_dimid(ncid, 'Time', DimID)
   call nc_check(io, routine,'inq_dimid','Time',filename)

   io = nf90_inquire_dimension(ncid, DimID, len=ntimes)
   call nc_check(io, routine,'inquire_dimension','Time',filename)

   io = nf90_inq_dimid(ncid, 'DateStrLen', DimID)
   call nc_check(io, routine,'inq_dimid','DateStrLen',filename)

   io = nf90_inquire_dimension(ncid, DimID, len=strlen)
   call nc_check(io, routine,'inquire_dimension','DateStrLen',filename)

   if (strlen /= STRINGLENGTH) then
      write(string1,*)'DatStrLen string length ',strlen,' /= ',STRINGLENGTH
      call error_handler(E_ERR,routine, string1, source )
   endif

   ! Get all the Time strings, use the last one.
   io = nf90_inq_varid(ncid, 'Times', VarID)
   call nc_check(io, routine, 'inq_varid','Times',filename)

   allocate(datestring(ntimes))

   io = nf90_get_var(ncid, VarID, datestring)
   call nc_check(io, routine, 'get_var','Times',filename)

else ! Get the time from the hydro or parameter file

   io = nf90_inquire_attribute(ncid, NF90_GLOBAL, 'Restart_Time', len=strlen)
   call nc_check(io, routine, 'inquire_attribute','Restart_Time', filename)

   io = nf90_get_att(ncid, NF90_GLOBAL, 'Restart_Time', datestring_scalar)
   call nc_check(io, routine, 'get_att','Restart_Time', filename)

   ntimes = 1
   allocate(datestring(ntimes))
   datestring(1) = datestring_scalar

endif

call nc_close_file(ncid, routine, filename)

read(datestring(ntimes),'(i4,5(1x,i2))') year, month, day, hour, minute, second

read_model_time = set_date(year, month, day, hour, minute, second)

if ( do_output() .and. debug > 0 ) then
   write(*,*)'routine: Last time string is '//trim(datestring(ntimes))
   call print_date(read_model_time,' valid date is ')
   call print_time(read_model_time,' valid time is ')
endif

deallocate(datestring)

end function read_model_time


!-----------------------------------------------------------------------
!> The LSM restart files have "time".
!> We are always using the 'most recent' which is, by defn, the last one.
!>
!> RESTART.2004010102_DOMAIN1 has
!>     Times = '2004-01-01_02:00:00' ;
!>
!> The data is valid @ 2004-01-01_02:00:00

subroutine write_model_time(ncid, dart_time)

integer,             intent(in) :: ncid
type(time_type),     intent(in) :: dart_time

character(len=*), parameter :: routine = 'write_model_time'

integer, parameter :: STRINGLENGTH = 19
character(len=STRINGLENGTH) datestring

integer :: year, month, day, hour, minute, second
integer :: VarID, strlen, ntimes
integer :: TimeDimID
integer :: DateStrLenDimID
integer :: ios

integer :: mystart(2), mycount(2)

call get_date(dart_time, year, month, day, hour, minute, second)
write(datestring,'(i4,"-",i2.2,"-",i2.2,"_",i2.2,":",i2.2,":",i2.2)') &
                  year, month, day, hour, minute, second

! TJH This has the assumption that Times(DateStrLen,Time)
! Get the dimensions for the strings of times
ios = nf90_inq_dimid(ncid, 'Time', TimeDimID)
call nc_check(ios, routine, 'inq_dimid', 'Time', ncid=ncid)

ios = nf90_inquire_dimension(ncid, TimeDimID, len=ntimes)
call nc_check(ios, routine, 'inquire_dimension', 'Time', ncid=ncid)

ios = nf90_inq_dimid(ncid, 'DateStrLen', DateStrLenDimID)
call nc_check(ios, routine, 'inq_dimid', 'DateStrLen', ncid=ncid)

ios = nf90_inquire_dimension(ncid, DateStrLenDimID, len=strlen)
call nc_check(ios, routine, 'inquire_dimension', 'DateStrLen', ncid=ncid)

if (strlen /= STRINGLENGTH) then
   write(string1,*)'DatStrLen string length ',strlen,' /= ',STRINGLENGTH
   call error_handler(E_ERR, routine, string1, source )
endif

ios = nf90_inq_varid(ncid, 'Times', VarID)
call nc_check(ios, routine, "inquire Times variable")

mystart = (/ ntimes, 1 /)
mycount = (/ 19, 1 /)

! The first time this routine is called the variable is empty.
if (ntimes == 0) mystart(1) = 1

ios = nf90_put_var(ncid, VarID, datestring, start=mystart, count=mycount)
call nc_check(ios, routine, 'put_var', 'Times', ncid=ncid)

if ( do_output() .and. debug > 0 ) then
   write(*,*)'write_model_time: Last time string is '//datestring
   call print_date(dart_time,' valid date is ')
   call print_time(dart_time,' valid time is ')
endif

end subroutine write_model_time


!------------------------------------------------------------------
!> Returns the smallest increment in time that the model is capable
!> of advancing the state in a given implementation, or the shortest
!> time you want the model to advance between assimilations.
!> This interface is required for all applications.

function shortest_time_between_assimilations()

type(time_type) :: shortest_time_between_assimilations

shortest_time_between_assimilations = &
      set_time(assimilation_period_seconds, assimilation_period_days)

end function shortest_time_between_assimilations


!------------------------------------------------------------------
!> Given a state handle, a location, and a model state variable type,
!> interpolates the state variable field to that location and returns
!> the value in expected_obs(:). The istatus variable should be returned as
!> 0 unless there is some problem in computing the interpolation in
!> which case an alternate value should be returned.
!>
!> istatus = 11 means the DART vector does not have the quantity of interest

subroutine model_interpolate(state_handle, ens_size, location, obs_type, &
                             expected_obs, istatus)

type(ensemble_type), intent(in) :: state_handle
integer,             intent(in) :: ens_size
type(location_type), intent(in) :: location
integer,             intent(in) :: obs_type
real(r8),           intent(out) :: expected_obs(ens_size)
integer,            intent(out) :: istatus(ens_size)

character(len=*), parameter :: routine = 'model_interpolate'

! This is the forward operator -- compute model state value for given point.
! Strategy is to do the horizontal interpolation on the levels above and below
! the observation location and then to do the vertical interpolation.

real(r8) :: loc_lon, loc_lat, loc_depth
real(r8) :: loc(3)
real(r8) :: xloc                       !< fractional longitudinal index
real(r8) :: yloc                       !< fractional latitudinal  index
real(r8) :: zloc                       !< fractional depth/height index
real(r8) :: dx, dy, dz, dxm, dym, dzm  !< fractional weights

integer  :: i, j, k, e

integer  :: domid, varid, obs_kind, num_dims
logical  :: is_lev0

integer, dimension(2) :: ll, lr, ul, ur
integer            :: rc
integer(i8)        :: ill, ilr, iul, iur

real(r8) :: fld(2,ens_size)
real(r8) :: x_iul(ens_size) ,x_iur(ens_size)
real(r8) :: x_ill(ens_size), x_ilr(ens_size)

real(r8) :: layer_midpoints(max(num_soil_layers,num_soil_nitrogen_layers))
integer  :: num_layers

if ( .not. module_initialized ) call static_init_model

! The return code for successful return should be 0.
! Any positive number is an error.
! Negative values are reserved for use by the DART framework.
! Using distinct positive values for different types of errors can be
! useful in diagnosing problems.

istatus      = 1
expected_obs = MISSING_R8
fld          = MISSING_R8

! rename for sanity - we can't change the argument names
! to this subroutine, but this really is a kind.
obs_kind = obs_type

! if observation is outside region encompassed in the history file - fail

loc       = get_location(location) ! loc is in DEGREES
loc_lon   = loc(1)
loc_lat   = loc(2)
loc_depth = loc(3)

! one use of model_interpolate is to allow other modules/routines
! the ability to 'count' the model levels. To do this, we can create
! locations with model levels and 'interpolate' them to
! QTY_GEOPOTENTIAL_HEIGHT

if ( (obs_type == QTY_GEOPOTENTIAL_HEIGHT) .and. is_vertical(location,'LEVEL') ) then
   if (nint(loc_depth) > num_soil_layers) then
      expected_obs = MISSING_R8
   else
      expected_obs = soil_layer_thickness(nint(loc_depth))
      istatus = 0
   endif
   return
endif

! need to know if the DART state has this quantity

DOMAIN: do domid = 1,domain_count
   varid = get_varid_from_kind(domid, obs_kind)
   if (varid < 0) cycle DOMAIN
enddo DOMAIN

if (varid < 0) then
   if (debug > 2) then
      call write_location(0,location,charstring=string3)
      write(string1,*) 'obs kind not found in state.  kind ', obs_kind, &
                       ' (',trim(get_name_for_quantity(obs_kind)),')'
      write(string2,*) ' loc: ', trim(string3)
      call error_handler(E_MSG,routine,string1,text2=string2)
   endif
   istatus(:) = 11
   return ! this kind isn't found in the state vector
endif

!--------------------------------------------------------
! 0. Prelude to Interpolation
!--------------------------------------------------------

! first obtain (innermost) domain ID, and mass points (i,j)
! xloc, yloc are fractional indices ...

call get_domain_info(loc(1), loc(2), domid, xloc, yloc)

if (domid < 1) then
   istatus(:) = 10
   return ! observation outside all domains
endif

if( debug > 99 ) then
   i = xloc
   j = yloc

   write(*,*)'domain found is ',domid
   write(*,*)'loc(1),loc(2),xloc,yloc,i,j',loc(1), loc(2), xloc,yloc,i,j
   write(*,*)'corners of lat:'
   write(*,*) xlat(i,j), xlat(i+1,j), xlat(i,j+1), xlat(i+1,j+1)
   write(*,*)'corners of lon:'
   write(*,*) xlong(i,j),xlong(i+1,j), xlong(i,j+1),xlong(i+1,j+1)
endif

! get integer (west/south) grid point and distances to neighboring grid points
! distances are used as weights to carry out horizontal interpolations

call toGrid(xloc, i, dx, dxm)
call toGrid(yloc, j, dy, dym)

num_dims = get_num_dims(domid, varid)

! Some variables have vertical levels of soil_layers, soil_nitrogen_layers, etc.
call get_vertical_array(domid, varid, num_dims, layer_midpoints, num_layers)

call height_to_zk(loc_depth, layer_midpoints, num_layers, zloc, is_lev0)

if( debug > 99 ) then
   print*,' obs is by height and zloc,lev0 =',zloc, is_lev0
   print*,'soil_depth ', soil_depth
   if (num_soil_nitrogen_layers > 0) &
   print*,'soil_nitrogen_depth', soil_nitrogen_depth
endif

if(zloc == missing_r8) istatus = 2

!----------------------------------
! 1. Horizontal Interpolation
!----------------------------------

! Strategy is to do the horizontal interpolation on two different levels
! and then to do the vertical interpolation afterwards.

call getCorners(i, j, domid, ll, ul, lr, ur, rc )

! Interpolation at level k
if (num_dims == 3) then
   ill = get_dart_vector_index(ll(1), int(zloc), ll(2), domid, varid)
   iul = get_dart_vector_index(ul(1), int(zloc), ul(2), domid, varid)
   ilr = get_dart_vector_index(lr(1), int(zloc), lr(2), domid, varid)
   iur = get_dart_vector_index(ur(1), int(zloc), ur(2), domid, varid)
else
   ill = get_dart_vector_index(ll(1), ll(2), int(zloc), domid, varid)
   iul = get_dart_vector_index(ul(1), ul(2), int(zloc), domid, varid)
   ilr = get_dart_vector_index(lr(1), lr(2), int(zloc), domid, varid)
   iur = get_dart_vector_index(ur(1), ur(2), int(zloc), domid, varid)
endif

x_ill = get_state(ill, state_handle)
x_iul = get_state(iul, state_handle)
x_ilr = get_state(ilr, state_handle)
x_iur = get_state(iur, state_handle)

if (debug > 99) then
   print*,'  ill,   iul,   ilr,   iur',   ill,   iul,   ilr,   iur
   print*,'x_ill, x_iul, x_ilr, x_iur', x_ill, x_iul, x_ilr, x_iur
   print*,'   dx,   dxm,    dy,   dym',    dx,   dxm,    dy,   dym
   print*,'corner weights            ',dym*dxm, dym*dx, dy*dxm, dy*dx
endif

fld(1,:) = dym*( dxm*x_ill + dx*x_ilr ) + dy*( dxm*x_iul + dx*x_iur )
if (debug > 0) print*,'vector comp:',fld(1,:)

where(x_ill == missing_r8) fld(1,:) = missing_r8
where(x_ilr == missing_r8) fld(1,:) = missing_r8
where(x_iul == missing_r8) fld(1,:) = missing_r8
where(x_iur == missing_r8) fld(1,:) = missing_r8

! Interpolation at level k+1

if (num_dims == 3) then
   ill = get_dart_vector_index(ll(1), int(zloc), ll(2), domid, varid)
   iul = get_dart_vector_index(ul(1), int(zloc), ul(2), domid, varid)
   ilr = get_dart_vector_index(lr(1), int(zloc), lr(2), domid, varid)
   iur = get_dart_vector_index(ur(1), int(zloc), ur(2), domid, varid)
else
   ill = get_dart_vector_index(ll(1), ll(2), int(zloc+1), domid, varid)
   iul = get_dart_vector_index(ul(1), ul(2), int(zloc+1), domid, varid)
   ilr = get_dart_vector_index(lr(1), lr(2), int(zloc+1), domid, varid)
   iur = get_dart_vector_index(ur(1), ur(2), int(zloc+1), domid, varid)
endif

x_ill = get_state(ill, state_handle)
x_ilr = get_state(ilr, state_handle)
x_iul = get_state(iul, state_handle)
x_iur = get_state(iur, state_handle)

fld(2,:) = dym*( dxm*x_ill + dx*x_ilr ) + dy*( dxm*x_iul + dx*x_iur )

where(x_ill == missing_r8) fld(2,:) = missing_r8
where(x_ilr == missing_r8) fld(2,:) = missing_r8
where(x_iul == missing_r8) fld(2,:) = missing_r8
where(x_iur == missing_r8) fld(2,:) = missing_r8

!----------------------------------
! 2. Vertical Interpolation
!----------------------------------

! The previous section (1. Horizontal Interpolation) has produced a variable
! called "fld", which nominally has two entries in it. 3D fields have
! hopefully produced 2 non-zero entries, whereas surface fields only have
! filled the first entry.  If a full 3D field, then do vertical interpolation
! between encompassing model levels (k and k+1).

! Do vertical interpolation -- at this point zloc is >= 1

do e = 1,ens_size
   if ( fld(1,e) == missing_r8 ) then
     expected_obs(e) = missing_r8
     istatus(e) = 12
   else
      call toGrid(zloc, k, dz, dzm)
      if (debug > 4) print*, 'zloc, k, dz, dzm = ', zloc, k, dz, dzm
      if (int(zloc) >= 1) then
         ! Linearly interpolate between grid points
         expected_obs(e) = dzm*fld(1,e) + dz*fld(2,e)
         if (debug > 4) print*, 'interpolated obs_val = ', expected_obs(e)
      else
         ! Extrapolate below first level.
         expected_obs(e) = fld(1,e) - (fld(2,e)-fld(1,e))*dzm
         if (debug > 4) print*, 'extrapolated obs_val = ', expected_obs(e)
      endif
      istatus(e) = 0
   endif
enddo

if(debug > 1) then
   do e = 1,ens_size
   write(     *     ,*) 'model_interpolate: member, kind, "obs", status = ', &
              e, obs_kind, expected_obs(e), istatus(e)
   write(logfileunit,*) 'model_interpolate: member, kind, "obs", status = ', &
              e, obs_kind, expected_obs(e), istatus(e)
   enddo
endif

!This is the check routine, skip the Non-uniform part, 
! If one or more of the ensemble members doesn't have snow,
! we must reject the observation because we do not have
! a full rank ensemble for the regression. We don't know
! how to make snow  
if ( (obs_kind == QTY_SNOWCOVER_FRAC) .and. &
     any(expected_obs <= 0.0) )  then
    expected_obs = missing_r8
    istatus      = 12
endif

end subroutine model_interpolate


!------------------------------------------------------------------
!> Given an integer index into the state vector structure, returns the
!> associated location. A second intent(out) optional argument kind
!> can be returned if the model has more than one type of field (for
!> instance temperature and zonal wind component). This interface is
!> required for all filter applications as it is required for computing
!> the distance between observations and state variables.

subroutine get_state_meta_data(index_in, location, var_type)

integer(i8),         intent(in)            :: index_in
type(location_type), intent(out)           :: location
integer,             intent(out), optional :: var_type

character(len=*), parameter :: routine = 'get_state_meta_data'
integer :: idim1, idim2, idim3, varid, domid, qtyid, num_dims
character(len=obstypelength) :: qtystring

character(len=NF90_MAX_NAME) :: varname

real(r8) :: layer_midpoints(max(num_soil_layers,num_soil_nitrogen_layers))
integer  :: num_layers

if ( .not. module_initialized ) call static_init_model

! Be careful with the indices from get_model_variable_indices.
! lat/lon variables are declared (west_east,south_north)
! 2D      variables are declared (west_east,south_north)
! 3D      variables are declared (west_east,soil_layers_stag,south_north)

call get_model_variable_indices(index_in, idim1, idim2, idim3, &
                                varid, domid, qtyid, qtystring)   

num_dims = get_num_dims(domid,varid)
varname  = get_variable_name(domid,varid)

!>@todo Support more than just soil layers.
!> Different vertical dimensions will have different ways of specifying depth.
!> Snow layers are still a problem.

call get_vertical_array(domid, varid, num_dims, layer_midpoints, num_layers)

if ( num_dims == 2 ) then

   location = set_location(xlong(idim1,idim2), xlat(idim1,idim2),  &
                  0.0_r8, VERTISSURFACE)

elseif ( num_dims == 3 ) then

   location = set_location(xlong(idim1,idim3), xlat(idim1,idim3),  &
                  layer_midpoints(idim2), VERTISHEIGHT)

else
   write(string1,*) 'unsupported number of dimensions (',num_dims,') for "',trim(varname),'"'
   call error_handler(E_ERR,'get_state_meta_data',string1,source)
endif

if (present(var_type)) var_type = qtyid

end subroutine get_state_meta_data


!------------------------------------------------------------------
!> Does any shutdown and clean-up needed for model.

subroutine end_model()

deallocate(xlong,xlat)
deallocate(domain_info)
deallocate(soil_depth)

if ( .not. module_initialized ) call static_init_model

end subroutine end_model


!------------------------------------------------------------------
!> This routine writes all the netCDF 'infrastructure' and sets up the
!> global attributes, dimensions, coordinate variables, and output variables.
!> The actuall filling of the output variables is done by
!> nc_write_model_vars() which can be called repeatedly for each
!> assimilation cycle.
!>
!> All errors are fatal, so the return code is always '0 == normal'
!>
!> assim_model_mod:init_diag_output uses information from the location_mod
!>     to define the location dimension and variable ID. All we need to do
!>     is query, verify, and fill ...

subroutine nc_write_model_atts( ncid, domain_id )

integer, intent(in)  :: ncid
integer, intent(in)  :: domain_id

character(len=*), parameter :: routine = 'nc_write_model_atts'

! variables for the geographic metadata

integer :: TimeDimID
integer :: DateStrLenDimID
integer :: weDimID
integer :: snDimID
integer :: nsoilDimID
integer :: nitsoilDimID
integer :: varid

! variables for the namelist output

character(len=129), allocatable, dimension(:) :: textblock
integer :: LineLenDimID, nlinesDimID, nmlVarID
integer :: nlines, linelen
logical :: has_lsm_namelist

integer :: io

if ( .not. module_initialized ) call static_init_model

call nc_begin_define_mode(ncid)

call write_noah_global_atts(ncid)

call nc_add_global_attribute(ncid, "model", "NOAHMP", routine)

! Determine shape of namelist.
! long lines are truncated when read into textblock

call find_textfile_dims(lsm_namelist_filename, nlines, linelen)
if (nlines > 0) then
   has_lsm_namelist = .true.
else
   has_lsm_namelist = .false.
endif

if (has_lsm_namelist) then
   allocate(textblock(nlines))
   textblock = ''

   io= nf90_def_dim(ncid, name='noahNMLnlines', len = nlines, dimid = nlinesDimID)
   call nc_check(io, routine, 'def_dim noahNMLnlines')

   io= nf90_def_dim(ncid, name='line_length', len = 129, dimid = linelenDimID)
   call nc_check(io, routine, 'def_dim line_length')

   io = nf90_def_var(ncid,name=trim(lsm_namelist_filename), xtype=nf90_char, &
          dimids = (/ linelenDimID, nlinesDimID /),  varid=nmlVarID)
   call nc_check(io, routine, 'def_var noah_namelist')

   io= nf90_put_att(ncid, nmlVarID, 'long_name', &
          'contents of '//trim(lsm_namelist_filename))
   call nc_check(io, routine, 'put_att noah_namelist')
endif

! We need to output the prognostic variables.
! Define the additional dimensions IDs
!>@todo multiple domains ... perhaps use the domain_id argument ...

io = nf90_def_dim(ncid, name='Time', len=nf90_unlimited, dimid=TimeDimID)
call nc_check(io, routine, 'def_dim DateStrLen')

io = nf90_def_dim(ncid, name='DateStrLen', len=19, dimid=DateStrLenDimID)
call nc_check(io, routine, 'def_dim DateStrLen')

io = nf90_def_dim(ncid, name='west_east', len=lsm%dom(1)%we, dimid=weDimID)
call nc_check(io, routine,'west_east def_dim')

io = nf90_def_dim(ncid, name='south_north', len=lsm%dom(1)%sn, dimid=snDimID)
call nc_check(io, routine, 'south_north def_dim')

io = nf90_def_dim(ncid, name='soil_layers_stag', len=num_soil_layers, dimid=nsoilDimID)
call nc_check(io, routine, 'def_dim soil_layers_stag')

if (num_soil_nitrogen_layers > 0 ) then
   io = nf90_def_dim(ncid, name='soil_nitrogen_layers_stag', &
             len=num_soil_nitrogen_layers, dimid=nitsoilDimID)
   call nc_check(io, routine, 'def_dim soil_nitrogen_layers_stag')
endif

! Create the (empty) Coordinate Variables and the Attributes

call nc_check(nf90_def_var(ncid, name='Times', xtype=nf90_char, &
              dimids=(/ DateStrLenDimID, TimeDimID /), varid=VarID),&
              routine, 'Times def_var')

! Grid Longitudes
call nc_check(nf90_def_var(ncid,name='XLONG', xtype=nf90_real, &
              dimids=(/ weDimID, snDimID /), varid=VarID),&
              routine, 'XLONG def_var')
call nc_check(nf90_put_att(ncid,  VarID, 'long_name', 'coordinate longitude'), &
              routine, 'XLONG long_name')
call nc_check(nf90_put_att(ncid,  VarID, 'coordinates', 'XLONG XLAT'),  &
              routine, 'XLONG coordinates')
call nc_check(nf90_put_att(ncid,  VarID, 'FieldType', 104),  &
              routine, 'XLONG FieldType')
call nc_check(nf90_put_att(ncid,  VarID, 'MemoryOrder', 'XY'),  &
              routine, 'XLONG MemoryOrder ')
call nc_check(nf90_put_att(ncid,  VarID, '_FillValue', -999.0 ), &
              routine, 'XLONG _FillValue')

! Grid Latitudes
call nc_check(nf90_def_var(ncid,name='XLAT', xtype=nf90_real, &
              dimids=(/ weDimID, snDimID /), varid=VarID),&
              routine, 'XLAT def_var ')
call nc_check(nf90_put_att(ncid,  VarID, 'long_name', 'coordinate latitude'), &
              routine, 'XLAT long_name ')
call nc_check(nf90_put_att(ncid,  VarID, 'coordinates', 'XLONG XLAT'),  &
              routine, 'XLAT coordinates ')
call nc_check(nf90_put_att(ncid,  VarID, 'FieldType', 104),  &
              routine, 'XLAT FieldType ')
call nc_check(nf90_put_att(ncid,  VarID, 'MemoryOrder', 'XY'),  &
              routine, 'XLAT MemoryOrder ')
call nc_check(nf90_put_att(ncid,  VarID, '_FillValue', -999.0 ), &
              routine, 'XLAT _FillValue ')

! subsurface levels
call nc_check(nf90_def_var(ncid,name='soil_layers_stag', xtype=nf90_real, &
              dimids=(/ nsoilDimID /), varid=VarID),&
              routine, 'soil_layers_stag def_var ')
call nc_check(nf90_put_att(ncid,  VarID, 'long_name', 'coordinate soil levels'), &
              routine, 'soil_layers_stag long_name ')
call nc_check(nf90_put_att(ncid,  VarID, 'units', 'm'),  &
              routine, 'soil_layers_stag units ')

! Finished with dimension/variable definitions, must end 'define' mode to fill.

call nc_check(nf90_enddef(ncid), routine, 'prognostic enddef')

! Fill the coordinate variables

call nc_check(nf90_inq_varid(ncid, 'XLONG', VarID), &
             routine, 'inq_varid XLONG ')
call nc_check(nf90_put_var(ncid, VarID, xlong ), &
             routine, 'put_var XLONG ')

call nc_check(nf90_inq_varid(ncid, 'XLAT', VarID), &
             routine, 'inq_varid XLAT ')
call nc_check(nf90_put_var(ncid, VarID, xlat ), &
             routine, 'put_var XLAT ')

call nc_check(nf90_inq_varid(ncid, 'soil_layers_stag', VarID), &
             routine, 'inq_varid soil_layers_stag ')
call nc_check(nf90_put_var(ncid, VarID, soil_layer_thickness(1:num_soil_layers)), &
             routine, 'put_var soil_layers_stag ')


! Fill the variables we can

if (has_lsm_namelist) then
   call file_to_text(lsm_namelist_filename, textblock)
   call nc_check(nf90_put_var(ncid, nmlVarID, textblock ), &
                 routine, 'put_var nmlVarID')
   deallocate(textblock)
endif

! Flush the buffer and leave netCDF file open

call nc_check(nf90_sync(ncid),routine, 'sync')

end subroutine nc_write_model_atts


!------------------------------------------------------------------
!> Perturbs a model state for generating initial ensembles.
!> The perturbed state is referenced by the state_ens_handle.
!> We do not know how to generate an ensemble for an LSM,
!> so this dies right away (for now).

subroutine pert_model_copies(state_ens_handle, ens_size, pert_amp, interf_provided)

type(ensemble_type), intent(inout) :: state_ens_handle
integer,             intent(in)    :: ens_size
real(r8),            intent(in)    :: pert_amp
logical,             intent(out)   :: interf_provided

character(len=*), parameter :: routine = 'pert_model_copies'

logical, save :: seed_unset = .true.
integer, save :: seed
type(random_seq_type) :: random_seq

real(r8) :: stddev, rng
real(r8) :: new_state
integer  :: j, k, copy, ivar
integer  :: start_ind, end_ind
logical  :: positive
integer, parameter :: NTRIES = 100

if ( .not. module_initialized ) call static_init_model

call error_handler(E_MSG,routine, &
                  'NOAH cannot be started from a single vector', source, &
                  text2='see comments in noah/model_mod.f90::pert_model_copies()',&
                  text3='or noah/model_mod.html#pert_model_copies')

interf_provided = .true.

! Generate a unique (but repeatable - if task count is same)
! seed for each ensemble member
if (seed_unset) then
   seed = (my_task_id()+1) * 1000
   seed_unset = .false.
endif

if ( debug > 99 ) then
   write(string1,*)'seed_unset is ',seed_unset,'; seed is ',seed
   call error_handler(E_MSG,routine,string1)
endif

call init_random_seq(random_seq, seed)
seed = seed + 1  ! next ensemble member gets a different seed

stddev = model_perturbation_amplitude

DOMAIN : do idom = 1, domain_count

   !>@todo skip variables and use lognormal for moisture, gaussians for temperatures ...
   do ivar = 1, get_num_variables(idom)

      start_ind = get_index_start(idom, ivar)
      end_ind   = get_index_end(  idom, ivar)

      MINE : do j = 1, state_ens_handle%my_num_vars

         if (state_ens_handle%my_vars(j) >= start_ind .and. &
             state_ens_handle%my_vars(j) <= end_ind   ) then

         COPIES : do copy = 1, ens_size

               if (state_ens_handle%copies(copy,j) == MISSING_R8) cycle COPIES

               if (trim(perturb_distribution) == 'lognormal') then
                  rng = random_gaussian(random_seq, 0.0_r8, 1.0_r8)
                  state_ens_handle%copies(copy,j) = state_ens_handle%copies(copy,j)*exp(stddev*rng)

               else !if it's not lognormal, then the only other option is Gaussian.

                  positive = .false.

                  !>@todo might want to make sure it is within the variable bounds from the namelist
                  POSDEF : do k=1,NTRIES ! prevent runoff from being negative
                     new_state = random_gaussian(random_seq, state_ens_handle%copies(copy,j), stddev)

                     if ( new_state >= 0.0_r8 ) then
                       positive = .true.
                       state_ens_handle%copies(copy,j) = new_state
                       exit POSDEF
                     endif
                  enddo POSDEF
                  if (.not. positive) then
                     write(string1,*)'tried ',NTRIES,' times to get something >= 0.0_r8 and failed'
                     write(string2,*)'state value ',state_ens_handle%copies(copy,j)
                     call error_handler(E_ERR, routine, string1, source, text2=string2)
                  endif
               endif
            enddo COPIES
         endif
      enddo MINE
   enddo
enddo DOMAIN

end subroutine pert_model_copies


!==================================================================
! PUBLIC interfaces that aren't required by the DART code but are
! generally useful for other related utility programs.
! (less necessary for small models; generally used for larger models
! with predefined file formats and control structures.)
!==================================================================


!------------------------------------------------------------------
!> given the list of variables and a filename, check user input
!> return the handle to the open netCDF file and the number of variables
!> in this 'domain'

subroutine verify_variables( variable_table, filename, ngood, &
                       var_names, var_qtys, var_ranges, var_update)

character(len=*), intent(in)  :: variable_table(:,:)
character(len=*), intent(in)  :: filename
integer,          intent(out) :: ngood
character(len=*), intent(out) :: var_names(:)
real(r8),         intent(out) :: var_ranges(:,:)
logical,          intent(out) :: var_update(:)
integer ,         intent(out) :: var_qtys(:)

character(len=*), parameter :: routine = 'verify_variables'

integer  :: io, i, quantity
real(r8) :: minvalue, maxvalue

character(len=NF90_MAX_NAME) :: varname
character(len=NF90_MAX_NAME) :: dartstr
character(len=NF90_MAX_NAME) :: minvalstring
character(len=NF90_MAX_NAME) :: maxvalstring
character(len=NF90_MAX_NAME) :: state_or_aux

ngood = 0
MyLoop : do i = 1, size(variable_table,2)

   varname      = variable_table(VT_VARNAMEINDX,i)
   dartstr      = variable_table(VT_KINDINDX   ,i)
   minvalstring = variable_table(VT_MINVALINDX ,i)
   maxvalstring = variable_table(VT_MAXVALINDX ,i)
   state_or_aux = variable_table(VT_STATEINDX  ,i)

   if ( varname == ' ' .and. dartstr == ' ' ) exit MyLoop ! Found end of list.

   if ( varname == ' ' .or.  dartstr == ' ' ) then
      string1 = 'model_nml: variable list not fully specified'
      string2 = 'reading from "'//trim(filename)//'"'
      call error_handler(E_ERR,routine, string1, source, text2=string2)
   endif

   ! The internal DART routines check if the variable name is valid.

   ! Make sure DART kind is valid
   quantity = get_index_for_quantity(dartstr)
   if( quantity < 0 ) then
      write(string1,'(''there is no obs_kind "'',a,''" in obs_kind_mod.f90'')') &
                    trim(dartstr)
      call error_handler(E_ERR,routine,string1,source)
   endif

   ! All good to here - fill the output variables

   ngood = ngood + 1
   var_names( ngood)   = varname
   var_qtys(  ngood)   = quantity
   var_ranges(ngood,:) = (/ MISSING_R8, MISSING_R8 /)
   var_update(ngood)   = .false.   ! at least initially

   ! convert the [min,max]valstrings to numeric values if possible
   read(minvalstring,*,iostat=io) minvalue
   if (io == 0) var_ranges(ngood,1) = minvalue

   read(maxvalstring,*,iostat=io) maxvalue
   if (io == 0) var_ranges(ngood,2) = maxvalue

   call to_upper(state_or_aux)
   if (state_or_aux == 'UPDATE') var_update(ngood) = .true.

enddo MyLoop

if (ngood == MAX_STATE_VARIABLES) then
   string1 = 'WARNING: you may need to increase "MAX_STATE_VARIABLES"'
   write(string2,'(''you have specified at least '',i4,'' perhaps more.'')') ngood
   call error_handler(E_MSG,routine,string1,source,text2=string2)
endif


end subroutine verify_variables


!-----------------------------------------------------------------------
!> Sets the location information arrays for each domain
!> Each location array is declared to be 3D to make it easy to use
!> the state_structure routines.

subroutine configure_domains()

! the xlong, xlat, west_east, south_north variables are needed later,
! so we keep them around in module storage
!>@todo should be part of the per-domain structure

character(len=*), parameter :: routine = 'configure_domains'

integer :: i, j

logical, save :: initialized = .false.

if (initialized) return ! can only call this function ONCE.

initialized = .true.

do idom = 1,domain_count

   ! determine the wrfinput file for this domain
   call get_lsm_domain_filename(idom, lsm%dom(idom)%filename)

   ! fills the lsm%dom(idom) structure
   call read_wrf_file_attributes(idom)

   !>@todo should have xlong,xlat one per domain ...
   allocate(xlong(lsm%dom(idom)%we, lsm%dom(idom)%sn))
   allocate( xlat(lsm%dom(idom)%we, lsm%dom(idom)%sn))

   !>@todo re-evlauate
   call get_lsm_domain_info(west_east,south_north,xlong,xlat)

   !>@todo do we even need a DART location type representation of xlong,xlat
   allocate( domain_info(idom)%location(west_east,south_north) )

   do j = 1,south_north
   do i = 1,west_east
       domain_info(idom)%location(i,j) = &
            set_location(xlong(i,j), xlat(i,j), 0.0_r8, VERTISSURFACE)
   enddo
   enddo

   ! Get model projection information
   call setup_map_projection(idom)

enddo

end subroutine configure_domains


!-----------------------------------------------------------------------
!> only support variables with a vertical coordinate of 'soil_layers_stag'
!>@todo this will have to change when we support snow layers ...

subroutine check_vertical_dimension(domainID)

integer, intent(in) :: domainID

character(len=*), parameter  :: routine = 'check_vertical_dimension'
character(len=obstypelength) :: varname, dimname
integer :: ivar,i

do ivar = 1,get_num_variables(domainID)

   varname = get_variable_name(domainID,ivar)

   do i = 1, get_num_dims(domainID, ivar)

      dimname = get_dim_name(domainID, ivar, i)

      if(dimname /= 'south_north' .and. &
         dimname /= 'west_east'   .and. &
         dimname /= 'soil_layers_stag' .and. &
         dimname /= 'soil_nitrogen_layers_stag') then
         write(string1,*) 'Only supporting variables with a vertical coordinate &
                           &of "soil_layers_stag" or "soil_nirogen_layers_stag"'
         write(string2,*) trim(varname), ' has "',trim(dimname),'"'
         call error_handler(E_ERR, routine, string1, source, text2=string2)
      endif
   enddo
enddo

end subroutine check_vertical_dimension


!-----------------------------------------------------------------------
!> Get the global projection attributes from a 'wrfinput' file.
!> This routine has an analogue in the wrf model_mod.f90

subroutine read_wrf_file_attributes(dom_id)

integer, intent(in) :: dom_id

character(len=*), parameter :: routine = 'read_wrf_file_attributes'
integer :: io, ncid, DimID
real(r8) :: dummy

! set the metadata that should be in the netCDF file but is being read from the module

lsm%dom(dom_id)%periodic_x = periodic_x
lsm%dom(dom_id)%periodic_y = periodic_y
lsm%dom(dom_id)%polar      = polar

! Only a subset of the wrf grid logic is implemented in this (noah) 
! If people are going across the prime meridian or ... that logic is not
! implemented yet. 

if (periodic_x .or. periodic_y .or. polar) then

   write(string1,*)'Only a subset of the wrf grid logic is supported in our noah implementation.'
   write(string2,*)'We do not support periodic boundary conditions. Unsupported namelist setting:'
   write(string3,*)'periodic_x=',periodic_x,'; periodic_y=',periodic_y,'; polar=',polar
   call error_handler(E_ERR, routine, string1, source, text2=string2, text3=string3)
endif

! get meta data and static data we need

write(string1,*)trim(routine),' "'//trim(lsm%dom(dom_id)%filename)//'"'
ncid = nc_open_file_readonly(lsm%dom(dom_id)%filename, string1)

io = nf90_inq_dimid(ncid, 'south_north', DimID)
call nc_check(io, routine, 'inq_dimid','south_north',ncid=ncid)
io = nf90_inquire_dimension(ncid, DimID, len=lsm%dom(dom_id)%sn)
call nc_check(io, routine, 'inquire_dimension','south_north',ncid=ncid)

io = nf90_inq_dimid(ncid, 'west_east', DimID)
call nc_check(io, routine, 'inq_dimid','west_east',ncid=ncid)
io = nf90_inquire_dimension(ncid, DimID, len=lsm%dom(dom_id)%we)
call nc_check(io, routine, 'inquire_dimension','west_east',ncid=ncid)

call nc_get_global_attribute(ncid, 'DX',        lsm%dom(dom_id)%dx, routine)
call nc_get_global_attribute(ncid, 'DY',        lsm%dom(dom_id)%dy, routine)
call nc_get_global_attribute(ncid, 'MAP_PROJ',  lsm%dom(dom_id)%map_proj, routine)
call nc_get_global_attribute(ncid, 'TRUELAT1',  lsm%dom(dom_id)%truelat1, routine)
call nc_get_global_attribute(ncid, 'TRUELAT2',  lsm%dom(dom_id)%truelat2, routine)
call nc_get_global_attribute(ncid, 'STAND_LON', lsm%dom(dom_id)%stand_lon, routine)
call nc_get_global_attribute(ncid, 'POLE_LAT',  dummy, routine)
call nc_get_global_attribute(ncid, 'POLE_LON',  dummy, routine)

if(debug > 0) then
   write(*,*) ' filename    "'//trim(lsm%dom(dom_id)%filename)//'"'
   write(*,*) ' south_north ',lsm%dom(dom_id)%sn
   write(*,*) ' west_east   ',lsm%dom(dom_id)%we
   write(*,*) ' dx          ',lsm%dom(dom_id)%dx
   write(*,*) ' dy          ',lsm%dom(dom_id)%dy
   write(*,*) ' map_proj    ',lsm%dom(dom_id)%map_proj
   write(*,*) ' truelat1    ',lsm%dom(dom_id)%truelat1
   write(*,*) ' truelat2    ',lsm%dom(dom_id)%truelat2
   write(*,*) ' stand_lon      ',lsm%dom(dom_id)%stand_lon
   write(*,*) ' pole_lat    ',lsm%dom(dom_id)%pole_lat
   write(*,*) ' pole_lon    ',lsm%dom(dom_id)%pole_lon

   write(logfileunit,*) ' filename    "'//trim(lsm%dom(dom_id)%filename)//'"'
   write(logfileunit,*) ' south_north ',lsm%dom(dom_id)%sn
   write(logfileunit,*) ' west_east   ',lsm%dom(dom_id)%we
   write(logfileunit,*) ' dx          ',lsm%dom(dom_id)%dx
   write(logfileunit,*) ' dy          ',lsm%dom(dom_id)%dy
   write(logfileunit,*) ' map_proj    ',lsm%dom(dom_id)%map_proj
   write(logfileunit,*) ' truelat1    ',lsm%dom(dom_id)%truelat1
   write(logfileunit,*) ' truelat2    ',lsm%dom(dom_id)%truelat2
   write(logfileunit,*) ' stand_lon      ',lsm%dom(dom_id)%stand_lon
   write(logfileunit,*) ' pole_lat    ',lsm%dom(dom_id)%pole_lat
   write(logfileunit,*) ' pole_lon    ',lsm%dom(dom_id)%pole_lon
endif

call nc_close_file(ncid, routine)

end subroutine read_wrf_file_attributes


!------------------------------------------------------------------------
!> Initialize the map projection structure
!> Populate the map projection structure with desired values

subroutine setup_map_projection(dom_id)

integer, intent(in) :: dom_id

character(len=*), parameter :: routine = 'setup_map_projection'

integer  :: proj_code
real(r8) :: latinc, loninc

call map_init(lsm%dom(dom_id)%proj)

if ( lsm%dom(dom_id)%scm ) then
   ! JPH -- set to zero which should cause the map utils to return NaNs if called
   latinc = 0.0_r8
   loninc = 0.0_r8
   write(string1,*)'Single Column Model not supported.'
   call error_handler(E_ERR, routine, string1, source)
else
   latinc = 180.0_r8/lsm%dom(dom_id)%sn
   loninc = 360.0_r8/lsm%dom(dom_id)%we
endif

if(lsm%dom(dom_id)%map_proj == map_latlon) then
   ! latinc and loninc should be the interval between two neighbors
   latinc   =  xlat(1,2) -  xlat(1,1)
   loninc   = xlong(2,1) - xlong(1,1)
   lsm%dom(dom_id)%truelat1 = latinc
   lsm%dom(dom_id)%stand_lon   = loninc
   proj_code = map_latlon
elseif(lsm%dom(dom_id)%map_proj == map_lambert) then
   proj_code = map_lambert
elseif(lsm%dom(dom_id)%map_proj == map_polar_stereo) then
   proj_code = map_polar_stereo
elseif(lsm%dom(dom_id)%map_proj == map_mercator) then
   proj_code = map_mercator
elseif(lsm%dom(dom_id)%map_proj == map_cylindrical) then
   proj_code = map_cylindrical
elseif(lsm%dom(dom_id)%map_proj == map_cassini) then
   proj_code = map_cassini
else
   write(string1,*)'Map projection ',lsm%dom(dom_id)%map_proj,' not supported.'
   write(string2,*)'Valid values are: ',map_latlon, map_lambert, map_polar_stereo
   write(string3,*)'                  ',map_mercator, map_cylindrical, map_cassini
   call error_handler(E_ERR, routine, string1, source, text2=string2, text3=string3)
endif

if (debug > 99) then
   print*,'projection code is ',proj_code
   print*,'xlat(1,1),xlong(1,1)',xlat(1,1),xlong(1,1)
   print*,'latinc,loninc', latinc, loninc
endif

! Throw an error if using a rotated pole. That option is not advised.
!>@todo include filename in message

if (lsm%dom(dom_id)%pole_lat /= 90.0_r8 .or. &
    lsm%dom(dom_id)%pole_lon /=  0.0_r8 ) then
   write(string1,*)'Rotated poles are not supported.'
   write(string2,*)'POLE_LAT must be 90.0'
   write(string3,*)'POLE_LON must be  0.0'
   call error_handler(E_ERR, routine, string1, source, text2=string2, text3=string3)
endif

call map_set( proj_code = proj_code, &
              proj      = lsm%dom(dom_id)%proj, &
              lat1      = xlat(1,1), &
              lon1      = xlong(1,1), &
              lat0      = lsm%dom(dom_id)%pole_lat, &
              lon0      = lsm%dom(dom_id)%pole_lon, &
              knowni    = 1.0_r8, &
              knownj    = 1.0_r8, &
              dx        = lsm%dom(dom_id)%dx, &
              latinc    = latinc, &
              loninc    = loninc, &
              stdlon    = lsm%dom(dom_id)%stand_lon, &
              truelat1  = lsm%dom(dom_id)%truelat1, &
              truelat2  = lsm%dom(dom_id)%truelat2 )

end subroutine setup_map_projection


!------------------------------------------------------------------------
!> given arbitrary lat and lon values, returns closest domain id and
!> horizontal mass point grid points (xloc,yloc)

subroutine get_domain_info(obslon, obslat, dom_id, iloc, jloc, domain_id_start)

real(r8), intent(in)           :: obslon, obslat
integer,  intent(out)          :: dom_id
real(r8), intent(out)          :: iloc, jloc
integer,  intent(in), optional :: domain_id_start

character(len=*), parameter :: routine = 'get_domain_info'
logical :: dom_found

! the default is to start at the innermost domain and stop when
! the location is found.  however if you want to start at a particular
! domain id number, pass it in as the last optional arg.

dom_id = domain_count
if (present(domain_id_start)) then
   if (domain_id_start < 1 .or. domain_id_start > domain_count) then
      write(string1,  '(A,I1)') 'bad domain_id_start: ', domain_id_start
      write(string2, '(A,I1)') 'must be between 1 and ', domain_count
      call error_handler(E_ERR, routine, string1, source, text2=string2)
   endif
   dom_id = domain_id_start
endif

dom_found = .false.

do while (.not. dom_found)

   ! Checking for exact equality on real variable types is generally a bad idea.

   if(  abs(obslat) > 90.0_r8 ) then

      ! catch latitudes that are out of range - ignore them but print out a warning.
      write(string1, *) 'obs with latitude out of range: ', obslat
      call error_handler(E_MSG, routine, string1)

   else

      !>@todo is the min(max()) bit necessary given the abs(oblat) just above
      call latlon_to_ij(lsm%dom(dom_id)%proj, &
               min(max(obslat,-89.9999999_r8),89.9999999_r8),obslon,iloc,jloc)

      ! Array bound checking depends on whether periodic or not -- these are
      !   real-valued indices here, so we cannot use boundsCheck  :(

      if ( lsm%dom(dom_id)%periodic_x .and. .not. lsm%dom(dom_id)%periodic_y  ) then
         if ( lsm%dom(dom_id)%polar ) then
            !   Periodic     X & M_grid ==> [1 we+1)
            !   Periodic     Y & M_grid ==> [0.5 sn+0.5]
            if ( iloc >= 1.0_r8 .and. iloc <  real(size(xlong(:,1)),r8)+1.0_r8 .and. &
                 jloc >= 0.5_r8 .and. jloc <= real(size(xlong(1,:)),r8)+0.5_r8 ) &
                 dom_found = .true.
         else
            !   Periodic     X & M_grid ==> [1 we+1)
            !   NOT Periodic Y & M_grid ==> [1 sn]
            if ( iloc >= 1.0_r8 .and. iloc <  real(size(xlong(:,1)),r8)+1.0_r8 .and. &
                 jloc >= 1.0_r8 .and. jloc <= real(size(xlong(1,:)),r8) ) &
                 dom_found = .true.
         endif

      elseif ( lsm%dom(dom_id)%periodic_x .and. lsm%dom(dom_id)%periodic_y ) then
            !   Periodic     X & M_grid ==> [1 we+1)
            !   Periodic     Y & M_grid ==> [1 sn+1]
            if ( iloc >= 1.0_r8 .and. iloc <  real(size(xlong(:,1)),r8)+1.0_r8 .and. &
                 jloc >= 1.0_r8 .and. jloc <= real(size(xlong(1,:)),r8)+1.0_r8 ) &
                 dom_found = .true.

      else
         if ( lsm%dom(dom_id)%polar ) then
            !   NOT Periodic X & M_grid ==> [1 we]
            !   Periodic     Y & M_grid ==> [0.5 sn+0.5]
            if ( iloc >= 1.0_r8 .and. iloc <= real(size(xlong(:,1)),r8) .and. &
                 jloc >= 0.5_r8 .and. jloc <= real(size(xlong(1,:)),r8)+0.5_r8 ) &
                 dom_found = .true.
         else
            !   NOT Periodic X & M_grid ==> [1 we]
            !   NOT Periodic Y & M_grid ==> [1 sn]
            if ( iloc >= 1.0_r8 .and. iloc <= real(size(xlong(:,1)),r8) .and. &
                 jloc >= 1.0_r8 .and. jloc <= real(size(xlong(1,:)),r8) ) &
                 dom_found = .true.
         endif
      endif

   endif

   if (.not. dom_found) then  ! check the next domain
      dom_id = dom_id - 1
      if (dom_id == 0) return
   endif

enddo

end subroutine get_domain_info

!------------------------------------------------------------------------
!> Transfer obs. x to grid j and calculate its
!> distance to grid j and j+1

subroutine toGrid (x, j, dx, dxm)

real(r8), intent(in)  :: x
real(r8), intent(out) :: dx
real(r8), intent(out) :: dxm
integer,  intent(out) :: j

j   = int(x)
dx  = x - real(j)
dxm = 1.0_r8 - dx

end subroutine toGrid


!------------------------------------------------------------------------
!>  Calculate the model level "zk" on half (mass) levels,
!>  corresponding to height "obs_v".

subroutine height_to_zk(obs_v, mdl_v, n3, zk, lev0)

  real(r8), intent(in)  :: obs_v
  integer,  intent(in)  :: n3
  real(r8), intent(in)  :: mdl_v(1:n3)
  real(r8), intent(out) :: zk
  logical,  intent(out) :: lev0

  integer   :: k

  zk = missing_r8
  lev0 = .false.

  ! if out of range completely, return missing_r8 and lev0 false
  if ( obs_v > mdl_v(n3)) return

  ! if above surface but below lowest 3-d height level, return the
  ! height value but set lev0 true
  if(obs_v >= 0.0_r8 .and. obs_v < mdl_v(1)) then
    lev0 = .true.
    zk = 1.0_r8
    return
  endif

  ! find the 2 height levels the value is between and return that
  ! as a real number, including the fraction between the levels.
  do k = 1,n3-1
     if(obs_v >= mdl_v(k) .and. obs_v <= mdl_v(k+1)) then
        zk = real(k) + (mdl_v(k) - obs_v)/(mdl_v(k) - mdl_v(k+1))
        exit
     endif
  enddo

end subroutine height_to_zk

!-------------------------------------------------------------------
!> takes in an i and j index, information about domain and grid staggering,
!> and then returns the four cornering gridpoints' 2-element integer indices.
!>
!> This routine is a derivative of a routine by the same name in the
!> wrf model_mod.f90

subroutine getCorners(i, j, id, ll, ul, lr, ur, rc)

integer, intent(in)  :: i, j, id
integer, intent(out) :: ll(2), ul(2), lr(2), ur(2)
integer, intent(out) :: rc

character(len=*), parameter :: routine = 'getCorners'
logical, parameter :: restrict_polar = .false.

! set return code to 0, and change this if necessary
rc = 0

!----------------
! LOWER LEFT i and j are the lower left (ll) corner already
!----------------

! NOTE :: once we allow for polar periodicity, the incoming j index could actually
!           be 0, which would imply a ll(2) value of 1, with a ll(1) value 180 degrees
!           of longitude away from the incoming i index!  But we have not included
!           this possibility yet.

! As of 22 Oct 2007, this option is not allowed!
!   Note that j = 0 can only happen if we are on the M (or U) wrt to latitude

if ( lsm%dom(id)%polar .and. j == 0 .and. .not. restrict_polar ) then
   ! j = 0 should be mapped to j = 1 (ll is on other side of globe)
   ll(2) = 1
   ! Need to map i index 180 degrees away
   ll(1) = i + size(xlong(1,:))/2
   ! Check validity of bounds & adjust by periodicity if necessary
   if ( ll(1) > size(xlong(:,1)) ) ll(1) = ll(1) - size(xlong(:,1))
else
   ll(1) = i
   ll(2) = j
endif

!----------------
! LOWER RIGHT
!----------------

if ( lsm%dom(id)%periodic_x ) then
   write(string1,*)'not supporting periodic x just yet'
   call error_handler(E_ERR,routine,string1,source)
endif

! Regardless of grid, NOT Periodic always has i+1

lr(1) = ll(1) + 1
lr(2) = ll(2)

!----------------
! UPPER LEFT
!----------------
! Regardless of grid, NOT Periodic always has j+1

ul(1) = ll(1)
ul(2) = ll(2) + 1

!----------------
! UPPER RIGHT
!----------------

ur(2) = ul(2)

! Need to check if ur(1) .ne. lr(1)
if ( lsm%dom(id)%polar .and. .not. restrict_polar ) then
   ! Only if j == 0 or j == sn
   if ( j == 0 .or. j == size(xlong(1,:)) ) then
      ! j == 0 means that ll(1) = i + 180 deg, so we cannot use lr(1) -- hence, we will
      !   add 1 to ul(1), unless doing so spans the longitude seam point.
      ! j == sn means that ul(1) = i + 180 deg.  Here we cannot use lr(1) either because
      !   it will be half a domain away from ul(1)+1.  Be careful of longitude seam point.

      !   Here we need to check longitude periodicity and the type of grid
   ! If not a special j value, then we are set for the ur(1) = lr(1)
   else
      ur(1) = lr(1)
   endif
! If not an unrestricted polar periodic domain, then we have nothing to worry about
else
   ur(1) = lr(1)
endif

end subroutine getCorners

!-------------------------------------------------------------------

subroutine get_vertical_array(domid, varid, num_dims, layer_midpoints, num_layers)

integer,          intent(in) :: domid
integer,          intent(in) :: varid
integer,          intent(in) :: num_dims
real(r8),         intent(out) :: layer_midpoints(:)
integer,          intent(out) :: num_layers

integer :: jdim
character(len=NF90_MAX_NAME) :: dimname

! Assume everything is at the surface until proven different

layer_midpoints = 0.0_r8
num_layers      = 1

VERTICAL : do jdim = 1, num_dims
   dimname = get_dim_name(domid, varid, jdim)
   select case (dimname)
      case ('soil_layers_stag')
         layer_midpoints = soil_depth
         num_layers      = num_soil_layers
         exit VERTICAL

      case ('soil_nitrogen_layers_stag')
         if (num_soil_nitrogen_layers > 0) then
            layer_midpoints = soil_nitrogen_depth
            num_layers      = num_soil_nitrogen_layers
         else
            write(string1,*)'variable uses soil_nitrogen_layers but'
            write(string2,*)'namelist.hrldas has no nitrogen layers (NITSOIL).'
            call error_handler(E_ERR,'get_vertical_array',string1, source, text2=string2)
         endif
         exit VERTICAL

   end select
enddo VERTICAL

end subroutine get_vertical_array


end module model_mod

