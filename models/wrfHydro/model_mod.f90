! DART software - Copyright UCAR. This open source software is provided
! by UCAR, "as is", without charge, subject to all terms of use at
! http://www.image.ucar.edu/DAReS/DART/DART_download
!
! DART $Id$

module model_mod

!> This is the interface for wrfHydro and the DART data assimilation infrastructure. 

use              types_mod, only : r8, i8, MISSING_R8, obstypelength

use      time_manager_mod, only : time_type, set_time, get_time, set_time_missing, &
                                  set_date, print_time, print_date, set_calendar_type

use          location_mod, only : location_type, get_close_type, get_dist, &
                                  loc_get_close_obs => get_close_obs, &
                                  loc_get_close_state => get_close_state, &
                                  convert_vertical_obs, convert_vertical_state, &
                                  set_location, set_location_missing, VERTISHEIGHT, &
                                  write_location

use         utilities_mod, only : register_module, error_handler, &
                                  E_ERR, E_MSG, file_exist, logfileunit, &
                                  nmlfileunit, do_output, do_nml_file, do_nml_term, &
                                  find_namelist_in_file, check_namelist_read, &
                                  to_upper

use  netcdf_utilities_mod, only : nc_check, nc_add_global_attribute, &
                                  nc_synchronize_file, nc_end_define_mode, &
                                  nc_add_global_creation_time, nc_begin_define_mode

use obs_def_utilities_mod, only : track_status

use         obs_kind_mod,  only : get_index_for_quantity, get_name_for_quantity, &
                                  QTY_BUCKET_MULTIPLIER, QTY_RUNOFF_MULTIPLIER, &
                                  SOIL_MOISTURE

use  ensemble_manager_mod, only : ensemble_type

use distributed_state_mod, only : get_state

use      dart_time_io_mod, only : write_model_time

use     default_model_mod, only : adv_1step, nc_write_model_vars

use        noah_hydro_mod, only : configure_lsm, configure_hydro, &
                                  n_link, linkLong, linkLat, linkAlt, get_link_tree, &
                                  full_to_connection, get_downstream_links

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

use        random_seq_mod, only : random_seq_type, init_random_seq, random_gaussian

use netcdf

implicit none
private

! Required by DART code - will be called from filter and other DART executables
! These interfaces are fixed and cannot be changed in any way.

public :: static_init_model,      &
          init_time,              &
          init_conditions,        &
          pert_model_copies,      &
          get_model_size,         &
          read_model_time,        &
          shortest_time_between_assimilations, &
          get_state_meta_data,    &
          model_interpolate,      &
          get_close_obs,          &
          get_close_state,        &
          nc_write_model_atts,    &
          end_model

! Also required, but the default routines are sufficient.

public :: nc_write_model_vars,    &
          adv_1step,              &
          convert_vertical_obs,   &
          convert_vertical_state, &
          write_model_time

! version controlled file description for error handling, do not edit
character(len=*), parameter :: source   = &
   "$URL$"
character(len=*), parameter :: revision = "$Revision$"
character(len=*), parameter :: revdate  = "$Date$"

logical, save :: module_initialized = .false.

character(len=512) :: string1, string2, string3

integer(i8) :: model_size
type(time_type) :: time_step

! Codes for interpreting the columns of the variable_table
integer, parameter :: VT_VARNAMEINDX  = 1 ! ... variable name
integer, parameter :: VT_KINDINDX     = 2 ! ... DART kind
integer, parameter :: VT_MINVALINDX   = 3 ! ... minimum value if any
integer, parameter :: VT_MAXVALINDX   = 4 ! ... maximum value if any
integer, parameter :: VT_STATEINDX    = 5 ! ... update (state) or not
integer, parameter :: MAX_STATE_VARIABLES     = 40
integer, parameter :: NUM_STATE_TABLE_COLUMNS = 5

integer :: domain_count
integer :: idom, idom_hydro = -1, idom_parameters = -1, idom_lsm = -1

! Model namelist declarations with defaults 
integer             :: assimilation_period_days     = 0
integer             :: assimilation_period_seconds  = 3600
character(len=128)  :: lsm_model_choice             = 'noahMP'
character(len=256)  :: domain_shapefiles(3)         = ''
integer             :: debug                        = 0
real(r8)            :: model_perturbation_amplitude = 0.002
real(r8)            :: max_link_distance            = 10000.0   ! 10 km
character(len=256)  :: perturb_distribution         = 'lognormal'

character(len=obstypelength) :: &
    lsm_variables(  NUM_STATE_TABLE_COLUMNS,MAX_STATE_VARIABLES) = '', &
    hydro_variables(NUM_STATE_TABLE_COLUMNS,MAX_STATE_VARIABLES) = '', &
    parameters(     NUM_STATE_TABLE_COLUMNS,MAX_STATE_VARIABLES) = ''

namelist /model_nml/ assimilation_period_days,     &
                     assimilation_period_seconds,  &
                     lsm_model_choice,             &
                     domain_shapefiles,            &
                     debug,                        & 
                     model_perturbation_amplitude, &
                     perturb_distribution,         &
                     max_link_distance,            &
                     lsm_variables,                &
                     hydro_variables,              &
                     parameters

type domain_locations
   private
   type(location_type), allocatable :: location(:,:,:)
end type domain_locations

type(domain_locations), allocatable :: domain_info(:)

contains

!------------------------------------------------------------------
!> One-time initialization of the model. This reads the namelists
!> and determines what 'domains' are being used and fills the
!> appropriate metadata for the module (location arrays, etc).

subroutine static_init_model()

character(len=*), parameter :: routine = 'static_init_model'

integer  :: iunit, io, domainID
integer  :: n_lsm_fields
integer  :: n_hydro_fields
integer  :: n_parameters

character(len=obstypelength) :: var_names(MAX_STATE_VARIABLES)
real(r8) :: var_ranges(MAX_STATE_VARIABLES,2)
logical  :: var_update(MAX_STATE_VARIABLES)
integer  :: var_qtys(  MAX_STATE_VARIABLES)

character(len=256) :: filename

if ( module_initialized ) return ! only need to do this once

! Print module information to log file and stdout.
call register_module(source, revision, revdate)

! Since this routine calls other routines that could call this routine
! we'll say we've been initialized pretty dang early.
module_initialized = .true.

! Read the DART namelist, and record.
call find_namelist_in_file("input.nml", "model_nml", iunit)
read(iunit, nml = model_nml, iostat = io)
call check_namelist_read(iunit, io, "model_nml")

if (do_nml_file()) write(nmlfileunit, nml=model_nml)
if (do_nml_term()) write(     *     , nml=model_nml)

! Set the time aspects of the model
call set_calendar_type('Gregorian')
time_step  = set_time(assimilation_period_seconds, assimilation_period_days)

! Determine the composition of the DART state vector.

if ( .not. file_exist(domain_shapefiles(1)) ) then 
   write(string1,*) 'Domain 1 shapefile "', trim(domain_shapefiles(1)),'" does not exist.'
   write(string2,*) 'There must be a domain 1.'
   call error_handler(E_ERR,routine,string1,source,revision,revdate,text2=string2)
endif

DOMAINS: do domainID = 1,size(domain_shapefiles)

   if (len_trim(domain_shapefiles(domainID)) == 0) exit DOMAINS

   filename = domain_shapefiles(domainID)
   call to_upper(filename)

   if ( .not. file_exist(domain_shapefiles(domainID)) ) then 
      write(string1,'("Domain shapefile ",i3, A," does not exist.")')  &
                     domainID, '"'//trim(domain_shapefiles(domainID))//'"'
      call error_handler(E_ERR,routine,string1,source,revision,revdate)
   endif

   if (index(filename,'HYDRO') > 0) then 

      call configure_hydro()
      call verify_variables(hydro_variables, domain_shapefiles(domainID), n_hydro_fields, &
                       var_names, var_qtys, var_ranges, var_update)
      idom_hydro = add_domain(domain_shapefiles(domainID), n_hydro_fields, var_names, &
                         kind_list=var_qtys, &
                        clamp_vals=var_ranges(1:n_hydro_fields,:), &
                       update_list=var_update)

      if (debug > 99) call state_structure_info(idom_hydro)

   elseif (index(filename,'PARAMETER') > 0) then 

      call verify_variables(parameters, domain_shapefiles(domainID), n_parameters, &
                       var_names, var_qtys, var_ranges, var_update)
      idom_parameters = add_domain(domain_shapefiles(domainID), n_parameters, &
                       var_names, kind_list=var_qtys, clamp_vals=var_ranges(1:n_parameters,:), &
                       update_list=var_update )
      if (debug > 99) call state_structure_info(idom_parameters)

   elseif (index(filename,'LSM') > 0) then 

      call configure_lsm()
      call verify_variables(lsm_variables, domain_shapefiles(domainID), n_lsm_fields, &
                       var_names, var_qtys, var_ranges, var_update)
      idom_lsm = add_domain(domain_shapefiles(domainID), n_lsm_fields, var_names, &
                         kind_list=var_qtys, &
                        clamp_vals=var_ranges(1:n_lsm_fields,:), &
                       update_list=var_update)
      if (debug > 99) call state_structure_info(idom_lsm)
   endif
enddo DOMAINS

if (idom_parameters > 0 .and. idom_hydro < 1) then
   write(string1,*) 'Since we are getting the parameter locations from the hydro domain'
   write(string2,*) 'the hydro shapefile must exist.'
   call error_handler(E_ERR,'init_time',string1,source,revision,revdate,text2=string2)
endif

! Now that we know the composition of the DART state, determine the model size.

domain_count = get_num_domains()
model_size = 0
do idom = 1,domain_count
   model_size = model_size + get_domain_size(idom)
enddo

! Each domain has its own array of locations. 

allocate(domain_info(domain_count))
call configure_domains()

end subroutine static_init_model


!------------------------------------------------------------------
!> Returns a time that is somehow appropriate for starting up a new
!> instance of the model. This is only used if the namelist parameter 
!> start_from_restart is set to .false. in the program perfect_model_obs.

subroutine init_time(time)

type(time_type), intent(out) :: time

! for now, just set to 0
time = set_time(0,0)

write(string1,*) 'no good way to specify initial time'
call error_handler(E_ERR,'init_time',string1,source,revision,revdate)

end subroutine init_time


!------------------------------------------------------------------
!> Companion to init_time(). Returns a model state vector, x, that is
!> an appropriate initial condition for starting an integration of the model.
!> At present, this is only used if the namelist parameter 
!> start_from_restart is set to .false. in the program perfect_model_obs.

subroutine init_conditions(x)

real(r8), intent(out) :: x(:)

x = MISSING_R8

write(string1,*) 'no good way to specify initial conditions'
call error_handler(E_ERR,'init_conditions',string1,source,revision,revdate)

end subroutine init_conditions


!------------------------------------------------------------------
!> Perturbs a model state for generating initial ensembles.
!> The perturbed state is referenced by the state_ens_handle.
!> We do not know how to generate an ensemble for an LSM,
!> so this will skip that domain.

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

if ( idom_lsm > 0 ) then
   string1 = 'LSM cannot be started from a single vector - not perturbing that domain.'
   string2 = 'LSM cannot be started from a single vector - not perturbing that domain.'
   string3 = 'LSM cannot be started from a single vector - not perturbing that domain.'
   call error_handler(E_MSG, routine, string1, text2=string2, text3=string3)
endif

! If you want the default model perturbations specified by the
! ensemble_manager (which uses input.nml:ensemble_manager_nml:perturbation_amplitude)
! interf_provided = .false. and return.

interf_provided = .true.   ! .true. means we have our own way of perturbing

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

   if (idom == idom_lsm) cycle DOMAIN

   do ivar = 1, get_num_variables(idom)

      ! This is where we can skip variables we do not want to perturb
      ! or change the stddev based on the variable dynamic range
      ! for the dynamic range option, wrf, mpas_atm have examples.

      start_ind = get_index_start(idom, ivar)
      end_ind   = get_index_end(  idom, ivar)

      do j = 1, state_ens_handle%my_num_vars

         if (state_ens_handle%my_vars(j) >= start_ind .and. &
             state_ens_handle%my_vars(j) <= end_ind   ) then 

            do copy = 1, ens_size
               if (trim(perturb_distribution) == 'lognormal') then 
                  rng = random_gaussian(random_seq, 0.0_r8, 1.0_r8)
                  state_ens_handle%copies(copy,j) = state_ens_handle%copies(copy,j)*exp(stddev*rng)
              else 
               !MOHA: if it's not lognormal, then the only 
               !current other option is Gaussian. We can add more later. 

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
                     call error_handler(E_ERR, routine, string1, &
                                source, revision, revdate, text2=string2)
                  endif
               endif
            enddo
         endif
      enddo
   enddo
enddo DOMAIN

end subroutine pert_model_copies


!------------------------------------------------------------------
!> Returns the number of items in the state vector as an integer. 

function get_model_size()

integer(i8) :: get_model_size

get_model_size = model_size

end function get_model_size


!-----------------------------------------------------------------------
!> The LSM restart files have "time".
!> We are always using the 'most recent' which is, by defn, the last one.
!> The time in the restart file USED to DIFFER from time at which the state was valid.
!> And therefore we used to shift the time one time step back
!> But in the recent development of the model, this has been fixed
!> Therefore, the time of the restart file is the valid time now.
!> Example:
!> the noah_timestep = 3600 seconds &
!> restart_frequency_hours = 1
!> LSM restart filename is RESTART.2004010102_DOMAIN1 has
!>     Times = '2004-01-01_02:00:00' ;
!> The data is valid @ 2004-01-01_02:00:00

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

io = nf90_open(filename, NF90_NOWRITE, ncid)
call nc_check(io,routine,'open',filename)

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
      call error_handler(E_ERR,routine, string1, source, revision, revdate)
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

io = nf90_close(ncid)
call nc_check(io, routine, 'close', filename)
 
read(datestring(ntimes),'(i4,5(1x,i2))') year, month, day, hour, minute, second

read_model_time = set_date(year, month, day, hour, minute, second)

if ( do_output() .and. debug > 0 ) write(*,*)'routine: Last time string is '//trim(datestring(ntimes))
if ( do_output() .and. debug > 0 ) call print_date(read_model_time,' valid time is ')

deallocate(datestring)
      
end function read_model_time
   
!------------------------------------------------------------------
!> Returns the smallest increment in time that the model is capable 
!> of advancing the state in a given implementation, or the shortest
!> time you want the model to advance between assimilations.
!> This interface is required for all applications.

function shortest_time_between_assimilations()

type(time_type) :: shortest_time_between_assimilations

shortest_time_between_assimilations = time_step

end function shortest_time_between_assimilations


!------------------------------------------------------------------
!> Given an integer index into the state vector structure, returns the
!> associated location. A second intent(out) optional argument 
!> returns the DART 'QUANTITY' (i.e. QTY_STREAM_FLOW) for this index.
!> This interface is required for computing distances.

subroutine get_state_meta_data(index_in, location, var_type)

integer(i8),         intent(in)  :: index_in
type(location_type), intent(out) :: location
integer,             intent(out), optional :: var_type

integer :: iloc, jloc, kloc, varid, domid

call get_model_variable_indices(index_in, iloc, jloc, kloc, varid, domid, var_type)

location = domain_info(domid)%location(iloc,jloc,kloc)

if (do_output() .and. debug > 99) then
   call write_location(0,location,charstring=string1)
   write(*,*)'gsmd index,i,j,k = ',index_in, iloc, jloc, kloc, trim(string1)
endif

end subroutine get_state_meta_data


!------------------------------------------------------------------
!> Given a state handle, a location, and a DART QUANTITY,
!> return the expected estimate of the model state quantity at that location.
!> This function returns the estimates for the entire ensemble in a single call.
!>
!> The return code for successful return must be 0, any positive number is 
!> an error, and the expected_obs array contains unreliable values.
!> Negative values are reserved for use by the DART framework.
!> Distinct positive values for different types of errors are used
!> to aid in diagnosing problems.

subroutine model_interpolate(state_handle, ens_size, location, obs_type, &
                             expected_obs, istatus)

type(ensemble_type), intent(in) :: state_handle
integer,             intent(in) :: ens_size
type(location_type), intent(in) :: location
integer,             intent(in) :: obs_type
real(r8),           intent(out) :: expected_obs(ens_size)
integer,            intent(out) :: istatus(ens_size)

character(len=*), parameter :: routine = 'model_interpolate'

real(r8) :: closest_dist, distance
integer  :: domid, varid, closest_i, closest_j, closest_k
integer  :: i,j,k
integer(i8) :: closest_index

expected_obs(:) = MISSING_R8
istatus(:)      = 1  ! presume failure

varid = -1 
domid = -1
! Use the obs_type to figure out which variable we need, and in what domain.
DOMLOOP : do idom = 1,domain_count
   varid = get_varid_from_kind(idom, obs_type)
   if (varid > 0) then
      domid = idom
      exit DOMLOOP
   endif
enddo DOMLOOP

if (varid < 0) then
   istatus = 2 ! variable is not present in DART state
   return
endif

! Now we know the domain and the variable ID
! but we want to know the closest location to the variable on its grid
!>@todo only do this for the smaller grids ... there has to be a better way
!> for the larger grids

!>@todo might also need to check if the closest is in the same/right watershed
!>      as the observation ... somehow

!>@todo If the parameter domain has both 1D and 2D parameters ... we have to revisit this.

closest_dist = huge(1.0_r8)

do k = 1,size(domain_info(domid)%location,3)
do j = 1,size(domain_info(domid)%location,2)
do i = 1,size(domain_info(domid)%location,1)
   distance = get_dist(location, domain_info(domid)%location(i,j,k), no_vert=.false.)
   if (distance < closest_dist) then
      closest_dist  = distance
      closest_i     = i
      closest_j     = j
      closest_k     = k
   endif
enddo
enddo
enddo

closest_index = get_dart_vector_index(closest_i, closest_j, closest_k, domid, varid)
expected_obs  = get_state(closest_index, state_handle)
istatus       = 0

end subroutine model_interpolate


!-----------------------------------------------------------------------
!> Given an observation and its TYPE (the 'base') and a list of observations
!> and their TYPEs, returns the subset close to the base, 
!> their indices, and their distances.

subroutine get_close_obs(gc, base_loc, base_type, locs, loc_qtys, loc_types, &
                         num_close, close_ind, dist, ens_handle)

type(get_close_type), intent(in)    :: gc           !< handle to a get_close structure
type(location_type),  intent(inout) :: base_loc     !< location of interest
type(location_type),  intent(inout) :: locs(:)      !< all possible locations
integer,              intent(in)    :: base_type    !< observation TYPE
integer,              intent(in)    :: loc_qtys(:)  !< QTYs for all locations
integer,              intent(in)    :: loc_types(:) !< TYPEs for all locations
integer,              intent(out)   :: num_close    !< how many are close
integer,              intent(out)   :: close_ind(:) !< indices in the the locs array
real(r8),             optional, intent(out) :: dist(:)      !< distances
type(ensemble_type),  optional, intent(in)  :: ens_handle

integer :: t_ind, k

! Initialize variables to missing status

num_close = 0
close_ind = -99
if (present(dist)) dist = huge(1.0_r8)  !something big and positive (far away)

! Get all the potentially close obs but no dist (optional argument dist(:)
! is not present) This way, we are decreasing the number of distance
! computations that will follow.  This is a horizontal-distance operation and
! we don't need to have the relevant vertical coordinate information yet
! (for obs).

! TJH When called with state; subset based on HUC, watershed, link ID, ... whatever ...

! this calls the 'traditional' get_close routine but does not calculate the distances
call loc_get_close_obs(gc, base_loc, base_type, locs, loc_qtys, loc_types, &
                         num_close, close_ind)

! Loop over potentially close subset of obs priors 

if ( .false. ) then ! this block demonstrates how to do something new
   if (present(dist)) then
      do k = 1, num_close

         t_ind = close_ind(k)
         ! if something, leave original HUGE value.  otherwise, compute real dist.
         if (loc_qtys(t_ind) /= SOIL_MOISTURE) then
            dist(k) = get_dist(base_loc, locs(t_ind), base_type, loc_qtys(t_ind))
         endif

      enddo
   endif
else  ! repeat default behavior
   if (present(dist)) then
      do k = 1, num_close
         t_ind   = close_ind(k)
         dist(k) = get_dist(base_loc, locs(t_ind), base_type, loc_qtys(t_ind))
      enddo
   endif
endif

end subroutine get_close_obs


!-----------------------------------------------------------------------
!> Given a location, return the subset of close states that are local to
!> the mpi task

subroutine get_close_state(gc, base_loc, base_type, locs, loc_qtys, loc_indx, &
                         num_close, close_ind, dist, ens_handle)

type(get_close_type), intent(in)    :: gc           !< handle to a get_close structure
type(location_type),  intent(inout) :: base_loc     !< location of interest
integer,              intent(in)    :: base_type    !< observation TYPE
type(location_type),  intent(inout) :: locs(:)      !< all possible locations
integer,              intent(in)    :: loc_qtys(:)  !< QTYs for all locations
integer(i8),          intent(in)    :: loc_indx(:)  !< index into DART state
integer,              intent(out)   :: num_close    !< how many are close
integer,              intent(out)   :: close_ind(:) !< indices in the the locs array
real(r8),             optional, intent(out) :: dist(:)      !< distances
type(ensemble_type),  optional, intent(in)  :: ens_handle


integer  :: nclose                         !< how many are close
integer(i8) :: close_indices(size(close_ind)) !< indices in the the locs array
real(r8) :: distances(size(close_ind))     !< distances

character(len=*),parameter :: routine = 'get_close_state'
integer :: t_ind, k, i
integer :: depth
logical :: in_list
integer(i8) :: full_index
real(r8) :: reach_length
integer :: connection_index

! Initialize variables to missing status
reach_length = 0.0_r8
nclose       = 0
depth        = 1

if (present(dist)) dist = huge(1.0_r8)  !something far away

if (base_type < 0) then ! we are working with identity observations
   full_index = abs(base_type)
else
   call error_handler(E_ERR,routine,'unexpected',source,revision,revdate)
endif

connection_index = full_to_connection(full_index)
if (connection_index < 0) then
   call error_handler(E_ERR,routine,'unable to relate',source,revision,revdate)
endif

call get_link_tree(connection_index, max_link_distance, depth, &
                   reach_length, nclose, close_indices, distances)

! get the downstream close bits and augment the list

call get_downstream_links(connection_index, max_link_distance, depth, &
                   nclose, close_indices, distances)

! Nancy: which indices are mine?

call get_my_close(nclose, close_indices, distances, loc_indx, &
            num_close, close_ind, dist)

if (debug > 99) then
   write(string1,'("PE ",I7)') my_task_id()
   write(*,*)trim(string1),' get_close_state:num_close ',num_close
   write(*,*)trim(string1),' get_close_state:close_ind ',close_ind(1:num_close)
   write(*,*)trim(string1),' get_close_state:dist      ',dist(1:num_close)

   do k = 1,num_close
      write(*,*)trim(string1),' get_close_state: state_vector_index ',k, &
                ' is ',loc_indx(close_ind(k))
   enddo
endif

! If distances are not needed, we are done.
if (.not. present(dist)) return

!do k = 1, num_close
!
!      t_ind = close_ind(k)
!
!      ! Specialized treatment of the global parameters.
!      ! By definition, they are close to everything.
!      !MOHA: Multipliers are made spatially variable, keeping this block for future use.
!      !if (loc_qtys(t_ind) == QTY_BUCKET_MULTIPLIER .or. &
!      !    loc_qtys(t_ind) == QTY_RUNOFF_MULTIPLIER) then
!      !   dist(k) = 0.0_r8
!      !else
!         ! Compute real distance
!
!         dist(k) = get_dist(base_loc, locs(t_ind), base_type, loc_qtys(t_ind))
!
!      !endif
!      
!enddo

end subroutine get_close_state


!------------------------------------------------------------------
!> Write attributes and metadata to the DART output and diagnostic files

subroutine nc_write_model_atts(ncid, domain_id)

integer, intent(in) :: ncid      ! netCDF file identifier
integer, intent(in) :: domain_id

call nc_begin_define_mode(ncid) ! put file into define mode.

call nc_add_global_creation_time(ncid)
call nc_add_global_attribute(ncid, "model_source",   source )
call nc_add_global_attribute(ncid, "model_revision", revision )
call nc_add_global_attribute(ncid, "model_revdate",  revdate )
call nc_add_global_attribute(ncid, "model",          "wrfHydro")

!>@todo add the metadata for each domain ... lat/lon/alt etc.

call nc_end_define_mode(ncid)

call nc_synchronize_file(ncid) ! Flush the buffer and leave netCDF file open

end subroutine nc_write_model_atts


!------------------------------------------------------------------
!> Shutdown and clean-up

subroutine end_model()

do idom = 1,domain_count
   deallocate(domain_info(idom)%location)
enddo
deallocate(domain_info)

end subroutine end_model

!=======================================================================
! End of the required interfaces
!=======================================================================

!-----------------------------------------------------------------------
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
      call error_handler(E_ERR,routine, string1, &
                 source, revision, revdate, text2=string2)
   endif

   ! The internal DART routines check if the variable name is valid.

   ! Make sure DART kind is valid
   quantity = get_index_for_quantity(dartstr)
   if( quantity < 0 ) then
      write(string1,'(''there is no obs_kind "'',a,''" in obs_kind_mod.f90'')') &
                    trim(dartstr)
      call error_handler(E_ERR,routine,string1,source,revision,revdate)
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
   call error_handler(E_MSG,routine,string1,source,revision,revdate,text2=string2)
endif

end subroutine verify_variables


!-----------------------------------------------------------------------
!> Sets the location information arrays for each domain
!> Each location array is declared to be 3D to make it easy to use
!> the state_structure routines.

subroutine configure_domains()

character(len=*), parameter :: routine = 'configure_domains'

integer :: i,j,k, var_size, num_vars, ivar

do idom = 1,domain_count

   if (idom == idom_hydro) then 
      j = 1
      k = 1
      allocate( domain_info(idom)%location(n_link,j,k) )
      do i = 1,n_link
         domain_info(idom)%location(i,j,k) = set_location(linkLong(i),linkLat(i),linkAlt(i),VERTISHEIGHT)
      enddo

   elseif (idom == idom_lsm) then
      call error_handler(E_ERR,routine,'LSM locations unsupported for now',source,revision,revdate)

   !MOHA: checking the shape of the parameters to decide on the location array.
   elseif (idom == idom_parameters) then
      num_vars = get_num_variables(idom)
      var_size = 0
      do ivar = 1, num_vars
         var_size = var_size + get_variable_size(idom,ivar)
      enddo
      if (var_size /= num_vars) then ! spatially varying parameters case
         j = 1
         k = 1
         allocate( domain_info(idom)%location(n_link,j,k) )
         do i = 1,n_link
            domain_info(idom)%location(i,j,k) = set_location(linkLong(i),linkLat(i),linkAlt(i),VERTISHEIGHT)
         enddo              
      else ! global parameters case  
         i = 1
         j = 1
         k = 1
         allocate( domain_info(idom)%location(i,j,k) )
         domain_info(idom)%location(i,j,k) = set_location(35.38_r8, 33.60_r8, linkAlt(i), VERTISHEIGHT)
      endif
   endif
enddo

end subroutine configure_domains


!-----------------------------------------------------------------------
!> determine which indices are on this task, return the ...

subroutine get_my_close(nclose, close_indices, distances, state_indices, &
            num_close, close_ind, dist)

integer,          intent(in)  :: nclose
integer(i8),      intent(in)  :: close_indices(:)
real(r8),         intent(in)  :: distances(:)
integer(i8),      intent(in)  :: state_indices(:)
integer,          intent(out) :: num_close
integer,          intent(out) :: close_ind(:)
real(r8),         intent(out) :: dist(:)

integer :: imine, iclose

num_close = 0

ALLPOSSIBLE : do imine = 1,size(state_indices)
   CLOSE : do iclose = 1,nclose

      if ( close_indices(iclose) == state_indices(imine) ) then
          num_close            = num_close + 1
          close_ind(num_close) = imine
          dist(num_close)      = distances(iclose)
      endif

   enddo CLOSE
enddo ALLPOSSIBLE

end subroutine get_my_close


!===================================================================
! End of model_mod
!===================================================================

end module model_mod

! <next few lines under version control, do not edit>
! $URL$
! $Id$
! $Revision$
! $Date$
