! DART software - Copyright UCAR. This open source software is provided
! by UCAR, "as is", without charge, subject to all terms of use at
! http://www.image.ucar.edu/DAReS/DART/DART_download
!

module model_mod

! This is the interface between pywatershed and DART

use types_mod,             only : r4, r8, i8, i4, MISSING_R8, vtablenamelength, &
                                  earth_radius, MISSING_I, MISSING_I8

use time_manager_mod,      only : time_type, set_time, set_calendar_type

use location_mod,          only : location_type, set_location, get_location,  &
                                  get_close_obs, get_close_state,             &
                                  convert_vertical_obs, convert_vertical_state

use utilities_mod,         only : do_nml_file, do_nml_term, E_ERR, E_MSG,        &
                                  nmlfileunit, find_namelist_in_file, to_upper,  &
                                  check_namelist_read, file_exist, error_handler

use location_io_mod,      only : nc_write_location_atts, nc_write_location

use mpi_utilities_mod,    only : my_task_id 

use netcdf_utilities_mod, only : nc_add_global_attribute, nc_synchronize_file,      &
                                 nc_add_global_creation_time, nc_begin_define_mode, &
                                 nc_end_define_mode, nc_open_file_readonly,         &
                                 nc_close_file, nc_get_global_attribute,            &
                                 nc_get_dimension_size, nc_get_variable, nc_check 

use         obs_kind_mod,  only : get_index_for_quantity, &
                                  get_name_for_quantity, &
                                  get_quantity_for_type_of_obs, &
                                  QTY_BUCKET_MULTIPLIER, &
                                  QTY_RUNOFF_MULTIPLIER, &
                                  QTY_STREAM_FLOW

use ensemble_manager_mod,  only : ensemble_type

use distributed_state_mod, only : get_state

use state_structure_mod,   only : add_domain, get_domain_size, get_index_start,      &
                                  get_index_end, get_num_domains, get_num_variables, &
                                  get_num_dims, get_dim_name, get_dim_length,        &
                                  get_variable_name, get_model_variable_indices,     &
                                  get_varid_from_kind, get_dart_vector_index,        &
                                  get_variable_size, state_structure_info

use default_model_mod,     only : adv_1step, end_model, pert_model_copies,              & 
                                  nc_write_model_vars, MAX_STATE_VARIABLE_FIELDS_CLAMP, &
                                  init_time => fail_init_time,                          & 
                                  init_conditions => fail_init_conditions,              &
                                  parse_variables_clamp 

use dart_time_io_mod,      only : read_model_time, write_model_time

use netcdf

implicit none
private

! required by DART code - will be called from filter and other
! DART executables. 
public :: get_model_size,         &
          get_state_meta_data,    &
          model_interpolate,      &
          static_init_model,      &
          init_conditions,        &
          init_time,              &
          nc_write_model_atts,    & 
          pert_model_copies,      &
          nc_write_model_vars,    &
          get_close_obs,          &
          get_close_state,        &
          end_model,              &
          adv_1step,              &
          convert_vertical_obs,   &
          convert_vertical_state, &
          read_model_time,        &
          write_model_time,       &
          shortest_time_between_assimilations

character(len=256), parameter :: source = "pywatershed/model_mod.f90"

integer, parameter :: IDSTRLEN = 15 ! Number of character for a USGS gauge ID

!TODO: Revisit for bigger (e.g., CONUS) domain
!For DRB the maximum number of contributing HRUs was 3
integer, parameter :: NUM_HRU  = 10 ! Number of contributing HRUs for each segment  

type(location_type), allocatable :: state_loc(:)  ! state locations, compute once and store for speed

! user-defined type to enable a 'linked list' of stream links
type link_relations
   private
   character(len=IDSTRLEN) :: gaugeName          = ''         ! USGS gauge ID
   integer                 :: segID              = MISSING_I  ! National Hydrologic Model (NHM) ID
   real(r4)                :: segLength          = 0.0_r8     ! seg length (meters)
   real(r4)                :: segSlope           = 0.0_r8     ! seg slope (decimal fraction)
   real(r4)                :: segWidth           = 0.0_r8     ! seg width (meters)
   real(r4)                :: segElevation       = 0.0_r8     ! seg elevation (meters)
   real(r4)                :: segManning         = 0.0_r8     ! seg manning's n in seconds / meter ** (1/3)
   integer                 :: segType            = MISSING_I  ! seg type(0: seg, 1: headwater, 2: lake, ...)
   integer(i8)             :: domain_offset      = MISSING_I8 ! into DART state vector
   integer                 :: downstream_segID   = MISSING_I
   integer                 :: downstream_index   = MISSING_I  ! into link_type structure
   integer, allocatable    :: upstream_segID(:)
   integer, allocatable    :: upstream_index(:)               ! into link_type structure
   integer, allocatable    :: avail_hru(:)                    ! into 
end type link_relations

type(link_relations), allocatable  :: connections(:)          ! Connections structure of the entire river network
integer, allocatable, dimension(:) :: num_up_links            ! Number of links upstream from each segment

integer, allocatable, dimension(:) :: GaugeID
integer, allocatable, dimension(:) :: segNHMid
integer, allocatable, dimension(:) :: segGauge
integer, allocatable, dimension(:) :: segTyp
integer, allocatable, dimension(:) :: hruMask

real(r8), allocatable, dimension(:) :: segLat
real(r8), allocatable, dimension(:) :: segLon
real(r8), allocatable, dimension(:) :: segLen
real(r8), allocatable, dimension(:) :: segSlp
real(r8), allocatable, dimension(:) :: segElv
real(r8), allocatable, dimension(:) :: segWid
real(r8), allocatable, dimension(:) :: segMan

integer            :: nseg, ngages, nhru                                      

type(time_type)    :: time_step
integer(i8)        :: model_size
logical, save      :: module_initialized = .false. 
character(len=256) :: msg1, msg2, msg3              ! Strings for error and warning messages
integer            :: idom                          ! domain counter
integer            :: idom_chn = -1                 ! channel domain number placeholder
integer            :: idom_hru = -1                 ! hru domain number placeholder
integer            :: idom_prm = -1                 ! parameters domain number placeholder

! Model namelist declarations with defaults
integer            :: assimilation_period_days     = -1           ! model instantaneous streamflow at the hour 23 of the day
integer            :: assimilation_period_seconds  = -1           ! NA
integer            :: debug                        = 0            ! debugging flag
real(r8)           :: model_perturbation_amplitude = 0.1          ! perturb parameter for initial ensemble generation
real(r8)           :: max_link_distance            = 10000.0      ! Max distance along the stream: 10 km
character(len=256) :: perturb_distribution         = 'lognormal'  ! distribution needed for initial ensemble generation
character(len=256) :: domain_order(3)              = ''           ! Domains: channel, hru, params??
character(len=256) :: domain_shapefiles(3)         = ''           ! template restart files for each domain
character(len=256) :: channel_config_file          = 'dis_seg.nc' ! file with stream connectivity information
character(len=256) :: hru_config_file              = 'PRMS.nc'    ! file relating segments to HRUs

character(len=vtablenamelength) :: channel_variables(MAX_STATE_VARIABLE_FIELDS_CLAMP) = '' ! channel state variables
character(len=vtablenamelength) :: hru_variables(MAX_STATE_VARIABLE_FIELDS_CLAMP)     = '' ! hru (hydrologic response unit) variables
character(len=vtablenamelength) :: parameters(MAX_STATE_VARIABLE_FIELDS_CLAMP)        = '' ! parameters

namelist /model_nml/ assimilation_period_days,     &
                     assimilation_period_seconds,  &
                     model_perturbation_amplitude, &
                     perturb_distribution,         &
                     max_link_distance,            &
                     domain_order,                 &
                     domain_shapefiles,            &
                     debug,                        &
                     channel_config_file,          &
                     hru_config_file,              &
                     channel_variables,            &
                     hru_variables,                & 
                     parameters

contains


!---------------------------------------------------------------------
! Initialize the model, read namelist, configure domains and variables

subroutine static_init_model()

integer  :: iunit, io
integer  :: domain_count

! only do this once
if (module_initialized) return 
module_initialized = .true. 

! Read the namelist 
call find_namelist_in_file("input.nml", "model_nml", iunit)
read(iunit, nml = model_nml, iostat = io)
call check_namelist_read(iunit, io, "model_nml")

! Output the namelist values if requested
if (do_nml_file()) write(nmlfileunit, nml=model_nml)
if (do_nml_term()) write(     *     , nml=model_nml)

! Time handling 
call set_calendar_type('gregorian')
time_step  = set_time(assimilation_period_seconds, assimilation_period_days)

! Add, configure and manage domains
call manage_domains()

! Determine the model size
domain_count = get_num_domains()
model_size = 0
do idom = 1,domain_count
   model_size = model_size + get_domain_size(idom)
enddo

end subroutine static_init_model


!----------------------------------------------------
! Add hydro + other domains and configure the network 
! Make sure we have consistent templates and domains

subroutine manage_domains()

character(len=*), parameter :: routine = 'manage_domains'

integer            :: adom, ndom 
character(len=256) :: domain_name

! Make sure there are enough domains 
! If no domains are entered, then error out
if (len_trim(domain_order(1)) == 0) then 
   write(msg1, *) 'The domain order is incorrectly specified.'
   write(msg2, *) 'The domain order cannot be empty, it should be something like:'
   write(msg3, *) 'domain_order = "channel", "hru", "parameters"'
   call error_handler(E_ERR, routine, msg1, source, text2=msg2, text3=msg3)
endif

! With a correct domain order, now look
! for the corresponding template files and 
! start processing each domain
ndom = size(domain_order)
adom = 0 ! We should have at least 1 active domain

do idom = 1, ndom
   ! Figure how many active domains are available 
   if (len_trim(domain_order(idom)) == 0) exit 
   adom = adom + 1
enddo
if (debug > 99) write(*, '(A, i3)') 'The number of active domains is:', adom

! Check to see if the template files
! for these active domains are correctly specified
do idom = 1, adom
   if (len_trim(domain_shapefiles(idom)) == 0) then 
      write(msg1, '(A, i3, 1x, A)') 'The template file for domain', & 
                  idom, 'is missing. Exiting ..' 
      call error_handler(E_ERR, routine, msg1, source)
   endif
   if (debug > 99) then 
      write(*, '(3A, i3)') 'Found shapefile "', &
               trim(domain_shapefiles(idom)), &
               '" for domain', idom
   endif
enddo

! Read and add variables in the active domains
DOMAINS: do idom = 1, adom
   ! Domain name and its template
   domain_name = domain_order(idom)
   call to_upper(domain_name)

   if (index(domain_name, 'CHANNEL') > 0) then 
      ! Read the stream network and channel attributes   
      call read_stream_network(channel_config_file)
      call read_channel_atts(domain_shapefiles(idom)) !TODO: may need to remove

      ! Add channel variables: streamflow
      idom_chn = add_domain(domain_shapefiles(idom), &
                 parse_variables_clamp(channel_variables)) 
      call verify_domain(idom_chn, nseg, domain_shapefiles(idom), &
                         channel_config_file)
  
   elseif (index(domain_name, 'HRU') > 0) then 
       ! Read the hru properties
       call read_hru_properties(hru_config_file)
 
       ! Add hru variables: groundwater
       idom_hru = add_domain(domain_shapefiles(idom), &
                  parse_variables_clamp(hru_variables))  
       call verify_domain(idom_hru, nhru, domain_shapefiles(idom), &
                          hru_config_file) 

   elseif (index(domain_name, 'PARAMETERS') > 0) then 
       ! To be added 
   endif

    
enddo DOMAINS 

end subroutine manage_domains


!--------------------------------------------------------------
! Read stream network variables and construct the connectivity 

subroutine read_stream_network(filename)

character(len=*), intent(in) :: filename
character(len=*), parameter  :: routine = 'read_stream_network'

integer              :: ncid, nindex, io, varid
integer, allocatable :: fromIndices(:)
integer, allocatable :: fromIndsStart(:)
integer, allocatable :: fromIndsEnd(:)
integer, allocatable :: toIndex(:)

ncid = nc_open_file_readonly(filename, routine)        

! Read in the dimensions
nseg   = nc_get_dimension_size(ncid, 'nsegment' , routine)
nindex = nc_get_dimension_size(ncid, 'index'    , routine)
ngages = nc_get_dimension_size(ncid, 'npoigages', routine)

allocate(segLat(nseg), segLon(nseg), &
         segSlp(nseg), segLen(nseg), &
         segWid(nseg), segElv(nseg), &
         segTyp(nseg), segNHMid(nseg))

! Read segment properties
call nc_get_variable(ncid, 'nhm_seg'     , segNHMid, routine)
call nc_get_variable(ncid, 'seg_lat'     , segLat  , routine)
call nc_get_variable(ncid, 'seg_length'  , segLen  , routine)
call nc_get_variable(ncid, 'seg_elev'    , segElv  , routine)
call nc_get_variable(ncid, 'seg_slope'   , segSlp  , routine)
call nc_get_variable(ncid, 'seg_width'   , segWid  , routine)
call nc_get_variable(ncid, 'segment_type', segTyp  , routine)

allocate(GaugeID(ngages), segGauge(ngages))

! Read gauge-related variables
call nc_get_variable(ncid, 'poi_gauges'      , GaugeID , routine)
call nc_get_variable(ncid, 'poi_gage_segment', segGauge, routine)

allocate(fromIndices(nindex))
allocate(fromIndsStart(nseg))
allocate(fromIndsEnd(nseg))
allocate(toIndex(nseg))
allocate(num_up_links(nseg))

! Read in temporary connection arrays
call nc_get_variable(ncid, 'fromIndices'  , fromIndices  , routine)
call nc_get_variable(ncid, 'fromIndsStart', fromIndsStart, routine)
call nc_get_variable(ncid, 'fromIndsEnd'  , fromIndsEnd  , routine)
call nc_get_variable(ncid, 'tosegment'    , toIndex      , routine)
call nc_get_variable(ncid, 'num_up_links' , num_up_links , routine)

! Fill the connections
call fill_connections(toIndex, fromIndices, fromIndsStart, fromIndsEnd)

deallocate(fromIndices, fromIndsStart, &
           fromIndsEnd, toIndex)

end subroutine read_stream_network


!-------------------------------------------------------------
! Populate the connections using the segment properties and 
! the offline network calculations

subroutine fill_connections(toIndex,fromIndices,fromIndsStart,fromIndsEnd)

integer, intent(in) :: toIndex(:)
integer, intent(in) :: fromIndices(:)
integer, intent(in) :: fromIndsStart(:)
integer, intent(in) :: fromIndsEnd(:)

integer             :: iseg

allocate(connections(nseg))

! segment properties
do iseg = 1, nseg

   ! naming
   connections(iseg)%segID         = segNHMid(iseg)
   connections(iseg)%gaugeName     = ''  

   ! geometry
   connections(iseg)%segLength     = segLen(iseg)
   connections(iseg)%segWidth      = segWid(iseg)
   connections(iseg)%segSlope      = segSlp(iseg)
   connections(iseg)%segElevation  = segElv(iseg)
   connections(iseg)%segType       = segTyp(iseg)

   ! dart-related
   connections(iseg)%domain_offset = iseg 
   
   ! downstream
   if (toIndex(iseg) == 0) then 
      connections(iseg)%downstream_segID = MISSING_I
      connections(iseg)%downstream_index = MISSING_I
   else
      connections(iseg)%downstream_segID = segNHMid(toIndex(iseg))
      connections(iseg)%downstream_index = toIndex(iseg) 
   endif

   ! upstream
   if (num_up_links(iseg) == 0) then 
      allocate(connections(iseg)%upstream_segID(1))
      allocate(connections(iseg)%upstream_index(1))

      connections(iseg)%upstream_segID(1) = MISSING_I 
      connections(iseg)%upstream_index(1) = MISSING_I
   else
      allocate(connections(iseg)%upstream_segID(num_up_links(iseg)))
      allocate(connections(iseg)%upstream_index(num_up_links(iseg)))

      connections(iseg)%upstream_segID(:) = segNHMid(fromIndices(fromIndsStart(iseg):fromIndsEnd(iseg)))
      connections(iseg)%upstream_index(:) = fromIndices(fromIndsStart(iseg):fromIndsEnd(iseg))
   endif
 
enddo

! Update USGS gauge IDs
do iseg = 1, ngages
   write(connections(segGauge(iseg))%gaugeName, '(i0)') GaugeID(iseg)
enddo

if (debug > 99) then 
   write(msg1, '("PE",i2)') my_task_id()
   do iseg = 1, nseg
      write(*, *) ''
      write(*, *) trim(msg1), ' Connectivity for segment: ', iseg
      write(*, *) trim(msg1), ' Segment NHM ID          : ', connections(iseg)%segID
      write(*, *) trim(msg1), ' Segment USGS Gauge      : ', connections(iseg)%gaugeName
      write(*, *) trim(msg1), ' Segment Length          : ', connections(iseg)%segLength
      write(*, *) trim(msg1), ' Segment Width           : ', connections(iseg)%segWidth
      write(*, *) trim(msg1), ' Segment Slope           : ', connections(iseg)%segSlope
      write(*, *) trim(msg1), ' Segment Elevation       : ', connections(iseg)%segElevation
      write(*, *) trim(msg1), ' Segment Type            : ', connections(iseg)%segType
      write(*, *) trim(msg1), ' Segment down segID      : ', connections(iseg)%downstream_segID
      write(*, *) trim(msg1), ' Segment down index      : ', connections(iseg)%downstream_index
      write(*, *) trim(msg1), ' Segment up segID        : ', connections(iseg)%upstream_segID
      write(*, *) trim(msg1), ' Segment up index        : ', connections(iseg)%upstream_index
   enddo
endif

end subroutine fill_connections


!-----------------------------------------------------------
! Make sure the domain is read properly and that shapefile
! and restart files are consistent

subroutine verify_domain(dom_id, nvars, rest_file, config_file)

character(len=*), parameter  :: routine = 'verify_domain'

integer, intent(in)          :: dom_id
integer, intent(in)          :: nvars
character(len=*), intent(in) :: rest_file, config_file

integer :: var_size

if (debug > 99) call state_structure_info(dom_id)

var_size = get_variable_size(dom_id, 1) 
if (var_size /= nvars) then 
   write(msg1, *) 'Restart file and shapefile are not consustent.'
   write(msg2, *) 'Number of variables in the restart "', & 
                  trim(rest_file), '" is ', var_size
   write(msg3, *) 'Number of variables in the shapefile "', &
                  trim(config_file), '" is ', nvars
   call error_handler(E_ERR, routine, msg1, source, text2=msg2, text3=msg3)
endif

end subroutine verify_domain


!--------------------------------------------------------------
! Read hru parameters and add to the segment the connectivity 

subroutine read_hru_properties(filename)

character(len=*), intent(in) :: filename
character(len=*), parameter  :: routine = 'read_hru_properties'

integer              :: ncid, iseg, ihru, count_hrus
integer, allocatable :: hru_seg(:), hru_inds(:)

ncid = nc_open_file_readonly(filename, routine)

! Read in the HRU dimension
nhru = nc_get_dimension_size(ncid, 'nhru' , routine)

allocate(hru_seg(nhru), segMan(nseg), hruMask(nseg))

! Read hru/segment relation
call nc_get_variable(ncid, 'hru_segment', hru_seg, routine)
call nc_get_variable(ncid, 'mann_n'     , segMan , routine)

! Additional segment properties
connections(:)%segManning = segMan

! Construct some sort of a bucket mask
! HRU index contributes (thru groundwater) to segment .. 
allocate(hru_inds(NUM_HRU))
do iseg = 1, nseg 

   ! Find which HRUs contribute to segement iseg
   hru_inds   = pack([(ihru, ihru = 1, nhru)], hru_seg == iseg)
   count_hrus = size(hru_inds)
   
   ! Fill connections with HRU information
   if (count_hrus == 0) then
      allocate(connections(iseg)%avail_hru(1))
      connections(iseg)%avail_hru(1) = MISSING_I
   else
      allocate(connections(iseg)%avail_hru(count_hrus))
      connections(iseg)%avail_hru(:) = hru_inds(1:count_hrus)
   endif

   if (debug > 99) then   
       write(msg1, '("PE",i2)') my_task_id()
       write(*, *) trim(msg1), ' segment: ', iseg, 'Contributing HRUs: ', connections(iseg)%avail_hru(:)
   endif  
enddo 

end subroutine read_hru_properties


!-------------------------------------------------------------
! Read channel global attributes from streamflow restart file

subroutine read_channel_atts(filename)

character(len=*), intent(in) :: filename

end subroutine read_channel_atts


!------------------------------------------------------------------
! Returns the number of items in the state vector as an integer. 
! This interface is required for all applications.

function get_model_size()

integer(i8) :: get_model_size

get_model_size = model_size

end function get_model_size


!------------------------------------------------------------------
! Returns the smallest increment in time that the model is capable 
! of advancing the state in a given implementation, or the shortest
! time you want the model to advance between assimilations.
! This interface is required for all applications.

function shortest_time_between_assimilations()

type(time_type) :: shortest_time_between_assimilations

shortest_time_between_assimilations = time_step

end function shortest_time_between_assimilations


!------------------------------------------------------------------
! Given a state handle, a location, and a model state variable quantity,
! interpolates the state variable fields to that location and returns
! the values in expected_obs. The istatus variables should be returned as
! 0 unless there is some problem in computing the interpolation in
! which case an alternate value should be returned. The itype variable
! is an integer that specifies the quantity of field (for
! instance temperature, zonal wind component, etc.). In low order
! models that have no notion of types of variables this argument can
! be ignored. For applications in which only perfect model experiments
! with identity observations (i.e. only the value of a particular
! state variable is observed), this can be a NULL INTERFACE.

subroutine model_interpolate(state_handle, ens_size, location, iqty, expected_obs, istatus)

type(ensemble_type),  intent(in) :: state_handle
integer,              intent(in) :: ens_size
type(location_type),  intent(in) :: location
integer,              intent(in) :: iqty
real(r8),            intent(out) :: expected_obs(ens_size) !< array of interpolated values
integer,             intent(out) :: istatus(ens_size)

! This should be the result of the interpolation of a
! given quantity (iqty) of variable at the given location.
expected_obs(:) = MISSING_R8

! The return code for successful return should be 0. 
! Any positive number is an error.
! Negative values are reserved for use by the DART framework.
! Using distinct positive values for different types of errors can be
! useful in diagnosing problems.
istatus(:) = 1

end subroutine model_interpolate



!------------------------------------------------------------------
! Given an integer index into the state vector structure, returns the
! associated location. A second intent(out) optional argument quantity
! (qty) can be returned if the model has more than one type of field
! (for instance temperature and zonal wind component). This interface is
! required for all filter applications as it is required for computing
! the distance between observations and state variables.

subroutine get_state_meta_data(index_in, location, qty_type)

integer(i8),         intent(in)  :: index_in
type(location_type), intent(out) :: location
integer,             intent(out), optional :: qty_type

! these should be set to the actual location and state quantity
location = state_loc(index_in)

end subroutine get_state_meta_data


!------------------------------------------------------------------
! Writes model-specific attributes to a netCDF file

subroutine nc_write_model_atts(ncid, domain_id)

integer, intent(in) :: ncid
integer, intent(in) :: domain_id

! put file into define mode.

integer :: msize

msize = int(model_size, i4)

call nc_begin_define_mode(ncid)

call nc_add_global_creation_time(ncid)

call nc_add_global_attribute(ncid, "model_source", source )

call nc_add_global_attribute(ncid, "model", "template")

call nc_write_location_atts(ncid, msize)
call nc_end_define_mode(ncid)
call nc_write_location(ncid, state_loc, msize)

! Flush the buffer and leave netCDF file open
call nc_synchronize_file(ncid)

end subroutine nc_write_model_atts

!===================================================================
! End of model_mod
!===================================================================
end module model_mod

