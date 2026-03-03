! DART software - Copyright UCAR. This open source software is provided
! by UCAR, "as is", without charge, subject to all terms of use at
! http://www.image.ucar.edu/DAReS/DART/DART_download
!

module model_mod

! This is the interface between pywatershed and DART

use types_mod,             only : r4, r8, i8, i4, MISSING_R8, vtablenamelength,      &
                                  earth_radius, MISSING_I, MISSING_I8
use time_manager_mod,      only : time_type, set_time, set_calendar_type,            &
                                  set_date, print_date, increment_time 
use utilities_mod,         only : do_nml_file, do_nml_term, E_ERR, E_MSG,            &
                                  nmlfileunit, find_namelist_in_file, to_upper,      &
                                  check_namelist_read, error_handler
use         location_mod,  only : location_type, get_close_type,                     &
                                  loc_get_close_obs => get_close_obs,                &
                                  loc_get_close_state => get_close_state,            &
                                  convert_vertical_obs, convert_vertical_state,      &
                                  set_location, VERTISHEIGHT,                        &
                                  write_location
use mpi_utilities_mod,     only : my_task_id 
use netcdf_utilities_mod,  only : nc_add_global_attribute, nc_synchronize_file,      &
                                  nc_add_global_creation_time, nc_begin_define_mode, &
                                  nc_end_define_mode, nc_open_file_readonly,         &
                                  nc_close_file,                                     &
                                  nc_get_dimension_size, nc_get_variable, nc_check,  &
                                  nc_get_attribute_from_variable, nc_synchronize_file 
use         obs_kind_mod,  only : QTY_STREAM_FLOW
use ensemble_manager_mod,  only : ensemble_type
use state_structure_mod,   only : add_domain, get_domain_size, get_index_start,      &
                                  get_index_end, get_num_domains, get_num_variables, &
                                  get_variable_name, get_model_variable_indices,     &
                                  get_variable_size, state_structure_info
use default_model_mod,     only : adv_1step, end_model, nc_write_model_vars,         & 
                                  MAX_STATE_VARIABLE_FIELDS_CLAMP,                   &
                                  init_time => fail_init_time,                       &  
                                  init_conditions => fail_init_conditions,           &
                                  parse_variables_clamp 
use dart_time_io_mod,      only : write_model_time
use random_seq_mod,        only : random_seq_type, init_random_seq,                  &
                                  random_gaussian
use netcdf

implicit none
private

! required by DART code - will be called from filter and other
! DART executables. 
public :: get_model_size,         &
          get_number_of_segments, &
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
integer, parameter :: JINDEX   = 1  ! 1D river network
integer, parameter :: KINDEX   = 1  ! 1D river network

! TODO: Revisit for bigger (e.g., CONUS) domain
! For DRB the maximum number of contributing HRUs 
! (hydrologic response unit) was 3
integer, parameter :: NUM_HRU      = 10 ! Number of contributing HRUs for each segment  
integer, parameter :: INACTIVE_HRU = 0  ! 0=inactive; 1=land; 2=lake; 3=swale
integer, parameter :: MASKED_HRU   = 0

type domain_locations
   private
   type(location_type), allocatable :: location(:,:,:)
end type domain_locations
type(domain_locations), allocatable :: domain_info(:)

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
   integer, allocatable    :: avail_hru(:)                    ! into hru_type structure 
end type link_relations

type(link_relations), allocatable  :: connections(:)          ! Connections structure of the entire river network
integer, allocatable, dimension(:) :: num_up_links            ! Number of links upstream from each segment

character(len=IDSTRLEN), allocatable, dimension(:) :: GaugeID
integer, allocatable, dimension(:) :: segNHMid
integer, allocatable, dimension(:) :: segGauge  
integer, allocatable, dimension(:) :: segTyp     ! Segment type; could be headwater, lake, ...
integer, allocatable, dimension(:) :: hruMask    ! Construct it as either 0 (no HRUs) or non-zero (# of contributing HRUs) 
integer, allocatable, dimension(:) :: hruTyp     ! HRU type (could be inactive)

real(r8), allocatable, dimension(:) :: segLat
real(r8), allocatable, dimension(:) :: segLon
real(r8), allocatable, dimension(:) :: segLen
real(r8), allocatable, dimension(:) :: segSlp
real(r8), allocatable, dimension(:) :: segElv
real(r8), allocatable, dimension(:) :: segWid
real(r8), allocatable, dimension(:) :: segMan

real(r8), allocatable, dimension(:) :: hruLat
real(r8), allocatable, dimension(:) :: hruLon
real(r8), allocatable, dimension(:) :: hruElv

integer            :: nseg, ngages, nhru                                      
integer            :: domain_count

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
character(len=256) :: hru_config_file              = 'dis_hru.nc' ! file with hru information
character(len=256) :: hydro_config_file            = 'PRMS.nc'    ! file relating segments to HRUs

character(len=vtablenamelength) :: channel_variables(MAX_STATE_VARIABLE_FIELDS_CLAMP) = '' ! channel state variables
character(len=vtablenamelength) :: hru_variables(MAX_STATE_VARIABLE_FIELDS_CLAMP)     = '' ! hru variables
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
                     hydro_config_file,            &
                     channel_variables,            &
                     hru_variables,                & 
                     parameters

contains


!---------------------------------------------------------------------
! Initialize the model, read namelist, configure domains and variables

subroutine static_init_model()

integer  :: iunit, io

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

! Locations for the domains
allocate(domain_info(domain_count))
call configure_domloc()

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

   !TODO: address variable units in the restart files 
   if (index(domain_name, 'CHANNEL') > 0) then 
      ! Read the stream network and channel attributes   
      call read_stream_network(channel_config_file)

      ! Add channel variables: streamflow (cfs)
      idom_chn = add_domain(domain_shapefiles(idom), &
                 parse_variables_clamp(channel_variables)) 
      call verify_domain(idom_chn, nseg, domain_shapefiles(idom), &
                         channel_config_file)
  
   elseif (index(domain_name, 'HRU') > 0) then 
       ! Read the hru properties
       call read_hru_properties(hru_config_file, hydro_config_file)
 
       ! Add hru variables: groundwater (inches)
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

integer              :: ncid, nindex, strLen, igage
integer, allocatable :: fromIndices(:)
integer, allocatable :: fromIndsStart(:)
integer, allocatable :: fromIndsEnd(:)
integer, allocatable :: toIndex(:)

logical, save :: network_read = .false. 

! Only do this once
if (network_read) return 
network_read = .true.

ncid = nc_open_file_readonly(filename, routine)        

! Read in the dimensions
nseg   = nc_get_dimension_size(ncid, 'nsegment' , routine)
nindex = nc_get_dimension_size(ncid, 'index'    , routine)
ngages = nc_get_dimension_size(ncid, 'npoigages', routine)
strLen = nc_get_dimension_size(ncid, 'IDLength' , routine)

allocate(segLat(nseg), segLon(nseg), &
         segSlp(nseg), segLen(nseg), &
         segWid(nseg), segElv(nseg), &
         segTyp(nseg), segNHMid(nseg))

! Read segment properties
call nc_get_variable(ncid, 'nhm_seg'     , segNHMid, routine)
call nc_get_variable(ncid, 'seg_lats'    , segLat  , routine)
call nc_get_variable(ncid, 'seg_lons'    , segLon  , routine)
call nc_get_variable(ncid, 'seg_length'  , segLen  , routine)
call nc_get_variable(ncid, 'seg_elev'    , segElv  , routine)
call nc_get_variable(ncid, 'seg_slope'   , segSlp  , routine)
call nc_get_variable(ncid, 'seg_width'   , segWid  , routine)
call nc_get_variable(ncid, 'segment_type', segTyp  , routine)

! DART-style longitude for the segments
where(segLon <    0.0_r8) segLon = segLon + 360.0_r8
where(segLon == 360.0_r8) segLon = 0.0_r8

allocate(GaugeID(ngages), segGauge(ngages))

! Read gauge-related variables
call get_string_array(ncid, 'IDLength', strLen, 'poi_gauges', GaugeID)
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
call nc_close_file(ncid, routine)

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
   connections(segGauge(iseg))%gaugeName = GaugeID(iseg)
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

subroutine read_hru_properties(hru_file, prms_file)

character(len=*), intent(in) :: hru_file, prms_file
character(len=*), parameter  :: routine = 'read_hru_properties'

integer              :: ncid, iseg, ihru, count_hrus
integer, allocatable :: hru_seg(:), hru_inds(:)
character(len=256)   :: fmt
logical, save        :: basin_read = .false. 

! Only do this once
if (basin_read) return 
basin_read = .true.

ncid = nc_open_file_readonly(prms_file, routine)

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
      connections(iseg)%avail_hru(1) = MISSING_I       ! Similar to "Masked Buckets"
   else
      allocate(connections(iseg)%avail_hru(count_hrus))
      connections(iseg)%avail_hru(:) = hru_inds
   endif
   hruMask(iseg) = count_hrus  ! Number of contributing HRUS >= 0

   if (debug > 99) then   
       write(msg1, '("PE", i2)') my_task_id()
       write(fmt, '(A, i2, A)') '(2A, i5, X, A, i3, X, A,', max(count_hrus, 1), 'i10)' 
       write(*, fmt) trim(msg1), ' segment: ', iseg, & 
                                 'HRU Count:', hruMask(iseg), &
                                 'Contributing HRUs:', connections(iseg)%avail_hru(:)
   endif  
enddo 
deallocate(hru_seg, hru_inds)

call nc_close_file(ncid, routine)

! Now, get HRU data
ncid = nc_open_file_readonly(hru_file, routine)

allocate(hruLat(nhru), hruLon(nhru), &
         hruElv(nhru), hruTyp(nhru))

! Read hru location variables
call nc_get_variable(ncid, 'hru_lat' , hruLat, routine)
call nc_get_variable(ncid, 'hru_lon' , hruLon, routine)
call nc_get_variable(ncid, 'hru_elev', hruElv, routine)
call nc_get_variable(ncid, 'hru_type', hruTyp, routine)

! DART-style longitude
where(hruLon <    0.0_r8) hruLon = hruLon + 360.0_r8
where(hruLon == 360.0_r8) hruLon = 0.0_r8

call nc_close_file(ncid, routine)

end subroutine read_hru_properties


!------------------------------------------------------------------
! Returns the number of items in the state vector as an integer. 
! This interface is required for all applications.

function get_model_size()

integer(i8) :: get_model_size

get_model_size = model_size

end function get_model_size


!------------------------------------------------------------------
! Read model time from restart file

function read_model_time(filename)

type(time_type) :: read_model_time, base_hydro_date
character(len=*), intent(in) :: filename
character(len=*), parameter  :: routine = 'read_model_time'

integer, parameter          :: STRINGLENGTH = 30
character(len=STRINGLENGTH) :: datestr
integer                     :: ncid, days, pos
integer                     :: year, month, day, hour, minute, second

! Need to get time from the streamflow restart file 
ncid = nc_open_file_readonly(filename, routine)

call nc_get_variable(ncid, 'time', days, routine)
call nc_get_attribute_from_variable(ncid, 'time', 'units', datestr)

pos = scan(datestr, "0123456789") 

! Read "Date since"
read(datestr(pos:STRINGLENGTH), '(i4, 5(1x,i2))') year, month, day, hour, minute, second

base_hydro_date = set_date(year, month, day, hour, minute, second)
read_model_time = increment_time(base_hydro_date, 0, days)

if (debug > 99) then 
   write(*, *) 'Time string from restart is "'//trim(datestr)//'"'
   call print_date(read_model_time, ' Valid time is ')
endif

call nc_close_file(ncid, routine)

end function read_model_time


!------------------------------------------------------------------
! write any additional attributes to the output and diagnostic files

subroutine nc_write_model_atts(ncid, domain_id)

integer, intent(in) :: ncid      ! netCDF file identifier
integer, intent(in) :: domain_id

if (.not. module_initialized) call static_init_model

! put file into define mode
call nc_begin_define_mode(ncid)

call nc_add_global_creation_time(ncid)

call nc_add_global_attribute(ncid, "model_source", source )
call nc_add_global_attribute(ncid, "model", "pywatershed")

call nc_end_define_mode(ncid)

! Flush the buffer and leave netCDF file open
call nc_synchronize_file(ncid)

end subroutine nc_write_model_atts

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
! Set the 3D location array for ach domain

subroutine configure_domloc()

character(len=*), parameter :: routine = 'configure_domloc'

integer :: idom, iseg, ihru

do idom = 1, domain_count
   if (idom == idom_chn) then 
      allocate(domain_info(idom)%location(nseg, JINDEX, KINDEX))
      do iseg = 1, nseg
         domain_info(idom)%location(iseg, JINDEX, KINDEX) = &
                                 set_location(segLon(iseg), &
                                              segLat(iseg), & 
                                              segElv(iseg), &
                                              VERTISHEIGHT)
      enddo
   elseif (idom == idom_hru) then
      allocate(domain_info(idom)%location(nhru, JINDEX, KINDEX))
      do ihru = 1, nhru
         domain_info(idom)%location(ihru, JINDEX, KINDEX) = &
                                 set_location(hruLon(ihru), &
                                              hruLat(ihru), & 
                                              hruElv(ihru), &
                                              VERTISHEIGHT)
      enddo
   elseif (idom == idom_prm) then 
      msg1 = 'Parameter locations are not supported for now.'
      call error_handler(E_ERR, routine, msg1, source)  
   endif
enddo

end subroutine configure_domloc


!------------------------------------------------------------------
! Given an integer index into the state vector structure, returns the
! associated location. A second intent(out) optional argument quantity
! (qty: QTY_STREAM_FLOW). This interface is
! required for all filter applications as it is required for computing
! the distance between observations and state variables.

subroutine get_state_meta_data(index_in, location, qty_type)

integer(i8),         intent(in)            :: index_in
type(location_type), intent(out)           :: location
integer,             intent(out), optional :: qty_type

integer :: iloc, jloc, kloc, varid, idom

if (.not. module_initialized) call static_init_model

call get_model_variable_indices(index_in, iloc, jloc, kloc, varid, idom, qty_type)

location = domain_info(idom)%location(iloc, jloc, kloc)

if (debug > 99) then
   call write_location(0, location, charstring=msg1)
   write(*, '(A, i5, i5, i3, i3, 2x, A)') 'gsmd index,i,j,k = ', index_in, iloc, jloc, kloc, trim(msg1)
endif

end subroutine get_state_meta_data


!------------------------------------------------------------------
! Starting from an initial state, generate perturbations to form an 
! ensemble. 

subroutine pert_model_copies(state_ens_handle, ens_size, pert_amp, interf_provided)

character(len=*), parameter :: routine = 'pert_model_copies'

type(ensemble_type), intent(inout) :: state_ens_handle
integer,             intent(in)    :: ens_size 
real(r8),            intent(in)    :: pert_amp
logical,             intent(out)   :: interf_provided

logical, save         :: seed_unset = .true.
integer, save         :: seed             
type(random_seq_type) :: random_seq

integer               :: idom, ind1, ind2, ivar, jvar, iens
real(r8)              :: pert

if (.not. module_initialized) call static_init_model

! let's do our own perturbation
interf_provided = .true. 

! Generate a unique (but repeatable - if task count is same)
! seed for each task
if (seed_unset) then
   seed = (my_task_id()+1) * 1000
   seed_unset = .false.
endif

call init_random_seq(random_seq, seed)
seed = seed + 1  

pert = model_perturbation_amplitude

DOMAINS: do idom = 1, domain_count
   ! Skip the parameters for now
   if (idom == idom_prm) cycle DOMAINS

   ! Both Channel & HRU domains have 1 variable in them for now
   ! So, get_num_variables(idom) should be just 1   
   do ivar = 1, get_num_variables(idom)
      ind1 = get_index_start(idom, ivar)
      ind2 = get_index_end(  idom, ivar)
      ! Could change perturbation size for different variables 
      ! in different domains
      do jvar = 1, state_ens_handle%my_num_vars
         if (state_ens_handle%my_vars(jvar) >= ind1 .and. &
             state_ens_handle%my_vars(jvar) <= ind2) then
            ! Generate the ensemble members (copies) one-by-one  
            do iens = 1, ens_size
               if (trim(perturb_distribution) == 'lognormal') then
                  state_ens_handle%copies(iens, jvar) = state_ens_handle%copies(iens, jvar) &
                                 * exp(pert * random_gaussian(random_seq, 0.0_r8, 1.0_r8))             
               elseif (trim(perturb_distribution) == 'normal') then 
                  state_ens_handle%copies(iens, jvar) =  random_gaussian(random_seq, &
                                   state_ens_handle%copies(iens, jvar), pert)
               else
                  write(msg1, *) 'pert_model_copies'
                  write(msg2, *) 'Selected distribution "', trim(perturb_distribution), &
                                 '" is not currently supported.'
                  call error_handler(E_ERR, routine, msg1, source, text2=msg2)
               endif
            enddo ! Ensemble Members
         endif
      enddo ! Ensemble Handle Variables
   enddo ! Domain Variables
enddo DOMAINS

end subroutine pert_model_copies


!------------------------------------------------------------------
! ATS Localization: Find close obs along the network
! Only supports identity streamflow obs

subroutine get_close_obs(gc, base_loc, base_type, locs, loc_qtys, & 
               loc_types, num_close, close_ind, dist, ens_handle)

type(get_close_type), intent(in)            :: gc           ! handle to a get_close structure
type(location_type),  intent(inout)         :: base_loc     ! location of interest
type(location_type),  intent(inout)         :: locs(:)      ! all locations on my task
integer,              intent(in)            :: base_type    ! observation TYPE
integer,              intent(in)            :: loc_qtys(:)  ! QTYs for all locations
integer,              intent(in)            :: loc_types(:) ! TYPEs for all locations
integer,              intent(out)           :: num_close    ! how many are close
integer,              intent(out)           :: close_ind(:) ! index on task of close locs
real(r8),             optional, intent(out) :: dist(:)      ! distances (in radians)
type(ensemble_type),  optional, intent(in)  :: ens_handle   

character(len=*), parameter :: routine = 'get_close_obs'

!==================== local variables ====================
integer             :: base_qty
integer(i8)         :: full_index       ! gsmd requires a larger integer into DART state index
type(location_type) :: location

integer             :: stream_nclose
integer(i8)         :: stream_indices(nseg), loc_indx(size(close_ind))
real(r8)            :: stream_dists(nseg)

! Only supporting identity streamflow obs for now!!
if (base_type > 0) then 
   write(msg1, *) 'Issue with the observation type:'
   write(msg2, *) 'Localization is only done for streamflow obs using ', &
                  'ATS localization strategy.'
   write(msg3, *) 'For now, anything that is not identity obs will fail.'
   call error_handler(E_ERR, routine, msg1, source, text2=msg2, text3=msg3)
endif

! Assert that this is a streamflow obs
full_index = abs(base_type)
call get_state_meta_data(full_index, location, base_qty)
if (base_qty /= QTY_STREAM_FLOW) then 
   write(msg1, *) 'Issue with the observation quantity:'
   write(msg2, *) 'The observation must be QTY_STREAM_FLOW: ', &
                  QTY_STREAM_FLOW
   write(msg3, *) 'The used obs QTY is: ', base_qty
   call error_handler(E_ERR, routine, msg1, source, text2=msg2, text3=msg3)  
endif

! We now have an identity streamflow observation
! Find close downstream+upstream streamflow obs
call close_streams(full_index, stream_nclose, stream_indices, stream_dists)  

!TODO: Should I subset further based on which streams have gauges (obs)?

! Determine which of the global indices are mine
loc_indx = abs(loc_types)
call get_my_close(stream_nclose, stream_indices, stream_dists, loc_indx, &
                  num_close, close_ind, dist)

if (debug > 99) then
   write(msg1,'("PE ", i3)') my_task_id()
   write(*, *) trim(msg1), ' get_close_obs: num_close ', num_close
   write(*, *) trim(msg1), ' get_close_obs: close_ind ', close_ind(1:num_close)
   write(*, *) trim(msg1), ' get_close_obs: dist      ', dist(1:num_close)
endif

end subroutine get_close_obs


!------------------------------------------------------------------
! Determine streams that are close to me. Return a list of indices 
! into the DART state vector. 

subroutine close_streams(full_index, num_close, close_indices, distances)

integer(i8), intent(in)  :: full_index       ! identity observation
integer,     intent(out) :: num_close        ! Number of close streams
integer(i8), intent(out) :: close_indices(:) ! full DART indices
real(r8),    intent(out) :: distances(:)     ! in radians

integer     :: depth, seg_index
real(r8)    :: reach_length

integer     :: stream_nclose
integer(i8) :: stream_indices(nseg)
real(r8)    :: stream_distances(nseg)

! Initialize variable
seg_index     = full_index
depth         = 1            ! starting value for recursive routine
reach_length  = 0.0_r8       ! initial length along the network
stream_nclose = 0            ! # of close streams
distances     = huge(1.0_r8) ! initially: entire array is far away

! Determine the upstream close streams
call link_tree(seg_index, max_link_distance, depth, &
                   reach_length, stream_nclose, stream_indices, stream_distances)

! Add the close downstreams ones to the list
call downstream_links(seg_index, max_link_distance, &
                   stream_nclose, stream_indices, stream_distances)

num_close                  = stream_nclose
close_indices(1:num_close) = stream_indices(1:stream_nclose)
distances(1:num_close)     = stream_distances(1:stream_nclose) / 1000.0_r8  ! to km
distances(1:num_close)     = distances(1:num_close) / earth_radius          ! to radians

if (debug > 99) then
   write(msg1, '("PE ",I3)') my_task_id()
   write(*, *) trim(msg1), ' gc_streams: num_close     ', num_close
   write(*, *) trim(msg1), ' gc_streams: close_indices ', close_indices(1:num_close)
   write(*, *) trim(msg1), ' gc_streams: distances     ', distances(1:num_close)
endif

end subroutine close_streams


!-----------------------------------------------------------------------
! Recurrsive routine to find upstream segments for ATS localization 

recursive subroutine link_tree(my_index, reach_cutoff, depth, &
                    reach_length, nclose, close_indices, distances)

integer,     intent(in)    :: my_index
real(r8),    intent(in)    :: reach_cutoff   ! meters
integer,     intent(in)    :: depth
real(r8),    intent(inout) :: reach_length
integer,     intent(inout) :: nclose
integer(i8), intent(inout) :: close_indices(:)
real(r8),    intent(inout) :: distances(:)

real(r8) :: direct_length
integer  :: iup

direct_length = reach_length

if (depth > 1) direct_length = direct_length + connections(my_index)%segLength

! reach_length may also need to include elevation change
if (direct_length >= reach_cutoff) return

nclose                = nclose + 1
close_indices(nclose) = my_index
distances(nclose)     = direct_length

if (debug > 99) then 
   write(*, '(A, i5, f10.2, i8, i4, 5i8)') 'depth, distance, node, # upstream, up nodes: ', &
                     depth, direct_length, my_index, num_up_links(my_index),                &
                     connections(my_index)%upstream_index(:)
endif

do iup = 1, num_up_links(my_index)
   call link_tree(connections(my_index)%upstream_index(iup), &
            reach_cutoff, depth+1, direct_length, nclose, close_indices, distances)
enddo

end subroutine link_tree


!-----------------------------------------------------------------------
! Routine to collect the close downstream segments from me

subroutine downstream_links(my_index, reach_cutoff, &
                    nclose, close_indices, distances)

integer,     intent(in)    :: my_index
real(r8),    intent(in)    :: reach_cutoff   ! meters
integer,     intent(inout) :: nclose
integer(i8), intent(inout) :: close_indices(:)
real(r8),    intent(inout) :: distances(:)

real(r8) :: direct_length
integer  :: idown

idown = my_index

direct_length = connections(idown)%segLength

do while(direct_length < reach_cutoff .and. &
         connections(idown)%downstream_index > 0)

   idown = connections(idown)%downstream_index

   nclose                = nclose + 1
   close_indices(nclose) = idown
   distances(nclose)     = direct_length

   direct_length = direct_length + connections(idown)%segLength
enddo

end subroutine downstream_links


!-----------------------------------------------------------------------
! determine the indices of the objects that are on this task
!
! num_superset         length of the superset of desired indices
! superset_indices     superset of desired indices
! superset_distances   superset of desired distances
! my_task_indices      all possible indices ON MY TASK (i.e. a local subset)
! num_close            length of desired elements ON MY TASK (i.e. intersection)
! close_ind(:)         the list of desired indices ON MY TASK
! dist(:)              the corresponding distances of the desired indices

subroutine get_my_close(num_superset, superset_indices, superset_distances, &
                        my_task_indices, num_close, close_ind, dist)

integer,          intent(in)  :: num_superset
integer(i8),      intent(in)  :: superset_indices(:)
real(r8),         intent(in)  :: superset_distances(:)
integer(i8),      intent(in)  :: my_task_indices(:)
integer,          intent(out) :: num_close
integer,          intent(out) :: close_ind(:)
real(r8),         intent(out) :: dist(:)

integer, dimension(:), allocatable :: index_map
integer :: i, idx, il, ir

num_close = 0

! Determine the range of my_task_indices
il = minval(my_task_indices)
ir = maxval(my_task_indices)

! Create a map for quick lookup
allocate(index_map(il:ir))
index_map = 0
do i = 1, num_superset
  idx = superset_indices(i)
  if (idx >= il .and. idx <= ir) then
    index_map(idx) = i
  end if
end do

! Loop over my_task_indices and find matches using the map
do i = 1, size(my_task_indices)
  idx = my_task_indices(i)
  if (index_map(idx) > 0) then
     num_close = num_close + 1
     close_ind(num_close) = i 
     dist(num_close) = superset_distances(index_map(idx))
  end if
end do

! Deallocate the map
deallocate(index_map)

end subroutine get_my_close


!------------------------------------------------------------------
! ATS Localization: Find close state along the network

subroutine get_close_state(gc, base_loc, base_type, locs, loc_qtys, & 
                  loc_indx, num_close, close_ind, dist, ens_handle)

type(get_close_type), intent(in)            :: gc           ! handle to a get_close structure
type(location_type),  intent(inout)         :: base_loc     ! location of interest
integer,              intent(in)            :: base_type    ! observation TYPE
type(location_type),  intent(inout)         :: locs(:)      ! locations on my task
integer,              intent(in)            :: loc_qtys(:)  ! QTYs for locations on my task
integer(i8),          intent(in)            :: loc_indx(:)  ! indices into DART state on my task
integer,              intent(out)           :: num_close    ! how many are close
integer,              intent(out)           :: close_ind(:) ! indices in the the locs array
real(r8),             optional, intent(out) :: dist(:)      ! distances
type(ensemble_type),  optional, intent(in)  :: ens_handle

character(len=*), parameter :: routine = 'get_close_state'

!====================== local variables ======================
integer     :: num_vars_chn, num_vars_hru, contributing_HRUs
integer     :: num_close_seg, iclose, start_index, close_index
integer     :: ihru, index_skip, close_hru_index, local_hru_index
integer(i8) :: full_index 

integer     :: stream_nclose
integer(i8) :: stream_indices(nseg) 
real(r8)    :: stream_dists(nseg)

full_index   = abs(base_type)
num_vars_chn = get_num_variables(idom_chn)

if (num_vars_chn > 1) then 
   write(msg1, *) 'Current implementation of ATS Localization ',   &
                  'only supports streamflow and groundwater.'
   write(msg2, *) 'If there are more variables in each of the ',   &
                  'Channel and HRU domains, filter will fail.'
   write(msg3, *) 'Localizing other variables require knowledge ', &
                  'of their respective grids.'
   call error_handler(E_ERR, routine, msg1, source, text2=msg2, text3=msg3) 
endif

! We now have an identity streamflow observation
! Find close downstream+upstream streamflow obs
call close_streams(full_index, stream_nclose, stream_indices, stream_dists)

! Add groundwater variables if we're updating them 
if (domain_count > 1) then 

   num_vars_hru = get_num_variables(idom_hru)
   if (num_vars_chn > 1 .or. & 
       get_variable_name(idom_hru, 1) /= 'gwres_stor') then 
      write(msg1, *) 'Only supports groundwater in HRU domain!'
      call error_handler(E_ERR, routine, msg1, source)
   endif

   ! Remove groundwater in HRUs that don't contribute to the 
   ! close segments and add the one(s) that do contribute to the flow 
   num_close_seg = stream_nclose
   start_index   = get_index_start(idom_hru, 1)
   index_skip    = start_index - 1
   
   GWLOOP: do iclose = 1, num_close_seg
      close_index = stream_indices(iclose)

      if (hruMask(close_index) == MASKED_HRU) cycle GWLOOP
      
      contributing_HRUs = hruMask(close_index)

      if (debug > 99) print *, 'iclose:', iclose,  &
                      'close_index:', close_index, & 
                      'HRU count:', contributing_HRUs

      ! Add them to the list
      HRULOOP : do ihru = 1, contributing_HRUs
         local_hru_index = connections(close_index)%avail_hru(ihru)
         close_hru_index = index_skip + local_hru_index

         ! Found HRU index, now make sure it's active 
         if (hruTyp(local_hru_index) == INACTIVE_HRU) cycle HRULOOP

         stream_nclose                 = stream_nclose + 1
         stream_indices(stream_nclose) = close_hru_index
         stream_dists(  stream_nclose) = stream_dists(iclose)            
      enddo HRULOOP       
   enddo GWLOOP
   
   if (debug > 99) write(*, '(A, i5, X, A)') &
                   'In total, we added', stream_nclose-num_close_seg, & 
                   'HRUs to the list of close streamflow segments'
endif

! Determine which of the global indices are mine
call get_my_close(stream_nclose, stream_indices, stream_dists, loc_indx, &
                  num_close, close_ind, dist)

if (debug > 99) then
   write(msg1, '("PE ", i3)') my_task_id()
   write(*, *) trim(msg1), ' get_close_state: num_close ', num_close
   write(*, *) trim(msg1), ' get_close_state: close_ind ', close_ind(1:num_close)
   write(*, *) trim(msg1), ' get_close_state: dist      ', dist(1:num_close)
endif

end subroutine get_close_state

!------------------------------------------------------------------
! Returns the number of links in the state vector

function get_number_of_segments()

integer(i8) :: get_number_of_segments

get_number_of_segments = nseg

end function get_number_of_segments


! -------------------------------------
! Read array of strings from input file

subroutine get_string_array(ncid, dimname, dimlen, varname, string_var)
   
integer,          intent(in)    :: ncid, dimlen
character(len=*), intent(in)    :: dimname, varname
character(len=*), intent(inout) :: string_var(:)
   
character(len=*), parameter  :: routine = 'get_string_array'  

integer :: i, io, varid
      
io = nf90_inq_varid(ncid, varname, varid)
call nc_check(io, routine, 'inq_varid "'//varname//'"')
  
io = nf90_get_var(ncid, varid, string_var)
call nc_check(io, routine, 'get_var "'//varname//'"')
   
do i = 1, dimlen
   string_var(i) = adjustl(string_var(i))
enddo

if (debug > 99) then
   do i = 1, dimlen
      write(msg1, *) i, ': "'//string_var(i)//'"'
      call error_handler(E_MSG, routine, msg1)
   enddo
endif

end subroutine get_string_array

!===================================================================
! End of model_mod
!===================================================================
end module model_mod

