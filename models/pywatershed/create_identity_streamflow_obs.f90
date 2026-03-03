! DART software - Copyright UCAR. This open source software is provided
! by UCAR, "as is", without charge, subject to all terms of use at
! http://www.image.ucar.edu/DAReS/DART/DART_download
!
! DART $Id$

program create_identity_streamflow_obs

!---------------------------------------------------------------------------
! create_identity_streamflow_obs - reads streamflow data and a
! pywatershed model state and writes a DART obs_seq file where the
! observations are actually identity observations.
!---------------------------------------------------------------------------

use types_mod,              only : r8, i8
use location_mod,           only : VERTISHEIGHT, location_type, get_location
use utilities_mod,          only : nmlfileunit, do_nml_file, do_nml_term,        &
                                   initialize_utilities, finalize_utilities,     &
                                   find_namelist_in_file, check_namelist_read,   &
                                   error_handler, E_ERR, E_MSG, open_file,       &
                                   find_textfile_dims, close_file 
use time_manager_mod,       only : time_type, set_calendar_type, GREGORIAN,      &
                                   set_time, get_time, set_date
use obs_sequence_mod,       only : obs_sequence_type, obs_type, read_obs_seq,    &
                                   static_init_obs_sequence, write_obs_seq,      &
                                   init_obs, init_obs_sequence, get_num_obs,     &
                                   set_copy_meta_data, set_qc_meta_data,         &
                                   get_num_copies, get_num_qc
use obs_kind_mod,           only : QTY_STREAM_FLOW
use obs_def_streamflow_mod, only : set_streamflow_metadata
use obs_utilities_mod,      only : add_obs_to_seq, create_3d_obs, getdimlen
use model_mod,              only : static_init_model, get_state_meta_data,       &
                                   get_number_of_segments
use sort_mod,               only : index_sort
use netcdf_utilities_mod,   only : nc_check, nc_get_dimension_size,              &
                                   nc_open_file_readonly, nc_get_variable,       &
                                   nc_close_file, nc_get_attribute_from_variable
use netcdf

implicit none

character(len=*), parameter :: source  = "create_identity_streamflow_obs.f90"

! Input file (USGS) has data quality control fields, whether to use or ignore them.
integer,  parameter :: NUM_COPIES      = 1            ! number of copies in sequence
integer,  parameter :: NUM_QC          = 1            ! number of QC entries
real(r8), parameter :: MIN_OBS_ERR_STD = 0.1_r8       ! ft^3/sec
real(r8), parameter :: MAX_OBS_ERR_STD = 1000000.0_r8 ! obs error shouldn't exceed
real(r8), parameter :: UNIT_CONVERSION = 35.3147_r8   ! convert from cms to cfs

! These should/must match what is in the netCDF files ... should error check
! or at least provide a decent error message.

integer, parameter  :: IDLEN        = 15
integer, parameter  :: STATIDSTRLEN = 15
integer, parameter  :: TIMESTRLEN   = 19

integer             :: nseg, ngages, nobs, nfiles 

character(len=512)  :: msg1, msg2, msg3
character(len=34)   :: prognml = 'create_identity_streamflow_obs_nml'
character(len=9)    :: dartnml = 'input.nml'

! Variables read from the parameters_dis_seg_app.nc file
real(r8),             allocatable :: lat(:)
real(r8),             allocatable :: lon(:)
integer,              allocatable :: seg(:)  ! NHM segment ID
real(r8),             allocatable :: elv(:)
character(len=IDLEN), allocatable :: gauge_strings(:)

! Variables read from observation data file
character(len=TIMESTRLEN),   allocatable :: time_str(:)
character(len=STATIDSTRLEN), allocatable :: stations(:)
real(r8),                    allocatable :: discharge(:)
character(len=IDLEN),        allocatable :: desired_gages(:)

! structure to hold the lat/lon/depth/indexing information
! for the DART state vector
type database
   integer(i8)              :: nsegments
   real(r8),    allocatable :: lat(:)
   real(r8),    allocatable :: lon(:)
   real(r8),    allocatable :: dep(:)
   integer(i8), allocatable :: sortedindex(:)
end type database
type(database) :: lookup_table

type gauge_lookup
   character(len=256), allocatable :: gauge_keys(:)
   integer, allocatable            :: gauge_indices(:)
end type gauge_lookup
type(gauge_lookup) :: gauge_map

type(obs_type)          :: obs 
type(obs_sequence_type) :: obs_seq
type(obs_type)          :: prev_obs
type(time_type)         :: time_obs, prev_time

! namelist variables
character(len=256) :: input_files            = ''            ! Input raw data file; timeslices from USGS
character(len=256) :: output_file            = 'obs_seq.out' ! Output DART-style obs-seq file
character(len=256) :: network_file           = 'dis_seg.nc'  ! File with geometry and network information
character(len=256) :: gages_list_file        = ''            ! List of gauges to be assimilated
real(r8)           :: obs_fraction_for_error = 0.01          ! Parameter used to parametrize obs error
logical            :: assimilate_all         = .false.       ! Flag to turn on assimilation of all available gauges
integer            :: debug                  = 0             ! Verbosity

namelist /create_identity_streamflow_obs_nml/ &
               input_files,                   &
               output_file,                   &
               network_file,                  &
               gages_list_file,               &
               obs_fraction_for_error,        &
               assimilate_all,                &
               debug

!-------------------------------------------------------------------------------

! Initialize and read namelist
call initialize_utilities(source)
call read_namelist() 

! Put the reference date into DART format
call set_calendar_type(GREGORIAN)

! Initialize the DART state vector so we can relate the observation location
! to the location in the DART state vector.
call static_init_model()

! Read location info from the network file 
call read_network_information()

! Prepare output obs seq file: create it || append to it
call prep_obs_seq() 

! Now, start reading the obs from slices and add
! them to the obs-seq file
call add_obs_from_slices() 

! end of main program

contains


!---------------------------------------------
! Read the namelist from DART's input.nml file

subroutine read_namelist()

integer :: iunit, io

! Read the DART namelist
call find_namelist_in_file(dartnml, prognml, iunit)
read(iunit, nml = create_identity_streamflow_obs_nml, iostat = io) 
call check_namelist_read(iunit, io, prognml)

! Print the content of the namelist for clarification
if (debug > 1) then 
   print *, 'list of input files is in : "'//trim(input_files)//'"'
   print *, 'output_file               : "'//trim(output_file)//'"'
   print *, 'network_file              : "'//trim(network_file)//'"'
   print *, 'gages_list_file           : "'//trim(gages_list_file)//'"'
   print *, 'obs_fraction_for_error    : ', obs_fraction_for_error
   print *, 'debug                     : ', debug
endif

! Record the DART namelist values used for the run ...
if (do_nml_file()) write(nmlfileunit, nml=create_identity_streamflow_obs_nml)
if (do_nml_term()) write(     *     , nml=create_identity_streamflow_obs_nml)

end subroutine read_namelist


!--------------------------------------------------------
! Open the network file and retrieve location information

subroutine read_network_information()

character(len=*), parameter  :: routine = 'read_network_information'

integer              :: ncid, igage, iseg
integer, allocatable :: segGauge(:)

character(len=STATIDSTRLEN), allocatable :: gID(:)

ncid = nc_open_file_readonly(network_file, routine)     

! Read in the dimensions
nseg   = nc_get_dimension_size(ncid, 'nsegment' , routine)
ngages = nc_get_dimension_size(ncid, 'npoigages', routine)

allocate(lat(nseg), lon(nseg), elv(nseg), &
         seg(nseg), gauge_strings(nseg),  &
         gID(ngages), segGauge(ngages))

! Read segment properties
call nc_get_variable(ncid, 'nhm_seg'         , seg, routine)
call nc_get_variable(ncid, 'seg_lats'        , lat, routine)
call nc_get_variable(ncid, 'seg_lons'        , lon, routine)
call nc_get_variable(ncid, 'seg_elev'        , elv, routine)
call nc_get_variable(ncid, 'poi_gage_segment', segGauge, routine)

call get_string_array(ncid, 'IDLength', ngages, 'poi_gauges', gID)

do igage = 1, ngages
   gauge_strings(segGauge(igage)) = gID(igage)
enddo

! Build gauges lookup table **once** before any searches
call build_gauge_lookup()

! DART-style longitude for the segments
where(lon <    0.0_r8) lon = lon + 360.0_r8
where(lon == 360.0_r8) lon = 0.0_r8

call nc_close_file(ncid, routine)

if (debug > 5) then 
   do iseg = 1, nseg
      write(*, '(A, i5, X, A, i6, X, A, f8.3, X, A, f8.3, X, A, f8.3, X, A, X, A)') &
               'num:', iseg, 'segID:', seg(iseg), 'lat:', lat(iseg), &
               'lon:', lon(iseg), 'elv:', elv(iseg), &
               'gauge:', adjustl(gauge_strings(iseg))
   enddo
endif

end subroutine read_network_information


!--------------------------------------------------------
! Open already available obs-seq file or create a new one 

subroutine prep_obs_seq()

character(len=*), parameter  :: routine = 'prep_obs_seq'

integer :: num_new_obs, n_desired_gages
logical :: file_exist

integer :: existing_num_copies, existing_num_qc

call static_init_obs_sequence()

call init_obs(obs,      num_copies=NUM_COPIES, num_qc=NUM_QC)
call init_obs(prev_obs, num_copies=NUM_COPIES, num_qc=NUM_QC)

! Check if the user has provided a list of gauage to assimilate 
n_desired_gages = set_desired_gages(gages_list_file)

! How many timeslices are there in my input list
call find_textfile_dims(input_files, nfiles)

! some estimate of the number of obs
num_new_obs = estimate_total_obs_count(input_files, nfiles)

! Begin here the check
inquire(file=output_file, exist=file_exist)

if (file_exist) then ! existing file found, append to it

   write(msg1, *) "Found existing obs_seq file, appending to ", trim(output_file)
   write(msg2, *) "Adding up to a maximum of ", num_new_obs, " new observations"
   call error_handler(E_MSG, routine, msg1, source, text2=msg2)

   call read_obs_seq(output_file, 0, 0, num_new_obs, obs_seq)

   ! check to see if existing file is compatible
   existing_num_copies = get_num_copies(obs_seq)
   existing_num_qc     = get_num_qc(obs_seq)

   if (existing_num_copies /= NUM_COPIES .or.  existing_num_qc /= NUM_QC) then
      write(msg1, *) 'Incompatible existing observation sequence file'
      write(msg2, '(A,i4,A,i4)') 'Expected ',NUM_COPIES, &
                                 ' obs copies got ',existing_num_copies
      write(msg3, '(A,i4,A,i4)') 'Expected ',NUM_QC, &
                                 ' QC  copies got ',existing_num_qc
      call error_handler(E_ERR, routine, msg1, source, text2=msg2, text3=msg3)
   endif

else ! create a new one ...

   call init_obs_sequence(obs_seq, NUM_COPIES, NUM_QC, num_new_obs)
   call set_copy_meta_data(obs_seq, NUM_COPIES, 'observation')
   call set_qc_meta_data(obs_seq, NUM_QC, 'Data QC')

   write(msg1, *) "No existing obs_seq file, creating ", trim(output_file)
   write(msg2, *) "with up to a maximum of ", num_new_obs, " observations"

   call error_handler(E_MSG, routine, msg1, source, text2=msg2)
endif

end subroutine prep_obs_seq


!-------------------------------------------------------------
! Go through the time slices and add the data to obs_seq.out

subroutine add_obs_from_slices()

integer                 :: ifile, ncid, iob, obind
integer                 :: iunit, io, dart_index, key
integer                 :: oday, osec
real(r8)                :: oerr, qc
logical                 :: first_obs
character(len=256)      :: Q_unit, input_file

character(len=*), parameter  :: routine = 'add_obs_from_slices'

iunit = open_file(input_files, form='formatted', action='read')

first_obs = .true. 

FILELOOP : do ifile = 1, nfiles 

   ! Read in the name of the timeslices file
   read(iunit, '(A)', iostat=io) input_file
   if (io /= 0 ) then
     write(msg1, *) 'Unable to read input file from "'//trim(input_files)//'"'
     write(msg2, *) 'file ', ifile
     call error_handler(E_ERR, routine, msg1, source, text2=msg2)
   endif
   
   ! For the current file, get the number of available obs
   ncid = nc_open_file_readonly(input_file, routine)
   nobs = nc_get_dimension_size(ncid, 'stationIdInd' , routine) 

   write(msg1, *) 'Reading time slice file "',trim(input_file)//'"'
   write(msg2, *) 'number of obs in the time slice is ', nobs
   call error_handler(E_MSG, routine, msg1, text2=msg2)

   allocate(discharge(nobs), stations(nobs), time_str(nobs))

   ! Read the discharge and the associated quality
   call nc_get_variable(ncid, 'discharge'  , discharge  , routine)
   call nc_get_attribute_from_variable(ncid, 'discharge', 'units', Q_unit, routine)

   if (Q_unit == 'm^3/s') then 
      write(msg1, *) "Model's streamflow is in cfs. Converting discharge from cms to cfs ..."
      call error_handler(E_MSG, routine, msg1)

      discharge = discharge * UNIT_CONVERSION
   endif 

   ! Read the time and station string arrays
   call get_string_array(ncid, 'stationIdStrLen', nobs, 'stationId', stations)
   call get_string_array(ncid, 'timeStrLen'     , nobs, 'time'     , time_str)

   call nc_close_file(ncid, routine) 

   ! Process the observations one by one
   OBSLOOP : do iob = 1, nobs
   
      ! Make sure discharge is physically meaningful 
      if (discharge(iob) < 0.0_r8 .or. &
          discharge(iob) /= discharge(iob)) cycle OBSLOOP

      obind = find_matching_gauge(iob) 
      if (obind == 0) cycle OBSLOOP

      ! Relate the physical location to the dart state vector index
      dart_index = linkloc_to_dart(lat(obind), lon(obind))

      ! Desired gauges get the provided obs error
      ! Remaining ones are for verification purposes
      if (ANY(desired_gages == stations(iob)) .or. assimilate_all) then
        oerr = max(discharge(iob)*obs_fraction_for_error, MIN_OBS_ERR_STD)
      else
        oerr = MAX_OBS_ERR_STD
      endif

      call convert_time_string(time_str(iob), oday, osec, iob)
      time_obs = set_time(osec,oday) 

      call set_streamflow_metadata(key, stations(iob), seg(obind))

      call create_3d_obs(lat(obind), lon(obind), elv(obind), VERTISHEIGHT, &
                         discharge(iob), dart_index, oerr, oday, osec, qc, obs, key)
      
      call add_obs_to_seq(obs_seq, obs, time_obs, prev_obs, prev_time, first_obs)

   enddo OBSLOOP
   
   ! Release memory associated with the current obs file
   deallocate(discharge, time_str, stations)
enddo FILELOOP

nobs = get_num_obs(obs_seq)

! If we added any obs to the sequence, write it now.
if (nobs > 0)  call write_obs_seq(obs_seq, output_file)

if (debug > 1) then
   write(msg1, *) 'Writing ', nobs, ' observations to file.'
   call error_handler(E_MSG, routine, msg1)
endif

call close_file(iunit)

deallocate(lat, lon, elv, seg, gauge_strings)
if (allocated(desired_gages)) deallocate(desired_gages)

call finalize_utilities()

end subroutine add_obs_from_slices

!-------------------------------------------------------------------------------
! Read the list of the gages to process.
! If the filename is empty, use all the gages.
! The list of desired gages is globally scoped:  desired_gages()

function set_desired_gages(filename) result(ngages)

character(len=*), intent(in) :: filename
integer                      :: ngages

character(len=*), parameter  :: routine = 'set_desired_gages'

integer :: nlines, i, iunit, io

ngages = 0

if (filename == '' .or. filename == 'null') then
   write(msg1, *) 'No gages_list_file supplied. Using all the gages.'
   call error_handler(E_MSG, routine, msg1)
   return
endif               
   
call find_textfile_dims(filename, nlines)

allocate(desired_gages(nlines))

iunit = open_file(filename, form='formatted', action='read')
do i = 1, nlines
   read(iunit, *, iostat=io)  desired_gages(i)
   if (io /= 0) then
      write(msg1, *) 'Unable to read gage from "'//trim(filename)//'"'
      write(msg2, *) 'line ',i
      call error_handler(E_ERR, routine, msg1, source, text2=msg2)
   endif
enddo
call close_file(iunit)

if (debug > 5) then
   print *, ''
   print *, 'List of the gauges to process:'
   do i = 1, nlines
      print *, i, desired_gages(i)
   enddo
endif

ngages = nlines

end function set_desired_gages


!---------------------------------------------------
! Rough estimate of the total number of observations

function estimate_total_obs_count(file_list, nfiles) result (num_obs)
  
character(len=*), intent(in) :: file_list
integer                      :: num_obs, nfiles
   
character(len=*), parameter :: routine = 'estimate_total_obs_count'

integer            :: ncid, no, iunit, io
character(len=256) :: input_file

! Use the first observation timeslice file in the list                     
iunit = open_file(file_list, form='formatted', action='read')
read(iunit, '(A)', iostat=io) input_file 
if (io /= 0 ) then
  write(msg1, *) 'Unable to read input file from "'//trim(file_list)//'"'
  call error_handler(E_ERR, routine, msg1, source)
endif
call close_file(iunit) 

! Need to know about how many observations are in each file.
io = nf90_open(input_file, nf90_nowrite, ncid) 
call nc_check(io, routine, 'opening data file "'//trim(input_file)//'"' )
call getdimlen(ncid, 'stationIdInd', no)
io = nf90_close(ncid)
call nc_check(io, routine, 'closing file "'//trim(input_file)//'"' )

! We need to know how many observations there may be.
! Specifying too many is not really a problem.
! I am multiplying by 10.

num_obs = 10 * no * nfiles

end function estimate_total_obs_count


! -------------------------------------
! Read array of strings from input file

subroutine get_string_array(ncid, dimname, dimlen, varname, string_var)

integer,          intent(in)    :: ncid, dimlen
character(len=*), intent(in)    :: dimname, varname
character(len=*), intent(inout) :: string_var(:)

character(len=*), parameter  :: routine = 'get_string_array'  
 
integer :: i, io, dim1, varid

call getdimlen(ncid, dimname, dim1)

io = nf90_inq_varid(ncid, varname, varid)
call nc_check(io, routine, 'inq_varid "'//varname//'"')
   
io = nf90_get_var(ncid, varid, string_var)
call nc_check(io, routine, 'get_var "'//varname//'"')

do i = 1, dimlen
   string_var(i) = adjustl(string_var(i))
enddo 
      
if (debug > 99) then
   do i = 1, nobs
      write(msg1, *) i, ': "'//string_var(i)//'"'
      call error_handler(E_MSG, routine, msg1)
   enddo
endif

end subroutine get_string_array


!-----------------------------------------------------------------------
! Match the observation stationId string to the Network gauge string to
! determine which gauge has the location matching the stationId.
! gauge_strings : gauges available in the network file
! stations      : gauges from USGS input (timeslice) file 

function find_matching_gauge(counter) result(id)

integer, intent(in) :: counter
integer :: id  

character(len=*), parameter  :: routine = 'find_matching_gauge'  

id = find_gauge_index(stations(counter))

if (debug > 5) then 
   if (id == 0) then
      write(msg1, *) 'Unable to match station id for obs #', counter, &
                  ' "', stations(counter), '"'
      call error_handler(E_MSG, routine, msg1)
   else
      write(msg1, *) 'Identified station: "', stations(counter), '"'
      call error_handler(E_MSG, routine, msg1)
   endif
endif

end function find_matching_gauge


!-----------------------------------------------
! Use gsmd to retreive state index from location

function linkloc_to_dart(lat, lon) result (dartindx)

real(r8), intent(in) :: lat
real(r8), intent(in) :: lon
integer(i8)          :: dartindx
integer(i8)          :: sorted_index, n
logical, save        :: sorted = .false.

character(len=*), parameter  :: routine = 'linkloc_to_dart'

if (.not. sorted) then
   call create_fast_lookup_table()
   sorted = .true.
endif

! Find the first longitude that matches
dartindx = 0

LONLOOP : do n = 1,lookup_table%nsegments

   sorted_index = lookup_table%sortedindex(n)

   if (abs(lon - lookup_table%lon(sorted_index)) < .0000001 .and. &
       abs(lat - lookup_table%lat(sorted_index)) < .0000001) then
      dartindx = sorted_index
      exit LONLOOP
   endif

enddo LONLOOP

if (dartindx == 0) then
   write(msg2, *) 'longitude is ', lon
   write(msg3, *) 'latitude  is ', lat
   call error_handler(E_ERR, routine, msg1, source, text2=msg2, text3=msg3)
endif

! to make it an identity observation, it has to be negative
dartindx = -1 * dartindx

if (debug > 20) then
   write(msg1, *) 'lon, lat ', lon, lat
   write(msg2, *) 'matched at DART location ', dartindx
   call error_handler(E_MSG, routine, msg1, source, text2=msg2)
endif

end function linkloc_to_dart


!----------------------------------------------
! Build a sorted look up table; longitude-based

subroutine create_fast_lookup_table()

integer(i8)         :: indx, n
type(location_type) :: location
real(r8)            :: loc_array(3)
integer             :: var_type

character(len=*), parameter  :: routine = 'create_fast_lookup_table'

n = get_number_of_segments()

lookup_table%nsegments = n

allocate(lookup_table%lat(n), &
         lookup_table%lon(n), &
         lookup_table%dep(n), &
         lookup_table%sortedindex(n))

TABLE : do indx = 1, lookup_table%nsegments

   call get_state_meta_data(indx, location, var_type)

   ! The thought is that the segments are the first domain
   if (var_type /= QTY_STREAM_FLOW ) then
      write(msg1, *) 'Model size is ', lookup_table%nsegments, &
                     'working on index', indx
      write(msg2, *) 'var_type   is ', var_type, ' not ', QTY_STREAM_FLOW
      call error_handler(E_ERR, routine, msg1, source, text2=msg2)
   endif

   loc_array              = get_location(location)
   lookup_table%lon(indx) = loc_array(1)
   lookup_table%lat(indx) = loc_array(2)
   lookup_table%dep(indx) = loc_array(3)

enddo TABLE

! Do the index sort on lookup_table%longitude
call index_sort(lookup_table%lon, &
                lookup_table%sortedindex, &
                lookup_table%nsegments)

if (debug > 20) then
   write(*, *) 'In original, then sorted order:'
   do indx = 1, n
      write(*, *) indx, lookup_table%sortedindex(indx), &
                        lookup_table%lon(indx), &
                        lookup_table%lon(lookup_table%sortedindex(indx))
   enddo
endif

end subroutine create_fast_lookup_table


!-------------------------------------------------------------------------------
! read the character matrix from the netCDF file and parse into
! useable strings
      
subroutine convert_time_string(string, days, seconds, n)
      
character(len=*), intent(in)  :: string
integer, intent(out)          :: days, seconds
integer, optional, intent(in) :: n
      
integer         :: year, month, day, hour, minute, second 
type(time_type) :: darttime

character(len=*), parameter  :: routine = 'convert_time_strings'   

read(string, '(i4,5(1x,i2.2))') year, month, day, hour, minute, second
                       
if (debug > 20) then
   if (present(n)) then
      write(msg1, *) ' ..  read "', trim(string), " from observation ", n
   else
      write(msg1, *) ' ..  read ', trim(string)
   endif
   write(msg2, *) ' interpreted as ', year, month, day, hour, minute, second
   call error_handler(E_MSG, routine, msg1, text2=msg2)
endif

darttime = set_date(year, month, day, hour, minute, second)
call get_time(darttime, seconds, days)

end subroutine convert_time_string


!---------------------------------
! Build a gauge map for searching 

subroutine build_gauge_lookup()
  
integer :: iseg

allocate(gauge_map%gauge_keys(nseg))
allocate(gauge_map%gauge_indices(nseg))

! Fill the lookup table
do iseg = 1, nseg
   gauge_map%gauge_keys(iseg)    = adjustl(gauge_strings(iseg))
   gauge_map%gauge_indices(iseg) = iseg
end do

end subroutine build_gauge_lookup


!----------------------------------------
! Check if gauge matches incoming station

function find_gauge_index(station_id) result(idx)

character(len=*), intent(in) :: station_id
integer :: iseg, idx

idx = 0 ! Default: Not found

do iseg = 1, nseg
   if (station_id == gauge_map%gauge_keys(iseg)) then
      idx = gauge_map%gauge_indices(iseg)
      return
   end if
end do

end function find_gauge_index


! ...


end program create_identity_streamflow_obs

! <next few lines under version control, do not edit>
! $URL$
! $Id$
! $Revision$
! $Date$
