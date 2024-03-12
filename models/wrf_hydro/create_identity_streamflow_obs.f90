! DART software - Copyright UCAR. This open source software is provided
! by UCAR, "as is", without charge, subject to all terms of use at
! http://www.image.ucar.edu/DAReS/DART/DART_download
!
! DART $Id$

program create_identity_streamflow_obs

!-------------------------------------------------------------------------------
!
! create_identity_streamflow_obs - reads streamflow data and a
!     channel-only model and writes a DART obs_seq file where the
!     observations are actually identity observations.
!
!-------------------------------------------------------------------------------

use         types_mod, only : r8, missing_r8, i8

use      location_mod, only : VERTISHEIGHT, location_type, get_location

use     utilities_mod, only : nmlfileunit, do_nml_file, do_nml_term, &
                              initialize_utilities, finalize_utilities, &
                              find_namelist_in_file, check_namelist_read, &
                              error_handler, E_ERR, E_MSG, &
                              find_textfile_dims, &
                              open_file, close_file

use  time_manager_mod, only : time_type, set_calendar_type, GREGORIAN, &
                              increment_time, set_time, get_time, set_date, &
                              operator(>=)

use  obs_sequence_mod, only : obs_sequence_type, obs_type, read_obs_seq, &
                              static_init_obs_sequence, write_obs_seq, &
                              init_obs, init_obs_sequence, get_num_obs, &
                              set_copy_meta_data, set_qc_meta_data, &
                              get_num_copies, get_num_qc

use      obs_kind_mod, only : QTY_STREAM_FLOW

use obs_def_streamflow_mod, only : set_streamflow_metadata

use obs_utilities_mod, only : getvar_real, getvar_int, get_or_fill_QC, &
                              add_obs_to_seq, create_3d_obs, getvar_int, &
                              getdimlen, set_missing_name

use         model_mod, only : static_init_model, get_state_meta_data, &
                              get_number_of_links

use          sort_mod, only : index_sort
use    netcdf_utilities_mod, only : nc_check
use netcdf

implicit none

! version controlled file description for error handling, do not edit
character(len=*), parameter :: source   = &
   "$URL$"
character(len=*), parameter :: revision = "$Revision$"
character(len=*), parameter :: revdate  = "$Date$"
character(len=*), parameter :: routine  = 'create_identity_streamflow_obs'

character(len=512) :: string1, string2, string3 ! strings for messages

! input file has data quality control fields, whether to use or ignore them.

integer,  parameter :: NUM_COPIES      = 1      ! number of copies in sequence
integer,  parameter :: NUM_QC          = 1      ! number of QC entries
real(r8), parameter :: MIN_OBS_ERR_STD = 0.1_r8 ! m^3/sec
real(r8), parameter :: MAX_OBS_ERR_STD = 1000000.0_r8 
real(r8), parameter :: NORMAL_FLOW     = 10.0_r8
real(r8), parameter :: contract        = 0.001_r8

integer :: existing_num_copies, existing_num_qc

! These should/must match what is in the netCDF files ... should error check
! or at least provide a decent error message.

integer, parameter :: IDLength = 15
integer, parameter :: stationIdStrLen = 15
integer, parameter :: timeStrLen = 19

integer :: num_new_obs, nobs, n, i, nlinks, indx, key
integer :: dart_index                 ! matching state vector index
integer :: ifile, iunit, io, nfiles
integer :: ncid, varid
logical :: file_exist, first_obs
character(len=256) :: input_file

type(obs_sequence_type) :: obs_seq
type(obs_type)          :: prev_obs
type(time_type)         :: time_obs, prev_time

! variables read from RouteLink file
real(r8),                allocatable :: lat(:)
real(r8),                allocatable :: lon(:)
integer,                 allocatable :: link(:)
real(r8),                allocatable :: altitude(:)
character(len=IDLength), allocatable :: gage_strings(:)

! variables read from observation data file
character(len=timeStrLen),      allocatable :: time_string(:)
character(len=stationIdStrLen), allocatable :: station_strings(:)
integer,                        allocatable :: discharge_quality(:)
real(r8),                       allocatable :: discharge(:)

character(len=IDLength),        allocatable :: desired_gages(:)
integer         :: n_wanted_gages, n_desired_gages
real(r8)        :: oerr, qc
integer         :: oday, osec
type(obs_type)  :: obs

! structure to hold the lat/lon/depth/indexing information
! for the DART state vector
type database
   integer(i8)              :: nlinks
   real(r8),    allocatable :: latitude(:)
   real(r8),    allocatable :: longitude(:)
   real(r8),    allocatable :: depth(:)
   integer(i8), allocatable :: sortedindex(:)
end type database
type(database) :: lookup_table


! namelist variables
character(len=256) :: input_files     = 'inputs.txt'
character(len=256) :: output_file     = 'obs_seq.out'
character(len=256) :: location_file   = 'location.nc'
character(len=256) :: gages_list_file = ''
real(r8)           :: obs_fraction_for_error = 0.01
logical            :: assimilate_all  = .false. 
integer            :: debug = 0

namelist / create_identity_streamflow_obs_nml / &
               input_files, &
               output_file, &
               location_file, &
               gages_list_file, &
               obs_fraction_for_error, &
               assimilate_all, &
               debug

!-------------------------------------------------------------------------------

call initialize_utilities(routine)

! Read the DART namelist
call find_namelist_in_file('input.nml', 'create_identity_streamflow_obs_nml', iunit)
read(iunit, nml = create_identity_streamflow_obs_nml, iostat = io)
call check_namelist_read(iunit, io, 'create_identity_streamflow_obs_nml')

! print the content of the namelist for clarification
print*, 'list of input files is in : "'//trim(input_files)//'"'
print*, 'output_file               : "'//trim(output_file)//'"'
print*, 'location_file             : "'//trim(location_file)//'"'
print*, 'obs_fraction_for_error    : ', obs_fraction_for_error
print*, 'debug                     : ', debug
print*, 'gages_list_file           : "'//trim(gages_list_file)//'"'

! Record the DART namelist values used for the run ...
if (do_nml_file()) write(nmlfileunit, nml=create_identity_streamflow_obs_nml)
if (do_nml_term()) write(     *     , nml=create_identity_streamflow_obs_nml)

! put the reference date into DART format
call set_calendar_type(GREGORIAN)

! Initialize the DART state vector so we can relate the observation location
! to the location in the DART state vector.

call static_init_model()

!-------------------------------------------------------------------------------
! READ the file with the location information
!-------------------------------------------------------------------------------

io =  nf90_open(location_file, nf90_nowrite, ncid)
call nc_check(io, routine, 'opening location file "'//trim(location_file)//'"' )

call getdimlen(ncid, 'feature_id', nlinks)

write(string1,*)'number of links is ',nlinks
call error_handler(E_MSG, routine, string1)

allocate(         lat(nlinks))
allocate(         lon(nlinks))
allocate(    altitude(nlinks))
allocate(        link(nlinks))
allocate(gage_strings(nlinks))

! read in the data arrays from the RouteLink file
call getvar_real(     ncid, 'lat',   lat       ) ! latitudes
call getvar_real(     ncid, 'lon',   lon       ) ! longitudes
call getvar_real(     ncid, 'alt',   altitude  ) ! elevation
call getvar_int(      ncid, 'link',  link      ) ! "Link ID (NHDFlowline_network COMID)"
call get_gage_strings(ncid, 'gages', gage_strings) ! "NHD Gage Event ID from SOURCE_FEA ... "

! convert all lons to [0,360], replace any bad lat lims
! check the lat/lon values to see if they are ok
! change lon from -180 to 180 into 0-360

where (lat >  90.0_r8) lat =  90.0_r8
where (lat < -90.0_r8) lat = -90.0_r8
where (lon <   0.0_r8) lon = lon + 360.0_r8

call nc_check(nf90_close(ncid), routine, 'closing file '//trim(location_file))

!-------------------------------------------------------------------------------
! prepare output observation sequence - append or creating new ...
!-------------------------------------------------------------------------------

call static_init_obs_sequence()
call init_obs(obs,      num_copies=NUM_COPIES, num_qc=NUM_QC)
call init_obs(prev_obs, num_copies=NUM_COPIES, num_qc=NUM_QC)

! Collect all the gauges: 
! - desired ones will have the provided obs_err_sd
! - remaining gauges are dummy with very large obs_err_sd

n_desired_gages = set_desired_gages(gages_list_file)
n_wanted_gages  = 0 !set_desired_gages(gages_list_file)
call find_textfile_dims(input_files, nfiles)

num_new_obs = estimate_total_obs_count(input_files, nfiles)

inquire(file=output_file, exist=file_exist)

if ( file_exist ) then ! existing file found, append to it

   write(string1,*) "found existing obs_seq file, appending to ", &
                    trim(output_file)
   write(string2,*) "adding up to a maximum of ", num_new_obs, &
                    " new observations"
   call error_handler(E_MSG, routine, string1, &
                      source, revision, revdate, text2=string2)

   call read_obs_seq(output_file, 0, 0, num_new_obs, obs_seq)

   ! check to see if existing file is compatible
   existing_num_copies = get_num_copies(obs_seq)
   existing_num_qc     = get_num_qc(obs_seq)

   if (existing_num_copies /= NUM_COPIES .or.  existing_num_qc /= NUM_QC) then
      write(string1,*)'incompatible existing observation sequence file'
      write(string2,'(A,i4,A,i4)')'expected ',NUM_COPIES, &
                           ' obs copies got ',existing_num_copies
      write(string3,'(A,i4,A,i4)')'expected ',NUM_QC, &
                           ' QC  copies got ',existing_num_qc
      call error_handler(E_ERR, routine, string1, &
                  source, revision, revdate, text2=string2, text3=string3)
   endif

else ! create a new one ...

   call init_obs_sequence(obs_seq, NUM_COPIES, NUM_QC, num_new_obs)

   do i=1,NUM_COPIES ! kinda silly ... only 1 type of observation
      call set_copy_meta_data(obs_seq, i, 'observation')
   enddo
   do i=1,NUM_QC ! kinda silly ... only 1 type of qc
      call set_qc_meta_data(obs_seq, i, 'Data QC')
   enddo

   write(string1,*) "no existing obs_seq file, creating ", trim(output_file)
   write(string2,*) "with up to a maximum of ", num_new_obs, " observations"
   call error_handler(E_MSG, routine, string1, &
                      source, revision, revdate, text2=string2)

endif

!-------------------------------------------------------------------------------
! Loop through the time slices and adding the data to obs_seq.out
! if the gage is in the gage_strings
!-------------------------------------------------------------------------------

iunit = open_file(input_files,form='formatted',action='read')

first_obs = .true.

FILELOOP : do ifile=1,nfiles

   read(iunit,'(A)', iostat=io) input_file
   if (io /= 0 ) then
     write(string1,*) 'Unable to read input file from "'//trim(input_files)//'"'
     write(string2,*) 'file ',ifile
     call error_handler(E_ERR,'create_identity_streamflow_obs',string1, &
                source, revision, revdate, text2=string2)
   endif

   io = nf90_open(input_file, nf90_nowrite, ncid)
   call nc_check(io, routine, 'opening data file "'//trim(input_file)//'"' )
   call getdimlen(ncid, 'stationIdInd', nobs)

   write(string1,*)'Reading time slice file "',trim(input_file)//'"'
   write(string2,*)'number of obs in the time slice is ',nobs
   call error_handler(E_MSG, routine, string1, text2=string2)

   allocate(         discharge(nobs))
   allocate( discharge_quality(nobs))
   allocate(   station_strings(nobs))
   allocate(       time_string(nobs))

   ! read in the data arrays
   call getvar_real(ncid, 'discharge', discharge )  ! streamflow
   call getvar_int( ncid, 'discharge_quality', discharge_quality)

   call get_station_strings(ncid, 'stationId' ) ! Name of the gage
   call get_time_strings(   ncid, 'time'      ) ! observation time

   io = nf90_close(ncid)
   call nc_check(io, routine, 'closing file '//trim(input_file) )

   ! Set the DART data quality control.  Be consistent with NCEP codes;
   ! 0 is 'must use', 1 is good, no reason not to use it.

   qc = 1.0_r8   ! modify based on discharge_quality ... perhaps

   OBSLOOP: do n = 1, nobs

      ! make sure discharge is physical 
      if ( discharge(n) < 0.0_r8 .or. discharge(n) /= discharge(n) ) cycle OBSLOOP

      ! relate the TimeSlice:station to the RouteLink:gage so we can
      ! determine the location
      indx = find_matching_gage_index(n)
      if (indx == 0) cycle OBSLOOP

      ! relate the physical location to the dart state vector index
      dart_index = linkloc_to_dart(lat(indx), lon(indx))

      ! desired gauges get the provided obs_err
      ! remaining ones are for verification purposes
      if (ANY(desired_gages == station_strings(n)) .or. assimilate_all) then  
        oerr = max(discharge(n)*obs_fraction_for_error, MIN_OBS_ERR_STD)
      else 
        oerr = MAX_OBS_ERR_STD
      endif
         ! don't correct that much, the gauge observations imply that the flow 
         ! in the stream is small. This is not a flood period. Streamflow values
         ! indicate a more or less lake situation rather than a strongly flowing stream. 
         ! For this, choose a large value for the observation error standard deviation 
         ! in order not to crush the ensemble spread. 
      !   oerr = MAX_OBS_ERR_STD
      !else 
         ! This is a more interesting scenario where the flow in teh stream 
         ! is big enough for DA to make sense. 
         
         ! NEW MOHA
      !   oerr = MIN_OBS_ERR_STD + contract*(discharge(n) - NORMAL_FLOW)**2
      !   oerr = max(discharge(n)*obs_fraction_for_error, MIN_OBS_ERR_STD)
      !endif 

      call convert_time_string(time_string(n),oday,osec,n)
      time_obs = set_time(osec,oday)  ! yes, seconds then days

      call set_streamflow_metadata(key, station_strings(n), link(indx))

      call create_3d_obs(lat(indx), lon(indx), altitude(indx), VERTISHEIGHT, &
                  discharge(n), dart_index, oerr, oday, osec, qc, obs, key)
      call add_obs_to_seq(obs_seq, obs, time_obs, prev_obs, prev_time, first_obs)

    enddo OBSLOOP

   deallocate(discharge, discharge_quality)
   deallocate(time_string, station_strings)

enddo FILELOOP

nobs = get_num_obs(obs_seq)

! If we added any obs to the sequence, write it now.
if ( nobs > 0 )  call write_obs_seq(obs_seq, output_file)

if (debug > 0) then
   write(string1,*)'Writing ',nobs,' observations to file.'
   call error_handler(E_MSG, routine, string1)
endif

call close_file(iunit)

deallocate(lat, lon, altitude, link, gage_strings)
if (allocated(desired_gages)) deallocate(desired_gages)
call finalize_utilities()

! end of main program

contains

!-------------------------------------------------------------------------------
!> Match the observation stationId string to the RouteLink gage string to
!> determine which gage has the location matching the stationId.
!>
!> If the gage is not one we desire, set the id to 0 to just skip it.
!>
!>@todo There are probably faster ways to do this with index sorting, etc.
!>      as the list of observations or the number of links increases, this
!>      will become more important.

function find_matching_gage_index(counter) result(id)

integer, intent(in) :: counter
integer :: id

integer :: i

id = 0  ! Indicate the station has no matching gage or is not wanted.

LINKS : do i = 1,nlinks
   if (station_strings(counter) == gage_strings(i)) then
      if (debug > 0) then
         write(string1,*)'identified station "',station_strings(counter),'"'
         call error_handler(E_MSG, 'find_matching_gage_index:', string1)
      endif
      id = i
      exit LINKS
   endif
enddo LINKS

if (id == 0 .and. debug > 1) then
   write(string1,*)'Unable to match station id for obs #',counter, &
                   ' "',station_strings(counter),'"'
   call error_handler(E_MSG, 'find_matching_gage_index', string1, &
              source, revision, revdate)
endif

if (id == 0) return  ! unable to find location information

! If no subsetting gage list is specified, we want them all.
! If the gage is desired, we want it.
! Just return with whatever id we have already found.
if ( n_wanted_gages == 0 ) then
   continue
else if( ANY(desired_gages == station_strings(counter)) ) then
   continue
else
   ! We do not want the gage
   id = 0
endif

end function find_matching_gage_index


!-------------------------------------------------------------------------------
!> read the character matrix from the netCDF file and parse into
!> useable strings


subroutine convert_time_string(string,days,seconds,n)

character(len=*), intent(in) :: string
integer, intent(out) :: days, seconds
integer, optional, intent(in) :: n

integer :: year, month, day, hour, minute, second
type(time_type) :: darttime

read(string,'(i4,5(1x,i2.2))') year, month, day, hour, minute, second

if (debug > 1) then
   if ( present(n) ) then
      write(string1,*) ' ..  read "',trim(string)," from observation ",n
   else
      write(string1,*) ' ..  read ',trim(string)
   endif
   write(string2,*)' interpreted as ', year, month, day, hour, minute, second
   call error_handler(E_MSG,'convert_time_string:',string1,text2=string2)
endif

darttime = set_date(year, month, day, hour, minute, second)
call get_time(darttime, seconds, days)

end subroutine convert_time_string


!-------------------------------------------------------------------------------
!> read the character matrix of UTC times from the data netCDF file and
!> parse into an array of strings.
!> dimensions:
!>         stationIdInd = UNLIMITED ; // (2 currently)
!>         timeStrLen = 19 ;
!> variables:
!>         char time(stationIdInd, timeStrLen) ;
!>              time:units = "UTC" ;
!>              time:long_name = "YYYY-MM-DD_HH:mm:ss UTC" ;

subroutine  get_time_strings(ncid, varname)

integer,          intent(in) :: ncid
character(len=*), intent(in) :: varname

integer :: io, dim1

call getdimlen(ncid, 'timeStrLen', dim1)

io = nf90_inq_varid(ncid,varname,varid)
call nc_check(io, 'get_time_strings', 'inq_varid "'//varname//'"')

io = nf90_get_var(ncid,varid,time_string)
call nc_check(io, 'get_time_strings', 'get_var "'//varname//'"')

if (debug > 2) then
   do i = 1,nobs
      write(string1,*) 'time ',i,' is "'//time_string(i)//'"'
      call error_handler(E_MSG, 'get_time_strings:', string1)
   enddo
endif

end subroutine  get_time_strings


!-------------------------------------------------------------------------------
!> read the character matrix of USGS station identifiers from the data
!> netCDF file and parse into an array of strings.
!> dimensions:
!>         stationIdInd = UNLIMITED ; // (2 currently)
!>         stationIdStrLen = 15 ;
!> variables:
!>         char stationId(stationIdInd, stationIdStrLen) ;
!>              stationId:long_name = "USGS station identifer of length 15" ;


subroutine  get_station_strings(ncid, varname)

integer,          intent(in) :: ncid
character(len=*), intent(in) :: varname

integer :: io, dim1

call getdimlen(ncid, 'stationIdStrLen', dim1)

io = nf90_inq_varid(ncid,varname,varid)
call nc_check(io, 'get_station_strings', 'inq_varid "'//varname//'"')

io = nf90_get_var(ncid,varid,station_strings)
call nc_check(io, 'get_station_strings', 'get_var "'//varname//'"')

do i = 1,nobs
   station_strings(i) = adjustl(station_strings(i))
enddo

if (debug > 2) then
   do i = 1,nobs
      write(string1,*) 'station ',i,' is "'//station_strings(i)//'"'
      call error_handler(E_MSG, 'get_station_strings:', string1)
   enddo
endif

end subroutine  get_station_strings


!-------------------------------------------------------------------------------
!> read the character matrix of NHD Gage Event IDs  from the metadata
!> netCDF file and parse into an array of strings.
!> dimensions:
!>     feature_id = 157 ;
!>     IDLength = 15 ;
!> variables:
!>     char gages(feature_id, IDLength) ;
!>          gages:long_name = "NHD Gage Event ID from SOURCE_FEA field in Gages feature class" ;
!>          gages:coordinates = "lat lon" ;

subroutine get_gage_strings(ncid, varname, gagelist)

integer,          intent(in)  :: ncid
character(len=*), intent(in)  :: varname
character(len=*), intent(out) :: gagelist(:)

integer :: io, dim1

call getdimlen(ncid, 'feature_id', dim1)

io = nf90_inq_varid(ncid,varname,varid)
call nc_check(io, 'get_gage_strings', 'inq_varid "'//varname//'"')

io = nf90_get_var(ncid,varid,gagelist)
call nc_check(io, 'get_gage_strings', 'get_var "'//varname//'"')

do i = 1,nlinks
   gagelist(i) = adjustl(gagelist(i))
enddo

if (debug > 2) then
   do i = 1,nlinks
      write(string1,*) 'gage_strings ',i,' is "'//gagelist(i)//'"'
      call error_handler(E_MSG, 'get_gage_strings:', string1)
   enddo
endif

end subroutine get_gage_strings


!-------------------------------------------------------------------------------
!> Read the list of the gages to process.
!> If the filename is empty, use all the gages.
!> the list of desired gages is globally scoped:  desired_gages()

function set_desired_gages(filename) result(ngages)
character(len=*), intent(in) :: filename
integer                      :: ngages

character(len=*), parameter :: routine = 'set_desired_gages'
integer :: nlines, i

ngages = 0

if (filename == '' .or. filename == 'null') then
   write(string1,*)'No gages_list_file supplied. Using all the gages.'
   call error_handler(E_MSG,routine,string1)
   return
endif

call find_textfile_dims(filename, nlines)

write(string1,*) 'number of gages to read from the time slices :', nlines
call error_handler(E_MSG,routine,string1)

allocate(desired_gages(nlines))

iunit = open_file(filename,form='formatted',action='read')
do i = 1, nlines
   read(iunit,*,iostat=io)  desired_gages(i)
   if (io /= 0) then
      write(string1,*) 'Unable to read gage from "'//trim(filename)//'"'
      write(string2,*) 'line ',i
      call error_handler(E_ERR,routine,string1, &
                 source, revision, revdate, text2=string2)
   endif
enddo
call close_file(iunit)

if (debug > 1) then
   print*, 'List of the gages to process:'
   do i = 1,nlines
      print*, i, desired_gages(i)
   enddo
endif

ngages = nlines

end function set_desired_gages


!-------------------------------------------------------------------------------
!>

function estimate_total_obs_count(file_list,nfiles) result (num_obs)

character(len=*), intent(in) :: file_list
integer,          intent(in) :: nfiles
integer                      :: num_obs

character(len=*), parameter :: routine = 'estimate_total_obs_count'
integer :: iunit, io, ncid, nobs
character(len=256) :: input_file

iunit = open_file(file_list,form='formatted',action='read')
read(iunit,'(A)', iostat=io) input_file
if (io /= 0 ) then
  write(string1,*) 'Unable to read input file from "'//trim(file_list)//'"'
  call error_handler(E_ERR,routine,string1,source,revision,revdate)
endif
call close_file(iunit)

! Need to know about how many observations are in each file.
io = nf90_open(input_file, nf90_nowrite, ncid)
call nc_check(io, routine, 'opening data file "'//trim(input_file)//'"' )
call getdimlen(ncid, 'stationIdInd', nobs)
io = nf90_close(ncid)
call nc_check(io, routine, 'closing file "'//trim(input_file)//'"' )

! We need to know how many observations there may be.
! Specifying too many is not really a problem.
! I am multiplying by 10.

num_obs = 10.0_r8 * nobs * nfiles

end function estimate_total_obs_count


!-------------------------------------------------------------------------------
!>

function linkloc_to_dart(lat, lon) result (dartindx)

real(r8), intent(in) :: lat
real(r8), intent(in) :: lon
integer(i8)          :: dartindx

integer(i8) :: sorted_index, n

logical, save :: sorted = .false.

if ( .not. sorted ) then
   call create_fast_lookup_table()
   sorted = .true.
endif

! Find the first longitude that matches

dartindx = 0
LONLOOP : do n = 1,lookup_table%nlinks

   sorted_index = lookup_table%sortedindex(n)

   if ( abs(lon - lookup_table%longitude(sorted_index)) < .0000001 .and. &
        abs(lat - lookup_table%latitude( sorted_index)) < .0000001      ) then
      dartindx = sorted_index
      exit LONLOOP
   endif

enddo LONLOOP

if (dartindx == 0) then
   write(string2,*)'longitude is ',lon
   write(string3,*)'latitude  is ',lat
   call error_handler(E_ERR,'linkloc_to_dart','no matching location', &
              source, revision, revdate, text2=string2, text3=string3)
endif

! to make it an identity observation, it has to be negative

dartindx = -1 * dartindx

if (debug > 99) then
   write(string1,*)'lon, lat ',lon,lat
   write(string2,*)'matched at DART location ',dartindx
   call error_handler(E_MSG,'linkloc_to_dart',string1, &
              source, revision, revdate, text2=string2)
endif

end function linkloc_to_dart


!-------------------------------------------------------------------------------
!>

subroutine create_fast_lookup_table()

integer(i8)         :: indx, n
type(location_type) :: location
real(r8)            :: loc_array(3)
integer             :: var_type

n = get_number_of_links()

lookup_table%nlinks = n

allocate(lookup_table%latitude(n), &
         lookup_table%longitude(n), &
         lookup_table%depth(n), &
         lookup_table%sortedindex(n))

TABLE : do indx = 1,lookup_table%nlinks

   call get_state_meta_data(indx, location, var_type)

   ! The thought is that the links are the first domain
   if (var_type /= QTY_STREAM_FLOW ) then
      write(string1,*)'model size is ',lookup_table%nlinks, &
                      'working on index', indx
      write(string2,*)'var_type   is ',var_type, ' not ', QTY_STREAM_FLOW
      call error_handler(E_ERR, 'create_fast_lookup_table', string1, &
                 source, revision, revdate, text2=string2)
   endif

   loc_array                    = get_location(location)
   lookup_table%longitude(indx) = loc_array(1)
   lookup_table%latitude( indx) = loc_array(2)
   lookup_table%depth(    indx) = loc_array(3)

enddo TABLE

!  do the index sort on lookup_table%longitude
call index_sort(lookup_table%longitude, &
                lookup_table%sortedindex, &
                lookup_table%nlinks )

if (debug > 99) then
   write(*,*)'In original, then sorted order:'
   do indx = 1,n
      write(*,*) indx, lookup_table%sortedindex(indx), &
                       lookup_table%longitude(indx), &
                       lookup_table%longitude(lookup_table%sortedindex(indx))
   enddo
endif

end subroutine create_fast_lookup_table


end program create_identity_streamflow_obs

! <next few lines under version control, do not edit>
! $URL$
! $Id$
! $Revision$
! $Date$
