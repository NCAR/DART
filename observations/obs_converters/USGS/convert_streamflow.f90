! DART software - Copyright UCAR. This open source software is provided
! by UCAR, "as is", without charge, subject to all terms of use at
! http://www.image.ucar.edu/DAReS/DART/DART_download
!
! DART $Id$

program convert_streamflow

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! convert_streamflow - reads streamflow data and writes a DART 
!                      obs_seq file using the DART library routines.
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

use         types_mod, only : r8, missing_r8

use      location_mod, only : VERTISHEIGHT

use     utilities_mod, only : initialize_utilities, finalize_utilities, &
                              nmlfileunit, do_nml_file, do_nml_term, &
                              find_namelist_in_file, check_namelist_read, &
                              error_handler, E_ERR, E_MSG, find_textfile_dims, &
                              open_file, close_file

use  time_manager_mod, only : time_type, set_calendar_type, GREGORIAN, &
                              increment_time, set_time, get_time, set_date, &
                              operator(>=)

use  obs_sequence_mod, only : obs_sequence_type, obs_type, read_obs_seq, &
                              static_init_obs_sequence, init_obs, write_obs_seq, &
                              init_obs_sequence, get_num_obs, &
                              set_copy_meta_data, set_qc_meta_data, &
                              get_num_copies, get_num_qc

use      obs_kind_mod, only : STREAM_FLOW

use obs_def_streamflow_mod, only : set_streamflow_metadata

use obs_utilities_mod, only : getvar_real, getvar_int, get_or_fill_QC, add_obs_to_seq, &
                              create_3d_obs, getvar_int, getdimlen, set_missing_name
use    netcdf_utilities_mod, only : nc_check
use netcdf

implicit none

! version controlled file description for error handling, do not edit
character(len=*), parameter :: source   = "USGS/convert_streamflow.f90"
character(len=*), parameter :: revision = "$Revision$"
character(len=*), parameter :: revdate  = "$Date$"
character(len=*), parameter :: routine  = 'convert_streamflow:'

character(len=512) :: string1, string2, string3 ! strings for messages

! input file has data quality control fields, whether to use or ignore them.

integer,  parameter :: NUM_COPIES      = 1      ! number of copies in sequence
integer,  parameter :: NUM_QC          = 1      ! number of QC entries

integer :: existing_num_copies, existing_num_qc

! These should/must match what is in the netCDF files ... should error check
! or at least provide a decent error message.

integer, parameter :: IDLength = 15
integer, parameter :: stationIdStrLen = 15
integer, parameter :: timeStrLen = 19

integer :: num_new_obs, nobs, n, i, nlinks, indx, key
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
integer         :: n_wanted_gages
real(r8)        :: oerr, qc
integer         :: oday, osec
type(obs_type)  :: obs

! namelist variables
character(len=256) :: input_files     = 'inputs.txt'
character(len=256) :: output_file     = 'obs_seq.out'
character(len=256) :: location_file   = 'location.nc'
character(len=256) :: gages_list_file = ''
real(r8)           :: obs_fraction_for_error = 0.01_r8
real(r8)           :: obs_min_err_std = 0.5_r8
integer            :: verbose = 0

namelist / convert_streamflow_nml / &
     input_files, output_file, location_file, &
     obs_fraction_for_error, obs_min_err_std, &
     verbose, gages_list_file

!*****************************************************************************************

call initialize_utilities(routine)

! Read the DART namelist
call find_namelist_in_file('input.nml', 'convert_streamflow_nml', iunit)
read(iunit, nml = convert_streamflow_nml, iostat = io)
call check_namelist_read(iunit, io, 'convert_streamflow_nml')

! print the content of the namelist for clarification 
print*, 'list of input files is in : "'//trim(input_files)//'"'
print*, 'output_file               : "'//trim(output_file)//'"'
print*, 'location_file             : "'//trim(location_file)//'"'
print*, 'obs_fraction_for_error    : ', obs_fraction_for_error
print*, 'verbose                   : ', verbose
print*, 'gages_list_file           : "'//trim(gages_list_file)//'"'

! Record the DART namelist values used for the run ...
if (do_nml_file()) write(nmlfileunit, nml=convert_streamflow_nml)
if (do_nml_term()) write(     *     , nml=convert_streamflow_nml)

! put the reference date into DART format
call set_calendar_type(GREGORIAN)

!*****************************************************************************************
! READ the file with the location information 
!*****************************************************************************************

io =  nf90_open(location_file, nf90_nowrite, ncid)
call nc_check(io, routine, 'opening location file "'//trim(location_file)//'"' )

call getdimlen(ncid, 'linkDim', nlinks)

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

!*****************************************************************************************

call static_init_obs_sequence()
call init_obs(obs,      num_copies=NUM_COPIES, num_qc=NUM_QC)
call init_obs(prev_obs, num_copies=NUM_COPIES, num_qc=NUM_QC)

n_wanted_gages = set_desired_gages(gages_list_file)
call find_textfile_dims(input_files, nfiles)

num_new_obs = estimate_total_obs_count(input_files, nfiles)

inquire(file=output_file, exist=file_exist)

if ( file_exist ) then ! existing file found, append to it

   write(string1,*) "found existing obs_seq file, appending to ", trim(output_file)
   write(string2,*) "adding up to a maximum of ", num_new_obs, " new observations"
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

!*****************************************************************************************
!   Loop through the time slices and adding the data to obs_seq.out 
!   if the gage is in the gage_strings
!*****************************************************************************************

iunit = open_file(input_files,form='formatted',action='read')

first_obs = .true.

FILELOOP : do ifile=1,nfiles

   read(iunit,'(A)', iostat=io) input_file
   if (io /= 0 ) then 
     write(string1,*) 'Unable to read input file from "'//trim(input_files)//'"'
     write(string2,*) 'file ',ifile
     call error_handler(E_ERR,'convert_streamflow',string1, &
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

      if ( discharge(n) < 0.0_r8 ) cycle OBSLOOP

      ! relate the TimeSlice:station to the RouteLink:gage
      indx = find_matching_gage_index(n)
      if (indx == 0) cycle OBSLOOP

      ! oerr is the observation error standard deviation in this application.
      ! The observation error variance encoded in the observation file
      ! will be oerr*oerr
      oerr = max(discharge(n)*obs_fraction_for_error, obs_min_err_std)

      call convert_time_string(time_string(n),oday,osec,n)
      time_obs = set_time(osec,oday)  ! yes, seconds then days

      call set_streamflow_metadata(key, station_strings(n), link(indx))

      call create_3d_obs(lat(indx), lon(indx), altitude(indx), VERTISHEIGHT, &
                         discharge(n), STREAM_FLOW, oerr, oday, osec, qc, obs, key)
      call add_obs_to_seq(obs_seq, obs, time_obs, prev_obs, prev_time, first_obs)
 
    enddo OBSLOOP

   deallocate(discharge, discharge_quality)
   deallocate(time_string, station_strings)

enddo FILELOOP

nobs = get_num_obs(obs_seq)

! If we added any obs to the sequence, write it now.
if ( nobs > 0 )  call write_obs_seq(obs_seq, output_file)

if (verbose > 0) then
   write(string1,*)'Writing ',nobs,' observations to file.'
   call error_handler(E_MSG, routine, string1)
endif

call close_file(iunit)

deallocate(lat, lon, altitude, link, gage_strings)
if (allocated(desired_gages)) deallocate(desired_gages)
call finalize_utilities()

! end of main program

contains

!-----------------------------------------------------------------------
!> Match the observation stationId string to the RouteLink gage string to 
!> determine which gage has the location matching the stationId.
!>
!> If the gage is not one we desire, set the id to 0 to just skip it.

function find_matching_gage_index(counter) result(id)

integer, intent(in) :: counter
integer :: id

integer :: i

id = 0  ! Indicate the station has no matching gage or is not wanted.

LINKS : do i = 1,nlinks
   if (station_strings(counter) == gage_strings(i)) then
      if (verbose > 0) then
         write(string1,*)'identified station "',station_strings(counter),'"'
         call error_handler(E_MSG, 'find_matching_gage_index:', string1)
      endif
      id = i
      exit LINKS
   endif
enddo LINKS

if (id == 0 .and. verbose > 1) then
   write(string1,*)'Unable to match station id for obs #',counter,' "',station_strings(counter),'"'
   call error_handler(E_MSG, 'find_matching_gage_index', string1, source, revision, revdate)
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


!-----------------------------------------------------------------------
!> read the character matrix from the netCDF file and parse into
!> useable strings


subroutine convert_time_string(string,days,seconds,n)

character(len=*), intent(in) :: string
integer, intent(out) :: days, seconds
integer, optional, intent(in) :: n

integer :: year, month, day, hour, minute, second
type(time_type) :: darttime

read(string,'(i4,5(1x,i2.2))') year, month, day, hour, minute, second

if (verbose > 1) then
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


!-----------------------------------------------------------------------
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

if (verbose > 2) then
   do i = 1,nobs
      write(string1,*) 'time ',i,' is "'//time_string(i)//'"'
      call error_handler(E_MSG, 'get_time_strings:', string1)
   enddo
endif

end subroutine  get_time_strings


!-----------------------------------------------------------------------
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

if (verbose > 2) then
   do i = 1,nobs
      write(string1,*) 'station ',i,' is "'//station_strings(i)//'"'
      call error_handler(E_MSG, 'get_station_strings:', string1)
   enddo
endif

end subroutine  get_station_strings


!-----------------------------------------------------------------------
!> read the character matrix of NHD Gage Event IDs  from the metadata
!> netCDF file and parse into an array of strings.
!> dimensions:
!>        linkDim = 157 ;
!>        IDLength = 15 ;
!> variables:
!>        char gages(linkDim, IDLength) ;
!>             gages:long_name = "NHD Gage Event ID from SOURCE_FEA field in Gages feature class" ;
!>             gages:coordinates = "lat lon" ;

subroutine get_gage_strings(ncid, varname, gagelist)

integer,          intent(in)  :: ncid
character(len=*), intent(in)  :: varname
character(len=*), intent(out) :: gagelist(:)

integer :: io, dim1

call getdimlen(ncid, 'linkDim', dim1)

io = nf90_inq_varid(ncid,varname,varid)
call nc_check(io, 'get_gage_strings', 'inq_varid "'//varname//'"')

io = nf90_get_var(ncid,varid,gagelist)
call nc_check(io, 'get_gage_strings', 'get_var "'//varname//'"')

do i = 1,nlinks
   gagelist(i) = adjustl(gagelist(i))
enddo

if (verbose > 2) then
   do i = 1,nlinks
      write(string1,*) 'gage_strings ',i,' is "'//gagelist(i)//'"'
      call error_handler(E_MSG, 'get_gage_strings:', string1)
   enddo
endif

end subroutine get_gage_strings


!-----------------------------------------------------------------------
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

if (verbose > 1) then
   print*, 'List of the gages to process:'
   do i = 1,nlines
      print*, i, desired_gages(i)
   enddo
endif

ngages = nlines

end function set_desired_gages


!-----------------------------------------------------------------------
!> 

function estimate_total_obs_count(file_list,nfiles) result (num_obs)

character(len=*), intent(in) :: file_list
integer,          intent(in) :: nfiles
integer                      :: num_obs

! Rather than opening up each file and getting the actual number of possible
! observations in each file, we are just going to extrapolate from a single
! file and add a fudge factor. There does not seem to be much variability in
! the number of possible observations from file-to-file.

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
! I am adding 20% 

num_obs = 1.2_r8 * nobs * nfiles

end function estimate_total_obs_count


end program convert_streamflow

! <next few lines under version control, do not edit>
! $URL$
! $Id$
! $Revision$
! $Date$
