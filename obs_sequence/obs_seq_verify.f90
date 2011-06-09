! DART software - Copyright 2004 - 2011 UCAR. This open source software is
! provided by UCAR, "as is", without charge, subject to all terms of use at
! http://www.image.ucar.edu/DAReS/DART/DART_download

program obs_seq_verify

! <next few lines under version control, do not edit>
! $URL$
! $Id$
! $Revision$
! $Date$

!-----------------------------------------------------------------------
!
! This program creates a netCDF file suitable for forecast evaluation.
!
! (copy, station, level, ensemble, analysisT, fcstlead)
!   |       |       |       |          |        |
!   |       |       |       |          |        +- forecast length : 0,3,6,9,12,...
!   |       |       |       |          |
!   |       |       |       |          +---------- analysis time/date
!   |       |       |       |
!   |       |       |       +--------------------- ensemble member index
!   |       |       |
!   |       |       +----------------------------- vertical level index
!   |       |
!   |       +------------------------------------- (horizontal) station index
!   |
!   +--------------------------------------------- obs value, prior, obs_err.
!
! I think the logic of the program should be as follows:
! 1) The list of 'stations' is read from a netCDF file -- the product of obs_seq_coverage.f90 
! 2) The set of forecast lead times must be determined. Having 8 might not be
!    specific enough ... we want 8 separated by 3 hours ... for example.
! 3) A series of obs_seq.fcst files (each from one analysisT and containing multiple
!    forecast lead times) is read and stuffed into an appropriate structure.
! 4) The structure is written into the netCDF file.
! 5) On to the next obs_seq.fcst file ... (step 3)
! 6) wrap up ... 
! 
! Soyoung's wish list:
! Now I'm done running filter for 24-hr forecast with 3-hrly observations as an
! evaluation mode only, and am ready to hand obs_seq.final over to you for the final
! conversion process.
! 
! be1005en.ucar.edu:/ptmp/syha/wrfruns/ENS_FCST/Verify>
! -rw-r--r--    1 syha     ncar     1006513852 Nov 19 12:24 Prior_Diag.nc
! -rw-r--r--    1 syha     ncar     1006513856 Nov 19 12:24 Posterior_Diag.nc
! -rw-r--r--    1 syha     ncar      111484368 Nov 19 12:24 prior_inflate_restart
! -rw-r--r--    1 syha     ncar      126848308 Nov 19 12:24 obs_seq.final
! 
! Ideally, in obs_seq_fcst.nc (if I can name it on my own), I would like to have a
! data structure of (copy, station, level, ensemble, date, time) for each variable
! and each obs type, where copy is (observation value, prior observation value
! corresponding to the observation, obs error standard deviation), date is the
! number of cycle, and time is the number of forecast lead times.
! So, the dimension "date" is 1 for my current single obs_seq.final, but reserved
! for multiple initialization cases (i.e., multiple obs_seq.final files) with the
! same length of forecast lead times. This dimension is the most critical one in the
! statistical verification since this dimension number is the actual sample size to
! tell if we have enough samples to make the verification result statistically
! significant or not.
!
!-----------------------------------------------------------------------

use        types_mod, only : r4, r8, digits12, MISSING_R8, MISSING_I, &
                             metadatalength, obstypelength
use obs_sequence_mod, only : read_obs_seq, obs_type, obs_sequence_type, get_first_obs, &
                             get_obs_def, get_copy_meta_data, get_obs_values, &
                             get_next_obs, init_obs, init_obs_sequence, &
                             assignment(=), get_num_copies, get_num_qc, get_qc, &
                             static_init_obs_sequence, destroy_obs_sequence, destroy_obs, &
                             read_obs_seq_header, get_qc_meta_data
use      obs_def_mod, only : obs_def_type, get_obs_def_time, get_obs_kind, write_obs_def, &
                             get_obs_def_location, set_obs_def_time, &
                             set_obs_def_location, set_obs_def_kind, &
                             set_obs_def_error_variance, get_obs_def_error_variance
use     obs_kind_mod, only : max_obs_kinds, get_obs_kind_name, get_obs_kind_index, &
                             write_obs_kind
use     location_mod, only : location_type, get_location, set_location_missing, &
                             write_location, operator(/=), operator(==), &
                             set_location, is_location_in_region, query_location, &
                             nc_write_location_atts, nc_get_location_varids, &
                             nc_write_location, LocationDims
use time_manager_mod, only : time_type, set_date, set_time, get_time, print_time, &
                             set_time_missing, print_date, set_calendar_type, &
                             operator(>), operator(<), operator(==), &
                             operator(<=), operator(-), operator(+), operator(/=)
use    utilities_mod, only : get_unit, close_file, register_module, &
                             file_exist, error_handler, E_ERR, E_WARN, E_MSG, &
                             initialize_utilities, nmlfileunit, timestamp, &
                             find_namelist_in_file, check_namelist_read, nc_check, &
                             next_file, get_next_filename, find_textfile_dims, &
                             file_to_text, do_nml_file, do_nml_term
use         sort_mod, only : sort

use typeSizes
use netcdf

implicit none

! version controlled file description for error handling, do not edit
character(len=128), parameter :: &
   source   = '$URL$', &
   revision = '$Revision$', &
   revdate  = '$Date$'

!---------------------------------------------------------------------
! An array of these structures will be filled for each obs sequence file.
! After each obs sequence file has been read, the information will be
! written to the appropriate slots in the netCDF file.
!---------------------------------------------------------------------

type station
   integer                  :: obs_type
   type(location_type)      :: location
   type(time_type)          :: first_time
   type(time_type)          :: last_time
   integer                  :: ntimes
   type(time_type), pointer :: times(:)
   integer                  :: analysisT
   integer,  pointer, dimension(  :) :: orgqc
   integer,  pointer, dimension(  :) :: dartqc
   real(r8), pointer, dimension(  :) :: observation       ! (fcstlead)
   real(r8), pointer, dimension(  :) :: obserror          ! (fcstlead)
   real(r8), pointer, dimension(:,:) :: forecast ! (enssize, fcstlead) priors
end type station

! (copy, station, level, ensemble, analysisT, fcstlead)

logical,       allocatable, dimension(:) :: DesiredStations
type(station), allocatable, dimension(:) :: stations
integer :: ensemble_size ! the # of ensemble members in the obs_seq file 
integer :: num_stations  ! This is the current number of unique locations

!---------------------------------------------------------------------
! variables associated with the observation
!---------------------------------------------------------------------

type(obs_sequence_type) :: seq
type(obs_type)          :: obs1, obs2
type(obs_def_type)      :: obs_def
type(location_type)     :: obs_loc

character(len = 129) :: obs_seq_in_file_name
character(len = 129), allocatable, dimension(:) :: obs_seq_filenames

integer :: flavor
integer :: num_copies, num_qc, num_obs, max_num_obs, obs_seq_file_id
integer :: obtype_integer

character(len=metadatalength), allocatable, dimension(:) :: qc_copy_names
real(r8),                      allocatable, dimension(:) :: qc_values
character(len=metadatalength), allocatable, dimension(:) :: obs_copy_names
real(r8),                      allocatable, dimension(:) :: copy_values

character(len=129) :: obs_seq_read_format
logical :: pre_I_format
logical :: last_ob_flag

!-----------------------------------------------------------------------
! Namelist with (some scalar) default values
!-----------------------------------------------------------------------

character(len = 129) :: obs_sequence_name = 'obs_seq.final'
character(len = 129) :: obs_sequence_list = ''
character(len = 129) :: station_template  = 'obsdef_mask.nc'
character(len = 129) :: netcdf_out        = 'obs_seq_fcst.nc'
character(len = 129) :: calendar          = 'gregorian'
character(len = obstypelength) :: obtype_string

logical :: verbose = .false.
logical :: debug   = .false.   ! undocumented ... on purpose

namelist /obs_seq_verify_nml/ obs_sequence_name, obs_sequence_list, &
                             station_template, netcdf_out, &
                             calendar, obtype_string, verbose, debug

!-----------------------------------------------------------------------
! Quantities of interest
!-----------------------------------------------------------------------

integer ::       qc_index   ! copy index of the original qc value
integer ::  dart_qc_index   ! copy index of the DART qc value
integer :: obs_copy_index   ! copy index of the observation
integer :: station_id
integer :: time_index       ! verification time/forecast lead time index
integer :: AnalysisIndex    ! index of the forecast experiment 
character(len=metadatalength), dimension(:), allocatable :: module_obs_copy_names
integer,                       dimension(:), allocatable :: copy_indices
integer,                       dimension(:), allocatable :: forecast_leads
real(digits12),              dimension(:,:), allocatable :: ExperimentTimes
type(time_type),             dimension(:,:), allocatable :: VerifyTimes

!-----------------------------------------------------------------------
! General purpose variables
!-----------------------------------------------------------------------

integer  :: ifile, nread, ngood
integer  :: i, io, ncunit

type(time_type) :: obs_time
type(time_type) :: analysisT ! valid time of analysis at start of forecast 

character(len = 129) :: ncName, string1, string2

! ~# of degrees for 1/2 meter at Earth equator 
! 360 deg-earth/(40000 km-earth * 1000m-km)
real(r8), parameter :: HALF_METER = 180.0_r8 / (40000.0_r8 * 1000.0_r8)

!=======================================================================
! Get the party started
!=======================================================================

call initialize_utilities('obs_seq_verify')
call register_module(source,revision,revdate) 
call static_init_obs_sequence()  ! Initialize the obs sequence module 

call init_obs(obs1, 0, 0)
call init_obs(obs2, 0, 0)
call init_obs_sequence(seq,0,0,0)

!----------------------------------------------------------------------
! Read the namelist
!----------------------------------------------------------------------

call find_namelist_in_file('input.nml', 'obs_seq_verify_nml', ncunit)
read(ncunit, nml = obs_seq_verify_nml, iostat = io)
call check_namelist_read(ncunit, io, 'obs_seq_verify_nml')

! Record the namelist values used for the run ...
if (do_nml_file()) write(nmlfileunit, nml=obs_seq_verify_nml)
if (do_nml_term()) write(    *      , nml=obs_seq_verify_nml)

string1 = adjustl(obtype_string)
string2 = adjustl(netcdf_out)
obtype_string  = trim(string1)
netcdf_out     = trim(string2)
obtype_integer = get_obs_kind_index(obtype_string)

if (obtype_integer < 1) then
   write(string1,*)'obtype_string ',trim(obtype_string),' is unknown. change input.nml'
   call error_handler(E_ERR,'obs_seq_verify',string1,source,revision,revdate)
endif

call set_calendar_type(calendar)

! Check the user input for sanity
if ((obs_sequence_name /= '') .and. (obs_sequence_list /= '')) then
   write(string1,*)'specify "obs_sequence_name" or "obs_sequence_list"'
   write(string2,*)'set other to an empty string ... i.e. ""'
   call error_handler(E_ERR, 'obs_seq_verify', string1, &
                   source, revision, revdate, text2=string2)
endif

! Determine the number of stations from the input netcdf file,
! Initialize the output netcdf file. Must know the ensemble size.

ensemble_size = find_ensemble_size() ! from the obs_seq copy metadata
num_stations  = fill_stations( station_template, ensemble_size )

allocate(copy_indices(ensemble_size), module_obs_copy_names(ensemble_size))

ncName = trim(obtype_string)//'_'//trim(netcdf_out)

ncunit = InitNetCDF(trim(ncName))

!----------------------------------------------------------------------
! Prepare the variables
!----------------------------------------------------------------------

allocate(obs_seq_filenames(1000))
obs_seq_filenames = 'null'

ObsFileLoop : do ifile=1, size(obs_seq_filenames)
!-----------------------------------------------------------------------

  ! Because of the ability to 'cycle' the ObsFileLoop, we need to
  ! destroy and deallocate at the top of the loop.

   call destroy_obs(obs1)
   call destroy_obs(obs2)
   call destroy_obs_sequence(seq)

   if (allocated(   copy_values)) deallocate(   copy_values)
   if (allocated(     qc_values)) deallocate(     qc_values)
   if (allocated( qc_copy_names)) deallocate( qc_copy_names)
   if (allocated(obs_copy_names)) deallocate(obs_copy_names)

   ! Determine the next input filename ... 

   if (obs_sequence_list == '') then
      obs_seq_in_file_name = trim(next_file(obs_sequence_name,ifile))
   else
      obs_seq_in_file_name = trim(get_next_filename(obs_sequence_list,ifile))
      if (obs_seq_in_file_name == '') exit ObsFileLoop
   endif

   if ( file_exist(trim(obs_seq_in_file_name)) ) then
      write(*,*) ! whitespace
      write(string1,*)'opening ', trim(obs_seq_in_file_name)
      call error_handler(E_MSG,'obs_seq_verify:',string1,source,revision,revdate)
   else
      write(string1,*)trim(obs_seq_in_file_name),&
                        ' does not exist. Finishing up.'
      call error_handler(E_MSG,'obs_seq_verify:',string1,source,revision,revdate)
      exit ObsFileLoop
   endif

   ! Read in information about observation sequence so we can allocate
   ! observations. We need info about how many copies, qc values, etc.

   obs_seq_filenames(ifile) = trim(obs_seq_in_file_name)

   ! Determine the analysis time from the file name

   call find_analysis_time(obs_seq_in_file_name, AnalysisIndex, analysisT)

   if (verbose) then
      write(string1,*)'analysis ',AnalysisIndex,' date is '
      call print_date(analysisT,string1)
   endif

   call read_obs_seq_header(obs_seq_in_file_name, &
             num_copies, num_qc, num_obs, max_num_obs, &
             obs_seq_file_id, obs_seq_read_format, pre_I_format, &
             close_the_file = .true.)

   ! Initialize some (individual) observation variables

   call init_obs(obs1, num_copies, num_qc) ! First obs in sequence
   call init_obs(obs2, num_copies, num_qc)

   if ((num_qc <= 0) .or. (num_copies <=0)) then
      write(string1,*)'need at least 1 qc and 1 observation copy'
      call error_handler(E_ERR,'obs_seq_verify',string1,source,revision,revdate)
   endif

   allocate(    copy_values(num_copies),     qc_values(num_qc))
   allocate( obs_copy_names(num_copies), qc_copy_names(num_qc))

   if ( debug ) then
      write(*,*)
      write(*,*)'num_copies          is ',num_copies
      write(*,*)'num_qc              is ',num_qc
      write(*,*)'num_obs             is ',num_obs
      write(*,*)'max_num_obs         is ',max_num_obs
      write(*,*)'obs_seq_read_format is ',trim(obs_seq_read_format)
      write(*,*)
   endif

   if (num_obs <= 1) then
      last_ob_flag = .TRUE.
   else
      last_ob_flag = .FALSE.
   endif

   if ( debug ) write(*,*)'num obs/last_ob_flag are ',num_obs, last_ob_flag

   !--------------------------------------------------------------------
   ! * Read the entire observation sequence - allocates 'seq' internally
   ! * Find the N copy indices we want
   ! * Do something with the QC values?
   !     - I suspect they will all be 'evaluate only'
   !     - what if the prior forward operator fails?
   !--------------------------------------------------------------------

   call read_obs_seq(obs_seq_in_file_name, 0, 0, 0, seq)

   do i=1, num_copies
         string1 = trim(get_copy_meta_data(seq,i))
         obs_copy_names(i) = adjustl(string1)
   enddo
   do i=1, num_qc
         string1 = trim(get_qc_meta_data(seq,i))
         qc_copy_names(i) = adjustl(string1)
   enddo

   ! Determine which qc copy is original qc and dart qc 

   call find_our_copies(seq, obs_copy_index, copy_indices, qc_index, dart_qc_index )

   ! The first trip through sets the module_obs_copy_names so we can
   ! be sure we are stuffing compatible objects into the same slots 
   if ( ifile == 1 ) then
      module_obs_copy_names = obs_copy_names
   else
      ! FIXME more robust checks
      do i = 1,ensemble_size
         if ( obs_copy_names(copy_indices(i)) /= module_obs_copy_names(i) ) then
         
            write(string1,*)'mismatch in observation copies ',&
                             trim(obs_copy_names(copy_indices(i)))
            string2 = trim(module_obs_copy_names(i))
            call error_handler(E_ERR,'obs_seq_verify', &
               string1,source,revision,revdate,text2=string2)
         endif
      enddo
   endif

   !--------------------------------------------------------------------
   ! Read the first observation in the sequence
   !--------------------------------------------------------------------

   ngood = 0

   if ( .not. get_first_obs(seq, obs1) )              &
      call error_handler(E_ERR,'obs_seq_verify', &
              'No first observation in sequence.',    &
              source,revision,revdate)

   !--------------------------------------------------------------------
   ObservationLoop : do nread = 1,num_obs
   !--------------------------------------------------------------------

      if ( verbose .and. (mod(nread,10000) == 0) ) &
         write(*,*)'Processing obs ',nread,' of ',num_obs

      call get_obs_def(obs1,          obs_def)
      flavor   = get_obs_kind(        obs_def)
      obs_time = get_obs_def_time(    obs_def)
      obs_loc  = get_obs_def_location(obs_def)

      if (debug .and. (nread == 1)) then
         call print_time(obs_time,'First observation time')
         call print_date(obs_time,'First observation date')
         call get_qc(        obs1,   qc_values)
         call get_obs_values(obs1, copy_values)
         write(*,*)'looking for observation kinds of ',obtype_integer, &
                                     get_obs_kind_name(obtype_integer)

         write(*,*)'first observation "kind" is ',flavor, &
                                get_obs_kind_name(flavor)
         write(*,*)'observation values are:'
         do i = 1,size(copy_values)
            write(*,*)i,copy_values(i)
         enddo
         write(*,*)'observation QC values are:'
         do i = 1,size(qc_values)
            write(*,*)i,qc_values(i)
         enddo

         call write_location(6,obs_loc,fform='ascii')
         call write_location(6,obs_loc,fform='ascii',charstring=string1)
         write(*,*)trim(string1)

         if (qc_index > 0) &
         write(*,*)'     qc_index and value is ',     qc_index, qc_values(     qc_index)
         if (dart_qc_index > 0) &
         write(*,*)'dart_qc_index and value is ',dart_qc_index, qc_values(dart_qc_index)

      endif

      ! 0) is it the right type of observation
      ! 1) what/if any station does it belong to?
      ! 2) Is it one of the times of interest?
      ! 3) stuff it into the appropriate station structure

      if ( obtype_integer /= flavor ) goto 100

      station_id = find_station_location(flavor, obs_loc, stations)
      if ( station_id < 1 ) goto 100

      if ( is_time_unwanted( obs_time, station_id, stations, time_index) ) &
             goto 100

      if (debug) write(*,*)'obs ',nread,' is station ',station_id,' at time ',time_index

      call get_qc(        obs1,   qc_values)
      call get_obs_values(obs1, copy_values)

      if (     qc_index > 0) stations(station_id)%orgqc( time_index) = qc_values(     qc_index)
      if (dart_qc_index > 0) stations(station_id)%dartqc(time_index) = qc_values(dart_qc_index)

      stations(station_id)%obserror(   time_index) = get_obs_def_error_variance(obs_def)
      stations(station_id)%observation(time_index) = copy_values(obs_copy_index)
      stations(station_id)%forecast(:, time_index) = copy_values(copy_indices)

      if (debug) then   ! GIGANTIC OUTPUT - BE CAREFUL
         do i = 1,size(copy_indices)
            write(*,'(''ob'',3(1x,i6),4(1x,e15.5))') &
            nread, time_index, copy_indices(i), &
            stations(station_id)%obserror(   time_index),&
            stations(station_id)%observation(time_index),&
            stations(station_id)%forecast(i, time_index),&
            copy_values(copy_indices(i))
         enddo
         write(*,*) ! a little whitespace
      endif

 100  continue

      call get_next_obs(seq, obs1, obs2, last_ob_flag)
      if (.not. last_ob_flag) obs1 = obs2

   !--------------------------------------------------------------------
   enddo ObservationLoop
   !--------------------------------------------------------------------
   ! So by now, all the observations have been stuffed into the 'station'
   ! array - not that we know if the station array is complete ... 
   ! Take what we have (i.e. for one analysis time) and push it into
   ! the output forecast netcdf file.

   call WriteNetCDF(ncunit, trim(ncName), stations)

enddo ObsFileLoop

call CloseNetCDF(ncunit, trim(ncName))

!-----------------------------------------------------------------------
! Really, really, done.
!-----------------------------------------------------------------------

call destroy_obs(obs1)
call destroy_obs(obs2)
call destroy_obs_sequence(seq)
call destroy_stations(stations)

if (allocated(qc_values))             deallocate(qc_values)
if (allocated(qc_copy_names))         deallocate(qc_copy_names)
if (allocated(copy_values))           deallocate(copy_values)
if (allocated(obs_copy_names))        deallocate(obs_copy_names)
if (allocated(obs_seq_filenames))     deallocate(obs_seq_filenames)
if (allocated(DesiredStations))       deallocate(DesiredStations)

call timestamp(source,revision,revdate,'end') ! That closes the log file, too.

!======================================================================
CONTAINS
!======================================================================


function find_station_location(ObsType, ObsLocation, stationlist) result(station_id)
!----------------------------------------------------------------------------
! Simply try to find a matching lat/lon for an observation type
! The lons/lats get yanked around "a lot" - being converted from ASCII radians
! to r8 degrees to r8 radians to r8 degrees and then checked for "equality".
! So - we're actually just checking to see if the lat/lon is within something
! like 500 cm either direction. Seems like a reasonable definition of 'match'.
!
integer,                     intent(in) :: ObsType
type(location_type),         intent(in) :: ObsLocation
type(station), dimension(:), intent(in) :: stationlist
integer                                 :: station_id

integer :: i
real(r8), dimension(3) :: obslocarray, stnlocarray
real(r8) :: londiff, latdiff

station_id = 0;

FindLoop : do i = 1,num_stations

   obslocarray = get_location(ObsLocation)  ! returns degrees
   stnlocarray = get_location(stationlist(i)%location)

   londiff = abs(obslocarray(1) - stnlocarray(1)) 
   latdiff = abs(obslocarray(2) - stnlocarray(2)) 

   if ( (londiff <= HALF_METER) .and. &
        (latdiff <= HALF_METER) .and. &
        (ObsType == stationlist(i)%obs_type) ) then

      station_id = i
      exit FindLoop
   endif

enddo FindLoop

end function find_station_location


!============================================================================


function fill_stations( filename, nensemble ) result(nstations)
!----------------------------------------------------------------------------
! Check to make sure the file exists.
! Read the 'station' variable (binary value for desired or not).
! count the number of non-zero entries ...
!-----------------------------------------------------------------------

character(len=*), intent(in) :: filename
integer,          intent(in) :: nensemble
integer                      :: nstations

integer :: DimID, VarID
integer ::  nmax, ntimes, nforecasts, nverify, locNdim, strlen
integer :: ncid, i, j, istation, ndims, mylen
integer :: my_type
integer :: seconds, days
type(time_type)     :: mytime
type(location_type) :: myloc

integer, dimension(nf90_max_var_dims) :: dimIDs

! These temporarily mimic the contents of the netCDF file
character(len=obstypelength), dimension(:), allocatable :: obs_types
integer,                      dimension(:), allocatable :: allstations
real(r8),                   dimension(:,:), allocatable :: locations
integer,                      dimension(:), allocatable :: which_vert
integer,                      dimension(:), allocatable :: ntimearray
real(digits12),               dimension(:), allocatable :: first_time
real(digits12),               dimension(:), allocatable :: last_time
real(digits12),             dimension(:,:), allocatable :: ReportTimes

!-----------------------------------------------------------------------
! Open the netCDF file

call nc_check(nf90_open(trim(filename), NF90_NOWRITE, &
         ncid = ncid), 'fill_stations', 'open '//trim(filename))

!-----------------------------------------------------------------------
! Determine the maximum number of stations from the station dimension,
! the number of times, and the number of dimensions in the location 

! used for integer array indicating yes/no - do we want the station
call nc_check(nf90_inq_dimid(ncid, 'station', dimid=DimID), &
           'fill_stations', 'inquire station dimid '//trim(filename))
call nc_check(nf90_inquire_dimension(ncid, DimID, len=nmax), &
           'fill_stations', 'inquire station nmax '//trim(filename))

! used for array of all possible verification times needed
call nc_check(nf90_inq_dimid(ncid, 'time', dimid=DimID), &
           'fill_stations', 'inquire time dimid '//trim(filename))
call nc_check(nf90_inquire_dimension(ncid, DimID, len=ntimes), &
           'fill_stations', 'inquire time len '//trim(filename))

! used for array of the analysis times / 'zero-length' forecast times  
call nc_check(nf90_inq_dimid(ncid, 'analysisT', dimid=DimID), &
           'fill_stations', 'inquire analysisT dimid '//trim(filename))
call nc_check(nf90_inquire_dimension(ncid, DimID, len=nforecasts), &
           'fill_stations', 'inquire analysisT len '//trim(filename))

! used for array of the verification times 
call nc_check(nf90_inq_dimid(ncid, 'forecast_lead', dimid=DimID), &
           'fill_stations', 'inquire forecast_lead dimid '//trim(filename))
call nc_check(nf90_inquire_dimension(ncid, DimID, len=nverify), &
           'fill_stations', 'inquire forecast_lead len '//trim(filename))

call nc_check(nf90_inq_dimid(ncid, 'stringlength', dimid=DimID), &
           'fill_stations', 'inquire stringlength dimid '//trim(filename))
call nc_check(nf90_inquire_dimension(ncid, DimID, len=strlen), &
           'fill_stations', 'inquire stringlength len '//trim(filename))

call nc_check(nf90_inq_dimid(ncid, 'location', dimid=DimID), &
           'fill_stations', 'inquire location dimid '//trim(filename))
call nc_check(nf90_inquire_dimension(ncid, DimID, len=locNdim), &
           'fill_stations', 'inquire location len '//trim(filename))

!-----------------------------------------------------------------------
! Read the "stations" 1D variable. Each non-zero entry indicates a
! 'station' we wanted to keep when we ran obs_seq_coverage.f90.

call nc_check(nf90_inq_varid(ncid, 'station', varid=VarID), &
          'fill_stations', 'inq_varid:station '//trim(filename))
call nc_check(nf90_inquire_variable(ncid, VarID, ndims=ndims, dimids=dimIDs), &
           'fill_stations', 'inquire station ndims '//trim(filename))

if (ndims /= 1) then
   write(string1,*)'station is expected to be a 1D array, it is ',ndims
   call error_handler(E_MSG,'fill_stations:',string1,source,revision,revdate)
endif

call nc_check(nf90_inquire_dimension(ncid, dimIDs(1), len=mylen), &
           'fill_stations', 'station inquire dimid(1) '//trim(filename))

if (mylen /= nmax) then
   write(string1,*)'station length expected to be ',nmax,' it is ',mylen
   call error_handler(E_MSG,'fill_stations:',string1,source,revision,revdate)
endif

allocate(allstations(nmax),  obs_types(nmax), which_vert(nmax), &
          ntimearray(nmax), first_time(nmax),  last_time(nmax)) 
allocate(locations(locNdim,nmax))
allocate(forecast_leads(nverify))       ! how far into the forecast (seconds)
allocate(ReportTimes(nverify,nmax))    ! observation times that are closest
allocate(ExperimentTimes(nverify,nforecasts)) ! verification times of interest 
allocate(VerifyTimes(nforecasts,nverify))     ! verification times of interest 

!-----------------------------------------------------------------------
! Read all the information from the station list netCDF file.

call nc_check(nf90_get_var(ncid, VarID, allstations), &
                   'fill_stations', 'get_var:allstations '//trim(filename))

call nc_check(nf90_inq_varid(ncid, 'obs_type', varid=VarID), &
        'fill_stations', 'inq_varid:obs_type '//trim(filename))
call nc_check(nf90_get_var(ncid, VarID, obs_types), &
        'fill_stations',   'get_var:obs_type '//trim(filename))

call nc_check(nf90_inq_varid(ncid, 'location', varid=VarID), &
        'fill_stations', 'inq_varid:location '//trim(filename))
call nc_check(nf90_get_var(ncid, VarID, locations), &
        'fill_stations',   'get_var:location '//trim(filename))

call nc_check(nf90_inq_varid(ncid, 'which_vert', varid=VarID), &
        'fill_stations', 'inq_varid:which_vert '//trim(filename))
call nc_check(nf90_get_var(ncid, VarID, which_vert), &
        'fill_stations',   'get_var:which_vert '//trim(filename))

call nc_check(nf90_inq_varid(ncid, 'ntimes', varid=VarID), &
        'fill_stations', 'inq_varid:ntimes '//trim(filename))
call nc_check(nf90_get_var(ncid, VarID, ntimearray), &
        'fill_stations',   'get_var:ntimes '//trim(filename)) ! unused

call nc_check(nf90_inq_varid(ncid, 'forecast_lead', varid=VarID), &
        'fill_stations', 'inq_varid:forecast_lead '//trim(filename))
call nc_check(nf90_get_var(ncid, VarID, forecast_leads), &
        'fill_stations',   'get_var:forecast_lead '//trim(filename))

call nc_check(nf90_inq_varid(ncid, 'verification_times', varid=VarID), &
        'fill_stations', 'inq_varid:verification_times '//trim(filename))
call nc_check(nf90_get_var(ncid, VarID, ExperimentTimes), &
        'fill_stations',   'get_var:verification_times '//trim(filename))

call nc_check(nf90_inq_varid(ncid, 'first_time', varid=VarID), &
        'fill_stations', 'inq_varid:first_time '//trim(filename))
call nc_check(nf90_get_var(ncid, VarID, first_time), &
        'fill_stations',   'get_var:first_time '//trim(filename))

call nc_check(nf90_inq_varid(ncid, 'last_time', varid=VarID), &
        'fill_stations', 'inq_varid:last_time '//trim(filename))
call nc_check(nf90_get_var(ncid, VarID, last_time), &
        'fill_stations',   'get_var:last_time '//trim(filename))

call nc_check(nf90_inq_varid(ncid, 'ReportTime', varid=VarID), &
        'fill_stations', 'inq_varid:ReportTime '//trim(filename))
call nc_check(nf90_get_var(ncid, VarID, ReportTimes), &
        'fill_stations',   'get_var:ReportTime '//trim(filename))

call nc_check(nf90_close(ncid), &
        'fill_stations', 'close '//trim(filename))

!-----------------------------------------------------------------------
! Convert all the Experiment Times to DART time types.
! While we're at it, let's reshape them to our liking:

if (debug) write(*,*)'size(VerifyTimes,1) is ',size(VerifyTimes,1),'nforecasts is ',nforecasts
if (debug) write(*,*)'size(VerifyTimes,2) is ',size(VerifyTimes,2),'nverify is',nverify

do j = 1,nforecasts
do i = 1,nverify

   days    = floor(ExperimentTimes(i,j))
   seconds = nint((ExperimentTimes(i,j) - days) * 86400.0_digits12)
   VerifyTimes(j,i) = set_time(seconds, days)
   
enddo
enddo

!-----------------------------------------------------------------------
! Allocate the desired number of stations. 
! The 'allstations' array contains either 0 (unwanted) or 1 (desired)
! from the temporal coverage standpoint. Combine that with namelist input
! to generate subset of stations we want.

GoodLoop : do istation = 1,nmax
   if (allstations(istation) < 1) cycle GoodLoop

   ! flag the station as uninteresting if wrong type
   my_type = get_obs_kind_index(obs_types(istation))
   if (my_type /= obtype_integer) allstations(istation) = 0

enddo GoodLoop

nstations = sum(allstations)

! check to make sure nstations is not zero

if (nstations < 1) then
   write(string1,*)'no valid stations of ',trim(obtype_string)
   call error_handler(E_ERR,'fill_stations',string1,source,revision,revdate)
endif

allocate(stations(nstations))   ! global scope

i = 0
StationLoop : do istation = 1,nmax

   if (allstations(istation) < 1) cycle StationLoop

   i = i + 1

   if ( i > nstations ) then
      write(string1,*)'count of "good" stations wrong - got',nstations
      call error_handler(E_ERR,'fill_stations',string1,source,revision,revdate)
   endif

   stations(i)%obs_type = get_obs_kind_index(obs_types(istation))

   myloc = set_location(locations(1,istation),  locations(2,istation), &
                        locations(3,istation), which_vert(istation) ) 
   stations(i)%location   = myloc

   days    = floor(first_time(istation))
   seconds = nint((first_time(istation) - days) * 86400.0_digits12)
   mytime  = set_time(seconds, days)
   stations(i)%first_time = mytime

   days    = floor(last_time(istation))
   seconds = nint((last_time(istation) - days) * 86400.0_digits12)
   mytime  = set_time(seconds, days)
   stations(i)%last_time  = mytime

   stations(i)%ntimes     = nverify

   allocate(stations(i)%times(      nverify))
   allocate(stations(i)%orgqc(      nverify))
   allocate(stations(i)%dartqc(     nverify))
   allocate(stations(i)%obserror(   nverify))
   allocate(stations(i)%observation(nverify))
   allocate(stations(i)%forecast(nensemble,nverify))

   do j = 1,nverify
      days    = floor(ReportTimes(j,istation))
      seconds = nint((ReportTimes(j,istation) - days) * 86400.0_digits12)
      mytime  = set_time(seconds, days)
      stations(i)%times(j) = mytime
   enddo
   
   stations(i)%orgqc       = MISSING_I
   stations(i)%dartqc      = MISSING_I
   stations(i)%obserror    = MISSING_R8
   stations(i)%observation = MISSING_R8
   stations(i)%forecast    = MISSING_R8
   ! The only thing left to fill is the analaysisT ... which
   ! comes from the obs_seq.fcst directory name or something. 
   ! It surely does not come from the station list netcdf file.

enddo StationLoop

!-----------------------------------------------------------------------
! Print a summary of all the station information if so desired.
! FIXME can think of lots of other things to check ... 

if (verbose) then
   write(*,*) ! whitespace
   write(string1,*)'There are ',nstations,' stations of interest,'
   write(string2,*)'and   ',nverify,  ' times    of interest.'
   call error_handler(E_MSG,'fill_stations:',string1,text2=string2)
endif

! housekeeping

deallocate(obs_types, allstations, locations, which_vert, &
           ntimearray, first_time,  last_time, ReportTimes )

end function fill_stations


!============================================================================


function is_time_unwanted(ObsTime, stationid, stationlist, timeindex)
!----------------------------------------------------------------------------
! Determine if the observation time is one
! of the times for the particular station.
! Due to roundoff, we will tolerate up to a 1 second difference.

type(time_type),             intent(in)  :: ObsTime
integer,                     intent(in)  :: stationid
type(station), dimension(:), intent(in)  :: stationlist
integer,                     intent(out) :: timeindex
logical                                  :: is_time_unwanted

integer :: i
type(time_type) :: one_second, dt


one_second = set_time(1,0)

timeindex        = 0
is_time_unwanted = .TRUE.

if ( stationlist(stationid)%ntimes == 0 ) return

TimeLoop : do i = 1,stationlist(stationid)%ntimes

   dt = stationlist(stationid)%times(i) - ObsTime

   if (dt <= one_second) then
      if (debug) write(string1,*)'desired ob at station ',stationid,' time ',i
      if (debug) call print_time(ObsTime,trim(string1))

      timeindex        = i
      is_time_unwanted = .FALSE.
      exit TimeLoop
   endif

enddo TimeLoop

end function is_time_unwanted


!============================================================================


function InitNetCDF(fname)
!----------------------------------------------------------------------------
character(len=*), intent(in) :: fname
integer                      :: InitNetCDF

integer :: ncid, i, nlines, linelen
integer :: AnalysisDimID,     CopyDimID, StationsDimID
integer ::   LevelsDimID, NmembersDimID, ForecastDimID
integer ::  LineLenDimID,   nlinesDimID,   stringDimID
integer :: VarID, LocationVarID, WhichVertVarID

character(len=8)      :: crdate      ! needed by F90 DATE_AND_TIME intrinsic
character(len=10)     :: crtime      ! needed by F90 DATE_AND_TIME intrinsic
character(len=5)      :: crzone      ! needed by F90 DATE_AND_TIME intrinsic
integer, dimension(8) :: values      ! needed by F90 DATE_AND_TIME intrinsic

character(len=129), allocatable, dimension(:) :: textblock

integer, parameter :: ncopies = 3 ! obs value, prior, obs error variance, orgqc, dartqc
integer, parameter :: nlevels = 1 ! hardwired for now - do no level subsetting anyway

character(len=metadatalength), dimension(ncopies) :: copy_metadata = &
     (/ 'observation               ', &
        'forecast                  ', &
        'observation error variance' /)
 
if(.not. byteSizesOK()) then
    call error_handler(E_ERR,'InitNetCDF', &
   'Compiler does not support required kinds of variables.',source,revision,revdate)
endif

InitNetCDF = 0

call nc_check(nf90_create(path = trim(fname), cmode = nf90_clobber, &
         ncid = ncid), 'obs_seq_verify:InitNetCDF', 'create '//trim(fname))

if (verbose) then
   write(string1,*)trim(ncName), ' is fortran unit ',ncid
   call error_handler(E_MSG,'InitNetCDF:',string1,source,revision,revdate)
endif

!----------------------------------------------------------------------------
! Write Global Attributes (mostly namelist input, that sort of thing)
!----------------------------------------------------------------------------

call DATE_AND_TIME(crdate,crtime,crzone,values)
write(string1,'(''YYYY MM DD HH MM SS = '',i4,5(1x,i2.2))') &
               values(1), values(2), values(3), values(5), values(6), values(7)
call nc_check(nf90_put_att(ncid, NF90_GLOBAL, 'creation_date', trim(string1) ), &
           'InitNetCDF', 'put_att creation_date '//trim(fname))

call nc_check(nf90_put_att(ncid, NF90_GLOBAL, 'source',   source   ), &
           'InitNetCDF', 'put_att   source '//trim(fname))
call nc_check(nf90_put_att(ncid, NF90_GLOBAL, 'revision', revision ), &
           'InitNetCDF', 'put_att revision '//trim(fname))
call nc_check(nf90_put_att(ncid, NF90_GLOBAL, 'revdate',  revdate  ), &
           'InitNetCDF', 'put_att  revdate '//trim(fname))

!----------------------------------------------------------------------------
! Define the dimensions
! Set nofill mode - supposed to be performance gain
!----------------------------------------------------------------------------

call nc_check(nf90_set_fill(ncid, NF90_NOFILL, i),  &
             'obs_seq_verify:InitNetCDF', 'set_nofill '//trim(fname))

call nc_check(nf90_def_dim(ncid=ncid, name='analysisT', len=NF90_UNLIMITED, &
               dimid = AnalysisDimID), 'InitNetCDF', 'def_dim:analysisT '//trim(fname))

call nc_check(nf90_def_dim(ncid=ncid, name='copy',      len=ncopies, &
               dimid =     CopyDimID), 'InitNetCDF', 'def_dim:copy '//trim(fname))

call nc_check(nf90_def_dim(ncid=ncid, name='station',  len=num_stations, &
               dimid = StationsDimID), 'InitNetCDF', 'def_dim:station '//trim(fname))

call nc_check(nf90_def_dim(ncid=ncid, name='level',    len=nlevels, &
               dimid =   LevelsDimID), 'InitNetCDF', 'def_dim:level '//trim(fname))

call nc_check(nf90_def_dim(ncid=ncid, name='ensemble',  len=ensemble_size, &
               dimid = NmembersDimID), 'InitNetCDF', 'def_dim:ensemble '//trim(fname))

call nc_check(nf90_def_dim(ncid=ncid, name='forecast_lead',  len=stations(1)%ntimes, &
               dimid = ForecastDimID), 'InitNetCDF', 'def_dim:forecast_lead '//trim(fname))

! namelist quantities

call find_textfile_dims('input.nml', nlines, linelen)
allocate(textblock(nlines))
textblock = ''

call nc_check(nf90_def_dim(ncid=ncid, name='linelen', len = len(textblock(1)), &
               dimid =  linelenDimID), 'InitNetCDF', 'def_dim:linelen '//'input.nml')

call nc_check(nf90_def_dim(ncid=ncid, name='nlines', len = nlines, &
               dimid =   nlinesDimID), 'InitNetCDF', 'def_dim:nlines '//'input.nml')

call nc_check(nf90_def_dim(ncid=ncid, name='stringlength', len = metadatalength, &
               dimid =  StringDimID),  'InitNetCDF', 'def_dim:stringlength '//trim(fname))

!----------------------------------------------------------------------------
! Define the coordinate variables
!----------------------------------------------------------------------------

call nc_check(nf90_def_var(ncid=ncid, name='analysisT', xtype=nf90_double, &
          dimids=(/ AnalysisDimID /), varid=VarID), &
          'InitNetCDF', 'analysisT:def_var')
call nc_check(nf90_put_att(ncid, VarID, 'long_name', 'time of analysis'), &
          'InitNetCDF', 'analysisT:put_att long_name')
call nc_check(nf90_put_att(ncid, VarID, 'units',     'days since 1601-1-1'), &
          'InitNetCDF', 'analysisT:put_att units')
call nc_check(nf90_put_att(ncid, VarID, 'calendar',  'Gregorian'), &
          'InitNetCDF', 'analysisT:put_att calendar')
call nc_check(nf90_put_att(ncid, VarID, 'missing_value', 0.0_digits12), &
          'InitNetCDF', 'analysisT:put_att missing')
call nc_check(nf90_put_att(ncid, VarID, '_FillValue',    0.0_digits12), &
          'InitNetCDF', 'analysisT:put_att fill_value')

call nc_check(nf90_def_var(ncid=ncid, name='copy', xtype=nf90_int, &
          dimids = (/ CopyDimID /), varid=VarID), &
          'InitNetCDF', 'copy:def_var')
call nc_check(nf90_put_att(ncid, VarID, 'long_name', 'observation copy'), &
          'InitNetCDF', 'copy:long_name')

do i = 1,ncopies
   write(string1,'(''note'',i1)')i
   write(string2,'(i1,'' == '',a)')i,trim(copy_metadata(i))
   call nc_check(nf90_put_att(ncid, VarID, trim(string1),trim(string2)), &
          'InitNetCDF', 'copy'//trim(string1))
enddo

call nc_check(nf90_put_att(ncid, VarID, 'explanation', 'see CopyMetaData variable'), &
          'InitNetCDF', 'copy:explanation')

call nc_check(nf90_def_var(ncid=ncid, name='station', xtype=nf90_int, &
          dimids = (/ StationsDimID /), varid=VarID), &
          'InitNetCDF', 'station:def_var')
call nc_check(nf90_put_att(ncid, VarID, 'long_name', 'station index'), &
          'InitNetCDF', 'station:long_name')

call nc_check(nf90_def_var(ncid=ncid, name='original_qc', xtype=nf90_int, &
          dimids = (/ ForecastDimID, StationsDimID, AnalysisDimID /), varid=VarID), &
          'InitNetCDF', 'original_qc:def_var')
call nc_check(nf90_put_att(ncid, VarID, 'long_name', 'original QC value'), &
          'InitNetCDF', 'original_qc:long_name')
call nc_check(nf90_put_att(ncid, VarID, 'missing_value', MISSING_I), &
          'InitNetCDF', 'original_qc:missing')
call nc_check(nf90_put_att(ncid, VarID, '_FillValue',    MISSING_I), &
          'InitNetCDF', 'original_qc:fill_value')

call nc_check(nf90_def_var(ncid=ncid, name='dart_qc', xtype=nf90_int, &
          dimids = (/ ForecastDimID, StationsDimID, AnalysisDimID /), varid=VarID), &
          'InitNetCDF', 'dart_qc:def_var')
call nc_check(nf90_put_att(ncid, VarID, 'long_name', 'DART QC value'), &
          'InitNetCDF', 'dart_qc:long_name')
call nc_check(nf90_put_att(ncid, VarID, 'explanation1', '1 == prior evaluated only'), &
          'InitNetCDF', 'dart_qc:explanation1')
call nc_check(nf90_put_att(ncid, VarID, 'explanation2', '4 == forward operator failed'), &
          'InitNetCDF', 'dart_qc:explanation1')
call nc_check(nf90_put_att(ncid, VarID, 'missing_value', MISSING_I), &
          'InitNetCDF', 'dart_qc:missing')
call nc_check(nf90_put_att(ncid, VarID, '_FillValue',    MISSING_I), &
          'InitNetCDF', 'dart_qc:fill_value')

call nc_check(nf90_def_var(ncid=ncid, name='level', xtype=nf90_double, &
          dimids = (/ LevelsDimID /), varid=VarID), &
          'InitNetCDF', 'level:def_var')
call nc_check(nf90_put_att(ncid, VarID, 'long_name', 'vertical level of observation'), &
          'InitNetCDF', 'level:long_name')

call nc_check(nf90_def_var(ncid=ncid, name='ensemble', xtype=nf90_int, &
          dimids = (/ NmembersDimID /), varid=VarID), &
          'InitNetCDF', 'ensemble:def_var')
call nc_check(nf90_put_att(ncid, VarID, 'long_name', 'ensemble member'), &
          'InitNetCDF', 'ensemble:long_name')

call nc_check(nf90_def_var(ncid=ncid, name='forecast_lead', xtype=nf90_int, &
          dimids = (/ ForecastDimID /), varid=VarID), &
          'InitNetCDF', 'forecast_lead:def_var')
call nc_check(nf90_put_att(ncid, VarID, 'long_name', 'forecast lead time'), &
          'InitNetCDF', 'forecast_lead:long_name')
call nc_check(nf90_put_att(ncid, VarID, 'units',     'seconds'), &
          'InitNetCDF', 'forecast_lead:units')

!----------------------------------------------------------------------------
! Define the RECORD variables
!----------------------------------------------------------------------------

! let the location module write what it needs to ...

if ( nc_write_location_atts( ncid, fname, StationsDimID ) /= 0 ) then
   write(string1,*)'problem initializing netCDF location attributes'
   call error_handler(E_ERR,'InitNetCDF',string1,source,revision,revdate)
endif

! Define the variable to record the input parameters ... the namelist

call nc_check(nf90_def_var(ncid=ncid, name='namelist', xtype=nf90_char, &
          dimids = (/ linelenDimID, nlinesDimID /), varid=VarID), &
          'InitNetCDF', 'namelist:def_var')
call nc_check(nf90_put_att(ncid, VarID, 'long_name', 'input.nml contents'), &
          'InitNetCDF', 'namelist:long_name')

! Define the variable for interpreting the 'copy' dimension

call nc_check(nf90_def_var(ncid=ncid, name='CopyMetaData', xtype=nf90_char, &
          dimids = (/ StringDimID, CopyDimID /), varid=VarID), &
          'InitNetCDF', 'copymeta:def_var')
call nc_check(nf90_put_att(ncid, VarID, 'long_name', 'copy quantity names'), &
          'InitNetCDF', 'copymeta:long_name')

!----------------------------------------------------------------------------
! And now for the main event. For each observation type, we need:
!
! (analysisT, station, level, copy, ensemble, forecast_lead)
!      |         |       |      |      |         |
!      |         |       |      |      |         +- forecast length : 0,3,6,9,12,...
!      |         |       |      |      |
!      |         |       |      |      +----------- ensemble member index
!      |         |       |      |
!      |         |       |      +------------------ obs value, prior, obs_err
!      |         |       |
!      |         |       +------------------------- vertical level index
!      |         |
!      |         +--------------------------------- (horizontal) station index
!      |
!      +------------------------------------------- analysis time/date
!
!----------------------------------------------------------------------------

call nc_check(nf90_def_var(ncid=ncid, name=trim(obtype_string), xtype=nf90_double, &
           dimids = (/ ForecastDimID, NmembersDimID,     CopyDimID, &
                         LevelsDimID, StationsDimID, AnalysisDimID /), &
           varid = VarID), 'InitNetCDF', 'obtype_string:def_var')

call nc_check(nf90_put_att(ncid, VarID, 'long_name', 'forecast variable quantities'), &
          'InitNetCDF', 'obtype_string:long_name')
call nc_check(nf90_put_att(ncid, VarID, 'missing_value', MISSING_R8), &
          'InitNetCDF', 'obtype_string:missing')
call nc_check(nf90_put_att(ncid, VarID, '_FillValue',    MISSING_R8), &
          'InitNetCDF', 'obtype_string:fill_value')

!----------------------------------------------------------------------------
! Leave define mode so we can fill
! The analysisT variable is filled as new obs_sequence files are read.
!----------------------------------------------------------------------------

call nc_check(nf90_enddef(ncid), 'InitNetCDF', 'enddef '//trim(fname))

call file_to_text('input.nml', textblock)

call nc_check(nf90_inq_varid(ncid, 'namelist', varid=VarID), &
          'InitNetCDF', 'inq_varid:namelist '//trim(fname))
call nc_check(nf90_put_var(ncid, VarID, textblock ), &
          'InitNetCDF', 'put_var:namelist')

deallocate(textblock)

call nc_check(nf90_inq_varid(ncid, 'CopyMetaData', varid=VarID), &
          'InitNetCDF', 'inq_varid:CopyMetaData '//trim(fname))
call nc_check(nf90_put_var(ncid, VarID, copy_metadata), &
          'InitNetCDF', 'put_var:CopyMetaData')

call nc_check(nf90_inq_varid(ncid, 'copy', varid=VarID), &
          'InitNetCDF', 'inq_varid:copy '//trim(fname))
call nc_check(nf90_put_var(ncid, VarID, (/ (i,i=1,ncopies) /) ), &
          'InitNetCDF', 'put_var:copy')

call nc_check(nf90_inq_varid(ncid, 'ensemble', varid=VarID), &
          'InitNetCDF', 'inq_varid:ensemble '//trim(fname))
call nc_check(nf90_put_var(ncid, VarID, (/ (i,i=1,ensemble_size) /) ), &
          'InitNetCDF', 'put_var:ensemble')

call nc_check(nf90_inq_varid(ncid, 'station', varid=VarID), &
          'InitNetCDF', 'inq_varid:station '//trim(fname))
call nc_check(nf90_put_var(ncid, VarID, (/ (i,i=1,num_stations) /) ), &
          'InitNetCDF', 'put_var:station')

call nc_check(nf90_inq_varid(ncid, 'level', varid=VarID), &
          'InitNetCDF', 'inq_varid:level '//trim(fname))
call nc_check(nf90_put_var(ncid, VarID, (/ -9999 /) ), &
          'InitNetCDF', 'put_var:level')

call nc_check(nf90_inq_varid(ncid, 'forecast_lead', varid=VarID), &
          'InitNetCDF', 'inq_varid:forecast_lead '//trim(fname))
call nc_check(nf90_put_var(ncid, VarID, forecast_leads ), &
          'InitNetCDF', 'put_var:forecast_lead')

! Record station locations using location_mod:nc_write_location()

call nc_get_location_varids(ncid, fname, LocationVarID, WhichVertVarID)

do i = 1,num_stations
   call nc_write_location(ncid, LocationVarId, stations(i)%location, &
             i, WhichVertVarId)
enddo

!----------------------------------------------------------------------------
! Finish up ...
!----------------------------------------------------------------------------

call nc_check(nf90_sync( ncid), 'InitNetCDF', 'sync '//trim(fname))  

InitNetCDF = ncid

end function InitNetCDF


!============================================================================


subroutine WriteNetCDF(ncid, fname, stations)
!----------------------------------------------------------------------------
integer,                     intent(in) :: ncid
character(len=*),            intent(in) :: fname
type(station), dimension(:), intent(in) :: stations

integer, dimension(1) :: istart, icount

integer :: ilev, obscopyindex, priorcopyindex, errorcopyindex, ensmem
integer :: stationindex, i, seconds, days

integer, dimension(nf90_max_var_dims) :: dimIDs, dimLengths
integer ::  TimeVarID, VarID, ndims
integer ::  QCVarID, DARTQCVarID

integer :: AnalysisDimID  ! dimension ID for the number of analysis times
integer ::     CopyDimID  ! dimension ID for the number of 'copies'
integer :: StationsDimID  ! dimension ID for each station/location
integer ::   LevelsDimID  ! dimension ID for each level
integer :: NmembersDimID  ! dimension ID for the number of ensemble members
integer :: ForecastDimID  ! dimension ID for each forecast step

integer :: AnalysisDimlen
integer ::     CopyDimlen
integer :: StationsDimlen
integer ::   LevelsDimlen
integer :: NmembersDimlen
integer :: ForecastDimlen

real(digits12) :: fdays
real(digits12), allocatable, dimension(:) :: mytimes
real(r8), allocatable, dimension(:,:,:,:) :: datmat

if (debug) write(*,*)'DEBUG --- entering WriteNetCDF'

!----------------------------------------------------------------------------
! Find the dimension ID's and lengths
!----------------------------------------------------------------------------

call nc_check(nf90_inq_dimid(ncid, 'copy', dimid=CopyDimID), &
           'WriteNetCDF', 'inquire copy dimid '//trim(fname))
call nc_check(nf90_inquire_dimension(ncid, CopyDimID, len=CopyDimlen), &
           'WriteNetCDF', 'inquire copy dimlen '//trim(fname))

call nc_check(nf90_inq_dimid(ncid, 'station', dimid=StationsDimID), &
           'WriteNetCDF', 'inquire station dimid '//trim(fname))
call nc_check(nf90_inquire_dimension(ncid, StationsDimID, len=StationsDimlen), &
           'WriteNetCDF', 'inquire station dimlen '//trim(fname))

call nc_check(nf90_inq_dimid(ncid, 'level', dimid=LevelsDimID), &
           'WriteNetCDF', 'inquire level dimid '//trim(fname))
call nc_check(nf90_inquire_dimension(ncid, LevelsDimID, len=LevelsDimlen), &
           'WriteNetCDF', 'inquire level dimlen '//trim(fname))

call nc_check(nf90_inq_dimid(ncid, 'ensemble', dimid=NmembersDimID), &
           'WriteNetCDF', 'inquire ensemble dimid '//trim(fname))
call nc_check(nf90_inquire_dimension(ncid, NmembersDimID, len=NmembersDimlen), &
           'WriteNetCDF', 'inquire ensemble dimlen '//trim(fname))

call nc_check(nf90_inq_dimid(ncid, 'forecast_lead', dimid=ForecastDimID), &
           'WriteNetCDF', 'inquire forecast_lead dimid '//trim(fname))
call nc_check(nf90_inquire_dimension(ncid, ForecastDimID, len=ForecastDimlen), &
           'WriteNetCDF', 'inquire forecast_lead dimlen '//trim(fname))

call nc_check(nf90_inq_dimid(ncid, 'analysisT', dimid=AnalysisDimID), &
           'WriteNetCDF', 'inquire analysisT dimid '//trim(fname))
call nc_check(nf90_inquire_dimension(ncid, AnalysisDimID, len=AnalysisDimlen), &
           'WriteNetCDF', 'inquire analysisT dimlen '//trim(fname))

call nc_check(nf90_inq_varid(ncid, 'analysisT', varid=TimeVarID), &
          'WriteNetCDF', 'inq_varid:analysisT '//trim(fname))

call nc_check(nf90_inq_varid(ncid, 'original_qc', varid=QCVarID), &
          'WriteNetCDF', 'inq_varid:original_qc '//trim(fname))

call nc_check(nf90_inq_varid(ncid, 'dart_qc', varid=DARTQCVarID), &
          'WriteNetCDF', 'inq_varid:dart_qc '//trim(fname))

! Increase the record dimension 

istart(1) = AnalysisDimlen + 1
icount(1) = 1

call get_time(analysisT,seconds,days)
fdays = days + seconds/86400.0_digits12

call nc_check(nf90_put_var(ncid, TimeVarId, (/ fdays /), &
            start=istart, count=icount ), 'WriteNetCDF', 'put_var:analysisT')

! Get shape of variable to ensure conformable

call nc_check(nf90_inq_varid(ncid, trim(obtype_string), varid=VarID), &
          'WriteNetCDF', 'inq_varid:obtype_string '//trim(fname))
call nc_check(nf90_inquire_variable(ncid, VarID, ndims=ndims, dimids=dimIDs), &
          'WriteNetCDF', 'inquire obtype_string dimids '//trim(fname))

if (verbose) then
   write(*,*)   ! a little whitespace
   do i = 1,ndims
      write(string1,*)'dimlen ',i,' of ',trim(obtype_string),trim(fname)
      call nc_check(nf90_inquire_dimension(ncid, dimIDs(i), len=dimLengths(i)), &
              'WriteNetCDF', string1)

      ! FIXME check the dimlength agains that expected ...  
      write(*,*)trim(obtype_string),' dimlen ',i,' is ',dimLengths(i)
   enddo
   write(*,*)   ! a little whitespace
endif

!----------------------------------------------------------------------------
! Parse the rest of the station information into arrays to stuff into netcdf
!----------------------------------------------------------------------------
! This routine gets called for every analysis time

allocate(mytimes(AnalysisDimlen))
allocate(datmat(ForecastDimlen, NmembersDimlen, CopyDimlen, LevelsDimlen))

! FIXME These are hardwired for now ... kills me
ilev           = 1
  obscopyindex = 1 
priorcopyindex = 2 
errorcopyindex = 3 

WriteObs : do stationindex = 1,num_stations

   if (debug) write(*,*)'writing station ',stationindex,' to output.'

   ! FIXME - some sort of error check on making sure the time of the station
   ! matches the time of the forecast verification 

   do ensmem = 1,NmembersDimlen
   do i = 1,ForecastDimlen

      ! observation value
      datmat(i,ensmem,  obscopyindex,ilev) = stations(stationindex)%observation(i) 
      ! ensemble prior observation value
      datmat(i,ensmem,priorcopyindex,ilev) = stations(stationindex)%forecast(ensmem,i) 
      ! observation error variance
      datmat(i,ensmem,errorcopyindex,ilev) = stations(stationindex)%obserror(i) 

   enddo
   enddo

   call nc_check(nf90_put_var(ncid, VarId, datmat, &
        start=(/ 1, 1, 1, ilev, stationindex, istart(1) /), &
        count=(/ ForecastDimlen, NmembersDimlen, CopyDimlen, LevelsDimlen, 1, 1 /) ), &
        'WriteNetCDF', 'put_var:times')

   ! Write QC values

   call nc_check(nf90_put_var(ncid, QCVarId, stations(stationindex)%orgqc, &
           start=(/ 1, stationindex, istart(1) /), count=(/ ForecastDimlen, 1, 1 /)), &
           'WriteNetCDF', 'put_var:orgqcmat')

   call nc_check(nf90_put_var(ncid, DARTQCVarId, stations(stationindex)%dartqc, &
           start=(/ 1, stationindex, istart(1) /), count=(/ ForecastDimlen, 1, 1 /)), &
           'WriteNetCDF', 'put_var:orgqcmat')

enddo WriteObs

deallocate(mytimes)
deallocate(datmat)

!----------------------------------------------------------------------------
! finished ...
!----------------------------------------------------------------------------

call nc_check(nf90_sync( ncid), 'WriteNetCDF', 'sync '//trim(fname))  

if (debug) write(*,*)'DEBUG --- leaving WriteNetCDF'

end subroutine WriteNetCDF


!============================================================================


subroutine CloseNetCDF(ncid, fname)
!----------------------------------------------------------------------------
integer,          intent(in) :: ncid
character(len=*), intent(in) :: fname

integer :: indx1

if ( debug ) write(*,*)'DEBUG --- Closing ',trim(fname)

! Enter define mode (again) to record the observation sequence
! files used. For expediency, should learn how to allocate some
! extra space for this ... possible FIXME

call nc_check(nf90_redef(ncid), 'CloseNetCDF', 'redef '//trim(fname))  

! write all observation sequence files used
FILEloop : do i = 1,SIZE(obs_seq_filenames)

  indx1 = index(obs_seq_filenames(i),'null')

  if (indx1 > 0) exit FILEloop

  write(string1,'(''obs_seq_file_'',i3.3)')i
  call nc_check(nf90_put_att(ncid, NF90_GLOBAL, &
         trim(string1), trim(obs_seq_filenames(i)) ), &
         'CloseNetCDF', 'put_att:filenames')

enddo FILEloop

call nc_check(nf90_enddef(ncid), 'CloseNetCDF', 'enddef '//trim(fname))  
call nc_check(nf90_sync(  ncid), 'CloseNetCDF', 'sync   '//trim(fname))  
call nc_check(nf90_close( ncid), 'CloseNetCDF', 'close  '//trim(fname))  

end subroutine CloseNetCDF


!============================================================================


subroutine destroy_stations(stations)
!----------------------------------------------------------------------------
type(station), allocatable, dimension(:), intent(inout) :: stations

integer :: i,N

N = size(stations)

do i = 1,N
   if (associated(stations(i)%times)) then
      deallocate( stations(i)%times )
      nullify(    stations(i)%times )
   endif
enddo

if (allocated(stations)) deallocate(stations)

end subroutine destroy_stations


!============================================================================


function find_ensemble_size()
!----------------------------------------------------------------------------
! The only way I know of to determine the ensemble size is to read the
! copy metadata from one of the obs_seq.xxx files. I'm just going to 
! open the first file, read the entire damn sequence, and count them up. 
!
! Then, I need to preserve all the indices for the observation value
! as well as the indices for the ensemble members. The observation error
! variance comes from get_obs_def_error_variance()

integer :: find_ensemble_size

integer :: filenum = 1

if (obs_sequence_list == '') then
   obs_seq_in_file_name = next_file(obs_sequence_name,filenum)
else
   obs_seq_in_file_name = get_next_filename(obs_sequence_list,filenum)
endif

if ( file_exist(trim(obs_seq_in_file_name)) ) then
   write(string1,*)'opening ', trim(obs_seq_in_file_name)
   call error_handler(E_MSG,'find_ensemble_size:',string1,source,revision,revdate)
else
   write(string1,*)trim(obs_seq_in_file_name), &
                   ' does not exist. Dying a dramatic death.'
   call error_handler(E_ERR,'find_ensemble_size:',string1,source,revision,revdate)
endif

!-----------------------------------------------------------------------
! Read the entire observation sequence - allocates 'seq' internally
! And then cruise throught the metadata 
!-----------------------------------------------------------------------

call read_obs_seq(obs_seq_in_file_name, 0, 0, 0, seq)

find_ensemble_size = 0

do i=1, get_num_copies(seq)
   if( index(get_copy_meta_data(seq,i), 'prior ensemble member') > 0) &
      find_ensemble_size = find_ensemble_size + 1
enddo

write(string1,'(''There are '',i4,'' ensemble members.'')') find_ensemble_size
if (verbose) call error_handler(E_MSG,'find_ensemble_size:',string1,source,revision,revdate)

if (find_ensemble_size < 1) then
   write(string1,*)'no ensemble member info in ', trim(obs_seq_in_file_name)
   write(string2,*)'Unable to continue.'
   call error_handler(E_ERR,'find_ensemble_size', string1, &
                 source, revision, revdate, text2=string2)
endif

! call destroy_obs_sequence(seq) gets destroyed at top of loop around input files

end function find_ensemble_size


!============================================================================


subroutine find_analysis_time(filename, ExpIndex, analysis_time)
!----------------------------------------------------------------------------
! Parse the filename into an integer string that defines the 
! analysis time.
!

character(len=*), intent(in)  :: filename
integer,          intent(out) :: ExpIndex
type(time_type),  intent(out) :: analysis_time

character(len=LEN_TRIM(filename)) :: lj_filename

integer :: iyear,imonth,iday,ihour
integer :: indx, j, mystat
type(time_type) :: mytime

ExpIndex    = 0
lj_filename = adjustl(filename)

! Find the start of the filename extension, the YYYYMMDDHH string

indx = index(lj_filename,'.',back=.TRUE.) + 1
string2 = 'The extension must be the analysis time - YYYYMMDDHH'

if (indx == 1) then
   write(string1,*)'Cannot find the YYYYMMDDHH extension on '//trim(lj_filename)
   call error_handler(E_ERR,'find_analysis_time', string1, &
                 source, revision, revdate, text2=string2)
endif

read(lj_filename(indx+0:indx+3),'(i4)',iostat=mystat)iyear
if (mystat /= 0) then
   write(string1,*)'Cannot find the year in the extension (YYYYMMDDHH) on '//trim(lj_filename)
   call error_handler(E_ERR,'find_analysis_time', string1, &
                 source, revision, revdate, text2=string2)
endif

read(lj_filename(indx+4:indx+5),'(i2)',iostat=mystat)imonth
if (mystat /= 0) then
   write(string1,*)'Cannot find the month in the extension (YYYYMMDDHH) on '//trim(lj_filename)
   call error_handler(E_ERR,'find_analysis_time', string1, &
                 source, revision, revdate, text2=string2)
endif

read(lj_filename(indx+6:indx+7),'(i2)',iostat=mystat)iday
if (mystat /= 0) then
   write(string1,*)'Cannot find the day in the extension (YYYYMMDDHH) on '//trim(lj_filename)
   call error_handler(E_ERR,'find_analysis_time', string1, &
                 source, revision, revdate, text2=string2)
endif

read(lj_filename(indx+8:indx+9),'(i2)',iostat=mystat)ihour
if (mystat /= 0) then
   write(string1,*)'Cannot find the hour in the extension (YYYYMMDDHH) on '//trim(lj_filename)
   call error_handler(E_ERR,'find_analysis_time', string1, &
                 source, revision, revdate, text2=string2)
endif

if (debug) then
   write(*,*)'iyear  is ',iyear
   write(*,*)'imonth is ',imonth
   write(*,*)'iday   is ',iday
   write(*,*)'ihour  is ',ihour
endif

mytime = set_date(iyear, imonth, iday, hours=ihour)

TimeLoop : do j = 1,size(VerifyTimes,1)
   if ( mytime == VerifyTimes(j,1) ) then
      analysis_time = mytime 
      ExpIndex      = j
      exit TimeLoop
   endif
enddo TimeLoop

if ( ExpIndex == 0 ) then
   write(string1,*)'Cannot find an analysis time that matches '//lj_filename(indx+0:indx+9)
   call error_handler(E_ERR,'find_analysis_time',string1,source,revision,revdate)
endif
 
end subroutine find_analysis_time


!============================================================================


subroutine find_our_copies(myseq, obsindex, indices, qcindex, dart_qcindex)
!----------------------------------------------------------------------------
!
type(obs_sequence_type), intent(in)  :: myseq
integer,                 intent(out) :: obsindex
integer, dimension(:),   intent(out) :: indices
integer,                 intent(out) :: qcindex
integer,                 intent(out) :: dart_qcindex

integer :: i, myindex

character(len=obstypelength) :: metadata

    obsindex = -1
     qcindex = -1
dart_qcindex = -1
indices      = -1

! Find the observation copy 
ObsDataLoop : do i=1, get_num_copies(myseq)

   metadata = get_copy_meta_data(myseq,i)
   if( index(metadata, 'observation') > 0) then  ! FIXME - not robust
      obsindex = i
      exit ObsDataLoop
   endif

enddo ObsDataLoop

if (obsindex < 0) then
   write(string1,*)'Could not find observation copy index in sequence.'
   call error_handler(E_ERR,'find_our_copies',string1,source,revision,revdate)
endif

myindex = 0
! All prior copies get tacked one after another ...
MetaDataLoop : do i=1, get_num_copies(myseq)

   metadata = get_copy_meta_data(myseq,i)

   if( index(metadata, 'prior ensemble member') > 0) then
      myindex = myindex + 1

      if ( myindex > size(indices) ) then
         write(string1,*)'Found too many prior copies in metadata.'
         call error_handler(E_ERR,'find_our_copies',string1,source,revision,revdate)
      endif
      indices(myindex) = i
   endif

enddo MetaDataLoop

if (myindex == 0) then
   write(string1,*)'Could not find any prior copy indices in sequence.'
   call error_handler(E_ERR,'find_our_copies',string1,source,revision,revdate)
endif

! Here is the sanity check to make sure the indices pick off just the
! prior ensemble members.

if (verbose) then
   write(*,*)   ! just a little whitespace
   do i = 1,size(indices)
      metadata = get_copy_meta_data(myseq,indices(i))
      write(*,'('' index '',i4,'' is '',a)')indices(i),trim(metadata)
   enddo
   write(*,*)   ! just a little whitespace
endif

! Sometimes the first QC field is 'Quality Control' or 'NCEP QC index'
! to me ... they mean the same thing.

QCMetaDataLoop : do i=1, get_num_qc(myseq)
   if(index(get_qc_meta_data(myseq,i),'Quality Control'     ) > 0)      qcindex = i
   if(index(get_qc_meta_data(myseq,i),'NCEP QC index'       ) > 0)      qcindex = i
   if(index(get_qc_meta_data(myseq,i),'DART quality control') > 0) dart_qcindex = i
enddo QCMetaDataLoop

if (      qcindex < 0 ) then 
   write(string1,*)'metadata:Quality Control copyindex not found' 
   call error_handler(E_MSG,'find_qc_indices:',string1,source,revision,revdate)
endif
if ( dart_qcindex < 0 ) then 
   write(string1,*)'metadata:DART quality control copyindex not found' 
   call error_handler(E_MSG,'find_qc_indices:',string1,source,revision,revdate)
endif

! Just echo what we know

if (verbose) then
   write(*,*)'QC index',     qcindex,' ',trim(get_qc_meta_data(myseq,     qcindex))
   write(*,*)'QC index',dart_qcindex,' ',trim(get_qc_meta_data(myseq,dart_qcindex))
   write(*,*)
endif

end subroutine find_our_copies


!============================================================================


end program obs_seq_verify

