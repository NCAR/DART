! DART software - Copyright UCAR. This open source software is provided
! by UCAR, "as is", without charge, subject to all terms of use at
! http://www.image.ucar.edu/DAReS/DART/DART_download

!> This program creates a netCDF file suitable for forecast evaluation.
!>
!> (analysisT, station, level, copy, ensemble, forecast_lead)
!>      |         |       |      |      |         |
!>      |         |       |      |      |         +- forecast length : 0,3,6,9,12,...
!>      |         |       |      |      |
!>      |         |       |      |      +----------- ensemble member index
!>      |         |       |      |
!>      |         |       |      +------------------ obs value, prior, obs_err
!>      |         |       |
!>      |         |       +------------------------- vertical level index
!>      |         |
!>      |         +--------------------------------- (horizontal) station index
!>      |
!>      +------------------------------------------- analysis time/date
!>
!> I think the logic of the program should be as follows:
!> 1) The list of 'stations' is read from a netCDF file -- the product of obs_seq_coverage.f90
!> 2) The set of forecast lead times must be determined. Having 8 might not be
!>    specific enough ... we want 8 separated by 3 hours ... for example.
!> 3) A series of obs_seq.fcst files (each from one analysisT and containing multiple
!>    forecast lead times) is read and stuffed into an appropriate structure.
!> 4) The structure is written into the netCDF file.
!> 5) On to the next obs_seq.fcst file ... (step 3)
!> 6) wrap up ...
!>
!> Soyoung's (original) wish list - which has been subsequently modified:
!> "Now I'm done running filter for 24-hr forecast with 3-hrly observations as an
!> evaluation mode only, and am ready to hand obs_seq.final over to you for the final
!> conversion process.
!>
!> Ideally, in obs_seq_fcst.nc (if I can name it on my own), I would like to have a
!> data structure of (copy, station, level, ensemble, date, time) for each variable
!> and each obs type, where copy is (observation value, prior observation value
!> corresponding to the observation, obs error standard deviation), date is the
!> number of cycle, and time is the number of forecast lead times.
!> So, the dimension "date" is 1 for my current single obs_seq.final, but reserved
!> for multiple initialization cases (i.e., multiple obs_seq.final files) with the
!> same length of forecast lead times. This dimension is the most critical one in the
!> statistical verification since this dimension number is the actual sample size to
!> tell if we have enough samples to make the verification result statistically
!> significant or not."

program obs_seq_verify

use        types_mod, only : r4, r8, digits12, MISSING_R8, MISSING_I, &
                             metadatalength, obstypelength

use obs_sequence_mod, only : read_obs_seq, obs_type, obs_sequence_type, get_first_obs, &
                             get_obs_def, get_copy_meta_data, get_obs_values, &
                             get_next_obs, init_obs, init_obs_sequence, &
                             assignment(=), get_num_copies, get_num_qc, get_qc, &
                             static_init_obs_sequence, destroy_obs_sequence, destroy_obs, &
                             read_obs_seq_header, get_qc_meta_data

use      obs_def_mod, only : obs_def_type, get_obs_def_time, get_obs_def_type_of_obs, write_obs_def, &
                             get_obs_def_location, set_obs_def_time, &
                             set_obs_def_location, set_obs_def_type_of_obs, &
                             set_obs_def_error_variance, get_obs_def_error_variance

use     obs_kind_mod, only : max_defined_types_of_obs, get_name_for_type_of_obs, get_index_for_type_of_obs, &
                             write_type_of_obs_table

use     location_mod, only : location_type, get_location, set_location_missing, &
                             write_location, operator(/=), operator(==), &
                             set_location, is_location_in_region, query_location, &
                             VERTISUNDEF, is_vertical

use  location_io_mod, only : nc_write_location_atts, nc_write_location

use time_manager_mod, only : time_type, set_date, set_time, get_time, print_time, &
                             set_time_missing, print_date, set_calendar_type, &
                             operator(>), operator(<), operator(==), &
                             operator(<=), operator(-), operator(+), operator(/=)

use    utilities_mod, only : get_unit, close_file, timestamp, &
                             file_exist, error_handler, E_ERR, E_WARN, E_MSG, &
                             initialize_utilities, finalize_utilities, nmlfileunit, &
                             find_namelist_in_file, check_namelist_read, &
                             next_file, set_filename_list, find_textfile_dims, &
                             file_to_text, do_nml_file, do_nml_term

use  netcdf_utilities_mod, only : nc_check

use typeSizes
use netcdf

implicit none

character(len=*), parameter :: source = 'obs_seq_verify.f90'

!---------------------------------------------------------------------
! An array of these structures will be filled for each obs sequence file.
! After each obs sequence file has been read, the information will be
! written to the appropriate slots in the netCDF file.
!---------------------------------------------------------------------

type voxel_type
   integer                  :: good
   integer                  :: obs_type
   type(time_type)          :: first_time
   type(time_type)          :: last_time
   type(time_type), pointer :: times(:)     ! actual obs time closest to nominal
   type(location_type)      :: location
   real(r8), dimension(3)   :: LonLatLev    ! easy access to location components
   integer                  :: station_id   ! index of unique horizontal location
   integer                  :: levelindex   ! index of mandatory level
end type voxel_type

type station_type
   integer                  :: obs_type
   type(location_type)      :: location
end type station_type

type nc6Dvar
   character(len=NF90_MAX_NAME) :: AnalysisDimName = 'analysisT'
   integer                      :: AnalysisDimID
   integer                      :: AnalysisDimLen
   character(len=NF90_MAX_NAME) :: CopyDimName = 'copy'
   integer                      :: CopyDimID
   integer                      :: CopyDimLen
   character(len=NF90_MAX_NAME) :: StationsDimName = 'station'
   integer                      :: StationsDimID
   integer                      :: StationsDimLen
   character(len=NF90_MAX_NAME) :: LevelsDimName = 'level'
   integer                      :: LevelsDimID
   integer                      :: LevelsDimLen
   character(len=NF90_MAX_NAME) :: EnsembleDimName = 'ensemble'
   integer                      :: EnsembleDimID
   integer                      :: EnsembleDimLen
   character(len=NF90_MAX_NAME) :: ForecastDimName = 'forecast_lead'
   integer                      :: ForecastDimID
   integer                      :: ForecastDimLen
end type nc6Dvar

logical,            allocatable, dimension(:) :: DesiredStations
type(station_type), allocatable, dimension(:) :: station
type(voxel_type),   allocatable, dimension(:) :: voxel
type(nc6Dvar) :: ncmeta

integer :: ensemble_size   ! the # of ensemble members in the obs_seq files
integer :: num_voxels      ! The number of unique locations
integer :: num_stations    ! The number of unique horizontal locations
integer :: num_levels      ! number of mandatory levels
integer :: num_verif_times ! number of all possible verification times.
integer :: num_forecasts   ! number of forecasts 
integer, parameter :: NCOPIES = 3 ! obs value, prior, obs error variance

character(len=metadatalength), dimension(NCOPIES) :: copy_metadata = &
     (/ 'observation               ', &
        'forecast                  ', &
        'observation error variance' /)

!---------------------------------------------------------------------
! variables associated with the observation
!---------------------------------------------------------------------

type(obs_sequence_type) :: seq
type(obs_type)          :: obs1, obs2
type(obs_def_type)      :: obs_def
type(location_type)     :: obs_loc

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

integer, parameter :: MAX_NUM_INPUT_FILES = 500
integer            :: num_input_files
character(len=256) :: obs_seq_in_file_name

!-----------------------------------------------------------------------
! Namelist with (some scalar) default values
!-----------------------------------------------------------------------

character(len=256) :: obs_sequences(MAX_NUM_INPUT_FILES) = ''
character(len=256) :: obs_sequence_list = ''
character(len=256) :: input_template    = 'obsdef_mask.nc'
character(len=256) :: netcdf_out        = 'forecast.nc'
character(len=129) :: calendar          = 'gregorian'
character(len=obstypelength) :: obtype_string

integer :: print_every = 10000
logical :: verbose     = .true.
logical :: debug       = .false.

namelist /obs_seq_verify_nml/ obs_sequences, obs_sequence_list, &
                             input_template, netcdf_out, &
                             obtype_string, calendar, &
                             print_every, verbose, debug

!-----------------------------------------------------------------------
! Quantities of interest
!-----------------------------------------------------------------------

integer ::       qc_index   ! copy index of the original qc value
integer ::  dart_qc_index   ! copy index of the DART qc value
integer :: obs_copy_index   ! copy index of the observation
integer :: voxel_id
integer :: fcst_lead_index  ! forecast lead time index
integer :: AnalysisIndex    ! index of the forecast experiment
type(time_type) :: analysisT ! valid time of analysis at start of forecast
character(len=metadatalength), dimension(:), allocatable :: module_obs_copy_names
integer,                       dimension(:), allocatable :: prior_copy_indices
integer,                       dimension(:), allocatable :: forecast_leads
real(r8),                      dimension(:), allocatable :: mandatory_level
type(time_type),             dimension(:,:), allocatable :: VerifyTimes
real(r8) :: obs_error_variance

!-----------------------------------------------------------------------
! General purpose variables
!-----------------------------------------------------------------------

integer  :: ifile, iobs, ngood
integer  :: i, io, ncunit
integer, dimension(6) :: mystart, mycount

type(time_type) :: obs_time

character(len=256) :: ncName
character(len=512) :: string1, string2, string3

! ~# of degrees for 1/2 meter at Earth equator
! 360 deg-earth/(40000 km-earth * 1000m-km)
real(r8), parameter :: HALF_METER = 180.0_r8 / (40000.0_r8 * 1000.0_r8)
real(r8), parameter :: OnePa = 1.0_r8

!=======================================================================
! Get the party started
!=======================================================================

call initialize_utilities('obs_seq_verify')
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

num_input_files =  set_filename_list(obs_sequences, obs_sequence_list,'obs_seq_verify')

if (debug) then
   write(*,*)'There are ',num_input_files,' input observation sequence files.'
   do ifile = 1,num_input_files
      write(*,*)'file      ',ifile,' is ['//trim(obs_sequences(ifile))//']'
   enddo 
endif

string1 = adjustl(obtype_string)
string2 = adjustl(netcdf_out)
obtype_string  = trim(string1)
netcdf_out     = trim(string2)
obtype_integer = get_index_for_type_of_obs(obtype_string)

if (obtype_integer < 1) then
   write(string1,*)'obtype_string ',trim(obtype_string),' is unknown. change input.nml'
   call error_handler(E_ERR,'obs_seq_verify:',string1,source)
endif

call set_calendar_type(calendar)

! Determine the number of voxels from the input netcdf file.
! Initialize the output netcdf file. Must know the ensemble size.
! Determine which voxels share the same horizontal location.
! Each voxel already has a vertical level index associated with it.

call find_ensemble_size() ! from the obs_seq copy metadata of the first file
num_voxels   = fill_voxels( input_template )
num_stations = relate_voxels_to_stations()

allocate(prior_copy_indices(ensemble_size), module_obs_copy_names(ensemble_size))

ncName = trim(obtype_string)//'_'//trim(netcdf_out)

ncunit = InitNetCDF(trim(ncName))

!----------------------------------------------------------------------
! Prepare the variables
!----------------------------------------------------------------------

ObsFileLoop : do ifile=1,num_input_files
!-----------------------------------------------------------------------
! By construction - 
! each file contains all the observations for an entire forecast.
! The file NAME specifies the analysis time from which the forecast is started.

   ! Because of the ability to 'cycle' the ObsFileLoop, we need to
   ! destroy and deallocate at the top of the loop.

   if (ifile /= 1) then
      call destroy_obs(obs1)
      call destroy_obs(obs2)
      call destroy_obs_sequence(seq)
   endif

   if (allocated(   copy_values)) deallocate(   copy_values)
   if (allocated(     qc_values)) deallocate(     qc_values)
   if (allocated( qc_copy_names)) deallocate( qc_copy_names)
   if (allocated(obs_copy_names)) deallocate(obs_copy_names)

   ! Determine the next input filename and check for existence.

   obs_seq_in_file_name = obs_sequences(ifile)

   if ( file_exist(obs_seq_in_file_name) ) then
      write(*,*) ! whitespace
      write(*,*) ! whitespace
      write(string1,*)'opening ', trim(obs_seq_in_file_name)
      call error_handler(E_MSG,'obs_seq_verify:',string1)
   else
      write(string1,*)'['//trim(obs_seq_in_file_name)//'] does not exist.'
      call error_handler(E_ERR,'obs_seq_verify:',string1,source)
      exit ObsFileLoop
   endif

   ! Determine the analysis time from the file name

   call find_analysis_time(obs_seq_in_file_name, AnalysisIndex, analysisT)

   ! Read in information about observation sequence so we can allocate
   ! observations. We need info about how many copies, qc values, etc.

   call read_obs_seq_header(obs_seq_in_file_name, &
             num_copies, num_qc, num_obs, max_num_obs, &
             obs_seq_file_id, obs_seq_read_format, pre_I_format, &
             close_the_file = .true.)

   if ((num_qc <= 0) .or. (num_copies <= 0)) then
      write(string1,*)'need at least 1 qc and 1 observation copy'
      call error_handler(E_ERR,'obs_seq_verify:',string1,source)
   endif

   ! Initialize some (individual) observation variables

   call init_obs(obs1, num_copies, num_qc) ! First obs in sequence
   call init_obs(obs2, num_copies, num_qc)

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

   !--------------------------------------------------------------------
   ! * Read the entire observation sequence - allocates 'seq' internally
   ! * Find the N copy indices we want
   ! * Do something with the QC values?
   !     - I suspect they will all be 'evaluate only'
   !     - what if the prior forward operator fails?
   !--------------------------------------------------------------------

   call read_obs_seq(obs_seq_in_file_name, 0, 0, 0, seq)

   do i=1, num_copies
      obs_copy_names(i) = get_copy_meta_data(seq,i)
   enddo
   do i=1, num_qc
      qc_copy_names(i) = get_qc_meta_data(seq,i)
   enddo

   ! Determine which qc copy is original qc and dart qc

   call find_our_copies(seq, obs_copy_index, prior_copy_indices)

   ! The first trip through sets the module_obs_copy_names so we can
   ! be sure we are stuffing compatible objects into the same slots
   if ( ifile == 1 ) then
      module_obs_copy_names = obs_copy_names(prior_copy_indices(1:ensemble_size))
   else
      ! Check to make sure the ensemble members are in the expected copies
      do i = 1,ensemble_size
         if ( obs_copy_names(prior_copy_indices(i)) /= module_obs_copy_names(i) ) then
            write(string1,'(''module has '',A)') trim(module_obs_copy_names(i))
            write(string2,'(A,'' has '',A)') trim(obs_seq_in_file_name), &
                                             trim(obs_copy_names(prior_copy_indices(i)))
            call error_handler(E_ERR,'obs_seq_verify:', &
               'mismatch in observation copies',source,text2=string1,text3=string2)
         endif
      enddo
   endif

   ngood = 0

   !--------------------------------------------------------------------
   ObservationLoop : do iobs = 1,num_obs
   !--------------------------------------------------------------------

      if (iobs == 1) then
         if ( .not. get_first_obs(seq, obs1) )           &
            call error_handler(E_ERR,'obs_seq_verify:',   &
              'No first observation in sequence.', source)
      else
         call get_next_obs(seq, obs1, obs2, last_ob_flag)
         obs1 = obs2
      endif

      if ( verbose .and. (mod(iobs,print_every) == 0) ) &
         write(*,*)'Processing obs ',iobs,' of ',num_obs

      call get_obs_def(obs1,          obs_def)
      flavor   = get_obs_def_type_of_obs(        obs_def)
      obs_time = get_obs_def_time(    obs_def)
      obs_loc  = get_obs_def_location(obs_def)
      call get_qc(obs1, qc_values)

   !  if (debug .and. any(iobs == (/1, 295, 1908, 6265/) )) then
      if (debug .and. (iobs == 1)) then
         write(*,*)'looking for observation types of ',obtype_integer, &
                                     get_name_for_type_of_obs(obtype_integer)

         write(*,*)'First observation "type" happens to be ',flavor, &
                                get_name_for_type_of_obs(flavor)
         call print_time(obs_time,'First observation time')
         call print_date(obs_time,'First observation date')
         call write_location(6,obs_loc,fform='ascii')
         call write_location(6,obs_loc,fform='ascii',charstring=string1)
         write(*,*)trim(string1)

         write(*,*)'     qc_index and value is ',     qc_index, qc_values(     qc_index)
         write(*,*)'dart_qc_index and value is ',dart_qc_index, qc_values(dart_qc_index)

         call get_obs_values(obs1, copy_values)
         write(*,*)'observation values are:'
         do i = 1,size(copy_values)
            write(*,*)i,copy_values(i)
         enddo
         write(*,*)'observation QC values are:'
         do i = 1,size(qc_values)
            write(*,*)i,qc_values(i)
         enddo

      endif

      ! 0) is it the right type of observation
      ! 1) what/if any station does it belong to?
      ! 2) Is it one of the times of interest?
      ! 3) stuff it into the appropriate station structure

      ! FIXME ... What do we do with DART QCs of 7 
      if ( flavor /= obtype_integer)     cycle ObservationLoop
      if (qc_values(dart_qc_index) >= 4) cycle ObservationLoop

      call match_obs_to_voxel(iobs, flavor, obs_loc, obs_time, voxel_id, fcst_lead_index)
      if ( fcst_lead_index < 1 )            cycle ObservationLoop

      call get_obs_values(obs1, copy_values)
      obs_error_variance = get_obs_def_error_variance(obs_def)

   !  if (debug) then ! GIGANTIC OUTPUT - BE CAREFUL
      if (debug .and. (voxel(voxel_id)%station_id == 1)) then ! LESS GIGANTIC OUTPUT
         do i = 1,size(prior_copy_indices)
            write(*,'('' DEBUG ob '',i10,3(1x,i6),3(1x,f17.7))') &
            iobs, i, fcst_lead_index, prior_copy_indices(i), &
            obs_error_variance, &
            copy_values(obs_copy_index), &
            copy_values(prior_copy_indices(i))
         enddo
         write(*,*) ! a little whitespace
      endif

      ! Determine where to stuff it into the output netCDF file.
      ! YES ... one observation at a time.
      ! FIXME allocate a potentially LARGE array instead of writing one-at-a-time
      ! FIXME ... one for each Analysis Time

      call determine_hyperslab_indices(AnalysisIndex, fcst_lead_index, &
              voxel_id, mystart, mycount)

      call WriteNetCDF(ncunit, trim(ncName), copy_values, obs_copy_index, &
              prior_copy_indices, obs_error_variance, qc_values(qc_index), &
              qc_values(dart_qc_index), mystart, mycount)

   !--------------------------------------------------------------------
   enddo ObservationLoop
   !--------------------------------------------------------------------

enddo ObsFileLoop

call CloseNetCDF(ncunit, ncName)

!-----------------------------------------------------------------------
! Really, really, done.
!-----------------------------------------------------------------------

call destroy_obs(obs1)
call destroy_obs(obs2)
call destroy_obs_sequence(seq)
! call destroy_stations(station)
call destroy_voxels()

if (allocated(qc_values))             deallocate(qc_values)
if (allocated(qc_copy_names))         deallocate(qc_copy_names)
if (allocated(copy_values))           deallocate(copy_values)
if (allocated(obs_copy_names))        deallocate(obs_copy_names)
if (allocated(DesiredStations))       deallocate(DesiredStations)

call error_handler(E_MSG,'obs_seq_verify:','Finished successfully.')
call finalize_utilities()



!======================================================================
CONTAINS
!======================================================================



function fill_voxels( filename ) result(nvoxels)
!----------------------------------------------------------------------------
! Check to make sure the file exists.
! Read the 'voxel' variable (binary value for desired or not).
! count the number of non-zero entries ...
!
! Since the 'voxels' coming in are assumed vertically independent, but the output
! structure 'station' has a vertical component ... we must determine which input
! 'voxels' are actually part of a station.
!-----------------------------------------------------------------------

character(len=*), intent(in) :: filename
integer                      :: nvoxels

integer :: DimID, VarID
integer :: nanalyses, locNdim, strlen
integer :: ncid, i, j, ivoxel, ndims, mylen
integer :: seconds, days
type(location_type) :: myloc

integer, dimension(nf90_max_var_dims) :: dimIDs

! These mimic the contents of the netCDF file
character(len=obstypelength), dimension(:), allocatable :: obs_types
integer,                      dimension(:), allocatable :: voxel_flag
integer,                      dimension(:), allocatable :: which_vert
integer,                      dimension(:), allocatable :: voxel_level_index
real(digits12),               dimension(:), allocatable :: first_time
real(digits12),               dimension(:), allocatable :: last_time
real(digits12),             dimension(:,:), allocatable :: ReportTimes
real(digits12),             dimension(:,:), allocatable :: ExperimentTimes
real(r8),                   dimension(:,:), allocatable :: locations

if (verbose) then
   call timestamp(string1='fill_voxels: reading the voxel information at ',pos='brief')
endif

!-----------------------------------------------------------------------
! Open the netCDF file

call nc_check(nf90_open(trim(filename), NF90_NOWRITE, &
         ncid = ncid), 'fill_voxels', 'open '//trim(filename))

!-----------------------------------------------------------------------
! Determine the maximum number of stations from the station dimension,
! the number of times, and the number of dimensions in the location

! used for integer array indicating yes/no - do we want the station
call nc_check(nf90_inq_dimid(ncid, 'voxel', dimid=DimID), &
           'fill_voxels', 'inquire voxel dimid '//trim(filename))
call nc_check(nf90_inquire_dimension(ncid, DimID, len=nvoxels), &
           'fill_voxels', 'inquire voxel nvoxels '//trim(filename))

! used for array of all possible verification times needed
call nc_check(nf90_inq_dimid(ncid, 'time', dimid=DimID), &
           'fill_voxels', 'inquire time dimid '//trim(filename))
call nc_check(nf90_inquire_dimension(ncid, DimID, len=num_verif_times), &
           'fill_voxels', 'inquire time len '//trim(filename))

! used for array of the analysis times / 'zero-length' forecast times
call nc_check(nf90_inq_dimid(ncid, 'analysisT', dimid=DimID), &
           'fill_voxels', 'inquire analysisT dimid '//trim(filename))
call nc_check(nf90_inquire_dimension(ncid, DimID, len=nanalyses), &
           'fill_voxels', 'inquire analysisT len '//trim(filename))

! used for array of the verification times
call nc_check(nf90_inq_dimid(ncid, 'forecast_lead', dimid=DimID), &
           'fill_voxels', 'inquire forecast_lead dimid '//trim(filename))
call nc_check(nf90_inquire_dimension(ncid, DimID, len=num_forecasts), &
           'fill_voxels', 'inquire forecast_lead len '//trim(filename))

call nc_check(nf90_inq_dimid(ncid, 'stringlength', dimid=DimID), &
           'fill_voxels', 'inquire stringlength dimid '//trim(filename))
call nc_check(nf90_inquire_dimension(ncid, DimID, len=strlen), &
           'fill_voxels', 'inquire stringlength len '//trim(filename))

call nc_check(nf90_inq_dimid(ncid, 'location', dimid=DimID), &
           'fill_voxels', 'inquire location dimid '//trim(filename))
call nc_check(nf90_inquire_dimension(ncid, DimID, len=locNdim), &
           'fill_voxels', 'inquire location len '//trim(filename))

call nc_check(nf90_inq_dimid(ncid, 'nlevels', dimid=DimID), &
           'fill_voxels', 'inquire nlevels dimid '//trim(filename))
call nc_check(nf90_inquire_dimension(ncid, DimID, len=num_levels), &
           'fill_voxels', 'inquire nlevels len '//trim(filename))

!-----------------------------------------------------------------------
! Read the "voxel" 1D variable. Each non-zero entry indicates a
! voxel we wanted to keep when we ran obs_seq_coverage.f90.

call nc_check(nf90_inq_varid(ncid, 'voxel', varid=VarID), &
          'fill_voxels', 'inq_varid:voxel '//trim(filename))
call nc_check(nf90_inquire_variable(ncid, VarID, ndims=ndims, dimids=dimIDs), &
          'fill_voxels', 'inquire voxel ndims '//trim(filename))

if (ndims /= 1) then
   write(string1,*)'voxel is expected to be a 1D array, it is ',ndims
   call error_handler(E_MSG,'fill_voxels:',string1)
endif

call nc_check(nf90_inquire_dimension(ncid, dimIDs(1), len=mylen), &
           'fill_voxels', 'voxel inquire dimid(1) '//trim(filename))

if (mylen /= nvoxels) then
   write(string1,*)'voxel length expected to be ',nvoxels,' it is ',mylen
   call error_handler(E_MSG,'fill_voxels:',string1)
endif

allocate( voxel_flag(nvoxels),  obs_types(nvoxels),        which_vert(nvoxels), &
          first_time(nvoxels),  last_time(nvoxels), voxel_level_index(nvoxels))
allocate(locations(locNdim,nvoxels))
allocate(ReportTimes(num_verif_times,nvoxels))     ! observation times that are closest
allocate(ExperimentTimes(num_forecasts,nanalyses)) ! verification times of interest

! The following are all GLOBAL scope

allocate(VerifyTimes(nanalyses,num_forecasts)) ! verification times of interest
allocate(mandatory_level(num_levels))          ! mandatory pressure levels
allocate(forecast_leads(num_forecasts))        ! how far into the forecast (seconds)
allocate(voxel(nvoxels))

!-----------------------------------------------------------------------
! Read all the information from the template netCDF file.

call nc_check(nf90_get_var(ncid, VarID, voxel_flag), &
         'fill_voxels',  'get_var:voxel_flag '//trim(filename))

call nc_check(nf90_inq_varid(ncid, 'obs_type', varid=VarID), &
        'fill_voxels', 'inq_varid:obs_type '//trim(filename))
call nc_check(nf90_get_var(ncid, VarID, obs_types), &
        'fill_voxels',   'get_var:obs_type '//trim(filename))

call nc_check(nf90_inq_varid(ncid, 'location', varid=VarID), &
        'fill_voxels', 'inq_varid:location '//trim(filename))
call nc_check(nf90_get_var(ncid, VarID, locations), &
        'fill_voxels',   'get_var:location '//trim(filename))

call nc_check(nf90_inq_varid(ncid, 'mandatory_level', varid=VarID), &
        'fill_voxels', 'inq_varid:mandatory_level '//trim(filename))
call nc_check(nf90_get_var(ncid, VarID, mandatory_level), &
        'fill_voxels',   'get_var:mandatory_level'//trim(filename))

call nc_check(nf90_inq_varid(ncid, 'which_vert', varid=VarID), &
        'fill_voxels', 'inq_varid:which_vert '//trim(filename))
call nc_check(nf90_get_var(ncid, VarID, which_vert), &
        'fill_voxels',   'get_var:which_vert '//trim(filename))

call nc_check(nf90_inq_varid(ncid, 'voxel_level_index', varid=VarID), &
        'fill_voxels', 'inq_varid:voxel_level_index '//trim(filename))
call nc_check(nf90_get_var(ncid, VarID, voxel_level_index), &
        'fill_voxels',   'get_var:voxel_level_index '//trim(filename))

call nc_check(nf90_inq_varid(ncid, 'forecast_lead', varid=VarID), &
        'fill_voxels', 'inq_varid:forecast_lead '//trim(filename))
call nc_check(nf90_get_var(ncid, VarID, forecast_leads), &
        'fill_voxels',   'get_var:forecast_lead '//trim(filename))

call nc_check(nf90_inq_varid(ncid, 'verification_times', varid=VarID), &
        'fill_voxels', 'inq_varid:verification_times '//trim(filename))
call nc_check(nf90_get_var(ncid, VarID, ExperimentTimes), &
        'fill_voxels',   'get_var:verification_times '//trim(filename))

call nc_check(nf90_inq_varid(ncid, 'first_time', varid=VarID), &
        'fill_voxels', 'inq_varid:first_time '//trim(filename))
call nc_check(nf90_get_var(ncid, VarID, first_time), &
        'fill_voxels',   'get_var:first_time '//trim(filename))

call nc_check(nf90_inq_varid(ncid, 'last_time', varid=VarID), &
        'fill_voxels', 'inq_varid:last_time '//trim(filename))
call nc_check(nf90_get_var(ncid, VarID, last_time), &
        'fill_voxels',   'get_var:last_time '//trim(filename))

call nc_check(nf90_inq_varid(ncid, 'ReportTime', varid=VarID), &
        'fill_voxels', 'inq_varid:ReportTime '//trim(filename))
call nc_check(nf90_get_var(ncid, VarID, ReportTimes), &
        'fill_voxels',   'get_var:ReportTime '//trim(filename))

call nc_check(nf90_close(ncid), &
        'fill_voxels', 'close '//trim(filename))

!-----------------------------------------------------------------------
! Convert all the Experiment Times to DART time types.
! While we're at it, let's reshape them to our liking:

if (debug) then
   write(*,*)'nanalyses/num_forecasts are :',nanalyses,num_forecasts
   write(*,*)'shape of VerifyTimes     is :',shape(VerifyTimes)
   write(*,*)'shape of ReportTimes     is :',shape(ReportTimes)
   write(*,*)'                   expected :',num_verif_times,nvoxels
   write(*,*)'ReportTimes for     voxel 1 :',ReportTimes(:,1)
endif

do j = 1,nanalyses
do i = 1,num_forecasts

   days    = floor(ExperimentTimes(i,j))
   seconds = nint((ExperimentTimes(i,j) - days) * 86400.0_digits12)
   VerifyTimes(j,i) = set_time(seconds, days)

   ! FIXME ... could probably write a summary table to stdout/log

enddo
enddo

!-----------------------------------------------------------------------
! The 'voxel_flag' array contains either 0 (unwanted) or 1 (desired)
! from the temporal coverage standpoint.
! Combine that with namelist input to generate subset of voxels we want.

do ivoxel = 1,nvoxels

   allocate(voxel(ivoxel)%times(num_verif_times))

   voxel(ivoxel)%times    = set_time(0,0)        ! SET ONE TIME ONLY
   voxel(ivoxel)%obs_type = get_index_for_type_of_obs(obs_types(ivoxel))

   ! flag the voxel as uninteresting if not desired in namelist
   if (voxel(ivoxel)%obs_type /= obtype_integer) voxel_flag(ivoxel) = 0
   voxel(ivoxel)%good     = voxel_flag(ivoxel)

   myloc = set_location(locations(1,ivoxel), &
                        locations(2,ivoxel), &
                        locations(3,ivoxel), which_vert(ivoxel) )
   voxel(ivoxel)%location    = myloc
   voxel(ivoxel)%LonLatLev   = locations(1:3,ivoxel)
   voxel(ivoxel)%levelindex  = voxel_level_index(ivoxel)
   voxel(ivoxel)%station_id  = 0

   days    = floor(first_time(ivoxel))
   seconds = nint((first_time(ivoxel) - days) * 86400.0_digits12)
   voxel(ivoxel)%first_time  = set_time(seconds, days)

   days    = floor(last_time(ivoxel))
   seconds = nint((last_time(ivoxel) - days) * 86400.0_digits12)
   voxel(ivoxel)%last_time   = set_time(seconds, days)

   do j = 1,num_verif_times
      days    = floor(ReportTimes(j,ivoxel))
      seconds = nint((ReportTimes(j,ivoxel) - days) * 86400.0_digits12)
      voxel(ivoxel)%times(j) = set_time(seconds, days)
   enddo

   if (verbose .and. (ivoxel == 1)) then
      write(*,*)
      write(*,*)'voxel 1 ReportTimes is ',ReportTimes(:,1)
      call print_time(voxel(1)%times(1),'fill_voxels ... first voxel,first time: ')
      call print_date(voxel(1)%times(1),'fill_voxels ... first voxel,first date: ')
      write(*,*)
   endif

enddo

! check to make sure nvoxels is not zero
if (sum(voxel_flag) < 1) then
   write(string1,*)'no valid voxels of ',trim(obtype_string)
   call error_handler(E_ERR,'fill_voxels:',string1,source)
endif

!-----------------------------------------------------------------------
! Summary of voxel information is printed in relate_voxels_to_stations()
!-----------------------------------------------------------------------

deallocate(voxel_flag, obs_types, which_vert, &
           first_time, last_time, voxel_level_index, &
           locations, ReportTimes, ExperimentTimes)

end function fill_voxels



!============================================================================



function relate_voxels_to_stations() result(nstations)

! Determine the number of unique horizontal locations ... AKA 'stations'
! Then, when we know the actual number, create and fill a station array.
!
! Each observation type has its own set of voxels and stations
! voxels   have a obs_type/lon/lat/vert
! stations have a obs_type/lon/lat      but no vert

integer :: nstations

! Local Variables

integer :: ivoxel, itime
integer :: iknown
integer :: istation
logical :: matched

real(r8), allocatable, dimension(:,:) :: unique_locations

if (verbose) then
   call timestamp('relate_voxels_to_stations: starting at ', pos='brief')
endif

! Determine a list of unique horizontal locations.

! temporary array of candidate stations
! worst case scenario - everyone unique.
! unique_locations(:,1) = longitudes
! unique_locations(:,2) = latitudes
! unique_locations(:,3) = observation type

allocate( unique_locations(num_voxels,3) )

! Find the first good voxel and save the location.
! So the list of known locations has something in it for
! future comparisons.

nstations = 0
FirstLocation : do ivoxel = 1,num_voxels
   if (voxel(ivoxel)%good == 1) then
      nstations = 1
      unique_locations(nstations,1) = voxel(1)%LonLatLev(1)
      unique_locations(nstations,2) = voxel(1)%LonLatLev(2)
      unique_locations(nstations,3) = voxel(1)%obs_type
      exit FirstLocation
   endif
enddo FirstLocation

if (nstations == 0) then
   write(string1,*)'Unable to find even 1 good voxel.'
   call error_handler(E_ERR, 'relate_voxels_to_stations:', string1, source)
endif

! Compare each good voxel to the list of known stations.
! If new, add to the list of known stations.

Horizontal : do ivoxel = 1,num_voxels

   if (voxel(ivoxel)%good /= 1) cycle Horizontal

   matched = .false.

   Known : do iknown = 1,nstations
      if ((unique_locations(iknown,1) == voxel(ivoxel)%LonLatLev(1)) .and. &
          (unique_locations(iknown,2) == voxel(ivoxel)%LonLatLev(2)) .and. &
          (unique_locations(iknown,3) == voxel(ivoxel)%obs_type    )) then
          matched = .true.
          voxel(ivoxel)%station_id = iknown
          exit Known
      endif
   enddo Known

   if (.not. matched) then
      nstations = nstations + 1
      unique_locations(nstations,1) = voxel(ivoxel)%LonLatLev(1)
      unique_locations(nstations,2) = voxel(ivoxel)%LonLatLev(2)
      unique_locations(nstations,3) = voxel(ivoxel)%obs_type
      voxel(ivoxel)%station_id = nstations
   endif

enddo Horizontal

! So now we know what the unique horizontal locations are, and how
! many of them there are. Save these to facilitate stowing them in
! the netCDF output file. 

allocate( station(nstations) )

do istation = 1,nstations
   station(istation)%obs_type = unique_locations(istation,3)
   station(istation)%location = set_location( unique_locations(istation,1), &
                                unique_locations(istation,2), 0.0_r8, VERTISUNDEF)
enddo

if (debug) then
   write(string1,*)'There are ',nstations,' unique horizontal locations - AKA "stations".'
   write(string2,*)'and   ', num_forecasts, ' forecast times of interest.'
   call error_handler(E_MSG,'relate_voxels_to_stations:',string1,text2=string2)
   do ivoxel = 1,num_voxels
      write(*,*)'-----------------'
      write(*,'(''Summary for voxel '',i8,'' good '',i1,'' obs_type '',i3, &
         &'' station_id '',i8,'' level index '',i2)') ivoxel, &
         voxel(ivoxel)%good,       voxel(ivoxel)%obs_type, &
         voxel(ivoxel)%station_id, voxel(ivoxel)%levelindex
      call print_date(voxel(ivoxel)%first_time,'first time')
      call print_date(voxel(ivoxel)%last_time, 'last  time')

      do itime = 1,num_verif_times
         write(string1,*)'verify time ',itime,' is '
         write(string2,*)'verify date ',itime,' is '
         call print_time(voxel(ivoxel)%times(itime),trim(string1))
         call print_date(voxel(ivoxel)%times(itime),trim(string2))
      enddo

      write(*,*)'LonLatLev      ',voxel(ivoxel)%LonLatLev
      call write_location(6,voxel(ivoxel)%location,fform='ascii',charstring=string1)
      write(*,*) trim(string1)
   enddo
   write(*,*)'-----------------'
endif

deallocate( unique_locations )

end function relate_voxels_to_stations



!============================================================================



function InitNetCDF(fname)
!----------------------------------------------------------------------------
! Initialize the output netcdf file and fill a local structure with all
! the coordinate variable IDs etc. so we are not incessantly querying them.

character(len=*), intent(in) :: fname
integer                      :: InitNetCDF

integer :: ncid, i, nlines, linelen, ndims
integer ::  LineLenDimID,   nlinesDimID,   stringDimID
integer :: VarID
character(len=nf90_max_name) :: dimName
integer, dimension(nf90_max_var_dims) :: dimIDs, dimLengths

character(len=8)      :: crdate      ! needed by F90 DATE_AND_TIME intrinsic
character(len=10)     :: crtime      ! needed by F90 DATE_AND_TIME intrinsic
character(len=5)      :: crzone      ! needed by F90 DATE_AND_TIME intrinsic
integer, dimension(8) :: values      ! needed by F90 DATE_AND_TIME intrinsic

character(len=256), allocatable, dimension(:) :: textblock

if (verbose) then
   call timestamp('Initializing output file '//trim(fname), pos='brief')
endif

if(.not. byteSizesOK()) then
    call error_handler(E_ERR,'InitNetCDF:', &
   'Compiler does not support required kinds of variables.',source)
endif

InitNetCDF = 0

call nc_check(nf90_create(path = trim(fname), cmode = nf90_clobber, &
         ncid = ncid), 'obs_seq_verify:InitNetCDF', 'create '//trim(fname))

if (debug) then
   write(string1,*)trim(ncName), ' is fortran unit ',ncid
   call error_handler(E_MSG,'InitNetCDF:',string1)
endif


!----------------------------------------------------------------------------
! Write Global Attributes (mostly namelist input, that sort of thing)
!----------------------------------------------------------------------------

call DATE_AND_TIME(crdate,crtime,crzone,values)
write(string1,'(''YYYY MM DD HH MM SS = '',i4,5(1x,i2.2))') &
               values(1), values(2), values(3), values(5), values(6), values(7)
call nc_check(nf90_put_att(ncid, NF90_GLOBAL, 'creation_date', trim(string1) ), &
           'InitNetCDF', 'put_att creation_date '//trim(fname))

call nc_check(nf90_put_att(ncid, NF90_GLOBAL, 'source',   source), &
           'InitNetCDF', 'put_att   source '//trim(fname))

! write all observation sequence files used
FILEloop : do i = 1,num_input_files

  write(string1,'(''obs_seq_file_'',i3.3)')i
  call nc_check(nf90_put_att(ncid, NF90_GLOBAL, &
             trim(string1), trim(obs_sequences(i)) ), &
             'InitNetCDF', 'put_att:filenames')

enddo FILEloop

!----------------------------------------------------------------------------
! Define the dimensions
! Set nofill mode - supposed to be performance gain
!----------------------------------------------------------------------------

ncmeta%AnalysisDimLen = NF90_UNLIMITED
ncmeta%CopyDimLen     = NCOPIES
ncmeta%StationsDimLen = num_stations
ncmeta%LevelsDimLen   = num_levels
ncmeta%EnsembleDimLen = ensemble_size
ncmeta%ForecastDimLen = size(forecast_leads)

call nc_check(nf90_set_fill(ncid, NF90_NOFILL, i),  &
              'obs_seq_verify:InitNetCDF', 'set_nofill '//trim(fname))

call nc_check(nf90_def_dim(ncid=ncid, name=ncmeta%AnalysisDimName, &
               len=ncmeta%AnalysisDimLen, dimid = ncmeta%AnalysisDimID), &
               'InitNetCDF', 'def_dim:analysisT '//trim(fname))

call nc_check(nf90_def_dim(ncid=ncid, name=ncmeta%CopyDimName, &
               len=ncmeta%CopyDimLen, dimid = ncmeta%CopyDimID), &
               'InitNetCDF', 'def_dim:copy '//trim(fname))

call nc_check(nf90_def_dim(ncid=ncid, name=ncmeta%StationsDimName, &
               len=ncmeta%StationsDimLen, dimid = ncmeta%StationsDimID), &
               'InitNetCDF', 'def_dim:station '//trim(fname))

call nc_check(nf90_def_dim(ncid=ncid, name=ncmeta%LevelsDimName, &
               len=ncmeta%LevelsDimLen, dimid = ncmeta%LevelsDimID), &
               'InitNetCDF', 'def_dim:level '//trim(fname))

call nc_check(nf90_def_dim(ncid=ncid, name=ncmeta%EnsembleDimName, &
               len=ncmeta%EnsembleDimLen, dimid = ncmeta%EnsembleDimID), &
               'InitNetCDF', 'def_dim:ensemble '//trim(fname))

call nc_check(nf90_def_dim(ncid=ncid, name=ncmeta%ForecastDimName, &
               len=ncmeta%ForecastDimLen, dimid = ncmeta%ForecastDimID), &
               'InitNetCDF', 'def_dim:forecast_lead '//trim(fname))

! namelist quantities

call find_textfile_dims('input.nml', nlines, linelen)
allocate(textblock(nlines))
textblock = ''

call nc_check(nf90_def_dim(ncid=ncid, name='linelen', len=len(textblock(1)), &
               dimid =  linelenDimID), 'InitNetCDF', 'def_dim:linelen '//'input.nml')

call nc_check(nf90_def_dim(ncid=ncid, name='nlines', len=nlines, &
               dimid =   nlinesDimID), 'InitNetCDF', 'def_dim:nlines '//'input.nml')

call nc_check(nf90_def_dim(ncid=ncid, name='stringlength', len=metadatalength, &
               dimid =  StringDimID),  'InitNetCDF', 'def_dim:stringlength '//trim(fname))

! Define the variable to record the input parameters ... the namelist

call nc_check(nf90_def_var(ncid=ncid, name='namelist', xtype=nf90_char, &
          dimids = (/ linelenDimID, nlinesDimID /), varid=VarID), &
          'InitNetCDF', 'namelist:def_var')
call nc_check(nf90_put_att(ncid, VarID, 'long_name', 'input.nml contents'), &
          'InitNetCDF', 'namelist:long_name')

! Define the variable for interpreting the 'copy' dimension

call nc_check(nf90_def_var(ncid=ncid, name='CopyMetaData', xtype=nf90_char, &
          dimids = (/ StringDimID, ncmeta%CopyDimID /), varid=VarID), &
          'InitNetCDF', 'copymeta:def_var')
call nc_check(nf90_put_att(ncid, VarID, 'long_name', 'copy quantity names'), &
          'InitNetCDF', 'copymeta:long_name')

!----------------------------------------------------------------------------
! Define the coordinate variables
!----------------------------------------------------------------------------

call nc_check(nf90_def_var(ncid=ncid, name=ncmeta%AnalysisDimName, xtype=nf90_double, &
          dimids=(/ ncmeta%AnalysisDimID /), varid=VarID), &
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

call nc_check(nf90_def_var(ncid=ncid, name=ncmeta%CopyDimName, xtype=nf90_int, &
          dimids = (/ ncmeta%CopyDimID /), varid=VarID), &
          'InitNetCDF', 'copy:def_var')
call nc_check(nf90_put_att(ncid, VarID, 'long_name', 'observation copy'), &
          'InitNetCDF', 'copy:long_name')

do i = 1,NCOPIES
   write(string1,'(''note'',i1)')i
   write(string2,'(i1,'' == '',a)')i,trim(copy_metadata(i))
   call nc_check(nf90_put_att(ncid, VarID, trim(string1),trim(string2)), &
          'InitNetCDF', 'copy'//trim(string1))
enddo

call nc_check(nf90_put_att(ncid, VarID, 'explanation', 'see CopyMetaData variable'), &
          'InitNetCDF', 'copy:explanation')

call nc_check(nf90_def_var(ncid=ncid, name=ncmeta%StationsDimName, xtype=nf90_int, &
          dimids = (/ ncmeta%StationsDimID /), varid=VarID), &
          'InitNetCDF', 'station:def_var')
call nc_check(nf90_put_att(ncid, VarID, 'long_name', 'station index'), &
          'InitNetCDF', 'station:long_name')

call nc_check(nf90_def_var(ncid=ncid, name=ncmeta%LevelsDimName, xtype=nf90_double, &
          dimids = (/ ncmeta%LevelsDimID /), varid=VarID), &
          'InitNetCDF', 'level:def_var')
call nc_check(nf90_put_att(ncid, VarID, 'long_name', 'vertical level of observation'), &
          'InitNetCDF', 'level:long_name')

call nc_check(nf90_def_var(ncid=ncid, name=ncmeta%EnsembleDimName, xtype=nf90_int, &
          dimids = (/ ncmeta%EnsembleDimID /), varid=VarID), &
          'InitNetCDF', 'ensemble:def_var')
call nc_check(nf90_put_att(ncid, VarID, 'long_name', 'ensemble member'), &
          'InitNetCDF', 'ensemble:long_name')

call nc_check(nf90_def_var(ncid=ncid, name=ncmeta%ForecastDimName, xtype=nf90_int, &
          dimids = (/ ncmeta%ForecastDimID /), varid=VarID), &
          'InitNetCDF', 'forecast_lead:def_var')
call nc_check(nf90_put_att(ncid, VarID, 'long_name', 'forecast lead time'), &
          'InitNetCDF', 'forecast_lead:long_name')
call nc_check(nf90_put_att(ncid, VarID, 'units',     'seconds'), &
          'InitNetCDF', 'forecast_lead:units')

! let the location module write what it needs to ...

call nc_write_location_atts(ncid, num_stations, ncmeta%StationsDimID)

!----------------------------------------------------------------------------
! Define the RECORD variables
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
        dimids = (/ ncmeta%ForecastDimID, ncmeta%EnsembleDimID,     ncmeta%CopyDimID, &
                      ncmeta%LevelsDimID, ncmeta%StationsDimID, ncmeta%AnalysisDimID /), &
        varid = VarID), 'InitNetCDF', 'obtype_string:def_var')

call nc_check(nf90_put_att(ncid, VarID, 'long_name', 'forecast variable quantities'), &
            'InitNetCDF', 'obtype_string:long_name')
call nc_check(nf90_put_att(ncid, VarID, 'missing_value', MISSING_R8), &
            'InitNetCDF', 'obtype_string:missing_value')
call nc_check(nf90_put_att(ncid, VarID, '_FillValue',    MISSING_R8), &
             'InitNetCDF', 'obtype_string:fill_value')

call nc_check(nf90_def_var(ncid=ncid, name='original_qc', xtype=nf90_int, &
          dimids = (/ ncmeta%ForecastDimID, ncmeta%LevelsDimID, ncmeta%StationsDimID, ncmeta%AnalysisDimID /), varid=VarID), &
          'InitNetCDF', 'original_qc:def_var')
call nc_check(nf90_put_att(ncid, VarID, 'long_name', 'original QC value'), &
          'InitNetCDF', 'original_qc:long_name')
call nc_check(nf90_put_att(ncid, VarID, 'missing_value', MISSING_I), &
          'InitNetCDF', 'original_qc:missing')
call nc_check(nf90_put_att(ncid, VarID, '_FillValue',    MISSING_I), &
          'InitNetCDF', 'original_qc:fill_value')

call nc_check(nf90_def_var(ncid=ncid, name='dart_qc', xtype=nf90_int, &
          dimids = (/ ncmeta%ForecastDimID, ncmeta%LevelsDimID, ncmeta%StationsDimID, ncmeta%AnalysisDimID /), varid=VarID), &
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

call nc_check(nf90_inq_varid(ncid, ncmeta%CopyDimName, varid=VarID), &
          'InitNetCDF', 'inq_varid:copy '//trim(fname))
call nc_check(nf90_put_var(ncid, VarID, (/ (i,i=1,NCOPIES) /) ), &
          'InitNetCDF', 'put_var:copy')

call nc_check(nf90_inq_varid(ncid, ncmeta%EnsembleDimName, varid=VarID), &
          'InitNetCDF', 'inq_varid:ensemble '//trim(fname))
call nc_check(nf90_put_var(ncid, VarID, (/ (i,i=1,ensemble_size) /) ), &
          'InitNetCDF', 'put_var:ensemble')

call nc_check(nf90_inq_varid(ncid, ncmeta%StationsDimName, varid=VarID), &
          'InitNetCDF', 'inq_varid:station '//trim(fname))
call nc_check(nf90_put_var(ncid, VarID, (/ (i,i=1,num_stations) /) ), &
          'InitNetCDF', 'put_var:station')

call nc_check(nf90_inq_varid(ncid, ncmeta%LevelsDimName, varid=VarID), &
          'InitNetCDF', 'inq_varid:level '//trim(fname))
call nc_check(nf90_put_var(ncid, VarID, mandatory_level ), &
          'InitNetCDF', 'put_var:level')

call nc_check(nf90_inq_varid(ncid, ncmeta%ForecastDimName, varid=VarID), &
          'InitNetCDF', 'inq_varid:forecast_lead '//trim(fname))
call nc_check(nf90_put_var(ncid, VarID, forecast_leads ), &
          'InitNetCDF', 'put_var:forecast_lead')

! Record station locations using location_mod:nc_write_location()

do i = 1,num_stations
   call nc_write_location(ncid, station(i)%location, i, do_vert=.true.)
enddo

!----------------------------------------------------------------------------
! Finish up ...
!----------------------------------------------------------------------------

call nc_check(nf90_sync( ncid), 'InitNetCDF', 'sync '//trim(fname))

InitNetCDF = ncid

if (verbose) then
   write(*,*)   ! a little whitespace

   call nc_check(nf90_inq_varid(ncid, trim(obtype_string), varid=VarID), &
          'InitNetCDF', 'inq_varid:summary '//trim(fname))
   call nc_check(nf90_inquire_variable(ncid, VarID, ndims=ndims, dimids=dimIDs), &
          'InitNetCDF', 'inquire obtype_string dimids '//trim(fname))

   do i = 1,ndims
      write(string1,*)'dimlen ',i,' of ',trim(obtype_string),trim(fname)
      call nc_check(nf90_inquire_dimension(ncid, dimIDs(i), name=dimName, &
               len=dimLengths(i)), 'InitNetCDF', string1)
      write(string2,*)trim(obtype_string), ' dimID ', i, ' length/name ', &
                      dimLengths(i),trim(dimName)
      call error_handler(E_MSG,'InitNetCDF:',string2)
   enddo
endif

end function InitNetCDF



!============================================================================



subroutine WriteNetCDF(ncid, fname, values, obs_index, &
              prior_indices, ob_err, org_qc, dart_qc, ncstart, nccount)
!----------------------------------------------------------------------------
! At this point, we have the (ncopies,ensemble_size) hyperslab of data, but
! it is not in the right shape to simply blast into the netCDF file.
 
integer,                     intent(in) :: ncid
character(len=*),            intent(in) :: fname
real(r8), dimension(:),      intent(in) :: values
integer,                     intent(in) :: obs_index
integer, dimension(:),       intent(in) :: prior_indices
real(r8),                    intent(in) :: ob_err
real(r8),                    intent(in) :: org_qc
real(r8),                    intent(in) :: dart_qc
integer, dimension(:),       intent(in) :: ncstart
integer, dimension(:),       intent(in) :: nccount

integer, dimension(1) :: istart, icount

integer :: obscopyindex, priorcopyindex, errorcopyindex
integer :: i, seconds, days

integer, dimension(nf90_max_var_dims) :: dimIDs

integer ::  TimeVarID
integer ::  VarID
integer ::  ndims
integer ::  QCVarID
integer ::  DARTQCVarID

real(digits12) :: fdays

real(r8), dimension(NCOPIES,ensemble_size) :: datmat

!----------------------------------------------------------------------------
! Use the dimension ID's and lengths from the 'ncmeta' variable
!----------------------------------------------------------------------------

call nc_check(nf90_inq_varid(ncid, 'analysisT', varid=TimeVarID), &
          'WriteNetCDF', 'inq_varid:analysisT '//trim(fname))

call nc_check(nf90_inq_varid(ncid, 'original_qc', varid=QCVarID), &
          'WriteNetCDF', 'inq_varid:original_qc '//trim(fname))

call nc_check(nf90_inq_varid(ncid, 'dart_qc', varid=DARTQCVarID), &
          'WriteNetCDF', 'inq_varid:dart_qc '//trim(fname))

! Increase the record dimension ... once for each input/analysis file

istart(1) = AnalysisIndex
icount(1) = 1

call get_time(analysisT,seconds,days)
fdays = days + seconds/86400.0_digits12

call nc_check(nf90_put_var(ncid, TimeVarId, (/ fdays /), &
            start=istart, count=icount ), 'WriteNetCDF', 'put_var:analysisT')

call nc_check(nf90_inq_varid(ncid, trim(obtype_string), varid=VarID), &
          'WriteNetCDF', 'inq_varid:obtype_string '//trim(fname))

! FIXME ... Get shape of variable to ensure conformable ... 
! no actual checking going on ... onerous if called a million times
! call nc_check(nf90_inquire_variable(ncid, VarID, ndims=ndims, dimids=dimIDs), &
!           'WriteNetCDF', 'inquire obtype_string dimids '//trim(fname))

datmat = MISSING_R8

  obscopyindex = 1
priorcopyindex = 2
errorcopyindex = 3

do i = 1,ensemble_size

      datmat(  obscopyindex,i) = values(obs_index)
      datmat(priorcopyindex,i) = values(prior_indices(i))
      datmat(errorcopyindex,i) = ob_err

enddo

call nc_check(nf90_put_var(ncid, VarId, datmat, start=ncstart, count=nccount), &
        'WriteNetCDF', 'put_var:datmat')

! original_qc(analysisT, station, level,                 forecast_lead) ;
!     dart_qc(analysisT, station, level,                 forecast_lead) ;
!     trouble(analysisT, station, level, copy, ensemble, forecast_lead) ;
!              dim6       dim5     dim4  dim3    dim2        dim1

! possible FIXME ... shape check ... implicit assumption about shape
call nc_check(nf90_inquire_variable(ncid, QCVarID, ndims=ndims, dimids=dimIDs), &
          'WriteNetCDF', 'inquire original_qc dimids '//trim(fname))
call nc_check(nf90_put_var(ncid, QCVarId, (/ org_qc /), &
          start=(/ ncstart(1), ncstart(4), ncstart(5), ncstart(6) /), &
          count=(/ nccount(1), nccount(4), nccount(5), nccount(6) /)), &
          'WriteNetCDF', 'put_var:orgqcmat')

! possible FIXME ... shape check ... implicit assumption about shape
call nc_check(nf90_inquire_variable(ncid, DARTQCVarID, ndims=ndims, dimids=dimIDs), &
          'WriteNetCDF', 'inquire dart_qc dimids '//trim(fname))
call nc_check(nf90_put_var(ncid, DARTQCVarId, (/ dart_qc /), &
          start=(/ ncstart(1), ncstart(4), ncstart(5), ncstart(6) /), &
          count=(/ nccount(1), nccount(4), nccount(5), nccount(6) /)), &
          'WriteNetCDF', 'put_var:orgqcmat')

!----------------------------------------------------------------------------
! finished ...
!----------------------------------------------------------------------------

call nc_check(nf90_sync( ncid), 'WriteNetCDF', 'sync '//trim(fname))

end subroutine WriteNetCDF



!============================================================================



subroutine CloseNetCDF(ncid, fname)
!----------------------------------------------------------------------------
integer,          intent(in) :: ncid
character(len=*), intent(in) :: fname

integer :: i, VarID, ndims
integer, dimension(nf90_max_var_dims) :: dimIDs, dimLengths
character(len=nf90_max_name) :: dimName

if (verbose) then
   write(*,*)   ! a little whitespace

   call nc_check(nf90_inq_varid(ncid, trim(obtype_string), varid=VarID), &
          'CloseNetCDF', 'inq_varid:summary '//trim(fname))
   call nc_check(nf90_inquire_variable(ncid, VarID, ndims=ndims, dimids=dimIDs), &
          'CloseNetCDF', 'inquire obtype_string dimids '//trim(fname))

   do i = 1,ndims
      write(string1,*)'dimlen ',i,' of ',trim(obtype_string),trim(fname)
      call nc_check(nf90_inquire_dimension(ncid, dimIDs(i), name=dimName, &
               len=dimLengths(i)), 'CloseNetCDF', string1)
      write(string2,*)trim(obtype_string), ' dimID ', i, ' length/name ', &
                      dimLengths(i),trim(dimName)
      call error_handler(E_MSG,'CloseNetCDF:',string2)
   enddo
endif

call nc_check(nf90_close(ncid), 'CloseNetCDF:', 'close '//trim(fname))

end subroutine CloseNetCDF



!============================================================================



subroutine destroy_voxels()
!----------------------------------------------------------------------------

integer :: i

do i = 1,num_voxels
  if (associated(voxel(i)%times)) then
     deallocate( voxel(i)%times )
     nullify(    voxel(i)%times )
  endif
enddo

if (allocated(voxel)) deallocate(voxel)

end subroutine destroy_voxels



!============================================================================



subroutine find_ensemble_size()
!----------------------------------------------------------------------------
! The only way I know of to determine the ensemble size is to read the
! copy metadata from one of the obs_seq.xxx files. I'm just going to
! open the first file, read the entire damn sequence, and count them up.
!
! Then, I need to preserve all the indices for the observation value
! as well as the indices for the ensemble members. The observation error
! variance comes from get_obs_def_error_variance()

integer :: filenum = 1

if ( file_exist(obs_sequences(filenum)) ) then
   write(string1,*)'opening ', trim(obs_sequences(filenum))
   call error_handler(E_MSG,'find_ensemble_size:',string1)
else
   write(string1,*)trim(obs_sequences(filenum)), &
                   ' does not exist. Dying a dramatic death.'
   call error_handler(E_ERR,'find_ensemble_size:',string1,source)
endif

!-----------------------------------------------------------------------
! Read the entire observation sequence - allocates 'seq' internally
! And then cruise throught the metadata
!-----------------------------------------------------------------------

call read_obs_seq(obs_sequences(filenum), 0, 0, 0, seq)

ensemble_size = 0

do i=1, get_num_copies(seq)
   if( index(get_copy_meta_data(seq,i), 'prior ensemble member') > 0) &
      ensemble_size = ensemble_size + 1
enddo

if (ensemble_size < 1) then
   write(string1,*)'no ensemble member info in ', trim(obs_sequences(filenum))
   write(string2,*)'looking for a prior ensemble mean'
   call error_handler(E_MSG,'find_ensemble_size:', string1, text2=string2)

   ! could be a 'deterministic' forecast from the ensemble mean
   do i=1, get_num_copies(seq)
      if( index(get_copy_meta_data(seq,i), 'prior ensemble mean') > 0) &
         ensemble_size = ensemble_size + 1
   enddo
endif

if (ensemble_size < 1) then
   write(string1,*)'no ensemble member info in ', trim(obs_sequences(filenum))
   write(string2,*)'cannot continue'
   call error_handler(E_ERR,'find_ensemble_size:', string1, source, text2=string2)
else
   call error_handler(E_MSG,'find_ensemble_size:', 'found a prior ensemble mean.')
endif

write(string1,'(''There are '',i4,'' ensemble members.'')') ensemble_size
if (verbose) call error_handler(E_MSG,'find_ensemble_size:',string1)

call destroy_obs_sequence(seq)

end subroutine find_ensemble_size



!============================================================================



subroutine find_analysis_time(filename, ExpIndex, analysis_time)
!----------------------------------------------------------------------------
! Parse the filename into an integer string that defines the
! analysis time.
!

character(len=*), intent(in)  :: filename
integer,          intent(out) :: ExpIndex
type(time_type),  intent(out) :: analysis_time

character(len=LEN(filename)) :: lj_filename

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
   call error_handler(E_ERR,'find_analysis_time:', string1, source, text2=string2)
endif

read(lj_filename(indx+0:indx+3),'(i4)',iostat=mystat)iyear
if (mystat /= 0) then
   write(string1,*)'Cannot find the year in the extension (YYYYMMDDHH) on '//trim(lj_filename)
   call error_handler(E_ERR,'find_analysis_time:', string1, source, text2=string2)
endif

read(lj_filename(indx+4:indx+5),'(i2)',iostat=mystat)imonth
if (mystat /= 0) then
   write(string1,*)'Cannot find the month in the extension (YYYYMMDDHH) on '//trim(lj_filename)
   call error_handler(E_ERR,'find_analysis_time:', string1, source, text2=string2)
endif

read(lj_filename(indx+6:indx+7),'(i2)',iostat=mystat)iday
if (mystat /= 0) then
   write(string1,*)'Cannot find the day in the extension (YYYYMMDDHH) on '//trim(lj_filename)
   call error_handler(E_ERR,'find_analysis_time:', string1, source, text2=string2)
endif

read(lj_filename(indx+8:indx+9),'(i2)',iostat=mystat)ihour
if (mystat /= 0) then
   write(string1,*)'Cannot find the hour in the extension (YYYYMMDDHH) on '//trim(lj_filename)
   call error_handler(E_ERR,'find_analysis_time:', string1, source, text2=string2)
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
   call error_handler(E_ERR,'find_analysis_time:',string1,source)
endif

if (debug) then
   write(*,*)'find_analysis_time iyear  is ',iyear
   write(*,*)'find_analysis_time imonth is ',imonth
   write(*,*)'find_analysis_time iday   is ',iday
   write(*,*)'find_analysis_time ihour  is ',ihour
   write(*,*)'find_analysis_time resulting AnalysisIndex is ',ExpIndex
endif

if (verbose) then
   write(string1,*)'analysis ',AnalysisIndex,' date is '
   call print_date(analysisT,string1)
endif

end subroutine find_analysis_time



!============================================================================



subroutine find_our_copies(myseq, obsindex, indices)
!----------------------------------------------------------------------------
!
type(obs_sequence_type), intent(in)  :: myseq
integer,                 intent(out) :: obsindex
integer, dimension(:),   intent(out) :: indices
! integer,               intent(out) :: qc_index      ! global variable
! integer,               intent(out) :: dart_qc_index ! global variable

integer :: i, myindex

character(len=obstypelength) :: metadata

obsindex      = -1
qc_index      = -1
dart_qc_index = -1
indices       = -1

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
   call error_handler(E_ERR,'find_our_copies:',string1,source)
endif

myindex = 0
! All prior copies get tacked one after another ...
MetaDataLoop : do i=1, get_num_copies(myseq)

   metadata = get_copy_meta_data(myseq,i)

   if( index(metadata, 'prior ensemble member') > 0) then
      myindex = myindex + 1
      if ( myindex > size(indices) ) then
         write(string1,*)'Found too many prior copies in metadata.'
         call error_handler(E_ERR,'find_our_copies:',string1,source)
      endif
      indices(myindex) = i
   endif

   ! Check for deterministic forecast, apparently
   if ( (ensemble_size == 1) .and. (index(metadata, 'prior ensemble mean') > 0) ) then
      myindex = myindex + 1
      if ( myindex > size(indices) ) then
         write(string1,*)'Found too many prior copies in metadata.'
         call error_handler(E_ERR,'find_our_copies:',string1,source)
      endif
      indices(myindex) = i
      exit MetaDataLoop
   endif

enddo MetaDataLoop

if (myindex == 0) then
   write(string1,*)'Could not find any prior copy indices in sequence.'
   call error_handler(E_ERR,'find_our_copies:',string1,source)
endif

! Here is the sanity check to make sure the indices pick off just the
! prior ensemble members.

if (debug) then
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
   if(index(get_qc_meta_data(myseq,i),'Quality Control'     ) > 0)      qc_index = i
   if(index(get_qc_meta_data(myseq,i),'NCEP QC index'       ) > 0)      qc_index = i
   if(index(get_qc_meta_data(myseq,i),'DART quality control') > 0) dart_qc_index = i
enddo QCMetaDataLoop

if (      qc_index < 0 ) then
   write(string1,*)'metadata:Quality Control copyindex not found'
   call error_handler(E_ERR,'find_qc_indices:',string1,source)
endif
if ( dart_qc_index < 0 ) then
   write(string1,*)'metadata:DART quality control copyindex not found'
   call error_handler(E_ERR,'find_qc_indices:',string1,source)
endif

! Just echo what we know

if (debug) then
   write(*,*)'QC index',     qc_index,' ',trim(get_qc_meta_data(myseq,     qc_index))
   write(*,*)'QC index',dart_qc_index,' ',trim(get_qc_meta_data(myseq,dart_qc_index))
   write(*,*)
endif

end subroutine find_our_copies



!============================================================================



subroutine match_obs_to_voxel(iobs, ob_type, ob_loc, ob_time, vid, forecast_lead_index)
! This routine determines if the observation is within a voxel. 
! The list of voxels should be considered a superset, in that there are some
! voxels for observation types not likely to be encountered. (see Note 1)
!
! There are also observations that are not in a voxel. Just skip them.
!
! The voxels are defined by a lat/lon/level/obs_type ... and a time.
! The time is designed to be _close_ to a forecast/verification time, but 
! sometimes the observations are time-tagged slightly before or after the 
! nominal time.
!
! Note 1: You can create a list of voxels for multiple observation types, 
! but obs_seq_verify() only works for one observation type at-a-time.
!
integer,             intent(in)  :: iobs
integer,             intent(in)  :: ob_type
type(location_type), intent(in)  :: ob_loc
type(time_type),     intent(in)  :: ob_time
integer,             intent(out) :: vid
integer,             intent(out) :: forecast_lead_index

! Local variables

integer                :: ilevel, ivoxel, vt1, vtN
integer,  dimension(1) :: iz
real(r8), dimension(3) :: obslocarray
real(r8), dimension(num_levels) :: distances

! Set these values right away because there are some 'early' returns.

vid = -1
forecast_lead_index = -1

! Find the level

obslocarray = get_location(ob_loc)

if ( is_vertical(ob_loc, "PRESSURE") ) then
   ! Determine which mandatory vertical level is closest

   distances = abs(obslocarray(3) - mandatory_level)
   iz        = minloc(distances)
   ilevel    = iz(1)

   if ( abs(obslocarray(3) - mandatory_level(ilevel)) < OnePa ) then
      obslocarray(3) = mandatory_level(ilevel)
   else
      ! not close enough
      return
   endif

elseif ( is_vertical(ob_loc, "UNDEFINED") ) then
   ilevel = 1

elseif ( is_vertical(ob_loc, "SURFACE") ) then
   ilevel = 1

else ! unwanted
   return
endif

! Find the voxel and time index
! Must assume there are multiple voxels with the same horizontal 
! location but different vertical values (eg. for each mandatory pressure level)
! voxels for observations with vertisundef or vertissurface may have garbage
! in their vertical level. Those voxels may also have indeterminate values.
! There are also some voxels that are not sampled in time sufficiently to be of interest.

VoxelLoop : do ivoxel = 1,num_voxels

   if (voxel(ivoxel)%good       /=       1 ) cycle VoxelLoop  ! bad temporal sampling
   if (voxel(ivoxel)%obs_type   /= ob_type ) cycle VoxelLoop
   if (voxel(ivoxel)%levelindex /=  ilevel ) cycle VoxelLoop

   if (abs(obslocarray(1) - voxel(ivoxel)%LonLatLev(1)) > HALF_METER) cycle VoxelLoop
   if (abs(obslocarray(2) - voxel(ivoxel)%LonLatLev(2)) > HALF_METER) cycle VoxelLoop

   vid = ivoxel

   ! By now, we know the horizontal and vertical are desired.
   ! Could still encounter an observation at a time we do not need.
   !
   ! The voxel times are the superset of all times.
   ! should only loop over times for this analysis cycle ...

   vt1 = AnalysisIndex
   vtN = AnalysisIndex + num_forecasts - 1

   TimeLoop : do i = vt1,vtN

      if  ( ob_time /= voxel(ivoxel)%times(i) ) cycle TimeLoop

      ! Found the matching observation time ... 
      ! record the corresponding Forecast lead index
      
      forecast_lead_index = i - vt1 + 1
      exit VoxelLoop

   enddo TimeLoop
enddo VoxelLoop

if (forecast_lead_index < 0) then
   ! Could not find a location with a desired time.
   ! Not an error, just return.
   vid = -1
   return
endif

! Summarize what happened to this observation

if (debug) then

   write(*,*)'-------------------'
   write(string1,'('' desired ob ('',i7,'') voxel='',i6, &
      &'' analysis index= '',i3, '' lead_index= '',i3,'' ob date='')') &
       iobs, vid, AnalysisIndex, forecast_lead_index

   call print_time(ob_time,trim(string1))
   call write_location(6,ob_loc,             fform='ascii',charstring=string1)
   call write_location(6,voxel(vid)%location,fform='ascii',charstring=string2)
   call write_location(6,station(voxel(vid)%station_id)%location, &
                                             fform='ascii',charstring=string3)
   write(*,*)'observation ',trim(string1)
   write(*,*)'voxel       ',trim(string2)
   write(*,*)'station     ',trim(string3)
   write(*,*)'station ID  ',voxel(vid)%station_id

endif

end subroutine match_obs_to_voxel



!============================================================================



subroutine determine_hyperslab_indices(analysisTindex, forecast_lead_index, &
              vid, ncstart, nccount)
! We have harvested all the information required to be able to stuff the observation
! into the proper slot in the output structure, i.e.:
! RADIOSONDE_TEMPERATURE(analysisT, station, level, copy, ensemble, forecast_lead) ;
! The job now is to construct the appropriate hyperslabbing given what we know.

integer,               intent(in)  :: analysisTindex
integer,               intent(in)  :: forecast_lead_index
integer,               intent(in)  :: vid
integer, dimension(:), intent(out) :: ncstart
integer, dimension(:), intent(out) :: nccount

ncstart(:) = 0
nccount(:) = 0

ncstart(6) = analysisTindex
ncstart(5) = voxel(vid)%station_id
ncstart(4) = voxel(vid)%levelindex
ncstart(3) = 1
ncstart(2) = 1
ncstart(1) = forecast_lead_index

nccount(6) = 1
nccount(5) = 1
nccount(4) = 1
nccount(3) = NCOPIES   ! observation, forecast, obs error variance
nccount(2) = ensemble_size
nccount(1) = 1

if ( debug ) then
   write(   *   ,'('' determine_hyperslab_indices           '' ,6(1x,a11))') &
               'lead_index','ensemble','copy','level','station','analysis'
   write(string1,*)'"ncstart"',ncstart
   write(string2,*)'"nccount"',nccount
   call error_handler(E_MSG,'determine_hyperslab_indices:',string1)
   call error_handler(E_MSG,'determine_hyperslab_indices:',string2)
endif

end subroutine determine_hyperslab_indices



!============================================================================



end program obs_seq_verify

