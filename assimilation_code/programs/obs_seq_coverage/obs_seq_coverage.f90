! DART software - Copyright UCAR. This open source software is provided
! by UCAR, "as is", without charge, subject to all terms of use at
! http://www.image.ucar.edu/DAReS/DART/DART_download

!> This program queries a bunch of obs_seq.xxxx files and tries to
!> figure out 'voxel coverage' ... what locations are consistently
!> reported through time. Absolutely a 'reverse-engineering exercise'.
!>
!> The observation sequence file only contains lat/lon/level/which_vert,
!> so this is all we have to work with.

program obs_seq_coverage

use        types_mod, only : r4, r8, digits12, MISSING_R8, MISSING_R4, MISSING_I, &
                             metadatalength, obstypelength

use obs_sequence_mod, only : read_obs_seq, obs_type, obs_sequence_type, get_first_obs, &
                             get_obs_def, get_copy_meta_data, &
                             get_next_obs, init_obs, init_obs_sequence, &
                             assignment(=), get_num_copies, get_num_qc, get_qc, &
                             static_init_obs_sequence, destroy_obs_sequence, destroy_obs, &
                             read_obs_seq_header, get_qc_meta_data

use      obs_def_mod, only : obs_def_type, get_obs_def_time, get_obs_def_type_of_obs, write_obs_def, &
                             get_obs_def_location, set_obs_def_time, &
                             set_obs_def_location, set_obs_def_type_of_obs, set_obs_def_error_variance

use     obs_kind_mod, only : max_defined_types_of_obs, get_name_for_type_of_obs, get_index_for_type_of_obs, &
                             write_type_of_obs_table

use     location_mod, only : location_type, get_location, set_location_missing, &
                             write_location, operator(/=), operator(==), &
                             set_location, is_location_in_region, query_location, &
                             is_vertical, VERTISUNDEF

use  location_io_mod, only : nc_write_location_atts, nc_write_location

use time_manager_mod, only : time_type, set_date, set_time, get_time, &
                             set_calendar_type, get_calendar_string, &
                             print_time, print_date, &
                             operator(+), operator(-), operator(<), operator(>), &
                             operator(==), operator(/=), operator(<=), operator(/), &
                             operator(>=), operator(*)

use    utilities_mod, only : get_unit, close_file, &
                             file_exist, error_handler, E_ERR, E_WARN, E_MSG, &
                             initialize_utilities, nmlfileunit, finalize_utilities, &
                             find_namelist_in_file, check_namelist_read, &
                             next_file, set_filename_list, find_textfile_dims, &
                             file_to_text, do_nml_file, do_nml_term

use  netcdf_utilities_mod, only : nc_check

use typeSizes
use netcdf

implicit none

character(len=*), parameter :: source = 'obs_seq_coverage.f90'

!---------------------------------------------------------------------
!---------------------------------------------------------------------

type voxel_type
   integer                  :: obs_type
   type(location_type)      :: location
   integer                  :: levelindx
   type(time_type)          :: first_time
   type(time_type)          :: last_time
   integer                  :: ntimes
   type(time_type), pointer :: times(:)
end type voxel_type

logical,          allocatable, dimension(:) :: Desiredvoxels
type(voxel_type), allocatable, dimension(:) :: voxels

integer :: num_voxels    ! This is the current number of unique locations
integer :: max_voxels    ! This is the largest possible number of uniq locs
integer :: voxel_id      ! the index (into voxels) of an existing location
integer :: timeindex     ! the index (into the time array of a voxel)
integer :: num_out_stat  ! total number of desired voxels found
integer :: num_out_total ! total number of desired locations * times found

integer,  parameter :: MAX_OBS_NAMES_IN_NAMELIST = 500  ! lazy, just going big
integer,  parameter :: NUM_MANDATORY_LEVELS  = 14  ! number of mandatory pressure levels
real(r8), parameter :: OnePa = 1.0_r8   ! 1 Pa ... tolerance for vertical comparisons

!---------------------------------------------------------------------
! variables associated with the observation
!---------------------------------------------------------------------

type(obs_sequence_type) :: seq
type(obs_type)          :: obs1, obs2
type(obs_def_type)      :: obs_def
type(location_type)     :: obs_loc, minl, maxl

integer :: flavor, flavor_of_interest
integer :: num_copies, num_qc, num_obs, max_num_obs, obs_seq_file_id

character(len=129) :: obs_seq_read_format
logical :: pre_I_format
logical :: last_ob_flag

character(len=256) :: obs_seq_in_file_name
integer, parameter :: MAX_INPUT_FILES = 500
integer            :: num_input_files

!-----------------------------------------------------------------------
! Namelist with (some scalar) default values
!-----------------------------------------------------------------------

character(len=256) :: obs_sequences(MAX_INPUT_FILES) = ''
character(len=256) :: obs_sequence_list = ''
character(len=256) :: textfile_out      = 'obsdef_mask.txt'
character(len=256) :: netcdf_out        = 'obsdef_mask.nc'
character(len=129) :: calendar          = 'Gregorian'
character(len=obstypelength) :: obs_of_interest(MAX_OBS_NAMES_IN_NAMELIST) = ''

integer, dimension(6) :: first_analysis = (/ 2003, 1, 1, 0, 0, 0 /)
integer, dimension(6) ::  last_analysis = (/ 2003, 1, 2, 0, 0, 0 /)
integer  :: forecast_length_days          = 1
integer  :: forecast_length_seconds       = 0
integer  :: verification_interval_seconds = 21600 ! 6 hours
real(r8) :: temporal_coverage_percent     = 100.0 ! all times required
real(r8) :: lonlim1 = MISSING_R8
real(r8) :: lonlim2 = MISSING_R8
real(r8) :: latlim1 = MISSING_R8
real(r8) :: latlim2 = MISSING_R8
logical  :: verbose = .false.
logical  :: debug   = .false.  

namelist /obs_seq_coverage_nml/ obs_sequences, obs_sequence_list, &
              obs_of_interest, textfile_out, netcdf_out, calendar, &
              first_analysis, last_analysis, forecast_length_days, &
              forecast_length_seconds, verification_interval_seconds, &
              temporal_coverage_percent, lonlim1, lonlim2, latlim1, latlim2, &
              verbose, debug

!-----------------------------------------------------------------------
! General purpose variables
!-----------------------------------------------------------------------

integer  :: ifile, iobs, ngood
integer  :: i, j, io, ncunit

type(time_type) :: obs_time, no_time, last_possible_time

character(len=256) :: ncName
character(len=512) :: string1, string2, string3

! ~# of degrees for 1/2 meter at Earth equator 
! 360 deg-earth/(40000 km-earth * 1000m-km)
real(r8), parameter :: HALF_METER = 180.0_r8 / (40000.0_r8 * 1000.0_r8)

!-----------------------------------------------------------------------
! Quantities of interest
!-----------------------------------------------------------------------

integer ::       qc_index   ! copy index of the original qc value
integer ::  dart_qc_index   ! copy index of the DART qc value

character(len=metadatalength), allocatable, dimension(:) :: module_obs_copy_names
character(len=metadatalength), allocatable, dimension(:) ::        obs_copy_names
character(len=metadatalength), allocatable, dimension(:) :: module_qc_copy_names
character(len=metadatalength), allocatable, dimension(:) ::        qc_copy_names

integer, dimension(max_defined_types_of_obs) :: obs_type_inds = 0
real(r8),        allocatable, dimension(:)   :: qc_values
type(time_type), allocatable, dimension(:)   :: all_verif_times
type(time_type), allocatable, dimension(:,:) :: verification_times
real(digits12),  allocatable, dimension(:,:) :: experiment_Tr8
type(time_type) :: verification_stride, half_stride

integer :: num_analyses             ! # of fcsts from first_analysis to last_analysis
integer :: num_verify_per_fcst
integer :: num_verification_times   ! number of verification times - total
integer :: nT_minimum               ! will settle for this many verif times - total
integer :: ilev                     ! index of mandatory level

real(r8), dimension(NUM_MANDATORY_LEVELS) :: mandatory_levels = MISSING_R8 ! pressure level (hPa)

!=======================================================================
! Get the party started
!=======================================================================

call initialize_utilities('obs_seq_coverage')
call static_init_obs_sequence()  ! Initialize the obs sequence module 

call init_obs(obs1, 0, 0)
call init_obs(obs2, 0, 0)
call init_obs_sequence(seq,0,0,0)
no_time = set_time(0,0)

!----------------------------------------------------------------------
! Read the namelist
!----------------------------------------------------------------------

call find_namelist_in_file('input.nml', 'obs_seq_coverage_nml', ncunit)
read(ncunit, nml = obs_seq_coverage_nml, iostat = io)
call check_namelist_read(ncunit, io, 'obs_seq_coverage_nml')

! Record the namelist values used for the run ...
if (do_nml_file()) write(nmlfileunit, nml=obs_seq_coverage_nml)
if (do_nml_term()) write(    *      , nml=obs_seq_coverage_nml)

! Check the user input for sanity
if (temporal_coverage_percent < 100.0_r8) then
   write(string1,*)'namelist: temporal_coverage_percent (',temporal_coverage_percent,&
                   ') must be == 100.0 for now.)' 
   call error_handler(E_MSG, 'obs_seq_coverage', string1, source)
endif

num_input_files =  set_filename_list(obs_sequences, obs_sequence_list,'obs_seq_coverage')

call set_calendar_type(calendar)
call get_calendar_string(calendar)
call SetPressureLevels()

minl = set_location(lonlim1, latlim1, 0.0_r8, VERTISUNDEF)
maxl = set_location(lonlim2, latlim2, 0.0_r8, VERTISUNDEF)

! Determine if the desired observation types exist

TypeLoop : do i = 1,MAX_OBS_NAMES_IN_NAMELIST

   if ( (len(obs_of_interest(i)) == 0) .or. &
            (obs_of_interest(i)  == "") ) exit TypeLoop

   string2 = adjustl(obs_of_interest(i))

   flavor_of_interest = get_index_for_type_of_obs(trim(string2))

   if (flavor_of_interest < 0) then
      write(string1,*)trim(string2),' is not a known observation type.'
      call error_handler(E_ERR, 'obs_seq_coverage', string1, source)
   endif

   if (verbose) write(*,*)trim(string2),' is type ',flavor_of_interest

   obs_type_inds(flavor_of_interest) = 1 ! indicate that we want this one

enddo TypeLoop

if (all(obs_type_inds < 1)) then
      write(string1,*)' No desired observation types of interest specified.'
      write(string2,*)' Check the value of "obs_of_interest".'
      call error_handler(E_ERR, 'obs_seq_coverage', string1, source, text2=string2)
endif

! Set the verification time array (global storage)

call set_required_times(first_analysis, last_analysis, &
          forecast_length_days, forecast_length_seconds, &
          verification_interval_seconds, temporal_coverage_percent)

write(*,*) ! whitespace
write(*,*)'The analysis times (the start of the forecasts) are:'
do i=1,size(verification_times,1)
   write(string1,*)'analysis # ',i,' at '
   call print_date(verification_times(i,1),trim(string1))
enddo
write(*,*) ! whitespace
write(*,*)'At least',nT_minimum,' observations are required from the following:'
do i=1,num_verification_times
   write(string1,*)'verification # ',i,' at '
   call print_date(all_verif_times(i),trim(string1))
   call print_time(all_verif_times(i),trim(string1))
enddo
write(*,*) ! whitespace

last_possible_time = all_verif_times(num_verification_times) + half_stride

! Allocate a hunk of voxels. If we fill this up, we will
! have to create temporary storage, copy, deallocate, reallocate  ...

num_voxels = 0
max_voxels = 4000
call initialize_voxels(max_voxels, voxels)

!====================================================================
!====================================================================


!----------------------------------------------------------------------
! Prepare the variables
!----------------------------------------------------------------------

ObsFileLoop : do ifile = 1, num_input_files
!-----------------------------------------------------------------------

  ! Because of the ability to 'cycle' the ObsFileLoop, we need to
  ! destroy and deallocate at the top of the loop.

   call destroy_obs(obs1)
   call destroy_obs(obs2)
   call destroy_obs_sequence(seq)

   if (allocated(qc_values))      deallocate(qc_values)
   if (allocated(qc_copy_names))  deallocate(qc_copy_names)
   if (allocated(obs_copy_names)) deallocate(obs_copy_names)

   ! Determine the next input filename and check for existence.

   obs_seq_in_file_name = obs_sequences(ifile)

   if ( file_exist(trim(obs_seq_in_file_name)) ) then
      write(string1,*)'opening ', trim(obs_seq_in_file_name)
      call error_handler(E_MSG,'obs_seq_coverage',string1,source)
   else
      write(string1,*)trim(obs_seq_in_file_name),&
                        ' does not exist. Finishing up.'
      call error_handler(E_MSG,'obs_seq_coverage',string1,source)
      write(*,*) ! whitespace
      exit ObsFileLoop
   endif

   ! Read in information about observation sequence so we can allocate
   ! observations. We need info about how many copies, qc values, etc.

   call read_obs_seq_header(obs_seq_in_file_name, &
             num_copies, num_qc, num_obs, max_num_obs, &
             obs_seq_file_id, obs_seq_read_format, pre_I_format, &
             close_the_file = .true.)

   ! Initialize some (individual) observation variables

   call init_obs(obs1, num_copies, num_qc) ! First obs in sequence
   call init_obs(obs2, num_copies, num_qc)

   ! I am taking the observational error variance and making it one of the copies

   if ((num_qc <= 0) .or. (num_copies <=0)) then
      write(string1,*)'need at least 1 qc and 1 observation copy'
      call error_handler(E_ERR,'obs_seq_coverage',string1,source)
   endif

   allocate( obs_copy_names(num_copies), qc_copy_names(num_qc), qc_values(num_qc))

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
   ! Read the entire observation sequence - allocates 'seq' internally
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

   call find_our_copies(seq, qc_index, dart_qc_index)

   if ( ifile == 1 ) then ! record the metadata for comparison

      allocate(module_obs_copy_names(num_copies), &
                module_qc_copy_names(num_qc) )

      do i=1, num_copies
         module_obs_copy_names(i) = obs_copy_names(i)
      enddo
      do i=1, num_qc
         module_qc_copy_names(i) = qc_copy_names(i)
      enddo

   else ! Compare all subsequent files' metadata to the first one

      if (num_copies /= size(module_obs_copy_names)) then
            write(string1,'(''num_copies '',i3,'' does not match '',i3)') &
                              num_copies, size(module_obs_copy_names) 
            call error_handler(E_ERR,'obs_seq_coverage',string1,source)
      endif

      do i = 1,num_copies
         if (trim(obs_copy_names(i)) /= trim(module_obs_copy_names(i))) then
            write(string1,'(''obs copy '',i3,'' from '',a)') i,trim(obs_seq_in_file_name)
            call error_handler(E_MSG,'obs_seq_coverage',string1,source)

            string1 = 'does not match the same observation copy from the first file.'
            write(string2,'(''obs copy >'',a,''<'')') trim(obs_copy_names(i))
            write(string3,'(''expected >'',a,''<'')') trim(module_obs_copy_names(i))
            call error_handler(E_ERR, 'obs_seq_coverage', string1, source, &
                    text2=string2,text3=string3)
         endif
      enddo

      do i = 1,num_qc
         if (trim(qc_copy_names(i)) /= trim(module_qc_copy_names(i))) then
            write(string1,'(''qc copy '',i3,'' from '',a)') i,trim(obs_seq_in_file_name)
            call error_handler(E_MSG,'obs_seq_coverage',string1,source)

            string1 = 'does not match the same qc copy from the first file.'
            write(string2,'(''qc  copy '',a)') trim(qc_copy_names(i))
            write(string3,'(''expected '',a)') trim(module_qc_copy_names(i))
            call error_handler(E_ERR, 'obs_seq_coverage', string1, source, &
                    text2=string2,text3=string3)
         endif
      enddo

   endif

   ngood = 0

   !--------------------------------------------------------------------
   ObservationLoop : do iobs = 1,num_obs
   !--------------------------------------------------------------------

      if (iobs == 1) then
         if ( .not. get_first_obs(seq, obs1) )           &
            call error_handler(E_ERR,'obs_seq_coverage', &
                    'No first observation in sequence.', source)
      else
         call get_next_obs(seq, obs1, obs2, last_ob_flag)
         obs1 = obs2
      endif

      if ( verbose .and. (mod(iobs,10000) == 0) ) &
         write(*,*)'Processing obs ',iobs,' of ',num_obs

      call get_obs_def(obs1,          obs_def)
      flavor   = get_obs_def_type_of_obs(        obs_def)
      obs_time = get_obs_def_time(    obs_def)
      obs_loc  = get_obs_def_location(obs_def)

      call get_qc( obs1, qc_values)

      if (verbose .and. (iobs == 1)) then
         call print_time(obs_time,'First observation time')
         call print_date(obs_time,'First observation date')
         write(*,*) ! whitespace
      endif

      if (obs_time > last_possible_time) exit ObservationLoop

      !-----------------------------------------------------------------
      ! * reject if dart_qc exists and QC is undesirable
      ! * reject if not a type we want [tracked in obs_type_inds(:)]
      ! * reject if not in desired region
      ! * reject if outside vertical region of interest
      !-----------------------------------------------------------------

      if ( dart_qc_index > 0 ) then
         if (qc_values(dart_qc_index) >= 4) cycle ObservationLoop
      endif 

      if ( obs_type_inds(flavor) <= 0 ) cycle ObservationLoop

      if ( .not. is_location_in_region(obs_loc,minl,maxl) ) cycle ObservationLoop

      if ( .not. vertically_desired(obs_loc, ilev) ) cycle ObservationLoop

      ngood = ngood + 1

      ! determine if obs is a new location or time at an existing loc

      voxel_id = find_voxel_location(flavor, obs_loc) 

      if ( voxel_id < 1 ) then
            voxel_id = add_new_voxel(flavor, obs_loc, ilev)
      endif

      if (debug) then
         call write_location(0,obs_loc,'ascii',string1)
         write(*,*)'observation',iobs,'is voxel',voxel_id,trim(string1)
      endif

      if ( time_is_wanted( obs_time, voxel_id, timeindex) ) &
         call update_time( obs_time, voxel_id, timeindex)

   !--------------------------------------------------------------------
   enddo ObservationLoop
   !--------------------------------------------------------------------

enddo ObsFileLoop

! Determine which voxels match the temporal selection requirements

allocate(Desiredvoxels(num_voxels))
Desiredvoxels = .FALSE.
num_out_stat = 0
num_out_total = 0

do i = 1,num_voxels

   voxels(i)%ntimes = 0

   do j = 1,num_verification_times
      if (voxels(i)%times(j) /= no_time) &
         voxels(i)%ntimes = voxels(i)%ntimes + 1
   enddo

   if (voxels(i)%ntimes >= nT_minimum) then
      Desiredvoxels(i) = .TRUE.
      num_out_stat  = num_out_stat + 1
      num_out_total = num_out_total + voxels(i)%ntimes
   endif

   if (debug) write(*,*) 'voxel ID ',i,' has ',voxels(i)%ntimes, ' reports.'

enddo

if (verbose) write(*,*)'There were ',num_out_stat,' voxels matching the input criterion.'
if (verbose) write(*,*)'There were ',num_out_total,' voxels*times matching the input criterion.'

! Output a netCDF file of 'all' observations locations and times.
! Used to explore what is available.

ncName = adjustl(netcdf_out)
ncunit = InitNetCDF(trim(ncName))
call WriteNetCDF(ncunit, trim(ncName), voxels)
call CloseNetCDF(ncunit, trim(ncName))

! if no voxels are selected, do something.

if (num_out_stat < 1) then
   write(string1,*)'No location had at least ',nT_minimum,' reporting times.'
   call error_handler(E_ERR, 'obs_seq_coverage', string1, source)
endif

! Output the file of desired observation locations and times.
! Used to subset the observation sequence files.
call write_obsdefs

!-----------------------------------------------------------------------
! Really, really, done.
!-----------------------------------------------------------------------

call destroy_obs(obs1)
call destroy_obs(obs2)
call destroy_obs_sequence(seq)
call destroy_voxels(voxels)

if (allocated(qc_values))             deallocate(qc_values)
if (allocated(qc_copy_names))         deallocate(qc_copy_names)
if (allocated(obs_copy_names))        deallocate(obs_copy_names)
if (allocated(module_obs_copy_names)) deallocate(module_obs_copy_names)
if (allocated(module_qc_copy_names )) deallocate(module_qc_copy_names )
if (allocated(Desiredvoxels))         deallocate(Desiredvoxels)

call error_handler(E_MSG,'obs_seq_coverage','Finished successfully.',source)
call finalize_utilities()


!======================================================================
CONTAINS
!======================================================================


function find_voxel_location(ObsType, ObsLocation) result(voxel_id)
! Simply try to find a matching lat/lon for an observation type
! The lons/lats get yanked around "a lot" - being converted from ASCII radians
! to r8 degrees to r8 radians to r8 degrees and then checked for "equality".
! So - we're actually just checking to see if the lat/lon is within something
! like 500 cm either direction. Seems like a reasonable definition of 'match'.
!
! By this point the ObsLocation has one of the mandatory vertical levels.

integer,             intent(in) :: ObsType
type(location_type), intent(in) :: ObsLocation
integer                         :: voxel_id

integer :: i
real(r8), dimension(3) :: obslocarray, stnlocarray
real(r8) :: londiff, latdiff, hgtdiff

voxel_id = 0

if (num_voxels == 0) return

obslocarray = get_location(ObsLocation)

FindLoop : do i = 1,num_voxels

   if (ObsType /= voxels(i)%obs_type) cycle FindLoop

   stnlocarray = get_location(voxels(i)%location)

   londiff = abs(obslocarray(1) - stnlocarray(1)) 
   latdiff = abs(obslocarray(2) - stnlocarray(2)) 

   if (is_vertical(ObsLocation, "PRESSURE") ) then
      hgtdiff = abs(obslocarray(3) - stnlocarray(3)) 
   else
      ! if the vertical is undefined or surface, then
      ! the vertical will always match.
      hgtdiff = 0.0_r8
   endif 

   if ( (londiff <= HALF_METER) .and. &
        (latdiff <= HALF_METER) .and. &
        (hgtdiff <  OnePa)    ) then ! within 1 Pa is close enough
      voxel_id = i
      exit FindLoop
   endif

enddo FindLoop

end function find_voxel_location


!============================================================================


function add_new_voxel(ObsType, ObsLocation, ilevel) result(voxel_id)

! Ugh ... if a new location is found, add it. If the voxellist does not have
! enough space, must copy the info to a temporary list, deallocate/reallocate
! copy the info back, and deallocate the temporary list. Ugh. 

integer,             intent(in) :: ObsType
type(location_type), intent(in) :: ObsLocation
integer,             intent(in) :: ilevel
integer                         :: voxel_id

type(voxel_type), allocatable, dimension(:) :: templist
integer :: i

if ( num_voxels >= max_voxels ) then  ! need to make room

   if (verbose) write(*,*)'Doubling number of possible voxels from ', &
                           & num_voxels,' to ',2*max_voxels

   ! Allocate temporary space; Copy. 
   ! Deallocate/nullify existing space.
   ! Double the size of the existing space.
   ! Copy the information back.
   ! Deallocate/nullify the temporary space.
   ! Actually add the new voxel information
   
   ! Allocate the temporary space.
   ! We'll worry about the number of time steps later.
   call initialize_voxels(num_voxels, templist)
   
   ! Copy the information to the temporary space.
   DupLoop1 : do i = 1,num_voxels
   
      templist(i)%obs_type   = voxels(i)%obs_type
      templist(i)%location   = voxels(i)%location
      templist(i)%levelindx  = voxels(i)%levelindx
      templist(i)%first_time = voxels(i)%first_time
      templist(i)%last_time  = voxels(i)%last_time
      templist(i)%ntimes     = voxels(i)%ntimes
   
      ! Make sure the time array is the right size, then copy.
      if (associated(templist(i)%times)) then
         deallocate( templist(i)%times )
         nullify(    templist(i)%times )
      endif
      allocate( templist(i)%times( size(voxels(i)%times) ) )
      templist(i)%times      = voxels(i)%times
   
   enddo DupLoop1
   
   ! Deallocate the voxels, double the array length,
   ! allocate the new space.
   call destroy_voxels(voxels)
   max_voxels = 2 * max_voxels
   call initialize_voxels(max_voxels, voxels)
   
   ! Copy the information BACK to the new space.
   DupLoop2 : do i = 1,num_voxels
   
      voxels(i)%obs_type   = templist(i)%obs_type
      voxels(i)%location   = templist(i)%location
      voxels(i)%levelindx  = templist(i)%levelindx
      voxels(i)%first_time = templist(i)%first_time
      voxels(i)%last_time  = templist(i)%last_time
      voxels(i)%ntimes     = templist(i)%ntimes
   
      if (associated(voxels(i)%times)) then
         deallocate( voxels(i)%times )
         nullify(    voxels(i)%times )
      endif
      allocate( voxels(i)%times( size(templist(i)%times) ) )
      voxels(i)%times = templist(i)%times
   
   enddo DupLoop2
   
   ! Remove the temporary space.
   call destroy_voxels(templist)

endif

! Add the new voxel information.
! Create voxel with nominal (mandatory) vertical value. 
! The vertically_desired() routine replaces the ob vertical value with
! the closest mandatory value, so we're good.

num_voxels = num_voxels + 1
voxel_id   = num_voxels

voxels(voxel_id)%obs_type  = ObsType
voxels(voxel_id)%location  = ObsLocation
voxels(voxel_id)%levelindx = ilevel

if (debug) then
   call write_location(0,ObsLocation,'ascii',string1)
   call write_location(0,voxels(voxel_id)%location,'ascii',string2)
   write(*,*)
   write(*,*)'Added voxel ',voxel_id,' for type ',ObsType
   write(*,*)'observation location', trim(string1)
   write(*,*)'voxel       location', trim(string2)
   write(*,*)'voxel          level', voxels(voxel_id)%levelindx
   write(*,*)
endif

end function add_new_voxel


!============================================================================


function time_is_wanted(ObsTime, voxelid, timeindex)

! The voxel has a list of the observation times closest to the
! verification times. Determine if the observation time is closer to
! the verification time than what we already have.

type(time_type), intent(in)  :: ObsTime
integer,         intent(in)  :: voxelid
integer,         intent(out) :: timeindex
logical                      :: time_is_wanted

type(time_type) :: stndelta, obdelta
integer :: i

timeindex = 0
time_is_wanted = .FALSE.

! the time_minus function always returns a positive difference

TimeLoop : do i = 1,num_verification_times

   obdelta = ObsTime - all_verif_times(i)

   ! If observation is not within half a verification step,
   ! try the next verification time. 
   if (obdelta >= half_stride) cycle TimeLoop 

   ! we must be close now ...
   stndelta = voxels(voxelid)%times(i) - all_verif_times(i)

   ! Check to see if the observation is closer to the verification time
   ! than the one we have.
   if (obdelta < stndelta) then
      if (debug) call print_time(voxels(voxelid)%times(i),'replacing ')
      if (debug) call print_time(ObsTime,'with this observation time')
      timeindex      = i
      time_is_wanted = .TRUE.
      exit TimeLoop
   endif

enddo TimeLoop

end function time_is_wanted


!============================================================================


subroutine update_time(ObsTime, voxelid, timeindex)

! The voxel has a list of the observation times closest to the
! verification times. 
! Add a new time to the voxel registry.

type(time_type), intent(in)  :: ObsTime
integer,         intent(in)  :: voxelid
integer,         intent(in)  :: timeindex

! Update stuff that seems like a good idea, 
! but I don't really know if I'll use it ...
if ( voxels(voxelid)%ntimes == 0 ) then
     voxels(voxelid)%first_time = ObsTime
     voxels(voxelid)%last_time  = ObsTime
endif

if ( voxels(voxelid)%first_time > ObsTime ) &
     voxels(voxelid)%first_time = ObsTime   
if ( voxels(voxelid)%last_time  < ObsTime ) &
     voxels(voxelid)%last_time  = ObsTime   

if (debug) write(*,*)'Stuffing time into voxel ',voxelid,' at timestep ', timeindex

! as long as ntimes /= 0 we are OK.
! When the voxels get written to the netCDF file, count the
! number of non-zero times in the times array for a real count.
voxels(voxelid)%ntimes = voxels(voxelid)%ntimes + 1 

! Stuff the time in the appropriate slot ... finally.
voxels(voxelid)%times(timeindex) = ObsTime

end subroutine update_time


!============================================================================


Function InitNetCDF(fname)
character(len=*), intent(in) :: fname
integer                      :: InitNetCDF

integer :: ncid, i, nlines, linelen
integer :: LineLenDimID, nlinesDimID, stringDimID
integer :: TimeDimID, voxelsDimID, FcstDimID, VerifyDimID
integer :: VarID, FcstVarID, VerifVarID, ExperimentVarID
integer :: nlevDimID, plevelVarID

character(len=8)      :: crdate      ! needed by F90 DATE_AND_TIME intrinsic
character(len=10)     :: crtime      ! needed by F90 DATE_AND_TIME intrinsic
character(len=5)      :: crzone      ! needed by F90 DATE_AND_TIME intrinsic
integer, dimension(8) :: values      ! needed by F90 DATE_AND_TIME intrinsic

character(len=256), allocatable, dimension(:) :: textblock
real(digits12),     allocatable, dimension(:) :: mytimes
integer,            allocatable, dimension(:) :: forecast_length

integer :: ntypes, secs, days, ndims, mylen

integer, dimension(nf90_max_var_dims) :: dimIDs

if(.not. byteSizesOK()) then
    call error_handler(E_ERR,'InitNetCDF', &
   'Compiler does not support required kinds of variables.',source)
endif

InitNetCDF = 0

call nc_check(nf90_create(path = trim(fname), cmode = nf90_clobber, &
         ncid = ncid), 'obs_seq_coverage:InitNetCDF', 'create '//trim(fname))

!----------------------------------------------------------------------------
! Write Global Attributes (mostly namelist input, that sort of thing)
!----------------------------------------------------------------------------

call DATE_AND_TIME(crdate,crtime,crzone,values)
write(string1,'(''YYYY MM DD HH MM SS = '',i4,5(1x,i2.2))') &
               values(1), values(2), values(3), values(5), values(6), values(7)
call nc_check(nf90_put_att(ncid, NF90_GLOBAL, 'creation_date', trim(string1) ), &
           'InitNetCDF', 'put_att creation_date '//trim(fname))

call nc_check(nf90_put_att(ncid, NF90_GLOBAL, 'obs_seq_coverage_source', source), &
           'InitNetCDF', 'put_att obs_seq_coverage_source '//trim(fname))
call nc_check(nf90_put_att(ncid, NF90_GLOBAL, 'min_steps_required', nT_minimum ), &
           'InitNetCDF', 'put_att min_steps_required '//trim(fname))
call nc_check(nf90_put_att(ncid, NF90_GLOBAL, 'forecast_length_days', &
        forecast_length_days ), 'InitNetCDF', 'put_att forecast days '//trim(fname))
call nc_check(nf90_put_att(ncid, NF90_GLOBAL, 'forecast_length_seconds', &
        forecast_length_seconds ), 'InitNetCDF', 'put_att forecast seconds '//trim(fname))
call nc_check(nf90_put_att(ncid, NF90_GLOBAL, 'verification_interval_seconds', &
        verification_interval_seconds ), 'InitNetCDF', 'put_att verif interval '//trim(fname))

! Write all desired observation types.
! As a sanity check - do it from our working array.
ntypes = 0
TYPELOOP : do i = 1,size(obs_type_inds)

   if (obs_type_inds(i) < 1) cycle TYPELOOP

   ntypes = ntypes + 1

   ! create unique netCDF attribute name
   write(string1,'(''obs_of_interest_'',i3.3)') ntypes

   ! decode the index into an observation type name
   string2 = adjustl(get_name_for_type_of_obs(i))

   call nc_check(nf90_put_att(ncid, NF90_GLOBAL, string1, string2 ), &
              'InitNetCDF', 'put_att obs_of_interest '//trim(fname))
enddo TYPELOOP

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
 
call nc_check(nf90_set_fill(ncid, NF90_NOFILL, i),  &
            'InitNetCDF', 'set_fill '//trim(fname))

! the number of voxels

call nc_check(nf90_def_dim(ncid=ncid, &
             name='voxel', len = NF90_UNLIMITED, dimid = voxelsDimID), &
             'InitNetCDF', 'def_dim:voxel '//trim(fname))

call nc_check(nf90_def_var(ncid=ncid, name='voxel', xtype=nf90_int, &
             dimids = (/ voxelsDimID /), varid=VarID), &
             'InitNetCDF', 'voxel:def_var')
call nc_check(nf90_put_att(ncid, VarID, 'long_name', 'desired voxel flag'), &
             'InitNetCDF', 'voxel:long_name')
call nc_check(nf90_put_att(ncid, VarID, 'description', '1 == good voxel'), &
             'InitNetCDF', 'voxel:description')

! the number of verification times

call nc_check(nf90_def_dim(ncid=ncid, &
              name='time', len = num_verification_times, dimid = TimeDimID), &
              'InitNetCDF', 'def_dim:time '//trim(fname))

call nc_check(nf90_def_var(ncid=ncid, name='time', xtype=nf90_double, &
             dimids = (/ TimeDimID /), varid=VarID), &
             'InitNetCDF', 'time:def_var')
call nc_check(nf90_put_att(ncid, VarID, 'long_name', 'verification time'), &
             'InitNetCDF', 'time:long_name')
call nc_check(nf90_put_att(ncid, VarID, 'units',     'days since 1601-1-1'), &
             'InitNetCDF', 'time:put_att units')
call nc_check(nf90_put_att(ncid, VarID, 'calendar',  trim(calendar)), &
             'InitNetCDF', 'time:put_att calendar')

! the number of supported forecasts

call nc_check(nf90_def_dim(ncid=ncid, &
              name='analysisT', len = num_analyses, dimid = FcstDimID), &
              'InitNetCDF', 'def_dim:analysisT '//trim(fname))
call nc_check(nf90_def_var(ncid=ncid, name='analysisT', xtype=nf90_double, &
             dimids = (/ FcstDimID /), varid=FcstVarID), &
             'InitNetCDF', 'analysisT:def_var')
call nc_check(nf90_put_att(ncid, FcstVarID, 'long_name', 'analysis (start) time of each forecast'), &
             'InitNetCDF', 'analysisT:long_name')
call nc_check(nf90_put_att(ncid, FcstVarID, 'units',     'days since 1601-1-1'), &
             'InitNetCDF', 'analysisT:put_att units')
call nc_check(nf90_put_att(ncid, FcstVarID, 'calendar',  trim(calendar)), &
             'InitNetCDF', 'analysisT:put_att calendar')

! the number of verification times per forecast

call nc_check(nf90_def_dim(ncid=ncid, &
              name='forecast_lead',len= num_verify_per_fcst,dimid=VerifyDimID), &
              'InitNetCDF', 'def_dim:forecast_lead '//trim(fname))
call nc_check(nf90_def_var(ncid=ncid, name='forecast_lead', xtype=nf90_int, &
             dimids = (/ VerifyDimID /), varid=VerifVarID), &
             'InitNetCDF', 'forecast_lead:def_var')
call nc_check(nf90_put_att(ncid, VerifVarID, 'long_name', 'current forecast length'), &
             'InitNetCDF', 'forecast_lead:long_name')
call nc_check(nf90_put_att(ncid, VerifVarID, 'units',     'seconds'), &
             'InitNetCDF', 'forecast_lead:put_att units')

! the verification times for each forecast

call nc_check(nf90_def_var(ncid=ncid, name='verification_times', xtype=nf90_double, &
             dimids = (/ VerifyDimID, FcstDimID /), varid=ExperimentVarID), &
             'InitNetCDF', 'experiment:def_var')
call nc_check(nf90_put_att(ncid, ExperimentVarID, 'long_name', 'verification times during each forecast run'), &
             'InitNetCDF', 'experiment:long_name')
call nc_check(nf90_put_att(ncid, ExperimentVarID, 'units',     'days since 1601-1-1'), &
             'InitNetCDF', 'experiment:put_att units')
call nc_check(nf90_put_att(ncid, ExperimentVarID, 'calendar',  trim(calendar)), &
             'InitNetCDF', 'experiment:put_att calendar')
call nc_check(nf90_put_att(ncid, ExperimentVarID, 'rows',  'each forecast'), &
             'InitNetCDF', 'experiment:put_att rows')
call nc_check(nf90_put_att(ncid, ExperimentVarID, 'cols',  'each verification time'), &
             'InitNetCDF', 'experiment:put_att cols')

! the number of levels

call nc_check(nf90_def_dim(ncid=ncid, &
              name='nlevels', len= NUM_MANDATORY_LEVELS, dimid=nlevDimID), &
              'InitNetCDF', 'def_dim:nlevels '//trim(fname))
call nc_check(nf90_def_var(ncid=ncid, name='mandatory_level', xtype=nf90_real, &
              dimids = (/ nlevDimID /), varid=plevelVarID), &
              'InitNetCDF', 'mandatory_level:def_var')
call nc_check(nf90_put_att(ncid, plevelVarID, 'long_name', 'mandatory pressure levels'), &
              'InitNetCDF', 'mandatory_level:long_name')
call nc_check(nf90_put_att(ncid, plevelVarID, 'units', 'Pa'), &
              'InitNetCDF', 'mandatory_level:units')

! write all namelist quantities

call find_textfile_dims('input.nml', nlines, linelen)
allocate(textblock(nlines))
textblock = ''

call nc_check(nf90_def_dim(ncid=ncid, &
              name='linelen', len = len(textblock(1)), dimid = linelenDimID), &
              'InitNetCDF', 'def_dim:linelen '//'input.nml')

call nc_check(nf90_def_dim(ncid=ncid, &
              name='nlines', len = nlines, dimid = nlinesDimID), &
              'InitNetCDF', 'def_dim:nlines '//'input.nml')

call nc_check(nf90_def_dim(ncid=ncid, &
              name='stringlength', len = obstypelength, dimid = StringDimID), &
              'InitNetCDF', 'def_dim:obstypelength '//trim(fname))

! Define the variable to record the input parameters ... the namelist

call nc_check(nf90_def_var(ncid=ncid, name='namelist', xtype=nf90_char, &
             dimids = (/ linelenDimID, nlinesDimID /), varid=VarID), &
             'InitNetCDF', 'namelist:def_var')
call nc_check(nf90_put_att(ncid, VarID, 'long_name', 'input.nml contents'), &
             'InitNetCDF', 'namelist:long_name')

!----------------------------------------------------------------------------
! Define the RECORD variables
!----------------------------------------------------------------------------

! Define the observation type

call nc_check(nf90_def_var(ncid=ncid, name='obs_type', xtype=nf90_char, &
          dimids=(/ StringDimID, voxelsDimID /), varid=VarID), &
          'InitNetCDF', 'obs_type:def_var')
call nc_check(nf90_put_att(ncid, VarID, 'long_name', &
          'observation type string at this voxel'), &
          'InitNetCDF', 'obs_type:put_att long_name')

! let the location module write what it needs to ...

! create a 'location' variable
call nc_write_location_atts(ncid, 0, voxelsDimID)

! Define the mandatory level corresponding to each voxel

call nc_check(nf90_def_var(ncid=ncid, name='voxel_level_index', xtype=nf90_int, &
          dimids=(/ voxelsDimID /), varid=VarID), &
          'InitNetCDF', 'voxel_level_index:def_var')
call nc_check(nf90_put_att(ncid, VarID, 'long_name', &
          'index of the mandatory level of this voxel'), &
          'InitNetCDF', 'voxel_level_index:put_att long_name')
call nc_check(nf90_put_att(ncid, VarID, 'valid_range', &
          (/ 1, NUM_MANDATORY_LEVELS /) ), &
          'InitNetCDF', 'voxel_level_index:put_att valid_range')

! Define the number of observation times

call nc_check(nf90_def_var(ncid=ncid, name='ntimes', xtype=nf90_int, &
          dimids=(/ voxelsDimID /), varid=VarID), &
          'InitNetCDF', 'ntimes:def_var')
call nc_check(nf90_put_att(ncid, VarID, 'long_name', &
          'number of observation times at this voxel'), &
          'InitNetCDF', 'ntimes:put_att long_name')

! Define the first valid observation time

call nc_check(nf90_def_var(ncid=ncid, name='first_time', xtype=nf90_double, &
          dimids=(/ voxelsDimID /), varid=VarID), &
          'InitNetCDF', 'first_time:def_var')
call nc_check(nf90_put_att(ncid, VarID, 'long_name', &
          'first valid observation time at this voxel'), &
          'InitNetCDF', 'first_time:put_att long_name')
call nc_check(nf90_put_att(ncid, VarID, 'units',     'days since 1601-1-1'), &
          'InitNetCDF', 'first_time:put_att units')
call nc_check(nf90_put_att(ncid, VarID, 'calendar',  trim(calendar)), &
          'InitNetCDF', 'first_time:put_att calendar')

! Define the last valid observation time

call nc_check(nf90_def_var(ncid=ncid, name='last_time', xtype=nf90_double, &
          dimids=(/ voxelsDimID /), varid=VarID), &
          'InitNetCDF', 'last_time:def_var')
call nc_check(nf90_put_att(ncid, VarID, 'long_name', &
          'last valid observation time at this voxel'), &
          'InitNetCDF', 'last_time:put_att long_name')
call nc_check(nf90_put_att(ncid, VarID, 'units',     'days since 1601-1-1'), &
          'InitNetCDF', 'last_time:put_att units')
call nc_check(nf90_put_att(ncid, VarID, 'calendar',  trim(calendar)), &
          'InitNetCDF', 'last_time:put_att calendar')

! Define the observation times

call nc_check(nf90_def_var(ncid=ncid, name='ReportTime', xtype=nf90_double, &
          dimids=(/ TimeDimID, voxelsDimID /), varid=VarID), &
          'InitNetCDF', 'ReportTime:def_var')
call nc_check(nf90_put_att(ncid, VarID, 'long_name', 'report time of observation'), &
          'InitNetCDF', 'ReportTime:put_att long_name')
call nc_check(nf90_put_att(ncid, VarID, 'units',     'days since 1601-1-1'), &
          'InitNetCDF', 'ReportTime:put_att units')
call nc_check(nf90_put_att(ncid, VarID, 'calendar',  trim(calendar)), &
          'InitNetCDF', 'ReportTime:put_att calendar')
call nc_check(nf90_put_att(ncid, VarID, 'missing_value', 0.0_digits12), &
          'InitNetCDF', 'ReportTime:put_att missing')
call nc_check(nf90_put_att(ncid, VarID, '_FillValue',    0.0_digits12), &
          'InitNetCDF', 'time:put_att fill_value')

!----------------------------------------------------------------------------
! Leave define mode so we can fill
!----------------------------------------------------------------------------
call nc_check(nf90_enddef(ncid), 'InitNetCDF', 'enddef '//trim(fname))

!----------------------------------------------------------------------------
! Fill the coordinate variables.
! The time variable is filled as time progresses.
!----------------------------------------------------------------------------

call file_to_text('input.nml', textblock)

call nc_check(nf90_inq_varid(ncid, 'namelist', varid=VarID), &
           'InitNetCDF', 'inq_varid:namelist '//trim(fname))
call nc_check(nf90_put_var(ncid, VarID, textblock ), &
           'InitNetCDF', 'put_var:namelist')

deallocate(textblock)

! Fill all possible verification times

call nc_check(nf90_inq_varid(ncid, 'time', varid=VarID), &
           'InitNetCDF', 'inq_varid:time '//trim(fname))

allocate(mytimes(num_verification_times))
mytimes = 0.0_digits12
do i = 1,num_verification_times
   call get_time(all_verif_times(i),secs,days)
   mytimes(i) = days + secs/(60.0_digits12 * 60.0_digits12 * 24.0_digits12)
enddo
call nc_check(nf90_put_var(ncid, VarID, mytimes ), &
           'InitNetCDF', 'put_var:all_verif_times')
deallocate(mytimes)

! Fill forecast start times
! Each of the first N verification times may be used to start a forecast.

allocate(mytimes(num_analyses))
mytimes = 0.0_digits12
do i = 1,num_analyses
   call get_time(all_verif_times(i),secs,days)
   mytimes(i) = days + secs/(60.0_digits12 * 60.0_digits12 * 24.0_digits12)
enddo
call nc_check(nf90_put_var(ncid, FcstVarID, mytimes ), &
           'InitNetCDF', 'put_var:forecast start times')
deallocate(mytimes)

! Fill verification_interval_seconds ... this is the forecast length (seconds) ...

allocate(forecast_length(num_verify_per_fcst))
forecast_length = 0
do i = 1,num_verify_per_fcst
   forecast_length(i) = (i-1)*verification_interval_seconds
enddo
call nc_check(nf90_put_var(ncid, VerifVarID, forecast_length ), &
           'InitNetCDF', 'put_var:forecast_length')
deallocate(forecast_length)

! Fill verification_times into netCDF.
! Each row is a separate forecast.
! Each column is the full date of the forecast/verification time.

call nc_check(nf90_inq_varid(ncid, 'verification_times', varid=VarID ), &
           'InitNetCDF', 'inq_varid:verification_times')

call nc_check(nf90_inquire_variable(ncid, VarID, ndims=ndims, dimids=dimIDs), &
           'InitNetCDF', 'inquire experiment ndims ')

if (debug) then
do i = 1,ndims
   write(string1,*)'verification_times inquire dimid(',i,')'
   call nc_check(nf90_inquire_dimension(ncid, dimIDs(i), len=mylen), &
           'InitNetCDF', string1)
   write(*,*)'verification_times dimension ',i,' has length ',mylen
enddo
endif

call nc_check(nf90_put_var(ncid, ExperimentVarID, experiment_Tr8 ), &
            'InitNetCDF', 'put_var:verification_times')

call nc_check(nf90_put_var(ncid, plevelVarID, mandatory_levels ), &
            'InitNetCDF', 'put_var:plevel')

!----------------------------------------------------------------------------
! Finish up ...
!----------------------------------------------------------------------------

call nc_check(nf90_sync( ncid), 'InitNetCDF', 'sync '//trim(fname))  

InitNetCDF = ncid

end Function InitNetCDF


!============================================================================


subroutine WriteNetCDF(ncid, fname, voxels)
integer,                        intent(in) :: ncid
character(len=*),               intent(in) :: fname
type(voxel_type), dimension(:), intent(in) :: voxels

integer :: DimID, ntimes, voxelindex, days, secs, i
integer, dimension(1) :: istart, icount

integer :: voxelVarID, TimeVarID, NTimesVarID, &
           T1VarID, TNVarID, ObsTypeVarID, IlevVarID

real(digits12), allocatable, dimension(:) :: mytimes
integer, dimension(num_voxels) :: gooduns    ! Cray compiler likes this better

character(len=obstypelength) :: stringOb(1) ! must be an array of strings

!----------------------------------------------------------------------------
! Find the current length of the unlimited dimension so we can add correctly.
!----------------------------------------------------------------------------

call nc_check(nf90_inq_varid(ncid, 'voxel', varid=voxelVarID), &
                   'WriteNetCDF', 'inq_varid:voxelindex '//trim(fname))

gooduns = 0
where(Desiredvoxels) gooduns = 1

call nc_check(nf90_put_var(ncid, voxelVarID, gooduns), &
                   'WriteNetCDF', 'put_var:gooduns '//trim(fname))

call nc_check(nf90_inq_varid(ncid, 'ReportTime', varid=TimeVarID), &
          'WriteNetCDF', 'inq_varid:time '//trim(fname))

call nc_check(nf90_inq_dimid(ncid, 'time', dimid=DimID), &
           'WriteNetCDF', 'inquire time dimid '//trim(fname))
call nc_check(nf90_inquire_dimension(ncid, DimID, len=ntimes), &
           'WriteNetCDF', 'inquire time ntimes '//trim(fname))

call nc_check(nf90_inq_varid(ncid, 'obs_type', varid=ObsTypeVarID), &
          'WriteNetCDF', 'inq_varid:obs_type '//trim(fname))

call nc_check(nf90_inq_varid(ncid, 'ntimes', varid=NTimesVarID), &
          'WriteNetCDF', 'inq_varid:ntimes '//trim(fname))

call nc_check(nf90_inq_varid(ncid, 'first_time', varid=T1VarID), &
          'WriteNetCDF', 'inq_varid:first_time '//trim(fname))

call nc_check(nf90_inq_varid(ncid, 'last_time', varid=TNVarID), &
          'WriteNetCDF', 'inq_varid:last_time '//trim(fname))

call nc_check(nf90_inq_varid(ncid, 'voxel_level_index', varid=IlevVarID), &
          'WriteNetCDF', 'inq_varid:voxel_level_index '//trim(fname))

allocate(mytimes(ntimes))

WriteObs : do voxelindex = 1,num_voxels

   istart(1) = voxelindex
   icount(1) = 1

   ! Must go through Herculean tasks to create 'blank-filled' strings
   ! for the obs_type variable. Dunno why this is so hard.

   string1     = ' '
   stringOb(1) = ' '
   string1 = get_name_for_type_of_obs(voxels(voxelindex)%obs_type)
   write(stringOb(1),'(A)') string1(1:obstypelength)

   call nc_check(nf90_put_var(ncid, ObsTypeVarId, stringOb, &
                start=(/ 1, voxelindex /), count=(/ obstypelength, 1 /) ), &
                'WriteNetCDF', 'put_var:obs_type_string')

   call get_time(voxels(voxelindex)%first_time, secs, days)
   mytimes(1) = days + secs/(60.0_digits12 * 60.0_digits12 * 24.0_digits12)
   call nc_check(nf90_put_var(ncid, T1VarId, (/ mytimes(1) /), &
                start=(/ voxelindex /), count=(/ 1 /) ), &
                'WriteNetCDF', 'put_var:first_time')

   call get_time(voxels(voxelindex)%last_time, secs, days)
   mytimes(1) = days + secs/(60.0_digits12 * 60.0_digits12 * 24.0_digits12)
   call nc_check(nf90_put_var(ncid, TNVarId, (/ mytimes(1) /), &
                start=(/ voxelindex /), count=(/ 1 /) ), &
                'WriteNetCDF', 'put_var:last_time')

   call nc_check(nf90_put_var(ncid, NTimesVarId, (/ voxels(voxelindex)%ntimes /), &
                start=istart, count=icount), 'WriteNetCDF', 'put_var:ntimes')

   call nc_check(nf90_put_var(ncid, IlevVarId, (/ voxels(voxelindex)%levelindx /), &
                start=istart, count=icount), 'WriteNetCDF', 'put_var:voxel_level_index')

   !----------------------------------------------------------------------------
   ! time : fill, write
   !----------------------------------------------------------------------------
   mytimes = 0.0_digits12
   do i = 1,ntimes
      call get_time(voxels(voxelindex)%times(i), secs, days)
      mytimes(i) = days + secs/(60.0_digits12 * 60.0_digits12 * 24.0_digits12)
   enddo

   call nc_check(nf90_put_var(ncid, TimeVarId, mytimes, &
                start=(/ 1, voxelindex /), count=(/ ntimes, 1 /) ), &
                'WriteNetCDF', 'put_var:times')

   !----------------------------------------------------------------------------
   ! Using the location_mod:nc_write_location() routine.
   !----------------------------------------------------------------------------
   call nc_write_location(ncid, voxels(voxelindex)%location, voxelindex, do_vert=.true.)

enddo WriteObs

deallocate(mytimes)

!----------------------------------------------------------------------------
! finished ...
!----------------------------------------------------------------------------

call nc_check(nf90_sync( ncid), 'WriteNetCDF', 'sync '//trim(fname))  

end subroutine WriteNetCDF


!============================================================================


subroutine CloseNetCDF(ncid, fname)
integer,          intent(in) :: ncid
character(len=*), intent(in) :: fname

if ( debug ) write(*,*)'DEBUG --- Closing ',trim(fname)

call nc_check(nf90_sync( ncid), 'CloseNetCDF', 'sync '//trim(fname))  
call nc_check(nf90_close(ncid), 'CloseNetCDF', 'close '//trim(fname))  

end subroutine CloseNetCDF


!============================================================================


subroutine initialize_voxels(Nvoxels, myvoxels)
integer,                                     intent(in)  :: Nvoxels
type(voxel_type), allocatable, dimension(:), intent(out) :: myvoxels

integer :: i

allocate(myvoxels(Nvoxels))

do i = 1,Nvoxels
   myvoxels(i)%obs_type   = 0
   myvoxels(i)%location   = set_location_missing()
   myvoxels(i)%ntimes     = 0
   allocate( myvoxels(i)%times( num_verification_times ) )
   myvoxels(i)%first_time = no_time
   myvoxels(i)%last_time  = no_time
   myvoxels(i)%times      = no_time
enddo

end subroutine initialize_voxels


!============================================================================


subroutine destroy_voxels(myvoxels)
type(voxel_type), allocatable, dimension(:), intent(inout) :: myvoxels

integer :: i,N

N = size(myvoxels)

do i = 1,N
   if (associated(myvoxels(i)%times)) then
      deallocate( myvoxels(i)%times )
      nullify(    myvoxels(i)%times )
   endif
enddo

if (allocated(myvoxels)) deallocate(myvoxels)

end subroutine destroy_voxels


!============================================================================


subroutine write_obsdefs

! Write out a file containing the observation definitions of the desired 
! observations. This file is used to subset the 'big' observation sequence 
! files to harvest the observations used to validate the forecast. 
! Also, print a summary of the voxels we found, etc.

integer :: sec1,secN,day1,dayN
type(obs_def_type) :: obs_def

integer :: iunit, i, j

iunit = get_unit()
open(iunit,file=trim(textfile_out), form='formatted', &
                action='write', position='rewind')

! num_out_total is the result of traversing the list of voxels and times
! and finding the intersection with the user input. How many voxels
! and times fit the requirements.
write(iunit,*)'num_definitions ',num_out_total

call write_type_of_obs_table(iunit, fform='formatted', use_list=obs_type_inds)
call set_obs_def_error_variance(obs_def, MISSING_R8)

write(*,*) ! whitespace
write(*,*)'There are ',num_out_total,' locs/times.'
write(*,*)'Only interested in locations with at least ', nT_minimum,' obs times.' 
write(*,*)'minlon/minlat ', lonlim1, latlim1 
write(*,*)'maxlon/maxlat ', lonlim2, latlim2
write(*,*) ! whitespace

if ( debug ) then
   TYPELOOP : do i = 1,size(obs_type_inds) 
      if (obs_type_inds(i) < 1) cycle TYPELOOP
      string2 = adjustl(get_name_for_type_of_obs(i))
      write(*,*)'i,obs_type_inds(i)',i,obs_type_inds(i),trim(string2)
   enddo TYPELOOP
   write(*,*)
endif

Selections : do i = 1,num_voxels

   if ( .not. Desiredvoxels(i) ) cycle Selections

   call set_obs_def_type_of_obs(    obs_def, voxels(i)%obs_type)
   call set_obs_def_location(obs_def, voxels(i)%location)

   TimeLoop : do j = 1,num_verification_times
      if (voxels(i)%times(j) /= no_time) then
         call set_obs_def_time( obs_def, voxels(i)%times(j))
         call write_obs_def(iunit, obs_def, i, 'formatted')
      endif
   enddo TimeLoop

   if (verbose) then
      call get_time(voxels(i)%first_time,sec1,day1)
      call get_time(voxels(i)%last_time, secN,dayN)
      write(*,'(''voxel '',i6,'' has '',i3,'' obs between ['',&
                  &i7,1x,i5,'' and '',i7,1x,i5,'']'')') &
       i,voxels(i)%ntimes,day1,sec1,dayN,secN
   endif

enddo Selections

close(iunit)

end subroutine write_obsdefs


!============================================================================


subroutine set_required_times(analysis1, analysisN, flen_days, &
          flen_seconds, v_int_seconds, coverage_pcnt )
!-------------------------------------------------------------------------------
! Setting the required times for the verification consists of knowing
! the time of the first analysis through the time of the last analysis plus
! the forecast length (from the last analysis) AND the frequency of the
! verification_interval_seconds. An example with three analysis and 
! 4 verify times per forecast. Results in 7 total times needed to verify.
!
! analysis1      o          o          o          o
!            analysis2      o          o          o          o
!                       analysis3      o          o          o          o
!    |
!    |__________________________________________________________________|
!    T1                                                                 TN
!    V1         V2         V3         V4         V5         V6          V7
!
! Global variables modified/set by this routine:
!
! num_verify_per_fcst      the number of verification times per forecast
! num_analyses             the number of forecast experiments
! num_verification_times   the total number of verification times for ALL experiments
! nT_minimum               the minimum number of verif times we can live with
! verification_times       a matrix of times: rows are forecast runs,
!                                             columns are verification times
! experiment_Tr8           same as verification_times ... just in digits12 format


integer, dimension(:), intent(in)  :: analysis1     ! time of first analysis
integer, dimension(:), intent(in)  :: analysisN     ! time of last  analysis
integer,               intent(in)  :: flen_days     ! forecast length (days)
integer,               intent(in)  :: flen_seconds  ! forecast length (seconds)
integer,               intent(in)  :: v_int_seconds ! verification interval (seconds)
real(r8),              intent(in)  :: coverage_pcnt ! temporal coverage percentage

! declare some local variables

type(time_type) :: flen    ! forecast length
type(time_type) :: T1, TN, TimeMax, thistime, nexttime

integer :: i, j, nsteps, seconds, days

verification_stride = set_time(v_int_seconds,0)         ! global variable
half_stride         = verification_stride / 2
flen                = set_time(flen_seconds, flen_days) ! set the forecast length

T1   = set_date(analysis1(1), analysis1(2), analysis1(3), &
                analysis1(4), analysis1(5), analysis1(6) )

TN   = set_date(analysisN(1), analysisN(2), analysisN(3), &
                analysisN(4), analysisN(5), analysisN(6) )

if ( TN < T1 ) then
   write(string1,*)'namelist: last_analysis must be >= first analysis'
   call error_handler(E_ERR,'set_required_times',string1,source)
endif

! Check to make sure the forecast length is a multiple
! of the verification interval.

nsteps   = flen / verification_stride   ! time_manager_mod:time_divide
nexttime = verification_stride * nsteps ! time_manager_mod:time_scalar_mult 

if (nexttime /= flen) then
   call print_time(verification_stride,'verification stride')
   call print_time(flen,    'original forecast length')
   call print_time(nexttime,'implied  forecast length')
   write(string1,*)'namelist: forecast length is not a multiple of the verification interval'
   write(string2,*)'check forecast_length_[days,seconds] and verification_interval_seconds'
   call error_handler(E_ERR, 'set_required_times', string1, source, text2=string2)
endif

num_verify_per_fcst = nsteps+1   ! SET GLOBAL VALUE

if (verbose) then
   write(*,*) ! a little whitespace
   write(*,*)'There are ',num_verify_per_fcst,' verification times per forecast.' 
endif

! Check to make sure the last analysis is a multiple number
! of verification intervals from the first analysis.
! If we were to launch a forecast at every possible verification 
! time (between analysis1 and analsyisN), how many would that be? 
! Generates a list of potential analysis times.
! FIXME - pathological cases ... 3 hour verification, but
! FIXME - analysis times separated by 1 hour ??? 

if (TN /= T1) then

   thistime = TN - T1
   nsteps   = thistime / verification_stride
   nexttime = verification_stride * nsteps

   if (nexttime /= thistime) then
      call print_time(verification_stride,'verification stride')
      call print_time(thistime,'original analysis duration')
      call print_time(nexttime,'implied  analysis duration')

      ! FIXME - offer a last analysis time ... T1 + nint(steps)*verification_stride

      write(string1,*)'namelist: last analysis time is not a multiple of the verification interval'
      write(string2,*)'check [first,last]_analysis and verification_interval_seconds'
      call error_handler(E_ERR, 'set_required_times', string1, source, text2=string2)
   endif

   num_analyses = nsteps+1   ! SET GLOBAL VALUE
else
   num_analyses = 1   ! SET GLOBAL VALUE
endif

if (verbose) write(*,*)'There are ',num_analyses,' supported forecasts.' 

! Now that we know the start/end/stride are correct, 
! figure out the total number of verification times we need.

TimeMax  = TN + flen ! The time of the last verification
thistime = TimeMax - T1
nsteps   = thistime / verification_stride
nexttime = verification_stride * nsteps

if (nexttime /= thistime) then
   call print_time(verification_stride,'verification stride')
   call print_time(thistime,'total time interval')
   call print_time(nexttime,'implied total time')
   write(string1,*)'bad logic on Tims part. Should not be able to get here.'
   call error_handler(E_ERR,'set_required_times',string1,source)
endif

num_verification_times = nsteps + 1                                ! SET GLOBAL VALUE
nT_minimum = nint(num_verification_times * coverage_pcnt/100.0_r8) ! SET GLOBAL VALUE

if (verbose) then
   write(string1,*)'Need ',num_verification_times, &
      ' verification times for 100.00% coverage.' 
   write(string2,'(''at least '',i5,''  verification times for'',&
      & f7.2,''% coverage.'')') nT_minimum, coverage_pcnt

   call error_handler(E_MSG,'set_required_times',string1,text2=string2)
endif

allocate(all_verif_times(num_verification_times))
allocate(verification_times(num_analyses,num_verify_per_fcst))
allocate(experiment_Tr8(num_verify_per_fcst,num_analyses)) ! TRANSPOSED for netCDF

! Fill the desired verification times
all_verif_times(1) = T1
do i=2,num_verification_times
   all_verif_times(i) = all_verif_times(i-1) + verification_stride
enddo

! Use those desired times to fill the times for each forecast experiment
do i = 1,num_analyses   ! SET GLOBAL VALUE
   verification_times(i, 1:num_verify_per_fcst) = &
      all_verif_times(i:(i+num_verify_per_fcst-1))
enddo

! Use those to create the matching array of real(digits12)

do j = 1,num_verify_per_fcst
do i = 1,num_analyses
   call get_time(verification_times(i,j),seconds,days)
   experiment_Tr8(j,i) = days + seconds/(86400.0_digits12)
enddo
enddo


end subroutine set_required_times


!============================================================================


subroutine find_our_copies(myseq, qcindex, dart_qcindex)
!----------------------------------------------------------------------------
!
type(obs_sequence_type), intent(in)  :: myseq
integer,                 intent(out) :: qcindex
integer,                 intent(out) :: dart_qcindex

integer :: i

     qcindex = -1
dart_qcindex = -1

! Sometimes the first QC field is 'Quality Control' or 'NCEP QC index'
! to me ... they mean the same thing.

QCMetaDataLoop : do i=1, get_num_qc(myseq)
   if(index(get_qc_meta_data(myseq,i),'Quality Control'     ) > 0)      qcindex = i
   if(index(get_qc_meta_data(myseq,i),'NCEP QC index'       ) > 0)      qcindex = i
   if(index(get_qc_meta_data(myseq,i),'DART quality control') > 0) dart_qcindex = i
enddo QCMetaDataLoop

if (      qcindex < 0 ) then
   write(string1,*)'metadata:source Quality Control copyindex not found'
   call error_handler(E_MSG,'find_qc_indices',string1,source)
endif
if ( dart_qcindex < 0 ) then
   write(string1,*)'metadata:DART   Quality Control copyindex not found'
   call error_handler(E_MSG,'find_qc_indices',string1,source)
endif

! Just echo what we know

if (verbose) then
   write(*,*)'QC index',     qcindex,' ',trim(get_qc_meta_data(myseq,     qcindex))
   write(*,*)'QC index',dart_qcindex,' ',trim(get_qc_meta_data(myseq,dart_qcindex))
   write(*,*)
endif

end subroutine find_our_copies


!============================================================================


subroutine setPressureLevels()

! We are using the mandatory levels as defined in:
! http://www.meteor.wisc.edu/~hopkins/aos100/raobdoc.htm
!
! There are 14 MANDATORY pressure levels for radiosonde observations: 
! 1000, 925, 850, 700, 500, 400, 300, 250, 200, 150, 100, 70, 50, and 10 (hPa)

mandatory_levels = (/ 1000.0_r8, 925.0_r8, 850.0_r8, 700.0_r8, 500.0_r8, &
                       400.0_r8, 300.0_r8, 250.0_r8, 200.0_r8, 150.0_r8, &
                       100.0_r8,  70.0_r8,  50.0_r8,  10.0_r8/)

mandatory_levels = mandatory_levels * 100.0_r8  ! convert hPa to Pa

end subroutine setPressureLevels


!======================================================================


function vertically_desired(location,ilevel)

! For obs on pressure levels ... how close to the mandatory levels is close
! enough? A quick examination of 3000+ radiosonde obs revealed that all the 
! pressures were recorded with a precision of 10 Pa (e.g. 110.0, 120.0, ...)
!
! So if the observation is within 1 Pa of a mandatory level, that's good enough.

type(location_type), intent(inout) :: location
integer,             intent(out)   :: ilevel
logical                            :: vertically_desired

integer :: iz
real(r8), dimension(4) :: obslocarray
real(r8) :: zdist

ilevel = MISSING_I

if (is_vertical(location, "UNDEFINED")  .or. &
    is_vertical(location, "SURFACE")) then

   vertically_desired = .true.
   ilevel = 1
   
elseif ( is_vertical(location, "PRESSURE") ) then

   obslocarray(1:3) = get_location(location)
   obslocarray( 4 ) = query_location(location,'which_vert') 

   iz = FindClosestPressureLevel(obslocarray(3))

   zdist = abs(obslocarray(3) - mandatory_levels(iz))

   if ( zdist < OnePa ) then 

      obslocarray(3)     = mandatory_levels(iz)
      location           = set_location( obslocarray )
      vertically_desired = .true.
      ilevel             = iz

   else ! not close enough

      vertically_desired = .false.

   endif

else

   ! do not support height, scale_height, or level
   vertically_desired = .false.

endif

end function vertically_desired


!======================================================================


function FindClosestPressureLevel(obslevel)
real(r8), intent(in) :: obslevel
integer              :: FindClosestPressureLevel

real(r8), dimension(NUM_MANDATORY_LEVELS) :: distances
integer, dimension(1) :: iz

distances = abs(obslevel - mandatory_levels)
iz        = minloc(distances)

FindClosestPressureLevel = iz(1)

end function FindClosestPressureLevel


!======================================================================


end program obs_seq_coverage

