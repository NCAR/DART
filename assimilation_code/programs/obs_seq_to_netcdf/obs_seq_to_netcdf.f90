! DART software - Copyright UCAR. This open source software is provided
! by UCAR, "as is", without charge, subject to all terms of use at
! http://www.image.ucar.edu/DAReS/DART/DART_download

!> converts an observation sequence file to a netCDF file but only retains
!> basic metadata. Metadata to specific observation types is not converted.
!> For example, the location of the radar and direction are not converted,
!> but the location of the observation IS converted.

program obs_seq_to_netcdf

!-----------------------------------------------------------------------
! The programs defines a series of epochs (periods of time) 
!
! All 'possible' obs_kinds are treated separately.
!-----------------------------------------------------------------------

use        types_mod, only : r8, digits12, MISSING_R8
use obs_sequence_mod, only : read_obs_seq, obs_type, obs_sequence_type, get_first_obs, &
                             get_obs_from_key, get_obs_def, get_copy_meta_data, &
                             get_obs_time_range, get_time_range_keys, &
                             get_obs_values, init_obs, assignment(=), static_init_obs_sequence, &
                             get_qc, destroy_obs_sequence, read_obs_seq_header, & 
                             get_last_obs, destroy_obs, get_qc_meta_data
use      obs_def_mod, only : obs_def_type, get_obs_def_error_variance, get_obs_def_time, &
                             get_obs_def_location,  get_obs_def_type_of_obs
use     obs_kind_mod, only : max_defined_types_of_obs, get_name_for_type_of_obs
use     location_mod, only : location_type, write_location, operator(/=), operator(==), &
                             set_location, is_location_in_region, query_location
use  location_io_mod, only : nc_write_location_atts, nc_write_location
use time_manager_mod, only : time_type, get_time, print_time, &
                             get_calendar_string, print_date, &
                             operator(>), operator(<) 
use     schedule_mod, only : schedule_type, set_regular_schedule, get_schedule_length, &
                             get_time_from_schedule
use    utilities_mod, only : file_exist, error_handler, E_ERR, E_MSG, &
                             initialize_utilities, finalize_utilities, nmlfileunit, &
                             find_namelist_in_file, check_namelist_read, &
                             next_file, get_next_filename, find_textfile_dims, &
                             file_to_text, do_nml_file, do_nml_term
use netcdf_utilities_mod, only : nc_check

use typeSizes
use netcdf

implicit none

character(len=*), parameter :: source = 'obs_seq_to_netcdf.f90'

!---------------------------------------------------------------------
!---------------------------------------------------------------------

integer, parameter :: stringlength = 32

!---------------------------------------------------------------------
! variables associated with the observation
!---------------------------------------------------------------------

type(obs_sequence_type) :: seq
type(obs_type)          :: observation, next_obs
type(obs_type)          :: obs1, obsN
type(obs_def_type)      :: obs_def
type(location_type)     :: obs_loc, minl, maxl

character(len=256) :: obs_seq_in_file_name
character(len=256), allocatable, dimension(:) :: obs_seq_filenames

real(r8)            :: obs_err_var

integer :: flavor ! THIS IS THE (global) 'KIND' in the obs_def_mod list. 
integer :: num_copies, num_qc, num_obs, max_num_obs, obs_seq_file_id

integer :: num_obs_kinds

character(len=stringlength) :: obs_seq_read_format
logical :: pre_I_format

logical :: out_of_range, is_there_one, keeper

!-----------------------------------------------------------------------
! Namelist with (some scalar) default values
!-----------------------------------------------------------------------

character(len=256) :: obs_sequence_name = 'obs_seq.final'
character(len=256) :: obs_sequence_list = ''

real(r8) :: lonlim1= MISSING_R8, lonlim2= MISSING_R8
real(r8) :: latlim1= MISSING_R8, latlim2= MISSING_R8 

logical :: debug = .false.   ! undocumented ... on purpose
logical :: verbose = .false.
logical :: append_to_netcdf = .false.

namelist /obs_seq_to_netcdf_nml/ obs_sequence_name, obs_sequence_list, &
                                 lonlim1, lonlim2, latlim1, latlim2, &
                                 verbose, append_to_netcdf, debug

!-----------------------------------------------------------------------
! Quantities of interest
!-----------------------------------------------------------------------

integer, parameter :: Ncopies = 1
integer :: allNcopies
character(len=stringlength), dimension(Ncopies) :: copy_names = &
   (/ 'observation error variance' /)

character(len=stringlength), allocatable, dimension(:) :: module_obs_copy_names
character(len=stringlength), allocatable, dimension(:) :: module_qc_copy_names
character(len=stringlength), allocatable, dimension(:) :: obs_copy_names, qc_copy_names
character(len=stringlength), dimension(max_defined_types_of_obs) :: my_obs_kind_names

real(r8), allocatable, dimension(:) :: qc, copyvals
real(r8), allocatable, dimension(:) :: obscopies

integer,        dimension(2) :: key_bounds
integer,        allocatable, dimension(:)   ::       keys
real(digits12), allocatable, dimension(:)   ::  obs_times
integer,        allocatable, dimension(:)   ::  obs_types, obs_keys
real(r8),       allocatable, dimension(:,:) :: obs_copies
integer,        allocatable, dimension(:,:) ::  qc_copies
type(location_type), allocatable, dimension(:) ::  locations
integer,        allocatable, dimension(:)   :: which_vert
integer,        allocatable, dimension(:)   :: total_obs_in_epoch

!-----------------------------------------------------------------------
! General purpose variables
!-----------------------------------------------------------------------

integer  :: iepoch, ifile, num_obs_in_epoch, ngood

integer  :: i, io, obsindex, ncunit
integer  :: Nepochs
! Keep track of epoch files created during this run of obs_seq_to_netcdf
logical, allocatable, dimension(:) :: epoch_file_created 

type(schedule_type) :: schedule
type(time_type) :: TimeMin, TimeMax    ! of the entire period of interest
type(time_type) :: beg_time, end_time  ! of the particular bin
type(time_type) :: seqT1, seqTN        ! first,last time in entire observation sequence
type(time_type) :: obs_time

real(digits12)  :: mytime
integer         :: seconds, days

character(len=stringlength) :: ncName, calendarstring
character(len=512) :: string1, string2, string3

!=======================================================================
! Get the party started
!=======================================================================

call initialize_utilities('obs_seq_to_netcdf')
call static_init_obs_sequence()  ! Initialize the obs sequence module 

call init_obs(       obs1, 0, 0) ! I am initialiazing these obs
call init_obs(       obsN, 0, 0) ! simply to make the logic in 
call init_obs(observation, 0, 0) ! ObsFileLoop simpler. This way we
call init_obs(   next_obs, 0, 0) ! can destroy at the top of the loop.

!----------------------------------------------------------------------
! FIXME : at some point, I'd like to restrict this to just the ones
!         in the obs_seq file ... not the ones supported. Problem is
!         that you have to read all the obs sequence files and process
!         the region information before you know if you have any to
!         output. Maybe allocate more space in the netCDF header and
!         re-enter DEFINE mode after all is said and done ...
!----------------------------------------------------------------------

num_obs_kinds = max_defined_types_of_obs

do i = 1,max_defined_types_of_obs
   my_obs_kind_names(i) = get_name_for_type_of_obs(i)
enddo

!----------------------------------------------------------------------
! Read the namelist
!----------------------------------------------------------------------

call find_namelist_in_file('input.nml', 'obs_seq_to_netcdf_nml', ncunit)
read(ncunit, nml = obs_seq_to_netcdf_nml, iostat = io)
call check_namelist_read(ncunit, io, 'obs_seq_to_netcdf_nml')

! Record the namelist values used for the run ...
if (do_nml_file()) write(nmlfileunit, nml=obs_seq_to_netcdf_nml)
if (do_nml_term()) write(    *      , nml=obs_seq_to_netcdf_nml)

if ((obs_sequence_name /= '') .and. (obs_sequence_list /= '')) then
   write(string1,*)'specify "obs_sequence_name" -OR- "obs_sequence_list"'
   write(string2,*)'set other to an empty string ... i.e. ""'
   call error_handler(E_ERR, 'obs_seq_to_netcdf', string1, source, text2=string2)
endif

!----------------------------------------------------------------------
! SetSchedule rectifies user input and the final binning sequence.
!----------------------------------------------------------------------

call set_regular_schedule(schedule) ! also sets calendar type

Nepochs = get_schedule_length(schedule)

allocate(total_obs_in_epoch(Nepochs))
total_obs_in_epoch = 0
allocate(epoch_file_created(Nepochs))
epoch_file_created(:) = .false.

call get_time_from_schedule(TimeMin, schedule,       1, 1)
call get_time_from_schedule(TimeMax, schedule, Nepochs, 2)
call get_calendar_string(calendarstring)

minl = set_location( (/ lonlim1, latlim1, -huge(0.0_r8), 1.0_r8 /))
maxl = set_location( (/ lonlim2, latlim2,  huge(0.0_r8), 1.0_r8 /))

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
   call destroy_obs(obsN)
   call destroy_obs(observation)
   call destroy_obs(next_obs)
   call destroy_obs_sequence(seq)  ! hopefully destroys OK without preallocate

   if (allocated(qc)) deallocate( qc, copyvals, obs_copy_names, qc_copy_names, obscopies )

   ! Determine the next input filename ... 

   if (obs_sequence_list == '') then
      obs_seq_in_file_name = next_file(obs_sequence_name,ifile)
   else
      obs_seq_in_file_name = get_next_filename(obs_sequence_list,ifile)
      if (obs_seq_in_file_name == '') exit ObsFileLoop
   endif

   if ( file_exist(trim(obs_seq_in_file_name)) ) then
      write(string1,*)'opening ', trim(obs_seq_in_file_name)
      call error_handler(E_MSG,'obs_seq_to_netcdf',string1,source)
   else
      write(string1,*)trim(obs_seq_in_file_name),&
                        ' does not exist. Finishing up.'
      call error_handler(E_MSG,'obs_seq_to_netcdf',string1,source)
      exit ObsFileLoop
   endif

   ! Read in information about observation sequence so we can allocate
   ! observations. We need info about how many copies, qc values, etc.

   obs_seq_in_file_name     = trim(obs_seq_in_file_name) ! Lahey requirement
   obs_seq_filenames(ifile) = trim(obs_seq_in_file_name)

   call read_obs_seq_header(obs_seq_in_file_name, &
             num_copies, num_qc, num_obs, max_num_obs, &
             obs_seq_file_id, obs_seq_read_format, pre_I_format, &
             close_the_file = .true.)

   ! Initialize some (individual) observation variables

   call init_obs(       obs1, num_copies, num_qc)   ! First obs in sequence
   call init_obs(       obsN, num_copies, num_qc)   ! Last  obs in sequence
   call init_obs(observation, num_copies, num_qc)   ! current obs
   call init_obs(   next_obs, num_copies, num_qc)   ! duh ...

   ! I am taking the observational error variance and making it one of the copies

   allNcopies = num_copies + Ncopies

   if ((num_qc <= 0) .or. (num_copies <=0)) then
      write(string1,*)'need at least 1 qc and 1 observation copy'
      call error_handler(E_ERR,'obs_seq_to_netcdf',string1,source)
   endif

   allocate( copyvals(allNcopies), &
       obs_copy_names(allNcopies), &
        qc_copy_names(num_qc),     &
                   qc(num_qc),     &
            obscopies(num_copies))

   if ( debug ) then
      write(*,*)
      write(*,*)'num_copies          is ',num_copies
      write(*,*)'num_qc              is ',num_qc
      write(*,*)'num_obs             is ',num_obs
      write(*,*)'max_num_obs         is ',max_num_obs
      write(*,*)'obs_seq_read_format is ',trim(obs_seq_read_format)
      write(*,*)'pre_I_format        is ',pre_I_format
      write(*,*)
   endif

   !--------------------------------------------------------------------
   ! Read the entire observation sequence - allocates 'seq' internally
   !--------------------------------------------------------------------

   call read_obs_seq(obs_seq_in_file_name, 0, 0, 0, seq)

   do i=1, num_copies
         string1 = trim(get_copy_meta_data(seq,i))//'                          '
         obs_copy_names(i) = string1(1:stringlength)
   enddo
   do i=1, Ncopies
         obs_copy_names(num_copies+i) = trim(copy_names(i))
   enddo
   do i=1, num_qc
         string1 = trim(get_qc_meta_data(seq,i))//'                          '
         qc_copy_names(i) = string1(1:stringlength)
   enddo

   if ( ifile == 1 ) then ! record the metadata for comparison

      allocate(module_obs_copy_names(allNcopies), &
                module_qc_copy_names(num_qc) )

      do i=1, num_copies
         string1 = trim(get_copy_meta_data(seq,i))//'                          '
         module_obs_copy_names(i) = string1(1:stringlength)
      enddo
      do i=1, Ncopies
         module_obs_copy_names(num_copies+i) = trim(copy_names(i))
      enddo
      do i=1, num_qc
         string1 = trim(get_qc_meta_data(seq,i))//'                          '
         module_qc_copy_names(i) = string1(1:stringlength)
      enddo

   else ! Compare all subsequent files' metadata to the first one

      do i = 1,allNcopies
         if (trim(obs_copy_names(i)) /= trim(module_obs_copy_names(i))) then
            write(string1,'(''obs copy '',i3,'' from '',a)') i,trim(obs_seq_in_file_name)
            string2 = 'does not match the same observation copy from the first file.'
            write(string3,'(''obs copy >'',a,''< expected >'',a,''<'')') &
                        trim(obs_copy_names(i)), trim(module_obs_copy_names(i))
            call error_handler(E_ERR,'obs_seq_to_netcdf',string1, &
                source,text2=string2,text3=string3)

         endif
      enddo

      do i = 1,num_qc
         if (trim(qc_copy_names(i)) /= trim(module_qc_copy_names(i))) then
            write(string1,'(''qc copy '',i3,'' from '',a)') i,trim(obs_seq_in_file_name)
            string2 = 'does not match the same qc copy from the first file.'
            write(string3,'(''qc  copy >'',a,''< expected >'',a,''<'')') &
                         trim(qc_copy_names(i)), trim(module_qc_copy_names(i))
            call error_handler(E_ERR,'obs_seq_to_netcdf',string1, &
                source,text2=string2,text3=string3)
         endif
      enddo

   endif

   !--------------------------------------------------------------------
   ! Determine the time encompassed in the observation sequence.
   !--------------------------------------------------------------------

   is_there_one = get_first_obs(seq, obs1)
   if ( .not. is_there_one ) then
      call error_handler(E_ERR,'obs_seq_to_netcdf','No first observation  in sequence.', &
      source)
   endif
   call get_obs_def(obs1,   obs_def)
   seqT1 = get_obs_def_time(obs_def)

   is_there_one = get_last_obs(seq, obsN)
   if ( .not. is_there_one ) then
      call error_handler(E_ERR,'obs_seq_to_netcdf','No last observation in sequence.', &
      source)
   endif
   call get_obs_def(obsN,   obs_def)
   seqTN = get_obs_def_time(obs_def)

   if ( verbose ) then
      call print_time(  seqT1,'First observation time')
      call print_time(TimeMin,'TimeMin     from input')
      call print_time(  seqTN,'Last  observation time')
      call print_time(TimeMax,'TimeMax     from input')
      write(*,*)''
      call print_date(  seqT1,'First observation date')
      call print_date(TimeMin,'DateMin     from input')
      call print_date(  seqTN,'Last  observation date')
      call print_date(TimeMax,'DateMax     from input')
      write(*,*)''
   endif

   !--------------------------------------------------------------------
   ! If the last observation is before the period of interest, move on.
   !--------------------------------------------------------------------

   if ( seqTN < TimeMin ) then
      if (verbose) write(*,*)'seqTN < TimeMin ... trying next file.'
      cycle ObsFileLoop
   else
      if (verbose) write(*,*)'seqTN > TimeMin ... using ', trim(obs_seq_in_file_name)
   endif

   !--------------------------------------------------------------------
   ! If the first observation is after the period of interest, finish.
   !--------------------------------------------------------------------

   if ( seqT1 > TimeMax ) then
      if (verbose) write(*,*)'seqT1 > TimeMax ... stopping.'
      exit ObsFileLoop
   else
      if (verbose) write(*,*)'seqT1 < TimeMax ... using ',trim(obs_seq_in_file_name)
   endif

   !====================================================================
   EpochLoop : do iepoch = 1, Nepochs
   !====================================================================

      call get_time_from_schedule(beg_time, schedule, iepoch, 1)
      call get_time_from_schedule(end_time, schedule, iepoch, 2)

      call get_obs_time_range(seq, beg_time, end_time, key_bounds, &
                  num_obs_in_epoch, out_of_range)

      if( num_obs_in_epoch == 0 ) then
         if (verbose) write(*,*)' No observations in epoch ',iepoch,' cycling ...'
         cycle EpochLoop
      endif

      write(*,*)'num_obs_in_epoch (', iepoch, ') = ', num_obs_in_epoch

      allocate(      keys(            num_obs_in_epoch), &
                obs_times(            num_obs_in_epoch), &
                obs_types(            num_obs_in_epoch), &
                 obs_keys(            num_obs_in_epoch), &
               which_vert(            num_obs_in_epoch), &
               obs_copies(allNcopies, num_obs_in_epoch), &
                qc_copies(    num_qc, num_obs_in_epoch), &
                locations(            num_obs_in_epoch))

      call get_time_range_keys(seq, key_bounds, num_obs_in_epoch, keys)

      !-----------------------------------------------------------------
      ! Append epoch number to output file names.
      ! Clear out any existing files unless we are TRYING to append.
      !-----------------------------------------------------------------

      write(ncName,'(''obs_epoch_'',i3.3,''.nc'')')iepoch

      if ( file_exist(ncName) .and. append_to_netcdf ) then
         ncunit = NC_Compatibility_Check(ncName, iepoch)
      else 
         if (.not. epoch_file_created(iepoch)) then
            write(string1,*) 'creating file ', ncName
            call error_handler(E_MSG,'obs_seq_to_netcdf',string1,source)
            ncunit = InitNetCDF(ncName, iepoch)
            epoch_file_created(iepoch) = .true.
         else
            ncunit = NC_Compatibility_Check(ncName, iepoch)
         endif
      endif

      ngood = 0
      !-----------------------------------------------------------------
      ObservationLoop : do obsindex = 1, num_obs_in_epoch
      !-----------------------------------------------------------------

         ! 'flavor' is from the 'master list' in the obs_kind_mod.f90
         ! each obs_seq.final file has their own private kind - which
         ! gets mapped to the 'master list', if you will.

         if ( verbose .and. (mod(obsindex,10000) == 0) ) then
            write(*,*)'Processing obs ',obsindex,' of ',num_obs_in_epoch
         endif

         call get_obs_from_key(seq, keys(obsindex), observation)
         call get_obs_values(observation, obscopies)
         call get_obs_def(observation, obs_def)
         call get_qc(observation, qc)

         flavor      = get_obs_def_type_of_obs(obs_def)
         obs_time    = get_obs_def_time(obs_def)
         obs_loc     = get_obs_def_location(obs_def)

         ! replace missing values with NetCDF missing value
         where (obscopies == MISSING_R8 ) obscopies = NF90_FILL_DOUBLE

         ! paste on the observational error variance
         obs_err_var = get_obs_def_error_variance(obs_def)

         copyvals = (/ obscopies, obs_err_var /)

         call get_time(obs_time,seconds,days)
         mytime   = days + seconds/86400.0_digits12

         !--------------------------------------------------------------
         ! We have one Region of interest
         !--------------------------------------------------------------

         keeper = is_location_in_region( obs_loc, minl, maxl ) 

         if ( .not. keeper ) cycle ObservationLoop

         ngood = ngood + 1

         !--------------------------------------------------------------
         ! Summary of observation knowledge at this point
         ! can replace the hardcoded 6 when the write_location function
         ! can write to a character string.
         !--------------------------------------------------------------

         if ( debug ) then
            write(6, *)'observation # ',obsindex
            write(6, *)'key           ',keys(obsindex)
            write(6, *)'obs_flavor    ',flavor
            call write_location(6, obs_loc, 'ascii')
            write(6, *)'copyvals      ',copyvals
            write(6, *)'qc            ',qc
         endif

         obs_copies(:,ngood) = copyvals
          qc_copies(:,ngood) = nint(qc)
          locations(  ngood) = obs_loc
          obs_times(  ngood) = mytime
          obs_types(  ngood) = flavor
           obs_keys(  ngood) = keys(obsindex)
         which_vert(  ngood) = nint(query_location(obs_loc))

      !-----------------------------------------------------------------
      enddo ObservationLoop
      !-----------------------------------------------------------------

      total_obs_in_epoch(iepoch) = total_obs_in_epoch(iepoch) + ngood

      if ( ngood > 0 ) call WriteNetCDF(ncunit, ncname, ngood, obs_copies, &
                       qc_copies, locations, obs_times, obs_types, obs_keys) 

      call CloseNetCDF(ncunit, ncname)

      deallocate(keys, obs_times, obs_types, obs_keys, which_vert, &
                 obs_copies, qc_copies, locations)

   enddo EpochLoop

   if (verbose) write(*,*)'End of EpochLoop for ',trim(obs_seq_in_file_name)

enddo ObsFileLoop

!-----------------------------------------------------------------------
! Block to check if any output file(s) was(were) created ... 
!-----------------------------------------------------------------------
   
do iepoch = 1, Nepochs

   if ( total_obs_in_epoch(iepoch) == 0 ) then
      write(string1,'(''epoch '',i6,'' has no observations. verbose=.TRUE. may help'')')iepoch
      write(string2,*)'schedule_nml might not match observation times in files.'
      call error_handler(E_MSG, 'obs_seq_to_netcdf', string1, &
              source, text2=string2, text3=' ')
   endif

enddo
  
! Nancy: "I had put the wrong year in the schedule namelist
! and so no obs were within the window.  but
! without the verbose flag set to true i didn't get any
! message about why.  it would have helped to see
! a terse line like 'no output file created.  no obs within window'
! (a message, not error) before it exits in this case.
!
! turning on the verbose flag gave me all kinds of
! other info and the problem was clear - but i didn't
! have a clue about why there was no output until then."

!-----------------------------------------------------------------------
! Really, really, done.
!-----------------------------------------------------------------------

call destroy_obs(obs1)
call destroy_obs(obsN)
call destroy_obs(observation)
call destroy_obs(next_obs)
call destroy_obs_sequence(seq)

if (allocated(qc)) deallocate( qc, copyvals, obs_copy_names, qc_copy_names, obscopies )

if (allocated(module_obs_copy_names)) &
   deallocate(module_obs_copy_names, module_qc_copy_names)

deallocate(obs_seq_filenames)
deallocate(epoch_file_created)

call error_handler(E_MSG,'obs_seq_to_netcdf','Finished successfully.',source)
call finalize_utilities()


!======================================================================
CONTAINS
!======================================================================

Function InitNetCDF(fname, ibin)
character(len=*), intent(in) :: fname
integer,          intent(in) :: ibin
integer                      :: InitNetCDF

integer :: ncid, i, indx1, nlines, linelen
integer :: LineLenDimID, nlinesDimID, stringDimID
integer :: ObsCopyDimID, QCCopyDimID
integer ::   TypesDimID
integer ::  ObsNumDimID
integer ::   VarID

character(len=8)      :: crdate      ! needed by F90 DATE_AND_TIME intrinsic
character(len=10)     :: crtime      ! needed by F90 DATE_AND_TIME intrinsic
character(len=5)      :: crzone      ! needed by F90 DATE_AND_TIME intrinsic
integer, dimension(8) :: values      ! needed by F90 DATE_AND_TIME intrinsic

character(len=256), allocatable, dimension(:) :: textblock

real(digits12)  :: epoch_edges(2)
integer         :: seconds, days
type(time_type) :: mytime

if(.not. byteSizesOK()) then
    call error_handler(E_ERR,'InitNetCDF', &
   'Compiler does not support required kinds of variables.',source)
endif

InitNetCDF = 0

call get_time_from_schedule(mytime,schedule,ibin,1)
call get_time(mytime,seconds,days)
epoch_edges(1) = days + seconds/86400.0_digits12

call get_time_from_schedule(mytime,schedule,ibin,2)
call get_time(mytime,seconds,days)
epoch_edges(2) = days + seconds/86400.0_digits12

call nc_check(nf90_create(path = trim(fname), cmode = nf90_share, &
         ncid = ncid), 'obs_seq_to_netcdf:InitNetCDF', 'create '//trim(fname))

write(string1,*)trim(ncName), ' is fortran unit ',ncid
call error_handler(E_MSG,'InitNetCDF',string1,source)

!----------------------------------------------------------------------------
! Write Global Attributes 
!----------------------------------------------------------------------------

call DATE_AND_TIME(crdate,crtime,crzone,values)
write(string1,'(''YYYY MM DD HH MM SS = '',i4,5(1x,i2.2))') &
               values(1), values(2), values(3), values(5), values(6), values(7)
call nc_check(nf90_put_att(ncid, NF90_GLOBAL, 'creation_date', trim(string1) ), &
           'InitNetCDF', 'put_att creation_date '//trim(fname))

call nc_check(nf90_put_att(ncid, NF90_GLOBAL, 'obs_seq_to_netcdf_source', source ), &
           'InitNetCDF', 'put_att obs_seq_to_netcdf_source '//trim(fname))

! write all observation sequence files used
FILEloop : do i = 1,SIZE(obs_seq_filenames)

  indx1 = index(obs_seq_filenames(i),'null')

  if (indx1 > 0) exit FILEloop

  write(string1,'(''obs_seq_file_'',i3.3)')i
  call nc_check(nf90_put_att(ncid, NF90_GLOBAL, &
         trim(string1), trim(obs_seq_filenames(i)) ), &
         'InitNetCDF', 'put_att:filenames')

enddo FILEloop

!----------------------------------------------------------------------------
! Define the dimensions
! Set nofill mode - supposed to be performance gain
!----------------------------------------------------------------------------
 
call nc_check(nf90_set_fill(ncid, NF90_NOFILL, i),  &
            'obs_seq_to_netcdf:InitNetCDF', 'set_nofill '//trim(fname))

! write all namelist quantities

call find_textfile_dims('input.nml', nlines, linelen)
allocate(textblock(nlines))
textblock = ''

call nc_check(nf90_def_dim(ncid=ncid, &
              name="linelen", len = len(textblock(1)), dimid = linelenDimID), &
              'InitNetCDF', 'def_dim:linelen '//'input.nml')

call nc_check(nf90_def_dim(ncid=ncid, &
              name="nlines", len = nlines, dimid = nlinesDimID), &
              'InitNetCDF', 'def_dim:nlines '//'input.nml')

call nc_check(nf90_def_dim(ncid=ncid, &
              name='stringlength', len = stringlength, dimid = StringDimID), &
              'InitNetCDF', 'def_dim:stringlength '//trim(fname))

call nc_check(nf90_def_dim(ncid=ncid, &
              name='copy', len = allNcopies, dimid = ObsCopyDimID), &
              'InitNetCDF', 'def_dim:copy '//trim(fname))

call nc_check(nf90_def_dim(ncid=ncid, &
              name='qc_copy', len = num_qc, dimid = QCCopyDimID), &
              'InitNetCDF', 'def_dim:qc_copy '//trim(fname))

call nc_check(nf90_def_dim(ncid=ncid, &
              name='ObsTypes', len = num_obs_kinds, dimid = TypesDimID), &
              'InitNetCDF', 'def_dim:ObsTypes '//trim(fname))

call nc_check(nf90_def_dim(ncid=ncid, &
              name='ObsIndex', len = NF90_UNLIMITED, dimid = ObsNumDimID), &
              'InitNetCDF', 'def_dim:ObsIndex '//trim(fname))

!----------------------------------------------------------------------------
! Define the static variables
!----------------------------------------------------------------------------

! Define the types of observation quantities 

call nc_check(nf90_def_var(ncid=ncid, name='copy', xtype=nf90_int, &
             dimids=ObsCopyDimID, varid=VarID), &
             'InitNetCDF', 'copy:def_var')
call nc_check(nf90_put_att(ncid, VarID, 'explanation', 'see CopyMetaData'), &
             'InitNetCDF', 'copy:explanation')

! Define the types of qc quantities

call nc_check(nf90_def_var(ncid=ncid, name='qc_copy', xtype=nf90_int, &
             dimids=QCCopyDimID, varid=VarID), &
             'InitNetCDF', 'qc_copy:def_var')
call nc_check(nf90_put_att(ncid, VarID, 'explanation', 'see QCMetaData'), &
             'InitNetCDF', 'qc_copy:explanation')

! Define the observation type

call nc_check(nf90_def_var(ncid=ncid, name='ObsTypes', xtype=nf90_int, &
          dimids=TypesDimID, varid=VarID), &
          'InitNetCDF', 'ObsTypes:def_var')
call nc_check(nf90_put_att(ncid, VarID, 'explanation', 'see ObsTypesMetaData'), &
          'InitNetCDF', 'ObsTypes:explanation')

! Define the character strings

call nc_check(nf90_def_var(ncid=ncid, name='ObsTypesMetaData', xtype=nf90_char, &
             dimids=(/ StringDimID, TypesDimID /), varid=VarID), &
             'InitNetCDF', 'typesmeta:def_var')
call nc_check(nf90_put_att(ncid, VarID, 'long_name', 'DART observation types'), &
             'InitNetCDF', 'typesmeta:long_name')
call nc_check(nf90_put_att(ncid, VarID, 'comment', &
         'table relating integer to observation type string'), &
             'InitNetCDF', 'typesmeta:comment')

! Define the character strings for the QC flags

call nc_check(nf90_def_var(ncid=ncid, name='QCMetaData', xtype=nf90_char, &
             dimids=(/ StringDimID, QCCopyDimID /), varid=VarID), &
             'InitNetCDF', 'qcmeta:def_var')
call nc_check(nf90_put_att(ncid, VarID, 'long_name', 'quantity names'), &
             'InitNetCDF', 'qcmeta:long_name')

! Define the character strings for the quantities recorded

call nc_check(nf90_def_var(ncid=ncid, name='CopyMetaData', xtype=nf90_char, &
             dimids=(/ StringDimID, ObsCopyDimID /), varid=VarID), &
             'InitNetCDF', 'copymeta:def_var')
call nc_check(nf90_put_att(ncid, VarID, 'long_name', 'quantity names'), &
             'InitNetCDF', 'copymeta:long_name')

! Define the variable to record the input parameters ... the namelist

call nc_check(nf90_def_var(ncid=ncid, name="namelist", xtype=nf90_char, &
             dimids = (/ linelenDimID, nlinesDimID /), varid=VarID), &
             'InitNetCDF', 'namelist:def_var')
call nc_check(nf90_put_att(ncid, VarID, "long_name", "input.nml contents"), &
             'InitNetCDF', 'namelist:long_name')

!----------------------------------------------------------------------------
! Define the RECORD variables
!----------------------------------------------------------------------------

! Define the observation number coordinate variable (UNLIMITED DIMENSION) 

call nc_check(nf90_def_var(ncid=ncid, name='ObsIndex', xtype=nf90_int, &
          dimids=(/ ObsNumDimID /), varid=VarID), &
            'InitNetCDF', 'obsindex:def_var')
call nc_check(nf90_put_att(ncid, VarID, 'long_name', 'observation index'), &
          'InitNetCDF', 'obsindex:long_name')
call nc_check(nf90_put_att(ncid, VarID, 'units',     'dimensionless'), &
          'InitNetCDF', 'obsindex:units')

! Define the observation time 

call nc_check(nf90_def_var(ncid=ncid, name='time', xtype=nf90_double, &
          dimids=(/ ObsNumDimID /), varid=VarID), &
            'InitNetCDF', 'time:def_var')
call nc_check(nf90_put_att(ncid, VarID, 'long_name', 'time of observation'), &
          'InitNetCDF', 'time:long_name')
call nc_check(nf90_put_att(ncid, VarID, 'units',     'days since 1601-1-1'), &
          'InitNetCDF', 'time:units')
call nc_check(nf90_put_att(ncid, VarID, 'calendar', trim(calendarstring)), &
          'InitNetCDF', 'time:calendar')
call nc_check(nf90_put_att(ncid, VarID, 'valid_range', &
          (/ epoch_edges(1), epoch_edges(2) /)), &
          'InitNetCDF', 'time:valid_range')

! Define the observation type  (obs_type integer)

call nc_check(nf90_def_var(ncid=ncid, name='obs_type', xtype=nf90_int, &
          dimids=(/ ObsNumDimID /), varid=VarID), &
            'InitNetCDF', 'obs_type:def_var')
call nc_check(nf90_put_att(ncid, VarID, 'long_name', 'DART observation type'), &
          'InitNetCDF', 'obs_type:long_name')
call nc_check(nf90_put_att(ncid, VarID, 'explanation', 'see ObsTypesMetaData'), &
          'InitNetCDF', 'obs_type:explanation')

! Define the observation key  (index into linked list in original file)

call nc_check(nf90_def_var(ncid=ncid, name='obs_keys', xtype=nf90_int, &
          dimids=(/ ObsNumDimID /), varid=VarID), &
            'InitNetCDF', 'obs_keys:def_var')
call nc_check(nf90_put_att(ncid, VarID, 'long_name', 'DART key in linked list'), &
          'InitNetCDF', 'obs_keys:long_name')

! Define the observation copies

call nc_check(nf90_def_var(ncid=ncid, name='observations', xtype=nf90_double, &
          dimids=(/ ObsCopyDimID, ObsNumDimID /), varid=VarID), &
            'InitNetCDF', 'observations:def_var')
call nc_check(nf90_put_att(ncid, VarID,'long_name','org observation, estimates, etc.'), &
          'InitNetCDF', 'observations:long_name')
call nc_check(nf90_put_att(ncid, VarID, 'explanation', 'see CopyMetaData'), &
          'InitNetCDF', 'observations:explanation')
call nc_check(nf90_put_att(ncid, VarID, 'missing_value', NF90_FILL_DOUBLE), &
          'InitNetCDF', 'observations:missing_value')

! Define the QC copies

call nc_check(nf90_def_var(ncid=ncid, name='qc', xtype=nf90_int, &
          dimids=(/ QCCopyDimID, ObsNumDimID /), varid=VarID), &
            'InitNetCDF', 'qc:def_var')
call nc_check(nf90_put_att(ncid, VarID, 'long_name', 'QC values'), &
          'InitNetCDF', 'qc:long_name')
call nc_check(nf90_put_att(ncid, VarID, 'explanation', 'see QCMetaData'), &
          'InitNetCDF', 'qc:explanation')

! let the location module write what it needs to ...

call nc_write_location_atts(ncid, 0, ObsNumDimID)

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

call nc_check(nf90_inq_varid(ncid, 'copy', varid=VarID), &
           'InitNetCDF', 'inq_varid:copy '//trim(fname))
call nc_check(nf90_put_var(ncid, VarId, (/ (i,i=1,allNcopies) /) ), &
           'InitNetCDF', 'put_var:copy')

call nc_check(nf90_inq_varid(ncid, 'CopyMetaData', varid=VarID), &
           'InitNetCDF', 'inq_varid:CopyMetaData '//trim(fname))
call nc_check(nf90_put_var(ncid, VarID, module_obs_copy_names), &
           'InitNetCDF', 'put_var:CopyMetaData')

call nc_check(nf90_inq_varid(ncid, 'ObsTypes', varid=VarID), &
           'InitNetCDF', 'inq_varid:ObsTypes '//trim(fname))
call nc_check(nf90_put_var(ncid, VarId, (/ (i,i=1,num_obs_kinds) /) ), &
           'InitNetCDF', 'put_var:ObsTypes')

call nc_check(nf90_inq_varid(ncid, 'ObsTypesMetaData', varid=VarID), &
           'InitNetCDF', 'inq_varid:ObsTypesmetaData '//trim(fname))
call nc_check(nf90_put_var(ncid, VarID, my_obs_kind_names(1:num_obs_kinds)), &
           'InitNetCDF', 'put_var:ObsTypesMetaData')

call nc_check(nf90_inq_varid(ncid, 'qc_copy', varid=VarID), &
           'InitNetCDF', 'inq_varid:qc_copy '//trim(fname))
call nc_check(nf90_put_var(ncid, VarId, (/ (i,i=1,num_qc) /) ), &
           'InitNetCDF', 'put_var:qc_copy')

call nc_check(nf90_inq_varid(ncid, 'QCMetaData', varid=VarID), &
           'InitNetCDF', 'inq_varid:QCMetaData '//trim(fname))
call nc_check(nf90_put_var(ncid, VarID, module_qc_copy_names), &
           'InitNetCDF', 'put_var:QCMetaData')

!----------------------------------------------------------------------------
! Finish up ...
!----------------------------------------------------------------------------

call nc_check(nf90_sync( ncid), 'InitNetCDF', 'sync '//trim(fname))  

InitNetCDF = ncid

end Function InitNetCDF



Subroutine WriteNetCDF(ncid, fname, ngood, obs_copies, qc_copies, &
                       locations, obs_times, obs_types, obs_keys ) 
!============================================================================
integer,                             intent(in) :: ncid
character(len=*),                    intent(in) :: fname
integer,                             intent(in) :: ngood
real(r8),            dimension(:,:), intent(in) :: obs_copies
integer,             dimension(:,:), intent(in) :: qc_copies
type(location_type), dimension(:),   intent(in) :: locations
real(digits12),      dimension(:),   intent(in) :: obs_times
integer,             dimension(:),   intent(in) :: obs_types
integer,             dimension(:),   intent(in) :: obs_keys

integer :: DimID, dimlen, obsindex, iobs
integer, dimension(1) :: istart, icount, intval

integer :: obsldimlen, qcldimlen

integer :: ObsIndexVarID, TimeVarID, ObsTypeVarID, &
           ObsVarID, QCVarID, ObsKeyVarID

!----------------------------------------------------------------------------
! Find the current length of the unlimited dimension so we can add correctly.
!----------------------------------------------------------------------------

if (debug) write(*,*)'DEBUG --- entering WriteNetCDF'

obsldimlen = size(obs_copies,1)
 qcldimlen = size( qc_copies,1)

call nc_check(nf90_inquire(ncid, UnlimitedDimID=DimID), &
           'WriteNetCDF', 'inquire unlimited '//trim(fname))

call nc_check(nf90_inquire_dimension(ncid, DimID, len=dimlen), &
           'WriteNetCDF', 'inquire unlimited dimlen '//trim(fname))

obsindex  = dimlen + 1
istart(1) = obsindex
icount(1) = ngood

if (debug) write(*,*)'DEBUG --- WriteNetCDF istart/icount ',istart(1), icount(1)

call nc_check(nf90_inq_varid(ncid, 'ObsIndex', varid=ObsIndexVarID), &
          'WriteNetCDF', 'inq_varid:ObsIndex '//trim(fname))

call nc_check(nf90_inq_varid(ncid, 'time', varid=TimeVarID), &
          'WriteNetCDF', 'inq_varid:time '//trim(fname))

call nc_check(nf90_inq_varid(ncid, 'obs_type', varid=ObsTypeVarID), &
          'WriteNetCDF', 'inq_varid:obs_type '//trim(fname))

call nc_check(nf90_inq_varid(ncid, 'obs_keys', varid=ObsKeyVarID), &
          'WriteNetCDF', 'inq_varid:obs_keys '//trim(fname))

call nc_check(nf90_inq_varid(ncid, 'observations', varid=ObsVarID), &
          'WriteNetCDF', 'inq_varid:observations '//trim(fname))

call nc_check(nf90_inq_varid(ncid, 'qc', varid=QCVarID), &
          'WriteNetCDF', 'inq_varid:qc '//trim(fname))

WriteObs : do iobs = 1,ngood

   obsindex  = dimlen + iobs
   istart(1) = obsindex
   icount(1) = 1
   
   ! Fill the unlimited dimension coordinate variable 

   intval = obsindex 
   call nc_check(nf90_put_var(ncid, ObsIndexVarId, intval, &
                start=istart, count=icount), 'WriteNetCDF', 'put_var:ObsIndex')
   
   call nc_check(nf90_put_var(ncid, TimeVarId, (/ obs_times(iobs) /), &
                 start=istart, count=icount), 'WriteNetCDF', 'put_var:time')
   
   intval = obs_types(iobs) 
   call nc_check(nf90_put_var(ncid, ObsTypeVarId, intval, &
                 start=istart, count=icount), 'WriteNetCDF', 'put_var:obs_type')
   
   intval = obs_keys(iobs) 
   call nc_check(nf90_put_var(ncid, ObsKeyVarId, intval, &
                 start=istart, count=icount), 'WriteNetCDF', 'put_var:obs_keys')

   call nc_write_location(ncid, locations(iobs), obsindex, do_vert=.true.)

   call nc_check(nf90_put_var(ncid, ObsVarId, obs_copies(:,iobs), &
            start=(/ 1, obsindex /), count=(/ obsldimlen, 1 /) ), &
              'WriteNetCDF', 'put_var:observations')
   
   call nc_check(nf90_put_var(ncid, QCVarID,  qc_copies(:,iobs), &
            start=(/ 1, obsindex /), count=(/ qcldimlen, 1 /) ), &
              'WriteNetCDF', 'put_var:observations')

enddo WriteObs

!----------------------------------------------------------------------------
! finished ...
!----------------------------------------------------------------------------

call nc_check(nf90_sync( ncid), 'WriteNetCDF', 'sync '//trim(fname))  

if (debug) write(*,*)'DEBUG --- leaving WriteNetCDF'

end Subroutine WriteNetCDF



Subroutine CloseNetCDF(ncid, fname)
integer,          intent(in) :: ncid
character(len=*), intent(in) :: fname

if ( debug ) write(*,*)'DEBUG --- Closing ',trim(fname)

call nc_check(nf90_sync( ncid), 'CloseNetCDF', 'sync '//trim(fname))  
call nc_check(nf90_close(ncid), 'init_diag_output', 'close '//trim(fname))  

end Subroutine CloseNetCDF



function NC_Compatibility_Check(fname,ibin)
!----------------------------------------------------------------------------
! If the file exists, check for compatibility; i.e. it must have the same 
! unlimited dimension, CopyMetaData, and ObsTypesMetaData and return.
! If the file does not exist, we get to make one.
!----------------------------------------------------------------------------
character(len=*), intent(in) :: fname
integer,          intent(in) :: ibin
integer                      :: NC_Compatibility_Check

real(digits12)  :: epoch_edges(2), validRange(2)
integer         :: seconds, days
type(time_type) :: mytime

character(len=nf90_max_name) :: dimname
integer                      :: dimlen
integer, dimension(nf90_max_var_dims) :: dimIDs

integer :: ncid, StringDimID, ObsCopyDimID, DimID, VarID

! We will open it in write mode, but only read from it now.
! Later on it will need to be written to ...

call nc_check(nf90_open(path=trim(fname), mode=NF90_WRITE, ncid=ncid), &
        'NC_Compatibility_Check', 'open '//trim(fname))

! Check unlimited dimension variable dimension name

call nc_check(nf90_inquire(ncid, UnlimitedDimID=DimID), &
        'NC_Compatibility_Check', 'inquire unlimited '//trim(fname))
call nc_check(nf90_inquire_dimension(ncid, DimID, name=dimname, len=dimlen), &
        'NC_Compatibility_Check', 'inquire unlimited dimname '//trim(fname))

if ( trim(dimname) /= 'ObsIndex' ) then
   write(string1,*)'problem with '//trim(dimname)
   call error_handler(E_ERR,'NC_Compatibility_Check',string1,source)
endif

! Check time period - do not paste in observations from wrong epoch 
! The current epoch times can fit _inside_ ... no problem.

call get_time_from_schedule(mytime,schedule,ibin,1)
call get_time(mytime,seconds,days)
epoch_edges(1) = days + seconds/86400.0_digits12

call get_time_from_schedule(mytime,schedule,ibin,2)
call get_time(mytime,seconds,days)
epoch_edges(2) = days + seconds/86400.0_digits12

call nc_check(nf90_inq_varid(ncid, 'time', varid=VarID), &
        'NC_Compatibility_Check', 'inq_varid:time '//trim(fname))
call nc_check(nf90_get_att(ncid, VarID, 'valid_range', validRange), &
        'NC_Compatibility_Check', 'get_att:timerange '//trim(fname))

if ( epoch_edges(1) < validRange(1) ) then
   write(string1,*)'problem : obs_time ',epoch_edges(1),' < ',validRange(1),' netcdf time'
   call error_handler(E_ERR,'NC_Compatibility_Check',string1,source)
endif

if ( epoch_edges(2) > validRange(2) ) then
   write(string1,*)'problem : obs_time ',epoch_edges(2),' > ',validRange(1),' netcdf time'
   call error_handler(E_ERR,'NC_Compatibility_Check',string1,source)
endif

! Check stringlength dimension

call nc_check(nf90_inq_dimid(ncid, 'stringlength', dimid=StringDimID), &
        'NC_Compatibility_Check', 'inquire stringdimid '//trim(fname))
call nc_check(nf90_inquire_dimension(ncid, StringDimID, name=dimname, len=dimlen), &
        'NC_Compatibility_Check', 'inquire_dimension stringdimid '//trim(fname))

if ( dimlen /= stringlength ) then
   write(string1,*)'stringlength problem ... ',dimlen,' /= ',stringlength
   call error_handler(E_ERR,'NC_Compatibility_Check',string1,source)
endif

! Check the number of copies

call nc_check(nf90_inq_dimid(ncid, 'copy', dimid=ObsCopyDimID), &
        'NC_Compatibility_Check', 'inq_dimid:copy  '//trim(fname))
call nc_check(nf90_inquire_dimension(ncid, ObsCopyDimID, name=dimname, len=dimlen), &
        'NC_Compatibility_Check', 'inquire_dimension:copy '//trim(fname))

if ( dimlen /= allNcopies ) then
   write(string1,*)'different number of copies ... ',dimlen,' /= ',allNcopies
   call error_handler(E_ERR,'NC_Compatibility_Check',string1,source)
endif

! must check shape and actual values of copy metadata

call nc_check(nf90_inq_varid(ncid, 'CopyMetaData', varid=VarID), &
        'NC_Compatibility_Check', 'inq_varid CopyMetaData '//trim(fname))

call nc_check(nf90_inquire_variable(ncid, VarID, dimids=dimIDs), &
        'NC_Compatibility_Check', 'inquire_variable CopyMetaData '//trim(fname))

call nc_check(nf90_inquire_dimension(ncid, dimIDs(1), len=dimlen), &
        'NC_Compatibility_Check', 'inquire_dimension CopyMetaData(1) '//trim(fname))

if ( dimlen /= stringlength ) then
   write(string1,*)'copymetadata dim1 ',dimlen,' /= ',stringlength
   call error_handler(E_ERR,'NC_Compatibility_Check',string1,source)
endif

call nc_check(nf90_inquire_dimension(ncid, dimIDs(2), len=dimlen), &
        'NC_Compatibility_Check', 'inquire_dimension CopyMetaData(2) '//trim(fname))

if ( dimlen /= allNcopies ) then
   write(string1,*)'copymetadata dim2 ',dimlen,' /= ',allNcopies
   call error_handler(E_ERR,'NC_Compatibility_Check',string1,source)
endif

do i = 1,allNcopies
   call nc_check(nf90_get_var(ncid, VarID, values=dimname, &
                    start = (/ 1, i /), count = (/ stringlength, 1 /)), &
        'NC_Compatibility_Check', 'get_var CopyMetaData '//trim(fname))

   if ( trim(dimname) /= trim(module_obs_copy_names(i)) ) then
      write(string1,*)'copymetadata(',i,') ',trim(dimname),' /= ',trim(module_obs_copy_names(i))
      call error_handler(E_ERR,'NC_Compatibility_Check',string1,source)
   endif
enddo

! must check shape and actual values of QC metadata

call nc_check(nf90_inq_varid(ncid, 'QCMetaData', varid=VarID), &
        'NC_Compatibility_Check', 'inq_varid QCMetaData '//trim(fname))

call nc_check(nf90_inquire_variable(ncid, VarID, dimids=dimIDs), &
        'NC_Compatibility_Check', 'inquire_variable QCMetaData '//trim(fname))

call nc_check(nf90_inquire_dimension(ncid, dimIDs(1), len=dimlen), &
        'NC_Compatibility_Check', 'inquire_dimension QCMetaData(1) '//trim(fname))

if ( dimlen /= stringlength ) then
   write(string1,*)'QCMetaData dim1 ',dimlen,' /= ',stringlength
   call error_handler(E_ERR,'NC_Compatibility_Check',string1,source)
endif

call nc_check(nf90_inquire_dimension(ncid, dimIDs(2), len=dimlen), &
        'NC_Compatibility_Check', 'inquire_dimension QCMetaData(2) '//trim(fname))

if ( dimlen /= num_qc ) then
   write(string1,*)'QCMetaData dim2 ',dimlen,' /= ',num_qc
   call error_handler(E_ERR,'NC_Compatibility_Check',string1,source)
endif

do i = 1,num_qc
   call nc_check(nf90_get_var(ncid, VarID, dimname, &
                    start = (/ 1, i /), count = (/ stringlength, 1 /)), &
        'NC_Compatibility_Check', 'get_var QCMetaData '//trim(fname))

   if ( trim(dimname) /= trim(module_qc_copy_names(i)) ) then
      write(string1,*)'QCMetaData ',i,trim(dimname),' /= ',trim(module_qc_copy_names(i))
      call error_handler(E_ERR,'NC_Compatibility_Check',string1,source)
   endif
enddo

! must check shape and actual values of observation types

call nc_check(nf90_inq_varid(ncid, 'ObsTypesMetaData', varid=VarID), &
        'NC_Compatibility_Check', 'inq_varid ObsTypesMetaData '//trim(fname))

call nc_check(nf90_inquire_variable(ncid, VarID, dimids=dimIDs), &
        'NC_Compatibility_Check', 'inquire_variable ObsTypesMetaData '//trim(fname))

call nc_check(nf90_inquire_dimension(ncid, dimIDs(1), len=dimlen), &
        'NC_Compatibility_Check', 'inquire_dimension ObsTypesMetaData(1) '//trim(fname))

if ( dimlen /= stringlength ) then
   write(string1,*)'ObsTypesMetaData dim1 ',dimlen,' /= ',stringlength
   call error_handler(E_ERR,'NC_Compatibility_Check',string1,source)
endif

call nc_check(nf90_inquire_dimension(ncid, dimIDs(2), len=dimlen), &
        'NC_Compatibility_Check', 'inquire_dimension ObsTypesMetaData(2) '//trim(fname))

if ( dimlen /= num_obs_kinds ) then
   write(string1,*)'copymetadata dim2 ',dimlen,' /= ',num_obs_kinds
   call error_handler(E_ERR,'NC_Compatibility_Check',string1,source)
endif

do i = 1,num_obs_kinds
   call nc_check(nf90_get_var(ncid, VarID, dimname, &
                    start = (/ 1, i /), count = (/ stringlength, 1 /)), &
        'NC_Compatibility_Check', 'get_var ObsTypesMetaData '//trim(fname))
   if ( trim(dimname) /= trim(my_obs_kind_names(i)) ) then
      write(string1,*)'typesmetavrid ',i,trim(dimname),' /= ',trim(my_obs_kind_names(i))
      call error_handler(E_ERR,'NC_Compatibility_Check',string1,source)
   endif
enddo

NC_Compatibility_Check = ncid

end function NC_Compatibility_Check



end program obs_seq_to_netcdf

