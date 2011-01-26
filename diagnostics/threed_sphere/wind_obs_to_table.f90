! DART software - Copyright © 2004 - 2010 UCAR. This open source software is
! provided by UCAR, "as is", without charge, subject to all terms of use at
! http://www.image.ucar.edu/DAReS/DART/DART_download

program wind_obs_to_table

! <next few lines under version control, do not edit>
! $URL$
! $Id$
! $Revision$
! $Date$

!-----------------------------------------------------------------------
! The programs defines a series of epochs (periods of time) 
!
! All 'possible' obs_kinds are treated separately.
!-----------------------------------------------------------------------

use        types_mod, only : r4, r8, digits12, MISSING_R8, MISSING_R4
use obs_sequence_mod, only : read_obs_seq, obs_type, obs_sequence_type, get_first_obs, &
                             get_obs_from_key, get_obs_def, get_copy_meta_data, &
                             get_obs_time_range, get_time_range_keys, get_num_obs, &
                             get_next_obs, get_num_times, get_obs_values, init_obs, &
                             assignment(=), get_num_copies, static_init_obs_sequence, &
                             get_qc, destroy_obs_sequence, read_obs_seq_header, & 
                             get_last_obs, destroy_obs, get_num_qc, get_qc_meta_data
use      obs_def_mod, only : obs_def_type, get_obs_def_error_variance, get_obs_def_time, &
                             get_obs_def_location,  get_obs_kind, get_obs_name
use     obs_kind_mod, only : max_obs_kinds, get_obs_kind_var_type, get_obs_kind_name, &
                             do_obs_form_pair, add_wind_names, &
                             KIND_U_WIND_COMPONENT, KIND_V_WIND_COMPONENT
use     location_mod, only : location_type, get_location, set_location_missing, &
                             write_location, operator(/=), operator(==), &
                             set_location, is_location_in_region, VERTISUNDEF
use time_manager_mod, only : time_type, set_date, set_time, get_time, print_time, &
                             set_calendar_type, get_calendar_string, print_date, &
                             operator(*), operator(+), operator(-), &
                             operator(>), operator(<), operator(/), &
                             operator(/=), operator(<=)
use     schedule_mod, only : schedule_type, set_regular_schedule, get_schedule_length, &
                             get_time_from_schedule
use    utilities_mod, only : open_file, close_file, register_module, &
                             file_exist, error_handler, E_ERR, E_WARN, E_MSG, &
                             initialize_utilities, nmlfileunit, timestamp, &
                             find_namelist_in_file, check_namelist_read, nc_check, &
                             next_file, find_textfile_dims, file_to_text

use typeSizes
use netcdf

implicit none

! version controlled file description for error handling, do not edit
character(len=128), parameter :: &
   source   = '$URL$', &
   revision = '$Revision$', &
   revdate  = '$Date$'

!---------------------------------------------------------------------
!---------------------------------------------------------------------

integer, parameter :: stringlength = 32
logical, parameter :: DEBUG = .false.

!---------------------------------------------------------------------
! variables associated with the observation
!---------------------------------------------------------------------

type(obs_sequence_type) :: seq
type(obs_type)          :: observation, next_obs
type(obs_type)          :: obs1, obsN
type(obs_def_type)      :: obs_def
type(location_type)     :: obs_loc, minl, maxl

character(len = 129) :: obs_seq_in_file_name
character(len = 129), allocatable, dimension(:) :: obs_seq_filenames

! We are treating winds as a vector pair, but we are handling the
! observations serially. Consequently, we exploit the fact that
! the U observations are _followed_ by the V observations.

real(r8)            :: U_obs         = 0.0_r8
real(r8)            :: U_obs_err_var = 0.0_r8
type(location_type) :: U_obs_loc
integer             :: U_flavor
integer             :: U_type        = KIND_V_WIND_COMPONENT ! intentional mismatch
integer             :: U_qc          = 0

integer :: obs_index
integer :: flavor, wflavor ! THIS IS THE (global) 'KIND' in the obs_def_mod list. 
integer :: num_copies, num_qc, num_obs, max_num_obs, obs_seq_file_id
integer :: num_obs_kinds
character(len=129) :: obs_seq_read_format
logical :: pre_I_format

integer,  dimension(2) :: key_bounds
real(r8), dimension(1) :: obs
real(r8) :: obs_err_var

integer,  allocatable, dimension(:) :: keys

logical :: out_of_range, is_there_one, keeper

!---------------------------------------------------------------------
! variables associated with quality control
!
! qc_index  reflects the 'original' QC value of the observation, if any.
!           Most frequently represents the value NCEP assigned to their
!           observations.
!
! dart_qc_index 
! 0     observation assimilated
! 1     observation evaluated only
!   --- everything above this means the prior and posterior are OK
! 2     assimilated, but the posterior forward operator failed
! 3     Evaluated only, but the posterior forward operator failed
!   --- everything above this means only the prior is OK
! 4     prior forward operator failed
! 5     not used
! 6     prior QC rejected
! 7     outlier rejected
! 8+    reserved for future use

integer             :: qc_index, dart_qc_index
integer             :: qc_integer
integer, parameter  :: QC_MAX = 7
integer, parameter  :: QC_MAX_PRIOR     = 3
integer, parameter  :: QC_MAX_POSTERIOR = 1
real(r8), allocatable, dimension(:) :: qc
real(r8), allocatable, dimension(:) :: copyvals

!-----------------------------------------------------------------------
! Namelist with (some scalar) default values
!-----------------------------------------------------------------------

character(len = 129) :: obs_sequence_name = 'obs_seq.final'

real(r8) :: lonlim1= MISSING_R8, lonlim2= MISSING_R8
real(r8) :: latlim1= MISSING_R8, latlim2= MISSING_R8 

logical :: verbose = .false.

namelist /wind_obs_to_table_nml/ obs_sequence_name, lonlim1, lonlim2, &
                            latlim1, latlim2, verbose

!-----------------------------------------------------------------------
! Variables used to accumulate the statistics.
!-----------------------------------------------------------------------

integer, parameter :: Ncopies = 10
character(len = stringlength), dimension(Ncopies) :: copy_names = &
   (/ 'Nposs      ', 'Nused      ', 'NbadQC     ', 'NbadIZ     ', 'NbadUV     ', &
      'NbadLV     ', 'rmse       ', 'bias       ', 'spread     ', 'totalspread' /)

!-----------------------------------------------------------------------
! General purpose variables
!-----------------------------------------------------------------------

integer  :: iepoch, ivar, ifile, num_obs_in_epoch
real(r8) :: obslon, obslat, obslevel, obsloc3(3)

integer  :: obsindex, iunit, io
integer  :: Nepochs

character(len = stringlength), pointer, dimension(:) :: my_obs_kind_names

type(schedule_type) :: schedule
type(time_type) :: TimeMin, TimeMax    ! of the entire period of interest
type(time_type) :: beg_time, end_time  ! of the particular bin
type(time_type) :: seqT1, seqTN        ! first,last time in entire observation sequence
type(time_type) :: obs_time

character(len = 129) :: ncName, windName, msgstring, calendarstring

!=======================================================================
! Get the party started
!=======================================================================

call initialize_utilities('wind_obs_to_table')
call register_module(source,revision,revdate) 
call static_init_obs_sequence()  ! Initialize the obs sequence module 

!----------------------------------------------------------------------
! Define/Append the 'horizontal wind' obs_kinds to supplant the list declared
! in obs_kind_mod.f90 i.e. if there is a RADIOSONDE_U_WIND_COMPONENT
! and a RADIOSONDE_V_WIND_COMPONENT, there must be a RADIOSONDE_HORIZONTAL_WIND
! Replace calls to 'get_obs_kind_name' with variable 'my_obs_kind_names'
!----------------------------------------------------------------------

num_obs_kinds = add_wind_names(my_obs_kind_names)

!----------------------------------------------------------------------
! Read the namelist
!----------------------------------------------------------------------

call find_namelist_in_file('input.nml', 'wind_obs_to_table_nml', iunit)
read(iunit, nml = wind_obs_to_table_nml, iostat = io)
call check_namelist_read(iunit, io, 'wind_obs_to_table_nml')

! Record the namelist values used for the run ...
write(nmlfileunit, nml=wind_obs_to_table_nml)
write(    *      , nml=wind_obs_to_table_nml)

!----------------------------------------------------------------------
! SetSchedule rectifies user input and the final binning sequence.
!----------------------------------------------------------------------

call set_regular_schedule(schedule) ! also sets calendar type

Nepochs = get_schedule_length(schedule)
call get_time_from_schedule(TimeMin, schedule,       1, 1)
call get_time_from_schedule(TimeMax, schedule, Nepochs, 2)
call get_calendar_string(calendarstring)

U_obs_loc = set_location_missing()
minl = set_location(lonlim1, latlim1, 0.0_r8, VERTISUNDEF) ! vertical unimportant
maxl = set_location(lonlim2, latlim2, 0.0_r8, VERTISUNDEF) ! vertical unimportant

!----------------------------------------------------------------------
! Prepare the variables
!----------------------------------------------------------------------

allocate(obs_seq_filenames(Nepochs*4))
obs_seq_filenames = 'null'

ObsFileLoop : do ifile=1, Nepochs*4
!-----------------------------------------------------------------------

   obs_seq_in_file_name = next_file(obs_sequence_name,ifile)

   if ( file_exist(trim(obs_seq_in_file_name)) ) then
      write(msgstring,*)'opening ', trim(obs_seq_in_file_name)
      call error_handler(E_MSG,'wind_obs_to_table',msgstring,source,revision,revdate)
   else
      write(msgstring,*)trim(obs_seq_in_file_name),&
                        ' does not exist. Finishing up.'
      call error_handler(E_MSG,'wind_obs_to_table',msgstring,source,revision,revdate)
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

   if (num_qc     > 0) allocate( qc(num_qc) )
   if (num_copies > 0) allocate( copyvals(num_copies) )

   if ( DEBUG ) then
      write(*,*)
      write(*,*)'num_copies          is ',num_copies
      write(*,*)'num_qc              is ',num_qc
      write(*,*)'num_obs             is ',num_obs
      write(*,*)'max_num_obs         is ',max_num_obs
      write(*,*)'obs_seq_read_format is ',trim(obs_seq_read_format)
      write(*,*)'pre_I_format        is ',pre_I_format
      write(*,*)
   endif

   ! Read in the entire observation sequence

   call read_obs_seq(obs_seq_in_file_name, 0, 0, 0, seq)

   !--------------------------------------------------------------------
   ! Determine the time encompassed in the observation sequence.
   !--------------------------------------------------------------------

   is_there_one = get_first_obs(seq, obs1)
   if ( .not. is_there_one ) then
      call error_handler(E_ERR,'wind_obs_to_table','No first observation  in sequence.', &
      source,revision,revdate)
   endif
   call get_obs_def(obs1,   obs_def)
   seqT1 = get_obs_def_time(obs_def)

   is_there_one = get_last_obs(seq, obsN)
   if ( .not. is_there_one ) then
      call error_handler(E_ERR,'wind_obs_to_table','No last observation in sequence.', &
      source,revision,revdate)
   endif
   call get_obs_def(obsN,   obs_def)
   seqTN = get_obs_def_time(obs_def)

   if ( verbose ) then
      call print_time(seqT1,'First observation time')
      call print_time(seqTN,'Last  observation time')
      call print_time(TimeMin,'TimeMin from input')
      call print_time(TimeMax,'TimeMax from input')
      call print_date(seqT1,'First observation date')
      call print_date(seqTN,'Last  observation date')
      call print_date(TimeMin,'DateMin from input')
      call print_date(TimeMax,'DateMax from input')
   endif

   !--------------------------------------------------------------------
   ! If the last observation is before the period of interest, move on.
   !--------------------------------------------------------------------

   if ( seqTN < TimeMin ) then
      if (verbose) write(*,*)'seqTN < TimeMin ... trying next file.'
      call destroy_obs(obs1)
      call destroy_obs(obsN)
      call destroy_obs(observation)
      call destroy_obs(next_obs)
      call destroy_obs_sequence(seq)
      if (allocated(qc)) deallocate( qc )
      if (allocated(copyvals)) deallocate( copyvals )
      cycle ObsFileLoop
   else
      if (verbose) write(*,*)'seqTN > TimeMin ... using ', trim(obs_seq_in_file_name)
   endif

   !--------------------------------------------------------------------
   ! If the first observation is after the period of interest, finish.
   !--------------------------------------------------------------------

   if ( seqT1 > TimeMax ) then
      if (verbose) write(*,*)'seqT1 > TimeMax ... stopping.'
      call destroy_obs(obs1)
      call destroy_obs(obsN)
      call destroy_obs(observation)
      call destroy_obs(next_obs)
      call destroy_obs_sequence(seq)
      if (allocated(qc)) deallocate( qc )
      if (allocated(copyvals)) deallocate( copyvals )
      exit ObsFileLoop
   else
      if (verbose) write(*,*)'seqT1 < TimeMax ... using ',trim(obs_seq_in_file_name)
   endif

   !--------------------------------------------------------------------
   ! Find the index of obs, ensemble mean, spread ... etc.
   !--------------------------------------------------------------------
   ! Only require obs_index to be present; this allows the program
   ! to be run on obs_seq.in files which have no means or spreads.
   ! You can still plot locations, but that's it.
   !--------------------------------------------------------------------

   call SetIndices( obs_index, qc_index, dart_qc_index )

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

      ! Append epoch number to output file name
      write(windName,'(a,i3.3,a)') 'wind_vectors.', iepoch, '.dat'
      iunit = open_file(windName, form='formatted', action='rewind')

      allocate(keys(num_obs_in_epoch))

      call get_time_range_keys(seq, key_bounds, num_obs_in_epoch, keys)

      !-----------------------------------------------------------------
      ObservationLoop : do obsindex = 1, num_obs_in_epoch
      !-----------------------------------------------------------------

         ! 'flavor' is from the 'master list' in the obs_kind_mod.f90
         ! each obs_seq.final file has their own private kind - which
         ! gets mapped to the 'master list', if you will.

         call get_obs_from_key(seq, keys(obsindex), observation)
         call get_obs_values(observation, obs, obs_index)
         call get_obs_def(observation, obs_def)
         call get_qc(observation, qc)

         flavor      = get_obs_kind(obs_def)
         obs_err_var = get_obs_def_error_variance(obs_def)
         obs_time    = get_obs_def_time(obs_def)
         obs_loc     = get_obs_def_location(obs_def)
         obsloc3     = get_location(obs_loc)

         obslon   = obsloc3(1) ! [  0, 360]
         obslat   = obsloc3(2) ! [-90,  90]
         obslevel = obsloc3(3) ! variable-dependent

         !--------------------------------------------------------------
         ! We have one Region of interest
         !--------------------------------------------------------------

         keeper = is_location_in_region( obs_loc, minl, maxl ) 

         if ( .not. keeper ) cycle ObservationLoop

         !--------------------------------------------------------------
         ! Convert the DART QC data
         !--------------------------------------------------------------

         if (dart_qc_index > 0) then
            qc_integer = min( nint(qc(dart_qc_index)), QC_MAX )
         else if (qc_index > 0) then ! If there is no dart_qc in obs_seq ...
            qc_integer = nint(qc(qc_index))
         else                        ! If there is no original QC value ...
            qc_integer = 999
         endif

         !--------------------------------------------------------------
         ! Summary of observation knowledge at this point
         !--------------------------------------------------------------

         if ( DEBUG ) then
            write(*,*)'observation # ',obsindex
            write(*,*)'key           ',keys(obsindex)
            write(*,*)'obs_flavor    ',flavor
            write(*,*)'obs_err_var   ',obs_err_var
            write(*,*)'lon/lat/level ',obslon,obslat,obslevel
            write(*,*)'obs(1)        ',obs(1)
            write(*,*)'qc            ',qc
         endif

         !--------------------------------------------------------------
         ! If it is a U wind component, all we need to do is save it.
         ! It will be matched up with the subsequent V component.
         ! At some point we have to remove the dependency that the 
         ! U component MUST preceed the V component.
         !--------------------------------------------------------------

         if ( get_obs_kind_var_type(flavor) == KIND_U_WIND_COMPONENT ) then

            U_obs         = obs(1)
            U_obs_err_var = obs_err_var
            U_obs_loc     = obs_loc
            U_flavor      = flavor
            U_type        = KIND_U_WIND_COMPONENT
            U_qc          = qc_integer

            cycle ObservationLoop

         endif


         !-----------------------------------------------------------
         ! Additional work for horizontal wind (given U,V)
         !-----------------------------------------------------------

         ObsIsWindCheck: if ( get_obs_kind_var_type(flavor) == KIND_V_WIND_COMPONENT ) then

            ! The big assumption is that the U wind component has
            ! immediately preceeded the V component and has been saved.
            ! 
            ! We check for observation compatibility and the index for 
            ! this wind 'kind' ... not originally in the max_obs_kind namelist.
            ! this will be the 'wflavor' (wind) flavor.

            if ((obs_loc == U_obs_loc) .and.   &
               do_obs_form_pair(flavor,U_flavor,keys(obsindex),my_obs_kind_names,wflavor)) then

               call WritePairs( iunit, U_flavor, flavor, &
                                obs_time, obslon, obslat, obslevel, &
                                U_qc, qc_integer, U_obs, obs(1))
            else
               write(*,*)'time series : V with no U obs index ', keys(obsindex)
               wflavor = -99
            endif

         endif ObsIsWindCheck

      !-----------------------------------------------------------------
      enddo ObservationLoop
      !-----------------------------------------------------------------

      deallocate(keys)

      close(iunit)

   enddo EpochLoop

   if (verbose) write(*,*)'End of EpochLoop for ',trim(obs_seq_in_file_name)

   call destroy_obs(obs1)
   call destroy_obs(obsN)
   call destroy_obs(observation)
   call destroy_obs(next_obs)
   call destroy_obs_sequence(seq)
   if (allocated(qc)) deallocate( qc )
   if (allocated(copyvals)) deallocate( copyvals )

enddo ObsFileLoop

!----------------------------------------------------------------------
! Open netCDF output file 
!----------------------------------------------------------------------

ncName = 'wind_obs_to_table_output.nc'
  
call WriteNetCDF(ncName)

!-----------------------------------------------------------------------
! Really, really, done.
!-----------------------------------------------------------------------

deallocate(obs_seq_filenames)


deallocate(my_obs_kind_names )

call timestamp(source,revision,revdate,'end') ! That closes the log file, too.

!======================================================================
CONTAINS
!======================================================================
! These routines use common variables from the scope of this file.
! If it's not in the argument list ... it's scoped within this file.
!======================================================================

   Subroutine WriteNetCDF(fname)
   character(len=129), intent(in) :: fname

   integer :: ncid, i, indx1, nlines, linelen
   integer ::    TimeDimID
   integer ::    CopyDimID,    CopyVarID,  CopyMetaVarID
   integer ::   TypesDimID,   TypesVarID, TypesMetaVarID
   integer ::  BoundsDimID,  BoundsVarID  
   integer ::  StringDimID
   integer ::  linelenDimID, nlinesDimID, nmlVarID

   integer :: TimeBoundsVarID

   character(len=8)      :: crdate      ! needed by F90 DATE_AND_TIME intrinsic
   character(len=10)     :: crtime      ! needed by F90 DATE_AND_TIME intrinsic
   character(len=5)      :: crzone      ! needed by F90 DATE_AND_TIME intrinsic
   integer, dimension(8) :: values      ! needed by F90 DATE_AND_TIME intrinsic

   character(len=129), allocatable, dimension(:) :: textblock

   real(digits12)  :: epoch_edges(2,Nepochs)
   integer         :: seconds, days
   type(time_type) :: mytime

   if(.not. byteSizesOK()) then
       call error_handler(E_ERR,'WriteNetCDF', &
      'Compiler does not support required kinds of variables.',source,revision,revdate)
   endif

   call nc_check(nf90_create(path = trim(fname), cmode = nf90_share, &
            ncid = ncid), 'wind_obs_to_table:WriteNetCDF', 'create '//trim(fname))

   write(msgstring,*)trim(ncName), ' is fortran unit ',ncid
   call error_handler(E_MSG,'WriteNetCDF',msgstring,source,revision,revdate)

   do i = 1,Nepochs
      call get_time_from_schedule(mytime,schedule,i,1)
      call get_time(mytime,seconds,days)
      epoch_edges(1,i) = days + seconds/86400.0_digits12

      call get_time_from_schedule(mytime,schedule,i,2)
      call get_time(mytime,seconds,days)
      epoch_edges(2,i) = days + seconds/86400.0_digits12
   enddo

   !----------------------------------------------------------------------------
   ! Write Global Attributes 
   !----------------------------------------------------------------------------

   call DATE_AND_TIME(crdate,crtime,crzone,values)
   write(msgstring,'(''YYYY MM DD HH MM SS = '',i4,5(1x,i2.2))') &
                  values(1), values(2), values(3), values(5), values(6), values(7)
   call nc_check(nf90_put_att(ncid, NF90_GLOBAL, 'creation_date', trim(msgstring) ), &
              'WriteNetCDF', 'put_att creation_date '//trim(fname))

   call nc_check(nf90_put_att(ncid, NF90_GLOBAL, 'wind_obs_to_table_source', source ), &
              'WriteNetCDF', 'put_att wind_obs_to_table_source '//trim(fname))
   call nc_check(nf90_put_att(ncid, NF90_GLOBAL, 'wind_obs_to_table_revision', revision ), &
              'WriteNetCDF', 'put_att wind_obs_to_table_revision '//trim(fname))
   call nc_check(nf90_put_att(ncid, NF90_GLOBAL, 'wind_obs_to_table_revdate', revdate ), &
              'WriteNetCDF', 'put_att wind_obs_to_table_revdate '//trim(fname))

   call nc_check(nf90_put_att(ncid, NF90_GLOBAL, 'horizontal_wind', &
              'vector wind derived from U,V components' ), &
              'WriteNetCDF', 'put_att wind '//trim(fname))


   ! write all observation sequence files used
   FILEloop : do i = 1,SIZE(obs_seq_filenames)

     indx1 = index(obs_seq_filenames(i),'null')

     if (indx1 > 0) exit FILEloop

     write(msgstring,'(''obs_seq_file_'',i3.3)')i
     call nc_check(nf90_put_att(ncid, NF90_GLOBAL, &
            trim(msgstring), trim(obs_seq_filenames(i)) ), &
            'WriteNetCDF', 'region_names:obs_kinds')

   enddo FILEloop

   ! write all 'known' observation types
   call nc_check(nf90_put_att(ncid, NF90_GLOBAL, 'comment', &
              'All known observation types follow. &
              &Also see ObservationTypes variable.' ), &
              'WriteNetCDF', 'put_att latlim2 '//trim(fname))
   do ivar = 1,max_obs_kinds
     call nc_check(nf90_put_att(ncid, NF90_GLOBAL, &
            trim(my_obs_kind_names(ivar)), ivar ), &
            'WriteNetCDF', 'region_names:obs_kinds')
   enddo

   !----------------------------------------------------------------------------
   ! Define the dimensions
   !----------------------------------------------------------------------------

   ! write all namelist quantities 
   call find_textfile_dims('input.nml', nlines, linelen)
   allocate(textblock(nlines))
   textblock = ''

   call nc_check(nf90_def_dim(ncid=ncid, &
              name="linelen", len = len(textblock(1)), dimid = linelenDimID), &
              'WriteNetCDF', 'linelen:def_dim '//'input.nml')

   call nc_check(nf90_def_dim(ncid=ncid, &
              name="nlines", len = nlines,           dimid = nlinesDimID), &
              'WriteNetCDF', 'nlines:def_dim '//'input.nml')

   call nc_check(nf90_def_dim(ncid=ncid, &
              name='copy', len = Ncopies,            dimid = CopyDimID), &
              'WriteNetCDF', 'copy:def_dim '//trim(fname))

   call nc_check(nf90_def_dim(ncid=ncid, &
              name='obstypes', len = max_obs_kinds,  dimid = TypesDimID), &
              'WriteNetCDF', 'types:def_dim '//trim(fname))

   call nc_check(nf90_def_dim(ncid=ncid, &
              name='time',   len = NF90_UNLIMITED,   dimid = TimeDimID), &
              'WriteNetCDF', 'time:def_dim '//trim(fname))
   call nc_check(nf90_def_dim(ncid=ncid, &
              name='bounds',   len = 2,  dimid = BoundsDimID), &
              'WriteNetCDF', 'bounds:def_dim '//trim(fname))

   call nc_check(nf90_def_dim(ncid=ncid, &
              name='stringlength', len = stringlength, dimid = StringDimID), &
              'WriteNetCDF', 'stringlength:def_dim '//trim(fname))

   !----------------------------------------------------------------------------
   ! Define the variables
   !----------------------------------------------------------------------------

   call nc_check(nf90_def_var(ncid=ncid, name="namelist", xtype=nf90_char, &
             dimids = (/ linelenDimID, nlinesDimID /), varid=nmlVarID), &
             'WriteNetCDF', 'namelist:def_var')
   call nc_check(nf90_put_att(ncid, nmlVarID, "long_name", "input.nml contents"), &
             'WriteNetCDF', 'namelist:long_name')

   ! Define the types of derived quantities - aka - 'copies'

   call nc_check(nf90_def_var(ncid=ncid, name='copy', xtype=nf90_int, &
             dimids=CopyDimID, varid=CopyVarID), &
             'WriteNetCDF', 'copy:def_var')
   call nc_check(nf90_put_att(ncid, CopyVarID, 'explanation', 'see CopyMetaData'), &
             'WriteNetCDF', 'copy:explanation')

   ! Define the observation types - needed to be a coordinate variable

   call nc_check(nf90_def_var(ncid=ncid, name='obstypes', xtype=nf90_int, &
             dimids=TypesDimID, varid=TypesVarID), &
             'WriteNetCDF', 'types:def_var')
   call nc_check(nf90_put_att(ncid, TypesVarID, 'explanation', 'see ObservationTypes'), &
             'WriteNetCDF', 'types:explanation')

   ! Define 'bounds' dimension

   call nc_check(nf90_def_var(ncid=ncid, name='bounds', xtype=nf90_int, &
             dimids=BoundsDimID, varid=BoundsVarID), 'WriteNetCDF', 'bounds:def_var')
   call nc_check(nf90_put_att(ncid, BoundsVarID, 'valid_range', (/1,2/)), &
             'WriteNetCDF', 'bounds:valid_range')

   ! Define the time edges coordinate variable and attributes

   call nc_check(nf90_def_var(ncid=ncid, name='time_bounds', xtype=nf90_double, &
             dimids=(/ BoundsDimID, TimeDimID/), varid=TimeBoundsVarID), &
               'WriteNetCDF', 'time_bounds:def_var')
   call nc_check(nf90_put_att(ncid, TimeBoundsVarID, 'long_name', 'temporal bin edges'), &
             'WriteNetCDF', 'time_bounds:long_name')
   call nc_check(nf90_put_att(ncid, TimeBoundsVarID, 'units',     'days since 1601-1-1'), &
             'WriteNetCDF', 'time_bounds:units')
   call nc_check(nf90_put_att(ncid, TimeBoundsVarID, 'calendar', trim(calendarstring)), &
             'WriteNetCDF', 'time_bounds:calendar')

   call nc_check(nf90_put_att(ncid, TimeBoundsVarID, 'valid_range', &
             (/ epoch_edges(1,1), epoch_edges(2,Nepochs) /)), &
             'WriteNetCDF', 'time_bounds:valid_range')

   call nc_check(nf90_def_var(ncid=ncid, name='CopyMetaData', xtype=nf90_char, &
             dimids=(/ StringDimID, CopyDimID /), varid=CopyMetaVarID), &
             'WriteNetCDF', 'copymeta:def_var')
   call nc_check(nf90_put_att(ncid, CopyMetaVarID, 'long_name', 'quantity names'), &
             'WriteNetCDF', 'copymeta:long_name')

   call nc_check(nf90_def_var(ncid=ncid, name='ObservationTypes', xtype=nf90_char, &
             dimids=(/ StringDimID, TypesDimID /), varid=TypesMetaVarID), &
             'WriteNetCDF', 'typesmeta:def_var')
   call nc_check(nf90_put_att(ncid, TypesMetaVarID, 'long_name', 'DART observation types'), &
             'WriteNetCDF', 'typesmeta:long_name')
   call nc_check(nf90_put_att(ncid, TypesMetaVarID, 'comment', &
         'table relating integer to observation type string'), &
             'WriteNetCDF', 'typesmeta:comment')

   ! Set nofill mode - supposed to be performance gain
 
   call nc_check(nf90_set_fill(ncid, NF90_NOFILL, i),  &
            'wind_obs_to_table:WriteNetCDF', 'set_nofill '//trim(fname))

   !----------------------------------------------------------------------------
   ! Leave define mode so we can fill
   !----------------------------------------------------------------------------
   call nc_check(nf90_enddef(ncid), 'WriteNetCDF', 'enddef '//trim(fname))

   !----------------------------------------------------------------------------
   ! Fill the coordinate variables.
   ! The time variable is filled as time progresses.
   !----------------------------------------------------------------------------

   call file_to_text('input.nml', textblock)
   call nc_check(nf90_put_var(ncid, nmlVarID, textblock ), &
              'WriteNetCDF', 'namelist:put_var')
   deallocate(textblock)

   call nc_check(nf90_put_var(ncid, CopyVarId, (/ (i,i=1,Ncopies) /) ), &
              'WriteNetCDF', 'copy:put_var')

   call nc_check(nf90_put_var(ncid, CopyMetaVarID, copy_names), &
              'WriteNetCDF', 'copymeta:put_var')

   call nc_check(nf90_put_var(ncid, TypesVarId, (/ (i,i=1,max_obs_kinds) /) ), &
              'WriteNetCDF', 'types:put_var')

   call nc_check(nf90_put_var(ncid, TypesMetaVarID, my_obs_kind_names(1:max_obs_kinds)), &
              'WriteNetCDF', 'typesmeta:put_var')

   call nc_check(nf90_put_var(ncid, BoundsVarID, (/ 1, 2 /)), &
              'WriteNetCDF', 'bounds:put_var')

   call nc_check(nf90_put_var(ncid, TimeBoundsVarID, epoch_edges ), &
              'WriteNetCDF', 'time_bounds:put_var')

   call nc_check(nf90_sync( ncid), 'WriteNetCDF', 'sync '//trim(fname))  

   !----------------------------------------------------------------------------
   ! finish ...
   !----------------------------------------------------------------------------

   call nc_check(nf90_sync( ncid), 'WriteNetCDF', 'sync '//trim(fname))  
   call nc_check(nf90_close(ncid), 'init_diag_output', 'close '//trim(fname))  

   end Subroutine WriteNetCDF



   Subroutine  SetIndices( obs_index, qc_index, dart_qc_index)

   integer, intent(out) :: obs_index, qc_index, dart_qc_index

   ! Using 'seq' from global scope

   integer :: i

   obs_index              = -1
   qc_index               = -1
   dart_qc_index          = -1

   MetaDataLoop : do i=1, get_num_copies(seq)
      if(index(get_copy_meta_data(seq,i), 'observation'              ) > 0) &
                          obs_index = i
   enddo MetaDataLoop

   ! Hmnnn ... sometimes the first QC field is 'Quality Control' or 'NCEP QC index'
   ! to me ... they mean the same thing.

   QCMetaDataLoop : do i=1, get_num_qc(seq)
      if(index(  get_qc_meta_data(seq,i), 'Quality Control'          ) > 0) &
                           qc_index = i
      if(index(  get_qc_meta_data(seq,i), 'NCEP QC index'            ) > 0) &
                           qc_index = i
      if(index(  get_qc_meta_data(seq,i), 'DART quality control'     ) > 0) &
                      dart_qc_index = i
   enddo QCMetaDataLoop

   !--------------------------------------------------------------------
   ! Make sure we find an index for each of them.
   !--------------------------------------------------------------------

   if (               qc_index < 0 ) then 
      write(msgstring,*)'metadata:Quality Control not found' 
      call error_handler(E_MSG,'wind_obs_to_table',msgstring,source,revision,revdate)
   endif
   if (          dart_qc_index < 0 ) then 
      write(msgstring,*)'metadata:DART quality control not found' 
      call error_handler(E_MSG,'wind_obs_to_table',msgstring,source,revision,revdate)
   endif

   ! Only require obs_index to be present; this allows the program
   ! to be run on obs_seq.in files which have no means or spread.  You get
   ! less info from them, but you can still plot locations, etc.

   if ( obs_index < 0 ) then
      write(msgstring,*)'metadata:observation not found'
      call error_handler(E_ERR,'wind_obs_to_table',msgstring,source,revision,revdate)
   endif

   !--------------------------------------------------------------------
   ! Echo what we found.
   !--------------------------------------------------------------------

   if ( verbose ) then
   write(msgstring,'(''observation      index '',i2,'' metadata '',a)') &
        obs_index, trim(get_copy_meta_data(seq,obs_index))
   call error_handler(E_MSG,'wind_obs_to_table',msgstring,source,revision,revdate)
   endif

   if (qc_index > 0 ) then
      write(msgstring,'(''Quality Control      index '',i2,'' metadata '',a)') &
           qc_index,      trim(get_qc_meta_data(seq,     qc_index))
      call error_handler(E_MSG,'wind_obs_to_table',msgstring,source,revision,revdate)
   endif

   if (dart_qc_index > 0 ) then
      write(msgstring,'(''DART quality control index '',i2,'' metadata '',a)') &
           dart_qc_index, trim(get_qc_meta_data(seq,dart_qc_index))
      call error_handler(E_MSG,'wind_obs_to_table',msgstring,source,revision,revdate)
   endif

   end Subroutine SetIndices



   Subroutine WritePairs( wunit, Uflavor, Vflavor, obs_time, obslon, obslat, &
                          obslevel, U_qc, V_qc, &
                          U_obs,     V_obs,        &
                          U_pr_mean, V_pr_mean,    &
                          U_po_mean, V_po_mean )

   integer,            intent(IN) :: wunit, Uflavor, Vflavor
   type(time_type),    intent(IN) :: obs_time
   real(r8),           intent(IN) :: obslon, obslat, obslevel
   integer,            intent(IN) :: U_qc,  V_qc
   real(r8),           intent(IN) :: U_obs, V_obs
   real(r8), optional, intent(IN) :: U_pr_mean, V_pr_mean
   real(r8), optional, intent(IN) :: U_po_mean, V_po_mean

   integer  :: seconds, days, flavor
   real(r8) :: U_pr, U_po, V_pr, V_po

   call get_time(obs_time,seconds,days)

   flavor = Uflavor*100 + Vflavor

   if (present(U_pr_mean) .and. present(V_pr_mean) .and. &
       present(U_po_mean) .and. present(V_po_mean)) then

      U_po = U_po_mean
      U_pr = U_pr_mean
      V_po = V_po_mean
      V_pr = V_pr_mean

      if ( U_qc > QC_MAX_PRIOR ) then
         U_pr = -999.0
         U_po = -999.0
      else if (U_qc > QC_MAX_POSTERIOR ) then
         U_po = -999.0
      endif

      if ( V_qc > QC_MAX_PRIOR ) then
         V_pr = -999.0
         V_po = -999.0
      else if (V_qc > QC_MAX_POSTERIOR ) then
         V_po = -999.0
      endif

      write(wunit,100) flavor, days, seconds, &
                    obslon, obslat, obslevel, &
                    U_qc, V_qc, U_obs, V_obs, &
                    U_pr, V_pr, U_po, V_po
   else

      write(wunit,110) flavor, days, seconds, &
                    obslon, obslat, obslevel, &
                    U_qc, V_qc, U_obs, V_obs
   endif
   
 100 format(i4,1x,i7,1x,i5,2(1x,f10.5),1x,e11.5,2(1x,i3),6(1x,e12.6))
 110 format(i4,1x,i7,1x,i5,2(1x,f10.5),1x,e11.5,2(1x,i3),2(1x,e12.6))

   end Subroutine WritePairs

end program wind_obs_to_table
