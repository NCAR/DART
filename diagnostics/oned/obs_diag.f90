! DART software - Copyright 2004 - 2013 UCAR. This open source software is
! provided by UCAR, "as is", without charge, subject to all terms of use at
! http://www.image.ucar.edu/DAReS/DART/DART_download
!
! $Id$

program obs_diag

!-----------------------------------------------------------------------
! The programs defines a series of epochs (periods of time) and geographic
! regions and accumulates statistics for these epochs and regions.
!-----------------------------------------------------------------------

! In Atmospheric Science, 'spread' has units of standard deviations ...
!
! I should rename some of the variables I use as variances to reflect this.
! 'priorspred' should really be 'priorvar' since you have to accumulate variances
! the math is correct as it is, but the variable names don't make it easy ...

use        types_mod, only : r4, r8, digits12, MISSING_R8, MISSING_R4, metadatalength
use obs_sequence_mod, only : read_obs_seq, obs_type, obs_sequence_type, get_first_obs, &
                             get_obs_from_key, get_obs_def, get_copy_meta_data, &
                             get_obs_time_range, get_time_range_keys, get_num_obs, &
                             get_next_obs, get_num_times, get_obs_values, init_obs, &
                             assignment(=), get_num_copies, static_init_obs_sequence, &
                             get_qc, destroy_obs_sequence, get_last_obs, get_num_qc, &
                             read_obs_seq_header, destroy_obs, get_qc_meta_data
use      obs_def_mod, only : obs_def_type, get_obs_def_error_variance, get_obs_def_time, &
                             get_obs_def_location, get_obs_kind, get_obs_name
use     obs_kind_mod, only : max_obs_kinds, get_obs_kind_name
use     location_mod, only : location_type, get_location, set_location_missing, &
                             write_location, operator(/=), LocationDims
use time_manager_mod, only : time_type, set_date, set_time, get_time, print_time, &
                             print_date, set_calendar_type, get_date, &
                             operator(*), operator(+), operator(-), &
                             operator(>), operator(<), operator(/), &
                             operator(/=), operator(<=), operator(>=)
use    utilities_mod, only : get_unit, open_file, close_file, register_module, &
                             file_exist, error_handler, E_ERR, E_WARN, E_MSG,  &
                             initialize_utilities, logfileunit, nmlfileunit,   &
                             find_namelist_in_file, check_namelist_read,       &
                             nc_check, do_nml_file, do_nml_term, finalize_utilities, &
                             next_file, get_next_filename, find_textfile_dims, &
                             file_to_text
use         sort_mod, only : sort
use   random_seq_mod, only : random_seq_type, init_random_seq, several_random_gaussians

use typeSizes
use netcdf

implicit none

! version controlled file description for error handling, do not edit
character(len=256), parameter :: source   = &
   "$URL$"
character(len=32 ), parameter :: revision = "$Revision$"
character(len=128), parameter :: revdate  = "$Date$"

!---------------------------------------------------------------------
!---------------------------------------------------------------------

integer, parameter :: MaxRegions = 4
integer, parameter :: MaxTrusted = 5
integer, parameter :: stringlength = 32

!---------------------------------------------------------------------
! variables associated with the observation
!---------------------------------------------------------------------

type(obs_sequence_type) :: seq
type(obs_type)          :: observation, next_obs
type(obs_type)          :: obs1, obsN
type(obs_def_type)      :: obs_def
type(location_type)     :: obs_loc

character(len = 129) :: obs_seq_in_file_name
character(len = 129), allocatable, dimension(:) :: obs_seq_filenames
character(len = stringlength), dimension(MaxTrusted) :: trusted_obsname = 'null'

! Storage with fixed size for observation space diagnostics
real(r8), dimension(1) :: prior_mean, posterior_mean, prior_spread, posterior_spread
real(r8) :: pr_mean, po_mean ! same as above, without useless dimension
real(r8) :: pr_sprd, po_sprd ! same as above, without useless dimension

integer :: obs_copy_index, prior_mean_index, posterior_mean_index
integer :: prior_spread_index, posterior_spread_index
integer :: flavor
integer :: num_copies, num_qc, num_obs, max_num_obs, obs_seq_file_id
integer :: num_obs_types

! variables used primarily/exclusively for the rank histogram
integer                :: ens_size, rank_histogram_bin
type(random_seq_type)  :: ran_seq
real(r8)               :: obs_err_var

character(len=129) :: obs_seq_read_format
logical :: pre_I_format

integer,  dimension(2) :: key_bounds
real(r8), dimension(1) :: obs

integer,  allocatable, dimension(:) :: keys
integer,  allocatable, dimension(:) :: ens_copy_index

logical :: out_of_range, is_there_one, keeper

!---------------------------------------------------------------------
! variables associated with quality control
!
! org_qc_index reflects the 'original' QC value of the observation, if any.
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

integer             :: org_qc_index, dart_qc_index
integer             :: qc_integer
integer, parameter  :: QC_MAX = 8
integer, parameter  :: QC_MAX_PRIOR     = 3
integer, parameter  :: QC_MAX_POSTERIOR = 1
integer, dimension(0:QC_MAX) :: qc_counter = 0
real(r8), allocatable, dimension(:) :: qc
real(r8), allocatable, dimension(:) :: copyvals

integer, parameter, dimension(5) :: hist_qcs = (/ 0, 1, 2, 3, 7 /)
integer, parameter, dimension(5) :: trusted_prior_qcs = (/ 0, 1, 2, 3, 7 /)
integer, parameter, dimension(3) :: trusted_poste_qcs = (/ 0, 1,       7 /)
integer :: numqcvals

!-----------------------------------------------------------------------
! Namelist with default values
!
character(len = 129) :: obs_sequence_name = "obs_seq.final"
character(len = 129) :: obs_sequence_list = ""

character(len = stringlength), dimension(MaxTrusted) :: trusted_obs = 'null'

integer :: max_num_bins       = 9999 ! maximum number of temporal bins to consider
integer :: bin_width_days     = -1   ! width of the assimilation bin - seconds
integer :: bin_width_seconds  = -1   ! width of the assimilation bin - days
integer :: init_skip_days     = 0
integer :: init_skip_seconds  = 0
logical :: verbose               = .false.
logical :: outliers_in_histogram = .true.
logical :: create_rank_histogram = .true.
logical :: use_zero_error_obs    = .false.

! index 1 == region 1 == [0.0, 1.0) i.e. Entire domain
! index 2 == region 2 == [0.0, 0.5)
! index 3 == region 3 == [0.5, 1.0)

integer :: Nregions = MaxRegions
real(r8), dimension(MaxRegions) :: lonlim1 = (/ 0.0_r8, 0.0_r8, 0.5_r8, -1.0_r8 /)
real(r8), dimension(MaxRegions) :: lonlim2 = (/ 1.0_r8, 0.5_r8, 1.0_r8, -1.0_r8 /)

character(len=6), dimension(MaxRegions) :: reg_names = &
                                   (/ 'whole ','yin   ','yang  ','bogus '/)

namelist /obs_diag_nml/ obs_sequence_name, obs_sequence_list,  &
                        bin_width_days, bin_width_seconds,     &
                        init_skip_days, init_skip_seconds, max_num_bins, &
                        Nregions, lonlim1, lonlim2, reg_names, &
                        verbose, outliers_in_histogram,        &
                        create_rank_histogram, trusted_obs, use_zero_error_obs

!-----------------------------------------------------------------------
! Variables used to accumulate the statistics.
!-----------------------------------------------------------------------

integer, parameter :: Ncopies = 18
character(len = stringlength), dimension(Ncopies) :: copy_names = &
   (/ 'Nposs      ', 'Nused      ',                               &
      'rmse       ', 'bias       ', 'spread     ', 'totalspread', &
      'NbadDARTQC ', 'observation', 'ens_mean   ',                &
      'N_DARTqc_0 ', 'N_DARTqc_1 ', 'N_DARTqc_2 ', 'N_DARTqc_3 ', &
      'N_DARTqc_4 ', 'N_DARTqc_5 ', 'N_DARTqc_6 ', 'N_DARTqc_7 ', 'N_trusted  ' /)

type TRV_type
   ! statistics by time-region-variable
   integer ::     time_dim = 1
   integer ::   region_dim = 2
   integer :: variable_dim = 3
   character(len=8) :: string
   integer :: num_times = 0, num_regions = 0, num_variables = 0
   integer,  dimension(:,:,:), pointer :: Nposs, Nused, Ntrusted
   real(r8), dimension(:,:,:), pointer :: rmse, bias, spread, totspread
   integer,  dimension(:,:,:), pointer :: NbadDartQC ! # bad DART QC values
   real(r8), dimension(:,:,:), pointer :: observation, ens_mean
   integer,  dimension(:,:,:), pointer :: NDartQC_0, NDartQC_1, NDartQC_2, NDartQC_3
   integer,  dimension(:,:,:), pointer :: NDartQC_4, NDartQC_5, NDartQC_6, NDartQC_7
   integer,  dimension(:,:,:,:), pointer :: hist_bin => NULL()
end type TRV_type

type(TRV_type) :: analy, guess

type(time_type), allocatable, dimension(:)   :: bincenter
type(time_type), allocatable, dimension(:,:) :: binedges
real(digits12),  allocatable, dimension(:)   :: epochcenter
real(digits12),  allocatable, dimension(:,:) :: epochedges
integer,         allocatable, dimension(:)   :: obs_used_in_epoch

!-----------------------------------------------------------------------
! General purpose variables
!-----------------------------------------------------------------------

integer  :: iregion, iepoch, ivar, ifile, num_obs_in_epoch
real(r8) :: rlocation

integer  :: obsindex, i, iunit, ierr, io, ireg
integer  :: seconds, days, Nepochs, Nfiles

integer  :: num_trusted
logical  :: trusted

! List of observations types
character(len = stringlength), pointer, dimension(:) :: obs_type_strings

! These pairs of variables are used when we diagnose which observations
! are far from the background.
integer, parameter :: MaxSigmaBins = 100
integer  :: nsigma(0:MaxSigmaBins) = 0
integer  :: nsigmaUnit, indx

real(r8) :: pr_zscore, po_zscore

type(time_type) :: TimeMin, TimeMax    ! of the entire period of interest
type(time_type) :: binwidth, halfbinwidth
type(time_type) :: seqT1, seqTN        ! first,last time of obs sequence
type(time_type) :: obsT1, obsTN        ! first,last time of all observations
type(time_type) :: obs_time, skip_time

character(len = 129) :: msgstring1, msgstring2
character(len = stringlength) :: str1, str2, str3

!-----------------------------------------------------------------------
! Some variables to keep track of who's rejected why ...
!-----------------------------------------------------------------------

integer :: Nidentity  = 0   ! identity observations are not appropriate.

!=======================================================================
! Get the party started
!=======================================================================

call initialize_utilities('obs_diag')
call register_module(source,revision,revdate)
call static_init_obs_sequence()

num_obs_types = max_obs_kinds ! for compatibility with 3D version
allocate(obs_type_strings(num_obs_types))
do ivar = 1,max_obs_kinds
   obs_type_strings(ivar) = get_obs_kind_name(ivar)
enddo

!----------------------------------------------------------------------
! Read the namelist
!----------------------------------------------------------------------

call find_namelist_in_file('input.nml', 'obs_diag_nml', iunit)
read(iunit, nml = obs_diag_nml, iostat = io)
call check_namelist_read(iunit, io, 'obs_diag_nml')

! Record the namelist values used for the run ...
if (do_nml_file()) write(nmlfileunit, nml=obs_diag_nml)
if (do_nml_term()) write(     *     , nml=obs_diag_nml)

if ((obs_sequence_name /= '') .and. (obs_sequence_list /= '')) then
   write(msgstring1,*)'specify "obs_sequence_name" or "obs_sequence_list"'
   write(msgstring2,*)'set other to an empty string ... i.e. ""'
   call error_handler(E_ERR, 'obs_diag', msgstring1, source, revision, revdate, text2=msgstring2)
endif

num_trusted = DefineTrustedObs()

!----------------------------------------------------------------------
! Check to see if we are including the outlier observations in the
! rank histogram calculation.
!----------------------------------------------------------------------

if ( outliers_in_histogram ) then
   numqcvals = size(hist_qcs)
else
   numqcvals = size(hist_qcs) - 1
endif

!----------------------------------------------------------------------
! Now that we have input, do some checking and setup
!----------------------------------------------------------------------

! Can lie about calendar for low-order models since its all hypothetical

call set_calendar_type('GREGORIAN')

!----------------------------------------------------------------------
! Determine temporal bin characteristics.
! Nepochs is the total number of time intervals of the period requested.
! if the namelist does not specify a start/stop time and binwidth, we will
! presume the first/last times in the input file(s) are to be used.
!----------------------------------------------------------------------

skip_time = set_time(init_skip_seconds, init_skip_days)

call DetermineFilenames(obsT1, obsTN, Nepochs, Nfiles) ! fills obs_seq_filenames array

call DefineTimes() ! Sets binwidth, halfbinwidth

allocate(  bincenter(Nepochs),   binedges(2,Nepochs)) ! time_type
allocate(epochcenter(Nepochs), epochedges(2,Nepochs)) ! 64bit reals for netCDF
allocate( obs_used_in_epoch(Nepochs) )

call SetSchedule(obsT1, Nepochs, binwidth, halfbinwidth, &
                 bincenter, binedges, epochcenter, epochedges)

TimeMin = binedges(1,      1) ! minimum time of interest
TimeMax = binedges(2,Nepochs) ! maximum time of interest
obs_used_in_epoch = 0

!----------------------------------------------------------------------
! Rectify the region namelist information
!----------------------------------------------------------------------

ireg = MaxRegions
Regions: do i = 1,MaxRegions
   if ((lonlim1(i) < 0.0_r8) .or. (lonlim2(i) < 0.0_r8) ) then
      exit Regions
   else
      ireg = i
   endif
enddo Regions
Nregions = min(Nregions, ireg)

if (verbose) then
   write(logfileunit,*)
   write(    *      ,*)
   do i = 1,Nregions
      write(*,'(''Region '',i02,1x,a32,'' : '',2(f10.4,1x))') &
             i, reg_names(i), lonlim1(i), lonlim2(i)
      write(logfileunit,'(''Region '',i02,1x,a32,'' : '',2(f10.4,1x))') &
             i, reg_names(i), lonlim1(i), lonlim2(i)
   enddo
endif


!----------------------------------------------------------------------
! Declares and initializes the guess and analy structures.
!----------------------------------------------------------------------

call PrepareVariables()

!----------------------------------------------------------------------
! Open file for histogram of innovations, as a function of standard deviation.
!----------------------------------------------------------------------

nsigmaUnit = open_file('LargeInnov.txt',form='formatted',action='rewind')
write(nsigmaUnit,'(a)')'Any observations flagged as bad are dumped into the last bin.'
write(nsigmaUnit,'(a)') '   day   secs    loc            obs         prior   zscore   key   kind'

ObsFileLoop : do ifile=1, Nfiles
!-----------------------------------------------------------------------

   write(*,*)'Reading file ',ifile, trim(obs_seq_filenames(ifile))

   ! Read in information about observation sequence so we can allocate
   ! observations. We need info about how many copies, qc values, etc.
   ! We have already read this sequence once, so no caution required.

   call read_obs_seq_header(obs_seq_filenames(ifile), &
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

   if ( verbose ) then
      write(logfileunit,*)
      write(logfileunit,*)'num_copies          is ',num_copies
      write(logfileunit,*)'num_qc              is ',num_qc
      write(logfileunit,*)'num_obs             is ',num_obs
      write(logfileunit,*)'max_num_obs         is ',max_num_obs
      write(logfileunit,*)'obs_seq_read_format is ',trim(obs_seq_read_format)
      write(    *      ,*)
      write(    *      ,*)'num_copies          is ',num_copies
      write(    *      ,*)'num_qc              is ',num_qc
      write(    *      ,*)'num_obs             is ',num_obs
      write(    *      ,*)'max_num_obs         is ',max_num_obs
      write(    *      ,*)'obs_seq_read_format is ',trim(obs_seq_read_format)
   endif

   ! Read in the entire observation sequence

   call read_obs_seq(obs_seq_filenames(ifile), 0, 0, 0, seq)

   !--------------------------------------------------------------------
   ! The observations for the low-order models are all exactly at
   ! the assimilation timestep. So we know the bin separation.
   !--------------------------------------------------------------------

   is_there_one = get_first_obs(seq, obs1)  ! already checked that this is true.
   call get_obs_def(obs1,     obs_def)
   seqT1        = get_obs_def_time(obs_def)

   is_there_one = get_last_obs(seq, obsN)  ! already checked that this is true.
   call get_obs_def(obsN,     obs_def)
   seqTN        = get_obs_def_time(obs_def)

   !--------------------------------------------------------------------
   ! If the last observation is before the period of interest, move on.
   !--------------------------------------------------------------------

   if ( seqTN < TimeMin ) then
      if (verbose) then
         write(logfileunit,*)'seqTN < TimeMin ... trying next file.'
         write(    *      ,*)'seqTN < TimeMin ... trying next file.'
      endif
      call destroy_obs(obs1)
      call destroy_obs(obsN)
      call destroy_obs(observation)
      call destroy_obs(next_obs)
      call destroy_obs_sequence(seq)
      if (allocated(qc)) deallocate( qc )
      if (allocated(copyvals)) deallocate( copyvals )
      cycle ObsFileLoop
   else
      if (verbose) then
         write(logfileunit,*)'seqTN > TimeMin ... using ', trim(obs_seq_filenames(ifile))
         write(    *      ,*)'seqTN > TimeMin ... using ', trim(obs_seq_filenames(ifile))
      endif
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
      if (verbose) then 
         write(logfileunit,*)'seqT1 < TimeMax ... using ', trim(obs_seq_filenames(ifile))
         write(    *      ,*)'seqT1 < TimeMax ... using ', trim(obs_seq_filenames(ifile))
      endif
   endif

   !--------------------------------------------------------------------
   ! Find the index of obs, ensemble mean, spread ... etc.
   ! Each observation sequence file can have its copies in any order.
   !--------------------------------------------------------------------

   ! FIXME : Make sure this observation sequence file has the same
   ! number of ensemble members as 'the last one' ...

   ens_size = GetEnsSize()

   if (ens_size > 0) then
      if (allocated(ens_copy_index)) then
         if (size(ens_copy_index) /= ens_size) then
            write(msgstring1,'(''expecting '',i3,'' ensemble members, got '',i3)') &
                                       size(ens_copy_index), ens_size
            call error_handler(E_ERR,'obs_diag',msgstring1,source,revision,revdate)
         endif
      else
         ! This should happen exactly once, if at all.
         allocate(guess%hist_bin( Nepochs, Nregions, num_obs_types, ens_size+1))
         allocate(ens_copy_index(ens_size))
         guess%hist_bin    = 0
         call init_random_seq(ran_seq, seed=23)
      endif
      if (verbose) then
         write(msgstring1,*) 'Creating rank histogram with ',ens_size+1,' bins.'
         call error_handler(E_MSG,'obs_diag',msgstring1)
      endif 
   else
      call error_handler(E_MSG,'obs_diag','Cannot create rank histogram.')
      create_rank_histogram = .false.
   endif

   call SetIndices()

   !====================================================================
   EpochLoop : do iepoch = 1, Nepochs
   !====================================================================

      call get_obs_time_range(seq, binedges(1,iepoch), binedges(2,iepoch), key_bounds, &
                  num_obs_in_epoch, out_of_range )

      if( num_obs_in_epoch == 0 ) then
         if ( verbose ) then
            call print_time(binedges(1,iepoch),' epoch  start ',logfileunit)
            call print_time(bincenter( iepoch),' epoch center ',logfileunit)
            call print_time(binedges(2,iepoch),' epoch    end ',logfileunit)
            call print_time(binedges(1,iepoch),' epoch  start ')
            call print_time(bincenter( iepoch),' epoch center ')
            call print_time(binedges(2,iepoch),' epoch    end ')
            write(msgstring1,*)' No observations in epoch ',iepoch,' cycling ...'
            call error_handler(E_MSG,'obs_diag',msgstring1)
         endif
         cycle EpochLoop
      endif

      if ( verbose ) then
         call print_time(binedges(1,iepoch),' epoch  start ',logfileunit)
         call print_time(bincenter( iepoch),' epoch center ',logfileunit)
         call print_time(binedges(2,iepoch),' epoch    end ',logfileunit)
         call print_time(binedges(1,iepoch),' epoch  start ')
         call print_time(bincenter( iepoch),' epoch center ')
         call print_time(binedges(2,iepoch),' epoch    end ')
      endif

      allocate(keys(num_obs_in_epoch))

      call get_time_range_keys(seq, key_bounds, num_obs_in_epoch, keys)

      !-----------------------------------------------------------------
      ObservationLoop : do obsindex = 1, num_obs_in_epoch
      !-----------------------------------------------------------------

         call get_obs_from_key(seq, keys(obsindex), observation)
         call get_obs_def(observation, obs_def)

         flavor      = get_obs_kind(obs_def) ! this is (almost) always [1,max_obs_kinds]
         obs_time    = get_obs_def_time(obs_def)
         obs_loc     = get_obs_def_location(obs_def)
         rlocation   = get_location(obs_loc)

         ! Check to make sure we are past the burn-in 
         if (obs_time < skip_time) cycle ObservationLoop

         ! Check to see if this is a trusted observation
         trusted = .false.
         if ( num_trusted > 0 ) then
            trusted = is_observation_trusted( get_obs_kind_name(flavor) )
         endif

         if ( use_zero_error_obs ) then
            obs_err_var = 0.0_r8
         else
            obs_err_var = get_obs_def_error_variance(obs_def)
         endif

         ! Check to see if it is an identity observation.
         ! If it is, we count them and skip them since they are better
         ! explored with the model-space diagnostics.
         if ( flavor < 0 ) then
         !  write(*,*)'obs ',obsindex,' is an identity observation - no fair.',obs_err_var
            Nidentity = Nidentity + 1
            cycle ObservationLoop
         endif

         !--------------------------------------------------------------
         ! retrieve observation prior and posterior means and spreads
         !--------------------------------------------------------------

         call get_obs_values(observation,              obs,         obs_copy_index)
         call get_obs_values(observation,       prior_mean,       prior_mean_index)
         call get_obs_values(observation,   posterior_mean,   posterior_mean_index)
         call get_obs_values(observation,     prior_spread,     prior_spread_index)
         call get_obs_values(observation, posterior_spread, posterior_spread_index)

         pr_mean =       prior_mean(1)
         po_mean =   posterior_mean(1)
         pr_sprd =     prior_spread(1)
         po_sprd = posterior_spread(1)

         !--------------------------------------------------------------
         ! Convert the DART QC data to an integer and create histogram
         !--------------------------------------------------------------

         call get_qc(observation, qc) ! populates 'qc' with ALL qc values

         if ( dart_qc_index > 0 ) then
            qc_integer = min( nint(qc(dart_qc_index)), QC_MAX )
            qc_counter(qc_integer) = qc_counter(qc_integer) + 1  ! histogram
         else
            ! Provide backwards compatibility. If no dart_qc in obs_seq,
            ! put qc_integer to 0 to replicate logic to be unable to treat
            ! prior and posterior separately.
            qc_integer = 0
         endif

         !--------------------------------------------------------------
         ! (DEBUG) Summary of observation knowledge at this point
         !--------------------------------------------------------------

         if ( 1 == 2 ) then
            call print_time(obs_time,'time is')
            call print_time(obs_time,'time is',logfileunit)
            write(*,*)'observation # ',obsindex
            write(*,*)'obs_flavor ',flavor
            write(*,*)'obs_err_var ',obs_err_var
            write(*,*)'qc ',qc
            write(*,*)'obs(1) ',obs(1)
            write(*,*)'pr_mean, po_mean ',pr_mean, po_mean
            write(*,*)'pr_sprd, po_sprd ',pr_sprd, po_sprd
         endif

         !--------------------------------------------------------------
         ! update the histogram of the magnitude of the innovation,
         ! where each bin is a single standard deviation. This is
         ! a one-sided histogram.
         !--------------------------------------------------------------

         pr_zscore = InnovZscore(obs(1), pr_mean, pr_sprd, obs_err_var, qc_integer, QC_MAX_PRIOR)
         po_zscore = InnovZscore(obs(1), po_mean, po_sprd, obs_err_var, qc_integer, QC_MAX_POSTERIOR)

         indx         = min(int(pr_zscore), MaxSigmaBins)
         nsigma(indx) = nsigma(indx) + 1

         ! Individual (valid) observations that are very far away get
         ! logged to a separate file.

         if( (pr_zscore > 3.0_r8) .and. (qc_integer <= QC_MAX_PRIOR) ) then
            call get_time(obs_time,seconds,days)

            write(nsigmaUnit,FMT='(i7,1x,i5,1x,f8.2,1x,2f13.2,f8.1,2i7)') &
                 days, seconds, rlocation, &
                 obs(1), pr_mean, pr_zscore, keys(obsindex), flavor
         endif

         obs_used_in_epoch(iepoch) = obs_used_in_epoch(iepoch) + 1

         !--------------------------------------------------------------
         ! If needed, calculate the rank histogram bin (once!) for 
         ! this observation - even if the QC value is bad.
         !--------------------------------------------------------------

         if ( create_rank_histogram ) then
            call get_obs_values(observation, copyvals)
            rank_histogram_bin = Rank_Histogram(copyvals, obs_copy_index, obs_err_var)
         endif

         !--------------------------------------------------------------
         ! We have Nregions of interest.
         ! FIXME: suppot if the region of interest is [0.8, 0.2]
         !--------------------------------------------------------------

         Areas : do iregion =1, Nregions

            keeper = InRegion( rlocation, lonlim1(iregion), lonlim2(iregion) )
            if ( .not. keeper ) cycle Areas

            !-----------------------------------------------------------
            ! Count DART QC values 
            !-----------------------------------------------------------

            call count_QC_values(qc_integer, iepoch, iregion, flavor)

            !-----------------------------------------------------------
            ! Do all the heavy lifting
            !-----------------------------------------------------------

            call Bin3D(qc_integer, iepoch, iregion, flavor, trusted, obs(1), &
                obs_err_var, pr_mean, pr_sprd, po_mean, po_sprd, rank_histogram_bin)

         enddo Areas

      !-----------------------------------------------------------------
      enddo ObservationLoop
      !-----------------------------------------------------------------

      deallocate(keys)

      if(verbose) then
         write(msgstring1,'(''num obs considered in epoch '',i4,'' = '',i8, &
                                  & '' out of '',i8,'' possible'')') &
                         iepoch, obs_used_in_epoch(iepoch), num_obs_in_epoch
         call error_handler(E_MSG,'obs_diag',msgstring1,source,revision,revdate)
         write(logfileunit,*)''
         write(     *     ,*)''
      endif

   enddo EpochLoop

   if (verbose) then
      write(logfileunit,*)'End of EpochLoop for ',trim(obs_seq_filenames(ifile))
      write(     *     ,*)'End of EpochLoop for ',trim(obs_seq_filenames(ifile))
   endif

   call destroy_obs(obs1)
   call destroy_obs(obsN)
   call destroy_obs(observation)
   call destroy_obs(next_obs)
   call destroy_obs_sequence(seq)
   if (allocated(qc))       deallocate( qc )
   if (allocated(copyvals)) deallocate( copyvals )

enddo ObsFileLoop

!-----------------------------------------------------------------------
! We have read all possible files, and stuffed the observations into the
! appropriate bins. Time to normalize the prior and posterior structures.
!-----------------------------------------------------------------------

call NormalizeTRV()

if (sum(obs_used_in_epoch) == 0 ) then
   call error_handler(E_ERR,'obs_diag','All identity observations. Stopping.', &
                     source, revision, revdate)
endif

!-----------------------------------------------------------------------
! Print final summary.
!-----------------------------------------------------------------------

write(*,*)
write(*,*) '# observations used  : ',sum(obs_used_in_epoch)
write(*,*) 'Count summary over all regions - obs may count for multiple regions:'
write(*,*) '# identity           : ',Nidentity
write(*,*) '# bad DART QC prior  : ',sum(guess%NbadDartQC)
write(*,*) '# bad DART QC post   : ',sum(analy%NbadDartQC)
write(*,*) '# TRUSTED            : ',sum(analy%Ntrusted)
write(*,*)
write(*,*) '# prior DART QC 0 : ',sum(guess%NDartQC_0)
write(*,*) '# prior DART QC 1 : ',sum(guess%NDartQC_1)
write(*,*) '# prior DART QC 2 : ',sum(guess%NDartQC_2)
write(*,*) '# prior DART QC 3 : ',sum(guess%NDartQC_3)
write(*,*) '# prior DART QC 4 : ',sum(guess%NDartQC_4)
write(*,*) '# prior DART QC 5 : ',sum(guess%NDartQC_5)
write(*,*) '# prior DART QC 6 : ',sum(guess%NDartQC_6)
write(*,*) '# prior DART QC 7 : ',sum(guess%NDartQC_7)
write(*,*)
write(*,*) '# poste DART QC 0 : ',sum(analy%NDartQC_0)
write(*,*) '# poste DART QC 1 : ',sum(analy%NDartQC_1)
write(*,*) '# poste DART QC 2 : ',sum(analy%NDartQC_2)
write(*,*) '# poste DART QC 3 : ',sum(analy%NDartQC_3)
write(*,*) '# poste DART QC 4 : ',sum(analy%NDartQC_4)
write(*,*) '# poste DART QC 5 : ',sum(analy%NDartQC_5)
write(*,*) '# poste DART QC 6 : ',sum(analy%NDartQC_6)
write(*,*) '# poste DART QC 7 : ',sum(analy%NDartQC_7)

write(logfileunit,*)
write(logfileunit,*) '# observations used  : ',sum(obs_used_in_epoch)
write(logfileunit,*) 'Count summary over all regions - obs may count for multiple regions:'
write(logfileunit,*) '# identity           : ',Nidentity
write(logfileunit,*) '# bad DART QC prior  : ',sum(guess%NbadDartQC)
write(logfileunit,*) '# bad DART QC post   : ',sum(analy%NbadDartQC)
write(logfileunit,*) '# TRUSTED            : ',sum(analy%Ntrusted)
write(logfileunit,*)
write(logfileunit,*) '# prior DART QC 0 : ',sum(guess%NDartQC_0)
write(logfileunit,*) '# prior DART QC 1 : ',sum(guess%NDartQC_1)
write(logfileunit,*) '# prior DART QC 2 : ',sum(guess%NDartQC_2)
write(logfileunit,*) '# prior DART QC 3 : ',sum(guess%NDartQC_3)
write(logfileunit,*) '# prior DART QC 4 : ',sum(guess%NDartQC_4)
write(logfileunit,*) '# prior DART QC 5 : ',sum(guess%NDartQC_5)
write(logfileunit,*) '# prior DART QC 6 : ',sum(guess%NDartQC_6)
write(logfileunit,*) '# prior DART QC 7 : ',sum(guess%NDartQC_7)
write(logfileunit,*)
write(logfileunit,*) '# poste DART QC 0 : ',sum(analy%NDartQC_0)
write(logfileunit,*) '# poste DART QC 1 : ',sum(analy%NDartQC_1)
write(logfileunit,*) '# poste DART QC 2 : ',sum(analy%NDartQC_2)
write(logfileunit,*) '# poste DART QC 3 : ',sum(analy%NDartQC_3)
write(logfileunit,*) '# poste DART QC 4 : ',sum(analy%NDartQC_4)
write(logfileunit,*) '# poste DART QC 5 : ',sum(analy%NDartQC_5)
write(logfileunit,*) '# poste DART QC 6 : ',sum(analy%NDartQC_6)
write(logfileunit,*) '# poste DART QC 7 : ',sum(analy%NDartQC_7)

! Print the histogram of innovations as a function of standard deviation. 
write(     *     ,*)
write(     *     ,*) 'Table that reflects the outlier_threshold -- '
write(     *     ,*) 'How are the (good) innovations distributed?'
write(logfileunit,*)
write(logfileunit,*) 'Table that reflects the outlier_threshold -- '
write(logfileunit,*) 'How are the (good) innovations distributed?'
do i=0,MaxSigmaBins
   if(nsigma(i) /= 0) then
      write(     *     ,*)'innovations within ',i+1,' stdev = ',nsigma(i)
      write(logfileunit,*)'innovations within ',i+1,' stdev = ',nsigma(i)
   endif
enddo
write(     *     ,*)
write(logfileunit,*)

!----------------------------------------------------------------------
! Open netCDF output file 
!----------------------------------------------------------------------

call WriteNetCDF('obs_diag_output.nc')

!-----------------------------------------------------------------------
! Really, really, done.
!-----------------------------------------------------------------------

call DestroyVariables()

call error_handler(E_MSG,'obs_diag','Finished successfully.')
call finalize_utilities()

!======================================================================
CONTAINS
!======================================================================
! These routines use common variables from the scope of this file.
! If it's not in the argument list ... it's scoped within this file.
!======================================================================


Subroutine PrepareVariables()

allocate(guess%rmse(       Nepochs, Nregions, num_obs_types), &
         guess%bias(       Nepochs, Nregions, num_obs_types), &
         guess%spread(     Nepochs, Nregions, num_obs_types), &
         guess%totspread(  Nepochs, Nregions, num_obs_types), &
         guess%observation(Nepochs, Nregions, num_obs_types), &
         guess%ens_mean(   Nepochs, Nregions, num_obs_types), &
         guess%Nposs(      Nepochs, Nregions, num_obs_types), &
         guess%Nused(      Nepochs, Nregions, num_obs_types), &
         guess%NbadDartQC( Nepochs, Nregions, num_obs_types), &
         guess%NDartQC_0(  Nepochs, Nregions, num_obs_types), &
         guess%NDartQC_1(  Nepochs, Nregions, num_obs_types), &
         guess%NDartQC_2(  Nepochs, Nregions, num_obs_types), &
         guess%NDartQC_3(  Nepochs, Nregions, num_obs_types), &
         guess%NDartQC_4(  Nepochs, Nregions, num_obs_types), &
         guess%NDartQC_5(  Nepochs, Nregions, num_obs_types), &
         guess%NDartQC_6(  Nepochs, Nregions, num_obs_types), &
         guess%NDartQC_7(  Nepochs, Nregions, num_obs_types), &
         guess%Ntrusted(   Nepochs, Nregions, num_obs_types)  )

guess%rmse        = 0.0_r8
guess%bias        = 0.0_r8
guess%spread      = 0.0_r8
guess%totspread   = 0.0_r8
guess%observation = 0.0_r8
guess%ens_mean    = 0.0_r8
guess%Nposs       = 0
guess%Nused       = 0
guess%NbadDartQC  = 0
guess%NDartQC_0   = 0
guess%NDartQC_1   = 0
guess%NDartQC_2   = 0
guess%NDartQC_3   = 0
guess%NDartQC_4   = 0
guess%NDartQC_5   = 0
guess%NDartQC_6   = 0
guess%NDartQC_7   = 0
guess%Ntrusted    = 0

guess%string        = 'guess'
guess%num_times     = Nepochs
guess%num_regions   = Nregions
guess%num_variables = num_obs_types

allocate(analy%rmse(       Nepochs, Nregions, num_obs_types), &
         analy%bias(       Nepochs, Nregions, num_obs_types), &
         analy%spread(     Nepochs, Nregions, num_obs_types), &
         analy%totspread(  Nepochs, Nregions, num_obs_types), &
         analy%observation(Nepochs, Nregions, num_obs_types), &
         analy%ens_mean(   Nepochs, Nregions, num_obs_types), &
         analy%Nposs(      Nepochs, Nregions, num_obs_types), &
         analy%Nused(      Nepochs, Nregions, num_obs_types), &
         analy%NbadDartQC( Nepochs, Nregions, num_obs_types), &
         analy%NDartQC_0(  Nepochs, Nregions, num_obs_types), &
         analy%NDartQC_1(  Nepochs, Nregions, num_obs_types), &
         analy%NDartQC_2(  Nepochs, Nregions, num_obs_types), &
         analy%NDartQC_3(  Nepochs, Nregions, num_obs_types), &
         analy%NDartQC_4(  Nepochs, Nregions, num_obs_types), &
         analy%NDartQC_5(  Nepochs, Nregions, num_obs_types), &
         analy%NDartQC_6(  Nepochs, Nregions, num_obs_types), &
         analy%NDartQC_7(  Nepochs, Nregions, num_obs_types), &
         analy%Ntrusted(   Nepochs, Nregions, num_obs_types)  )

analy%rmse        = 0.0_r8
analy%bias        = 0.0_r8
analy%spread      = 0.0_r8
analy%totspread   = 0.0_r8
analy%observation = 0.0_r8
analy%ens_mean    = 0.0_r8
analy%Nposs       = 0
analy%Nused       = 0
analy%NbadDartQC  = 0
analy%NDartQC_0   = 0
analy%NDartQC_1   = 0
analy%NDartQC_2   = 0
analy%NDartQC_3   = 0
analy%NDartQC_4   = 0
analy%NDartQC_5   = 0
analy%NDartQC_6   = 0
analy%NDartQC_7   = 0
analy%Ntrusted    = 0

analy%string        = 'analy'
analy%num_times     = Nepochs
analy%num_regions   = Nregions
analy%num_variables = num_obs_types

end Subroutine PrepareVariables


!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


Subroutine DestroyVariables()

deallocate(obs_seq_filenames)

deallocate(guess%rmse,        guess%bias,      guess%spread,    guess%totspread, &
           guess%observation, guess%ens_mean,  guess%Nposs,     guess%Nused,     &
                                               guess%NbadDartQC,guess%Ntrusted,  &
           guess%NDartQC_0,   guess%NDartQC_1, guess%NDartQC_2, guess%NDartQC_3, &
           guess%NDartQC_4,   guess%NDartQC_5, guess%NDartQC_6, guess%NDartQC_7)

if (associated(guess%hist_bin)) deallocate(guess%hist_bin)

if (allocated(ens_copy_index)) deallocate(ens_copy_index)

deallocate(analy%rmse,        analy%bias,      analy%spread,    analy%totspread, &
           analy%observation, analy%ens_mean,  analy%Nposs,     analy%Nused,     &
                                               analy%NbadDartQC,analy%Ntrusted,  &
           analy%NDartQC_0,   analy%NDartQC_1, analy%NDartQC_2, analy%NDartQC_3, &
           analy%NDartQC_4,   analy%NDartQC_5, analy%NDartQC_6, analy%NDartQC_7)

deallocate(epochcenter, epochedges, bincenter, binedges, obs_used_in_epoch)
deallocate(obs_type_strings)

end Subroutine DestroyVariables



!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


Function InnovZscore(obsval, prmean, prspred, errvar, qcval, qcmaxval)

! This function tries to get a handle on the magnitude of the innovations.
! If the ratio of the observation to the prior mean is 'big', it is an outlier.
! If the prior mean cannot be calculated (i.e. is missing) we put it in the
! last 'bin' of the crude histogram. This is pretty much a 'z' score in the
! statistical sense.

real(r8)             :: InnovZscore
real(r8), intent(in) :: obsval, prmean, prspred, errvar
integer,  intent(in) :: qcval, qcmaxval

real(r8) :: numer, denom

if ( qcval <= qcmaxval ) then ! QC indicates a valid obs
   numer = abs(prmean - obsval)
   denom = sqrt( prspred**2 + errvar )
   InnovZscore = numer / denom
else
   InnovZscore = real(MaxSigmaBins,r8)
endif

end Function InnovZscore


!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


Function InRegion( lon, lon1, lon2 ) result( keeper )
! FIXME ... this does not wrap around the origin ... 
!       ... are all 1D locations periodic
real(r8), intent(in) :: lon, lon1, lon2
logical :: keeper

keeper = .false.

if( (lon .ge. lon1) .and. (lon .lt. lon2) ) keeper = .true.

end Function InRegion



!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


Subroutine DetermineFilenames( time1, timeN, nsteps, numfiles )
! This routine reads through all the observation sequence files to determine the
! first and last times in the files. Can only do this by reading each and every
! observation sequence file.

type(time_type), intent(out) :: time1    ! first observation time
type(time_type), intent(out) :: timeN    ! last observation time
integer,         intent(out) :: nsteps   ! number of unique observation times in all files
integer,         intent(out) :: numfiles ! number of observation sequence files

integer :: seqNsteps

numfiles = 0
nsteps   = 0

allocate(obs_seq_filenames(1000))  ! GLOBAL scope
obs_seq_filenames = 'null'

TimeLoop : do ifile = 1, size(obs_seq_filenames)

   if (obs_sequence_list == '') then ! try to increment filename
      obs_seq_in_file_name = next_file(obs_sequence_name,ifile)
   else
      obs_seq_in_file_name = get_next_filename(obs_sequence_list,ifile)
      if (obs_seq_in_file_name == '') exit TimeLoop
   endif

   if ( file_exist(trim(obs_seq_in_file_name)) ) then
      write(msgstring1,*)'opening ', trim(obs_seq_in_file_name)
      call error_handler(E_MSG,'DetermineFilenames',msgstring1,source,revision,revdate)
   else
      if (numfiles < 1) then
         write(msgstring1,*)trim(obs_seq_in_file_name),&
                        ' does not exist. No observation files.'
         call error_handler(E_ERR,'DetermineFilenames',msgstring1,source,revision,revdate)
      else
         write(msgstring1,*)trim(obs_seq_in_file_name),&
                        ' does not exist. No more observation files.'
         call error_handler(E_MSG,'DetermineFilenames',msgstring1,source,revision,revdate)
      endif
      exit TimeLoop
   endif

   ! save a copy of desired filenames
   obs_seq_filenames(ifile) = trim(obs_seq_in_file_name)
   numfiles                 = ifile

   ! Read in information about observation sequence so we can allocate
   ! observations. We need info about how many copies, qc values, etc.

   call read_obs_seq_header(obs_seq_in_file_name, &
             num_copies, num_qc, num_obs, max_num_obs, &
             obs_seq_file_id, obs_seq_read_format, pre_I_format, &
             close_the_file = .true.)

   ! Initialize some (individual) observation variables
   ! Read in the entire observation sequence

   call init_obs( obs1, num_copies, num_qc)
   call init_obs( obsN, num_copies, num_qc)
   call read_obs_seq(obs_seq_in_file_name, 0, 0, 0, seq)

   ! Determine the time encompassed in the observation sequence.

   is_there_one = get_first_obs(seq, obs1)
   if ( .not. is_there_one ) then
      call error_handler(E_ERR,'DetermineFilenames','No first observation in sequence.', &
      source,revision,revdate,text2=obs_seq_in_file_name)
   endif
   call get_obs_def(obs1,   obs_def)
   seqT1 = get_obs_def_time(obs_def)

   is_there_one = get_last_obs(seq, obsN)
   if ( .not. is_there_one ) then
      call error_handler(E_ERR,'DetermineFilenames','No last observation in sequence.', &
      source,revision,revdate,text2=obs_seq_in_file_name)
   endif
   call get_obs_def(obsN,   obs_def)
   seqTN = get_obs_def_time(obs_def)

   seqNsteps = get_num_times(seq)
   nsteps    = nsteps + seqNsteps

   if (ifile == 1) then
      time1  = seqT1
      timeN  = seqTN
   endif

   if (seqT1 < time1) time1 = seqT1
   if (seqTN > timeN) timeN = seqTN

   if ( verbose ) then

      write(logfileunit,*)trim(obs_seq_in_file_name),' has ',seqNsteps,' unique times.'
      write(    *      ,*)trim(obs_seq_in_file_name),' has ',seqNsteps,' unique times.'

      call print_date(seqT1, ' DetermineFilenames: observation 1 date', logfileunit)
      call print_date(seqTN, ' DetermineFilenames: observation N date', logfileunit)
      call print_time(seqT1, ' DetermineFilenames: observation 1 time', logfileunit)
      call print_time(seqTN, ' DetermineFilenames: observation N time', logfileunit)

      call print_date(seqT1, ' DetermineFilenames: observation 1 date')
      call print_date(seqTN, ' DetermineFilenames: observation N date')
      call print_time(seqT1, ' DetermineFilenames: observation 1 time')
      call print_time(seqTN, ' DetermineFilenames: observation N time')

   endif

   call destroy_obs(obs1)
   call destroy_obs(obsN)
   call destroy_obs_sequence(seq)

enddo TimeLoop

if (nsteps < 1) then
   write(msgstring1,*)'cannot find any times in the ',numfiles,' input files.'
   call error_handler(E_ERR,'obs_diag:DetermineFilenames',msgstring1,source,revision,revdate)
endif

if ( verbose ) then

   write(logfileunit,*)
   write(    *      ,*)

   write(logfileunit,*)'Have ',nsteps,' unique times.'
   write(    *      ,*)'Have ',nsteps,' unique times.'

   call print_date(time1, ' DetermineFilenames: first bincenter date',logfileunit)
   call print_date(timeN, ' DetermineFilenames: last  bincenter date',logfileunit)
   call print_time(time1, ' DetermineFilenames: first bincenter time',logfileunit)
   call print_time(timeN, ' DetermineFilenames: last  bincenter time',logfileunit)

   call print_date(time1, ' DetermineFilenames: first bincenter date')
   call print_date(timeN, ' DetermineFilenames: last  bincenter date')
   call print_time(time1, ' DetermineFilenames: first bincenter time')
   call print_time(timeN, ' DetermineFilenames: last  bincenter time')

endif

end Subroutine DetermineFilenames


!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


Subroutine DefineTimes()

!  Sets the binwidth, halfbinwidth

! These are variables that can be modified by this routine
!  integer,         intent(inout) GLOBAL :: Nepochs
!  type(time_type), intent(out)   GLOBAL :: binwidth     ! period of interest around center
!  type(time_type), intent(out)   GLOBAL :: halfbinwidth ! half that period
!  bin_width_days, bin_width_seconds  GLOBAL from namelist

integer :: nbins
type(time_type) :: test_time

! do some error-checking first

if ( (bin_width_days < 0) .and. (bin_width_seconds >= 0) ) then

   write(msgstring1,*)'bin_width_[days,seconds] must be non-negative, they are ', &
   bin_width_days, bin_width_seconds
   call error_handler(E_ERR,'DefineTimes',msgstring1,source,revision,revdate, &
          text2='namelist parameter out-of-bounds. Fix and try again.')

elseif ( (bin_width_days >= 0) .and. (bin_width_seconds < 0) ) then

   write(msgstring1,*)'bin_width_[days,seconds] must be non-negative, they are ', &
   bin_width_days, bin_width_seconds
   call error_handler(E_ERR,'DefineTimes',msgstring1,source,revision,revdate, &
          text2='namelist parameter out-of-bounds. Fix and try again.')

elseif ( (bin_width_days <= 0) .and. (bin_width_seconds <= 0) ) then
   
   ! This is the 'default' case ... use all possible, up to "max_num_bins".
   ! 'space-filling' strategy: bin width and bin separation are same.
   ! Using Nepochs that comes from the number of unique times in the files.

   binwidth  = (obsTN - obsT1) / (Nepochs - 1)
   if (Nepochs > max_num_bins) then
      write(msgstring1,*)'default calculation results in ',Nepochs,' time bins.'
      write(msgstring2,*)'namelist "max_num_bins" requests ',max_num_bins,'. Using this value.'
      call error_handler(E_MSG,'DefineTimes',msgstring1,source,revision,revdate,text2=msgstring2)
      Nepochs = max_num_bins 
      obsTN   = obsT1 + (Nepochs-1)*binwidth
   endif

else
   ! honor the user input
   binwidth  = set_time(bin_width_seconds, bin_width_days)
   test_time = obsT1
   nbins = 0
   COUNTBINS : do i = 1,max_num_bins
      if (test_time >  obsTN) exit COUNTBINS
      test_time = test_time + binwidth
      nbins     = nbins + 1
   enddo COUNTBINS

   ! Warn about falling off end ...
   if ((nbins == max_num_bins) .and. verbose) then 
      write(msgstring1,*)'namelist "max_num_bins" requests ',max_num_bins,'. Using this value.'
      call error_handler(E_MSG,'DefineTimes',msgstring1,source,revision,revdate)
   endif

   Nepochs = nbins 
   obsTN   = obsT1 + (Nepochs-1)*binwidth
endif

halfbinwidth = binwidth / 2

if ( verbose ) then
   write(logfileunit,*)
   write(     *     ,*)

   call print_date(       obsT1,' DefineTimes: start             date',logfileunit)
   call print_date(       obsTN,' DefineTimes: end               date',logfileunit)
   call print_time(       obsT1,' DefineTimes: start             time',logfileunit)
   call print_time(       obsTN,' DefineTimes: end               time',logfileunit)
   call print_time(    binwidth,' DefineTimes: requested     binwidth',logfileunit)
   call print_time(halfbinwidth,' DefineTimes: implied   halfbinwidth',logfileunit)

   call print_date(       obsT1,' DefineTimes: start             date')
   call print_date(       obsTN,' DefineTimes: end               date')
   call print_time(       obsT1,' DefineTimes: start             time')
   call print_time(       obsTN,' DefineTimes: end               time')
   call print_time(    binwidth,' DefineTimes: requested     binwidth')
   call print_time(halfbinwidth,' DefineTimes: implied   halfbinwidth')
endif

end subroutine DefineTimes


!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


! TJH FIXME SetSchedule() should really come from the schedule module

Subroutine SetSchedule(bin1time, num_epochs, fullwidth, halfwidth, &
   bin_center, bin_edges, epoch_center, epoch_edges)

type(time_type),                 intent(in)  :: bin1time
integer,                         intent(in)  :: num_epochs
type(time_type),                 intent(in)  :: fullwidth
type(time_type),                 intent(in)  :: halfwidth
type(time_type), dimension(  :), intent(out) :: bin_center
type(time_type), dimension(:,:), intent(out) :: bin_edges
real(digits12),  dimension(  :), intent(out) :: epoch_center
real(digits12),  dimension(:,:), intent(out) :: epoch_edges

! Explicitly set first epoch - the rest will be predicated on this
! There is a chance the first bin center is at T=0, which is nigh
! impossible to subtract half a bin width from ...
! Since both the bin center and bin width have to be positive, the
! end of the bin MUST be a positive time, so it really pertains
! only to the first bin leading edge.

iepoch = 1
bin_center( iepoch)    = bin1time
if (bin1time <= halfwidth ) then
   bin_edges(1,iepoch) = set_time(0,0)
else
   bin_edges(1,iepoch) = bin1time - halfwidth + set_time(1,0)
endif
bin_edges(2,iepoch)    = bin1time + halfwidth

call get_time(bin_center(iepoch),seconds,days)
epoch_center(iepoch) = days + seconds/86400.0_digits12

call get_time(bin_edges(1,iepoch),seconds,days)
epoch_edges(1,iepoch) = days + seconds/86400.0_digits12

call get_time(bin_edges(2,iepoch),seconds,days)
epoch_edges(2,iepoch) = days + seconds/86400.0_digits12

! Now that we have the first bin center and extent defined ... we roll ...

BinLoop : do iepoch = 2,num_epochs

   bin_center( iepoch) = bin_center(iepoch-1) + fullwidth
   bin_edges(1,iepoch) = bin_center(iepoch) - halfwidth + set_time(1,0)
   bin_edges(2,iepoch) = bin_center(iepoch) + halfwidth

   call get_time(bin_center(iepoch),seconds,days)
   epoch_center(iepoch) = days + seconds/86400.0_digits12

   call get_time(bin_edges(1,iepoch),seconds,days)
   epoch_edges(1,iepoch) = days + seconds/86400.0_digits12

   call get_time(bin_edges(2,iepoch),seconds,days)
   epoch_edges(2,iepoch) = days + seconds/86400.0_digits12

enddo BinLoop

if ( verbose ) then
   do iepoch = 1,num_epochs
      write(logfileunit,*)
      write(     *     ,*)
      write(str1,'(''epoch '',i6,''  start'')')iepoch
      write(str2,'(''epoch '',i6,'' center'')')iepoch
      write(str3,'(''epoch '',i6,''    end'')')iepoch
      call print_time(bin_edges(1,iepoch),str1,logfileunit)
      call print_time(bin_center( iepoch),str2,logfileunit)
      call print_time(bin_edges(2,iepoch),str3,logfileunit)
      call print_time(bin_edges(1,iepoch),str1)
      call print_time(bin_center( iepoch),str2)
      call print_time(bin_edges(2,iepoch),str3)

      call print_date(bin_edges(1,iepoch),str1,logfileunit)
      call print_date(bin_center( iepoch),str2,logfileunit)
      call print_date(bin_edges(2,iepoch),str3,logfileunit)
      call print_date(bin_edges(1,iepoch),str1)
      call print_date(bin_center( iepoch),str2)
      call print_date(bin_edges(2,iepoch),str3)
   enddo
   write(logfileunit,*)
   write(     *     ,*)
endif

end Subroutine SetSchedule


!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


Function GetEnsSize()
!
!  Loop over all the metadata to count the number of ensemble members
!  available in the observation sequence file. We need this count to 
!  allocate space for the rank histogram information. Since the rank
!  histogram will be created for the priors only ...
!

integer :: GetEnsSize

! Using 'seq' from global scope

integer :: i

GetEnsSize = 0

MetaDataLoop : do i=1, get_num_copies(seq)
   if(index(get_copy_meta_data(seq,i), 'prior ensemble member') > 0) &
                   GetEnsSize = GetEnsSize + 1
enddo MetaDataLoop

write(msgstring1,'(''There are '',i4,'' ensemble members.'')') GetEnsSize
call error_handler(E_MSG,'GetEnsSize',msgstring1,source,revision,revdate)

end Function GetEnsSize


!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


Subroutine  SetIndices()

! integer, intent(out) :: obs_copy_index, org_qc_index, dart_qc_index, &
!                         prior_mean_index,   posterior_mean_index,    &
!                         prior_spread_index, posterior_spread_index

! Using 'seq' and 'ens_size' from global scope

integer :: i, ens_count
character(len=metadatalength) :: metadata

obs_copy_index         = -1
org_qc_index           = -1
dart_qc_index          = -1
prior_mean_index       = -1
posterior_mean_index   = -1
prior_spread_index     = -1
posterior_spread_index = -1

ens_count = 0

MetaDataLoop : do i=1, get_num_copies(seq)

   metadata = get_copy_meta_data(seq,i)

   if ( use_zero_error_obs ) then
      if(index(metadata, 'truth'       ) > 0) obs_copy_index = i
   else
      if(index(metadata, 'observation' ) > 0) obs_copy_index = i
   endif

   if(index(metadata, 'prior ensemble mean'      ) > 0) prior_mean_index = i
   if(index(metadata, 'posterior ensemble mean'  ) > 0) posterior_mean_index = i
   if(index(metadata, 'prior ensemble spread'    ) > 0) prior_spread_index = i
   if(index(metadata, 'posterior ensemble spread') > 0) posterior_spread_index = i

   if(index(metadata, 'prior ensemble member'    ) > 0 .and. &
      create_rank_histogram ) then
         ens_count  = ens_count + 1
         ens_copy_index(ens_count) = i
   endif

enddo MetaDataLoop

QCMetaDataLoop : do i=1, get_num_qc(seq)
   metadata = get_qc_meta_data(seq,i)
   if(index(metadata, 'Quality Control'      ) > 0)  org_qc_index = i
   if(index(metadata, 'DART quality control' ) > 0) dart_qc_index = i
enddo QCMetaDataLoop

!--------------------------------------------------------------------
! Make sure we find an index for each of them.
!--------------------------------------------------------------------

if (       prior_mean_index < 0 ) then
   write(msgstring1,*)'metadata:prior ensemble mean not found'
   call error_handler(E_MSG,'SetIndices',msgstring1)
endif
if (   posterior_mean_index < 0 ) then
   write(msgstring1,*)'metadata:posterior ensemble mean not found'
   call error_handler(E_MSG,'SetIndices',msgstring1)
endif
if (     prior_spread_index < 0 ) then
   write(msgstring1,*)'metadata:prior ensemble spread not found'
   call error_handler(E_MSG,'SetIndices',msgstring1)
endif
if ( posterior_spread_index < 0 ) then
   write(msgstring1,*)'metadata:posterior ensemble spread not found'
   call error_handler(E_MSG,'SetIndices',msgstring1)
endif
if (          org_qc_index < 0 ) then
   write(msgstring1,*)'metadata:Quality Control not found'
   call error_handler(E_MSG,'SetIndices',msgstring1)
endif
if (         dart_qc_index < 0 ) then
   write(msgstring1,*)'metadata:DART quality control not found'
   call error_handler(E_MSG,'SetIndices',msgstring1)
endif

! Only require obs_index to be present; this allows the program
! to be run on obs_seq.in files which have no means or spread.

if ( obs_copy_index < 0 ) then
   if ( use_zero_error_obs ) then
      write(msgstring1,*)'metadata:truth       not found'
   else
      write(msgstring1,*)'metadata:observation not found'
   endif
   call error_handler(E_MSG,'SetIndices',msgstring1)
endif

!--------------------------------------------------------------------
! Echo what we found. If we want to.
!--------------------------------------------------------------------
if (verbose) then
   if ( use_zero_error_obs ) then
      write(msgstring1,'(''truth                index '',i2,'' metadata '',a)') &
        obs_copy_index, trim(adjustl(get_copy_meta_data(seq,obs_copy_index)))
   else
      write(msgstring1,'(''observation          index '',i2,'' metadata '',a)') &
        obs_copy_index, trim(adjustl(get_copy_meta_data(seq,obs_copy_index)))
   endif
   call error_handler(E_MSG,'SetIndices',msgstring1,source,revision,revdate)
   
   write(msgstring1,'(''prior mean           index '',i2,'' metadata '',a)') &
        prior_mean_index, trim(adjustl(get_copy_meta_data(seq,prior_mean_index)))
   call error_handler(E_MSG,'SetIndices',msgstring1,source,revision,revdate)
   
   write(msgstring1,'(''posterior mean       index '',i2,'' metadata '',a)') &
        posterior_mean_index, trim(adjustl(get_copy_meta_data(seq,posterior_mean_index)))
   call error_handler(E_MSG,'SetIndices',msgstring1,source,revision,revdate)
   
   write(msgstring1,'(''prior spread         index '',i2,'' metadata '',a)') &
        prior_spread_index, trim(adjustl(get_copy_meta_data(seq,prior_spread_index)))
   call error_handler(E_MSG,'SetIndices',msgstring1,source,revision,revdate)
   
   write(msgstring1,'(''posterior spread     index '',i2,'' metadata '',a)') &
        posterior_spread_index, trim(adjustl(get_copy_meta_data(seq,posterior_spread_index)))
   call error_handler(E_MSG,'SetIndices',msgstring1,source,revision,revdate)
   
   write(msgstring1,'(''Quality Control      index '',i2,'' metadata '',a)') &
        org_qc_index, trim(adjustl(get_qc_meta_data(seq,org_qc_index)))
   call error_handler(E_MSG,'SetIndices',msgstring1,source,revision,revdate)
   
   if (          dart_qc_index > 0 ) then
   write(msgstring1,'(''DART quality control index '',i2,'' metadata '',a)') &
        dart_qc_index, trim(adjustl(get_qc_meta_data(seq,dart_qc_index)))
   call error_handler(E_MSG,'SetIndices',msgstring1,source,revision,revdate)
   endif
endif

end Subroutine SetIndices


!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


Function Rank_Histogram(copyvalues, obs_index, error_variance ) result(rank)

! Calculates the bin/rank
! We don't care about the QC value. If the ob wasn't assimilated
! the bin is meaningless.

real(r8),dimension(:), intent(in)  :: copyvalues
integer,               intent(in)  :: obs_index
real(r8),              intent(in)  :: error_variance
integer                            :: rank

! Local Variables

real(r8)                      :: obsvalue, mean, stddev
real(r8), dimension(ens_size) :: ensemble_values
real(r8), dimension(ens_size) :: sampling_noise

! Grab the observation value from the myriad copy values.

obsvalue = copyvalues(obs_index)
mean     = 0.0_r8
stddev   = sqrt(error_variance)

! Grab the ensemble member values from the myriad copy values.
! It should be noted that the ensemble values here WILL INCLUDE the 
! application of any prior inflation value ...

do i = 1,ens_size
   ensemble_values(i) = copyvalues(ens_copy_index(i))
enddo

! Create an array of 'observation error' ~ N(0,sigma) 
! ran_seq is 'global' and has been initialized already  

call several_random_gaussians(ran_seq, mean, stddev, ens_size, sampling_noise)

! Add the 'sampling noise' to the ensemble values, 
! sort, and find the rank of the observation. Simple.

ensemble_values = sort(ensemble_values + sampling_noise)
 ! ensemble_values = sort(ensemble_values) ! DEBUG - U-shaped valley, here we come.

rank = 0
RankLoop : do i = 1,ens_size

   if ( obsvalue <= ensemble_values(i) ) then
      rank = i
      exit RankLoop
   endif

enddo RankLoop

if (rank == 0) then ! ob is larger than largest ensemble member.
  rank = ens_size + 1
endif


if ( 2 == 1 )  then ! DEBUG block
   write(*,*)'observation error variance is ',error_variance
   write(*,*)'observation          value is ',obsvalue
   write(*,*)'observation           rank is ',rank
   write(*,*)'noisy ensemble values are '
   write(*,*)ensemble_values
   write(*,*)
endif

end Function Rank_Histogram


!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


   Subroutine RPE(x,y)
      real(r8), intent(inout) :: x
      real(r8), intent(in) :: y
      x = x + y
   end Subroutine RPE
   Subroutine IPE(x,y)
      integer, intent(inout) :: x
      integer, intent(in) :: y
      x = x + y
   end Subroutine IPE


!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


  Subroutine Bin3D(iqc, iepoch, iregion, flavor, trusted, &
                obsval,  obserrvar,  prmean,  prsprd,  pomean,  posprd, rank)
   !----------------------------------------------------------------------
   ! The 'guess' and 'analysis' structures are globally scoped.
   ! This function simply accumulates the appropriate sums. 
   ! The normalization occurrs after all the data has been read, naturally.
   !
   ! ... spread is computed via sqrt(ensemble_spread**2 + observation_error**2).  
   ! however, dart stores the variance, not the error, so we do not need
   ! to square it here.
   ! If you are verifying the ensemble against imperfect (real) observations, 
   ! it is necessary to account for the observation error when computing the 
   ! spread.  Since the observation error is not included as output from 
   ! obs_diag, it is not possible to compute this quantity now.

   integer,  intent(in) :: iqc, iepoch, iregion, flavor
   logical,  intent(in) :: trusted
   real(r8), intent(in) :: obsval,  obserrvar,  prmean,  prsprd,  pomean,  posprd
   integer,  intent(in) :: rank

   real(r8) :: priorsqerr      ! PRIOR     Squared Error
   real(r8) :: priorbias       ! PRIOR     simple bias
   real(r8) :: priorspred      ! PRIOR     (spread,variance)
   real(r8) :: priorspredplus  ! PRIOR     (spread,variance**)
   real(r8) :: postsqerr       ! POSTERIOR Squared Error
   real(r8) :: postbias        ! POSTERIOR simple bias
   real(r8) :: postspred       ! POSTERIOR (spread,variance)
   real(r8) :: postspredplus   ! POSTERIOR (spread,variance**)

   real(r8) :: priormean, postmean, obsmean

   priorsqerr     = (prmean - obsval)**2
   postsqerr      = (pomean - obsval)**2
   priorbias      =  prmean - obsval
   postbias       =  pomean - obsval
   priorspred     = prsprd**2
   postspred      = posprd**2
   priorspredplus = prsprd**2 + obserrvar
   postspredplus  = posprd**2 + obserrvar
   obsmean        = obsval
   priormean      = prmean
   postmean       = pomean

   !----------------------------------------------------------------------
   ! Track the number of possible observations
   !----------------------------------------------------------------------

   call IPE(guess%Nposs(iepoch,iregion,flavor), 1)
   call IPE(analy%Nposs(iepoch,iregion,flavor), 1)

   !----------------------------------------------------------------------
   ! Enforce the use of trusted observations.
   !----------------------------------------------------------------------

   if ( trusted ) then

      call IPE(guess%Ntrusted(iepoch,iregion,flavor), 1)
      call IPE(analy%Ntrusted(iepoch,iregion,flavor), 1)

      ! Accrue the PRIOR quantities
      if ( any(iqc == trusted_prior_qcs) ) then
         call IPE(guess%Nused(      iepoch,iregion,flavor),      1    )
         call RPE(guess%observation(iepoch,iregion,flavor), obsmean   )
         call RPE(guess%ens_mean(   iepoch,iregion,flavor), priormean )
         call RPE(guess%bias(       iepoch,iregion,flavor), priorbias )
         call RPE(guess%rmse(       iepoch,iregion,flavor), priorsqerr)
         call RPE(guess%spread(     iepoch,iregion,flavor), priorspred)
         call RPE(guess%totspread(  iepoch,iregion,flavor), priorspredplus)
      else
         call IPE(guess%NbadDartQC(iepoch,iregion,flavor),      1     )
      endif

      ! Accrue the POSTERIOR quantities
      if ( any(iqc == trusted_poste_qcs) ) then
         call IPE(analy%Nused(      iepoch,iregion,flavor),     1    )
         call RPE(analy%observation(iepoch,iregion,flavor), obsmean  )
         call RPE(analy%ens_mean(   iepoch,iregion,flavor), postmean )
         call RPE(analy%bias(       iepoch,iregion,flavor), postbias )
         call RPE(analy%rmse(       iepoch,iregion,flavor), postsqerr)
         call RPE(analy%spread(     iepoch,iregion,flavor), postspred)
         call RPE(analy%totspread(  iepoch,iregion,flavor), postspredplus)
      else
         call IPE(analy%NbadDartQC(iepoch,iregion,flavor),      1    )
      endif

      return  ! EXIT THE BINNING ROUTINE
   endif

   !----------------------------------------------------------------------
   ! Proceed 'as normal'.
   !----------------------------------------------------------------------

   if ( iqc > QC_MAX_PRIOR ) then  ! prior and posterior failed

      call IPE(guess%NbadDartQC(iepoch,iregion,flavor),      1    )
      call IPE(analy%NbadDartQC(iepoch,iregion,flavor),      1    )

   else if ( iqc > QC_MAX_POSTERIOR ) then

      ! Then at least the prior (A.K.A. guess) is good
      call IPE(guess%Nused(      iepoch,iregion,flavor),      1    )
      call RPE(guess%observation(iepoch,iregion,flavor), obsmean   )
      call RPE(guess%ens_mean(   iepoch,iregion,flavor), priormean )
      call RPE(guess%bias(       iepoch,iregion,flavor), priorbias )
      call RPE(guess%rmse(       iepoch,iregion,flavor), priorsqerr)
      call RPE(guess%spread(     iepoch,iregion,flavor), priorspred)
      call RPE(guess%totspread(  iepoch,iregion,flavor), priorspredplus)

      ! However, the posterior is bad
      call IPE(analy%NbadDartQC(iepoch,iregion,flavor),      1    )

   else

      ! The prior is good
      call IPE(guess%Nused(      iepoch,iregion,flavor),      1    )
      call RPE(guess%observation(iepoch,iregion,flavor), obsmean   )
      call RPE(guess%ens_mean(   iepoch,iregion,flavor), priormean )
      call RPE(guess%bias(       iepoch,iregion,flavor), priorbias )
      call RPE(guess%rmse(       iepoch,iregion,flavor), priorsqerr)
      call RPE(guess%spread(     iepoch,iregion,flavor), priorspred)
      call RPE(guess%totspread(  iepoch,iregion,flavor), priorspredplus)

      ! The posterior is good
      call IPE(analy%Nused(      iepoch,iregion,flavor),      1   )
      call RPE(analy%observation(iepoch,iregion,flavor), obsmean  )
      call RPE(analy%ens_mean(   iepoch,iregion,flavor), postmean )
      call RPE(analy%bias(       iepoch,iregion,flavor), postbias )
      call RPE(analy%rmse(       iepoch,iregion,flavor), postsqerr)
      call RPE(analy%spread(     iepoch,iregion,flavor), postspred)
      call RPE(analy%totspread(  iepoch,iregion,flavor), postspredplus)

   endif

   ! The rank histogram binning is a bit of a peculiar situation.
   ! Only the prior is of interest ... so DART QCs of 0 1 2 3 are 'good'.
   ! There is some debate about whether we should be considering the 
   ! 'outlier' observations (DART QC == 7), so that is namelist controlled.

   if (     (rank > 0) .and. create_rank_histogram ) then
      if ( any(iqc == hist_qcs(1:numqcvals) ) )  &
         call IPE(guess%hist_bin(iepoch,iregion,flavor,rank), 1)
   endif

   end Subroutine Bin3D


!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


   Subroutine WriteNetCDF(fname)
   character(len=*), intent(in) :: fname

   integer :: ncid, i, indx1, nobs, typesdimlen
   integer ::  RegionDimID,  RegionVarID
   integer ::    TimeDimID,    TimeVarID
   integer ::    CopyDimID,    CopyVarID,  CopyMetaVarID
   integer ::   TypesDimID,   TypesVarID, TypesMetaVarID
   integer ::  BoundsDimID,  BoundsVarID  
   integer ::    RankDimID,    RankVarID  
   integer ::  StringDimID

   integer :: TimeBoundsVarID, RegionNamesVarID

   character(len=8)      :: crdate      ! needed by F90 DATE_AND_TIME intrinsic
   character(len=10)     :: crtime      ! needed by F90 DATE_AND_TIME intrinsic
   character(len=5)      :: crzone      ! needed by F90 DATE_AND_TIME intrinsic
   integer, dimension(8) :: values      ! needed by F90 DATE_AND_TIME intrinsic

   integer :: year, month, day, hour, minute, second

   if(.not. byteSizesOK()) then
       call error_handler(E_ERR,'WriteNetCDF', &
      'Compiler does not support required kinds of variables.',source,revision,revdate)
   endif

   call nc_check(nf90_create(path = trim(fname), cmode = nf90_share, &
            ncid = ncid), 'obs_diag:WriteNetCDF', 'create '//trim(fname))

   !----------------------------------------------------------------------------
   ! Write Global Attributes 
   !----------------------------------------------------------------------------

   call DATE_AND_TIME(crdate,crtime,crzone,values)
   write(msgstring1,'(''YYYY MM DD HH MM SS = '',i4,5(1x,i2.2))') &
                  values(1), values(2), values(3), values(5), values(6), values(7)
   call nc_check(nf90_put_att(ncid, NF90_GLOBAL, 'creation_date', trim(msgstring1) ), &
              'WriteNetCDF', 'put_att creation_date '//trim(fname))
   call nc_check(nf90_put_att(ncid, NF90_GLOBAL, 'obs_diag_source', source ), &
              'WriteNetCDF', 'put_att obs_diag_source '//trim(fname))
   call nc_check(nf90_put_att(ncid, NF90_GLOBAL, 'obs_diag_revision', revision ), &
              'WriteNetCDF', 'put_att obs_diag_revision '//trim(fname))
   call nc_check(nf90_put_att(ncid, NF90_GLOBAL, 'obs_diag_revdate', revdate ), &
              'WriteNetCDF', 'put_att obs_diag_revdate '//trim(fname))
   call nc_check(nf90_put_att(ncid, NF90_GLOBAL, 'LocationRank', LocationDims ), &
              'WriteNetCDF', 'put_att LocationRank '//trim(fname))

   call nc_check(nf90_put_att(ncid, NF90_GLOBAL, 'bias_convention', &
              'model - observation' ), 'WriteNetCDF', 'put_att bias '//trim(fname))

   if ( create_rank_histogram ) then
      call nc_check(nf90_put_att(ncid, NF90_GLOBAL, 'DART_QCs_in_histogram', &
              hist_qcs(1:numqcvals) ), &
              'WriteNetCDF', 'put_att qcs in histogram '//trim(fname))

      if ( outliers_in_histogram ) then
      call nc_check(nf90_put_att(ncid, NF90_GLOBAL, 'outliers_in_histogram', &
              'TRUE' ), 'WriteNetCDF', 'put_att outliers histogram '//trim(fname))
      else
      call nc_check(nf90_put_att(ncid, NF90_GLOBAL, 'outliers_in_histogram', &
              'FALSE' ), 'WriteNetCDF', 'put_att outliers histogram '//trim(fname))
      endif

   endif

   !----------------------------------------------------------------------------
   ! write all namelist quantities - and some that simply match the 3D version.
   ! That makes it easier to write programs that will handle both.
   !----------------------------------------------------------------------------

   call get_date(bincenter(1), year, month, day, hour, minute, second)
   call nc_check(nf90_put_att(ncid, NF90_GLOBAL, 'first_bin_center', &
                      (/ year, month, day, hour, minute, second /)), &
              'WriteNetCDF', 'put_att first_bin_center '//trim(fname))

   call get_date(bincenter(Nepochs), year, month, day, hour, minute, second)
   call nc_check(nf90_put_att(ncid, NF90_GLOBAL, 'last_bin_center', &
                     (/ year, month, day, hour, minute, second /)), &
              'WriteNetCDF', 'put_att last_bin_center '//trim(fname))

   ! The "space-filling" design forces the bin_separation and the bin_width to be the same
   call get_time(binwidth, seconds, days)
   call nc_check(nf90_put_att(ncid, NF90_GLOBAL, 'bin_width', (/days, seconds/)), &
              'WriteNetCDF', 'put_att bin_width '//trim(fname))
   call nc_check(nf90_put_att(ncid, NF90_GLOBAL, 'bin_separation', (/days, seconds/)), &
              'WriteNetCDF', 'put_att bin_separation '//trim(fname))

   call nc_check(nf90_put_att(ncid, NF90_GLOBAL, 'bin_width_days', days ), &
              'WriteNetCDF', 'put_att bin_width_days '//trim(fname))
   call nc_check(nf90_put_att(ncid, NF90_GLOBAL, 'bin_width_seconds', seconds ), &
              'WriteNetCDF', 'put_att bin_width_seconds '//trim(fname))
   call nc_check(nf90_put_att(ncid, NF90_GLOBAL, 'time_to_skip', &
              (/ 0, 0, init_skip_days, 0, 0, init_skip_seconds /) ), &
              'WriteNetCDF', 'put_att time_to_skip '//trim(fname))
   call nc_check(nf90_put_att(ncid, NF90_GLOBAL, 'max_num_bins', Nepochs ), &
              'WriteNetCDF', 'put_att max_num_bins '//trim(fname))
   call nc_check(nf90_put_att(ncid, NF90_GLOBAL, 'Nregions', Nregions ), &
              'WriteNetCDF', 'put_att Nregions '//trim(fname))
   call nc_check(nf90_put_att(ncid, NF90_GLOBAL, 'lonlim1', lonlim1(1:Nregions) ), &
              'WriteNetCDF', 'put_att lonlim1 '//trim(fname))
   call nc_check(nf90_put_att(ncid, NF90_GLOBAL, 'lonlim2', lonlim2(1:Nregions) ), &
              'WriteNetCDF', 'put_att lonlim2 '//trim(fname))

   !----------------------------------------------------------------------------
   ! write all observation sequence files used
   !----------------------------------------------------------------------------

   FILEloop : do i = 1,SIZE(obs_seq_filenames)

     indx1 = index(obs_seq_filenames(i),'null')

     if (indx1 > 0) exit FILEloop

     write(msgstring1,'(''obs_seq_file_'',i3.3)')i
     call nc_check(nf90_put_att(ncid, NF90_GLOBAL, &
            trim(msgstring1), trim(obs_seq_filenames(i)) ), &
            'WriteNetCDF', 'region_names:obs_kinds')

   enddo FILEloop

   call nc_check(nf90_put_att(ncid, NF90_GLOBAL, 'NumIdentityObs', Nidentity ), &
              'WriteNetCDF', 'put_att identity '//trim(fname))

   !----------------------------------------------------------------------------
   ! Write all observation types that are used. Requires counting how many
   ! observations for each observation type.
   !----------------------------------------------------------------------------

   call nc_check(nf90_put_att(ncid, NF90_GLOBAL, 'comment', &
           'All used observation types follow. &
           &ObservationTypes variable has all types known.' ), &
           'WriteNetCDF', 'put_att obstypes comment '//trim(fname))

   typesdimlen = 0
   do ivar = 1,max_obs_kinds

      nobs = sum(analy%Nposs(:,:,ivar))

      if (nobs > 0) then
         typesdimlen = typesdimlen + 1

         call nc_check(nf90_put_att(ncid, NF90_GLOBAL, &
            trim(obs_type_strings(ivar)), ivar ), &
            'WriteNetCDF', 'region_names:obs_kinds')
      endif
   enddo

   if (typesdimlen < 1) then
      call nc_check(nf90_put_att(ncid, NF90_GLOBAL, 'comment', &
              'NO OBSERVATIONS. Check input time window &
              &against observation times.' ), &
              'WriteNetCDF', 'put_att empty file comment '//trim(fname))
      call nc_check(nf90_put_att(ncid, NF90_GLOBAL, 'continued', &
              'NO OBSERVATIONS. Expected if using an obs_seq.out'), &
              'WriteNetCDF', 'put_att empty file comment '//trim(fname))
   endif

   !----------------------------------------------------------------------------
   ! Define the dimensions
   !----------------------------------------------------------------------------

   call nc_check(nf90_def_dim(ncid=ncid, &
              name='time',   len = NF90_UNLIMITED,   dimid = TimeDimID), &
              'WriteNetCDF', 'time:def_dim '//trim(fname))
   call nc_check(nf90_def_dim(ncid=ncid, &
              name='bounds',   len = 2,  dimid = BoundsDimID), &
              'WriteNetCDF', 'bounds:def_dim '//trim(fname))

   call nc_check(nf90_def_dim(ncid=ncid, &
              name='copy', len = Ncopies,            dimid = CopyDimID), &
              'WriteNetCDF', 'copy:def_dim '//trim(fname))

   call nc_check(nf90_def_dim(ncid=ncid, &
              name='obstypes', len = max_obs_kinds,    dimid = TypesDimID), &
              'WriteNetCDF', 'types:def_dim '//trim(fname))

   call nc_check(nf90_def_dim(ncid=ncid, &
              name='region', len = Nregions,         dimid = RegionDimID), &
              'WriteNetCDF', 'region:def_dim '//trim(fname))

   call nc_check(nf90_def_dim(ncid=ncid, &
              name='stringlength', len = stringlength, dimid = StringDimID), &
              'WriteNetCDF', 'stringlength:def_dim '//trim(fname))

   if (create_rank_histogram) then
   call nc_check(nf90_def_dim(ncid=ncid, &
              name='rank_bins', len = ens_size+1,  dimid = RankDimID), &
              'WriteNetCDF', 'rank_bins:def_dim '//trim(fname))
   endif

   !----------------------------------------------------------------------------
   ! Define the types of derived quantities - aka - 'copies'
   !----------------------------------------------------------------------------

   call nc_check(nf90_def_var(ncid=ncid, name='copy', xtype=nf90_int, &
             dimids=CopyDimID, varid=CopyVarID), &
             'WriteNetCDF', 'copy:def_var')
   call nc_check(nf90_put_att(ncid, CopyVarID, 'explanation', 'see CopyMetaData'), &
             'WriteNetCDF', 'copy:explanation')

   !----------------------------------------------------------------------------
   ! Define the observation types - needed to be a coordinate variable
   !----------------------------------------------------------------------------

   call nc_check(nf90_def_var(ncid=ncid, name='obstypes', xtype=nf90_int, &
             dimids=TypesDimID, varid=TypesVarID), &
             'WriteNetCDF', 'types:def_var')
   call nc_check(nf90_put_att(ncid, TypesVarID, 'explanation', 'see ObservationTypes'), &
             'WriteNetCDF', 'types:explanation')

   !----------------------------------------------------------------------------
   ! Define the regions coordinate variable and attributes
   !----------------------------------------------------------------------------

   call nc_check(nf90_def_var(ncid=ncid, name='region', xtype=nf90_int, &
             dimids=RegionDimID, varid=RegionVarID), 'WriteNetCDF', 'region:def_var')
   call nc_check(nf90_put_att(ncid, RegionVarID, 'long_name', 'model region'), &
             'WriteNetCDF', 'region:long_name')
   call nc_check(nf90_put_att(ncid, RegionVarID, 'units',     'nondimensional'), &
             'WriteNetCDF', 'region:units')
   call nc_check(nf90_put_att(ncid, RegionVarID, 'valid_range', (/1,Nregions/)), &
             'WriteNetCDF', 'region:valid_range')
   call nc_check(nf90_put_att(ncid, RegionVarID, 'explanation', 'see region_names'), &
             'WriteNetCDF', 'types:explanation')

   !----------------------------------------------------------------------------
   ! Define 'bounds' dimension
   !----------------------------------------------------------------------------

   call nc_check(nf90_def_var(ncid=ncid, name='bounds', xtype=nf90_int, &
             dimids=BoundsDimID, varid=BoundsVarID), 'WriteNetCDF', 'bounds:def_var')
   call nc_check(nf90_put_att(ncid, BoundsVarID, 'valid_range', (/1,2/)), &
             'WriteNetCDF', 'bounds:valid_range')

   !----------------------------------------------------------------------------
   ! Define the time coordinate variable and attributes
   !----------------------------------------------------------------------------

   call nc_check(nf90_def_var(ncid=ncid, name='time', xtype=nf90_double, &
             dimids=TimeDimID, varid=TimeVarID), 'WriteNetCDF', 'time:def_var')
   call nc_check(nf90_put_att(ncid, TimeVarID, 'standard_name',    'time'), &
             'WriteNetCDF', 'time:standard_name')
   call nc_check(nf90_put_att(ncid, TimeVarID, 'long_name', 'temporal bin midpoints'), &
             'WriteNetCDF', 'time:long_name')
   call nc_check(nf90_put_att(ncid, TimeVarID, 'units',     'days since 1601-01-01'), &
             'WriteNetCDF', 'time:units')
   call nc_check(nf90_put_att(ncid, TimeVarID, 'calendar',    'Gregorian'), &
             'WriteNetCDF', 'time:calendar')
   call nc_check(nf90_put_att(ncid, TimeVarID, 'axis',    'T'), &
             'WriteNetCDF', 'time:axis')
   call nc_check(nf90_put_att(ncid, TimeVarID, 'bounds',    'time_bounds'), &
             'WriteNetCDF', 'time:bounds')
   call nc_check(nf90_put_att(ncid, TimeVarID, 'valid_range', &
             (/ epochcenter(1), epochcenter(Nepochs) /)), &
             'WriteNetCDF', 'time:valid_range')

   !----------------------------------------------------------------------------
   ! Define the time edges coordinate variable and attributes
   !----------------------------------------------------------------------------

   call nc_check(nf90_def_var(ncid=ncid, name='time_bounds', xtype=nf90_double, &
             dimids=(/ BoundsDimID, TimeDimID/), varid=TimeBoundsVarID), &
               'WriteNetCDF', 'time_bounds:def_var')
   call nc_check(nf90_put_att(ncid, TimeBoundsVarID, 'long_name', 'temporal bin edges'), &
             'WriteNetCDF', 'time_bounds:long_name')
   call nc_check(nf90_put_att(ncid, TimeBoundsVarID, 'units',     'days since 1601-01-01'), &
             'WriteNetCDF', 'time_bounds:units')
   call nc_check(nf90_put_att(ncid, TimeBoundsVarID, 'calendar',  'Gregorian'), &
             'WriteNetCDF', 'time_bounds:calendar')
   call nc_check(nf90_put_att(ncid, TimeBoundsVarID, 'valid_range', &
             (/ epochedges(1,1), epochedges(2,Nepochs) /)), &
             'WriteNetCDF', 'time_bounds:valid_range')

   !----------------------------------------------------------------------------
   ! Define the unusual coordinate variables
   !----------------------------------------------------------------------------

   call nc_check(nf90_def_var(ncid=ncid, name='region_names', xtype=nf90_char, &
             dimids=(/ StringDimID, RegionDimID /), varid=RegionNamesVarID), &
             'WriteNetCDF', 'region:def_var')
   call nc_check(nf90_put_att(ncid, RegionNamesVarID, 'long_name', 'region names'), &
             'WriteNetCDF', 'region:long_name')

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

   if ( create_rank_histogram ) then
   call nc_check(nf90_def_var(ncid=ncid, name='rank_bins', xtype=nf90_int, &
             dimids=(/ RankDimID /), varid=RankVarID), &
             'WriteNetCDF', 'rank_bins:def_var')
   call nc_check(nf90_put_att(ncid, RankVarID, 'long_name', 'rank histogram bins'), &
             'WriteNetCDF', 'rank_bins:long_name')
   call nc_check(nf90_put_att(ncid, RankVarID, 'comment', &
         'position of the observation among the sorted noisy ensemble members'), &
             'WriteNetCDF', 'rank_bins:comment')
   endif

   !----------------------------------------------------------------------------
   ! Set nofill mode - supposed to be performance gain
   !----------------------------------------------------------------------------

   call nc_check(nf90_set_fill(ncid, NF90_NOFILL, i),  &
            'obs_diag:WriteNetCDF', 'set_nofill '//trim(fname))

   !----------------------------------------------------------------------------
   ! Leave define mode so we can fill
   !----------------------------------------------------------------------------

   call nc_check(nf90_enddef(ncid), 'WriteNetCDF', 'enddef '//trim(fname))

   !----------------------------------------------------------------------------
   ! Fill the coordinate variables.
   ! The time variable is filled as time progresses.
   !----------------------------------------------------------------------------

   call nc_check(nf90_put_var(ncid, CopyVarId, (/ (i,i=1,Ncopies) /) ), &
              'WriteNetCDF', 'copy:put_var')

   call nc_check(nf90_put_var(ncid, CopyMetaVarID, copy_names), &
              'WriteNetCDF', 'copymeta:put_var')

   call nc_check(nf90_put_var(ncid, TypesVarId, (/ (i,i=1,max_obs_kinds) /) ), &
              'WriteNetCDF', 'types:put_var')

   call nc_check(nf90_put_var(ncid, TypesMetaVarID, obs_type_strings(1:max_obs_kinds)), &
              'WriteNetCDF', 'typesmeta:put_var')

   call nc_check(nf90_put_var(ncid, RegionVarID, (/ (i,i=1,Nregions) /) ), &
              'WriteNetCDF', 'region:put_var')

   call nc_check(nf90_put_var(ncid, BoundsVarID, (/ 1, 2 /)), &
              'WriteNetCDF', 'bounds:put_var')

   call nc_check(nf90_put_var(ncid, TimeVarID, epochcenter ), &
              'WriteNetCDF', 'time:put_var')

   call nc_check(nf90_put_var(ncid, TimeBoundsVarID, epochedges ), &
              'WriteNetCDF', 'time_bounds:put_var')

   call nc_check(nf90_put_var(ncid, RegionNamesVarID, reg_names(1:Nregions)), &
             'WriteNetCDF', 'region_names:put_var')

   call nc_check(nf90_sync( ncid), 'WriteNetCDF', 'sync '//trim(fname))  

   !----------------------------------------------------------------------------
   ! write the data we took such pains to collate ...
   !----------------------------------------------------------------------------

   if (verbose) write(*,*)'summary for Priors of time-region vars' 
   if ( create_rank_histogram ) then
      ierr = WriteTRV(ncid, guess, TimeDimID, CopyDimID, RegionDimID, RankDimID)
   else
      ierr = WriteTRV(ncid, guess, TimeDimID, CopyDimID, RegionDimID)
   endif
   if (verbose) write(*,*) ! a little whitespace
   if (verbose) write(*,*)'summary for Posteriors of time-region vars' 
   ierr = WriteTRV(ncid, analy,    TimeDimID, CopyDimID, RegionDimID)
   if (verbose) write(*,*) ! a little whitespace

   !----------------------------------------------------------------------------
   ! finish ...
   !----------------------------------------------------------------------------

   call nc_check(nf90_sync( ncid), 'WriteNetCDF', 'sync '//trim(fname))  
   call nc_check(nf90_close(ncid), 'init_diag_output', 'close '//trim(fname))  

   end Subroutine WriteNetCDF




   Function WriteTRV(ncid, vrbl, TimeDimID, CopyDimID, RegionDimID, RankDimID)
   integer,           intent(in) :: ncid
   type(TRV_type),    intent(in) :: vrbl
   integer,           intent(in) :: TimeDimID, CopyDimID, RegionDimID
   integer, optional, intent(in) :: RankDimID
   integer :: WriteTRV

   integer :: nobs, ivar, itime, iregion
   integer :: Nbins, irank, ndata
   character(len=NF90_MAX_NAME) :: string1

   integer :: VarID, VarID2, oldmode
   real(r4), allocatable, dimension(:,:,:) :: rchunk
   integer,  allocatable, dimension(:,:,:) :: ichunk

   FLAVORS : do ivar = 1,num_obs_types

      nobs = sum(vrbl%Nposs(:,:,ivar))
      if (nobs < 1) cycle FLAVORS

      if (verbose) then
         write(*,'(i4,1x,A,1x,i8)') ivar, obs_type_strings(ivar), nobs
      endif

      allocate(rchunk(Nregions,Ncopies,Nepochs))
      rchunk = MISSING_R4

      do itime   = 1,Nepochs
      do iregion = 1,Nregions

         rchunk(iregion, 1,itime) = vrbl%Nposs(      itime,iregion,ivar)
         rchunk(iregion, 2,itime) = vrbl%Nused(      itime,iregion,ivar)
         rchunk(iregion, 3,itime) = vrbl%rmse(       itime,iregion,ivar)
         rchunk(iregion, 4,itime) = vrbl%bias(       itime,iregion,ivar)
         rchunk(iregion, 5,itime) = vrbl%spread(     itime,iregion,ivar)
         rchunk(iregion, 6,itime) = vrbl%totspread(  itime,iregion,ivar)
         rchunk(iregion, 7,itime) = vrbl%NbadDartQC( itime,iregion,ivar)
         rchunk(iregion, 8,itime) = vrbl%observation(itime,iregion,ivar)
         rchunk(iregion, 9,itime) = vrbl%ens_mean(   itime,iregion,ivar)
         rchunk(iregion,10,itime) = vrbl%NDartQC_0(  itime,iregion,ivar)
         rchunk(iregion,11,itime) = vrbl%NDartQC_1(  itime,iregion,ivar)
         rchunk(iregion,12,itime) = vrbl%NDartQC_2(  itime,iregion,ivar)
         rchunk(iregion,13,itime) = vrbl%NDartQC_3(  itime,iregion,ivar)
         rchunk(iregion,14,itime) = vrbl%NDartQC_4(  itime,iregion,ivar)
         rchunk(iregion,15,itime) = vrbl%NDartQC_5(  itime,iregion,ivar)
         rchunk(iregion,16,itime) = vrbl%NDartQC_6(  itime,iregion,ivar)
         rchunk(iregion,17,itime) = vrbl%NDartQC_7(  itime,iregion,ivar)
         rchunk(iregion,18,itime) = vrbl%Ntrusted(   itime,iregion,ivar)

      enddo
      enddo

      call nc_check(nf90_redef(ncid), 'WriteTRV', 'redef')  

      ! Create netCDF variable name
      
      str1 = obs_type_strings(ivar)
      string1 = trim(str1)//'_'//adjustl(vrbl%string)

      call nc_check(nf90_def_var(ncid, name=string1, xtype=nf90_real, &
             dimids=(/ RegionDimID, CopyDimID, TimeDimID /), &
             varid=VarID), 'WriteTRV', 'region:def_var')
      call nc_check(nf90_put_att(ncid, VarID, '_FillValue',    MISSING_R4), &
              'WriteTRV','put_att:fillvalue')
      call nc_check(nf90_put_att(ncid, VarID, 'missing_value', MISSING_R4), &
              'WriteTRV','put_att:missing')

      call nc_check(nf90_set_fill(ncid, NF90_NOFILL, oldmode),  &
              'WriteTRV', 'set_nofill '//trim(vrbl%string))

      ! The rank histogram has no 'copy' dimension, so it must be handled differently.

      ndata = 0
      if (present(RankDimID)) then

         string1 = trim(string1)//'_RankHist'
         Nbins   = size(vrbl%hist_bin,4)
         ndata   = sum(vrbl%hist_bin(:,:,ivar,:))

         if ( ndata > 0 ) then

            allocate(ichunk(Nregions,Nbins,Nepochs))
            ichunk = 0

            do itime   = 1,Nepochs
            do iregion = 1,Nregions
            do irank   = 1,Nbins
   
            ichunk(iregion,irank,itime) = vrbl%hist_bin(itime,iregion,ivar,irank)

            enddo
            enddo
            enddo

            call nc_check(nf90_def_var(ncid, name=string1, xtype=nf90_int, &
                dimids=(/ RegionDimID, RankDimID, TimeDimID /), &
                varid=VarID2), 'WriteTRV', 'rank_hist:def_var')
         else
            write(logfileunit,*)string1//' has ',ndata,'"rank"able observations.'
            write(     *     ,*)string1//' has ',ndata,'"rank"able observations.'
         endif

      endif

      call nc_check(nf90_enddef(ncid), 'WriteTRV', 'enddef ')

      call nc_check(nf90_put_var(ncid, VarID, rchunk ), &
              'WriteTRV', 'realchunk:put_var')
      deallocate(rchunk)

      if (present(RankDimID) .and. (ndata > 0) ) then
         call nc_check(nf90_put_var(ncid, VarID2, ichunk ), &
                 'WriteTRV', 'intchunk:put_var')
         deallocate(ichunk)
      endif

   enddo FLAVORS

   WriteTRV = 0

   end Function WriteTRV



!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~



   Function DefineTrustedObs()

   ! Count up the number of 'trusted' observations
   ! Check to make sure trusted observations desired are supported.
   ! The namelist specifies 'trusted_obs(:)' ... we create a local list that
   ! is the intersection of the trusted_obs desired and those on hand ...

   integer :: DefineTrustedObs

   ! character(len=*), intent(out) :: trusted_obsname    the list of trusted observation types
   ! integer,          intent(out) :: DefineTrustedObs   the number of trusted observations types

   integer :: i, num_trusted, ikind
   logical :: matched
   character(len=NF90_MAX_NAME) :: string1

   num_trusted = 0
   
   CountTrusted : do i = 1,MaxTrusted
      if (trim(trusted_obs(i)) == 'null') exit CountTrusted
   
      matched = .false.
      VerifyTrusted : do ikind = 1,max_obs_kinds
         if (trim(trusted_obs(i)) == trim( get_obs_kind_name(ikind) )) then
            matched = .true.
            exit VerifyTrusted
         endif
      enddo VerifyTrusted
   
      if (matched) then
          num_trusted = num_trusted + 1
          trusted_obsname(num_trusted) = trim(trusted_obs(i))
      else
         write(string1,*)'trusted_obs "',trim(trusted_obs(i)),'" is not a supported observation type.'
         call error_handler(E_WARN, 'DefineTrustedObs', trim(string1), source, revision, revdate)
      endif
   enddo CountTrusted
   
   if (num_trusted == MaxTrusted) then
      write(string1,*)'There are ',num_trusted,' "trusted" observation types.'
      call error_handler(E_WARN, 'DefineTrustedObs', string1, &
           text2='This is the maximum allowed unless you increase "MaxTrusted" and recompile.')
   endif
   
   if (num_trusted > 0 ) then
      write(string1,*)'There are ',num_trusted,' "trusted" observation types, they are:'
      call error_handler(E_MSG, 'DefineTrustedObs', string1)
      do i = 1,num_trusted
         call error_handler(E_MSG, 'DefineTrustedObs', trim(trusted_obsname(i)) )
      enddo
   else
      write(string1,*)'There are no "trusted" observation types.'
      call error_handler(E_MSG, 'DefineTrustedObs', string1)
   endif

   DefineTrustedObs = num_trusted

   end Function DefineTrustedObs


!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


   Function is_observation_trusted(obsname)

   ! Is the observation one that we 'trust'.
   ! If so, disregard the DART QC ==7 (outlier rejection) and use it to
   ! in the statistics calculations.
   ! Since each obs_sequence file can have its own header/table, the safest
   ! way is to compare the string to a list of trusted (string) observation types.

   character(len=*), intent(in) :: obsname
   logical                      :: is_observation_trusted

   is_observation_trusted = .false.
   rUtrusted : do i = 1,num_trusted
      if ( trim(obsname) == trim(trusted_obsname(i)) ) then
         is_observation_trusted = .true.
         exit rUtrusted
      endif
   enddo rUtrusted

   end Function is_observation_trusted



!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

   Subroutine count_QC_values(iqc,itime,ireg,iflav)

   integer, intent(in) :: iqc, itime, ireg, iflav

   if (        iqc == 0 ) then
      call IPE(guess%NDartQC_0(itime,ireg,iflav), 1)
      call IPE(analy%NDartQC_0(itime,ireg,iflav), 1)

   elseif (    iqc == 1 ) then
      call IPE(guess%NDartQC_1(itime,ireg,iflav), 1)
      call IPE(analy%NDartQC_1(itime,ireg,iflav), 1)

   elseif (    iqc == 2 ) then
      call IPE(guess%NDartQC_2(itime,ireg,iflav), 1)
      call IPE(analy%NDartQC_2(itime,ireg,iflav), 1)

   elseif (    iqc == 3 ) then
      call IPE(guess%NDartQC_3(itime,ireg,iflav), 1)
      call IPE(analy%NDartQC_3(itime,ireg,iflav), 1)

   elseif (    iqc == 4 ) then
      call IPE(guess%NDartQC_4(itime,ireg,iflav), 1)
      call IPE(analy%NDartQC_4(itime,ireg,iflav), 1)

   elseif (    iqc == 5 ) then
      call IPE(guess%NDartQC_5(itime,ireg,iflav), 1)
      call IPE(analy%NDartQC_5(itime,ireg,iflav), 1)

   elseif (    iqc == 6 ) then
      call IPE(guess%NDartQC_6(itime,ireg,iflav), 1)
      call IPE(analy%NDartQC_6(itime,ireg,iflav), 1)

   elseif (    iqc == 7 ) then
      call IPE(guess%NDartQC_7(itime,ireg,iflav), 1)
      call IPE(analy%NDartQC_7(itime,ireg,iflav), 1)

   endif

   end Subroutine count_QC_values



!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


Subroutine NormalizeTRV()

integer :: ivar, iregion, iepoch

if (verbose) then
   write(logfileunit,*)'Normalizing time-region-variable quantities.'
   write(     *     ,*)'Normalizing time-region-variable quantities.'
endif

do ivar   = 1,num_obs_types
do iregion= 1,Nregions
do iepoch = 1,Nepochs

   ! The Priors

   if (  guess%Nused(      iepoch, iregion, ivar) == 0) then
         guess%observation(iepoch, iregion, ivar) = MISSING_R4
         guess%ens_mean(   iepoch, iregion, ivar) = MISSING_R4
         guess%bias(       iepoch, iregion, ivar) = MISSING_R4
         guess%rmse(       iepoch, iregion, ivar) = MISSING_R4
         guess%spread(     iepoch, iregion, ivar) = MISSING_R4
         guess%totspread(  iepoch, iregion, ivar) = MISSING_R4
   else
         guess%observation(iepoch, iregion, ivar) = &
         guess%observation(iepoch, iregion, ivar) / &
         guess%Nused(      iepoch, iregion, ivar)

         guess%ens_mean(   iepoch, iregion, ivar) = &
         guess%ens_mean(   iepoch, iregion, ivar) / &
         guess%Nused(      iepoch, iregion, ivar)

         guess%bias(       iepoch, iregion, ivar) = &
         guess%bias(       iepoch, iregion, ivar) / &
         guess%Nused(      iepoch, iregion, ivar)

         guess%rmse(       iepoch, iregion, ivar) = &
    sqrt(guess%rmse(       iepoch, iregion, ivar) / &
         guess%Nused(      iepoch, iregion, ivar) )

         guess%spread(     iepoch, iregion, ivar) = &
    sqrt(guess%spread(     iepoch, iregion, ivar) / &
         guess%Nused(      iepoch, iregion, ivar) )

         guess%totspread(  iepoch, iregion, ivar) = &
    sqrt(guess%totspread(  iepoch, iregion, ivar) / &
         guess%Nused(      iepoch, iregion, ivar) )

   endif

   ! The Posteriors

   if (  analy%Nused(      iepoch, iregion, ivar) == 0) then
         analy%observation(iepoch, iregion, ivar) = MISSING_R4
         analy%ens_mean(   iepoch, iregion, ivar) = MISSING_R4
         analy%bias(       iepoch, iregion, ivar) = MISSING_R4
         analy%rmse(       iepoch, iregion, ivar) = MISSING_R4
         analy%spread(     iepoch, iregion, ivar) = MISSING_R4
         analy%totspread(  iepoch, iregion, ivar) = MISSING_R4
   else
         analy%observation(iepoch, iregion, ivar) = &
         analy%observation(iepoch, iregion, ivar) / &
         analy%Nused(      iepoch, iregion, ivar)

         analy%ens_mean(   iepoch, iregion, ivar) = &
         analy%ens_mean(   iepoch, iregion, ivar) / &
         analy%Nused(      iepoch, iregion, ivar)

         analy%bias(       iepoch, iregion, ivar) = &
         analy%bias(       iepoch, iregion, ivar) / &
         analy%Nused(      iepoch, iregion, ivar)

         analy%rmse(       iepoch, iregion, ivar) = &
    sqrt(analy%rmse(       iepoch, iregion, ivar) / &
         analy%Nused(      iepoch, iregion, ivar) )

         analy%spread(     iepoch, iregion, ivar) = &
    sqrt(analy%spread(     iepoch, iregion, ivar) / &
         analy%Nused(      iepoch, iregion, ivar) )

         analy%totspread(  iepoch, iregion, ivar) = &
    sqrt(analy%totspread(  iepoch, iregion, ivar) / &
         analy%Nused(      iepoch, iregion, ivar) )

   endif
enddo
enddo
enddo

end Subroutine NormalizeTRV


!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

end program obs_diag

! <next few lines under version control, do not edit>
! $URL$
! $Id$
! $Revision$
! $Date$
