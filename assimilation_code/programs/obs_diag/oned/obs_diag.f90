! DART software - Copyright UCAR. This open source software is provided
! by UCAR, "as is", without charge, subject to all terms of use at
! http://www.image.ucar.edu/DAReS/DART/DART_download

!> The programs defines a series of epochs (periods of time) and geographic
!> regions and accumulates statistics for these epochs and regions.
!> All 'possible' observation types are treated separately.
!> The results are written to a netCDF file.
!> If the rank histogram is requested (and if the data is available),
!> only the PRIOR rank is calculated.

program obs_diag

! In Atmospheric Science, 'spread' has units of standard deviation ...
! In filter:obs_space_diagnostics() the 'spread' copies are converted to
! standard deviations.
!
! I should rename some of the variables I use as variances to reflect this.
! 'priorspred' should really be 'priorvar' since you have to accumulate variances
! the math is correct as it is, but the variable names don't make it easy ...

use        types_mod, only : r4, r8, digits12, MISSING_R4, MISSING_R8, &
                             metadatalength

use obs_sequence_mod, only : read_obs_seq, obs_type, obs_sequence_type, get_first_obs, &
                             get_obs_from_key, get_obs_def, get_copy_meta_data, &
                             get_obs_time_range, get_time_range_keys, &
                             get_num_times, get_obs_values, init_obs, &
                             assignment(=), get_num_copies, static_init_obs_sequence, &
                             get_qc, destroy_obs_sequence, get_last_obs, get_num_qc, &
                             read_obs_seq_header, destroy_obs, get_qc_meta_data

use      obs_def_mod, only : obs_def_type, get_obs_def_error_variance, get_obs_def_time, &
                             get_obs_def_location, get_obs_def_type_of_obs

use     obs_kind_mod, only : max_defined_types_of_obs, get_name_for_type_of_obs, &
                             RAW_STATE_VARIABLE

use     location_mod, only : location_type, get_location, operator(/=), LocationDims
use time_manager_mod, only : time_type, set_time, get_time, print_time, &
                             print_date, set_calendar_type, get_date, &
                             operator(*), operator(+), operator(-), &
                             operator(>), operator(<), operator(/), &
                             operator(/=), operator(<=), operator(>=)
use    utilities_mod, only : open_file, close_file, &
                             file_exist, error_handler, E_ERR, E_WARN, E_MSG,  &
                             initialize_utilities, logfileunit, nmlfileunit,   &
                             find_namelist_in_file, check_namelist_read,       &
                             do_nml_file, do_nml_term, finalize_utilities,     &
                             set_filename_list

use netcdf_utilities_mod, only : nc_check
                             
use         sort_mod, only : sort

use   random_seq_mod, only : random_seq_type, init_random_seq, &
                             several_random_gaussians

use typeSizes
use netcdf

implicit none

character(len=*), parameter :: source = 'oned/obs_diag.f90'

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

character(len=stringlength), dimension(MaxTrusted) :: trusted_list = 'null'

! Storage with fixed size for observation space diagnostics
real(r8), dimension(1) :: prior_mean, posterior_mean, prior_spread, posterior_spread
real(r8) :: pr_mean, po_mean ! same as above, without useless dimension
real(r8) :: pr_sprd, po_sprd ! same as above, without useless dimension

integer :: obs_index, prior_mean_index, posterior_mean_index
integer :: prior_spread_index, posterior_spread_index
integer :: flavor
integer :: num_copies, num_qc, num_obs, max_num_obs, obs_seq_file_id
integer :: num_obs_types

! variables used primarily/exclusively for the rank histogram
integer               :: ens_size, rank_histogram_bin
type(random_seq_type) :: ran_seq
real(r8)              :: obs_error_variance

character(len=stringlength) :: obs_seq_read_format
logical :: pre_I_format

integer,  dimension(2) :: key_bounds
real(r8), dimension(1) :: obs

integer,  allocatable, dimension(:) :: keys
integer,  allocatable, dimension(:) :: ens_copy_index

logical :: out_of_range, is_there_one, keeper

! Filter has an option to compute the posterior values or not; consequently
! the observation sequences may not have posterior values at all. This
! has a profound impact on the logic of 'obs_diag'
logical :: has_posteriors = .true.

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
! 5     not assimilated or evaluated because of namelist specification of
!         input.nml:obs_kind_nml:[assimilate,evaluate]_these_obs_types
! 6     prior QC rejected
! 7     outlier rejected
! 8     failed vertical conversion
! 9+    reserved for future use
!
! Some DART QC == 4 have meaningful posterior mean/spread (i.e. not MISSING)
! Anything with a DART QC == 5 has MISSING values for all DART copies
! Anything with a DART QC == 6 has MISSING values for all DART copies
! Anything with a DART QC == 7 has 'good' values for all DART copies, EXCEPT
! ambiguous case:
! prior rejected (7) ... posterior fails (should be 7 & 2)
!
! FIXME can there be a case where the prior is evaluated and the posterior QC is wrong
! FIXME ... there are cases where the prior fails but the posterior works ...

integer             :: org_qc_index, dart_qc_index, qc_value
integer, parameter  :: QC_MAX_PRIOR     = 3
integer, parameter  :: QC_MAX_POSTERIOR = 1
integer, parameter  :: QC_OUTLIER       = 7
integer, parameter  :: QC_PO_FOP_FAIL   = 2

real(r8), allocatable, dimension(:) :: qc
real(r8), allocatable, dimension(:) :: copyvals

integer, parameter, dimension(5) ::          hist_qcs = (/ 0, 1, 2, 3, 7 /)
integer, parameter, dimension(5) :: trusted_prior_qcs = (/ 0, 1, 2, 3, 7 /)
integer, parameter, dimension(3) :: trusted_poste_qcs = (/ 0, 1,       7 /)
integer, parameter, dimension(4) ::    good_prior_qcs = (/ 0, 1, 2, 3 /)
integer, parameter, dimension(2) ::    good_poste_qcs = (/ 0, 1       /)
integer :: numqcvals

integer, parameter :: max_num_input_files = 10000

!-----------------------------------------------------------------------
! Namelist with default values

character(len=256) :: obs_sequence_name(max_num_input_files) = ''
character(len=256) :: obs_sequence_list = ''
integer :: bin_width_days     = -1   ! width of the assimilation bin - seconds
integer :: bin_width_seconds  = -1   ! width of the assimilation bin - days
integer :: init_skip_days     = 0
integer :: init_skip_seconds  = 0
integer :: max_num_bins       = 9999 ! maximum number of temporal bins to consider

! index 1 == region 1 == [0.0, 1.0) i.e. Entire domain
! index 2 == region 2 == [0.0, 0.5)
! index 3 == region 3 == [0.5, 1.0)

integer :: Nregions = MaxRegions
real(r8), dimension(MaxRegions) :: lonlim1 = (/ 0.0_r8, 0.0_r8, 0.5_r8, -1.0_r8 /)
real(r8), dimension(MaxRegions) :: lonlim2 = (/ 1.0_r8, 0.5_r8, 1.0_r8, -1.0_r8 /)
character(len=6), dimension(MaxRegions) :: reg_names = &
                                   (/ 'whole ','yin   ','yang  ','bogus '/)

character(len=stringlength), dimension(MaxTrusted) :: trusted_obs = 'null'

logical :: verbose               = .false.
logical :: outliers_in_histogram = .true.
logical :: create_rank_histogram = .true.
logical :: use_zero_error_obs    = .false.

namelist /obs_diag_nml/ obs_sequence_name, obs_sequence_list,  &
                        bin_width_days, bin_width_seconds,     &
                        init_skip_days, init_skip_seconds, max_num_bins, &
                        Nregions, lonlim1, lonlim2, reg_names, &
                        verbose, outliers_in_histogram,        &
                        create_rank_histogram, trusted_obs, use_zero_error_obs

!-----------------------------------------------------------------------
! Variables used to accumulate the statistics.
!-----------------------------------------------------------------------

!>@todo must be a more clever way to relate the copy_names to the components

integer, parameter :: Ncopies = 19
character(len=stringlength), dimension(Ncopies) :: copy_names =                  &
   (/ 'Nposs      ', 'Nused      ',                                              &
      'rmse       ', 'bias       ', 'spread     ', 'totalspread',                &
      'NbadDARTQC ', 'observation', 'ens_mean   ', 'N_trusted  ',                &
      'N_DARTqc_0 ', 'N_DARTqc_1 ', 'N_DARTqc_2 ', 'N_DARTqc_3 ', 'N_DARTqc_4 ', &
      'N_DARTqc_5 ', 'N_DARTqc_6 ', 'N_DARTqc_7 ', 'N_DARTqc_8 '                /)

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
   integer,  dimension(:,:,:), pointer :: NDartQC_8
   integer,  dimension(:,:,:,:), pointer :: hist_bin => NULL()
end type TRV_type

type(TRV_type) :: prior, poste

type(time_type), allocatable, dimension(:)   :: bincenter
type(time_type), allocatable, dimension(:,:) :: binedges
real(digits12),  allocatable, dimension(:)   :: epoch_center
real(digits12),  allocatable, dimension(:,:) :: epoch_edges
integer,         allocatable, dimension(:)   :: obs_used_in_epoch

!-----------------------------------------------------------------------
! General purpose variables
!-----------------------------------------------------------------------

integer  :: iregion, iepoch, ivar, ifile, num_obs_in_epoch
real(r8) :: rlocation

integer  :: obsindex, i, iunit, ierr, io, ireg
integer  :: seconds, days, Nepochs, num_input_files

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

character(len=512) :: string1, string2, string3
character(len=stringlength) :: obsname

integer :: Nidentity = 0
integer :: num_ambiguous = 0   ! prior QC 7, posterior mean MISSING_R8

!=======================================================================
! Get the party started
!=======================================================================

call initialize_utilities('obs_diag')
call static_init_obs_sequence()

num_obs_types = max_defined_types_of_obs ! for compatibility with 3D version
allocate(obs_type_strings(num_obs_types))
do ivar = 1,max_defined_types_of_obs
   obs_type_strings(ivar) = get_name_for_type_of_obs(ivar)
enddo

! Read the namelist

call find_namelist_in_file('input.nml', 'obs_diag_nml', iunit)
read(iunit, nml = obs_diag_nml, iostat = io)
call check_namelist_read(iunit, io, 'obs_diag_nml')

! Record the namelist values used for the run ...
if (do_nml_file()) write(nmlfileunit, nml=obs_diag_nml)
if (do_nml_term()) write(    *      , nml=obs_diag_nml)

num_input_files = set_filename_list(obs_sequence_name, obs_sequence_list, 'obs_diag')

num_trusted = DefineTrustedObs()

! Check to see if we are including the outlier observations in the
! rank histogram calculation.

if ( outliers_in_histogram ) then
   numqcvals = size(hist_qcs)
else
   numqcvals = size(hist_qcs) - 1
endif

! Can lie about calendar for low-order models since its all hypothetical

call set_calendar_type('GREGORIAN')

! Determine temporal bin characteristics.
! Nepochs is the total number of time intervals of the period requested.
! if the namelist does not specify a start/stop time and binwidth, we will
! presume the first/last times in the input file(s) are to be used.

call DetermineNumEpochs(obsT1, obsTN, Nepochs, skip_time)
call DefineTemporalBins()

allocate(   bincenter(Nepochs),    binedges(2,Nepochs)) ! time_type
allocate(epoch_center(Nepochs), epoch_edges(2,Nepochs)) ! 64bit reals for netCDF
allocate( obs_used_in_epoch(Nepochs) )

call SetSchedule(obsT1, Nepochs, binwidth, halfbinwidth, &
                 bincenter, binedges, epoch_center, epoch_edges)

TimeMin = binedges(1,      1) ! minimum time of interest
TimeMax = binedges(2,Nepochs) ! maximum time of interest
obs_used_in_epoch = 0

! Rectify the region namelist information

ireg = MaxRegions
Regions: do i = 1,MaxRegions
   if ((lonlim1(i) < 0.0_r8) .or. (lonlim2(i) < 0.0_r8) ) then
      exit Regions
   else
      ireg = i
   endif
enddo Regions
Nregions = min(Nregions, ireg)

if ( verbose ) then
   do i = 1,Nregions
      write(string1,'(''Region '',i02,1x,a32,'' : '',2(f10.4,1x))') &
             i, reg_names(i), lonlim1(i), lonlim2(i)
      call error_handler(E_MSG,'obs_diag',string1)
   enddo
endif

call InitializeVariables()

! Open file for histogram of innovations, as a function of standard deviation.

nsigmaUnit = open_file('LargeInnov.txt',form='formatted',action='write')
write(nsigmaUnit,'(a)')'Any observations flagged as bad are dumped into the last bin.'
write(nsigmaUnit,'(a)') '   day   secs    loc            obs         prior   zscore   key   kind'

!-----------------------------------------------------------------------

ObsFileLoop : do ifile=1, num_input_files

   write(string1,*)'Reading file # ',ifile, ' of ',num_input_files, &
                   ' "'//trim(obs_sequence_name(ifile))//'"'
   call error_handler(E_MSG,'obs_diag',string1)

   ! Read in information about observation sequence so we can allocate
   ! observations. We need info about how many copies, qc values, etc.
   ! We have already read this sequence once, so no caution required.

   call read_obs_seq_header(obs_sequence_name(ifile), &
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

   call read_obs_seq(obs_sequence_name(ifile), 0, 0, 0, seq)

   ! Determine the time encompassed in the observation sequence.
   ! also compare to first/last times of ALL sequences

   call GetFirstLastObs(obs_sequence_name(ifile), seq, obs1, obsN, seqT1, seqTN)

   if (No_Time_Intersection(obs_sequence_name(ifile),seqT1,seqTN,TimeMin,TimeMax)) then
      call destroy_obs(obs1)
      call destroy_obs(obsN)
      call destroy_obs(observation)
      call destroy_obs(next_obs)
      call destroy_obs_sequence(seq)
      if (allocated(qc)) deallocate( qc )
      if (allocated(copyvals)) deallocate( copyvals )
      cycle ObsFileLoop
   endif

   ! Prepare some variables for the rank histogram.
   ! FIXME : Make sure this observation sequence file has the same number of
   ! ensemble members as 'the first one' (which defines the bins) ...

   ens_size = GetEnsSize()

   if ((ens_size == 0) .and. create_rank_histogram) then
      call error_handler(E_MSG,'obs_diag', &
                 'Cannot create rank histogram. Zero ensemble members.')
      create_rank_histogram = .false.

   elseif ((ens_size > 0) .and. create_rank_histogram ) then
      if (allocated(ens_copy_index)) then
         if (size(ens_copy_index) /= ens_size) then
            write(string1,'(''expecting '',i3,'' ensemble members, got '',i3)') &
                                       size(ens_copy_index), ens_size
            call error_handler(E_ERR,'obs_diag',string1,source)
         endif
      else
         ! This should happen exactly once, if at all.
         allocate(prior%hist_bin( Nepochs, Nregions, num_obs_types, ens_size+1))
         allocate(ens_copy_index(ens_size))
         prior%hist_bin    = 0
         call init_random_seq(ran_seq, seed=23)
      endif
      if ( verbose ) then
         write(string1,*) 'Creating rank histogram with ',ens_size+1,' bins.'
         call error_handler(E_MSG,'obs_diag',string1)
      endif
   endif

   ! Find the index of obs, ensemble mean, spread ... etc.
   !
   ! Only require obs_index to be present; this allows the program
   ! to be run on obs_seq.[in,out] files which have no means or spreads.
   ! You can still plot obs count, incoming QC, obs values ...
   !
   ! Each observation sequence file can have its copies in any order.

   call SetIndices()

   ! Loop over all potential time periods ... the observation sequence
   ! files are not required to be in any particular order.

   EpochLoop : do iepoch = 1, Nepochs

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
            write(string1,*)' No observations in epoch ',iepoch,' cycling ...'
            call error_handler(E_MSG,'obs_diag',string1)
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

      ObservationLoop : do obsindex = 1, num_obs_in_epoch

         ! 'flavor' is from the 'master list' in the obs_kind_mod.f90

         call get_obs_from_key(seq, keys(obsindex), observation)
         call get_obs_def(observation, obs_def)

         flavor    = get_obs_def_type_of_obs(obs_def)

         ! Check to see if it is an identity observation.
         ! Redefine identity observations as flavor = RAW_STATE_VARIABLE
         !> Still have a problem determining what state type best relates
         !> to the observation kind - but it would allow us to
         !> do this for all models, regardless of dimensionality.

         if ( flavor < 0 ) then
            Nidentity = Nidentity + 1
            flavor = RAW_STATE_VARIABLE
         endif

         obsname   = get_name_for_type_of_obs(flavor)
         obs_time  = get_obs_def_time(obs_def)
         obs_loc   = get_obs_def_location(obs_def)
         rlocation = get_location(obs_loc)

         ! Check to make sure we are past the burn-in
         if (obs_time < skip_time) cycle ObservationLoop

         ! Check to see if this is a trusted observation
         if ( num_trusted > 0 ) then
            trusted = is_observation_trusted( obsname )
         else
            trusted = .false.
         endif

         if ( use_zero_error_obs ) then
            obs_error_variance = 0.0_r8
         else
            obs_error_variance = get_obs_def_error_variance(obs_def)
         endif
         ! retrieve observation prior and posterior means and spreads

         prior_mean(1)       = 0.0_r8
         posterior_mean(1)   = 0.0_r8
         prior_spread(1)     = 0.0_r8
         posterior_spread(1) = 0.0_r8

            call get_obs_values(observation,              obs,              obs_index)
         if (prior_mean_index > 0) &
            call get_obs_values(observation,       prior_mean,       prior_mean_index)
         if (posterior_mean_index > 0) &
            call get_obs_values(observation,   posterior_mean,   posterior_mean_index)
         if (prior_spread_index > 0) &
            call get_obs_values(observation,     prior_spread,     prior_spread_index)
         if (posterior_spread_index > 0) &
            call get_obs_values(observation, posterior_spread, posterior_spread_index)

         pr_mean =       prior_mean(1)
         po_mean =   posterior_mean(1)
         pr_sprd =     prior_spread(1)
         po_sprd = posterior_spread(1)

         call get_qc(observation, qc)

         if ( dart_qc_index > 0 ) then
            qc_value = qc(dart_qc_index)
         else
            ! If there is no dart_qc, this must be a case where we 
            ! are interested only in getting the location information.
            qc_value = 0
         endif

         ! There is an ambiguous case wherein the prior is rejected (DART QC == 7)
         ! and the posterior forward operator fails (DART QC == 2). In this case,
         ! the DART_QC only reflects the fact the prior was rejected - HOWEVER -
         ! the posterior mean,spread are set to MISSING.
         !
         ! If it is your intent to compare identical prior and posterior TRUSTED
         ! observations, then you should enable the following few lines of code.
         ! and realize that the number of observations rejected because of the
         ! outlier threshold will be wrong.
         !
         ! This is the only block of code you should need to change.

         if (qc_value == QC_OUTLIER .and. posterior_mean(1) == MISSING_R8) then
            write(string1,*)'WARNING ambiguous case for obs index ',obsindex
            string2 = 'obs failed outlier threshhold AND posterior operator failed.'
            string3 = 'Counting as a Prior QC == 7, Posterior QC == 2.'
            if (trusted) then
! COMMENT      string3 = 'WARNING changing DART QC from 7 to 2'
! COMMENT      qc_value = 2
            endif
            call error_handler(E_MSG,'obs_diag',string1,text2=string2,text3=string3)
            num_ambiguous = num_ambiguous + 1
         endif

         ! (DEBUG) Summary of observation knowledge at this point

         if ( .false. ) then
            write(*,*)
            write(*,*)'observation #,flavor ', obsindex, flavor
            write(*,*)'obs(1), qc ', obs(1), qc
            write(*,*)'obs_error_variance ', obs_error_variance
            call print_time(obs_time,'time is')
            write(*,*)'pr_mean, po_mean ', pr_mean, po_mean
            write(*,*)'pr_sprd, po_sprd ', pr_sprd, po_sprd
         endif

         ! update the histogram of the magnitude of the innovation,
         ! where each bin is a single standard deviation.
         ! This is a one-sided histogram.

         pr_zscore = InnovZscore(obs(1), pr_mean, pr_sprd, obs_error_variance, &
                                 qc_value, QC_MAX_PRIOR)

         if (has_posteriors) po_zscore = InnovZscore(obs(1), po_mean, po_sprd, &
                                    obs_error_variance, qc_value, QC_MAX_POSTERIOR)

         indx         = min(int(pr_zscore), MaxSigmaBins)
         nsigma(indx) = nsigma(indx) + 1

         ! Individual (valid) observations that are very far away get
         ! logged to a separate file.

         if( (pr_zscore > 3.0_r8) .and. (qc_value <= QC_MAX_PRIOR) ) then
            call get_time(obs_time,seconds,days)

            write(nsigmaUnit,'(i7,1x,i5,1x,f8.2,1x,2f13.2,f8.1,2i7)') &
                 days, seconds, rlocation, &
                 obs(1), pr_mean, pr_zscore, keys(obsindex), flavor
         endif

         ! At this point, the observation has passed all checks.

         obs_used_in_epoch(iepoch) = obs_used_in_epoch(iepoch) + 1

         ! If needed, calculate the rank histogram bin (once!) for
         ! this observation - even if the QC value is bad.

         if ( create_rank_histogram ) then
            call get_obs_values(observation, copyvals)
            rank_histogram_bin = Rank_Histogram(copyvals, obs_index, &
                 obs_error_variance)
         endif

         ! We have Nregions of interest.

         Areas : do iregion =1, Nregions

            keeper = is_location_in_region( rlocation, lonlim1(iregion), lonlim2(iregion) )
            if ( .not. keeper ) cycle Areas

            call count_QC_values(qc_value, iepoch, iregion, flavor)

            call Bin3D(qc_value, iepoch, iregion, flavor, trusted, obs(1), &
                obs_error_variance, pr_mean, pr_sprd, po_mean, po_sprd, rank_histogram_bin)

         enddo Areas

      enddo ObservationLoop

      deallocate(keys)

      if( verbose ) then
         write(string1,'(''num obs considered in epoch '',i4,'' = '',i8, &
                                  & '' out of '',i8,'' possible'')') &
                         iepoch, obs_used_in_epoch(iepoch), num_obs_in_epoch
         call error_handler(E_MSG,'obs_diag',string1)
         write(logfileunit,*)''
         write(     *     ,*)''
      endif

   enddo EpochLoop

   if ( verbose ) then
      write(string1,*)'Finished reading "',trim(obs_sequence_name(ifile))//'"'
      call error_handler(E_MSG,'obs_diag',string1)
   endif

   call destroy_obs(obs1)
   call destroy_obs(obsN)
   call destroy_obs(observation)
   call destroy_obs(next_obs)
   call destroy_obs_sequence(seq)
   if (allocated(qc))       deallocate( qc )
   if (allocated(copyvals)) deallocate( copyvals )

enddo ObsFileLoop

! We have read all possible files, and stuffed the observations into the
! appropriate bins. Time to normalize.

call NormalizeTRV()

if (sum(obs_used_in_epoch) == 0 ) then
   call error_handler(E_ERR,'obs_diag','All identity observations. Stopping.', source)
endif

! Print final summary.

write(*,*)
write(*,*) '# observations used  : ',sum(obs_used_in_epoch)
write(*,*) 'Count summary over all regions - obs may count for multiple regions:'
write(*,*) '# identity           : ',Nidentity
write(*,*) '# bad DART QC prior  : ',sum(prior%NbadDartQC)
if (has_posteriors) write(*,*) '# bad DART QC post   : ',sum(poste%NbadDartQC)
write(*,*) '# priorQC 7 postQC 2 : ',num_ambiguous
write(*,*)
write(*,*) '# trusted prior   : ',sum(prior%Ntrusted)
write(*,*) '# prior DART QC 0 : ',sum(prior%NDartQC_0)
write(*,*) '# prior DART QC 1 : ',sum(prior%NDartQC_1)
write(*,*) '# prior DART QC 2 : ',sum(prior%NDartQC_2)
write(*,*) '# prior DART QC 3 : ',sum(prior%NDartQC_3)
write(*,*) '# prior DART QC 4 : ',sum(prior%NDartQC_4)
write(*,*) '# prior DART QC 5 : ',sum(prior%NDartQC_5)
write(*,*) '# prior DART QC 6 : ',sum(prior%NDartQC_6)
write(*,*) '# prior DART QC 7 : ',sum(prior%NDartQC_7)
write(*,*) '# prior DART QC 8 : ',sum(prior%NDartQC_8)
write(*,*)

if (has_posteriors) then
   write(*,*) '# trusted poste   : ',sum(poste%Ntrusted)
   write(*,*) '# poste DART QC 0 : ',sum(poste%NDartQC_0)
   write(*,*) '# poste DART QC 1 : ',sum(poste%NDartQC_1)
   write(*,*) '# poste DART QC 2 : ',sum(poste%NDartQC_2)
   write(*,*) '# poste DART QC 3 : ',sum(poste%NDartQC_3)
   write(*,*) '# poste DART QC 4 : ',sum(poste%NDartQC_4)
   write(*,*) '# poste DART QC 5 : ',sum(poste%NDartQC_5)
   write(*,*) '# poste DART QC 6 : ',sum(poste%NDartQC_6)
   write(*,*) '# poste DART QC 7 : ',sum(poste%NDartQC_7)
   write(*,*) '# poste DART QC 8 : ',sum(poste%NDartQC_8)
   write(*,*)
endif

write(logfileunit,*)
write(logfileunit,*) '# observations used  : ',sum(obs_used_in_epoch)
write(logfileunit,*) 'Count summary over all regions - obs may count for multiple regions:'
write(logfileunit,*) '# identity           : ',Nidentity
write(logfileunit,*) '# bad DART QC prior  : ',sum(prior%NbadDartQC)
if (has_posteriors) write(logfileunit,*) '# bad DART QC post   : ',sum(poste%NbadDartQC)
write(logfileunit,*)
write(logfileunit,*) '# trusted prior   : ',sum(prior%Ntrusted)
write(logfileunit,*) '# prior DART QC 0 : ',sum(prior%NDartQC_0)
write(logfileunit,*) '# prior DART QC 1 : ',sum(prior%NDartQC_1)
write(logfileunit,*) '# prior DART QC 2 : ',sum(prior%NDartQC_2)
write(logfileunit,*) '# prior DART QC 3 : ',sum(prior%NDartQC_3)
write(logfileunit,*) '# prior DART QC 4 : ',sum(prior%NDartQC_4)
write(logfileunit,*) '# prior DART QC 5 : ',sum(prior%NDartQC_5)
write(logfileunit,*) '# prior DART QC 6 : ',sum(prior%NDartQC_6)
write(logfileunit,*) '# prior DART QC 7 : ',sum(prior%NDartQC_7)
write(logfileunit,*) '# prior DART QC 8 : ',sum(prior%NDartQC_8)
write(logfileunit,*)

if (has_posteriors) then
   write(logfileunit,*) '# trusted poste   : ',sum(poste%Ntrusted)
   write(logfileunit,*) '# poste DART QC 0 : ',sum(poste%NDartQC_0)
   write(logfileunit,*) '# poste DART QC 1 : ',sum(poste%NDartQC_1)
   write(logfileunit,*) '# poste DART QC 2 : ',sum(poste%NDartQC_2)
   write(logfileunit,*) '# poste DART QC 3 : ',sum(poste%NDartQC_3)
   write(logfileunit,*) '# poste DART QC 4 : ',sum(poste%NDartQC_4)
   write(logfileunit,*) '# poste DART QC 5 : ',sum(poste%NDartQC_5)
   write(logfileunit,*) '# poste DART QC 6 : ',sum(poste%NDartQC_6)
   write(logfileunit,*) '# poste DART QC 7 : ',sum(poste%NDartQC_7)
   write(logfileunit,*) '# poste DART QC 8 : ',sum(poste%NDartQC_8)
   write(logfileunit,*)
endif

! Print the histogram of innovations as a function of standard deviation.
if ( verbose ) then
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
endif

! Open netCDF output file

call WriteNetCDF('obs_diag_output.nc')

call close_file(nsigmaUnit)
call DestroyVariables()
call error_handler(E_MSG,'obs_diag','Finished successfully.')
call finalize_utilities()


CONTAINS

!======================================================================
! These routines use common variables from the scope of this file.
! If it's not in the argument list ... it's scoped within this file.
!======================================================================


subroutine InitializeVariables()

allocate(prior%rmse(       Nepochs, Nregions, num_obs_types), &
         prior%bias(       Nepochs, Nregions, num_obs_types), &
         prior%spread(     Nepochs, Nregions, num_obs_types), &
         prior%totspread(  Nepochs, Nregions, num_obs_types), &
         prior%observation(Nepochs, Nregions, num_obs_types), &
         prior%ens_mean(   Nepochs, Nregions, num_obs_types), &
         prior%Nposs(      Nepochs, Nregions, num_obs_types), &
         prior%Nused(      Nepochs, Nregions, num_obs_types), &
         prior%NbadDartQC( Nepochs, Nregions, num_obs_types), &
         prior%Ntrusted(   Nepochs, Nregions, num_obs_types), &
         prior%NDartQC_0(  Nepochs, Nregions, num_obs_types), &
         prior%NDartQC_1(  Nepochs, Nregions, num_obs_types), &
         prior%NDartQC_2(  Nepochs, Nregions, num_obs_types), &
         prior%NDartQC_3(  Nepochs, Nregions, num_obs_types), &
         prior%NDartQC_4(  Nepochs, Nregions, num_obs_types), &
         prior%NDartQC_5(  Nepochs, Nregions, num_obs_types), &
         prior%NDartQC_6(  Nepochs, Nregions, num_obs_types), &
         prior%NDartQC_7(  Nepochs, Nregions, num_obs_types), &
         prior%NDartQC_8(  Nepochs, Nregions, num_obs_types))

prior%rmse        = 0.0_r8
prior%bias        = 0.0_r8
prior%spread      = 0.0_r8
prior%totspread   = 0.0_r8
prior%observation = 0.0_r8
prior%ens_mean    = 0.0_r8
prior%Nposs       = 0
prior%Nused       = 0
prior%NbadDartQC  = 0
prior%Ntrusted    = 0
prior%NDartQC_0   = 0
prior%NDartQC_1   = 0
prior%NDartQC_2   = 0
prior%NDartQC_3   = 0
prior%NDartQC_4   = 0
prior%NDartQC_5   = 0
prior%NDartQC_6   = 0
prior%NDartQC_7   = 0
prior%NDartQC_8   = 0

prior%string        = 'guess'
prior%num_times     = Nepochs
prior%num_regions   = Nregions
prior%num_variables = num_obs_types

allocate(poste%rmse(       Nepochs, Nregions, num_obs_types), &
         poste%bias(       Nepochs, Nregions, num_obs_types), &
         poste%spread(     Nepochs, Nregions, num_obs_types), &
         poste%totspread(  Nepochs, Nregions, num_obs_types), &
         poste%observation(Nepochs, Nregions, num_obs_types), &
         poste%ens_mean(   Nepochs, Nregions, num_obs_types), &
         poste%Nposs(      Nepochs, Nregions, num_obs_types), &
         poste%Nused(      Nepochs, Nregions, num_obs_types), &
         poste%NbadDartQC( Nepochs, Nregions, num_obs_types), &
         poste%Ntrusted(   Nepochs, Nregions, num_obs_types), &
         poste%NDartQC_0(  Nepochs, Nregions, num_obs_types), &
         poste%NDartQC_1(  Nepochs, Nregions, num_obs_types), &
         poste%NDartQC_2(  Nepochs, Nregions, num_obs_types), &
         poste%NDartQC_3(  Nepochs, Nregions, num_obs_types), &
         poste%NDartQC_4(  Nepochs, Nregions, num_obs_types), &
         poste%NDartQC_5(  Nepochs, Nregions, num_obs_types), &
         poste%NDartQC_6(  Nepochs, Nregions, num_obs_types), &
         poste%NDartQC_7(  Nepochs, Nregions, num_obs_types), &
         poste%NDartQC_8(  Nepochs, Nregions, num_obs_types))

poste%rmse        = 0.0_r8
poste%bias        = 0.0_r8
poste%spread      = 0.0_r8
poste%totspread   = 0.0_r8
poste%observation = 0.0_r8
poste%ens_mean    = 0.0_r8
poste%Nposs       = 0
poste%Nused       = 0
poste%NbadDartQC  = 0
poste%Ntrusted    = 0
poste%NDartQC_0   = 0
poste%NDartQC_1   = 0
poste%NDartQC_2   = 0
poste%NDartQC_3   = 0
poste%NDartQC_4   = 0
poste%NDartQC_5   = 0
poste%NDartQC_6   = 0
poste%NDartQC_7   = 0
poste%NDartQC_8   = 0

poste%string        = 'analy'
poste%num_times     = Nepochs
poste%num_regions   = Nregions
poste%num_variables = num_obs_types

end subroutine InitializeVariables


!======================================================================


subroutine DestroyVariables()

if (associated(prior%hist_bin)) deallocate(prior%hist_bin)
if (allocated(ens_copy_index))  deallocate(ens_copy_index)

deallocate(prior%rmse,        prior%bias,      prior%spread,    prior%totspread, &
           prior%observation, prior%ens_mean,  prior%Nposs,     prior%Nused,     &
                                               prior%NbadDartQC,prior%Ntrusted,  &
           prior%NDartQC_0,   prior%NDartQC_1, prior%NDartQC_2, prior%NDartQC_3, &
           prior%NDartQC_4,   prior%NDartQC_5, prior%NDartQC_6, prior%NDartQC_7, &
           prior%NDartQC_8)

deallocate(poste%rmse,        poste%bias,      poste%spread,    poste%totspread, &
           poste%observation, poste%ens_mean,  poste%Nposs,     poste%Nused,     &
                                               poste%NbadDartQC,poste%Ntrusted,  &
           poste%NDartQC_0,   poste%NDartQC_1, poste%NDartQC_2, poste%NDartQC_3, &
           poste%NDartQC_4,   poste%NDartQC_5, poste%NDartQC_6, poste%NDartQC_7, &
           poste%NDartQC_8)

deallocate(epoch_center, epoch_edges, bincenter, binedges, obs_used_in_epoch)
deallocate(obs_type_strings)

end subroutine DestroyVariables


!======================================================================


function InnovZscore(obsval, prmean, prspred, errvar, qcval, qcmaxval)

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

end function InnovZscore


!======================================================================


function is_location_in_region( lon, lon1, lon2 ) result( keeper )

!>@todo support if the region of interest is [0.8, 0.2]
!       ... are all 1D locations periodic

real(r8), intent(in) :: lon, lon1, lon2
logical :: keeper

keeper = .false.

if( (lon .ge. lon1) .and. (lon .lt. lon2) ) keeper = .true.

end function is_location_in_region


!======================================================================
!> This routine reads through all the observation sequence files to determine
!> the first and last observation times and how many distinct times are spanned.
!> Can only do this by reading every observation sequence file.
!> Does not automatically discount identical times, so if identical times
!> are specified in multiple files, the bin width must be manually specified.
!> Alternatively, use the obs_sequence_tool to concatenate the input files and
!> feed the result to obs_diag.

subroutine DetermineNumEpochs( time1, timeN, nsteps, time_to_skip )

type(time_type), intent(out) :: time1    ! first observation time
type(time_type), intent(out) :: timeN    ! last observation time
integer,         intent(out) :: nsteps   ! number of distinct obs times in all files
type(time_type), intent(out) :: time_to_skip

character(len=*), parameter :: routine = 'DetermineNumEpochs'
integer :: seqNsteps

nsteps   = 0

do ifile = 1, num_input_files

   if ( file_exist(trim(obs_sequence_name(ifile))) ) then
      write(string1,*)'opening "'//trim(obs_sequence_name(ifile))//'"'
      call error_handler(E_MSG,routine,string1)
   else
      write(string1,*)'input observation file does not exist'
      write(string2,*)'looking for "'//trim(obs_sequence_name(ifile))//'"'
      call error_handler(E_ERR, routine, string1, source, text2=string2)
   endif

   ! Read in information about observation sequence so we can allocate
   ! observations. We need info about how many copies, qc values, etc.

   call read_obs_seq_header(obs_sequence_name(ifile), &
             num_copies, num_qc, num_obs, max_num_obs, &
             obs_seq_file_id, obs_seq_read_format, pre_I_format, &
             close_the_file = .true.)

   ! Initialize some (individual) observation variables
   ! Read in the entire observation sequence

   call init_obs( obs1, num_copies, num_qc)
   call init_obs( obsN, num_copies, num_qc)
   call read_obs_seq(obs_sequence_name(ifile), 0, 0, 0, seq)

   ! Determine the time encompassed in the observation sequence.

   !>@todo  replace with GetFirstLastObs()

   is_there_one = get_first_obs(seq, obs1)
   if ( .not. is_there_one ) then
      call error_handler(E_ERR,routine,'No first observation in sequence.', &
      source,text2=obs_sequence_name(ifile))
   endif
   call get_obs_def(obs1,   obs_def)
   seqT1 = get_obs_def_time(obs_def)

   is_there_one = get_last_obs(seq, obsN)
   if ( .not. is_there_one ) then
      call error_handler(E_ERR,routine,'No last observation in sequence.', &
      source,text2=obs_sequence_name(ifile))
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

      write(logfileunit,*)trim(obs_sequence_name(ifile)),' has ', &
                          seqNsteps,' distinct times.'
      write(    *      ,*)trim(obs_sequence_name(ifile)),' has ', &
                          seqNsteps,' distinct times.'

      call print_date(seqT1, ' DetermineNumEpochs: observation 1 date', logfileunit)
      call print_date(seqTN, ' DetermineNumEpochs: observation N date', logfileunit)
      call print_time(seqT1, ' DetermineNumEpochs: observation 1 time', logfileunit)
      call print_time(seqTN, ' DetermineNumEpochs: observation N time', logfileunit)

      call print_date(seqT1, ' DetermineNumEpochs: observation 1 date')
      call print_date(seqTN, ' DetermineNumEpochs: observation N date')
      call print_time(seqT1, ' DetermineNumEpochs: observation 1 time')
      call print_time(seqTN, ' DetermineNumEpochs: observation N time')

   endif

   call destroy_obs(obs1)
   call destroy_obs(obsN)
   call destroy_obs_sequence(seq)

enddo

if (nsteps < 1) then
   write(string1,*)'cannot find any times in the ',num_input_files,' input files.'
   call error_handler(E_ERR,'DetermineNumEpochs',string1,source)
endif

! Summarize first and last observations, first and last bins
write(logfileunit,*)
write(    *      ,*)

write(logfileunit,*)'Have ',nsteps,' distinct times.'
write(    *      ,*)'Have ',nsteps,' distinct times.'

call print_date(time1, ' first bincenter date',logfileunit)
call print_date(timeN, ' last  bincenter date',logfileunit)
call print_time(time1, ' first bincenter time',logfileunit)
call print_time(timeN, ' last  bincenter time',logfileunit)

call print_date(time1, ' first bincenter date')
call print_date(timeN, ' last  bincenter date')
call print_time(time1, ' first bincenter time')
call print_time(timeN, ' last  bincenter time')

write(logfileunit,*)
write(    *      ,*)

time_to_skip = set_time(init_skip_seconds, init_skip_days)

if (timeN < time_to_skip) then
   call print_date(time_to_skip, ' implied skip-to-date',logfileunit)
   call print_time(time_to_skip, ' implied skip-to-time',logfileunit)
   call print_date(time_to_skip, ' implied skip-to-date')
   call print_time(time_to_skip, ' implied skip-to-time')
   write(string1,*)'Namelist set to skip beyond the last observation time.'
   call error_handler(E_ERR,routine,string1,source)
endif

end subroutine DetermineNumEpochs


!======================================================================
!>  Sets the binwidth, halfbinwidth


subroutine DefineTemporalBins()

! These are variables that can be modified by this routine
!  integer,         intent(inout) GLOBAL :: Nepochs
!  type(time_type), intent(out)   GLOBAL :: binwidth     ! period of interest around center
!  type(time_type), intent(out)   GLOBAL :: halfbinwidth ! half that period
!  bin_width_days, bin_width_seconds  GLOBAL from namelist

integer :: nbins
type(time_type) :: test_time

if (Nepochs == 1) then

   ! If there is only 1 time in the file, then just do the right thing 
   ! and ignore the user input

   binwidth = set_time(60,0)  ! one minute, more than enough

elseif ( (bin_width_days <  0) .and. (bin_width_seconds >= 0) .or. &
         (bin_width_days >= 0) .and. (bin_width_seconds <  0) ) then

   write(string1,*)'bin_width_[days,seconds] must be non-negative, they are ', &
   bin_width_days, bin_width_seconds
   call error_handler(E_ERR,'DefineTemporalBins',string1,source, &
          text2='namelist parameter out-of-bounds. Fix and try again.')

elseif ( (bin_width_days <= 0) .and. (bin_width_seconds <= 0) ) then

   ! This is the 'default' case ... use all possible, up to "max_num_bins".
   ! 'space-filling' strategy: bin width and bin separation are same.
   ! Using Nepochs that comes from the number of distinct times in the files.

   binwidth  = (obsTN - obsT1) / (Nepochs - 1)
   if (Nepochs > max_num_bins) then
      write(string1,*)'default calculation results in ',Nepochs,' time bins.'
      write(string2,*)'namelist "max_num_bins" requests ',max_num_bins,'. Using this value.'
      call error_handler(E_MSG,'DefineTemporalBins',string1,text2=string2)
      Nepochs = max_num_bins
   endif
   obsTN = obsT1 + (Nepochs-1)*binwidth

else ! honor the user input

   binwidth  = set_time(bin_width_seconds, bin_width_days)
   test_time = obsT1
   nbins = 0
   COUNTBINS : do i = 1,max_num_bins
      if (test_time >  obsTN) exit COUNTBINS
      test_time = test_time + binwidth
      nbins     = nbins + 1
   enddo COUNTBINS

   ! Warn about falling off end ...
   if (nbins == max_num_bins) then
      write(string1,*)'namelist "max_num_bins" requests ',max_num_bins,'. Using this value.'
      call error_handler(E_MSG,'DefineTemporalBins',string1)
      Nepochs = nbins

   elseif (nbins == 0) then
      write(string1,*)'namelist settings for bin width results in no useful bins.'
      write(string2,*)'Stopping.'
      call error_handler(E_MSG,'DefineTemporalBins',string1,source,text2=string2)
      Nepochs = nbins
   endif
   obsTN = obsT1 + (Nepochs-1)*binwidth

endif

halfbinwidth = binwidth / 2

if ( verbose ) then
   call print_date(       obsT1,' DefineTemporalBins: start             date',logfileunit)
   call print_date(       obsTN,' DefineTemporalBins: end               date',logfileunit)
   call print_time(       obsT1,' DefineTemporalBins: start             time',logfileunit)
   call print_time(       obsTN,' DefineTemporalBins: end               time',logfileunit)
   call print_time(    binwidth,' DefineTemporalBins: requested     binwidth',logfileunit)
   call print_time(halfbinwidth,' DefineTemporalBins: implied   halfbinwidth',logfileunit)

   call print_date(       obsT1,' DefineTemporalBins: start             date')
   call print_date(       obsTN,' DefineTemporalBins: end               date')
   call print_time(       obsT1,' DefineTemporalBins: start             time')
   call print_time(       obsTN,' DefineTemporalBins: end               time')
   call print_time(    binwidth,' DefineTemporalBins: requested     binwidth')
   call print_time(halfbinwidth,' DefineTemporalBins: implied   halfbinwidth')
endif

end subroutine DefineTemporalBins


!======================================================================


subroutine SetSchedule(bin1time, num_epochs, fullwidth, halfwidth, &
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
      write(string1,'(''epoch '',i6,''  start'')')iepoch
      write(string2,'(''epoch '',i6,'' center'')')iepoch
      write(string3,'(''epoch '',i6,''    end'')')iepoch
      call print_time(bin_edges(1,iepoch),string1,logfileunit)
      call print_time(bin_center( iepoch),string2,logfileunit)
      call print_time(bin_edges(2,iepoch),string3,logfileunit)
      call print_time(bin_edges(1,iepoch),string1)
      call print_time(bin_center( iepoch),string2)
      call print_time(bin_edges(2,iepoch),string3)

      call print_date(bin_edges(1,iepoch),string1,logfileunit)
      call print_date(bin_center( iepoch),string2,logfileunit)
      call print_date(bin_edges(2,iepoch),string3,logfileunit)
      call print_date(bin_edges(1,iepoch),string1)
      call print_date(bin_center( iepoch),string2)
      call print_date(bin_edges(2,iepoch),string3)
   enddo
   write(logfileunit,*)
   write(     *     ,*)
endif

end subroutine SetSchedule


!======================================================================


function GetEnsSize()

!  Loop over all the metadata to count the number of ensemble members
!  available in the observation sequence file. We need this count to
!  allocate space for the rank histogram information. Since the rank
!  histogram will be created for the priors only ...

integer :: GetEnsSize

! Using 'seq' from global scope

integer :: i

GetEnsSize = 0

MetaDataLoop : do i=1, get_num_copies(seq)
   if(index(get_copy_meta_data(seq,i), 'prior ensemble member') > 0) &
                   GetEnsSize = GetEnsSize + 1
enddo MetaDataLoop

write(string1,'(''There are '',i4,'' ensemble members.'')') GetEnsSize
call error_handler(E_MSG,'GetEnsSize',string1)

end function GetEnsSize


!======================================================================


subroutine  SetIndices()

! integer, intent(out) :: obs_index, org_qc_index, dart_qc_index, &
!                         prior_mean_index,   posterior_mean_index,    &
!                         prior_spread_index, posterior_spread_index

! Using 'seq' and 'ens_size' from global scope

integer :: i, ens_count
character(len=metadatalength) :: metadata

obs_index              = -1
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
      if(index(metadata, 'truth'      ) > 0) obs_index = i
   else
      if(index(metadata, 'observation') > 0) obs_index = i
   endif

   if(index(metadata, 'prior ensemble mean'      ) > 0)       prior_mean_index = i
   if(index(metadata, 'posterior ensemble mean'  ) > 0)   posterior_mean_index = i
   if(index(metadata, 'prior ensemble spread'    ) > 0)     prior_spread_index = i
   if(index(metadata, 'posterior ensemble spread') > 0) posterior_spread_index = i

   if(index(metadata, 'prior ensemble member') > 0 .and. &
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

if ( prior_mean_index       < 0 ) then
   write(string1,*)'metadata:"prior ensemble mean" not found'
   call error_handler(E_MSG,'SetIndices',string1)
endif
if ( posterior_mean_index   < 0 ) then
   write(string1,*)'metadata:"posterior ensemble mean" not found'
   call error_handler(E_MSG,'SetIndices',string1)
endif
if ( prior_spread_index     < 0 ) then
   write(string1,*)'metadata:"prior ensemble spread" not found'
   call error_handler(E_MSG,'SetIndices',string1)
endif
if ( posterior_spread_index < 0 ) then
   write(string1,*)'metadata:"posterior ensemble spread" not found'
   call error_handler(E_MSG,'SetIndices',string1)
endif
if (           org_qc_index < 0 ) then
   write(string1,*)'metadata:"Quality Control" not found'
   call error_handler(E_MSG,'SetIndices',string1)
endif
if (          dart_qc_index < 0 ) then
   write(string1,*)'metadata:"DART quality control" not found'
   call error_handler(E_MSG,'SetIndices',string1)
endif

! Only require obs_index to be present; this allows the program
! to be run on obs_seq.[in,out] files which have no means or spread.
! Can still count number of obs, observation mean, ...

if ( obs_index < 0 ) then
   if ( use_zero_error_obs ) then
      write(string1,*)'metadata:"truth"       not found'
   else
      write(string1,*)'metadata:"observation" not found'
   endif
   call error_handler(E_ERR,'SetIndices',string1,source)
endif

!--------------------------------------------------------------------
! Echo what we found.
!--------------------------------------------------------------------

if ( use_zero_error_obs ) then
   write(string1,'(''"truth"                index '',i2,'' metadata '',a)') &
     obs_index, trim(get_copy_meta_data(seq,obs_index))
else
   write(string1,'(''"observation"          index '',i2,'' metadata '',a)') &
     obs_index, trim(get_copy_meta_data(seq,obs_index))
endif
call error_handler(E_MSG,'SetIndices',string1)

if (prior_mean_index > 0 ) then
   write(string1,'(''"prior mean"           index '',i2,'' metadata '',a)') &
        prior_mean_index, trim(get_copy_meta_data(seq,prior_mean_index))
   call error_handler(E_MSG,'SetIndices',string1)
endif

if (posterior_mean_index > 0 ) then
   write(string1,'(''"posterior mean"       index '',i2,'' metadata '',a)') &
        posterior_mean_index, trim(get_copy_meta_data(seq,posterior_mean_index))
   call error_handler(E_MSG,'SetIndices',string1)
endif

if (prior_spread_index > 0 ) then
   write(string1,'(''"prior spread"         index '',i2,'' metadata '',a)') &
        prior_spread_index, trim(get_copy_meta_data(seq,prior_spread_index))
   call error_handler(E_MSG,'SetIndices',string1)
endif

if (posterior_spread_index > 0 ) then
   write(string1,'(''"posterior spread"     index '',i2,'' metadata '',a)') &
        posterior_spread_index, trim(get_copy_meta_data(seq,posterior_spread_index))
   call error_handler(E_MSG,'SetIndices',string1)
endif

if (org_qc_index > 0 ) then
   write(string1,'(''"Quality Control"      index '',i2,'' metadata '',a)') &
         org_qc_index, trim(get_qc_meta_data(seq,org_qc_index))
   call error_handler(E_MSG,'SetIndices',string1)
endif

if (dart_qc_index > 0 ) then
   write(string1,'(''"DART quality control" index '',i2,'' metadata '',a)') &
        dart_qc_index, trim(get_qc_meta_data(seq,dart_qc_index))
   call error_handler(E_MSG,'SetIndices',string1)
endif

if ( any( (/ prior_mean_index, prior_spread_index/) < 0) ) then
   string1 = 'Observation sequence has no prior information.'
   string2 = 'You will still get a count, maybe observation value, incoming qc, ...'
   string3 = 'For simple information, you may want to use "obs_seq_to_netcdf" instead.'
   call error_handler(E_MSG, 'SetIndices', string1, source, text2=string2, text3=string3)
endif

has_posteriors = .true.
if ( any( (/ posterior_mean_index, posterior_spread_index /) < 0) ) then
   has_posteriors = .false.
   string1 = 'Observation sequence has no posterior information,'
   string2 = 'therefore - posterior diagnostics are not possible.'
   call error_handler(E_WARN, 'SetIndices', string1, source, text2=string2)
endif

end subroutine SetIndices


!======================================================================


function Rank_Histogram(copyvalues, observation_index, error_variance ) result(rank)

! Calculates the bin/rank
! We don't care about the QC value. If the ob wasn't assimilated
! the bin is meaningless.

real(r8),dimension(:), intent(in)  :: copyvalues
integer,               intent(in)  :: observation_index
real(r8),              intent(in)  :: error_variance
integer                            :: rank

! Local Variables

real(r8)                      :: obsvalue, mean, stddev
real(r8), dimension(ens_size) :: ensemble_values
real(r8), dimension(ens_size) :: sampling_noise

! Grab the observation value from the myriad copy values.

obsvalue = copyvalues(observation_index)
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


if ( .false. )  then ! DEBUG block
   write(*,*)'observation error variance is ',error_variance
   write(*,*)'observation          value is ',obsvalue
   write(*,*)'observation           rank is ',rank
   write(*,*)'noisy ensemble values are '
   write(*,*)ensemble_values
   write(*,*)
endif

end function Rank_Histogram


!======================================================================


   subroutine RPE(x,y)
      real(r8), intent(inout) :: x
      real(r8), intent(in) :: y
      x = x + y
   end subroutine RPE
   subroutine IPE(x,y)
      integer, intent(inout) :: x
      integer, intent(in) :: y
      x = x + y
   end subroutine IPE


!======================================================================
!> This function simply accumulates the appropriate sums.
!> The normalization occurs after all the data has been read, naturally.
!> ... spread is computed via sqrt(ensemble_spread**2 + observation_error**2).

  subroutine Bin3D(iqc, iepoch, iregion, flavor, trusted, &
                obsval,  obserrvar,  prmean,  prsprd,  pomean,  posprd, rank)

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

   call IPE(prior%Nposs(iepoch,iregion,flavor), 1)
   if (has_posteriors) &
      call IPE(poste%Nposs(iepoch,iregion,flavor), 1)

   !----------------------------------------------------------------------
   ! Enforce the use of trusted observations.
   !----------------------------------------------------------------------

   if ( trusted ) then

      call IPE(prior%Ntrusted(iepoch,iregion,flavor), 1)
      if (has_posteriors) &
         call IPE(poste%Ntrusted(iepoch,iregion,flavor), 1)

      ! Accrue the PRIOR quantities
      if ( any(iqc == trusted_prior_qcs) ) then
         call IPE(prior%Nused(      iepoch,iregion,flavor),      1    )
         call RPE(prior%observation(iepoch,iregion,flavor), obsmean   )
         call RPE(prior%ens_mean(   iepoch,iregion,flavor), priormean )
         call RPE(prior%bias(       iepoch,iregion,flavor), priorbias )
         call RPE(prior%rmse(       iepoch,iregion,flavor), priorsqerr)
         call RPE(prior%spread(     iepoch,iregion,flavor), priorspred)
         call RPE(prior%totspread(  iepoch,iregion,flavor), priorspredplus)
      else
         call IPE(prior%NbadDartQC(iepoch,iregion,flavor),      1     )
      endif

      ! Accrue the POSTERIOR quantities
      if (has_posteriors) then
         if ( any(iqc == trusted_poste_qcs) ) then
            call IPE(poste%Nused(      iepoch,iregion,flavor),     1    )
            call RPE(poste%observation(iepoch,iregion,flavor), obsmean  )
            call RPE(poste%ens_mean(   iepoch,iregion,flavor), postmean )
            call RPE(poste%bias(       iepoch,iregion,flavor), postbias )
            call RPE(poste%rmse(       iepoch,iregion,flavor), postsqerr)
            call RPE(poste%spread(     iepoch,iregion,flavor), postspred)
            call RPE(poste%totspread(  iepoch,iregion,flavor), postspredplus)
         else
            call IPE(poste%NbadDartQC(iepoch,iregion,flavor),      1    )
         endif
      endif

      return  ! EXIT THE BINNING ROUTINE
   endif

   !----------------------------------------------------------------------
   ! Proceed 'as normal'.
   !----------------------------------------------------------------------

   if ( iqc > QC_MAX_PRIOR ) then  ! prior and posterior failed

      call IPE(prior%NbadDartQC(iepoch,iregion,flavor),      1    )
      if (has_posteriors) &
         call IPE(poste%NbadDartQC(iepoch,iregion,flavor),      1    )

   else if ( iqc > QC_MAX_POSTERIOR ) then

      ! Then at least the prior (A.K.A. prior) is good
      call IPE(prior%Nused(      iepoch,iregion,flavor),      1    )
      call RPE(prior%observation(iepoch,iregion,flavor), obsmean   )
      call RPE(prior%ens_mean(   iepoch,iregion,flavor), priormean )
      call RPE(prior%bias(       iepoch,iregion,flavor), priorbias )
      call RPE(prior%rmse(       iepoch,iregion,flavor), priorsqerr)
      call RPE(prior%spread(     iepoch,iregion,flavor), priorspred)
      call RPE(prior%totspread(  iepoch,iregion,flavor), priorspredplus)

      ! However, the posterior is bad
      if (has_posteriors) &
         call IPE(poste%NbadDartQC(iepoch,iregion,flavor),      1    )

   else

      ! The prior is good
      call IPE(prior%Nused(      iepoch,iregion,flavor),      1    )
      call RPE(prior%observation(iepoch,iregion,flavor), obsmean   )
      call RPE(prior%ens_mean(   iepoch,iregion,flavor), priormean )
      call RPE(prior%bias(       iepoch,iregion,flavor), priorbias )
      call RPE(prior%rmse(       iepoch,iregion,flavor), priorsqerr)
      call RPE(prior%spread(     iepoch,iregion,flavor), priorspred)
      call RPE(prior%totspread(  iepoch,iregion,flavor), priorspredplus)

      ! The posterior is good
      if (has_posteriors) then
         call IPE(poste%Nused(      iepoch,iregion,flavor),      1   )
         call RPE(poste%observation(iepoch,iregion,flavor), obsmean  )
         call RPE(poste%ens_mean(   iepoch,iregion,flavor), postmean )
         call RPE(poste%bias(       iepoch,iregion,flavor), postbias )
         call RPE(poste%rmse(       iepoch,iregion,flavor), postsqerr)
         call RPE(poste%spread(     iepoch,iregion,flavor), postspred)
         call RPE(poste%totspread(  iepoch,iregion,flavor), postspredplus)
      endif

   endif

   ! The rank histogram binning is a bit of a peculiar situation.
   ! Only the prior is of interest ... so DART QCs of 0 1 2 3 are 'good'.
   ! There is some debate about whether we should be considering the
   ! 'outlier' observations (DART QC == 7), so that is namelist controlled.

   if (     (rank > 0) .and. create_rank_histogram ) then
      if ( any(iqc == hist_qcs(1:numqcvals) ) ) then 
         call IPE(prior%hist_bin(iepoch,iregion,flavor,rank), 1)
      endif
   endif

   end subroutine Bin3D


!======================================================================


   subroutine WriteNetCDF(fname)
   character(len=*), intent(in) :: fname

   integer :: ncid, i, nobs, typesdimlen
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
      'Compiler does not support required kinds of variables.',source)
   endif

   call nc_check(nf90_create(path = trim(fname), cmode = nf90_share, &
            ncid = ncid), 'obs_diag:WriteNetCDF', 'create '//trim(fname))

   !----------------------------------------------------------------------------
   ! Write Global Attributes
   !----------------------------------------------------------------------------

   call DATE_AND_TIME(crdate,crtime,crzone,values)
   write(string1,'(''YYYY MM DD HH MM SS = '',i4,5(1x,i2.2))') &
                  values(1), values(2), values(3), values(5), values(6), values(7)
   call nc_check(nf90_put_att(ncid, NF90_GLOBAL, 'creation_date', trim(string1) ), &
              'WriteNetCDF', 'put_att creation_date '//trim(fname))
   call nc_check(nf90_put_att(ncid, NF90_GLOBAL, 'obs_diag_source', source), &
              'WriteNetCDF', 'put_att obs_diag_source '//trim(fname))
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

   FILEloop : do i = 1,num_input_files

     write(string1,'(''obs_seq_file_'',i3.3)')i
     call nc_check(nf90_put_att(ncid, NF90_GLOBAL, &
            trim(string1), trim(obs_sequence_name(i)) ), &
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
   do ivar = 1,max_defined_types_of_obs

      nobs = sum(prior%Nposs(:,:,ivar))

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
              name='obstypes', len = max_defined_types_of_obs,    dimid = TypesDimID), &
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
             (/ epoch_center(1), epoch_center(Nepochs) /)), &
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
             (/ epoch_edges(1,1), epoch_edges(2,Nepochs) /)), &
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

   call nc_check(nf90_put_var(ncid, TypesVarId, (/ (i,i=1,max_defined_types_of_obs) /) ), &
              'WriteNetCDF', 'types:put_var')

   call nc_check(nf90_put_var(ncid, TypesMetaVarID, obs_type_strings(1:max_defined_types_of_obs)), &
              'WriteNetCDF', 'typesmeta:put_var')

   call nc_check(nf90_put_var(ncid, RegionVarID, (/ (i,i=1,Nregions) /) ), &
              'WriteNetCDF', 'region:put_var')

   call nc_check(nf90_put_var(ncid, BoundsVarID, (/ 1, 2 /)), &
              'WriteNetCDF', 'bounds:put_var')

   call nc_check(nf90_put_var(ncid, TimeVarID, epoch_center ), &
              'WriteNetCDF', 'time:put_var')

   call nc_check(nf90_put_var(ncid, TimeBoundsVarID, epoch_edges ), &
              'WriteNetCDF', 'time_bounds:put_var')

   call nc_check(nf90_put_var(ncid, RegionNamesVarID, reg_names(1:Nregions)), &
             'WriteNetCDF', 'region_names:put_var')

   call nc_check(nf90_sync( ncid), 'WriteNetCDF', 'sync '//trim(fname))

   !----------------------------------------------------------------------------
   ! write the data we took such pains to collate ...
   !----------------------------------------------------------------------------

   if ( verbose ) write(*,*)'summary for Priors of time-region vars'
   if ( create_rank_histogram ) then
      ierr = WriteTRV(ncid, prior, TimeDimID, CopyDimID, RegionDimID, RankDimID)
   else
      ierr = WriteTRV(ncid, prior, TimeDimID, CopyDimID, RegionDimID)
   endif
   if ( verbose ) write(*,*)
   if ( verbose ) write(*,*)'summary for Posteriors of time-region vars'
   ierr = WriteTRV(ncid, poste,    TimeDimID, CopyDimID, RegionDimID)
   if ( verbose ) write(*,*)

   !----------------------------------------------------------------------------
   ! finish ...
   !----------------------------------------------------------------------------

   call nc_check(nf90_sync( ncid), 'WriteNetCDF', 'sync '//trim(fname))
   call nc_check(nf90_close(ncid), 'init_diag_output', 'close '//trim(fname))

   end subroutine WriteNetCDF


!======================================================================


   function WriteTRV(ncid, vrbl, TimeDimID, CopyDimID, RegionDimID, RankDimID)
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

      if ( verbose ) then
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
         rchunk(iregion,10,itime) = vrbl%Ntrusted(   itime,iregion,ivar)
         rchunk(iregion,11,itime) = vrbl%NDartQC_0(  itime,iregion,ivar)
         rchunk(iregion,12,itime) = vrbl%NDartQC_1(  itime,iregion,ivar)
         rchunk(iregion,13,itime) = vrbl%NDartQC_2(  itime,iregion,ivar)
         rchunk(iregion,14,itime) = vrbl%NDartQC_3(  itime,iregion,ivar)
         rchunk(iregion,15,itime) = vrbl%NDartQC_4(  itime,iregion,ivar)
         rchunk(iregion,16,itime) = vrbl%NDartQC_5(  itime,iregion,ivar)
         rchunk(iregion,17,itime) = vrbl%NDartQC_6(  itime,iregion,ivar)
         rchunk(iregion,18,itime) = vrbl%NDartQC_7(  itime,iregion,ivar)
         rchunk(iregion,19,itime) = vrbl%NDartQC_8(  itime,iregion,ivar)

      enddo
      enddo

      call nc_check(nf90_redef(ncid), 'WriteTRV', 'redef')

      ! Create netCDF variable name

      obsname = obs_type_strings(ivar)
      string1 = trim(obsname)//'_'//adjustl(vrbl%string)

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
            write(logfileunit,*)trim(string1)//' has ',ndata,'"rank"able observations.'
            write(     *     ,*)trim(string1)//' has ',ndata,'"rank"able observations.'
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

   end function WriteTRV


!======================================================================
!> We need to know the time of the first and last observations in the sequence,
!> primarily just to see if they intersect the desired Epoch window.
!> We also record these times so we can report the first/last times of all
!> observations in all the obs_seq files.

subroutine GetFirstLastObs(my_fname, my_sequence, my_obs1, my_obsN, &
                   my_seqT1, my_seqTN, my_AllseqT1, my_AllseqTN)

character(len=*),          intent(in)    :: my_fname
type(obs_sequence_type),   intent(in)    :: my_sequence
type(obs_type),            intent(out)   :: my_obs1
type(obs_type),            intent(out)   :: my_obsN
type(time_type),           intent(out)   :: my_seqT1
type(time_type),           intent(out)   :: my_seqTN
type(time_type), optional, intent(inout) :: my_AllseqT1  ! ALL observation sequences
type(time_type), optional, intent(inout) :: my_AllseqTN  ! ALL observation sequences

type(obs_def_type) :: obs_def

logical,         SAVE :: first_time = .true.
type(time_type), SAVE :: absolute_first, absolute_last

if ( .not. get_first_obs(my_sequence, my_obs1) ) then
   call error_handler(E_ERR,'obs_diag','No first observation in '//trim(my_fname),source)
endif
call get_obs_def(my_obs1,   obs_def)
my_seqT1 = get_obs_def_time(obs_def)

if ( .not. get_last_obs(my_sequence, my_obsN) ) then
   call error_handler(E_ERR,'obs_diag','No last observation in '//trim(my_fname),source)
endif
call get_obs_def(my_obsN,   obs_def)
my_seqTN = get_obs_def_time(obs_def)

! Capture a little information to assist in an error message if the
! namelist input does not intersect the observation sequence file.

if ( first_time ) then
   absolute_first = my_seqT1
   absolute_last  = my_seqTN
   first_time     = .false.
else
   if (my_seqT1 < absolute_first) absolute_first = my_seqT1
   if (my_seqTN > absolute_last ) absolute_last  = my_seqTN
endif

! these are always logged
call print_time(my_seqT1,'First observation time',logfileunit)
call print_time(my_seqTN,'Last  observation time',logfileunit)
call print_date(my_seqT1,'First observation date',logfileunit)
call print_date(my_seqTN,'Last  observation date',logfileunit)

if ( verbose ) then
   call print_time(my_seqT1,'First observation time')
   call print_time(my_seqTN,'Last  observation time')
   call print_date(my_seqT1,'First observation date')
   call print_date(my_seqTN,'Last  observation date')
endif

write(logfileunit,*)
write(*,*)

if (present(my_allseqT1)) my_AllseqT1 = absolute_first
if (present(my_allseqTN)) my_AllseqTN = absolute_last

end subroutine GetFirstLastObs


!======================================================================
!> Function to determine if the observation sequence file has any
!> observations in the desired time window.

function No_Time_Intersection(filename,sequence_T1,sequence_TN,first_time,last_time)

character(len=*), intent(in) :: filename
type(time_type),  intent(in) :: sequence_T1 !< first ob time in sequence
type(time_type),  intent(in) :: sequence_TN !< last ob time in sequence
type(time_type),  intent(in) :: first_time  !< first time of interest
type(time_type),  intent(in) :: last_time   !< last time of interest
logical                      :: No_Time_Intersection

character(len=*), parameter :: routine = 'No_Time_Intersection'

if ( sequence_T1 >= first_time .and. sequence_TN <= last_time ) then
   No_Time_Intersection = .false.
else
   No_Time_Intersection = .true.
endif

if ( sequence_TN < first_time ) then
   if ( verbose ) then
      string1 = '"'//trim(filename)//'" last obs before first time ... trying next file.'
      call error_handler(E_MSG, routine, string1)
   endif
endif

if ( sequence_T1 > last_time ) then
   if ( verbose ) then
      string1 = '"'//trim(filename)//'" first obs after last_time ... trying next file.'
      call error_handler(E_MSG, routine, string1)
   endif
endif

end function No_Time_Intersection


!======================================================================


   function DefineTrustedObs()

   ! Count up the number of 'trusted' observations
   ! Check to make sure trusted observations desired are supported.
   ! The namelist specifies 'trusted_obs(:)' ... we create a local list that
   ! is the intersection of the trusted_obs desired and those on hand ...

   integer :: DefineTrustedObs

   ! character(len=*), intent(out) :: trusted_list    the list of trusted observation types
   ! integer,          intent(out) :: DefineTrustedObs   the number of trusted observations types

   integer :: i, num_trusted, ikind
   logical :: matched
   character(len=NF90_MAX_NAME) :: string1

   num_trusted = 0

   CountTrusted : do i = 1,MaxTrusted
      if (trim(trusted_obs(i)) == 'null') exit CountTrusted

      matched = .false.
      VerifyTrusted : do ikind = 1,max_defined_types_of_obs
         if (trim(trusted_obs(i)) == trim( get_name_for_type_of_obs(ikind) )) then
            matched = .true.
            exit VerifyTrusted
         endif
      enddo VerifyTrusted

      if (matched) then
          num_trusted = num_trusted + 1
          trusted_list(num_trusted) = trim(trusted_obs(i))
      else
         write(string1,*)'trusted_obs "',trim(trusted_obs(i)),'" is not a supported observation type.'
         call error_handler(E_WARN, 'DefineTrustedObs', trim(string1), source)
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
         call error_handler(E_MSG, 'DefineTrustedObs', trim(trusted_list(i)) )
      enddo
   else
      write(string1,*)'There are no "trusted" observation types.'
      call error_handler(E_MSG, 'DefineTrustedObs', string1)
   endif

   DefineTrustedObs = num_trusted

   end function DefineTrustedObs


!======================================================================


function is_observation_trusted(obsname)

! Is the observation one that we 'trust'.
! If so, disregard the DART QC ==7 (outlier rejection) and use it to
! in the statistics calculations.
! Since each obs_sequence file can have its own header/table, the safest
! way is to compare the string to a list of trusted (string) observation types.

character(len=*), intent(in) :: obsname
logical                      :: is_observation_trusted

is_observation_trusted = .false.
rUtrusted : do i = 1,num_trusted
   if ( trim(obsname) == trim(trusted_list(i)) ) then
      is_observation_trusted = .true.
      exit rUtrusted
   endif
enddo rUtrusted

end function is_observation_trusted


!======================================================================


   subroutine count_QC_values(iqc,itime,ireg,iflav)

   integer, intent(in) :: iqc, itime, ireg, iflav

   if (        iqc == 0 ) then
      call IPE(prior%NDartQC_0(itime,ireg,iflav), 1)
      if (has_posteriors) &
         call IPE(poste%NDartQC_0(itime,ireg,iflav), 1)

   elseif (    iqc == 1 ) then
      call IPE(prior%NDartQC_1(itime,ireg,iflav), 1)
      if (has_posteriors) &
         call IPE(poste%NDartQC_1(itime,ireg,iflav), 1)

   elseif (    iqc == 2 ) then
      call IPE(prior%NDartQC_2(itime,ireg,iflav), 1)
      if (has_posteriors) &
         call IPE(poste%NDartQC_2(itime,ireg,iflav), 1)

   elseif (    iqc == 3 ) then
      call IPE(prior%NDartQC_3(itime,ireg,iflav), 1)
      if (has_posteriors) &
         call IPE(poste%NDartQC_3(itime,ireg,iflav), 1)

   elseif (    iqc == 4 ) then
      call IPE(prior%NDartQC_4(itime,ireg,iflav), 1)
      if (has_posteriors) &
         call IPE(poste%NDartQC_4(itime,ireg,iflav), 1)

   elseif (    iqc == 5 ) then
      call IPE(prior%NDartQC_5(itime,ireg,iflav), 1)
      if (has_posteriors) &
         call IPE(poste%NDartQC_5(itime,ireg,iflav), 1)

   elseif (    iqc == 6 ) then
      call IPE(prior%NDartQC_6(itime,ireg,iflav), 1)
      if (has_posteriors) &
         call IPE(poste%NDartQC_6(itime,ireg,iflav), 1)

   elseif (    iqc == 7 ) then
      call IPE(prior%NDartQC_7(itime,ireg,iflav), 1)
      if (has_posteriors) &
         call IPE(poste%NDartQC_7(itime,ireg,iflav), 1)

   elseif (    iqc == 8 ) then
      call IPE(prior%NDartQC_8(itime,ireg,iflav), 1)
      if (has_posteriors) &
         call IPE(poste%NDartQC_8(itime,ireg,iflav), 1)

   endif

   end subroutine count_QC_values


!======================================================================


subroutine NormalizeTRV()

integer :: ivar, iregion, iepoch

if ( verbose ) then
   write(logfileunit,*)'Normalizing time-region-variable quantities.'
   write(     *     ,*)'Normalizing time-region-variable quantities.'
endif

do ivar   = 1,num_obs_types
do iregion= 1,Nregions
do iepoch = 1,Nepochs

   ! The Priors

   if (  prior%Nused(      iepoch, iregion, ivar) == 0) then
         prior%observation(iepoch, iregion, ivar) = MISSING_R4
         prior%ens_mean(   iepoch, iregion, ivar) = MISSING_R4
         prior%bias(       iepoch, iregion, ivar) = MISSING_R4
         prior%rmse(       iepoch, iregion, ivar) = MISSING_R4
         prior%spread(     iepoch, iregion, ivar) = MISSING_R4
         prior%totspread(  iepoch, iregion, ivar) = MISSING_R4
   else
         prior%observation(iepoch, iregion, ivar) = &
         prior%observation(iepoch, iregion, ivar) / &
         prior%Nused(      iepoch, iregion, ivar)

         prior%ens_mean(   iepoch, iregion, ivar) = &
         prior%ens_mean(   iepoch, iregion, ivar) / &
         prior%Nused(      iepoch, iregion, ivar)

         prior%bias(       iepoch, iregion, ivar) = &
         prior%bias(       iepoch, iregion, ivar) / &
         prior%Nused(      iepoch, iregion, ivar)

         prior%rmse(       iepoch, iregion, ivar) = &
    sqrt(prior%rmse(       iepoch, iregion, ivar) / &
         prior%Nused(      iepoch, iregion, ivar) )

         prior%spread(     iepoch, iregion, ivar) = &
    sqrt(prior%spread(     iepoch, iregion, ivar) / &
         prior%Nused(      iepoch, iregion, ivar) )

         prior%totspread(  iepoch, iregion, ivar) = &
    sqrt(prior%totspread(  iepoch, iregion, ivar) / &
         prior%Nused(      iepoch, iregion, ivar) )

   endif

   ! The Posteriors aka analy
   if (has_posteriors) then

      if (  poste%Nused(      iepoch, iregion, ivar) == 0) then
            poste%observation(iepoch, iregion, ivar) = MISSING_R4
            poste%ens_mean(   iepoch, iregion, ivar) = MISSING_R4
            poste%bias(       iepoch, iregion, ivar) = MISSING_R4
            poste%rmse(       iepoch, iregion, ivar) = MISSING_R4
            poste%spread(     iepoch, iregion, ivar) = MISSING_R4
            poste%totspread(  iepoch, iregion, ivar) = MISSING_R4
      else
            poste%observation(iepoch, iregion, ivar) = &
            poste%observation(iepoch, iregion, ivar) / &
            poste%Nused(      iepoch, iregion, ivar)

            poste%ens_mean(   iepoch, iregion, ivar) = &
            poste%ens_mean(   iepoch, iregion, ivar) / &
            poste%Nused(      iepoch, iregion, ivar)

            poste%bias(       iepoch, iregion, ivar) = &
            poste%bias(       iepoch, iregion, ivar) / &
            poste%Nused(      iepoch, iregion, ivar)

            poste%rmse(       iepoch, iregion, ivar) = &
       sqrt(poste%rmse(       iepoch, iregion, ivar) / &
            poste%Nused(      iepoch, iregion, ivar) )

            poste%spread(     iepoch, iregion, ivar) = &
       sqrt(poste%spread(     iepoch, iregion, ivar) / &
            poste%Nused(      iepoch, iregion, ivar) )

            poste%totspread(  iepoch, iregion, ivar) = &
       sqrt(poste%totspread(  iepoch, iregion, ivar) / &
            poste%Nused(      iepoch, iregion, ivar) )

      endif
   endif
enddo
enddo
enddo

end subroutine NormalizeTRV


!======================================================================


end program obs_diag

