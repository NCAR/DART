! DART software - Copyright UCAR. This open source software is provided
! by UCAR, "as is", without charge, subject to all terms of use at
! http://www.image.ucar.edu/DAReS/DART/DART_download

!> The programs defines a series of epochs (periods of time) and geographic
!> regions and accumulates statistics for these epochs and regions.
!> All 'possible' observation types are treated separately.
!> The results are written to a netCDF file.
!> This version is specifically for identity streamflow observations, i.e.
!> as defined by the wrf-Hydro support.

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

use        types_mod, only : r4, r8, digits12, i8, metadatalength, &
                             MISSING_R8, MISSING_R4, MISSING_I 

use obs_sequence_mod, only : read_obs_seq, obs_type, obs_sequence_type,      &
                             get_first_obs, get_obs_from_key, get_obs_def,   &
                             get_copy_meta_data, get_obs_time_range,         &
                             get_time_range_keys, get_obs_values, init_obs,  &
                             assignment(=), get_num_copies, get_qc,          &
                             static_init_obs_sequence, destroy_obs_sequence, &
                             read_obs_seq_header, get_last_obs, destroy_obs, &
                             get_num_qc, get_qc_meta_data

use      obs_def_mod, only : obs_def_type, get_obs_def_error_variance, &
                             get_obs_def_time, get_obs_def_location, &
                             get_obs_def_type_of_obs

use     obs_kind_mod, only : get_quantity_for_type_of_obs, &
                             get_name_for_type_of_obs

use     location_mod, only : location_type, get_location, set_location_missing,   &
                             write_location, operator(/=), is_location_in_region, &
                             set_location, query_location, LocationDims,          &
                             is_vertical, VERTISUNDEF, VERTISSURFACE,             &
                             VERTISLEVEL, VERTISPRESSURE, VERTISHEIGHT,           &
                             VERTISSCALEHEIGHT

use time_manager_mod, only : time_type, set_date, set_time, get_time, print_time, &
                             set_calendar_type, print_date, GREGORIAN,            &
                             operator(*),  operator(+),  operator(-),             &
                             operator(>),  operator(<),  operator(/),             &
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

character(len=*), parameter :: source = 'streamflow_obs_diag.f90'

!---------------------------------------------------------------------

integer, parameter :: MaxGages   = 200
integer, parameter :: MaxLevels  = 50
integer, parameter :: MaxRegions = 50
integer, parameter :: MaxTrusted = 50
integer, parameter :: stringlength = 32

integer :: max_defined_types_of_obs

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
integer(i8) :: flavor 
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

logical :: out_of_range, keeper

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

!>@todo remove after verifying NbiqQC, NbadIZ not used in plotting scripts
real(r8):: rat_cri               = 5000.0_r8 ! QC ratio
real(r8):: input_qc_threshold    = 3.0_r8    ! maximum NCEP QC factor

!-----------------------------------------------------------------------
! Namelist with default values

character(len=256) :: obs_sequence_name(max_num_input_files) = ''
character(len=256) :: obs_sequence_list = ''
integer, dimension(6) :: first_bin_center = (/ 2003, 1, 1, 0, 0, 0 /)
integer, dimension(6) :: last_bin_center  = (/ 2003, 1, 2, 0, 0, 0 /)
integer, dimension(6) :: bin_separation   = (/    0, 0, 0, 6, 0, 0 /)
integer, dimension(6) :: bin_width        = (/    0, 0, 0, 6, 0, 0 /)
integer, dimension(6) :: time_to_skip     = (/    0, 0, 1, 0, 0, 0 /)
integer :: max_num_bins = 1000       ! maximum number of temporal bins to consider

real(r8), dimension(MaxLevels)   :: plevel       = MISSING_R8 ! pressure level (hPa)
real(r8), dimension(MaxLevels)   :: hlevel       = MISSING_R8 ! height (meters)
real(r8), dimension(MaxLevels)   :: mlevel       = MISSING_R8 ! model level (integer index)
real(r8), dimension(MaxLevels+1) :: plevel_edges = MISSING_R8 ! pressure level (hPa)
real(r8), dimension(MaxLevels+1) :: hlevel_edges = MISSING_R8 ! height (meters)
real(r8), dimension(MaxLevels+1) :: mlevel_edges = MISSING_R8 ! model levels (nondimensional)

integer :: Nregions = 0
real(r8), dimension(MaxRegions) :: lonlim1= MISSING_R8, lonlim2= MISSING_R8
real(r8), dimension(MaxRegions) :: latlim1= MISSING_R8, latlim2= MISSING_R8
character(len=stringlength), dimension(MaxRegions) :: reg_names = 'null'
type(location_type), dimension(MaxRegions) :: min_loc, max_loc

character(len=stringlength), dimension(MaxTrusted) :: trusted_obs = 'null'

logical :: print_mismatched_locs = .false.
logical :: print_obs_locations   = .false.
logical :: verbose               = .false.
logical :: outliers_in_histogram = .false.
logical :: create_rank_histogram = .true.
logical :: use_zero_error_obs    = .false.

namelist /obs_diag_nml/ obs_sequence_name, obs_sequence_list,                 &
                       first_bin_center, last_bin_center,                     &
                       bin_separation, bin_width, time_to_skip, max_num_bins, &
                       plevel, hlevel, mlevel, &
                       Nregions, lonlim1, lonlim2, latlim1, latlim2,          &
                       reg_names, print_mismatched_locs, print_obs_locations, &
                       create_rank_histogram, outliers_in_histogram,          &
                       plevel_edges, hlevel_edges, mlevel_edges,              &
                       verbose, trusted_obs, use_zero_error_obs

!-----------------------------------------------------------------------
! Variables used to accumulate the statistics.
!-----------------------------------------------------------------------

!>@todo must be a more clever way to relate the copy_names to the components

integer, parameter :: Ncopies = 23
character(len=stringlength), dimension(Ncopies) :: copy_names =                  &
   (/ 'Nposs      ', 'Nused      ', 'NbigQC     ', 'NbadIZ     ', 'NbadUV     ', &
      'NbadLV     ', 'rmse       ', 'bias       ', 'spread     ', 'totalspread', &
      'NbadDARTQC ', 'observation', 'ens_mean   ', 'N_trusted  ',                &
      'N_DARTqc_0 ', 'N_DARTqc_1 ', 'N_DARTqc_2 ', 'N_DARTqc_3 ', 'N_DARTqc_4 ', &
      'N_DARTqc_5 ', 'N_DARTqc_6 ', 'N_DARTqc_7 ', 'N_DARTqc_8 '                /)

type TLRV_type
   ! statistics by time-level-region-variable
   integer ::     time_dim = 1
   integer ::    level_dim = 2
   integer ::   region_dim = 3
   integer :: variable_dim = 4
   character(len=8) :: string
   integer :: num_times = 0, num_levels = 0, num_regions = 0, num_variables = 0
   integer,  dimension(:,:,:,:), pointer :: Nposs, Nused, Ntrusted
   integer,  dimension(:,:,:,:), pointer :: NbigQC ! # original QC values >= input_qc_threshold
   integer,  dimension(:,:,:,:), pointer :: NbadIZ ! # bad (ie huge) Innovation Zscore
   integer,  dimension(:,:,:,:), pointer :: NbadUV ! # unmatched U/V wind pairs
   integer,  dimension(:,:,:,:), pointer :: NbadLV ! # obs above/below top/bottom
   real(r8), dimension(:,:,:,:), pointer :: rmse, bias, spread, totspread
   integer,  dimension(:,:,:,:), pointer :: NbadDartQC ! # bad DART QC values
   real(r8), dimension(:,:,:,:), pointer :: observation, ens_mean
   integer,  dimension(:,:,:,:), pointer :: NDartQC_0, NDartQC_1, NDartQC_2, NDartQC_3
   integer,  dimension(:,:,:,:), pointer :: NDartQC_4, NDartQC_5, NDartQC_6, NDartQC_7
   integer,  dimension(:,:,:,:), pointer :: NDartQC_8
   integer,  dimension(:,:,:,:,:), pointer :: hist_bin => NULL()
end type TLRV_type

type(TLRV_type) :: poste,    prior

type(time_type), allocatable, dimension(:)   :: bin_center
type(time_type), allocatable, dimension(:,:) :: bin_edges
real(digits12),  allocatable, dimension(:)   :: epoch_center
real(digits12),  allocatable, dimension(:,:) :: epoch_edges
integer,         allocatable, dimension(:)   :: obs_used_in_epoch

!-----------------------------------------------------------------------
! General purpose variables
!-----------------------------------------------------------------------

integer  :: Nplevels=0, Nhlevels=0, Nmlevels=0
integer  :: iregion, iepoch, ivar, ifile, num_obs_in_epoch
real(r8) :: obslon, obslat, obslevel, obsloc3(3)

integer  :: obsindex, i, iunit, io, ikind
integer  :: seconds, days, Nepochs, num_input_files

integer  :: num_trusted
logical  :: trusted

integer  :: ivert
integer  :: level_index
integer  :: Nlevels, ilev

integer,  allocatable :: which_vert(:) ! relates kind of level for each obs kind
real(r8), allocatable :: scale_factor(:) ! to convert to plotting units
integer,  allocatable :: ob_defining_vert(:) ! obs index defining vert coord type
character(len=stringlength), allocatable :: obs_type_strings(:)

! These pairs of variables are used when we diagnose which observations
! are far from the background.
integer, parameter :: MaxSigmaBins = 100
integer  :: nsigma(0:MaxSigmaBins) = 0
integer  :: nsigmaUnit, indx

real(r8) :: pr_zscore, po_zscore

type(time_type) :: TimeMin, TimeMax    ! of the entire period of interest
type(time_type) :: beg_time, end_time  ! of the particular bin
type(time_type) :: binsep, binwidth, halfbinwidth
type(time_type) :: seqT1, seqTN        ! first,last time in entire observation sequence
type(time_type) :: AllseqT1, AllseqTN  ! first,last time in ALL observation sequences
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

! This is the part that makes identity observations hard.
! Rather than knowing exactly how many observation types we are supporting,
! we have to have one for each state vector element - each identity obs.

num_obs_types = grok_observation_names()

allocate( which_vert(num_obs_types), &
        scale_factor(num_obs_types), &
    obs_type_strings(num_obs_types), &
    ob_defining_vert(num_obs_types))

which_vert       = VERTISUNDEF
scale_factor     = 1.0_r8
ob_defining_vert = -1

! Read the namelist

call find_namelist_in_file('input.nml', 'obs_diag_nml', iunit)
read(iunit, nml = obs_diag_nml, iostat = io)
call check_namelist_read(iunit, io, 'obs_diag_nml')

! Record the namelist values used for the run ...
if (do_nml_file()) write(nmlfileunit, nml=obs_diag_nml)
if (do_nml_term()) write(    *      , nml=obs_diag_nml)

num_input_files = set_filename_list(obs_sequence_name, obs_sequence_list, 'obs_diag')

! Check to see if we are including the outlier observations in the
! rank histogram calculation.

if ( outliers_in_histogram ) then
   numqcvals = size(hist_qcs)
else
   numqcvals = size(hist_qcs) - 1
endif

! Now that we have input, do some checking and setup

call set_calendar_type(GREGORIAN)
call Namelist2Times()    ! convert namelist times to DART times
call DefineTimeBins()    ! Set the temporal binning variables
call SetPressureLevels(plevel, plevel_edges, Nplevels)
call SetHeightLevels(  hlevel, hlevel_edges, Nhlevels)
call SetModelLevels(   mlevel, mlevel_edges, Nmlevels)
call DefineRegions()
call CountTrustedObsTypes()
call SetScaleFactors() ! for plotting purposes

Nlevels = maxval((/ Nplevels, Nhlevels, Nmlevels /))
call InitializeVariables( Nepochs, Nlevels, Nregions, num_obs_types)

! Open file for histogram of innovations, as a function of standard deviation.

nsigmaUnit = open_file('LargeInnov.txt',form='formatted',action='write')
write(nsigmaUnit,'(a)')'Any observations flagged as bad are dumped into the last bin.'
write(nsigmaUnit,'(a)')'   day   secs    lon      lat    level         obs    prior&
                       &   zscore   key   kind'

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

   call GetFirstLastObs(obs_sequence_name(ifile), seq, obs1, obsN, seqT1, seqTN, &
                        AllseqT1, AllseqTN)

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
         allocate(prior%hist_bin( Nepochs, Nlevels, Nregions, num_obs_types, ens_size+1))
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

   call SetIndices( obs_index, org_qc_index, dart_qc_index, &
            prior_mean_index,   posterior_mean_index,   &
            prior_spread_index, posterior_spread_index)

   ! Loop over all potential time periods ... the observation sequence
   ! files are not required to be in any particular order.

   EpochLoop : do iepoch = 1, Nepochs

      beg_time = bin_edges(1,iepoch)
      end_time = bin_edges(2,iepoch)

      ! Using linked list information in the observation sequence,
      ! find the number of observations that are within this epoch.

      call get_obs_time_range(seq, beg_time, end_time, key_bounds, &
                  num_obs_in_epoch, out_of_range)

      if( num_obs_in_epoch == 0 ) then
         if ( verbose ) then
            write(logfileunit,*)' No observations in epoch ',iepoch,' cycling ...'
            write(     *     ,*)' No observations in epoch ',iepoch,' cycling ...'
         endif
         cycle EpochLoop
      endif

      write(logfileunit, *) 'num_obs_in_epoch (', iepoch, ') = ', num_obs_in_epoch
      write(     *     , *) 'num_obs_in_epoch (', iepoch, ') = ', num_obs_in_epoch

      ! actually get the indices (keys) to the observations of interest

      allocate(keys(num_obs_in_epoch))

      call get_time_range_keys(seq, key_bounds, num_obs_in_epoch, keys)

      ObservationLoop : do obsindex = 1, num_obs_in_epoch

         ! 'flavor' is from the 'master list' in the obs_kind_mod.f90
         ! each obs_seq.final file has their own private kind - which
         ! gets mapped to the 'master list', if you will.

         call get_obs_from_key(seq, keys(obsindex), observation)
         call get_obs_def(observation, obs_def)

         flavor    = get_obs_def_type_of_obs(obs_def)
         obs_time  = get_obs_def_time(obs_def)
         obs_loc   = get_obs_def_location(obs_def)
         obsloc3   = get_location(obs_loc)
         ivert     = nint(query_location(obs_loc))

         obslon    = obsloc3(1) ! [  0, 360]
         obslat    = obsloc3(2) ! [-90,  90]
         obslevel  = obsloc3(3) ! variable-dependent

         trusted = .false.

         ! Check to see if it is an identity observation.
         ! All streamflow observations are identity observations at this point
         if (flavor < 0) then

            ! Not every model location has a streamflow observation, but the
            ! location can be anywhere in the model state, so the 'flavor' can
            ! be quite large. We cannot allocate that much space, so we create
            ! a lookup array  ....

            max_defined_types_of_obs = lookup_flavor(flavor)

            if (flavor > num_obs_types) then
               write(string1,*)'unsupported identity observation ',-1*flavor
               write(string2,*)'Increase "num_obs_types" and recompile.'
               call error_handler(E_ERR, 'obs_diag', string1, source, text2=string2)
            endif

            obsname = 'STREAM_FLOW'
            Nidentity = Nidentity + 1
         endif

         ! same sort of thing for the scale factors
         if ( use_zero_error_obs ) then
            obs_error_variance = 0.0_r8
         else
            obs_error_variance = get_obs_def_error_variance(obs_def)
         endif
         obs_error_variance = obs_error_variance * &
                       scale_factor(flavor) * scale_factor(flavor)

         ! Check consistency of the vertical coordinate system
         ! Sometimes observations of the same flavor from different sources
         ! are defined on different vertical coordinates.

         call CheckVertical(obs_loc, flavor)

         ! Figure out which level the observation relates to ...

         level_index = ParseLevel(obs_loc, obslevel, flavor)

         if ( 1 == 0 ) then
            write(*,*)'obsindx ',obsindex, keys(obsindex), obsloc3(3), level_index
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

         call get_qc(observation, qc)

         if ( dart_qc_index > 0 ) then
            qc_value = qc(dart_qc_index)
         else
            ! If there is no dart_qc, this must be a case where we 
            ! are interested only in getting the location information.
            qc_value = 0
         endif

         ! Check to see if there are any observations with wild values
         ! and a DART QC flag that is inconsistent. I checked it once
         ! with qc_value < 4 and found that ONLY the posterior values
         ! were bad. While the underlying bug is being fixed, the workaround
         ! is to simply manually set the DART QC value such that the
         ! posterior is flagged as bad.  I cannot acutally tell if the
         ! observation is supposed to have been assimilated or just evaluated,
         ! so I erred on the conservative side. TJH 24 Aug 2007

         ! DEBUG block for strange looking observations:
         if ( .false. ) then
            call get_obs_values(observation, copyvals)
            ! add your own if test here and enable logic
        !   if ( obsindex == 311 ) then
            if (any(copyvals(2:size(copyvals)) /= -888888.0) .and. &
                  (qc_value == 4)) then
        !   if (abs(prior_mean(1)) > 1000 .and. qc_value < 4) then
        !   if (.true.) then
              write(*,*)
              write(*,*)'Observation index is ',keys(obsindex),' and has:'
              write(*,*)'observation type       is',' '//trim(obsname)
              write(*,*)'flavor                 is',flavor
              write(*,*)'obs              value is',obs(1)
              write(*,*)'prior_mean       value is',prior_mean(1)
              write(*,*)'posterior_mean   value is',posterior_mean(1)
              write(*,*)'prior_spread     value is',prior_spread(1)
              write(*,*)'posterior_spread value is',posterior_spread(1)
              write(*,*)'DART QC          value is',qc_value
              write(*,*)'observation  is   trusted',trusted
              do i= 1,num_copies
                 write(*,*)copyvals(i),trim(get_copy_meta_data(seq,i))
              enddo
              write(*,*)
         endif
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

         ! Scale the quantities so they plot sensibly.

         obs(1)  = obs(1)             *scale_factor(flavor)
         pr_mean = prior_mean(1)      *scale_factor(flavor)
         po_mean = posterior_mean(1)  *scale_factor(flavor)
         pr_sprd = prior_spread(1)    *scale_factor(flavor)  ! standard deviations
         po_sprd = posterior_spread(1)*scale_factor(flavor)  ! standard deviations

         ! DEBUG block Summary of observation knowledge at this point

         if ( .false. ) then
            write(*,*)
            write(*,*)'observation #,flavor ', obsindex, flavor
            write(*,*)'obs(1), qc ', obs(1), qc
            write(*,*)'obs_error_variance ', obs_error_variance
            call print_time(obs_time,'time is')
            write(*,*)'pr_mean, po_mean ', pr_mean, po_mean
            write(*,*)'pr_sprd, po_sprd ', pr_sprd, po_sprd
            write(*,*)'obslon/obslat ', obslon, obslat
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
            ivert = nint(query_location(obs_loc))

            write(nsigmaUnit,'(i7,1x,i5,1x,2f8.2,i7,1x,2f13.2,f8.1,2i7)') &
                 days, seconds, obslon, obslat, ivert, &
                 obs(1), pr_mean, pr_zscore, keys(obsindex), flavor
         endif

         ! At this point, the observation has passed all checks.

         obs_used_in_epoch(iepoch) = obs_used_in_epoch(iepoch) + 1

         ! If needed, calculate the rank histogram bin (once!) for
         ! this observation - even if the QC value is bad.

         if ( create_rank_histogram ) then
            call get_obs_values(observation, copyvals)
            rank_histogram_bin = Rank_Histogram(copyvals, obs_index, obs_error_variance )
         endif

         ! We have Nregions of interest.

         Areas : do iregion =1, Nregions

            keeper = is_location_in_region( obs_loc, min_loc(iregion), max_loc(iregion) )
            if ( .not. keeper ) cycle Areas

            !===========================================================
            ! temporal statistics part ... all the 'evolution' variables
            !===========================================================

            ! Reject observations too high or too low without counting it
            ! as a possible observation for this region.

            if ((level_index < 1) .or. (level_index > Nlevels)) then
               prior%NbadLV(iepoch,:,iregion,flavor) = &
               prior%NbadLV(iepoch,:,iregion,flavor) + 1
               if (has_posteriors) then
                  poste%NbadLV(iepoch,:,iregion,flavor) = &
                  poste%NbadLV(iepoch,:,iregion,flavor) + 1
               endif
               cycle Areas
            endif

            ! Count original QC values 'of interest' ...
            !>@todo remove after verifying NbiqQC not used in plotting scripts

            if (      org_qc_index  > 0 ) then
               if (qc(org_qc_index) > input_qc_threshold ) then
               call IPE(prior%NbigQC(iepoch,level_index,iregion,flavor), 1)
                  if (has_posteriors) &
                     call IPE(poste%NbigQC(iepoch,level_index,iregion,flavor), 1)
               endif
            endif

            ! Count DART QC values

            call CountDartQC_4D(qc_value, iepoch, level_index, iregion, &
                    flavor, prior, poste, posterior_mean=po_mean)

            ! Count 'large' innovations
            !>@todo remove after verifying NbiqIZ not used in plotting scripts

            if( pr_zscore > rat_cri ) then
               call IPE(prior%NbadIZ(iepoch,level_index,iregion,flavor), 1)
            endif

            if(po_zscore > rat_cri .and. has_posteriors) then
               call IPE(poste%NbadIZ(iepoch,level_index,iregion,flavor), 1)
            endif

            ! Do all the heavy lifting

            call Bin4D(qc_value, iepoch, level_index, iregion, flavor, trusted, &
                  obs(1), obs_error_variance, pr_mean, pr_sprd, po_mean, po_sprd, &
                  rank_histogram_bin)

            !===========================================================
            ! end of time series statistics
            !===========================================================

            if ( obs_time < skip_time ) cycle Areas

         enddo Areas

      enddo ObservationLoop

      deallocate(keys)

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

call Normalize4Dvars()

! Print final summary.

write(*,*)
write(*,*) '# observations used  : ',sum(obs_used_in_epoch)
write(*,*) 'Count summary over all regions - obs may count for multiple regions:'
write(*,*) '# identity           : ',Nidentity
write(*,*) '# bad Innov  (ratio) : ',sum(prior%NbadIZ)
write(*,*) '# bad UV (wind pairs): ',sum(prior%NbadUV)
write(*,*) '# bad Level          : ',sum(prior%NbadLV(:,1,:,:))
write(*,*) '# big (original) QC  : ',sum(prior%NbigQC)
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
write(logfileunit,*) '# bad Innov  (ratio) : ',sum(prior%NbadIZ)
write(logfileunit,*) '# bad UV (wind pairs): ',sum(prior%NbadUV)
write(logfileunit,*) '# bad Level          : ',sum(prior%NbadLV(:,1,:,:))
write(logfileunit,*) '# big (original) QC  : ',sum(prior%NbigQC)
write(logfileunit,*) '# bad DART QC prior  : ',sum(prior%NbadDartQC)
if (has_posteriors) write(logfileunit,*) '# bad DART QC post   : ',sum(poste%NbadDartQC)
write(logfileunit,*) '# priorQC 7 postQC 2 : ',num_ambiguous
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

! If the namelist input does not result in some observations, print
! a little summary that may result in better user input.

if ( sum(obs_used_in_epoch) == 0 ) then
   write(    *      ,*)
   write(logfileunit,*)
   call print_date(AllseqT1,' First observation date',logfileunit)
   call print_date(AllseqTN,' Last  observation date',logfileunit)
   call print_date( TimeMin,' First REQUESTED   date',logfileunit)
   call print_date( TimeMax,' Last  REQUESTED   date',logfileunit)
   call print_date(AllseqT1,' First observation date')
   call print_date(AllseqTN,' Last  observation date')
   call print_date( TimeMin,' First REQUESTED   date')
   call print_date( TimeMax,' Last  REQUESTED   date')
   write(string1,*)'NO OBSERVATIONS in requested time bins.'
   call error_handler(E_ERR,'obs_diag',string1,source)
endif

call WriteNetCDF('obs_diag_output.nc')

call DestroyVariables()
call error_handler(E_MSG,'obs_diag','Finished successfully.')
call finalize_utilities()


CONTAINS

!======================================================================
! These routines use common variables from the scope of this file.
! If it's not in the argument list ... it's scoped within this file.
!======================================================================


subroutine InitializeVariables( ntimes, nlevs, nareas, ntypes )

! Global variables set in this routine:
! type(TLRV_type), intent(out) :: poste,    prior

integer, intent(in)  :: ntimes
integer, intent(in)  :: nlevs
integer, intent(in)  :: nareas
integer, intent(in)  :: ntypes

allocate(prior%rmse(       ntimes, nlevs, nareas, ntypes), &
         prior%bias(       ntimes, nlevs, nareas, ntypes), &
         prior%spread(     ntimes, nlevs, nareas, ntypes), &
         prior%totspread(  ntimes, nlevs, nareas, ntypes), &
         prior%observation(ntimes, nlevs, nareas, ntypes), &
         prior%ens_mean(   ntimes, nlevs, nareas, ntypes), &
         prior%Nposs(      ntimes, nlevs, nareas, ntypes), &
         prior%Nused(      ntimes, nlevs, nareas, ntypes), &
         prior%NbigQC(     ntimes, nlevs, nareas, ntypes), &
         prior%NbadIZ(     ntimes, nlevs, nareas, ntypes), &
         prior%NbadUV(     ntimes, nlevs, nareas, ntypes), &
         prior%NbadLV(     ntimes, nlevs, nareas, ntypes), &
         prior%NbadDartQC( ntimes, nlevs, nareas, ntypes), &
         prior%Ntrusted(   ntimes, nlevs, nareas, ntypes), &
         prior%NDartQC_0(  ntimes, nlevs, nareas, ntypes), &
         prior%NDartQC_1(  ntimes, nlevs, nareas, ntypes), &
         prior%NDartQC_2(  ntimes, nlevs, nareas, ntypes), &
         prior%NDartQC_3(  ntimes, nlevs, nareas, ntypes), &
         prior%NDartQC_4(  ntimes, nlevs, nareas, ntypes), &
         prior%NDartQC_5(  ntimes, nlevs, nareas, ntypes), &
         prior%NDartQC_6(  ntimes, nlevs, nareas, ntypes), &
         prior%NDartQC_7(  ntimes, nlevs, nareas, ntypes), &
         prior%NDartQC_8(  ntimes, nlevs, nareas, ntypes))

prior%rmse        = 0.0_r8
prior%bias        = 0.0_r8
prior%spread      = 0.0_r8
prior%totspread   = 0.0_r8
prior%observation = 0.0_r8
prior%ens_mean    = 0.0_r8
prior%Nposs       = 0
prior%Nused       = 0
prior%NbigQC      = 0
prior%NbadIZ      = 0
prior%NbadUV      = 0
prior%NbadLV      = 0
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
prior%num_times     = ntimes
prior%num_levels    = nlevs
prior%num_regions   = nareas
prior%num_variables = ntypes

allocate(poste%rmse(       ntimes, nlevs, nareas, ntypes), &
         poste%bias(       ntimes, nlevs, nareas, ntypes), &
         poste%spread(     ntimes, nlevs, nareas, ntypes), &
         poste%totspread(  ntimes, nlevs, nareas, ntypes), &
         poste%observation(ntimes, nlevs, nareas, ntypes), &
         poste%ens_mean(   ntimes, nlevs, nareas, ntypes), &
         poste%Nposs(      ntimes, nlevs, nareas, ntypes), &
         poste%Nused(      ntimes, nlevs, nareas, ntypes), &
         poste%NbigQC(     ntimes, nlevs, nareas, ntypes), &
         poste%NbadIZ(     ntimes, nlevs, nareas, ntypes), &
         poste%NbadUV(     ntimes, nlevs, nareas, ntypes), &
         poste%NbadLV(     ntimes, nlevs, nareas, ntypes), &
         poste%NbadDartQC( ntimes, nlevs, nareas, ntypes), &
         poste%Ntrusted(   ntimes, nlevs, nareas, ntypes), &
         poste%NDartQC_0(  ntimes, nlevs, nareas, ntypes), &
         poste%NDartQC_1(  ntimes, nlevs, nareas, ntypes), &
         poste%NDartQC_2(  ntimes, nlevs, nareas, ntypes), &
         poste%NDartQC_3(  ntimes, nlevs, nareas, ntypes), &
         poste%NDartQC_4(  ntimes, nlevs, nareas, ntypes), &
         poste%NDartQC_5(  ntimes, nlevs, nareas, ntypes), &
         poste%NDartQC_6(  ntimes, nlevs, nareas, ntypes), &
         poste%NDartQC_7(  ntimes, nlevs, nareas, ntypes), &
         poste%NDartQC_8(  ntimes, nlevs, nareas, ntypes))

poste%rmse        = 0.0_r8
poste%bias        = 0.0_r8
poste%spread      = 0.0_r8
poste%totspread   = 0.0_r8
poste%observation = 0.0_r8
poste%ens_mean    = 0.0_r8
poste%Nposs       = 0
poste%Nused       = 0
poste%NbigQC      = 0
poste%NbadIZ      = 0
poste%NbadUV      = 0
poste%NbadLV      = 0
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
poste%num_times     = ntimes
poste%num_levels    = nlevs
poste%num_regions   = nareas
poste%num_variables = ntypes

end subroutine InitializeVariables


!======================================================================


subroutine DestroyVariables()

if (associated(prior%hist_bin)) deallocate(prior%hist_bin)
if (allocated(ens_copy_index))  deallocate(ens_copy_index)

deallocate(prior%rmse,        prior%bias,      prior%spread,    prior%totspread, &
           prior%observation, prior%ens_mean,  prior%Nposs,     prior%Nused,     &
           prior%NbigQC,      prior%NbadIZ,    prior%NbadUV,    prior%NbadLV,    &
           prior%NbadDartQC,  prior%Ntrusted)

deallocate(prior%NDartQC_0,   prior%NDartQC_1, prior%NDartQC_2, prior%NDartQC_3, &
           prior%NDartQC_4,   prior%NDartQC_5, prior%NDartQC_6, prior%NDartQC_7, &
           prior%NDartQC_8)

deallocate(poste%rmse,        poste%bias,      poste%spread,    poste%totspread, &
           poste%observation, poste%ens_mean,  poste%Nposs,     poste%Nused,     &
           poste%NbigQC,      poste%NbadIZ,    poste%NbadUV,    poste%NbadLV,    &
           poste%NbadDartQC,  poste%Ntrusted)

deallocate(poste%NDartQC_0,   poste%NDartQC_1, poste%NDartQC_2, poste%NDartQC_3, &
           poste%NDartQC_4,   poste%NDartQC_5, poste%NDartQC_6, poste%NDartQC_7, &
           poste%NDartQC_8)

deallocate(epoch_center, epoch_edges, bin_center, obs_used_in_epoch)

deallocate(obs_type_strings, which_vert, scale_factor)

end subroutine DestroyVariables


!======================================================================


function grok_observation_names()
!----------------------------------------------------------------------
! All identity observations are STREAM_FLOW
!----------------------------------------------------------------------

integer :: grok_observation_names

grok_observation_names = MaxGages

end function grok_observation_names


!======================================================================


subroutine Namelist2Times( )

! Global-scope variables read in this routine:
! first_bin_center   comes from namelist
!  last_bin_center   comes from namelist
!   bin_separation   comes from namelist
!        bin_width   comes from namelist
!     time_to_skip   comes from namelist
!
! Global-scope variables set in this routine:
! type(time_type), intent(out) :: beg_time     ! first_bin_center
! type(time_type), intent(out) :: end_time     ! last_bin_center
! type(time_type), intent(out) :: skip_time    ! time AFTER first bin leading edge
! type(time_type), intent(out) :: binsep       ! time between bin centers
! type(time_type), intent(out) :: binwidth     ! period of interest around center
! type(time_type), intent(out) :: halfbinwidth ! half that period

! We are using bin_separation and bin_width as offsets relative to the
! first time. to ensure this, the year and month must be zero.

logical :: error_out = .false.
integer :: seconds

! do some error-checking first

if ( (bin_separation(1) /= 0) .or. (bin_separation(2) /= 0) ) then
   write(string1,*)'bin_separation:year,month must both be zero, they are ', &
   bin_separation(1),bin_separation(2)
   call error_handler(E_WARN,'Namelist2Times',string1,source)
   error_out = .true.
endif

if ( (bin_width(1) /= 0) .or. (bin_width(2) /= 0) ) then
   write(string1,*)'bin_width:year,month must both be zero, they are ', &
   bin_width(1),bin_width(2)
   call error_handler(E_WARN,'Namelist2Times',string1,source)
   error_out = .true.
endif

if ( (time_to_skip(1) /= 0) .or. (time_to_skip(2) /= 0) ) then
   write(string1,*)'time_to_skip:year,month must both be zero, they are ', &
   time_to_skip(1),time_to_skip(2)
   call error_handler(E_WARN,'Namelist2Times',string1,source)
   error_out = .true.
endif

if ( error_out ) call error_handler(E_ERR,'Namelist2Times', &
    'namelist parameter out-of-bounds. Fix and try again.',source)

! Set time of first bin center
beg_time   = set_date(first_bin_center(1), first_bin_center(2), &
                      first_bin_center(3), first_bin_center(4), &
                      first_bin_center(5), first_bin_center(6) )

! Set time of _intended_ last bin center
end_time   = set_date( last_bin_center(1),  last_bin_center(2), &
                       last_bin_center(3),  last_bin_center(4), &
                       last_bin_center(5),  last_bin_center(6) )

seconds  = bin_separation(4)*3600 + bin_separation(5)*60 + bin_separation(6)
binsep   = set_time(seconds, bin_separation(3))

seconds  = bin_width(     4)*3600 + bin_width(     5)*60 + bin_width(     6)
binwidth = set_time(seconds, bin_width(     3))

halfbinwidth = binwidth / 2

! Set time that defines the end of the burn-in period.
! Anything at or after this time will be used to accumulate time-averaged quantities.
seconds   = time_to_skip(4)*3600 + time_to_skip(5)*60 + time_to_skip(6)

if ( (beg_time + set_time(seconds,time_to_skip(3))) <= halfbinwidth) then
   skip_time = set_time(0,0)
else
   skip_time = beg_time - halfbinwidth + set_time(seconds, time_to_skip(3))
endif

if ( verbose ) then
   call print_date(    beg_time,'requested first bincenter date',logfileunit)
   call print_date(    beg_time,'requested first bincenter date')
   call print_date(    end_time,'requested last  bincenter date',logfileunit)
   call print_date(    end_time,'requested last  bincenter date')
   call print_date(   skip_time,'implied   skip-to         date',logfileunit)
   call print_date(   skip_time,'implied   skip-to         date')
   write(logfileunit,*)
   write(*,*)
   call print_time(    beg_time,'requested first bincenter time',logfileunit)
   call print_time(    beg_time,'requested first bincenter time')
   call print_time(    end_time,'requested last  bincenter time',logfileunit)
   call print_time(    end_time,'requested last  bincenter time')
   call print_time(   skip_time,'implied   skip-to         time',logfileunit)
   call print_time(   skip_time,'implied   skip-to         time')
   write(logfileunit,*)
   write(*,*)
   call print_time(      binsep,'requested bin separation',logfileunit)
   call print_time(      binsep,'requested bin separation')
   call print_time(    binwidth,'requested bin      width',logfileunit)
   call print_time(    binwidth,'requested bin      width')
   call print_time(halfbinwidth,'implied     halfbinwidth',logfileunit)
   call print_time(halfbinwidth,'implied     halfbinwidth')
endif

end subroutine Namelist2Times


!======================================================================


subroutine DefineTimeBins()
! Determine temporal bin characteristics.
! The user input is not guaranteed to align on bin centers.
! So -- we will assume the start time is correct and take strides till we
! get past the last time of interest.
! Nepochs will be the total number of time intervals of the period requested.

! Global variables read in this routine:
! integer,         intent(in)  :: max_num_bins
! type(time_type), intent(in)  :: beg_time, end_time  ! of the particular bin
! type(time_type), intent(in)  :: binsep, halfbinwidth

! Global variables set in this routine:
! type(time_type), intent(out) :: TimeMin, TimeMax
! integer,         intent(out) :: Nepochs
! type(time_type), intent(out), dimension(  :) :: bin_center
! type(time_type), intent(out), dimension(:,:) :: bin_edges
! real(digits12),  intent(out), dimension(  :) :: epoch_center
! real(digits12),  intent(out), dimension(:,:) :: epoch_edges
! integer,         intent(out), dimension(  :) :: obs_used_in_epoch

integer :: iepoch, seconds, days

TimeMin  = beg_time
NepochLoop : do iepoch = 1,max_num_bins
   Nepochs = iepoch
   TimeMax = TimeMin + binsep
   if ( TimeMax > end_time ) exit NepochLoop
   TimeMin = TimeMax
enddo NepochLoop

write(string1,*)'Requesting ',Nepochs,' assimilation periods.'
call error_handler(E_MSG,'DefineTimeBins',string1)

allocate(   bin_center(Nepochs),    bin_edges(2,Nepochs), &
         epoch_center(Nepochs), epoch_edges(2,Nepochs), &
    obs_used_in_epoch(Nepochs) )

obs_used_in_epoch = 0

! Explicitly set first epoch - the rest will be predicated on this
! There is a chance the first bin center is at T=0, which is nigh
! impossible to subtract half a bin width from ...
! Since both the bin center and bin width have to be positive, the
! end of the bin MUST be a positive time, so it really pertains
! only to the first bin leading edge.

iepoch = 1
bin_center( iepoch)    = beg_time
bin_edges(2,iepoch)    = beg_time + halfbinwidth

if (beg_time <= halfbinwidth ) then
   bin_edges(1,iepoch) = set_time(0,0)
else
   bin_edges(1,iepoch) = beg_time - halfbinwidth + set_time(1,0)
endif

call get_time(bin_center(iepoch),seconds,days)
epoch_center(iepoch) = days + seconds/86400.0_digits12

call get_time(bin_edges(1,iepoch),seconds,days)
epoch_edges(1,iepoch) = days + seconds/86400.0_digits12

call get_time(bin_edges(2,iepoch),seconds,days)
epoch_edges(2,iepoch) = days + seconds/86400.0_digits12

BinLoop : do iepoch = 2,Nepochs

   bin_center( iepoch) = bin_center(iepoch-1) + binsep
   bin_edges(1,iepoch) = bin_center(iepoch) - halfbinwidth + set_time(1,0)
   bin_edges(2,iepoch) = bin_center(iepoch) + halfbinwidth

   call get_time(bin_center(iepoch),seconds,days)
   epoch_center(iepoch) = days + seconds/86400.0_digits12

   call get_time(bin_edges(1,iepoch),seconds,days)
   epoch_edges(1,iepoch) = days + seconds/86400.0_digits12

   call get_time(bin_edges(2,iepoch),seconds,days)
   epoch_edges(2,iepoch) = days + seconds/86400.0_digits12

enddo BinLoop

if ( verbose ) then
   do iepoch = 1,Nepochs
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

TimeMin = bin_edges(1,      1) ! minimum time of interest
TimeMax = bin_edges(2,Nepochs) ! maximum time of interest

if ( verbose ) then
   call print_time(TimeMin,'First time of interest',logfileunit)
   call print_time(TimeMax,'Last  time of interest',logfileunit)
   call print_time(TimeMin,'First time of interest')
   call print_time(TimeMax,'Last  time of interest')
   write(logfileunit,*)
   write(     *     ,*)
   call print_date(TimeMin,'First date of interest',logfileunit)
   call print_date(TimeMax,'Last  date of interest',logfileunit)
   call print_date(TimeMin,'First date of interest')
   call print_date(TimeMax,'Last  date of interest')
endif

end subroutine DefineTimeBins


!======================================================================

!> Supported scenarios:
!
!> 1) they are taking the default - 11 standard layerMiddles
!> 2) they are specifying their own layerMiddles
!> 3) they are specifying their own layer interfaces (layerEdges)

subroutine setHeightLevels(layerMiddles, layerEdges, nLayers)

! Compatibility issues:
! Historically, it was only possible to specify the layer middles.

real(r8), dimension(:), intent(inout) :: layerMiddles
real(r8), dimension(:), intent(inout) :: layerEdges
integer,                intent(  out) :: nLayers

! FIXME ... do these NEED to be sorted ...

if (any(layerMiddles /= MISSING_R8) .and. any(layerEdges /= MISSING_R8)) then

   ! user has specified both ... do we want to support that ... no.

   write(string1,*)'Specify "hlevel" or "hlevel_edges" but not both.'
   write(string2,*)'Check your input.nml obs_diag_nml settings.'
   call error_handler(E_ERR, 'setHeightLevels', string1, source, text2=string2)

elseif (all(layerMiddles == MISSING_R8) .and. all(layerEdges == MISSING_R8)) then

   ! default case -- specify the standard layer midpoints and define
   ! the layer edges to be halfway between each midpoint. The top
   ! and bottom edges are inferred.

   layerMiddles(1:11) = (/ 1000, 2000, 3000,  4000,  5000, 6000, &
                         7000, 8000, 9000, 10000, 11000        /)
   nLayers = Rlevels2edges(layerMiddles, layerEdges)
   layerEdges(1:nLayers+1) = sort(layerEdges(1:nLayers+1))

elseif ( any(layerEdges /= MISSING_R8) ) then

   ! user specified edges and we have to find the midpoints

   NLayers = Redges2levels(layerEdges, layerMiddles)

else

   ! user specified midpoints and we have to find the edges

   nLayers = Rlevels2edges(layerMiddles, layerEdges)
   layerEdges(1:nLayers+1) = sort(layerEdges(1:nLayers+1))

endif

end subroutine setHeightLevels


!======================================================================


subroutine setPressureLevels(layerMiddles, layerEdges, nLayers)

! Supported scenarios:
!
! 1) they are taking the default - 11 standard layerMiddles
! 2) they are specifying their own layerMiddles
! 3) they are specifying their own layer interfaces (layerEdges)
!
! Compatibility issues:
! Historically, it was only possible to specify the layer middles.
!

real(r8), dimension(:), intent(inout) :: layerMiddles
real(r8), dimension(:), intent(inout) :: layerEdges
integer,                intent(  out) :: nLayers

if (any(layerMiddles /= MISSING_R8) .and. any(layerEdges /= MISSING_R8)) then

   ! user has specified both ... do we want to support that ... no.

   write(string1,*)'Specify "plevel" or "plevel_edges" but not both.'
   write(string2,*)'Check your input.nml obs_diag_nml settings.'
   call error_handler(E_ERR, 'setPressureLevels', string1, source, text2=string2)

elseif (all(layerMiddles == MISSING_R8) .and. all(layerEdges == MISSING_R8)) then

   ! default case -- specify the standard layer midpoints and define
   ! the layer edges to be halfway between each midpoint. The top
   ! and bottom edges are inferred.

   layerMiddles(1:11) = (/ 1000, 925, 850, 700, 500, 400, 300, 250, 200, 150, 100 /)
   nLayers = Rlevels2edges(layerMiddles, layerEdges)

elseif ( any(layerEdges /= MISSING_R8) ) then

   ! user specified edges and we have to find the midpoints

   NLayers = Redges2levels(layerEdges, layerMiddles)

else

   ! user specified midpoints and we have to find the edges

   nLayers = Rlevels2edges(layerMiddles, layerEdges)

endif

end subroutine setPressureLevels


!======================================================================


subroutine setModelLevels(layerMiddles, layerEdges, nLayers)

! Supported scenarios:
!
! 1) they are taking the default - 11 standard layerMiddles
! 2) they are specifying their own layerMiddles
! 3) they are specifying their own layer interfaces (layerEdges)
!
! Compatibility issues:
! Historically, it was only possible to specify the layer middles.
!

real(r8), dimension(:), intent(inout) :: layerMiddles
real(r8), dimension(:), intent(inout) :: layerEdges
integer,                intent(  out) :: nLayers

! FIXME ... do these NEED to be sorted ...

if (any(layerMiddles /= MISSING_R8) .and. any(layerEdges /= MISSING_R8)) then

   ! user has specified both ... do we want to support that ... no.

   write(string1,*)'Specify "mlevel" or "mlevel_edges" but not both.'
   write(string2,*)'Check your input.nml obs_diag_nml settings.'
   call error_handler(E_ERR, 'setModelLevels', string1, source, text2=string2)

elseif (all(layerMiddles == MISSING_R8) .and. all(layerEdges == MISSING_R8)) then

   ! default case -- specify the standard layer midpoints and define
   ! the layer edges to be halfway between each midpoint. The top
   ! and bottom edges are inferred.

   layerMiddles(1:11) = (/ (i,i=1,11) /)   ! set model levels to indices

   nLayers = Rlevels2edges(layerMiddles, layerEdges)

elseif ( any(layerEdges /= MISSING_R8) ) then

   ! user specified edges and we have to find the midpoints

   NLayers = Redges2levels(layerEdges, layerMiddles)

else

   ! user specified midpoints and we have to find the edges

   nLayers = Rlevels2edges(layerMiddles, layerEdges)

endif

end subroutine setModelLevels


!======================================================================

!> Set the min and max location_types for each region

subroutine DefineRegions()

! globally-scoped variables set by this routine:
! xlim1, xlim2
! ylim1, ylim2
! reg_names
! Nregions   may be modified
!----------------------------------------------------------------------

integer  :: i
real(r8) :: lon

! Provide some level of backward compatiblilty ... if Nregions is more
! than 1, the assumption is that everything is set and no action needs
! to be taken. If it is less than 1, we are trying to provide some
! 'the old' defaults.

! Set the min and max location_types for each region

if (Nregions < 1) then
   Nregions     = 4
   lonlim1(1:4) = (/   0.0_r8,   0.0_r8,   0.0_r8, 235.0_r8 /)
   lonlim2(1:4) = (/ 360.0_r8, 360.0_r8, 360.0_r8, 295.0_r8 /)
   latlim1(1:4) = (/  20.0_r8, -80.0_r8, -20.0_r8,  25.0_r8 /)
   latlim2(1:4) = (/  80.0_r8, -20.0_r8,  20.0_r8,  55.0_r8 /)
   reg_names(1) = 'Northern Hemisphere '
   reg_names(2) = 'Southern Hemisphere '
   reg_names(3) = 'Tropics             '
   reg_names(4) = 'North America       '
else

   ! If lonlim2 < lonlim1 ... we add 360.0 to lonlim2 in anticipation
   ! that you are going 'the long way around'
   ! e.g.    300E to  25E  is a perfectly valid domain ... which
   ! is also 300E to 385E

   do i = 1,Nregions
      if (lonlim2(i) < lonlim1(i)) lonlim2(i) = lonlim2(i) + 360.0_r8
   enddo

endif

! only set the horizontal limits; vertical will be ignored.
! modulo() return is always positive, even if input is negative.
! this is not true for mod()
do i = 1, Nregions
   lon = modulo(lonlim1(i), 360.0_r8)
   min_loc(i) = set_location(lon, latlim1(i), 0.0_r8, VERTISUNDEF)
   lon = modulo(lonlim2(i), 360.0_r8)
   max_loc(i) = set_location(lon, latlim2(i), 0.0_r8, VERTISUNDEF)
enddo

! Print a summary if desired. There can be a lot of levels so I cannot
! predict if they will fit in a string to be sent to the error_handler
if ( verbose ) then

   write(     *     ,*)'pressure levels     = ',plevel(      1:Nplevels)
   write(     *     ,*)'pressure interfaces = ',plevel_edges(1:Nplevels+1)
   write(     *     ,*)'height   levels     = ',hlevel(      1:Nhlevels)
   write(     *     ,*)'height   interfaces = ',hlevel_edges(1:Nhlevels+1)
   write(     *     ,*)'model    levels     = ',mlevel(      1:Nmlevels)
   write(     *     ,*)'model    interfaces = ',mlevel_edges(1:Nmlevels+1)

   write(logfileunit,*)'pressure levels     = ',plevel(      1:Nplevels)
   write(logfileunit,*)'pressure interfaces = ',plevel_edges(1:Nplevels+1)
   write(logfileunit,*)'height   levels     = ',hlevel(      1:Nhlevels)
   write(logfileunit,*)'height   interfaces = ',hlevel_edges(1:Nhlevels+1)
   write(logfileunit,*)'model    levels     = ',mlevel(      1:Nmlevels)
   write(logfileunit,*)'model    interfaces = ',mlevel_edges(1:Nmlevels+1)

   do i = 1,Nregions
      write(     *     ,'(''Region '',i02,1x,a32,'' (WESN): '',4(f10.4,1x))') &
             i, reg_names(i), lonlim1(i),lonlim2(i),latlim1(i),latlim2(i)
      write(logfileunit,'(''Region '',i02,1x,a32,'' (WESN): '',4(f10.4,1x))') &
             i, reg_names(i), lonlim1(i),lonlim2(i),latlim1(i),latlim2(i)
   enddo
endif

end subroutine DefineRegions


!======================================================================


subroutine CountTrustedObsTypes()

! Count up the number of 'trusted' observations types from those
! specified in the input namelist 'trusted_obs' argument.
! Check to make sure trusted observations desired are supported.
!
! trusted_obs  is a namelist variable ... global storage
! trusted_list is in global storage and is SET by this routine.
! num_trusted  is in global storage and is SET by this routine.

integer :: i
logical :: matched
character(len=stringlength) :: possible_obs_type

! Loop over all user input candidates for 'trusted' observations.
! Check each candidate against list of known observation names.

num_trusted = 0
InputList : do i = 1,MaxTrusted

   if (trim(trusted_obs(i)) == 'null') exit InputList

   matched = .false.

   PossibleList : do ikind = 1,num_obs_types

      possible_obs_type = get_name_for_type_of_obs(ikind)
      if (trim(trusted_obs(i)) == trim(possible_obs_type)) then
         matched = .true.
         exit PossibleList
   endif

   enddo PossibleList

   ! At this point, we either have a trusted observation name or
   ! user is trying to trust an observation we do not know about ...
   ! perhaps a mis-spelling in the user input list?
   ! If the TERMLEVEL is set to 1 ... the WARNING will be fatal.

   if (matched) then
       num_trusted = num_trusted + 1
       trusted_list(num_trusted) = trim(trusted_obs(i))
      else
      write(string1,*)'trusted_obs "',trim(trusted_obs(i)), &
                    & '" is not a supported observation type.'
      write(string2,*)'to get past this warning, set utilities_nml:TERMLEVEL > 1 and rerun.'
      call error_handler(E_WARN,'CountTrustedObsTypes', string1, source, text2=string2)
      endif

enddo InputList

! There is a (remote) possibility someone wants to trust more than MaxTrusted
! different observation types. We will issue a WARNING, but carry on.
! If the TERMLEVEL is set to 1 ... the WARNING will be fatal.

if (num_trusted == MaxTrusted) then
   write(string1,*)'There are ',num_trusted,' "trusted" observation types.'
   write(string2,*)'This is the maximum allowed unless you increase "MaxTrusted" and recompile.'
   write(string3,*)'to get past this warning, set utilities_nml:TERMLEVEL > 1 and rerun.'
   call error_handler(E_WARN, 'CountTrustedObsTypes', string1, &
                         source, text2=string2, text3=string3)
endif

! If there are some trusted observation types, list them.
! If not, just say so and keep going.

if (num_trusted > 0 ) then
   write(string1,*)'There are ',num_trusted, &
                   ' "trusted" observation types, they are:'
   call error_handler(E_MSG, 'CountTrustedObsTypes', string1)
   do i = 1,num_trusted
      call error_handler(E_MSG, 'CountTrustedObsTypes', trusted_list(i))
   enddo
else
   write(string1,*)'There are no "trusted" observation types.'
   call error_handler(E_MSG, 'CountTrustedObsTypes', string1)
endif

end subroutine CountTrustedObsTypes


!======================================================================


subroutine  SetScaleFactors()

! The surface pressure in the obs_sequence is in Pa, we want to convert
! from Pa to hPa for plotting. The specific humidity is a similar thing.
! In the obs_sequence file, the units are kg/kg, we want to plot
! in the g/kg world...
!
! Gloval-scope variables used in this routine:
! real(r8), dimension(:), intent(inout) :: scale_factor
! integer,                intent(in)    :: logfileunit

! If kind_surface_pressure or ... does not exist, we are in trouble here.
! the scale_factor should be defined to reflect the type, which are not
! guaranteed to be numbered sequentially ... vortices 81, for example

character(len=stringlength) :: obs_string
integer :: ivar

scale_factor = 1.0_r8

! The scale_factor array is dimensioned from obs_kind_mod:num_obs_types

do ivar = 1,SIZE(scale_factor)

   obs_string = obs_type_strings(ivar)

   if ( index(obs_string,'SURFACE_PRESSURE') > 0 ) &
          scale_factor(ivar) = 0.01_r8

   if ( index(obs_string,'SPECIFIC_HUMIDITY') > 0 ) &
          scale_factor(ivar) = 1000.0_r8

   ! Somehow, we should plot statistics on the dBZ scale for these ...
   ! scale_factor(QTY_RADAR_REFLECTIVITY) = 10log10(z)

   ! The scaling is summarized in WriteNetCDF if the verbose option is chosen.

enddo

end subroutine SetScaleFactors


!======================================================================
!> We need to know the time of the first and last observations in the sequence,
!> primarily just to see if they intersect the desired Epoch window.
!> We also record these times so we can report the first/last times of all
!> observations in all the obs_seq files.

subroutine GetFirstLastObs(my_fname, my_sequence, my_obs1, my_obsN, &
                   my_seqT1, my_seqTN, my_AllseqT1, my_AllseqTN)

character(len=*),        intent(in)    :: my_fname
type(obs_sequence_type), intent(in)    :: my_sequence
type(obs_type),          intent(out)   :: my_obs1
type(obs_type),          intent(out)   :: my_obsN
type(time_type),         intent(out)   :: my_seqT1
type(time_type),         intent(out)   :: my_seqTN
type(time_type),         intent(inout) :: my_AllseqT1 ! ALL observation sequences
type(time_type),         intent(inout) :: my_AllseqTN ! ALL observation sequences

type(obs_def_type) :: obs_def

logical, SAVE :: first_time = .true.

if ( .not. get_first_obs(my_sequence, my_obs1) ) then
   call error_handler(E_ERR,'GetFirstLastObs','No first observation in '//trim(my_fname), &
   source)
endif
call get_obs_def(my_obs1,   obs_def)
my_seqT1 = get_obs_def_time(obs_def)

if ( .not. get_last_obs(my_sequence, my_obsN) ) then
   call error_handler(E_ERR,'GetFirstLastObs','No last observation in '//trim(my_fname), &
   source)
endif
call get_obs_def(my_obsN,   obs_def)
my_seqTN = get_obs_def_time(obs_def)

! Capture a little information to assist in an error message if the
! namelist input does not intersect the observation sequence file.

if ( first_time ) then
   my_AllseqT1 = my_seqT1
   my_AllseqTN = my_seqTN
   first_time  = .false.
else
   if (my_seqT1 < my_AllseqT1) my_AllseqT1 = my_seqT1
   if (my_seqTN > my_AllseqTN) my_AllseqTN = my_seqTN
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

! If the last observation is before the period of interest, move on.

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


function GetEnsSize()
!-----------------------------------------------------------------------
!  Loop over all the metadata to count the number of ensemble members
!  available in the observation sequence file. We need this count to
!  allocate space for the rank histogram information. Since the rank
!  histogram will be created for the priors only ...
!-----------------------------------------------------------------------

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


subroutine SetIndices( obs_index, org_qc_index, dart_qc_index,     &
                       prior_mean_index,   posterior_mean_index,   &
                       prior_spread_index, posterior_spread_index)

! There are many 'copy' indices that need to be set from the obs_sequence
! metadata. Some are required, some are optional.

integer, intent(out) :: obs_index
integer, intent(out) :: org_qc_index
integer, intent(out) :: dart_qc_index
integer, intent(out) :: prior_mean_index
integer, intent(out) :: posterior_mean_index
integer, intent(out) :: prior_spread_index
integer, intent(out) :: posterior_spread_index

! Using 'seq', 'ens_size', and ens_copy_index from global scope

integer :: i, ens_count
character(len=metadatalength) :: metadata

obs_index              = -1
prior_mean_index       = -1
posterior_mean_index   = -1
prior_spread_index     = -1
posterior_spread_index = -1
org_qc_index           = -1
dart_qc_index          = -1

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

! Hmnnn ... sometimes the first QC field is 'Quality Control' or 'NCEP QC index'
! to me ... they mean the same thing.

QCMetaDataLoop : do i=1, get_num_qc(seq)
   if(index(  get_qc_meta_data(seq,i), 'Quality Control'          ) > 0) &
                    org_qc_index = i
   if(index(  get_qc_meta_data(seq,i), 'NCEP QC index'            ) > 0) &
                    org_qc_index = i
   if(index(  get_qc_meta_data(seq,i), 'DART quality control'     ) > 0) &
                   dart_qc_index = i
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
   call error_handler(E_MSG, 'SetIndices', string1, &
              source, text2=string2, text3=string3)
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


subroutine CheckVertical(obslocation, flavor)

! It is not allowed to have one observation flavor exist on multiple
! types of vertical coordinate systems. If the flavor has been
! initialized (ob_defining_vert(flavor) > 0) then it is an
! error to change the type of vertical coordinate system
! {tracked in which_vert()}.

type(location_type),   intent(in) :: obslocation
integer(i8),           intent(in) :: flavor

!  integer, dimension(:), intent(in) :: which_vert       GLOBAL SCOPE
!  integer, dimension(:), intent(in) :: ob_defining_vert GLOBAL SCOPE

if (is_vertical(obslocation, "HEIGHT")        .and. &
   (      which_vert(flavor) /= VERTISHEIGHT) .and. &
   (ob_defining_vert(flavor)  >     0       ) ) then
        write(string1,'(''obs '', i8, '' type '', i3, &
        & '' changing from '', i2, '' to height - def by obs '',i8)') &
        keys(obsindex), flavor, which_vert(flavor), ob_defining_vert(flavor)
   call error_handler(E_WARN,'CheckVertical',string1,source)
endif

end subroutine CheckVertical


!======================================================================


function ParseLevel(obslocation, obslevel, flavor)
! returns the index of the level closest to 'obslevel'
! There are three independent 'level' arrays
type(location_type),   intent(in   ) :: obslocation
real(r8),              intent(inout) :: obslevel
integer(i8),           intent(in   ) :: flavor
integer :: ParseLevel

!  integer, dimension(:), intent(inout) :: which_vert        GLOBAL SCOPE
!  logical, dimension(:), intent(inout) :: ob_defining_vert  GLOBAL SCOPE

if(is_vertical(obslocation, "SURFACE")) then
   obslevel           = 1
   ParseLevel         = 1
   which_vert(flavor) = VERTISSURFACE

elseif(is_vertical(obslocation, "LEVEL")) then
   ! at one point, negative (model) levels were used to indicate
   ! surface observations - ancient history ...
   if (obslevel < 1.0_r8 ) obslevel = 1.0_r8  ! TJH negative model level
   ParseLevel         = find_my_level(obslevel, VERTISLEVEL)
   which_vert(flavor) = VERTISLEVEL

   if ( (1 == 2) .and. (ParseLevel > 0) .and. (ParseLevel < MaxLevels)) &
         write(*,301) keys(obsindex), obsindex, &
         mlevel_edges(ParseLevel), obslevel, mlevel_edges(ParseLevel+1), &
         mlevel(ParseLevel), ParseLevel, 'level'

elseif(is_vertical(obslocation, "PRESSURE")) then
   ParseLevel         = find_my_level(obslevel, VERTISPRESSURE)
   which_vert(flavor) = VERTISPRESSURE

   if ( (1 == 2) .and. (ParseLevel > 0) .and. (ParseLevel < MaxLevels)) &
         write(*,300) keys(obsindex), obsindex, &
         plevel_edges(ParseLevel), obslevel, plevel_edges(ParseLevel+1), &
         plevel(ParseLevel), ParseLevel, 'pressure'

elseif(is_vertical(obslocation, "HEIGHT")) then
   ParseLevel         = 1
   which_vert(flavor) = VERTISHEIGHT

   if ( (1 == 2) .and. (ParseLevel > 0) .and. (ParseLevel < MaxLevels)) &
         write(*,300) keys(obsindex), obsindex, &
         hlevel_edges(ParseLevel), obslevel, hlevel_edges(ParseLevel+1), &
         hlevel(ParseLevel), ParseLevel, 'height'

elseif(is_vertical(obslocation, "SCALE_HEIGHT")) then
   ParseLevel         = find_my_level(obslevel, VERTISSCALEHEIGHT)
   which_vert(flavor) = VERTISSCALEHEIGHT

   if ( (1 == 2) .and. (ParseLevel > 0) .and. (ParseLevel < MaxLevels)) &
         write(*,300) keys(obsindex), obsindex, &
         hlevel_edges(ParseLevel), obslevel, hlevel_edges(ParseLevel+1), &
         hlevel(ParseLevel), ParseLevel, 'height'

elseif(is_vertical(obslocation, "UNDEFINED")) then
   obslevel           = 1
   ParseLevel         = 1
   which_vert(flavor) = VERTISUNDEF

else
   call error_handler(E_ERR,'ParseLevel','Vertical coordinate not recognized', source)
endif

ob_defining_vert(flavor) = keys(obsindex)

 300 format('obskey ',i7,' # ',i7,' [',3(f10.4,1x),'] mid == ',f10.4,' level ',i3,1x,a)
 301 format('obskey ',i7,' # ',i7,' [',3(f10.4,1x),'] mid == ',i10  ,' level ',i3,1x,a)

end function ParseLevel


!======================================================================


function InnovZscore(obsval, prmean, prspred, errvar, qcval, qcmaxval)

! This function tries to get a handle on the magnitude of the innovations.
! If the ratio of the observation to the prior mean is 'big', it is an outlier.
! If the prior mean cannot be calculated (i.e. is missing) we put it in the
! last 'bin' of the crude histogram. This is pretty much a 'z' score in the
! statistical sense.
!
! In Jan of 2014 I ran into a special case. We are performing a perfect model
! experiment - perturbing a single state and then assimilating. We wanted
! to compare against the true observation value (errvar = 0.0) and for the
! first time step, there is no ensemble spread either. In this case the denom
! was really zero and the calculation blew up. Since we only use this to track
! how far apart these things are, we can restrict the distance to the worst-case
! scenario ... the last bin.

real(r8)             :: InnovZscore
real(r8), intent(in) :: obsval, prmean, prspred, errvar
integer,  intent(in) :: qcval, qcmaxval

real(r8) :: numer, denom

InnovZscore = real(MaxSigmaBins,r8) ! worst-case ... really far apart

if ( qcval <= qcmaxval ) then ! QC indicates a valid obs
   numer = abs(prmean - obsval)
   denom = sqrt( prspred**2 + errvar )

   ! At worst, the InnovZscore can be 'MaxSigmaBins' i.e. 100
   ! protect against dividing by a really small number
   ! if numer/denom < 100 then go ahead and calculate
   ! if numer < 100 * denom ...

   if ( numer < InnovZscore*denom ) then
      InnovZscore = numer / denom
   endif

endif

end function InnovZscore


!======================================================================


function Rank_Histogram(copyvalues, obs_index, &
                    error_variance ) result(rank)

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


subroutine CountDartQC_4D(inqc, iepoch, ilevel, iregion, itype, prior, poste, &
                          posterior_mean, uqc)

integer,         intent(in)    :: inqc
integer,         intent(in)    :: iepoch
integer,         intent(in)    :: ilevel
integer,         intent(in)    :: iregion
integer(i8),     intent(in)    :: itype
type(TLRV_type), intent(inout) :: prior
type(TLRV_type), intent(inout) :: poste
real(r8),        intent(in), optional :: posterior_mean
integer,         intent(in), optional :: uqc

integer :: myqc

if (present(uqc)) then
   myqc = maxval((/ inqc, uqc /))
else
   myqc = inqc
endif

if (        myqc == 0 ) then
   call IPE(prior%NDartQC_0(iepoch,ilevel,iregion,itype), 1)
   if (has_posteriors) &
      call IPE(poste%NDartQC_0(iepoch,ilevel,iregion,itype), 1)

elseif (    myqc == 1 ) then
   call IPE(prior%NDartQC_1(iepoch,ilevel,iregion,itype), 1)
   if (has_posteriors) &
      call IPE(poste%NDartQC_1(iepoch,ilevel,iregion,itype), 1)

elseif (    myqc == 2 ) then
   call IPE(prior%NDartQC_2(iepoch,ilevel,iregion,itype), 1)
   if (has_posteriors) &
      call IPE(poste%NDartQC_2(iepoch,ilevel,iregion,itype), 1)

elseif (    myqc == 3 ) then
   call IPE(prior%NDartQC_3(iepoch,ilevel,iregion,itype), 1)
   if (has_posteriors) &
      call IPE(poste%NDartQC_3(iepoch,ilevel,iregion,itype), 1)

elseif (    myqc == 4 ) then
   call IPE(prior%NDartQC_4(iepoch,ilevel,iregion,itype), 1)
   if (has_posteriors) &
      call IPE(poste%NDartQC_4(iepoch,ilevel,iregion,itype), 1)

elseif (    myqc == 5 ) then
   call IPE(prior%NDartQC_5(iepoch,ilevel,iregion,itype), 1)
   if (has_posteriors) &
      call IPE(poste%NDartQC_5(iepoch,ilevel,iregion,itype), 1)

elseif (    myqc == 6 ) then
   call IPE(prior%NDartQC_6(iepoch,ilevel,iregion,itype), 1)
   if (has_posteriors) &
      call IPE(poste%NDartQC_6(iepoch,ilevel,iregion,itype), 1)

elseif (    myqc == QC_OUTLIER ) then
   call IPE(prior%NDartQC_7(iepoch,ilevel,iregion,itype), 1)

   if (present(posterior_mean) .and. has_posteriors) then
      if (posterior_mean == MISSING_R8) then
         ! ACTUALLY A FAILED POSTERIOR FORWARD OPERATOR - ambiguous case
         call IPE(poste%NDartQC_2(iepoch,ilevel,iregion,itype), 1)
      else
         call IPE(poste%NDartQC_7(iepoch,ilevel,iregion,itype), 1)
      endif
   endif

elseif (    myqc == 8 ) then
   call IPE(prior%NDartQC_8(iepoch,ilevel,iregion,itype), 1)
   if (has_posteriors) &
      call IPE(poste%NDartQC_8(iepoch,ilevel,iregion,itype), 1)

endif

end subroutine CountDartQC_4D


!----------------------------------------------------------------------
!> This function simply accumulates the appropriate sums.
!> The normalization occurrs after all the data has been read, naturally.

subroutine Bin4D(iqc, iepoch, ilevel, iregion, flavor, trusted, &
             obsval,  obserrvar,  prmean,  prsprd,  pomean,  posprd, rank, &
               uobs, uobserrvar, uprmean, uprsprd, upomean, uposprd, uqc)

! The 'prior' and 'poste' structures are globally scoped.
!
! Wind measurements are vector quantities - so we are collapsing them to
! scalar speed for the bias. The optional arguments specify the U components
! while the mandatory arguments specify the V components.
! Its an 'all-or-nothing' optional argument situation.

! If you are verifying the ensemble against imperfect (real) observations,
! it is necessary to account for the observation error when computing the
! spread.  Since the observation error is not included as output from
! obs_diag, it is not possible to compute this quantity now.

integer,  intent(in)           :: iqc, iepoch, ilevel, iregion
integer(i8), intent(in)        :: flavor
logical,  intent(in)           :: trusted
real(r8), intent(in)           :: obsval,  obserrvar,  prmean,  prsprd,  pomean,  posprd
integer,  intent(in)           :: rank
real(r8), intent(in), optional ::   uobs, uobserrvar, uprmean, uprsprd, upomean, uposprd
integer,  intent(in), optional :: uqc

real(r8) :: priorsqerr      ! PRIOR     Squared Error
real(r8) :: priorbias       ! PRIOR     simple bias
real(r8) :: postsqerr       ! POSTERIOR Squared Error
real(r8) :: postbias        ! POSTERIOR simple bias

real(r8) :: prior_variance
real(r8) :: prior_varianceplus
real(r8) :: posterior_variance
real(r8) :: posterior_varianceplus

real(r8) :: priormean, postmean, obsmean
integer  :: myrank, prior_qc, posterior_qc

logical, dimension(7) :: optionals

prior_qc     = iqc
posterior_qc = iqc

! There is an ambiguous case wherein the prior is rejected (DART QC == 7)
! and the posterior forward operator fails (DART QC == 2). In this case,
! the DART_QC reflects the fact the prior was rejected - HOWEVER -
! the posterior mean,spread are set to MISSING.

if (prior_qc == QC_OUTLIER .and. pomean == MISSING_R8) posterior_qc = QC_PO_FOP_FAIL

! Check to see if we are creating wind speeds from U,V components

optionals = (/ present(uobs), present(uobserrvar), present(uprmean), &
               present(uprsprd), present(upomean), present(uposprd), present(uqc) /)

if ( all(optionals) ) then

   ! the wind QC is only as good as the worst of the U,V QCs
   prior_qc     = maxval( (/ iqc, uqc /) )

   ! If either the U or V is ambiguous, the wind is ambiguous
   if     (uqc == QC_OUTLIER .and. upomean == MISSING_R8) then
      posterior_qc = QC_PO_FOP_FAIL
   elseif (iqc == QC_OUTLIER .and.  pomean == MISSING_R8) then
      posterior_qc = QC_PO_FOP_FAIL
   else
      posterior_qc = maxval( (/ iqc, uqc /) )
   endif

   priorsqerr     = (prmean - obsval)**2 + (uprmean - uobs)**2
   postsqerr      = (pomean - obsval)**2 + (upomean - uobs)**2

   ! This calculation is the bias in the wind speed (F-O)
   obsmean        = sqrt(obsval**2 + uobs**2)
   priormean      = sqrt(prmean**2 + uprmean**2)
   postmean       = sqrt(pomean**2 + upomean**2)
   priorbias      = priormean - obsmean
   postbias       = postmean  - obsmean

   ! convert standard deviations to variances and add
   prior_variance          = prsprd**2 + uprsprd**2
   posterior_variance      = posprd**2 + uposprd**2
   prior_varianceplus      = prsprd**2 + obserrvar + uprsprd**2 + uobserrvar
   posterior_varianceplus  = posprd**2 + obserrvar + uposprd**2 + uobserrvar

   ! If we are working with 'horizontal winds', we do not have enough
   ! information to recreate the appropriate rank histogram. We will
   ! set the rank to a bogus value which will ensure it does not get
   ! counted.
   myrank = -99

elseif ( any(optionals) ) then
   call error_handler(E_ERR,'Bin4D','wrong number of optional arguments', source)
else

   priorsqerr     = (prmean - obsval)**2
   postsqerr      = (pomean - obsval)**2

   obsmean        = obsval
   priormean      = prmean
   postmean       = pomean
   priorbias      = prmean - obsval
   postbias       = pomean - obsval

   ! convert standard deviations to variances and add
   prior_variance          = prsprd**2
   posterior_variance      = posprd**2
   prior_varianceplus      = prsprd**2 + obserrvar
   posterior_varianceplus  = posprd**2 + obserrvar

   myrank = rank
endif

!----------------------------------------------------------------------
! The rank histogram binning is a bit of a peculiar situation.
! Only the prior is of interest ... so DART QCs of 0 1 2 3 are 'good'.
! There is some debate about whether we should be considering the
! 'outlier' observations (DART QC == 7), so that is namelist controlled.
!----------------------------------------------------------------------

if (     (myrank > 0) .and. create_rank_histogram ) then
   if ( any(prior_qc == hist_qcs(1:numqcvals) ) )  &
      call IPE(prior%hist_bin(iepoch,ilevel,iregion,flavor,myrank), 1)
endif

!----------------------------------------------------------------------
! Track the number of possible observations
!----------------------------------------------------------------------

call IPE(prior%Nposs(iepoch,ilevel,iregion,flavor), 1)
if (has_posteriors) &
   call IPE(poste%Nposs(iepoch,ilevel,iregion,flavor), 1)

!----------------------------------------------------------------------
! Select which set of qcs are valid and accrue everything
!----------------------------------------------------------------------

if ( trusted ) then
   call IPE(prior%Ntrusted(iepoch,ilevel,iregion,flavor), 1)
   if (has_posteriors) &
      call IPE(poste%Ntrusted(iepoch,ilevel,iregion,flavor), 1)
endif

! Accrue the PRIOR quantities
if ((      trusted .and.  any(trusted_prior_qcs == prior_qc)) .or. &
    (.not. trusted .and.  any(   good_prior_qcs == prior_qc))) then
   call IPE(prior%Nused(      iepoch,ilevel,iregion,flavor),      1    )
   call RPE(prior%observation(iepoch,ilevel,iregion,flavor), obsmean   )
   call RPE(prior%ens_mean(   iepoch,ilevel,iregion,flavor), priormean )
   call RPE(prior%bias(       iepoch,ilevel,iregion,flavor), priorbias )
   call RPE(prior%rmse(       iepoch,ilevel,iregion,flavor), priorsqerr)
   call RPE(prior%spread(     iepoch,ilevel,iregion,flavor), prior_variance)
   call RPE(prior%totspread(  iepoch,ilevel,iregion,flavor), prior_varianceplus)
else
   call IPE(prior%NbadDartQC(iepoch,ilevel,iregion,flavor),       1    )
endif

! Accrue the POSTERIOR quantities
if (has_posteriors) then
   if ((      trusted .and.  any(trusted_poste_qcs == posterior_qc)) .or. &
       (.not. trusted .and.  any(   good_poste_qcs == posterior_qc))) then
      call IPE(poste%Nused(      iepoch,ilevel,iregion,flavor),      1   )
      call RPE(poste%observation(iepoch,ilevel,iregion,flavor), obsmean  )
      call RPE(poste%ens_mean(   iepoch,ilevel,iregion,flavor), postmean )
      call RPE(poste%bias(       iepoch,ilevel,iregion,flavor), postbias )
      call RPE(poste%rmse(       iepoch,ilevel,iregion,flavor), postsqerr)
      call RPE(poste%spread(     iepoch,ilevel,iregion,flavor), posterior_variance)
      call RPE(poste%totspread(  iepoch,ilevel,iregion,flavor), posterior_varianceplus)
   else
      call IPE(poste%NbadDartQC(iepoch,ilevel,iregion,flavor),       1    )
   endif
endif

end subroutine Bin4D


!======================================================================


subroutine Normalize4Dvars()

! The vertical coordinate system definition is summarized in WriteNetCDF
! if the verbose option is chosen.

if ( verbose ) then
   write(logfileunit,*)'Normalizing time-level-region-variable quantities.'
   write(     *     ,*)'Normalizing time-level-region-variable quantities.'
endif

do ivar   = 1,num_obs_types
do iregion= 1,Nregions
do ilev   = 1,Nlevels
do iepoch = 1,Nepochs

   ! Normalize the priors

   if (  prior%Nused(      iepoch, ilev, iregion, ivar) == 0) then
         prior%observation(iepoch, ilev, iregion, ivar) = MISSING_R4
         prior%ens_mean(   iepoch, ilev, iregion, ivar) = MISSING_R4
         prior%bias(       iepoch, ilev, iregion, ivar) = MISSING_R4
         prior%rmse(       iepoch, ilev, iregion, ivar) = MISSING_R4
         prior%spread(     iepoch, ilev, iregion, ivar) = MISSING_R4
         prior%totspread(  iepoch, ilev, iregion, ivar) = MISSING_R4
   else
         prior%observation(iepoch, ilev, iregion, ivar) = &
         prior%observation(iepoch, ilev, iregion, ivar) / &
         prior%Nused(      iepoch, ilev, iregion, ivar)

         prior%ens_mean(   iepoch, ilev, iregion, ivar) = &
         prior%ens_mean(   iepoch, ilev, iregion, ivar) / &
         prior%Nused(      iepoch, ilev, iregion, ivar)

         prior%bias(       iepoch, ilev, iregion, ivar) = &
         prior%bias(       iepoch, ilev, iregion, ivar) / &
         prior%Nused(      iepoch, ilev, iregion, ivar)

         prior%rmse(       iepoch, ilev, iregion, ivar) = &
    sqrt(prior%rmse(       iepoch, ilev, iregion, ivar) / &
         prior%Nused(      iepoch, ilev, iregion, ivar) )

    ! convert the (pooled) variances back to standard deviations AKA 'spread'
         prior%spread(     iepoch, ilev, iregion, ivar) = &
    sqrt(prior%spread(     iepoch, ilev, iregion, ivar) / &
         prior%Nused(      iepoch, ilev, iregion, ivar) )

    ! convert the (pooled) variances back to standard deviations AKA 'spread'
         prior%totspread(  iepoch, ilev, iregion, ivar) = &
    sqrt(prior%totspread(  iepoch, ilev, iregion, ivar) / &
         prior%Nused(      iepoch, ilev, iregion, ivar) )

   endif

   ! Same thing for the posteriors

   if (  poste%Nused(      iepoch, ilev, iregion, ivar) == 0) then
         poste%observation(iepoch, ilev, iregion, ivar) = MISSING_R4
         poste%ens_mean(   iepoch, ilev, iregion, ivar) = MISSING_R4
         poste%bias(       iepoch, ilev, iregion, ivar) = MISSING_R4
         poste%rmse(       iepoch, ilev, iregion, ivar) = MISSING_R4
         poste%spread(     iepoch, ilev, iregion, ivar) = MISSING_R4
         poste%totspread(  iepoch, ilev, iregion, ivar) = MISSING_R4
   else
         poste%observation(iepoch, ilev, iregion, ivar) = &
         poste%observation(iepoch, ilev, iregion, ivar) / &
         poste%Nused(      iepoch, ilev, iregion, ivar)

         poste%ens_mean( iepoch, ilev, iregion, ivar) = &
         poste%ens_mean( iepoch, ilev, iregion, ivar) / &
         poste%Nused(      iepoch, ilev, iregion, ivar)

         poste%bias(       iepoch, ilev, iregion, ivar) = &
         poste%bias(       iepoch, ilev, iregion, ivar) / &
         poste%Nused(      iepoch, ilev, iregion, ivar)

         poste%rmse(       iepoch, ilev, iregion, ivar) = &
    sqrt(poste%rmse(       iepoch, ilev, iregion, ivar) / &
         poste%Nused(      iepoch, ilev, iregion, ivar) )

    ! convert the (pooled) variances back to standard deviations AKA 'spread'
         poste%spread(     iepoch, ilev, iregion, ivar) = &
    sqrt(poste%spread(     iepoch, ilev, iregion, ivar) / &
         poste%Nused(      iepoch, ilev, iregion, ivar) )

    ! convert the (pooled) variances back to standard deviations AKA 'spread'
         poste%totspread(  iepoch, ilev, iregion, ivar) = &
    sqrt(poste%totspread(  iepoch, ilev, iregion, ivar) / &
         poste%Nused(      iepoch, ilev, iregion, ivar) )

   endif
enddo
enddo
enddo
enddo

! Actually print the histogram of innovations as a function of standard deviation.
write(     *    ,*)
write(nsigmaUnit,*)
write(     *    ,'(''last bin contains all (flagged) bad observations'')')
write(nsigmaUnit,'(''last bin contains all (flagged) bad observations'')')
do i=0,MaxSigmaBins
   if(nsigma(i) /= 0) then
      write(     *    ,'(''(prior) innovations in stdev bin '',i3,'' = '',i10)')i+1,nsigma(i)
      write(nsigmaUnit,'(''(prior) innovations in stdev bin '',i3,'' = '',i10)')i+1,nsigma(i)
   endif
enddo
call close_file(nsigmaUnit)

end subroutine Normalize4Dvars


!======================================================================


subroutine WriteNetCDF(fname)

character(len=*), intent(in) :: fname

integer :: ncid, i, nobs, typesdimlen, io
integer ::  RegionDimID,  RegionVarID
integer ::  MlevelDimID,  MlevelVarID
integer ::  PlevelDimID,  PlevelVarID
integer ::  HlevelDimID,  HlevelVarID
integer ::  SlevelDimID,  UlevelDimID
integer ::    TimeDimID,    TimeVarID
integer ::    CopyDimID,    CopyVarID,  CopyMetaVarID
integer ::   TypesDimID,   TypesVarID, TypesMetaVarID
integer :: PlevIntDimID, PlevIntVarID
integer :: HlevIntDimID, HlevIntVarID
integer :: MlevIntDimID, MlevIntVarID
integer ::  BoundsDimID,  BoundsVarID
integer ::    RankDimID,    RankVarID
integer ::  StringDimID

integer :: TimeBoundsVarID, RegionNamesVarID

character(len=8)      :: crdate      ! needed by F90 DATE_AND_TIME intrinsic
character(len=10)     :: crtime      ! needed by F90 DATE_AND_TIME intrinsic
character(len=5)      :: crzone      ! needed by F90 DATE_AND_TIME intrinsic
integer, dimension(8) :: values      ! needed by F90 DATE_AND_TIME intrinsic

if(.not. byteSizesOK()) then
    call error_handler(E_ERR,'WriteNetCDF', &
   'Compiler does not support required kinds of variables.',source)
endif

call nc_check(nf90_create(path = trim(fname), cmode = nf90_share, &
         ncid = ncid), 'WriteNetCDF', 'create "'//trim(fname)//'"')

!----------------------------------------------------------------------------
! Write Global Attributes
!----------------------------------------------------------------------------

call DATE_AND_TIME(crdate,crtime,crzone,values)
write(string1,'(''YYYY MM DD HH MM SS = '',i4,5(1x,i2.2))') &
               values(1), values(2), values(3), values(5), values(6), values(7)
call nc_check(nf90_put_att(ncid, NF90_GLOBAL, 'creation_date', trim(string1) ), &
           'WriteNetCDF', 'put_att creation_date '//trim(fname))

!  call nc_check(nf90_put_att(ncid, NF90_GLOBAL, 'title', global_meta_data), &
!             'WriteNetCDF', 'put_att title '//trim(fname))
call nc_check(nf90_put_att(ncid, NF90_GLOBAL, 'obs_diag_source', source ), &
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
! write all namelist quantities
!----------------------------------------------------------------------------

call nc_check(nf90_put_att(ncid, NF90_GLOBAL, 'first_bin_center', first_bin_center ), &
           'WriteNetCDF', 'put_att first_bin_center '//trim(fname))
call nc_check(nf90_put_att(ncid, NF90_GLOBAL, 'last_bin_center', last_bin_center ), &
           'WriteNetCDF', 'put_att last_bin_center '//trim(fname))
call nc_check(nf90_put_att(ncid, NF90_GLOBAL, 'bin_separation', bin_separation ), &
           'WriteNetCDF', 'put_att bin_separation '//trim(fname))
call nc_check(nf90_put_att(ncid, NF90_GLOBAL, 'bin_width', bin_width ), &
           'WriteNetCDF', 'put_att bin_width '//trim(fname))
call nc_check(nf90_put_att(ncid, NF90_GLOBAL, 'time_to_skip', time_to_skip ), &
           'WriteNetCDF', 'put_att time_to_skip '//trim(fname))
call nc_check(nf90_put_att(ncid, NF90_GLOBAL, 'max_num_bins', max_num_bins ), &
           'WriteNetCDF', 'put_att max_num_bins '//trim(fname))
call nc_check(nf90_put_att(ncid, NF90_GLOBAL, 'rat_cri', rat_cri ), &
           'WriteNetCDF', 'put_att rat_cri '//trim(fname))
call nc_check(nf90_put_att(ncid, NF90_GLOBAL, 'input_qc_threshold', input_qc_threshold ), &
           'WriteNetCDF', 'put_att input_qc_threshold '//trim(fname))
call nc_check(nf90_put_att(ncid, NF90_GLOBAL, 'lonlim1', lonlim1(1:Nregions) ), &
           'WriteNetCDF', 'put_att lonlim1 '//trim(fname))
call nc_check(nf90_put_att(ncid, NF90_GLOBAL, 'lonlim2', lonlim2(1:Nregions) ), &
           'WriteNetCDF', 'put_att lonlim2 '//trim(fname))
call nc_check(nf90_put_att(ncid, NF90_GLOBAL, 'latlim1', latlim1(1:Nregions) ), &
           'WriteNetCDF', 'put_att latlim1 '//trim(fname))
call nc_check(nf90_put_att(ncid, NF90_GLOBAL, 'latlim2', latlim2(1:Nregions) ), &
           'WriteNetCDF', 'put_att latlim2 '//trim(fname))

do i = 1,num_trusted
   write(string1,'(''trusted_obs_'',i2.2)') i
   call nc_check(nf90_put_att(ncid, NF90_GLOBAL, trim(string1), trim(trusted_list(i))), &
        'WriteNetCDF', 'put_att trusted_list '//trim(fname))
enddo

!----------------------------------------------------------------------------
! write all observation sequence files used
!----------------------------------------------------------------------------

FILEloop : do i = 1, num_input_files

  write(string1,'(''obs_seq_file_'',i5.5)')i
  io = nf90_put_att(ncid, NF90_GLOBAL, trim(string1), trim(obs_sequence_name(i)))
  call nc_check(io, 'WriteNetCDF', 'put_att input file names')

enddo FILEloop

io = nf90_put_att(ncid, NF90_GLOBAL, 'NumIdentityObs', Nidentity)
call nc_check(io, 'WriteNetCDF', 'put_att identity '//trim(fname))

!----------------------------------------------------------------------------
! Write all observation types that are used. Requires counting how many
! observations for each observation type.
!----------------------------------------------------------------------------

call nc_check(nf90_put_att(ncid, NF90_GLOBAL, 'comment', &
        'All used observation types follow. &
        &ObservationTypes variable has all types known.' ), &
        'WriteNetCDF', 'put_att obstypes comment '//trim(fname))

if ( verbose ) then
   ! print a banner to help identify the columns - the whitespace makes it work.
   write(string1,*)'                                        # obs      vertical'
   write(string2,*)'observation type                 possible        system    scaling'
   call error_handler(E_MSG,'WriteNetCDF',string1,text2=string2)
endif

typesdimlen = 0
do ivar = 1,max_defined_types_of_obs

   nobs = sum(prior%Nposs(:,:,:,ivar))

   if ( verbose ) then
      write(string1,'(i4,1x,(a32),1x,i8,1x,'' obs@vert '',i3,f11.3)') ivar, &
         obs_type_strings(ivar), nobs, which_vert(ivar), scale_factor(ivar)
      call error_handler(E_MSG,'WriteNetCDF',string1)
   endif

   if (nobs > 0) then
      typesdimlen = typesdimlen + 1

      if ( is_observation_trusted(obs_type_strings(ivar)) ) then
         string1 = trim(obs_type_strings(ivar))//'--TRUSTED'
      else
         string1 = obs_type_strings(ivar)
      endif

      call nc_check(nf90_put_att(ncid, NF90_GLOBAL, string1, ivar), &
         'WriteNetCDF', 'put_att:obs_type_string '//trim(obs_type_strings(ivar)))
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
          name='time', len=NF90_UNLIMITED, dimid=TimeDimID), &
          'WriteNetCDF', 'time:def_dim '//trim(fname))
call nc_check(nf90_def_dim(ncid=ncid, &
          name='bounds', len=2, dimid=BoundsDimID), &
          'WriteNetCDF', 'bounds:def_dim '//trim(fname))

call nc_check(nf90_def_dim(ncid=ncid, &
          name='copy', len=Ncopies, dimid=CopyDimID), &
          'WriteNetCDF', 'copy:def_dim '//trim(fname))

call nc_check(nf90_def_dim(ncid=ncid, &
          name='obstypes', len=max_defined_types_of_obs, dimid=TypesDimID), &
          'WriteNetCDF', 'types:def_dim '//trim(fname))

call nc_check(nf90_def_dim(ncid=ncid, &
          name='region', len=Nregions, dimid=RegionDimID), &
          'WriteNetCDF', 'region:def_dim '//trim(fname))

call nc_check(nf90_def_dim(ncid=ncid, &
          name='surface', len=1, dimid=SlevelDimID), &
          'WriteNetCDF', 'slevel:def_dim '//trim(fname))

call nc_check(nf90_def_dim(ncid=ncid, &
          name='undef', len=1, dimid=UlevelDimID), &
          'WriteNetCDF', 'ulevel:def_dim '//trim(fname))

call nc_check(nf90_def_dim(ncid=ncid, &
          name='mlevel', len=Nmlevels, dimid=MlevelDimID), &
          'WriteNetCDF', 'mlevel:def_dim '//trim(fname))
call nc_check(nf90_def_dim(ncid=ncid, &
          name='mlevel_edges', len=Nmlevels+1, dimid=MlevIntDimID), &
          'WriteNetCDF', 'mlevel_edges:def_dim '//trim(fname))

call nc_check(nf90_def_dim(ncid=ncid, &
          name='plevel', len=Nplevels, dimid=PlevelDimID), &
          'WriteNetCDF', 'plevel:def_dim '//trim(fname))
call nc_check(nf90_def_dim(ncid=ncid, &
          name='plevel_edges', len=Nplevels+1, dimid=PlevIntDimID), &
          'WriteNetCDF', 'plevel_edges:def_dim '//trim(fname))

call nc_check(nf90_def_dim(ncid=ncid, &
          name='hlevel', len=Nhlevels, dimid=HlevelDimID), &
          'WriteNetCDF', 'hlevel:def_dim '//trim(fname))
call nc_check(nf90_def_dim(ncid=ncid, &
          name='hlevel_edges', len=Nhlevels+1, dimid=HlevIntDimID), &
          'WriteNetCDF', 'hlevel_edges:def_dim '//trim(fname))

call nc_check(nf90_def_dim(ncid=ncid, &
          name='stringlength', len=stringlength, dimid=StringDimID), &
          'WriteNetCDF', 'stringlength:def_dim '//trim(fname))

if (create_rank_histogram) then
call nc_check(nf90_def_dim(ncid=ncid, &
          name='rank_bins', len=ens_size+1, dimid=RankDimID), &
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

!----------------------------------------------------------------------------
! Define the model level coordinate variable and attributes
!----------------------------------------------------------------------------

call nc_check(nf90_def_var(ncid=ncid, name='mlevel', xtype=nf90_int, &
          dimids=MlevelDimID, varid=MlevelVarID), 'WriteNetCDF', 'mlevel:def_var')
call nc_check(nf90_put_att(ncid, MlevelVarID, 'long_name', 'model level'), &
          'WriteNetCDF', 'mlevel:long_name')
call nc_check(nf90_put_att(ncid, MlevelVarID, 'units',     'model level'), &
          'WriteNetCDF', 'mlevel:units')
call nc_check(nf90_put_att(ncid, MlevelVarID, 'axis',     'Z'), &
          'WriteNetCDF', 'mlevel:axis')
call nc_check(nf90_put_att(ncid, MlevelVarID, 'valid_range',  &
           (/ minval(mlevel(1:Nmlevels)), maxval(mlevel(1:Nmlevels)) /)), &
          'WriteNetCDF', 'mlevel:valid_range')

!----------------------------------------------------------------------------
! Define the model level interface coordinate variable and attributes
!----------------------------------------------------------------------------

call nc_check(nf90_def_var(ncid=ncid, name='mlevel_edges', xtype=nf90_real, &
          dimids=MlevIntDimID, varid=MlevIntVarID), 'WriteNetCDF', 'mlevel_edges:def_var')
call nc_check(nf90_put_att(ncid, MlevIntVarID, 'long_name', 'model level edges'), &
          'WriteNetCDF', 'mlevel_edges:long_name')
call nc_check(nf90_put_att(ncid, MlevIntVarID, 'units',     'model level'), &
          'WriteNetCDF', 'mlevel_edges:units')
call nc_check(nf90_put_att(ncid, MlevIntVarID, 'axis',     'Z'), &
          'WriteNetCDF', 'mlevel_edges:axis')
call nc_check(nf90_put_att(ncid, MlevIntVarID, 'valid_range', &
   (/ minval(mlevel_edges(1:Nmlevels+1)), maxval(mlevel_edges(1:Nmlevels+1)) /)), &
          'WriteNetCDF', 'mlevel_edges:valid_range')

!----------------------------------------------------------------------------
! Define the pressure level coordinate variable and attributes
!----------------------------------------------------------------------------

call nc_check(nf90_def_var(ncid=ncid, name='plevel', xtype=nf90_real, &
          dimids=PlevelDimID, varid=PlevelVarID), 'WriteNetCDF', 'plevel:def_var')
call nc_check(nf90_put_att(ncid, PlevelVarID, 'long_name', 'pressure bin midpoints'), &
          'WriteNetCDF', 'plevel:long_name')
call nc_check(nf90_put_att(ncid, PlevelVarID, 'units',     'hPa'), &
          'WriteNetCDF', 'plevel:units')
call nc_check(nf90_put_att(ncid, PlevelVarID, 'axis',     'Z'), &
          'WriteNetCDF', 'plevel:axis')
call nc_check(nf90_put_att(ncid, PlevelVarID, 'valid_range', &
   (/ minval(plevel(1:Nplevels)), maxval(plevel(1:Nplevels)) /)), &
          'WriteNetCDF', 'plevel:valid_range')

!----------------------------------------------------------------------------
! Define the pressure level interface coordinate variable and attributes
!----------------------------------------------------------------------------

call nc_check(nf90_def_var(ncid=ncid, name='plevel_edges', xtype=nf90_real, &
          dimids=PlevIntDimID, varid=PlevIntVarID), 'WriteNetCDF', 'plevel_edges:def_var')
call nc_check(nf90_put_att(ncid, PlevIntVarID, 'long_name', 'pressure bin edges'), &
          'WriteNetCDF', 'plevel_edges:long_name')
call nc_check(nf90_put_att(ncid, PlevIntVarID, 'units',     'hPa'), &
          'WriteNetCDF', 'plevel_edges:units')
call nc_check(nf90_put_att(ncid, PlevIntVarID, 'axis',     'Z'), &
          'WriteNetCDF', 'plevel_edges:axis')
call nc_check(nf90_put_att(ncid, PlevIntVarID, 'valid_range', &
   (/ minval(plevel_edges(1:Nplevels+1)), maxval(plevel_edges(1:Nplevels+1)) /)), &
          'WriteNetCDF', 'plevel_edges:valid_range')

!----------------------------------------------------------------------------
! Define the height level coordinate variable and attributes
!----------------------------------------------------------------------------

call nc_check(nf90_def_var(ncid=ncid, name='hlevel', xtype=nf90_real, &
          dimids=HlevelDimID, varid=HlevelVarID), 'WriteNetCDF', 'hlevel:def_var')
call nc_check(nf90_put_att(ncid, HlevelVarID, 'long_name', 'height bin midpoints'), &
          'WriteNetCDF', 'hlevel:long_name')
call nc_check(nf90_put_att(ncid, HlevelVarID, 'units',     'm'), &
          'WriteNetCDF', 'hlevel:units')
call nc_check(nf90_put_att(ncid, HlevelVarID, 'axis',     'Z'), &
          'WriteNetCDF', 'hlevel:axis')
call nc_check(nf90_put_att(ncid, HlevelVarID, 'valid_range', &
          (/ minval(hlevel(1:Nhlevels)), maxval(hlevel(1:Nhlevels)) /)), &
          'WriteNetCDF', 'hlevel:valid_range')

!----------------------------------------------------------------------------
! Define the height level interface coordinate variable and attributes
!----------------------------------------------------------------------------

call nc_check(nf90_def_var(ncid=ncid, name='hlevel_edges', xtype=nf90_real, &
          dimids=HlevIntDimID, varid=HlevIntVarID), 'WriteNetCDF', 'hlevel_edges:def_var')
call nc_check(nf90_put_att(ncid, HlevIntVarID, 'long_name', 'height bin edges'), &
          'WriteNetCDF', 'hlevel_edges:long_name')
call nc_check(nf90_put_att(ncid, HlevIntVarID, 'units',     'm'), &
          'WriteNetCDF', 'hlevel_edges:units')
call nc_check(nf90_put_att(ncid, HlevIntVarID, 'axis',     'Z'), &
          'WriteNetCDF', 'hlevel_edges:axis')
call nc_check(nf90_put_att(ncid, HlevIntVarID, 'valid_range', &
          (/ minval(hlevel_edges(1:Nhlevels+1)), maxval(hlevel_edges(1:Nhlevels+1)) /)), &
          'WriteNetCDF', 'hlevel_edges:valid_range')

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
call nc_check(nf90_put_att(ncid, TimeVarID, 'units',     'days since 1601-1-1'), &
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
call nc_check(nf90_put_att(ncid, TimeBoundsVarID, 'units',     'days since 1601-1-1'), &
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
          'WriteNetCDF', 'set_nofill '//trim(fname))

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

call nc_check(nf90_put_var(ncid, MlevelVarID, mlevel(1:Nmlevels) ), &
          'WriteNetCDF', 'mlevel:put_var')
call nc_check(nf90_put_var(ncid, MlevIntVarID, mlevel_edges(1:Nmlevels+1)), &
          'WriteNetCDF', 'mlevel_edges:put_var')

call nc_check(nf90_put_var(ncid, PlevelVarID, plevel(1:Nplevels)), &
          'WriteNetCDF', 'plevel:put_var')
call nc_check(nf90_put_var(ncid, PlevIntVarID, plevel_edges(1:Nplevels+1)), &
          'WriteNetCDF', 'plevel_edges:put_var')

call nc_check(nf90_put_var(ncid, HlevelVarID, hlevel(1:Nhlevels)), &
          'WriteNetCDF', 'hlevel:put_var')
call nc_check(nf90_put_var(ncid, HlevIntVarID, hlevel_edges(1:Nhlevels+1)), &
          'WriteNetCDF', 'hlevel_edges:put_var')

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
! The priors always have values. It is possible that there are no posteriors.
!----------------------------------------------------------------------------

write(*,*)'Finished writing metadata ... defining output variables ...'
write(*,*)'(This may take a couple minutes with a thousand timesteps and hundreds of gages.)'

if ( create_rank_histogram ) then
   call Define_Prior_and_Posterior(ncid, prior, poste, &
                    TimeDimID, CopyDimID, RegionDimID, RankDimID)
   write(*,*)'Finished defining output variables ... writing priors ...'
   call WriteTLRV(ncid, prior, RankDimID)
else
   call Define_Prior_and_Posterior(ncid, prior, poste, &
                    TimeDimID, CopyDimID, RegionDimID)
   write(*,*)'Finished defining output variables ... writing priors ...'
   call WriteTLRV(ncid, prior)
endif

write(*,*)'Finished writing priors ... writing posteriors.'
if (has_posteriors) call WriteTLRV(ncid, poste)

!----------------------------------------------------------------------------
! finish ...
!----------------------------------------------------------------------------

write(*,*)'Finishing up.'
call nc_check(nf90_sync( ncid), 'WriteNetCDF', 'sync '//trim(fname))
call nc_check(nf90_close(ncid), 'init_diag_output', 'close '//trim(fname))

end subroutine WriteNetCDF


!======================================================================


function Rlevels2edges(level_middle, level_edge)
! determine layer edges ... ascending vs. descending complications
! pressures should be descending.
! heights   should be ascending, except when they are depths ...
integer                               :: Rlevels2edges
real(r8), dimension(:), intent(inout) :: level_middle
real(r8), dimension(:), intent(  out) :: level_edge

logical :: increasing
integer :: n, i

Rlevels2edges = SIZE(level_middle,1)
if (Rlevels2edges > MaxLevels) then
   write(string1,*)Rlevels2edges,' is too many levels - max is ',MaxLevels,' (MaxLevels)'
   call error_handler(E_ERR,'Rlevels2edges', string1,source)
endif

! find length of useful portion of level_middle array ... among other things.
n = RRemoveDuplicates(level_middle)
Rlevels2edges = n

! single-level case works for pressures, must change sign for heights
! but there is no way to know if heights vs. depths right now ...
if ( n == 1 ) then
   level_edge(1) = level_middle(1) + 0.1 * level_middle(1)
   level_edge(2) = level_middle(1) - 0.1 * level_middle(1)
   return
endif

! ensure monotonicity - up or down -
if (level_middle(1) < level_middle(n)) then
   increasing        = .true.
   level_middle(1:n) = sort(level_middle(1:n))
else
   increasing        = .false.
   level_edge(1:n)   = sort(level_middle(1:n)) ! use level_edge as a temporary array
   level_middle(1:n) = level_edge(n:1:-1)
endif

level_edge(  1) = level_middle(1) - (level_middle(2) - level_middle(  1))/2.0_r8
level_edge(n+1) = level_middle(n) + (level_middle(n) - level_middle(n-1))/2.0_r8
do i = 2,n
  level_edge(i) = level_middle(i-1) + (level_middle(i) - level_middle(i-1))/2.0_r8
enddo

if ( increasing ) then ! must be heights or model levels
!  level_edge(  1) = max(   0.0_r8, level_edge(  1))
else                   ! must be pressures
   level_edge(  1) = min(1025.0_r8, level_edge(  1))
   level_edge(n+1) = max(   0.0_r8, level_edge(n+1))
endif

end function Rlevels2edges


!======================================================================


function Redges2levels(level_edge, level_middle)
! determine layer edge ... ascending vs. descending complications
integer :: Redges2levels
real(r8), dimension(:), intent(inout) :: level_edge
real(r8), dimension(:), intent(  out) :: level_middle

logical :: increasing
integer :: n, i

Redges2levels = SIZE(level_edge,1) - 1
if (Redges2levels > MaxLevels) then
   write(string1,*)Redges2levels,' is too many level_middle - max is ',MaxLevels,' (MaxLevels)'
   call error_handler(E_ERR,'Redges2levels', string1,source)
endif

! find the number of edges the user specified
n = RRemoveDuplicates(level_edge)
if ( n == 1 ) then
   write(string1,*)'obs_diag_nml: need at least two entries in hlevel_edges. Only got 1.'
   call error_handler(E_ERR,'Redges2levels', string1,source)
   return
endif

Redges2levels = n - 1  ! returns the number of layers, not the number of edges

! resolve [ascending/descending] monotonicity
if (level_edge(1) < level_edge(n) ) then
   increasing      = .true.
   level_edge(1:n) = sort(level_edge(1:n))
else
   increasing      = .false.
   level_edge(1:n) = sort(level_edge(1:n))
   level_edge(1:n) = level_edge(n:1:-1)
endif

do i = 1,n-1
  level_middle(i) = (level_edge(i) + level_edge(i+1))/2.0_r8
enddo

end function Redges2levels


!======================================================================


function find_my_level(level_in, leveltype)
! The levels, intervals are ordered :
! closest to the center of the earth == 1
! outer space                        == N
!
! nlev, levels and levels_int are scoped in this unit.
!
! The complexity arises now that the strategy that uses minloc() does not
! work for level midpoints that were the result of being specified by the
! level edges.
!
! Scenario 1: layer midpoints are defined. layer edges are derived.
! This case is the classic one. The layer edges are COMPUTED to be halfway
! between the layer midpoints.
!
! defined:    |  100  |  200  |  500  |  1000  |
! derived:   50      150     350     750      1250
!
! You can always bin by detemining the closest midpoint. By construction.
!
! Scenario 2: layer edges are defined. layer midpoints are derived.
! Now the closest midpoint is NOT definitive.
!
! defined:   50      150     350     750      1250
! derived:    |  100  |  250  |  550  |  1000  |
!
! HOWEVER: A level of 760 is only 210 away from the 550,
!          and 240 away from the proper bin center of 1000.

real(r8), intent(in) :: level_in         ! target level
integer,  intent(in) :: leveltype        ! indicator for surface, height, etc.
integer              :: find_my_level

integer  :: i
real(r8) :: surface, top

find_my_level = -1 ! set this to a 'bad' value.

if (leveltype == VERTISHEIGHT .or. &
    leveltype == VERTISSCALEHEIGHT) then ! we have heights levels

   if      ( level_in < hlevel_edges(1) ) then
      find_my_level = -1
   elseif ( level_in > hlevel_edges(Nhlevels+1) ) then
      find_my_level = 100 + Nhlevels
   elseif ( level_in == hlevel_edges(1) ) then
      find_my_level = 1
   else
      ! Starting close to the center of the earth ...
      HEIGHTLOOP : do i = 1,Nhlevels
         if (level_in >  hlevel_edges(i) .and. &
             level_in <= hlevel_edges(i+1)) then
            find_my_level = i
            exit HEIGHTLOOP
         endif
      enddo HEIGHTLOOP
    endif

elseif (leveltype == VERTISLEVEL) then   ! we have model levels

   surface = minval(mlevel_edges(1:Nmlevels+1))
   top     = maxval(mlevel_edges(1:Nmlevels+1))

   if (      level_in < surface ) then
      find_my_level = -1
   elseif ( level_in > top     ) then
      find_my_level = 100 + Nmlevels
   elseif ( level_in == mlevel_edges(1) ) then
      find_my_level = 1
   else
      ! Starting close to the center of the earth ...
      LEVELLOOP : do i = 1,Nmlevels
         if (level_in >  mlevel_edges(i) .and. &
             level_in <= mlevel_edges(i+1)) then
            find_my_level = i
            exit LEVELLOOP
         endif
      enddo LEVELLOOP
   endif

elseif (leveltype == VERTISUNDEF) then

      find_my_level = VERTISUNDEF    ! VERTISUNDEF == -2, a good bad value

else  ! we have pressure levels (or surface obs ... also in pressure units)

   surface = plevel_edges(1)
   top     = plevel_edges(Nplevels+1)

   if (      level_in > surface ) then ! below the lowest layer
      find_my_level = -1
   elseif ( level_in < top     ) then ! above the highest layer
      find_my_level = 100 + Nplevels
   elseif ( level_in == plevel_edges(1) ) then
      find_my_level = 1
   else
      ! Starting close to the center of the earth ...
      ! pressure greatest at surface.
      PRESSURELOOP : do i = 1,Nplevels
         if (level_in <  plevel_edges(i) .and. &
             level_in >= plevel_edges(i+1)) then
            find_my_level = i
            exit PRESSURELOOP
         endif
      enddo PRESSURELOOP
   endif

endif

end function find_my_level


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

subroutine Define_Prior_and_Posterior(ncid, prior, poste, TimeDimID, CopyDimID, RegionDimID, RankDimID)

integer,           intent(in) :: ncid
type(TLRV_type),   intent(in) :: prior
type(TLRV_type),   intent(in) :: poste
integer,           intent(in) :: TimeDimID, CopyDimID, RegionDimID
integer, optional, intent(in) :: RankDimID

character(len=*), parameter :: routine = 'Define_Prior_and_Posterior'

integer :: nobs, Nlevels, ivar
integer :: ndata
character(len=NF90_MAX_NAME) :: string1, string2

integer :: VarID, LevelDimID, oldmode

call nc_check(nf90_redef(ncid), routine, 'redef')

call nc_check(nf90_set_fill(ncid, NF90_NOFILL, oldmode),  &
           routine, 'set_nofill')

DEFINE : do ivar = 1,num_obs_types

   nobs = sum(prior%Nposs(:,:,:,ivar))
   if (nobs < 1) cycle DEFINE

   ! Create netCDF variable name

   string3 = obs_type_strings(ivar)
   string1 = trim(string3)//'_'//adjustl(prior%string)
   string2 = trim(string3)//'_'//adjustl(poste%string)

   ! Determine what kind of levels to use (LevelDimID) ... models, pressure, height ...

   Nlevels = FindVertical(ncid, ivar, LevelDimID)

   ! Define the prior variable and its attributes.

   call nc_check(nf90_def_var(ncid, name=string1, xtype=nf90_real, &
          dimids=(/ RegionDimID, LevelDimID, CopyDimID, TimeDimID /), &
          varid=VarID), routine, 'region:def_var '//trim(string1))
   call nc_check(nf90_put_att(ncid, VarID, '_FillValue',    MISSING_R4), &
           routine,'put_att:fillvalue '//trim(string1))
   call nc_check(nf90_put_att(ncid, VarID, 'missing_value', MISSING_R4), &
           routine,'put_att:missing '//trim(string1))

   if ( is_observation_trusted(obs_type_strings(ivar)) ) then
      call nc_check(nf90_put_att(ncid, VarID, 'TRUSTED', 'TRUE'), &
           routine,'put_att:trusted '//trim(string1))
      call error_handler(E_MSG,routine,string1,text2='is TRUSTED.')
   endif

   ! Define the poste variable and its attributes.

   if (has_posteriors) then

      call nc_check(nf90_def_var(ncid, name=string2, xtype=nf90_real, &
             dimids=(/ RegionDimID, LevelDimID, CopyDimID, TimeDimID /), &
             varid=VarID), routine, 'region:def_var '//trim(string2))
      call nc_check(nf90_put_att(ncid, VarID, '_FillValue',    MISSING_R4), &
              routine,'put_att:fillvalue '//trim(string2))
      call nc_check(nf90_put_att(ncid, VarID, 'missing_value', MISSING_R4), &
              routine,'put_att:missing '//trim(string2))
   
      if ( is_observation_trusted(obs_type_strings(ivar)) ) then
         call nc_check(nf90_put_att(ncid, VarID, 'TRUSTED', 'TRUE'), &
              routine,'put_att:trusted '//trim(string2))
         call error_handler(E_MSG,routine,string2,text2='is TRUSTED.')
      endif

   endif

   ! The rank histogram has no 'copy' dimension, so it must be handled differently.

   if (present(RankDimID)) then

      string2 = trim(string1)//'_RankHist'
      ndata   = sum(prior%hist_bin(:,:,:,ivar,:))

      if ( ndata > 0 ) then
         call nc_check(nf90_def_var(ncid, name=string2, xtype=nf90_int, &
             dimids=(/ RegionDimID, LevelDimID, RankDimID, TimeDimID /), &
             varid=VarID), routine, 'rank_hist:def_var '//trim(string2))
      endif

      if (verbose) then
         write(string3,*)ndata,' observations '//trim(string2)
         call error_handler(E_MSG,routine,string3)
      endif

   endif

enddo DEFINE

call nc_check(nf90_enddef(ncid), routine, 'enddef ')

end subroutine Define_Prior_and_Posterior


!======================================================================


subroutine WriteTLRV(ncid, vrbl, RankDimID)
integer,           intent(in) :: ncid
type(TLRV_type),   intent(in) :: vrbl
integer, optional, intent(in) :: RankDimID

integer :: nobs, Nlevels, ivar, itime, ilevel, iregion
integer :: Nbins, irank, ndata
character(len=NF90_MAX_NAME) :: string1, string2

integer :: VarID, LevelDimID
real(r4), allocatable, dimension(:,:,:,:) :: rchunk
integer,  allocatable, dimension(:,:,:,:) :: ichunk

FILL : do ivar = 1,num_obs_types

   nobs = sum(vrbl%Nposs(:,:,:,ivar))
   if (nobs < 1) cycle FILL

   ! Create netCDF variable name

   string2 = obs_type_strings(ivar)
   string1 = trim(string2)//'_'//adjustl(vrbl%string)

   ! determine what kind of levels to use ... models, pressure, height ...

   Nlevels = FindVertical(ncid, ivar, LevelDimID)

   allocate(rchunk(Nregions,Nlevels,Ncopies,Nepochs))
   rchunk = MISSING_R4

   do itime   = 1,Nepochs
   do ilevel  = 1,Nlevels
   do iregion = 1,Nregions

      rchunk(iregion,ilevel, 1,itime) = vrbl%Nposs(      itime,ilevel,iregion,ivar)
      rchunk(iregion,ilevel, 2,itime) = vrbl%Nused(      itime,ilevel,iregion,ivar)
      rchunk(iregion,ilevel, 3,itime) = vrbl%NbigQC(     itime,ilevel,iregion,ivar)
      rchunk(iregion,ilevel, 4,itime) = vrbl%NbadIZ(     itime,ilevel,iregion,ivar)
      rchunk(iregion,ilevel, 5,itime) = vrbl%NbadUV(     itime,ilevel,iregion,ivar)
      rchunk(iregion,ilevel, 6,itime) = vrbl%NbadLV(     itime,ilevel,iregion,ivar)
      rchunk(iregion,ilevel, 7,itime) = vrbl%rmse(       itime,ilevel,iregion,ivar)
      rchunk(iregion,ilevel, 8,itime) = vrbl%bias(       itime,ilevel,iregion,ivar)
      rchunk(iregion,ilevel, 9,itime) = vrbl%spread(     itime,ilevel,iregion,ivar)
      rchunk(iregion,ilevel,10,itime) = vrbl%totspread(  itime,ilevel,iregion,ivar)
      rchunk(iregion,ilevel,11,itime) = vrbl%NbadDartQC( itime,ilevel,iregion,ivar)
      rchunk(iregion,ilevel,12,itime) = vrbl%observation(itime,ilevel,iregion,ivar)
      rchunk(iregion,ilevel,13,itime) = vrbl%ens_mean(   itime,ilevel,iregion,ivar)
      rchunk(iregion,ilevel,14,itime) = vrbl%Ntrusted(   itime,ilevel,iregion,ivar)
      rchunk(iregion,ilevel,15,itime) = vrbl%NDartQC_0(  itime,ilevel,iregion,ivar)
      rchunk(iregion,ilevel,16,itime) = vrbl%NDartQC_1(  itime,ilevel,iregion,ivar)
      rchunk(iregion,ilevel,17,itime) = vrbl%NDartQC_2(  itime,ilevel,iregion,ivar)
      rchunk(iregion,ilevel,18,itime) = vrbl%NDartQC_3(  itime,ilevel,iregion,ivar)
      rchunk(iregion,ilevel,19,itime) = vrbl%NDartQC_4(  itime,ilevel,iregion,ivar)
      rchunk(iregion,ilevel,20,itime) = vrbl%NDartQC_5(  itime,ilevel,iregion,ivar)
      rchunk(iregion,ilevel,21,itime) = vrbl%NDartQC_6(  itime,ilevel,iregion,ivar)
      rchunk(iregion,ilevel,22,itime) = vrbl%NDartQC_7(  itime,ilevel,iregion,ivar)
      rchunk(iregion,ilevel,23,itime) = vrbl%NDartQC_8(  itime,ilevel,iregion,ivar)

   enddo
   enddo
   enddo

   call nc_check(nf90_inq_varid(ncid, string1, VarID), &
           'WriteTLRV', 'region:inq_varid '//trim(string1))
   call nc_check(nf90_put_var(ncid, VarID, rchunk ), &
           'WriteTLRV', 'realchunk:put_var '//trim(string1))
   deallocate(rchunk)

   ! The rank histogram has no 'copy' dimension, so it must be handled differently.

   if ( present(RankDimID) ) then

      string2 = trim(string1)//'_RankHist'
      Nbins   = size(vrbl%hist_bin,5)
      ndata   = sum(vrbl%hist_bin(:,:,:,ivar,:))

      if ( ndata > 0 ) then

         allocate(ichunk(Nregions,Nlevels,Nbins,Nepochs))
         ichunk = 0

         do itime   = 1,Nepochs
         do ilevel  = 1,Nlevels
         do iregion = 1,Nregions
         do irank   = 1,Nbins

         ichunk(iregion,ilevel,irank,itime) = vrbl%hist_bin(itime,ilevel,iregion,ivar,irank)

         enddo
         enddo
         enddo
         enddo

         call nc_check(nf90_inq_varid(ncid, string2, VarID), &
                 'WriteTLRV', 'rank_hist:inq_varid '//trim(string2))
         call nc_check(nf90_put_var(ncid, VarID, ichunk ), &
                 'WriteTLRV', 'intchunk:put_var '//trim(string2))

         deallocate(ichunk)

      endif
   endif

enddo FILL

end subroutine WriteTLRV


!======================================================================


function FindVertical(ncid, flav, dimid)
integer, intent(in)  :: ncid, flav
integer, intent(out) :: dimid
integer              :: FindVertical

if      ( which_vert(flav) == VERTISSURFACE ) then

   FindVertical = 1
   call nc_check(nf90_inq_dimid(ncid, 'surface', dimid), &
                                 'FindVertical', 'vertissurface')

elseif ( which_vert(flav) == VERTISUNDEF   ) then
   FindVertical = 1
   call nc_check(nf90_inq_dimid(ncid, 'undef', dimid), &
                                 'FindVertical', 'vertisundef')

elseif ( which_vert(flav) == VERTISLEVEL   ) then

   FindVertical = Nmlevels
   call nc_check(nf90_inq_dimid(ncid, 'mlevel', dimid), &
                                 'FindVertical', 'vertislevel')

elseif ( which_vert(flav) == VERTISPRESSURE) then

   FindVertical = Nplevels
   call nc_check(nf90_inq_dimid(ncid, 'plevel', dimid), &
                                 'FindVertical', 'vertispressure')

elseif ( which_vert(flav) == VERTISHEIGHT .or. &
          which_vert(flav) == VERTISSCALEHEIGHT  ) then

   FindVertical = Nhlevels
   call nc_check(nf90_inq_dimid(ncid, 'hlevel', dimid), &
                                 'FindVertical', 'vertisheight')

else
   call error_handler(E_ERR,'FindVertical','unknown vertical', source)
endif

end function FindVertical


!======================================================================


function RRemoveDuplicates(mydata)

! remove duplicate levels from user input.
! This does not sort first, probably could be faster if it did.
! however, the levels arrays are length about a dozen or so.
integer                               :: RRemoveDuplicates
real(r8), dimension(:), intent(inout) :: mydata

real(r8), dimension(SIZE(mydata)) :: temparray
integer :: n

if (mydata(1) == MISSING_R8) then
   RRemoveDuplicates = 0
   return
endif

n = 1
temparray(1) = mydata(1)

do i = 2,SIZE(mydata)
   if (mydata(i) == MISSING_R8) then
      ! Do nothing.
   elseif (any(mydata(i) == temparray(1:n))) then
      write(string1,*)'   WARNING ... duplicate namelist level found ',mydata(i)
      write(string2,*)'WARNING ... duplicate namelist level found ',mydata(i)
      write(string3,*)'removing duplicate value.'
      call error_handler(E_MSG,'RRemoveDuplicates', &
                 string1,text2=string2,text3=string3)
   else
     n = n + 1
     temparray(n) = mydata(i)
   endif
enddo

mydata(1:n) = temparray(1:n)
RRemoveDuplicates = n

end function RRemoveDuplicates


!======================================================================


function IRemoveDuplicates(mydata)

! remove duplicate levels from user input.
! This does not sort first, probably could be faster if it did.
! however, the levels arrays are length about a dozen or so.
integer                              :: IRemoveDuplicates
integer, dimension(:), intent(inout) :: mydata

integer, dimension(SIZE(mydata)) :: temparray
integer :: n

if (mydata(1) == MISSING_I) then
   IRemoveDuplicates = 0
   return
endif

n = 1
temparray(1) = mydata(1)

do i = 2,SIZE(mydata)
   if (mydata(i) == MISSING_I) then
      ! Do nothing.
   elseif (any(mydata(i) == temparray(1:n))) then
      write(string1,*)'warning ... duplicate namelist level found ',mydata(i)
      write(string2,*)'removing duplicate value.'
      call error_handler(E_MSG,'IRemoveDuplicates',string1,text2=string2)
   else
     n = n + 1
     temparray(n) = mydata(i)
   endif
enddo

mydata(1:n) = temparray(1:n)
IRemoveDuplicates = n

end function IRemoveDuplicates

!======================================================================
!> convert the state vector index into the index of the array normally used
!> to hold the integer codes of the observation kinds

function lookup_flavor(indx) result(nflavors)

integer(i8), intent(inout) :: indx
integer                    :: nflavors

logical, save :: initialized = .false.
integer(i8), allocatable, save :: identity_locations(:)
integer, save :: maxflavors 

integer :: iflavor

if (.not. initialized) then
   allocate(identity_locations(MaxGages))
   maxflavors    = 0
   initialized = .true.
endif

indx = abs(indx)

REMAP : do iflavor = 1,maxflavors
   if (indx == identity_locations(iflavor)) then
      indx = iflavor
      nflavors = maxflavors
      return
   endif
enddo REMAP

! If we make it to here, we have a new identity observation 
! to add to the list.

if (maxflavors < MaxGages) then
   maxflavors = maxflavors + 1
   identity_locations(maxflavors) = indx

   write(obs_type_strings(maxflavors),'("STREAM_FLOW_",i6.6)') indx

   if (verbose) then
      write(*,*)'making ',trim(obs_type_strings(maxflavors))
      write(string1,*)'Identity observation at ',indx, &
                      ' being added to list as ',maxflavors
      call error_handler(E_MSG, 'lookup_flavor', string1,  source )
   endif

   indx = maxflavors
   nflavors = maxflavors

else
   write(string1,*)'Identity observation at ',indx,' exceeds capacity.'
   write(string2,*)'Increase "MaxGages", recompile, and rerun.'
   write(string3,*)'"MaxGages" is currently ',MaxGages
   call error_handler(E_ERR, 'lookup_flavor', string1, & 
              source, text2=string2, text3=string3)
endif

end function lookup_flavor

!======================================================================


end program obs_diag

