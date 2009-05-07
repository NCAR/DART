! Data Assimilation Research Testbed -- DART
! Copyright 2004-2007, Data Assimilation Research Section
! University Corporation for Atmospheric Research
! Licensed under the GPL -- www.gpl.org/licenses/gpl.html

program obs_diag

! <next few lines under version control, do not edit>
! $URL$
! $Id$
! $Revision$
! $Date$

!-----------------------------------------------------------------------
! The programs defines a series of epochs (periods of time) and geographic
! regions and accumulates statistics for these epochs and regions.
!-----------------------------------------------------------------------

use        types_mod, only : r8, digits12
use obs_sequence_mod, only : read_obs_seq, obs_type, obs_sequence_type, get_first_obs, &
                             get_obs_from_key, get_obs_def, get_copy_meta_data, &
                             get_obs_time_range, get_time_range_keys, get_num_obs, &
                             get_next_obs, get_num_times, get_obs_values, init_obs, &
                             assignment(=), get_num_copies, static_init_obs_sequence, &
                             get_qc, destroy_obs_sequence, get_last_obs, get_num_qc, &
                             read_obs_seq_header, destroy_obs, get_qc_meta_data
use      obs_def_mod, only : obs_def_type, get_obs_def_error_variance, get_obs_def_time, &
                             get_obs_def_location, get_obs_kind, get_obs_name
use     obs_kind_mod, only : max_obs_kinds
use     location_mod, only : location_type, get_location, set_location_missing, &
                             write_location, operator(/=)
use time_manager_mod, only : time_type, set_date, set_time, get_time, print_time, &
                             print_date, &
                             operator(*), operator(+), operator(-), &
                             operator(>), operator(<), operator(/), &
                             operator(/=)
use    utilities_mod, only : get_unit, open_file, close_file, register_module, &
                             file_exist, error_handler, E_ERR, E_MSG,          &
                             initialize_utilities, logfileunit, nmlfileunit,   &
                             find_namelist_in_file, check_namelist_read,       &
                             do_output, do_nml_file, do_nml_term, timestamp

implicit none

! version controlled file description for error handling, do not edit
character(len=128), parameter :: &
   source   = "$URL$", &
   revision = "$Revision$", &
   revdate  = "$Date$"

type(obs_sequence_type) :: seq
type(obs_type)          :: observation, next_obs
type(obs_type)          :: obs1, obsN
type(obs_def_type)      :: obs_def
type(location_type)     :: obs_loc


!---------------------
integer :: obsindex, i, j, iunit, io, ivarcount
character(len=129) :: obs_seq_in_file_name

! Storage with fixed size for observation space diagnostics
real(r8), dimension(1) :: prior_mean, posterior_mean, prior_spread, posterior_spread
real(r8) :: pr_mean, po_mean ! same as above, without useless dimension 
real(r8) :: pr_sprd, po_sprd ! same as above, without useless dimension

integer :: qc_index
integer :: obs_copy_index, prior_mean_index, posterior_mean_index
integer :: prior_spread_index, posterior_spread_index
integer :: key_bounds(2), flavor
integer :: num_copies, num_qc, num_obs, max_num_obs, obs_seq_file_id
character(len=129) :: obs_seq_read_format
logical :: pre_I_format

real(r8), dimension(1) :: obs
real(r8), dimension(2) :: qc
real(r8) :: obs_err_var

integer,  allocatable :: keys(:)

logical :: out_of_range, is_there_one, is_this_last, keeper

integer, parameter :: MaxRegions = 4 

!-----------------------------------------------------------------------
! Namelist with default values
!
character(len = 129) :: obs_sequence_name = "obs_seq.final"
integer :: iskip_days = 0        ! skip the first 'iskip' days
integer :: obs_select = 1        ! obs type selection: 1=all, 2 =RAonly, 3=noRA
real(r8):: rat_cri    = 4.0      ! QC ratio
real(r8):: input_qc_threshold = 4.0    ! input obs QC values >= this not used
integer :: bin_width_seconds  = 0   ! width of the bin seconds
logical :: verbose = .false.

! index 1 == region 1 == [0.0, 1.0) i.e. Entire domain
! index 2 == region 2 == [0.0, 0.5)
! index 3 == region 3 == [0.5, 1.0)

real(r8), dimension(MaxRegions) :: lonlim1 = (/ 0.0_r8, 0.0_r8, 0.5_r8, -1.0_r8 /)
real(r8), dimension(MaxRegions) :: lonlim2 = (/ 1.0_r8, 0.5_r8, 1.0_r8, -1.0_r8 /)

character(len=6), dimension(MaxRegions) :: reg_names = &
                                   (/ 'whole ','yin   ','yang  ','bogus '/)

namelist /obs_diag_nml/ obs_sequence_name, &
                       iskip_days, obs_select, rat_cri, &
                       input_qc_threshold, bin_width_seconds, &
                       lonlim1, lonlim2, reg_names, verbose

!-----------------------------------------------------------------------
! Variables used to accumulate the statistics.
! Dimension 1 is temporal, actually - these are time-by-region-by-var
!-----------------------------------------------------------------------

integer,  allocatable, dimension(:,:,:) :: ges_level_N, anl_level_N

! statistics by time, for a particular level:  time-region-variable
real(r8), allocatable, dimension(:,:,:) :: rms_ges_mean   ! prior mean - obs
real(r8), allocatable, dimension(:,:,:) :: rms_ges_spread ! prior spread     
real(r8), allocatable, dimension(:,:,:) :: rms_anl_mean   ! posterior mean - obs
real(r8), allocatable, dimension(:,:,:) :: rms_anl_spread ! posterior spread

type(time_type), allocatable, dimension(:) :: bincenter
real(digits12),  allocatable, dimension(:) :: epoch_center
integer,         allocatable, dimension(:) :: obs_used_in_epoch

!-----------------------------------------------------------------------
! General purpose variables
!-----------------------------------------------------------------------

integer  :: seconds, days
integer  :: Nepochs
integer  :: gesUnit, anlUnit

! These pairs of variables are used when we diagnose which observations 
! are far from the background.
integer, parameter :: MaxSigmaBins = 100  
integer  :: nsigma(0:MaxSigmaBins) = 0
integer  :: indx

real(r8) :: ratio

type(time_type) :: beg_time, end_time  ! of the particular bin
type(time_type) :: seqT1, seqTN        ! first,last time in entire observation sequence
type(time_type) :: binsep, binwidth, halfbinwidth 
type(time_type) :: obs_time, skip_time

character(len = 129) :: gesName, anlName, msgstring

integer  :: num_obs_in_epoch 
integer  :: Nregions, iregion, iepoch, ivar
real(r8) :: rlocation

!-----------------------------------------------------------------------
! If the observation is a 'prior' or a 'posterior',
! we can perform different levels of quality control checking.
! 'DART quality control' value and meaning:
!
! 7 == outlier rejected
! 6 == prior qc rejected
! 5 == not used
! 4 == prior forward operator failed.
!      ----- the prior observation cannot be used -----
! 3 == Evaluated only BUT posterior forward operator failed.
! 2 == O.K. BUT posterior forward operator failed.
!      ----- the posterior observation cannot be used -----
! 1 == Evaluated only
! 0 == all O.K.
!      ----- all can be used -----
!-----------------------------------------------------------------------

integer             :: dart_qc_index, qc_integer
integer, parameter  :: QC_MAX = 7
integer, parameter  :: QC_MAX_PRIOR     = 3
integer, parameter  :: QC_MAX_POSTERIOR = 1
integer, dimension(0:QC_MAX) :: qc_counter = 0

!-----------------------------------------------------------------------
! Some variables to keep track of who's rejected why ...
!-----------------------------------------------------------------------

integer :: Nidentity  = 0   ! identity observations are not appropriate.
integer :: NwrongType = 0   ! namelist discrimination
integer :: NbadLevel  = 0   ! out-of-range pressures
integer :: NbadQC     = 0   ! ratio

! track rejection by time-region-variable  [prior/posterior]
integer, allocatable, dimension(:,:,:)   :: Nrejected
integer, allocatable, dimension(:,:,:,:) :: NbadQClevel

!=======================================================================
! Get the party started
!=======================================================================

call initialize_utilities('obs_diag')
call register_module(source,revision,revdate) 
call static_init_obs_sequence()

prior_mean(1)       = 0.0_r8
prior_spread(1)     = 0.0_r8
posterior_mean(1)   = 0.0_r8
posterior_spread(1) = 0.0_r8

!----------------------------------------------------------------------
! Read the namelist
!----------------------------------------------------------------------

call find_namelist_in_file("input.nml", "obs_diag_nml", iunit)
read(iunit, nml = obs_diag_nml, iostat = io)
call check_namelist_read(iunit, io, "obs_diag_nml")

! Record the namelist values used for the run ...
if (do_nml_file()) write(nmlfileunit, nml=obs_diag_nml)
if (do_nml_term()) write(     *     , nml=obs_diag_nml)

!----------------------------------------------------------------------
! Now that we have input, do some checking and setup
!----------------------------------------------------------------------

! Determine the number of regions from namelist input
Nregions = 0
FindNumRegions : do iregion = 1,MaxRegions

   if (( lonlim1(iregion)  <   0      ) .or. &
       ( lonlim2(iregion)  <   0      ) .or. &
       (reg_names(iregion) == 'bogus ') ) exit FindNumRegions

   Nregions = Nregions + 1

enddo FindNumRegions


!----------------------------------------------------------------------
! Determine temporal bin characteristics.
! Nepochs is the total number of time intervals of the period requested.
!----------------------------------------------------------------------
! The low order models are fundamentally different than the high-order 
! models. The low-order models assimilate at every time step and have
! a much different connotation of 'time' ... i.e.  no real calendar.
!----------------------------------------------------------------------

!ObsFileLoop : do ifile=1, tot_days*4
!-----------------------------------------------------------------------

   ! Read in information about observation sequence so we can allocate
   ! observations. We need info about how many copies, qc values, etc.

   obs_seq_in_file_name = trim(adjustl(obs_sequence_name)) ! Lahey requirement

   call read_obs_seq_header(obs_seq_in_file_name, &
             num_copies, num_qc, num_obs, max_num_obs, &
             obs_seq_file_id, obs_seq_read_format, pre_I_format, &
             close_the_file = .true.)

   if ( verbose ) then
      write(logfileunit,*)'num_copies          is ',num_copies
      write(logfileunit,*)'num_qc              is ',num_qc
      write(logfileunit,*)'num_obs             is ',num_obs
      write(logfileunit,*)'max_num_obs         is ',max_num_obs
      write(logfileunit,*)'obs_seq_read_format is ',trim(adjustl(obs_seq_read_format))
      write(logfileunit,*)'pre_I_format        is ',pre_I_format
      write(    *      ,*)'num_copies          is ',num_copies
      write(    *      ,*)'num_qc              is ',num_qc
      write(    *      ,*)'num_obs             is ',num_obs
      write(    *      ,*)'max_num_obs         is ',max_num_obs
      write(    *      ,*)'obs_seq_read_format is ',trim(adjustl(obs_seq_read_format))
      write(    *      ,*)'pre_I_format        is ',pre_I_format
   endif

   ! Initialize some (individual) observation variables

   call init_obs(       obs1, num_copies, num_qc)   ! First obs in sequence
   call init_obs(       obsN, num_copies, num_qc)   ! Last  obs in sequence
   call init_obs(observation, num_copies, num_qc)   ! current obs
   call init_obs(   next_obs, num_copies, num_qc)   ! duh ...

   ! Read in the entire observation sequence

   call read_obs_seq(obs_sequence_name, 0, 0, 0, seq)

   !--------------------------------------------------------------------
   ! The observations for the low-order models are all exactly at 
   ! the assimilation timestep. So we know the bin separation. 
   !--------------------------------------------------------------------

   is_there_one = get_first_obs(seq, obs1)
   if ( .not. is_there_one ) then
         call error_handler(E_ERR,'obs_diag','No first observation in sequence.', &
         source,revision,revdate)
   endif
   call get_obs_def(obs1,     obs_def)
   seqT1   = get_obs_def_time(obs_def)

   is_there_one = get_last_obs(seq, obsN)
   if ( .not. is_there_one ) then
         call error_handler(E_ERR,'obs_diag','No last observation in sequence.', &
         source,revision,revdate)
   endif
   call get_obs_def(obsN,     obs_def)
   seqTN   = get_obs_def_time(obs_def)

   call print_time(seqT1,'First observation time',logfileunit)
   call print_time(seqT1,'First observation time')
   call print_time(seqTN,'Last  observation time',logfileunit)
   call print_time(seqTN,'Last  observation time')

   !--------------------------------------------------------------------
   ! If the last observation is before the period of interest, move on.
   !--------------------------------------------------------------------

   ! we always process the entire file

   !--------------------------------------------------------------------
   ! If the first observation is after the period of interest, finish.
   !--------------------------------------------------------------------

   ! we always process the entire file

   !--------------------------------------------------------------------
   ! Get the number of different times in the sequence.
   ! The low order models assimilate at every time step.
   !--------------------------------------------------------------------

   beg_time  = seqT1
   end_time  = set_time(0, iskip_days)  ! convert to a time_type object
   skip_time = beg_time + end_time  ! the start time for the accumulation of vertical statistics 

   Nepochs = get_num_times(seq)
   write(*,*)'Sequence has ',Nepochs,' different times in it.'

   allocate(bincenter(Nepochs), epoch_center(Nepochs))  ! time_type and 'real'
   allocate(obs_used_in_epoch(Nepochs))
   obs_used_in_epoch = 0
   epoch_center      = -999.0_digits12

   allocate(rms_ges_mean(  Nepochs, Nregions, max_obs_kinds), &
            rms_ges_spread(Nepochs, Nregions, max_obs_kinds), &
            rms_anl_mean(  Nepochs, Nregions, max_obs_kinds), &
            rms_anl_spread(Nepochs, Nregions, max_obs_kinds))
   allocate( ges_level_N(  Nepochs, Nregions, max_obs_kinds), &
             anl_level_N(  Nepochs, Nregions, max_obs_kinds), &
               Nrejected(  Nepochs, Nregions, max_obs_kinds), &
             NbadQClevel(  Nepochs, Nregions, max_obs_kinds, 2))

   rms_ges_mean   = 0
   rms_ges_spread = 0
   rms_anl_mean   = 0
   rms_anl_spread = 0
   Nrejected      = 0
   NbadQClevel    = 0
   ges_level_N    = 0
   anl_level_N    = 0

   !--------------------------------------------------------------------
   ! Navigate the entire obs sequence once to determine the unique times
   !--------------------------------------------------------------------
   observation       = obs1
   iepoch            = 1
   bincenter(iepoch) = seqT1
   call get_time(bincenter(iepoch),seconds,days)
   epoch_center(iepoch) = days + seconds/86400.0_digits12

   write(msgstring,'(''epoch '',i4,'' center '')')iepoch
   call print_time(bincenter(iepoch),trim(adjustl(msgstring)),logfileunit)

   BinLoop : do obsindex = 1,max_num_obs

      call get_next_obs(seq, observation, next_obs, is_this_last)
      if ( is_this_last) exit BinLoop

      call get_obs_def(next_obs, obs_def)
      obs_time = get_obs_def_time(obs_def)

      if ( obs_time /= bincenter(iepoch) ) then
         if (iepoch == Nepochs) then
            ! Just making sure we stay within array bounds.
            write(msgstring,*)'no room for another bin center'
            call error_handler(E_ERR,'obs_diag',msgstring,source,revision,revdate)
         endif

         iepoch = iepoch + 1
         bincenter(iepoch) = obs_time
         call get_time(bincenter(iepoch),seconds,days)
         epoch_center(iepoch) = days + seconds/86400.0_digits12

         write(msgstring,'(''epoch '',i4,'' center '')')iepoch
         call print_time(bincenter(iepoch),trim(adjustl(msgstring)),logfileunit)
!        call print_date(bincenter(iepoch),trim(adjustl(msgstring)),logfileunit)
      endif

      observation = next_obs

   enddo BinLoop

   if (iepoch /= Nepochs) then
      write(msgstring,*)'Nepochs is ',Nepochs,' but only found ',&
                         iepoch,' unique times in sequence.'
      call error_handler(E_ERR,'obs_diag',msgstring,source,revision,revdate)
   endif

   if (bin_width_seconds <= 0) then  ! no user guidance.
      binsep       = bincenter(2) - bincenter(1)
      binwidth     = binsep
      halfbinwidth = binwidth / 2
   else
      binsep       = bincenter(2) - bincenter(1)
      binwidth     = set_time(bin_width_seconds,0)
      halfbinwidth = binwidth / 2
      call get_time(binwidth,bin_width_seconds,days) 
   endif

   call print_time(      binsep,'binsep       ')
   call print_time(halfbinwidth,'halfbinwidth ')
   call print_time(      binsep,'binsep       ',logfileunit)
   call print_time(halfbinwidth,'halfbinwidth ',logfileunit)

   !--------------------------------------------------------------------
   ! Find the index of obs, ensemble mean, spread ... etc.
   !--------------------------------------------------------------------

   obs_copy_index         = -1
   prior_mean_index       = -1
   posterior_mean_index   = -1
   prior_spread_index     = -1
   posterior_spread_index = -1
   qc_index               = -1
   dart_qc_index          = -1

   MetaDataLoop : do i=1, get_num_copies(seq)
      if(index(get_copy_meta_data(seq,i), 'observation'              ) > 0) &
                          obs_copy_index = i
      if(index(get_copy_meta_data(seq,i), 'prior ensemble mean'      ) > 0) &
                   prior_mean_index = i
      if(index(get_copy_meta_data(seq,i), 'posterior ensemble mean'  ) > 0) &
               posterior_mean_index = i
      if(index(get_copy_meta_data(seq,i), 'prior ensemble spread'    ) > 0) &
                 prior_spread_index = i
      if(index(get_copy_meta_data(seq,i), 'posterior ensemble spread') > 0) &
             posterior_spread_index = i
   enddo MetaDataLoop

   QCMetaDataLoop : do i=1, get_num_qc(seq)
      if(index(  get_qc_meta_data(seq,i), 'Quality Control'          ) > 0) &
                           qc_index = i
      if(index(  get_qc_meta_data(seq,i), 'DART quality control'     ) > 0) &
                      dart_qc_index = i
   enddo QCMetaDataLoop

   !--------------------------------------------------------------------
   ! Make sure we find an index for each of them.
   !--------------------------------------------------------------------

   if (         obs_copy_index < 0 ) then
      write(msgstring,*)'metadata:observation not found'
      call error_handler(E_MSG,'obs_diag',msgstring,source,revision,revdate)
   endif
   if (       prior_mean_index < 0 ) then
      write(msgstring,*)'metadata:prior ensemble mean not found'
      call error_handler(E_MSG,'obs_diag',msgstring,source,revision,revdate)
   endif
   if (   posterior_mean_index < 0 ) then 
      write(msgstring,*)'metadata:posterior ensemble mean not found' 
      call error_handler(E_MSG,'obs_diag',msgstring,source,revision,revdate) 
   endif
   if (     prior_spread_index < 0 ) then 
      write(msgstring,*)'metadata:prior ensemble spread not found'
      call error_handler(E_MSG,'obs_diag',msgstring,source,revision,revdate) 
   endif
   if ( posterior_spread_index < 0 ) then 
      write(msgstring,*)'metadata:posterior ensemble spread not found' 
      call error_handler(E_MSG,'obs_diag',msgstring,source,revision,revdate) 
   endif
   if (               qc_index < 0 ) then 
      write(msgstring,*)'metadata:Quality Control not found' 
      call error_handler(E_MSG,'obs_diag',msgstring,source,revision,revdate) 
   endif
   if ( any( (/obs_copy_index, prior_mean_index, posterior_mean_index, & 
               prior_spread_index, posterior_spread_index, &
               qc_index /) < 0) ) then
      write(msgstring,*)'metadata incomplete'
      call error_handler(E_ERR,'obs_diag',msgstring,source,revision,revdate)
   endif

   if (          dart_qc_index < 0 ) then 
      write(msgstring,*)'metadata:DART quality control not found' 
      call error_handler(E_MSG,'obs_diag',msgstring,source,revision,revdate) 
   endif

   !--------------------------------------------------------------------
   ! Echo what we found.
   !--------------------------------------------------------------------

   write(msgstring,'(''observation          index '',i2,'' metadata '',a)') &
        obs_copy_index, trim(adjustl(get_copy_meta_data(seq,obs_copy_index)))
   call error_handler(E_MSG,'obs_diag',msgstring,source,revision,revdate)

   write(msgstring,'(''prior mean           index '',i2,'' metadata '',a)') &
        prior_mean_index, trim(adjustl(get_copy_meta_data(seq,prior_mean_index)))
   call error_handler(E_MSG,'obs_diag',msgstring,source,revision,revdate)

   write(msgstring,'(''posterior mean       index '',i2,'' metadata '',a)') &
        posterior_mean_index, trim(adjustl(get_copy_meta_data(seq,posterior_mean_index)))
   call error_handler(E_MSG,'obs_diag',msgstring,source,revision,revdate) 

   write(msgstring,'(''prior spread         index '',i2,'' metadata '',a)') &
        prior_spread_index, trim(adjustl(get_copy_meta_data(seq,prior_spread_index)))
   call error_handler(E_MSG,'obs_diag',msgstring,source,revision,revdate)

   write(msgstring,'(''posterior spread     index '',i2,'' metadata '',a)') &
        posterior_spread_index, trim(adjustl(get_copy_meta_data(seq,posterior_spread_index)))
   call error_handler(E_MSG,'obs_diag',msgstring,source,revision,revdate)

   write(msgstring,'(''Quality Control      index '',i2,'' metadata '',a)') &
        qc_index,      trim(adjustl(get_qc_meta_data(seq,     qc_index)))
   call error_handler(E_MSG,'obs_diag',msgstring,source,revision,revdate)

   if (          dart_qc_index > 0 ) then 
   write(msgstring,'(''DART quality control index '',i2,'' metadata '',a)') &
        dart_qc_index, trim(adjustl(get_qc_meta_data(seq,dart_qc_index)))
   call error_handler(E_MSG,'obs_diag',msgstring,source,revision,revdate)
   endif

   !====================================================================
   EpochLoop : do iepoch = 1, Nepochs
   !====================================================================

      ! If the start time was 0,0 ... subtracting the halfbinwidth
      ! would result in a negative time, which is not allowed. Bad
      ! things happened. Since, by definition, the first observation
      ! is at the earliest time, we lose nothing by using it for the 
      ! start time of the first epoch.
      if (iepoch == 1) then
         beg_time = seqT1
      else
         beg_time = bincenter(iepoch) - halfbinwidth + set_time(1,0)
      endif
      end_time     = bincenter(iepoch) + halfbinwidth

      call get_obs_time_range(seq, beg_time, end_time, key_bounds, &
                  num_obs_in_epoch, out_of_range )

      if( num_obs_in_epoch == 0 ) then
         if ( verbose ) then
            call print_time(         beg_time,' epoch  start ',logfileunit)
            call print_time(         beg_time,' epoch  start ')
            call print_time(bincenter(iepoch),' epoch center ',logfileunit)
            call print_time(bincenter(iepoch),' epoch center ')
            call print_time(         end_time,' epoch    end ',logfileunit)
            call print_time(         end_time,' epoch    end ')
            write(logfileunit,*)' No observations in epoch ',iepoch,' cycling ...'
            write(     *     ,*)' No observations in epoch ',iepoch,' cycling ...'
         endif
         cycle EpochLoop
      endif

      if ( verbose ) then
         call print_time(         beg_time,' epoch  start ',logfileunit)
         call print_time(bincenter(iepoch),' epoch center ',logfileunit)
         call print_time(         end_time,' epoch    end ',logfileunit)
         write(logfileunit, *)'num_obs_in_epoch (', iepoch, ') = ', num_obs_in_epoch

         call print_time(         beg_time,' epoch  start ')
         call print_time(bincenter(iepoch),' epoch center ')
         call print_time(         end_time,' epoch    end ')
         write(     *     , *)'num_obs_in_epoch (', iepoch, ') = ', num_obs_in_epoch
      endif

      allocate(keys(num_obs_in_epoch))

      call get_time_range_keys(seq, key_bounds, num_obs_in_epoch, keys)

      !-----------------------------------------------------------------
      ObservationLoop : do obsindex = 1, num_obs_in_epoch
      !-----------------------------------------------------------------

         call get_obs_from_key(seq, keys(obsindex), observation) 
         call get_obs_def(observation, obs_def)
         obs_time    = get_obs_def_time(obs_def)
         flavor      = get_obs_kind(obs_def) ! this is (almost) always [1,max_obs_kinds]
         obs_err_var = get_obs_def_error_variance(obs_def) 
         obs_loc     = get_obs_def_location(obs_def)
         rlocation   = get_location(obs_loc) 

         call get_obs_values(observation,              obs,         obs_copy_index)
         call get_obs_values(observation,       prior_mean,       prior_mean_index)
         call get_obs_values(observation,   posterior_mean,   posterior_mean_index)
         call get_obs_values(observation,     prior_spread,     prior_spread_index)
         call get_obs_values(observation, posterior_spread, posterior_spread_index)

         pr_mean = prior_mean(1)
         po_mean = posterior_mean(1)
         pr_sprd = prior_spread(1)
         po_sprd = posterior_spread(1)

         ! Convert the DART QC data to an integer and create histogram 

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
         ! A Whole bunch of reasons to be rejected
         !--------------------------------------------------------------

         if ( flavor < 0 ) then
         !  write(*,*)'obs ',obsindex,' is an identity observation - no fair.',obs_err_var
            Nidentity = Nidentity + 1
            cycle ObservationLoop
         endif

         keeper = .true.
!        keeper = CheckObsType(obs_select, obs_err_var, flavor, ipressure)
         if ( .not. keeper ) then
            write(*,*)'obs ',obsindex,' rejected by CheckObsType ',obs_err_var
            NwrongType = NwrongType + 1
            cycle ObservationLoop
         endif

         if( qc(qc_index) >= input_qc_threshold ) then
         !  write(*,*)'obs ',obsindex,' rejected by qc ',qc
            NbadQC = NbadQC + 1
            cycle ObservationLoop 
         endif

         ! update the histogram of the magnitude of the innovation,
         ! where each bin is a single standard deviation. This is 
         ! a one-sided histogram.

         ratio = GetRatio(obs(1), pr_mean, pr_sprd, obs_err_var, qc_integer)
         indx         = min(int(ratio), MaxSigmaBins)
         nsigma(indx) = nsigma(indx) + 1

         obs_used_in_epoch(iepoch) = obs_used_in_epoch(iepoch) + 1

         !--------------------------------------------------------------
         ! We have four regions of interest
         !--------------------------------------------------------------

         Areas : do iregion =1, Nregions

            keeper = InRegion( rlocation, lonlim1(iregion), lonlim2(iregion) )
            if ( .not. keeper ) cycle Areas

            !-----------------------------------------------------------
            ! Check for desired type ... joke
            !-----------------------------------------------------------

!           if (flavor > 0 ) then  ! keep all types, for now ...

               ratio = GetRatio(obs(1), pr_mean, pr_sprd, obs_err_var, qc_integer)

               if ( ratio > rat_cri ) then
                  if (verbose) then
                     write(logfileunit,*)'obsindex ',obsindex,' ratio ', ratio
                     write(logfileunit,*)'val,prm,prs,var ',obs(1), &
                                          pr_mean, pr_sprd, obs_err_var
                  endif
                  Nrejected(iepoch, iregion, flavor) = Nrejected(iepoch, iregion, flavor) + 1 
                  cycle Areas
               endif

               if ( qc_integer > QC_MAX_PRIOR ) then  ! prior and posterior failed

                  NbadQClevel(iepoch, iregion, flavor, 1) = &
                  NbadQClevel(iepoch, iregion, flavor, 1) + 1 
                  NbadQClevel(iepoch, iregion, flavor, 2) = &
                  NbadQClevel(iepoch, iregion, flavor, 2) + 1 

               else if ( qc_integer > QC_MAX_POSTERIOR ) then

                  ! Then at least the prior (A.K.A. guess) is good

                  ges_level_N(   iepoch, iregion, flavor) = &
                  ges_level_N(   iepoch, iregion, flavor) + 1

                  rms_ges_mean(  iepoch, iregion, flavor) = &
                  rms_ges_mean(  iepoch, iregion, flavor) + (pr_mean - obs(1))**2

                  rms_ges_spread(iepoch, iregion, flavor) = &
                  rms_ges_spread(iepoch, iregion, flavor) + pr_sprd**2

                  ! However, the posterior is bad

                  NbadQClevel(iepoch, iregion, flavor, 2) = &
                  NbadQClevel(iepoch, iregion, flavor, 2) + 1 

               else

                  ! The prior is good

                  ges_level_N(   iepoch, iregion, flavor) = &
                  ges_level_N(   iepoch, iregion, flavor) + 1

                  rms_ges_mean(  iepoch, iregion, flavor) = &
                  rms_ges_mean(  iepoch, iregion, flavor) + (pr_mean - obs(1))**2

                  rms_ges_spread(iepoch, iregion, flavor) = &
                  rms_ges_spread(iepoch, iregion, flavor) + pr_sprd**2

                  ! The posterior is good

                  anl_level_N(   iepoch, iregion, flavor) = &
                  anl_level_N(   iepoch, iregion, flavor) + 1

                  rms_anl_mean(  iepoch, iregion, flavor) = &
                  rms_anl_mean(  iepoch, iregion, flavor) + (po_mean - obs(1))**2

                  rms_anl_spread(iepoch, iregion, flavor) = &
                  rms_anl_spread(iepoch, iregion, flavor) + po_sprd**2

               endif

        !   endif

         enddo Areas

      !-----------------------------------------------------------------
      enddo ObservationLoop
      !-----------------------------------------------------------------

      deallocate(keys)

!     if(verbose) then
!        write(msgstring,'(''num obs considered in epoch '',i4,'' = '',i8)') &
!                        iepoch, obs_used_in_epoch(iepoch)
!        call error_handler(E_MSG,'obs_diag',msgstring,source,revision,revdate)
!        write(logfileunit,*)''
!        write(     *     ,*)''
!     endif

   enddo EpochLoop

   if (verbose) then
      write(logfileunit,*)'End of EpochLoop for ',trim(adjustl(obs_seq_in_file_name))
      write(     *     ,*)'End of EpochLoop for ',trim(adjustl(obs_seq_in_file_name))
   endif

   call destroy_obs(obs1)
   call destroy_obs(obsN)
   call destroy_obs(observation)
   call destroy_obs(next_obs)
   call destroy_obs_sequence(seq)

!enddo ObsFileLoop

!-----------------------------------------------------------------------
! We have read all possible files, and stuffed the observations into the
! appropriate bins. Time to normalize and finish up. Almost.
!-----------------------------------------------------------------------
! First, echo attributes to a file to facilitate plotting.
! This file is also a matlab function ... the variables are
! loaded by just typing the name at the matlab prompt.
! The output file names are included in this master file, 
! but they are created on-the-fly, so we must leave this open till then.
!-----------------------------------------------------------------------

if (sum(obs_used_in_epoch) == 0 ) then
   call error_handler(E_ERR,'obs_diag','All identity observations. Stopping.', &
                     source, revision, revdate)
endif

iunit = open_file('ObsDiagAtts.m',form='formatted',action='rewind')
write(iunit,'(''iskip_days     = '',i6,'';'')')iskip_days
write(iunit,'(''obs_select     = '',i6,'';'')')obs_select
write(iunit,'(''rat_cri        = '',f9.2,'';'')')rat_cri
write(iunit,'(''qc_threshold   = '',f9.2,'';'')')input_qc_threshold
write(iunit,'(''bin_width_seconds = '',i5,'';'')')bin_width_seconds
write(iunit,'(''t1             = '',f20.6,'';'')')epoch_center(1)
write(iunit,'(''tN             = '',f20.6,'';'')')epoch_center(Nepochs)
do iregion=1,Nregions
   write(iunit,'(''lonlim1('',i3,'') = '',f9.2,'';'')')iregion, lonlim1(iregion)
   write(iunit,'(''lonlim2('',i3,'') = '',f9.2,'';'')')iregion, lonlim2(iregion)
   write(iunit,47)iregion,trim(adjustl(reg_names(iregion)))
enddo
47 format('Regions(',i3,') = {''',a,'''};')

write(*,*)''
do iepoch = 1, Nepochs
   write(msgstring,'(''num obs used in epoch '',i5,'' = '',i8)') iepoch, obs_used_in_epoch(iepoch)
   call error_handler(E_MSG,'obs_diag',msgstring,source,revision,revdate)
enddo

if (verbose) then
   write(logfileunit,*)'Normalizing time-region-variable quantities for desired level.'
   write(    *      ,*)'Normalizing time-region-variable quantities for desired level.'
endif

ivarcount = 0

OneLevel : do ivar=1,max_obs_kinds
   do iregion=1, Nregions
   do iepoch=1, Nepochs

      ! prior ...

      if ( ges_level_N(iepoch, iregion, ivar) == 0) then
          rms_ges_mean(iepoch, iregion, ivar) = -99.0_r8
        rms_ges_spread(iepoch, iregion, ivar) = -99.0_r8
      else
          rms_ges_mean(iepoch, iregion, ivar) = sqrt(   rms_ges_mean(iepoch, iregion, ivar) / &
                                                         ges_level_N(iepoch, iregion, ivar) )
        rms_ges_spread(iepoch, iregion, ivar) = sqrt( rms_ges_spread(iepoch, iregion, ivar) / &
                                                         ges_level_N(iepoch, iregion, ivar) )
      endif

      ! posterior ...

      if ( anl_level_N(iepoch, iregion, ivar) == 0) then
          rms_anl_mean(iepoch, iregion, ivar) = -99.0_r8
        rms_anl_spread(iepoch, iregion, ivar) = -99.0_r8
      else
          rms_anl_mean(iepoch, iregion, ivar) = sqrt(   rms_anl_mean(iepoch, iregion, ivar) / &
                                                         anl_level_N(iepoch, iregion, ivar) )
        rms_anl_spread(iepoch, iregion, ivar) = sqrt( rms_anl_spread(iepoch, iregion, ivar) / &
                                                         anl_level_N(iepoch, iregion, ivar) )
      endif

   enddo
   enddo

   !--------------------------------------------------------------------
   ! Create data files for all observation kinds we have used
   !--------------------------------------------------------------------

   if( any (ges_level_N(:, :, ivar) > 0) ) then

      ivarcount = ivarcount + 1

      write(gesName,'(a,''_ges_times.dat'')') trim(adjustl(get_obs_name(ivar)))
      gesUnit = open_file(trim(adjustl(gesName)),form='formatted',action='rewind')

      if (verbose) then
         write(logfileunit,*)'Creating '//trim(adjustl(gesName))
         write(     *     ,*)'Creating '//trim(adjustl(gesName))
      endif

      do i=1, Nepochs
         if( any (ges_level_N(i, :, ivar) /= 0) ) then
            call get_time(bincenter(i),seconds,days)
            write(gesUnit,91) days, seconds, &
                 (rms_ges_mean(i,j,ivar),rms_ges_spread(i,j,ivar),ges_level_N(i,j,ivar),j=1,Nregions)
         endif
      enddo
      close(gesUnit)

      write(iunit,95) ivarcount, trim(adjustl(get_obs_name(ivar)))

   endif


   if( any (anl_level_N(:, :, ivar) > 0) ) then

      write(anlName,'(a,''_anl_times.dat'')') trim(adjustl(get_obs_name(ivar)))
      anlUnit = open_file(trim(adjustl(anlName)),form='formatted',action='rewind')

      if (verbose) then
         write(logfileunit,*)'Creating '//trim(adjustl(anlName))
         write(     *     ,*)'Creating '//trim(adjustl(anlName))
      endif

      do i=1, Nepochs
         if( any (anl_level_N(i, :, ivar) /= 0) ) then
            call get_time(bincenter(i),seconds,days)
            write(anlUnit,91) days, seconds, &
                 (rms_anl_mean(i,j,ivar),rms_anl_spread(i,j,ivar),anl_level_N(i,j,ivar),j=1,Nregions)
         endif
      enddo
      close(anlUnit)

      write(iunit,95) ivarcount, trim(adjustl(get_obs_name(ivar)))

   endif

enddo OneLevel

91 format(i7,1x,i5,4(1x,2f7.2,1x,i8))
95 format('One_Level_Varnames(',i3,') = {''',a,'''};')

! Print the histogram of innovations as a function of standard deviation. 
write(*,*) ''
write(*,*) 'Table to indicate how to choose rat_cri -- '
write(*,*) 'How are the innovations distributed?'
do i=0,100
   if(nsigma(i) /= 0) write(*,*)'innovations within ',i+1,' stdev = ',nsigma(i)
enddo
write(*,*) ''

! Print the histogram of DART QC values. 
do i=0,7
   write(*,*)'DART QC value of ',i,' N = ',qc_counter(i)
enddo
write(*,*)'Desired obs with DART QC values, N = ',sum(qc_counter)
write(*,*) ''

!-----------------------------------------------------------------------
close(iunit)   ! Finally close the 'master' matlab diagnostic file.
!-----------------------------------------------------------------------
! Print final rejection summary.
!-----------------------------------------------------------------------

write(*,*) ''
write(*,*) '# observations used  : ',sum(obs_used_in_epoch)
write(*,*) 'Rejected Observations summary.'
write(*,*) '# NwrongType         : ',NwrongType
write(*,*) '# NbadQC             : ',NbadQC
write(*,*) '# NIdentityObs       : ',Nidentity
write(*,*) '# NbadQClevel (prior): ',sum(NbadQClevel(:,:,:,1))
write(*,*) '# NbadQClevel (post) : ',sum(NbadQClevel(:,:,:,2))
write(*,*) '# NbadLevel          : ',NbadLevel
write(*,'(a)')'Table of observations rejected by region:'
write(*,'(5a)')'                                                   ',reg_names(1:Nregions)
do ivar=1,max_obs_kinds
   write(*,'(3a,4i8)') ' # ',get_obs_name(ivar),' (ratio)  : ',sum(Nrejected(:,:,ivar),1)
   write(*,'(3a,4i8)') ' # ',get_obs_name(ivar),' (prior)  : ',sum(NbadQClevel(:,:,ivar,1),1)
   write(*,'(3a,4i8)') ' # ',get_obs_name(ivar),' (post)   : ',sum(NbadQClevel(:,:,ivar,2),1)
enddo
write(*,*)''

write(logfileunit,*) ''
write(logfileunit,*) '# observations used  : ',sum(obs_used_in_epoch)
write(logfileunit,*) 'Rejected Observations summary.'
write(logfileunit,*) '# NwrongType         : ',NwrongType
write(logfileunit,*) '# NbadQC             : ',NbadQC
write(logfileunit,*) '# NIdentityObs       : ',Nidentity
write(logfileunit,*) '# NbadQClevel (prior): ',sum(NbadQClevel(:,:,:,1),1)
write(logfileunit,*) '# NbadQClevel (post) : ',sum(NbadQClevel(:,:,:,2),1)
write(logfileunit,*) '# NbadLevel          : ',NbadLevel
write(logfileunit,'(a)')'Table of observations rejected by region:'
write(logfileunit,'(5a)')' ',reg_names(1:Nregions)
do ivar=1,max_obs_kinds
   write(logfileunit,'(3a,4i8)') ' # ',get_obs_name(ivar),' (ratio)  : ',sum(Nrejected(:,:,ivar),1)
   write(logfileunit,'(3a,4i8)') ' # ',get_obs_name(ivar),' (prior)  : ',sum(NbadQClevel(:,:,ivar,1),1)
   write(logfileunit,'(3a,4i8)') ' # ',get_obs_name(ivar),' (post)   : ',sum(NbadQClevel(:,:,ivar,2),1)
enddo

!-----------------------------------------------------------------------
! Really, really, done.
!-----------------------------------------------------------------------

deallocate(rms_ges_mean, rms_ges_spread, &
           rms_anl_mean, rms_anl_spread, &
           ges_level_N, anL_level_N, Nrejected, NbadQClevel, &
           obs_used_in_epoch, bincenter, epoch_center)

call timestamp(source,revision,revdate,'end') ! That closes the log file, too.

contains

   Function GetRatio(obsval, prmean, prspred, errcov, qcval ) result (ratio)

   ! This function tries to get a handle on the magnitude of the innovations.
   ! If the ratio of the observation to the prior mean is 'big', it is an outlier. 
   ! If the prior mean cannot be calculated (i.e. is missing) we give it a HUGE
   ! innovation -- sufficiently large to put it in the last 'bin' of the crude
   ! histogram.

   real(r8), intent(in) :: obsval, prmean, prspred, errcov
   integer,  intent(in) :: qcval
   real(r8)             :: ratio

   real(r8) :: numer, denom

   if (qcval <= QC_MAX_PRIOR ) then
      numer = abs(prmean- obsval) 
      denom = sqrt( prspred**2 + errcov )
      ratio = numer / denom
   else
      ratio = real(MaxSigmaBins,r8)
   endif

   end Function GetRatio



   Function InRegion( lon, lon1, lon2 ) result( keeper )
   real(r8), intent(in) :: lon, lon1, lon2
   logical :: keeper

   keeper = .false.

   if( (lon .ge. lon1) .and. (lon .lt. lon2) ) keeper = .true.

   end Function InRegion


   Function CheckObsType(obs_select, obs_err_var, flavor, ipressure) result(keeper)
   ! Since the observation kind does not have platform information, Hui
   ! has determined an ad-hoc set of rules to determine the origin of the
   ! observation. Hence the 'magic' numbers.

   integer,  intent(in) :: obs_select 
   real(r8), intent(in) :: obs_err_var 
   integer,  intent(in) :: flavor, ipressure
   logical              :: keeper

   keeper = .true. ! Innocent till proven guilty ...

   end Function CheckObsType

end program obs_diag
