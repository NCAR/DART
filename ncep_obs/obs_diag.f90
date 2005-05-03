! Data Assimilation Research Testbed -- DART
! Copyright 2004, 2005, Data Assimilation Initiative, University Corporation for Atmospheric Research
! Licensed under the GPL -- www.gpl.org/licenses/gpl.html

program obs_diag

! <next five lines automatically updated by CVS, do not edit>
! $Source$
! $Revision$
! $Date$
! $Author$
! $Name$

! At present (April, 2005) this program assumes that one observation sequence
! file contains a DAY of observations. This must certainly change.
!
! The programs defines a series of epochs (periods of time) and geographic
! regions and accumulates statistics for these epochs and regions.

use        types_mod, only : r8
use obs_sequence_mod, only : read_obs_seq, obs_type, obs_sequence_type, get_first_obs, &
                             get_obs_from_key, get_obs_def, get_copy_meta_data, &
                             get_obs_time_range, get_time_range_keys, get_num_obs, &
                             get_next_obs, get_num_times, get_obs_values, init_obs, &
                             assignment(=), get_num_copies, static_init_obs_sequence, &
                             get_qc, destroy_obs_sequence 
use      obs_def_mod, only : obs_def_type, get_obs_def_error_variance, get_obs_def_time, &
                             get_obs_def_location,  get_obs_def_kind 
use     location_mod, only : location_type, get_location, set_location_missing, &
                             write_location, operator(/=)
use     obs_kind_mod, only : KIND_U, KIND_V, KIND_PS, KIND_T, KIND_QV, &
                             obs_kind_type, get_obs_kind 
use time_manager_mod, only : time_type, set_date, set_time, get_time, print_time, &
                             set_calendar_type, operator(*), &
                             operator(+), operator(-), operator(/=), operator(>)
use    utilities_mod, only : get_unit, open_file, close_file, register_module, &
                             check_nml_error, file_exist, error_handler, E_ERR, E_MSG, &
                             initialize_utilities, logfileunit, timestamp

implicit none

! CVS Generated file description for error handling, do not edit
character(len=128) :: &
source   = "$Source$", &
revision = "$Revision$", &
revdate  = "$Date$"

type(obs_sequence_type) :: seq
type(obs_type)          :: observation, next_obs
type(obs_def_type)      :: obs_def
type(obs_kind_type)     :: obs_kind
type(location_type)     :: obs_loc
type(time_type)         :: next_time

!---------------------
integer :: ibintoday, obsindex, i, iunit, ierr, io
integer :: num_obs_in_set, obs_used
integer :: NbadUV = 0

! Storage with fixed size for observation space diagnostics
real(r8) :: prior_mean(1), posterior_mean(1)
real(r8) :: prior_spread(1), posterior_spread(1)

! We are treating winds as a vector pair, but we are handling the
! observations serially. Consequenctly, we exploit the fact that
! the U observations are followed by the V observations. If we
! keep the 'last' obs -- when we hit a V obs, the U is still available.
type(location_type) :: U_obs_loc
integer  :: U_flavor
real(r8) :: U_prior_mean(1), U_posterior_mean(1)
real(r8) :: U_prior_spread(1), U_posterior_spread(1)
real(r8) :: U_obs, U_obs_err_cov

integer :: obs_index, prior_mean_index, posterior_mean_index
integer :: prior_spread_index, posterior_spread_index
integer :: key_bounds(2), flavor

real(r8), allocatable :: obs_err_cov(:), obs(:), qc(:)
integer,  allocatable :: keys(:)

logical :: out_of_range, is_there_one, is_this_last, keeper

!----------------------------------------------------------------
! Namelists with default values
!
real(r8) :: outlier_threshold = -1.0_r8
character(len = 129) :: obs_sequence_name = "obs_seq.final"

namelist /filter_nml/ obs_sequence_name, outlier_threshold

integer :: obs_year   = 2003     ! the first date of the diagnostics
integer :: obs_month  = 1
integer :: obs_day    = 1
integer :: tot_days   = 1        ! total days
integer :: iskip      = 0        ! skip the first 'iskip' days
integer :: level      = 500      ! level (hPa)
integer :: obs_select = 1        ! obs type selection: 1=all, 2 =RAonly, 3=noRA
real(r8):: rat_cri    = 5.0      ! QC ratio
real(r8):: bin_separation = 6.0  ! Bins every so often (hours)
real(r8):: bin_width = 6.0       ! width of the bin (hour)

namelist /obsdiag_nml/ obs_year, obs_month, obs_day, tot_days, &
                       iskip, level, obs_select, rat_cri, &
                       bin_separation, bin_width

!--------------------------------------------------------------------------
! Spatial
! Each observation kind gets its own mean, spread, for Guess/Analysis
!--------------------------------------------------------------------------
! index 1 == region 1 == Northern Hemisphere
! index 2 == region 2 == Southern Hemisphere
! index 3 == region 3 == Tropics
! index 4 == region 4 == North America
! TJH - some kind of crazy nomenclature that South Pole = lat 0?

integer, parameter :: Nregions = 4 

real(r8) :: lonlim1(Nregions), lonlim2(Nregions), &
            latlim1(Nregions), latlim2(Nregions)

data lonlim1 /   0.0_r8,   0.0_r8,   0.0_r8, 235.0_r8 /
data lonlim2 / 360.0_r8, 360.0_r8, 360.0_r8, 295.0_r8 /
data latlim1 / 110.0_r8,  10.0_r8,  70.0_r8, 115.0_r8 /
data latlim2 / 170.0_r8,  70.0_r8, 110.0_r8, 145.0_r8 /

integer  :: iregion, iepoch, iday
real(r8) :: lon0, lat0, obsloc3(3)
real(r8) :: speed_obs2, speed_ges2, speed_anl2

!---------------------------------------------------------
! Vertical
!---------------------------------------------------------

integer, parameter :: nlev=11
integer  :: ipressure, plev(nlev), pint(nlev+1)

real(r8) :: rms_ges_ver_W(nlev, Nregions),  rms_anl_ver_W(nlev, Nregions)
real(r8) :: rms_ges_ver_T(nlev, Nregions),  rms_anl_ver_T(nlev, Nregions)
real(r8) :: rms_ges_ver_Q(nlev, Nregions),  rms_anl_ver_Q(nlev, Nregions)

real(r8) :: bias_ges_ver_W(nlev, Nregions), bias_anl_ver_W(nlev, Nregions)
real(r8) :: bias_ges_ver_T(nlev, Nregions), bias_anl_ver_T(nlev, Nregions)
real(r8) :: bias_ges_ver_Q(nlev, Nregions), bias_anl_ver_Q(nlev, Nregions)
integer  :: num_ver_W(nlev, Nregions), num_ver_T(nlev, Nregions), num_ver_Q(nlev, Nregions)

integer  :: ilev          ! counter
integer  :: level_index   ! index of pressure level closest to input 
data plev / 1000, 925, 850, 700, 500, 400, 300, 250, 200, 150, 100/
data pint / 1025, 950, 900, 800, 600, 450, 350, 275, 225, 175, 125, 75/

!-------------------------------------------------------------------------
! Spatio-Temporal Variables
! Dimension 1 is temporal, actually - these are time-by-region- 
!-------------------------------------------------------------------------

integer,  allocatable, dimension(:,:) :: num_in_level_W, &
                                         num_in_level_T, &
                                         num_in_level_Q, &
                                         num_in_level_P

real(r8), allocatable, dimension(:,:) :: rms_ges_mean_W, rms_ges_spread_W, &
                                         rms_anl_mean_W, rms_anl_spread_W, &
                                         rms_ges_mean_T, rms_ges_spread_T, &
                                         rms_anl_mean_T, rms_anl_spread_T, &
                                         rms_ges_mean_Q, rms_ges_spread_Q, &
                                         rms_anl_mean_Q, rms_anl_spread_Q, &
                                         rms_ges_mean_P, rms_ges_spread_P, &
                                         rms_anl_mean_P, rms_anl_spread_P 

real(r8), allocatable, dimension(:) :: epoch_center

!-------------------------------------------------------------------------
! General purpose variables
!-------------------------------------------------------------------------

integer  :: seconds, days, DayOne
integer  :: k0, kkk, NBinsPerDay, Nepochs
integer  :: calendar_type
integer  :: WgesUnit, WanlUnit, TgesUnit, TanlUnit
integer  :: QgesUnit, QanlUnit, PgesUnit, PanlUnit

real(r8) :: numer, denom, ratio

type(time_type) :: beg_time, end_time
type(time_type) :: bincenter, binwidth, halfbinwidth 

character(len =   6) :: day_num 
character(len = 129) :: WgesName, WanlName, TgesName, TanlName, msgstring
character(len = 129) :: QgesName, QanlName, PgesName, PanlName

!----------------------------------------------------------------

call initialize_utilities('obs_diag')
call register_module(source,revision,revdate) 
call static_init_obs_sequence()  ! Initialize the obs sequence module 
call init_obs(observation, 0, 0) ! Initialize the observation type variables
call init_obs(   next_obs, 0, 0)

! Begin by reading the namelist input for obs_diag

if(file_exist('input.nml')) then
   iunit = open_file('input.nml', action = 'read')
   ierr = 1

   do while(ierr /= 0)
      read(iunit, nml = obsdiag_nml, iostat = io, end = 11)
      ierr = check_nml_error(io, 'obsdiag_nml')
   enddo
 11 continue
   call close_file(iunit)
endif
write(logfileunit,nml=obsdiag_nml)

! Reading the filter namelist for obs sequence filename.

if(file_exist('input.nml')) then
   iunit = open_file('input.nml', action = 'read')
   ierr = 1
   do while(ierr /= 0)
        read(iunit, nml = filter_nml, iostat = io, end = 12)
        ierr = check_nml_error(io, 'filter_nml')
   enddo
12 continue
   call close_file(iunit)
endif
write(logfileunit,nml=filter_nml)

! Now that we have input, do some checking and setup

level_index = GetClosestLevel(level) 

calendar_type = 3
call set_calendar_type(calendar_type)

NBinsPerDay  = nint( 24.0_r8 / bin_separation )
    binwidth = set_time(nint(bin_width * 3600.0_r8), 0)    ! full bin width 
halfbinwidth = set_time(nint(bin_width * 1800.0_r8), 0)    ! half bin width 

! Echo to a resource file to facilitate plotting

iunit = open_file('Tanl_times_level.dat',form='formatted',action='rewind')
write(iunit,*) plev(level_index), tot_days, iskip
close(iunit)

print*,'tot_days,iskip,level= ',tot_days,iskip,plev(level_index),'mb'
print*,'OBS type, QC ratio = ', obs_select, rat_cri

! Initialize.

prior_mean(1)       = 0.0_r8
prior_spread(1)     = 0.0_r8
posterior_mean(1)   = 0.0_r8
posterior_spread(1) = 0.0_r8

rms_ges_ver_W = 0.0_r8
rms_anl_ver_W = 0.0_r8
rms_ges_ver_T = 0.0_r8
rms_anl_ver_T = 0.0_r8
rms_ges_ver_Q = 0.0_r8
rms_anl_ver_Q = 0.0_r8

bias_ges_ver_W = 0.0_r8
bias_anl_ver_W = 0.0_r8
bias_ges_ver_T = 0.0_r8
bias_anl_ver_T = 0.0_r8
bias_ges_ver_Q = 0.0_r8
bias_anl_ver_Q = 0.0_r8

num_ver_W = 0
num_ver_T = 0
num_ver_Q = 0

!----------------------------------------------------------------------------------
! Nepochs is the total number of time intervals of the period requested
!----------------------------------------------------------------------------------
Nepochs = NBinsPerDay*tot_days

allocate(rms_ges_mean_W(Nepochs, Nregions), rms_ges_spread_W(Nepochs, Nregions), &
         rms_anl_mean_W(Nepochs, Nregions), rms_anl_spread_W(Nepochs, Nregions), &
         rms_ges_mean_T(Nepochs, Nregions), rms_ges_spread_T(Nepochs, Nregions), &
         rms_anl_mean_T(Nepochs, Nregions), rms_anl_spread_T(Nepochs, Nregions), &
         rms_ges_mean_Q(Nepochs, Nregions), rms_ges_spread_Q(Nepochs, Nregions), &
         rms_anl_mean_Q(Nepochs, Nregions), rms_anl_spread_Q(Nepochs, Nregions), &
         rms_ges_mean_P(Nepochs, Nregions), rms_ges_spread_P(Nepochs, Nregions), &
         rms_anl_mean_P(Nepochs, Nregions), rms_anl_spread_P(Nepochs, Nregions))

allocate(num_in_level_W(Nepochs, Nregions), &
         num_in_level_T(Nepochs, Nregions), &
         num_in_level_Q(Nepochs, Nregions), &
         num_in_level_P(Nepochs, Nregions))

allocate(epoch_center(Nepochs))

rms_ges_mean_W   = 0.0_r8
rms_anl_mean_W   = 0.0_r8
rms_ges_mean_T   = 0.0_r8
rms_anl_mean_T   = 0.0_r8
rms_ges_mean_Q   = 0.0_r8
rms_anl_mean_Q   = 0.0_r8
rms_ges_mean_P   = 0.0_r8
rms_anl_mean_P   = 0.0_r8

rms_ges_spread_W = 0.0_r8
rms_anl_spread_W = 0.0_r8
rms_ges_spread_T = 0.0_r8
rms_anl_spread_T = 0.0_r8
rms_ges_spread_Q = 0.0_r8
rms_anl_spread_Q = 0.0_r8
rms_ges_spread_P = 0.0_r8
rms_anl_spread_P = 0.0_r8

num_in_level_W   = 0
num_in_level_T   = 0
num_in_level_Q   = 0
num_in_level_P   = 0

iepoch = 0  ! epoch counter

!--------------------------------------------------------------
DayLoop : do iday=1, tot_days
!--------------------------------------------------------------

   ! Directory/file names are similar to    01_03/obs_seq.final

   write(day_num, '(i2.2,''_'',i2.2,''/'')') obs_month, obs_day + (iday-1) 
   write(msgstring,*)'opening ', day_num, trim(obs_sequence_name)
   call error_handler(E_MSG,'obs_diag',msgstring,source,revision,revdate)

   ! Read in with enough space for diagnostic output values

   call read_obs_seq(day_num//obs_sequence_name, 0, 0, 0, seq)

   print*,'Number of times = ',get_num_times(seq)

   write(msgstring,*)'get_num_copies = ', get_num_copies(seq)
   call error_handler(E_MSG,'obs_diag',msgstring,source,revision,revdate)

   ! Find the index of obs, ensemble mean, spread ... etc.

   obs_index              = -1
   prior_mean_index       = -1
   posterior_mean_index   = -1
   prior_spread_index     = -1
   posterior_spread_index = -1

   MetaDataLoop : do i=1, get_num_copies(seq)
      if(index(get_copy_meta_data(seq,i), 'observation'              ) > 0) &
                          obs_index = i
      if(index(get_copy_meta_data(seq,i), 'prior ensemble mean'      ) > 0) &
                   prior_mean_index = i
      if(index(get_copy_meta_data(seq,i), 'posterior ensemble mean'  ) > 0) &
               posterior_mean_index = i
      if(index(get_copy_meta_data(seq,i), 'prior ensemble spread'    ) > 0) &
                 prior_spread_index = i
      if(index(get_copy_meta_data(seq,i), 'posterior ensemble spread') > 0) &
             posterior_spread_index = i
   enddo MetaDataLoop

   ! Make sure we find an index for each of them.

   if ( obs_index              < 0 ) then
      write(msgstring,*)'metadata:observation not found'
      call error_handler(E_MSG,'obs_diag',msgstring,source,revision,revdate)
   endif
   if ( prior_mean_index       < 0 ) then
      write(msgstring,*)'metadata:prior ensemble mean not found'
      call error_handler(E_MSG,'obs_diag',msgstring,source,revision,revdate)
   endif
   if ( posterior_mean_index   < 0 ) then 
      write(msgstring,*)'metadata:posterior ensemble mean not found' 
      call error_handler(E_MSG,'obs_diag',msgstring,source,revision,revdate) 
   endif
   if ( prior_spread_index     < 0 ) then 
      write(msgstring,*)'metadata:prior ensemble spread not found'
      call error_handler(E_MSG,'obs_diag',msgstring,source,revision,revdate) 
   endif
   if ( posterior_spread_index < 0 ) then 
      write(msgstring,*)'metadata:posterior ensemble spread not found' 
      call error_handler(E_MSG,'obs_diag',msgstring,source,revision,revdate) 
   endif
   if ( any( (/obs_index, prior_mean_index, posterior_mean_index, & 
               prior_spread_index, posterior_spread_index /) < 0) ) then
      write(msgstring,*)'metadata incomplete'
      call error_handler(E_ERR,'obs_diag',msgstring,source,revision,revdate)
   endif

   ! Echo what we found.

   write(msgstring,'(''observation      index '',i2,'' metadata '',a)') &
        obs_index, trim(adjustl(get_copy_meta_data(seq,obs_index)))
   call error_handler(E_MSG,'obs_diag',msgstring,source,revision,revdate)

   write(msgstring,'(''prior mean       index '',i2,'' metadata '',a)') &
        prior_mean_index, trim(adjustl(get_copy_meta_data(seq,prior_mean_index)))
   call error_handler(E_MSG,'obs_diag',msgstring,source,revision,revdate)

   write(msgstring,'(''posterior mean   index '',i2,'' metadata '',a)') &
        posterior_mean_index, trim(adjustl(get_copy_meta_data(seq,posterior_mean_index)))
   call error_handler(E_MSG,'obs_diag',msgstring,source,revision,revdate) 

   write(msgstring,'(''prior spread     index '',i2,'' metadata '',a)') &
        prior_spread_index, trim(adjustl(get_copy_meta_data(seq,prior_spread_index)))
   call error_handler(E_MSG,'obs_diag',msgstring,source,revision,revdate)

   write(msgstring,'(''posterior spread index '',i2,'' metadata '',a)') &
        posterior_spread_index, trim(adjustl(get_copy_meta_data(seq,posterior_spread_index)))
   call error_handler(E_MSG,'obs_diag',msgstring,source,revision,revdate)

   ! Get the time of the first observation in the sequence
   is_there_one = get_first_obs(seq, observation)
   if ( is_there_one /= .TRUE. ) then
      call error_handler(E_ERR,'obs_diag','No Observations in sequence.', &
      source,revision,revdate)
   endif

   call get_obs_def(observation, obs_def)
   next_time    = get_obs_def_time(obs_def)
   call get_time(next_time, seconds, DayOne)   ! ultimately want to throw away seconds
   bincenter = set_time(0,DayOne) + binwidth   ! 06Z, 12Z, 18Z ... for example

   !====================================================================
   Advancesets : do ibintoday = 1, NBinsPerDay
   !====================================================================

      iepoch = iepoch + 1

      if ( ibintoday > 1 ) then ! Get the next time (if any) in the obs sequence
         call get_next_obs(seq, observation, next_obs, is_this_last)

         if( is_this_last ) exit Advancesets

         call get_obs_def(next_obs, obs_def)
         next_time = get_obs_def_time(obs_def)
      endif

      ! set bin begin and end time 

      bincenter = set_time(0,DayOne) + ibintoday * binwidth
      beg_time  = bincenter - halfbinwidth
      end_time  = bincenter + halfbinwidth

      call get_time(bincenter,seconds,days)
      epoch_center(iepoch) = days + seconds/86400.0_r8

      call get_obs_time_range(seq, beg_time, end_time, key_bounds, &
                  num_obs_in_set, out_of_range, observation)

      obs_used = 0
      call print_time( beg_time,'bin  start ',logfileunit)
      call print_time(bincenter,'bin center ',logfileunit)
      call print_time( end_time,'bin    end ',logfileunit)
      write(logfileunit, *) 'num_obs_in_bin (', iepoch, ') = ', num_obs_in_set
      write(     *     , *) 'num_obs_in_bin (', iepoch, ') = ', num_obs_in_set

      allocate(keys(num_obs_in_set), obs_err_cov(num_obs_in_set), &
                obs(num_obs_in_set), qc(num_obs_in_set))

      call get_time_range_keys(seq, key_bounds, num_obs_in_set, keys)
      call print_time(next_time)

      obs(1)         = 0.0_r8
      obs_err_cov(1) = 0.0_r8
      obs_loc        = set_location_missing()
      flavor         = KIND_U

      ! Get the observation value and interpolated ones at obs locations
      !-----------------------------------------------------------------
      ObservationLoop : do obsindex = 1, num_obs_in_set
      !-----------------------------------------------------------------
         ! The U, V components are ASSUMED to come in adjacent pairs.
         ! So we capture what _might_ be the U component and when we
         ! know we have a V component, we can check these to make sure
         ! they are the matching U.  They are SUPPOSED to come U ... V. 
         ! Since we are using 'cycle observationloop' we have to capture
         ! the previous obs before reading a new one ... even for the 
         ! first trip through ...

         U_obs                 = obs(max(obsindex,1))
         U_obs_err_cov         = obs_err_cov(max(obsindex,1))
         U_obs_loc             = obs_loc
         U_flavor              = flavor
         U_prior_mean(1)       = prior_mean(1)
         U_prior_spread(1)     = prior_spread(1)
         U_posterior_mean(1)   = posterior_mean(1)
         U_posterior_spread(1) = posterior_spread(1)

         call get_obs_from_key(seq, keys(obsindex), observation)

         call get_obs_def(observation, obs_def)

         obs_loc   = get_obs_def_location(obs_def)
         obsloc3   = get_location(obs_loc) 
         lon0      = obsloc3(1) + 1.01_r8                            ! 0-360
         lat0      = obsloc3(2) + 1.01_r8 + 90.0_r8                  ! 0-180
         ipressure = 0.01_r8 * obsloc3(3)                       ! mb TJH ...  rounding?

         obs_err_cov(obsindex) = get_obs_def_error_variance(obs_def)
         call get_qc(observation, qc(obsindex:obsindex), 1)

         obs_kind = get_obs_def_kind(obs_def)
         flavor   = get_obs_kind(obs_kind)

         call get_obs_values(observation, obs(obsindex:obsindex), obs_index)

         ! get interpolated values of prior and posterior ensemble mean
         call get_obs_values(observation,       prior_mean,       prior_mean_index)
         call get_obs_values(observation,   posterior_mean,   posterior_mean_index)
         call get_obs_values(observation,     prior_spread,     prior_spread_index)
         call get_obs_values(observation, posterior_spread, posterior_spread_index)
         ! if(obsindex==1)  print*, 'prior_mean_index= ', prior_mean_index, posterior_mean_index
         ! if(obsindex==1)  print*, 'prior_spread_index= ', prior_spread_index, posterior_spread_index

         !--------------------------------------------------------
         ! A Whole bunch of reasons to be rejected
         !--------------------------------------------------------
         if( qc(obsindex) >= 4.0_r8 )     cycle ObservationLoop 
         if (prior_mean(1) == 0.0_r8)     cycle ObservationLoop   ! test for default value
         keeper = CheckObsType(obs_select, obs_err_cov(obsindex), flavor, ipressure)
         if ( .not. keeper )              cycle ObservationLoop

         obs_used = obs_used + 1

         Areas : do iregion =1, Nregions

            !  Set the area of interest

            keeper = InRegion( lon0, lat0, lonlim1(iregion), lonlim2(iregion), &
                                           latlim1(iregion), latlim2(iregion))
            if ( .not. keeper ) cycle Areas

            !-------------------------------------------------------------------
            ! Surface pressure
            !-------------------------------------------------------------------

            if (flavor == KIND_PS ) then

               ratio = GetRatio(obs(obsindex), prior_mean(1), prior_spread(1), &
                        obs_err_cov(obsindex))

               if ( ratio <= rat_cri ) then
                    num_in_level_P(iepoch, iregion) =   num_in_level_P(iepoch, iregion) + 1
                    rms_ges_mean_P(iepoch, iregion) =   rms_ges_mean_P(iepoch, iregion) + &
                                                       (    prior_mean(1)-obs(obsindex))**2
                    rms_anl_mean_P(iepoch, iregion) =   rms_anl_mean_P(iepoch, iregion) + &
                                                       (posterior_mean(1)-obs(obsindex))**2
                  rms_ges_spread_P(iepoch, iregion) = rms_ges_spread_P(iepoch, iregion) + &
                                                      prior_spread(1)**2
                  rms_anl_spread_P(iepoch, iregion) = rms_anl_spread_P(iepoch, iregion) + &
                                                      posterior_spread(1)**2
               endif

               cycle Areas

            endif

            !-------------------------------------------------
            ! Time series statistics of the selected layer
            !-------------------------------------------------
            ! If the observation is within our layer, great -- if
            ! not ... we still need to do vertical stats.

            condition3:  if(ipressure .le. pint(level_index) .and. &
                            ipressure .gt. pint(level_index+1))    then

               select case ( flavor )
                  case ( KIND_V ) !! Wind component

                     ratio = abs( prior_mean(1)-obs(obsindex) ) / &
                             sqrt(prior_spread(1)**2 + obs_err_cov(obsindex))

                     ! We are accumulating BOTH the U and V components in
                     ! the same storage ... must think about normalizing later ...  

                     if( ratio <= rat_cri )  then
                          num_in_level_W(iepoch,iregion) =   num_in_level_W(iepoch,iregion) + 1
                          rms_ges_mean_W(iepoch,iregion) =   rms_ges_mean_W(iepoch,iregion) + &
                                                               (prior_mean(1)-obs(obsindex))**2
                          rms_anl_mean_W(iepoch,iregion) =   rms_anl_mean_W(iepoch,iregion) + &
                                                           (posterior_mean(1)-obs(obsindex))**2
                        rms_ges_spread_W(iepoch,iregion) = rms_ges_spread_W(iepoch,iregion) + &
                                                             prior_spread(1)**2 
                        rms_anl_spread_W(iepoch,iregion) = rms_anl_spread_W(iepoch,iregion) + &
                                                         posterior_spread(1)**2 
                     endif

                  case ( KIND_T ) !! Temperature

                     ratio = GetRatio(obs(obsindex), prior_mean(1), prior_spread(1), &
                              obs_err_cov(obsindex))
                     if(ratio <= rat_cri ) then
                          num_in_level_T(iepoch,iregion) =   num_in_level_T(iepoch,iregion) + 1
                          rms_ges_mean_T(iepoch,iregion) =   rms_ges_mean_T(iepoch,iregion) + &
                                                  (prior_mean(1)-obs(obsindex))**2
                          rms_anl_mean_T(iepoch,iregion) =   rms_anl_mean_T(iepoch,iregion) + &
                                              (posterior_mean(1)-obs(obsindex))**2
                        rms_ges_spread_T(iepoch,iregion) = rms_ges_spread_T(iepoch,iregion) + &
                                                             prior_spread(1)**2
                        rms_anl_spread_T(iepoch,iregion) = rms_anl_spread_T(iepoch,iregion) + &
                                                         posterior_spread(1)**2
                     endif

                  case ( KIND_QV ) !! Moisture Q

                     ratio = GetRatio(obs(obsindex), prior_mean(1), prior_spread(1), &
                              obs_err_cov(obsindex))
                     if(ratio <= rat_cri ) then
                          num_in_level_Q(iepoch,iregion) =   num_in_level_Q(iepoch,iregion) + 1
                          rms_ges_mean_Q(iepoch,iregion) =   rms_ges_mean_Q(iepoch,iregion) + &
                                                  (prior_mean(1)- obs(obsindex))**2
                          rms_anl_mean_Q(iepoch,iregion) =   rms_anl_mean_Q(iepoch,iregion) + &
                                              (posterior_mean(1)- obs(obsindex))**2
                        rms_ges_spread_Q(iepoch,iregion) = rms_ges_spread_Q(iepoch,iregion) + &
                                                             prior_spread(1)**2
                        rms_anl_spread_Q(iepoch,iregion) = rms_anl_spread_Q(iepoch,iregion) + &
                                                         posterior_spread(1)**2
                     endif

               end select

            endif condition3

            !-------------------------------------------------------------------
            ! end of time series statistics
            !-------------------------------------------------------------------
            ! vertical statistical part
            !-------------------------------------------------------------------

            if(iday <= iskip)             cycle Areas
            if(ipressure > pint(1) )      cycle Areas   ! pressure more than 1025
            if(ipressure < pint(nlev+1) ) cycle Areas   ! pressure less than   75

            k0 = GetClosestLevel(ipressure)

            select case (flavor)

               case ( KIND_V ) ! Wind v-component

                  ierr = CheckMate(flavor, U_flavor, obs_loc, U_obs_loc) 
                  if ( ierr /= 0 ) then
                     write(*,*)'Aiiiii ... V with no matching U ...'
                     NbadUV = NbadUV + 1
                     cycle ObservationLoop
                  endif

                  ! Since we are treating wind as a vector quantity, the ratio
                  ! calculation is a bit different.

                  numer = (prior_mean(1)-obs(obsindex))**2 + (U_prior_mean(1)-U_obs)**2
                  denom = prior_spread(1)**2 + U_prior_spread(1)**2 + &
                          obs_err_cov(obsindex) +  U_obs_err_cov 
                  ratio = sqrt( numer / denom )

                  if(ratio <= rat_cri ) then
                     speed_ges2 = sqrt( U_prior_mean(1)**2 + prior_mean(1)**2 )
                     speed_anl2 = sqrt( U_posterior_mean(1)**2 + posterior_mean(1)**2 )
                     speed_obs2 = sqrt( U_obs**2 + obs(obsindex)**2 )

                          num_ver_W(k0,iregion) =      num_ver_W(k0,iregion) + 1
                      rms_ges_ver_W(k0,iregion) =  rms_ges_ver_W(k0,iregion) + &
                                                   (prior_mean(1)-obs(obsindex))**2 + &
                                                   (U_prior_mean(1)- U_obs)**2
                      rms_anl_ver_W(k0,iregion) =  rms_anl_ver_W(k0,iregion) + &
                                                   (posterior_mean(1)-obs(obsindex))**2 + &
                                                   (U_posterior_mean(1)- U_obs)**2
                     bias_ges_ver_W(k0,iregion) = bias_ges_ver_W(k0,iregion) + &
                                                  speed_ges2 - speed_obs2
                     bias_anl_ver_W(k0,iregion) = bias_anl_ver_W(k0,iregion) + &
                                                  speed_anl2 - speed_obs2
                  endif

               case ( KIND_T ) ! Temperature

                  ratio = GetRatio(obs(obsindex), prior_mean(1), prior_spread(1), &
                           obs_err_cov(obsindex))

                  if(ratio <= rat_cri )  then
                          num_ver_T(k0,iregion) =      num_ver_T(k0,iregion) + 1
                      rms_ges_ver_T(k0,iregion) =  rms_ges_ver_T(k0,iregion) + &
                                                  (prior_mean(1)    - obs(obsindex))**2
                      rms_anl_ver_T(k0,iregion) =  rms_anl_ver_T(k0,iregion) + &
                                                  (posterior_mean(1)- obs(obsindex))**2
                     bias_ges_ver_T(k0,iregion) = bias_ges_ver_T(k0,iregion) + &
                                                  prior_mean(1)     - obs(obsindex)
                     bias_anl_ver_T(k0,iregion) = bias_anl_ver_T(k0,iregion) + &
                                                  posterior_mean(1) - obs(obsindex)
                  endif

               case ( KIND_QV ) ! Moisture

                  ratio = GetRatio(obs(obsindex), prior_mean(1), prior_spread(1), &
                           obs_err_cov(obsindex))

                  if(ratio <= rat_cri )  then
                          num_ver_Q(k0,iregion) =      num_ver_Q(k0,iregion) + 1  
                      rms_ges_ver_Q(k0,iregion) =  rms_ges_ver_Q(k0,iregion) + &
                                                  (prior_mean(1)    - obs(obsindex))**2
                      rms_anl_ver_Q(k0,iregion) =  rms_anl_ver_Q(k0,iregion) + &
                                                  (posterior_mean(1)- obs(obsindex))**2
                     bias_ges_ver_Q(k0,iregion) = bias_ges_ver_Q(k0,iregion) + &
                                                 prior_mean(1)     - obs(obsindex)
                     bias_anl_ver_Q(k0,iregion) = bias_anl_ver_Q(k0,iregion) + &
                                                 posterior_mean(1) - obs(obsindex)
                  endif

            end select

            !-------------------------------------------------------------------
            !  end of vertical statistics
            !-------------------------------------------------------------------

         enddo Areas

      !-------------------------------------------------------- 
      end do ObservationLoop
      !-------------------------------------------------------- 

      deallocate(keys, obs,  obs_err_cov, qc)

      do iregion=1, Nregions
         if ( num_in_level_W(iepoch, iregion) .gt. 0) then
              rms_ges_mean_W(iepoch, iregion) = sqrt(   rms_ges_mean_W(iepoch, iregion) / &
                                                        num_in_level_W(iepoch, iregion) )
              rms_anl_mean_W(iepoch, iregion) = sqrt(   rms_anl_mean_W(iepoch, iregion) / &
                                                        num_in_level_W(iepoch, iregion) )
            rms_ges_spread_W(iepoch, iregion) = sqrt( rms_ges_spread_W(iepoch, iregion) / &
                                                        num_in_level_W(iepoch, iregion) )
            rms_anl_spread_W(iepoch, iregion) = sqrt( rms_anl_spread_W(iepoch, iregion) / &
                                                        num_in_level_W(iepoch, iregion) )
         else
              rms_ges_mean_W(iepoch, iregion) = 0.0_r8
              rms_anl_mean_W(iepoch, iregion) = 0.0_r8
            rms_ges_spread_W(iepoch, iregion) = 0.0_r8
            rms_anl_spread_W(iepoch, iregion) = 0.0_r8
         endif

         if ( num_in_level_T(iepoch, iregion) .gt. 0) then
              rms_ges_mean_T(iepoch, iregion) = sqrt(   rms_ges_mean_T(iepoch, iregion) / &
                                                        num_in_level_T(iepoch, iregion) )
              rms_anl_mean_T(iepoch, iregion) = sqrt(   rms_anl_mean_T(iepoch, iregion) / &
                                                        num_in_level_T(iepoch, iregion) )
            rms_ges_spread_T(iepoch, iregion) = sqrt( rms_ges_spread_T(iepoch, iregion) / &
                                                        num_in_level_T(iepoch, iregion) )
            rms_anl_spread_T(iepoch, iregion) = sqrt( rms_anl_spread_T(iepoch, iregion) / &
                                                        num_in_level_T(iepoch, iregion) )
         else
              rms_ges_mean_T(iepoch, iregion) = 0.0_r8
              rms_anl_mean_T(iepoch, iregion) = 0.0_r8
            rms_ges_spread_T(iepoch, iregion) = 0.0_r8
            rms_anl_spread_T(iepoch, iregion) = 0.0_r8
         endif

         if ( num_in_level_Q(iepoch, iregion) .gt. 0) then
              rms_ges_mean_Q(iepoch, iregion) = sqrt(   rms_ges_mean_Q(iepoch, iregion) / &
                                                     num_in_level_Q(iepoch,iregion) )*1000.0_r8
              rms_anl_mean_Q(iepoch, iregion) = sqrt(   rms_anl_mean_Q(iepoch, iregion) / &
                                                     num_in_level_Q(iepoch,iregion) )*1000.0_r8
            rms_ges_spread_Q(iepoch, iregion) = sqrt( rms_ges_spread_Q(iepoch, iregion) / &
                                                     num_in_level_Q(iepoch,iregion) )*1000.0_r8
            rms_anl_spread_Q(iepoch, iregion) = sqrt( rms_anl_spread_Q(iepoch, iregion) / &
                                                     num_in_level_Q(iepoch,iregion) )*1000.0_r8
         else
              rms_ges_mean_Q(iepoch, iregion) = 0.0_r8
              rms_anl_mean_Q(iepoch, iregion) = 0.0_r8
            rms_ges_spread_Q(iepoch, iregion) = 0.0_r8
            rms_anl_spread_Q(iepoch, iregion) = 0.0_r8
         endif

         if ( num_in_level_P(iepoch, iregion) .gt. 0) then
              rms_ges_mean_P(iepoch, iregion) = sqrt(   rms_ges_mean_P(iepoch, iregion) / &
                                                        num_in_level_P(iepoch,iregion) )
              rms_anl_mean_P(iepoch, iregion) = sqrt(   rms_anl_mean_P(iepoch, iregion) / &
                                                        num_in_level_P(iepoch,iregion) )
            rms_ges_spread_P(iepoch, iregion) = sqrt( rms_ges_spread_P(iepoch, iregion) / &
                                                      num_in_level_P(iepoch,iregion) )
            rms_anl_spread_P(iepoch, iregion) = sqrt( rms_anl_spread_P(iepoch, iregion) / &
                                                      num_in_level_P(iepoch,iregion) )
         else
              rms_ges_mean_P(iepoch, iregion) = 0.0_r8
              rms_anl_mean_P(iepoch, iregion) = 0.0_r8
            rms_ges_spread_P(iepoch, iregion) = 0.0_r8
            rms_anl_spread_P(iepoch, iregion) = 0.0_r8
         endif
      enddo

   end do Advancesets

   call destroy_obs_sequence(seq)

end do Dayloop

!-----------------------------------------------------------------------
write(WgesName,'(''Wges_times_'',i4.4,''mb.dat'')') plev(level_index)
write(WanlName,'(''Wanl_times_'',i4.4,''mb.dat'')') plev(level_index)
write(TgesName,'(''Tges_times_'',i4.4,''mb.dat'')') plev(level_index)
write(TanlName,'(''Tanl_times_'',i4.4,''mb.dat'')') plev(level_index)
write(QgesName,'(''Qges_times_'',i4.4,''mb.dat'')') plev(level_index)
write(QanlName,'(''Qanl_times_'',i4.4,''mb.dat'')') plev(level_index)
write(PgesName,'(''Pges_times.dat'')')
write(PanlName,'(''Panl_times.dat'')')

WgesUnit = get_unit()
OPEN(WgesUnit,FILE=trim(adjustl(WgesName)),FORM='formatted')
WanlUnit = get_unit()
OPEN(WanlUnit,FILE=trim(adjustl(WanlName)),FORM='formatted')

TgesUnit = get_unit()
OPEN(TgesUnit,FILE=trim(adjustl(TgesName)),FORM='formatted')
TanlUnit = get_unit()
OPEN(TanlUnit,FILE=trim(adjustl(TanlName)),FORM='formatted')

QgesUnit = get_unit()
OPEN(QgesUnit,FILE=trim(adjustl(QgesName)),FORM='formatted')
QanlUnit = get_unit()
OPEN(QanlUnit,FILE=trim(adjustl(QanlName)),FORM='formatted')

PgesUnit = get_unit()
OPEN(PgesUnit,FILE=trim(adjustl(PgesName)),FORM='formatted')
PanlUnit = get_unit()
OPEN(PanlUnit,FILE=trim(adjustl(PanlName)),FORM='formatted')

do i=1, iepoch
   write(WgesUnit,91) epoch_center(i),rms_ges_mean_W(i,:)
   write(WanlUnit,91) epoch_center(i),rms_anl_mean_W(i,:)

   write(TgesUnit,91) epoch_center(i),rms_ges_mean_T(i,:)
   write(TanlUnit,91) epoch_center(i),rms_anl_mean_T(i,:)

   write(QgesUnit,91) epoch_center(i),rms_ges_mean_Q(i,:)
   write(QanlUnit,91) epoch_center(i),rms_anl_mean_Q(i,:)

   write(PgesUnit,91) epoch_center(i),rms_ges_mean_P(i,:)
   write(PanlUnit,91) epoch_center(i),rms_anl_mean_P(i,:)
enddo

91 format(f6.0, 4(2f7.2, 8x) )

close(WgesUnit)
close(WanlUnit)
close(TgesUnit)
close(TanlUnit)
close(QgesUnit)
close(QanlUnit)
close(PgesUnit)
close(PanlUnit)

deallocate(rms_ges_mean_W, rms_ges_spread_W, &
           rms_anl_mean_W, rms_anl_spread_W, &
           rms_ges_mean_T, rms_ges_spread_T, &
           rms_anl_mean_T, rms_anl_spread_T, &
           rms_ges_mean_Q, rms_ges_spread_Q, &
           rms_anl_mean_Q, rms_anl_spread_Q, &
           rms_ges_mean_P, rms_ges_spread_P, &
           rms_anl_mean_P, rms_anl_spread_P)

deallocate(num_in_level_W, &
           num_in_level_T, &
           num_in_level_Q, &
           num_in_level_P)

deallocate( epoch_center)

!-------------------------------------------------
!  do all day average of the vertical statistics
!-------------------------------------------------
do iregion=1, Nregions
   do ilev=1, nlev
      if(     num_ver_W(ilev,iregion) .ne.0) then
          rms_ges_ver_W(ilev,iregion) = sqrt( rms_ges_ver_W(ilev,iregion) / &
                                                  num_ver_W(ilev,iregion))
          rms_anl_ver_W(ilev,iregion) = sqrt( rms_anl_ver_W(ilev,iregion) / &
                                                  num_ver_W(ilev,iregion))
         bias_ges_ver_W(ilev,iregion) =      bias_ges_ver_W(ilev,iregion) / &
                                                  num_ver_W(ilev,iregion) 
         bias_anl_ver_W(ilev,iregion) =      bias_anl_ver_W(ilev,iregion) / &
                                                  num_ver_W(ilev,iregion)
      endif

      if(     num_ver_T(ilev,iregion) .ne.0) then
          rms_ges_ver_T(ilev,iregion) = sqrt( rms_ges_ver_T(ilev,iregion) / &
                                                  num_ver_T(ilev,iregion))
          rms_anl_ver_T(ilev,iregion) = sqrt( rms_anl_ver_T(ilev,iregion) / &
                                                  num_ver_T(ilev,iregion))
         bias_ges_ver_T(ilev,iregion) =      bias_ges_ver_T(ilev,iregion) / &
                                                  num_ver_T(ilev,iregion)
         bias_anl_ver_T(ilev,iregion) =      bias_anl_ver_T(ilev,iregion) / &
                                                  num_ver_T(ilev,iregion)
      endif

      if(     num_ver_Q(ilev,iregion) .ne.0) then
          rms_ges_ver_Q(ilev,iregion) = sqrt( rms_ges_ver_Q(ilev,iregion) / &
                                                  num_ver_Q(ilev,iregion))*1000.0_r8
          rms_anl_ver_Q(ilev,iregion) = sqrt( rms_anl_ver_Q(ilev,iregion) / &
                                                  num_ver_Q(ilev,iregion))*1000.0_r8
         bias_ges_ver_Q(ilev,iregion) =      bias_ges_ver_Q(ilev,iregion) / &
                                                  num_ver_Q(ilev,iregion)*1000.0_r8
         bias_anl_ver_Q(ilev,iregion) =      bias_anl_ver_Q(ilev,iregion) / &
                                                  num_ver_Q(ilev,iregion)*1000.0_r8
      endif
   enddo
enddo

OPEN(185,FILE='Wges_ver_ave.dat',FORM='FORMATTED')
OPEN(186,FILE='Wanl_ver_ave.dat',FORM='FORMATTED')
OPEN(285,FILE='Tges_ver_ave.dat',FORM='FORMATTED')
OPEN(286,FILE='Tanl_ver_ave.dat',FORM='FORMATTED')
OPEN(385,FILE='Qges_ver_ave.dat',FORM='FORMATTED')
OPEN(386,FILE='Qanl_ver_ave.dat',FORM='FORMATTED')

OPEN(195,FILE='Wges_ver_ave_bias.dat',FORM='FORMATTED')
OPEN(196,FILE='Wanl_ver_ave_bias.dat',FORM='FORMATTED')
OPEN(295,FILE='Tges_ver_ave_bias.dat',FORM='FORMATTED')
OPEN(296,FILE='Tanl_ver_ave_bias.dat',FORM='FORMATTED')
OPEN(395,FILE='Qges_ver_ave_bias.dat',FORM='FORMATTED')
OPEN(396,FILE='Qanl_ver_ave_bias.dat',FORM='FORMATTED')

do ilev = nlev, 1, -1
   write(185, 610) plev(ilev), ( rms_ges_ver_W(ilev,iregion), num_ver_W(ilev,iregion), iregion=1, Nregions)
   write(186, 610) plev(ilev), ( rms_anl_ver_W(ilev,iregion), num_ver_W(ilev,iregion), iregion=1, Nregions) 
   write(285, 610) plev(ilev), ( rms_ges_ver_T(ilev,iregion), num_ver_T(ilev,iregion), iregion=1, Nregions)
   write(286, 610) plev(ilev), ( rms_anl_ver_T(ilev,iregion), num_ver_T(ilev,iregion), iregion=1, Nregions) 
   write(385, 610) plev(ilev), ( rms_ges_ver_Q(ilev,iregion), num_ver_Q(ilev,iregion), iregion=1, Nregions)
   write(386, 610) plev(ilev), ( rms_anl_ver_Q(ilev,iregion), num_ver_Q(ilev,iregion), iregion=1, Nregions) 
   write(195, 610) plev(ilev), (bias_ges_ver_W(ilev,iregion), num_ver_W(ilev,iregion), iregion=1, Nregions)
   write(196, 610) plev(ilev), (bias_anl_ver_W(ilev,iregion), num_ver_W(ilev,iregion), iregion=1, Nregions) 
   write(295, 610) plev(ilev), (bias_ges_ver_T(ilev,iregion), num_ver_T(ilev,iregion), iregion=1, Nregions)
   write(296, 610) plev(ilev), (bias_anl_ver_T(ilev,iregion), num_ver_T(ilev,iregion), iregion=1, Nregions) 
   write(395, 610) plev(ilev), (bias_ges_ver_Q(ilev,iregion), num_ver_Q(ilev,iregion), iregion=1, Nregions)
   write(396, 610) plev(ilev), (bias_anl_ver_Q(ilev,iregion), num_ver_Q(ilev,iregion), iregion=1, Nregions)
enddo
610 format(i5, 4(f8.3, i8) )

close(185)
close(186)
close(285)
close(286)
close(385)
close(386)

close(195)
close(196)
close(295)
close(296)
close(395)
close(396)

call timestamp(source,revision,revdate,'end') ! That closes the log file, too.

contains

   Function GetRatio(obsval, prmean, prspred, errcov) result (ratio)
   real(r8), intent(in) :: obsval, prmean, prspred, errcov
   real(r8)             :: ratio

   real(r8) :: numer, denom

   numer = abs(prmean- obsval) 
   denom = sqrt( prspred**2 + errcov )
   ratio = numer / denom

   end Function GetRatio



   Function CheckMate(flavor1, flavor2, obsloc1, obsloc2) result (ierr)
   integer,             intent(in) :: flavor1, flavor2
   type(location_type), intent(in) :: obsloc1, obsloc2
   integer                         :: ierr

   ierr = -1 ! Assume no match ... till proven otherwise

   ! flavor 1 has to be either U or V, flavor 2 has to be the complement
   if ( ((flavor1 + flavor2) /= (KIND_U + KIND_V)) .and. &
        ((flavor1 == KIND_U) .or. (flavor1 == KIND_V)) ) then
      write(*,*) 'flavors not appropriate ...',flavor1, flavor2
      return
   endif 

   if ( obsloc1 /= obsloc2 ) then
      write(*,*) 'locations no match ...'
      call write_location(6,obsloc1,'FORMATTED')
      call write_location(6,obsloc2,'FORMATTED')
      return
   endif

   ierr = 0

   end Function CheckMate



   Function GetClosestLevel(pressure) result (level_index)
   ! The pressure intervals are ordered  surface == 1
   ! So if we start at the surface and work upwards the surface, ...
   !
   ! We are using nlev and the pressure arrays from global storage
   !
   integer, intent(in) :: pressure
   integer             :: level_index

   integer, dimension(nlev) :: dx
   integer, dimension(nlev), save :: inds = (/ (i,i=1,nlev) /)
   integer :: i

   dx = abs(pressure - plev)    ! whole array
 
   ! level_index = minloc(dx)   why this does not work is beyond me
   level_index = minval( inds , mask=(dx == minval(dx)))

   end Function GetClosestLevel


   Function InRegion( lon, lat, lon1, lon2, lat1, lat2 ) result( keeper )
   real(r8), intent(in) :: lon, lat, lon1, lon2, lat1, lat2
   logical :: keeper

   keeper = .false.

   if( (lon .ge. lon1) .and. (lon .le. lon2) .and. &
       (lat .ge. lat1) .and. (lat .le. lat2) ) keeper = .true.

   end Function InRegion


   Function CheckObsType(obs_select, obs_err_cov, flavor, ipressure) result(keeper)
   ! Since the observation kind does not have platform information, Hui
   ! has determined an ad-hoc set of rules to determine the origin of the
   ! observation. Hence the 'magic' numbers.

   integer,  intent(in) :: obs_select 
   real(r8), intent(in) :: obs_err_cov 
   integer,  intent(in) :: flavor, ipressure
   logical              :: keeper

   keeper = .true. ! Innocent till proven guilty ...

   select case ( flavor )
      case ( KIND_V ) !! Wind component

         if(obs_select == 2) then ! keep RA only and skip ACARS and SATWND data
            if( abs(sqrt( obs_err_cov) - 2.5_r8 ) <= 0.1_r8)               keeper = .false.
            if(     sqrt( obs_err_cov) > 3.3_r8)                           keeper = .false.
            if(     sqrt( obs_err_cov) > 2.5_r8 .and. ipressure .gt. 420 ) keeper = .false.
            if(     sqrt( obs_err_cov) > 1.7_r8 .and. ipressure .gt. 680 ) keeper = .false.
         endif

         if(obs_select == 3) then ! keep only ACARS and SATWND data
            if( sqrt( obs_err_cov) <  2.5_r8 .and. ipressure < 400 ) keeper = .false.
            if( sqrt( obs_err_cov) >  2.5_r8 .and. &
                sqrt( obs_err_cov) <  3.5_r8 .and. ipressure < 400 ) keeper = .false.
            if( sqrt( obs_err_cov) <= 1.7_r8)                        keeper = .false.
         endif

      case ( KIND_T ) !! Temperature

         if(obs_select == 2) then ! temporarily keep RA only 
            if( abs(sqrt(obs_err_cov) - 1.0_r8) <= 0.1_r8) keeper = .false.
         endif

         if(obs_select == 3) then ! temporarily keep only ACARS and SATWND data
            if( abs(sqrt(obs_err_cov) - 1.0_r8) >  0.1_r8) keeper = .false.
         endif

      case default

         ! Moisture (Q) and Surface Pressure (P) are always accepted
         keeper = .true.

   end select
   end Function CheckObsType

end program obs_diag
