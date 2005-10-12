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

!-----------------------------------------------------------------------
! The programs defines a series of epochs (periods of time) and geographic
! regions and accumulates statistics for these epochs and regions.
!
! All 'possible' obs_kinds are treated separately.
!-----------------------------------------------------------------------

use        types_mod, only : r8, digits12
use obs_sequence_mod, only : read_obs_seq, obs_type, obs_sequence_type, get_first_obs, &
                             get_obs_from_key, get_obs_def, get_copy_meta_data, &
                             get_obs_time_range, get_time_range_keys, get_num_obs, &
                             get_next_obs, get_num_times, get_obs_values, init_obs, &
                             assignment(=), get_num_copies, static_init_obs_sequence, &
                             get_qc, destroy_obs_sequence, read_obs_seq_header, & 
                             get_last_obs, destroy_obs
use      obs_def_mod, only : obs_def_type, get_obs_def_error_variance, get_obs_def_time, &
                             get_obs_def_location,  get_obs_kind, get_obs_name
use     obs_kind_mod, only : max_obs_kinds, get_obs_kind_var_type, get_obs_kind_name, &
                             KIND_U_WIND_COMPONENT, KIND_V_WIND_COMPONENT, &
                             KIND_SURFACE_PRESSURE, KIND_SPECIFIC_HUMIDITY
use     location_mod, only : location_type, get_location, set_location_missing, &
                             write_location, operator(/=), vert_is_surface, vert_is_pressure, &
                             vert_is_height
use time_manager_mod, only : time_type, set_date, set_time, get_time, print_time, &
                             set_calendar_type, print_date, &
                             operator(*), operator(+), operator(-), &
                             operator(>), operator(<), &
                             operator(/=), operator(<=)
use    utilities_mod, only : get_unit, open_file, close_file, register_module, &
                             file_exist, error_handler, E_ERR, E_MSG, &
                             initialize_utilities, logfileunit, timestamp, &
                             find_namelist_in_file, check_namelist_read

implicit none

! CVS Generated file description for error handling, do not edit
character(len=128) :: &
source   = "$Source$", &
revision = "$Revision$", &
revdate  = "$Date$"

type(obs_sequence_type) :: seq
type(obs_type)          :: observation, next_obs
type(obs_type)          :: obs1, obsN
type(obs_def_type)      :: obs_def
type(location_type)     :: obs_loc

!---------------------
integer :: obsindex, i, j, iunit, ierr, io
character(len = 129) :: obs_seq_in_file_name

! Storage with fixed size for observation space diagnostics
real(r8), dimension(1) :: prior_mean, posterior_mean, prior_spread, posterior_spread
real(r8) :: pr_mean, po_mean ! same as above, without useless dimension 
real(r8) :: pr_sprd, po_sprd ! same as above, without useless dimension

!-----------------------------------------------------------------------
! We are treating winds as a vector pair, but we are handling the
! observations serially. Consequently, we exploit the fact that
! the U observations are _followed_ by the V observations.

real(r8)            :: U_obs         = 0.0_r8
real(r8)            :: U_obs_err_var = 0.0_r8
type(location_type) :: U_obs_loc
integer             :: U_type        = KIND_V_WIND_COMPONENT ! intentional mismatch
real(r8)            :: U_pr_mean     = 0.0_r8
real(r8)            :: U_pr_sprd     = 0.0_r8
real(r8)            :: U_po_mean     = 0.0_r8
real(r8)            :: U_po_sprd     = 0.0_r8

integer :: obs_index, prior_mean_index, posterior_mean_index
integer :: prior_spread_index, posterior_spread_index
integer :: key_bounds(2), flavor
integer :: num_copies, num_qc, num_obs, max_num_obs, obs_seq_file_id
character(len=129) :: obs_seq_read_format
logical :: pre_I_format

real(r8), dimension(1) :: obs, qc
real(r8) :: obs_err_var

integer,  allocatable :: keys(:)

logical :: out_of_range, is_there_one, is_this_last, keeper

!-----------------------------------------------------------------------
! Namelist with default values
!
character(len = 129) :: obs_sequence_name = "obs_seq.final"
integer :: obs_year       = 2003     ! the first date of the diagnostics
integer :: obs_month      = 1
integer :: obs_day        = 1
integer :: tot_days       = 1        ! total days
integer :: iskip          = 0        ! skip the first 'iskip' days
integer :: plevel         = 500      ! pressure level (hPa)
integer :: hlevel         = 5000     ! height (meters)
integer :: obs_select     = 1        ! obs type selection: 1=all, 2 =RAonly, 3=noRA
integer :: Nregions       = 4
real(r8):: rat_cri        = 3.0_r8   ! QC ratio
real(r8):: qc_threshold   = 4.0_r8   ! maximum NCEP QC factor
real(r8):: bin_separation = 6.0_r8   ! Bins every so often (hours)
real(r8):: bin_width      = 6.0_r8   ! width of the bin (hours)
logical :: print_mismatched_locs = .false.
logical :: verbose = .false.

! TJH - some kind of crazy nomenclature that South Pole = lat 0?

real(r8), dimension(4) :: lonlim1 = (/   0.0_r8,   0.0_r8,   0.0_r8, 235.0_r8 /)
real(r8), dimension(4) :: lonlim2 = (/ 360.0_r8, 360.0_r8, 360.0_r8, 295.0_r8 /)
real(r8), dimension(4) :: latlim1 = (/ 110.0_r8,  10.0_r8,  70.0_r8, 115.0_r8 /)
real(r8), dimension(4) :: latlim2 = (/ 170.0_r8,  70.0_r8, 110.0_r8, 145.0_r8 /)

character(len = 20), dimension(4) :: reg_names = (/ 'Northern Hemisphere ', &
                                                    'Southern Hemisphere ', &
                                                    'Tropics             ', &
                                                    'North America       ' /)

namelist /obsdiag_nml/ obs_sequence_name, obs_year, obs_month, obs_day, &
                       tot_days, iskip, plevel, hlevel, obs_select, Nregions, rat_cri, &
                       qc_threshold, bin_separation, bin_width, &
                       lonlim1, lonlim2, latlim1, latlim2, reg_names, &
                       print_mismatched_locs, verbose

integer  :: iregion, iepoch, iday, ivar, ifile, num_obs_in_epoch
real(r8) :: lon0, lat0, obsloc3(3)
real(r8) :: speed_obs2, speed_ges2, speed_anl2

!-----------------------------------------------------------------------
! There are (at least) two types of vertical coordinates.
! Observations may be on pressure levels, height levels, etc.
! Since we calculate separate statistics for each observation type,
! we really don't care what kind of level they are on. BUT, since we
! are using one multidimensional array to accomodate all levels, each
! array of levels must be the same length. 
!
! The 'levels' array is a tricky beast.
! Since which_vert() is 2 for pressure and 3 for height, we 
! are designing an array to exploit that.
!-----------------------------------------------------------------------

integer, parameter :: nlev = 11
integer  :: levels(nlev,3)            ! set to 'mean' surface pressure (hPa)
integer  :: levels_int(nlev+1,3)      ! set to a bogus value ... never needed
integer  :: ivert
integer  :: which_vert(max_obs_kinds) ! need to know which kind of level for each obs kind

real(r8) :: scale_factor(max_obs_kinds)

integer  :: ilev          ! counter
integer  :: plevel_index  ! index of pressure level closest to input
integer  :: hlevel_index  ! index of height   level closest to input
integer  :: level_index   ! generic level index

data levels(:,1)     / 1013, 1013, 1013, 1013, 1013, 1013, 1013, 1013, 1013,  1013,  1013/ 
data levels(:,2)     / 1000,  925,  850,  700,  500,  400,  300,  250,  200,   150,   100/ 
data levels(:,3)     / 1000, 2000, 3000, 4000, 5000, 6000, 7000, 8000, 9000, 10000, 11000/

data levels_int(:,1) / 1025, 1025, 1025, 1025, 1025, 1025, 1000, 1000, 1000,  1000, 1000,  1000/ 
data levels_int(:,2) / 1025,  950,  900,  800,  600,  450,  350,  275,  225,  175,   125,    75/
data levels_int(:,3) /    0, 1500, 2500, 3500, 4500, 5500, 6500, 7500, 8500, 9500, 10500, 11500/

!-----------------------------------------------------------------------
! Variables used to accumulate the statistics.
! Dimension 1 is temporal, actually - these are time-by-region-by-var
!-----------------------------------------------------------------------

integer,  allocatable, dimension(:,:,:) :: num_in_level, num_ver

! statistics by time, for a particular level:  time-region-variable
real(r8), allocatable, dimension(:,:,:) :: rms_ges_mean   ! prior mean - obs
real(r8), allocatable, dimension(:,:,:) :: rms_ges_spread ! prior spread     
real(r8), allocatable, dimension(:,:,:) :: rms_anl_mean   ! posterior mean - obs
real(r8), allocatable, dimension(:,:,:) :: rms_anl_spread ! posterior spread

! statistics by level, averaged over time :  level-region-variable
real(r8), allocatable, dimension(:,:,:) ::  rms_ges_ver   ! prior mean - obs
real(r8), allocatable, dimension(:,:,:) :: bias_ges_ver   ! prior mean - obs (L1 norm)
real(r8), allocatable, dimension(:,:,:) ::  rms_anl_ver   ! posterior mean - obs
real(r8), allocatable, dimension(:,:,:) :: bias_anl_ver   ! posterior mean - obs (L1 norm)
                                           
type(time_type), allocatable, dimension(:) :: bincenter
real(digits12),  allocatable, dimension(:) :: epoch_center
integer,         allocatable, dimension(:) :: obs_used_in_epoch

!-----------------------------------------------------------------------
! General purpose variables
!-----------------------------------------------------------------------

integer  :: seconds, days
integer  :: obslevel, NBinsPerDay, Nepochs
integer  :: calendar_type
integer  :: gesUnit, anlUnit

integer  :: nsigma(0:100), nsigmaUnit

real(r8) :: ratio, ratioU
real(digits12) :: t1, tN

type(time_type) :: TimeMin, TimeMax    ! of the entire period of interest
type(time_type) :: beg_time, end_time  ! of the particular bin
type(time_type) :: binsep, binwidth, halfbinwidth 
type(time_type) :: seqT1, seqTN        ! first,last time in entire observation sequence
type(time_type) :: obs_time, skip_time

character(len =   6) :: day_num 
character(len = 129) :: gesName, anlName, msgstring, nml_string

!-----------------------------------------------------------------------
! Some variables to keep track of who's rejected why ...
!-----------------------------------------------------------------------

integer                      :: NwrongType = 0   ! namelist discrimination
integer                      :: NbadQC     = 0   ! out-of-range QC values
integer                      :: NbadLevel  = 0   ! out-of-range pressures

integer, allocatable, dimension(:)   :: NbadW      ! V with no U, all days, regions
integer, allocatable, dimension(:)   :: NbadWvert  ! V with no U, select days, all levels
integer, allocatable, dimension(:,:) :: Nrejected  ! all days, each region, one layer

!=======================================================================
! Get the party started
!=======================================================================

call initialize_utilities('obs_diag')
call register_module(source,revision,revdate) 
call static_init_obs_sequence()  ! Initialize the obs sequence module 


write(logfileunit,*)levels
write(logfileunit,*)levels_int


scale_factor = 1.0_r8
which_vert   = -2     ! set to 'undefined' value 

scale_factor(KIND_SURFACE_PRESSURE)  =    0.01_r8
scale_factor(KIND_SPECIFIC_HUMIDITY) = 1000.0_r8

prior_mean(1)       = 0.0_r8
prior_spread(1)     = 0.0_r8
posterior_mean(1)   = 0.0_r8
posterior_spread(1) = 0.0_r8

U_obs_loc = set_location_missing()

calendar_type = 3   ! HARDWIRED VALUE TO MATCH TIME_MANAGER_MOD 'GREGORIAN'
call set_calendar_type(calendar_type)

!----------------------------------------------------------------------
! Read the namelist
!----------------------------------------------------------------------

call find_namelist_in_file("input.nml", "obsdiag_nml", iunit)
read(iunit, nml = obsdiag_nml, iostat = io)
call check_namelist_read(iunit, io, "obsdiag_nml")

! Record the namelist values used for the run ...
call error_handler(E_MSG,'obs_diag','obsdiag_nml values are',' ',' ',' ')
write(logfileunit,nml=obsdiag_nml)
write(    *      ,nml=obsdiag_nml)

!----------------------------------------------------------------------
! Now that we have input, do some checking and setup
!----------------------------------------------------------------------

allocate(rms_ges_ver(nlev, Nregions, max_obs_kinds), &
         rms_anl_ver(nlev, Nregions, max_obs_kinds), &
        bias_ges_ver(nlev, Nregions, max_obs_kinds), &
        bias_anl_ver(nlev, Nregions, max_obs_kinds), &
             num_ver(nlev, Nregions, max_obs_kinds))

 rms_ges_ver = 0.0_r8
 rms_anl_ver = 0.0_r8
bias_ges_ver = 0.0_r8
bias_anl_ver = 0.0_r8
num_ver      = 0

allocate(NbadW(Nregions), NbadWvert(Nregions), Nrejected(Nregions, max_obs_kinds))

NbadW      = 0
NbadWvert  = 0
Nrejected  = 0

plevel_index = GetClosestLevel(plevel, 2)  ! Set desired input level to nearest actual
hlevel_index = GetClosestLevel(hlevel, 3)  ! ditto for heights.

nsigma = 0
nsigmaUnit = get_unit()
OPEN(nsigmaUnit,FILE='nsigma.dat',FORM='FORMATTED')
write(nsigmaUnit,*) '  day   secs    lon    lat   level    obs   guess  ratio    key   kind'

! Determine temporal bin characteristics.
! Nepochs is the total number of time intervals of the period requested.

NBinsPerDay  = nint( 24.0_r8 / bin_separation )
    binsep   = set_time(nint(bin_separation * 3600.0_r8), 0)
    binwidth = set_time(nint(bin_width      * 3600.0_r8), 0) ! full bin width 
halfbinwidth = set_time(nint(bin_width      * 1800.0_r8), 0) ! half bin width 
beg_time     = set_date(obs_year, obs_month, obs_day, 0, 0, 0)
end_time     = set_time(0, iskip)   ! convert days to skip to a time_type object
skip_time    = beg_time + end_time  ! the start time for the accumulation of vertical statistics 
Nepochs      = NBinsPerDay*tot_days

write(*,*)'Requesting ',Nepochs,' assimilation periods.'

allocate(rms_ges_mean(  Nepochs, Nregions, max_obs_kinds), &
         rms_ges_spread(Nepochs, Nregions, max_obs_kinds), &
         rms_anl_mean(  Nepochs, Nregions, max_obs_kinds), &
         rms_anl_spread(Nepochs, Nregions, max_obs_kinds))
allocate(num_in_level(  Nepochs, Nregions, max_obs_kinds))

rms_ges_mean   = 0.0_r8
rms_anl_mean   = 0.0_r8
rms_ges_spread = 0.0_r8
rms_anl_spread = 0.0_r8
num_in_level   = 0

! Define the centers of the temporal bins.

allocate(bincenter(Nepochs), epoch_center(Nepochs))  ! time_type and 'real'
allocate(obs_used_in_epoch(Nepochs))
obs_used_in_epoch = 0

BinLoop : do iepoch = 1,Nepochs
      bincenter(iepoch) = beg_time + (iepoch-1) * binsep
      call get_time(bincenter(iepoch),seconds,days)
      epoch_center(iepoch) = days + seconds/86400.0_r8

      write(msgstring,'(''epoch '',i4,'' center '')')iepoch
      call print_time(bincenter(iepoch),trim(adjustl(msgstring)),logfileunit)
      call print_date(bincenter(iepoch),trim(adjustl(msgstring)),logfileunit)
      write(logfileunit,*)''
enddo BinLoop

TimeMin = bincenter(   1   ) - halfbinwidth   ! minimum time of interest
TimeMax = bincenter(Nepochs) + halfbinwidth   ! maximum time of interest

if (verbose) call print_time(TimeMin,'minimum time of interest',logfileunit)
if (verbose) call print_time(TimeMax,'maximum time of interest',logfileunit)


!-----------------------------------------------------------------------
! We must assume the observation sequence files span an unknown amount
! of time. We must make some sort of assumption about the naming structure
! of these files. Each file name is the same, but they live in sequentially-
! numbered directories. At one point, the first node in the directory name
! referred to 'month', so we will continue to interpret it that way.
! The last part of the directory name will be incremented ad infinitum.
!
! Directory/file names are similar to    01_03/obs_seq.final
!
! We will stop when we encounter an observation with a time
! stamp out of the desired range. 
!-----------------------------------------------------------------------

! Open and prepare the first observation sequence.
! This entails creating the file name, 

!call GetNextObsSequence(obs_sequence_name, obs_month, obs1, obsN, observation, &
!        next_obs, seq, seqT1, seqTN, obs_index, &
!        prior_mean_index, posterior_mean_index, prior_spread_index, posterior_spread_index)
!
!EpochLoop : do iepoch = 1, Nepochs
!
!
!
!enddo EpochLoop


! The strategy at this point is to open WAY too many files and 
! check the observation sequences against ALL of the temporal bins.
! If the sequence is completely past the time period of interest, we stop.
! If the sequence is completely before the time period of interest, we skip.

ObsFileLoop : do ifile=1, tot_days*4
!-----------------------------------------------------------------------

   write(day_num, '(i2.2,''_'',i2.2,''/'')') obs_month, obs_day + (ifile-1) 
   write(obs_seq_in_file_name,*)trim(adjustl(day_num))//trim(adjustl(obs_sequence_name))

   if ( file_exist(trim(adjustl(obs_seq_in_file_name))) ) then
      write(msgstring,*)'opening ', trim(adjustl(obs_seq_in_file_name))
      call error_handler(E_MSG,'obs_diag',msgstring,source,revision,revdate)
   else
      write(msgstring,*)trim(adjustl(obs_seq_in_file_name)),&
                        ' does not exist. Finishing up.'
      call error_handler(E_MSG,'obs_diag',msgstring,source,revision,revdate)
      exit ObsFileLoop
   endif

   ! Read in information about observation sequence so we can allocate
   ! observations. We need info about how many copies, qc values, etc.

   obs_seq_in_file_name = trim(adjustl(obs_seq_in_file_name)) ! Lahey requirement

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

   call read_obs_seq(obs_seq_in_file_name, 0, 0, 0, seq)

   !--------------------------------------------------------------------
   ! Determine the time encompassed in the observation sequence.
   !--------------------------------------------------------------------

   is_there_one = get_first_obs(seq, obs1)
   if ( .not. is_there_one ) then
      call error_handler(E_ERR,'obs_diag','No first observation  in sequence.', &
      source,revision,revdate)
   endif
   call get_obs_def(obs1,   obs_def)
   seqT1 = get_obs_def_time(obs_def)

   is_there_one = get_last_obs(seq, obsN)
   if ( .not. is_there_one ) then
      call error_handler(E_ERR,'obs_diag','No last observation in sequence.', &
      source,revision,revdate)
   endif
   call get_obs_def(obsN,   obs_def)
   seqTN = get_obs_def_time(obs_def)

   call print_time(seqT1,'First observation time',logfileunit)
   call print_time(seqTN,'Last  observation time',logfileunit)

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
      cycle ObsFileLoop
   else
      if (verbose) write(*,*)'seqTN > TimeMin ... using ',trim(adjustl(obs_seq_in_file_name))
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
      exit ObsFileLoop
   else
      if (verbose) write(*,*)'seqT1 < TimeMax ... using ',trim(adjustl(obs_seq_in_file_name))
   endif

   !--------------------------------------------------------------------
   ! Find the index of obs, ensemble mean, spread ... etc.
   !--------------------------------------------------------------------

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

   !--------------------------------------------------------------------
   ! Make sure we find an index for each of them.
   !--------------------------------------------------------------------

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

   !--------------------------------------------------------------------
   ! Echo what we found.
   !--------------------------------------------------------------------

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

   !====================================================================
   EpochLoop : do iepoch = 1, Nepochs
   !====================================================================

      beg_time = bincenter(iepoch) - halfbinwidth + set_time(1, 0)
      end_time = bincenter(iepoch) + halfbinwidth

      if ( verbose ) then
         call print_time(         beg_time,'epoch  start ',logfileunit)
         call print_time(bincenter(iepoch),'epoch center ',logfileunit)
         call print_time(         end_time,'epoch    end ',logfileunit)
         call print_time(         beg_time,'epoch  start ')
         call print_time(bincenter(iepoch),'epoch center ')
         call print_time(         end_time,'epoch    end ')
      endif

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

      allocate(keys(num_obs_in_epoch))

      call get_time_range_keys(seq, key_bounds, num_obs_in_epoch, keys)

      !-----------------------------------------------------------------
      ObservationLoop : do obsindex = 1, num_obs_in_epoch
      !-----------------------------------------------------------------

         call get_obs_from_key(seq, keys(obsindex), observation)
         call get_obs_def(observation, obs_def)
         obs_time    = get_obs_def_time(obs_def)
         flavor      = get_obs_kind(obs_def)
         obs_err_var = get_obs_def_error_variance(obs_def) * &
                          scale_factor(get_obs_kind_var_type(flavor)) * &
                          scale_factor(get_obs_kind_var_type(flavor))

         obs_loc = get_obs_def_location(obs_def)
         obsloc3 = get_location(obs_loc)
         lon0    = obsloc3(1) + 1.01_r8             ! 0-360
         lat0    = obsloc3(2) + 1.01_r8 + 90.0_r8   ! 0-180

         if(vert_is_surface(obs_loc)) then    ! use closest height equivalent 
            level_index        = 1            ! use 'surface' level
            ivert              = 1
            which_vert(flavor) = 1
            obslevel           = 1
         elseif(vert_is_pressure(obs_loc)) then
            level_index        = plevel_index    ! compare to desired pressure level
            ivert              = nint(0.01_r8 * obsloc3(3))  ! cvrt to hPa
            which_vert(flavor) = 2               ! use pressure column of 'levels'
            obslevel           = GetClosestLevel(ivert, which_vert(flavor))
         elseif(vert_is_height(obs_loc)) then
            level_index        = hlevel_index    ! compare to desired height level
            ivert              = nint(obsloc3(3)) 
            which_vert(flavor) = 3               ! use height column of 'levels'
            obslevel           = GetClosestLevel(ivert, which_vert(flavor))
         else
            call error_handler(E_ERR,'obs_diag','Vertical coordinate not recognized', &
                 source,revision,revdate)
         endif

         call get_qc(observation, qc, 1)
         call get_obs_values(observation, obs, obs_index)
         obs(1) = obs(1)*scale_factor(get_obs_kind_var_type(flavor))

         ! get interpolated values of prior and posterior ensemble mean

         call get_obs_values(observation,       prior_mean,       prior_mean_index)
         call get_obs_values(observation,   posterior_mean,   posterior_mean_index)
         call get_obs_values(observation,     prior_spread,     prior_spread_index)
         call get_obs_values(observation, posterior_spread, posterior_spread_index)

         pr_mean = prior_mean(1)      *scale_factor(get_obs_kind_var_type(flavor))
         po_mean = posterior_mean(1)  *scale_factor(get_obs_kind_var_type(flavor))
         pr_sprd = prior_spread(1)    *scale_factor(get_obs_kind_var_type(flavor))
         po_sprd = posterior_spread(1)*scale_factor(get_obs_kind_var_type(flavor))

         !--------------------------------------------------------------
         ! (DEBUG) Summary of observation knowledge at this point
         !--------------------------------------------------------------

         ! write(*,*)'observation # ',obsindex
         ! write(*,*)'obs_flavor ',flavor
         ! write(*,*)'obs_err_var ',obs_err_var
         ! write(*,*)'lon0/lat0 ',lon0,lat0
         ! write(*,*)'ivert,which_vert,closestlevel ',ivert,which_vert(flavor),obslevel
         ! write(*,*)'qc ',qc
         ! write(*,*)'obs(1) ',obs(1)
         ! write(*,*)'pr_mean,po_mean ',pr_mean,po_mean
         ! write(*,*)'pr_sprd,po_sprd ',pr_sprd,po_sprd

         !--------------------------------------------------------------
         ! A Whole bunch of reasons to be rejected
         !--------------------------------------------------------------

         keeper = CheckObsType(obs_select, flavor)
         if ( .not. keeper ) then
            write(*,*)'obs ',obsindex,' rejected by CheckObsType ',flavor
            NwrongType = NwrongType + 1
            cycle ObservationLoop
         endif

         if( qc(1) >= qc_threshold ) then
         !  write(*,*)'obs ',obsindex,' rejected by qc ',qc(1)
            NbadQC = NbadQC + 1
            cycle ObservationLoop
         endif

         if ( obslevel < 1 .or. obslevel > nlev )   then
         !  write(*,*)'obs ',obsindex,' rejected. Uninteresting level ',obslevel
         !  write(*,*)'obs ',obsindex,' rejected. ivert was ',ivert,&
         !                            ' for which_vert ',which_vert(flavor)
            NbadLevel = NbadLevel + 1
            cycle ObservationLoop
         endif

         ratio = GetRatio(obs(1), pr_mean, pr_sprd, obs_err_var)
         nsigma(int(ratio)) = nsigma(int(ratio)) + 1
  !      if(ratio > 10.0_r8) then
  !         call get_time(obs_time,seconds,days)
  !         write(nsigmaUnit,FMT='(i7,1x,i5,1x,2f7.2,i6,2f8.2,f7.1,2i7)') &
  !              days, seconds, lon0,lat0,ivert, &
  !              obs(1), pr_mean,ratio,obsindex,flavor
  !      endif

         obs_used_in_epoch(iepoch) = obs_used_in_epoch(iepoch) + 1
  !      write(*,*)'incrementing obs in this epoch to ',obs_used_in_epoch(iepoch)

         !--------------------------------------------------------------
         ! If it is a U wind component, all we need to do is save it.
         ! It will be matched up with the subsequent V component.
         ! At some point we have to remove the dependency that the 
         ! U component MUST preceed the V component.
         !--------------------------------------------------------------

         if ( get_obs_kind_var_type(flavor) == KIND_U_WIND_COMPONENT ) then

           !if (verbose) then
           !   write(logfileunit,*)'obs ',obsindex,' is a U wind'
           !   write(     *     ,*)'obs ',obsindex,' is a U wind'
           !endif

            U_obs         = obs(1)
            U_obs_err_var = obs_err_var
            U_obs_loc     = obs_loc
            U_type        = KIND_U_WIND_COMPONENT
            U_pr_mean     = pr_mean
            U_pr_sprd     = pr_sprd
            U_po_mean     = po_mean
            U_po_sprd     = po_sprd

            cycle ObservationLoop

         endif

         !--------------------------------------------------------------
         ! We have Nregions of interest
         !--------------------------------------------------------------

         Areas : do iregion =1, Nregions

            keeper = InRegion( lon0, lat0, lonlim1(iregion), lonlim2(iregion), &
                                           latlim1(iregion), latlim2(iregion))
            if ( .not. keeper ) then
            !  if (verbose) then
            !     write(logfileunit,*)'lon/lat ',lon0,lat0,' not in region ',iregion
            !     write(     *     ,*)'lon/lat ',lon0,lat0,' not in region ',iregion
            !  endif
               cycle Areas
            endif

            !-----------------------------------------------------------
            ! Time series statistics of the selected layer and surface observations.
            !-----------------------------------------------------------
            ! If the observation is within our layer, great -- if
            ! not ... we still need to do vertical stats.

            DesiredLevel: if ( obslevel == level_index ) then

               if ( get_obs_kind_var_type(flavor) == KIND_V_WIND_COMPONENT ) then
                  ! The big assumption is that the U wind component has
                  ! immediately preceeded the V component and has been saved.
                  ! We check for compatibility and proceed.  

                  ierr = CheckMate(KIND_V_WIND_COMPONENT, U_type, obs_loc, U_obs_loc) 
                  if ( ierr /= 0 ) then
                     !  write(*,*)'time series ... V with no matching U ...'
                     NbadW(iregion) = NbadW(iregion) + 1
                     cycle ObservationLoop
                  endif

                  ! since we don't have the necessary covariance between U,V
                  ! we will reject if either univariate z score is bad 

                  ratioU = GetRatio(U_obs, U_pr_mean, U_pr_sprd, U_obs_err_var)

                  if( (ratio <= rat_cri) .and. (ratioU <= rat_cri) )  then
                     num_in_level(iepoch,iregion,flavor) = &
                     num_in_level(iepoch,iregion,flavor) + 1

                     rms_ges_mean(iepoch,iregion,flavor) = &
                     rms_ges_mean(iepoch,iregion,flavor) + &
                     (pr_mean - obs(1))**2 + (U_pr_mean - U_obs)**2

                     rms_anl_mean(iepoch,iregion,flavor) = &
                     rms_anl_mean(iepoch,iregion,flavor) + &
                     (po_mean - obs(1))**2 + (U_po_mean - U_obs)**2

                     ! these are wrong, but we don't have the off-diagonal
                     ! terms of the covariance matrix ... 
                     rms_ges_spread(iepoch,iregion,flavor) = &
                     rms_ges_spread(iepoch,iregion,flavor) + pr_sprd**2 + U_pr_sprd**2

                     rms_anl_spread(iepoch,iregion,flavor) = &
                     rms_anl_spread(iepoch,iregion,flavor) + po_sprd**2 + U_po_sprd**2
                  else
                     Nrejected(iregion,flavor) = Nrejected(iregion,flavor) + 1
                  endif

               else

                  if(ratio <= rat_cri ) then
                     num_in_level(iepoch,iregion,flavor) = &
                     num_in_level(iepoch,iregion,flavor) + 1

                     rms_ges_mean(iepoch,iregion,flavor) = &
                     rms_ges_mean(iepoch,iregion,flavor) + (pr_mean-obs(1))**2

                     rms_anl_mean(iepoch,iregion,flavor) = &
                     rms_anl_mean(iepoch,iregion,flavor) + (po_mean-obs(1))**2

                     rms_ges_spread(iepoch,iregion,flavor) = &
                     rms_ges_spread(iepoch,iregion,flavor) + pr_sprd**2

                     rms_anl_spread(iepoch,iregion,flavor) = &
                     rms_anl_spread(iepoch,iregion,flavor) + po_sprd**2
                  else
                     Nrejected(iregion,flavor) = Nrejected(iregion,flavor) + 1
                  endif

               endif

            endif DesiredLevel

            if (vert_is_surface(obs_loc)) cycle Areas

            !-----------------------------------------------------------
            ! end of time series statistics
            !-----------------------------------------------------------
            ! vertical statistical part
            !-----------------------------------------------------------

            if ( obs_time <= skip_time) cycle Areas

            if ( get_obs_kind_var_type(flavor) == KIND_V_WIND_COMPONENT ) then

               ierr = CheckMate(KIND_V_WIND_COMPONENT, U_type, obs_loc, U_obs_loc)

               if ( ierr /= 0 ) then

                  !  write(*,*)'vertical ... V with no matching U ...'

                  NbadWvert(iregion) = NbadWvert(iregion) + 1
                  cycle ObservationLoop
               endif

               ! Since we are treating wind as a vector quantity, the ratio
               ! calculation is a bit different.

               ratioU = GetRatio(U_obs, U_pr_mean, U_pr_sprd, U_obs_err_var)

               if( (ratio <= rat_cri) .and. (ratioU <= rat_cri) )  then
                  speed_ges2 = sqrt( U_pr_mean**2 + pr_mean**2 )
                  speed_anl2 = sqrt( U_po_mean**2 + po_mean**2 )
                  speed_obs2 = sqrt( U_obs**2 + obs(1)**2 )

                  num_ver(obslevel,iregion,flavor) =      num_ver(obslevel,iregion,flavor) + 1

                  rms_ges_ver(obslevel,iregion,flavor) =  rms_ges_ver(obslevel,iregion,flavor) + &
                       (pr_mean-obs(1))**2 + (U_pr_mean- U_obs)**2

                  rms_anl_ver(obslevel,iregion,flavor) =  rms_anl_ver(obslevel,iregion,flavor) + &
                       (po_mean-obs(1))**2 + (U_po_mean- U_obs)**2

                  bias_ges_ver(obslevel,iregion,flavor) = bias_ges_ver(obslevel,iregion,flavor) + &
                       speed_ges2 - speed_obs2

                  bias_anl_ver(obslevel,iregion,flavor) = bias_anl_ver(obslevel,iregion,flavor) + &
                       speed_anl2 - speed_obs2
               endif

            else

               if(ratio <= rat_cri )  then
                  num_ver(obslevel,iregion,flavor) =      num_ver(obslevel,iregion,flavor) + 1
                  rms_ges_ver(obslevel,iregion,flavor) =  rms_ges_ver(obslevel,iregion,flavor) + &
                       (pr_mean - obs(1))**2
                  rms_anl_ver(obslevel,iregion,flavor) =  rms_anl_ver(obslevel,iregion,flavor) + &
                       (po_mean - obs(1))**2
                  bias_ges_ver(obslevel,iregion,flavor) = bias_ges_ver(obslevel,iregion,flavor) + &
                       pr_mean - obs(1)
                  bias_anl_ver(obslevel,iregion,flavor) = bias_anl_ver(obslevel,iregion,flavor) + &
                       po_mean - obs(1)
               endif

            endif

            !-----------------------------------------------------------
            !  end of vertical statistics
            !-----------------------------------------------------------

         enddo Areas

      !-----------------------------------------------------------------
      enddo ObservationLoop
      !-----------------------------------------------------------------

      deallocate(keys)

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

enddo ObsFileLoop

!-----------------------------------------------------------------------
! We have read all possible files, and stuffed the observations into the
! appropriate bins. Time to normalize and finish up.
!-----------------------------------------------------------------------

do iepoch = 1, Nepochs
   write(msgstring,'(''num obs used in epoch '',i3,'' = '',i8)') iepoch, obs_used_in_epoch(iepoch)
   call error_handler(E_MSG,'obs_diag',msgstring,source,revision,revdate)
enddo

if (verbose) then
   write(logfileunit,*)'Normalizing time-region-variable quantities for desired level.'
   write(     *     ,*)'Normalizing time-region-variable quantities for desired level.'
endif

OurLevel : do ivar=1,max_obs_kinds

   do iregion=1, Nregions
      do iepoch = 1, Nepochs
         if ( num_in_level(iepoch, iregion, ivar) == 0) then
              rms_ges_mean(iepoch, iregion, ivar) = -99.0_r8
              rms_anl_mean(iepoch, iregion, ivar) = -99.0_r8
            rms_ges_spread(iepoch, iregion, ivar) = -99.0_r8
            rms_anl_spread(iepoch, iregion, ivar) = -99.0_r8
         else

         !  write(logfileunit,*)'e,r,v #',iepoch,iregion,ivar,num_in_level(iepoch,iregion,ivar)
         !  write(     *     ,*)'e,r,v #',iepoch,iregion,ivar,num_in_level(iepoch,iregion,ivar)

              rms_ges_mean(iepoch, iregion, ivar) =   sqrt(rms_ges_mean(iepoch, iregion, ivar) / &
                                                           num_in_level(iepoch, iregion, ivar) )
              rms_anl_mean(iepoch, iregion, ivar) =   sqrt(rms_anl_mean(iepoch, iregion, ivar) / &
                                                           num_in_level(iepoch, iregion, ivar) )
            rms_ges_spread(iepoch, iregion, ivar) = sqrt(rms_ges_spread(iepoch, iregion, ivar) / &
                                                           num_in_level(iepoch, iregion, ivar) )
            rms_anl_spread(iepoch, iregion, ivar) = sqrt(rms_anl_spread(iepoch, iregion, ivar) / &
                                                           num_in_level(iepoch, iregion, ivar) )
         endif
      enddo
   enddo

   !--------------------------------------------------------------------
   ! Create data files for all observation kinds we have used
   !--------------------------------------------------------------------

   if( all (num_in_level(:, :, ivar) == 0) ) then
   !  write(logfileunit,*)'No valid data for variable number ',ivar
   !  write(     *     ,*)'No valid data for variable number ',ivar
      cycle OurLevel
   endif

   write(gesName,'(a,''_ges_times.dat'')') trim(adjustl(get_obs_name(ivar)))
   write(anlName,'(a,''_anl_times.dat'')') trim(adjustl(get_obs_name(ivar)))
   gesUnit = open_file(trim(adjustl(gesName)),form='formatted',action='rewind')
   anlUnit = open_file(trim(adjustl(anlName)),form='formatted',action='rewind')

   if (verbose) then
      write(logfileunit,*)'Creating '//trim(adjustl(anlName))
      write(     *     ,*)'Creating '//trim(adjustl(gesName))
   endif

   do i=1, Nepochs
      if( any (num_in_level(i, :, ivar) /= 0) ) then
         call get_time(bincenter(i),seconds,days)
         write(gesUnit,91) days, seconds, &
              (rms_ges_mean(i,j,ivar),rms_ges_spread(i,j,ivar),num_in_level(i,j,ivar),j=1,Nregions)
         write(anlUnit,91) days, seconds, &
              (rms_anl_mean(i,j,ivar),rms_anl_spread(i,j,ivar),num_in_level(i,j,ivar),j=1,Nregions)
      endif
   enddo
   close(gesUnit)
   close(anlUnit)

enddo OurLevel

91 format(i7,1x,i5,4(1x,2f7.2,1x,i8))


close(nsigmaUnit)
!do i=0,100
!   if(nsigma(i) /= 0) print*,i,nsigma(i)
!enddo

deallocate(rms_ges_mean, rms_ges_spread, &
           rms_anl_mean, rms_anl_spread)

deallocate(num_in_level)

!-----------------------------------------------------------------------
! temporal average of the vertical statistics
!-----------------------------------------------------------------------
if (verbose) then
   write(logfileunit,*)'Normalize quantities for all levels.'
   write(     *     ,*)'Normalize quantities for all levels.'
endif

AllLevels : do ivar=1,max_obs_kinds

   do iregion=1, Nregions
      do ilev=1, nlev
         if(     num_ver(ilev,iregion,ivar) == 0) then
             rms_ges_ver(ilev,iregion,ivar) = -99.0_r8
             rms_anl_ver(ilev,iregion,ivar) = -99.0_r8
            bias_ges_ver(ilev,iregion,ivar) = -99.0_r8
            bias_anl_ver(ilev,iregion,ivar) = -99.0_r8
         else

         !  write(logfileunit,*)'l,r,v #',ilev,iregion,ivar,num_ver(ilev,iregion,ivar)
         !  write(     *     ,*)'l,r,v #',ilev,iregion,ivar,num_ver(ilev,iregion,ivar)

             rms_ges_ver(ilev,iregion,ivar) = sqrt(  rms_ges_ver(ilev,iregion,ivar) / &
                                                         num_ver(ilev,iregion,ivar))
             rms_anl_ver(ilev,iregion,ivar) = sqrt(  rms_anl_ver(ilev,iregion,ivar) / &
                                                         num_ver(ilev,iregion,ivar))
            bias_ges_ver(ilev,iregion,ivar) =       bias_ges_ver(ilev,iregion,ivar) / &
                                                         num_ver(ilev,iregion,ivar)
            bias_anl_ver(ilev,iregion,ivar) =       bias_anl_ver(ilev,iregion,ivar) / &
                                                         num_ver(ilev,iregion,ivar)
         endif
      enddo
   enddo

   !--------------------------------------------------------------------
   ! Create data files for all observation kinds we have used
   !--------------------------------------------------------------------

   if( all (num_ver(:,:,ivar) == 0) ) cycle AllLevels

   write(gesName,'(a,''_ges_ver_ave.dat'')') trim(adjustl(get_obs_name(ivar)))
   write(anlName,'(a,''_anl_ver_ave.dat'')') trim(adjustl(get_obs_name(ivar)))
   gesUnit = open_file(trim(adjustl(gesName)),form='formatted',action='rewind')
   anlUnit = open_file(trim(adjustl(anlName)),form='formatted',action='rewind')

   if (verbose) then
      write(logfileunit,*)'Creating '//trim(adjustl(anlName))
      write(     *     ,*)'Creating '//trim(adjustl(gesName))
   endif

   do ilev = nlev, 1, -1
      write(gesUnit, 610) levels(ilev,which_vert(ivar)), &
           (rms_ges_ver(ilev,iregion,ivar), num_ver(ilev,iregion,ivar), iregion=1, Nregions)
      write(anlUnit, 610) levels(ilev,which_vert(ivar)), &
           (rms_anl_ver(ilev,iregion,ivar), num_ver(ilev,iregion,ivar), iregion=1, Nregions) 
   enddo
   close(gesUnit)
   close(anlUnit)

   write(gesName,'(a,''_ges_ver_ave_bias.dat'')') trim(adjustl(get_obs_name(ivar)))
   write(anlName,'(a,''_anl_ver_ave_bias.dat'')') trim(adjustl(get_obs_name(ivar)))
   gesUnit = open_file(trim(adjustl(gesName)),form='formatted',action='rewind')
   anlUnit = open_file(trim(adjustl(anlName)),form='formatted',action='rewind')

   if (verbose) then
      write(logfileunit,*)'Creating '//trim(adjustl(anlName))
      write(     *     ,*)'Creating '//trim(adjustl(gesName))
   endif

   do ilev = nlev, 1, -1
      write(gesUnit, 610) levels(ilev,which_vert(ivar)), &
           (bias_ges_ver(ilev,iregion,ivar), num_ver(ilev,iregion,ivar), iregion=1, Nregions)
      write(anlUnit, 610) levels(ilev,which_vert(ivar)), &
           (bias_anl_ver(ilev,iregion,ivar), num_ver(ilev,iregion,ivar), iregion=1, Nregions) 
   enddo
   close(gesUnit)
   close(anlUnit)

enddo AllLevels

610 format(i5, 4(f8.3, i8) )

write(*,*) ''
write(*,*) '# observations used  : ',sum(obs_used_in_epoch)
write(*,*) 'Rejected Observations summary.'
write(*,*) '# NwrongType         : ',NwrongType
write(*,*) '# NbadQC             : ',NbadQC
write(*,*) '# NbadLevel          : ',NbadLevel
write(*,'(a)')'Table of observations rejected by region for specified level'
write(*,'(5a)')'                                       ',reg_names(1:Nregions)
do ivar=1,max_obs_kinds
   write(*,'(3a,4i8)') ' # ',get_obs_name(ivar),'         : ',Nrejected(:,ivar)
enddo
write(*,'(a,4i8)') ' # bad winds per region for specified level : ',NbadW
write(*,'(a,4i8)') ' # bad winds per region for all levels      : ',NbadWvert
write(*,*)''

write(logfileunit,*) ''
write(logfileunit,*) '# observations used  : ',sum(obs_used_in_epoch)
write(logfileunit,*) 'Rejected Observations summary.'
write(logfileunit,*) '# NwrongType         : ',NwrongType
write(logfileunit,*) '# NbadQC             : ',NbadQC
write(logfileunit,*) '# NbadLevel          : ',NbadLevel
write(logfileunit,'(a)')'Table of observations rejected by region for specified level'
write(logfileunit,'(5a)')'                                       ',reg_names(1:Nregions)
do ivar=1,max_obs_kinds
   write(logfileunit,'(3a,4i8)') ' # ',get_obs_name(ivar),'         : ',Nrejected(:,ivar)
enddo
write(logfileunit,'(a,4i8)') ' # bad winds per region for specified level : ',NbadW
write(logfileunit,'(a,4i8)') ' # bad winds per region for all levels      : ',NbadWvert

!-----------------------------------------------------------------------
! Echo attributes to a file to facilitate plotting.
! This file is also a matlab function ... the variables are
! loaded by just typing the name at the matlab prompt.
!-----------------------------------------------------------------------

iunit = open_file('ObsDiagAtts.m',form='formatted',action='rewind')
write(iunit,'(''obs_year       = '',i6,'';'')')obs_year
write(iunit,'(''obs_month      = '',i6,'';'')')obs_month
write(iunit,'(''obs_day        = '',i6,'';'')')obs_day
write(iunit,'(''tot_days       = '',i6,'';'')')tot_days
write(iunit,'(''iskip          = '',i6,'';'')')iskip
write(iunit,'(''plevel         = '',i6,'';'')')plevel
write(iunit,'(''psurface       = '',i6,'';'')')levels_int(1,2)
write(iunit,'(''ptop           = '',i6,'';'')')levels_int(nlev+1,2)
write(iunit,'(''hlevel         = '',i6,'';'')')hlevel
write(iunit,'(''obs_select     = '',i6,'';'')')obs_select
write(iunit,'(''rat_cri        = '',f9.2,'';'')')rat_cri
write(iunit,'(''qc_threshold   = '',f9.2,'';'')')qc_threshold
write(iunit,'(''bin_width      = '',f9.2,'';'')')bin_width
write(iunit,'(''bin_separation = '',f9.2,'';'')')bin_separation 
write(iunit,'(''t1             = '',f20.6,'';'')')epoch_center(1)
write(iunit,'(''tN             = '',f20.6,'';'')')epoch_center(Nepochs)
write(iunit,'(''plev    = ['',11(1x,i5),''];'')')levels(:,2)
write(iunit,'(''lonlim1 = ['',4(1x,f9.2),''];'')')lonlim1
write(iunit,'(''lonlim2 = ['',4(1x,f9.2),''];'')')lonlim2
write(iunit,'(''latlim1 = ['',4(1x,f9.2),''];'')')latlim1
write(iunit,'(''latlim2 = ['',4(1x,f9.2),''];'')')latlim2
close(iunit)


deallocate(epoch_center, bincenter)
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
   if ( .not.((flavor1 == KIND_U_WIND_COMPONENT .and. flavor2 == KIND_V_WIND_COMPONENT) .or. &
              (flavor2 == KIND_U_WIND_COMPONENT .and. flavor1 == KIND_V_WIND_COMPONENT)) ) then
      write(*,*) 'flavors not complementary ...',flavor1, flavor2
      return
   endif 

   if ( obsloc1 /= obsloc2 ) then
      if ( print_mismatched_locs ) then
         write(*,*) 'locations do not match ...'
         call write_location(6,obsloc1,'FORMATTED')
         call write_location(6,obsloc2,'FORMATTED')
      endif
      return
   endif

   ierr = 0

   end Function CheckMate



   Function GetClosestLevel(level_in, levind) result (level_index_out)
   ! The levels, intervals are ordered  surface == 1
   !
   ! nlev, levels and levels_int are scoped in this unit.
   !
   integer, intent(in) :: level_in         ! target level
   integer, intent(in) :: levind
   integer             :: level_index_out, a(1)

   integer, dimension(nlev) :: dx

   if (levind == 3) then   ! we have heights levels

      if (level_in < levels_int(1,levind) ) then             ! below surface
         level_index_out = -1
      else if ( level_in > levels_int(nlev+1,levind) ) then  ! outer space
         level_index_out = 100 + nlev
      else
         dx = abs(level_in - levels(:,levind))               ! whole array
         a  = minloc(dx)
         level_index_out = a(1)
      endif

      write(*,*)'level_in/column ',level_in,levind,' results in ',level_index_out

   else  ! we have pressure levels (or surface obs ... also in pressure units)

      if (level_in > levels_int(1,levind) )             then ! greater than surface
         level_index_out = -1
      else if ( level_in <= levels_int(nlev+1,levind) ) then ! outer space
         level_index_out = 100 + nlev
      else
         dx = abs(level_in - levels(:,levind))               ! whole array
         a  = minloc(dx)
         level_index_out = a(1)
      endif

   endif


   end Function GetClosestLevel



   Function InRegion( lon, lat, lon1, lon2, lat1, lat2 ) result( keeper )
   real(r8), intent(in) :: lon, lat, lon1, lon2, lat1, lat2
   logical :: keeper

   keeper = .false.

   if( (lon .ge. lon1) .and. (lon .le. lon2) .and. &
       (lat .ge. lat1) .and. (lat .le. lat2) ) keeper = .true.

   end Function InRegion


   Function CheckObsType(obs_select, flavor) result(keeper)

   integer,  intent(in) :: obs_select 
   integer,  intent(in) :: flavor
   logical              :: keeper

   character(len = 10) :: platform

   keeper = .true. ! Innocent till proven guilty ...

   platform = get_obs_kind_name(flavor)
   if( (obs_select == 2) .and. (platform /= 'RADIOSONDE') ) keeper = .false.
   if( (obs_select == 3) .and. (platform == 'RADIOSONDE') ) keeper = .false.

   end Function CheckObsType

end program obs_diag
