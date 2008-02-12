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
                             KIND_U_WIND_COMPONENT, KIND_V_WIND_COMPONENT, &
                             KIND_SURFACE_PRESSURE, KIND_SPECIFIC_HUMIDITY
use     location_mod, only : location_type, get_location, set_location_missing, &
                             write_location, operator(/=),     &
                             vert_is_undef,    VERTISUNDEF,    &
                             vert_is_surface,  VERTISSURFACE,  &
                             vert_is_level,    VERTISLEVEL,    &
                             vert_is_pressure, VERTISPRESSURE, &
                             vert_is_height,   VERTISHEIGHT
use time_manager_mod, only : time_type, set_date, set_time, get_time, print_time, &
                             set_calendar_type, print_date, GREGORIAN, &
                             operator(*), operator(+), operator(-), &
                             operator(>), operator(<), operator(/), &
                             operator(/=), operator(<=)
use    utilities_mod, only : get_unit, open_file, close_file, register_module, &
                             file_exist, error_handler, E_ERR, E_MSG, &
                             initialize_utilities, logfileunit, timestamp, &
                             find_namelist_in_file, check_namelist_read, nc_check
use         sort_mod, only : sort

use typeSizes
use netcdf

implicit none

! version controlled file description for error handling, do not edit
character(len=128), parameter :: &
   source   = "$URL$", &
   revision = "$Revision$", &
   revdate  = "$Date$"

!---------------------------------------------------------------------
!---------------------------------------------------------------------

integer, parameter :: MaxLevels  = 50
integer, parameter :: MaxRegions = 50
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

! Storage with fixed size for observation space diagnostics
real(r8), dimension(1) :: prior_mean, posterior_mean, prior_spread, posterior_spread
real(r8) :: pr_mean, po_mean ! same as above, without useless dimension 
real(r8) :: pr_sprd, po_sprd ! same as above, without useless dimension

! We are treating winds as a vector pair, but we are handling the
! observations serially. Consequently, we exploit the fact that
! the U observations are _followed_ by the V observations.

real(r8)            :: U_obs         = 0.0_r8
real(r8)            :: U_obs_err_var = 0.0_r8
type(location_type) :: U_obs_loc
integer             :: U_flavor
integer             :: U_type        = KIND_V_WIND_COMPONENT ! intentional mismatch
real(r8)            :: U_pr_mean     = 0.0_r8
real(r8)            :: U_pr_sprd     = 0.0_r8
real(r8)            :: U_po_mean     = 0.0_r8
real(r8)            :: U_po_sprd     = 0.0_r8
integer             :: U_qc          = 0

integer :: obs_index, prior_mean_index, posterior_mean_index
integer :: prior_spread_index, posterior_spread_index
integer :: flavor, wflavor ! THIS IS THE (global) 'KIND' in the obs_def_mod list. 
integer :: num_copies, num_qc, num_obs, max_num_obs, obs_seq_file_id
integer :: num_obs_kinds
character(len=129) :: obs_seq_read_format
logical :: pre_I_format
logical :: only_print_locations = .false.

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
integer             :: qc_integer, my_qc_integer
integer, parameter  :: QC_MAX = 7
integer, parameter  :: QC_MAX_PRIOR     = 3
integer, parameter  :: QC_MAX_POSTERIOR = 1
integer, dimension(0:QC_MAX) :: qc_counter = 0
real(r8), allocatable, dimension(:) :: qc
real(r8), allocatable, dimension(:) :: copyvals

!-----------------------------------------------------------------------
! Namelist with (some scalar) default values
!-----------------------------------------------------------------------

character(len = 129) :: obs_sequence_name = "obs_seq.final"
integer, dimension(6) :: first_bin_center = (/ 2003, 1, 1, 0, 0, 0 /)
integer, dimension(6) :: last_bin_center  = (/ 2003, 1, 2, 0, 0, 0 /)
integer, dimension(6) :: bin_separation   = (/    0, 0, 0, 6, 0, 0 /)
integer, dimension(6) :: bin_width        = (/    0, 0, 0, 6, 0, 0 /)
integer, dimension(6) :: time_to_skip     = (/    0, 0, 1, 0, 0, 0 /)
integer :: max_num_bins = 1000       ! maximum number of temporal bins to consider

real(r8), dimension(MaxLevels) :: plevel = MISSING_R8 ! pressure level (hPa)
real(r8), dimension(MaxLevels) :: hlevel = MISSING_R8 ! height (meters)
integer,  dimension(MaxLevels) :: mlevel = MISSING_R8 ! model level (integer index)

integer :: Nregions = 0
real(r8), dimension(MaxRegions) :: lonlim1= MISSING_R8, lonlim2= MISSING_R8
real(r8), dimension(MaxRegions) :: latlim1= MISSING_R8, latlim2= MISSING_R8 
character(len = stringlength), dimension(MaxRegions) :: reg_names = 'null'

real(r8):: rat_cri               = 5000.0_r8 ! QC ratio
real(r8):: input_qc_threshold    = 4.0_r8    ! maximum NCEP QC factor
logical :: print_mismatched_locs = .false.
logical :: print_obs_locations   = .false.
logical :: verbose               = .false.

namelist /obs_diag_nml/ obs_sequence_name, first_bin_center, last_bin_center, &
                       bin_separation, bin_width, time_to_skip, max_num_bins, &
                       plevel, hlevel, mlevel, rat_cri, input_qc_threshold, &
                       Nregions, lonlim1, lonlim2, latlim1, latlim2, &
                       reg_names, print_mismatched_locs, print_obs_locations, verbose

!-----------------------------------------------------------------------
! Variables used to accumulate the statistics.
!-----------------------------------------------------------------------

integer, parameter :: Ncopies = 10
character(len = stringlength), dimension(Ncopies) :: copy_names = &
   (/ 'Nposs      ', 'Nused      ', 'NbadQC     ', 'NbadIZ     ', 'NbadUV     ', &
      'NbadLV     ', 'rmse       ', 'bias       ', 'spread     ', 'totalspread' /)

type TLRV_type
   ! statistics by time-level-region-variable
   integer ::     time_dim = 1
   integer ::    level_dim = 2
   integer ::   region_dim = 3
   integer :: variable_dim = 4
   character(len=8) :: string
   integer :: num_times = 0, num_levels = 0, num_regions = 0, num_variables = 0
   integer,  dimension(:,:,:,:), pointer :: Nposs, Nused
   integer,  dimension(:,:,:,:), pointer :: NbadQC ! # bad (original) QC values
   integer,  dimension(:,:,:,:), pointer :: NbadIZ ! # bad (ie huge) Innovation Zscore
   integer,  dimension(:,:,:,:), pointer :: NbadUV ! # unmatched U/V wind pairs
   integer,  dimension(:,:,:,:), pointer :: NbadLV ! # obs above/below top/bottom
   real(r8), dimension(:,:,:,:), pointer :: rmse, bias, spread, totspread
end type TLRV_type

type LRV_type
   ! statistics (averaged over time) level-region-variable
   integer ::    level_dim = 1
   integer ::   region_dim = 2
   integer :: variable_dim = 3
   character(len=8) :: string
   integer :: num_levels = 0, num_regions = 0, num_variables = 0
   integer,  dimension(:,:,:), pointer :: Nposs, Nused
   integer,  dimension(:,:,:), pointer :: NbadQC ! # bad (original) QC values
   integer,  dimension(:,:,:), pointer :: NbadIZ ! # bad (ie huge) Innovation Zscore
   integer,  dimension(:,:,:), pointer :: NbadUV ! # unmatched U/V wind pairs
   integer,  dimension(:,:,:), pointer :: NbadLV ! # obs above/below top/bottom
   real(r8), dimension(:,:,:), pointer :: rmse, bias, spread, totspread
end type LRV_type

type(TLRV_type) :: analy,    guess
type( LRV_type) :: analyAVG, guessAVG

type(time_type), pointer, dimension(:)   :: bincenter
type(time_type), pointer, dimension(:,:) :: binedges

real(digits12),  pointer, dimension(:)   :: epoch_center
real(digits12),  pointer, dimension(:,:) :: epoch_edges

integer,         pointer, dimension(:) :: obs_used_in_epoch

!-----------------------------------------------------------------------
! General purpose variables
!-----------------------------------------------------------------------

integer  :: Nplevels=0, Nhlevels=0, Nmlevels=0
integer  :: iregion, iepoch, ivar, ifile, num_obs_in_epoch
real(r8) :: obslon, obslat, obslevel, obsloc3(3)

real(r8), dimension(MaxLevels+1) :: plevel_edges = MISSING_R8 ! pressure level (hPa)
real(r8), dimension(MaxLevels+1) :: hlevel_edges = MISSING_R8 ! height (meters)

integer  :: obsindex, i, iunit, lunit, ierr, io

integer  :: ivert
integer  :: level_index
integer  :: Nlevels, ilev   ! counters
integer  :: seconds, days, Nepochs

integer,  allocatable, dimension(:) :: which_vert ! relates kind of level for each obs kind
real(r8), allocatable, dimension(:) :: scale_factor ! to convert to plotting units

character(len = stringlength), pointer, dimension(:) :: my_obs_kind_names

! These pairs of variables are used when we diagnose which observations 
! are far from the background.
integer, parameter :: MaxSigmaBins = 100
integer  :: nsigma(0:MaxSigmaBins) = 0
integer  :: nsigmaUnit, indx

real(r8) :: pr_zscore, po_zscore, zscoreU

type(time_type) :: TimeMin, TimeMax    ! of the entire period of interest
type(time_type) :: beg_time, end_time  ! of the particular bin
type(time_type) :: binsep, binwidth, halfbinwidth 
type(time_type) :: seqT1, seqTN        ! first,last time in entire observation sequence
type(time_type) :: obs_time, skip_time

character(len = 129) :: ncName, locName, msgstring
character(len = stringlength) :: str1, str2, str3

integer  :: Nidentity  = 0   ! identity observations

!=======================================================================
! Get the party started
!=======================================================================

call initialize_utilities('obs_diag')
call register_module(source,revision,revdate) 
call static_init_obs_sequence()  ! Initialize the obs sequence module 

!----------------------------------------------------------------------
! Define/Append the 'wind velocity' obs_kinds to supplant the list declared
! in obs_kind_mod.f90 i.e. if there is a RADIOSONDE_U_WIND_COMPONENT
! and a RADIOSONDE_V_WIND_COMPONENT, there must be a RADIOSONDE_WIND_VELOCITY
! Replace calls to 'get_obs_kind_name' with variable 'my_obs_kind_names'
!----------------------------------------------------------------------

num_obs_kinds = grok_observation_names(my_obs_kind_names)

allocate(which_vert(num_obs_kinds), scale_factor(num_obs_kinds))
which_vert = VERTISUNDEF

!----------------------------------------------------------------------
! Read the namelist
!----------------------------------------------------------------------

call find_namelist_in_file("input.nml", "obs_diag_nml", iunit)
read(iunit, nml = obs_diag_nml, iostat = io)
call check_namelist_read(iunit, io, "obs_diag_nml")

! Record the namelist values used for the run ...
call error_handler(E_MSG,'obs_diag','obs_diag_nml values are',' ',' ',' ')
write(logfileunit,nml=obs_diag_nml)
write(    *      ,nml=obs_diag_nml)

!----------------------------------------------------------------------
! Now that we have input, do some checking and setup
!----------------------------------------------------------------------

call ObsLocationsExist( print_obs_locations )
call set_calendar_type(GREGORIAN)
call Convert2Time(beg_time, end_time, skip_time, binsep, binwidth, halfbinwidth)
call ActOnNamelist( Nregions )

! layer boundaries - heights get sorted to solve single-layer case

Nplevels = Rmidpoints2edges(plevel, plevel_edges)
Nhlevels = Rmidpoints2edges(hlevel, hlevel_edges)
Nmlevels = Imidpoints2edges(mlevel)

hlevel_edges(1:Nhlevels+1) = sort(hlevel_edges(1:Nhlevels+1))

!----------------------------------------------------------------------
! SetTime rectifies user input and the final binning sequence.
!----------------------------------------------------------------------

call SetScaleFactors(scale_factor, logfileunit) ! for plotting purposes

if (verbose) then
   write(*,*)'pressure levels     = ',plevel(    1:Nplevels)
   write(*,*)'pressure interfaces = ',plevel_edges(1:Nplevels+1)
   write(*,*)'height   levels     = ',hlevel(    1:Nhlevels)
   write(*,*)'height   interfaces = ',hlevel_edges(1:Nhlevels+1)
   write(*,*)'model    levels     = ',mlevel(    1:Nmlevels)
   do i = 1,Nregions
      write(*,'(''Region '',i02,1x,a32,'' (WESN): '',4(f10.4,1x))') i, &
             reg_names(i), lonlim1(i),lonlim2(i),latlim1(i),latlim2(i)
   enddo
endif

call SetTime(beg_time, end_time, binsep, binwidth, halfbinwidth, &
     TimeMin, TimeMax, Nepochs, bincenter, binedges, epoch_center, epoch_edges, &
     obs_used_in_epoch)

prior_mean(1)       = 0.0_r8
prior_spread(1)     = 0.0_r8
posterior_mean(1)   = 0.0_r8
posterior_spread(1) = 0.0_r8

U_obs_loc = set_location_missing()

!----------------------------------------------------------------------
! Prepare the variables
!----------------------------------------------------------------------

Nlevels = maxval((/ Nplevels, Nhlevels, Nmlevels /))

write(*,*)'level dimension can be ',Nlevels

allocate(obs_seq_filenames(Nepochs*4))
obs_seq_filenames = 'null'

allocate(guess%rmse(  Nepochs, Nlevels, Nregions, num_obs_kinds), &
         guess%bias(  Nepochs, Nlevels, Nregions, num_obs_kinds), &
         guess%spread(Nepochs, Nlevels, Nregions, num_obs_kinds), &
      guess%totspread(Nepochs, Nlevels, Nregions, num_obs_kinds), &
         guess%Nposs( Nepochs, Nlevels, Nregions, num_obs_kinds), &
         guess%Nused( Nepochs, Nlevels, Nregions, num_obs_kinds), &
         guess%NbadQC(Nepochs, Nlevels, Nregions, num_obs_kinds), &
         guess%NbadIZ(Nepochs, Nlevels, Nregions, num_obs_kinds), &
         guess%NbadUV(Nepochs, Nlevels, Nregions, num_obs_kinds), &
         guess%NbadLV(Nepochs, Nlevels, Nregions, num_obs_kinds), &
         analy%rmse(  Nepochs, Nlevels, Nregions, num_obs_kinds), &
         analy%bias(  Nepochs, Nlevels, Nregions, num_obs_kinds), &
         analy%spread(Nepochs, Nlevels, Nregions, num_obs_kinds), &
      analy%totspread(Nepochs, Nlevels, Nregions, num_obs_kinds), &
         analy%Nposs( Nepochs, Nlevels, Nregions, num_obs_kinds), &
         analy%Nused( Nepochs, Nlevels, Nregions, num_obs_kinds), &
         analy%NbadQC(Nepochs, Nlevels, Nregions, num_obs_kinds), &
         analy%NbadIZ(Nepochs, Nlevels, Nregions, num_obs_kinds), &
         analy%NbadUV(Nepochs, Nlevels, Nregions, num_obs_kinds), &
         analy%NbadLV(Nepochs, Nlevels, Nregions, num_obs_kinds))

guess%rmse   = 0.0_r8
guess%bias   = 0.0_r8
guess%spread = 0.0_r8
guess%totspread = 0.0_r8
guess%Nposs  = 0
guess%Nused  = 0
guess%NbadQC = 0
guess%NbadIZ = 0
guess%NbadUV = 0
guess%NbadLV = 0

guess%string        = "guess"
guess%num_times     = Nepochs
guess%num_levels    = Nlevels
guess%num_regions   = Nregions
guess%num_variables = num_obs_kinds

analy%rmse   = 0.0_r8
analy%bias   = 0.0_r8
analy%spread = 0.0_r8
analy%totspread = 0.0_r8
analy%Nposs  = 0
analy%Nused  = 0
analy%NbadQC = 0
analy%NbadIZ = 0
analy%NbadUV = 0
analy%NbadLV = 0

analy%string        = "analy"
analy%num_times     = Nepochs
analy%num_levels    = Nlevels
analy%num_regions   = Nregions
analy%num_variables = num_obs_kinds

allocate(guessAVG%rmse(  Nlevels, Nregions, num_obs_kinds), &
         guessAVG%bias(  Nlevels, Nregions, num_obs_kinds), &
         guessAVG%spread(Nlevels, Nregions, num_obs_kinds), &
      guessAVG%totspread(Nlevels, Nregions, num_obs_kinds), &
         guessAVG%Nposs( Nlevels, Nregions, num_obs_kinds), &
         guessAVG%Nused( Nlevels, Nregions, num_obs_kinds), &
         guessAVG%NbadQC(Nlevels, Nregions, num_obs_kinds), &
         guessAVG%NbadIZ(Nlevels, Nregions, num_obs_kinds), &
         guessAVG%NbadUV(Nlevels, Nregions, num_obs_kinds), &
         guessAVG%NbadLV(Nlevels, Nregions, num_obs_kinds), &
         analyAVG%rmse(  Nlevels, Nregions, num_obs_kinds), &
         analyAVG%bias(  Nlevels, Nregions, num_obs_kinds), &
         analyAVG%spread(Nlevels, Nregions, num_obs_kinds), &
      analyAVG%totspread(Nlevels, Nregions, num_obs_kinds), &
         analyAVG%Nposs( Nlevels, Nregions, num_obs_kinds), &
         analyAVG%Nused( Nlevels, Nregions, num_obs_kinds), &
         analyAVG%NbadQC(Nlevels, Nregions, num_obs_kinds), &
         analyAVG%NbadIZ(Nlevels, Nregions, num_obs_kinds), &
         analyAVG%NbadUV(Nlevels, Nregions, num_obs_kinds), &
         analyAVG%NbadLV(Nlevels, Nregions, num_obs_kinds))

guessAVG%rmse   = 0.0_r8
guessAVG%bias   = 0.0_r8
guessAVG%spread = 0.0_r8
guessAVG%totspread = 0.0_r8
guessAVG%Nposs  = 0
guessAVG%Nused  = 0
guessAVG%NbadQC = 0
guessAVG%NbadIZ = 0
guessAVG%NbadUV = 0
guessAVG%NbadLV = 0

guessAVG%string        = "VPguess"
guessAVG%num_levels    = Nlevels
guessAVG%num_regions   = Nregions
guessAVG%num_variables = num_obs_kinds

analyAVG%rmse   = 0.0_r8
analyAVG%bias   = 0.0_r8
analyAVG%spread = 0.0_r8
analyAVG%totspread = 0.0_r8
analyAVG%Nposs  = 0
analyAVG%Nused  = 0
analyAVG%NbadQC = 0
analyAVG%NbadIZ = 0
analyAVG%NbadUV = 0
analyAVG%NbadLV = 0

analyAVG%string        = "VPanaly"
analyAVG%num_levels    = Nlevels
analyAVG%num_regions   = Nregions
analyAVG%num_variables = num_obs_kinds

!----------------------------------------------------------------------
! Open file for histogram of innovations, as a function of standard deviation.
!----------------------------------------------------------------------
nsigmaUnit = open_file('LargeInnov.txt',form='formatted',action='rewind')
write(nsigmaUnit,'(a)')'Any observations flagged as bad are dumped into the last bin.'
write(nsigmaUnit,'(a)') '   day   secs    lon      lat    level         obs    guess   zscore   key   kind'
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
!-----------------------------------------------------------------------
! The strategy at this point is to open WAY too many files and 
! check the observation sequences against ALL of the temporal bins.
! If the sequence is completely before the time period of interest, we skip.
! If the sequence is completely past the time period of interest, we stop.

ObsFileLoop : do ifile=1, Nepochs*4
!-----------------------------------------------------------------------

   ! TODO: NextFile only handles relative filenames (as of 12/20/2007).
   ! The assumption is that the file are organized  obs_0001/obs_seq.final   
   ! NextFile finds the first slash instead of the last slash. If everything
   ! were slaved off the last slash I suspect relative and absolute would work.
   ! TJH 2007/12/20

   obs_seq_in_file_name = NextFile(obs_sequence_name,ifile)

   if ( file_exist(trim(obs_seq_in_file_name)) ) then
      write(msgstring,*)'opening ', trim(obs_seq_in_file_name)
      call error_handler(E_MSG,'obs_diag',msgstring,source,revision,revdate)
   else
      write(msgstring,*)trim(obs_seq_in_file_name),&
                        ' does not exist. Finishing up.'
      call error_handler(E_MSG,'obs_diag',msgstring,source,revision,revdate)
      exit ObsFileLoop
   endif

   ! Read in information about observation sequence so we can allocate
   ! observations. We need info about how many copies, qc values, etc.

   obs_seq_in_file_name     = trim(adjustl(obs_seq_in_file_name)) ! Lahey requirement
   obs_seq_filenames(ifile) = trim(adjustl(obs_seq_in_file_name))

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

   if ( verbose ) then
      write(logfileunit,*)
      write(logfileunit,*)'num_copies          is ',num_copies
      write(logfileunit,*)'num_qc              is ',num_qc
      write(logfileunit,*)'num_obs             is ',num_obs
      write(logfileunit,*)'max_num_obs         is ',max_num_obs
      write(logfileunit,*)'obs_seq_read_format is ',trim(adjustl(obs_seq_read_format))
      write(logfileunit,*)'pre_I_format        is ',pre_I_format
      write(    *      ,*)
      write(    *      ,*)'num_copies          is ',num_copies
      write(    *      ,*)'num_qc              is ',num_qc
      write(    *      ,*)'num_obs             is ',num_obs
      write(    *      ,*)'max_num_obs         is ',max_num_obs
      write(    *      ,*)'obs_seq_read_format is ',trim(adjustl(obs_seq_read_format))
      write(    *      ,*)'pre_I_format        is ',pre_I_format
   endif

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
   call print_date(seqT1,'First observation date',logfileunit)
   call print_date(seqTN,'Last  observation date',logfileunit)
   if ( verbose ) then
      call print_time(seqT1,'First observation time')
      call print_time(seqTN,'Last  observation time')
      call print_time(TimeMin,'TimeMin ')
      call print_time(TimeMax,'TimeMax ')
      call print_date(seqT1,'First observation date')
      call print_date(seqTN,'Last  observation date')
      call print_date(TimeMin,'TimeMin ')
      call print_date(TimeMax,'TimeMax ')
   endif

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
         write(logfileunit,*)'seqTN > TimeMin ... using ', &
                             trim(adjustl(obs_seq_in_file_name))
         write(    *      ,*)'seqTN > TimeMin ... using ', &
                             trim(adjustl(obs_seq_in_file_name))
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
      if (verbose) write(*,*)'seqT1 < TimeMax ... using ',trim(adjustl(obs_seq_in_file_name))
   endif

   !--------------------------------------------------------------------
   ! Find the index of obs, ensemble mean, spread ... etc.
   !--------------------------------------------------------------------
   ! Only require obs_index to be present; this allows the program
   ! to be run on obs_seq.in files which have no means or spreads.
   ! You can still plot locations, but that's it.
   !--------------------------------------------------------------------

   call SetIndices( obs_index, qc_index, dart_qc_index, &
            prior_mean_index,   posterior_mean_index,   &
            prior_spread_index, posterior_spread_index  )

   if ( any( (/ prior_mean_index,     prior_spread_index, &
            posterior_mean_index, posterior_spread_index /) < 0) ) then
      only_print_locations = .true.
      msgstring = 'observation sequence has no prior/posterior information'
      call error_handler(E_MSG,'obs_diag',msgstring,source,revision,revdate)
   else
      only_print_locations = .false.
   endif

   !====================================================================
   EpochLoop : do iepoch = 1, Nepochs
   !====================================================================

      beg_time = binedges(1,iepoch)
      end_time = binedges(2,iepoch) 

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

      ! Open each epoch file here if writing out auxiliary location files
      if (print_obs_locations) then
          ! Append epoch number to name
          write(locName,'(a,i3.3,a)') 'observation_locations.', iepoch, '.dat'
          if (file_exist(locName)) then
             lunit = open_file(trim(adjustl(locName)),form='formatted',action='append')
          else
             lunit = open_file(trim(adjustl(locName)),form='formatted',action='rewind')
             write(lunit, '(a)') '     lon         lat      lev       kind     key   QCval'
          endif
      endif

      allocate(keys(num_obs_in_epoch))

      call get_time_range_keys(seq, key_bounds, num_obs_in_epoch, keys)

      !-----------------------------------------------------------------
      ObservationLoop : do obsindex = 1, num_obs_in_epoch
      !-----------------------------------------------------------------

         ! 'flavor' is from the 'master list' in the obs_kind_mod.f90
         ! each obs_seq.final file has their own private kind - which
         ! gets mapped to the 'master list', if you will.

         call get_obs_from_key(seq, keys(obsindex), observation)
         call get_obs_def(observation, obs_def)

         flavor   = get_obs_kind(obs_def)
         obs_time = get_obs_def_time(obs_def)
         obs_loc  = get_obs_def_location(obs_def)
         obsloc3  = get_location(obs_loc)

         obslon   = obsloc3(1) ! [  0, 360]
         obslat   = obsloc3(2) ! [-90,  90]
         obslevel = obsloc3(3) ! variable-dependent

         ! Check to see if it is an identity observation.
         ! If it is, we count them and skip them since they are better
         ! explored with the model-space diagnostics.
         if (flavor < 0) then
            write(*,*)'skipping identity obs at ',obslon,obslat,obslevel
            Nidentity = Nidentity + 1
            cycle ObservationLoop
         endif

         ! address change of units ... DART tracks around Pascals,
         ! however, we like to plot hPa ... that sort of thing.
         ! cvrt to hPa
         if(vert_is_pressure(obs_loc)) obslevel = 0.01_r8 * obsloc3(3)

         ! same sort of thing for the scale factors 
         obs_err_var = get_obs_def_error_variance(obs_def) * &
                       scale_factor(flavor) * &
                       scale_factor(flavor)

         !--------------------------------------------------------------
         ! Figure out which level the observation relates to ...
         !--------------------------------------------------------------

         level_index = ParseLevel(obs_loc, obslevel, flavor, which_vert)

         if ( 1 == 2 ) then
            write(8,*)'obsindx ',obsindex, keys(obsindex), obsloc3(3), level_index
         endif

         !--------------------------------------------------------------
         ! Convert the DART QC data to an integer and create histogram 
         !--------------------------------------------------------------

         call get_qc(observation, qc)

         if (dart_qc_index > 0) then
            qc_integer = min( nint(qc(dart_qc_index)), QC_MAX )
            qc_counter(qc_integer) = qc_counter(qc_integer) + 1  ! histogram
         else
            ! If there is no dart_qc in obs_seq, make sure the observation
            ! is never used. This must be a case where we are interested
            ! only in getting the location information.
            qc_integer = QC_MAX + 9999
         endif

         !--------------------------------------------------------------
         ! Write location of observation if namelist item is true
         !--------------------------------------------------------------

         if (print_obs_locations) then

            if (dart_qc_index > 0) then 
               my_qc_integer =      nint(qc(dart_qc_index))
            elseif  (qc_index > 0) then
               my_qc_integer = -1 * nint(qc(     qc_index))
            else
               my_qc_integer = -99
            endif

            write(lunit, FMT='(3(f10.2,1x),3(i7,1x))') &
               obslon, obslat, obslevel, flavor, keys(obsindex), my_qc_integer
         endif

         !--------------------------------------------------------------
         ! Early exit from the observation loop if the observation 
         ! does not have all the required copies (attributes).
         !--------------------------------------------------------------

         if (only_print_locations) cycle ObservationLoop

         !--------------------------------------------------------------
         ! retrieve observation prior and posterior means and spreads
         !--------------------------------------------------------------

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

         !--------------------------------------------------------------
         ! Check to see if there are any observations with wild values
         ! and a DART QC flag that is inconsistent. I checked it once
         ! with qc_integer < 4 and found that ONLY the posterior values
         ! were bad. While the underlying bug is being fixed, the workaround
         ! is to simply manually set the DART QC value such that the
         ! posterior is flagged as bad.  I cannot acutally tell if the
         ! observation is supposed to have been assimilated or just evaluated,
         ! so I erred on the conservative side. TJH 24 Aug 2007
         !--------------------------------------------------------------

         ! integer, parameter  :: QC_MAX = 7
         ! integer, parameter  :: QC_MAX_PRIOR     = 3
         ! integer, parameter  :: QC_MAX_POSTERIOR = 1

         if ( 1 == 2 ) then
         call get_obs_values(observation, copyvals)
     !   if (any(copyvals < -87.0) .and. ( qc_integer < (QC_MAX_PRIOR+1) ) ) then
         if ( obsindex == 311 ) then
              write(*,*)
              write(*,*)'Observation index 311 is ',keys(obsindex),' and has:'
              write(*,*)'flavor                 is',flavor
              write(*,*)'obs              value is',obs(1)
              write(*,*)'prior_mean       value is',prior_mean(1)
              write(*,*)'posterior_mean   value is',posterior_mean(1)
              write(*,*)'prior_spread     value is',prior_spread(1)
              write(*,*)'posterior_spread value is',posterior_spread(1)
              write(*,*)'DART QC          value is',qc_integer
              do i= 1,num_copies 
                 write(*,*)copyvals(i),trim(get_copy_meta_data(seq,i))
              enddo
              write(*,*)
         endif
         endif

         !--------------------------------------------------------------
         ! Scale the quantities so they plot sensibly.
         !--------------------------------------------------------------

         obs(1)  = obs(1)             *scale_factor(flavor)
         pr_mean = prior_mean(1)      *scale_factor(flavor)
         po_mean = posterior_mean(1)  *scale_factor(flavor)
         pr_sprd = prior_spread(1)    *scale_factor(flavor)
         po_sprd = posterior_spread(1)*scale_factor(flavor)

         !--------------------------------------------------------------
         ! (DEBUG) Summary of observation knowledge at this point
         !--------------------------------------------------------------

         if ( 1 == 2 ) then
            write(*,*)'observation # ',obsindex
            write(*,*)'obs_flavor ',flavor
            write(*,*)'obs_err_var ',obs_err_var
            write(*,*)'obslon/obslat ',obslon,obslat
            write(*,*)'qc ',qc
            write(*,*)'obs(1) ',obs(1)
            write(*,*)'pr_mean,po_mean ',pr_mean,po_mean
            write(*,*)'pr_sprd,po_sprd ',pr_sprd,po_sprd
         endif

         !--------------------------------------------------------------
         ! update the histogram of the magnitude of the innovation,
         ! where each bin is a single standard deviation. 
         ! This is a one-sided histogram. 
         !--------------------------------------------------------------

         pr_zscore = InnovZscore(obs(1), pr_mean, pr_sprd, obs_err_var, qc_integer, QC_MAX_PRIOR)
         po_zscore = InnovZscore(obs(1), po_mean, po_sprd, obs_err_var, qc_integer, QC_MAX_POSTERIOR)

         indx         = min(int(pr_zscore), MaxSigmaBins)
         nsigma(indx) = nsigma(indx) + 1

         ! Individual (valid) observations that are very far away get
         ! logged to a separate file.

         if( (pr_zscore > 3.0_r8) .and. (qc_integer <= QC_MAX_PRIOR) ) then
            call get_time(obs_time,seconds,days)

            write(nsigmaUnit,FMT='(i7,1x,i5,1x,2f8.2,i7,1x,2f13.2,f8.1,2i7)') &
                 days, seconds, obslon, obslat, ivert, &
                 obs(1), pr_mean, pr_zscore, keys(obsindex), flavor
         endif

         obs_used_in_epoch(iepoch) = obs_used_in_epoch(iepoch) + 1

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
            U_pr_mean     = pr_mean
            U_pr_sprd     = pr_sprd
            U_po_mean     = po_mean
            U_po_sprd     = po_sprd
            U_qc          = qc_integer

         endif

         !--------------------------------------------------------------
         ! We have Nregions of interest
         !--------------------------------------------------------------

         Areas : do iregion =1, Nregions

            keeper = InRegion( obslon, obslat, lonlim1(iregion), lonlim2(iregion), &
                                           latlim1(iregion), latlim2(iregion))
            if ( .not. keeper ) cycle Areas

            !-----------------------------------------------------------
            ! Reject observations too high or too low.
            !-----------------------------------------------------------

            if ( level_index < 1 .or. level_index > Nlevels )   then
               guess%NbadLV(iepoch,:,iregion,flavor) = &
               guess%NbadLV(iepoch,:,iregion,flavor) + 1
               analy%NbadLV(iepoch,:,iregion,flavor) = &
               analy%NbadLV(iepoch,:,iregion,flavor) + 1
               cycle Areas
            endif

            !-----------------------------------------------------------
            ! Count original (bad) QC values
            !-----------------------------------------------------------

            if( qc_index > 0 ) then
               if (qc(qc_index) >= input_qc_threshold ) then
               call IPE(guess%NbadQC(iepoch,level_index,iregion,flavor), 1)
               call IPE(analy%NbadQC(iepoch,level_index,iregion,flavor), 1)
               endif
            endif

            !-----------------------------------------------------------
            ! Count 'large' innovations
            !-----------------------------------------------------------

            if( pr_zscore > rat_cri ) then
               call IPE(guess%NbadIZ(iepoch,level_index,iregion,flavor), 1)
            endif

            if( po_zscore > rat_cri ) then
               call IPE(analy%NbadIZ(iepoch,level_index,iregion,flavor), 1)
            endif

            !-----------------------------------------------------------
            ! Do all the heavy lifting
            !-----------------------------------------------------------

            call Bin4D(qc_integer, iepoch, level_index, iregion, flavor, &
                     obs(1), obs_err_var, pr_mean, pr_sprd, po_mean, po_sprd)

            !-----------------------------------------------------------
            ! Additional work for wind speed (given U,V)
            !-----------------------------------------------------------

            ObsIsWindCheck: if ( get_obs_kind_var_type(flavor) == KIND_V_WIND_COMPONENT ) then

               ! The big assumption is that the U wind component has
               ! immediately preceeded the V component and has been saved.
               ! 
               ! We check for observation compatibility and the index for 
               ! this wind 'kind' ... not originally in the max_obs_kind namelist.
               ! this will be the 'wflavor' (wind) flavor.

               ierr = CheckMate(flavor, U_flavor, obs_loc, U_obs_loc, wflavor ) 

               if ( ierr /= 0 ) then
                  write(*,*)'time series ... V with no matching U ...'
                  call IPE(guess%NbadUV(iepoch,level_index,iregion,wflavor), 1)
                  call IPE(analy%NbadUV(iepoch,level_index,iregion,wflavor), 1)
               else

                  ! The next big assumption is that the 'wind speed' flavors
                  ! need to have their which_vert explicitly set so the netcdf
                  ! files know what sort of vertical coordinate to use. The only
                  ! way I know of is to use the observation of that type to tell us.
                  ! (over and over and over ... last one wins!)

                  ierr = ParseLevel(obs_loc, obslevel, wflavor, which_vert)

                  ! since we don't have the necessary covariance between U,V
                  ! we will reject if either univariate z score is bad 

                  zscoreU = InnovZscore(U_obs, U_pr_mean, U_pr_sprd, U_obs_err_var, U_qc, QC_MAX_PRIOR)
                  if( (pr_zscore > rat_cri) .or. (zscoreU > rat_cri) )  then
                     call IPE(guess%NbadIZ(iepoch,level_index,iregion,wflavor), 1)
                  endif

                  zscoreU = InnovZscore(U_obs, U_po_mean, U_po_sprd, U_obs_err_var, U_qc, QC_MAX_POSTERIOR)
                  if( (po_zscore > rat_cri) .or. (zscoreU > rat_cri) )  then
                     call IPE(analy%NbadIZ(iepoch,level_index,iregion,wflavor), 1)
                  endif

                  call Bin4D(qc_integer, iepoch, level_index, iregion, wflavor, &
                       obs(1),  obs_err_var,   pr_mean,   pr_sprd,   po_mean,   po_sprd,  &
                       U_obs, U_obs_err_var, U_pr_mean, U_pr_sprd, U_po_mean, U_po_sprd   )
               endif

            endif ObsIsWindCheck

            !-----------------------------------------------------------
            ! end of time series statistics
            !-----------------------------------------------------------

            if (which_vert(flavor) == VERTISSURFACE) cycle Areas

            if ( obs_time < skip_time ) cycle Areas

            !-----------------------------------------------------------
            ! vertical statistical part
            !-----------------------------------------------------------

            if( qc_index > 0 ) then
               if (qc(qc_index) >= input_qc_threshold ) then
               call IPE(guessAVG%NbadQC(level_index,iregion,flavor), 1)
               call IPE(analyAVG%NbadQC(level_index,iregion,flavor), 1)
               endif
            endif

            if(pr_zscore > rat_cri )  then
               call IPE(guessAVG%NbadIZ(level_index,iregion,flavor), 1)
            endif

            if(po_zscore > rat_cri )  then
               call IPE(analyAVG%NbadIZ(level_index,iregion,flavor), 1)
            endif

            call Bin3D(qc_integer, level_index,   iregion,    flavor, &
                   obs(1),  obs_err_var,  pr_mean,   pr_sprd,   po_mean,   po_sprd  )

            ! Handle wind speed given U,V components 

            if ( get_obs_kind_var_type(flavor) == KIND_V_WIND_COMPONENT ) then

               ierr = CheckMate(flavor, U_flavor, obs_loc, U_obs_loc, wflavor)
               if ( ierr /= 0 ) then
                  call IPE(guessAVG%NbadUV(level_index,iregion,wflavor), 1)
                  call IPE(analyAVG%NbadUV(level_index,iregion,wflavor), 1)
               else

                  ierr = ParseLevel(obs_loc, obslevel, wflavor, which_vert)

                  zscoreU = InnovZscore(U_obs, U_pr_mean, U_pr_sprd, U_obs_err_var, U_qc, QC_MAX_PRIOR)
                  if( (pr_zscore > rat_cri) .or. (zscoreU > rat_cri) )  then
                     call IPE(guessAVG%NbadIZ(level_index,iregion,wflavor), 1)
                  endif

                  zscoreU = InnovZscore(U_obs, U_po_mean, U_po_sprd, U_obs_err_var, U_qc, QC_MAX_POSTERIOR)
                  if( (po_zscore > rat_cri) .or. (zscoreU > rat_cri) )  then
                     call IPE(analyAVG%NbadIZ(level_index,iregion,wflavor), 1)
                  endif

                  call Bin3D(qc_integer, level_index,   iregion,    wflavor, &
                      obs(1),  obs_err_var,   pr_mean,   pr_sprd,   po_mean,   po_sprd, &
                      U_obs, U_obs_err_var, U_pr_mean, U_pr_sprd, U_po_mean, U_po_sprd  ) 
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

      if (print_obs_locations) close(lunit)

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
   if (allocated(qc)) deallocate( qc )
   if (allocated(copyvals)) deallocate( copyvals )

enddo ObsFileLoop

!-----------------------------------------------------------------------
! We have read all possible files, and stuffed the observations into the
! appropriate bins. Time to normalize.
!-----------------------------------------------------------------------

if (verbose) then
   do ivar   = 1,SIZE(which_vert)
      write(logfileunit,*)'which_vert(',ivar,' of ',num_obs_kinds,') = ',which_vert(ivar)
      write(     *     ,*)'which_vert(',ivar,' of ',num_obs_kinds,') = ',which_vert(ivar)
   enddo
endif

if (verbose) then
   write(logfileunit,*)'Normalizing time-level-region-variable quantities.'
   write(     *     ,*)'Normalizing time-level-region-variable quantities.'
endif

do ivar   = 1,num_obs_kinds
do iregion= 1,Nregions
do ilev   = 1,Nlevels
do iepoch = 1,Nepochs
   if (  guess%Nused( iepoch, ilev, iregion, ivar) == 0) then
         guess%rmse(  iepoch, ilev, iregion, ivar) = MISSING_R4
         guess%bias(  iepoch, ilev, iregion, ivar) = MISSING_R4
         guess%spread(iepoch, ilev, iregion, ivar) = MISSING_R4
      guess%totspread(iepoch, ilev, iregion, ivar) = MISSING_R4
   else
         guess%rmse(  iepoch, ilev, iregion, ivar) = &
    sqrt(guess%rmse(  iepoch, ilev, iregion, ivar) / &
         guess%Nused( iepoch, ilev, iregion, ivar) )

         guess%bias(  iepoch, ilev, iregion, ivar) = &
         guess%bias(  iepoch, ilev, iregion, ivar) / &
         guess%Nused( iepoch, ilev, iregion, ivar)  

         guess%spread(iepoch, ilev, iregion, ivar) = &
    sqrt(guess%spread(iepoch, ilev, iregion, ivar) / &
         guess%Nused( iepoch, ilev, iregion, ivar) )

         guess%totspread(iepoch, ilev, iregion, ivar) = &
    sqrt(guess%totspread(iepoch, ilev, iregion, ivar) / &
         guess%Nused(    iepoch, ilev, iregion, ivar) )

   endif

   if (  analy%Nused( iepoch, ilev, iregion, ivar) == 0) then
         analy%rmse(  iepoch, ilev, iregion, ivar) = MISSING_R4
         analy%bias(  iepoch, ilev, iregion, ivar) = MISSING_R4
         analy%spread(iepoch, ilev, iregion, ivar) = MISSING_R4
      analy%totspread(iepoch, ilev, iregion, ivar) = MISSING_R4
   else
         analy%rmse(  iepoch, ilev, iregion, ivar) = &
    sqrt(analy%rmse(  iepoch, ilev, iregion, ivar) / &
         analy%Nused( iepoch, ilev, iregion, ivar) )

         analy%bias(  iepoch, ilev, iregion, ivar) = &
         analy%bias(  iepoch, ilev, iregion, ivar) / &
         analy%Nused( iepoch, ilev, iregion, ivar)  

         analy%spread(iepoch, ilev, iregion, ivar) = &
    sqrt(analy%spread(iepoch, ilev, iregion, ivar) / &
         analy%Nused( iepoch, ilev, iregion, ivar) )

         analy%totspread(iepoch, ilev, iregion, ivar) = &
    sqrt(analy%totspread(iepoch, ilev, iregion, ivar) / &
         analy%Nused(    iepoch, ilev, iregion, ivar) )

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
close(nsigmaUnit)

!-----------------------------------------------------------------------
! temporal average of the vertical statistics
!-----------------------------------------------------------------------
if (verbose) then
   write(logfileunit,*)'Normalize quantities for all levels.'
   write(     *     ,*)'Normalize quantities for all levels.'
endif

do ivar=1,num_obs_kinds
do iregion=1, Nregions
do ilev=1, Nlevels

   if (    guessAVG%Nused( ilev, iregion, ivar) == 0) then
           guessAVG%rmse(  ilev, iregion, ivar) = MISSING_R4
           guessAVG%bias(  ilev, iregion, ivar) = MISSING_R4
           guessAVG%spread(ilev, iregion, ivar) = MISSING_R4
        guessAVG%totspread(ilev, iregion, ivar) = MISSING_R4

   else
           guessAVG%rmse(  ilev, iregion, ivar) = &
      sqrt(guessAVG%rmse(  ilev, iregion, ivar) / &
           guessAVG%Nused( ilev, iregion, ivar) )

           guessAVG%bias(  ilev, iregion, ivar) = &
           guessAVG%bias(  ilev, iregion, ivar) / &
           guessAVG%Nused( ilev, iregion, ivar)

           guessAVG%spread(ilev, iregion, ivar) = &
      sqrt(guessAVG%spread(ilev, iregion, ivar) / &
           guessAVG%Nused( ilev, iregion, ivar) )

           guessAVG%totspread(ilev, iregion, ivar) = &
      sqrt(guessAVG%totspread(ilev, iregion, ivar) / &
           guessAVG%Nused(    ilev, iregion, ivar) )

   endif

   if (    analyAVG%Nused( ilev, iregion, ivar) == 0) then
           analyAVG%rmse(  ilev, iregion, ivar) = MISSING_R4
           analyAVG%bias(  ilev, iregion, ivar) = MISSING_R4
           analyAVG%spread(ilev, iregion, ivar) = MISSING_R4
        analyAVG%totspread(ilev, iregion, ivar) = MISSING_R4

   else
           analyAVG%rmse(  ilev, iregion, ivar) = &
      sqrt(analyAVG%rmse(  ilev, iregion, ivar) / &
           analyAVG%Nused( ilev, iregion, ivar) )

           analyAVG%bias(  ilev, iregion, ivar) = &
           analyAVG%bias(  ilev, iregion, ivar) / &
           analyAVG%Nused( ilev, iregion, ivar)

           analyAVG%spread(ilev, iregion, ivar) = &
      sqrt(analyAVG%spread(ilev, iregion, ivar) / &
           analyAVG%Nused( ilev, iregion, ivar) )

           analyAVG%totspread(ilev, iregion, ivar) = &
      sqrt(analyAVG%totspread(ilev, iregion, ivar) / &
           analyAVG%Nused(    ilev, iregion, ivar) )

   endif
enddo
enddo
enddo

!-----------------------------------------------------------------------
! Print final rejection summary.
!-----------------------------------------------------------------------

write(*,*)
write(*,*) '# observations used  : ',sum(obs_used_in_epoch)
write(*,*) 'Observation summary:'
write(*,*) '# identity           : ',Nidentity
write(*,*) '# bad Innov  (ratio) : ',sum(guess%NbadIZ)
write(*,*) '# bad UV (wind pairs): ',sum(guess%NbadUV)
write(*,*) '# bad Level          : ',sum(guess%NbadLV(:,1,:,:))
write(*,*) '# bad DART QC        : ',sum(guess%NbadQC)
write(*,*)

write(logfileunit,*)
write(logfileunit,*) '# observations used  : ',sum(obs_used_in_epoch)
write(logfileunit,*) 'Observation summary:'
write(logfileunit,*) '# identity           : ',Nidentity
write(logfileunit,*) '# bad Innov  (ratio) : ',sum(guess%NbadIZ)
write(logfileunit,*) '# bad UV (wind pairs): ',sum(guess%NbadUV)
write(logfileunit,*) '# bad Level          : ',sum(guess%NbadLV(:,1,:,:))
write(logfileunit,*) '# bad DART QC        : ',sum(guess%NbadQC)

if (Nidentity > 0) then
   write(*,*)'There were identity observations in this observation sequence file.'
   write(*,*)'At present, obs_diag does not support the analysis of identity '
   write(*,*)'observations. In general, identity observation are explored with'
   write(*,*)'state space diagnostics, i.e. take a peek in the matlab directory. '

   write(logfileunit,*)'There were identity observations in this observation sequence file.'
   write(logfileunit,*)'At present, obs_diag does not support the analysis of identity '
   write(logfileunit,*)'observations. In general, identity observation are explored with'
   write(logfileunit,*)'state space diagnostics, i.e. take a peek in the matlab directory'
endif

!----------------------------------------------------------------------
! Open netCDF output file 
!----------------------------------------------------------------------

ncName = 'obs_diag_output.nc'
  
call WriteNetCDF(ncName)

!-----------------------------------------------------------------------
! Really, really, done.
!-----------------------------------------------------------------------

deallocate(obs_seq_filenames)

deallocate(guess%NbadQC, guess%NbadIZ, guess%NbadUV, guess%NbadLV, &
           analy%NbadQC, analy%NbadIZ, analy%NbadUV, analy%NbadLV, &
           guess%rmse, guess%bias, guess%spread, guess%totspread, guess%Nposs, guess%Nused,  &
           analy%rmse, analy%bias, analy%spread, analy%totspread, analy%Nposs, analy%Nused)

deallocate(guessAVG%NbadQC, guessAVG%NbadIZ, guessAVG%NbadUV, guessAVG%NbadLV, & 
           analyAVG%NbadQC, analyAVG%NbadIZ, analyAVG%NbadUV, analyAVG%NbadLV, &
           guessAVG%rmse, guessAVG%bias, guessAVG%spread, guessAVG%totspread, &
           analyAVG%rmse, analyAVG%bias, analyAVG%spread, analyAVG%totspread, &
           guessAVG%Nposs, guessAVG%Nused,  &
           analyAVG%Nposs, analyAVG%Nused)

deallocate(epoch_center, epoch_edges, bincenter, obs_used_in_epoch)

deallocate(my_obs_kind_names, which_vert, scale_factor)

call timestamp(source,revision,revdate,'end') ! That closes the log file, too.

!======================================================================
CONTAINS
!======================================================================
! These routines use common variables from the scope of this file.
! If it's not in the argument list ... it's scoped within this file.
!======================================================================

   Subroutine  SetScaleFactors(scale_factor, logfileunit)

      real(r8), dimension(:), intent(inout) :: scale_factor
      integer,                intent(in)    :: logfileunit

      ! The surface pressure in the obs_sequence is in Pa, we want to convert
      ! from Pa to hPa for plotting. The specific humidity is a similar thing.
      ! In the obs_sequence file, the units are kg/kg, we want to plot
      ! in the g/kg world...

      ! If kind_surface_pressure or ... does not exist, we are in trouble here.
      ! the scale_factor should be defined to reflect the type, which are not
      ! guaranteed to be numbered sequentially ... vortices 81, for example

      character(len = stringlength) :: obs_kind_name
      integer :: ivar

      scale_factor = 1.0_r8
   
      ! The scale_factor array is dimensioned from obs_kind_mod:num_obs_kinds

      do ivar = 1,SIZE(scale_factor)

         obs_kind_name = my_obs_kind_names(ivar)

         if ( index(obs_kind_name,'SURFACE_PRESSURE') > 0 ) &
                scale_factor(ivar) = 0.01_r8

         if ( index(obs_kind_name,'SPECIFIC_HUMIDITY') > 0 ) &
                scale_factor(ivar) = 1000.0_r8

         ! Somehow, we should plot statistics on the dBZ scale for these ...
         ! scale_factor(KIND_RADAR_REFLECTIVITY) = 10log10(z)

         write(     *     ,*)'scaling of ',scale_factor(ivar),obs_kind_name
         write(logfileunit,*)'scaling of ',scale_factor(ivar),obs_kind_name

      enddo

   end Subroutine SetScaleFactors



   Subroutine SetTime(beg_time, end_time, binsep, binwidth, halfbinwidth, &
        TimeMin, TimeMax, Nepochs, bincenter, binedges, epoch_center, epoch_edges, &
        obs_used_in_epoch)

      type(time_type), intent(in)  :: beg_time, end_time  ! of the particular bin
      type(time_type), intent(in)  :: binsep, binwidth, halfbinwidth
      type(time_type), intent(out) :: TimeMin, TimeMax
      integer,         intent(out) :: Nepochs

      type(time_type), pointer, dimension(  :) :: bincenter
      type(time_type), pointer, dimension(:,:) :: binedges
      real(digits12),  pointer, dimension(  :) :: epoch_center
      real(digits12),  pointer, dimension(:,:) :: epoch_edges
      integer,         pointer, dimension(  :) :: obs_used_in_epoch

      ! Global-scope variables in use in this routine:
      ! max_num_bins    comes from namelist
      ! verbose         comes from namelist
      ! msgstring
      ! logfileunit

      ! Determine temporal bin characteristics.
      ! The user input is not guaranteed to align on bin centers. 
      ! So -- we will assume the start time is correct and take strides till we
      ! get past the last time of interest. 
      ! Nepochs will be the total number of time intervals of the period requested.

      integer :: iepoch, seconds, days
 
      TimeMin  = beg_time
      NepochLoop : do iepoch = 1,max_num_bins
         Nepochs = iepoch
         TimeMax = TimeMin + binsep
         if ( TimeMax > end_time ) exit NepochLoop
         TimeMin = TimeMax
      enddo NepochLoop

      write(*,*)'Requesting ',Nepochs,' assimilation periods.'

      allocate(   bincenter(Nepochs),    binedges(2,Nepochs), &
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
      bincenter( iepoch)    = beg_time
      binedges(2,iepoch)    = beg_time + halfbinwidth

      if (beg_time <= halfbinwidth ) then
         binedges(1,iepoch) = set_time(0,0)
      else
         binedges(1,iepoch) = beg_time - halfbinwidth + set_time(1,0)
      endif

      call get_time(bincenter(iepoch),seconds,days)
      epoch_center(iepoch) = days + seconds/86400.0_digits12

      call get_time(binedges(1,iepoch),seconds,days)
      epoch_edges(1,iepoch) = days + seconds/86400.0_digits12

      call get_time(binedges(2,iepoch),seconds,days)
      epoch_edges(2,iepoch) = days + seconds/86400.0_digits12

      BinLoop : do iepoch = 2,Nepochs

         bincenter( iepoch) = bincenter(iepoch-1) + binsep
         binedges(1,iepoch) = bincenter(iepoch) - halfbinwidth + set_time(1,0)
         binedges(2,iepoch) = bincenter(iepoch) + halfbinwidth

         call get_time(bincenter(iepoch),seconds,days)
         epoch_center(iepoch) = days + seconds/86400.0_digits12

         call get_time(binedges(1,iepoch),seconds,days)
         epoch_edges(1,iepoch) = days + seconds/86400.0_digits12

         call get_time(binedges(2,iepoch),seconds,days)
         epoch_edges(2,iepoch) = days + seconds/86400.0_digits12

      enddo BinLoop

      if ( verbose ) then
      do iepoch = 1,Nepochs
         write(logfileunit,*)
         write(     *     ,*)
         write(str1,'(''epoch '',i6,''  start'')')iepoch
         write(str2,'(''epoch '',i6,'' center'')')iepoch
         write(str3,'(''epoch '',i6,''    end'')')iepoch
         call print_time(binedges(1,iepoch),str1,logfileunit)
         call print_time(bincenter( iepoch),str2,logfileunit)
         call print_time(binedges(2,iepoch),str3,logfileunit)
         call print_time(binedges(1,iepoch),str1)
         call print_time(bincenter( iepoch),str2)
         call print_time(binedges(2,iepoch),str3)

         call print_date(binedges(1,iepoch),str1,logfileunit)
         call print_date(bincenter( iepoch),str2,logfileunit)
         call print_date(binedges(2,iepoch),str3,logfileunit)
         call print_date(binedges(1,iepoch),str1)
         call print_date(bincenter( iepoch),str2)
         call print_date(binedges(2,iepoch),str3)
      enddo
      write(logfileunit,*)
      write(     *     ,*)
      endif

      TimeMin = binedges(1,      1) ! minimum time of interest
      TimeMax = binedges(2,Nepochs) ! maximum time of interest

   end Subroutine SetTime



   Subroutine Convert2Time( beg_time, end_time, skip_time, binsep, &
              binwidth, halfbinwidth)

   ! Global-scope variables in use in this routine:
   ! first_bin_center   comes from namelist
   !  last_bin_center   comes from namelist
   !   bin_separation   comes from namelist
   !        bin_width   comes from namelist
   !     time_to_skip   comes from namelist

   ! We are using bin_separation and bin_width as offsets relative to the
   ! first time. to ensure this, the year and month must be zero.

   type(time_type), intent(out) :: beg_time     ! first_bin_center
   type(time_type), intent(out) :: end_time     ! last_bin_center
   type(time_type), intent(out) :: skip_time    ! time AFTER first_bin_center
   type(time_type), intent(out) :: binsep       ! time between bin centers
   type(time_type), intent(out) :: binwidth     ! period of interest around center
   type(time_type), intent(out) :: halfbinwidth ! half that period

   logical :: error_out = .false.
   integer :: seconds

   ! do some error-checking first

   if ( (bin_separation(1) /= 0) .or. (bin_separation(2) /= 0) ) then
      write(msgstring,*)'bin_separation:year,month must both be zero, they are ', &
      bin_separation(1),bin_separation(2)
      call error_handler(E_MSG,'obs_diag:Convert2Time',msgstring,source,revision,revdate)
      error_out = .true.
   endif

   if ( (bin_width(1) /= 0) .or. (bin_width(2) /= 0) ) then
      write(msgstring,*)'bin_width:year,month must both be zero, they are ', &
      bin_width(1),bin_width(2)
      call error_handler(E_MSG,'obs_diag:Convert2Time',msgstring,source,revision,revdate)
      error_out = .true.
   endif

   if ( (time_to_skip(1) /= 0) .or. (time_to_skip(2) /= 0) ) then
      write(msgstring,*)'time_to_skip:year,month must both be zero, they are ', &
      time_to_skip(1),time_to_skip(2)
      call error_handler(E_MSG,'obs_diag:Convert2Time',msgstring,source,revision,revdate)
      error_out = .true.
   endif

   if ( error_out ) call error_handler(E_ERR,'obs_diag:Convert2Time', &
       'namelist parameter out-of-bounds. Fix and try again.',source,revision,revdate)

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

   seconds  = bin_width(4)*3600 + bin_width(5)*60 + bin_width(6)
   binwidth = set_time(seconds, bin_width(3))

   halfbinwidth = binwidth / 2


   seconds   = time_to_skip(4)*3600 + time_to_skip(5)*60 + time_to_skip(6)

   if ( (beg_time + set_time(seconds,time_to_skip(3))) <= halfbinwidth) then
      skip_time = set_time(0,0)
   else
      skip_time = beg_time - halfbinwidth + set_time(seconds, time_to_skip(3)) 
   endif

   if ( verbose ) then    
      call print_date(    beg_time,'requested first bincenter date',logfileunit)
      call print_date(    end_time,'requested last  bincenter date',logfileunit)
      call print_date(   skip_time,'implied   skip-to         date',logfileunit)
      call print_time(    beg_time,'requested first bincenter time',logfileunit)
      call print_time(    end_time,'requested last  bincenter time',logfileunit)
      call print_time(   skip_time,'implied   skip-to         time',logfileunit)
      call print_time(      binsep,'requested bin separation',logfileunit)
      call print_time(    binwidth,'requested bin      width',logfileunit)
      call print_time(halfbinwidth,'implied     halfbinwidth',logfileunit)

      call print_date(    beg_time,'requested first bincenter date')
      call print_date(    end_time,'requested last  bincenter date')
      call print_date(   skip_time,'implied skip-to           date')
      call print_time(    beg_time,'requested first bincenter time')
      call print_time(    end_time,'requested last  bincenter time')
      call print_time(   skip_time,'implied skip-to           time')
      call print_time(      binsep,'requested bin separation')
      call print_time(    binwidth,'requested bin      width')
      call print_time(halfbinwidth,'implied     halfbinwidth')
   endif

   end Subroutine Convert2Time


   Subroutine ActOnNamelist(Nregions)

   ! Many of the namelist variables are set to 'difficult' defaults.
   ! This routine attempts to divine whether or not the defaults 
   ! have been changed by the namelist. 
   !
   ! globally-scoped variables used 
   ! plevel
   ! hlevel
   ! mlevel
   ! lonlim1, lonlim2
   ! latlim1, latlim2
   ! reg_names

   !----------------------------------------------------------------------
   ! Provide some level of backward compatiblilty ... if Nregions is more
   ! than 1, the assumption is that everything is set and no action needs
   ! to be taken. If it is less than 1, we are trying to provide some
   ! 'the old' defaults.
   !----------------------------------------------------------------------

   integer,                intent(inout) :: Nregions

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

   ! Pressure levels

   if (all(plevel <  1.0)) then
       plevel(1:11) = (/ 1000,  925,  850,   700,   500,  400, &
                          300,  250,  200,   150,   100        /)
   endif

   ! Height levels

   if (all(hlevel < 1.0)) then
       hlevel(1:11) = (/ 1000, 2000, 3000,  4000,  5000, 6000, &
                         7000, 8000, 9000, 10000, 11000        /)
   endif

   ! Model levels

   if (all(mlevel < 1)) then
       mlevel(1:11) = (/ (i,i=1,11) /)   ! set model levels to indices
   endif

   end Subroutine ActOnNamelist



   Subroutine  SetIndices( obs_index, qc_index, dart_qc_index,        &
                           prior_mean_index,   posterior_mean_index,  &
                           prior_spread_index, posterior_spread_index ) 

   integer, intent(out) :: obs_index, qc_index, dart_qc_index, &
                           prior_mean_index,   posterior_mean_index, &
                           prior_spread_index, posterior_spread_index

   ! Using 'seq' from global scope

   integer :: i

   obs_index              = -1
   prior_mean_index       = -1
   posterior_mean_index   = -1
   prior_spread_index     = -1
   posterior_spread_index = -1
   qc_index               = -1
   dart_qc_index          = -1

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
   if (               qc_index < 0 ) then 
      write(msgstring,*)'metadata:Quality Control not found' 
      call error_handler(E_MSG,'obs_diag',msgstring,source,revision,revdate)
   endif
   if (          dart_qc_index < 0 ) then 
      write(msgstring,*)'metadata:DART quality control not found' 
      call error_handler(E_MSG,'obs_diag',msgstring,source,revision,revdate)
   endif

   ! Only require obs_index to be present; this allows the program
   ! to be run on obs_seq.in files which have no means or spread.  You get
   ! less info from them, but you can still plot locations, etc.

   if ( obs_index < 0 ) then
      write(msgstring,*)'metadata:observation not found'
      call error_handler(E_ERR,'obs_diag',msgstring,source,revision,revdate)
   endif

   !--------------------------------------------------------------------
   ! Echo what we found.
   !--------------------------------------------------------------------

   write(msgstring,'(''observation      index '',i2,'' metadata '',a)') &
        obs_index, trim(adjustl(get_copy_meta_data(seq,obs_index)))
   call error_handler(E_MSG,'obs_diag',msgstring,source,revision,revdate)

   if (prior_mean_index > 0 ) then
      write(msgstring,'(''prior mean       index '',i2,'' metadata '',a)') &
           prior_mean_index, trim(adjustl(get_copy_meta_data(seq,prior_mean_index)))
      call error_handler(E_MSG,'obs_diag',msgstring,source,revision,revdate)
   endif

   if (posterior_mean_index > 0 ) then
      write(msgstring,'(''posterior mean   index '',i2,'' metadata '',a)') &
           posterior_mean_index, trim(adjustl(get_copy_meta_data(seq,posterior_mean_index)))
      call error_handler(E_MSG,'obs_diag',msgstring,source,revision,revdate) 
   endif

   if (prior_spread_index > 0 ) then
      write(msgstring,'(''prior spread     index '',i2,'' metadata '',a)') &
           prior_spread_index, trim(adjustl(get_copy_meta_data(seq,prior_spread_index)))
      call error_handler(E_MSG,'obs_diag',msgstring,source,revision,revdate)
   endif

   if (posterior_spread_index > 0 ) then
      write(msgstring,'(''posterior spread index '',i2,'' metadata '',a)') &
           posterior_spread_index, trim(adjustl(get_copy_meta_data(seq,posterior_spread_index)))
      call error_handler(E_MSG,'obs_diag',msgstring,source,revision,revdate)
   endif

   if (qc_index > 0 ) then
      write(msgstring,'(''Quality Control      index '',i2,'' metadata '',a)') &
           qc_index,      trim(adjustl(get_qc_meta_data(seq,     qc_index)))
      call error_handler(E_MSG,'obs_diag',msgstring,source,revision,revdate)
   endif

   if (dart_qc_index > 0 ) then
      write(msgstring,'(''DART quality control index '',i2,'' metadata '',a)') &
           dart_qc_index, trim(adjustl(get_qc_meta_data(seq,dart_qc_index)))
      call error_handler(E_MSG,'obs_diag',msgstring,source,revision,revdate)
   endif

   end Subroutine SetIndices



   Function Rmidpoints2edges(levels, midpoints)
   ! determine layer midpoints ... ascending vs. descending complications
   integer :: Rmidpoints2edges
   real(r8), dimension(:), intent(inout) :: levels
   real(r8), dimension(:), intent(  out) :: midpoints

   logical :: increasing
   integer :: n, i
   integer :: dummyindxs(1)
   integer, dimension(SIZE(levels)) :: indxs
   logical, dimension(SIZE(levels)) :: logicarray

   Rmidpoints2edges = SIZE(levels,1)
   if (Rmidpoints2edges > MaxLevels) then
      write(msgstring,*)Rmidpoints2edges,' is too many levels - max is ',Maxlevels,' (Maxlevels)'
      call error_handler(E_ERR,'obs_diag:Rmidpoints2edges', msgstring,source,revision,revdate)
   endif

   ! find length of useful portion of levels array
   indxs      = (/ (i,i=1,SIZE(levels)) /) 
   logicarray = levels > MISSING_R8
   dummyindxs = maxloc(indxs, logicarray)
   n          = dummyindxs(1)

   Rmidpoints2edges = n

   ! single-level case works for pressures, must change sign for heights
   if ( n == 1 ) then
      midpoints(1) = levels(1) + 0.1 * levels(1) 
      midpoints(2) = levels(1) - 0.1 * levels(1) 
      return
   endif

   ! ensure monotonicity - up or down - 
   if (levels(1) < levels(n) ) then 
      increasing  = .true.
      levels(1:n) = sort(levels(1:n))
   else
      increasing     = .false.
      midpoints(1:n) = sort(levels(1:n))
      levels(1:n)    = midpoints(n:1:-1)
   endif

   midpoints(  1) = levels(1) - (levels(2) - levels(  1))/2.0_r8
   midpoints(n+1) = levels(n) + (levels(n) - levels(n-1))/2.0_r8
   do i = 2,n
     midpoints(i) = levels(i-1) + (levels(i) - levels(i-1))/2.0_r8
   enddo

   if ( increasing ) then ! must be heights
      midpoints(  1) = max(   0.0_r8, midpoints(  1))
   else                   ! must be pressures
      midpoints(  1) = min(1025.0_r8, midpoints(  1))
      midpoints(n+1) = max(   0.0_r8, midpoints(n+1))
   endif

   end Function Rmidpoints2edges




   Function Imidpoints2edges(levels)
   ! determine layer midpoints ... sorta.
   ! Mostly for compatibility with FindRMidpoints.

   integer :: Imidpoints2edges
   integer, dimension(:), intent(inout) :: levels

   logical :: increasing
   integer :: n, i
   integer :: dummyindxs(1)
   integer,  dimension(SIZE(levels)) :: indxs
   logical,  dimension(SIZE(levels)) :: logicarray
   real(r8), dimension(SIZE(levels)) :: midpoints ! used for the 'sort' routine

   Imidpoints2edges = SIZE(levels,1)
   if (Imidpoints2edges > MaxLevels) then
      write(msgstring,*)Imidpoints2edges,' is too many levels - max is ',Maxlevels,' (Maxlevels)'
      call error_handler(E_ERR,'obs_diag:Imidpoints2edges', msgstring,source,revision,revdate)
   endif

   midpoints = levels

   ! find length of useful portion of levels array
   indxs      = (/ (i,i=1,SIZE(levels)) /) 
   logicarray = levels > nint(MISSING_R8)
   dummyindxs = maxloc(indxs, logicarray)
   n          = dummyindxs(1)
   Imidpoints2edges = n

   ! single-level case 
   if ( n == 1 ) then
      midpoints(1) = levels(1)
      midpoints(2) = levels(1)
      return
   endif

   ! ensure monotonicity - up or down - 
   if (levels(1) < levels(n) ) then 
      increasing  = .true.
      levels(1:n) = nint(sort(midpoints(1:n)))
   else
      increasing     = .false.
      midpoints(1:n) = nint(sort(midpoints(1:n)))
      levels(1:n)    = midpoints(n:1:-1)
   endif

   end Function Imidpoints2edges


 
   Function ParseLevel(obslocation, obslevel, flavor, whichvert )

   type(location_type),   intent(in   ) :: obslocation
   real(r8),              intent(inout) :: obslevel
   integer,               intent(in   ) :: flavor
   integer, dimension(:), intent(inout) :: whichvert
   integer :: ParseLevel
   character(len=8) :: bob

   if(vert_is_surface(obslocation)) then
      obslevel = 1
      bob = 'surface'
      ParseLevel        = 1
      whichvert(flavor) = VERTISSURFACE

   elseif(vert_is_level(obslocation)) then
      ! at one point, negative (model) levels were used to indicate
      ! surface observations - ancient history ...
      bob = 'level'
      if (obslevel < 1.0_r8 ) obslevel = 1.0_r8  ! TJH negative model level
      ParseLevel        = ClosestLevel(obslevel, VERTISLEVEL)
      whichvert(flavor) = VERTISLEVEL

   elseif(vert_is_pressure(obslocation)) then
      bob = 'pressure'
      ParseLevel         = ClosestLevel(obslevel, VERTISPRESSURE)
      whichvert(flavor) = VERTISPRESSURE

   elseif(vert_is_height(obslocation)) then
      bob = 'height'
      ParseLevel         = ClosestLevel(obslevel, VERTISHEIGHT)
      whichvert(flavor) = VERTISHEIGHT
   else
      call error_handler(E_ERR,'obs_diag','Vertical coordinate not recognized', &
           source,revision,revdate)
   endif

   if ( 1 == 2 ) then ! TJH debug statement
      write(*,'(i6,'' obslevel '',f10.4,'', level index '',i3,1x,(a8))') &
         obsindex, obslevel, ParseLevel, bob
   endif

   end Function ParseLevel



   Function ClosestLevel(level_in, levind)

   ! The levels, intervals are ordered  surface == 1
   !
   ! nlev, levels and levels_int are scoped in this unit.
   !
   real(r8), intent(in) :: level_in         ! target level
   integer,  intent(in) :: levind           ! indicator for surface, height, etc.
   integer              :: ClosestLevel

   integer :: a(1)

   real(r8), dimension(MaxLevels) :: dx
   real(r8) :: surface, top

   if (levind == VERTISHEIGHT) then   ! we have heights levels

      surface = hlevel_edges(1)
      top     = hlevel_edges(Nhlevels+1)

      if      ( level_in < surface ) then
         ClosestLevel = -1
      else if ( level_in > top     ) then
         ClosestLevel = 100 + Nhlevels
      else
         dx(1:Nhlevels) = abs(level_in - hlevel(1:Nhlevels))
         a              = minloc(dx(1:Nhlevels))
         ClosestLevel   = a(1)
      endif

   else if (levind == VERTISLEVEL) then   ! we have model levels

      surface = mlevel(1)
      top     = mlevel(Nmlevels)

      if (      level_in < surface ) then
         write(*,*)'obs index ',obsindex,' below surface.'
         ClosestLevel = -1
      else if ( level_in > top     ) then
         ClosestLevel = 100 + Nmlevels
      else
         ClosestLevel = level_in
      endif

   else if (levind == VERTISUNDEF) then

         ClosestLevel = VERTISUNDEF    ! VERTISUNDEF == -2, a good bad value

   else  ! we have pressure levels (or surface obs ... also in pressure units)

      surface = plevel_edges(1)
      top     = plevel_edges(Nplevels+1)

      if (      level_in > surface ) then 
         ClosestLevel = -1
      else if ( level_in < top     ) then
         ClosestLevel = 100 + Nplevels
      else
         dx(1:Nplevels) = abs(level_in - plevel(1:Nplevels))
         a              = minloc(dx(1:Nplevels))
         ClosestLevel   = a(1)
      endif

   endif

   end Function ClosestLevel



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



   Function InRegion( lon, lat, lon1, lon2, lat1, lat2 )
   ! InRegion handles regions that 'wrap' in longitude
   ! if the easternmost longitude of the region is > 360.0
   ! For example; lon1 == 320, lon = 10, lon2 == 380 --> .true.
   ! 
   logical :: InRegion
   real(r8), intent(in) :: lon, lat, lon1, lon2, lat1, lat2
   real(r8) :: mylon

   InRegion = .false.

   if( (lon >= lon1) .and. (lon <= lon2) .and. &
       (lat >= lat1) .and. (lat <= lat2) ) InRegion = .true.

   if( lon2 > 360.0_r8 ) then ! checking the wraparound case
      mylon = lon + 360.0_r8
      if( (mylon >= lon1) .and. (mylon <= lon2) .and. &
          (  lat >= lat1) .and. (  lat <= lat2) ) then
         InRegion = .true.
      endif
   endif

   end Function InRegion



   Function CheckMate(flavor1, flavor2, obsloc1, obsloc2, flavor )

   ! This routine ensures that the U,V components of wind
   ! are from the same observation location so we can convert
   ! them to wind speed. I suppose I _could_ also check the time,
   ! but a mismatch there is supremely unlikely.

   integer,             intent(in)  :: flavor1, flavor2
   type(location_type), intent(in)  :: obsloc1, obsloc2
   integer,             intent(out) :: flavor
   integer                         :: CheckMate

   character(len=stringlength) :: str1, str2

   integer :: mykind1, mykind2
   integer :: indx1,indx2

   CheckMate = -1 ! Assume no match ... till proven otherwise

   mykind1 = get_obs_kind_var_type(flavor1)
   mykind2 = get_obs_kind_var_type(flavor2)

   ! flavor 1 has to be either U or V, flavor 2 has to be the complement
   if ( .not.((mykind2 == KIND_U_WIND_COMPONENT .and. &
               mykind1 == KIND_V_WIND_COMPONENT) .or. &
              (mykind1 == KIND_U_WIND_COMPONENT .and. &
               mykind2 == KIND_V_WIND_COMPONENT)) ) then
      write(*,*) 'flavors not complementary ...',flavor1, flavor2
      flavor = -flavor1
      return
   endif

   if ( obsloc1 /= obsloc2 ) then
      if ( print_mismatched_locs ) then
         write(*,*) 'locations do not match ...'
         call write_location(6,obsloc1,'FORMATTED')
         call write_location(6,obsloc2,'FORMATTED')
         flavor = -flavor1
      endif
      return
   endif

   ! must be matching wind components ... 
   ! must figure out what the obs_kind_name should be ...
   ! There are only two viable wind component strings:
   ! '_?_WIND_COMPONENT' and '_?_10_METER_WIND_SPEED'

   str1  = get_obs_kind_name(flavor1)
   str2  = get_obs_kind_name(flavor2)
   indx1 = index(str1,'_WIND_COMPONENT') - 3
   indx2 = index(str1, '_10_METER_WIND') - 3

   if ( (indx1 < 1) .and. (indx2 < 1) ) then
      write(*,*)'IMPOSSIBLE ... to be CheckMated. har har'
      flavor = -flavor2
   endif

   if (indx1 > 0) then ! must be _?_WIND_COMPONENT
      str3 = str1(1:indx1)//'_WIND_VELOCITY'
   else                ! must be _?_10_METER_WIND_SPEED
      str3 = str1(1:indx2)//'_WIND_VELOCITY'
      indx1 = indx2
   endif

   if (index(str1, str2(1:indx1)) < 1) then
      msgstring = 'observation types not compatible.'
      call error_handler(E_MSG,'obs_diag',msgstring,source,revision,revdate)
   endif

   MyKind : do ivar = 1,num_obs_kinds
      indx1 = index(str3,my_obs_kind_names(ivar))
      if (indx1 > 0) then
         flavor = ivar
         exit MyKind
      endif
   enddo MyKind

   CheckMate = 0

   end Function CheckMate



   Subroutine Bin4D(iqc, iepoch, ilevel, iregion, flavor, &
                obsval,  obserrvar,  prmean,  prsprd,  pomean,  posprd, &
                  uobs, uobserrvar, uprmean, uprsprd, upomean, uposprd)
   !----------------------------------------------------------------------
   ! The 'guess' and 'analysis' structures are globally scoped.
   ! This function simply accumulates the appropriate sums. 
   ! The normalization occurrs after all the data has been read, naturally.
   !
   ! Wind measurements are bivariate - so we are collapsing them to 
   ! scalar speed. The optional arguments specify the U components
   ! while the mandatory arguments specify the V components.
   ! Its an 'all-or-nothing' optional argument situation.
   !----------------------------------------------------------------------

   ! ... spread is computed via sqrt(ensemble_spread**2 + observation_error**2).  
   ! however, dart stores the variance, not the error, so we do not need
   ! to square it here.
   ! If you are verifying the ensemble against imperfect (real) observations, 
   ! it is necessary to account for the observation error when computing the 
   ! spread.  Since the observation error is not included as output from 
   ! obs_diag, it is not possible to compute this quantity now.

   integer,  intent(in) :: iqc, iepoch, ilevel, iregion, flavor
   real(r8), intent(in)           :: obsval,  obserrvar,  prmean,  prsprd,  pomean,  posprd
   real(r8), intent(in), optional ::   uobs, uobserrvar, uprmean, uprsprd, upomean, uposprd

   real(r8) :: priorsqerr      ! PRIOR     Squared Error
   real(r8) :: priorbias       ! PRIOR     simple bias
   real(r8) :: priorspred      ! PRIOR     (spread,variance)
   real(r8) :: priorspredplus  ! PRIOR     (spread,variance**)
   real(r8) :: postsqerr       ! POSTERIOR Squared Error
   real(r8) :: postbias        ! POSTERIOR simple bias
   real(r8) :: postspred       ! POSTERIOR (spread,variance)
   real(r8) :: postspredplus   ! POSTERIOR (spread,variance**)
   logical, dimension(6) :: optionals

   optionals = (/ present(uobs), present(uobserrvar), present(uprmean), &
                  present(uprsprd), present(upomean), present(uposprd) /)

   if ( all(optionals) ) then

      priorsqerr     = (prmean - obsval)**2 + (uprmean - uobs)**2
      postsqerr      = (pomean - obsval)**2 + (upomean - uobs)**2
      priorbias      = (prmean - obsval)    + (uprmean - uobs) 
      postbias       = (pomean - obsval)    + (upomean - uobs)
      priorspred     = prsprd**2 + uprsprd**2
      postspred      = posprd**2 + uposprd**2
      priorspredplus = prsprd**2 + obserrvar + uprsprd**2 + uobserrvar
      postspredplus  = posprd**2 + obserrvar + uposprd**2 + uobserrvar

   else if ( any(optionals) ) then
      call error_handler(E_ERR,'Bin4D','wrong number of optional arguments', &
                         source,revision,revdate)
   else
      priorsqerr     = (prmean - obsval)**2
      postsqerr      = (pomean - obsval)**2
      priorbias      =  prmean - obsval
      postbias       =  pomean - obsval
      priorspred     = prsprd**2
      postspred      = posprd**2
      priorspredplus = prsprd**2 + obserrvar
      postspredplus  = posprd**2 + obserrvar
   endif

   !----------------------------------------------------------------------
   ! Track the number of possible observations
   !----------------------------------------------------------------------

   call IPE(guess%Nposs(iepoch,ilevel,iregion,flavor), 1)
   call IPE(analy%Nposs(iepoch,ilevel,iregion,flavor), 1)

   if ( iqc > QC_MAX_PRIOR ) then  ! prior and posterior failed

      call IPE(guess%NbadQC(iepoch,ilevel,iregion,flavor),      1    )
      call IPE(analy%NbadQC(iepoch,ilevel,iregion,flavor),      1    )

   else if ( iqc > QC_MAX_POSTERIOR ) then

      ! Then at least the prior (A.K.A. guess) is good
      call IPE(guess%Nused(    iepoch,ilevel,iregion,flavor),      1    )
      call RPE(guess%rmse(     iepoch,ilevel,iregion,flavor), priorsqerr)
      call RPE(guess%bias(     iepoch,ilevel,iregion,flavor), priorbias )
      call RPE(guess%spread(   iepoch,ilevel,iregion,flavor), priorspred)
      call RPE(guess%totspread(iepoch,ilevel,iregion,flavor), priorspredplus)

      ! However, the posterior is bad
      call IPE(analy%NbadQC(iepoch,ilevel,iregion,flavor),      1    )

   else

      ! The prior is good
      call IPE(guess%Nused(    iepoch,ilevel,iregion,flavor),      1    ) 
      call RPE(guess%rmse(     iepoch,ilevel,iregion,flavor), priorsqerr)
      call RPE(guess%bias(     iepoch,ilevel,iregion,flavor), priorbias )
      call RPE(guess%spread(   iepoch,ilevel,iregion,flavor), priorspred)
      call RPE(guess%totspread(iepoch,ilevel,iregion,flavor), priorspredplus)

      ! The posterior is good
      call IPE(analy%Nused(    iepoch,ilevel,iregion,flavor),      1    )
      call RPE(analy%rmse(     iepoch,ilevel,iregion,flavor), postsqerr )
      call RPE(analy%bias(     iepoch,ilevel,iregion,flavor), postbias  )
      call RPE(analy%spread(   iepoch,ilevel,iregion,flavor), postspred )
      call RPE(analy%totspread(iepoch,ilevel,iregion,flavor), postspredplus)

   endif

   end Subroutine Bin4D



   Subroutine Bin3D(iqc, ilevel, iregion, flavor, &
                obsval,  obserrvar,  prmean,  prsprd,  pomean,  posprd, &
                  uobs, uobserrvar, uprmean, uprsprd, upomean, uposprd  )
   !----------------------------------------------------------------------
   ! The 'guess' and 'analysis' structures are globally scoped.
   ! This function simply accumulates the appropriate sums. 
   ! The normalization occurrs after all the data has been read, naturally.
   !
   ! Wind measurements are bivariate - so we are collapsing them to 
   ! scalar speed. The optional arguments specify the U components
   ! while the mandatory arguments specify the V components.
   ! Its an 'all-or-nothing' optional argument situation.
   !----------------------------------------------------------------------
   integer,  intent(in)           :: iqc, ilevel, iregion, flavor
   real(r8), intent(in)           :: obsval,  obserrvar,  prmean,  prsprd,  pomean,  posprd
   real(r8), intent(in), optional ::   uobs, uobserrvar, uprmean, uprsprd, upomean, uposprd

   real(r8) :: obsspeed       ! PRIOR     Squared Error
   real(r8) :: priorsqerr     ! PRIOR     Squared Error
   real(r8) :: priorbias      ! PRIOR     simple bias
   real(r8) :: priorspred     ! PRIOR     (spread,variance)
   real(r8) :: priorspredplus ! PRIOR     (spread,variance**)
   real(r8) :: postsqerr      ! POSTERIOR Squared Error
   real(r8) :: postbias       ! POSTERIOR simple bias
   real(r8) :: postspred      ! POSTERIOR (spread,variance)
   real(r8) :: postspredplus  ! POSTERIOR (spread,variance**)
   logical, dimension(6) :: optionals

   optionals = (/ present(uobs), present(uobserrvar), present(uprmean), &
                  present(uprsprd), present(upomean), present(uposprd) /)

   if ( all(optionals) ) then
      priorsqerr     = (prmean - obsval)**2 + (uprmean - uobs)**2
      postsqerr      = (pomean - obsval)**2 + (upomean - uobs)**2

      ! This calculation is the bias in the wind vector  
      priorbias      = (prmean - obsval)    + (uprmean - uobs) 
      postbias       = (pomean - obsval)    + (upomean - uobs)

      ! This calculation is the bias in the wind speed
      obsspeed       = sqrt(uobs**2 + obsval**2)
      priorbias      = obsspeed - sqrt(prmean**2 + uprmean**2) 
      postbias       = obsspeed - sqrt(pomean**2 + upomean**2) 

      priorspred     = prsprd**2 + uprsprd**2
      postspred      = posprd**2 + uposprd**2
      priorspredplus = prsprd**2 + obserrvar + uprsprd**2 + uobserrvar
      postspredplus  = posprd**2 + obserrvar + uposprd**2 + uobserrvar

   else if ( any(optionals) ) then
      call error_handler(E_ERR,'Bin3D','wrong number of optional arguments', &
                         source,revision,revdate)
   else
      priorsqerr     = (prmean - obsval)**2
      postsqerr      = (pomean - obsval)**2
      priorbias      =  prmean - obsval
      postbias       =  pomean - obsval
      priorspred     = prsprd**2
      postspred      = posprd**2
      priorspredplus = prsprd**2 + obserrvar
      postspredplus  = posprd**2 + obserrvar

   endif

   !----------------------------------------------------------------------
   ! Track the number of possible observations
   !----------------------------------------------------------------------
   call IPE(guessAVG%Nposs(ilevel,iregion,flavor), 1)
   call IPE(analyAVG%Nposs(ilevel,iregion,flavor), 1)

   if ( iqc > QC_MAX_PRIOR ) then  ! prior and posterior failed

      call IPE(guessAVG%NbadQC(ilevel,iregion,flavor),      1    )
      call IPE(analyAVG%NbadQC(ilevel,iregion,flavor),      1    )

   else if ( iqc > QC_MAX_POSTERIOR ) then

      ! Then at least the prior (A.K.A. guess) is good
      call IPE(guessAVG%Nused(    ilevel,iregion,flavor),      1    )
      call RPE(guessAVG%rmse(     ilevel,iregion,flavor), priorsqerr)
      call RPE(guessAVG%bias(     ilevel,iregion,flavor), priorbias )
      call RPE(guessAVG%spread(   ilevel,iregion,flavor), priorspred)
      call RPE(guessAVG%totspread(ilevel,iregion,flavor), priorspredplus)

      ! However, the posterior is bad
      call IPE(analyAVG%NbadQC(ilevel,iregion,flavor),      1    )

   else

      ! The prior is good
      call IPE(guessAVG%Nused(    ilevel,iregion,flavor),      1    ) 
      call RPE(guessAVG%rmse(     ilevel,iregion,flavor), priorsqerr)
      call RPE(guessAVG%bias(     ilevel,iregion,flavor), priorbias )
      call RPE(guessAVG%spread(   ilevel,iregion,flavor), priorspred)
      call RPE(guessAVG%totspread(ilevel,iregion,flavor), priorspredplus)

      ! The posterior is good
      call IPE(analyAVG%Nused(    ilevel,iregion,flavor),      1    )
      call RPE(analyAVG%rmse(     ilevel,iregion,flavor), postsqerr )
      call RPE(analyAVG%bias(     ilevel,iregion,flavor), postbias  )
      call RPE(analyAVG%spread(   ilevel,iregion,flavor), postspred)
      call RPE(analyAVG%totspread(ilevel,iregion,flavor), postspredplus)

   endif

   end Subroutine Bin3D



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



   Subroutine WriteNetCDF(fname)
   character(len=129), intent(in) :: fname

   integer :: ncid, i, indx1
   integer ::  RegionDimID,  RegionVarID
   integer ::  MlevelDimID,  MlevelVarID
   integer ::  PlevelDimID,  PlevelVarID
   integer ::  HlevelDimID,  HlevelVarID
   integer ::  SlevelDimID
   integer ::    TimeDimID,    TimeVarID
   integer ::    CopyDimID,    CopyVarID,  CopyMetaVarID
   integer ::   TypesDimID,   TypesVarID, TypesMetaVarID
   integer :: PlevIntDimID, PlevIntVarID
   integer :: HlevIntDimID, HlevIntVarID
   integer ::  BoundsDimID,  BoundsVarID  
   integer ::  StringDimID

   integer :: TimeBoundsVarID, RegionNamesVarID

   character(len=8)      :: crdate      ! needed by F90 DATE_AND_TIME intrinsic
   character(len=10)     :: crtime      ! needed by F90 DATE_AND_TIME intrinsic
   character(len=5)      :: crzone      ! needed by F90 DATE_AND_TIME intrinsic
   integer, dimension(8) :: values      ! needed by F90 DATE_AND_TIME intrinsic

   if(.not. byteSizesOK()) then
       call error_handler(E_ERR,'WriteNetCDF', &
      'Compiler does not support required kinds of variables.',source,revision,revdate)
   endif

   call nc_check(nf90_create(path = trim(fname), cmode = nf90_share, &
            ncid = ncid), 'obs_diag:WriteNetCDF', 'create '//trim(fname))

   write(msgstring,*)trim(ncName), ' is fortran unit ',ncid
   call error_handler(E_MSG,'WriteNetCDF',msgstring,source,revision,revdate)

   !----------------------------------------------------------------------------
   ! Write Global Attributes 
   !----------------------------------------------------------------------------

   call DATE_AND_TIME(crdate,crtime,crzone,values)
   write(msgstring,'(''YYYY MM DD HH MM SS = '',i4,5(1x,i2.2))') &
                  values(1), values(2), values(3), values(5), values(6), values(7)
   call nc_check(nf90_put_att(ncid, NF90_GLOBAL, "creation_date", trim(msgstring) ), &
              'WriteNetCDF', 'put_att creation_date '//trim(fname))

!  call nc_check(nf90_put_att(ncid, NF90_GLOBAL, "title", global_meta_data), &
!             'WriteNetCDF', 'put_att title '//trim(fname))
   call nc_check(nf90_put_att(ncid, NF90_GLOBAL, "obs_diag_source", source ), &
              'WriteNetCDF', 'put_att obs_diag_source '//trim(fname))
   call nc_check(nf90_put_att(ncid, NF90_GLOBAL, "obs_diag_revision", revision ), &
              'WriteNetCDF', 'put_att obs_diag_revision '//trim(fname))
   call nc_check(nf90_put_att(ncid, NF90_GLOBAL, "obs_diag_revdate", revdate ), &
              'WriteNetCDF', 'put_att obs_diag_revdate '//trim(fname))
   call nc_check(nf90_put_att(ncid, NF90_GLOBAL, "bias_convention", "forecast - observation" ), &
              'WriteNetCDF', 'put_att bias '//trim(fname))

   ! write all namelist quantities 

   call nc_check(nf90_put_att(ncid, NF90_GLOBAL, "first_bin_center", first_bin_center ), &
              'WriteNetCDF', 'put_att first_bin_center '//trim(fname))
   call nc_check(nf90_put_att(ncid, NF90_GLOBAL, "last_bin_center", last_bin_center ), &
              'WriteNetCDF', 'put_att last_bin_center '//trim(fname))
   call nc_check(nf90_put_att(ncid, NF90_GLOBAL, "bin_separation", bin_separation ), &
              'WriteNetCDF', 'put_att bin_separation '//trim(fname))
   call nc_check(nf90_put_att(ncid, NF90_GLOBAL, "bin_width", bin_width ), &
              'WriteNetCDF', 'put_att bin_width '//trim(fname))
   call nc_check(nf90_put_att(ncid, NF90_GLOBAL, "time_to_skip", time_to_skip ), &
              'WriteNetCDF', 'put_att time_to_skip '//trim(fname))
   call nc_check(nf90_put_att(ncid, NF90_GLOBAL, "max_num_bins", max_num_bins ), &
              'WriteNetCDF', 'put_att max_num_bins '//trim(fname))
   call nc_check(nf90_put_att(ncid, NF90_GLOBAL, "rat_cri", rat_cri ), &
              'WriteNetCDF', 'put_att rat_cri '//trim(fname))
   call nc_check(nf90_put_att(ncid, NF90_GLOBAL, "input_qc_threshold", input_qc_threshold ), &
              'WriteNetCDF', 'put_att input_qc_threshold '//trim(fname))
   call nc_check(nf90_put_att(ncid, NF90_GLOBAL, "lonlim1", lonlim1(1:Nregions) ), &
              'WriteNetCDF', 'put_att lonlim1 '//trim(fname))
   call nc_check(nf90_put_att(ncid, NF90_GLOBAL, "lonlim2", lonlim2(1:Nregions) ), &
              'WriteNetCDF', 'put_att lonlim2 '//trim(fname))
   call nc_check(nf90_put_att(ncid, NF90_GLOBAL, "latlim1", latlim1(1:Nregions) ), &
              'WriteNetCDF', 'put_att latlim1 '//trim(fname))
   call nc_check(nf90_put_att(ncid, NF90_GLOBAL, "latlim2", latlim2(1:Nregions) ), &
              'WriteNetCDF', 'put_att latlim2 '//trim(fname))

   ! write all observation sequence files used
   FILEloop : do i = 1,SIZE(obs_seq_filenames)

     indx1 = index(obs_seq_filenames(i),'null')

     if (indx1 > 0) exit FILEloop

     write(msgstring,'(''obs_seq_file_'',i3.3)')i
     call nc_check(nf90_put_att(ncid, NF90_GLOBAL, &
            trim(msgstring), trim(obs_seq_filenames(i)) ), &
            'WriteNetCDF', 'region_names:obs_kinds')

   enddo FILEloop

   call nc_check(nf90_put_att(ncid, NF90_GLOBAL, "NumIdentityObs", Nidentity ), &
              'WriteNetCDF', 'put_att identity '//trim(fname))

   ! write all 'known' observation types
   call nc_check(nf90_put_att(ncid, NF90_GLOBAL, "comment", &
              'All known observation types follow. &
              &Also see ObservationTypes variable.' ), &
              'WriteNetCDF', 'put_att latlim2 '//trim(fname))
   do ivar = 1,max_obs_kinds
     call nc_check(nf90_put_att(ncid, NF90_GLOBAL, &
            trim(adjustl(my_obs_kind_names(ivar))), ivar ), &
            'WriteNetCDF', 'region_names:obs_kinds')
   enddo

   !----------------------------------------------------------------------------
   ! Define the dimensions
   !----------------------------------------------------------------------------

   call nc_check(nf90_def_dim(ncid=ncid, &
              name="copy", len = Ncopies,            dimid = CopyDimID), &
              'WriteNetCDF', 'copy:def_dim '//trim(fname))

   call nc_check(nf90_def_dim(ncid=ncid, &
              name="obstypes", len = max_obs_kinds,  dimid = TypesDimID), &
              'WriteNetCDF', 'types:def_dim '//trim(fname))

   call nc_check(nf90_def_dim(ncid=ncid, &
              name="region", len = Nregions,         dimid = RegionDimID), &
              'WriteNetCDF', 'region:def_dim '//trim(fname))

   call nc_check(nf90_def_dim(ncid=ncid, &
              name="surface", len = 1,               dimid = SlevelDimID), &
              'WriteNetCDF', 'slevel:def_dim '//trim(fname))

   call nc_check(nf90_def_dim(ncid=ncid, &
              name="mlevel", len = Nmlevels,         dimid = MlevelDimID), &
              'WriteNetCDF', 'mlevel:def_dim '//trim(fname))

   call nc_check(nf90_def_dim(ncid=ncid, &
              name="plevel", len = Nplevels,         dimid = PlevelDimID), &
              'WriteNetCDF', 'plevel:def_dim '//trim(fname))
   call nc_check(nf90_def_dim(ncid=ncid, &
              name="plevel_edges", len = Nplevels+1, dimid = PlevIntDimID), &
              'WriteNetCDF', 'plevel_edges:def_dim '//trim(fname))

   call nc_check(nf90_def_dim(ncid=ncid, &
              name="hlevel", len = Nhlevels,         dimid = HlevelDimID), &
              'WriteNetCDF', 'hlevel:def_dim '//trim(fname))
   call nc_check(nf90_def_dim(ncid=ncid, &
              name="hlevel_edges", len = Nhlevels+1, dimid = HlevIntDimID), &
              'WriteNetCDF', 'hlevel_edges:def_dim '//trim(fname))

   call nc_check(nf90_def_dim(ncid=ncid, &
              name="time",   len = NF90_UNLIMITED,   dimid = TimeDimID), &
              'WriteNetCDF', 'time:def_dim '//trim(fname))
   call nc_check(nf90_def_dim(ncid=ncid, &
              name="bounds",   len = 2,  dimid = BoundsDimID), &
              'WriteNetCDF', 'bounds:def_dim '//trim(fname))

   call nc_check(nf90_def_dim(ncid=ncid, &
              name="stringlength", len = stringlength, dimid = StringDimID), &
              'WriteNetCDF', 'stringlength:def_dim '//trim(fname))

   ! Define the types of derived quantities - aka - 'copies'

   call nc_check(nf90_def_var(ncid=ncid, name="copy", xtype=nf90_int, &
             dimids=CopyDimID, varid=CopyVarID), &
             'WriteNetCDF', 'copy:def_var')
   call nc_check(nf90_put_att(ncid, CopyVarID, "explanation", 'see CopyMetaData'), &
             'WriteNetCDF', 'copy:explanation')

   ! Define the observation types - needed to be a coordinate variable

   call nc_check(nf90_def_var(ncid=ncid, name="obstypes", xtype=nf90_int, &
             dimids=TypesDimID, varid=TypesVarID), &
             'WriteNetCDF', 'types:def_var')
   call nc_check(nf90_put_att(ncid, TypesVarID, "explanation", 'see ObservationTypes'), &
             'WriteNetCDF', 'types:explanation')

   ! Define the regions coordinate variable and attributes

   call nc_check(nf90_def_var(ncid=ncid, name="region", xtype=nf90_int, &
             dimids=RegionDimID, varid=RegionVarID), 'WriteNetCDF', 'region:def_var')
   call nc_check(nf90_put_att(ncid, RegionVarID, "long_name", "model region"), &
             'WriteNetCDF', 'region:long_name')
   call nc_check(nf90_put_att(ncid, RegionVarID, "units",     "nondimensional"), &
             'WriteNetCDF', 'region:units')
   call nc_check(nf90_put_att(ncid, RegionVarID, "valid_range", (/1,Nregions/)), &
             'WriteNetCDF', 'region:valid_range')

   ! Define the model level coordinate variable and attributes

   call nc_check(nf90_def_var(ncid=ncid, name="mlevel", xtype=nf90_int, &
             dimids=MlevelDimID, varid=MlevelVarID), 'WriteNetCDF', 'mlevel:def_var')
   call nc_check(nf90_put_att(ncid, MlevelVarID, "long_name", "model level"), &
             'WriteNetCDF', 'mlevel:long_name')
   call nc_check(nf90_put_att(ncid, MlevelVarID, "units",     "nondimensional"), &
             'WriteNetCDF', 'mlevel:units')
   call nc_check(nf90_put_att(ncid, MlevelVarID, "valid_range", (/1,Nmlevels/)), &
             'WriteNetCDF', 'mlevel:valid_range')

   ! Define the pressure level coordinate variable and attributes

   call nc_check(nf90_def_var(ncid=ncid, name="plevel", xtype=nf90_real, &
             dimids=PlevelDimID, varid=PlevelVarID), 'WriteNetCDF', 'plevel:def_var')
   call nc_check(nf90_put_att(ncid, PlevelVarID, "long_name", "pressure bin midpoints"), &
             'WriteNetCDF', 'plevel:long_name')
   call nc_check(nf90_put_att(ncid, PlevelVarID, "units",     "hPa"), &
             'WriteNetCDF', 'plevel:units')
   call nc_check(nf90_put_att(ncid, PlevelVarID, "valid_range", &
      (/ minval(plevel(1:Nplevels)), maxval(plevel(1:Nplevels)) /)), &
             'WriteNetCDF', 'plevel:valid_range')

   ! Define the pressure level interface coordinate variable and attributes

   call nc_check(nf90_def_var(ncid=ncid, name="plevel_edges", xtype=nf90_real, &
             dimids=PlevIntDimID, varid=PlevIntVarID), 'WriteNetCDF', 'plevel_edges:def_var')
   call nc_check(nf90_put_att(ncid, PlevIntVarID, "long_name", "pressure bin edges"), &
             'WriteNetCDF', 'plevel_edges:long_name')
   call nc_check(nf90_put_att(ncid, PlevIntVarID, "units",     "hPa"), &
             'WriteNetCDF', 'plevel_edges:units')
   call nc_check(nf90_put_att(ncid, PlevIntVarID, "valid_range", &
      (/ minval(plevel_edges(1:Nplevels+1)), maxval(plevel_edges(1:Nplevels+1)) /)), &
             'WriteNetCDF', 'plevel_edges:valid_range')

   ! Define the height level coordinate variable and attributes

   call nc_check(nf90_def_var(ncid=ncid, name="hlevel", xtype=nf90_real, &
           dimids=HlevelDimID, varid=HlevelVarID), 'WriteNetCDF', 'hlevel:def_var')
   call nc_check(nf90_put_att(ncid, HlevelVarID, "long_name", "height bin midpoints"), &
             'WriteNetCDF', 'hlevel:long_name')
   call nc_check(nf90_put_att(ncid, HlevelVarID, "units",     "m"), &
             'WriteNetCDF', 'hlevel:units')
   call nc_check(nf90_put_att(ncid, HlevelVarID, "valid_range", &
      (/ minval(hlevel(1:Nhlevels)), maxval(hlevel(1:Nhlevels)) /)), &
             'WriteNetCDF', 'hlevel:valid_range')

   ! Define the height level interface coordinate variable and attributes

   call nc_check(nf90_def_var(ncid=ncid, name="hlevel_edges", xtype=nf90_real, &
           dimids=HlevIntDimID, varid=HlevIntVarID), 'WriteNetCDF', 'hlevel_edges:def_var')
   call nc_check(nf90_put_att(ncid, HlevIntVarID, "long_name", "height bin edges"), &
             'WriteNetCDF', 'hlevel_edges:long_name')
   call nc_check(nf90_put_att(ncid, HlevIntVarID, "units",     "m"), &
             'WriteNetCDF', 'hlevel_edges:units')
   call nc_check(nf90_put_att(ncid, HlevIntVarID, "valid_range", &
      (/ minval(hlevel_edges(1:Nhlevels+1)), maxval(hlevel_edges(1:Nhlevels+1)) /)), &
             'WriteNetCDF', 'hlevel_edges:valid_range')

   ! Define "bounds" dimension

   call nc_check(nf90_def_var(ncid=ncid, name="bounds", xtype=nf90_int, &
             dimids=BoundsDimID, varid=BoundsVarID), 'WriteNetCDF', 'bounds:def_var')
   call nc_check(nf90_put_att(ncid, BoundsVarID, "valid_range", (/1,2/)), &
             'WriteNetCDF', 'bounds:valid_range')

   ! Define the time coordinate variable and attributes

   call nc_check(nf90_def_var(ncid=ncid, name="time", xtype=nf90_double, &
             dimids=TimeDimID, varid=TimeVarID), 'WriteNetCDF', 'time:def_var')
   call nc_check(nf90_put_att(ncid, TimeVarID, "standard_name",    "time"), &
             'WriteNetCDF', 'time:standard_name')
   call nc_check(nf90_put_att(ncid, TimeVarID, "long_name", "temporal bin midpoints"), &
             'WriteNetCDF', 'time:long_name')
   call nc_check(nf90_put_att(ncid, TimeVarID, "units",     "days since 1601-1-1"), &
             'WriteNetCDF', 'time:units')
   call nc_check(nf90_put_att(ncid, TimeVarID, "calendar",    "Gregorian"), &
             'WriteNetCDF', 'time:calendar')
   call nc_check(nf90_put_att(ncid, TimeVarID, "axis",    "T"), &
             'WriteNetCDF', 'time:axis')
   call nc_check(nf90_put_att(ncid, TimeVarID, "bounds",    "time_bounds"), &
             'WriteNetCDF', 'time:bounds')
   call nc_check(nf90_put_att(ncid, TimeVarID, "valid_range", &
             (/ epoch_center(1), epoch_center(Nepochs) /)), &
             'WriteNetCDF', 'time:valid_range')

   ! Define the time edges coordinate variable and attributes

   call nc_check(nf90_def_var(ncid=ncid, name="time_bounds", xtype=nf90_double, &
             dimids=(/ BoundsDimID, TimeDimID/), varid=TimeBoundsVarID), &
               'WriteNetCDF', 'time_bounds:def_var')
   call nc_check(nf90_put_att(ncid, TimeBoundsVarID, "long_name", "temporal bin edges"), &
             'WriteNetCDF', 'time_bounds:long_name')
   call nc_check(nf90_put_att(ncid, TimeBoundsVarID, "units",     "days since 1601-1-1"), &
             'WriteNetCDF', 'time_bounds:units')
   call nc_check(nf90_put_att(ncid, TimeBoundsVarID, "calendar",  "Gregorian"), &
             'WriteNetCDF', 'time_bounds:calendar')
   call nc_check(nf90_put_att(ncid, TimeBoundsVarID, "valid_range", &
             (/ epoch_edges(1,1), epoch_edges(2,Nepochs) /)), &
             'WriteNetCDF', 'time_bounds:valid_range')

   ! Define the regions variable ... an unusual coordinate variable IMHO

   call nc_check(nf90_def_var(ncid=ncid, name="region_names", xtype=nf90_char, &
             dimids=(/ StringDimID, RegionDimID /), varid=RegionNamesVarID), &
             'WriteNetCDF', 'region:def_var')
   call nc_check(nf90_put_att(ncid, RegionNamesVarID, "long_name", "region names"), &
             'WriteNetCDF', 'region:long_name')

   call nc_check(nf90_def_var(ncid=ncid, name="CopyMetaData", xtype=nf90_char, &
             dimids=(/ StringDimID, CopyDimID /), varid=CopyMetaVarID), &
             'WriteNetCDF', 'copymeta:def_var')
   call nc_check(nf90_put_att(ncid, CopyMetaVarID, "long_name", "quantity names"), &
             'WriteNetCDF', 'copymeta:long_name')

   call nc_check(nf90_def_var(ncid=ncid, name="ObservationTypes", xtype=nf90_char, &
             dimids=(/ StringDimID, TypesDimID /), varid=TypesMetaVarID), &
             'WriteNetCDF', 'typesmeta:def_var')
   call nc_check(nf90_put_att(ncid, TypesMetaVarID, "long_name", "DART observation types"), &
             'WriteNetCDF', 'typesmeta:long_name')
   call nc_check(nf90_put_att(ncid, TypesMetaVarID, "comment", &
         "table relating integer to observation type string"), &
             'WriteNetCDF', 'typesmeta:comment')

   ! Set nofill mode - supposed to be performance gain
 
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

   call nc_check(nf90_put_var(ncid, TypesMetaVarID, my_obs_kind_names(1:max_obs_kinds)), &
              'WriteNetCDF', 'typesmeta:put_var')

   call nc_check(nf90_put_var(ncid, RegionVarID, (/ (i,i=1,Nregions) /) ), &
              'WriteNetCDF', 'region:put_var')

   call nc_check(nf90_put_var(ncid, MlevelVarID, (/ (i,i=1,Nmlevels) /) ), &
              'WriteNetCDF', 'mlevel:put_var')

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
   !----------------------------------------------------------------------------

   if (verbose) write(*,*)'summary for Priors of time-level-region vars' 
   ierr = WriteTLRV(ncid, guess,    TimeDimID, CopyDimID, RegionDimID)
   if (verbose) write(*,*)'summary for Posteriors of time-level-region vars' 
   ierr = WriteTLRV(ncid, analy,    TimeDimID, CopyDimID, RegionDimID)
   ierr = WriteLRV( ncid, guessAVG,            CopyDimID, RegionDimID)
   ierr = WriteLRV( ncid, analyAVG,            CopyDimID, RegionDimID)

   !----------------------------------------------------------------------------
   ! finish ...
   !----------------------------------------------------------------------------

   call nc_check(nf90_sync( ncid), 'WriteNetCDF', 'sync '//trim(fname))  
   call nc_check(nf90_close(ncid), 'init_diag_output', 'close '//trim(fname))  

   end Subroutine WriteNetCDF


   Function WriteTLRV(ncid, vrbl, TimeDimID, CopyDimID, RegionDimID)
   integer,         intent(in) :: ncid
   type(TLRV_type), intent(in) :: vrbl
   integer,         intent(in) :: TimeDimID, CopyDimID, RegionDimID
   integer :: WriteTLRV

   integer :: nobs, Nlevels, ivar, itime, ilevel, iregion
   character(len=40) :: string1

   integer :: VarID, LevelDimID, oldmode
   real(r4), allocatable, dimension(:,:,:,:) :: chunk

   FLAVORS : do ivar = 1,num_obs_kinds

      nobs = sum(vrbl%Nposs(:,:,:,ivar))

      if (verbose) then
         write(*,'(i4,1x,(a32),1x,i8,1x,'' obs@vert '',i3,f11.3)') ivar, &
          my_obs_kind_names(ivar), nobs, which_vert(ivar), scale_factor(ivar)
      endif

      if (nobs < 1) cycle FLAVORS

      ! determine what kind of levels to use ... models, pressure, height ...

      Nlevels = FindVertical(ncid, ivar, LevelDimID)

      allocate(chunk(Nregions,Nlevels,Ncopies,Nepochs))
      chunk = MISSING_R4

      do itime   = 1,Nepochs
      do ilevel  = 1,Nlevels
      do iregion = 1,Nregions

         chunk(iregion,ilevel, 1,itime) = vrbl%Nposs(    itime,ilevel,iregion,ivar)
         chunk(iregion,ilevel, 2,itime) = vrbl%Nused(    itime,ilevel,iregion,ivar)
         chunk(iregion,ilevel, 3,itime) = vrbl%NbadQC(   itime,ilevel,iregion,ivar)
         chunk(iregion,ilevel, 4,itime) = vrbl%NbadIZ(   itime,ilevel,iregion,ivar)
         chunk(iregion,ilevel, 5,itime) = vrbl%NbadUV(   itime,ilevel,iregion,ivar)
         chunk(iregion,ilevel, 6,itime) = vrbl%NbadLV(   itime,ilevel,iregion,ivar)
         chunk(iregion,ilevel, 7,itime) = vrbl%rmse(     itime,ilevel,iregion,ivar)
         chunk(iregion,ilevel, 8,itime) = vrbl%bias(     itime,ilevel,iregion,ivar)
         chunk(iregion,ilevel, 9,itime) = vrbl%spread(   itime,ilevel,iregion,ivar)
         chunk(iregion,ilevel,10,itime) = vrbl%totspread(itime,ilevel,iregion,ivar)

      enddo
      enddo
      enddo

      call nc_check(nf90_redef(ncid), 'WriteTLRV', 'redef')  

      ! Create netCDF variable name
      
      str1 = my_obs_kind_names(ivar)
      string1 = trim(adjustl(str1))//'_'//trim(adjustl(vrbl%string))

      call nc_check(nf90_def_var(ncid, name=string1, xtype=nf90_real, &
             dimids=(/ RegionDimID, LevelDimID, CopyDimID, TimeDimID /), &
             varid=VarID), 'WriteTLRV', 'region:def_var')
      call nc_check(nf90_put_att(ncid, VarID, '_FillValue',    MISSING_R4), &
              'WriteTLRV','put_att:fillvalue')
      call nc_check(nf90_put_att(ncid, VarID, 'missing_value', MISSING_R4), &
              'WriteTLRV','put_att:missing')

      call nc_check(nf90_set_fill(ncid, NF90_NOFILL, oldmode),  &
              'WriteTLRV', 'set_nofill '//trim(vrbl%string))

      call nc_check(nf90_enddef(ncid), 'WriteTLRV', 'enddef ')

      call nc_check(nf90_put_var(ncid, VarID, chunk ), &
              'WriteTLRV', 'time_bounds:put_var')

      deallocate(chunk)

   enddo FLAVORS

   WriteTLRV = 0

   end Function WriteTLRV



   Function WriteLRV(ncid, vrbl, CopyDimID, RegionDimID)
   integer,         intent(in) :: ncid
   type(LRV_type),  intent(in) :: vrbl
   integer,         intent(in) :: CopyDimID, RegionDimID
   integer :: WriteLRV

   integer :: nobs, Nlevels, ivar, ilevel, iregion
   character(len=40) :: string1

   integer :: VarID, LevelDimID, oldmode
   real(r4), allocatable, dimension(:,:,:) :: chunk

   FLAVORS : do ivar = 1,num_obs_kinds

      nobs = sum(vrbl%Nposs(:,:,ivar))
      if (nobs < 1) cycle FLAVORS

      ! determine what kind of levels to use ... models, pressure, height ...

      Nlevels = FindVertical(ncid, ivar, LevelDimID)

      allocate(chunk(Nregions,Nlevels,Ncopies))
      chunk = MISSING_R4

      do ilevel  = 1,Nlevels
      do iregion = 1,Nregions

         chunk(iregion,ilevel, 1) = vrbl%Nposs(    ilevel,iregion,ivar)
         chunk(iregion,ilevel, 2) = vrbl%Nused(    ilevel,iregion,ivar)
         chunk(iregion,ilevel, 3) = vrbl%NbadQC(   ilevel,iregion,ivar)
         chunk(iregion,ilevel, 4) = vrbl%NbadIZ(   ilevel,iregion,ivar)
         chunk(iregion,ilevel, 5) = vrbl%NbadUV(   ilevel,iregion,ivar)
         chunk(iregion,ilevel, 6) = vrbl%NbadLV(   ilevel,iregion,ivar)
         chunk(iregion,ilevel, 7) = vrbl%rmse(     ilevel,iregion,ivar)
         chunk(iregion,ilevel, 8) = vrbl%bias(     ilevel,iregion,ivar)
         chunk(iregion,ilevel, 9) = vrbl%spread(   ilevel,iregion,ivar)
         chunk(iregion,ilevel,10) = vrbl%totspread(ilevel,iregion,ivar)

      enddo
      enddo

      call nc_check(nf90_redef(ncid), 'WriteLRV', 'redef')  

      ! Create netCDF variable name
      
      str1 = my_obs_kind_names(ivar)
      string1 = trim(adjustl(str1))//'_'//trim(adjustl(vrbl%string))

      call nc_check(nf90_def_var(ncid, name=string1, xtype=nf90_real, &
             dimids=(/ RegionDimID, LevelDimID, CopyDimID /), &
             varid=VarID), 'WriteLRV', 'region:def_var')
      call nc_check(nf90_put_att(ncid, VarID, '_FillValue',    MISSING_R4), &
              'WriteLRV','put_att:fillvalue')
      call nc_check(nf90_put_att(ncid, VarID, 'missing_value', MISSING_R4), &
              'WriteLRV','put_att:missing')

      call nc_check(nf90_set_fill(ncid, NF90_NOFILL, oldmode),  &
              'WriteLRV', 'set_nofill '//trim(vrbl%string))

      call nc_check(nf90_enddef(ncid), 'WriteLRV', 'enddef ')

      call nc_check(nf90_put_var(ncid, VarID, chunk ), &
              'WriteLRV', 'time_bounds:put_var')

      deallocate(chunk)

   enddo FLAVORS

   WriteLRV = 0

   end Function WriteLRV



   Function FindVertical(ncid, flav, dimid)
      integer, intent(in)  :: ncid, flav
      integer, intent(out) :: dimid
      integer              :: FindVertical

      if      ( which_vert(flav) == VERTISSURFACE ) then

         FindVertical = 1
         call nc_check(nf90_inq_dimid(ncid, 'surface', dimid), &
                                       'FindVertical', 'vertissurface')

      else if ( which_vert(flav) == VERTISLEVEL   ) then

         FindVertical = Nmlevels
         call nc_check(nf90_inq_dimid(ncid, 'mlevel', dimid), &
                                       'FindVertical', 'vertislevel')

      else if ( which_vert(flav) == VERTISPRESSURE) then

         FindVertical = Nplevels
         call nc_check(nf90_inq_dimid(ncid, 'plevel', dimid), &
                                       'FindVertical', 'vertispressure')

      else if ( which_vert(flav) == VERTISHEIGHT  ) then

         FindVertical = Nhlevels
         call nc_check(nf90_inq_dimid(ncid, 'hlevel', dimid), &
                                       'FindVertical', 'vertisheight')

      else if ( which_vert(flav) == VERTISUNDEF   ) then
         write(msgstring,*)flav,trim(my_obs_kind_names(flav)),' has undefined vertical.'
         call error_handler(E_MSG,'FindVertical',msgstring,source,revision,revdate)
         call error_handler(E_ERR,'FindVertical','vertical undefined not implemented', &
                    source,revision,revdate)
      else 
         call error_handler(E_ERR,'FindVertical','unknown vertical', &
                    source,revision,revdate)
      endif

   end Function FindVertical



   Function grok_observation_names(my_names)
   !----------------------------------------------------------------------
   ! Define/Append the 'wind velocity' obs_kinds to supplant the list declared
   ! in obs_kind_mod.f90 i.e. if there is a RADIOSONDE_U_WIND_COMPONENT
   ! and a RADIOSONDE_V_WIND_COMPONENT, there must be a RADIOSONDE_WIND_VELOCITY
   ! Replace calls to 'get_obs_kind_name' with variable 'my_obs_kind_names'
   !----------------------------------------------------------------------

   character(len=stringlength), pointer :: my_names(:) ! INTENT OUT, btw
   integer :: grok_observation_names

   integer :: ivar, nwinds
   character(len=stringlength) :: str1, str2, str3
   character(len=stringlength), dimension(2*max_obs_kinds) :: names

   ! Initially, the array of obs_kind_names is exactly 'max_num_obs' in length.
   ! This block finds the U,V wind pairs and searches for pre-existing wind_velocity
   ! equivalents. Depending on the number of unique wind pairs - we can allocate
   ! space, copy the existing names into that array, and append the new unique ones.
   ! easy ...

   integer :: indx1, indx2
   integer :: indx1N,indx2N,indxN

   nwinds = 0

   ! Copy all the known obs kinds to a local list that is SURELY too BIG
   ! as we find wind pairs, we insert the new name at the end of the known 
   ! names.

   do ivar = 1,max_obs_kinds
      names(ivar) = get_obs_kind_name(ivar)
   enddo

   ! Search through the obs_kind_name list for matching U,V components.
   ! The U component always comes before the V component, so we exploit that.
   ! Once we have counted the pairs - we know how far to expand the obs_kind list.

   do ivar = 2,max_obs_kinds

      str1   = names(ivar-1)
      indx1  = index(str1,'_U_WIND_COMPONENT') - 1
      indx1N = len_trim(str1)

      str2   = names(ivar)
      indx2  = index(str2,'_V_WIND_COMPONENT') - 1
      indx2N = len_trim(str2)

   !  write(*,*)'Checking ',ivar, indx1, indx2, trim(adjustl(str2))

      if ((indx1 > 0) .and. (indx2 > 0)) then             ! we know we have u,v wind components

         indxN = index(str1(1:indx1),str2(1:indx2))

      !  write(*,*)' have u,v components at ',ivar,indxN

         if (indxN > 0) then ! we know they are matching kinds
            nwinds = nwinds + 1
            str3   = str1(1:indx2)//'_WIND_VELOCITY'
            names(max_obs_kinds + nwinds) = str3

         !  write(*,*)'Seems like ',str1(1:indx1N),' matches ',str2(1:indx2N)
         !  write(*,*)'results in ',str3
         endif
      endif

   enddo

   ! Turns out there is also a [U,V]_10_METER_WIND
   ! Need to find and count them, too.

   do ivar = 2,max_obs_kinds

      str1   = get_obs_kind_name(ivar-1)
      indx1  = index(str1,'_U_10_METER_WIND') - 1
      indx1N = len_trim(str1)

      str2   = get_obs_kind_name(ivar)
      indx2  = index(str2,'_V_10_METER_WIND') - 1
      indx2N = len_trim(str2)

      if ((indx1 > 0) .and. (indx2 > 0)) then             ! we know we have u,v wind components
         indxN = index(str1(1:indx1),str2(1:indx2))
         if (indxN > 0) then ! we know they are matching kinds
            nwinds = nwinds + 1
            str3   = str1(1:indx2)//'_10_METER_WIND_VELOCITY'
            names(max_obs_kinds + nwinds) = str3
         endif
      endif
   enddo

   ! write(*,*)'There are ',nwinds,' pairs of winds.'

   ! Now that we know how many wind pairs there are, we return the
   ! exact number and new array of observation kind names  

   grok_observation_names = max_obs_kinds + nwinds

   allocate(my_names(grok_observation_names))

   do ivar = 1,grok_observation_names
      my_names(ivar) = names(ivar)
   enddo

   end Function grok_observation_names



   Function NextFile(fname,ifile)
   !----------------------------------------------------------------------
   ! The file name can take one of two forms:
   ! obs_0001/obs_seq.final     -or- 
   ! obs_seq.final
   !
   ! If there is a '/' in the file name, we grab the portion before the
   ! slash and look for an underscore. Anything following the underscore
   ! is presumed to be the portion to increment.
   !
   ! If there is no slash AND ifile is > 1 ... we have already read
   ! the 'one and only' obs_seq.final file and we return 'done'.
   !----------------------------------------------------------------------

   character(len=129), intent(in) :: fname
   integer,            intent(in) :: ifile
   character(len=129)             :: NextFile

   integer,                     SAVE :: filenum = 0
   integer,                     SAVE :: dir_prec = 0
   character(len=129),          SAVE :: dir_base
   character(len=129),          SAVE :: filename
   character(len=stringlength), SAVE :: dir_ext

   character(len=129) :: dir_name
   integer :: slashindex, splitindex

   if (ifile == 1) then ! First time through ... find things.

      filename   = adjustl(fname)
      NextFile   = trim(filename)
      slashindex = index(filename,'/')

      if (slashindex > 0) then ! we have a directory structure

         dir_name   = trim(fname(1:slashindex-1))
         filename   = trim(fname(slashindex+1:129))
         splitindex = index(dir_name,'_')

         if (splitindex <= 0) then
            write(msgstring,*)'cannot find _ in ', trim(dir_name)
            call error_handler(E_ERR,'NextFile',msgstring,source,revision,revdate)
         endif

         dir_base   = dir_name(1:splitindex-1)
         dir_ext    = dir_name(splitindex+1:slashindex-1)
         dir_prec   = slashindex - splitindex - 1

         ! write(*,*)'base is ',trim(dir_base)
         ! write(*,*)'ext  is ',trim(dir_ext)
         ! write(*,*)'prec is ',dir_prec

         read(dir_ext,*) filenum   ! TJH - needs error handling

      else ! we have one single file - on the first trip through

         filenum  = -1 ! indicates no next file

      endif 

   else

      if (filenum < 0) then
         NextFile = 'doneDONEdoneDONE'
      else

         filenum = filenum + 1
         if (dir_prec == 2) then
         write(NextFile,'(a,''_'',i2.2,''/'',a)') trim(dir_base),filenum,trim(filename)
         elseif (dir_prec == 3) then
         write(NextFile,'(a,''_'',i3.3,''/'',a)') trim(dir_base),filenum,trim(filename)
         else
         write(NextFile,'(a,''_'',i4.4,''/'',a)') trim(dir_base),filenum,trim(filename)
         endif

      !  write(*,*)'base is ',trim(dir_base)
      !  write(*,*)'ext  is ',trim(dir_ext)
      !  write(*,*)'prec is ',dir_prec
      !  write(*,*)'ext  is ',filenum
      !  write(*,*)'newfile ',trim(NextFile)

      endif

   endif

   end Function NextFile

  
 
   Subroutine ObsLocationsExist( printswitch )

   ! This routine checks for the existence of observation location files.
   ! Each epoch writes out its own observation location file - and since
   ! it is possible that multiple observation sequence files contribute
   ! to the same epoch ... opening and appending is not a well-posed 
   ! strategy. Furthermore, some people do not start close enough to the
   ! beginning of their epoch definition to ensure that epoch 1 exists.
   ! So ... we check for any of observation_locations.00[1-4].dat
   ! Completely arbitrary.

   logical :: printswitch
   integer :: i

   ! locname, msgstring, source, revision, and revdate are globally-scoped

   if (printswitch) then

      do i = 1,4
      write(locName,'(a,i3.3,a)') 'observation_locations.', i, '.dat'
   
      if (file_exist(locName)) then
         write(msgstring,*)'please remove file(s) like ', trim(locName)
         call error_handler(E_MSG,'ObsLocationsExist',msgstring,source,revision,revdate)
         write(msgstring,*)'Cannot have pre-existing obs location output files. Stopping.'
         call error_handler(E_ERR,'ObsLocationsExist',msgstring,source,revision,revdate)
      endif
      enddo

   endif

   end Subroutine ObsLocationsExist

end program obs_diag
