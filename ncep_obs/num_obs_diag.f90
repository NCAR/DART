! Data Assimilation Research Testbed -- DART
! Copyright 2004, 2005, Data Assimilation Initiative, University Corporation for Atmospheric Research
! Licensed under the GPL -- www.gpl.org/licenses/gpl.html

program num_obs_diag

! <next five lines automatically updated by CVS, do not edit>
! $Source$
! $Revision$
! $Date$
! $Author$
! $Name$

!-----------------------------------------------------------------------
! This just echoes the number of observations in a set of observation
! windows as defined by bin_separation and bin_width.
!
! The programs defines a series of epochs (periods of time) and geographic
! regions and accumulates statistics for these epochs and regions.
!-----------------------------------------------------------------------

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
integer :: ibintoday, obsindex, i, j, iunit, ierr, io
integer :: num_obs_in_set, obs_used_in_set, obs_used = 0

integer :: obs_index
integer :: key_bounds(2), flavor

logical :: out_of_range, is_there_one, is_this_last, keeper

!-----------------------------------------------------------------------
! Namelist with default values
!

integer, parameter :: tot_days = 1

character(len = 129) :: obs_sequence_name = "obs_seq.final"
integer :: obs_year   = 2003     ! the first date of the diagnostics
integer :: obs_month  = 1
integer :: obs_day    = 1
integer :: level = 500
integer :: obs_select = 1        ! obs type selection: 1=all, 2 =RAonly, 3=noRA
real(r8):: bin_separation = 6.0  ! Bins every so often (hours)
real(r8):: bin_width = 6.0       ! width of the bin (hour)

namelist /numobsdiag_nml/ obs_sequence_name, obs_year, obs_month, obs_day, &
                       level, obs_select, bin_separation, bin_width

!-----------------------------------------------------------------------
! Spatial
! Each observation kind gets its own mean, spread, for Guess/Analysis
!-----------------------------------------------------------------------
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

!-----------------------------------------------------------------------
! Vertical
!-----------------------------------------------------------------------

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

!-----------------------------------------------------------------------
! Temporal Variables
!-----------------------------------------------------------------------

real(r8),        allocatable, dimension(:) :: epoch_center
type(time_type), allocatable, dimension(:) :: bincenter

!-----------------------------------------------------------------------
! General purpose variables
!-----------------------------------------------------------------------

integer  :: seconds, days, DayOne
integer  :: obslevel, k1, kkk, NBinsPerDay, Nepochs
integer  :: calendar_type

real(r8) :: numer, denom, ratio, ratioU

type(time_type) :: beg_time, end_time
type(time_type) :: binsep, binwidth, halfbinwidth 

character(len =   6) :: day_num 
character(len = 129) :: msgstring

!-----------------------------------------------------------------------

call initialize_utilities('num_obs_diag')
call register_module(source,revision,revdate) 
call static_init_obs_sequence()  ! Initialize the obs sequence module 
call init_obs(observation, 0, 0) ! Initialize the observation type variables
call init_obs(   next_obs, 0, 0)

! Begin by reading the namelist input for numobsdiag

if(file_exist('input.nml')) then
   iunit = open_file('input.nml', action = 'read')
   read(iunit, nml = numobsdiag_nml, iostat = io )
   if ( io /= 0 ) then                                                          
      write(msgstring,*)'numobsdiag_nml read error ',io                            
      call error_handler(E_ERR,'num_obs_diag',msgstring,source,revision,revdate)    
   endif
   call close_file(iunit)
endif
call error_handler(E_MSG,'num_obs_diag','numobsdiag_nml values are',' ',' ',' ')
write(logfileunit,nml=numobsdiag_nml)
write(    *      ,nml=numobsdiag_nml)

! Now that we have input, do some checking and setup

level_index = GetClosestLevel(level) 

calendar_type = 3
call set_calendar_type(calendar_type)

NBinsPerDay  = nint( 24.0_r8 / bin_separation )
    binsep   = set_time(nint(bin_separation * 3600.0_r8), 0)
    binwidth = set_time(nint(bin_width      * 3600.0_r8), 0) ! full bin width 
halfbinwidth = set_time(nint(bin_width      * 1800.0_r8), 0) ! half bin width 

!-----------------------------------------------------------------------
! Nepochs is the total number of time intervals of the period requested
!-----------------------------------------------------------------------
Nepochs = NBinsPerDay*tot_days
write(*,*)'Requesting ',Nepochs,' assimilation periods.' 

allocate(bincenter(Nepochs), epoch_center(Nepochs))

iepoch = 0  ! epoch counter

!-----------------------------------------------------------------------
DayLoop : do iday=1, tot_days
!-----------------------------------------------------------------------

   ! Directory/file names are similar to    01_03/obs_seq.final

   write(day_num, '(i2.2,''_'',i2.2,''/'')') obs_month, obs_day + (iday-1) 
   write(msgstring,*)'opening ', day_num, trim(obs_sequence_name)
   call error_handler(E_MSG,'num_obs_diag',msgstring,source,revision,revdate)

   ! Read in with enough space for diagnostic output values

   call read_obs_seq(day_num//obs_sequence_name, 0, 0, 0, seq)

   !--------------------------------------------------------------------
   ! Get the time of the first observation in the sequence.
   ! We slave the epochs to the day of the first observation.
   ! We completely ignore the seconds of the observation. 
   !--------------------------------------------------------------------

   is_there_one = get_first_obs(seq, observation)
   if ( is_there_one /= .TRUE. ) then
      call error_handler(E_ERR,'num_obs_diag','No Observations in sequence.', &
      source,revision,revdate)
   endif
   call get_obs_def(observation, obs_def)
   next_time    = get_obs_def_time(obs_def)
   call get_time(next_time, seconds, DayOne)

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

      bincenter(iepoch) = set_time(0,DayOne) + ibintoday * binsep
      beg_time          = bincenter(iepoch) - halfbinwidth
      end_time          = bincenter(iepoch) + halfbinwidth

      call get_time(bincenter(iepoch),seconds,days)
      epoch_center(iepoch) = days + seconds/86400.0_r8

      call get_obs_time_range(seq, beg_time, end_time, key_bounds, &
                  num_obs_in_set, out_of_range, observation)

      obs_used_in_set = 0
      call print_time(         beg_time,'bin  start ',logfileunit)
      call print_time(bincenter(iepoch),'bin center ',logfileunit)
      call print_time(         end_time,'bin    end ',logfileunit)
      write(logfileunit, *) 'num_obs_in_bin (', iepoch, ') = ', num_obs_in_set
      call print_time(         beg_time,'bin  start ')
      call print_time(bincenter(iepoch),'bin center ')
      call print_time(         end_time,'bin    end ')
      write(     *     , *) 'num_obs_in_bin (', iepoch, ') = ', num_obs_in_set

   enddo Advancesets

   call destroy_obs_sequence(seq)

enddo Dayloop

deallocate(epoch_center, bincenter)

call timestamp(source,revision,revdate,'end') ! That closes the log file, too.

contains


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

   if ( pressure > pint(1) ) then ! pressure greater than 1025
      level_index = -1
      return
   else if ( pressure <= pint(nlev+1) ) then ! pressure less than 75
      level_index = 100 + nlev
      return
   else
      dx = abs(pressure - plev)    ! whole array
      level_index = minval( inds , mask=(dx == minval(dx)))
   endif

   ! This defines bins relative to the pressure interval array, which
   ! was somewhat arbitrarily defined by Hui. There are slight
   ! differences between this method and the 'GetClosestLevel' method
   ! which defines bins relative to the strict midpoint of the plev
   ! array. 
   !PressureLoop: do kkk=1, nlev
   !   if(ipressure .le. pint(kkk) .and. ipressure .gt. pint(kkk+1) ) then
   !   level_index = kkk
   !   exit PressureLoop
   !   endif
   !enddo PressureLoop

   end Function GetClosestLevel

end program num_obs_diag
