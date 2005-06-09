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
!-----------------------------------------------------------------------

use        types_mod, only : r8
use obs_sequence_mod, only : read_obs_seq, obs_type, obs_sequence_type, get_first_obs, &
                             get_obs_from_key, get_obs_def, get_copy_meta_data, &
                             get_obs_time_range, get_time_range_keys, get_num_obs, &
                             get_next_obs, get_num_times, get_obs_values, init_obs, &
                             assignment(=), get_num_copies, static_init_obs_sequence, &
                             get_qc, destroy_obs_sequence, get_last_obs, get_num_qc
use      obs_def_mod, only : obs_def_type, get_obs_def_error_variance, get_obs_def_time, &
                             get_obs_def_location, get_obs_kind
use     location_mod, only : location_type, get_location, set_location_missing, &
                             write_location, operator(/=)
use time_manager_mod, only : time_type, set_date, set_time, get_time, print_time, &
                             operator(*), &
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
!type(obs_kind_type)     :: obs_kind
type(location_type)     :: obs_loc
type(time_type)         :: first_time, last_time, this_time

!---------------------
integer :: obsindex, i, j, iunit, io
integer :: num_obs_in_set, obs_used_in_set, obs_used = 0

! Storage with fixed size for observation space diagnostics
real(r8) :: prior_mean(1), posterior_mean(1)
real(r8) :: prior_spread(1), posterior_spread(1)
real(r8) :: pr_mean, po_mean ! same as above, without useless dimension 
real(r8) :: pr_sprd, po_sprd ! same as above, without useless dimension

integer :: obs_index, prior_mean_index, posterior_mean_index
integer :: prior_spread_index, posterior_spread_index
integer :: key_bounds(2), flavor

real(r8), allocatable :: obs_err_var(:), obs(:), qc(:)
integer,  allocatable :: keys(:)

logical :: out_of_range, is_there_one, keeper
logical :: is_this_last = .false.

!-----------------------------------------------------------------------
! Namelist with default values
!
character(len = 129) :: obs_sequence_name = "obs_seq.final"
integer :: obs_year   = 0        ! the first date of the diagnostics
integer :: obs_month  = 1
integer :: obs_day    = 1
integer :: tot_days   = 1        ! total days
integer :: iskip      = 0        ! skip the first 'iskip' days
integer :: level      = 1013     ! artificial level for plotting ease.
integer :: obs_select = 1        ! obs type selection: 1=all, 2 =RAonly, 3=noRA
real(r8):: rat_cri    = 3.0      ! QC ratio
real(r8):: qc_threshold = 4.0    ! maximum NCEP QC factor
integer :: bin_separation = 3600 ! Bins every so often seconds
integer :: bin_width = 0         ! width of the bin seconds

namelist /obsdiag_nml/ obs_sequence_name, obs_year, obs_month, obs_day, &
                       tot_days, iskip, level, obs_select, rat_cri, &
                       qc_threshold, bin_separation, bin_width

!-----------------------------------------------------------------------
! Spatial
! Each observation kind gets its own mean, spread, for Guess/Analysis
!-----------------------------------------------------------------------
! index 1 == region 1 == [0.0, 1.0) i.e. Entire domain
! index 2 == region 2 == [0.0, 0.5)
! index 3 == region 3 == [0.5, 1.0)

integer, parameter :: Nregions = 3 

real(r8) :: lonlim1(Nregions), lonlim2(Nregions)

data lonlim1 / 0.0_r8, 0.0_r8, 0.5_r8 /
data lonlim2 / 1.0_r8, 0.5_r8, 1.0_r8 /

integer  :: iregion, iepoch
real(r8) :: rlocation

character(len=3), parameter :: varnames = 'raw'
character(len=5), dimension(Nregions), parameter :: Regions = &
   (/ 'whole','ying ','yang ' /)

!-----------------------------------------------------------------------
! Spatio-Temporal Variables
! Dimension 1 is temporal, actually - these are time-by-region- 
!-----------------------------------------------------------------------

real(r8), allocatable, dimension(:,:) :: rms_ges_mean, rms_ges_spread, &
                                         rms_anl_mean, rms_anl_spread
integer,  allocatable, dimension(:,:) :: num_in_region, num_rejected

real(r8),        allocatable, dimension(:) :: epoch_center
type(time_type), allocatable, dimension(:) :: bincenter

!-----------------------------------------------------------------------
! General purpose variables
!-----------------------------------------------------------------------

integer  :: seconds, days
integer  :: Nepochs, num_copies, num_qc
integer  :: gesUnit, anlUnit

real(r8) :: ratio

type(time_type) :: beg_time, end_time
type(time_type) :: binsep, halfbinwidth 

character(len = 129) :: gesName, anlName, msgstring

!-----------------------------------------------------------------------
! Some variables to keep track of who's rejected why ...
!-----------------------------------------------------------------------

integer :: NwrongType = 0   ! namelist discrimination
integer :: NbadQC     = 0   ! out-of-range QC values
integer :: NbadLevel  = 0   ! out-of-range pressures

!-----------------------------------------------------------------------

call initialize_utilities('obs_diag')
call register_module(source,revision,revdate) 

! Begin by reading the namelist input for obs_diag

if(file_exist('input.nml')) then
   iunit = open_file('input.nml', action = 'read')
   read(iunit, nml = obsdiag_nml, iostat = io)
   if ( io /= 0 ) then
      write(msgstring,*)'obsdiag_nml read error ',io
      call error_handler(E_ERR,'obs_diag',msgstring,source,revision,revdate)
   endif
   call close_file(iunit)
endif
call error_handler(E_MSG,'obs_diag','obsdiag_nml values are',' ',' ',' ') 
write(logfileunit,nml=obsdiag_nml)
write(    *      ,nml=obsdiag_nml)

! Initialize the obs sequence module.
! Read in with enough space for diagnostic output values.

call static_init_obs_sequence()
call read_obs_seq(obs_sequence_name, 0, 0, 0, seq)
num_copies = get_num_copies(seq)
num_qc     = get_num_qc(seq)

call init_obs(observation, num_copies, num_qc)
call init_obs(   next_obs, num_copies, num_qc)

!--------------------------------------------------------------------
! The observations for the low-order models are all exactly at 
! the assimilation timestep. So if we get the first and second
! observation sequence, we know the bin separation. 
!--------------------------------------------------------------------
! Get the time of the last observation in the sequence. For Fun.

is_there_one = get_last_obs(seq, observation)
!if ( is_there_one /= .TRUE. ) then
if ( .not. is_there_one ) then
      call error_handler(E_ERR,'obs_diag','No "last" observation in sequence.', &
      source,revision,revdate)
endif
call get_obs_def(observation, obs_def)
last_time   = get_obs_def_time(obs_def)

! Get the time of the first observation in the sequence.

is_there_one = get_first_obs(seq, observation)
!if ( is_there_one /= .TRUE. ) then
if ( .not. is_there_one ) then
      call error_handler(E_ERR,'obs_diag','No Observations in sequence.', &
      source,revision,revdate)
endif
call get_obs_def(observation, obs_def)
first_time   = get_obs_def_time(obs_def)

! Get the time of the second observation in the sequence.

!call get_next_obs(seq, observation, next_obs, is_this_last)
!call get_obs_def(next_obs, obs_def)
!this_time = get_obs_def_time(obs_def)

! Log what we have and determine the bin separation 

call print_time(first_time,  'first time  ')
!call print_time( this_time,  'second time ')
call print_time( last_time,  'last time   ')
call print_time(first_time,  'first time  ',logfileunit)
!call print_time( this_time,  'second time ',logfileunit)
call print_time( last_time,  'last time   ',logfileunit)

! Then we have to reset to the first observation for processing

!is_there_one = get_first_obs(seq, observation)
!if ( is_there_one /= .TRUE. ) then
!      call error_handler(E_ERR,'obs_diag','No Observations in sequence.', &
!      source,revision,revdate)
!endif
!call get_obs_def(observation, obs_def)

! Get the number of different times in the sequence.
! The low order models assimilate at every time step.

Nepochs = get_num_times(seq)
allocate(    bincenter(Nepochs),           epoch_center(Nepochs))
allocate( rms_ges_mean(Nepochs, Nregions), rms_ges_spread(Nepochs, Nregions), &
          rms_anl_mean(Nepochs, Nregions), rms_anl_spread(Nepochs, Nregions), &
         num_in_region(Nepochs, Nregions),   num_rejected(Nepochs, Nregions) )

write(*,*)'Sequence has ',Nepochs,' different times in it.'

binsep       = set_time(bin_separation, 0)     ! bin separation
halfbinwidth = set_time(bin_width/2, 0)        ! half bin width 

call print_time(      binsep,'binsep       ')
call print_time(halfbinwidth,'halfbinwidth ')
call print_time(      binsep,'binsep       ',logfileunit)
call print_time(halfbinwidth,'halfbinwidth ',logfileunit)

bincenter(1) = first_time
beg_time     = bincenter(1) - halfbinwidth
end_time     = bincenter(1) + halfbinwidth

rms_ges_mean   = 0.0_r8
rms_anl_mean   = 0.0_r8
rms_ges_spread = 0.0_r8
rms_anl_spread = 0.0_r8

!-----------------------------------------------------------------------
! Initialize.
!-----------------------------------------------------------------------

prior_mean(1)       = 0.0_r8
prior_spread(1)     = 0.0_r8
posterior_mean(1)   = 0.0_r8
posterior_spread(1) = 0.0_r8

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
   Advancesets : do iepoch = 1, Nepochs
   !====================================================================

      if ( iepoch > 1 ) then ! Get the next time (if any) in the obs sequence

         if( is_this_last ) exit Advancesets

         call get_next_obs(seq, observation, next_obs, is_this_last)
         call get_obs_def(next_obs, obs_def)
         this_time = get_obs_def_time(obs_def)

         bincenter(iepoch) = bincenter(iepoch-1) + binsep
         beg_time          = bincenter(iepoch) - halfbinwidth
         end_time          = bincenter(iepoch) + halfbinwidth

      endif

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

      allocate(keys(num_obs_in_set), obs_err_var(num_obs_in_set), &
                obs(num_obs_in_set), qc(num_obs_in_set))

      call get_time_range_keys(seq, key_bounds, num_obs_in_set, keys)
      call print_time(this_time,'time is')
      call print_time(this_time,'time is',logfileunit)

      !-----------------------------------------------------------------
      ObservationLoop : do obsindex = 1, num_obs_in_set
      !-----------------------------------------------------------------

         call get_obs_from_key(seq, keys(obsindex), observation) 
         call get_obs_def(observation, obs_def)

         obs_err_var(obsindex) = get_obs_def_error_variance(obs_def) 
!        obs_kind              = get_obs_def_kind(obs_def)
         flavor                = get_obs_kind(obs_def)
         obs_loc               = get_obs_def_location(obs_def)
         rlocation             = get_location(obs_loc) 

         call get_qc(observation,          qc(obsindex:obsindex),         1)
         call get_obs_values(observation, obs(obsindex:obsindex), obs_index)

         ! get interpolated values of prior and posterior ensemble mean

         call get_obs_values(observation,       prior_mean,       prior_mean_index)
         call get_obs_values(observation,   posterior_mean,   posterior_mean_index)
         call get_obs_values(observation,     prior_spread,     prior_spread_index)
         call get_obs_values(observation, posterior_spread, posterior_spread_index)

         pr_mean = prior_mean(1)
         po_mean = posterior_mean(1)
         pr_sprd = prior_spread(1)
         po_sprd = posterior_spread(1)

         !--------------------------------------------------------------
         ! A Whole bunch of reasons to be rejected
         !--------------------------------------------------------------

         keeper = .true.
!        keeper = CheckObsType(obs_select, obs_err_var(obsindex), flavor, ipressure)
         if ( .not. keeper ) then
!        !  write(*,*)'obs ',obsindex,' rejected by CheckObsType ',obs_err_var(obsindex)
            NwrongType = NwrongType + 1
!           cycle ObservationLoop
         endif

         if( qc(obsindex) >= qc_threshold ) then
         !  write(*,*)'obs ',obsindex,' rejected by qc ',qc(obsindex)
            NbadQC = NbadQC + 1
            cycle ObservationLoop 
         endif

         !--------------------------------------------------------------
         ! 
         !--------------------------------------------------------------

         obs_used_in_set = obs_used_in_set + 1
         obs_used        = obs_used        + 1

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

               ratio = GetRatio(obs(obsindex), pr_mean, pr_sprd, &
                        obs_err_var(obsindex))

               if ( ratio <= rat_cri ) then
                  num_in_region(iepoch, iregion) = &
                  num_in_region(iepoch, iregion) + 1

                  rms_ges_mean(iepoch, iregion) = &
                  rms_ges_mean(iepoch, iregion) + (pr_mean - obs(obsindex))**2

                  rms_anl_mean(iepoch, iregion) = &
                  rms_anl_mean(iepoch, iregion) + (po_mean - obs(obsindex))**2

                  rms_ges_spread(iepoch, iregion) = &
                  rms_ges_spread(iepoch, iregion) + pr_sprd**2

                  rms_anl_spread(iepoch, iregion) = &
                  rms_anl_spread(iepoch, iregion) + po_sprd**2
               else
                  num_rejected(iepoch, iregion) = num_rejected(iepoch, iregion) + 1 
               endif

        !   endif

         enddo Areas

      !-----------------------------------------------------------------
      enddo ObservationLoop
      !-----------------------------------------------------------------

      deallocate(keys, obs,  obs_err_var, qc)

      do iregion=1, Nregions
         if (num_in_region(iepoch, iregion) .gt. 0) then
             rms_ges_mean(iepoch, iregion) = sqrt(   rms_ges_mean(iepoch, iregion) / &
                                                    num_in_region(iepoch, iregion) )
             rms_anl_mean(iepoch, iregion) = sqrt(   rms_anl_mean(iepoch, iregion) / &
                                                    num_in_region(iepoch, iregion) )
           rms_ges_spread(iepoch, iregion) = sqrt( rms_ges_spread(iepoch, iregion) / &
                                                    num_in_region(iepoch, iregion) )
           rms_anl_spread(iepoch, iregion) = sqrt( rms_anl_spread(iepoch, iregion) / &
                                                    num_in_region(iepoch, iregion) )
         else
              rms_ges_mean(iepoch, iregion) = -99.0_r8
              rms_anl_mean(iepoch, iregion) = -99.0_r8
            rms_ges_spread(iepoch, iregion) = -99.0_r8
            rms_anl_spread(iepoch, iregion) = -99.0_r8
         endif

      enddo

      write(msgstring,'(''num obs considered in epoch '',i4,'' = '',i8)') iepoch, obs_used_in_set
      call error_handler(E_MSG,'obs_diag',msgstring,source,revision,revdate)

   enddo Advancesets

   call destroy_obs_sequence(seq)

!enddo Dayloop

!-----------------------------------------------------------------------
write(gesName,'(''rawges_times_'',i4.4,''mb.dat'')')level
write(anlName,'(''rawanl_times_'',i4.4,''mb.dat'')')level

gesUnit = get_unit()
OPEN(gesUnit,FILE=trim(adjustl(gesName)),FORM='formatted')
anlUnit = get_unit()
OPEN(anlUnit,FILE=trim(adjustl(anlName)),FORM='formatted')

do i=1, Nepochs

   call get_time(bincenter(i),seconds,days)
!  call print_time(bincenter(i),'bin center ')
!  write(*,*)'bin center (',i,') is ',days, seconds

   write(gesUnit,91) days, seconds, &
       (rms_ges_mean(i,j),rms_ges_spread(i,j),num_in_region(i,j),j=1,Nregions)
   write(anlUnit,91) days, seconds, &
       (rms_anl_mean(i,j),rms_anl_spread(i,j),num_in_region(i,j),j=1,Nregions)

enddo

91 format(i6,1x,i5,3(1x,2f8.3,1x,i4))

close(gesUnit)
close(anlUnit)

write(*,*)''
write(*,*)'# NwrongType            : ',NwrongType
write(*,'('' # QC > '',f9.4,''        :   '',i10)')qc_threshold, NbadQC
write(*,*)'--------------------------------------'
write(*,*)'Rejected Observations   : ',NwrongType+NbadQC
write(*,*)'Considered Observations : ',obs_used

write(logfileunit,*)''
write(logfileunit,*)'# NwrongType            : ',NwrongType
write(logfileunit,'('' # QC > '',f9.4,''        :   '',i10)')qc_threshold, NbadQC
write(logfileunit,*)'--------------------------------------'
write(logfileunit,*)'Rejected Observations   : ',NwrongType+NbadQC
write(logfileunit,*)'Considered Observations : ',obs_used

!-----------------------------------------------------------------------
! Echo attributes to a file to facilitate plotting.
! This file is also a matlab function ... the variables are
! loaded by just typing the name at the matlab prompt.
!-----------------------------------------------------------------------
!  varnames = {'T','W','Q'};
!  varnames = {'T','W','Q','P'};

!  Regions = {'Northern Hemisphere', ...
!             'Southern Hemisphere', ...
!             'Tropics', 'North America'};

iunit = open_file('ObsDiagAtts.m',form='formatted',action='rewind')

!write(iunit,*)'Regions(1) = ',trim(adjustl(Regions(1))),';'
!write(iunit,*)'Regions(2) = ',trim(adjustl(Regions(2))),';'
!write(iunit,*)'Regions(3) = ',trim(adjustl(Regions(3))),';'
write(iunit,'(''varnames  = { ... '')')
write(iunit,'("''raw''",''};'')')
write(iunit,'(''Regions  = { ... '')')
write(iunit,'("''whole''",'',...'')')
write(iunit,'("''ying''",'',...'')')
write(iunit,'("''yang''",''};'')')
do i = 1,Nregions
!   write(iunit,*)trim(adjustl(Regions(i)))//', ...'
!   write(iunit,'('',a(8),'',...'')')trim(adjustl(Regions(i)))
!   write(iunit,100)trim(adjustl(Regions(i)))
enddo
100  format('',a,'',', ...')
write(iunit,'(''obs_year       = '',i,'';'')')obs_year
write(iunit,'(''obs_month      = '',i,'';'')')obs_month
write(iunit,'(''obs_day        = '',i,'';'')')obs_day
write(iunit,'(''tot_days       = '',i,'';'')')tot_days
write(iunit,'(''iskip          = '',i,'';'')')iskip
write(iunit,'(''level          = '',i,'';'')')level
write(iunit,'(''obs_select     = '',i,'';'')')obs_select
write(iunit,'(''rat_cri        = '',f9.2,'';'')')rat_cri
write(iunit,'(''qc_threshold   = '',f9.2,'';'')')qc_threshold
write(iunit,'(''bin_width      = '',i,'';'')')bin_width
write(iunit,'(''bin_separation = '',i,'';'')')bin_separation 
write(iunit,'(''t1             = '',f20.6,'';'')')epoch_center(1)
write(iunit,'(''tN             = '',f20.6,'';'')')epoch_center(Nepochs)
write(iunit,'(''lonlim1 = ['',3(1x,f9.2),''];'')')lonlim1
write(iunit,'(''lonlim2 = ['',3(1x,f9.2),''];'')')lonlim2
close(iunit)

deallocate( epoch_center, bincenter)
deallocate( rms_ges_mean, rms_ges_spread, &
            rms_anl_mean, rms_anl_spread, &
           num_in_region, num_rejected )

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
