! Data Assimilation Research Testbed -- DART
! Copyright 2004, Data Assimilation Initiative, University Corporation for Atmospheric Research
! Licensed under the GPL -- www.gpl.org/licenses/gpl.html

program obs_diag

! <next four lines automatically updated by CVS, do not edit>
! $Source$
! $Revision$
! $Date$
! $Author$
!
use        types_mod, only: r8, pi
use obs_sequence_mod, only: read_obs_seq, obs_type, obs_sequence_type, get_first_obs, &
                get_obs_from_key, set_copy_meta_data, get_copy_meta_data, get_obs_def, &
                   get_obs_time_range, get_time_range_keys, get_num_obs, &
                 get_next_obs, get_num_times, get_obs_values, init_obs, assignment(=), &
                get_num_copies, static_init_obs_sequence, get_qc, destroy_obs_sequence

use obs_def_mod, only : obs_def_type, get_obs_def_error_variance, get_obs_def_time, &
                        get_obs_def_location,  get_obs_def_kind

use location_mod, only : location_type, get_location
use obs_kind_mod, only : obs_kind_type, get_obs_kind
use real_obs_mod, only : obs_model_type

use time_manager_mod, only : time_type, set_time, print_time, operator(/=), operator(>)
use    utilities_mod, only : get_unit, open_file, close_file, register_module, &
                             check_nml_error, file_exist, error_handler, E_ERR, E_MSG, &
                             initialize_utilities, finalize_utilities, logfileunit

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
integer :: i, j, k, ind, iunit, io
integer :: num_obs_in_set, ierr, obs_used
integer :: num_obs_sets, in_obs_copy

! Storage with fixed size for observation space diagnostics
real(r8), allocatable ::  prior_ens_obs(:), posterior_ens_obs(:)
real(r8) :: prior_mean(1), posterior_mean(1)
real(r8) :: prior_spread(1), posterior_spread(1)
real(r8) :: pre_prior_mean(1), pre_posterior_mean(1)
real(r8) :: pre_prior_spread(1), pre_posterior_spread(1)
real(r8) :: pre_obs

integer :: prior_mean_obs_index, posterior_mean_obs_index
integer :: prior_spread_obs_index, posterior_spread_obs_index
integer :: key_bounds(2), kind, model_type, ncep_type
integer :: num_state_copies, num_obs_copies

real(r8), allocatable  :: obs_err_cov(:), obs(:), qc(:)
integer, allocatable :: keys(:)

character(len=129), allocatable :: prior_copy_meta_data(:), posterior_copy_meta_data(:)
logical :: out_of_range, is_there_one, is_this_last

!----------------------------------------------------------------
! Namelist input with default values
!
integer  :: ens_size = 20
logical  :: output_obs_ens_mean   = .true., output_obs_ens_spread   = .true.
integer  :: num_output_obs_members   = 0

character(len = 129) :: obs_sequence_in_name = "obs_seq.final"

namelist /filter_nml/obs_sequence_in_name, output_obs_ens_mean, &
                     output_obs_ens_spread, num_output_obs_members, ens_size

!--------------------------------------------------------------------------
integer, parameter :: max_sets=999, narea = 4       ! /'NH', 'SH', 'TR', 'NA'/
integer, parameter :: nlev=11, nlon=361, nlat=181

integer :: tot_sets, iskip, ktype, n
integer :: ipressure, plev(nlev), pint(nlev+1)
real(r8) :: lon0, lat0, obsloc3(3)
real(r8) :: speed_obs, speed_ges, speed_anl
real(r8) :: speed_obs2, speed_ges2, speed_anl2

integer  :: level, iday, tot_days
real(r8) :: lonlim1(narea), lonlim2(narea), latlim1(narea), latlim2(narea)
real(r8) ::   rms_ges_mean_W(narea),   rms_anl_mean_W(narea), &
            rms_ges_spread_W(narea), rms_anl_spread_W(narea)
integer  ::   num_in_level_W(narea)
real(r8) ::   rms_ges_mean_T(narea),   rms_anl_mean_T(narea), &
            rms_ges_spread_T(narea), rms_anl_spread_T(narea)
integer  ::   num_in_level_T(narea)

real(r8) :: rms_ges_mean_W1(max_sets, narea),   rms_anl_mean_W1(max_sets, narea),  &
          rms_ges_spread_W1(max_sets, narea), rms_anl_spread_W1(max_sets, narea)
integer  :: num_in_level_W1(max_sets, narea)

real(r8) :: rms_ges_mean_T1(max_sets, narea),   rms_anl_mean_T1(max_sets, narea),  &
          rms_ges_spread_T1(max_sets, narea), rms_anl_spread_T1(max_sets, narea)
integer  :: num_in_level_T1(max_sets, narea)

real(r8) ::  rms_ges_ver_W(nlev, narea),  rms_anl_ver_W(nlev, narea)
real(r8) ::  rms_ges_ver_T(nlev, narea),  rms_anl_ver_T(nlev, narea)
real(r8) :: bias_ges_ver_W(nlev, narea), bias_anl_ver_W(nlev, narea)
real(r8) :: bias_ges_ver_T(nlev, narea), bias_anl_ver_T(nlev, narea)
integer  :: num_ver_W(nlev, narea),  num_ver_T(nlev, narea)

real(r8) :: alon(nlon), alat(nlat)
integer  :: k0, kkk, level_index(5)
integer  :: WgesUnit, WanlUnit, TgesUnit, TanlUnit
real(r8) :: rat_cri, ratio
character(len = 6)  day_num(999)          !! total days = 999
character(len = 129) :: WgesName, WanlName, TgesName, TanlName, msgstring

data level_index / 1, 3, 4, 5, 9/
!   level: 1=1000,2=850,3=700,4=500,5=200
data plev / 1000, 925, 850, 700, 500, 400, 300, 250, 200, 150, 100/
data pint / 1025, 950, 900, 800, 600, 450, 350, 275, 225, 175, 125, 75/
!----------------------------------------------------------------

  call initialize_utilities
  call register_module(source,revision,revdate)

! Initialize the obs sequence module
  call static_init_obs_sequence()

! Initialize the observation type variables
  call init_obs(observation, 0, 0) 
  call init_obs(next_obs, 0, 0)

! Begin by reading the namelist input
  if(file_exist('input.nml')) then
   iunit = open_file('input.nml', action = 'read')
   ierr = 1
   do while(ierr /= 0)
      read(iunit, nml = filter_nml, iostat = io, end = 11)
      ierr = check_nml_error(io, 'filter_nml')
   enddo
   11 continue
   call close_file(iunit)
  endif
  write(logfileunit,nml=filter_nml) ! echo namelist params to log file.

  allocate ( prior_ens_obs(ens_size), posterior_ens_obs(ens_size) )

! Determine the number of output obs space fields
    num_obs_copies = 2 * num_output_obs_members
    if(output_obs_ens_mean) num_obs_copies = num_obs_copies + 2
    if(output_obs_ens_spread) num_obs_copies = num_obs_copies + 2

  iunit = get_unit()
  open(iunit, file='obs_diag.in', form='formatted')
  read(iunit,*) tot_days, iskip
  read(iunit,*) level
  read(iunit,*) ktype
  read(iunit,*) rat_cri

  do i=1, tot_days
  read(iunit,FMT='(a6)') day_num(i)
  enddo
  close(iunit)

  print*,'tot_days,iskip,level= ',tot_days,iskip,plev(level_index(level)),'mb'
  print*,'OBS type, QC ratio = ', ktype, rat_cri

  open(iunit, file='Tanl_times_level.dat', form='formatted')
  write(iunit,*) plev(level_index(level)), tot_days, iskip
  close(iunit)

     do i=1, nlon
        alon(i) = i-1
     enddo
     do j=1, nlat
        alat(j) = j-1
     enddo

     do n=1, narea
     do k=1, nlev
     rms_ges_ver_W(k,n) = 0.0
     rms_anl_ver_W(k,n) = 0.0
     rms_ges_ver_T(k,n) = 0.0
     rms_anl_ver_T(k,n) = 0.0

     bias_ges_ver_W(k,n) = 0.0
     bias_anl_ver_W(k,n) = 0.0
     bias_ges_ver_T(k,n) = 0.0
     bias_anl_ver_T(k,n) = 0.0

     num_ver_W(k,n) = 0
     num_ver_T(k,n) = 0
     enddo
     enddo

! set up the areas' limits
  lonlim1(1)= 0.0        ! NH
  lonlim2(1)= 360.0
  latlim1(1)= 110.0
  latlim2(1)= 170.0

  lonlim1(2)= 0.0        ! SH
  lonlim2(2)= 360.0
  latlim1(2)= 10.0
  latlim2(2)= 70.0

  lonlim1(3)= 0.0        ! TR
  lonlim2(3)= 360.0
  latlim1(3)= 70.0
  latlim2(3)= 110.0

  lonlim1(4)= 235.0      ! NA
  lonlim2(4)= 295.0
  latlim1(4)= 115.0
  latlim2(4)= 145.0

  tot_sets = 0
!-------------------------------
  DayLoop : do iday=1, tot_days
!-------------------------------

     write(msgstring,*)'opened ', day_num(iday), trim(obs_sequence_in_name)
     call error_handler(E_MSG,'new_obs_diag',msgstring,source,revision,revdate)

! Read in with enough space for diagnostic output values
    call read_obs_seq(day_num(iday)//obs_sequence_in_name, 0, 0, 0, seq)

     write(msgstring,*)'get_num_copies = ', get_num_copies(seq)
     call error_handler(E_MSG,'new_obs_diag',msgstring,source,revision,revdate)

    if(iday == 1) in_obs_copy = get_num_copies(seq) - num_obs_copies   !! NCEP obs copy =1

! Count of number of sets in the sequence
    num_obs_sets = get_num_times(seq)

     write(msgstring,*)'num_obs_sets =',   num_obs_sets 
     call error_handler(E_MSG,'new_obs_diag',msgstring,source,revision,revdate)

! Initialize the output sequences and state files and set their meta data
  if(iday ==1 ) then
call gen_copy_meta_data(output_obs_ens_mean, output_obs_ens_spread, num_output_obs_members,&
                        num_obs_copies, prior_mean_obs_index, posterior_mean_obs_index, &
                       prior_spread_obs_index, posterior_spread_obs_index, in_obs_copy, seq)
!  print*, 'index= ', prior_mean_obs_index, posterior_mean_obs_index, &
!          prior_spread_obs_index, posterior_spread_obs_index
  endif

! Get the time of the first observation in the sequence
    is_there_one = get_first_obs(seq, observation)
    call get_obs_def(observation, obs_def)
    next_time = get_obs_def_time(obs_def)

!====================================================
Advancesets : do i = 1, num_obs_sets
!====================================================

    do n=1, narea
     rms_ges_mean_W(n)   = 0.0
     rms_anl_mean_W(n)   = 0.0
     rms_ges_mean_T(n)   = 0.0
     rms_anl_mean_T(n)   = 0.0

     rms_ges_spread_W(n) = 0.0
     rms_anl_spread_W(n) = 0.0
     rms_ges_spread_T(n) = 0.0
     rms_anl_spread_T(n) = 0.0

     num_in_level_W(n)   = 0
     num_in_level_T(n)   = 0
    enddo

   call get_obs_time_range(seq, next_time, next_time, key_bounds, num_obs_in_set, &
                           out_of_range, observation)

    obs_used = 0
     write(*, *) 'num_obs_in_set of ', i, ' = ', num_obs_in_set

   allocate(keys(num_obs_in_set), obs_err_cov(num_obs_in_set), obs(num_obs_in_set), &
            qc(num_obs_in_set))

   call get_time_range_keys(seq, key_bounds, num_obs_in_set, keys)
   call print_time(next_time)

! Allocate storage for the ensemble priors for this number of observations

   ! Get the observation value and interpolated ones at obs locations
!--------------------------------------------------------
  ObservationLoop : do j = 1, num_obs_in_set
!--------------------------------------------------------
      call get_obs_from_key(seq, keys(j), observation)

      call get_obs_def(observation, obs_def)

      obs_loc   = get_obs_def_location(obs_def)
      obsloc3 = get_location(obs_loc)

      lon0 = obsloc3(1) + 1.01_r8                            ! 0-360
      lat0 = obsloc3(2) + 1.01_r8 + 90.0_r8                  ! 0-180
      ipressure = 0.01 * obsloc3(3)                          ! mb

      obs_err_cov(j) = get_obs_def_error_variance(obs_def)
      call get_qc(observation, qc(j:j), 1)

      obs_kind  = get_obs_def_kind(obs_def)
      kind = get_obs_kind(obs_kind)
       model_type = kind
!     call obs_model_type(kind, model_type, ncep_type)

!if(mod(j,800) == 0 ) then
! print*, 'kind, type= ', j, kind
! print*, 'loc= ', j, obsloc3(1)*3.14159/180.0, obsloc3(2)*3.14159/180.0, obsloc3(3)
! print*, 'cov, qc= ', j, obs_err_cov(j), qc(j)
!endif

      call get_obs_values(observation, obs(j:j), 1)

!   get interpolated values of prior and posterior ensembles 
     do k = 1, num_output_obs_members
      call get_obs_values(observation,     prior_ens_obs(k:k), in_obs_copy + 2*k-1)
      call get_obs_values(observation, posterior_ens_obs(k:k), in_obs_copy + 2*k  )
     end do

!   get interpolated values of prior and posterior ensemble mean
      call get_obs_values(observation,     prior_mean,     prior_mean_obs_index)
      call get_obs_values(observation, posterior_mean, posterior_mean_obs_index)

!    print*, 'prior_mean_obs_index= ', prior_mean_obs_index, posterior_mean_obs_index

!   get interpolated values of prior and posterior ensemble spread
      call get_obs_values(observation,     prior_spread,     prior_spread_obs_index)
      call get_obs_values(observation, posterior_spread, posterior_spread_obs_index)
!    print*, 'prior_spread_obs_index= ', prior_spread_obs_index, posterior_spread_obs_index
!if(mod(j,800) == 0 ) then
!  print*, 'obs values= ', j, obs(j),prior_mean, posterior_mean
!  print*, 'obs sp values= ', j, prior_spread, posterior_spread
!endif

!--------------------------------------------------------
condition1: if( qc(j) .lt. 4.0 ) then

   obs_used = obs_used + 1

areas : do n =1, narea

!   select the area for statistics
condition2:  if(lon0 .ge.lonlim1(n) .and. lon0  .le. lonlim2(n) .and. &
                lat0 .ge.latlim1(n) .and. lat0  .le. latlim2(n))         then

!--------------------------------------------------------------
!  time series statistics of one specific layer
!--------------------------------------------------------------

!   select the layer interested
condition3:  if(ipressure .le. pint(level_index(level)) .and.    &
                  ipressure .gt. pint(level_index(level)+1))             then
    
variable: if(model_type == 2 ) then               !! Wind v-component

! if(ktype == 2 .and. ncep_type >= 230) go to 250
! if(ktype == 3 .and. ncep_type == 220) go to 250

  if(ktype == 2) then
!   temporary keep RA only skip ACARS and SATWND data 
     if( abs(sqrt( obs_err_cov(j))-2.5) <= 0.1)              go to 250           !!  ACARS wind
     if( sqrt( obs_err_cov(j)) > 3.3)                        go to 250
     if( sqrt( obs_err_cov(j)) > 2.5 .and. ipressure .gt. 420 ) go to 250
     if( sqrt( obs_err_cov(j)) > 1.7 .and. ipressure .gt. 680 ) go to 250 
  endif

  if(ktype == 3) then
!   temporary to keep ACARS and SATWND data  only 
     if( sqrt( obs_err_cov(j)) < 2.5 .and. ipressure < 400 )    go to 250
     if( sqrt( obs_err_cov(j)) > 2.5 .and. sqrt( obs_err_cov(j)) < 3.5 .and. ipressure < 400 ) go to 250
     if( sqrt( obs_err_cov(j)) <= 1.7)                          go to 250
  endif


    ratio = sqrt( (prior_mean(1)-obs(j))**2 + (pre_prior_mean(1)- pre_obs)**2 ) /     &
            sqrt( prior_spread(1)**2 + pre_prior_spread(1)**2 + obs_err_cov(j) )

  if( ratio <= rat_cri )  then
  num_in_level_W(n) = num_in_level_W(n) + 1
  rms_ges_mean_W(n) = rms_ges_mean_W(n) + &
                      (prior_mean(1)-obs(j))**2 + (pre_prior_mean(1)-pre_obs)**2
  rms_anl_mean_W(n) = rms_anl_mean_W(n) + &
                      (posterior_mean(1)-obs(j))**2 + (pre_posterior_mean(1)-pre_obs)**2

  rms_ges_spread_W(n)= rms_ges_spread_W(n)+ &
                       prior_spread(1)**2 + pre_prior_spread(1)**2
  rms_anl_spread_W(n)= rms_anl_spread_W(n)+ &
                       posterior_spread(1)**2 + pre_posterior_spread(1)**2
  endif
 250 continue

elseif (model_type == 4 ) then               !! Temperature

! if(ktype == 2 .and. ncep_type >= 130) go to 252
! if(ktype == 3 .and. ncep_type == 120) go to 252

  if(ktype == 2) then
  !   temporary keep RA only 
     if( abs(sqrt( obs_err_cov(j))-1.0) <= 0.1) go to 252
  endif

  if(ktype == 3) then
!   temporary keep only ACARS and SATWND data  
     if( abs(sqrt( obs_err_cov(j))-1.0) > 0.1) go to 252
  endif

    ratio = abs(prior_mean(1)- obs(j)) / sqrt( prior_spread(1)**2 + obs_err_cov(j) )

    if(ratio <= rat_cri ) then
     num_in_level_T(n)   = num_in_level_T(n) + 1
     rms_ges_mean_T(n)   = rms_ges_mean_T(n) + (prior_mean(1)- obs(j))**2
     rms_anl_mean_T(n)   = rms_anl_mean_T(n) + (posterior_mean(1)- obs(j))**2
     rms_ges_spread_T(n) = rms_ges_spread_T(n) + prior_spread(1)**2
     rms_anl_spread_T(n) = rms_anl_spread_T(n) + posterior_spread(1)**2
    endif
 252 continue

endif variable

endif condition3
!  end of time series statistics


!-------------------------------------------------------------------
!   vertical statistical part
!-------------------------------------------------------------------

condition4: if(iday > iskip)     then
condition5: if(ipressure .le. pint(1) .and. ipressure .gt. pint(nlev+1) ) then

     do kkk=1, nlev
      if(ipressure .le. pint(kkk) .and. ipressure .gt. pint(kkk+1) ) then
      k0 = kkk
      go to 222
      endif
     enddo
  222 continue 
      if(k0 .lt. 1) k0 =1
      if(k0 .gt. nlev) k0 = nlev

variable2: if(model_type == 2 ) then               !! Wind v-component

! if(ktype == 2 .and. ncep_type >= 230) go to 254
! if(ktype == 3 .and. ncep_type == 220) go to 254

!------------------------------------------------
  if(ktype == 2) then
!   temporary keep RA only skip ACARS and SATWND data 
     if( abs(sqrt( obs_err_cov(j))-2.5) <= 0.1)                  go to 254  !!  ACARS wind
     if( sqrt( obs_err_cov(j)) > 3.3)                            go to 254
     if( sqrt( obs_err_cov(j)) > 2.5 .and. ipressure .gt. 420 )  go to 254
     if( sqrt( obs_err_cov(j)) > 1.7 .and. ipressure .gt. 680 )  go to 254 
  endif

  if(ktype == 3) then
!   temporary to keep ACARS and SATWND data  only 
     if( sqrt( obs_err_cov(j)) < 2.5 .and. ipressure < 400 )     go to 254
     if( sqrt( obs_err_cov(j)) > 2.5 .and. sqrt( obs_err_cov(j)) < 3.5 .and. ipressure < 400 ) go to 254
     if( sqrt( obs_err_cov(j)) <= 1.7)                           go to 254
  endif
!------------------------------------------------

    ratio = sqrt( (prior_mean(1)-obs(j))**2 + (pre_prior_mean(1)- pre_obs)**2 ) /     &
            sqrt( prior_spread(1)**2 + pre_prior_spread(1)**2 + obs_err_cov(j) )

  if(ratio <= rat_cri ) then

     speed_ges2 = sqrt ( pre_prior_mean(1)**2 + prior_mean(1)**2 )
     speed_anl2 = sqrt ( pre_posterior_mean(1)**2 + posterior_mean(1)**2 )
     speed_obs2 = sqrt ( pre_obs**2 + obs(j)**2 )

        num_ver_W(k0,n)=      num_ver_W(k0,n)+ 1   
    rms_ges_ver_W(k0,n)=  rms_ges_ver_W(k0,n)+  &
                          (prior_mean(1)-obs(j))**2 +(pre_prior_mean(1)- pre_obs)**2
    rms_anl_ver_W(k0,n)=  rms_anl_ver_W(k0,n)+  &
                          (posterior_mean(1)-obs(j))**2 +(pre_posterior_mean(1)- pre_obs)**2
   bias_ges_ver_W(k0,n)= bias_ges_ver_W(k0,n)+ speed_ges2 - speed_obs2
   bias_anl_ver_W(k0,n)= bias_anl_ver_W(k0,n)+ speed_anl2 - speed_obs2
  endif

254 continue

else if (model_type == 4 ) then               !! Temperature

! if(ktype == 2 .and. ncep_type >= 130) go to 256
! if(ktype == 3 .and. ncep_type == 120) go to 256

  if(ktype == 2) then
  !   temporary keep RA only 
     if( abs(sqrt( obs_err_cov(j))-1.0) <= 0.1) go to 256
  endif

  if(ktype == 3) then
!   temporary keep only ACARS and SATWND data  
     if( abs(sqrt( obs_err_cov(j))-1.0) > 0.1) go to 256
  endif


    ratio = abs(prior_mean(1)- obs(j)) / sqrt( prior_spread(1)**2 + obs_err_cov(j) )

  if(ratio <= rat_cri )  then

         num_ver_T(k0,n) =     num_ver_T(k0,n) + 1  
     rms_ges_ver_T(k0,n) = rms_ges_ver_T(k0,n) + (prior_mean(1)- obs(j))**2
     rms_anl_ver_T(k0,n) = rms_anl_ver_T(k0,n) + (posterior_mean(1)- obs(j))**2
    bias_ges_ver_T(k0,n) = bias_ges_ver_T(k0,n) + prior_mean(1)- obs(j)
    bias_anl_ver_T(k0,n) = bias_anl_ver_T(k0,n) + posterior_mean(1)- obs(j)
  endif

256 continue

endif variable2

endif condition5
endif condition4
!  end of vertical statistics


endif condition2
enddo areas 
endif condition1

    pre_obs = obs(j)
    pre_prior_mean(1) = prior_mean(1)
    pre_prior_spread(1) = prior_spread(1)
    pre_posterior_mean(1) = posterior_mean(1)
    pre_posterior_spread(1) = posterior_spread(1)

 end do ObservationLoop 
!--------------------------------------------------------

! Deallocate storage used for each set
    deallocate(keys, obs,  obs_err_cov, qc)

    do n=1, narea
      if ( num_in_level_W(n) .gt. 1) then
       rms_ges_mean_W(n) =   sqrt( rms_ges_mean_W(n) / num_in_level_W(n) )
       rms_anl_mean_W(n) =   sqrt( rms_anl_mean_W(n) / num_in_level_W(n) )
     rms_ges_spread_W(n) = sqrt( rms_ges_spread_W(n) / num_in_level_W(n) )
     rms_anl_spread_W(n) = sqrt( rms_anl_spread_W(n) / num_in_level_W(n) )
      endif

      if ( num_in_level_T(n) .gt. 1) then
       rms_ges_mean_T(n) =   sqrt( rms_ges_mean_T(n) / num_in_level_T(n) )
       rms_anl_mean_T(n) =   sqrt( rms_anl_mean_T(n) / num_in_level_T(n) )
     rms_ges_spread_T(n) = sqrt( rms_ges_spread_T(n) / num_in_level_T(n) )
     rms_anl_spread_T(n) = sqrt( rms_anl_spread_T(n) / num_in_level_T(n) )
      endif
    enddo

      tot_sets = tot_sets + 1
      print*, '    sets = ', tot_sets
!  print*, 'obs_used = ', obs_used

    do n=1, narea
        rms_ges_mean_W1(tot_sets, n) =   rms_ges_mean_W(n)
        rms_anl_mean_W1(tot_sets, n) =   rms_anl_mean_W(n)
      rms_ges_spread_W1(tot_sets, n) = rms_ges_spread_W(n)
      rms_anl_spread_W1(tot_sets, n) = rms_anl_spread_W(n)

        rms_ges_mean_T1(tot_sets, n) =   rms_ges_mean_T(n)
        rms_anl_mean_T1(tot_sets, n) =   rms_anl_mean_T(n)
      rms_ges_spread_T1(tot_sets, n) = rms_ges_spread_T(n)
      rms_anl_spread_T1(tot_sets, n) = rms_anl_spread_T(n)

        num_in_level_W1(tot_sets, n) =   num_in_level_W(n)
        num_in_level_T1(tot_sets, n) =   num_in_level_T(n)
    enddo

! Get the next time (if any) in the sequence
   call get_next_obs(seq, observation, next_obs, is_this_last)
   if(is_this_last) exit

   call get_obs_def(next_obs, obs_def)
   next_time = get_obs_def_time(obs_def)

end do Advancesets

  call destroy_obs_sequence(seq)

end do Dayloop

!-------------------------------------------------
    write(WgesName,'(''Wges_times_'',i4.4,''mb.dat'')') plev(level_index(level))
    write(WanlName,'(''Wanl_times_'',i4.4,''mb.dat'')') plev(level_index(level))
    write(TgesName,'(''Tges_times_'',i4.4,''mb.dat'')') plev(level_index(level))
    write(TanlName,'(''Tanl_times_'',i4.4,''mb.dat'')') plev(level_index(level))

   WgesUnit = get_unit()
   OPEN(WgesUnit,FILE=trim(adjustl(WgesName)),FORM='formatted')
   WanlUnit = get_unit()
   OPEN(WanlUnit,FILE=trim(adjustl(WanlName)),FORM='formatted')
   TgesUnit = get_unit()
   OPEN(TgesUnit,FILE=trim(adjustl(TgesName)),FORM='formatted')
   TanlUnit = get_unit()
   OPEN(TanlUnit,FILE=trim(adjustl(TanlName)),FORM='formatted')

do i=1, tot_sets
 write(WgesUnit,91) i,(rms_ges_mean_W1(i,n),rms_ges_spread_W1(i,n),num_in_level_W1(i,n),n=1,narea)
 write(WanlUnit,91) i,(rms_anl_mean_W1(i,n),rms_anl_spread_W1(i,n),num_in_level_W1(i,n),n=1,narea)
 write(TgesUnit,91) i,(rms_ges_mean_T1(i,n),rms_ges_spread_T1(i,n),num_in_level_T1(i,n),n=1,narea)
 write(TanlUnit,91) i,(rms_anl_mean_T1(i,n),rms_anl_spread_T1(i,n),num_in_level_T1(i,n),n=1,narea)
enddo

91 format(i2, 4(2f7.2, i8) )

   close(WgesUnit)
   close(WanlUnit)
   close(TgesUnit)
   close(TanlUnit)

!-------------------------------------------------
!  do all day average of the vertical statistics
!-------------------------------------------------
  do n=1, narea
    do k=1, nlev
     if(num_ver_W(k,n) .ne.0) then
       rms_ges_ver_W(k,n) = sqrt( rms_ges_ver_W(k,n) / num_ver_W(k,n) )
       rms_anl_ver_W(k,n) = sqrt( rms_anl_ver_W(k,n) / num_ver_W(k,n) )
      bias_ges_ver_W(k,n) =  bias_ges_ver_W(k,n) / num_ver_W(k,n) 
      bias_anl_ver_W(k,n) =  bias_anl_ver_W(k,n) / num_ver_W(k,n)
     endif

     if(num_ver_T(k,n) .ne.0) then
       rms_ges_ver_T(k,n) = sqrt( rms_ges_ver_T(k,n) / num_ver_T(k,n) )
       rms_anl_ver_T(k,n) = sqrt( rms_anl_ver_T(k,n) / num_ver_T(k,n) )
      bias_ges_ver_T(k,n) =  bias_ges_ver_T(k,n) / num_ver_T(k,n)
      bias_anl_ver_T(k,n) =  bias_anl_ver_T(k,n) / num_ver_T(k,n)
     endif
    enddo
  enddo

     OPEN(185,FILE='Wges_ver_ave.dat',FORM='FORMATTED')
     OPEN(186,FILE='Wanl_ver_ave.dat',FORM='FORMATTED')
     OPEN(285,FILE='Tges_ver_ave.dat',FORM='FORMATTED')
     OPEN(286,FILE='Tanl_ver_ave.dat',FORM='FORMATTED')

     OPEN(195,FILE='Wges_ver_ave_bias.dat',FORM='FORMATTED')
     OPEN(196,FILE='Wanl_ver_ave_bias.dat',FORM='FORMATTED')
     OPEN(295,FILE='Tges_ver_ave_bias.dat',FORM='FORMATTED')
     OPEN(296,FILE='Tanl_ver_ave_bias.dat',FORM='FORMATTED')

     do k=nlev, 1, -1
      write(185, 610) plev(k), (rms_ges_ver_W(k,n), num_ver_W(k,n), n=1, narea)
      write(186, 610) plev(k), (rms_anl_ver_W(k,n), num_ver_W(k,n), n=1, narea)
      write(285, 610) plev(k), (rms_ges_ver_T(k,n), num_ver_T(k,n), n=1, narea)
      write(286, 610) plev(k), (rms_anl_ver_T(k,n), num_ver_T(k,n), n=1, narea)
     enddo

     do k=nlev, 1, -1
      write(195, 610) plev(k), (bias_ges_ver_W(k,n), num_ver_W(k,n), n=1, narea)
      write(196, 610) plev(k), (bias_anl_ver_W(k,n), num_ver_W(k,n), n=1, narea)
      write(295, 610) plev(k), (bias_ges_ver_T(k,n), num_ver_T(k,n), n=1, narea)
      write(296, 610) plev(k), (bias_anl_ver_T(k,n), num_ver_T(k,n), n=1, narea)
     enddo
610 format(i5, 4(f8.2, i8) )

   close(185)
   close(186)
   close(285)
   close(286)

   close(195)
   close(196)
   close(295)
   close(296)

  deallocate ( prior_ens_obs, posterior_ens_obs )

contains

!-----------------------------------------------------------
subroutine gen_copy_meta_data(output_obs_ens_mean, output_obs_ens_spread, &
                          num_output_obs_members, num_obs_copies,  &
                          prior_mean_obs_index, posterior_mean_obs_index, &
           prior_spread_obs_index, posterior_spread_obs_index, in_obs_copy, seq)

implicit none

! Figures out the strings describing the output copies for the three output files.
! THese are the prior and posterior state output files and the observation sequence
! output file which contains both prior and posterior data.

logical, intent(in) :: output_obs_ens_mean, output_obs_ens_spread
integer, intent(in) :: num_output_obs_members, in_obs_copy
integer, intent(out) :: prior_mean_obs_index, posterior_mean_obs_index
integer, intent(out) :: prior_spread_obs_index, posterior_spread_obs_index
integer, intent(out) :: num_obs_copies
type(obs_sequence_type), intent(inout) :: seq

character(len=129) :: prior_meta_data, posterior_meta_data
integer :: i

! Set up the metadata for the output ensemble observations space file
do i = 1, num_output_obs_members
   if(i < 10000) then
      write(prior_meta_data, '(a21, 1x, i6)') 'prior ensemble member', i
      write(posterior_meta_data, '(a25, 1x, i6)') 'posterior ensemble member', i
   else
      write(*, *) 'output metadata in filter needs ensemble size < 10000'
      stop
   endif
!  call set_copy_meta_data(seq, in_obs_copy + 2*i - 1, prior_meta_data)
!  call set_copy_meta_data(seq, in_obs_copy + 2*i, posterior_meta_data)
end do

num_obs_copies = in_obs_copy + 2 * num_output_obs_members
if(output_obs_ens_mean) then
   num_obs_copies = num_obs_copies + 1
   prior_meta_data = 'prior ensemble mean'
!  call set_copy_meta_data(seq, num_obs_copies, prior_meta_data)
   prior_mean_obs_index = num_obs_copies

   num_obs_copies = num_obs_copies + 1
   posterior_meta_data = 'posterior ensemble mean'
!  call set_copy_meta_data(seq, num_obs_copies, posterior_meta_data)
   posterior_mean_obs_index = num_obs_copies 
endif

if(output_obs_ens_spread) then
   num_obs_copies = num_obs_copies + 1
   prior_meta_data = 'prior ensemble spread'
!  call set_copy_meta_data(seq, num_obs_copies, prior_meta_data)
   prior_spread_obs_index = num_obs_copies

   num_obs_copies = num_obs_copies + 1
   posterior_meta_data = 'posterior ensemble spread'
!  call set_copy_meta_data(seq, num_obs_copies, posterior_meta_data)
   posterior_spread_obs_index = num_obs_copies
endif

end subroutine gen_copy_meta_data
!-------------------------------------------------------------------------

end program obs_diag
