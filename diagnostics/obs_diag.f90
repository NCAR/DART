! Data Assimilation Research Testbed -- DART
! Copyright 2004, Data Assimilation Initiative, University Corporation for Atmospheric Research
! Licensed under the GPL -- www.gpl.org/licenses/gpl.html

program obs_diag

! <next four lines automatically updated by CVS, do not edit>
! $Source$
! $Revision$
! $Date$
! $Author$

use        types_mod, only : r8
use obs_sequence_mod, only : obs_sequence_type, &
   read_obs_sequence, get_num_obs_sets, get_obs_sequence_time, &
   get_num_obs_in_set, get_expected_obs, get_obs_values, &
   get_obs_location1, get_obs_kind1, get_single_obs_value
use time_manager_mod, only : time_type, set_time, print_time
use    utilities_mod, only :  get_unit, open_file, close_file, &
   check_nml_error, file_exist, error_handler, FATAL
use  assim_model_mod, only : assim_model_type
use typesizes            ! from netCDF F90 interface

implicit none

! CVS Generated file description for error handling, do not edit
character(len=128) :: &
source   = "$Source$", &
revision = "$Revision$", &
revdate  = "$Date$"

! Define a type for doing direct access to ensemble state vectors

type(obs_sequence_type) :: seq, prior_seq, posterior_seq
type(time_type)         :: time, time2

integer :: i, j, k, ind, iunit, prior_obs_unit, posterior_obs_unit, io
integer :: ierr
integer :: num_obs_sets, obs_type

! Storage for direct access to ensemble state vectors

real(r8), allocatable :: obs(:)
integer               :: ifirst = 1
integer, allocatable  :: num_obs_in_set(:)
real(r8), allocatable :: ges(:, :), anl(:, :)
real(r8), allocatable :: obsloc(:, :), obskind(:)
real(r8)              :: pressure
integer               :: lon, lat

!----------------------------------------------------------------
! Namelist input with default values
!
integer  :: async = 0, ens_size    = 20
real(r8) :: cutoff      = 200.0_r8
real(r8) :: cov_inflate = 1.0_r8
logical  :: start_from_restart = .false., output_restart = .false.
integer  :: init_time_days    = -1
integer  :: init_time_seconds = -1
! Control diagnostic output for state variables
logical :: output_state_ens_mean = .true., output_state_ens_spread = .true.
integer :: num_output_ens_members = 0
integer :: output_interval = 1

integer, parameter :: nlev = 55, nlon=361, nlat=181 

integer :: ens_size_plus , level
real(r8), allocatable :: rms_ges_mean(:), rms_anl_mean(:)
real(r8), allocatable :: rms_ges_spread(:), rms_anl_spread(:)
real(r8), allocatable :: rms_ges(:,:), rms_anl(:,:)

real(r8) :: plev(nlev), level_value(4)
real(r8) :: rms_ges_ver(nlev), rms_anl_ver(nlev)                                    
real(r8) :: rms_ges_hori(nlon, nlat), rms_anl_hori(nlon, nlat)                
real(r8) :: alon(nlon), alat(nlat)

integer, allocatable :: num_in_level(:)
integer :: k0, kkk,  num_ver(nlev) , num_hori(nlon, nlat), iy,im,id,ih, ista
character(len = 7) :: staid

data  level_value / 850, 700, 500, 200/

character(len = 129) :: obs_sequence_file_name = "obs_sequence", &
                        restart_in_file_name = 'filter_restart_in', &
                        restart_out_file_name = 'filter_restart_out'

namelist /filter_nml/async, ens_size, cutoff, cov_inflate, &
   start_from_restart, output_restart, &
   obs_sequence_file_name, restart_in_file_name, restart_out_file_name, &
   init_time_days, init_time_seconds, output_state_ens_mean, &
   output_state_ens_spread, num_output_ens_members, output_interval
!----------------------------------------------------------------

   ens_size_plus = ens_size + 1

   ! Begin by reading the namelist input
   if(file_exist('input.nml')) then
      iunit = open_file(file = 'input.nml', action = 'read')
      ierr = 1
      do while(ierr /= 0)
         read(iunit, nml = filter_nml, iostat = io, end = 11)
         ierr = check_nml_error(io, 'filter_nml')
      enddo
 11 continue
      call close_file(iunit)
   endif

   write(*, *) 'the ensemble size is ', ens_size
   write(*,*) 'input observation type for verification, 1=u,2=v,3=ps,4=t'
   read(*,*) obs_type

   ! Input the obs_sequence
   iunit = get_unit()
   open(unit = iunit, file = obs_sequence_file_name)
   seq = read_obs_sequence(iunit)
   close(iunit)

   ! Count of number of sets in the sequence
   num_obs_sets = get_num_obs_sets(seq)

!---------------------------------------------------------------
!  input the sequences of guess and analyses at observation space 
!---------------------------------------------------------------

   prior_obs_unit = get_unit()
   open(unit = prior_obs_unit, file = 'prior_obs_diagnostics')

   prior_seq =  read_obs_sequence(prior_obs_unit)
   close(prior_obs_unit)

   posterior_obs_unit = get_unit()
   open(unit = posterior_obs_unit, file = 'posterior_obs_diagnostics')

   posterior_seq = read_obs_sequence(posterior_obs_unit)
   close(posterior_obs_unit)


   allocate(num_obs_in_set(num_obs_sets), num_in_level(num_obs_sets) ) 
   allocate(rms_ges_mean(num_obs_sets), rms_anl_mean(num_obs_sets) ) 
   allocate(rms_ges_spread(num_obs_sets), rms_anl_spread(num_obs_sets) ) 
   allocate(rms_ges(num_obs_sets, ens_size_plus), rms_anl(num_obs_sets, ens_size_plus)) 

!  start the rms errors statistics.
   do k=1, nlev
      rms_ges_ver(k) = 0.0_r8
      rms_anl_ver(k) = 0.0_r8
      num_ver(k) = 0
   enddo

   do i=1, nlon
   do j=1, nlat
      rms_ges_hori(i,j) = 0.0_r8
      rms_anl_hori(i,j) = 0.0_r8
      num_hori(i,j) = 0
   enddo
   enddo

   do i = 1, num_obs_sets
      rms_ges_mean(i)   = 0.0_r8
      rms_ges_spread(i) = 0.0_r8

      rms_anl_mean(i)   = 0.0_r8
      rms_anl_spread(i) = 0.0_r8

      num_in_level(i) = 0

      do k=1, ens_size
         rms_ges(i,k) = 0.0_r8
         rms_anl(i,k) = 0.0_r8
      enddo
   enddo

   do k=1, nlev                                                                          
      plev(k) = 20.0_r8*(k-1)                                                                  
   enddo          

   do i=1, nlon
      alon(i) = i-1
   enddo
   do j=1, nlat
      alat(j) = j-1
   enddo

!  Loop through the time intervals
!-----------------------------------
AdvanceTime : do i = 1, num_obs_sets  
!-----------------------------------

   call get_obs_sequence_time(seq, i, time)
!  call print_time(time)

!  get the number of the observations in this set
   num_obs_in_set(i) = get_num_obs_in_set(seq, i)

!  Allocate storage for the ensemble priors for this number of observations
   allocate(obs(num_obs_in_set(i)), obsloc(num_obs_in_set(i), 3)) 
   allocate( obskind(num_obs_in_set(i)) ) 

   allocate(ges(num_obs_in_set(i),ens_size_plus), anl(num_obs_in_set(i),ens_size_plus)) 

!  Read in the observation values
   call get_obs_values(seq, i, obs, 1)

!  Get the observation locations & kinds
   call get_obs_location1(seq, i, obsloc)
   call get_obs_kind1(seq, i, obskind)

!  Read in the ensemble guess and analyses prior and posterior from diagnostic files
   do j = 1, num_obs_in_set(i)

      do k = 1, ens_size_plus
       call get_single_obs_value(prior_seq, i, j, ges(j,k), k)
       call get_single_obs_value(posterior_seq, i, j, anl(j,k), k)
      end do

!     Write(*,*) i,obs(j),ges(j,1),anl(j,1),obsloc(j,1),obsloc(j,2),obsloc(j,3),obskind(j)

      lon = obsloc(j,1) * 180.0_r8/pi + 0.5_r8             !! in degree
      lat = obsloc(j,2) * 180.0_r8/pi + 0.5_r8 + 90.0_r8   !! in degree from 0 - 180
      pressure = obsloc(j,3) * 0.01_r8                     !! in mb
!     print*, 'lon= ', lon, obsloc(j,1) * 180.0/pi, lat, obsloc(j,2) * 180.0/pi+90.0

      if (obskind(j) .ne. obs_type ) go to 4440

variable: if(obskind(j) == 3 ) then     !! for Ps

      obs(j) = obs(j) * 0.01_r8          !!  convert to mb
      do k = 1, ens_size_plus
         ges(j,k) = ges(j,k) * 0.01_r8   !!  convert to mb
         anl(j,k) = anl(j,k) * 0.01_r8   !!  convert to mb
      enddo

      rms_ges_mean(i) = rms_ges_mean(i) + (ges(j, ens_size_plus)- obs(j))**2
      rms_anl_mean(i) = rms_anl_mean(i) + (anl(j, ens_size_plus)- obs(j))**2
      num_in_level(i) = num_in_level(i) + 1

      do k = 1, ens_size
         rms_ges(i,k) = rms_ges(i,k) + (ges(j, k)- obs(j))**2
         rms_anl(i,k) = rms_anl(i,k) + (anl(j, k)- obs(j))**2
      enddo

!     horizontal distribution of the error for ensemble mean
      num_hori(lon,lat)     = num_hori(lon,lat) + 1                                                          
      rms_ges_hori(lon,lat) = rms_ges_hori(lon,lat) + (ges(j, ens_size_plus)- obs(j))**2
      rms_anl_hori(lon,lat) = rms_anl_hori(lon,lat) + (anl(j, ens_size_plus)- obs(j))**2

else   ! for T & wind components

      if(ifirst == 1) then
         Write(*,*) 'input vertical level for time series; 1=850,2=700,3=500,4=200, 5=all levels'
         read(*,*) level
         ifirst = 0
      endif
    
!  for ensemble mean and every ensemble member

      if( level .lt. 5 ) then
         if( abs( level_value(level) - pressure) .lt. 25.0_r8 ) then

            rms_ges_mean(i) = rms_ges_mean(i) + (ges(j, ens_size_plus)- obs(j))**2
            rms_anl_mean(i) = rms_anl_mean(i) + (anl(j, ens_size_plus)- obs(j))**2
            num_in_level(i) = num_in_level(i) + 1

            do k = 1, ens_size
               rms_ges(i,k) = rms_ges(i,k) + (ges(j, k)- obs(j))**2
               rms_anl(i,k) = rms_anl(i,k) + (anl(j, k)- obs(j))**2
            enddo

!           horizontal distribution of the error for ensemble mean
            num_hori(lon,lat)     = num_hori(lon,lat) + 1                                                          
            rms_ges_hori(lon,lat) = rms_ges_hori(lon,lat) + (ges(j, ens_size_plus)- obs(j))**2
            rms_anl_hori(lon,lat) = rms_anl_hori(lon,lat) + (anl(j, ens_size_plus)- obs(j))**2

         endif
      else
         rms_ges_mean(i) = rms_ges_mean(i) + (ges(j, ens_size_plus)- obs(j))**2
         rms_anl_mean(i) = rms_anl_mean(i) + (anl(j, ens_size_plus)- obs(j))**2
         num_in_level(i) = num_in_level(i) + 1

         do k = 1, ens_size
            rms_ges(i,k) = rms_ges(i,k) + (ges(j, k)- obs(j))**2
            rms_anl(i,k) = rms_anl(i,k) + (anl(j, k)- obs(j))**2
         enddo

!        horizontal distribution of the error for ensemble mean
         num_hori(lon,lat) = num_hori(lon,lat) + 1                                                          
         rms_ges_hori(lon,lat) = rms_ges_hori(lon,lat) + (ges(j, ens_size_plus)- obs(j))**2
         rms_anl_hori(lon,lat) = rms_anl_hori(lon,lat) + (anl(j, ens_size_plus)- obs(j))**2

      endif

!     averaged vertical profile of the error for ensemble mean
      k0 = ifix((pressure+10.0_r8)/20 ) + 1            

      num_ver(k0)     = num_ver(k0) + 1                                                          
      rms_ges_ver(k0) = rms_ges_ver(k0) + (ges(j, ens_size_plus)- obs(j))**2
      rms_anl_ver(k0) = rms_anl_ver(k0) + (anl(j, ens_size_plus)- obs(j))**2

endif variable
 4440 continue

    end do

    deallocate( obs, obsloc, obskind, ges, anl)
!------------------------------------------------------------------------

end do AdvanceTime

   do i = 1, num_obs_sets    !  time interval

      if ( num_in_level(i) .gt. 1) then
       rms_ges_mean(i) = sqrt( rms_ges_mean(i) / (num_in_level(i)-1) )
       rms_anl_mean(i) = sqrt( rms_anl_mean(i) / (num_in_level(i)-1) )

      do k = 1, ens_size
       rms_ges(i,k) = sqrt( rms_ges(i,k) / (num_in_level(i)-1) )
       rms_anl(i,k) = sqrt( rms_anl(i,k) / (num_in_level(i)-1) )
      enddo
     endif

   enddo

     do k=1, nlev
     if(num_ver(k) .ne.0) then
     rms_ges_ver(k) = sqrt( rms_ges_ver(k) / num_ver(k) )
     rms_anl_ver(k) = sqrt( rms_anl_ver(k) / num_ver(k) )
     endif
     enddo

     do i=1, nlon
     do j=1, nlat
     if(num_hori(i,j) .ne.0) then
     rms_ges_hori(i,j) = sqrt( rms_ges_hori(i,j) / num_hori(i,j) )
     rms_anl_hori(i,j) = sqrt( rms_anl_hori(i,j) / num_hori(i,j) )
     endif
     enddo
     enddo

!   for ensemble spread over observation space
     do i = 1, num_obs_sets    !  time interval
       rms_ges_spread(i) = (rms_ges(i,1) - rms_ges_mean(i) )**2
       rms_anl_spread(i) = (rms_anl(i,1) - rms_anl_mean(i) )**2
      do k=2, ens_size
       rms_ges_spread(i) = rms_ges_spread(i) + (rms_ges(i,k) - rms_ges_mean(i) )**2
       rms_anl_spread(i) = rms_anl_spread(i) + (rms_anl(i,k) - rms_anl_mean(i) )**2
      enddo
       rms_ges_spread(i) = sqrt( rms_ges_spread(i) / (ens_size -1) )
       rms_anl_spread(i) = sqrt( rms_anl_spread(i) / (ens_size -1) )
     enddo


!----------------------------------------
!    write out the statistics for ploting
!----------------------------------------

     OPEN(285,FILE='ges_fitobs_ver.dat',FORM='FORMATTED')
     OPEN(286,FILE='anl_fitobs_ver.dat',FORM='FORMATTED')
     do k=nlev-4, 6, -1
     if(num_ver(k) .ne.0) then
      write(285, *) rms_ges_ver(k), plev(k), num_ver(k)
      write(286, *) rms_anl_ver(k), plev(k), num_ver(k)
     endif
     enddo
     close(285)
     close(286)

     OPEN(285,FILE='ges_fitobs.dat',FORM='FORMATTED')
     OPEN(286,FILE='anl_fitobs.dat',FORM='FORMATTED')
     do i=1, num_obs_sets
      write(285, *) i, rms_ges_mean(i), rms_ges_spread(i), num_in_level(i)
      write(286, *) i, rms_anl_mean(i), rms_anl_spread(i), num_in_level(i)
     enddo
     close(285)
     close(286)

deallocate(rms_ges_mean, rms_anl_mean, rms_ges_spread, rms_anl_spread, num_in_level)

!------------------------------------------------------------------------------
!   write irregular horizontal distribution of the error in GRADS station format
!------------------------------------------------------------------------------
      open(80,file='rms_hori.dat',form='formatted')
!   no need to change the date here, only for output for consistent purpose later.
      iy = 2003                                              
      im = 7
      id = 23                                            
      ih = 0                                                

      do  j=1, nlat                                                
      do  i=1, nlon                                                   
       ista = (j-1)*nlon + i                                                    
       write(staid,310) ista           
       if(num_hori(i,j).ne.0) &
       write(80,888) iy,im,id,ih,staid,alat(j)-90.0,alon(i),rms_ges_hori(i,j)  
      enddo      
      enddo      

      ih = 1                                                  
      do  j=1, nlat                                                
      do  i=1, nlon                                                   
       ista = (j-1)*nlon + i                                                    
       write(staid,310) ista           
       if(num_hori(i,j).ne.0) &
       write(80,888) iy,im,id,ih,staid,alat(j)-90.0,alon(i),rms_anl_hori(i,j)  
      enddo      
      enddo      

  310 format(i7)                                            
  888 format(4i4,4x,a7,1x,3f8.2)                                               
      close(80)


end program obs_diag
