! DART software - Copyright 2004 - 2011 UCAR. This open source software
! is
! provided by UCAR, "as is", without charge, subject to all terms of use
! at
! http://www.image.ucar.edu/DAReS/DART/DART_download

program test_ccl

! ALICIA K. January 26th, 2016
use        types_mod
use    utilities_mod
use     location_mod
use      obs_def_mod
use     obs_kind_mod
use time_manager_mod
use cross_comp_localization_mod
use cov_cutoff_mod, only: comp_cov_factor


implicit none

!integer, parameter :: VERTISUNDEF       = -2 ! has no vertical location(undefined)
!integer, parameter :: VERTISSURFACE     = -1 ! surface value
!integer, parameter :: VERTISLEVEL       =  1 ! by level
!integer, parameter :: VERTISPRESSURE    =  2 ! by pressure
!integer, parameter :: VERTISHEIGHT      =  3 ! by height
!integer, parameter :: VERTISSCALEHEIGHT =  4 ! by scale height


type(location_type)      :: loc_obs,loc_state
integer                  :: kind_obs,kind_state
real(r8),dimension(3)    :: this_location_obs, this_location_state
real(r8)                 :: alpha,d1,d2,d3,c,aa,d 

print *,'THIS IS THE TEST_CCL PROGRAM'
!lon,lat,lev

!pressure here is in Pascals.
this_location_obs = (/250.0 , -20.,925/.01/) !pressure in mb/.01 to get Pa
this_location_state = (/240.0 ,-20., 50.0/) !depth in m

loc_obs = set_location(this_location_obs(1),this_location_obs(2),this_location_obs(3),VERTISPRESSURE)
loc_state =set_location(this_location_state(1),this_location_state(2),this_location_state(3),VERTISHEIGHT)
kind_obs = get_obs_kind_index('KIND_TEMPERATURE')
kind_state = get_obs_kind_index('KIND_TEMPERATURE')

!this cutoff should be read from the namelist
c=0.2

d1= get_xcomp_dist(loc_state,loc_obs,kind_state,kind_obs)
aa=comp_cov_factor(d1,c)
print*,'factor based on finalroutine = ' , aa


!alpha = get_vert_alpha(loc_obs,loc_state,kind_obs,kind_state)
!print*,'factor based on vertical only = ' , alpha
!d1 = get_vert_dist(alpha)
!print*,'vert dist in radians = ', d1
!alpha = get_vert_alpha(loc_state,loc_obs,kind_state,kind_obs)
!print*,'factor based on vertical only = ' , alpha
!d2 = get_vert_dist(alpha)
!print*,'vert dist in radians = ', d2


!d3=get_dist(loc_obs,loc_state,0,0,.true.) 
!print*,'horizontal distance in radians = ', d3
!aa = comp_cov_factor(d3,c)
!print*,'factor based on horizontal only = ', aa
!d=sqrt(d1**2 + d3**2);
!aa = comp_cov_factor(d,c)
!print*,'1: factor based on horizontal + vert = ', aa
!d=sqrt(d2**2 + d3**2);
!aa = comp_cov_factor(d,c)
!print*,'2: factor based on horizontal + vert = ', aa
end program test_ccl



