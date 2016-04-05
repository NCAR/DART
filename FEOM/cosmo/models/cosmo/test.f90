PROGRAM test

! <next few lines under version control, do not edit>
! $URL$
! $Id$
! $Revision$
! $Date$

!----------------------------------------------------------------------
! purpose: test routines
!----------------------------------------------------------------------

  use        types_mod, only : r4, r8, digits12, SECPERDAY, MISSING_R8,          &
                               rad2deg, deg2rad, PI

  use time_manager_mod, only : time_type, set_time, set_date, get_date, get_time,&
                               print_time, print_date, set_calendar_type,        &
                               operator(*),  operator(+), operator(-),           &
                               operator(>),  operator(<), operator(/),           &
                               operator(/=), operator(<=)

  use        model_mod, only : static_init_model,get_model_size,get_state_meta_data,&
                               model_interpolate, state_vector,grib_to_sv

  use     location_mod, only : location_type, get_dist, query_location,          &
                               get_close_maxdist_init, get_close_type,           &
                               set_location, get_location, horiz_dist_only,      & 
                               vert_is_undef,    VERTISUNDEF,                    &
                               vert_is_surface,  VERTISSURFACE,                  &
                               vert_is_level,    VERTISLEVEL,                    &
                               vert_is_pressure, VERTISPRESSURE,                 &
                               vert_is_height,   VERTISHEIGHT,                   &
                               get_close_obs_init, loc_get_close_obs => get_close_obs

  IMPLICIT NONE

! version controlled file description for error handling, do not edit
character(len=128), parameter :: &
   source   = "$URL$", &
   revision = "$Revision$", &
   revdate  = "$Date$"

  INTEGER :: var_type,index,istat,n,i
  integer,allocatable :: seed(:)
  TYPE(location_type) :: loc
  REAL(r8) :: lo,la,h,val
  REAL(r8),allocatable :: x(:)
  REAL(r8) :: lop,lap,x1,lo1,la1,x2,lo2,la2,xo,loo,lao
  TYPE(time_type) :: time

!  CALL static_init_model()
!  print*,get_model_size()

!  index=1940880
!  CALL get_state_meta_data(index,loc,var_type)
!  print*,get_location(loc)
!  print*,var_type


!  n=get_model_size()
!  allocate(x(1:n))
!  allocate(seed(1:n))
!  seed(1:n)=1
!
!  CALL RANDOM_SEED (PUT = seed (1:n))
!  CALL RANDOM_NUMBER(x)
!  x=x*10-5

!  h=70500.
!  lo=0.8
!  la=54.5
!  loc=set_location(lo,la,h,VERTISPRESSURE)
!  call model_interpolate(x,loc,1,val,istat)
!  print*,istat

STOP
  h=10.
  lo=22.9
  la=53.2
  loc=set_location(lo,la,h,VERTISLEVEL)

  call model_interpolate(x,loc,1,val,istat)
  print*,istat

  h=10.
  lo=22.95
  la=53.2
  loc=set_location(lo,la,h,VERTISLEVEL)

  call model_interpolate(x,loc,1,val,istat)
  print*,istat

  h=10.
  lo=22.97
  la=53.2
  loc=set_location(lo,la,h,VERTISLEVEL)

  call model_interpolate(x,loc,1,val,istat)
  print*,istat

  h=10.
  lo=23.0
  la=53.2
  loc=set_location(lo,la,h,VERTISLEVEL)

  call model_interpolate(x,loc,1,val,istat)
  print*,istat

END PROGRAM test
