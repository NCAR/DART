! Data Assimilation Research Testbed -- DART
! Copyright 2004-2009, Data Assimilation Research Section
! University Corporation for Atmospheric Research
! Licensed under the GPL -- www.gpl.org/licenses/gpl.html

module airs_obs_mod

! <next few lines under version control, do not edit>
! $URL: https://subversion.ucar.edu/DAReS/DART/trunk/observations/quikscat/quikscat_JPL_mod.f90 $
! $Id: quikscat_JPL_mod.f90 3809 2009-04-13 16:21:33Z nancy $
! $Revision: 3809 $
! $Date: 2009-04-13 10:21:33 -0600 (Mon, 13 Apr 2009) $

use types_mod,        only : r4, r8, digits12, deg2rad, rad2deg

use obs_def_mod,      only : obs_def_type, get_obs_def_time, read_obs_def, &
                             write_obs_def, destroy_obs_def, interactive_obs_def, &
                             copy_obs_def, set_obs_def_time, set_obs_def_kind, &
                             set_obs_def_error_variance, set_obs_def_location

use time_manager_mod, only : time_type, get_date, set_date, get_time, set_time, &
                             set_calendar_type, GREGORIAN, print_date, print_time, &
                             operator(+), operator(>=)

use    utilities_mod, only : get_unit, open_file, close_file, file_exist, &
                             register_module, error_handler, &
                             E_ERR, E_MSG, timestamp, is_longitude_between

use     location_mod, only : location_type, set_location, VERTISPRESSURE

use     obs_kind_mod, only : get_obs_kind_index, &
                             KIND_TEMPERATURE, KIND_SPECIFIC_HUMIDITY

use obs_kind_mod,     only : AIRS_TEMPERATURE, AIRS_SPECIFIC_HUMIDITY

use obs_sequence_mod, only : init_obs_sequence, init_obs, insert_obs_in_seq, &
                             set_obs_values, set_qc, obs_sequence_type, obs_type, &
                             copy_obs, set_copy_meta_data, set_qc_meta_data, set_obs_def, &
                             get_first_obs, get_last_obs, get_obs_def

use airs_JPL_mod   ! need ', only' list here

implicit none
private

public :: real_obs_sequence, create_output_filename

! version controlled file description for error handling, do not edit
character(len=128), parameter :: &
   source   = "$URL: $", &
   revision = "$Revision: 3809 $", &
   revdate  = "$Date: 2009-04-13 10:21:33 -0600 (Mon, 13 Apr 2009) $"

logical, save :: module_initialized = .false.
character(len=129) :: msgstring

logical :: DEBUG = .false.

! fixed size storage
real ::   T    (AIRS_RET_STDPRESSURELEV, AIRS_RET_GEOXTRACK, AIRS_RET_GEOTRACK)
real ::   T_err(AIRS_RET_STDPRESSURELEV, AIRS_RET_GEOXTRACK, AIRS_RET_GEOTRACK)
real :: MMR    (AIRS_RET_H2OPRESSURELEV, AIRS_RET_GEOXTRACK, AIRS_RET_GEOTRACK)
real :: MMR_err(AIRS_RET_H2OPRESSURELEV, AIRS_RET_GEOXTRACK, AIRS_RET_GEOTRACK)
real ::   Q    (AIRS_RET_H2OPRESSURELEV, AIRS_RET_GEOXTRACK, AIRS_RET_GEOTRACK)
real ::   Q_err(AIRS_RET_H2OPRESSURELEV, AIRS_RET_GEOXTRACK, AIRS_RET_GEOTRACK)
real ::       PBest(                     AIRS_RET_GEOXTRACK, AIRS_RET_GEOTRACK)
real ::       PGood(                     AIRS_RET_GEOXTRACK, AIRS_RET_GEOTRACK)
integer :: nBestStd(                     AIRS_RET_GEOXTRACK, AIRS_RET_GEOTRACK)
integer :: nGoodStd(                     AIRS_RET_GEOXTRACK, AIRS_RET_GEOTRACK)
real ::         lat(                     AIRS_RET_GEOXTRACK, AIRS_RET_GEOTRACK)
real ::         lon(                     AIRS_RET_GEOXTRACK, AIRS_RET_GEOTRACK)
real ::         tim(                     AIRS_RET_GEOXTRACK, AIRS_RET_GEOTRACK)
integer :: qual_h2o(                     AIRS_RET_GEOXTRACK, AIRS_RET_GEOTRACK)

contains


subroutine initialize_module
!-------------------------------------------------
call register_module(source, revision, revdate)

call set_calendar_type(GREGORIAN)

module_initialized = .true.

end subroutine initialize_module


subroutine create_output_filename(l2name, ofname)
!-------------------------------------------------
! The L2 filenames have a very long extension that
! records when the data was published - not very interesting
! for our purposes. replace with something DART-y.
character(len=*), intent(IN)  :: l2name
character(len=*), intent(OUT) :: ofname

integer :: i, basestart, extstart, strlen

! hardcoded and brittle, but for now...  the first 19 chars
! of the input filename have the date & granule number, which
! seems like the bulk of the useful info.  find the last / and
! copy from there to +19 chars.

strlen = len_trim(l2name)

basestart = 1
slashloop : do i = strlen-1,1,-1
   if (l2name(i:i) == '/' ) then
      basestart = i+1
      exit slashloop
   endif
enddo slashloop

extstart = basestart+19-1

ofname = l2name(basestart:extstart)//'.out'
print *, 'output filename = ', ofname

end subroutine create_output_filename



function real_obs_sequence ( granule, lon1, lon2, lat1, lat2 )
!------------------------------------------------------------------------------
!  extract the temperature and humidity observations from a granule
!  and convert to DART observation format.  allow caller to specify
!  a bounding box and only extract data within that region.

type(airs_granule_type), intent(in) :: granule
real(r8), intent(in) :: lon1, lon2, lat1, lat2

! max possible obs from this one granule.
integer :: max_num=  AIRS_RET_STDPRESSURELEV * AIRS_RET_GEOXTRACK * AIRS_RET_GEOTRACK + &
                     AIRS_RET_H2OPRESSURELEV * AIRS_RET_GEOXTRACK * AIRS_RET_GEOTRACK 

type(obs_sequence_type) :: real_obs_sequence
type(obs_def_type)      :: obs_def
type(obs_type)          :: obs, prev_obs

integer :: i, irow, icol, ivert, num_copies, num_qc
integer :: days, seconds
integer :: obs_num
integer :: which_vert, tobstype, qobstype

real(r4) :: speed, dir
real(r8) :: olon, olat, vloc
real(r8) :: obs_value, obs_var
real(r8) :: tqc, qqc
real(r8) :: sintheta, costheta, dirvar, speedvar

type(time_type) :: obs_time, base_time, pre_time, time

character(len = 129) :: meta_data

if ( .not. module_initialized ) call initialize_module

num_copies  = 1
num_qc      = 1

! Initialize an obs_sequence 
call init_obs_sequence(real_obs_sequence, num_copies, num_qc, max_num)

! set meta data of obs_seq
do i = 1, num_copies
   meta_data = 'observation'  
   call set_copy_meta_data(real_obs_sequence, i, meta_data)
end do

do i = 1, num_qc
   meta_data = 'AIRS QC'
   call set_qc_meta_data(real_obs_sequence, i, meta_data)
end do

! Initialize the obs variables
call init_obs(     obs, num_copies, num_qc)
call init_obs(prev_obs, num_copies, num_qc)

! assign each observation the correct observation type
tobstype = get_obs_kind_index('AIRS_TEMPERATURE')
qobstype = get_obs_kind_index('AIRS_SPECIFIC_HUMIDITY')
if ((tobstype < 1) .or. (qobstype < 1)) then
   msgstring = 'unknown observation type [AIRS_TEMPERATURE, AIRS_SPECIFIC_HUMIDITY]'
   call error_handler(E_ERR,'real_obs_sequence',msgstring,source,revision,revdate)
endif



!filename = '../data/200711/20071101/AIRS.2007.11.01.001.L2.RetStd.v5.2.2.0.G08078150655.hdf'


base_time = set_date(1993, 1, 1, 0, 0, 0)   ! reference date: jan 1st, 1993


print *, 'size TairStd', size(granule%TairStd)
T = reshape(granule%TairStd, shape(T))
print *, 'first  row T', T(:, 1, 1)
print *, 'second row T', T(:, 2, 1)
print *, 'second col T', T(:, 1, 2)

print *, 'size TairStdErr', size(granule%TairStdErr)
T_err = reshape(granule%TairStdErr, shape(T))
print *, 'first  row T_err', T_err(:, 1, 1)

print *, 'PBest', size(granule%PBest)
PBest = reshape(granule%PBest, shape(PBest))
print *, 'first  row PBest', PBest(:, 1)

print *, 'PGood', size(granule%PGood)
PGood = reshape(granule%PGood, shape(PGood))
print *, 'first  row PGood', PGood(:, 1)

print *, 'nBestStd', size(granule%nBestStd)
nBestStd = reshape(granule%nBestStd, shape(nBestStd))
print *, 'first  row nBestStd', nBestStd(:, 1)

print *, 'nGoodStd', size(granule%nGoodStd)
nGoodStd = reshape(granule%nGoodStd, shape(nGoodStd))
print *, 'first  row nGoodStd', nGoodStd(:, 1)

print *, 'size water vapor mass mixing ratio', size(granule%H2OMMRStd)
MMR = reshape(granule%H2OMMRStd, shape(MMR))
print *, 'first  row MMR', MMR(:, 1, 1)
print *, 'second col MMR', MMR(:, 2, 1)
print *, 'second col MMR', MMR(:, 1, 2)

print *, 'size water vapor mass mixing ratio err', size(granule%H2OMMRStdErr)
MMR_err = reshape(granule%H2OMMRStdErr, shape(MMR_err))
print *, 'first  row MMR_err', MMR_err(:, 1, 1)

print *, 'water vapor mass mixing ratio quality', size(granule%Qual_H2o)
qual_h2o = reshape(granule%Qual_H2O, shape(qual_h2o))
print *, 'first  row qual_h2o', qual_h2o(:, 1)

Q = (MMR / 1000.) / ( 1.0 + (MMR/1000.0))
Q_err = (MMR_err / 1000.) / ( 1.0 + (MMR_err/1000.0))
print *, 'first  row Q', Q(:, 1, 1)
print *, 'second col Q', Q(:, 2, 1)
print *, 'second col Q', Q(:, 1, 2)

print *, 'first  row Q_err', Q_err(:, 1, 1)

lat = reshape(granule%Latitude, shape(lat))
lon = reshape(granule%Longitude, shape(lon))
tim = reshape(granule%Time, shape(tim))

print *, '1st lat, lon, time', lat(1,1), lon(1,1), tim(1,1)
print *, 'r 2 lat, lon, time', lat(2,1), lon(2,1), tim(2,1)
print *, 'c 2 lat, lon, time', lat(1,2), lon(1,2), tim(1,2)

! compute time of this obs
print *, 'int time ', int(tim(1, 1))
obs_time = base_time + set_time(int(tim(1, 1)))
call print_date(obs_time)
print *, 'int time ', int(tim(AIRS_RET_GEOXTRACK, AIRS_RET_GEOTRACK))
obs_time = base_time + set_time(int(tim(AIRS_RET_GEOXTRACK, AIRS_RET_GEOTRACK)))
call print_date(obs_time)

print *, 'pressStd', granule%pressStd
print *, 'pressH2O', granule%pressH2O

!  loop over all observations within the file
!------------------------------------------------------------------------------

obs_num    = 1
which_vert = VERTISPRESSURE

rowloop:  do irow=1,AIRS_RET_GEOTRACK

   colloop:  do icol=1,AIRS_RET_GEOXTRACK

      olat  = lat(icol,irow) ! valid range [ -90.00,  90.00]
      olon  = lon(icol,irow) ! valid range [-180.00, 180.00]

      ! verify the location is not outside valid limits
      if((olon > 180.0_r8) .or. (olon < -180.0_r8) .or.  &
         (olat >  90.0_r8) .or. (olat <  -90.0_r8)) then
         write(*,*)'WARNING : invalid location.  col,row,lon,lat = ', icol,irow,olon,olat
         cycle colloop
      endif

      ! reject observations outside the bounding box (allowing wrapping)
      if(( olat < lat1) .or. ( olat > lat2 ) .or. &
         (.not. is_longitude_between(olon, lon1, lon2))) cycle colloop

      ! make sure lon is between 0 and 360
      if (olon < 0) olon = olon + 360.0_r8

      obs_time = base_time + set_time(int(tim(icol, irow)))
      call get_time(obs_time, seconds, days)
     

      vert_T_loop:  do ivert=nBestStd(icol, irow), AIRS_RET_STDPRESSURELEV

         tqc = 0   ! FIXME

         ! create the obs_def for this observation, add to sequence
         !---------------------------------------------------------------
      
         vloc = granule%pressStd(ivert)
   
         obs_value = T(ivert, icol, irow)
         obs_var = T_err(ivert, icol, irow) * T_err(ivert, icol, irow)
         call real_obs(num_copies, num_qc, obs, olon, olat, vloc, obs_value, &
                       obs_var, tqc, AIRS_TEMPERATURE, which_vert, seconds, days)
      
         if(obs_num == 1) then ! for the first observation 
            call insert_obs_in_seq(real_obs_sequence, obs)
         else  !  not the first observation 
            if(time >= pre_time) then  ! same time or later than previous obs
               call insert_obs_in_seq(real_obs_sequence, obs, prev_obs)
            else  ! earlier
               call insert_obs_in_seq(real_obs_sequence, obs)
            endif
         endif
         prev_obs = obs
         pre_time = time
   
         if ( DEBUG ) then
            write(*,*)irow, icol, ivert, olat, olon
         endif

         obs_num = obs_num + 1
 
      enddo vert_T_loop

      vert_Q_loop:  do ivert=1,AIRS_RET_H2OPRESSURELEV

         if (qual_h2o(icol, irow) > 0) exit vert_Q_loop

         qqc = 0   ! FIXME

         ! create the obs_def for this observation, add to sequence
         !---------------------------------------------------------------
      
         vloc = granule%pressH2O(ivert)
   
         obs_value = Q(ivert, icol, irow)
         obs_var = Q_err(ivert, icol, irow) * Q_err(ivert, icol, irow)
         call real_obs(num_copies, num_qc, obs, olon, olat, vloc, obs_value, &
                       obs_var, qqc, AIRS_SPECIFIC_HUMIDITY, which_vert, seconds, days)
      
         if(obs_num == 1) then ! for the first observation 
            call insert_obs_in_seq(real_obs_sequence, obs)
         else  !  not the first observation 
            if(time >= pre_time) then  ! same time or later than previous obs
               call insert_obs_in_seq(real_obs_sequence, obs, prev_obs)
            else  ! earlier
               call insert_obs_in_seq(real_obs_sequence, obs)
            endif
         endif
         prev_obs = obs
         pre_time = time
   
         if ( DEBUG ) then
            write(*,*)irow, icol, ivert, olat, olon
         endif

         obs_num = obs_num + 1
 
      enddo vert_Q_loop

   enddo colloop
enddo rowloop

! Print a little summary
print*, 'obs used = ', obs_num, ' obs skipped = ', max_num - obs_num

if ( get_first_obs(real_obs_sequence, obs) ) then
   call get_obs_def(obs, obs_def)
   pre_time = get_obs_def_time(obs_def)
   call print_time(pre_time,' first time in sequence is ')
   call print_date(pre_time,' first date in sequence is ')
endif
if( get_last_obs(real_obs_sequence, obs)) then
   call get_obs_def(obs, obs_def)
   time = get_obs_def_time(obs_def)
   call print_time(time,' last  time in sequence is ')
   call print_date(time,' last  date in sequence is ')
endif
print*, ''

end function real_obs_sequence



subroutine real_obs(num_copies, num_qc, obs, lon, lat, vloc, obs_value, &
                      var2, aqc, obs_kind, which_vert, seconds, days)
!------------------------------------------------------------------------------
integer,        intent(in)    :: num_copies, num_qc
type(obs_type), intent(inout) :: obs
real(r8),       intent(in)    :: lon, lat, vloc, obs_value, var2, aqc
integer,        intent(in)    :: obs_kind, which_vert, seconds, days

integer            :: i
real(r8)           :: aqc01(1), obs_value01(1)
type(obs_def_type) :: obsdef0

if ( .not. module_initialized ) call initialize_module

! Does real initialization of an observation type

call real_obs_def(obsdef0, lon, lat, vloc, &
                    var2, obs_kind, which_vert, seconds, days)
call set_obs_def(obs, obsdef0)

do i = 1, num_copies
   obs_value01(1) = obs_value
   call set_obs_values(obs, obs_value01(1:1) )
end do

do i = 1, num_qc
   aqc01(1) = aqc
   call set_qc(obs, aqc01(1:1))
end do

end subroutine real_obs



subroutine real_obs_def(obs_def, lon, lat, vloc, &
                        var2, obs_kind, which_vert, seconds, days)
!----------------------------------------------------------------------
type(obs_def_type), intent(inout) :: obs_def
real(r8),intent(in) :: lon, lat, vloc, var2
integer, intent(in) :: obs_kind, which_vert, seconds, days

type(location_type) :: loc0

if ( .not. module_initialized ) call initialize_module

! set obs location
loc0 = set_location(lon, lat, vloc, which_vert )
call set_obs_def_location(obs_def, loc0)

! set obs kind
call set_obs_def_kind(obs_def, obs_kind)

call set_obs_def_time(obs_def, set_time(seconds, days) )
call set_obs_def_error_variance(obs_def, var2)

end subroutine real_obs_def


end module airs_obs_mod
