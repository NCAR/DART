program create_ncep_obs_set_def

!
! <next four lines automatically updated by CVS, do not edit>
! $Source$
! $Revision$
! $Date$
! $Author$
!

! Main program to create a set_def_list description file for a 
! paticular obs_kind set. This is a prototype for more
! user friendly GUI driven methods of creating set_def_lists. For now, there
! is no support for nested observation definition subsets, but that will be
! needed in the long run.

use        types_mod, only : r8
use    utilities_mod, only : open_file
use      obs_def_mod, only : obs_def_type, init_obs_def, interactive_obs_def, &
                             read_ncep_obs_def
use  obs_set_def_mod, only : obs_set_def_type, init_obs_set_def, add_obs
use set_def_list_mod, only : set_def_list_type, init_set_def_list, &
                             add_to_list, write_set_def_list, read_set_def_list

implicit none

! CVS Generated file description for error handling, do not edit
character(len=128) :: &
source   = "$Source$", &
revision = "$Revision$", &
revdate  = "$Date$"

type(obs_def_type) :: obs_def
type(set_def_list_type) :: set_def_list
type(obs_set_def_type) :: obs_set_def

integer :: max_sets, num_obs
integer :: i, j, obs_set_def_index, unit
character(len = 129) :: file_name

integer :: obs_year, obs_month, obs_day
integer :: obsunit, obsunit2
character(len = 2) :: obstime
character(len = 8) :: obsdate
character(len = 80) :: obsfile, obsfile2



file_name = 'set_def.out'
write(*, *) 'filename for output of observation = ', file_name

open(800, file='/home/hliu/DART/ncep_obs/ncepobs.input', form='formatted')  
read(800, *) obs_year, obs_month, obs_day                                 
read(800, *) max_sets                                                  
close(800)                                                   
                                                               
write(*, *) 'number of unique observation sets is set to', max_sets 
                                                              
if(obs_month .lt. 10) then                                   
   write(obsdate, '(i4,A1,i1,i2)') obs_year, '0',obs_month, obs_day   
else                                                             
   write(obsdate, '(i4,i2,i2)') obs_year, obs_month, obs_day           
endif                                                  
print*, 'ncep obsdate = ', obsdate                         

 set_def_list = init_set_def_list(max_sets)

!  Loop to get definitions for each set of the 1-day period
 maxset: do i = 1, max_sets

     if(i.lt.10) then  
     write(obstime, '(A1, i1)') '0', i  
     else     
     write(obstime, '( i2)') i                                                           
     endif                                           

!  open NCEP observation data file ( hourly interval files)
    obsunit = 800 + i
!   obsfile='/home/hliu/ncepobs/all/obs.'//obstime//'.'//obsdate
    obsfile='/home/hliu/ncepobs/all/temp_obs.'//obstime//'.'//obsdate

    open (obsunit , file=obsfile , form='formatted')
    print*, 'obs file opened = ', obsfile

     call num_ncep_obs(obsunit, num_obs)
     write(*, *) 'number of observations in set ', i, ' is ', num_obs
     if(num_obs .ge. 900000) print*, 'num_obs is too large !!!'

!  Initialize the obs_set_def
      obs_set_def = init_obs_set_def(num_obs)

   rewind(obsunit)

   do j = 1, num_obs
    obs_def = read_ncep_obs_def(obsunit)
!   Insert this obs_def into an obs_set
    call add_obs(obs_set_def, obs_def)
   end do
   close (obsunit)

!  Insert this obs_set_def into the list
    obs_set_def_index = add_to_list(set_def_list, obs_set_def)

 end do maxset

!  Output the set_def_list
  unit = open_file(file_name, action = 'write')
  call write_set_def_list(unit, set_def_list)
  close(unit)

 write(*, *)trim(adjustl(file_name)),' successfully created.'


contains

 subroutine num_ncep_obs(obsunit, num_obs)

 integer, intent(in) :: obsunit
 integer, intent(out) :: num_obs
 integer :: obs_prof, obsindex
 real (r8) :: var, lon,lat,lev,zob,dummy,time,type,count, aaa
 integer :: ktot

! read and count the observation number in each hourly interval
    rewind (obsunit)

    num_obs = 0
    do ktot=1, 9000000
     read(obsunit,880,end=200) var, lon, lat, lev, zob, dummy,count,time,type

      num_obs = num_obs + 1

!     print out for test only
!     write(obsunit+100, 880) var, lon, lat, lev, zob, dummy,count,time,type
     aaa= lon*180.0/pi
     if(aaa.gt.360.0) print*, 'wr lon = ', aaa,lon
    enddo

  200 continue
  666 format(4e20.13)
  880 format(f4.2, 2f7.3, f7.1, f7.2, f7.2, f9.0, f7.3, f5.0)

end subroutine num_ncep_obs
!========================================

end program create_ncep_obs_set_def
