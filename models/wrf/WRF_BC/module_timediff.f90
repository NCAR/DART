! DART software - Copyright UCAR. This open source software is provided
! by UCAR, "as is", without charge, subject to all terms of use at
! http://www.image.ucar.edu/DAReS/DART/DART_download
!
! $Id$

MODULE module_timediff

  use        types_mod, only : r8, i8
  use time_manager_mod, only : time_type, set_date_gregorian, get_time

  implicit none
  private

public :: time_diff, find_time_index

CONTAINS

SUBROUTINE time_diff(stime,etime,diff)
!#######################################
! computes the difference in seconds between stime and etime, where stime
! and etime are character*19 as: YYYY-MM-DD_hh:mm:ss

  implicit none

  character(len=19), intent(in)  :: stime, etime
  real(r8),          intent(out) :: diff

  type(time_type) :: sdate, edate
  integer         :: syy,smm,sdd,shr,smi,sss
  integer         :: eyy,emm,edd,ehr,emi,ess
  integer(i8)     :: s_secs, e_secs

  read(stime(1:4),*) syy
  read(stime(6:7),*) smm
  read(stime(9:10),*) sdd
  read(stime(12:13),*) shr
  read(stime(15:16),*) smi
  read(stime(18:19),*) sss
  read(etime(1:4),*) eyy
  read(etime(6:7),*) emm
  read(etime(9:10),*) edd
  read(etime(12:13),*) ehr
  read(etime(15:16),*) emi
  read(etime(18:19),*) ess

  sdate = set_date_gregorian(syy,smm,sdd,shr,smi,sss)
  edate = set_date_gregorian(eyy,emm,edd,ehr,emi,ess)

!!$  s_secs = shr*3600 + smi*60 +sss
!!$  e_secs = ehr*3600 + emi*60 +ess
  call get_time(sdate, sss, sdd)
  call get_time(edate, ess, edd)
  s_secs = sss + sdd*86400_i8
  e_secs = ess + edd*86400_i8

  diff = e_secs - s_secs

  if ( diff < 0 ) then
     print*,sss,ess,sdd,edd,e_secs,s_secs,diff
     print*, "your time difference is negative - aborting:"
     stop 'time_diff'
  endif

END SUBROUTINE time_diff

SUBROUTINE find_time_index(timelist, time_to_find, ntimes, itime)
!################################################
! finds first index in timelist such that time_to_find is later than
! timelist(itime)
!##################################################
  implicit none

  character(len=19), dimension(:), intent(in)  :: timelist
  character(len=19),               intent(in)  :: time_to_find
  integer,                         intent(in)  :: ntimes
  integer,                         intent(out) :: itime

  integer          :: it
  integer(kind=i8) :: tfyr, tfmo, tfdy, tfhr, tfmm, tfss, tf
  integer(kind=i8) :: lyr, lmo, ldy, lhr, lmm, lss, l
  
  read(time_to_find(1:4),*) tfyr
  read(time_to_find(6:7),*) tfmo
  read(time_to_find(9:10),*) tfdy
  read(time_to_find(12:13),*) tfhr
  read(time_to_find(15:16),*) tfmm
  read(time_to_find(18:19),*) tfss
  
  tf = (tfyr*10000000000_i8)+(tfmo*100000000_i8) + (tfdy*1000000_i8) + &
       (tfhr*10000_i8) + (tfmm*100_i8) + tfss

  itime = ntimes
  do it = ntimes,1,-1
    read(timelist(it)(1:4),*) lyr
    read(timelist(it)(6:7),*) lmo
    read(timelist(it)(9:10),*) ldy
    read(timelist(it)(12:13),*) lhr
    read(timelist(it)(15:16),*) lmm
    read(timelist(it)(18:19),*) lss
    l = (lyr*10000000000_i8)+(lmo*100000000_i8) + (ldy*1000000_i8) + &
         (lhr*10000_i8) + (lmm*100_i8) + lss
    if ( l <= tf ) then
      itime = it 
      exit
    endif
  enddo 
   
END SUBROUTINE find_time_index

END MODULE module_timediff

! <next few lines under version control, do not edit>
! $URL$
! $Id$
! $Revision$
! $Date$
