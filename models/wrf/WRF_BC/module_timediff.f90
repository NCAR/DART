MODULE module_timediff

CONTAINS

SUBROUTINE time_diff(stime,etime,diff)
!#######################################
! computes the difference in seconds between stime and etime, where stime
! and etime are character*80 as: YYYY-MM-DD_hh:mm:ss
! NOTE: assumes that the difference is less than a day (should be fine for this)

  implicit none

  character(len=80), intent(in) :: stime, etime
  real, intent(out) :: diff

  integer :: shr,smm,sss, s_secs
  integer :: ehr,emm,ess, e_secs

  read(stime(12:13),*) shr
  read(stime(15:16),*) smm
  read(stime(18:19),*) sss
  read(etime(12:13),*) ehr
  read(etime(15:16),*) emm
  read(etime(18:19),*) ess

  s_secs = shr*3600 + smm*60 +sss
  e_secs = ehr*3600 + emm*60 +ess

  diff = e_secs - s_secs

  if ( diff < 0 ) then
    print*, "your time difference is negative - aborting:"
    stop 'time_diff'
  endif

END SUBROUTINE time_diff

SUBROUTINE find_time_index(timelist, time_to_find, ntimes, itime)
!################################################
! finds first index in timelist such that tim_to_find is later than
! timelist(itime)
!##################################################
  implicit none

  character(len=80), dimension(:), intent(in) :: timelist
  character(len=80), intent(in) :: time_to_find
  integer, intent(in) :: ntimes
  integer, intent(out) :: itime

  integer :: it
  integer*8 :: tfyr, tfmo, tfdy, tfhr, tfmm, tfss, tf
  integer*8 :: lyr, lmo, ldy, lhr, lmm, lss, l
  
  read(time_to_find(1:4),*) tfyr
  read(time_to_find(6:7),*) tfmo
  read(time_to_find(9:10),*) tfdy
  read(time_to_find(12:13),*) tfhr
  read(time_to_find(15:16),*) tfmm
  read(time_to_find(18:19),*) tfss
  
  tf = (tfyr*10000000000)+(tfmo*100000000) + (tfdy*1000000) + (tfhr*10000) + (tfmm*100) + tfss

  itime = ntimes
  do it = ntimes,1,-1
    read(timelist(it)(1:4),*) lyr
    read(timelist(it)(6:7),*) lmo
    read(timelist(it)(9:10),*) ldy
    read(timelist(it)(12:13),*) lhr
    read(timelist(it)(15:16),*) lmm
    read(timelist(it)(18:19),*) lss
    l = (lyr*10000000000)+(lmo*100000000) + (ldy*1000000) + (lhr*10000) + (lmm*100) + lss
    if ( l <= tf ) then
      itime = it 
      exit
    endif
  enddo 
   
END SUBROUTINE find_time_index

END MODULE module_timediff
