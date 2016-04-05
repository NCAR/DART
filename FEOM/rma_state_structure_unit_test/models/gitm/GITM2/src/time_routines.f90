! This code may (or may not) be part of the GITM distribution,
! So it is not protected by the DART copyright agreement.
!
! DART $Id$

integer function jday(year, mon, day) result(Julian_Day)

  implicit none

  integer :: i
  integer, dimension(1:12) :: dayofmon
  integer :: year, mon, day

  dayofmon(1) = 31
  dayofmon(2) = 28
  dayofmon(3) = 31
  dayofmon(4) = 30
  dayofmon(5) = 31
  dayofmon(6) = 30
  dayofmon(7) = 31
  dayofmon(8) = 31
  dayofmon(9) = 30
  dayofmon(10) = 31
  dayofmon(11) = 30
  dayofmon(12) = 31

  if (mod(year,4).eq.0) dayofmon(2) = dayofmon(1) + 1
  Julian_Day = 0
  do i = 1, mon-1
     Julian_Day = Julian_Day + dayofmon(i)
  enddo
  Julian_Day = Julian_Day + day

end

subroutine time_int_to_real(itime, timereal)

  implicit none

  integer, dimension(1:12) :: dayofmon
  integer, dimension(1:7) :: itime
  double precision :: timereal
  integer :: nyear, nleap, nmonth, nday, nhour, nmin, nsec, i

  dayofmon(1) = 31
  dayofmon(2) = 28
  dayofmon(3) = 31
  dayofmon(4) = 30
  dayofmon(5) = 31
  dayofmon(6) = 30
  dayofmon(7) = 31
  dayofmon(8) = 31
  dayofmon(9) = 30
  dayofmon(10) = 31
  dayofmon(11) = 30
  dayofmon(12) = 31

  if (mod(itime(1),4).eq.0) dayofmon(2) = 29

  timereal = 0.0;

  if (itime(1) > 1900) then
     nyear = itime(1) - 1965
  else 
     if (itime(1) > 65) then
        nyear = itime(1) - 65 
     else 
        nyear = itime(1) + 100 - 65
     endif
  endif
  nleap = nyear/4

  nmonth = itime(2) - 1

  nday = 0

  do i=1, nmonth
     nday = nday + dayofmon(i)
  enddo

  nday = nday + itime(3) - 1
  nhour = itime(4)
  nmin = itime(5)
  nsec = itime(6)

  timereal = (dble(nsec) * dble(1.0)) +                  &
       (dble(nmin) * dble(60.0)) +                       &
       (dble(nhour) * dble(60.0*60.0)) +                 &
       (dble(nday) * dble(24.0*60.0*60.0)) +             &
       (dble(nleap) * dble(24.0*60.0*60.0)) +            &
       (dble(nyear) * dble(365.0*24.0*60.0*60.0)) +      &
       itime(7)/1000.0

end subroutine time_int_to_real



subroutine time_real_to_int(timereal, itime)

  implicit none

  integer, dimension(1:12) :: dayofmon
  integer, dimension(1:7) :: itime
  double precision :: timereal
  integer :: nyear, nleap, nmonth, nday, nhour, nmin, nsec
  double precision :: speryear = 31536000.0
  double precision :: sperday = 86400.0
  double precision :: sperhour = 3600.0
  double precision :: spermin = 60.0
  double precision :: timeleft

  dayofmon(1) = 31
  dayofmon(2) = 28
  dayofmon(3) = 31
  dayofmon(4) = 30
  dayofmon(5) = 31
  dayofmon(6) = 30
  dayofmon(7) = 31
  dayofmon(8) = 31
  dayofmon(9) = 30
  dayofmon(10) = 31
  dayofmon(11) = 30
  dayofmon(12) = 31

  nyear = int(timereal/speryear)
  nleap = nyear/4
  nday = int((timereal - (dble(nyear)*speryear))/sperday)

  if (nday.le.nleap) then
     nyear = int((timereal - (dble(nleap)*sperday))/speryear)
     nleap = nyear/4
     nday = int((timereal - (dble(nyear)*speryear))/sperday)
     if (nday.le.nleap) then
        nyear = int((timereal - (dble(nleap)*sperday))/speryear)
        nleap = nyear/4
        nday = int((timereal - (dble(nyear)*speryear))/sperday)
     endif
  endif

  if (mod((nyear+65),4).eq.0) dayofmon(2) = dayofmon(2) + 1

  nday = nday - nleap

  timeleft = timereal - dble(nyear)*speryear
  timeleft = timeleft - dble(nday+nleap)*sperday

  nhour = int(timeleft/sperhour)
  timeleft = timeleft - dble(nhour)*sperhour

  nmin = int(timeleft/spermin)
  timeleft = timeleft - dble(nmin)*spermin

  nsec = int(timeleft)

  nmonth = 1;

  do while (nday.ge.dayofmon(nmonth))
     nday = nday - dayofmon(nmonth)
     nmonth = nmonth + 1
  end do

  itime(1) = nyear + 1965
  itime(2) = nmonth
  itime(3) = nday + 1
  itime(4) = nhour
  itime(5) = nmin
  itime(6) = nsec
  itime(7) = (timeleft - nsec) * 1000

end subroutine time_real_to_int

! <next few lines under version control, do not edit>
! $URL$
! $Id$
! $Revision$
! $Date$
