! DART software - Copyright UCAR. This open source software is provided
! by UCAR, "as is", without charge, subject to all terms of use at
! http://www.image.ucar.edu/DAReS/DART/DART_download
!
! $Id$
 
program advance_cymdh

  implicit none

! version controlled file description for error handling, do not edit
character(len=256), parameter :: source   = &
   "$URL$"
character(len=32 ), parameter :: revision = "$Revision$"
character(len=128), parameter :: revdate  = "$Date$"

  integer :: ccyy, mm, dd, hh, dh

  integer :: nargum, i, n, sign

  character(len=80), dimension(2) :: argum

  character(len=10) :: ccyymmddhh

  nargum=COMMAND_ARGUMENT_COUNT()

  if(nargum /= 2) then
     write(unit=*, fmt='(a)') &
          'Usage: advance_cymdh ccyymmddhh dh'
     stop 'try again.'
  endif

  do i=1,nargum
     do n=1,80
        argum(i)(n:n)=' '
     enddo
     call GET_COMMAND_ARGUMENT(i,argum(i))
  enddo

  ccyymmddhh = trim(argum(1))

  read(ccyymmddhh(1:10), fmt='(i4, 3i2)')  ccyy, mm, dd, hh

  sign = 1

  dh = 0

  do n=1,len_trim(argum(2))
     if(argum(2)(n:n) == '-') then
        sign = -1
        cycle
     else
        read(argum(2)(n:n), fmt='(i1)') i
        dh=10*dh + i
     end if
  enddo

  dh = sign * dh

  hh = hh + dh

  do while (hh < 0) 
     hh = hh + 24
     call change_date ( ccyy, mm, dd, -1 )
  end do

  do while (hh > 23) 
     hh = hh - 24
     call change_date ( ccyy, mm, dd, 1 )
  end do

  write(ccyymmddhh(1:10), fmt='(i4, 3i2.2)')  ccyy, mm, dd, hh
  write(unit=*, fmt='(a)') ccyymmddhh

contains

  subroutine change_date( ccyy, mm, dd, delta )

    implicit none

    integer, intent(inout) :: ccyy, mm, dd
    integer, intent(in)    :: delta

    integer, dimension(12) :: mmday

    mmday = (/31,28,31,30,31,30,31,31,30,31,30,31/)

    if (mod(ccyy,4) == 0) then
       mmday(2) = 29

       if ( mod(ccyy,100) == 0) then
          mmday(2) = 28
       endif

       if(mod(ccyy,400) == 0) then
          mmday(2) = 29
       end if
    endif

    dd = dd + delta

    if(dd == 0) then
       mm = mm - 1

       if(mm == 0) then
          mm = 12
          ccyy = ccyy - 1
       endif

       dd = mmday(mm)
    elseif ( dd .gt. mmday(mm) ) then
       dd = 1
       mm = mm + 1
       if(mm > 12 ) then
          mm = 1
          ccyy = ccyy + 1
       end if
    end if

  end subroutine change_date

end program advance_cymdh

! <next few lines under version control, do not edit>
! $URL$
! $Id$
! $Revision$
! $Date$
