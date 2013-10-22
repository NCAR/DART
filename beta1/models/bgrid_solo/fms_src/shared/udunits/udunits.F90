!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!                                                                   !!
!!                   GNU General Public License                      !!
!!                                                                   !!
!! This file is part of the Flexible Modeling System (FMS).          !!
!!                                                                   !!
!! FMS is free software; you can redistribute it and/or modify       !!
!! it and are expected to follow the terms of the GNU General Public !!
!! License as published by the Free Software Foundation.             !!
!!                                                                   !!
!! FMS is distributed in the hope that it will be useful,            !!
!! but WITHOUT ANY WARRANTY; without even the implied warranty of    !!
!! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the     !!
!! GNU General Public License for more details.                      !!
!!                                                                   !!
!! You should have received a copy of the GNU General Public License !!
!! along with FMS; if not, write to:                                 !!
!!          Free Software Foundation, Inc.                           !!
!!          59 Temple Place, Suite 330                               !!
!!          Boston, MA  02111-1307  USA                              !!
!! or see:                                                           !!
!!          http://www.gnu.org/licenses/gpl.txt                      !!
!!                                                                   !!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
#include <os.h>

module udunits_mod
!
!<CONTACT EMAIL="mjh@gfdl.noaa.gov">M.J. Harrison</CONTACT>
!
!<REVIEWER EMAIL="wfc@gfdl.noaa.gov">Will Cooke</REVIEWER>
!

!</OVERVIEW>
! Interfaces to udunits library (http://www.unidata.ucar.edu/packages/udunits)
!</OVERVIEW>

!<DESCRIPTION>
!
! subroutine udunits_init : Initialize udunits module
! subroutine get_cal_time : Retrieve calandar time based on units string
! subroutine convert_units : convert units <to_units> = <from_units>*slope + intercept
! subroutinen udunits_exit : Exit udunits module 
!</DESCRIPTION>
!
!<NOTES>
! tested succesfully with sgi mipspro compiler
! pgf90: problem with -I compilation flag,  not able to open udunits 
!        database without specifying full path
!</NOTES>
!
#ifdef __aix 
#define utopen utopen_ 
#define utdec utdec_ 
#define uttime  uttime_ 
#define utorigin utorigin_ 
#define utcaltime utcaltime_ 
#define utmake utmake_
#define utcvt utcvt_
#define utcls utcls_
#endif

  !use mpp_mod, only : mpp_error, FATAL, WARNING
  use fms_mod, only : error_mesg, FATAL, WARNING

  implicit none

  character(len=128) :: c_utname = '/usr/local/etc/udunits.dat'
  integer :: utmake, utopen, utdec,utorigin, utcaltime, uttime, utcvt
  external utmake, utopen, utdec,utorigin, utcaltime, uttime, utcvt

  logical :: udunits_initialized = .false.

  public udunits_init, udunits_exit 
  

  contains

    subroutine udunits_init()

      integer(LONG_KIND) :: iret

      if (udunits_initialized) return
      iret = utopen(c_utname) ! should probably open and close udunits file somewhere else
      if (iret.ne.0) call error_mesg('in udunits_mod udunints_init','cant open udunits package', FATAL)
      udunits_initialized = .true.

      return
    end subroutine udunits_init


    subroutine udunits_exit()

      call utcls()

      return
      
    end subroutine udunits_exit

  end module udunits_mod
  
#ifdef test_udunits

program test

use udunits_mod

implicit none

character(len=32) :: time_units='days since 1900-01-01 00:00:00'

real :: slope, intercept, time=36524.

integer :: yr, mon, day, hr, min, sec

call udunits_init()

write(*,*) 'Opened udunits module successfully! ...'

call convert_units('deg_C','deg_K',slope, intercept)

write(*,'(a,f7.3,a,f7.3)') '<deg_K> = ',slope,' * <deg_C> + ',intercept

call get_cal_time(time,time_units, yr, mon, day, hr, min, sec)

write(*,'(f7.1,1x,a,a,i4,a,i2,a,i2,1x,i2,a,i2,a,i2)') time,trim(time_units),' = ',yr,'-',mon,'-',day,hr,':',min,':',sec

call udunits_exit()

stop

end program test

#endif

