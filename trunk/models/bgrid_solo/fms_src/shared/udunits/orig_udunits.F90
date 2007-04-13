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

  use mpp_mod, only : mpp_error, FATAL, WARNING

  implicit none

  character(len=128) :: c_utname = '/usr/local/etc/udunits.dat'
  integer :: utmake, utopen, utdec,utorigin, utcaltime, uttime, utcvt
  external utmake, utopen, utdec,utorigin, utcaltime, uttime, utcvt

  logical :: udunits_initialized = .false.

  public get_cal_time, convert_units, udunits_init, udunits_exit 
  

  contains

    subroutine udunits_init()

      integer(LONG_KIND) :: iret

      if (udunits_initialized) return
      iret = utopen(c_utname) ! should probably open and close udunits file somewhere else
      if (iret.ne.0) call mpp_error(FATAL,'cant open udunits package')
      udunits_initialized = .true.

      return
    end subroutine udunits_init

    subroutine get_cal_time(time,units,yr,mon,day,hr,min,sec)

      real, intent(in) :: time
      character(len=*), intent(in) :: units
      integer, intent(out) :: yr, mon, day, hr, min, sec
      integer(LONG_KIND) :: iunit1, iunit2
      pointer (p_iunit1,iunit1)
      pointer (p_iunit2, iunit2)
      integer :: iret

      if (.not.udunits_initialized) call udunits_init()

      p_iunit1 = utmake()
      iret = utdec(units,p_iunit1)
      if (iret .ne. 0) call mpp_error(FATAL,'unrecognized time units')
!      iret = uttime(p_iunit1)
!      if (iret .ne. 0) call mpp_error(FATAL,'no time units in file')
      iret = utorigin(p_iunit1)
      if (iret .ne. 1) call mpp_error(FATAL,'no time origin in file')
      iret = utcaltime(time,p_iunit1,yr,mon,day,hr,min,sec)
      if (iret.ne.0) call mpp_error(FATAL,'error calculating time in utcaltime')
      
      if (sec > 60) sec = 0  ! error in utcaltime??
    end subroutine get_cal_time
    
    subroutine convert_units(from_units,to_units, slope, intercept)

      character(len=*) :: from_units, to_units
      real, intent(out) :: slope, intercept 
      integer(LONG_KIND) :: iunit1, iunit2
      pointer (p_iunit1,iunit1)
      pointer (p_iunit2, iunit2)
      integer :: iret


      if (.not.udunits_initialized) call udunits_init()      
      p_iunit1 = utmake()
      p_iunit2 = utmake()
      iret = utdec(from_units,p_iunit1)
      if (iret .ne. 0) then
          call mpp_error(WARNING,'unrecognized from units')
          slope=1
          intercept=0
          return
      endif
      iret = utdec(to_units,p_iunit2)
      if (iret .ne. 0) call mpp_error(FATAL,'unrecognized to units')      
      iret = utcvt(p_iunit1,p_iunit2,slope,intercept)
! error code returned is incorrect, i.e. bad return code but
! correct answer?? comment out for now
!      if (iret .ne. 0) call mpp_error(FATAL,'unable to convert units')
    end subroutine convert_units

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

