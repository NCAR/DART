module radinfo
! adapted from GSI/src/main/radinfo.f90 to include only minimum code for reading channel info

use         kinds, only : i_kind
use utilities_mod, only : open_file

implicit none

! all variables and subroutines are public

! public variables
integer(i_kind),allocatable,dimension(:) :: iuse_rad ! use to turn off satellite radiance data
character(len=20),allocatable,dimension(:):: nusis   ! sensor/instrument/satellite indicator
integer(i_kind) :: jpch_rad ! number of channels*sat
integer(i_kind),allocatable,dimension(:):: nuchan ! satellite channel

! namelist variables
integer(i_kind) :: npred = 12       ! number of radiance biases predictors
logical :: adp_anglebc   = .false.  ! logical to turn off or on the variational radiance angle bias correction
logical :: emiss_bc      = .false.  ! logical to turn off or on the emissivity predictor

contains

  subroutine radinfo_read

    implicit none

    integer(i_kind) i,j,k,lunin,nlines
    integer(i_kind) istat,n
    character(len=1):: cflg
    character(len=120) crecord

    data lunin / 49 /

!============================================================================

!   Determine number of entries in satellite information file
    lunin = open_file('satinfo',form='formatted',action='read')

    j=0
    nlines=0
    read1:  do
       read(lunin,100,iostat=istat) cflg,crecord
       if (istat /= 0) exit read1
       nlines=nlines+1
       if (cflg == '!') cycle read1
       j=j+1
    end do read1
    if (istat>0) then
       close(lunin)
       write(6,*)'RADINFO_READ:  ***ERROR*** error reading radinfo, istat=',istat
       write(6,*)'RADINFO_READ:  stop program execution'
       stop
    endif
    jpch_rad = j
    if(jpch_rad == 0)then
      close(lunin)
       write(6,*)'RADINFO_READ:  no channels from satinfo'
       write(6,*)'RADINFO_READ:  no channels from satinfo'
      return
    end if

!   Allocate arrays to hold radiance information
!     nuchan    - channel number
!     nusis     - sensor/instrument/satellite
!     iuse_rad  - use parameter

    allocate(nuchan(jpch_rad),nusis(jpch_rad),iuse_rad(0:jpch_rad))
    iuse_rad(0)=-999

    rewind(lunin)
    j=0
    do k=1,nlines
       read(lunin,100) cflg,crecord
       if (cflg == '!') cycle
       j=j+1
       read(crecord,*,iostat=istat) nusis(j),nuchan(j),iuse_rad(j)
       if(istat/=0) then
          write(6,*)'RADINFO_READ:  ***ERROR*** error reading satinfo, crecord = ',trim(crecord)
          write(6,*)'RADINFO_READ:  ***ERROR*** error reading satinfo, istat   = ',istat
          write(6,*)'RADINFO_READ:  stop program execution'
          stop
       endif
    end do
    close(lunin)
100 format(a1,a120)

    return

  end subroutine radinfo_read

  subroutine radinfo_clean
     if( allocated(nuchan))    deallocate(nuchan)
     if( allocated(nusis))     deallocate(nusis)
     if( allocated(iuse_rad))  deallocate(iuse_rad)
  end subroutine radinfo_clean

end module radinfo
