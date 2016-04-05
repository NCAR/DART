! This code is not protected by the DART copyright agreement.
! DART $Id$

! This is the only file from the rose source distribution needed
! for the DART interfaces. It is included only so that we can build
! the DART routines as part of our code maintenance. - TJH

      module params

      implicit none

!... general model 

      integer, parameter :: nx = 32, ny = 36, nz = 38
      integer, parameter :: nyp1=37, nzp1=nz+1, nzp2=nz+2
      integer, parameter :: nxhlf=nx/2

!... radiation

      integer, parameter :: nztrop=7, nzz=45

!... chemistry

      integer, parameter :: nbcon = 26
      integer, parameter :: nphot = 33 

!... dynamics

      integer, parameter :: ndyn = 4 
  
      integer, parameter :: nvars = ndyn + nbcon

!... variables for program control

      integer :: ntime                    ! number of steps per hour

      integer :: year0                    ! start year
      integer :: day0                     ! start day of year (day 1 = Jan. 1)
      integer :: ut0                      ! start universal time (sec) 

      integer :: nstep, nend, nsave, nseg, nstart 
      integer :: nout, nouttid, noutdiag

      save

      end module params 

! <next few lines under version control, do not edit>
! $URL$
! $Id$
! $Revision$
! $Date$
