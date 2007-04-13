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

module fft99_mod

use types_mod, only : r8
use constants_mod, only: pi

implicit none
private

public :: fft99, fft991, set99

contains

!##########################################################################

    subroutine fft99 (a,work,trigs,ifax,inc,jump,n,lot,isign) 

! purpose      performs multiple fast fourier transforms.  this package
!              will perform a number of simultaneous real/half-complex
!              periodic fourier transforms or corresponding inverse
!              transforms, i.e.  given a set of real data vectors, the
!              package returns a set of 'half-complex' fourier
!              coefficient vectors, or vice versa.  the length of the
!              transforms must be an even number greater than 4 that has
!              no other factors except possibly powers of 2, 3, and 5.
!              this is an all-fortran version of a optimized routine
!              fft99 written for xmp/ymps by dr. clive temperton of
!              ecmwf.
!
!              the package fft99f contains several user-level routines:
!
!            subroutine set99
!                an initialization routine that must be called once
!                before a sequence of calls to the fft routines
!                (provided that n is not changed).
!
!            subroutines fft99 and fft991
!                two fft routines that return slightly different
!                arrangements of the data in gridpoint space.
!
! usage        let n be of the form 2**p * 3**q * 5**r, where p .ge. 1,
!              q .ge. 0, and r .ge. 0.  then a typical sequence of
!              calls to transform a given set of real vectors of length
!              n to a set of 'half-complex' fourier coefficient vectors
!              of length n is
!
!                   dimension ifax(13),trigs(3*n/2+1),a(m*(n+2)),
!                  +          work(m*(n+1))
!
!                   call set99 (trigs, ifax, n)
!                   call fft99 (a,work,trigs,ifax,inc,jump,n,m,isign)
!
!              see the individual write-ups for set99, fft99, and
!              fft991 below, for a detailed description of the
!              arguments.
!
! history      the package was written by clive temperton at ecmwf in
!              november, 1978.  it was modified, documented, and tested
!              for ncar by russ rew in september, 1980.
!
!-----------------------------------------------------------------------
!
! subroutine set99 (trigs, ifax, n)
!
! purpose      a set-up routine for fft99 and fft991.  it need only be
!              called once before a sequence of calls to the fft
!              routines (provided that n is not changed).
!
! argument     ifax(13),trigs(3*n/2+1)
! dimensions
!
! arguments
!
! on input     trigs
!               a floating point array of dimension 3*n/2 if n/2 is
!               even, or 3*n/2+1 if n/2 is odd.
!
!              ifax
!               an integer array.  the number of elements actually used
!               will depend on the factorization of n.  dimensioning
!               ifax for 13 suffices for all n less than a million.
!
!              n
!               an even number greater than 4 that has no prime factor
!               greater than 5.  n is the length of the transforms (see
!               the documentation for fft99 and fft991 for the
!               definitions of the transforms).
!
! on output    ifax
!               contains the factorization of n/2.  ifax(1) is the
!               number of factors, and the factors themselves are stored
!               in ifax(2),ifax(3),...  if set99 is called with n odd,
!               or if n has any prime factors greater than 5, ifax(1)
!               is set to -99.
!
!              trigs
!               an array of trigonometric function values subsequently
!               used by the fft routines.
!
!-----------------------------------------------------------------------
!
! subroutine fft991 (a,work,trigs,ifax,inc,jump,n,m,isign)
!                       and
! subroutine fft99 (a,work,trigs,ifax,inc,jump,n,m,isign)
!
! purpose      perform a number of simultaneous real/half-complex
!              periodic fourier transforms or corresponding inverse
!              transforms, using ordinary spatial order of gridpoint
!              values (fft991) or explicit cyclic continuity in the
!              gridpoint values (fft99).  given a set
!              of real data vectors, the package returns a set of
!              'half-complex' fourier coefficient vectors, or vice
!              versa.  the length of the transforms must be an even
!              number that has no other factors except possibly powers
!              of 2, 3, and 5.  this is an all-fortran version of 
!              optimized routine fft991 written for xmp/ymps by
!              dr. clive temperton of ecmwf.
!
! argument     a(m*(n+2)), work(m*(n+1)), trigs(3*n/2+1), ifax(13)
! dimensions
!
! arguments
!
! on input     a
!               an array of length m*(n+2) containing the input data
!               or coefficient vectors.  this array is overwritten by
!               the results.
!
!              work
!               a work array of dimension m*(n+1)
!
!              trigs
!               an array set up by set99, which must be called first.
!
!              ifax
!               an array set up by set99, which must be called first.
!
!              inc
!               the increment (in words) between successive elements of
!               each data or coefficient vector (e.g.  inc=1 for
!               consecutively stored data).
!
!              jump
!               the increment (in words) between the first elements of
!               successive data or coefficient vectors.  on crays, 
!               try to arrange data so that jump is not a multiple of 8
!               (to avoid memory bank conflicts).  for clarification of
!               inc and jump, see the examples below.
!
!              n
!               the length of each transform (see definition of
!               transforms, below).
!
!              m
!               the number of transforms to be done simultaneously.
!
!              isign
!               = +1 for a transform from fourier coefficients to
!                    gridpoint values.
!               = -1 for a transform from gridpoint values to fourier
!                    coefficients.
!
! on output    a
!               if isign = +1, and m coefficient vectors are supplied
!               each containing the sequence:
!
!               a(0),b(0),a(1),b(1),...,a(n/2),b(n/2)  (n+2 values)
!
!               then the result consists of m data vectors each
!               containing the corresponding n+2 gridpoint values:
!
!               for fft991, x(0), x(1), x(2),...,x(n-1),0,0.
!               for fft99, x(n-1),x(0),x(1),x(2),...,x(n-1),x(0).
!                   (explicit cyclic continuity)
!
!               when isign = +1, the transform is defined by:
!                 x(j)=sum(k=0,...,n-1)(c(k)*exp(2*i*j*k*pi/n))
!                 where c(k)=a(k)+i*b(k) and c(n-k)=a(k)-i*b(k)
!                 and i=sqrt (-1)
!
!               if isign = -1, and m data vectors are supplied each
!               containing a sequence of gridpoint values x(j) as
!               defined above, then the result consists of m vectors
!               each containing the corresponding fourier cofficients
!               a(k), b(k), 0 .le. k .le n/2.
!
!               when isign = -1, the inverse transform is defined by:
!                 c(k)=(1/n)*sum(j=0,...,n-1)(x(j)*exp(-2*i*j*k*pi/n))
!                 where c(k)=a(k)+i*b(k) and i=sqrt(-1)
!
!               a call with isign=+1 followed by a call with isign=-1
!               (or vice versa) returns the original data.
!
!               note: the fact that the gridpoint values x(j) are real
!               implies that b(0)=b(n/2)=0.  for a call with isign=+1,
!               it is not actually necessary to supply these zeros.
!
! examples      given 19 data vectors each of length 64 (+2 for explicit
!               cyclic continuity), compute the corresponding vectors of
!               fourier coefficients.  the data may, for example, be
!               arranged like this:
!
! first data   a(1)=    . . .                a(66)=             a(70)
! vector       x(63) x(0) x(1) x(2) ... x(63) x(0)  (4 empty locations)
!
! second data  a(71)=   . . .                                  a(140)
! vector       x(63) x(0) x(1) x(2) ... x(63) x(0)  (4 empty locations)
!
!               and so on.  here inc=1, jump=70, n=64, m=19, isign=-1,
!               and fft99 should be used (because of the explicit cyclic
!               continuity).
!
!               alternatively the data may be arranged like this:
!
!                first         second                          last
!                data          data                            data
!                vector        vector                          vector
!
!                 a(1)=         a(2)=                           a(19)=
!
!                 x(63)         x(63)       . . .               x(63)
!        a(20)=   x(0)          x(0)        . . .               x(0)
!        a(39)=   x(1)          x(1)        . . .               x(1)
!                  .             .                               .
!                  .             .                               .
!                  .             .                               .
!
!               in which case we have inc=19, jump=1, and the remaining
!               parameters are the same as before.  in either case, each
!               coefficient vector overwrites the corresponding input
!               data vector.
!
!-----------------------------------------------------------------------
    integer, intent(in)    :: inc,jump,n,lot,isign
    integer, intent(inout) :: ifax(:)
    real(r8),    intent(in)    :: trigs(:)
    real(r8),    intent(inout) :: a(*),work(*)

!     dimension a(n),work(n),trigs(n),ifax(1)
!
!     subroutine "fft99" - multiple fast real periodic transform
!     corresponding to old scalar routine fft9
!     procedure used to convert to half-length complex transform
!     is given by cooley, lewis and welch (j. sound vib., vol. 12
!     (1970), 315-337)
!
!     a is the array containing input and output data
!     work is an area of size (n+1)*lot
!     trigs is a previously prepared list of trig function values
!     ifax is a previously prepared list of factors of n/2
!     inc is the increment within each data 'vector'
!         (e.g. inc=1 for consecutively stored data)
!     jump is the increment between the start of each data vector
!     n is the length of the data vectors
!     lot is the number of data vectors
!     isign = +1 for transform from spectral to gridpoint
!           = -1 for transform from gridpoint to spectral
!
!     ordering of coefficients:
!         a(0),b(0),a(1),b(1),a(2),b(2),...,a(n/2),b(n/2)
!         where b(0)=b(n/2)=0; (n+2) locations required
!
!     ordering of data:
!         x(n-1),x(0),x(1),x(2),...,x(n),x(0)
!         i.e. explicit cyclic continuity; (n+2) locations required
!
!     vectorization is achieved on cray by doing the transforms in
!     parallel
!
!     *** n.b. n is assumed to be an even number
!
!     definition of transforms:
!     -------------------------
!
!     isign=+1: x(j)=sum(k=0,...,n-1)(c(k)*exp(2*i*j*k*pi/n))
!         where c(k)=a(k)+i*b(k) and c(n-k)=a(k)-i*b(k)
!
!     isign=-1: a(k)=(1/n)*sum(j=0,...,n-1)(x(j)*cos(2*j*k*pi/n))
!               b(k)=-(1/n)*sum(j=0,...,n-1)(x(j)*sin(2*j*k*pi/n))
!

    integer :: nfax, nx, nh, ink, igo, ibase, jbase
    integer :: i, j, k, L, m, ia, la, ib

      nfax=ifax(1)
      nx=n+1
      nh=n/2
      ink=inc+inc
      if (isign.eq.+1) go to 30

!   if necessary, transfer data to work area
      igo=50
      if (mod(nfax,2).eq.1) goto 40
      ibase=inc+1
      jbase=1
    do L=1,lot
      i=ibase
      j=jbase
!dir$ ivdep
      do m=1,n
        work(j)=a(i)
        i=i+inc
        j=j+1
      enddo
      ibase=ibase+jump
      jbase=jbase+nx
    enddo

      igo=60
      go to 40

!   preprocessing (isign=+1)
!   ------------------------

   30 continue
      call fft99a(a,work,trigs,inc,jump,n,lot)
      igo=60

!   complex transform
!   -----------------

   40 continue
      ia=inc+1
      la=1
      do 80 k=1,nfax
      if (igo.eq.60) go to 60
   50 continue
      call vpassm (a(ia),a(ia+inc),work(1),work(2),trigs, &
                   ink,2,jump,nx,lot,nh,ifax(k+1),la)
      igo=60
      go to 70
   60 continue
      call vpassm (work(1),work(2),a(ia),a(ia+inc),trigs, &
                   2,ink,nx,jump,lot,nh,ifax(k+1),la)
      igo=50
   70 continue
      la=la*ifax(k+1)
   80 continue

    if (isign.eq.-1) go to 130

! if necessary, transfer data from work area

    if (mod(nfax,2).ne.1) then
      ibase=1
      jbase=ia
      do L=1,lot
        i=ibase
        j=jbase
!dir$ ivdep
        do m=1,n
          a(j)=work(i)
          i=i+1
          j=j+inc
        enddo
        ibase=ibase+nx
        jbase=jbase+jump
      enddo
    endif

!   fill in cyclic boundary points
      ia=1
      ib=n*inc+1
!dir$ ivdep
      do L=1,lot
        a(ia)=a(ib)
        a(ib+inc)=a(ia+inc)
        ia=ia+jump
        ib=ib+jump
      enddo
      go to 140

!   postprocessing (isign=-1):
!   --------------------------

  130 continue
      call fft99b(work,a,trigs,inc,jump,n,lot)

  140 continue

    end subroutine fft99

!##########################################################################

    subroutine fft99a (a,work,trigs,inc,jump,n,lot)
    integer, intent(in)    :: inc,jump,n,lot
    real(r8),    intent(in)    :: trigs(:)
    real(r8),    intent(inout) :: a(*),work(*)

!     dimension a(n),work(n),trigs(n)
!
!     subroutine fft99a - preprocessing step for fft99, isign=+1
!     (spectral to gridpoint transform)

    integer :: nh, nx, ink, k, L
    integer :: ia, ib, ja, jb, iabase, ibbase, jabase, jbbase
    real(r8)    :: c, s

      nh=n/2
      nx=n+1
      ink=inc+inc

!   a(0) and a(n/2)
      ia=1
      ib=n*inc+1
      ja=1
      jb=2
!dir$ ivdep
    do L=1,lot
      work(ja)=a(ia)+a(ib)
      work(jb)=a(ia)-a(ib)
      ia=ia+jump
      ib=ib+jump
      ja=ja+nx
      jb=jb+nx
    enddo
 
!   remaining wavenumbers
      iabase=2*inc+1
      ibbase=(n-2)*inc+1
      jabase=3
      jbbase=n-1

    do k=3,nh,2
      ia=iabase
      ib=ibbase
      ja=jabase
      jb=jbbase
      c=trigs(n+k)
      s=trigs(n+k+1)
!dir$ ivdep
      do L=1,lot
        work(ja)=(a(ia)+a(ib))- &
            (s*(a(ia)-a(ib))+c*(a(ia+inc)+a(ib+inc)))
        work(jb)=(a(ia)+a(ib))+ &
            (s*(a(ia)-a(ib))+c*(a(ia+inc)+a(ib+inc)))
        work(ja+1)=(c*(a(ia)-a(ib))-s*(a(ia+inc)+a(ib+inc)))+ &
            (a(ia+inc)-a(ib+inc))
        work(jb+1)=(c*(a(ia)-a(ib))-s*(a(ia+inc)+a(ib+inc)))- &
            (a(ia+inc)-a(ib+inc))
        ia=ia+jump
        ib=ib+jump
        ja=ja+nx
        jb=jb+nx
      enddo
      iabase=iabase+ink
      ibbase=ibbase-ink
      jabase=jabase+2
      jbbase=jbbase-2
    enddo

!   wavenumber n/4 (if it exists)
    if (iabase.eq.ibbase) then
      ia=iabase
      ja=jabase
!dir$ ivdep
      do L=1,lot
        work(ja)=2.0*a(ia)
        work(ja+1)=-2.0*a(ia+inc)
        ia=ia+jump
        ja=ja+nx
      enddo
    endif

    end subroutine fft99a

!##########################################################################

    subroutine fft99b (work,a,trigs,inc,jump,n,lot)
    integer, intent(in)    :: inc,jump,n,lot
    real(r8),    intent(in)    :: trigs(:)
    real(r8),    intent(inout) :: a(*),work(*)

!     dimension work(n),a(n),trigs(n)
!
!     subroutine fft99b - postprocessing step for fft99, isign=-1
!     (gridpoint to spectral transform)

    integer :: nh, nx, ink, k, L
    integer :: ia, ib, ja, jb, iabase, ibbase, jabase, jbbase
    real(r8)    :: rscale, c, s

      nh=n/2
      nx=n+1
      ink=inc+inc

!   a(0) and a(n/2)
      rscale=1.0/(1.0 * n)
      ia=1
      ib=2
      ja=1
      jb=n*inc+1
!dir$ ivdep
    do L=1,lot
      a(ja)=rscale*(work(ia)+work(ib))
      a(jb)=rscale*(work(ia)-work(ib))
      a(ja+inc)=0.0
      a(jb+inc)=0.0
      ia=ia+nx
      ib=ib+nx
      ja=ja+jump
      jb=jb+jump
    enddo

!   remaining wavenumbers
      rscale=0.5*rscale
      iabase=3
      ibbase=n-1
      jabase=2*inc+1
      jbbase=(n-2)*inc+1

    do k=3,nh,2
      ia=iabase
      ib=ibbase
      ja=jabase
      jb=jbbase
      c=trigs(n+k)
      s=trigs(n+k+1)
!dir$ ivdep
      do L=1,lot
        a(ja)=rscale*((work(ia)+work(ib)) &
           +(c*(work(ia+1)+work(ib+1))+s*(work(ia)-work(ib))))
        a(jb)=rscale*((work(ia)+work(ib)) &
           -(c*(work(ia+1)+work(ib+1))+s*(work(ia)-work(ib))))
        a(ja+inc)=rscale*((c*(work(ia)-work(ib))-s*(work(ia+1)+work(ib+1))) &
            +(work(ib+1)-work(ia+1)))
        a(jb+inc)=rscale*((c*(work(ia)-work(ib))-s*(work(ia+1)+work(ib+1))) &
            -(work(ib+1)-work(ia+1)))
        ia=ia+nx
        ib=ib+nx
        ja=ja+jump
        jb=jb+jump
      enddo
      iabase=iabase+2
      ibbase=ibbase-2
      jabase=jabase+ink
      jbbase=jbbase-ink
    enddo

!   wavenumber n/4 (if it exists)
    if (iabase.eq.ibbase) then
      ia=iabase
      ja=jabase
      rscale=2.0*rscale
!dir$ ivdep
      do L=1,lot
        a(ja)=rscale*work(ia)
        a(ja+inc)=-rscale*work(ia+1)
        ia=ia+nx
        ja=ja+jump
      enddo
    endif

    end subroutine fft99b

!##########################################################################

    subroutine fft991(a,work,trigs,ifax,inc,jump,n,lot,isign)
    integer, intent(in)    :: inc,jump,n,lot,isign
    integer, intent(inout) :: ifax(:)
    real(r8),    intent(in)    :: trigs(:)
    real(r8),    intent(inout) :: a(*),work(*)

!     dimension a(n),work(n),trigs(n),ifax(1)
!
!     subroutine "fft991" - multiple real/half-complex periodic
!     fast fourier transform
!
!     same as fft99 except that ordering of data corresponds to
!     that in mrfft2
!
!     procedure used to convert to half-length complex transform
!     is given by cooley, lewis and welch (j. sound vib., vol. 12
!     (1970), 315-337)
!
!     a is the array containing input and output data
!     work is an area of size (n+1)*lot
!     trigs is a previously prepared list of trig function values
!     ifax is a previously prepared list of factors of n/2
!     inc is the increment within each data 'vector'
!         (e.g. inc=1 for consecutively stored data)
!     jump is the increment between the start of each data vector
!     n is the length of the data vectors
!     lot is the number of data vectors
!     isign = +1 for transform from spectral to gridpoint
!           = -1 for transform from gridpoint to spectral
!
!     ordering of coefficients:
!         a(0),b(0),a(1),b(1),a(2),b(2),...,a(n/2),b(n/2)
!         where b(0)=b(n/2)=0; (n+2) locations required
!
!     ordering of data:
!         x(0),x(1),x(2),...,x(n-1)
!
!     vectorization is achieved on cray by doing the transforms in
!     parallel
!
!     *** n.b. n is assumed to be an even number
!
!     definition of transforms:
!     -------------------------
!
!     isign=+1: x(j)=sum(k=0,...,n-1)(c(k)*exp(2*i*j*k*pi/n))
!         where c(k)=a(k)+i*b(k) and c(n-k)=a(k)-i*b(k)
!
!     isign=-1: a(k)=(1/n)*sum(j=0,...,n-1)(x(j)*cos(2*j*k*pi/n))
!               b(k)=-(1/n)*sum(j=0,...,n-1)(x(j)*sin(2*j*k*pi/n))
!

    integer :: nfax, nx, nh, ink, igo, ibase, jbase
    integer :: i, j, k, L, m, ia, la, ib

      nfax=ifax(1)
      nx=n+1
      nh=n/2
      ink=inc+inc
      if (isign.eq.+1) go to 30

!     if necessary, transfer data to work area
      igo=50
      if (mod(nfax,2).eq.1) goto 40
      ibase=1
      jbase=1
    do L=1,lot
      i=ibase
      j=jbase
!dir$ ivdep
      do m=1,n
        work(j)=a(i)
        i=i+inc
        j=j+1
      enddo
      ibase=ibase+jump
      jbase=jbase+nx
    enddo
!
      igo=60
      go to 40
!
!   preprocessing (isign=+1)
!   ------------------------
!
   30 continue
      call fft99a(a,work,trigs,inc,jump,n,lot)
      igo=60
!
!   complex transform
!   -----------------
!
   40 continue
      ia=1
      la=1
    do k=1,nfax
      if (igo.eq.60) go to 60
   50 continue
        call vpassm (a(ia),a(ia+inc),work(1),work(2),trigs, &
                     ink,2,jump,nx,lot,nh,ifax(k+1),la)
      igo=60
      go to 70
   60 continue
        call vpassm (work(1),work(2),a(ia),a(ia+inc),trigs, &
                     2,ink,nx,jump,lot,nh,ifax(k+1),la)
      igo=50
   70 continue
      la=la*ifax(k+1)
    enddo

    if (isign.eq.-1) go to 130

! if necessary, transfer data from work area
    if (mod(nfax,2).ne.1) then
      ibase=1
      jbase=1
      do L=1,lot
        i=ibase
        j=jbase
!dir$ ivdep
        do m=1,n
          a(j)=work(i)
          i=i+1
          j=j+inc
        enddo
        ibase=ibase+nx
        jbase=jbase+jump
      enddo
    endif

!   fill in zeros at end
      ib=n*inc+1
!dir$ ivdep
      do L=1,lot
        a(ib)=0.0
        a(ib+inc)=0.0
        ib=ib+jump
      enddo
      go to 140

!     postprocessing (isign=-1):
!     --------------------------

  130 continue
      call fft99b (work,a,trigs,inc,jump,n,lot)

  140 continue

    end subroutine fft991

!##########################################################################

    subroutine set99 (trigs, ifax, n)
    integer, intent(in)  :: n
    integer, intent(out) :: ifax(:)
    real(r8),    intent(out) :: trigs(:)

!     dimension ifax(13),trigs(1)
!
! mode 3 is used for real/half-complex transforms.  it is possible
! to do complex/complex transforms with other values of mode, but
! documentation of the details were not available when this routine
! was written.
!
    integer :: mode = 3
    integer :: i

      call fax (ifax, n, mode)
      i = ifax(1)
      if (ifax(i+1) .gt. 5 .or. n .le. 4) ifax(1) = -99
      if (ifax(1) .le. 0 ) then 
        write(6,*) ' set99 -- invalid n'
        stop
      endif
      call fftrig (trigs, n, mode)

    end subroutine set99

!##########################################################################

    subroutine fax (ifax,n,mode)
    integer, intent(out) :: ifax(:)
    integer, intent(in)  :: n, mode

    integer :: nn, k, L, inc, nfax, ii, istop, i, item

      nn=n
      if (iabs(mode).eq.1) go to 10
      if (iabs(mode).eq.8) go to 10
      nn=n/2
      if ((nn+nn).eq.n) go to 10
      ifax(1)=-99
      return
   10 k=1
!     test for factors of 4
   20 if (mod(nn,4).ne.0) go to 30
      k=k+1
      ifax(k)=4
      nn=nn/4
      if (nn.eq.1) go to 80
      go to 20
!     test for extra factor of 2
   30 if (mod(nn,2).ne.0) go to 40
      k=k+1
      ifax(k)=2
      nn=nn/2
      if (nn.eq.1) go to 80
!     test for factors of 3
   40 if (mod(nn,3).ne.0) go to 50
      k=k+1
      ifax(k)=3
      nn=nn/3
      if (nn.eq.1) go to 80
      go to 40
!     now find remaining factors
   50 L=5
      inc=2
!     inc alternately takes on values 2 and 4
   60 if (mod(nn,L).ne.0) go to 70
      k=k+1
      ifax(k)=L
      nn=nn/L
      if (nn.eq.1) go to 80
      go to 60
   70 L=L+inc
      inc=6-inc
      go to 60
   80 ifax(1)=k-1
!     ifax(1) contains number of factors
      nfax=ifax(1)
!     sort factors into ascending order
      if (nfax.eq.1) go to 110
      do 100 ii=2,nfax
      istop=nfax+2-ii
      do 90 i=2,istop
      if (ifax(i+1).ge.ifax(i)) go to 90
      item=ifax(i)
      ifax(i)=ifax(i+1)
      ifax(i+1)=item
   90 continue
  100 continue
  110 continue

    end subroutine fax

!##########################################################################

    subroutine fftrig (trigs,n,mode)
    real(r8),    intent(out) :: trigs(:)
    integer, intent(in)  :: n, mode

    real(r8)    :: del, angle
    integer :: imode, nn, nh, i, L, la

      imode=iabs(mode)
      nn=n
      if (imode.gt.1.and.imode.lt.6) nn=n/2
      del=(pi+pi)/(1.0 * nn)
      L=nn+nn
      do i=1,L,2
        angle=0.5*(i-1.0)*del
        trigs(i)=cos(angle)
        trigs(i+1)=sin(angle)
      enddo
      if (imode.eq.1) return
      if (imode.eq.8) return

      del=0.5*del
      nh=(nn+1)/2
      L=nh+nh
      la=nn+nn
      do i=1,L,2
        angle=0.5*(i-1.0)*del
        trigs(la+i)=cos(angle)
        trigs(la+i+1)=sin(angle)
      enddo
      if (imode.le.3) return

      del=0.5*del
      la=la+nn
    if (mode.ne.5) then
      do i=2,nn
        angle=(i-1.0)*del
        trigs(la+i)=2.0*sin(angle)
      enddo
      return
    endif

      del=0.5*del
      do i=2,n
        angle=(i-1.0)*del
        trigs(la+i)=sin(angle)
      enddo

    end subroutine fftrig

!##########################################################################

    subroutine vpassm (a,b,c,d,trigs,inc1,inc2,inc3,inc4,lot,n,ifac,la)
    integer, intent(in)  :: inc1, inc2, inc3, inc4, lot, n, ifac, la
    real(r8),    intent(in)  :: a(*),b(*),trigs(*)
    real(r8),    intent(out) :: c(*),d(*)
!
!     subroutine "vpassm" - multiple version of "vpassa"
!     performs one pass through data
!     as part of multiple complex fft routine
!     a is first real input vector
!     b is first imaginary input vector
!     c is first real output vector
!     d is first imaginary output vector
!     trigs is precalculated table of sines " cosines
!     inc1 is addressing increment for a and b
!     inc2 is addressing increment for c and d
!     inc3 is addressing increment between a"s & b"s
!     inc4 is addressing increment between c"s & d"s
!     lot is the number of vectors
!     n is length of vectors
!     ifac is current factor of n
!     la is product of previous factors
!

    real(r8) :: sin36=0.587785252292473
    real(r8) :: cos36=0.809016994374947
    real(r8) :: sin72=0.951056516295154
    real(r8) :: cos72=0.309016994374947
    real(r8) :: sin60=0.866025403784437

    integer :: i, j, k, L, m, iink, jink, jump, ibase, jbase, igo, ijk, la1
    integer :: ia, ja, ib, jb, kb, ic, jc, kc, id, jd, kd, ie, je, ke
    real(r8)    :: c1, s1, c2, s2, c3, s3, c4, s4
      m=n/ifac
      iink=m*inc1
      jink=la*inc2
      jump=(ifac-1)*jink
      ibase=0
      jbase=0
      igo=ifac-1
      if (igo.gt.4) return
!del  go to (10,50,90,130),igo

  select case (igo)

!   coding for factor 2

    case (1)
   10 ia=1
      ja=1
      ib=ia+iink
      jb=ja+jink
      do 20 L=1,la
      i=ibase
      j=jbase
!dir$ ivdep
      do 15 ijk=1,lot
!!!if(ib + 1 > size(a)) goto 15
      c(ja+j)=a(ia+i)+a(ib+i)
      d(ja+j)=b(ia+i)+b(ib+i)
      c(jb+j)=a(ia+i)-a(ib+i)
      d(jb+j)=b(ia+i)-b(ib+i)
      i=i+inc3
      j=j+inc4
   15 continue
      ibase=ibase+inc1
      jbase=jbase+inc2
   20 continue
      if (la.eq.m) return
      la1=la+1
      jbase=jbase+jump
      do 40 k=la1,m,la
      kb=k+k-2
      c1=trigs(kb+1)
      s1=trigs(kb+2)
      do 30 L=1,la
      i=ibase
      j=jbase
!dir$ ivdep
      do 25 ijk=1,lot
      c(ja+j)=a(ia+i)+a(ib+i)
      d(ja+j)=b(ia+i)+b(ib+i)
      c(jb+j)=c1*(a(ia+i)-a(ib+i))-s1*(b(ia+i)-b(ib+i))
      d(jb+j)=s1*(a(ia+i)-a(ib+i))+c1*(b(ia+i)-b(ib+i))
      i=i+inc3
      j=j+inc4
   25 continue
      ibase=ibase+inc1
      jbase=jbase+inc2
   30 continue
      jbase=jbase+jump
   40 continue
!     return

!   coding for factor 3

    case (2)
   50 ia=1
      ja=1
      ib=ia+iink
      jb=ja+jink
      ic=ib+iink
      jc=jb+jink
      do 60 L=1,la
      i=ibase
      j=jbase
!dir$ ivdep
      do 55 ijk=1,lot
      c(ja+j)=a(ia+i)+(a(ib+i)+a(ic+i))
      d(ja+j)=b(ia+i)+(b(ib+i)+b(ic+i))
      c(jb+j)=(a(ia+i)-0.5*(a(ib+i)+a(ic+i)))-(sin60*(b(ib+i)-b(ic+i)))
      c(jc+j)=(a(ia+i)-0.5*(a(ib+i)+a(ic+i)))+(sin60*(b(ib+i)-b(ic+i)))
      d(jb+j)=(b(ia+i)-0.5*(b(ib+i)+b(ic+i)))+(sin60*(a(ib+i)-a(ic+i)))
      d(jc+j)=(b(ia+i)-0.5*(b(ib+i)+b(ic+i)))-(sin60*(a(ib+i)-a(ic+i)))
      i=i+inc3
      j=j+inc4
   55 continue
      ibase=ibase+inc1
      jbase=jbase+inc2
   60 continue
      if (la.eq.m) return
      la1=la+1
      jbase=jbase+jump
      do 80 k=la1,m,la
      kb=k+k-2
      kc=kb+kb
      c1=trigs(kb+1)
      s1=trigs(kb+2)
      c2=trigs(kc+1)
      s2=trigs(kc+2)
      do 70 L=1,la
      i=ibase
      j=jbase
!dir$ ivdep
      do 65 ijk=1,lot
      c(ja+j)=a(ia+i)+(a(ib+i)+a(ic+i))
      d(ja+j)=b(ia+i)+(b(ib+i)+b(ic+i))
      c(jb+j)=                                                           &
          c1*((a(ia+i)-0.5*(a(ib+i)+a(ic+i)))-(sin60*(b(ib+i)-b(ic+i)))) &
         -s1*((b(ia+i)-0.5*(b(ib+i)+b(ic+i)))+(sin60*(a(ib+i)-a(ic+i))))
      d(jb+j)=                                                           &
          s1*((a(ia+i)-0.5*(a(ib+i)+a(ic+i)))-(sin60*(b(ib+i)-b(ic+i)))) &
         +c1*((b(ia+i)-0.5*(b(ib+i)+b(ic+i)))+(sin60*(a(ib+i)-a(ic+i))))
      c(jc+j)=                                                           &
          c2*((a(ia+i)-0.5*(a(ib+i)+a(ic+i)))+(sin60*(b(ib+i)-b(ic+i)))) &
         -s2*((b(ia+i)-0.5*(b(ib+i)+b(ic+i)))-(sin60*(a(ib+i)-a(ic+i))))
      d(jc+j)=                                                           &
          s2*((a(ia+i)-0.5*(a(ib+i)+a(ic+i)))+(sin60*(b(ib+i)-b(ic+i)))) &
         +c2*((b(ia+i)-0.5*(b(ib+i)+b(ic+i)))-(sin60*(a(ib+i)-a(ic+i))))
      i=i+inc3
      j=j+inc4
   65 continue
      ibase=ibase+inc1
      jbase=jbase+inc2
   70 continue
      jbase=jbase+jump
   80 continue
!     return

!   coding for factor 4

    case (3)
   90 ia=1
      ja=1
      ib=ia+iink
      jb=ja+jink
      ic=ib+iink
      jc=jb+jink
      id=ic+iink
      jd=jc+jink
      do 100 L=1,la
      i=ibase
      j=jbase
!dir$ ivdep
      do 95 ijk=1,lot
      c(ja+j)=(a(ia+i)+a(ic+i))+(a(ib+i)+a(id+i))
      c(jc+j)=(a(ia+i)+a(ic+i))-(a(ib+i)+a(id+i))
      d(ja+j)=(b(ia+i)+b(ic+i))+(b(ib+i)+b(id+i))
      d(jc+j)=(b(ia+i)+b(ic+i))-(b(ib+i)+b(id+i))
      c(jb+j)=(a(ia+i)-a(ic+i))-(b(ib+i)-b(id+i))
      c(jd+j)=(a(ia+i)-a(ic+i))+(b(ib+i)-b(id+i))
      d(jb+j)=(b(ia+i)-b(ic+i))+(a(ib+i)-a(id+i))
      d(jd+j)=(b(ia+i)-b(ic+i))-(a(ib+i)-a(id+i))
      i=i+inc3
      j=j+inc4
   95 continue
      ibase=ibase+inc1
      jbase=jbase+inc2
  100 continue
      if (la.eq.m) return
      la1=la+1
      jbase=jbase+jump
      do 120 k=la1,m,la
      kb=k+k-2
      kc=kb+kb
      kd=kc+kb
      c1=trigs(kb+1)
      s1=trigs(kb+2)
      c2=trigs(kc+1)
      s2=trigs(kc+2)
      c3=trigs(kd+1)
      s3=trigs(kd+2)
      do 110 L=1,la
      i=ibase
      j=jbase
!dir$ ivdep
      do 105 ijk=1,lot
      c(ja+j)=(a(ia+i)+a(ic+i))+(a(ib+i)+a(id+i))
      d(ja+j)=(b(ia+i)+b(ic+i))+(b(ib+i)+b(id+i))
      c(jc+j)=                                     &
          c2*((a(ia+i)+a(ic+i))-(a(ib+i)+a(id+i))) &
         -s2*((b(ia+i)+b(ic+i))-(b(ib+i)+b(id+i)))
      d(jc+j)=                                     &
          s2*((a(ia+i)+a(ic+i))-(a(ib+i)+a(id+i))) &
         +c2*((b(ia+i)+b(ic+i))-(b(ib+i)+b(id+i)))
      c(jb+j)=                                     &
          c1*((a(ia+i)-a(ic+i))-(b(ib+i)-b(id+i))) &
         -s1*((b(ia+i)-b(ic+i))+(a(ib+i)-a(id+i)))
      d(jb+j)=                                     &
          s1*((a(ia+i)-a(ic+i))-(b(ib+i)-b(id+i))) &
         +c1*((b(ia+i)-b(ic+i))+(a(ib+i)-a(id+i)))
      c(jd+j)=                                     &
          c3*((a(ia+i)-a(ic+i))+(b(ib+i)-b(id+i))) &
         -s3*((b(ia+i)-b(ic+i))-(a(ib+i)-a(id+i)))
      d(jd+j)=                                     &
          s3*((a(ia+i)-a(ic+i))+(b(ib+i)-b(id+i))) &
         +c3*((b(ia+i)-b(ic+i))-(a(ib+i)-a(id+i)))
      i=i+inc3
      j=j+inc4
  105 continue
      ibase=ibase+inc1
      jbase=jbase+inc2
  110 continue
      jbase=jbase+jump
  120 continue
!     return

!   coding for factor 5

    case (4)
  130 ia=1
      ja=1
      ib=ia+iink
      jb=ja+jink
      ic=ib+iink
      jc=jb+jink
      id=ic+iink
      jd=jc+jink
      ie=id+iink
      je=jd+jink
      do 140 L=1,la
      i=ibase
      j=jbase
!dir$ ivdep
      do 135 ijk=1,lot
      c(ja+j)=a(ia+i)+(a(ib+i)+a(ie+i))+(a(ic+i)+a(id+i))
      d(ja+j)=b(ia+i)+(b(ib+i)+b(ie+i))+(b(ic+i)+b(id+i))
      c(jb+j)=(a(ia+i)+cos72*(a(ib+i)+a(ie+i))-cos36*(a(ic+i)+a(id+i))) &
        -(sin72*(b(ib+i)-b(ie+i))+sin36*(b(ic+i)-b(id+i)))
      c(je+j)=(a(ia+i)+cos72*(a(ib+i)+a(ie+i))-cos36*(a(ic+i)+a(id+i))) &
        +(sin72*(b(ib+i)-b(ie+i))+sin36*(b(ic+i)-b(id+i)))
      d(jb+j)=(b(ia+i)+cos72*(b(ib+i)+b(ie+i))-cos36*(b(ic+i)+b(id+i))) &
        +(sin72*(a(ib+i)-a(ie+i))+sin36*(a(ic+i)-a(id+i)))
      d(je+j)=(b(ia+i)+cos72*(b(ib+i)+b(ie+i))-cos36*(b(ic+i)+b(id+i))) &
        -(sin72*(a(ib+i)-a(ie+i))+sin36*(a(ic+i)-a(id+i)))
      c(jc+j)=(a(ia+i)-cos36*(a(ib+i)+a(ie+i))+cos72*(a(ic+i)+a(id+i))) &
        -(sin36*(b(ib+i)-b(ie+i))-sin72*(b(ic+i)-b(id+i)))
      c(jd+j)=(a(ia+i)-cos36*(a(ib+i)+a(ie+i))+cos72*(a(ic+i)+a(id+i))) &
        +(sin36*(b(ib+i)-b(ie+i))-sin72*(b(ic+i)-b(id+i)))
      d(jc+j)=(b(ia+i)-cos36*(b(ib+i)+b(ie+i))+cos72*(b(ic+i)+b(id+i))) &
        +(sin36*(a(ib+i)-a(ie+i))-sin72*(a(ic+i)-a(id+i)))
      d(jd+j)=(b(ia+i)-cos36*(b(ib+i)+b(ie+i))+cos72*(b(ic+i)+b(id+i))) &
        -(sin36*(a(ib+i)-a(ie+i))-sin72*(a(ic+i)-a(id+i)))
      i=i+inc3
      j=j+inc4
  135 continue
      ibase=ibase+inc1
      jbase=jbase+inc2
  140 continue
      if (la.eq.m) return
      la1=la+1
      jbase=jbase+jump
      do 160 k=la1,m,la
      kb=k+k-2
      kc=kb+kb
      kd=kc+kb
      ke=kd+kb
      c1=trigs(kb+1)
      s1=trigs(kb+2)
      c2=trigs(kc+1)
      s2=trigs(kc+2)
      c3=trigs(kd+1)
      s3=trigs(kd+2)
      c4=trigs(ke+1)
      s4=trigs(ke+2)
      do 150 L=1,la
      i=ibase
      j=jbase
!dir$ ivdep
      do 145 ijk=1,lot
      c(ja+j)=a(ia+i)+(a(ib+i)+a(ie+i))+(a(ic+i)+a(id+i))
      d(ja+j)=b(ia+i)+(b(ib+i)+b(ie+i))+(b(ic+i)+b(id+i))
      c(jb+j)=                                                          &
          c1*((a(ia+i)+cos72*(a(ib+i)+a(ie+i))-cos36*(a(ic+i)+a(id+i))) &
            -(sin72*(b(ib+i)-b(ie+i))+sin36*(b(ic+i)-b(id+i))))         &
         -s1*((b(ia+i)+cos72*(b(ib+i)+b(ie+i))-cos36*(b(ic+i)+b(id+i))) &
            +(sin72*(a(ib+i)-a(ie+i))+sin36*(a(ic+i)-a(id+i))))
      d(jb+j)=                                                          &
          s1*((a(ia+i)+cos72*(a(ib+i)+a(ie+i))-cos36*(a(ic+i)+a(id+i))) &
            -(sin72*(b(ib+i)-b(ie+i))+sin36*(b(ic+i)-b(id+i))))         &
         +c1*((b(ia+i)+cos72*(b(ib+i)+b(ie+i))-cos36*(b(ic+i)+b(id+i))) &
            +(sin72*(a(ib+i)-a(ie+i))+sin36*(a(ic+i)-a(id+i))))
      c(je+j)=                                                          &
          c4*((a(ia+i)+cos72*(a(ib+i)+a(ie+i))-cos36*(a(ic+i)+a(id+i))) &
            +(sin72*(b(ib+i)-b(ie+i))+sin36*(b(ic+i)-b(id+i))))         &
         -s4*((b(ia+i)+cos72*(b(ib+i)+b(ie+i))-cos36*(b(ic+i)+b(id+i))) &
            -(sin72*(a(ib+i)-a(ie+i))+sin36*(a(ic+i)-a(id+i))))
      d(je+j)=                                                          &
          s4*((a(ia+i)+cos72*(a(ib+i)+a(ie+i))-cos36*(a(ic+i)+a(id+i))) &
            +(sin72*(b(ib+i)-b(ie+i))+sin36*(b(ic+i)-b(id+i))))         &
         +c4*((b(ia+i)+cos72*(b(ib+i)+b(ie+i))-cos36*(b(ic+i)+b(id+i))) &
            -(sin72*(a(ib+i)-a(ie+i))+sin36*(a(ic+i)-a(id+i))))
      c(jc+j)=                                                          &
          c2*((a(ia+i)-cos36*(a(ib+i)+a(ie+i))+cos72*(a(ic+i)+a(id+i))) &
            -(sin36*(b(ib+i)-b(ie+i))-sin72*(b(ic+i)-b(id+i))))         &
         -s2*((b(ia+i)-cos36*(b(ib+i)+b(ie+i))+cos72*(b(ic+i)+b(id+i))) &
            +(sin36*(a(ib+i)-a(ie+i))-sin72*(a(ic+i)-a(id+i))))
      d(jc+j)=                                                          &
          s2*((a(ia+i)-cos36*(a(ib+i)+a(ie+i))+cos72*(a(ic+i)+a(id+i))) &
            -(sin36*(b(ib+i)-b(ie+i))-sin72*(b(ic+i)-b(id+i))))         &
         +c2*((b(ia+i)-cos36*(b(ib+i)+b(ie+i))+cos72*(b(ic+i)+b(id+i))) &
            +(sin36*(a(ib+i)-a(ie+i))-sin72*(a(ic+i)-a(id+i))))
      c(jd+j)=                                                          &
          c3*((a(ia+i)-cos36*(a(ib+i)+a(ie+i))+cos72*(a(ic+i)+a(id+i))) &
            +(sin36*(b(ib+i)-b(ie+i))-sin72*(b(ic+i)-b(id+i))))         &
         -s3*((b(ia+i)-cos36*(b(ib+i)+b(ie+i))+cos72*(b(ic+i)+b(id+i))) &
            -(sin36*(a(ib+i)-a(ie+i))-sin72*(a(ic+i)-a(id+i))))
      d(jd+j)=                                                          &
          s3*((a(ia+i)-cos36*(a(ib+i)+a(ie+i))+cos72*(a(ic+i)+a(id+i))) &
            +(sin36*(b(ib+i)-b(ie+i))-sin72*(b(ic+i)-b(id+i))))         &
         +c3*((b(ia+i)-cos36*(b(ib+i)+b(ie+i))+cos72*(b(ic+i)+b(id+i))) &
            -(sin36*(a(ib+i)-a(ie+i))-sin72*(a(ic+i)-a(id+i))))
      i=i+inc3
      j=j+inc4
  145 continue
      ibase=ibase+inc1
      jbase=jbase+inc2
  150 continue
      jbase=jbase+jump
  160 continue

  end select

    end subroutine vpassm

!##########################################################################

end module fft99_mod

