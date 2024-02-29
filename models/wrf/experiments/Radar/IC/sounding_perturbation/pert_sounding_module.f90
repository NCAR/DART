! DART software - Copyright UCAR. This open source software is provided
! by UCAR, "as is", without charge, subject to all terms of use at
! http://www.image.ucar.edu/DAReS/DART/DART_download
!
! $Id$

module pert_sounding_mod

! This program adds perturbations to a sounding
! based on an initial sounding that's input from
! the file "input_sounding" - this is considered
! the base sounding
!
! Written by Altug Aksoy for WRF/DART 11/09/2006
! Based on code by David Dowell
!>@todo FIXME remove the numerical recipes routines

implicit none

real, parameter  :: PI = 3.1415926535897932346
integer          :: iseed1, iseed2


contains


! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


subroutine rand_init()

! Initializes seeds for the Gaussian and uniform random number generators
! This is globally done (iseed1 and iseed2 defined in module) so that one
! single sequence is followed every time the program is run
! iseed1 and iseed2 should not be changed anywhere else in the code

implicit none

integer                :: values(8)
character(len=10)      :: date, time, zone

call date_and_time(date,time,zone,values)
iseed1 = -sum(values)
call date_and_time(date,time,zone,values)
iseed2 = -sum(values)

end subroutine rand_init


! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


subroutine read_numlevels(nz)

! Determines number of levels from the file "input_sounding"

implicit none

integer, intent(out)    :: nz
real                    :: ps, ts, qvs
character(len=10)       :: dummy
logical                 :: end_file
integer                 :: k

open(unit=10,file='input_sounding',form='formatted',status='old')
rewind(10)
read(10,*) ps, ts, qvs

end_file = .false.
k = 0
do while (.not. end_file)
   read(10,*,end=1) dummy
   k = k + 1
   goto 2
1  end_file = .true.
2  continue
enddo

close(unit=10)

nz = k

return

end subroutine read_numlevels


! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


subroutine read_sounding(nz,s1d,v1d)

! Reads in base sounding data from the file "input_sounding"
! Based on WRF subroutine of same name in WRF/dyn_em/module_initialize_quater_ss.F

implicit none

integer, intent(in)     :: nz
real, intent(inout)     :: s1d(4), v1d(5,nz)

integer                 :: i, iz

open(unit=10,file='input_sounding',form='formatted',status='old')
rewind(10)
read(10,*) (s1d(i),i=1,4)

do iz = 1, nz
   read(10,*) (v1d(i,iz),i=1,5)
enddo

close(unit=10,status = 'keep')

return

end subroutine read_sounding


! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


subroutine level1d(nz,h,hgt,lvl)

! This subroutine determines the level index of sounding height indices that is
! closest to the specified height

implicit none

integer,intent(in)             :: nz
real, dimension(:)             :: h
real, intent(in)               :: hgt
integer, intent(inout)         :: lvl

real                           :: hmin, hh(nz)
integer                        :: k

hh = h
hh(1:nz) = abs(hh(1:nz) - hgt)
hmin    = minval(hh(1:nz))

! Determine the lowest level at which height=minimum
! This is actually not necessary for a monotonic variable like height
do k = nz,1,-1
   if (hh(k).eq.hmin) lvl=k
enddo

return

end subroutine level1d


! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


subroutine pert1d(pkind,q,mag,h,nz,kmodes,nmodes)

! This routine computes perturbations for a sounding with random phase/amplitude

implicit none

integer,intent(in)     :: pkind,nz,nmodes
real, intent(in)       :: mag
integer, dimension(:)  :: kmodes
real, dimension(:)     :: h
real, dimension(:)     :: q

real, parameter        :: PI = 3.1415926535897932346
real                   :: amp(nmodes), phs(nmodes)
real                   :: gaussdev, ran2
integer                :: k, m


if (pkind.eq.1) then
   ! Determine random amplitude/phase for each mode
   do m=1, nmodes
      amp(m) = mag*gaussdev(iseed1)
      if (kmodes(m).eq.0) then
         phs(m) = 0
      else
         phs(m) = 2.0*PI*(2.0*ran2(iseed2)-1.0)
      endif
   enddo
   ! Add cosine perturbations at each level
   do k=1, nz
      do m=1, nmodes
         q(k) = q(k) + amp(m)*cos(phs(m)+kmodes(m)*2.0*PI*(h(k)-h(1))/(h(nz)-h(1)))
      enddo
   enddo

elseif (pkind.eq.2) then
   ! Add uncorrelated Gaussian perturbations at each level
   do k=1, nz
      q(k) = mag*gaussdev(iseed1)
   enddo

endif

return

end subroutine pert1d


! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


subroutine norm1d(q,stdbot,stdtop,nm,ns,h,nz,nens)

! This subroutine normalizes the spread of a given perturbed profile so that
! low-level spread is constant between lowest level and mid-level "nm" at
! magnitude "stdbot" (low-level spread zone), varies linearly between "nm"
! and "ns" from "stdbot" to "stdtop" (transition zone), and is again constant
! from "ns" to the highest level at magnitude "stdtop" (upper-level spread zone)

! Standard deviation computed according to Numerical Recipes 14.1

implicit none

integer,intent(in)     :: nz, nm, ns, nens
real, intent(in)       :: stdbot, stdtop
real, dimension(:)     :: h
real, dimension(:,:)   :: q

real                   :: std(nz)
real                   :: mm, erv, ep, ss, sprd
integer                :: i, k, ntop

! Compute standard deviation at each level
! See Num. Recipes 14.1 for this variant
! (corrects for bias due to sampling)
std = 0.0
do k= 1, nz
   mm  = sum(q(:,k))/nens
   ep  = 0.0
   std(k) = 0.0
   do i = 1, nens
      ss = q(i,k) - mm
      ep = ep + ss
      std(k) = std(k) + ss*ss
   enddo
   std(k) = sqrt((std(k)-ep**2/nens)/(nens-1))
enddo
! Normalize according to specified std deviation bounds
if ((stdbot.eq.stdtop).or.(nm.eq.ns).or.(nm.ge.nz)) then
   ! If bounds are equal, no transition zone specified, or
   ! transition zone above domain top, constant spread
   ! all across the domain
   do k = 1, nz
      do i = 1, nens
         q(i,k) = q(i,k)/std(k)*stdbot
      enddo
   enddo
else
   ! Constant spread from domain bottom to nm
   do k = 1, nm
      do i = 1, nens
         q(i,k) = q(i,k)/std(k)*stdbot
      enddo
   enddo
   ! Linearly varying spread in transition zone
   ntop = ns
   if (ns.ge.nz) ntop = nz
   do k = nm+1, ntop
      sprd = stdbot + (stdtop-stdbot)*(h(k)-h(nm+1))/(h(ntop)-h(nm+1))
      do i = 1, nens
         q(i,k) = q(i,k)/std(k)*sprd
      enddo
   enddo
   ! If specified top height lower than domain top, rest of domain set
   ! to top-level specified spread
   if (ntop.ne.nz) then
      do k = ntop+1,nz
         do i = 1, nens
            q(i,k) = q(i,k)/std(k)*stdtop
         enddo
      enddo
   endif
endif
return

end subroutine norm1d


! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


subroutine write_sounding_ens( nz,s1d,v1d,imem )

! Writes (perturbed) sounding data to file "input_soundingx" with member index attached

implicit none

integer, intent(in)     :: nz, imem
real, dimension(:)      :: s1d
real, dimension(:,:)    :: v1d

integer                 :: i, iz
character(len=3)        :: mem_num
character(len=20)       :: filename
real                    :: ps, ts, qs, hs
real                    :: hk, tk, qk, uk, vk

if (imem.le.9) then
   write(mem_num,'(i1.0,2x)') imem
elseif (imem.le.99) then
   write(mem_num,'(i2.0,1x)') imem
else
   write(mem_num,'(i3.0)') imem
endif

filename = 'input_sounding'//mem_num

open(unit=10,file=trim(adjustl(filename)),form='formatted',status='replace')
ps = s1d(1)
ts = s1d(2)
qs = s1d(3)
hs = s1d(4)
write(10,'(f9.2,1x,f7.2,1x,f7.3,1x,f9.2)') ps, ts, qs, hs

do iz = 1, nz
   hk = v1d(1,iz)
   tk = v1d(2,iz)
   qk = v1d(3,iz)
   uk = v1d(4,iz)
   vk = v1d(5,iz)
   write(10,'(f9.2,1x,f7.2,1x,f7.3,1x,f13.8,1x,f13.8)') hk, tk, qk, uk, vk
enddo

close(unit=10,status = 'keep')

return

end subroutine write_sounding_ens


end module pert_sounding_mod

! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


function gaussdev(idum)

! Returns a normally distributed deviate with 0 mean and unit variance,
! using ran2(idum) as the source of uniform deviates; see Num'l Recipes,
! ch 7.2.

! Note: Need to initialize with a negative idum and shouldn't change idum
! during successive deviate calls

integer ::  idum
real    ::  gaussdev

integer ::  iset
real    ::  fac,gset,rsq,v1,v2,ran2

save    ::  iset,gset

data iset/0/

if (idum.lt.0) iset=0
if (iset.eq.0) then
1  v1 = 2.*ran2(idum)-1.
   v2 = 2.*ran2(idum)-1.
   rsq = v1**2 + v2**2
   if (rsq.ge.1. .or. rsq.eq.0.) goto 1
   fac = sqrt( -2.*log(rsq)/rsq )
   gset = v1*fac
   gaussdev = v2*fac
   iset= 1
else
   gaussdev = gset
   iset=0
endif

return

end function gaussdev


!------------------------------------------------------------------------


function ran2(idum)

! This function produces a random number in (0,1) from a uniform dist.
! Taken from Numerical Recipes 7.1

! Note: Need to initialize with a negative idum and shouldn't change idum
! during successive deviate calls

integer ::  idum,IM1,IM2,IMM1,IA1,IA2,IQ1,IQ2,IR1,IR2,NTAB,NDIV
real    ::  ran2,AM,EPS,RNMX

parameter (IM1=2147483563,IM2=2147483399,AM=1./IM1,IMM1=IM1-1, &
           IA1=40014,IA2=40692,IQ1=53668,IQ2=52774,IR1=12211,  &
           IR2=3791,NTAB=32,NDIV=1+IMM1/NTAB,EPS=1.2e-7,RNMX=1.-EPS)

integer ::  idum2,j,k,iv(NTAB),iy
save iv,iy,idum2
data idum2/123456789/, iv/NTAB*0/, iy/0/

if (idum.le.0) then
   idum=max(-idum,1)
   idum2=idum
   do j=NTAB+8,1,-1
      k=idum/IQ1
      idum=IA1*(idum-k*IQ1)-k*IR1
      if (idum.lt.0) idum=idum+IM1
      if (j.le.NTAB) iv(j)=idum
   enddo
   iy=iv(1)
endif

k=idum/IQ1
idum=IA1*(idum-k*IQ1)-k*IR1
if (idum.lt.0) idum=idum+IM1
k=idum2/IQ2
idum2=IA2*(idum2-k*IQ2)-k*IR2
if (idum2.lt.0) idum2=idum2+IM2
j=1+iy/NDIV
iy=iv(j)-idum2
iv(j)=idum
if (iy.lt.1) iy=iy+IMM1
ran2=min(AM*iy,RNMX)

return

end function ran2

! <next few lines under version control, do not edit>
! $URL$
! $Id$
! $Revision$
! $Date$
