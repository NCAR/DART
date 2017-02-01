! DART software - Copyright UCAR. This open source software is provided
! by UCAR, "as is", without charge, subject to all terms of use at
! http://www.image.ucar.edu/DAReS/DART/DART_download
!
! $Id$

program prog_dobsdx

!-----------------------------------------------------------
! Numerically differentiate the observation values for small
! change in input vector component. Used to interpolate 
! observation.

implicit none
real(8)        :: deltax = 0.1_8, x
real(8),allocatable :: obsx(:), obsxpd(:), dobsdx(:), locx(:)
character(128) :: file_in ='obsval.dat.bck'
character(128) :: file_out ='obsval.dat'
character(128) :: file_x ='model.in.bck'
integer :: nobs,i,j

nobs = 0
open(11,file=trim(file_in),status='old')
do
  read(11,'(1x)',end=123)
  nobs = nobs + 1
enddo
123 rewind(11)
allocate(obsx(nobs),obsxpd(nobs),dobsdx(nobs),locx(nobs))
open(12,file=trim(file_out),status='old')
do i = 1, nobs
   read(11,*)locx(i),j,obsx(i)
   read(12,*)locx(i),j,obsxpd(i)
   dobsdx(i) = (obsxpd(i) - obsx(i))/deltax
end do
close(11)
rewind(12)
open(11,file=trim(file_x),status='old')
read(11,*)
read(11,*)
read(11,*) x
close(11)
do i = 1, nobs
   write(12,*)locx(i),x,obsx(i),dobsdx(i)
enddo
close(12)
end program prog_dobsdx

! <next few lines under version control, do not edit>
! $URL$
! $Id$
! $Revision$
! $Date$
