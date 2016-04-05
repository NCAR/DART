! DART software - Copyright 2004 - 2011 UCAR. This open source software is
! provided by UCAR, "as is", without charge, subject to all terms of use at
! http://www.image.ucar.edu/DAReS/DART/DART_download

program tester

! <next few lines under version control, do not edit>
! $URL$
! $Id$
! $Revision$
! $Date$

  implicit none
  real*8, dimension(:,:,:), allocatable  :: val, rval
  integer :: i, j, k

  allocate(val(2,3,4))

  do i = 1, 2
     do j = 1, 3
        do k = 1, 4
           val(i,j,k) = 2*i+3*j+4*k
        end do
     end do
  end do

  print *, val(1,:,:)
  print *, val(2,:,:)

  allocate(rval(4,2,3))

  rval = reshape(val, (/4, 2, 3/), order = (/3,1,2/))

  print *, rval(:,1,:)
  print *, rval(:,2,:)

end program tester
