program tester

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
