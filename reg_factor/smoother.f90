! Used to smooth regression confidence factor functions
! to get rid of sampling bumps

program smoother

implicit none

integer :: i, j , n, width
double precision, allocatable :: x(:), y(:), garb(:)
character (len = 30) :: in_file, out_file

write(*, *) 'input file name '
read(*, *) in_file

open(unit = 47, file = in_file)

write(*, *) 'output file name'
read(*, *) out_file

open(unit = 48, file = out_file)

write(*, *) 'input number of elements'
read(*, *) n

write(*, *) 'input smoothing half-width '
read(*, *) width

allocate(x(n), y(n), garb(n))
do i = 1, n
   read(47, *) garb(i), x(i)
end do


! Do the smoothing
y = x
do i = 1 + width, n - width
   y(i) = 0.0
   do j = i - width, i + width
      y(i) = y(i) + x(j)
   end do 
   y(i) = y(i) / (2.0 * width + 1.0)
end do

do i = 1, n
   write(48, *) garb(i), y(i)
end do

end program smoother

