program main

implicit none

integer :: i
integer, parameter :: interval = 7

write(*, *) 'set_def.out'
write(*, *) 1
write(*, *) 1800

do i = 1, 1800
   write(*, *) 100.0
   write(*, *) 1 + (i - 1) * 6
end do

end program main
