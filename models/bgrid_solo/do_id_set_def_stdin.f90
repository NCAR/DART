program main

implicit none

integer :: i
integer, parameter :: interval = 7

write(*, *) 'set_def.out'
write(*, *) 1
write(*, *) 28200 / interval

do i = 1, 28200, interval
   write(*, *) 1.0
   write(*, *) i
end do

end program main
