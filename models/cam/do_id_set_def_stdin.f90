program main

implicit none

integer :: i
integer, parameter :: interval = 105

write(*, *) 'set_def.out'
write(*, *) 1
write(*, *) 13440 / interval

do i = 1, 13440, interval
   write(*, *) 100.0
   write(*, *) i
end do

end program main
