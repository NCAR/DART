
program driver

use model_mod

double precision :: x(model_size)

integer :: i, j

call init_conditions(x)

do i = 1, 100
   call adv_1step(x)
   call output(real(x), real(i))
end do

end program driver
