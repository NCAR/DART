program test_random_nr


use random_numerical_recipes_mod, only : random_seq_type, init_ran1, ran1, gasdev

implicit none

type (random_seq_type) :: r
integer :: i, n
double precision :: r1, dist, mean_dist

mean_dist = 0.0
call init_ran1(r, -5)

n = 10000000

do i = 1, n
   r1 = gasdev(r)
   dist = dabs(r1)
   mean_dist = mean_dist + dist
end do

write(*, *) 'sd is ', mean_dist / n

end program test_random_nr
