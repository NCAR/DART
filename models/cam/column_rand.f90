program column_rand

! Allows creation of input file for generating a set of randomly located
! observation stations with full column of obs for CAM model. 

use random_seq_mod, only : random_seq_type, init_random_seq, random_uniform

implicit none

integer :: num_sets, level, num_cols, num_levs, i, j
integer, allocatable :: levels(:)
real :: lat, lon, t_err_var, uv_err_var, ps_err_var, q_err_var
real, parameter :: pi = 3.14159
type(random_seq_type) :: r

! Initialize the random sequence
call init_random_seq(r)

! Set up constants
write(*, *) 'Input number of sets'
read(*, *) num_sets

! Open an output file and write header info
open(unit = 20, file = 'cam_column_rand.out')
write(20, *) 'set_def.out'
write(20, *) num_sets

write(*, *) 'input the number of columns per set'
read(*, *) num_cols

write(*, *) 'input the number of model levels in column'
read(*, *) num_levs

allocate(levels(num_levs))
do i = 1, num_levs
   write(*, *) 'Input vertical level ', i
   read(*, *) levels(i)
end do

! First get error variance for surface pressure
write(*, *) 'Input error VARIANCE for surface pressure obs'
read(*, *) ps_err_var

! Get error variance for t, and u and v
write(*, *) 'Input error VARIANCE for T obs'
read(*, *) t_err_var
write(*, *) 'Input error VARIANCE for U and V obs'
read(*, *) uv_err_var
write(*, *) 'Input error VARIANCE for Q obs'
read(*, *) q_err_var

! Loop through each set
do j = 1, num_sets

! Output the total number of obs in set; Q is being observed, too
   write(*, *) 'total num is ', num_cols * (num_levs * 4 + 1)
   write(20, *) num_cols * (num_levs * 4 + 1)

! Loop through each column
   do i = 1, num_cols
! Get a random lon lat location for this column
   ! Longitude is random from 0 to 360
      lon = random_uniform(r) * 360.0

   ! Latitude must be area weighted
      lat = asin(random_uniform(r) * 2.0 - 1.0)
   ! Now convert from radians to degrees latitude
      lat = lat * 360.0 / (2.0 * pi)

! Do ps ob
      write(20, *) ps_err_var
      write(20, *) -1
! Level is -1 for ps
      write(20, *) -1
      write(20, *) lon
      write(20, *) lat
! Kind for surface pressure is 3
      write(20, *) 3

! Loop through each observation in the column
      do level = 1, num_levs
! Write out the t observation
         write(20, *) t_err_var
         write(20, *) -1
         write(20, *) levels(level)
         write(20, *) lon
         write(20, *) lat
! Kind for t is 4
         write(20, *) 4
! Write out the u observation
         write(20, *) uv_err_var
         write(20, *) -1
         write(20, *) levels(level)
         write(20, *) lon
         write(20, *) lat
! Kind for u is 1
         write(20, *) 1
! Write out the v observation
         write(20, *) uv_err_var
         write(20, *) -1
         write(20, *) levels(level)
         write(20, *) lon
         write(20, *) lat
! Kind for v is 2
         write(20, *) 2
! Write out the q observation
         write(20, *) q_err_var
         write(20, *) -1
         write(20, *) levels(level)
         write(20, *) lon
         write(20, *) lat
! Kind for q is 5???
         write(20, *) 5
      end do
   end do
end do

end program column_rand
