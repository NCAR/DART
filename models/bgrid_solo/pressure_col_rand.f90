program pressure_col_rand

!  $Source$
!  $Revision$
!  $Date$

! Allows creation of input file for generating a set of randomly located
! observation stations with  column of obs for b-grid model. Observations
! are placed at 800, 600, 400 and 200 mb in each column. PS is also 
! observed.

use random_seq_mod, only : random_seq_type, init_random_seq, random_uniform

implicit none

! let CVS fill strings ... DO NOT EDIT ...
character(len=128) :: &
   source   = "$Source$", &
   revision = "$Revision$", &
   revdate  = "$Date$"

integer :: num_sets, level, num_cols, i
real :: lat, lon, t_err_var, uv_err_var, ps_err_var
real, parameter :: pi = 3.14159
type(random_seq_type) :: r

integer, parameter :: num_levs = 4
real :: pressure(num_levs)

! Can't remember array constructor syntax
do i = 1, num_levs
   pressure(i) = 100000.0 - 20000.0 * i
end do

! Initialize the random sequence
call init_random_seq(r)

! Set up constants
num_sets = 1

! Open an output file and write header info
open(unit = 20, file = 'pressure_col_rand.out')
write(20, *) 'set_def.out'
write(20, *) num_sets

write(*, *) 'input the number of columns'
read(*, *) num_cols

! Output the total number of obs
write(*, *) 'total num is ', num_cols * (num_levs * 3 + 1)
write(20, *) num_cols * (num_levs * 3 + 1)

! First get error variance for surface pressure
write(*, *) 'Input error VARIANCE for surface pressure obs'
read(*, *) ps_err_var

! Get error variance for t, and u and v
write(*, *) 'Input error VARIANCE for T obs'
read(*, *) t_err_var
write(*, *) 'Input error VARIANCE for U and V obs'
read(*, *) uv_err_var


! Loop through each column
do i = 1, num_cols
! Get a random lon lat location for this column
   ! Longitude is random from 0 to 360
   lon = random_uniform(r) * 360.0
! Watch out for truncation range errors
   if(lon < 0.0) lon = 0.0
   if(lon > 360.0) lon = 360.0

   ! Latitude must be area weighted
   lat = asin(random_uniform(r) * 2.0 - 1.0)
   ! Now convert from radians to degrees latitude
   lat = lat * 360.0 / (2.0 * pi)
! Watch for truncation range errors
   if(lat < -90.0) lat = -90.0
   if(lat > 90.0) lat = 90.0

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
      write(20, *) -99
      write(20, *) pressure(level) 
      write(20, *) lon
      write(20, *) lat
! Kind for t is 4
      write(20, *) 4
! Write out the u observation
      write(20, *) uv_err_var
      write(20, *) -1
      write(20, *) -99
      write(20, *) pressure(level)
      write(20, *) lon
      write(20, *) lat
! Kind for u is 1
      write(20, *) 1
! Write out the t observation
      write(20, *) uv_err_var
      write(20, *) -1
      write(20, *) -99
      write(20, *) pressure(level)
      write(20, *) lon
      write(20, *) lat
! Kind for v is 2
      write(20, *) 2
   end do
end do

end program pressure_col_rand
