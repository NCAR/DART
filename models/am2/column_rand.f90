! DART software - Copyright UCAR. This open source software is provided
! by UCAR, "as is", without charge, subject to all terms of use at
! http://www.image.ucar.edu/DAReS/DART/DART_download
!
! $Id$
 
program column_rand

! Allows creation of input file for generating a set of randomly located
! observation stations with full column of obs for CAM model. 

use      types_mod, only : r8, PI
use  utilities_mod, only : get_unit, initialize_utilities, finalize_utilities
use random_seq_mod, only : random_seq_type, init_random_seq, random_uniform

implicit none

! version controlled file description for error handling, do not edit
character(len=256), parameter :: source   = &
   "$URL$"
character(len=32 ), parameter :: revision = "$Revision$"
character(len=128), parameter :: revdate  = "$Date$"

integer, allocatable :: levels(:)
integer  :: level, num_cols, num_levs, i, iunit
real(r8) :: lat, lon, t_err_var, uv_err_var, ps_err_var, q_err_var
type(random_seq_type) :: r

! Initialize the utilities
call initialize_utilities('Column_rand')

! Initialize the random sequence
call init_random_seq(r)

! Open an output file and write header info
iunit = get_unit()
open(unit = iunit, file = 'cam_column_rand.out')

write(*, *) 'input the number of columns per set'
read(*, *) num_cols

write(*, *) 'input the number of model levels in column'
read(*, *) num_levs

allocate(levels(num_levs))
do i = 1, num_levs
   write(*, *) 'Input vertical level ', i
   read(*, *) levels(i)
end do

! Output the total number of obs in set; Q is being observed, too
write(*, *) 'total num is ', num_cols * (num_levs * 4 + 1)
write(iunit, *) num_cols * (num_levs * 4 + 1)

! No copies or qc
write(iunit, *) 0
write(iunit, *) 0

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


! Loop through each column
do i = 1, num_cols

   ! Get a random lon lat location for this column
   ! Longitude is random from 0 to 360
   lon = random_uniform(r) * 360.0_r8

   ! Latitude must be area weighted
   lat = asin(random_uniform(r) * 2.0_r8 - 1.0_r8)

   ! Now convert from radians to degrees latitude
   lat = lat * 360.0_r8 / (2.0_r8 * PI)

   ! Do ps ob
   write(iunit, *) 0
   ! Kind for surface pressure is 3
   write(iunit, *) 3
   write(iunit, *) 1
   ! Level is -1 for ps
   write(iunit, *) -1
   write(iunit, *) lon
   write(iunit, *) lat
   write(iunit, *) 0, 0
   write(iunit, *) ps_err_var

   ! Loop through each observation in the column
   do level = 1, num_levs

      ! Write out the t observation
      write(iunit, *) 0
      ! Kind for t is 4
      write(iunit, *) 4
      write(iunit, *) 1
      write(iunit, *) levels(level)
      write(iunit, *) lon
      write(iunit, *) lat
      write(iunit, *) 0, 0
      write(iunit, *) t_err_var


      ! Write out the u observation
      write(iunit, *) 0
      ! Kind for u is 1
      write(iunit, *) 1
      write(iunit, *) 1
      write(iunit, *) levels(level)
      write(iunit, *) lon
      write(iunit, *) lat
      write(iunit, *) 0, 0
      write(iunit, *) uv_err_var


      ! Write out the v observation
      write(iunit, *) 0
      ! Kind for v is 2
      write(iunit, *) 2
      write(iunit, *) 1
      write(iunit, *) levels(level)
      write(iunit, *) lon
      write(iunit, *) lat
      write(iunit, *) 0, 0
      write(iunit, *) uv_err_var


      ! Write out the q observation
      write(iunit, *) 0
      ! Kind for q is 5
      write(iunit, *) 5
      write(iunit, *) 1
      write(iunit, *) levels(level)
      write(iunit, *) lon
      write(iunit, *) lat
      write(iunit, *) 0, 0
      write(iunit, *) q_err_var

   end do
end do

write(iunit, *) 'set_def.out'

! Shut down the utilities cleanly
call finalize_utilities()

end program column_rand

! <next few lines under version control, do not edit>
! $URL$
! $Id$
! $Revision$
! $Date$
