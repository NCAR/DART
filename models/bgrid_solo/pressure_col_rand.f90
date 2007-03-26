! Data Assimilation Research Testbed -- DART
! Copyright 2004-2007, Data Assimilation Research Section
! University Corporation for Atmospheric Research
! Licensed under the GPL -- www.gpl.org/licenses/gpl.html
 
program pressure_col_rand

! <next few lines under version control, do not edit>
! $URL$
! $Id$
! $Revision$
! $Date$

! Allows creation of input file for generating a set of randomly located
! observation stations with  column of obs for b-grid model. Observations
! are placed at 800, 600, 400 and 200 mb in each column. PS is also 
! observed.

use      types_mod, only : r8, PI
use  utilities_mod, only : get_unit
use random_seq_mod, only : random_seq_type, init_random_seq, random_uniform

implicit none

! version controlled file description for error handling, do not edit
character(len=128), parameter :: &
   source   = "$URL$", &
   revision = "$Revision$", &
   revdate  = "$Date$"

integer  :: num_sets, level, num_cols, i, iunit
real(r8) :: lat, lon, t_err_var, uv_err_var, ps_err_var
type(random_seq_type) :: r

integer, parameter :: num_levs = 4
real(r8) :: pressure(num_levs)

! Can't remember array constructor syntax
do i = 1, num_levs
   pressure(i) = 1000.0_r8 - 200.0_r8 * i
end do

! Initialize the random sequence
call init_random_seq(r)

! Set up constants
num_sets = 1

! Open an output file
iunit = get_unit()
open(unit = iunit, file = 'pressure_col_rand.out')

write(*, *) 'input the number of columns'
read(*, *) num_cols

! Output the total number of obs
write(*, *) 'total num is ', num_cols * (num_levs * 3 + 1)
write(iunit, *) num_cols * (num_levs * 3 + 1)

! First get error variance for surface pressure
write(*, *) 'Input error VARIANCE for surface pressure obs'
read(*, *) ps_err_var

! Get error variance for t, and u and v
write(*, *) 'Input error VARIANCE for T obs'
read(*, *) t_err_var
write(*, *) 'Input error VARIANCE for U and V obs'
read(*, *) uv_err_var

! No values of qc
write(iunit, *) 0
write(iunit, *) 0

! Loop through each column
do i = 1, num_cols
   ! More obs coming 
   write(iunit, *) 0

   ! Get a random lon lat location for this column

   ! Longitude is random from 0 to 360
   lon = random_uniform(r) * 360.0_r8

   ! Watch out for truncation range errors
   if(lon <   0.0_r8) lon =   0.0_r8
   if(lon > 360.0_r8) lon = 360.0_r8

   ! Latitude must be area weighted
   lat = asin(random_uniform(r) * 2.0_r8 - 1.0_r8)

   ! Now convert from radians to degrees latitude
   lat = lat * 360.0_r8 / (2.0_r8 * PI)

   ! Watch for truncation range errors
   if(lat < -90.0_r8) lat = -90.0_r8
   if(lat >  90.0_r8) lat =  90.0_r8

   ! Do ps ob
   write(iunit, *) 'RADIOSONDE_SURFACE_PRESSURE'
   ! MOdel level of -1 for ps
   write(iunit, *) 1
   write(iunit, *) -1
   write(iunit, *) lon
   write(iunit, *) lat
   write(iunit, *) 0, 0
   write(iunit, *) ps_err_var

   ! Loop through each observation in the column
   do level = 1, num_levs

      write(iunit, *) 0
      ! Kind for t is 
      write(iunit, *) 'RADIOSONDE_TEMPERATURE'
      ! On pressure level
      write(iunit, *) 2
      write(iunit, *) pressure(level) 
      write(iunit, *) lon
      write(iunit, *) lat
      write(iunit, *) 0, 0
      write(iunit, *) t_err_var

      write(iunit, *) 0
      ! Write out the u observation
      ! Kind for u is
      write(iunit, *) 'RADIOSONDE_U_WIND_COMPONENT'
      write(iunit, *) 2
      write(iunit, *) pressure(level)
      write(iunit, *) lon
      write(iunit, *) lat
      write(iunit, *) 0, 0
      write(iunit, *) uv_err_var

      write(iunit, *) 0
      ! Write out the v observation
      ! Kind for v is 
      write(iunit, *) 'RADIOSONDE_V_WIND_COMPONENT'
      write(iunit, *) 2
      write(iunit, *) pressure(level)
      write(iunit, *) lon
      write(iunit, *) lat
      write(iunit, *) 0, 0
      write(iunit, *) uv_err_var

   end do
end do

write(iunit, *) 'set_def.out'

close(iunit)

end program pressure_col_rand
