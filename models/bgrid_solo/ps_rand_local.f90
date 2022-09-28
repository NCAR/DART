! DART software - Copyright UCAR. This open source software is provided
! by UCAR, "as is", without charge, subject to all terms of use at
! http://www.image.ucar.edu/DAReS/DART/DART_download
!
! $Id$

program ps_rand_local

use      types_mod, only : r8, PI
use  utilities_mod, only : get_unit, error_handler, E_ERR, initialize_utilities
use random_seq_mod, only : random_seq_type, init_random_seq, random_uniform
use   location_mod, only : VERTISSURFACE

implicit none

! version controlled file description for error handling, do not edit
character(len=256), parameter :: source   = &
   "$URL$"
character(len=32 ), parameter :: revision = "$Revision$"
character(len=128), parameter :: revdate  = "$Date$"

integer  :: num_sets, level, num, num_done, iunit
real(r8) :: err_var, bot_lat, top_lat, bot_lon, top_lon
type(random_seq_type) :: r

! Initialize the random sequence
call init_random_seq(r)

! Set up constants
num_sets =  1
level    = -1

! Initializer to allow for use of functions and subroutines from utilities_mod
call initialize_utilities('ps_rand_local')

! Open an output file and write header info
iunit = get_unit()
open(unit = iunit, file = 'ps_rand.out')

write(*, *) 'input the number of observations'
read(*, *) num

write(*, *) 'input the obs error variance'
read(*, *) err_var

write(*, *) 'input a lower bound on latitude -90 to 90'
read(*, *) bot_lat
write(*, *) 'input an upper bound on latitude -90 to 90'
read(*, *) top_lat
write(*, *) 'input a lower bound on longitude: no wraparounds for now '
read(*, *) bot_lon
write(*, *) 'input an upper bound on longitude '
read(*, *) top_lon

! Simple error check to let people know limits
if(top_lat <= bot_lat .or. top_lon <= bot_lon) then
   call error_handler(E_ERR,'ps_rand_local', 'lat lon range error', source, revision, revdate)
endif

! Input number of obs
write(iunit, '(I5)') num
! No obs values or qc
write(iunit, '(I1)') 0
write(iunit, '(I1)') 0

num_done = 0
do while(num_done < num)
   ! There are more obs
   write(iunit, '(I1)') 0

   ! Kind is ps
   write(iunit, '(A)') 'RADIOSONDE_SURFACE_PRESSURE'

   ! Put this on model level -1
   write(iunit, '(I2)') VERTISSURFACE
   write(iunit, '(I2)') level

   ! Want randomly located in horizontal
   write(iunit, '(I2)') -1

   ! Input longitude and latitude bounds
   write(iunit, '(F6.1)') bot_lon
   write(iunit, '(F6.1)') top_lon
   write(iunit, '(F6.1)') bot_lat
   write(iunit, '(F6.1)') top_lat

   ! Time is 0 days and 0 seconds for create_obs_sequence base
   write(iunit, *) 0, 0

   ! Error variance
   write(iunit, '(F5.1)') err_var

   num_done = num_done + 1

end do

! File name default is set_def.out
write(iunit, *) 'set_def.out'

end program ps_rand_local

! <next few lines under version control, do not edit>
! $URL$
! $Id$
! $Revision$
! $Date$
