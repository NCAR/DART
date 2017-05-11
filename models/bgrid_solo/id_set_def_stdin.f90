! DART software - Copyright UCAR. This open source software is provided
! by UCAR, "as is", without charge, subject to all terms of use at
! http://www.image.ucar.edu/DAReS/DART/DART_download
!
! $Id$
 
program id_set_def_stdin

use     types_mod, only : r8, i8
use  location_mod, only : location_type
use utilities_mod, only : get_unit
use  obs_kind_mod, only : QTY_SURFACE_PRESSURE
use     model_mod, only : static_init_model, get_model_size, &
                          get_state_meta_data

implicit none

! version controlled file description for error handling, do not edit
character(len=256), parameter :: source   = &
   "$URL$"
character(len=32 ), parameter :: revision = "$Revision$"
character(len=128), parameter :: revdate  = "$Date$"

integer :: i, model_size, var_type, iunit
type(location_type) :: location

! Write to file
iunit  = get_unit()
open(unit = iunit, file = 'id_set_def_stdin.out')

! Get the model size
call static_init_model()
model_size = get_model_size()

! Set the number of state variables, all observed
write(iunit, *) model_size

! No values or qc
write(iunit, *) 0
write(iunit, *) 0

! Loop through all the state variables, set the obs variance accordingly
do i = 1, model_size
   ! There are more obs
   write(iunit, *) 0

   ! Identity obs
   write(iunit, *) -1 * i

   ! Time is 0 days 0 seconds for create obs sequence
   write(iunit, *) 0, 0

   call get_state_meta_data(int(i,i8), location, var_type)

   ! Output the appropriate observational error variance
   if(var_type == QTY_SURFACE_PRESSURE) then
   write(iunit, *)  10000.0_r8
   else
      write(iunit, *) 1.0_r8
   endif

end do

! Output the default set_def.out file name
write(iunit, *) 'set_def.out'

end program id_set_def_stdin

! <next few lines under version control, do not edit>
! $URL$
! $Id$
! $Revision$
! $Date$
