program id_set_def_stdin

use model_mod, only : static_init_model, get_model_size, &
   TYPE_PS, get_state_meta_data
use location_mod, only : location_type

implicit none

integer :: i, model_size, var_type
type(location_type) :: location

! Write to file
open(unit = 20, file = 'id_set_def_stdin.out')


! Get the model size
call static_init_model()
model_size = get_model_size()

! Set the output file name and a single set
write(20, *) 'set_def.out'
write(20, *) 1

! Set the number of state variables, all observed
write(20, *) model_size

! Loop through all the state variables, set the obs variance accordingly
do i = 1, model_size
   call get_state_meta_data(i, location, var_type)
! Output the appropriate observational error variance
   if(var_type == TYPE_PS) then
   write(20, *)  0.01
   else
      write(20, *) 0.00001
   endif
! Output the variable index
   write(20, *) i
end do

end program id_set_def_stdin
