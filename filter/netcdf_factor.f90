program netcdf_factor

!
! <next four lines automatically updated by CVS, do not edit>
! $Source$
! $Revision$
! $Date$
! $Author$
!

! Creates a netcdf file with the factors for all state variables and
! a given observation for the bgrid model. Can then plot with standard
! based matlab routines (plot_bgrid_factor.m).

use        types_mod, only : r8
use time_manager_mod, only : time_type, set_time
use    utilities_mod, only : get_unit, error_handler, FATAL
use  assim_model_mod, only : assim_model_type, static_init_assim_model, &
   init_diag_output, get_model_size, output_diagnostics, &
   init_assim_model, get_state_vector_ptr, set_model_time

use netcdf

implicit none
! Define a type for doing direct access to ensemble state vectors

! CVS Generated file description for error handling, do not edit
character(len=128) :: &
source   = "$Source$", &
revision = "$Revision$", &
revdate  = "$Date$"

type model_state_ptr_type
   private
   real(r8), pointer :: state(:)
end type model_state_ptr_type

type(time_type) :: time

integer :: prior_state_unit, meta_data_size, model_size, ierr
integer :: PriorStateUnit, i, time_index, obs_index, state_index
integer :: j, num_times, curr_time, ioerr, iunit
real(r8) :: factor_in
character(len = 129) :: copy_meta_data(1)
type(assim_model_type) :: factor
type(model_state_ptr_type) :: factor_ptr

! Need to get the model size, etc.
call static_init_assim_model()
model_size = get_model_size()

! Initialize the assim_model_data structure
call init_assim_model(factor)
factor_ptr%state => get_state_vector_ptr(factor)

meta_data_size = 1
copy_meta_data(1) = 'only copy'


PriorStateUnit = init_diag_output('Factor', &
   'factor', meta_data_size, copy_meta_data)


! Open the factor file for reading
! Read in the first time of the factor file
factor_ptr%state(:) = 0.0
! Read in the first entry
iunit = get_unit()
read(iunit, *) time_index, obs_index, state_index, factor_in
curr_time = time_index

do    
      factor_ptr%state(state_index) = factor_in
      read(iunit, *, IOSTAT = ioerr) time_index, obs_index, state_index, factor_in 
      if(time_index > curr_time .or. ioerr /= 0) then
         time = set_time(curr_time, 0)
         call set_model_time(factor, time)
         call output_diagnostics(PriorStateUnit, factor, 1)
         write(*, *) 'writing output'
         write(*, *) 'current and time ', curr_time, time_index
         factor_ptr%state(:) = 0.0
         curr_time = time_index
      endif
      if(ioerr /= 0) goto 111
end do

111 continue
ierr = NF90_close(PriorStateUnit)

end program netcdf_factor
