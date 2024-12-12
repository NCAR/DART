! DART software - Copyright UCAR. This open source software is provided
! by UCAR, "as is", without charge, subject to all terms of use at
! http://www.image.ucar.edu/DAReS/DART/DART_download

program test_get_state_variables

use        utilities_mod, only : find_namelist_in_file, check_namelist_read 
use    mpi_utilities_mod, only : initialize_mpi_utilities, finalize_mpi_utilities
use            types_mod, only : vtablenamelength, MISSING_R8
use    default_model_mod, only : get_state_variables, state_var_type

use test ! fortran-testanything

implicit none

integer :: iunit, io
logical :: use_clamping
type(state_var_type) :: state_vars

! DART state vector contents are specified in the input.nml:&model_nml namelist.
integer, parameter :: MAX_STATE_VARIABLES = 20
integer, parameter :: NUM_STATE_TABLE_COLUMNS = 3
character(len=vtablenamelength) :: state_variables(MAX_STATE_VARIABLES * NUM_STATE_TABLE_COLUMNS ) = ' '

namelist /model_nml/  &
   state_variables

call initialize_mpi_utilities('test_get_state_variables')

call plan(12)

! Using namelist file WITHOUT clamping values
use_clamping = .false.

call find_namelist_in_file('input.nml', 'model_nml', iunit)
read(iunit, nml = model_nml, iostat = io)
call check_namelist_read(iunit, io, 'model_nml')

call get_state_variables(state_variables, MAX_STATE_VARIABLES, use_clamping, state_vars)

write(*,*) 'state_vars%nvars: ', state_vars%nvars
write(*,*) 'state_vars%netcdf_var_names: ', state_vars%netcdf_var_names
write(*,*) 'state_vars%qtys: ', state_vars%qtys
write(*,*) 'state_vars%updates: ', state_vars%updates

call ok(state_vars%nvars == 5) !OK
call ok(state_vars%netcdf_var_names(2) == '') !NOT OK
call ok(state_vars%qtys(3) == 20) !NOT OK
call ok(state_vars%clamp_values(1,1) == MISSING_R8) !OK
call ok(state_vars%clamp_values(5,2) == 0.4) !NOT OK
call ok(state_vars%updates(4) == .true.) !OK

! Using namelist file WITH clamping values
use_clamping = .true.

call find_namelist_in_file('input.nml.clamp', 'model_nml', iunit)
read(iunit, nml = model_nml, iostat = io)
call check_namelist_read(iunit, io, 'model_nml')

call get_state_variables(state_variables, MAX_STATE_VARIABLES, use_clamping, state_vars)

write(*,*) 'state_vars%nvars: ', state_vars%nvars
write(*,*) 'state_vars%netcdf_var_names: ', state_vars%netcdf_var_names
write(*,*) 'state_vars%qtys: ', state_vars%qtys
write(*,*) 'state_vars%clamp_values(:,1): ', state_vars%clamp_values(:, 1)
write(*,*) 'state_vars%clamp_values(:,2): ', state_vars%clamp_values(:, 2)
write(*,*) 'state_vars%updates: ', state_vars%updates

call ok(state_vars%nvars == 7) !NOT OK
call ok(state_vars%netcdf_var_names(2) == 'TEMP_CUR') !OK
call ok(state_vars%qtys(3) == 356) !OK
call ok(state_vars%clamp_values(1,1) == 99) !NOT OK
call ok(state_vars%clamp_values(5,2) == 0.0) !OK
call ok(state_vars%updates(4) == .false.) !NOT OK

call finalize_mpi_utilities()

end program test_get_state_variables
