! DART software - Copyright UCAR. This open source software is provided
! by UCAR, "as is", without charge, subject to all terms of use at
! http://www.image.ucar.edu/DAReS/DART/DART_download

program test_read_variable_namelist

use        utilities_mod, only : find_namelist_in_file, check_namelist_read 
use    mpi_utilities_mod, only : initialize_mpi_utilities,                &
                                 finalize_mpi_utilities
use            types_mod, only : vtablenamelength
use    default_model_mod, only : verify_state_variables  

!use test

implicit none

integer :: iunit, io

! DART state vector contents are specified in the input.nml:&model_nml namelist.
integer, parameter :: MAX_STATE_VARIABLES = 20
integer, parameter :: NUM_STATE_TABLE_COLUMNS = 3

character(len=vtablenamelength) :: state_variables(MAX_STATE_VARIABLES * NUM_STATE_TABLE_COLUMNS ) = ' '
integer :: nfields ! number of fields in the state vector
character(len=vtablenamelength) :: variable_table(MAX_STATE_VARIABLES, NUM_STATE_TABLE_COLUMNS)
integer :: state_qty_list(MAX_STATE_VARIABLES)
logical :: update_var_list(MAX_STATE_VARIABLES)

namelist /model_nml/  &
   state_variables

call initialize_mpi_utilities('test_read_variable_namelist')

call find_namelist_in_file('input.nml', 'model_nml', iunit)
read(iunit, nml = model_nml, iostat = io)
call check_namelist_read(iunit, io, 'model_nml')

call verify_state_variables(state_variables, nfields, variable_table, state_qty_list, update_var_list)

call finalize_mpi_utilities()

end program test_read_variable_namelist
