! DART software - Copyright UCAR. This open source software is provided
! by UCAR, "as is", without charge, subject to all terms of use at
! http://www.image.ucar.edu/DAReS/DART/DART_download

program test_parse_variables

use        utilities_mod, only : find_namelist_in_file, check_namelist_read, &
                                 initialize_utilities, finalize_utilities
use            types_mod, only : vtablenamelength, MISSING_R8, r8
use    default_model_mod, only : parse_variables_clamp, parse_variables, &
                                 state_var_type, &
                                 MAX_STATE_VARIABLE_FIELDS_CLAMP, &
                                 MAX_STATE_VARIABLE_FIELDS
use         obs_kind_mod, only : QTY_SALINITY, QTY_POTENTIAL_TEMPERATURE, &
                                 QTY_U_CURRENT_COMPONENT

use test ! fortran-testanything

implicit none

integer :: iunit, io
type(state_var_type) :: state_vars, state_vars_clamp

character(len=vtablenamelength) :: state_variables(MAX_STATE_VARIABLE_FIELDS) = ' '
character(len=vtablenamelength) :: state_variables_clamp(MAX_STATE_VARIABLE_FIELDS_CLAMP) = ' '

namelist /model_nml/  &
   state_variables

namelist /model_nml_clamp/  &
   state_variables_clamp

call initialize_utilities('test_parse_variables')

call plan(28)

! Using namelist entry WITHOUT clamping values

call find_namelist_in_file('input.nml', 'model_nml', iunit)
read(iunit, nml = model_nml, iostat = io)
call check_namelist_read(iunit, io, 'model_nml')

state_vars = parse_variables(state_variables)

call ok(state_vars%nvars == 3)
call ok(state_vars%netcdf_var_names(1) == 'SALT_CUR')
call ok(state_vars%netcdf_var_names(2) == 'TEMP_CUR')
call ok(state_vars%netcdf_var_names(3) == 'UVEL_CUR')
call ok(state_vars%qtys(1) == QTY_SALINITY)
call ok(state_vars%qtys(2) == QTY_POTENTIAL_TEMPERATURE)
call ok(state_vars%qtys(3) == QTY_U_CURRENT_COMPONENT)
call ok(allocated(state_vars_clamp%clamp_values) .eqv. .false.)
call ok(state_vars%updates(1) .eqv. .true.)
call ok(state_vars%updates(2) .eqv. .false.)
call ok(state_vars%updates(3) .eqv. .true.)

! Using namelist entry WITH clamping values

call find_namelist_in_file('input.nml', 'model_nml_clamp', iunit)
read(iunit, nml = model_nml_clamp, iostat = io)
call check_namelist_read(iunit, io, 'model_nml_clamp')

state_vars_clamp = parse_variables_clamp(state_variables_clamp)

call ok(state_vars_clamp%nvars == 3)
call ok(state_vars_clamp%netcdf_var_names(1) == 'SALT_CUR')
call ok(state_vars_clamp%netcdf_var_names(2) == 'TEMP_CUR')
call ok(state_vars_clamp%netcdf_var_names(3) == 'UVEL_CUR')
call ok(state_vars_clamp%qtys(1) == QTY_SALINITY)
call ok(state_vars_clamp%qtys(2) == QTY_POTENTIAL_TEMPERATURE)
call ok(state_vars_clamp%qtys(3) == QTY_U_CURRENT_COMPONENT)
call ok(allocated(state_vars_clamp%clamp_values) .eqv. .true.)
call ok(state_vars_clamp%clamp_values(1,1) == 0.0_r8)
call ok(state_vars_clamp%clamp_values(1,2) == 0.0_r8)
call ok(state_vars_clamp%clamp_values(2,1) == 1.0_r8, "1.0")
call ok(state_vars_clamp%clamp_values(2,2) == 2.0_r8, "integer '2'")
call ok(state_vars_clamp%clamp_values(3,1) == -4.0_r8, "integer '-4'")
call ok(state_vars_clamp%clamp_values(3,2) == 3.0_r8, "integer '3'")
call ok(state_vars_clamp%updates(1) .eqv. .true.)
call ok(state_vars_clamp%updates(2) .eqv. .false.)
call ok(state_vars_clamp%updates(3) .eqv. .true.)

call finalize_utilities()

end program test_parse_variables
