program dart_to_aether

use       utilities_mod, only : initialize_utilities, finalize_utilities, &
                                find_namelist_in_file, check_namelist_read

use transform_state_mod, only : initialize_transform_state_mod, dart_to_model

implicit none

character(len=256) :: dart_file_directory, aether_file_directory

namelist /dart_to_aether_nml / dart_file_directory, aether_file_directory

integer :: iunit, io

!----------------------------------------------------------------

call initialize_utilities(progname='dart_to_aether')

! Read the namelist
call find_namelist_in_file('input.nml', 'dart_to_aether_nml', iunit)
read(iunit, nml = dart_to_aether_nml, iostat = io)
call check_namelist_read(iunit, io, 'dart_to_aether_nml')

call initialize_transform_state_mod()

call dart_to_model(dart_file_directory, aether_file_directory)

call finalize_utilities('dart_to_aether')

end program dart_to_aether
