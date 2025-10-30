! DART software - Copyright UCAR. This open source software is provided
! by UCAR, "as is", without charge, subject to all terms of use at
! http://www.image.ucar.edu/DAReS/DART/DART_download

! Convert DART filter files to Aether block restart files

program dart_to_aether

use       utilities_mod, only : initialize_utilities, finalize_utilities, &
                                find_namelist_in_file, check_namelist_read

use transform_state_mod, only : initialize_transform_state_mod, dart_to_model, &
                                get_ensemble_range_from_command_line

implicit none

character(len=256) :: dart_file_directory, aether_file_directory

namelist /dart_to_aether_nml / dart_file_directory, aether_file_directory

integer :: iunit, io
integer :: ens, start_ens, end_ens

!----------------------------------------------------------------

call initialize_utilities(progname='dart_to_aether')

! Read the namelist
call find_namelist_in_file('input.nml', 'dart_to_aether_nml', iunit)
read(iunit, nml = dart_to_aether_nml, iostat = io)
call check_namelist_read(iunit, io, 'dart_to_aether_nml')

call initialize_transform_state_mod()

! Convert for a range of ensemble members
call get_ensemble_range_from_command_line(start_ens, end_ens)

! The DART SE team has pointed out concerns about having the loop in the program
! Loop through the ensemble members and transform each
do ens = start_ens, end_ens
   call dart_to_model(dart_file_directory, aether_file_directory, ens)
end do

call finalize_utilities('dart_to_aether')

end program dart_to_aether
