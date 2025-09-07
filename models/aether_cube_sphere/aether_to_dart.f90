program aether_to_dart

use       utilities_mod, only : initialize_utilities, finalize_utilities, &
                                find_namelist_in_file, check_namelist_read

use transform_state_mod, only : initialize_transform_state_mod, model_to_dart, &
                                get_ensemble_range_from_command_line

implicit none

character(len=256) :: aether_file_directory, dart_file_directory

namelist /aether_to_dart_nml / aether_file_directory, dart_file_directory

integer :: iunit, io
integer :: ens, start_ens, end_ens

!----------------------------------------------------------------

call initialize_utilities(progname='aether_to_dart')

! Read the namelist
call find_namelist_in_file('input.nml', 'aether_to_dart_nml', iunit)
read(iunit, nml = aether_to_dart_nml, iostat = io)
call check_namelist_read(iunit, io, 'aether_to_dart_nml')

call initialize_transform_state_mod()

call get_ensemble_range_from_command_line(start_ens, end_ens)

! The DART SE team has pointed out concerns about having the loop in the program
! In the long-term, this may need to be moved back to be a command line argument
! and the program will only tranform a single file. 
! Loop through the ensemble members and transform each
do ens = start_ens, end_ens
   call model_to_dart(aether_file_directory, dart_file_directory, ens)
end do

call finalize_utilities('aether_to_dart')

end program aether_to_dart
