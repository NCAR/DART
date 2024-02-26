program aether_to_dart

use       utilities_mod, only : initialize_utilities, finalize_utilities
use transform_state_mod, only : read_namelist_and_command_line_arguments

implicit none

call initialize_utilities(progname='aether_to_dart')

call read_namelist_and_command_line_arguments()

call finalize_utilities('aether_to_dart')

end program aether_to_dart
