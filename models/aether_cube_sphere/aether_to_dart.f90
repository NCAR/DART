program aether_to_dart

use       utilities_mod, only : initialize_utilities, finalize_utilities
use transform_state_mod, only : initialize_transform_state_mod, model_to_dart, finalize_transform_state_mod

implicit none

call initialize_utilities(progname='aether_to_dart')

call initialize_transform_state_mod()

call model_to_dart()

call finalize_transform_state_mod()

call finalize_utilities('aether_to_dart')

end program aether_to_dart
