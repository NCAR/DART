program dart_to_aether

use       utilities_mod, only : initialize_utilities, finalize_utilities
use transform_state_mod, only : initialize_transform_state_mod, dart_to_model, finalize_transform_state_mod

implicit none

call initialize_utilities(progname='dart_to_aether')

call initialize_transform_state_mod()

call dart_to_model()

call finalize_transform_state_mod()

call finalize_utilities('dart_to_aether')

end program dart_to_aether
