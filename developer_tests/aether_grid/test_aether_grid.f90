program test_aether_grid

use mpi_utilities_mod, only : initialize_mpi_utilities, finalize_mpi_utilities
use model_mod,         only : test_grid_box
use assim_model_mod,   only : static_init_assim_model

call initialize_mpi_utilities('test_aether_grid')

call static_init_assim_model()

call test_grid_box

call finalize_mpi_utilities

end program test_aether_grid
