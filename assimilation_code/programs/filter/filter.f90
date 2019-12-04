! DART software - Copyright UCAR. This open source software is provided
! by UCAR, "as is", without charge, subject to all terms of use at
! http://www.image.ucar.edu/DAReS/DART/DART_download
!
! $Id$

!> Main DART Ensemble Filtering Program

program filter

use mpi_utilities_mod, only : initialize_mpi_utilities, finalize_mpi_utilities
use        filter_mod, only : filter_main

implicit none

!----------------------------------------------------------------

call initialize_mpi_utilities('Filter')

call filter_main()

call finalize_mpi_utilities()

end program filter

