! DART software - Copyright 2004 - 2013 UCAR. This open source software is
! provided by UCAR, "as is", without charge, subject to all terms of use at
! http://www.image.ucar.edu/DAReS/DART/DART_download
!
! $Id$

!> \dir filter  Main program contained here
!> \file filter.f90 Main program

program filter

!> \mainpage filter Main DART Ensemble Filtering Program
!> @{ \brief routine to perform ensemble filtering
!>

use mpi_utilities_mod, only : initialize_mpi_utilities, finalize_mpi_utilities
use        filter_mod, only : filter_main

implicit none

!----------------------------------------------------------------

call initialize_mpi_utilities('Filter')

call filter_main()

call finalize_mpi_utilities()

!> @}

end program filter

! <next few lines under version control, do not edit>
! $URL$
! $Id$
! $Revision$
! $Date$
