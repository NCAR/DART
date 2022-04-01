! DART software - Copyright UCAR. This open source software is provided
! by UCAR, "as is", without charge, subject to all terms of use at
! http://www.image.ucar.edu/DAReS/DART/DART_download

! Converts from MITgcm to DART-netcdf :
!  Inputs:
!       data        !GRID INFORMATION IS READ FROM HERE
!  15 variables >>
!              PSAL.data
!              PTMP.data
!              UVEL.data
!              VVEL.data
!              ETA.data
!              DIC.data
!              ALK.data
!              O2.data
!              NO3.data
!              PO4.data
!              FET.data
!              DON.data
!              DOP.data
!              PHY.data
!              CHL.data
!  Outputs:
!       OUTPUT.nc  ! Dart state netcdf file

program mit_to_dart

use    utilities_mod,     only : initialize_utilities, finalize_utilities
use    trans_mitdart_mod, only : static_init_trans, mit2dart

implicit none

call initialize_utilities(progname='mit_to_dart')

call static_init_trans()

call mit2dart()

call finalize_utilities('mit_to_dart')

end program mit_to_dart
