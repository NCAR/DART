! DART software - Copyright UCAR. This open source software is provided
! by UCAR, "as is", without charge, subject to all terms of use at
! http://www.image.ucar.edu/DAReS/DART/DART_download

! Converts from DART-netcdf to MITgcm:
!  Inputs:
!     data        ! GRID INFORMATION IS READ FROM HERE
!     INPUT.nc    ! PSAL,PTMP,UVEL,VVEL, and ETA, NO3, dart state is read from here
!
!  Outputs: 14 variables >>
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


program dart_to_mit

use    utilities_mod,     only : initialize_utilities, finalize_utilities
use    trans_mitdart_mod, only : static_init_trans, dart2mit

implicit none


call initialize_utilities(progname='dart_to_mit')

call static_init_trans()

call dart2mit()

call finalize_utilities('dart_to_mit')

end program dart_to_mit

