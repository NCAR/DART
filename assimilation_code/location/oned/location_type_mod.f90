! DART software - Copyright UCAR. This open source software is provided
! by UCAR, "as is", without charge, subject to all terms of use at
! http://www.image.ucar.edu/DAReS/DART/DART_download

module location_type_mod

use types_mod, only : r8

type location_type
   real(r8) :: x
end type location_type

public :: location_type

end module

