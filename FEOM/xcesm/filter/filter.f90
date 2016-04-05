! DART software - Copyright 2004 - 2013 UCAR. This open source software is
! provided by UCAR, "as is", without charge, subject to all terms of use at
! http://www.image.ucar.edu/DAReS/DART/DART_download
!
! $Id$

!> contents moved to a module for more flexible use of the
!> filter routines.

program filter

use filter_mod, only : filter_main

call filter_main()

end program

