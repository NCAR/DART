! DART software - Copyright UCAR. This open source software is provided
! by UCAR, "as is", without charge, subject to all terms of use at
! http://www.image.ucar.edu/DAReS/DART/DART_download

!> Module to provide version information for DART programs
!>
!> This module stores the version information that is captured at build time
!> from git describe and provides a function to access it.

module version_mod

implicit none
private 

! Assign a preprocessor variable for the version, which we'll get using git
#ifndef DART_VERSION
#define DART_VERSION 'unknown'
#endif

character(len=*), parameter :: dart_version_string = DART_VERSION

public :: get_dart_version

contains

! Return DART version in string
function get_dart_version() result(version)

character(len=len(dart_version_string)) :: version

version = dart_version_string

end function get_dart_version


end module version_mod


