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

! Version string captured at build time from git describe --tags --dirty
! This will be set by the build system via preprocessor definition
#ifndef DART_VERSION
#define DART_VERSION "unknown"
#endif

character(len=*), parameter :: dart_version_string = DART_VERSION

public :: get_dart_version

contains

!-----------------------------------------------------------------------
!> Return the DART version string captured at build time
!>
!> @return character string containing git version info like "v11.14.2-1-g375abbe" 
!>         or "v11.14.2-1-g375abbe-dirty" if there were uncommitted changes

function get_dart_version() result(version)

character(len=len(dart_version_string)) :: version

version = dart_version_string

end function get_dart_version

!-----------------------------------------------------------------------

end module version_mod