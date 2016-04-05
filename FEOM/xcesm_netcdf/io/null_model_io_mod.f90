! DART software - Copyright 2004 - 2013 UCAR. This open source software is
! provided by UCAR, "as is", without charge, subject to all terms of use at
! http://www.image.ucar.edu/DAReS/DART/DART_download
!
! $Id$

module null_model_io_mod

use time_manager_mod, only : time_type

implicit none

contains

end function get_info_file_name

!--------------------------------------------------------------------
!> read the time from the input file
function get_model_time_from_file(filename)

character(len=1024), intent(in) :: filename

type(time_type) :: get_model_time_from_file

end function get_model_time_from_file

!--------------------------------------------------------------------
!> construct restart file name for reading
function construct_file_name_in(stub, domain, copy)

character(len=512), intent(in) :: stub
integer,            intent(in) :: domain
integer,            intent(in) :: copy
character(len=1024)            :: construct_file_name_in

end function construct_file_name_in


!--------------------------------------------------------------------

end module null_model_io_mod
