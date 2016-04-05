! DART software - Copyright 2004 - 2013 UCAR. This open source software is
! provided by UCAR, "as is", without charge, subject to all terms of use at
! http://www.image.ucar.edu/DAReS/DART/DART_download
!
! $Id: io_null_netcdf_mod.f90 7314 2014-12-23 20:11:42Z hkershaw $

module null_model_io_mod

use time_manager_mod, only : time_type

implicit none

contains

!--------------------------------------------------------------------
!> read namelist and set up filename arrays
function fill_variable_list(num_variables_in_state)

integer             :: num_variables_in_state
character(len=256)  :: fill_variable_list(num_variables_in_state)


end function fill_variable_list

!--------------------------------------------------------------------
!> get the input file name
function get_info_file_name(domain)

integer, intent(in) :: domain
character(len=265)  :: get_info_file_name

end function get_info_file_name

!--------------------------------------------------------------------
!> read the time from the input file
function get_model_time_from_file(filename)

character(len=1024), intent(in) :: filename

type(time_type) :: get_model_time_from_file

end function get_model_time_from_file

!--------------------------------------------------------------------
!> pass number of variables in the state out to filter 
subroutine variables_domains(num_variables_in_state, num_doms)

integer, intent(out) :: num_variables_in_state
integer, intent(out) :: num_doms !< number of domains

end subroutine variables_domains

!--------------------------------------------------------------------
!> construct info filename for get_state_variable_info
function info_file_name(domain)

integer, intent(in) :: domain
character(len=256)  :: info_file_name

end function info_file_name

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
