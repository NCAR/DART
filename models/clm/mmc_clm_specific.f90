! DART software - Copyright UCAR. This open source software is provided
! by UCAR, "as is", without charge, subject to all terms of use at
! http://www.image.ucar.edu/DAReS/DART/DART_download
!
! $Id$

module mmc_clm_specific

!> the basic model_mod_check functions are in a common version now.
!> this module has only the clm specific parts which haven't been
!> incorporated into mmc yet.


use        types_mod, only : r8, digits12, metadatalength
use    utilities_mod, only : initialize_utilities, nc_check, &
                             open_file, close_file, find_namelist_in_file, &
                             check_namelist_read, finalize_utilities, &
                             error_handler, E_MSG
use     location_mod, only : location_type, set_location, write_location, get_dist, &
                             query_location, LocationDims, get_location, VERTISHEIGHT
use     obs_kind_mod, only : get_name_for_quantity, get_index_for_quantity, &
                             QTY_SNOWCOVER_FRAC, QTY_SOIL_TEMPERATURE
use  assim_model_mod, only : open_restart_read, open_restart_write, close_restart, &
                             aread_state_restart, awrite_state_restart, &
                             netcdf_file_type, aoutput_diagnostics, &
                             init_diag_output, finalize_diag_output
use time_manager_mod, only : time_type, set_calendar_type, GREGORIAN, &
                             read_time, get_time, set_time,  &
                             print_date, get_date, &
                             print_time, write_time, &
                             operator(-)
use        model_mod, only : static_init_model, get_model_size, get_state_meta_data, &
                             compute_gridcell_value, gridcell_components, &
                             model_interpolate, DART_get_var, get_grid_vertval

implicit none
private

! version controlled file description for error handling, do not edit
character(len=256), parameter :: source   = &
   "$URL$"
character(len=32 ), parameter :: revision = "$Revision$"
character(len=128), parameter :: revdate  = "$Date$"

!------------------------------------------------------------------
! The namelist variables
!------------------------------------------------------------------

public :: clm_mmc_checks

contains

!----------------------------------------------------------------------
!----------------------------------------------------------------------

subroutine clm_mmc_checks

if (test1thru > 7) then
   call gridcell_components( kind_of_interest )
endif

!----------------------------------------------------------------------
! Checking if the compute_gridcell_value works
!----------------------------------------------------------------------

if (test1thru > 8) then
   write(*,*)
   write(*,*)'Testing compute_gridcell_value() with QTY_SNOWCOVER_FRAC ...'

   loc = set_location(loc_of_interest(1), loc_of_interest(2), loc_of_interest(3), VERTISHEIGHT)

   call compute_gridcell_value(statevector, loc, QTY_SNOWCOVER_FRAC, interp_val, ios_out)

   if ( ios_out == 0 ) then
      write(*,*)'compute_gridcell_value : value is ',interp_val
   else
      write(*,*)'compute_gridcell_value : value is ',interp_val,'with error code',ios_out
   endif


   write(*,*)
   write(*,*)'Testing get_grid_vertval() with QTY_SOIL_TEMPERATURE ...'

   loc = set_location(loc_of_interest(1), loc_of_interest(2), loc_of_interest(3), VERTISHEIGHT)

   call get_grid_vertval(statevector, loc, QTY_SOIL_TEMPERATURE, interp_val, ios_out)

   if ( ios_out == 0 ) then
      write(*,*)'get_grid_vertval : value is ',interp_val
   else
      write(*,*)'get_grid_vertval : value is ',interp_val,'with error code',ios_out
   endif

endif

end subroutine clm_mmc_checks

end module mmc_clm_specific

! <next few lines under version control, do not edit>
! $URL$
! $Id$
! $Revision$
! $Date$
