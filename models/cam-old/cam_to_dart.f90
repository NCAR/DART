! DART software - Copyright UCAR. This open source software is provided
! by UCAR, "as is", without charge, subject to all terms of use at
! http://www.image.ucar.edu/DAReS/DART/DART_download
!
! $Id$

program cam_to_dart

!----------------------------------------------------------------------
! purpose: interface between CAM and DART
!
! method: Read CAM 'initial' file (netCDF format) for model state and time.
!         Reform fields into a DART state vector.
!         Write out state vector in "proprietary" format for DART
!
! author: Tim Hoar 4/4/2011
!
!----------------------------------------------------------------------

use        types_mod,    only : r8
use    utilities_mod,    only : initialize_utilities, finalize_utilities,   &
                                check_namelist_read, find_namelist_in_file, &
                                nmlfileunit, do_nml_file, do_nml_term
use        model_mod,    only : model_type, init_model_instance, &
                                end_model_instance, &
                                prog_var_to_vector, read_cam_init
use  assim_model_mod,    only : static_init_assim_model, get_model_size

use state_vector_io_mod, only : open_restart_write, awrite_state_restart, &
                                close_restart
use time_manager_mod,    only : time_type

implicit none

! version controlled file description for error handling, do not edit
character(len=256), parameter :: source   = &
   "$URL$"
character(len=32 ), parameter :: revision = "$Revision$"
character(len=128), parameter :: revdate  = "$Date$"


!-----------------------------------------------------------------------
! namelist parameters with default values.
!-----------------------------------------------------------------------

character (len = 128) :: cam_to_dart_input_file  = 'caminput.nc'
character (len = 128) :: cam_to_dart_output_file = 'dart_ics'

namelist /cam_to_dart_nml/ cam_to_dart_input_file, cam_to_dart_output_file

! allocatable storage to read in a native format for cam state
real(r8), allocatable  :: statevector(:)
type(model_type)       :: var
type(time_type)        :: model_time
integer                :: iunit, x_size, io

call initialize_utilities('cam_to_dart')

! Read the namelist entry
call find_namelist_in_file("input.nml", "cam_to_dart_nml", iunit)
read(iunit, nml = cam_to_dart_nml, iostat = io)
call check_namelist_read(iunit, io, "cam_to_dart_nml")

! Record the namelist values 
if (do_nml_file()) write(nmlfileunit, nml=cam_to_dart_nml)
if (do_nml_term()) write(     *     , nml=cam_to_dart_nml)

! Static init assim model sets the output file format (binary/ascii)
! and calls static_init_model
call static_init_assim_model()

! Allocate the local state vector
x_size = get_model_size()
allocate(statevector(x_size))

! Allocate the instance of the cam model type for storage
! I'll just point to the space I need, not;    
call init_model_instance(var)

! Read the file cam state fragments into var;
! transform fields into state vector for DART

call read_cam_init(cam_to_dart_input_file, var, model_time)

call prog_var_to_vector(var, statevector)
call end_model_instance(var)

! write out state vector in "proprietary" format
iunit = open_restart_write(cam_to_dart_output_file)
call awrite_state_restart(model_time, statevector, iunit)
call close_restart(iunit)

call finalize_utilities('cam_to_dart')

end program cam_to_dart

! <next few lines under version control, do not edit>
! $URL$
! $Id$
! $Revision$
! $Date$
