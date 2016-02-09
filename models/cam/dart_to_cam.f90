! DART software - Copyright 2004 - 2013 UCAR. This open source software is
! provided by UCAR, "as is", without charge, subject to all terms of use at
! http://www.image.ucar.edu/DAReS/DART/DART_download
!
! $Id$

program dart_to_cam

!----------------------------------------------------------------------
! purpose: interface between CAM and DART
!
! method: Read DART state vector ("proprietary" format), including time(s).
!         Reform state vector back into CAM fields.
!         Replace those fields on the CAM initial file with the new values,
!         preserving all other information on the file.
!         Write a 'times' file containing the time information needed to
!         update the CAM namelist start/stop times.
!
!----------------------------------------------------------------------

use       types_mod,     only : r8
use   utilities_mod,     only : open_file, close_file, &
                                initialize_utilities, finalize_utilities, &
                                logfileunit, nmlfileunit, do_nml_file, do_nml_term, &
                                check_namelist_read, find_namelist_in_file
use       model_mod,     only : model_type, init_model_instance, write_cam_init, &
                                vector_to_prog_var, static_init_model, get_model_size, &
                                write_cam_times, end_model_instance
use state_vector_io_mod, only : aread_state_restart, open_restart_read, close_restart
use time_manager_mod,    only : time_type, print_time, print_date

implicit none

! version controlled file description for error handling, do not edit
character(len=256), parameter :: source   = &
   "$URL$"
character(len=32 ), parameter :: revision = "$Revision$"
character(len=128), parameter :: revdate  = "$Date$"

!------------------------------------------------------------------
! The namelist variables
!------------------------------------------------------------------

character (len = 128) :: dart_to_cam_input_file  = 'dart_ics'
character (len = 128) :: dart_to_cam_output_file = 'caminput.nc'
logical               :: advance_time_present    = .true.

namelist /dart_to_cam_nml/ dart_to_cam_input_file,  &
                           dart_to_cam_output_file, &
                           advance_time_present

!----------------------------------------------------------------------

type(model_type)       :: var
type(time_type)        :: model_time, adv_to_time
real(r8), allocatable  :: statevector(:)
integer                :: file_unit, vecsize, iunit, io

!-----------------------------------------------------------------------------
! start of program
!-----------------------------------------------------------------------------

call initialize_utilities('dart_to_cam')

! Read the namelist to get the input filename and whether
! a simple restart/ic file or a model advance file.
call find_namelist_in_file("input.nml", "dart_to_cam_nml", iunit)
read(iunit, nml = dart_to_cam_nml, iostat = io)
call check_namelist_read(iunit, io, "dart_to_cam_nml")

! Record the namelist values used for the run 
if (do_nml_file()) write(nmlfileunit, nml=dart_to_cam_nml)
if (do_nml_term()) write(     *     , nml=dart_to_cam_nml)

! initialize model
call static_init_model()

vecsize = get_model_size()
allocate(statevector(vecsize))

! Allocate the instance of the cam model type for storage
call init_model_instance(var)

! Get file for DART vector input
file_unit = open_restart_read(dart_to_cam_input_file)
! read in model time and state vector from DART
if (advance_time_present) then
   call aread_state_restart(model_time, statevector, file_unit, adv_to_time)
else
   call aread_state_restart(model_time, statevector, file_unit)
endif
call close_restart(file_unit)

! decompose vector back into CAM fields
call vector_to_prog_var (statevector, var)
deallocate (statevector)

! write fields to the netCDF initial file.
call write_cam_init(dart_to_cam_output_file, model_time, var)
call end_model_instance(var)

! write cam times to a separate 'times' support file.  used to
! update cam namelist start/stop times.
if (advance_time_present) call write_cam_times(model_time, adv_to_time)

call print_date( model_time,'dart_to_cam:CAM  model date')
call print_time( model_time,'dart_to_cam:DART model time')
call print_date( model_time,'dart_to_cam:CAM  model date',logfileunit)
call print_time( model_time,'dart_to_cam:DART model time',logfileunit)

if ( advance_time_present ) then
call print_time(adv_to_time,'dart_to_cam:advance_to time')
call print_date(adv_to_time,'dart_to_cam:advance_to date')
call print_time(adv_to_time,'dart_to_cam:advance_to time',logfileunit)
call print_date(adv_to_time,'dart_to_cam:advance_to date',logfileunit)
endif


call finalize_utilities('dart_to_cam')

end program dart_to_cam

! <next few lines under version control, do not edit>
! $URL$
! $Id$
! $Revision$
! $Date$
