! DART software - Copyright 2004 - 2013 UCAR. This open source software is
! provided by UCAR, "as is", without charge, subject to all terms of use at
! http://www.image.ucar.edu/DAReS/DART/DART_download
!
! $Id$

!----------------------------------------------------------------------

!> purpose: interface between TIEGCM and DART
!>
!> method: Read TIEGCM restart file (netCDF format).
!>         Reform fields into a DART state vector.
!>         Write out state vector in DART format.
!>

program model_to_dart

use        types_mod, only : r8
use    utilities_mod, only : get_unit, initialize_utilities, finalize_utilities, &
                             error_handler, E_MSG, &
                             find_namelist_in_file, check_namelist_read
use        model_mod, only : static_init_model, get_model_size, &
                             read_TIEGCM_restart, get_restart_file_name
!                            clamp_bounded_variables
use  assim_model_mod, only : open_restart_write, awrite_state_restart, close_restart
use time_manager_mod, only : time_type, print_date, print_time

implicit none

! version controlled file description for error handling, do not edit
character(len=256), parameter :: source   = &
   "$URL$"
character(len=32 ), parameter :: revision = "$Revision$"
character(len=128), parameter :: revdate  = "$Date$"

!-----------------------------------------------------------------------
! namelist parameters with default values.
!-----------------------------------------------------------------------

character(len=256) :: file_out = 'dart_ics'

namelist /model_to_dart_nml/ file_out

!-----------------------------------------------------------------------
! global storage
!-----------------------------------------------------------------------

character(len=256)    :: filename
type(time_type)       :: model_time
integer               :: iunit, io, x_size
real(r8), allocatable :: x_state(:)

!=======================================================================
! Start the program.
!=======================================================================

call initialize_utilities(progname='model_to_dart', output_flag=.true.)

! read the namelist to get output file name
call find_namelist_in_file("input.nml", "model_to_dart_nml", iunit)
read(iunit, nml = model_to_dart_nml, iostat = io)
call check_namelist_read(iunit, io, "model_to_dart_nml") ! closes, too.

! static_init_model sets the geometry, model size, etc.
call static_init_model()

x_size = get_model_size()
allocate(x_state(x_size))

! Read the TIEGCM state variables into var and set the model_time
! to reflect the valid time of the TIEGCM state.

filename = get_restart_file_name()

call read_TIEGCM_restart(trim(filename), x_state, model_time)

! write out state vector in DART format
iunit = open_restart_write(trim(file_out))
call awrite_state_restart(model_time, x_state, iunit)
call close_restart(iunit)
deallocate(x_state)

! write a little summary
call print_date(model_time, str='model_to_dart: tiegcm date')
call print_time(model_time, str='model_to_dart: DART   time')

call error_handler(E_MSG,'model_to_dart','finished successfully.',source,revision,revdate)
call finalize_utilities()

end program model_to_dart

! <next few lines under version control, do not edit>
! $URL$
! $Id$
! $Revision$
! $Date$
