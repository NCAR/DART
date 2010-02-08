! DART software - Copyright © 2004 - 2010 UCAR. This open source software is
! provided by UCAR, "as is", without charge, subject to all terms of use at
! http://www.image.ucar.edu/DAReS/DART/DART_download

program trans_sv_pv

! <next few lines under version control, do not edit>
! $URL$
! $Id$
! $Revision$
! $Date$

!----------------------------------------------------------------------
! purpose: interface between CAM and DART
!
! method: Read DART state vector ("proprietary" format), but not time(s).
!         Reform state vector back into CAM fields.
!         Replace those fields on the CAM initial file with the new values,
!         preserving all other information on the file.
!
! author: Kevin Raeder 2/21/03
!         based on prog_var_to_vector and vector_to_prog_var by Jeff Anderson
! mod:    to read temp_ic (assim_model_state_ic; 2 times) or temp_ud (1 time) and put
!         the fields into the CAM initial file
!
!----------------------------------------------------------------------

use       types_mod, only : r8
use   utilities_mod, only : get_unit, file_exist, open_file, &
                            initialize_utilities, finalize_utilities
use       model_mod, only : model_type, init_model_instance, write_cam_init, &
   vector_to_prog_var 
use assim_model_mod, only : assim_model_type, static_init_assim_model, &
   init_assim_model, get_model_size, get_model_state_vector, read_state_restart, &
   open_restart_read, close_restart
use time_manager_mod, only : time_type, read_time

implicit none

! version controlled file description for error handling, do not edit
character(len=128), parameter :: &
   source   = "$URL$", &
   revision = "$Revision$", &
   revdate  = "$Date$"

type(assim_model_type) :: x
type(model_type)       :: var
! Guam clean out advance_model
type(time_type)        :: adv_to_time
real(r8), allocatable  :: x_state(:)
integer                :: file_unit, mem_unit, x_size
character (len = 128)  :: file_name, file_in 
logical                :: do_output = .false.

call initialize_utilities('Trans_sv_pv')

if(file_exist('element1')) do_output = .true.

! Static init assim model calls static_init_model
if (do_output) then
   WRITE(*,'(////A)') '========================================================================='
   PRINT*,'static_init_assim_model in trans_sv_pv'
endif
call static_init_assim_model()
call init_assim_model(x)

! Allocate the instance of the cam model type for storage
call init_model_instance(var)

if (file_exist( 'temp_ic' )) then
   file_in = 'temp_ic'
   file_name = 'caminput.nc'
   ! Get file for DART vector input
   file_unit = open_restart_read(file_in)
   ! read in target time and state vector from DART and throwing away the time(s)
   ! since those are handled by trans_time.f90
   call read_state_restart(x, file_unit, adv_to_time)
else if (file_exist( 'member' )) then
   mem_unit = open_file ('member')
   read(mem_unit,'(A)') file_in
   read(mem_unit,'(A)') file_name
   PRINT*,' file_in = ',file_in
   file_unit = open_restart_read(file_in)
   ! read state vector from DART and throw away the time
   call read_state_restart(x, file_unit)
endif
call close_restart(file_unit)

! Get the state part of the assim_model type x
x_size = get_model_size()
allocate(x_state(x_size))
x_state = get_model_state_vector(x)

! decompose vector back into CAM fields
call vector_to_prog_var (x_state, var)
deallocate (x_state)

! write fields to the netCDF initial file
! merge/MPI; this requires no change; a CAM state will exist in model_mod,
!            but this will ignore it and write out *this* CAM state.
call write_cam_init(file_name, var)

call finalize_utilities()

end program trans_sv_pv
