! DART software - Copyright UCAR. This open source software is provided
! by UCAR, "as is", without charge, subject to all terms of use at
! http://www.image.ucar.edu/DAReS/DART/DART_download
!
! $Id$

program dart_to_model

use        types_mod, only : r8
use    utilities_mod, only : initialize_utilities
use        model_mod, only : get_model_size, static_init_model
use  assim_model_mod, only : assim_model_type, aread_state_restart, &
                             open_restart_read, close_restart
use time_manager_mod, only : time_type, get_time

implicit none

! version controlled file description for error handling, do not edit
character(len=256), parameter :: source   = &
   "$URL$"
character(len=32 ), parameter :: revision = "$Revision$"
character(len=128), parameter :: revdate  = "$Date$"

type(time_type)        :: model_time, target_time
real(r8), allocatable  :: x_state(:)
integer                :: file_unit, x_size, ens_member, i, day, sec
character (len = 128)  :: file_out = 'model.in', file_in = 'temp_ic'

call initialize_utilities(progname='dart_to_model', output_flag=.true.)

call static_init_model()        ! reads input.nml, etc., sets the table 
x_size = get_model_size()       ! now that we know how big state vector is ...
allocate(x_state(x_size))       ! allocate space for the (empty) state vector

file_unit = open_restart_read(file_in)
call aread_state_restart(model_time, x_state, file_unit, target_time)
call close_restart(file_unit)

call get_time( model_time, sec, day )
open(11,file=trim(file_out),status='new')
write(11,*) sec, day
call get_time( target_time, sec, day )
write(11,*) sec, day
do i = 1, x_size
   write(11,*)x_state(i)
enddo
close(11)
deallocate(x_state)

end program dart_to_model

! <next few lines under version control, do not edit>
! $URL$
! $Id$
! $Revision$
! $Date$
