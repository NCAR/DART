! DART software - Copyright UCAR. This open source software is provided
! by UCAR, "as is", without charge, subject to all terms of use at
! http://www.image.ucar.edu/DAReS/DART/DART_download
!
! $Id$

program model_to_dart

use        types_mod, only : r8
use    utilities_mod, only : initialize_utilities
use        model_mod, only : get_model_size, static_init_model
use  assim_model_mod, only : assim_model_type, awrite_state_restart, &
                             open_restart_write, close_restart
use time_manager_mod, only : time_type, set_time


implicit none

! version controlled file description for error handling, do not edit
character(len=256), parameter :: source   = &
   "$URL$"
character(len=32 ), parameter :: revision = "$Revision$"
character(len=128), parameter :: revdate  = "$Date$"

type(time_type)        :: model_time, target_time
real(r8), allocatable  :: x_state(:)
integer                :: file_unit, x_size, ens_member, i, day, sec
character (len = 128)  :: file_in = 'model.out', file_out = 'temp_ud'

!----------------------------------------------------------------------
!----------------------------------------------------------------------
call initialize_utilities(progname='model_to_dart', output_flag=.true.)

call static_init_model()        ! reads input.nml, etc., sets the table 
x_size = get_model_size()       ! now that we know how big state vector is ...
allocate(x_state(x_size))       ! allocate space for the (empty) state vector

open(11,file=trim(file_in),status='old')
read(11,*) sec, day
do i = 1, x_size
   read(11,*)x_state(i)
enddo
close(11)
model_time = set_time( sec, day )

file_unit = open_restart_write(file_out,"unformatted")
call awrite_state_restart(model_time, x_state, file_unit)
call close_restart(file_unit)

deallocate(x_state)

end program model_to_dart

! <next few lines under version control, do not edit>
! $URL$
! $Id$
! $Revision$
! $Date$
