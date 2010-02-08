! DART software - Copyright © 2004 - 2010 UCAR. This open source software is
! provided by UCAR, "as is", without charge, subject to all terms of use at
! http://www.image.ucar.edu/DAReS/DART/DART_download

program trans_pv_sv

! <next few lines under version control, do not edit>
! $URL$
! $Id$
! $Revision$
! $Date$

!----------------------------------------------------------------------
! purpose: interface between AM2 and DART
!
! method: Read AM2 'initial' file for model state, but not time (netCDF format).
!         Get target time from assim_model_state_ic (temp_ic).
!         Reform fields into a state vector.
!         Write out state vector in "proprietary" format for DART
!
! author: Patrick Hofmann, updated: 6/2/2008
!         based on prog_var_to_vector and vector_to_prog_var by Jeff Anderson
!
!----------------------------------------------------------------------

use        types_mod, only : r8
use    utilities_mod, only : get_unit, file_exist, &
                             initialize_utilities, finalize_utilities
use        model_mod, only : model_type, init_model_instance, end_model_instance, &
                             prog_var_to_vector, read_model_init
use  assim_model_mod, only : assim_model_type, static_init_assim_model, &
     init_assim_model, get_model_size , set_model_state_vector, write_state_restart, &
     set_model_time, open_restart_read, open_restart_write, close_restart, &
     aread_state_restart
use time_manager_mod, only : time_type, read_time, set_time, set_date

implicit none

integer, external :: iargc

character (len = 128) :: dartSVout, RstFileIn, TrcFileIn
character (len = 256) :: string1, string2

!----------------------------------------------------------------------

! Temporary allocatable storage to read in a native format for cam state
type(assim_model_type) :: x
type(model_type)       :: var
type(time_type)        :: model_time, adv_to_time
real(r8), allocatable  :: x_state(:)
integer                :: file_unit, x_size, big_cld_iw, small_trcs
integer                :: year, month, day, hour, minute, second
logical                :: do_output = .false.

if(iargc()  == 0) stop "You must specify State Vector and input AM2 files"
call getarg(1, dartSVout)
call getarg(2, RstFileIn)
call getarg(3, TrcFileIn)

call initialize_utilities('Trans_pv_sv')

if(file_exist('element1')) do_output = .true.

! Static init assim model calls static_init_model
! which now (merge/MPI) calls read_model_init)
call static_init_assim_model()

! Initialize the assim_model instance
call init_assim_model(x)

! Allocate the local state vector
x_size = get_model_size()
allocate(x_state(x_size))

! Allocate the instance of the AM2 model type for storage  
call init_model_instance(var)

! Read the file AM2 state fragments into var, but not time
call read_model_init(RstFileIn, TrcFileIn, var)

! Ensure that all tracers that are <=1e-10 get set to zero.
! Further, ensure that CF is <=1 and exit if CLW or CIW are >1e-1
! Lastly, output error message that says how many values are getting adjusted.
!if (any(var%tracers(:,:,:,1:2) > 1e-1)) then
!   big_cld_iw = count(var%tracers(:,:,:,1:2) > 1e-1)
!   print*, 'Stopping due to ', big_cld_iw, ' values of cloud ice and water > 1e-1'
   !stop
!endif
!small_trcs = count(var%tracers < 1e-10)
!print*, 'Number of tracer values < 1e-10 ', small_trcs

!where(var%tracers < 1e-10) var%tracers = 0
!where(var%tracers(:,:,:,3) > 1) var%tracers(:,:,:,3) = 1

! transform fields into state vector for DART
call prog_var_to_vector(var, x_state)

call end_model_instance(var)

! Put this in the structure
call set_model_state_vector(x, x_state)

! Get current model time from line 3 of coupler.res
open(50, file = 'coupler.res',form = 'formatted', action = 'read')
read(50,*) string1
!print*, string1
read(50,*) string2
!print*, string2
read(50,*) year, month, day, hour, minute, second 
close(50)

! Set model_time
model_time = set_date(year, month, day, hour, minute, second)
call set_model_time(x, model_time)

!call close_restart(file_unit)

! Get channel for output,
! write out state vector in "proprietary" format
file_unit = open_restart_write(dartSVout)
call write_state_restart(x, file_unit)
call close_restart(file_unit)

call finalize_utilities()

end program trans_pv_sv
