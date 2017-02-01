! DART software - Copyright UCAR. This open source software is provided
! by UCAR, "as is", without charge, subject to all terms of use at
! http://www.image.ucar.edu/DAReS/DART/DART_download
!
! $Id$

program am2_to_dart

!----------------------------------------------------------------------
! purpose: interface between AM2 and DART
!
! method: Read a namelist for run-time information.
!         Read AM2 'restart' file for model state
!         Read AM2 coupler.res file for model time (for some reason)
!         Reform fields into a DART state vector.
!         Write out state vector in "proprietary" format for DART
!
! author: Tim Hoar 31 May 2013
!         based on previous trans_pv_sv.f90 which used non-portable
!         command-line arguments. 
!
!----------------------------------------------------------------------

use        types_mod, only : r8
use    utilities_mod, only : initialize_utilities, finalize_utilities, &
                             nmlfileunit, do_nml_file, do_nml_term, &
                             find_namelist_in_file, check_namelist_read, &
                             open_file, close_file, logfileunit
use        model_mod, only : model_type, static_init_model, init_model_instance, &
                             end_model_instance, prog_var_to_vector, read_model_init
use  assim_model_mod, only : get_model_size, awrite_state_restart, &
                             open_restart_write, close_restart
use time_manager_mod, only : time_type, set_time, set_date, &
                             print_time, print_date

implicit none

! version controlled file description for error handling, do not edit
character(len=256), parameter :: source   = &
   "$URL$"
character(len=32 ), parameter :: revision = "$Revision$"
character(len=128), parameter :: revdate  = "$Date$"

!-----------------------------------------------------------------------
! namelist parameters with default values.
!-----------------------------------------------------------------------

character(len=128) :: restart_file = 'fv_rst.res.nc'
character(len=128) :: tracer_file  = 'atmos_tracers.res.nc'
character(len=128) :: am2_to_dart_output_file  = 'dart_ics'

namelist /am2_to_dart_nml/ restart_file, tracer_file, am2_to_dart_output_file

!----------------------------------------------------------------------
! global storage
!----------------------------------------------------------------------

character(len=256)     :: string1
type(model_type)       :: var
type(time_type)        :: model_time
real(r8), allocatable  :: statevector(:)
integer                :: iunit, io, x_size, big_cld_iw, small_trcs
integer                :: year, month, day, hour, minute, second

!------------------------------------------------------------------

call initialize_utilities('am2_to_dart')

! Read the namelist information
call find_namelist_in_file("input.nml", "am2_to_dart_nml", iunit)
read(iunit, nml = am2_to_dart_nml, iostat = io)
call check_namelist_read(iunit, io, "am2_to_dart_nml")

! Record the namelist values 
if (do_nml_file()) write(nmlfileunit, nml=am2_to_dart_nml)
if (do_nml_term()) write(     *     , nml=am2_to_dart_nml)

write(*,*)
write(*,'(''am2_to_dart:converting am2 restart file '',A, &
      &'' to DART file '',A)') &
       trim(restart_file), trim(am2_to_dart_output_file)

!----------------------------------------------------------------------
! Get to work
!----------------------------------------------------------------------

call static_init_model()
x_size = get_model_size()
allocate(statevector(x_size))

! Allocate the instance of the AM2 model type for storage  
call init_model_instance(var)

! Read the file AM2 state fragments into var, but not time
call read_model_init(restart_file, tracer_file, var)

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
call prog_var_to_vector(var, statevector)
call end_model_instance(var)

! Get current model time from line 3 of coupler.res
iunit = open_file('coupler.res',form = 'formatted', action = 'read')
read(iunit,*) string1
read(iunit,*) string1
read(iunit,*) year, month, day, hour, minute, second 
call close_file(iunit)

! Set model_time
model_time = set_date(year, month, day, hour, minute, second)
call print_date(model_time,'am2_to_dart: AM2 model date')
call print_date(model_time,'am2_to_dart: AM2 model date',logfileunit)

! write out state vector in "proprietary" format
iunit = open_restart_write(am2_to_dart_output_file)
call awrite_state_restart(model_time, statevector, iunit)
call close_restart(iunit)

call finalize_utilities('am2_to_dart')

end program am2_to_dart

! <next few lines under version control, do not edit>
! $URL$
! $Id$
! $Revision$
! $Date$
