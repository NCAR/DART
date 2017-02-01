! DART software - Copyright UCAR. This open source software is provided
! by UCAR, "as is", without charge, subject to all terms of use at
! http://www.image.ucar.edu/DAReS/DART/DART_download
!
! $Id$

program dart_to_model

!----------------------------------------------------------------------
! purpose: interface between DART and the model model
!
! method: Read DART state vector and overwrite values in a model analysis file.
!         If the DART state vector has an 'advance_to_time' present, a
!         file called model_in.DART is created with a time_manager_nml namelist 
!         appropriate to advance model to the requested time.
!
!         The dart_to_model_nml namelist setting for advance_time_present 
!         determines whether or not the input file has an 'advance_to_time'.
!         Typically, only temporary files like 'assim_model_state_ic' have
!         an 'advance_to_time'.
!
! author: Tim Hoar 25 Jun 09, revised 12 July 2010
!----------------------------------------------------------------------

use        types_mod, only : r8
use    utilities_mod, only : initialize_utilities, finalize_utilities, &
                             find_namelist_in_file, check_namelist_read, &
                             logfileunit, open_file, close_file
use  assim_model_mod, only : open_restart_read, aread_state_restart, close_restart
use time_manager_mod, only : time_type, print_time, print_date, operator(-), &
                             get_time, get_date
use        model_mod, only : static_init_model, statevector_to_analysis_file, &
                             get_model_size, get_model_analysis_filename,     &
                             write_model_time, print_variable_ranges

implicit none

! version controlled file description for error handling, do not edit
character(len=256), parameter :: source   = &
   "$URL$"
character(len=32 ), parameter :: revision = "$Revision$"
character(len=128), parameter :: revdate  = "$Date$"

!------------------------------------------------------------------
! The namelist variables
!------------------------------------------------------------------

character(len=256)  :: dart_to_model_input_file = 'dart_restart'
logical             :: advance_time_present     = .false.
character(len=256)  :: time_filename            = 'mpas_time'
logical             :: print_data_ranges        = .true.

namelist /dart_to_model_nml/ dart_to_model_input_file, &
                             advance_time_present,     &
                             time_filename,            &
                             print_data_ranges

!----------------------------------------------------------------------

integer               :: iunit, io, x_size
type(time_type)       :: model_time, adv_to_time
real(r8), allocatable :: statevector(:)
character(len=256)    :: model_analysis_filename

!----------------------------------------------------------------------

call initialize_utilities(progname='dart_to_model')

! Read the namelist to get the input filename. 

call find_namelist_in_file("input.nml", "dart_to_model_nml", iunit)
read(iunit, nml = dart_to_model_nml, iostat = io)
call check_namelist_read(iunit, io, "dart_to_model_nml")


!----------------------------------------------------------------------
! Call model_mod:static_init_model() which reads the model namelists
! to set grid sizes, etc.
!----------------------------------------------------------------------

call static_init_model()

call get_model_analysis_filename(model_analysis_filename)

x_size = get_model_size()
allocate(statevector(x_size))

write(*,*)
write(*,*) 'dart_to_model: converting DART file ', "'"//trim(dart_to_model_input_file)//"'"
write(*,*) 'to model analysis file ', "'"//trim(model_analysis_filename)//"'" 

!----------------------------------------------------------------------
! Reads the valid time, the state, and the target time.
!----------------------------------------------------------------------

iunit = open_restart_read(dart_to_model_input_file)

if ( advance_time_present ) then
   call aread_state_restart(model_time, statevector, iunit, adv_to_time)
else
   call aread_state_restart(model_time, statevector, iunit)
endif
call close_restart(iunit)

!----------------------------------------------------------------------
! if requested, print out the data ranges variable by variable
! (note if we are clamping data values, that happens in the
! conversion routine and these values are before the clamping happens.)
!----------------------------------------------------------------------
if (print_data_ranges) then
    call print_variable_ranges(statevector)
endif


!----------------------------------------------------------------------
! update the current model state vector
! Convey the amount of time to integrate the model ...
! time_manager_nml: stop_option, stop_count increments
!----------------------------------------------------------------------

call statevector_to_analysis_file(statevector, model_analysis_filename, model_time)

! write time into in text format (YYYY-MM-DD_hh:mm:ss) into a file.
! if advance time is there, write the current time then advance time.
! otherwise just write current time.
if ( advance_time_present ) then
   call write_model_time(time_filename, model_time, adv_to_time)
else
   call write_model_time(time_filename, model_time)
endif

!----------------------------------------------------------------------
! Log what we think we're doing, and exit.
!----------------------------------------------------------------------

call print_date( model_time,'dart_to_model:model date')
call print_time( model_time,'dart_to_model:model time')
call print_date( model_time,'dart_to_model:model date',logfileunit)
call print_time( model_time,'dart_to_model:model time',logfileunit)

if ( advance_time_present ) then
call print_time(adv_to_time,'dart_to_model:advance_to time')
call print_date(adv_to_time,'dart_to_model:advance_to date')
call print_time(adv_to_time,'dart_to_model:advance_to time',logfileunit)
call print_date(adv_to_time,'dart_to_model:advance_to date',logfileunit)
endif

call finalize_utilities()

end program dart_to_model

! <next few lines under version control, do not edit>
! $URL$
! $Id$
! $Revision$
! $Date$
