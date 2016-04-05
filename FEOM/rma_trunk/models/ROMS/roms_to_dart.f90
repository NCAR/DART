! DART software - Copyright 2004 - 2013 UCAR. This open source software is
! provided by UCAR, "as is", without charge, subject to all terms of use at
! http://www.image.ucar.edu/DAReS/DART/DART_download
!
! $Id$

!-----------------------------------------------------------------------
!> Converts a subset of variables from a ROMS restart file to a DART file.
!>
!> Reads a ROMS "restart/analysis" file and
!> reforms them into a DART state vector (control vector).
!> Write out state vector in "proprietary" format for DART.
!> The output is a "DART restart format" file.
!> 
!> The ROMS filename is read from the input.nmnl:model_mod namelist.
!> The (output) DART filename is read from the 
!> input.nml:roms_to_dart_nml namelist.
!>
!> author: PENG XIU 12/2013 @ University of Maine
!>         peng.xiu@maine.edu
!>
!> subsequently modified by TJH 1/2015

program roms_to_dart

use        types_mod, only : r8
use    utilities_mod, only : initialize_utilities, finalize_utilities, &
                             find_namelist_in_file, check_namelist_read, &
                             logfileunit, error_handler, E_MSG
use        model_mod, only : restart_file_to_sv, static_init_model, &
                             get_model_size, get_model_restart_filename, &
                             get_time_from_namelist, pert_model_state,   &
                             print_variable_ranges
use  assim_model_mod, only : awrite_state_restart, open_restart_write, close_restart
use time_manager_mod, only : time_type, print_time, print_date

implicit none

! version controlled file description for error handling, do not edit
character(len=256), parameter :: source   = &
   "$URL$"
character(len=32 ), parameter :: revision = "$Revision$"
character(len=128), parameter :: revdate  = "$Date$"

!-----------------------------------------------------------------------
! namelist parameters with default values.
!-----------------------------------------------------------------------

character(len=256) :: roms_to_dart_output_file = 'dart.ics'
logical            :: perturb_state            = .false.
logical            :: print_data_ranges        = .true.
namelist /roms_to_dart_nml/ roms_to_dart_output_file, &
                            perturb_state, print_data_ranges

!----------------------------------------------------------------------
! global storage
!----------------------------------------------------------------------

integer               :: io, iunit, x_size
type(time_type)       :: model_time, last_file_time
real(r8), allocatable :: statevector(:)
real(r8), allocatable :: pert_state(:)
character(len=256)    :: roms_restart_filename
logical               :: interf_provided
character(len=256)    :: string1, string2

call initialize_utilities(progname='roms_to_dart')

!----------------------------------------------------------------------
! Read the namelist to get the output filename
!----------------------------------------------------------------------

call find_namelist_in_file("input.nml", "roms_to_dart_nml", iunit)
read(iunit, nml = roms_to_dart_nml, iostat = io)
call check_namelist_read(iunit, io, "roms_to_dart_nml") ! closes, too.

!----------------------------------------------------------------------
! Call model_mod:static_init_model() which reads the model namelists
! to set grid sizes, etc.  The 'ocean_time' times in the ROMS file
! cannot (at present) be correctly decoded to give dates for which we
! have real observations. Either the test files I have are from AD 223,
! or I am doing it wrong. Since I cannot decode them sensibly now,
! I am setting the time from the namelist value.
!----------------------------------------------------------------------

call static_init_model()

x_size = get_model_size()
allocate(statevector(x_size))

call get_model_restart_filename( roms_restart_filename )
model_time = get_time_from_namelist()

write(string1,*) '..  converting ROMS file <'//trim(roms_restart_filename)//'>'
write(string2,*) 'to     DART     file <'//trim(roms_to_dart_output_file)//'>'
call error_handler(E_MSG,'roms_to_dart:',string1,text2=string2)

iunit = open_restart_write(roms_to_dart_output_file)

!----------------------------------------------------------------------
! Now that we know the names, get to work.
!----------------------------------------------------------------------

call restart_file_to_sv(roms_restart_filename, statevector, last_file_time) 
if (print_data_ranges) call print_variable_ranges(statevector)

if(perturb_state) then

   call error_handler(E_MSG,'roms_to_dart:','Perturbing state.')

   allocate(pert_state(x_size))
   call pert_model_state(statevector, pert_state, interf_provided)
   call awrite_state_restart(model_time, pert_state, iunit)
   if (print_data_ranges) call print_variable_ranges(pert_state)
   deallocate(pert_state)

else

   call awrite_state_restart(model_time, statevector, iunit)

endif

deallocate(statevector)

call close_restart(iunit)

!----------------------------------------------------------------------
! Log what we think we're doing, and exit.
!----------------------------------------------------------------------

call print_date(    model_time, str='roms_to_dart:ROMS  model date')
call print_date(    model_time, str='roms_to_dart:ROMS  model date',iunit=logfileunit)
call print_time(    model_time, str='roms_to_dart:DART  model time')
call print_time(    model_time, str='roms_to_dart:DART  model time',iunit=logfileunit)

call print_date(last_file_time, str='roms_to_dart:ROMS   file date')
call print_date(last_file_time, str='roms_to_dart:ROMS   file date',iunit=logfileunit)
call print_time(last_file_time, str='roms_to_dart:DART   file time')
call print_time(last_file_time, str='roms_to_dart:DART   file time',iunit=logfileunit)

call finalize_utilities()

end program roms_to_dart

! <next few lines under version control, do not edit>
! $URL$
! $Id$
! $Revision$
! $Date$

