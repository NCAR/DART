! DART software - Copyright UCAR. This open source software is provided
! by UCAR, "as is", without charge, subject to all terms of use at
! http://www.image.ucar.edu/DAReS/DART/DART_download
!
! $Id$

program aether_to_dart

!----------------------------------------------------------------------
! purpose: interface between the GITM model and DART
!
! method: Read aether "restart" files of model state (multiple files, 
!         one block per aether mpi task)
!         Reform fields into a DART netcdf file
! TODO: Should this be an MPI program so that all members can be done at once?
!       Get the ensemble size from input.nml:filter_nml.
!       Can I send each member to a different node, so that the restart files
!       could all be read at once on separate processors, and still be local 
!       to the member's filter_input.nc?
!
! USAGE:  The aether restart dirname and output filename are read from 
!         the aether_to_dart_nml namelist.
!         <edit input.nml:aether_to_dart_nml>
!----------------------------------------------------------------------

use        types_mod, only : r8

use    utilities_mod, only : initialize_utilities, finalize_utilities, &
                             error_handler, E_MSG

use        model_mod, only : restart_files_to_netcdf

implicit none

! version controlled file description for error handling, do not edit
character(len=256), parameter :: source   = &
   "$URL$"
character(len=32 ), parameter :: revision = "$Revision$"
character(len=128), parameter :: revdate  = "$Date$"

character(len=*), parameter :: program_name = 'aether_to_dart'

! !-----------------------------------------------------------------------
! ! namelist parameters with default values.
! !-----------------------------------------------------------------------
! 
! character(len=256) :: aether_restart_input_dirname    = 'none'
! ! TODO: the calling script will need to move this to a name with $member in it,
! !       or use filter_nml:input_state_file_list
! ! TODO: Create the filter filename from filter_root, as in dart_to_aether.
! character(len=256) :: aether_to_dart_output_file    = 'filter_input.nc'
! 
! namelist /aether_to_dart_nml/ aether_restart_input_dirname,        &
!                               aether_to_dart_output_file, variables
! 
!----------------------------------------------------------------------
! global storage
!----------------------------------------------------------------------

integer :: member

!======================================================================

call initialize_utilities(program_name)

!----------------------------------------------------------------------
! Get the ensemble member
! TODO: The script must echo the member number to the aether_to_dart.
!----------------------------------------------------------------------
member = -88
read '(I3)', member
print*,'aether_to_dart: member = ',member

!----------------------------------------------------------------------
! Convert the files
!----------------------------------------------------------------------

! call error_handler(E_MSG, '', '')
! write(string1,*) 'converting aether restart files in directory ', &
!                  "'"//trim(aether_restart_input_dirname)//"'"
! write(string2,*) ' to the NetCDF file ', "'"//trim(aether_to_dart_output_file)//"'"
! call error_handler(E_MSG, program_name, string1, text2=string2)
! call error_handler(E_MSG, '', '')

call restart_files_to_netcdf(member)

! call error_handler(E_MSG, '', '')
! write(string1,*) 'Successfully converted the GITM restart files to ', &
!                  "'"//trim(aether_to_dart_output_file)//"'"
! call error_handler(E_MSG, program_name, string1)
! call error_handler(E_MSG, '', '')

!----------------------------------------------------------------------
! Finish up
!----------------------------------------------------------------------

! end - close the log, etc
call finalize_utilities()

end program aether_to_dart

! <next few lines under version control, do not edit>
! $URL$
! $Id$
! $Revision$
! $Date$
