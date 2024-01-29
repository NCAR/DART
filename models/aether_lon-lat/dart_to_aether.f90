! DART software - Copyright UCAR. This open source software is provided
! by UCAR, "as is", without charge, subject to all terms of use at
! http://www.image.ucar.edu/DAReS/DART/DART_download
!
! $Id$

program dart_to_aether

!----------------------------------------------------------------------
! purpose: interface between DART and the Aether model
!
! method: Read DART state netcdf files and overwrite values in Aether restart files.
!
! this version assumes that the DART grid is global and the data needs to be
! blocked into one block per Aether mpi task.  there is a different converter
! for when Aether only needs a single input/output file.
!
!----------------------------------------------------------------------

use    utilities_mod, only : initialize_utilities, finalize_utilities,   &
                             find_namelist_in_file, check_namelist_read, &
                             E_MSG, error_handler

use        model_mod, only : netcdf_to_restart_files

use time_manager_mod, only : operator(-)

implicit none

! version controlled file description for error handling, do not edit
character(len=256), parameter :: source   = "$URL$"
character(len=32 ), parameter :: revision = "$Revision$"
character(len=128), parameter :: revdate  = "$Date$"

character(len=*),   parameter :: progname = 'dart_to_aether'

!----------------------------------------------------------------------
! global storage
!----------------------------------------------------------------------

integer               :: member

!----------------------------------------------------------------------
! Get the ensemble member
! TODO: The script must echo the member number to the dart_to_aether.
!             There may be a mismatch between member numbers in DART and Aether; F or C indexing.
!----------------------------------------------------------------------
member = -88
read '(I3)', member
print*,'dart_to_aether: member = ',member

!======================================================================

call initialize_utilities(progname=progname)

!----------------------------------------------------------------------
! Reads the valid time, the state, and the target time.
!----------------------------------------------------------------------

! TODO: netcdf_to_restart_files; need all these file and dir names?
call netcdf_to_restart_files(member)

! end - close the log, etc
call finalize_utilities()

end program dart_to_aether

! <next few lines under version control, do not edit>
! $URL$
! $Id$
! $Revision$
! $Date$
