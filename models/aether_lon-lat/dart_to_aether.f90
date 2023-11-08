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

!-----------------------------------------------------------------------
! namelist parameters with default values.
!-----------------------------------------------------------------------

character (len = 256) :: aether_restart_dirname  = 'none'
character (len = 64) :: filter_io_root = 'filter_output'
character (len = 64) :: filter_io_name 

namelist /dart_to_aether_nml/   &
     aether_restart_dirname,    &
     filter_io_root

!----------------------------------------------------------------------
! global storage
!----------------------------------------------------------------------

integer               :: iunit, io, member
character(len=512)    :: string1, string2

!----------------------------------------------------------------------
! Get the ensemble member
! TODO: The script must echo the member number to the dart_to_aether.
!             There may be a mismatch between member numbers in DART and Aether; F or C indexing.
!----------------------------------------------------------------------
member = -88
read '(I3)', member
print*,'dart_to_aether: member = ',member

write(filter_io_name,'(2A,I0.4,A3)') trim(filter_io_root),'_',member,'.nc'

!======================================================================

call initialize_utilities(progname=progname)

!----------------------------------------------------------------------
! Read the namelist
!----------------------------------------------------------------------

call find_namelist_in_file("input.nml", "dart_to_aether_nml", iunit)
read(iunit, nml = dart_to_aether_nml, iostat = io)
call check_namelist_read(iunit, io, "dart_to_aether_nml") ! closes, too.

print*,'After namelist; aether_restart_dirname = ',aether_restart_dirname

call error_handler(E_MSG,progname,'','',revision,revdate)
write(string1,*) 'Extracting fields from DART file ', "'"//trim(filter_io_name)//"'"
write(string2,*) 'into Aether restart files in directory ', "'"//trim(aether_restart_dirname)//"'"
call error_handler(E_MSG,progname,string1,source,revision,revdate,text2=string2)

!----------------------------------------------------------------------
! Reads the valid time, the state, and the target time.
!----------------------------------------------------------------------

! TODO: netcdf_to_restart_files; need all these file and dir names?
call netcdf_to_restart_files(filter_io_name, member, aether_restart_dirname)

!----------------------------------------------------------------------
! Log what we think we're doing, and exit.
!----------------------------------------------------------------------
call error_handler(E_MSG,progname,'','',revision,revdate)
write(string1,*) 'Successfully converted to the Aether restart files in directory'
write(string2,*) "'"//trim(aether_restart_dirname)//"'"
call error_handler(E_MSG,progname,string1,source,revision,revdate,text2=string2)

! end - close the log, etc
call finalize_utilities()

end program dart_to_aether

! <next few lines under version control, do not edit>
! $URL$
! $Id$
! $Revision$
! $Date$
