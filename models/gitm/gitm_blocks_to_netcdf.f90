! DART software - Copyright UCAR. This open source software is provided
! by UCAR, "as is", without charge, subject to all terms of use at
! http://www.image.ucar.edu/DAReS/DART/DART_download
!
! $Id$

program gitm_blocks_to_netcdf

!----------------------------------------------------------------------
! purpose: interface between the GITM model and DART
!
! method: Read gitm "restart" files of model state (multiple files, one
!         block per gitm mpi task)
!         Reform fields into a DART netcdf file
!
! USAGE:  The gitm dirname is read from the gitm_in namelist
!         <edit gitm_blocks_to_netcdf_output_file in input.nml:gitm_blocks_to_netcdf_nml>
!         gitm_blocks_to_netcdf
!----------------------------------------------------------------------

use        types_mod, only : r8

use    utilities_mod, only : initialize_utilities, finalize_utilities, &
                             find_namelist_in_file, check_namelist_read

use        model_mod, only : static_init_model, get_model_size, &
                             get_gitm_restart_dirname !! , restart_file_to_statevector

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

character(len=128) :: gitm_blocks_to_netcdf_output_file  = 'filter_input.nc'

namelist /gitm_blocks_to_netcdf_nml/ gitm_blocks_to_netcdf_output_file

!----------------------------------------------------------------------
! global storage
!----------------------------------------------------------------------

integer               :: iunit, io, x_size
type(time_type)       :: model_time
character(len=256)    :: gitm_restart_dirname  = 'none'

!======================================================================

call initialize_utilities(progname='gitm_blocks_to_netcdf')

!----------------------------------------------------------------------
! Read the namelist.
!----------------------------------------------------------------------

call find_namelist_in_file("input.nml", "gitm_blocks_to_netcdf_nml", iunit)
read(iunit, nml = gitm_blocks_to_netcdf_nml, iostat = io)
call check_namelist_read(iunit, io, "gitm_blocks_to_netcdf_nml") ! closes, too.

!----------------------------------------------------------------------
! Call model_mod:static_init_model() which reads the gitm namelists
! to set grid sizes, etc.
!----------------------------------------------------------------------

call static_init_model()
call get_gitm_restart_dirname(gitm_restart_dirname)

write(*,*)
write(*,*) 'gitm_blocks_to_netcdf: converting gitm restart files in directory ', &
           "'"//trim(gitm_restart_dirname)//"'"
write(*,*) ' to DART file ', "'"//trim(gitm_blocks_to_netcdf_output_file)//"'"

x_size = get_model_size()
!allocate(statevector(x_size))

!call restart_file_to_statevector(gitm_restart_dirname, statevector, model_time)

!iunit = open_restart_write(gitm_blocks_to_netcdf_output_file)
!call awrite_state_restart(model_time, statevector, iunit)
!call close_restart(iunit)

!----------------------------------------------------------------------
! finish up
!----------------------------------------------------------------------

call print_date(model_time, str='gitm_blocks_to_netcdf:gitm model date')
call print_time(model_time, str='gitm_blocks_to_netcdf:DART model time')

! end - close the log, etc
call finalize_utilities()

end program gitm_blocks_to_netcdf

! <next few lines under version control, do not edit>
! $URL$
! $Id$
! $Revision$
! $Date$
