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
!         <edit gitm_to_netcdf_output_file in input.nml:gitm_blocks_to_netcdf_nml>
!         gitm_blocks_to_netcdf
!----------------------------------------------------------------------

use        types_mod, only : r8

use    utilities_mod, only : initialize_utilities, finalize_utilities, &
                             find_namelist_in_file, check_namelist_read, &
                             error_handler, E_MSG

use        model_mod, only : restart_files_to_netcdf

use time_manager_mod, only : time_type, print_time, print_date

implicit none

! version controlled file description for error handling, do not edit
character(len=256), parameter :: source   = &
   "$URL$"
character(len=32 ), parameter :: revision = "$Revision$"
character(len=128), parameter :: revdate  = "$Date$"

character(len=512) :: string1, string2
character(len=*), parameter :: program_name = 'gitm_blocks_to_netcdf'

!-----------------------------------------------------------------------
! namelist parameters with default values.
!-----------------------------------------------------------------------

character(len=256) :: gitm_restart_input_dirname    = 'none'
character(len=256) :: gitm_to_netcdf_output_file    = 'filter_input.nc'

namelist /gitm_blocks_to_netcdf_nml/ gitm_restart_input_dirname,        &
                                     gitm_to_netcdf_output_file

!----------------------------------------------------------------------
! global storage
!----------------------------------------------------------------------

integer               :: iunit, io

!======================================================================

call initialize_utilities(program_name)

!----------------------------------------------------------------------
! Read the namelist
!----------------------------------------------------------------------

call find_namelist_in_file("input.nml", "gitm_blocks_to_netcdf_nml", iunit)
read(iunit, nml = gitm_blocks_to_netcdf_nml, iostat = io)
call check_namelist_read(iunit, io, "gitm_blocks_to_netcdf_nml") ! closes, too.

!----------------------------------------------------------------------
! Convert the files
!----------------------------------------------------------------------

call error_handler(E_MSG, '', '')
write(string1,*) 'converting gitm restart files in directory ', &
                 "'"//trim(gitm_restart_input_dirname)//"'"
write(string2,*) ' to the NetCDF file ', "'"//trim(gitm_to_netcdf_output_file)//"'"
call error_handler(E_MSG, program_name, string1, text2=string2)
call error_handler(E_MSG, '', '')

call restart_files_to_netcdf(gitm_restart_input_dirname, gitm_to_netcdf_output_file)

call error_handler(E_MSG, '', '')
write(string1,*) 'Successfully converted the GITM restart files to ', &
                 "'"//trim(gitm_to_netcdf_output_file)//"'"
call error_handler(E_MSG, program_name, string1)
call error_handler(E_MSG, '', '')

!----------------------------------------------------------------------
! Finish up
!----------------------------------------------------------------------

! end - close the log, etc
call finalize_utilities()

end program gitm_blocks_to_netcdf

! <next few lines under version control, do not edit>
! $URL$
! $Id$
! $Revision$
! $Date$
