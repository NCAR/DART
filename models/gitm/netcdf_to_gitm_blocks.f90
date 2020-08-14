! DART software - Copyright UCAR. This open source software is provided
! by UCAR, "as is", without charge, subject to all terms of use at
! http://www.image.ucar.edu/DAReS/DART/DART_download
!
! $Id$

program netcdf_to_gitm_blocks

!----------------------------------------------------------------------
! purpose: interface between DART and the GITM model
!
! method: Read DART state netcdf files and overwrite values in a gitm restart file.
!
! NOT FINISHED!  HARDLY STARTED.
! this version assumes that the grid is global and the data needs to be
! blocked into one block per gitm mpi task.  there is a different converter
! for when gitm only needs a single input/output file.
!
!----------------------------------------------------------------------

use        types_mod, only : r8

use    utilities_mod, only : initialize_utilities, finalize_utilities, &
                             find_namelist_in_file, check_namelist_read, &
                             logfileunit, open_file, close_file

use        model_mod, only : static_init_model, get_model_size, &
                             get_gitm_restart_dirname !! , statevector_to_restart_file

use time_manager_mod, only : time_type, print_time, print_date, operator(-), &
                             get_time, get_date

implicit none

! version controlled file description for error handling, do not edit
character(len=256), parameter :: source   = &
   "$URL$"
character(len=32 ), parameter :: revision = "$Revision$"
character(len=128), parameter :: revdate  = "$Date$"

!-----------------------------------------------------------------------
! namelist parameters with default values.
!-----------------------------------------------------------------------

character (len = 256) :: netcdf_to_gitm_blocks_input_file = 'filter_restart.nc'

namelist /netcdf_to_gitm_blocks_nml/ netcdf_to_gitm_blocks_input_file

!----------------------------------------------------------------------
! global storage
!----------------------------------------------------------------------

integer               :: iunit, io, x_size
type(time_type)       :: model_time, adv_to_time
character(len=256)    :: gitm_restart_dirname  = 'none'

!======================================================================

call initialize_utilities(progname='netcdf_to_gitm_blocks')

!----------------------------------------------------------------------
! Read the namelist.
!----------------------------------------------------------------------

call find_namelist_in_file("input.nml", "netcdf_to_gitm_blocks_nml", iunit)
read(iunit, nml = netcdf_to_gitm_blocks_nml, iostat = io)
call check_namelist_read(iunit, io, "netcdf_to_gitm_blocks_nml")

!----------------------------------------------------------------------
! Call model_mod:static_init_model() which reads the gitm namelists
! to set grid sizes, etc.
!----------------------------------------------------------------------

call static_init_model()
call get_gitm_restart_dirname(gitm_restart_dirname)

write(*,*)
write(*,*) 'netcdf_to_gitm_blocks: converting DART file ', "'"//trim(netcdf_to_gitm_blocks_input_file)//"'"
write(*,*) 'to gitm restart files in directory ', "'"//trim(gitm_restart_dirname)//"'"

x_size = get_model_size()
!allocate(statevector(x_size))

!----------------------------------------------------------------------
! Reads the valid time, the state, and the target time.
!----------------------------------------------------------------------

!iunit = open_restart_read(netcdf_to_gitm_blocks_input_file)

!if ( advance_time_present ) then
!   call aread_state_restart(model_time, statevector, iunit, adv_to_time)
!else
!   call aread_state_restart(model_time, statevector, iunit)
!endif
!call close_restart(iunit)

!----------------------------------------------------------------------
! update the current gitm state vector
! Convey the amount of time to integrate the model ...
! time_manager_nml: stop_option, stop_count increments
!----------------------------------------------------------------------

!call statevector_to_restart_file(statevector, gitm_restart_dirname, model_time)

!if ( advance_time_present ) then
!   call write_gitm_time_control(model_time, adv_to_time)
!endif

!----------------------------------------------------------------------
! Log what we think we're doing, and exit.
!----------------------------------------------------------------------

call print_date( model_time,'netcdf_to_gitm_blocks:gitm model date')
call print_time( model_time,'netcdf_to_gitm_blocks:DART model time')
call print_date( model_time,'netcdf_to_gitm_blocks:gitm model date',logfileunit)
call print_time( model_time,'netcdf_to_gitm_blocks:DART model time',logfileunit)

! end - close the log, etc
call finalize_utilities()

!======================================================================
contains
!======================================================================

subroutine write_gitm_time_control(model_time, adv_to_time)
! The idea is to write a text file with the following structure:
!
!#TIMESTART
!2003            year
!06              month
!21              day
!00              hour
!00              minute
!00              second
!
!#TIMEEND
!2003            year
!07              month
!21              day
!00              hour
!00              minute
!00              second
!

type(time_type), intent(in) :: model_time, adv_to_time
integer :: iyear,imonth,iday,ihour,imin,isec

iunit = open_file('DART_GITM_time_control.txt', action='write')
write(iunit,*)

! the end time comes first.

call get_date(adv_to_time,iyear,imonth,iday,ihour,imin,isec)
write(iunit,'(''#TIMEEND'')')
write(iunit,'(i4.4,10x,''year''  )')iyear
write(iunit,'(i2.2,12x,''month'' )')imonth
write(iunit,'(i2.2,12x,''day''   )')iday
write(iunit,'(i2.2,12x,''hour''  )')ihour
write(iunit,'(i2.2,12x,''minute'')')imin
write(iunit,'(i2.2,12x,''second'')')isec
write(iunit,*)

call get_date(model_time,iyear,imonth,iday,ihour,imin,isec)
write(iunit,'(''#TIMESTART'')')
write(iunit,'(i4.4,10x,''year''  )')iyear
write(iunit,'(i2.2,12x,''month'' )')imonth
write(iunit,'(i2.2,12x,''day''   )')iday
write(iunit,'(i2.2,12x,''hour''  )')ihour
write(iunit,'(i2.2,12x,''minute'')')imin
write(iunit,'(i2.2,12x,''second'')')isec
write(iunit,*)

call close_file(iunit)
end subroutine write_gitm_time_control



end program netcdf_to_gitm_blocks

! <next few lines under version control, do not edit>
! $URL$
! $Id$
! $Revision$
! $Date$
