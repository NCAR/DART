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
! this version assumes that the grid is global and the data needs to be
! blocked into one block per gitm mpi task.  there is a different converter
! for when gitm only needs a single input/output file.
!
!----------------------------------------------------------------------

use        types_mod, only : r8

use    utilities_mod, only : initialize_utilities, finalize_utilities,   &
                             find_namelist_in_file, check_namelist_read, &
                             open_file, close_file, E_MSG, error_handler

use        model_mod, only : netcdf_to_restart_files

use time_manager_mod, only : time_type, print_time, print_date, operator(-), &
                             get_time, get_date

implicit none

! version controlled file description for error handling, do not edit
character(len=256), parameter :: source   = &
   "$URL$"
character(len=32 ), parameter :: revision = "$Revision$"
character(len=128), parameter :: revdate  = "$Date$"

character(len=*),   parameter :: progname = 'netcdf_to_gitm_blocks'

!-----------------------------------------------------------------------
! namelist parameters with default values.
!-----------------------------------------------------------------------

character (len = 256) :: gitm_restart_input_dirname  = 'none'
character (len = 256) :: gitm_restart_output_dirname = 'none'
character (len = 256) :: netcdf_to_gitm_blocks_input_file = 'filter_restart.nc'

namelist /netcdf_to_gitm_blocks_nml/   &
     gitm_restart_input_dirname,       &
     gitm_restart_output_dirname,      &
     netcdf_to_gitm_blocks_input_file         

!----------------------------------------------------------------------
! global storage
!----------------------------------------------------------------------

integer               :: iunit, io
character(len=512)    :: string1, string2, string3

!======================================================================

call initialize_utilities(progname=progname)

!----------------------------------------------------------------------
! Read the namelist.
!----------------------------------------------------------------------

call find_namelist_in_file("input.nml", "netcdf_to_gitm_blocks_nml", iunit)
read(iunit, nml = netcdf_to_gitm_blocks_nml, iostat = io)
call check_namelist_read(iunit, io, "netcdf_to_gitm_blocks_nml")

call error_handler(E_MSG,progname,'','',revision,revdate)
write(string1,*) 'converting DART file ', "'"//trim(netcdf_to_gitm_blocks_input_file)//"'"
write(string2,*) 'to gitm restart files in directory ', "'"//trim(gitm_restart_output_dirname)//"'"
write(string3,*) 'using the restart files in directory ', "'"//trim(gitm_restart_input_dirname)//"' as a template"
call error_handler(E_MSG,progname,string1,source,revision,revdate,text2=string2,text3=string3)

!----------------------------------------------------------------------
! Reads the valid time, the state, and the target time.
!----------------------------------------------------------------------

call netcdf_to_restart_files(netcdf_to_gitm_blocks_input_file,gitm_restart_output_dirname,&
                             gitm_restart_input_dirname)

!----------------------------------------------------------------------
! Log what we think we're doing, and exit.
!----------------------------------------------------------------------
call error_handler(E_MSG,progname,'','',revision,revdate)
call error_handler(E_MSG,progname,'','',revision,revdate)
write(string1,*) 'Successfully converted to the gitm restart files in directory'
write(string2,*) "'"//trim(gitm_restart_output_dirname)//"'"
call error_handler(E_MSG,progname,string1,source,revision,revdate,text2=string2)

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
