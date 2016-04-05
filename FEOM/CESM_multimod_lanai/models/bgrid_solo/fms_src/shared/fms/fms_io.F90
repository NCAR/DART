!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!                                                                   !!
!!                   GNU General Public License                      !!
!!                                                                   !!
!! This file is part of the Flexible Modeling System (FMS).          !!
!!                                                                   !!
!! FMS is free software; you can redistribute it and/or modify       !!
!! it and are expected to follow the terms of the GNU General Public !!
!! License as published by the Free Software Foundation.             !!
!!                                                                   !!
!! FMS is distributed in the hope that it will be useful,            !!
!! but WITHOUT ANY WARRANTY; without even the implied warranty of    !!
!! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the     !!
!! GNU General Public License for more details.                      !!
!!                                                                   !!
!! You should have received a copy of the GNU General Public License !!
!! along with FMS; if not, write to:                                 !!
!!          Free Software Foundation, Inc.                           !!
!!          59 Temple Place, Suite 330                               !!
!!          Boston, MA  02111-1307  USA                              !!
!! or see:                                                           !!
!!          http://www.gnu.org/licenses/gpl.txt                      !!
!!                                                                   !!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
module fms_io_mod

use utilities_mod, only : open_file, get_unit
  implicit none
  private

  character(len=32), private :: threading_read, fileset_read, threading_write,&
       fileset_write, format ! global i/o settings

  integer, private :: thread_r, thread_w, fset_r, fset_w, form

  logical, private :: module_is_initialized = .FALSE.
  logical, private :: read_all_pe = .TRUE.
  character(len=64) :: iospec_ieee32 = '-N ieee_32'
  
  public :: fms_io_init, fms_io_exit
  public :: open_namelist_file, open_restart_file, open_ieee32_file, close_file 

!  ---- version number -----

  character(len=128) :: version = '$Revision$'
  character(len=128) :: tagname = '$Id$'

contains

  subroutine fms_io_init()
    ! assign default values for default_file
    IMPLICIT NONE

    integer  :: i,j, unit, io_status
    logical :: file_exist

    namelist /fms_io_nml/ threading_read, fileset_read, threading_write,&
         fileset_write, format, read_all_pe, iospec_ieee32


  end subroutine fms_io_init






  subroutine fms_io_exit()
    ! <OVERVIEW>
    ! This routine is called after all fields have been written to temporary files
    ! The netcdf files are created here
    ! </OVERVIEW>
    IMPLICIT NONE


  end subroutine fms_io_exit


  ! ...




!#######################################################################
!
! routines for opening specific types of files:
!
!                       form        action 
! open_namelist_file  MPP_ASCII   MPP_RDONLY  
! open restart_file   MPP_NATIVE
! open_ieee32_file    MPP_IEEE32
!
! all have: access=MPP_SEQUENTIAL, threading=MPP_SINGLE, nohdrs=.true.
! use the close_file interface to close these files
!
! if other types of files need to be opened the mpp_open and
! mpp_close routines in the mpp_io_mod should be used
!
!#######################################################################
! opens single namelist file for reading only by all PEs
! the default file opened is called "input.nml"

  function open_namelist_file (file) result (unit)
    character(len=*), intent(in), optional :: file
    integer :: unit

    if (.not.module_is_initialized) call fms_io_init ( )

    if (present(file)) then
       unit = open_file(file, form = 'FORMATTED', action = 'READ')
       !call mpp_open ( unit, file, form=MPP_ASCII, action=MPP_RDONLY, &
       !     access=MPP_SEQUENTIAL, threading=MPP_SINGLE )
    else
       unit = open_file('input.nml', form = 'FORMATTED', action = 'READ')
       !call mpp_open ( unit, 'input.nml', form=MPP_ASCII, action=MPP_RDONLY, &
       !     access=MPP_SEQUENTIAL, threading=MPP_SINGLE )
    endif

  end function open_namelist_file

  !#######################################################################
  ! opens single restart file for reading by all PEs or
  ! writing by root PE only
  ! the file has native format and no mpp header records

  function open_restart_file (file, action) result (unit)
    character(len=*), intent(in) :: file, action
    integer :: unit

    integer :: mpp_action

    if (.not.module_is_initialized) call fms_io_init ( )

    !   --- action (read,write) ---

    unit = open_file(file, form = 'UNFORMATTED', action = action)

    !call mpp_open ( unit, file, form=MPP_NATIVE, action=mpp_action, &
    !     access=MPP_SEQUENTIAL, threading=MPP_SINGLE, nohdrs=.true. )

  end function open_restart_file

  !#######################################################################
  ! opens single 32-bit ieee file for reading by all PEs or 
  ! writing by root PE only (writing is not recommended)
  ! the file has no mpp header records

  function open_ieee32_file (file, action) result (unit)
    character(len=*), intent(in) :: file, action
    integer :: unit

    integer :: mpp_action

    if (.not.module_is_initialized) call fms_io_init ( )

    !   --- action (read,write) ---

    !select case (lowercase(trim(action)))
    !case ('read')
    !   mpp_action = MPP_RDONLY
    !case ('write')
    !   mpp_action = MPP_OVERWR
    !case default
    !   call mpp_error (FATAL,'invalid option for argument action')
    !end select

! Open the file
   unit = get_unit()

    if (iospec_ieee32(1:1) == ' ') then
       open(unit=unit, file=file, form = 'UNFORMATTED', action = action)
       !call mpp_open ( unit, file, form=MPP_IEEE32, action=mpp_action, &
       !     access=MPP_SEQUENTIAL, threading=MPP_SINGLE,    &
       !     nohdrs=.true. )
    else
       open(unit=unit, file=file, form = 'UNFORMATTED', action = action)
       !call mpp_open ( unit, file, form=MPP_IEEE32, action=mpp_action, &
       !     access=MPP_SEQUENTIAL, threading=MPP_SINGLE,    &
       !     nohdrs=.true., iospec=iospec_ieee32 )
    endif

  end function open_ieee32_file

  !#######################################################################

  subroutine close_file (unit, status)
    integer,          intent(in)           :: unit
    character(len=*), intent(in), optional :: status

    if (.not.module_is_initialized) call fms_io_init ( )

!    if (unit == stdlog()) return

    close(unit)

    !if (present(status)) then
    !   if (lowercase(trim(status)) == 'delete') then
    !      call mpp_close (unit, action=MPP_DELETE)
    !   else
    !      call mpp_error(FATAL,'invalid value for status')
    !   endif
    !else
    !   call mpp_close (unit)
    !endif

  end subroutine close_file

  !#######################################################################

end module fms_io_mod


