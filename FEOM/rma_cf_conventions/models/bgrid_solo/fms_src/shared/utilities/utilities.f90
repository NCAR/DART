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

module utilities_mod

!-----------------------------------------------------------------------
!                 wrapper around module fms_mod
!
!  this module provides a transition from utilities_mod to fms_mod
!  users should use fms_mod as soon as possible
!
!-----------------------------------------------------------------------

use    fms_mod, only: file_exist, stdlog, open_ieee32_file,  &
                      read_data, write_data, set_domain,     &
                      check_nml_error, get_my_pe=>mpp_pe,    &
                      get_root_pe=>mpp_root_pe, get_num_pes=>mpp_npes, &
                      fms_write_version_number=>write_version_number,  &
                      error_mesg, utilities_init=>fms_init,            &
                      utilities_end=>fms_end, close_file,              &
                      NOTE, WARNING, FATAL, lowercase, uppercase,      &
                      mpp_clock_begin, mpp_clock_end,                  &
                      fms_mpp_clock_init=>mpp_clock_init,              &
                      check_sum=>mpp_chksum, get_domain_decomp
use    mpp_mod, only: sync_all_pes=>mpp_sync, &
                      mpp_sum, mpp_min, mpp_max
use mpp_io_mod, only: mpp_open, MPP_RDONLY, MPP_OVERWR, MPP_APPEND,  &
                      MPP_ASCII, MPP_NATIVE, MPP_IEEE32, MPP_NETCDF, &
                      MPP_SEQUENTIAL, MPP_DIRECT, MPP_SINGLE, MPP_MULTI


! public routines in fez release of utilities

public :: file_exist, open_file, read_data, write_data,  &
          set_domain, get_domain_decomp,                 &
          check_nml_error, write_version_number,         &
          error_mesg, NOTE, WARNING, FATAL,              &
          get_my_pe, get_num_pes, sync_all_pes,          &
          get_root_pe, mpp_sum, mpp_min, mpp_max,        &
          utilities_init, utilities_end, close_file,     &
          set_system_clock, check_system_clock, check_sum, &
          mpp_clock_begin, mpp_clock_end, mpp_clock_init,  &
          lowercase, uppercase
public :: print_version_number ! to be removed

logical :: do_init = .true.
                   
 contains
!#######################################################################

 function open_file ( file, form, action, access, threading, recl ) &
             result ( unit )

 character(len=*), intent(in) :: file 
 character(len=*), intent(in), optional :: form, action, access, threading
 integer         , intent(in), optional :: recl 
 integer  :: unit 

 character(len=32) :: form_local, action_local, access_local, thread_local
 character(len=32) :: action_ieee32
 logical :: open, no_headers, do_ieee32
 integer :: mpp_format, mpp_action, mpp_access, mpp_thread
!-----------------------------------------------------------------------

   if ( do_init ) then
        call utilities_init ( )
        do_init = .false.
   endif

!   ---- return stdlog if this is the logfile ----

    if (trim(file) == 'logfile.out') then
       unit = stdlog()
       return
    endif

!   ---- is this file open and connected to a unit ?? ---- 

   inquire (file=trim(file), opened=open, number=unit)

!  cannot open a file that is already open
!  except for the log file

   if ( open .and. unit >= 0 ) then
      call error_mesg ('open_file in fms_mod', &
                       'file '//trim(file)//' is already open', FATAL)
   endif   

!  --- defaults ---

   form_local   = 'formatted';  if (present(form))      form_local   = form
   access_local = 'sequential'; if (present(access))    access_local = access
   thread_local = 'single';     if (present(threading)) thread_local = threading
   no_headers   = .true.
   do_ieee32    = .false.

   if (present(action)) then    ! must be present
      action_local = action
   else
      call error_mesg ('open_file in fms_mod', 'argument action not present', FATAL)
   endif


!   --- file format ---

    select case (lowercase(trim(form_local)))
       case ('formatted')
           mpp_format = MPP_ASCII
       case ('ascii')
           mpp_format = MPP_ASCII
       case ('unformatted')
           mpp_format = MPP_NATIVE
       case ('native')
           mpp_format = MPP_NATIVE
       case ('ieee32')
           do_ieee32 = .true.
       case ('netcdf')
           mpp_format = MPP_NETCDF
       case default
           call error_mesg ('open_file in fms_mod', &
                            'invalid option for argument form', FATAL)
    end select

!   --- action (read,write,append) ---

    select case (lowercase(trim(action_local)))
       case ('read')
           mpp_action = MPP_RDONLY
       case ('write')
           mpp_action = MPP_OVERWR
       case ('append')
           mpp_action = MPP_APPEND
       case default
           call error_mesg ('open_file in fms_mod', &
                            'invalid option for argument action', FATAL)
    end select

!   --- file access (sequential,direct) ---

    select case (lowercase(trim(access_local)))
       case ('sequential')
           mpp_access = MPP_SEQUENTIAL
       case ('direct')
           mpp_access = MPP_DIRECT
       case default
           call error_mesg ('open_file in fms_mod', &
                            'invalid option for argument access', FATAL)
    end select

!   --- threading (single,multi) ---

    select case (lowercase(trim(thread_local)))
       case ('single')
           mpp_thread = MPP_SINGLE
       case ('multi')
           mpp_thread = MPP_MULTI
       case default
           call error_mesg ('open_file in fms_mod', &
                            'invalid option for argument thread', FATAL)
           if (trim(file) /= '_read_error.nml') no_headers = .false.
    end select

!   ---------------- open file -----------------------

    if ( .not.do_ieee32 ) then
       call mpp_open ( unit, file, form=mpp_format, action=mpp_action, &
                       access=mpp_access, threading=mpp_thread,        &
                       nohdrs=no_headers, recl=recl )
    else
     ! special open for ieee32 file
     ! fms_mod has iospec value
     ! pass local action flag to open changing append to write
       action_ieee32 = action_local
       if (lowercase(trim(action_ieee32)) == 'append') action_ieee32 = 'write'
       unit = open_ieee32_file ( file, action_ieee32 )
    endif

!-----------------------------------------------------------------------

 end function open_file

!#######################################################################
! do nothing routines to mimic old interfaces
! use fms_mod instead
 
 subroutine set_system_clock
 end subroutine set_system_clock

 subroutine check_system_clock (string)
 character(len=*), intent(in), optional :: string
 end subroutine check_system_clock

!#######################################################################
! several wrappers for old form of interfaces

 subroutine write_version_number (unit, version, tag)
 integer,          intent(in) :: unit
 character(len=*), intent(in) :: version, tag
   call fms_write_version_number ( version, tag, unit )
 end subroutine write_version_number

 function mpp_clock_init (str1, str2, level)
 character(len=*), intent(in) :: str1, str2
 integer,          intent(in) :: level
 integer :: mpp_clock_init
   mpp_clock_init = fms_mpp_clock_init ( trim(str1)//', '//trim(str2), level)
 end function mpp_clock_init

!#######################################################################

end module utilities_mod

