! Data Assimilation Research Testbed -- DART
! Copyright 2004-2007, Data Assimilation Research Section
! University Corporation for Atmospheric Research
! Licensed under the GPL -- www.gpl.org/licenses/gpl.html
 
program trans_pv_sv

!----------------------------------------------------------------------
! purpose: interface between MITgcm_ocean and DART
!
! method: Read MITgcm_ocean "snapshot" files of model state
!         Reform fields into a DART state vector (control vector).
!         Write out state vector in "proprietary" format for DART.
!         The output is a "DART restart file" format.
!
! author: Tim Hoar 3/13/08
! more mangling:  nancy collins 14mar08
!
!----------------------------------------------------------------------

! <next few lines under version control, do not edit>
! $URL: http://subversion.ucar.edu/DAReS/DART/trunk/models/MITgcm_ocean/trans_pv_sv.f90 $
! $Id: trans_pv_sv.f90 3258 2008-03-14 15:58:36Z thoar $
! $Revision: 3258 $
! $Date: 2008-03-14 09:58:36 -0600 (Fri, 14 Mar 2008) $

use        types_mod, only : r4, r8
use    utilities_mod, only : get_unit, file_exist, E_ERR, E_WARN, E_MSG, &
                             initialize_utilities, finalize_utilities, &
                             error_handler
use        model_mod, only : MIT_meta_type, read_meta, read_snapshot, &
                             prog_var_to_vector, static_init_model, &
                             get_model_size
use  assim_model_mod, only : assim_model_type, static_init_assim_model, &
                             init_assim_model, get_model_size, set_model_state_vector, &
                             write_state_restart, set_model_time, open_restart_read, &
                             open_restart_write, close_restart, aread_state_restart
use time_manager_mod, only : time_type, read_time

implicit none

! version controlled file description for error handling, do not edit
character(len=128), parameter :: &
   source   = "$URL: http://subversion.ucar.edu/DAReS/DART/trunk/models/MITgcm_ocean/trans_pv_sv.f90 $", &
   revision = "$Revision: 3258 $", &
   revdate  = "$Date: 2008-03-14 09:58:36 -0600 (Fri, 14 Mar 2008) $"

character (len = 128) :: msgstring

! eg. [S,T,U,V,SSH].0000040992.[data,meta]
!
! The '.meta' file contains matlab-format information about shapes, etc.

character (len = 128) :: file_base = '0000000672'
character (len = 128) :: S_filename
character (len = 128) :: T_filename
character (len = 128) :: U_filename
character (len = 128) :: V_filename
character (len = 128) :: SSH_filename

character (len = 128) :: file_out = 'temp_ud', file_time = 'temp_ic'

! Temporary allocatable storage to read in a native format for cam state
type(assim_model_type) :: x
!type(model_type)       :: var
type(time_type)        :: model_time, adv_to_time

real(r4), allocatable  :: S(:,:,:)
real(r4), allocatable  :: T(:,:,:)
real(r4), allocatable  :: U(:,:,:)
real(r4), allocatable  :: V(:,:,:)
real(r4), allocatable  :: SSH(:,:)

real(r8), allocatable  :: x_state(:)

integer                :: file_unit, x_size, timestep

type(MIT_meta_type) :: mitmeta

!
! end of variable decls
! 
!----------------------------------------------------------------------
!
! code start
!

! The strategy here is to use the 'write_state_restart()' routine - 
! which requires an 'assim_model_type' variable.  
! Read the individual variables and pack them
! into one big state vector, add a time type and make nn

call initialize_utilities('trans_pv_sv')

call read_snapshot('XC',SSH,timestep)

! Static init assim model calls static_init_model
call static_init_assim_model()
call init_assim_model(x)

write(*,*)'dimensions are ',shape(SSH)
write(*,*)'model size is ',get_model_size()

stop

! Read the [meta,data] files 

mitmeta = read_meta(file_base,'U')
write(*,*)'timestep is ',timestep

call read_snapshot(file_base,S,timestep,'S')
call read_snapshot(file_base,T,timestep,'T')
call read_snapshot(file_base,U,timestep,'U')
call read_snapshot(file_base,V,timestep,'V')
call read_snapshot(file_base,SSH,timestep,'SSH')

! matlab debug messages
! write(8)mitmeta%dimList
! write(8)S
! close(8)

! Allocate the local state vector and fill it up.

x_size = size(S) + size(T) + size(U) + size(V) + size(SSH)

if ( x_size /= get_model_size() ) then
   write(msgstring,*)'data size ',x_size,' /= ',get_model_size(),' model size'
   call error_handler(E_ERR,'trans_pv_sv',msgstring,source,revision,revdate)
endif

allocate(x_state(x_size))

! transform fields into state vector for DART
call prog_var_to_vector(S,T,U,V,SSH,x_state)

! use the MIT namelist and timestepcount in the meta file to construct
! the current time.  we do not have a restart file to read; this program
! is creating one.

! We're done with x_state, so it can be uselessly filled in aread_state_restart,
! while getting model_time.
!call aread_state_restart(model_time, x_state, file_unit, adv_to_time)
!call set_model_time (x, adv_to_time)
!call close_restart(file_unit)

! Get channel for output,
! write out state vector in "proprietary" format
!file_unit = open_restart_write(file_out)
!call write_state_restart(x, file_unit)
!call close_restart(file_unit)

call finalize_utilities()

end program trans_pv_sv
