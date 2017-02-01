! DART software - Copyright UCAR. This open source software is provided
! by UCAR, "as is", without charge, subject to all terms of use at
! http://www.image.ucar.edu/DAReS/DART/DART_download
!
! $Id$
 
program trans_pv_sv

!----------------------------------------------------------------------
! purpose: interface between MITgcm_ocean and DART
!
! method: Read MITgcm_ocean "snapshot" files of model state
!         Reform fields into a DART state vector (control vector).
!         Write out state vector in "proprietary" format for DART.
!         The output is a "DART restart file" format.
! 
! USAGE: to read all the restartfiles known to model_mod:progvarnames()
!        for timestepcount 40992 (i.e.  S.0000040992.data ) you must 
!        redirect the timestepcount to the program. If the read fails,
!        the default timestepcount is 'zero' ...
!
! trans_pv_sv < 40992
!
! author: Tim Hoar 3/13/08
!
!----------------------------------------------------------------------

use        types_mod, only : r4, r8
use    utilities_mod, only : E_ERR, E_WARN, E_MSG, error_handler, logfileunit, &
                             initialize_utilities, finalize_utilities
use        model_mod, only : snapshot_files_to_sv, static_init_model, &
                             get_model_size, timestep_to_DARTtime
use  assim_model_mod, only : awrite_state_restart, open_restart_write, close_restart
use time_manager_mod, only : time_type, print_time, print_date

implicit none

! version controlled file description for error handling, do not edit
character(len=256), parameter :: source   = &
   "$URL$"
character(len=32 ), parameter :: revision = "$Revision$"
character(len=128), parameter :: revdate  = "$Date$"

character (len = 128) :: msgstring

! Reading from stdin 
! eg. [S,T,U,V,Eta].0000040992.[data,meta]
integer :: timestep = 40992
character (len = 128) :: file_base
character (len = 128) :: file_out  = 'assim_model_state_ud'

integer                :: io, iunit, x_size
type(time_type)        :: model_time
real(r8), allocatable  :: statevector(:)

!----------------------------------------------------------------------

call initialize_utilities('trans_pv_sv')

! Call model_mod:static_init_model(), which reads the namelists
! to set calendar type, starting date, deltaT, etc.

call static_init_model()

! Here's the tricky part ... the only piece of information we need
! is the timestepcount. Given that, we can construct the filenames, etc.
! So ... rather than use a namelist (which we'd need to rewrite every timestep)
! I am going to read it from stdin. Opens up the possibility of abuse if
! someone calls this program without a syntax like   trans_pv_sv < 40992
write(     *     ,*)'Waiting for timestep input ...'
write(logfileunit,*)'Waiting for timestep input ...'
read(*,*,iostat=io)timestep
if (io /= 0) then
   write(*,*)'ERROR trans_pv_sv - unable to read timestep from stdin.'
   msgstring = 'unable to read timestep from stdin'
   call error_handler(E_ERR,"trans_pv_sv", msgstring, source, revision, revdate)
endif

write(file_base,'(i10.10)')timestep

write(*,*)'Trying to read files like yyy.'//trim(file_base)//'.data'

! Use the MIT namelist and timestepcount in the meta file to construct
! the current time, allocate the local state vector, and fill

!model_time = timestep_to_DARTtime(0)
!call print_time(model_time,'time for timestep 0')
!call print_date(model_time,'time for timestep 0')

model_time = timestep_to_DARTtime(timestep)
!model_time = timestep_to_DARTtime(0)   ! ignore ... to generate restarts

call print_time(model_time,'time for '//file_base)
call print_date(model_time,'date for '//file_base)


x_size = get_model_size()
allocate(statevector(x_size))
call snapshot_files_to_sv(timestep, statevector) ! model_mod() knows all this

! could also compare the timestep from the snapshot file to model_time ...
! extra layer of bulletproofing.

iunit = open_restart_write(file_out)
call awrite_state_restart(model_time, statevector, iunit)
call close_restart(iunit)
call finalize_utilities()

end program trans_pv_sv

! <next few lines under version control, do not edit>
! $URL$
! $Id$
! $Revision$
! $Date$
