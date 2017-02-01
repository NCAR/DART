! DART software - Copyright UCAR. This open source software is provided
! by UCAR, "as is", without charge, subject to all terms of use at
! http://www.image.ucar.edu/DAReS/DART/DART_download
!
! $Id$

program dart_to_lmdz

!----------------------------------------------------------------------
! purpose: interface between LMDZ and DART
!
! method: Read DART state vector ("proprietary" format), including time(s).
!         Reform state vector back into LMDZ fields.
!         Replace those fields on the LMDZ initial file with the new values,
!         preserving all other information on the file.
!----------------------------------------------------------------------

use       types_mod, only : r8
use   utilities_mod, only : open_file, close_file, &
                            initialize_utilities, finalize_utilities, &
                            logfileunit, nmlfileunit, do_nml_file, do_nml_term, &
                            check_namelist_read, find_namelist_in_file
use       model_mod, only : data_2d_type,data_3d_type, init_model_instance, write_lmdz_init, &
                            vector_to_prog_var, static_init_model, get_model_size
use assim_model_mod, only : aread_state_restart, open_restart_read, close_restart
use time_manager_mod, only : time_type, print_time, print_date

implicit none

! version controlled file description for error handling, do not edit
character(len=128), parameter :: &
   source   = "$URL$", &
   revision = "$Revision$", &
   revdate  = "$Date$"

!------------------------------------------------------------------
! The namelist variables
!------------------------------------------------------------------

character (len = 128) :: dart_to_lmdz_input_file  = 'dart_ics'
character (len = 128) :: dart_to_lmdz_output_file = 'start.nc'
logical               :: advance_time_present    = .true.

namelist /dart_to_lmdz_nml/ dart_to_lmdz_input_file,  &
                           dart_to_lmdz_output_file, &
                           advance_time_present

!----------------------------------------------------------------------

type(data_2d_type)     :: PS_local
type(data_3d_type)     :: T_local,U_local,V_local,Q_local,CLDLIQ_local
type(time_type)        :: model_time, adv_to_time
real(r8), allocatable  :: statevector(:)
integer                :: file_unit, vecsize, iunit, io

!-----------------------------------------------------------------------------
! start of program
!-----------------------------------------------------------------------------

call initialize_utilities('dart_to_lmdz')

! Read the namelist to get the input filename and whether
! a simple restart/ic file or a model advance file.
call find_namelist_in_file("input.nml", "dart_to_lmdz_nml", iunit)
read(iunit, nml = dart_to_lmdz_nml, iostat = io)
call check_namelist_read(iunit, io, "dart_to_lmdz_nml")

! Record the namelist values used for the run 
if (do_nml_file()) write(nmlfileunit, nml=dart_to_lmdz_nml)
if (do_nml_term()) write(     *     , nml=dart_to_lmdz_nml)

! initialize model
call static_init_model()

vecsize = get_model_size()
allocate(statevector(vecsize))

! Allocate the instance of the lmdz model type for storage
call init_model_instance(PS_local,T_local,U_local,V_local,Q_local,CLDLIQ_local)

! Get file for DART vector input
file_unit = open_restart_read(dart_to_lmdz_input_file)


! read in model time and state vector from DART
if (advance_time_present) then
   call aread_state_restart(model_time, statevector, file_unit, adv_to_time)
else
   call aread_state_restart(model_time, statevector, file_unit)
endif
call close_restart(file_unit)

! decompose vector back into LMDZ fields
call vector_to_prog_var (statevector, PS_local,T_local,U_local,V_local, &
                                                   Q_local,CLDLIQ_local)

deallocate (statevector)

! write fields to the netCDF initial file.
call write_lmdz_init(dart_to_lmdz_output_file, PS_local,T_local,U_local, &
                               V_local,Q_local,CLDLIQ_local, model_time)

call finalize_utilities()

end program dart_to_lmdz

! <next few lines under version control, do not edit>
! $URL$
! $Id$
! $Revision$
! $Date$
