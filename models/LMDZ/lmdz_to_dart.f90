! DART software - Copyright UCAR. This open source software is provided
! by UCAR, "as is", without charge, subject to all terms of use at
! http://www.image.ucar.edu/DAReS/DART/DART_download
!
! $Id$

program lmdz_to_dart

!----------------------------------------------------------------------
! purpose: interface between LMDZ and DART
!
! method: Read LMDZ 'initial' file (netCDF format) for model state and time.
!         Reform fields into a DART state vector.
!         Write out state vector in "proprietary" format for DART
!----------------------------------------------------------------------

use        types_mod, only : r8
use    utilities_mod, only : initialize_utilities, finalize_utilities, do_output,     &
                             check_namelist_read, find_namelist_in_file, nmlfileunit, &
                             do_nml_file, do_nml_term
use        model_mod, only : data_2d_type,data_3d_type, init_model_instance, end_model_instance, &
                             prog_var_to_vector, read_lmdz_init, static_init_model,PS,T,U,V,Q,CLDLIQ
use  assim_model_mod, only : static_init_assim_model, get_model_size, &
                             open_restart_write, awrite_state_restart, close_restart
use time_manager_mod, only : time_type

implicit none

! version controlled file description for error handling, do not edit
character(len=128), parameter :: &
   source   = "$URL$", &
   revision = "$Revision$", &
   revdate  = "$Date$"

!-----------------------------------------------------------------------
! namelist parameters with default values.
!-----------------------------------------------------------------------

character (len = 128) :: lmdz_to_dart_input_file  = 'start.nc'
character (len = 128) :: lmdz_to_dart_output_file = 'dart_ics'

namelist /lmdz_to_dart_nml/ lmdz_to_dart_input_file, lmdz_to_dart_output_file

! allocatable storage to read in a native format for lmdz state

real(r8), allocatable  :: statevector(:)
type(time_type)        :: model_time
integer                :: iunit, x_size, io

call initialize_utilities('lmdz_to_dart')

! Read the namelist entry
call find_namelist_in_file("input.nml", "lmdz_to_dart_nml", iunit)
read(iunit, nml = lmdz_to_dart_nml, iostat = io)
call check_namelist_read(iunit, io, "lmdz_to_dart_nml")

! Record the namelist values

! do_nml_file () :  *** return whether nml should be written to nml file 

if (do_nml_file()) write(nmlfileunit, nml=lmdz_to_dart_nml)
if (do_nml_term()) write(     *     , nml=lmdz_to_dart_nml)

! Static init assim model sets the output file format (binary/ascii)
! call static_init_model() 

call static_init_assim_model()

! Allocate the local state vector
x_size = get_model_size()
allocate(statevector(x_size))


! Allocate the instance of the lmdz model type for storage
call init_model_instance(PS, T, U, V, Q, CLDLIQ)

! Read the file lmdz state fragments into var;
! transform fields into state vector for DART


call read_lmdz_init(lmdz_to_dart_input_file, model_time)

call prog_var_to_vector(statevector, PS, T, U, V, Q, CLDLIQ)

call end_model_instance(PS, T, U, V, Q, CLDLIQ)

! write out state vector in "proprietary" format
iunit = open_restart_write(lmdz_to_dart_output_file)
call awrite_state_restart(model_time, statevector, iunit)
call close_restart(iunit)

call finalize_utilities()

end program lmdz_to_dart

! <next few lines under version control, do not edit>
! $URL$
! $Id$
! $Revision$
! $Date$
