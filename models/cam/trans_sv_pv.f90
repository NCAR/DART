ns_sv_pv

!----------------------------------------------------------------------
! purpose: interface between CAM and DART
!
! method: Read DART state vector ("proprietary" format)
!         Reform state vector back into CAM fields.
!         Replace those fields on the CAM initial file with the new values,
!         preserving all other information on the file.
!
! author: Kevin Raeder 2/21/03
!         based on prog_var_to_vector and vector_to_prog_var by Jeff Anderson
!
!----------------------------------------------------------------------

use model_mod

! read in state vector from DART
call read_dart_vector

! decompose vector back into CAM fields
call vector_to_prog_var (x,siz,vars)
deallocate (x)

! optionally modify vector (into a recognizable pattern for testing)

! write fields to the netCDF initial file
call write_cam_init
deallocate (vars)

end program trans_sv_pv

