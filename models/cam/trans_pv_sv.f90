program trans_pv_sv

!----------------------------------------------------------------------
! purpose: interface between CAM and DART
!
! method: Read CAM 'initial' file (netCDF format).
!         Reform fields into a state vector.
!         Write out state vector in "proprietary" format for DART
!
! author: Kevin Raeder 2/21/03
!         based on prog_var_to_vector and vector_to_prog_var by Jeff Anderson
!
!----------------------------------------------------------------------

use model_mod, only : read_cam_init_size

integer :: num_lons, num_lats, num_levs
character (len = 128) :: file_name = 'H12-24icl.nc'

! read in field values from CAM in initial file
call read_cam_init_size(file_name, num_lons, num_lats, num_levs)

write(*, *) 'lons, lats, levs', num_lons, num_lats, num_levs

! transform fields into state vector for DART
!call prog_var_to_vector(vars,x,siz)
!deallocate (vars)

! write out state vector in "proprietary" format
!call write_dart_vector
!deallocate (x)

end program trans_pv_sv

