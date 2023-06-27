! DART software - Copyright UCAR. This open source software is provided
! by UCAR, "as is", without charge, subject to all terms of use at
! http://www.image.ucar.edu/DAReS/DART/DART_download
!

program create_X_Y_Z_interp_grid

!-------------------------------------------------------------------------------
! create a cartesian X, Y, Z grid and interpolate model data onto it.
! output a netcdf file with a 3d array of values for plotting with ncview.
! for a lat/lon version see the companion program create_lat_lon_interp_grid
!-------------------------------------------------------------------------------

use             types_mod, only : r8, i8, MISSING_R8, metadatalength

use         utilities_mod, only : register_module, error_handler, E_MSG, E_ERR, &
                                  initialize_utilities, finalize_utilities,     &
                                  find_namelist_in_file, check_namelist_read,   &
                                  E_MSG, open_file, close_file, do_output

use  netcdf_utilities_mod, only : nc_check, nc_create_file, nc_close_file, &
                                  nc_end_define_mode, nc_add_attribute_to_variable, &
                                  nc_define_dimension, nc_define_double_variable, &
                                  nc_put_variable

use          location_mod, only : location_type, set_location, write_location,  &
                                  get_dist, get_location, LocationDims

use          obs_kind_mod, only : get_name_for_quantity, get_index_for_quantity

use  ensemble_manager_mod, only : ensemble_type

use             model_mod, only : get_model_size, &
                                  get_state_meta_data, &
                                  model_interpolate

use netcdf

implicit none
private


character(len=*), parameter :: source   = "create_X_Y_Z_interp_grid"

! for messages
character(len=512) :: string1, string2

contains

!-------------------------------------------------------------------------------
!> Interpolate over a range of x, y, and z values.
!> Exercises model_mod:model_interpolate().
!> This will result in a netCDF file with all salient metadata.

subroutine interp_X_Y_Z_grid( ens_handle,            &
                              ens_size,              &
                              interp_test_dx,        &
                              interp_test_dy,        &
                              interp_test_dz,        &
                              interp_test_xrange,    &
                              interp_test_yrange,    &
                              interp_test_zrange,    &
                              quantity_string,       &
                              verbose )

type(ensemble_type)   , intent(inout) :: ens_handle
integer               , intent(in)    :: ens_size
real(r8)              , intent(in)    :: interp_test_dx
real(r8)              , intent(in)    :: interp_test_dy
real(r8)              , intent(in)    :: interp_test_dz
real(r8), dimension(2), intent(in)    :: interp_test_xrange
real(r8), dimension(2), intent(in)    :: interp_test_yrange
real(r8), dimension(2), intent(in)    :: interp_test_zrange
character(len=*),       intent(in)    :: quantity_string
logical               , intent(in)    :: verbose

integer :: test_create_X_Y_Z_grid

! Local variables

character(len=*), parameter :: routine = 'interp_X_Y_Z_grid'

real(r8), allocatable :: X(:), Y(:), Z(:)
real(r8), allocatable :: field(:,:,:,:)
integer               :: nx, ny, nz
integer               :: i, j, k, nfailed

character(len=*)    :: ncfilename = 'interpolation_results.nc'
character(len=32)   :: field_name(ens_size)
type(location_type) :: loc
integer             :: iunit, ios_out(ens_size), imem
integer             :: quantity_index


quantity_index = get_index_for_quantity(quantity_string)

if (ens_size == 1) then
   field_name(1) = "field"
else
   do imem = 1, ens_size
      write(field_name(imem),'(A,I3)') "field_",imem
   enddo
endif


! round down to avoid exceeding the specified range
nx = aint((interp_test_xrange(2) - interp_test_xrange(1))/interp_test_dx) + 1
ny = aint((interp_test_yrange(2) - interp_test_yrange(1))/interp_test_dy) + 1
nz = aint((interp_test_zrange(2) - interp_test_zrange(1))/interp_test_dz) + 1

allocate(X(nx), Y(ny), Z(nz), field(nx,ny,nz,ens_size))

nfailed = 0

do i = 1, nx
   X(i) = interp_test_xrange(1) + real(i-1,r8) * interp_test_dx
   do j = 1, ny
      Y(j) = interp_test_yrange(1) + real(j-1,r8) * interp_test_dy
      do k = 1, nz
         Z(k) = interp_test_zrange(1) + real(k-1,r8) * interp_test_dz

         loc  = set_location(X(i), Y(j), Z(k))

         call model_interpolate(ens_handle, ens_size, loc, quantity_index, &
                                field(i,j,k,:), ios_out)

         if (any(ios_out(:) /= 0)) then

            nfailed    = nfailed + 1

            if (verbose) then
               write(string1,*) 'interpolation return code was', ios_out
               write(string2,'(''i,j,k,X,Y,Z'',3(1x,i6),3(1x,f14.6))') i,j,k,X(i),Y(j),Z(k)
               call error_handler(E_MSG, routine, string1, &
                                  source, text2=string2)
            endif
         endif
      enddo
   enddo
enddo

if ( do_output() ) then
   write(*,'(A)')     '-------------------------------------------------------------'
   write(*,'(A,I10)') 'total  interpolations : ', nx*ny*nz
   write(*,'(A,I10)') 'failed interpolations : ', nfailed
   write(*,'(A)')     '-------------------------------------------------------------'
endif


! Write out the netCDF file for easy exploration.

ncid = nc_create_file(ncfilename,'create_X_Y_Z_interp_grid')

call nc_add_global_creation_time(ncid)

! Define dimensions

call nc_define_dimension(ncid, 'X', nx)
call nc_define_dimension(ncid, 'Y', ny)
call nc_define_dimension(ncid, 'Z', nz)

! Define variables

call nc_define_double_scalar(ncid, 'X', 'X')
call nc_define_double_scalar(ncid, 'Y', 'Y')
call nc_define_double_scalar(ncid, 'Z', 'Z')

call nc_add_attribute_to_variable(ncid, 'X', 'range', interp_test_xrange)
call nc_add_attribute_to_variable(ncid, 'Y', 'range', interp_test_yrange)
call nc_add_attribute_to_variable(ncid, 'Z', 'range', interp_test_zrange)

call nc_add_attribute_to_variable(ncid, 'X', 'cartesian_axis', 'X cartesian axis')
call nc_add_attribute_to_variable(ncid, 'Y', 'cartesian_axis', 'Y cartesian axis')
call nc_add_attribute_to_variable(ncid, 'Z', 'cartesian_axis', 'Z cartesian axis')

! loop over ensemble members

do imem = 1, ens_size
   call nc_define_double_variable(ncid, field_name(imem), (/ 'Z', 'Y', 'X' /)

   call nc_add_attribute_to_variable(ncid, field_name(imem), 'long_name',     quantity_string)
   call nc_add_attribute_to_variable(ncid, field_name(imem), '_FillValue',    MISSING_R8)
   call nc_add_attribute_to_variable(ncid, field_name(imem), 'missing_value', MISSING_R8)
enddo

! Leave define mode so we can fill the variables.
call nc_end_define_mode(ncid, routine)

! Fill the variables
call nc_put_variable(ncid, 'X', X)
call nc_put_variable(ncid, 'Y', Y)
call nc_put_variable(ncid, 'Z', Z)

do imem = 1, ens_size
   call nc_put_variable(ncid, field_name(imem), field(:,:,:,imem)
enddo

! tidy up
call nc_close_file(ncid, routine)

deallocate(X, Y, Z, field)


end subroutine interp_X_Y_Z_grid



end program create_X_Y_Z_grid

