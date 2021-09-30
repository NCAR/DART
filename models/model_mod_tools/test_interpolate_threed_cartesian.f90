! DART software - Copyright UCAR. This open source software is provided
! by UCAR, "as is", without charge, subject to all terms of use at
! http://www.image.ucar.edu/DAReS/DART/DART_download
!

module test_interpolate_mod

!-------------------------------------------------------------------------------
! interpolation test routines for threed cartesian locations.
!-------------------------------------------------------------------------------

use             types_mod, only : r8, i8, MISSING_R8

use         utilities_mod, only : error_handler, E_MSG, E_ERR, &
                                  E_MSG, open_file, close_file, do_output

use  netcdf_utilities_mod, only :  nc_create_file, nc_close_file, &
                                   nc_define_dimension, nc_define_double_variable, &
                                   nc_end_define_mode, nc_add_global_creation_time, &
                                   nc_add_attribute_to_variable, nc_put_variable

use          location_mod, only : location_type, set_location
                                 
use          obs_kind_mod, only : get_index_for_quantity

use  ensemble_manager_mod, only : ensemble_type

use model_check_utilities_mod, only : count_error_codes, &
                                      verify_consistent_istatus

use             model_mod, only : model_interpolate

use netcdf

implicit none
private

public :: setup_location, &
          test_interpolate_range

! filenames are 256, messages are 512
integer, parameter :: MAX_FILENAME_LEN = 256
integer, parameter :: MAX_MSG_LEN = 512

! for messages
character(len=MAX_MSG_LEN) :: string1, string2

contains

!-------------------------------------------------------------------------------
! convert an array of reals into a location for this type

function setup_location(vals, vertcoord_string)

real(r8),         intent(in) :: vals(:)
character(len=*), intent(in) :: vertcoord_string
type(location_type)          :: setup_location

setup_location = set_location(vals(1), vals(2), vals(3))

end function setup_location

!-------------------------------------------------------------------------------
!> Interpolate over a range of x, y, and z values.
!> Returns the number of failures.
!> Exercises model_mod:model_interpolate().

function test_interpolate_range( ens_handle,            &
                                 ens_size,              &
                                 interp_test_dx,        &
                                 interp_test_dy,        &
                                 interp_test_dz,        &
                                 interp_test_vertcoord, &
                                 interp_test_xrange,    &
                                 interp_test_yrange,    &
                                 interp_test_zrange,    &
                                 quantity_string,       &
                                 verbose )

! The arguments to this function must be the same across all types of location dimensions.

type(ensemble_type)   , intent(inout) :: ens_handle
integer               , intent(in)    :: ens_size
real(r8)              , intent(in)    :: interp_test_dx
real(r8)              , intent(in)    :: interp_test_dy
real(r8)              , intent(in)    :: interp_test_dz
character(len=*)      , intent(in)    :: interp_test_vertcoord
real(r8), dimension(2), intent(in)    :: interp_test_xrange
real(r8), dimension(2), intent(in)    :: interp_test_yrange
real(r8), dimension(2), intent(in)    :: interp_test_zrange
character(len=*),       intent(in)    :: quantity_string
logical               , intent(in)    :: verbose

integer :: test_interpolate_range

! Local variables

character(len=*), parameter :: routine = 'test_interpolate_range'

real(r8), allocatable :: X(:), Y(:), Z(:)
real(r8), allocatable :: field(:,:,:,:)
integer,  allocatable :: all_ios_out(:,:)
integer               :: nx, ny, nz
integer               :: i, j, k, nfailed
character(len=MAX_FILENAME_LEN) :: ncfilename, matlabfilename

integer :: ncid

type(location_type) :: loc
integer :: matunit, ios_out(ens_size)
integer :: quantity_index

test_interpolate_range = 0
quantity_index = get_index_for_quantity(quantity_string)

! round down to avoid exceeding the specified range
nx = aint((interp_test_xrange(2) - interp_test_xrange(1))/interp_test_dx) + 1
ny = aint((interp_test_yrange(2) - interp_test_yrange(1))/interp_test_dy) + 1
nz = aint((interp_test_zrange(2) - interp_test_zrange(1))/interp_test_dz) + 1


! netcdf and matlab output
! filenames are hardcoded and constructed in this routine

call create_output_files(ncfilename, matlabfilename, ncid, matunit)

call write_matlab_header(matunit, nx, ny, nz, ens_size)


allocate(X(nx), Y(ny), Z(nz), field(nx,ny,nz,ens_size))
allocate(all_ios_out(nx*ny*nz,ens_size))

all_ios_out = 0 ! assume successful interpolation for every grid location, all members.
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

         call verify_consistent_istatus(ens_size, field(i,j,k,:), ios_out)

         if (any(ios_out /= 0)) then

            nfailed = nfailed + 1
            ! don't really care which location was causing the failure
            all_ios_out(nfailed,:) = ios_out

            if (verbose) then
               write(string1,*) 'interpolation return code was', ios_out
               write(string2,'(''i,j,k,X,Y,Z'',3(1x,i6),3(1x,f14.6))') &
				 i,j,k,X(i),Y(j),Z(k)
               call error_handler(E_MSG, routine, string1, text2=string2)
            endif
         endif
      enddo
   enddo
enddo

call write_matlab_data(matunit, nx, ny, nz, field)
call write_matlab_footer_and_close(matunit)


! Write out the netCDF file for easy exploration.
call write_netcdf_output_and_close(ncid, ens_size, nx, ny, nz, &
                                   interp_test_xrange, &
                                   interp_test_yrange, &
                                   interp_test_zrange, &
                                   quantity_string, X, Y, Z, field)


if ( do_output() ) then
   write(*,'(A)')     '-------------------------------------------------------------'
   write(*,'(A,I10)') 'total  interpolations : ', nx*ny*nz
   write(*,'(A,I10)') 'failed interpolations : ', nfailed
   write(*,'(A)')     '-------------------------------------------------------------'
endif

call count_error_codes(all_ios_out(1:nfailed,:))

deallocate(X, Y, Z, field)
deallocate(all_ios_out)

test_interpolate_range = nfailed

end function test_interpolate_range


!-------------------------------------------------------------------------------
! THIS ROUTINE COULD MOVE INTO THE COMMON CODE FILE
! construct filenames and open output files for test_interpolate_range

subroutine create_output_files(ncfilename, matlabfilename, &
                               ncid, matid)

character(len=MAX_FILENAME_LEN), intent(out) :: ncfilename, matlabfilename
integer, intent(out) :: ncid, matid

character(len=*), parameter :: output_file_prefix = 'check_me'
character(len=*), parameter :: output_file_suffix = '_interptest'

character(len=*), parameter :: routine = 'create_output_files'

if (len(output_file_prefix) + len(output_file_suffix) + 3 > MAX_FILENAME_LEN) then
   write(string1,*) 'unexpected error, constructed filename too long'
   write(string2,*) 'contact DART support'
   call error_handler(E_ERR, routine, string1, text2=string2)
endif

ncfilename = trim(output_file_prefix)//output_file_suffix//'.nc'
ncid = nc_create_file(ncfilename, 'test_interpolate_range')

matlabfilename = trim(output_file_prefix)//output_file_suffix//'.m'
matid = open_file(matlabfilename, action='write')


end subroutine create_output_files

!-------------------------------------------------------------------------------
! header depends on number of dims - could be in common code if #dims passed

subroutine write_matlab_header(matunit, nx, ny, nz, ens_size)

integer, intent(in) :: matunit, nx, ny, nz, ens_size

write(matunit,'(''missingvals = '',f12.4,'';'')')MISSING_R8
write(matunit,'(''nx = '',i8,'';'')')nx
write(matunit,'(''ny = '',i8,'';'')')ny
write(matunit,'(''nz = '',i8,'';'')')nz
write(matunit,'(''nens = '',i8,'';'')')ens_size
write(matunit,'(''interptest = [ ... '')')

end subroutine write_matlab_header

!-------------------------------------------------------------------------------
! data array depends on dims

subroutine write_matlab_data(matunit, nx, ny, nz, field)

integer, intent(in) :: matunit, nx, ny, nz
real(r8), intent(in) :: field(:,:,:,:)

integer :: i, j, k

! write all ensemble values for each item
do i = 1, nx
   do j = 1, ny
      do k = 1, nz
         write(matunit,*) field(i,j,k,:)
      enddo
   enddo
enddo

end subroutine write_matlab_data

!-------------------------------------------------------------------------------

subroutine write_matlab_footer_and_close(matunit)

integer, intent(in) :: matunit

write(matunit,'(''];'')')
write(matunit,'(''datmat = reshape(interptest,nz,ny,nx,nens);'')')
write(matunit,'(''datmat = permute(datmat,[4,1,2,3]);'')')
write(matunit,'(''datmat(datmat == missingvals) = NaN;'')')

call close_file(matunit)

end subroutine write_matlab_footer_and_close

!-------------------------------------------------------------------------------

subroutine write_netcdf_output_and_close(ncid, ens_size, nx, ny, nz, &
                                         interp_test_xrange, &
                                         interp_test_yrange, &
                                         interp_test_zrange, &
                                         quantity_string, X, Y, Z, field)

integer,  intent(in) :: ncid
integer,  intent(in) :: ens_size, nx, ny, nz
real(r8), intent(in) :: interp_test_xrange(2)
real(r8), intent(in) :: interp_test_yrange(2)
real(r8), intent(in) :: interp_test_zrange(2)
character(len=*), intent(in) :: quantity_string
real(r8), intent(in) :: X(:), Y(:), Z(:)
real(r8), intent(in) :: field(:,:,:,:)

integer :: i, imem
character(len=32) :: field_name
character(len=*), parameter :: routine = 'test_interpolate_range'  ! most helpful context?


call nc_add_global_creation_time(ncid, routine)

call nc_define_dimension(ncid, 'X', nx, routine)
call nc_define_dimension(ncid, 'Y', ny, routine)
call nc_define_dimension(ncid, 'Z', nz, routine)

call nc_define_double_variable(ncid, 'X', 'X', routine)
call nc_add_attribute_to_variable(ncid, 'X', 'range', interp_test_xrange, routine)
call nc_add_attribute_to_variable(ncid, 'X', 'cartesian_axis', 'X', routine)

call nc_define_double_variable(ncid, 'Y', 'Y', routine)
call nc_add_attribute_to_variable(ncid, 'Y', 'range', interp_test_yrange, routine)
call nc_add_attribute_to_variable(ncid, 'Y', 'cartesian_axis', 'Y', routine)

call nc_define_double_variable(ncid, 'Z', 'Z', routine)
call nc_add_attribute_to_variable(ncid, 'Z', 'range', interp_test_zrange, routine)
call nc_add_attribute_to_variable(ncid, 'Z', 'cartesian_axis', 'Z', routine)


! create an output variable for each ensemble member

! default for a single ensemble member
field_name = "field"

do imem = 1, ens_size
   if (ens_size > 1) write(field_name,'(A,I2)') "field_",imem

   call nc_define_double_variable(ncid, field_name, (/ 'X', 'Y', 'Z' /) , routine)
   call nc_add_attribute_to_variable(ncid, field_name, 'long_name', quantity_string, routine)
   call nc_add_attribute_to_variable(ncid, field_name, '_FillValue', MISSING_R8, routine)
   call nc_add_attribute_to_variable(ncid, field_name, 'missing_value', MISSING_R8, routine)

enddo

! Leave define mode so we can fill the variables.
call nc_end_define_mode(ncid, routine)

! Fill the dimension variable, then the data field
call nc_put_variable(ncid, 'X', X, routine)
call nc_put_variable(ncid, 'Y', Y, routine)
call nc_put_variable(ncid, 'Z', Z, routine)


do imem = 1, ens_size
   if (ens_size > 1) write(field_name,'(A,I2)') "field_",imem

   call nc_put_variable(ncid, field_name, field(:,:,:,imem), routine)
enddo

! close and return
call nc_close_file(ncid, routine)

end subroutine write_netcdf_output_and_close

!-------------------------------------------------------------------------------
! End of test_interpolate_mod
!-------------------------------------------------------------------------------

end module test_interpolate_mod

