! DART software - Copyright UCAR. This open source software is provided
! by UCAR, "as is", without charge, subject to all terms of use at
! http://www.image.ucar.edu/DAReS/DART/DART_download
!

module test_interpolate_mod

!-------------------------------------------------------------------------------
! interpolation test routines for threed sphere locations.
!-------------------------------------------------------------------------------

use             types_mod, only : r8, i8, MISSING_R8

use         utilities_mod, only : error_handler, E_MSG, E_ERR, &
                                  E_MSG, open_file, close_file, do_output

use  netcdf_utilities_mod, only :  nc_create_file, nc_close_file, &
                                   nc_define_dimension, nc_define_double_variable, &
                                   nc_end_define_mode, nc_add_global_creation_time, &
                                   nc_add_attribute_to_variable, nc_put_variable

use          location_mod, only : location_type, set_location, &
                                  VERTISUNDEF, VERTISSURFACE, VERTISLEVEL, &
                                  VERTISPRESSURE, VERTISHEIGHT, VERTISSCALEHEIGHT, &
                                  query_location

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

integer :: vert_type

vert_type = get_location_index(vertcoord_string)
setup_location  = set_location((/ vals, real(vert_type,r8) /))

end function setup_location


!-------------------------------------------------------------------------------
!> Interpolate over a range of lat, lon, and vert values.
!> Returns the number of failures.
!> Exercises model_mod:model_interpolate().

function test_interpolate_range( ens_handle,            &
                                 ens_size,              &
                                 interp_test_dlon,      &
                                 interp_test_dlat,      &
                                 interp_test_dvert,     &
                                 interp_test_vertcoord, &
                                 interp_test_lonrange,  &
                                 interp_test_latrange,  &
                                 interp_test_vertrange, &
                                 quantity_string,       &
                                 verbose )

! The arguments to this function must be the same across all types of location dimensions.

type(ensemble_type)   , intent(inout) :: ens_handle
integer               , intent(in)    :: ens_size
real(r8)              , intent(in)    :: interp_test_dlon
real(r8)              , intent(in)    :: interp_test_dlat
real(r8)              , intent(in)    :: interp_test_dvert
character(len=*)      , intent(in)    :: interp_test_vertcoord
real(r8), dimension(2), intent(in)    :: interp_test_latrange
real(r8), dimension(2), intent(in)    :: interp_test_lonrange
real(r8), dimension(2), intent(in)    :: interp_test_vertrange
character(len=*),       intent(in)    :: quantity_string
logical               , intent(in)    :: verbose

integer :: test_interpolate_range

! Local variables

character(len=*), parameter :: routine = 'test_interpolate_range'

real(r8), allocatable :: lon(:), lat(:), vert(:)
real(r8), allocatable :: field(:,:,:,:)
integer,  allocatable :: all_ios_out(:,:)
real(r8) :: lonrange_top
integer :: nlon, nlat, nvert
integer :: ilon, jlat, kvert, nfailed
character(len=MAX_FILENAME_LEN) :: ncfilename, matlabfilename

integer :: ncid

type(location_type) :: loc
integer :: matunit, ios_out(ens_size)
integer :: quantity_index, vertcoord

test_interpolate_range = 0

if ((interp_test_dlon < 0.0_r8) .or. (interp_test_dlat < 0.0_r8)) then
   if ( do_output() ) then
      write(*,'(A)')    'Skipping the rigorous interpolation test because one of'
      write(*,'(A)')    'interp_test_dlon,interp_test_dlat are < 0.0'
      write(*,'(A,I2)') 'interp_test_dlon  = ',interp_test_dlon
      write(*,'(A,I2)') 'interp_test_dlat  = ',interp_test_dlat
      write(*,'(A,I2)') 'interp_test_dvert = ',interp_test_dvert
   endif
   return
endif

vertcoord = get_location_index(interp_test_vertcoord)
quantity_index = get_index_for_quantity(quantity_string)

! for longitude, allow wrap.
lonrange_top = interp_test_lonrange(2)
if (interp_test_lonrange(2) < interp_test_lonrange(1)) &
   lonrange_top = interp_test_lonrange(2) + 360.0_r8

! round down to avoid exceeding the specified range
nlat  = aint(( interp_test_latrange(2) - interp_test_latrange(1))  / interp_test_dlat) + 1
nlon  = aint((            lonrange_top - interp_test_lonrange(1))  / interp_test_dlon) + 1
nvert = aint((interp_test_vertrange(2) - interp_test_vertrange(1)) / interp_test_dvert) + 1


! netcdf and matlab output
! filenames are hardcoded and constructed in this routine

call create_output_files(ncfilename, matlabfilename, ncid, matunit)

call write_matlab_header(matunit, nlon, nlat, nvert, ens_size)


allocate(lon(nlon), lat(nlat), vert(nvert), field(nlon,nlat,nvert,ens_size))
allocate(all_ios_out(nlon*nlat*nvert,ens_size))

all_ios_out = 0 ! assume successful interpolation for every grid location, all members.
nfailed = 0

do ilon = 1, nlon
   lon(ilon) = interp_test_lonrange(1) + real(ilon-1,r8) * interp_test_dlon
   if (lon(ilon) >= 360.0_r8) lon(ilon) = lon(ilon) - 360.0_r8
   if (lon(ilon) <    0.0_r8) lon(ilon) = lon(ilon) + 360.0_r8
   do jlat = 1, nlat
      lat(jlat) = interp_test_latrange(1) + real(jlat-1,r8) * interp_test_dlat
      do kvert = 1, nvert
         vert(kvert) = interp_test_vertrange(1) + real(kvert-1,r8) * interp_test_dvert
         loc = set_location(lon(ilon), lat(jlat), vert(kvert), vertcoord)

         call model_interpolate(ens_handle, ens_size, loc, quantity_index, &
                                field(ilon,jlat,kvert,:), ios_out)

         call verify_consistent_istatus(ens_size, field(ilon,jlat,kvert,:), ios_out)

         if (any(ios_out /= 0)) then

	    nfailed = nfailed + 1
            ! don't really care which location was causing the failure
            all_ios_out(nfailed,:) = ios_out

            if (verbose) then
               write(string1,*) 'interpolation return code was', ios_out
               write(string2,'(''ilon,jlat,kvert,lon,lat,vert'',3(1x,i6),3(1x,f14.6))') &
                                 ilon,jlat,kvert,lon(ilon),lat(jlat),vert(kvert)
               call error_handler(E_MSG, routine, string1, text2=string2)
            endif
 
         endif
      enddo
   enddo
enddo

call write_matlab_data(matunit, nlon, nlat, nvert, field)
call write_matlab_footer_and_close(matunit)


! Write out the netCDF file for easy exploration.
call write_netcdf_output_and_close(ncid, ens_size, nlon, nlat, nvert, &
                                   interp_test_lonrange, &
                                   interp_test_latrange, &
                                   interp_test_vertrange, &
                                   interp_test_vertcoord, &
                                   quantity_string, lon, lat, vert, field)


if ( do_output() ) then
   write(*,'(A)')     '-------------------------------------------------------------'
   write(*,'(A,I10)') 'total  interpolations : ', nlon*nlat*nvert
   write(*,'(A,I10)') 'failed interpolations : ', nfailed
   write(*,'(A)')     '-------------------------------------------------------------'
endif

call count_error_codes(all_ios_out(1:nfailed,:))

deallocate(lon, lat, vert, field)
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

subroutine write_matlab_header(matunit, nlon, nlat, nvert, ens_size)

integer, intent(in) :: matunit, nlon, nlat, nvert, ens_size

write(matunit,'(''missingvals = '',f12.4,'';'')')MISSING_R8
write(matunit,'(''nlon = '',i8,'';'')')nlon
write(matunit,'(''nlat = '',i8,'';'')')nlat
write(matunit,'(''nvert = '',i8,'';'')')nvert
write(matunit,'(''nens = '',i8,'';'')')ens_size
write(matunit,'(''interptest = [ ... '')')

end subroutine write_matlab_header

!-------------------------------------------------------------------------------
! data array depends on dims

subroutine write_matlab_data(matunit, nlon, nlat, nvert, field)

integer, intent(in) :: matunit, nlon, nlat, nvert
real(r8), intent(in) :: field(:,:,:,:)

integer :: ilon, jlat, kvert

! write all ensemble values for each item
do ilon = 1, nlon
   do jlat = 1, nlat
      do kvert = 1, nvert
         write(matunit,*) field(ilon,jlat,kvert,:)
      enddo
   enddo
enddo

end subroutine write_matlab_data

!-------------------------------------------------------------------------------

subroutine write_matlab_footer_and_close(matunit)

integer, intent(in) :: matunit

write(matunit,'(''];'')')
write(matunit,'(''datmat = reshape(interptest,nvert,nlat,nlon,nens);'')')
write(matunit,'(''datmat = permute(datmat,[4,1,2,3]);'')')
write(matunit,'(''datmat(datmat == missingvals) = NaN;'')')

call close_file(matunit)

end subroutine write_matlab_footer_and_close

!-------------------------------------------------------------------------------

subroutine write_netcdf_output_and_close(ncid, ens_size, nlon, nlat, nvert, &
                                         interp_test_lonrange, &
                                         interp_test_latrange, &
                                         interp_test_vertrange, &
                                         interp_test_vertcoord, &
                                         quantity_string, lon, lat, vert, field)

integer,  intent(in) :: ncid
integer,  intent(in) :: ens_size, nlon, nlat, nvert
real(r8), intent(in) :: interp_test_lonrange(2)
real(r8), intent(in) :: interp_test_latrange(2)
real(r8), intent(in) :: interp_test_vertrange(2)
character(len=*), intent(in) :: interp_test_vertcoord
character(len=*), intent(in) :: quantity_string
real(r8), intent(in) :: lon(:), lat(:), vert(:)
real(r8), intent(in) :: field(:,:,:,:)

integer :: i, imem
character(len=32) :: field_name
character(len=*), parameter :: routine = 'test_interpolate_range'  ! most helpful context?


call nc_add_global_creation_time(ncid, routine)

call nc_define_dimension(ncid, 'lon',  nlon, routine)
call nc_define_dimension(ncid, 'lat',  nlat, routine)
call nc_define_dimension(ncid, 'vert', nvert, routine)

call nc_define_double_variable(ncid, 'lon', 'lon', routine)
call nc_add_attribute_to_variable(ncid, 'lon', 'range', interp_test_lonrange, routine)
call nc_add_attribute_to_variable(ncid, 'lon', 'cartesian_axis', 'X', routine)

call nc_define_double_variable(ncid, 'lat', 'lat', routine)
call nc_add_attribute_to_variable(ncid, 'lat', 'range', interp_test_latrange, routine)
call nc_add_attribute_to_variable(ncid, 'lat', 'cartesian_axis', 'Y', routine)

call nc_define_double_variable(ncid, 'vert', 'vert', routine)
call nc_add_attribute_to_variable(ncid, 'vert', 'range', interp_test_vertrange, routine)
call nc_add_attribute_to_variable(ncid, 'vert', 'cartesian_axis', 'Z', routine)


! create an output variable for each ensemble member

! default for a single ensemble member
field_name = "field"

do imem = 1, ens_size
   if (ens_size > 1) write(field_name,'(A,I2)') "field_",imem

   call nc_define_double_variable(ncid, field_name, (/ 'lon', 'lat', 'vert' /) , routine)
   call nc_add_attribute_to_variable(ncid, field_name, 'long_name', quantity_string, routine)
   call nc_add_attribute_to_variable(ncid, field_name, '_FillValue', MISSING_R8, routine)
   call nc_add_attribute_to_variable(ncid, field_name, 'missing_value', MISSING_R8, routine)
   call nc_add_attribute_to_variable(ncid, field_name, 'interp_test_vertcoord', routine)

enddo

! Leave define mode so we can fill the variables.
call nc_end_define_mode(ncid, routine)

! Fill the dimension variable, then the data field
call nc_put_variable(ncid, 'lon', lon, routine)
call nc_put_variable(ncid, 'lat', lat, routine)
call nc_put_variable(ncid, 'vert', vert, routine)


do imem = 1, ens_size
   if (ens_size > 1) write(field_name,'(A,I2)') "field_",imem

   call nc_put_variable(ncid, field_name, field(:,:,:,imem), routine)
enddo

! close and return
call nc_close_file(ncid, routine)

end subroutine write_netcdf_output_and_close

!-------------------------------------------------------------------------------
!> need to convert the character string for the test vertical coordinate into
!> the corresponding dart index.
!> TODO: FIXME - there should be a function in the locations mod for this.

function  get_location_index(test_vertcoord)

character(len=*) , intent(in) :: test_vertcoord

integer :: get_location_index

select case (test_vertcoord)
   case ('VERTISUNDEF')
      get_location_index = VERTISUNDEF
   case ('VERTISSURFACE')
      get_location_index = VERTISSURFACE
   case ('VERTISLEVEL')
      get_location_index = VERTISLEVEL
   case ('VERTISPRESSURE')
      get_location_index = VERTISPRESSURE
   case ('VERTISHEIGHT')
      get_location_index = VERTISHEIGHT
   case ('VERTISSCALEHEIGHT')
      get_location_index = VERTISSCALEHEIGHT
   case default
      get_location_index = VERTISUNDEF
end select

end function  get_location_index

!-------------------------------------------------------------------------------
! End of test_interpolate_mod
!-------------------------------------------------------------------------------

end module test_interpolate_mod

