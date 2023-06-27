! DART software - Copyright UCAR. This open source software is provided
! by UCAR, "as is", without charge, subject to all terms of use at
! http://www.image.ucar.edu/DAReS/DART/DART_download
!

program create_lat_lon_interp_grid

!-------------------------------------------------------------------------------
! create a lat/lon/vert grid and interpolate model data onto it.
! output a netcdf file with a 3d array of values for plotting with ncview.
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
                                  get_dist, get_location, LocationDims, &
                                  VERTISUNDEF, VERTISSURFACE, VERTISLEVEL, &
                                  VERTISPRESSURE, VERTISHEIGHT, VERTISSCALEHEIGHT, &
                                  query_location

use          obs_kind_mod, only : get_name_for_quantity, get_index_for_quantity

use  ensemble_manager_mod, only : ensemble_type

use             model_mod, only : get_model_size, &
                                  get_state_meta_data, &
                                  model_interpolate

use netcdf

implicit none
private


character(len=*), parameter :: source   = "create_lat_lon_interp_grid"

! for messages
character(len=512) :: string1, string2

contains

!-------------------------------------------------------------------------------
!> Interpolate over a range of x, y, and z values.
!> Exercises model_mod:model_interpolate().
!> This will result in a netCDF file with all salient metadata.

subroutine interp_lat_lon_grid( ens_handle,            &
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

integer :: test_create_lat_lon_grid

! Local variables

character(len=*), parameter :: routine = 'interp_lat_lon_grid'

real(r8), allocatable :: lon(:), lat(:), vert(:)
real(r8), allocatable :: field(:,:,:,:)
real(r8)              :: lonrange_top
integer               :: nlon, nlat, nvert

character(len=*)    :: ncfilename = 'interpolation_results.nc'
character(len=32)   :: field_name(ens_size)
type(location_type) :: loc
integer             :: iunit, ios_out(ens_size), imem
integer             :: quantity_index, vertcoord


quantity_index = get_index_for_quantity(quantity_string)
vertcoord = get_location_index(interp_test_vertcoord)

if (ens_size == 1) then
   field_name(1) = "field"
else
   do imem = 1, ens_size
      write(field_name(:),'(A,I3)') "field_",imem
   enddo
endif

! for longitude, allow wrap.
lonrange_top = interp_test_lonrange(2)
if (interp_test_lonrange(2) < interp_test_lonrange(1)) &
   lonrange_top = interp_test_lonrange(2) + 360.0_r8


! round down to avoid exceeding the specified range
nlat  = aint(( interp_test_latrange(2) - interp_test_latrange(1))  / interp_test_dlat) + 1
nlon  = aint((            lonrange_top - interp_test_lonrange(1))  / interp_test_dlon) + 1
nvert = aint((interp_test_vertrange(2) - interp_test_vertrange(1)) / interp_test_dvert) + 1


allocate(lon(nlon), lat(nlat), vert(nvert), field(nlon,nlat,nvert,ens_size))
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

         if (any(ios_out(:) /= 0)) then

            nfailed    = nfailed + 1

            if (verbose) then
               write(string1,*) 'interpolation return code was', ios_out
               write(string2,'(''ilon,jlat,kvert,lon,lat,vert'',3(1x,i6),3(1x,f14.6))') &
                                 ilon,jlat,kvert,lon(ilon),lat(jlat),vert(kvert)
               call error_handler(E_MSG, routine, string1, &
                                  source, text2=string2)
            endif
         endif
      enddo
   enddo
enddo

! Write out the netCDF file for easy exploration.

ncid = nc_create_file(ncfilename,'create_lat_lon_interp_grid')

call nc_add_global_creation_time(ncid)

! Define dimensions

call nc_define_dimension(ncid, 'lon', nlon)
call nc_define_dimension(ncid, 'lat', nlat)
call nc_define_dimension(ncid, 'vert, nvert)

! Define variables

call nc_define_double_scalar(ncid, 'lon',  'lon')
call nc_define_double_scalar(ncid, 'lat',  'lat')
call nc_define_double_scalar(ncid, 'vert', 'vert')

call nc_add_attribute_to_variable(ncid, 'lon',  'range', interp_test_lonrange)
call nc_add_attribute_to_variable(ncid, 'lat',  'range', interp_test_latrange)
call nc_add_attribute_to_variable(ncid, 'vert', 'range', interp_test_vertrange)

call nc_add_attribute_to_variable(ncid, 'lon',  'cartesian_axis', 'X')
call nc_add_attribute_to_variable(ncid, 'lat',  'cartesian_axis', 'Y')
call nc_add_attribute_to_variable(ncid, 'vert', 'cartesian_axis', 'Z')

! loop over ensemble members

do imem = 1, ens_size

   call nc_define_double_variable(ncid, field_name(imem), (/ 'vert', 'lat', 'lon' /)

   call nc_add_attribute_to_variable(ncid, field_name(imem), 'long_name',     quantity_string)
   call nc_add_attribute_to_variable(ncid, field_name(imem), '_FillValue',    MISSING_R8)
   call nc_add_attribute_to_variable(ncid, field_name(imem), 'missing_value', MISSING_R8)

   call nc_add_attribute_to_variable(ncid, field_name(imem), 'interp_test_vertcoord', interp_test_vertcoord)
enddo

! Leave define mode so we can fill the variables.
call nc_end_define_mode(ncid, routine)

! Fill the variables
call nc_put_variable(ncid, 'lon',  lon)
call nc_put_variable(ncid, 'lat',  lat)
call nc_put_variable(ncid, 'vert', vert)

do imem = 1, ens_size
   call nc_put_variable(ncid, field_name(imem), field(:,:,:,imem)
enddo

! tidy up
call nc_close_file(ncid, routine)

deallocate(lon, lat, vert, field)


end subroutine interp_lat_lon_grid



end program create_lat_lon_grid

