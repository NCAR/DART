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
 
use     mpi_utilities_mod, only : my_task_id, task_count, send_to, receive_from

use  netcdf_utilities_mod, only : nc_check, nc_create_file, nc_close_file, &
                                  nc_end_define_mode

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

! for messages
character(len=512) :: string1, string2, string3

contains

!-------------------------------------------------------------------------------
!> Interpolate over a range of lat, lon, and vert values.
!> Returns the number of failures.
!> Exercises model_mod:model_interpolate().
!> This will result in a netCDF file with all salient metadata.

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
character(len=128) :: ncfilename, txtfilename

character(len=8)      :: crdate      ! needed by F90 DATE_AND_TIME intrinsic
character(len=10)     :: crtime      ! needed by F90 DATE_AND_TIME intrinsic
character(len=5)      :: crzone      ! needed by F90 DATE_AND_TIME intrinsic
integer, dimension(8) :: values      ! needed by F90 DATE_AND_TIME intrinsic

integer :: ncid, nlonDimID, nlatDimID, nvertDimID
integer :: VarID(ens_size), lonVarID, latVarID, vertVarID

character(len=256)  :: output_file = 'check_me'
character(len=32)   :: field_name
type(location_type) :: loc
integer :: iunit, ios_out(ens_size), imem
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

write( ncfilename,'(a,a)')trim(output_file),'_interptest.nc'
write(txtfilename,'(a,a)')trim(output_file),'_interptest.m'

! for longitude, allow wrap.
lonrange_top = interp_test_lonrange(2)
if (interp_test_lonrange(2) < interp_test_lonrange(1)) &
   lonrange_top = interp_test_lonrange(2) + 360.0_r8

! round down to avoid exceeding the specified range
nlat  = aint(( interp_test_latrange(2) - interp_test_latrange(1))  / interp_test_dlat) + 1
nlon  = aint((            lonrange_top - interp_test_lonrange(1))  / interp_test_dlon) + 1
nvert = aint((interp_test_vertrange(2) - interp_test_vertrange(1)) / interp_test_dvert) + 1

iunit = open_file(trim(txtfilename), action='write')
write(iunit,'(''missingvals = '',f12.4,'';'')')MISSING_R8
write(iunit,'(''nlon = '',i8,'';'')')nlon
write(iunit,'(''nlat = '',i8,'';'')')nlat
write(iunit,'(''nvert = '',i8,'';'')')nvert
write(iunit,'(''nens = '',i8,'';'')')ens_size
write(iunit,'(''interptest = [ ... '')')

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

         write(iunit,*) field(ilon,jlat,kvert,:)
         if (any(ios_out /= 0)) then

            nfailed    = nfailed + 1
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

write(iunit,'(''];'')')
write(iunit,'(''datmat = reshape(interptest,nvert,nlat,nlon,nens);'')')
write(iunit,'(''datmat = permute(datmat,[4,1,2,3]);'')')
write(iunit,'(''datmat(datmat == missingvals) = NaN;'')')
call close_file(iunit)

if ( do_output() ) then
   write(*,'(A)')     '-------------------------------------------------------------'
   write(*,'(A,I10)') 'total  interpolations : ', nlon*nlat*nvert
   write(*,'(A,I10)') 'failed interpolations : ', nfailed
   write(*,'(A)')     '-------------------------------------------------------------'
endif

call count_error_codes(all_ios_out(1:nfailed,:))

! Write out the netCDF file for easy exploration.

call DATE_AND_TIME(crdate,crtime,crzone,values)
write(string1,'(''YYYY MM DD HH MM SS = '',i4,5(1x,i2.2))') &
                  values(1), values(2), values(3), values(5), values(6), values(7)

ncid = nc_create_file(ncfilename,'test_interpolate_threed_sphere')

call nc_check( nf90_put_att(ncid, NF90_GLOBAL, 'creation_date' ,trim(string1) ), &
                  routine, 'creation put '//trim(ncfilename))

! Define dimensions

call nc_check(nf90_def_dim(ncid=ncid, name='lon', len=nlon, &
        dimid = nlonDimID),routine, 'nlon def_dim '//trim(ncfilename))

call nc_check(nf90_def_dim(ncid=ncid, name='lat', len=nlat, &
        dimid = nlatDimID),routine, 'nlat def_dim '//trim(ncfilename))

call nc_check(nf90_def_dim(ncid=ncid, name='vert', len=nvert, &
        dimid = nvertDimID),routine, 'nvert def_dim '//trim(ncfilename))

! Define variables

call nc_check(nf90_def_var(ncid=ncid, name='lon', xtype=nf90_double, &
        dimids=nlonDimID, varid=lonVarID), routine, &
                 'lon def_var '//trim(ncfilename))
call nc_check(nf90_put_att(ncid, lonVarID, 'range', interp_test_lonrange), &
           routine, 'put_att lonrange '//trim(ncfilename))
call nc_check(nf90_put_att(ncid, lonVarID, 'cartesian_axis', 'X'),   &
           routine, 'lon cartesian_axis '//trim(ncfilename))

call nc_check(nf90_def_var(ncid=ncid, name='lat', xtype=nf90_double, &
        dimids=nlatDimID, varid=latVarID), routine, &
                 'lat def_var '//trim(ncfilename))
call nc_check(nf90_put_att(ncid, latVarID, 'range', interp_test_latrange), &
           routine, 'put_att latrange '//trim(ncfilename))
call nc_check(nf90_put_att(ncid, latVarID, 'cartesian_axis', 'Y'),   &
           routine, 'lat cartesian_axis '//trim(ncfilename))

call nc_check(nf90_def_var(ncid=ncid, name='vert', xtype=nf90_double, &
        dimids=nvertDimID, varid=vertVarID), routine, &
                 'vert def_var '//trim(ncfilename))
call nc_check(nf90_put_att(ncid, vertVarID, 'range', interp_test_vertcoord), &
           routine, 'put_att vertrange '//trim(ncfilename))
call nc_check(nf90_put_att(ncid, vertVarID, 'cartesian_axis', 'Z'),   &
           routine, 'vert cartesian_axis '//trim(ncfilename))

! loop over ensemble members
do imem = 1, ens_size
   if ( ens_size > 1) then
      write(field_name,'(A,I2)') "field_",imem
   else
      field_name = "field"
   endif
   call nc_check(nf90_def_var(ncid=ncid, name=field_name, xtype=nf90_double, &
           dimids=(/ nlonDimID, nlatDimID, nvertDimID /), varid=VarID(imem)), routine, &
                    'field def_var '//trim(ncfilename))
   call nc_check(nf90_put_att(ncid, VarID(imem), 'long_name', quantity_string), &
              routine, 'put_att field long_name '//trim(ncfilename))
   call nc_check(nf90_put_att(ncid, VarID(imem), '_FillValue', MISSING_R8), &
              routine, 'put_att field FillValue '//trim(ncfilename))
   call nc_check(nf90_put_att(ncid, VarID(imem), 'missing_value', MISSING_R8), &
              routine, 'put_att field missing_value '//trim(ncfilename))
   call nc_check(nf90_put_att(ncid, VarID(imem), 'interp_test_vertcoord', interp_test_vertcoord ), &
              routine, 'put_att field interp_test_vertcoord '//trim(ncfilename))
enddo

! Leave define mode so we can fill the variables.
call nc_end_define_mode(ncid, routine)

! Fill the variables
call nc_check(nf90_put_var(ncid, lonVarID, lon), &
              routine,'lon put_var '//trim(ncfilename))
call nc_check(nf90_put_var(ncid, latVarID, lat), &
              routine,'lat put_var '//trim(ncfilename))
call nc_check(nf90_put_var(ncid, vertVarID, vert), &
              routine,'vert put_var '//trim(ncfilename))

do imem = 1, ens_size
   call nc_check(nf90_put_var(ncid, VarID(imem), field(:,:,:,imem)), &
                 routine,'field put_var '//trim(ncfilename))
enddo

! tidy up
call nc_close_file(ncid, routine)

deallocate(lon, lat, vert, field)
deallocate(all_ios_out)

test_interpolate_range = nfailed

end function test_interpolate_range


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

