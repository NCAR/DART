! DART software - Copyright UCAR. This open source software is provided
! by UCAR, "as is", without charge, subject to all terms of use at
! http://www.image.ucar.edu/DAReS/DART/DART_download
!
! $Id$

module test_interpolate_mod

!-------------------------------------------------------------------------------
! interpolation test routines for threed cartesian locations.
!-------------------------------------------------------------------------------

use             types_mod, only : r8, i8, MISSING_R8, metadatalength

use         utilities_mod, only : register_module, error_handler, E_MSG, E_ERR, &
                                  initialize_utilities, finalize_utilities,     &
                                  find_namelist_in_file, check_namelist_read,   &
                                  E_MSG, open_file, close_file, do_output

use  netcdf_utilities_mod, only : nc_check, nc_create_file, nc_close_file, &
                                  nc_end_define_mode

use          location_mod, only : location_type, set_location, write_location,  &
                                  get_dist, get_location, LocationDims

use          obs_kind_mod, only : get_name_for_quantity, get_index_for_quantity

use  ensemble_manager_mod, only : ensemble_type

use model_check_utilities_mod, only : test_single_interpolation, &
                                      find_closest_gridpoint, &
                                      count_error_codes, &
                                      verify_consistent_istatus

use             model_mod, only : get_model_size, &
                                  get_state_meta_data, &
                                  model_interpolate

use netcdf

implicit none
private

public :: test_interpolate_single, &
          test_interpolate_range, &
          find_closest_gridpoint

! version controlled file description for error handling, do not edit
character(len=*), parameter :: source   = "$URL$"
character(len=*), parameter :: revision = "$Revision$"
character(len=*), parameter :: revdate  = "$Date$"

! for messages
character(len=512) :: string1, string2

contains

!-------------------------------------------------------------------------------
!> Interpolate over a range of x, y, and z values.
!> Returns the number of failures.
!> Exercises model_mod:model_interpolate().
!> This will result in a netCDF file with all salient metadata.

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
character(len=128)    :: ncfilename, txtfilename

character(len=8)      :: crdate      ! needed by F90 DATE_AND_TIME intrinsic
character(len=10)     :: crtime      ! needed by F90 DATE_AND_TIME intrinsic
character(len=5)      :: crzone      ! needed by F90 DATE_AND_TIME intrinsic
integer, dimension(8) :: values      ! needed by F90 DATE_AND_TIME intrinsic

integer :: ncid, nxDimID, nyDimID, nzDimID
integer :: VarID(ens_size), XVarID, YVarID, ZVarID

character(len=256)  :: output_file = 'check_me'
character(len=32)   :: field_name
type(location_type) :: loc
integer :: iunit, ios_out(ens_size), imem
integer :: quantity_index

test_interpolate_range = 0
quantity_index = get_index_for_quantity(quantity_string)

write( ncfilename,'(a,a)')trim(output_file),'_interptest.nc'
write(txtfilename,'(a,a)')trim(output_file),'_interptest.m'

! round down to avoid exceeding the specified range
nx = aint((interp_test_xrange(2) - interp_test_xrange(1))/interp_test_dx) + 1
ny = aint((interp_test_yrange(2) - interp_test_yrange(1))/interp_test_dy) + 1
nz = aint((interp_test_zrange(2) - interp_test_zrange(1))/interp_test_dz) + 1

iunit = open_file(trim(txtfilename), action='write')
write(iunit,'(''missingvals = '',f12.4,'';'')')MISSING_R8
write(iunit,'(''nx = '',i8,'';'')')nx
write(iunit,'(''ny = '',i8,'';'')')ny
write(iunit,'(''nz = '',i8,'';'')')nz
write(iunit,'(''nens = '',i8,'';'')')ens_size
write(iunit,'(''interptest = [ ... '')')

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

         write(iunit,*) field(i,j,k,:)
         if (any(ios_out(:) /= 0)) then

            nfailed    = nfailed + 1
            ! don't really care which location was causing the failure
            all_ios_out(nfailed,:) = ios_out

            if (verbose) then
               write(string1,*) 'interpolation return code was', ios_out
               write(string2,'(''i,j,k,X,Y,Z'',3(1x,i6),3(1x,f14.6))') i,j,k,X(i),Y(j),Z(k)
               call error_handler(E_MSG, routine, string1, &
                                  source, revision, revdate, text2=string2)
            endif
         endif
      enddo
   enddo
enddo

write(iunit,'(''];'')')
write(iunit,'(''datmat = reshape(interptest,nz,ny,nx,nens);'')')
write(iunit,'(''datmat = permute(datmat,[4,1,2,3]);'')')
write(iunit,'(''datmat(datmat == missingvals) = NaN;'')')
call close_file(iunit)

if ( do_output() ) then
   write(*,'(A)')     '-------------------------------------------------------------'
   write(*,'(A,I10)') 'total  interpolations : ', nx*ny*nz
   write(*,'(A,I10)') 'failed interpolations : ', nfailed
   write(*,'(A)')     '-------------------------------------------------------------'
endif

call count_error_codes(all_ios_out(1:nfailed,:))

! Write out the netCDF file for easy exploration.

call DATE_AND_TIME(crdate,crtime,crzone,values)
write(string1,'(''YYYY MM DD HH MM SS = '',i4,5(1x,i2.2))') &
                  values(1), values(2), values(3), values(5), values(6), values(7)

ncid = nc_create_file(ncfilename,'test_interpolate_threed_cartesian')

call nc_check( nf90_put_att(ncid, NF90_GLOBAL, 'creation_date' ,trim(string1) ), &
                  routine, 'creation put '//trim(ncfilename))

! Define dimensions

call nc_check(nf90_def_dim(ncid=ncid, name='X', len=nx, &
        dimid = nxDimID),routine, 'nx def_dim '//trim(ncfilename))

call nc_check(nf90_def_dim(ncid=ncid, name='Y', len=ny, &
        dimid = nyDimID),routine, 'ny def_dim '//trim(ncfilename))

call nc_check(nf90_def_dim(ncid=ncid, name='Z', len=nz, &
        dimid = nzDimID),routine, 'nz def_dim '//trim(ncfilename))

! Define variables

call nc_check(nf90_def_var(ncid=ncid, name='X', xtype=nf90_double, &
        dimids=nxDimID, varid=XVarID), routine, &
                 'X def_var '//trim(ncfilename))
call nc_check(nf90_put_att(ncid, XVarID, 'range', interp_test_xrange), &
           routine, 'put_att xrange '//trim(ncfilename))
call nc_check(nf90_put_att(ncid, XVarID, 'cartesian_axis', 'X'),   &
           routine, 'X cartesian_axis '//trim(ncfilename))

call nc_check(nf90_def_var(ncid=ncid, name='Y', xtype=nf90_double, &
        dimids=nyDimID, varid=YVarID), routine, &
                 'Y def_var '//trim(ncfilename))
call nc_check(nf90_put_att(ncid, YVarID, 'range', interp_test_yrange), &
           routine, 'put_att yrange '//trim(ncfilename))
call nc_check(nf90_put_att(ncid, YVarID, 'cartesian_axis', 'Y'),   &
           routine, 'Y cartesian_axis '//trim(ncfilename))

call nc_check(nf90_def_var(ncid=ncid, name='Z', xtype=nf90_double, &
        dimids=nzDimID, varid=ZVarID), routine, &
                 'Z def_var '//trim(ncfilename))
call nc_check(nf90_put_att(ncid, ZVarID, 'range', interp_test_zrange), &
           routine, 'put_att zrange '//trim(ncfilename))
call nc_check(nf90_put_att(ncid, ZVarID, 'cartesian_axis', 'Z'),   &
           routine, 'Z cartesian_axis '//trim(ncfilename))

! loop over ensemble members
do imem = 1, ens_size
   if ( ens_size > 1) then
      write(field_name,'(A,I2)') "field_",imem
   else
      field_name = "field"
   endif
   call nc_check(nf90_def_var(ncid=ncid, name=field_name, xtype=nf90_double, &
           dimids=(/ nxDimID, nyDimID, nzDimID /), varid=VarID(imem)), routine, &
                    'field def_var '//trim(ncfilename))
   call nc_check(nf90_put_att(ncid, VarID(imem), 'long_name', quantity_string), &
              routine, 'put_att field long_name '//trim(ncfilename))
   call nc_check(nf90_put_att(ncid, VarID(imem), '_FillValue', MISSING_R8), &
              routine, 'put_att field FillValue '//trim(ncfilename))
   call nc_check(nf90_put_att(ncid, VarID(imem), 'missing_value', MISSING_R8), &
              routine, 'put_att field missing_value '//trim(ncfilename))
enddo

! Leave define mode so we can fill the variables.
call nc_end_define_mode(ncid, routine)

! Fill the variables
call nc_check(nf90_put_var(ncid, XVarID, X), &
              routine,'X put_var '//trim(ncfilename))
call nc_check(nf90_put_var(ncid, YVarID, Y), &
              routine,'Y put_var '//trim(ncfilename))
call nc_check(nf90_put_var(ncid, ZVarID, Z), &
              routine,'Z put_var '//trim(ncfilename))

do imem = 1, ens_size
   call nc_check(nf90_put_var(ncid, VarID(imem), field(:,:,:,imem)), &
                 routine,'field put_var '//trim(ncfilename))
enddo

! tidy up
call nc_close_file(ncid, routine)

deallocate(X, Y, Z, field)
deallocate(all_ios_out)

test_interpolate_range = nfailed

end function test_interpolate_range


!-------------------------------------------------------------------------------
!> Do a single interpolation on a given location and kind.
!> Returns the interpolated values and ios_out.
!> Returns the number of ensemble members that passed.

function test_interpolate_single( ens_handle,       &
                                  ens_size,         &
                                  vertcoord_string, &
                                  xval,             &
                                  yval,             &
                                  zval,             &
                                  quantity_string,  &
                                  interp_vals,      &
                                  ios_out)

type(ensemble_type),intent(inout) :: ens_handle
integer            ,intent(in)    :: ens_size
character(len=*)   ,intent(in)    :: vertcoord_string
real(r8)           ,intent(in)    :: xval
real(r8)           ,intent(in)    :: yval
real(r8)           ,intent(in)    :: zval
character(len=*)   ,intent(in)    :: quantity_string
real(r8)           ,intent(out)   :: interp_vals(ens_size)
integer            ,intent(out)   :: ios_out(ens_size)

integer :: test_interpolate_single

type(location_type) :: location

location = set_location(xval, yval, zval)

test_interpolate_single = test_single_interpolation(ens_handle, ens_size, &
                               location, vertcoord_string, &
                               quantity_string, interp_vals, ios_out)

end function test_interpolate_single

!-------------------------------------------------------------------------------
! End of test_interpolate_mod
!-------------------------------------------------------------------------------

end module test_interpolate_mod

