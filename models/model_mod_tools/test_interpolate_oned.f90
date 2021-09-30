! DART software - Copyright UCAR. This open source software is provided
! by UCAR, "as is", without charge, subject to all terms of use at
! http://www.image.ucar.edu/DAReS/DART/DART_download
!

module test_interpolate_mod

!-------------------------------------------------------------------------------
! interpolation test routines for oned locations.
!-------------------------------------------------------------------------------

use             types_mod, only : r8, i8, MISSING_R8

use         utilities_mod, only : error_handler, E_MSG, E_ERR, &
                                  E_MSG, open_file, close_file, do_output

use  netcdf_utilities_mod, only : nc_check, nc_create_file, nc_close_file, &
                                  nc_end_define_mode

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

! for messages
character(len=512) :: string1, string2

contains

!-------------------------------------------------------------------------------
!> Interpolate over a range of x values.
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

! The arguments to this function must be the same across all types of location dimensions.
! Consequently, there are some unused variables. Could make them optional, but this is only
! a test.

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

real(r8), allocatable :: X(:)
real(r8), allocatable :: field(:,:)
integer,  allocatable :: all_ios_out(:,:)
integer               :: nx
integer               :: i, nfailed
character(len=128)    :: ncfilename, txtfilename

character(len=8)      :: crdate      ! needed by F90 DATE_AND_TIME intrinsic
character(len=10)     :: crtime      ! needed by F90 DATE_AND_TIME intrinsic
character(len=5)      :: crzone      ! needed by F90 DATE_AND_TIME intrinsic
integer, dimension(8) :: values      ! needed by F90 DATE_AND_TIME intrinsic

integer :: ncid, nxDimID
integer :: VarID(ens_size), XVarID

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

iunit = open_file(trim(txtfilename), action='write')
write(iunit,'(''missingvals = '',f12.4,'';'')')MISSING_R8
write(iunit,'(''nx = '',i8,'';'')')nx
write(iunit,'(''nens = '',i8,'';'')')ens_size
write(iunit,'(''interptest = [ ... '')')

allocate(X(nx), field(nx,ens_size))
allocate(all_ios_out(nx,ens_size))

all_ios_out = 0 ! assume successful interpolation for every grid location, all members.
nfailed = 0

do i = 1, nx
   X(i) = interp_test_xrange(1) + real(i-1,r8) * interp_test_dx
   loc  = set_location(X(i))

   call model_interpolate(ens_handle, ens_size, loc, quantity_index, field(i,:), ios_out)

   call verify_consistent_istatus(ens_size, field(i,:), ios_out)

   write(iunit,*) field(i,:)
   if (any(ios_out /= 0)) then

      nfailed = nfailed + 1
      ! don't really care which location was causing the failure
      all_ios_out(nfailed,:) = ios_out

      if (verbose) then
         write(string1,*) 'interpolation return code was', ios_out
         write(string2,'(''i,X'',1(1x,i6),1(1x,f14.6))') i,X(i)
         call error_handler(E_MSG, routine, string1, text2=string2)
      endif
   endif
enddo

write(iunit,'(''];'')')
write(iunit,'(''datmat = reshape(interptest,nx,nens);'')')
write(iunit,'(''datmat = permute(datmat,[2,1]);'')')
write(iunit,'(''datmat(datmat == missingvals) = NaN;'')')
call close_file(iunit)

if ( do_output() ) then
   write(*,'(A)')     '-------------------------------------------------------------'
   write(*,'(A,I10)') 'total  interpolations : ', nx
   write(*,'(A,I10)') 'failed interpolations : ', nfailed
   write(*,'(A)')     '-------------------------------------------------------------'
endif

call count_error_codes(all_ios_out(1:nfailed,:))

! Write out the netCDF file for easy exploration.

call DATE_AND_TIME(crdate,crtime,crzone,values)
write(string1,'(''YYYY MM DD HH MM SS = '',i4,5(1x,i2.2))') &
                  values(1), values(2), values(3), values(5), values(6), values(7)

ncid = nc_create_file(ncfilename,'test_interpolate_oned')

call nc_check( nf90_put_att(ncid, NF90_GLOBAL, 'creation_date' ,trim(string1) ), &
                  routine, 'creation put '//trim(ncfilename))

! Define dimensions

call nc_check(nf90_def_dim(ncid=ncid, name='X', len=nx, &
        dimid = nxDimID),routine, 'nx def_dim '//trim(ncfilename))

! Define variables

call nc_check(nf90_def_var(ncid=ncid, name='X', xtype=nf90_double, &
        dimids=nxDimID, varid=XVarID), routine, &
                 'X def_var '//trim(ncfilename))
call nc_check(nf90_put_att(ncid, XVarID, 'range', interp_test_xrange), &
           routine, 'put_att xrange '//trim(ncfilename))
call nc_check(nf90_put_att(ncid, XVarID, 'cartesian_axis', 'X'),   &
           routine, 'X cartesian_axis '//trim(ncfilename))

! loop over ensemble members
do imem = 1, ens_size
   if ( ens_size > 1) then
      write(field_name,'(A,I2)') "field_",imem
   else
      field_name = "field"
   endif
   call nc_check(nf90_def_var(ncid=ncid, name=field_name, xtype=nf90_double, &
           dimids=(/ nxDimID /), varid=VarID(imem)), routine, &
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

do imem = 1, ens_size
   call nc_check(nf90_put_var(ncid, VarID(imem), field(:,imem)), &
                 routine,'field put_var '//trim(ncfilename))
enddo

! tidy up
call nc_close_file(ncid, routine)

deallocate(X, field)
deallocate(all_ios_out)

test_interpolate_range = nfailed

end function test_interpolate_range


!-------------------------------------------------------------------------------
! convert an array of reals into a location for this type

function setup_location(vals, vertcoord_string)

real(r8),         intent(in) :: vals(:)
character(len=*), intent(in) :: vertcoord_string
type(location_type)          :: setup_location

setup_location = set_location(vals(1))

end function setup_location

!-------------------------------------------------------------------------------
! End of test_interpolate_mod
!-------------------------------------------------------------------------------

end module test_interpolate_mod

