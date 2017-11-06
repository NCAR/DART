! DART software - Copyright UCAR. This open source software is provided
! by UCAR, "as is", without charge, subject to all terms of use at
! http://www.image.ucar.edu/DAReS/DART/DART_download
!
! $Id$

module test_interpolate_mod

!-------------------------------------------------------------------------------
! interpolation test routines for oned locations.
!-------------------------------------------------------------------------------

use             types_mod, only : r8, i8, missing_r8, metadatalength

use         utilities_mod, only : register_module, error_handler, E_MSG, E_ERR, &
                                  initialize_utilities, finalize_utilities,     &
                                  find_namelist_in_file, check_namelist_read,   &
                                  nc_check, E_MSG, open_file, close_file, do_output

use          location_mod, only : location_type, set_location, write_location,  &
                                  get_dist, get_location, LocationDims

use          obs_kind_mod, only : get_name_for_quantity, get_index_for_quantity

use  ensemble_manager_mod, only : ensemble_type

use             model_mod, only : get_model_size, &
                                  get_state_meta_data, &
                                  model_interpolate

use netcdf

implicit none

public :: test_interpolate_single, &
          test_interpolate_range, &
          find_closest_gridpoint

! version controlled file description for error handling, do not edit
character(len=256), parameter :: source   = &
   "$URL$"
character(len=32 ), parameter :: revision = "$Revision$"
character(len=128), parameter :: revdate  = "$Date$"

! for messages
character(len=512) :: string1, string2, string3

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
nx  = aint(( interp_test_xrange(2) -  interp_test_xrange(1))/interp_test_dx) + 1

iunit = open_file(trim(txtfilename), action='write')
write(iunit,'(''missingvals = '',f12.4,'';'')')MISSING_R8
write(iunit,'(''nx = '',i8,'';'')')nx
write(iunit,'(''nens = '',i8,'';'')')ens_size
write(iunit,'(''interptest = [ ... '')')

allocate(X(nx), field(nx,ens_size))
allocate(all_ios_out(nx,ens_size))
nfailed = 0

do i = 1, nx
   X(i) = interp_test_xrange(1) + real(i-1,r8) * interp_test_dx
   loc  = set_location(X(i))
   call model_interpolate(ens_handle, ens_size, loc, mykindindex, field(i,:), ios_out)
   write(iunit,*) field(i,:)
   if (any(ios_out /= 0)) then
      if (verbose) then
         write(string1,*) 'interpolation return code was', ios_out
         write(string2,'(''i,X'',1(1x,i6),1(1x,f14.6))') i,X(i)
         call error_handler(E_MSG, routine, string1, &
                            source, revision, revdate, text2=string2)
      endif
      nfailed = nfailed + 1
      all_ios_out(nfailed,:) = ios_out
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

call count_error_codes(all_ios_out, nfailed)

! Write out the netCDF file for easy exploration.

call DATE_AND_TIME(crdate,crtime,crzone,values)
write(string1,'(''YYYY MM DD HH MM SS = '',i4,5(1x,i2.2))') &
                  values(1), values(2), values(3), values(5), values(6), values(7)

call nc_check( nf90_create(path=trim(ncfilename), cmode=NF90_clobber, ncid=ncid), &
                  routine, 'open '//trim(ncfilename))
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
call nc_check(nf90_enddef(ncid), &
              routine,'field enddef '//trim(ncfilename))

! Fill the variables
call nc_check(nf90_put_var(ncid, XVarID, X), &
              routine,'X put_var '//trim(ncfilename))

do imem = 1, ens_size
   call nc_check(nf90_put_var(ncid, VarID(imem), field(:,imem)), &
                 routine,'field put_var '//trim(ncfilename))
enddo

! tidy up
call nc_check(nf90_close(ncid), routine,'close '//trim(ncfilename))

deallocate(X, field)

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

type(ensemble_type)   , intent(inout) :: ens_handle
integer               , intent(in)    :: ens_size
character(len=*)      , intent(in)    :: vertcoord_string
real(r8)              , intent(in)    :: xval
real(r8)              , intent(in)    :: yval
real(r8)              , intent(in)    :: zval
character(len=*)      , intent(in)    :: quantity_string
real(r8)              , intent(out)   :: interp_vals(ens_size)
integer               , intent(out)   :: ios_out(ens_size)

integer :: test_interpolate_single

type(location_type) :: loc
integer :: imem, num_passed
character(len=128) :: my_location
integer :: quantity_index

quantity_index = get_index_for_quantity(quantity_string)

loc = set_location(xval)

if ( do_output() ) then
   call write_location(0, loc, charstring=my_location)
   write(*,'(A)') ''
   write(*,'(A)') '-------------------------------------------------------------'
   write(*,'("interpolating at ",A)') trim(my_location)//' for "'//trim(quantity_string)//'"'
   write(*,'(A)') '-------------------------------------------------------------'
   write(*,'(A)') ''
endif

call model_interpolate(ens_handle, ens_size, loc, quantity_index, interp_vals, ios_out)

num_passed = 0
do imem = 1, ens_size
   if (ios_out(imem) == 0 ) then
      if (do_output()) then
         write(string1,*)'model_interpolate SUCCESS with value    :: ', interp_vals(imem)
         write(*,'(A,I5,A,A)')'member ',imem,',',trim(string1)
         num_passed = num_passed + 1
      endif
   else
      if (do_output()) then
         write(string1,*)'model_interpolate ERROR with error code :: ', ios_out(imem)
         write(*,'(A,I5,A,A)')'member ',imem,',',trim(string1)
      endif
   endif
enddo

if ( do_output() ) write(*,'(A)') ''

test_interpolate_single = num_passed

end function test_interpolate_single


!-------------------------------------------------------------------------------
!> Count the number of different error codes and output the results.
!> This is just a helper function for test_interpolate_range.
!> Only sums error codes for the first ensemble member.

subroutine count_error_codes(error_codes, num_failed)

integer, intent(in) :: error_codes(:,:)
integer, intent(in) :: num_failed

integer :: i, count_errors, results

count_errors = 1

i = 1
do while (count_errors < num_failed)
   results = count(error_codes(:,1) == i)
   if (results /= 0) then
      if ( do_output() ) &
         write(*,'(i10, a, i3)') results + 1, " failed with ios_out ", i
      count_errors = count_errors + results
   endif
   i = i+1
enddo

end subroutine count_error_codes


!-----------------------------------------------------------------------
!> Expensive exhaustive search to find the indices into the
!> state vector of a particular lon/lat/vert. At present, only for a
!> single variable - could be extended to identify the closest location
!> for every variable in each domain. This could help ensure grid
!> staggering is being handled correctly.

subroutine find_closest_gridpoint(loc_of_interest, quantity_string)

real(r8),         intent(in) :: loc_of_interest(:)
character(len=*), intent(in) :: quantity_string

character(len=*), parameter :: routine = 'find_closest_gridpoint'

type(location_type)   :: loc0, loc1
integer(i8)           :: i
integer               :: quantity_index, var_type
real(r8)              :: closest, x
logical               :: matched
real(r8), allocatable :: thisdist(:)
real(r8),   parameter :: FARAWAY = huge(r8)
character(len=metadatalength) :: myquantity

!>@todo there should be arrays of length state_structure_mod:get_num_variables(domid)
!>      get_num_domains(), get_num_variables() ...

allocate( thisdist(get_model_size()) )
thisdist  = FARAWAY
matched   = .false.

! Trying to support the ability to specify matching a particular QUANTITY.
! With staggered grids, the closest gridpoint might not be of the quantity
! of interest.

quantity_index = get_index_for_quantity(quantity_string)
x = loc_of_interest(1)
loc0 = set_location(x)

write(string1,*)'Checking for the indices into the state vector that are close to'
call write_location(0, loc0, charstring=string2)
write(string3,*)'for "',trim(quantity_string),'"'
call error_handler(E_MSG,routine,string1,text2=string2,text3=string3)

! Since there can be/will be multiple variables with
! identical distances, we will just cruise once through
! the array and come back to find all the 'identical' values.

DISTANCE : do i = 1,get_model_size()

   call get_state_meta_data(i, loc1, var_type)

   if (var_type .ne. quantity_index) cycle DISTANCE

   thisdist(i) = get_dist(loc1, loc0)
   matched     = .true.

enddo DISTANCE

if (.not. matched) then
   write(string1,*)'No state vector elements of type "'//trim(quantity_string)//'"'
   call error_handler(E_MSG, routine, string1)
   deallocate( thisdist )
   return
endif

closest = minval(thisdist)

! Now that we know the distances ... report
! If more than one quantity has the same distance, report all.
! Be aware that if 'approximate_distance' is .true., everything
! in the box has a single location.

REPORT: do i = 1,get_model_size()

   if ( thisdist(i) == closest ) then
      call get_state_meta_data(i, loc1, var_type)
      myquantity = get_name_for_quantity(var_type)

      call write_location(0, loc1, charstring=string2)
      write(string1,'(A,I12,A)')trim(string2)//' is index ',i,' ('//trim(myquantity)//')'
      call error_handler(E_MSG, routine, string1)
   endif

enddo REPORT

deallocate( thisdist )

end subroutine find_closest_gridpoint

!-------------------------------------------------------------------------------
! End of test_interpolate_mod
!-------------------------------------------------------------------------------

end module test_interpolate_mod

! <next few lines under version control, do not edit>
! $URL$
! $Id$
! $Revision$
! $Date$
