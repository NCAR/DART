! DART software - Copyright UCAR. This open source software is provided
! by UCAR, "as is", without charge, subject to all terms of use at
! http://www.image.ucar.edu/DAReS/DART/DART_download
!
! $Id$

module test_roms_interpolate_mod

!-------------------------------------------------------------------------------
! interpolation test routines for threed sphere locations.
!-------------------------------------------------------------------------------

use             types_mod, only : r8, i8, MISSING_R8, metadatalength

use         utilities_mod, only : register_module, error_handler, E_MSG, E_ERR, &
                                  initialize_utilities, finalize_utilities,     &
                                  find_namelist_in_file, check_namelist_read,   &
                                  E_MSG, open_file, close_file, do_output

use  netcdf_utilities_mod, only : nc_check

use          location_mod, only : location_type, set_location, write_location,  &
                                  get_dist, VERTISUNDEF, VERTISSURFACE,         &
                                  VERTISLEVEL, VERTISPRESSURE, VERTISHEIGHT,    &
                                  VERTISSCALEHEIGHT, get_location

use          obs_kind_mod, only : get_name_for_quantity

use  ensemble_manager_mod, only : ensemble_type

use             model_mod, only : model_interpolate, get_location_from_ijk

use netcdf

implicit none

public :: test_interpolate_range, test_interpolate_single

! version controlled file description for error handling, do not edit
character(len=256), parameter :: source   = &
   "$URL$"
character(len=32 ), parameter :: revision = "$Revision$"
character(len=128), parameter :: revdate  = "$Date$"

contains

!-------------------------------------------------------------------------------
! Do a interpolation on a range of dj, di, dk values.  Returns the
! number of failures.
!-------------------------------------------------------------------------------
function test_interpolate_range( ens_handle,            &
                                 interp_test_di,        &
                                 interp_test_dj,        &
                                 interp_test_dk,        &
                                 interp_test_vertcoord, &
                                 interp_test_irange,    &
                                 interp_test_jrange,    &
                                 interp_test_krange,    &
                                 mykindindex,           &
                                 verbose )

type(ensemble_type)   , intent(inout) :: ens_handle
real(r8)              , intent(in)    :: interp_test_di
real(r8)              , intent(in)    :: interp_test_dj
real(r8)              , intent(in)    :: interp_test_dk
character(len=*)      , intent(in)    :: interp_test_vertcoord
real(r8), dimension(2), intent(in)    :: interp_test_irange
real(r8), dimension(2), intent(in)    :: interp_test_jrange
real(r8), dimension(2), intent(in)    :: interp_test_krange
integer               , intent(in)    :: mykindindex
logical               , intent(in)    :: verbose

! function to exercise the model_mod:model_interpolate() function
! This will result in a netCDF file with all salient metadata
integer :: test_interpolate_range

character(len=metadatalength) :: kind_of_interest

! Local variables

real(r8), allocatable :: di(:), dj(:), dk(:)
real(r8), allocatable :: fieldJ(:,:,:)
real(r8), allocatable :: fieldI(:,:,:)
real(r8), allocatable :: fieldK(:,:,:)
integer :: NI, NJ, NK
integer :: i, j, k, nfailed
character(len=128) :: ncfilename, txtfilename(3)

character(len=8)      :: crdate      ! needed by F90 DATE_AND_TIME intrinsic
character(len=10)     :: crtime      ! needed by F90 DATE_AND_TIME intrinsic
character(len=5)      :: crzone      ! needed by F90 DATE_AND_TIME intrinsic
integer, dimension(8) :: values      ! needed by F90 DATE_AND_TIME intrinsic

integer :: ncid, NIDimID, NJDimID, NKDimID
integer :: VarID(3), IVarID, JVarID, KVarID

character(len=256) :: output_file = 'check_me'

! for message strings
character(len=512) :: string1, string2

character(len=32)  :: field_name(3)
type(location_type) :: loc
integer :: iunit(3), ios_out, ifield, vertcoord
integer, allocatable :: all_ios_out(:)

real(r8) :: my_loc(3)

test_interpolate_range = 0

if ((interp_test_di < 0.0_r8) .or. (interp_test_dj < 0.0_r8)) then
   if ( do_output() ) then
      write(*,*)'Skipping the rigorous interpolation test because one of'
      write(*,*)'interp_test_di,interp_test_dj are < 0.0'
      write(*,*)'interp_test_di  = ',interp_test_di
      write(*,*)'interp_test_dj  = ',interp_test_dj
      write(*,*)'interp_test_dk = ',interp_test_dk
  endif
   return
endif

vertcoord = get_location_index(interp_test_vertcoord)

write( ncfilename,'(a,a)')trim(output_file),'_interptest.nc'
write(txtfilename(1),'(a,a,a)')trim(output_file),'_lon_interptest.m'
write(txtfilename(2),'(a,a,a)')trim(output_file),'_lat_interptest.m'
write(txtfilename(3),'(a,a,a)')trim(output_file),'_vrt_interptest.m'

! round down to avoid exceeding the specified range
NJ  = aint(( interp_test_irange(2) -  interp_test_irange(1))/interp_test_dj) + 1
NI  = aint(( interp_test_jrange(2) -  interp_test_jrange(1))/interp_test_di) + 1
NK  = aint(( interp_test_krange(2) -  interp_test_krange(1))/interp_test_dk) + 1

do i = 1,3
   iunit(i) = open_file(trim(txtfilename(i)), action='write')
   write(iunit(i),'(''missingvals = '',f12.4,'';'')')MISSING_R8
   write(iunit(i),'(''NI = '',i8,'';'')')NI
   write(iunit(i),'(''NJ = '',i8,'';'')')NJ
   write(iunit(i),'(''NK = '',i8,'';'')')NK
   write(iunit(i),'(''interptest = [ ... '')')
enddo

allocate(di(NI), dj(NJ), dk(NK))
allocate(fieldI(NI,NJ,NK))
allocate(fieldJ(NI,NJ,NK))
allocate(fieldK(NI,NJ,NK))
allocate(all_ios_out(NI*NJ*NK))

all_ios_out = 0 ! bad interpolation
nfailed = 0

do i = 1, NI
   di(i) = interp_test_irange(1) + real(i-1,r8) * interp_test_di
   do j = 1, NJ
      dj(j) = interp_test_jrange(1) + real(j-1,r8) * interp_test_dj
      do k = 1, NK
         dk(k) = interp_test_krange(1) + real(k-1,r8) * interp_test_dk
         
         ios_out = get_location_from_ijk(di(i), dj(j), dk(k), mykindindex, loc)

         my_loc = get_location(loc) 

         fieldI(i,j,k) = my_loc(1)
         fieldJ(i,j,k) = my_loc(2)
         fieldK(i,j,k) = my_loc(3)

         write(iunit(1),*) fieldI(i,j,k)
         write(iunit(2),*) fieldJ(i,j,k)
         write(iunit(3),*) fieldK(i,j,k)

         if ( ios_out /= 0 ) then
           if (verbose) then
              write(string2,'(''i,j,k,di,dj,dk'',3(1x,i6),3(1x,f14.6))') &
                          i,j,k,di(i),dj(j),dk(k)
              write(string1,*) 'interpolation return code was', ios_out
              call error_handler(E_MSG,'test_interpolate_range',string1,source,revision,revdate,text2=string2)
           endif
           nfailed = nfailed + 1
           all_ios_out(nfailed) = ios_out
         endif

      enddo
   end do
end do

do i=1,3
   write(iunit(i),'(''];'')')
   write(iunit(i),'(''datmat = reshape(interptest,NK,NJ,NI);'')')
   write(iunit(i),'(''datmat = permute(datmat,[4,1,2]);'')')
   write(iunit(i),'(''datmat(datmat == missingvals) = NaN;'')')
   call close_file(iunit(i))
enddo

if ( do_output() ) then
   write(*,*) 'total interpolations  : ', NI*NJ*NK
   write(*,*) 'failed interpolations : ', nfailed
endif

call count_error_codes(all_ios_out, nfailed)


! Write out the netCDF file for easy exploration.

call DATE_AND_TIME(crdate,crtime,crzone,values)
write(string1,'(''YYYY MM DD HH MM SS = '',i4,5(1x,i2.2))') &
                  values(1), values(2), values(3), values(5), values(6), values(7)

call nc_check( nf90_create(path=trim(ncfilename), cmode=NF90_clobber, ncid=ncid), &
                  'test_interpolate_range', 'open '//trim(ncfilename))
call nc_check( nf90_put_att(ncid, NF90_GLOBAL, 'creation_date' ,trim(string1) ), &
                  'test_interpolate_range', 'creation put '//trim(ncfilename))

! Define dimensions

call nc_check(nf90_def_dim(ncid=ncid, name='di', len=NI, &
        dimid = NIDimID),'test_interpolate_range', 'NI def_dim '//trim(ncfilename))

call nc_check(nf90_def_dim(ncid=ncid, name='dj', len=NJ, &
        dimid = NJDimID),'test_interpolate_range', 'NJ def_dim '//trim(ncfilename))

call nc_check(nf90_def_dim(ncid=ncid, name='dk', len=NK, &
        dimid = NKDimID),'test_interpolate_range', 'NK def_dim '//trim(ncfilename))

! Define variables

call nc_check(nf90_def_var(ncid=ncid, name='di', xtype=nf90_double, &
        dimids=NIDimID, varid=IVarID), 'test_interpolate_range', &
                 'di def_var '//trim(ncfilename))
call nc_check(nf90_put_att(ncid, IVarID, 'range', interp_test_jrange), &
           'test_interpolate_range', 'put_att lonrange '//trim(ncfilename))
call nc_check(nf90_put_att(ncid, IVarID, 'cartesian_axis', 'X'),   &
           'test_interpolate_range', 'di cartesian_axis '//trim(ncfilename))


call nc_check(nf90_def_var(ncid=ncid, name='dj', xtype=nf90_double, &
        dimids=NJDimID, varid=JVarID), 'test_interpolate_range', &
                 'dj def_var '//trim(ncfilename))
call nc_check(nf90_put_att(ncid, JVarID, 'range', interp_test_krange), &
           'test_interpolate_range', 'put_att latrange '//trim(ncfilename))
call nc_check(nf90_put_att(ncid, JVarID, 'cartesian_axis', 'Y'),   &
           'test_interpolate_range', 'dj cartesian_axis '//trim(ncfilename))

call nc_check(nf90_def_var(ncid=ncid, name='dk', xtype=nf90_double, &
        dimids=NKDimID, varid=KVarID), 'test_interpolate_range', &
                 'dk def_var '//trim(ncfilename))
call nc_check(nf90_put_att(ncid, KVarID, 'range', interp_test_vertcoord), &
           'test_interpolate_range', 'put_att vertrange '//trim(ncfilename))
call nc_check(nf90_put_att(ncid, KVarID, 'cartesian_axis', 'Z'),   &
           'test_interpolate_range', 'dk cartesian_axis '//trim(ncfilename))

! loop over ensemble members
field_name(1) = "field_lon"
field_name(2) = "field_lat"
field_name(3) = "field_vrt"
do ifield = 1, 3
   call nc_check(nf90_def_var(ncid=ncid, name=field_name(ifield), xtype=nf90_double, &
           dimids=(/ NIDimID, NJDimID, NKDimID /), varid=VarID(ifield)), 'test_interpolate_range', &
                    'field def_var '//trim(ncfilename))
   kind_of_interest = get_name_for_quantity(mykindindex)
   call nc_check(nf90_put_att(ncid, VarID(ifield), 'long_name', kind_of_interest), &
              'test_interpolate_range', 'put_att field long_name '//trim(ncfilename))
   call nc_check(nf90_put_att(ncid, VarID(ifield), '_FillValue', MISSING_R8), &
              'test_interpolate_range', 'put_att field FillValue '//trim(ncfilename))
   call nc_check(nf90_put_att(ncid, VarID(ifield), 'missing_value', MISSING_R8), &
              'test_interpolate_range', 'put_att field missing_value '//trim(ncfilename))
   call nc_check(nf90_put_att(ncid, VarID(ifield), 'interp_test_vertcoord', interp_test_vertcoord ), &
              'test_interpolate_range', 'put_att field interp_test_vertcoord '//trim(ncfilename))
enddo

! Leave define mode so we can fill the variables.
call nc_check(nf90_enddef(ncid), &
              'test_interpolate_range','field enddef '//trim(ncfilename))

! Fill the variables
call nc_check(nf90_put_var(ncid, IVarID, di), &
              'test_interpolate_range','di put_var '//trim(ncfilename))
call nc_check(nf90_put_var(ncid, JVarID, dj), &
              'test_interpolate_range','dj put_var '//trim(ncfilename))
call nc_check(nf90_put_var(ncid, KVarID, dk), &
              'test_interpolate_range','dk put_var '//trim(ncfilename))

call nc_check(nf90_put_var(ncid, VarID(1), fieldI(:,:,:)), &
              'test_interpolate_range','fieldI put_var '//trim(ncfilename))

call nc_check(nf90_put_var(ncid, VarID(2), fieldJ(:,:,:)), &
              'test_interpolate_range','fieldJ put_var '//trim(ncfilename))

call nc_check(nf90_put_var(ncid, VarID(3), fieldK(:,:,:)), &
              'test_interpolate_range','fieldvert put_var '//trim(ncfilename))

! tidy up
call nc_check(nf90_close(ncid), &
             'test_interpolate_range','close '//trim(ncfilename))

deallocate(di, dj, dk, fieldI, fieldJ, fieldK)

test_interpolate_range = nfailed

end function test_interpolate_range


!-------------------------------------------------------------------------------
! Do a single interpolation on a given location and kind.  Returns the
! interpolated values and ios_out. Returns the number of ensemble members that
! passed
!-------------------------------------------------------------------------------
function test_interpolate_single( ens_handle,       &
                                  vertcoord_string, &
                                  ival,             &
                                  jval,             &
                                  kval,             &
                                  mykindindex,      &
                                  value,            &
                                  ios_out) result(num_passed)

type(ensemble_type)   , intent(inout) :: ens_handle
character(len=*)      , intent(in)    :: vertcoord_string
real(r8)              , intent(in)    :: ival
real(r8)              , intent(in)    :: jval
real(r8)              , intent(in)    :: kval
integer               , intent(in)    :: mykindindex
real(r8)              , intent(out)   :: value
integer               , intent(out)   :: ios_out
integer :: num_passed

type(location_type) :: loc
integer :: vertcoord
real(r8) :: my_loc(3)

num_passed = 0

vertcoord = get_location_index(vertcoord_string)
ios_out   = get_location_from_ijk(ival, jval, kval, mykindindex, loc)
my_loc    = get_location(loc) 

! where is the code that actually does the interpolate??  nsc.
value = MISSING_R8

if (ios_out == 0 ) then
   if (do_output() ) &
      write(*,*) 'test_interpolate_single SUCCESS lat,lon,lev', my_loc(1),my_loc(2),my_loc(3)
   num_passed = num_passed + 1
else
   if (do_output() ) &
      write(*,*) 'test_interpolate_single ERROR with error code : ', ios_out
endif

end function test_interpolate_single

!-------------------------------------------------------------------------------
! Count the number of different error codes and output the results.  This
! is just a helper function for test_interpolate_range. Only sums error codes
! for the first ensemble member
!-------------------------------------------------------------------------------
subroutine count_error_codes(error_codes, num_failed)

integer, intent(in) :: error_codes(:)
integer, intent(in) :: num_failed

integer :: i, count_errors, results

count_errors = 1

i = 1
do while (count_errors < num_failed)
   results = count(error_codes(:) == i)
   if (results /= 0) then
      if ( do_output() ) &
         write(*,'(i10, a, i3)') results + 1, " failed with ios_out ", i
      count_errors = count_errors + results
   endif
   i = i+1
enddo

end subroutine count_error_codes

!-------------------------------------------------------------------------------
! need to convert the character string for the test vertical coordinate into 
! the corresponding dart index.
!-------------------------------------------------------------------------------
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

end module test_roms_interpolate_mod

! <next few lines under version control, do not edit>
! $URL$
! $Id$
! $Revision$
! $Date$
