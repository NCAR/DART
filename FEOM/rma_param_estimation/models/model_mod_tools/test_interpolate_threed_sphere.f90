! DART software - Copyright 2004 - 2013 UCAR. This open source software is
! provided by UCAR, "as is", without charge, subject to all terms of use at
! http://www.image.ucar.edu/DAReS/DART/DART_download
!
! $Id: model_mod_check.f90 6739 2014-01-15 20:44:54Z hkershaw $

module test_interpolate_mod

!-------------------------------------------------------------------------------
! interpolation test routines for threed sphere locations.
!-------------------------------------------------------------------------------

use             types_mod, only : r8, i8, missing_r8, metadatalength

use         utilities_mod, only : register_module, error_handler, E_MSG, E_ERR, &
                                  initialize_utilities, finalize_utilities,     &
                                  find_namelist_in_file, check_namelist_read,   &
                                  nc_check, E_MSG, open_file, close_file

use          location_mod, only : location_type, set_location, write_location,  &
                                  get_dist, VERTISUNDEF, VERTISSURFACE,         &
                                  VERTISLEVEL, VERTISPRESSURE, VERTISHEIGHT,    &
                                  VERTISSCALEHEIGHT

use          obs_kind_mod, only : get_raw_obs_kind_name

use  ensemble_manager_mod, only : ensemble_type

use             model_mod, only : model_interpolate

use netcdf

implicit none

public :: test_interpolate_range, test_interpolate_single

! version controlled file description for error handling, do not edit
character(len=256), parameter :: source   = &
   "$URL: https://proxy.subversion.ucar.edu/DAReS/DART/branches/rma_model_mod_check/models/template/model_mod_check.f90 $"
character(len=32 ), parameter :: revision = "$Revision: 6739 $"
character(len=128), parameter :: revdate  = "$Date: 2014-01-15 13:44:54 -0700 (Wed, 15 Jan 2014) $"

contains

!-------------------------------------------------------------------------------
! Do a interpolation on a range of lat, lon, vert values.  Returns the
! number of failures.
!-------------------------------------------------------------------------------
function test_interpolate_range( ens_handle,            &
                                 ens_size,              &
                                 interp_test_dlon,      &
                                 interp_test_dlat,      &
                                 interp_test_dvert,     &
                                 interp_test_vertcoord, &
                                 interp_test_lonrange,  &
                                 interp_test_latrange,  &
                                 interp_test_vertrange, &
                                 mykindindex,           &
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
integer               , intent(in)    :: mykindindex
logical               , intent(in)    :: verbose

! function to exercise the model_mod:model_interpolate() function
! This will result in a netCDF file with all salient metadata
integer :: test_interpolate_range

character(len=metadatalength) :: kind_of_interest

! Local variables

real(r8), allocatable :: lon(:), lat(:), vert(:)
real(r8), allocatable :: field(:,:,:,:)
integer :: nlon, nlat, nvert
integer :: ilon, jlat, kvert, nfailed
character(len=128) :: ncfilename, txtfilename

character(len=8)      :: crdate      ! needed by F90 DATE_AND_TIME intrinsic
character(len=10)     :: crtime      ! needed by F90 DATE_AND_TIME intrinsic
character(len=5)      :: crzone      ! needed by F90 DATE_AND_TIME intrinsic
integer, dimension(8) :: values      ! needed by F90 DATE_AND_TIME intrinsic

integer :: ncid, nlonDimID, nlatDimID, nvertDimID
integer :: VarID(ens_size), lonVarID, latVarID, vertVarID

character(len=256) :: output_file = 'check_me'

! for message strings
character(len=512) :: string1, string2

character(len=32)  :: field_name
type(location_type) :: loc
integer :: iunit, ios_out(ens_size), imem, vertcoord
integer, allocatable :: all_ios_out(:,:)

test_interpolate_range = 0

if ((interp_test_dlon < 0.0_r8) .or. (interp_test_dlat < 0.0_r8)) then
   write(*,*)'Skipping the rigorous interpolation test because one of'
   write(*,*)'interp_test_dlon,interp_test_dlat are < 0.0'
   write(*,*)'interp_test_dlon  = ',interp_test_dlon
   write(*,*)'interp_test_dlat  = ',interp_test_dlat
   write(*,*)'interp_test_dvert = ',interp_test_dvert
   return
endif

vertcoord = get_location_index(interp_test_vertcoord)

write( ncfilename,'(a,a)')trim(output_file),'_interptest.nc'
write(txtfilename,'(a,a)')trim(output_file),'_interptest.m'

! round down to avoid exceeding the specified range
nlat  = aint(( interp_test_latrange(2) -  interp_test_latrange(1))/interp_test_dlat) + 1
nlon  = aint(( interp_test_lonrange(2) -  interp_test_lonrange(1))/interp_test_dlon) + 1
nvert = aint((interp_test_vertrange(2) -  interp_test_vertrange(1))/interp_test_dvert) + 1

iunit = open_file(trim(txtfilename), action='write')
write(iunit,'(''missingvals = '',f12.4,'';'')')MISSING_R8
write(iunit,'(''nlon = '',i8,'';'')')nlon
write(iunit,'(''nlat = '',i8,'';'')')nlat
write(iunit,'(''nvert = '',i8,'';'')')nvert
write(iunit,'(''nens = '',i8,'';'')')ens_size
write(iunit,'(''interptest = [ ... '')')

allocate(lon(nlon), lat(nlat), vert(nvert), field(nlon,nlat,nvert,ens_size))
allocate(all_ios_out(nlon*nlat*nvert,ens_size))
nfailed = 0

do ilon = 1, nlon
   lon(ilon) = interp_test_lonrange(1) + real(ilon-1,r8) * interp_test_dlon
   do jlat = 1, nlat
      lat(jlat) = interp_test_latrange(1) + real(jlat-1,r8) * interp_test_dlat
      do kvert = 1, nvert
         vert(kvert) = interp_test_vertrange(1) + real(kvert-1,r8) * interp_test_dvert

         loc = set_location(lon(ilon), lat(jlat), vert(kvert), vertcoord)

         call model_interpolate(ens_handle, ens_size, loc, mykindindex, field(ilon,jlat,kvert,:), ios_out)

         write(iunit,*) field(ilon,jlat,kvert,:)
         if (any(ios_out /= 0)) then
           if (verbose) then
              write(string2,'(''ilon,jlat,kvert,lon,lat,vert'',3(1x,i6),3(1x,f14.6))') &
                          ilon,jlat,kvert,lon(ilon),lat(jlat),vert(kvert)
              write(string1,*) 'interpolation return code was', ios_out
              call error_handler(E_MSG,'test_interpolate_range',string1,source,revision,revdate,text2=string2)
           endif
           all_ios_out(nfailed,:) = ios_out
           nfailed = nfailed + 1
         endif

      enddo
   end do
end do

write(iunit,'(''];'')')
write(iunit,'(''datmat = reshape(interptest,nvert,nlat,nlon,nens);'')')
write(iunit,'(''datmat = permute(datmat,[4,1,2,3]);'')')
write(iunit,'(''datmat(datmat == missingvals) = NaN;'')')
call close_file(iunit)

write(*,*) 'total interpolations  : ', nlon*nlat*nvert
write(*,*) 'failed interpolations : ', nfailed

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

call nc_check(nf90_def_dim(ncid=ncid, name='lon', len=nlon, &
        dimid = nlonDimID),'test_interpolate_range', 'nlon def_dim '//trim(ncfilename))

call nc_check(nf90_def_dim(ncid=ncid, name='lat', len=nlat, &
        dimid = nlatDimID),'test_interpolate_range', 'nlat def_dim '//trim(ncfilename))

call nc_check(nf90_def_dim(ncid=ncid, name='vert', len=nvert, &
        dimid = nvertDimID),'test_interpolate_range', 'nvert def_dim '//trim(ncfilename))

! Define variables

call nc_check(nf90_def_var(ncid=ncid, name='lon', xtype=nf90_double, &
        dimids=nlonDimID, varid=lonVarID), 'test_interpolate_range', &
                 'lon def_var '//trim(ncfilename))
call nc_check(nf90_put_att(ncid, lonVarID, 'range', interp_test_lonrange), &
           'test_interpolate_range', 'put_att lonrange '//trim(ncfilename))
call nc_check(nf90_put_att(ncid, lonVarID, 'cartesian_axis', 'X'),   &
           'test_interpolate_range', 'lon cartesian_axis '//trim(ncfilename))


call nc_check(nf90_def_var(ncid=ncid, name='lat', xtype=nf90_double, &
        dimids=nlatDimID, varid=latVarID), 'test_interpolate_range', &
                 'lat def_var '//trim(ncfilename))
call nc_check(nf90_put_att(ncid, latVarID, 'range', interp_test_latrange), &
           'test_interpolate_range', 'put_att latrange '//trim(ncfilename))
call nc_check(nf90_put_att(ncid, latVarID, 'cartesian_axis', 'Y'),   &
           'test_interpolate_range', 'lat cartesian_axis '//trim(ncfilename))

call nc_check(nf90_def_var(ncid=ncid, name='vert', xtype=nf90_double, &
        dimids=nvertDimID, varid=vertVarID), 'test_interpolate_range', &
                 'vert def_var '//trim(ncfilename))
call nc_check(nf90_put_att(ncid, vertVarID, 'range', interp_test_vertcoord), &
           'test_interpolate_range', 'put_att vertrange '//trim(ncfilename))
call nc_check(nf90_put_att(ncid, vertVarID, 'cartesian_axis', 'Z'),   &
           'test_interpolate_range', 'vert cartesian_axis '//trim(ncfilename))

! loop over ensemble members
do imem = 1, ens_size
   if ( ens_size > 1) then
      write(field_name,'(A,I2)') "field_",imem
   else
      field_name = "field"
   endif
   call nc_check(nf90_def_var(ncid=ncid, name=field_name, xtype=nf90_double, &
           dimids=(/ nlonDimID, nlatDimID, nvertDimID /), varid=VarID(imem)), 'test_interpolate_range', &
                    'field def_var '//trim(ncfilename))
   kind_of_interest = get_raw_obs_kind_name(mykindindex)
   call nc_check(nf90_put_att(ncid, VarID(imem), 'long_name', kind_of_interest), &
              'test_interpolate_range', 'put_att field long_name '//trim(ncfilename))
   call nc_check(nf90_put_att(ncid, VarID(imem), '_FillValue', MISSING_R8), &
              'test_interpolate_range', 'put_att field FillValue '//trim(ncfilename))
   call nc_check(nf90_put_att(ncid, VarID(imem), 'missing_value', MISSING_R8), &
              'test_interpolate_range', 'put_att field missing_value '//trim(ncfilename))
   call nc_check(nf90_put_att(ncid, VarID(imem), 'interp_test_vertcoord', interp_test_vertcoord ), &
              'test_interpolate_range', 'put_att field interp_test_vertcoord '//trim(ncfilename))
enddo

! Leave define mode so we can fill the variables.
call nc_check(nf90_enddef(ncid), &
              'test_interpolate_range','field enddef '//trim(ncfilename))

! Fill the variables
call nc_check(nf90_put_var(ncid, lonVarID, lon), &
              'test_interpolate_range','lon put_var '//trim(ncfilename))
call nc_check(nf90_put_var(ncid, latVarID, lat), &
              'test_interpolate_range','lat put_var '//trim(ncfilename))
call nc_check(nf90_put_var(ncid, vertVarID, vert), &
              'test_interpolate_range','vert put_var '//trim(ncfilename))

do imem = 1, ens_size
   call nc_check(nf90_put_var(ncid, VarID(imem), field(:,:,:,imem)), &
                 'test_interpolate_range','field put_var '//trim(ncfilename))
enddo

! tidy up
call nc_check(nf90_close(ncid), &
             'test_interpolate_range','close '//trim(ncfilename))

deallocate(lon, lat, vert, field)

test_interpolate_range = nfailed

end function test_interpolate_range


!-------------------------------------------------------------------------------
! Do a single interpolation on a given location and kind.  Returns the
! interpolated values and ios_out. Returns the number of ensemble members that
! passed
!-------------------------------------------------------------------------------
function test_interpolate_single( ens_handle,       &
                                  ens_size,         &
                                  vertcoord_string, &
                                  lonval,           &
                                  latval,           &
                                  vertval,          &
                                  mykindindex,      &
                                  interp_vals,      &
                                  ios_out)

type(ensemble_type)   , intent(inout) :: ens_handle
integer               , intent(in)    :: ens_size
character(len=*)      , intent(in)    :: vertcoord_string
real(r8)              , intent(in)    :: lonval
real(r8)              , intent(in)    :: latval
real(r8)              , intent(in)    :: vertval
integer               , intent(in)    :: mykindindex
real(r8)              , intent(out)   :: interp_vals(ens_size)
integer               , intent(out)   :: ios_out(ens_size)

integer :: test_interpolate_single

type(location_type) :: loc
integer :: imem, num_passed, vertcoord

num_passed = 0

vertcoord = get_location_index(vertcoord_string)

loc = set_location(lonval, latval, vertval, vertcoord)

call model_interpolate(ens_handle, ens_size, loc, mykindindex, interp_vals, ios_out)

do imem = 1, ens_size
   if (ios_out(imem) == 0 ) then
      write(*,*) 'member ', imem, 'model_interpolate SUCCESS with value', interp_vals(imem)
      num_passed = num_passed + 1
   else
      write(*,*) 'member ', imem, 'model_interpolate ERROR with error code', ios_out(imem)
      test_interpolate_single = test_interpolate_single
   endif
enddo

test_interpolate_single = num_passed

end function test_interpolate_single

!-------------------------------------------------------------------------------
! Count the number of different error codes and output the results.  This
! is just a helper function for test_interpolate_range. Only sums error codes
! for the first ensemble member
!-------------------------------------------------------------------------------
subroutine count_error_codes(error_codes, num_failed)

integer, intent(in) :: error_codes(:,:)
integer, intent(in) :: num_failed

integer :: i, count_errors, results

count_errors = 1

i = 1
do while (count_errors < num_failed)
   results = count(error_codes(:,1) == i)
   if (results /= 0) then
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

end module test_interpolate_mod
