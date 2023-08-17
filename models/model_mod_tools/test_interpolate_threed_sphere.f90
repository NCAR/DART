! DART software - Copyright UCAR. This open source software is provided
! by UCAR, "as is", without charge, subject to all terms of use at
! http://www.image.ucar.edu/DAReS/DART/DART_download
!
! $Id$

module test_interpolate_mod

!-------------------------------------------------------------------------------
! interpolation test routines for threed sphere locations.
!-------------------------------------------------------------------------------

use             types_mod, only : r8, i8, MISSING_R8, metadatalength

use         utilities_mod, only : register_module, error_handler, E_MSG, E_ERR, &
                                  initialize_utilities, finalize_utilities,     &
                                  find_namelist_in_file, check_namelist_read,   &
                                  E_MSG, open_file, close_file, do_output
 
use  netcdf_utilities_mod, only : nc_check, nc_create_file, nc_close_file, &
                                  nc_end_define_mode

use          location_mod, only : location_type, set_location, write_location,  &
                                  get_dist, get_location, LocationDims, &
                                  VERTISUNDEF, VERTISSURFACE, VERTISLEVEL, &
                                  VERTISPRESSURE, VERTISHEIGHT, VERTISSCALEHEIGHT, &
                                  query_location

use          obs_kind_mod, only : get_name_for_quantity, get_index_for_quantity

use  ensemble_manager_mod, only : ensemble_type

use model_check_utilities_mod, only : test_single_interpolation, &
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
      write(*,*) 'interp_test_dlon  = ',interp_test_dlon
      write(*,*) 'interp_test_dlat  = ',interp_test_dlat
      write(*,*) 'interp_test_dvert = ',interp_test_dvert
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
               call error_handler(E_MSG, routine, string1, &
                                  source, revision, revdate, text2=string2)
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
!> Do a single interpolation on a given location and kind.
!> Returns the interpolated values and ios_out.
!> Returns the number of ensemble members that passed.

function test_interpolate_single( ens_handle,       &
                                  ens_size,         &
                                  vertcoord_string, &
                                  lonval,           &
                                  latval,           &
                                  vertval,          &
                                  quantity_string,  &
                                  interp_vals,      &
                                  ios_out)

type(ensemble_type),intent(inout) :: ens_handle
integer            ,intent(in)    :: ens_size
character(len=*)   ,intent(in)    :: vertcoord_string
real(r8)           ,intent(in)    :: lonval
real(r8)           ,intent(in)    :: latval
real(r8)           ,intent(in)    :: vertval
character(len=*)   ,intent(in)    :: quantity_string
real(r8)           ,intent(out)   :: interp_vals(ens_size)
integer            ,intent(out)   :: ios_out(ens_size)

integer :: test_interpolate_single

type(location_type) :: location
integer :: vertcoord

vertcoord = get_location_index(vertcoord_string)
location = set_location(lonval, latval, vertval, vertcoord)

test_interpolate_single = test_single_interpolation(ens_handle, ens_size, &
                               location, vertcoord_string, &
                               quantity_string, interp_vals, ios_out)

end function test_interpolate_single


!-------------------------------------------------------------------------------
!> need to convert the character string for the test vertical coordinate into
!> the corresponding dart index.

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

!-----------------------------------------------------------------------
!> Expensive exhaustive search to find the indices into the
!> state vector of a particular lon/lat/vert. At present, only for a
!> single variable - could be extended to identify the closest location
!> for every variable in each domain. This could help ensure grid
!> staggering is being handled correctly.
!> This differs from model_mod_tools:find_closest_gridpoint in that
!> this version has to do something with the vertical coordinate system.
!> Right now it is not perfect because it will indicate that something
!> at the same horizontal location and pressure = 1000 is at exactly
!> the same location as something with height = 1000 ... ugh.

subroutine find_closest_gridpoint(loc_of_interest, vertcoord_string, quantity_string)

real(r8),         intent(in) :: loc_of_interest(:)
character(len=*), intent(in) :: vertcoord_string
character(len=*), intent(in) :: quantity_string

character(len=*), parameter :: routine = ''   ! name not important in context

type(location_type)   :: location, loc0, loc1
integer(i8)           :: i
integer               :: quantity_index, var_type, vert_type
real(r8)              :: closest, rlon, rlat, rvert
logical               :: matched
real(r8), allocatable :: thisdist(:)
real(r8),   parameter :: FARAWAY = huge(r8)
character(len=metadatalength) :: myquantity

!>@todo there should be arrays of length state_structure_mod:get_num_variables(domid)
!>      get_num_domains(), get_num_variables() ...

vert_type = get_location_index(vertcoord_string)
location  = set_location((/ loc_of_interest, real(vert_type,r8) /))

allocate( thisdist(get_model_size()) )
thisdist  = FARAWAY
matched   = .false.

! Trying to support the ability to specify matching a particular QUANTITY.
! With staggered grids, the closest gridpoint might not be of the quantity
! of interest.

quantity_index  = get_index_for_quantity(quantity_string)
rlon  = loc_of_interest(1)
rlat  = loc_of_interest(2)
rvert = loc_of_interest(3)

write(string1,'("Checking for the indices into the state vector that are close to")')
call write_location(0, location, charstring=string2)
write(string3,'("for (",A,") variables.")')trim(quantity_string)
call error_handler(E_MSG,routine,string1,text2=string2,text3=string3)
call error_handler(E_MSG,routine,'')

! Since there can be/will be multiple variables with
! identical distances, we will just cruise once through
! the array and come back to find all the 'identical' values.

DISTANCE : do i = 1,get_model_size()

   call get_state_meta_data(i, loc1, var_type)

   if (var_type .ne. quantity_index) cycle DISTANCE

   ! Grab the vert_type from the grid and
   ! set out target location to have the same.
   ! Compute the distance.

   vert_type   = nint(query_location(loc1))
   loc0        = set_location(rlon, rlat, rvert, vert_type)
   thisdist(i) = get_dist( loc1, loc0, no_vert=.false.)
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

      call write_location(0, loc1, charstring=string1)
      write(string2,'(A,I12,A)')' is index ',i,' ('//trim(myquantity)//')'
      call error_handler(E_MSG, routine, string1,text2=string2)
   endif

enddo REPORT

deallocate( thisdist )

end subroutine find_closest_gridpoint

!-------------------------------------------------------------------------------
! End of test_interpolate_mod
!-------------------------------------------------------------------------------

end module test_interpolate_mod

