! DART software - Copyright UCAR. This open source software is provided
! by UCAR, "as is", without charge, subject to all terms of use at
! http://www.image.ucar.edu/DAReS/DART/DART_download
!
! $Id$

program exhaustion

!----------------------------------------------------------------------
! purpose: test routine for the new interpolation.  intended to
!   exhaustively test a large number of close locations and look
!   for unexpectedly large differences.
!
! right now this program still makes a rectangular grid and then
! tests a point a random distance from each existing point.  future
! versions will generate completely random pairs of points so the
! size of any arrays won't be an issue.
!
!----------------------------------------------------------------------

use        types_mod, only : r8, digits12, metadatalength, MISSING_R8, rad2deg
use    utilities_mod, only : initialize_utilities, finalize_utilities, nc_check, &
                             open_file, close_file, find_namelist_in_file, &
                             check_namelist_read, nmlfileunit, do_nml_file, do_nml_term, &
                             E_MSG, E_ERR, error_handler, get_unit
use   random_seq_mod, only : random_seq_type, init_random_seq, random_uniform
use     location_mod, only : location_type, set_location, write_location, get_dist, &
                             query_location, LocationDims, get_location, &
                             VERTISUNDEF, VERTISSURFACE, VERTISLEVEL, VERTISPRESSURE, &
                             VERTISHEIGHT, VERTISSCALEHEIGHT
use     obs_kind_mod, only : get_name_for_quantity, get_index_for_quantity, &
                             QTY_U_WIND_COMPONENT, QTY_V_WIND_COMPONENT
use  assim_model_mod, only : open_restart_read, open_restart_write, close_restart, &
                             aread_state_restart, awrite_state_restart, &
use time_manager_mod, only : time_type, set_calendar_type, GREGORIAN, &
                             read_time, get_time, set_time,  &
                             print_date, get_date, &
                             print_time, write_time, &
                             operator(-)
use        model_mod, only : static_init_model, get_model_size, get_state_meta_data, &
                             model_interpolate, get_analysis_time, &
                             get_model_analysis_filename, analysis_file_to_statevector, &
                             statevector_to_analysis_file, get_analysis_time,            &
                             write_model_time, get_grid_dims

use netcdf
use typesizes

implicit none

! version controlled file description for error handling, do not edit
character(len=256), parameter :: source   = &
   "$URL$"
character(len=32 ), parameter :: revision = "$Revision$"
character(len=128), parameter :: revdate  = "$Date$"

character(len=256) :: string1, string2

!------------------------------------------------------------------
! The namelist variables
!------------------------------------------------------------------

character (len = 129)  :: destroy_file         = 'temp_analysis_file.nc'
character (len = 129)  :: dart_input_file      = 'dart_ics'
character (len = 129)  :: output_file          = 'check_me'
logical                :: advance_time_present = .FALSE.
logical                :: verbose              = .TRUE.
logical                :: dojitter             = .TRUE.
logical                :: matlab_out           = .FALSE.
logical                :: netcdf_out           = .TRUE.
character (len = 129)  :: kind_of_interest     = 'QTY_POTENTIAL_TEMPERATURE'
real(r8)               :: interp_test_dlon     = 1.0
real(r8)               :: interp_test_dlat     = 1.0
real(r8)               :: interp_test_dvert    = 1.0
real(r8), dimension(2) :: interp_test_latrange = (/ -90.0,  90.0 /)
real(r8), dimension(2) :: interp_test_lonrange = (/   0.0, 360.0 /)
real(r8), dimension(2) :: interp_test_vertrange = (/  1000.0, 30000.0 /)
character(len=metadatalength) :: interp_test_vertcoord = 'VERTISHEIGHT'
real(r8)               :: hscale               = 100.0_r8
real(r8)               :: diff_threshold       = 100.0_r8
integer                :: pointcount           = 10000

namelist /exhaustion_nml/ dart_input_file, output_file, &
                        advance_time_present, verbose, &
                        matlab_out, netcdf_out, kind_of_interest, &
                        interp_test_dlon, interp_test_lonrange, &
                        interp_test_dlat, interp_test_latrange, &
                        interp_test_dvert, interp_test_vertrange, &
                        interp_test_vertcoord, destroy_file, hscale, &
                        diff_threshold, pointcount, dojitter

!----------------------------------------------------------------------
! other variables
!----------------------------------------------------------------------

integer :: ios_out, iunit, io, i
integer :: x_size, skip
integer :: mykindindex, vertcoord
real(r8) :: vmin, vmax, vscale

! lat/lon arrays for testing
real(r8), allocatable :: lon(:), lat(:), vert(:)
real(r8), allocatable :: full_lon(:,:,:), full_lat(:,:,:), full_vert(:,:,:)
real(r8), allocatable :: field(:,:,:), jit_field(:,:,:), field_diff(:,:,:)
integer :: nlon, nlat, nvert

type (random_seq_type) :: r
type(time_type)        :: model_time, adv_to_time
real(r8), allocatable  :: statevector(:)

character(len=129)  :: mpas_input_file  ! set with get_model_analysis_filename() if needed
type(location_type) :: loc


!----------------------------------------------------------------------
! program start
!----------------------------------------------------------------------

call setup_interpolate()

ios_out = test_interpolate(field)

if ( ios_out /= 0 ) then 
   print *, 'test interpolate had ', ios_out, ' failures.'
endif

if (dojitter) then
   call jitter_grid()

   ios_out = test_interpolate(jit_field, full_lon, full_lat, full_vert)
   
   call check_diffs(field, jit_field)
endif

call output_interpolate()

!call interpolate_centers()

!call interpolate_vertices()
!call interpolate_edges()

call random_test(pointcount)

call takedown_interpolate()

call finalize_utilities()

!----------------------------------------------------------------------
!----------------------------------------------------------------------

contains


!----------------------------------------------------------------------

subroutine setup_interpolate()
! allocate space, open files, etc.

! figure out what has to be global

integer :: ilon, jlat, kvert

call initialize_utilities(progname='exhaustion')
call set_calendar_type(GREGORIAN)

! Initialize repeatable random sequence
call init_random_seq(r)

print *, ''
print *, 'Reading the namelist.'

call find_namelist_in_file("input.nml", "exhaustion_nml", iunit)
read(iunit, nml = exhaustion_nml, iostat = io)
call check_namelist_read(iunit, io, "exhaustion_nml")

! Record the namelist values used for the run
if (do_nml_file()) write(nmlfileunit, nml=exhaustion_nml)
if (do_nml_term()) write(     *     , nml=exhaustion_nml)

call static_init_model()

x_size = get_model_size()

allocate(statevector(x_size))

print *, ' '
print *, 'Reading '//trim(dart_input_file)

iunit = open_restart_read(dart_input_file)
if ( advance_time_present ) then
   call aread_state_restart(model_time, statevector, iunit, adv_to_time)
else
   call aread_state_restart(model_time, statevector, iunit)
endif
call close_restart(iunit)

call print_date( model_time,'exhaustion:model date')
call print_time( model_time,'exhaustion:model time')


if ((interp_test_dlon < 0.0_r8) .or. (interp_test_dlat < 0.0_r8)) then
   write(*,*)'Error in interpolation test because one of'
   write(*,*)'interp_test_dlon,interp_test_dlat are < 0.0'
   write(*,*)'interp_test_dlon  = ',interp_test_dlon
   write(*,*)'interp_test_dlat  = ',interp_test_dlat
   write(*,*)'interp_test_dvert = ',interp_test_dvert
   stop
endif

! round down to avoid exceeding the specified range
nlat  = aint((interp_test_latrange(2) -  interp_test_latrange(1))/interp_test_dlat) + 1
nlon  = aint((interp_test_lonrange(2) -  interp_test_lonrange(1))/interp_test_dlon) + 1
nvert = aint((interp_test_vertrange(2) - interp_test_vertrange(1))/interp_test_dvert) + 1

allocate(lon(nlon), lat(nlat), vert(nvert))
allocate(full_lon(nlon,nlat,nvert), full_lat(nlon,nlat,nvert), full_vert(nlon,nlat,nvert))
allocate(field(nlon,nlat,nvert), jit_field(nlon,nlat,nvert), field_diff(nlon,nlat,nvert))
  
print *, 'interpolating to a grid ranging from: '
print *, nlon,  '  lons: ', interp_test_lonrange,  ' by ', interp_test_dlon
print *, nlat,  '  lats: ', interp_test_latrange,  ' by ', interp_test_dlat
print *, nvert, ' verts: ', interp_test_vertrange, ' by ', interp_test_dvert

! initialize the arrays

field = MISSING_R8
jit_field = MISSING_R8

do ilon = 1, nlon
   lon(ilon) = interp_test_lonrange(1) + real(ilon-1,r8) * interp_test_dlon
   do jlat = 1, nlat
      lat(jlat) = interp_test_latrange(1) + real(jlat-1,r8) * interp_test_dlat
      do kvert = 1, nvert
         vert(kvert) = interp_test_vertrange(1) + real(kvert-1,r8) * interp_test_dvert
      enddo
   end do 
end do


print *, 'vertical coord is ', interp_test_vertcoord
select case(trim(interp_test_vertcoord))
   case ('VERTISUNDEF')
      vertcoord = VERTISUNDEF
      vmin   = 0.0_r8
      vmax   = 1.0_r8
      vscale = 1.0_r8
   case ('VERTISSURFACE')
      vertcoord = VERTISSURFACE
      vmin   = 0.0_r8
      vmax   = 1.0_r8
      vscale = 1.0_r8
   case ('VERTISLEVEL')
      vertcoord = VERTISLEVEL
      vmin   =  1.0_r8
      vmax   = 41.0_r8
      vscale = 0.01_r8
   case ('VERTISPRESSURE')
      vertcoord = VERTISPRESSURE
      vmin   =  30000.0_r8
      vmax   = 100000.0_r8
      vscale =    100.0_r8
   case ('VERTISHEIGHT')
      vertcoord = VERTISHEIGHT
      vmin   =     0.0_r8
      vmax   = 30000.0_r8
      vscale =   100.0_r8
   case ('VERTISSCALEHEIGHT')
      vertcoord = VERTISSCALEHEIGHT
      vmin   = 0.0_r8
      vmax   = 5.0_r8
      vscale = 0.01_r8
   case default
      write(string1,*) 'unknown vertcoord ', trim(interp_test_vertcoord)
      call error_handler(E_ERR,'test_interpolate',string1,source,revision,revdate)
end select

print *, 'interpolating kind: ', kind_of_interest
mykindindex = get_index_for_quantity(kind_of_interest)

end subroutine setup_interpolate

!----------------------------------------------------------------------

subroutine takedown_interpolate()

deallocate(statevector)
deallocate(lon, lat, vert)
deallocate(field, jit_field, field_diff)
deallocate(full_lon, full_lat, full_vert)

end subroutine takedown_interpolate

!----------------------------------------------------------------------

function jitter_location(inloc)
 type(location_type), intent(in)  :: inloc
 type(location_type) :: jitter_location

! take existing lat/lon/vert location and add small 
! noise to each part of the coordinate.

real(r8) :: x(3), lon, lat, vert
real(r8) :: dellon, dellat, delvert, del


x = get_location(inloc)
lon = x(1)
lat = x(2)
vert = x(3)

call jitter_point(lon, lat, vert, dellon, dellat, delvert)

jitter_location = set_location(dellon, dellat, delvert, vertcoord)

end function jitter_location

!----------------------------------------------------------------------

subroutine jitter_point(lon, lat, vert, dellon, dellat, delvert)
 real(r8), intent(in ) :: lon, lat, vert
 real(r8), intent(out) :: dellon, dellat, delvert

! take existing lat/lon/vert location and add small 
! noise to each part of the coordinate.

real(r8) :: x(3), del


del = (random_uniform(r) - 0.5_r8)  / hscale
dellon = lon + del
if (dellon <   0.0_r8) dellon = -dellon
if (dellon > 360.0_r8) dellon = 360.0_r8 - (dellon - 360.0_r8)

del = (random_uniform(r) - 0.5_r8)  / hscale
dellat = lat + del
if (dellat < -90.0_r8) dellat = -90.0_r8 - (dellat + 90.0_r8)
if (dellat >  90.0_r8) dellat =  90.0_r8 - (dellat - 90.0_r8)

del = (random_uniform(r) - 0.5_r8) * vscale
delvert = vert + del
if (delvert < vmin) delvert = vert + (delvert - vmin)
if (delvert > vmax) delvert = vert - (delvert - vmax)

end subroutine jitter_point

!----------------------------------------------------------------------

subroutine jitter_grid()

! take existing lat/lon/vert array and add tiny noise to location
! should call random number gen

integer  :: ilon, jlat, kvert

print *, 'jittering coordinates with +/- hscale, vscale: ', 0.5_r8 / hscale, 0.5_r8 / vscale

do ilon = 1, nlon
   do jlat = 1, nlat
      do kvert = 1, nvert

         call jitter_point(lon(ilon), lat(jlat), vert(kvert), &
                           full_lon(ilon, jlat, kvert), &
                           full_lat(ilon, jlat, kvert), &
                           full_vert(ilon, jlat, kvert))

      enddo
   end do 
end do

print *, 'jittered field, min/max  lon: ', minval(full_lon),  maxval(full_lon)
print *, 'jittered field, min/max  lat: ', minval(full_lat),  maxval(full_lat)
print *, 'jittered field, min/max vert: ', minval(full_vert), maxval(full_vert)

end subroutine jitter_grid

!----------------------------------------------------------------------

function test_interpolate(f, flon, flat, fvert)
 real(r8), intent(inout) :: f(:,:,:)
 real(r8), intent(inout), optional :: flon(:,:,:)
 real(r8), intent(inout), optional :: flat(:,:,:)
 real(r8), intent(inout), optional :: fvert(:,:,:)

! function to exercise the model_mod:model_interpolate() function
integer :: test_interpolate

! Local variables

integer :: ilon, jlat, kvert, nfailed
logical :: compute_coords

compute_coords = .true.
if (present(flon) .and. present(flat) .and. present(fvert)) then
   compute_coords = .false.
endif

test_interpolate = 0   ! normal termination
nfailed = 0

do ilon = 1, nlon
   do jlat = 1, nlat
      do kvert = 1, nvert

         if (compute_coords) then
            loc = set_location(lon(ilon), lat(jlat), vert(kvert), vertcoord)
         else
            loc = set_location(flon(ilon,jlat,kvert), flat(ilon,jlat,kvert), &
                               fvert(ilon,jlat,kvert), vertcoord)
         endif

         call model_interpolate(statevector, loc, mykindindex, f(ilon,jlat,kvert), ios_out)

         if (ios_out /= 0) call report_loc_badinterp(loc, ios_out, nfailed)

      enddo
   end do 
end do

print *, 'total interpolations, num failed: ', nlon*nlat*nvert, nfailed
print *, ' '

test_interpolate = nfailed

end function test_interpolate

!----------------------------------------------------------------------

subroutine check_diffs(f1, f2)
 real(r8), intent(in) :: f1(:,:,:), f2(:,:,:)

integer  :: ilon, jlat, kvert
real(r8) :: v1, v2


print *, ' '
print *, 'checking and reporting differences larger than ', diff_threshold
print *, ' '

field_diff = f1 - f2

do ilon = 1, nlon
   do jlat = 1, nlat
      do kvert = 1, nvert
         v1 = f1(ilon,jlat,kvert)
         v2 = f2(ilon,jlat,kvert)
         if (abs(v1 - v2) > diff_threshold) then
            call report_diff(lon(ilon), lat(jlat), vert(kvert), &
                             full_lon(ilon,jlat,kvert), full_lat(ilon,jlat,kvert), &
                             full_vert(ilon,jlat,kvert), v1, v2)
         endif
         if (v1 == MISSING_R8 .or. v2 == MISSING_R8) then
            field_diff(ilon,jlat,kvert) = MISSING_R8
         endif
      enddo
   end do 
end do


print *, 'min, max diffs in data: ', minval(field_diff, mask= field_diff /= MISSING_R8), &
                                     maxval(field_diff)


end subroutine check_diffs

!----------------------------------------------------------------------

subroutine random_test(n)
 integer, intent(in) :: n

integer :: i

real(r8) :: x(3), lon, lat, vert
real(r8) :: dellon, dellat, delvert, del
real(r8) :: lonspan, latspan, vertspan
real(r8) :: oval, jval
type(location_type) :: original, jittered
integer  :: oout, jout

print *, ' '
print *, 'starting random point test for ', n, ' points, '
print *, 'checking and reporting differences larger than ', diff_threshold
print *, ' '

lonspan = interp_test_lonrange(2) - interp_test_lonrange(1)
latspan = interp_test_latrange(2) - interp_test_latrange(1)
vertspan = interp_test_vertrange(2) - interp_test_vertrange(1)

do i=1, n

   ! generate a random lon, lat, vert in range
   
   lon = (random_uniform(r) * lonspan) + interp_test_lonrange(1)
   lat = (random_uniform(r) * latspan) + interp_test_latrange(1)
   vert = (random_uniform(r) * vertspan) + interp_test_vertrange(1)
   
   original = set_location(lon, lat, vert, vertcoord)

   
   ! jitter it
   
   jittered = jitter_location(original)
   
   
   ! interp both and compute difference
   
   call model_interpolate(statevector, original, mykindindex, oval, oout)
   if (oout /= 0) then
      call report_loc_badinterp(original, oout)
      cycle
   endif
   
   call model_interpolate(statevector, jittered, mykindindex, jval, jout)
   if (jout /= 0) then
      call report_loc_badinterp(jittered, jout)
      cycle
   endif
   
   ! if del > threshold, print all info
   if (abs(oval - jval) > diff_threshold) call report_loc_diff(original, jittered, oval, jval)
   
   ! repeat
!if (mod(n,1000) == 0) print *, 'random test, iteration ', i
enddo

end subroutine random_test

!----------------------------------------------------------------------

subroutine report_loc_badinterp(loc, rc, failcount)
 type(location_type), intent(in)  :: loc
 integer,             intent(in)  :: rc
 integer, intent(inout), optional :: failcount

real(r8) :: l(3)

l = get_location(loc)

call report_badinterp(l(1), l(2), l(3), rc, failcount)

end subroutine report_loc_badinterp

!----------------------------------------------------------------------

subroutine report_badinterp(lon, lat, vert, rc, failcount)
 real(r8), intent(in)              :: lon, lat, vert
 integer,  intent(in)              :: rc
 integer,  intent(inout), optional :: failcount

if (verbose) then
   write(string2,'(A,3(1x,f14.6))') 'lon,lat,vert ', lon, lat, vert
   write(string1,*) 'interpolation return code was', rc
   call error_handler(E_MSG,'test_interpolate',string1,source,revision,revdate,text2=string2)
endif


if (present(failcount)) then
   failcount = failcount + 1
endif

end subroutine report_badinterp

!----------------------------------------------------------------------

subroutine report_loc_diff(loc1, loc2, v1, v2)
 type(location_type), intent(in) :: loc1, loc2
 real(r8),            intent(in) :: v1, v2

real(r8) :: l1(3), l2(3)

l1 = get_location(loc1)
l2 = get_location(loc2)

call report_diff(l1(1), l1(2), l1(3), &
                 l2(1), l2(2), l2(3), &
                 v1, v2)

end subroutine report_loc_diff

!----------------------------------------------------------------------

subroutine report_diff(lon1, lat1, vert1, lon2, lat2, vert2, v1, v2)
 real(r8), intent(in) :: lon1, lat1, vert1, lon2, lat2, vert2
 real(r8), intent(in) :: v1, v2

print *, lon1, lat1, vert1, v1, ' vs '

print *, lon2, lat2, vert2, v2, ' = ', abs(v1-v2)

end subroutine report_diff

!----------------------------------------------------------------------

subroutine output_interpolate()

! matlab, netcdf, or both

! This will result in a netCDF file with all salient metadata 
integer :: test_interpolate

! Local variables

integer :: ilon, jlat, kvert
character(len=128) :: ncfilename,txtfilename

character(len=8)      :: crdate      ! needed by F90 DATE_AND_TIME intrinsic
character(len=10)     :: crtime      ! needed by F90 DATE_AND_TIME intrinsic
character(len=5)      :: crzone      ! needed by F90 DATE_AND_TIME intrinsic
integer, dimension(8) :: values      ! needed by F90 DATE_AND_TIME intrinsic

integer :: ncid, nlonDimID, nlatDimID, nvertDimID
integer :: VarID, VarID2, lonVarID, latVarID, vertVarID


! matlab section

if (matlab_out) then
   write(txtfilename,'(a,a)')trim(output_file),'_test.m'
   
   iunit = open_file(trim(txtfilename), action='write')
   write(iunit,'(''missingvals = '',f12.4,'';'')')MISSING_R8
   write(iunit,'(''nlon = '',i8,'';'')')nlon
   write(iunit,'(''nlat = '',i8,'';'')')nlat
   write(iunit,'(''nvert = '',i8,'';'')')nvert
   write(iunit,'(''interptest = [ ... '')')
   
   do ilon = 1, nlon
      do jlat = 1, nlat
         do kvert = 1, nvert
            write(iunit,*) field(ilon,jlat,kvert)
         enddo
      end do 
   end do
   
   write(iunit,'(''];'')')
   write(iunit,'(''datmat = reshape(interptest,nvert,nlat,nlon);'')')
   write(iunit,'(''datmat(datmat == missingvals) = NaN;'')')
   call close_file(iunit)
endif

! Write out the netCDF file for easy exploration.

if (netcdf_out) then
   write( ncfilename,'(a,a)')trim(output_file),'_test.nc'
   
   call DATE_AND_TIME(crdate,crtime,crzone,values)
   write(string1,'(''YYYY MM DD HH MM SS = '',i4,5(1x,i2.2))') &
                     values(1), values(2), values(3), values(5), values(6), values(7)
   
   call nc_check( nf90_create(path=trim(ncfilename), cmode=NF90_clobber, ncid=ncid), &
                     'test_interpolate', 'open '//trim(ncfilename))
   call nc_check( nf90_put_att(ncid, NF90_GLOBAL, 'creation_date' ,trim(string1) ), &
                     'test_interpolate', 'creation put '//trim(ncfilename))
   call nc_check( nf90_put_att(ncid, NF90_GLOBAL, 'parent_file', dart_input_file ), &
                     'test_interpolate', 'put_att filename '//trim(ncfilename))
   
   ! Define dimensions
   
   call nc_check(nf90_def_dim(ncid=ncid, name='lon', len=nlon, &
           dimid = nlonDimID),'test_interpolate', 'nlon def_dim '//trim(ncfilename))
   
   call nc_check(nf90_def_dim(ncid=ncid, name='lat', len=nlat, &
           dimid = nlatDimID),'test_interpolate', 'nlat def_dim '//trim(ncfilename))
   
   call nc_check(nf90_def_dim(ncid=ncid, name='vert', len=nvert, &
           dimid = nvertDimID),'test_interpolate', 'nvert def_dim '//trim(ncfilename))
   
   ! Define variables
   
   call nc_check(nf90_def_var(ncid=ncid, name='lon', xtype=nf90_double, &
           dimids=nlonDimID, varid=lonVarID), 'test_interpolate', &
                    'lon def_var '//trim(ncfilename))
   call nc_check(nf90_put_att(ncid, lonVarID, 'range', interp_test_lonrange), &
              'test_interpolate', 'put_att lonrange '//trim(ncfilename))
   call nc_check(nf90_put_att(ncid, lonVarID, 'cartesian_axis', 'X'),   &
              'test_interpolate', 'lon cartesian_axis '//trim(ncfilename))
   
   
   call nc_check(nf90_def_var(ncid=ncid, name='lat', xtype=nf90_double, &
           dimids=nlatDimID, varid=latVarID), 'test_interpolate', &
                    'lat def_var '//trim(ncfilename))
   call nc_check(nf90_put_att(ncid, latVarID, 'range', interp_test_latrange), &
              'test_interpolate', 'put_att latrange '//trim(ncfilename))
   call nc_check(nf90_put_att(ncid, latVarID, 'cartesian_axis', 'Y'),   &
              'test_interpolate', 'lat cartesian_axis '//trim(ncfilename))
   
   call nc_check(nf90_def_var(ncid=ncid, name='vert', xtype=nf90_double, &
           dimids=nvertDimID, varid=vertVarID), 'test_interpolate', &
                    'vert def_var '//trim(ncfilename))
   call nc_check(nf90_put_att(ncid, vertVarID, 'range', interp_test_vertrange), &
              'test_interpolate', 'put_att vertrange '//trim(ncfilename))
   call nc_check(nf90_put_att(ncid, vertVarID, 'cartesian_axis', 'Z'),   &
              'test_interpolate', 'vert cartesian_axis '//trim(ncfilename))
   
   call nc_check(nf90_def_var(ncid=ncid, name='field', xtype=nf90_double, &
           dimids=(/ nlonDimID, nlatDimID, nvertDimID /), varid=VarID), 'test_interpolate', &
                    'field def_var '//trim(ncfilename))
   call nc_check(nf90_put_att(ncid, VarID, 'long_name', kind_of_interest), &
              'test_interpolate', 'put_att field long_name '//trim(ncfilename))
   call nc_check(nf90_put_att(ncid, VarID, '_FillValue', MISSING_R8), &
              'test_interpolate', 'put_att field FillValue '//trim(ncfilename))
   call nc_check(nf90_put_att(ncid, VarID, 'missing_value', MISSING_R8), &
              'test_interpolate', 'put_att field missing_value '//trim(ncfilename))
   call nc_check(nf90_put_att(ncid, VarID, 'vertcoord_string', interp_test_vertcoord ), &
              'test_interpolate', 'put_att field vertcoord_string '//trim(ncfilename))
   call nc_check(nf90_put_att(ncid, VarID, 'vertcoord', vertcoord ), &
              'test_interpolate', 'put_att field vertcoord '//trim(ncfilename))
   
if (dojitter) then
   call nc_check(nf90_def_var(ncid=ncid, name='field_diff', xtype=nf90_double, &
           dimids=(/ nlonDimID, nlatDimID, nvertDimID /), varid=VarID2), 'test_interpolate', &
                    'field def_var '//trim(ncfilename))
   call nc_check(nf90_put_att(ncid, VarID2, 'long_name', kind_of_interest), &
              'test_interpolate', 'put_att field long_name '//trim(ncfilename))
   call nc_check(nf90_put_att(ncid, VarID2, '_FillValue', MISSING_R8), &
              'test_interpolate', 'put_att field FillValue '//trim(ncfilename))
   call nc_check(nf90_put_att(ncid, VarID2, 'missing_value', MISSING_R8), &
              'test_interpolate', 'put_att field missing_value '//trim(ncfilename))
   call nc_check(nf90_put_att(ncid, VarID2, 'vertcoord_string', interp_test_vertcoord ), &
              'test_interpolate', 'put_att field vertcoord_string '//trim(ncfilename))
   call nc_check(nf90_put_att(ncid, VarID2, 'vertcoord', vertcoord ), &
              'test_interpolate', 'put_att field vertcoord '//trim(ncfilename))
endif
   
   ! Leave define mode so we can fill the variables.
   call nc_check(nf90_enddef(ncid), &
                 'test_interpolate','field enddef '//trim(ncfilename))
   
   ! Fill the variables
   call nc_check(nf90_put_var(ncid, lonVarID, lon), &
                 'test_interpolate','lon put_var '//trim(ncfilename))
   call nc_check(nf90_put_var(ncid, latVarID, lat), &
                 'test_interpolate','lat put_var '//trim(ncfilename))
   call nc_check(nf90_put_var(ncid, vertVarID, vert), &
                 'test_interpolate','vert put_var '//trim(ncfilename))
   call nc_check(nf90_put_var(ncid, VarID, field), &
                 'test_interpolate','field put_var '//trim(ncfilename))
if (dojitter) then
   call nc_check(nf90_put_var(ncid, VarID2, field_diff), &
                 'test_interpolate','field put_var '//trim(ncfilename))
endif
   
   ! tidy up
   call nc_check(nf90_close(ncid), &
                'test_interpolate','close '//trim(ncfilename))
endif

end subroutine output_interpolate


!----------------------------------------------------------------------
!----------------------------------------------------------------------
!----------------------------------------------------------------------
! unused (for the moment) below here.
!----------------------------------------------------------------------


! FIXME: repurpose this to simply interpolate to centers, edges, verts, etc
! and diff them with the original file to see if they are goofy.

subroutine interpolate_centers()

! The u(Time, nEdges, nVertLevels) variable should already be in the DART state
! vector. We are reading in the locations of the cell centers and predicting the
! double uReconstructMeridional(Time, nCells, nVertLevels) ;
! double uReconstructZonal(Time, nCells, nVertLevels) ;
! values and OVERWRITING the values in the 'destroy_file' ...
! This should allow us to easily create a difference file and plot that ...
!

integer :: ncid, dimid, VarID, uVarID, vVarID 
integer :: nCells, nVertLevelsP1, nVertLevels, xloc, zloc
integer :: nFailedU, nFailedV
real(r8), allocatable, dimension(:)   :: lonCell, latCell
real(r8), allocatable, dimension(:,:) :: zGridFace, zGridCenter, uhat, vhat

character(len=8)      :: crdate      ! needed by F90 DATE_AND_TIME intrinsic
character(len=10)     :: crtime      ! needed by F90 DATE_AND_TIME intrinsic
character(len=5)      :: crzone      ! needed by F90 DATE_AND_TIME intrinsic
integer, dimension(8) :: values      ! needed by F90 DATE_AND_TIME intrinsic

call nc_check( nf90_open(trim(destroy_file), NF90_WRITE, ncid), 'interpolate_centers', 'open '//trim(destroy_file))

! Get required dimensions

call nc_check(nf90_inq_dimid(ncid, 'nCells', dimid),           &
              'interpolate_centers','inq_dimid nCells '//trim(destroy_file))
call nc_check(nf90_inquire_dimension(ncid, dimid, len=nCells), &
              'interpolate_centers','inquire_dimension nCells '//trim(destroy_file))

call nc_check(nf90_inq_dimid(ncid, 'nVertLevelsP1', dimid),           &
              'interpolate_centers','inq_dimid nVertLevelsP1 '//trim(destroy_file))
call nc_check(nf90_inquire_dimension(ncid, dimid, len=nVertLevelsP1), &
              'interpolate_centers','inquire_dimension nVertLevelsP1 '//trim(destroy_file))

call nc_check(nf90_inq_dimid(ncid, 'nVertLevels', dimid),           &
              'interpolate_centers','inq_dimid NVertLevels '//trim(destroy_file))
call nc_check(nf90_inquire_dimension(ncid, dimid, len=nVertLevels), &
              'interpolate_centers','inquire_dimension NVertLevels '//trim(destroy_file))

allocate(lonCell(nCells), latCell(nCells))
allocate(zGridFace(nVertLevelsP1, nCells))
allocate(zGridCenter(nVertLevels, nCells))
allocate(       uhat(nVertLevels, nCells))
allocate(       vhat(nVertLevels, nCells))

call nc_check(nf90_inq_varid(ncid, "lonCell", VarID), &
              'interpolate_centers', 'inq_varid lonCell'//trim(destroy_file))
call nc_check(nf90_get_var(  ncid, VarID, lonCell),   &
              'interpolate_centers', 'get_var   lonCell'//trim(destroy_file))
call nc_check(nf90_inq_varid(ncid, "latCell", VarID), &
              'interpolate_centers', 'inq_varid latCell'//trim(destroy_file))
call nc_check(nf90_get_var(  ncid, VarID, latCell),   &
              'interpolate_centers', 'get_var   latCell'//trim(destroy_file))
call nc_check(nf90_inq_varid(ncid, 'zgrid', VarID),   &
              'interpolate_centers', 'inq_varid zgrid '//trim(destroy_file))
call nc_check(nf90_get_var(  ncid, VarID, zGridFace), &
              'interpolate_centers', 'get_var   zgrid '//trim(destroy_file))

! compute vertical center locations
do xloc=1, nCells
do zloc=1, nVertLevels
   zGridCenter(zloc,xloc) = (zGridFace(zloc,xloc) + zGridFace(zloc+1,xloc))*0.5_r8
enddo
enddo

latCell(:) = rad2deg*latCell
lonCell(:) = rad2deg*lonCell
uhat(:,:)  = MISSING_R8
vhat(:,:)  = MISSING_R8

nFailedU = 0
nFailedV = 0

! FIXME : make sure zgrid is VERTISHEIGHT
do xloc = 1, nCells
do zloc = 1, nVertLevels
   loc = set_location(lonCell(xloc),latCell(xloc), zGridCenter(zloc,xloc), VERTISHEIGHT)
   call model_interpolate(statevector, loc, QTY_U_WIND_COMPONENT, uhat(zloc,xloc), ios_out)
   if (ios_out /= 0) nFailedU = nFailedU + 1
   call model_interpolate(statevector, loc, QTY_V_WIND_COMPONENT, vhat(zloc,xloc), ios_out)
   if (ios_out /= 0) nFailedV = nFailedV + 1
enddo
if (mod(xloc, 100) == 0) then
   print *, 'finished interpolating ', xloc, ' of ', nCells, ' cells'
endif
enddo

write(*,*)'uReconstructed interpolations ',nCells*nVertLevels,'possible.'
write(*,*)'                       zonal: ',nFailedU,'failures'
write(*,*)'                  meridional: ',nFailedV,'failures'

! Grab the variable ID and replace ...

write(*,*)'Overwriting uReconstructZonal,uReconstructMeridional in '//trim(destroy_file) 
call nc_check(nf90_redef(ncid),'interpolate_centers','redef '//trim(destroy_file))

call DATE_AND_TIME(crdate,crtime,crzone,values)
write(string1,'(''YYYY MM DD HH MM SS = '',i4,5(1x,i2.2))') &
                  values(1), values(2), values(3), values(5), values(6), values(7)

call nc_check( nf90_put_att(ncid, NF90_GLOBAL, 'DART_time' ,trim(string1) ), &
                  'interpolate_centers', 'creation put '//trim(destroy_file))

call nc_check(nf90_inq_varid(ncid,               'uReconstructZonal', uVarID), &
              'interpolate_centers', 'inq_varid uReconstructZonal '//trim(destroy_file))
call nc_check(nf90_put_att(ncid, uVarID, 'history', 'interpolated by DART'), &
              'interpolate_centers', 'uReconstructZonal put_att history '//trim(destroy_file))
call nc_check(nf90_put_att(ncid, uVarID, 'parent_file', mpas_input_file ), &
              'interpolate_centers', 'uReconstructZonal put_att parent_file '//trim(destroy_file))
call nc_check(nf90_put_att(ncid, uVarID, '_FillValue', MISSING_R8), &
              'interpolate_centers', 'uReconstructZonal put_att FillValue '//trim(destroy_file))
call nc_check(nf90_put_att(ncid, uVarID, 'missing_value', MISSING_R8), &
           '   interpolate_centers', 'uReconstructZonal put_att missing_value '//trim(destroy_file))

call nc_check(nf90_inq_varid(ncid,               'uReconstructMeridional', vVarID), &
              'interpolate_centers', 'inq_varid uReconstructMeridional '//trim(destroy_file))
call nc_check(nf90_put_att(ncid, vVarID, 'history', 'interpolated by DART'), &
              'interpolate_centers', 'uReconstructMeridional put_att history '//trim(destroy_file))
call nc_check(nf90_put_att(ncid, vVarID, 'parent_file', mpas_input_file ), &
              'interpolate_centers', 'uReconstructMeridional put_att parent_file '//trim(destroy_file))
call nc_check(nf90_put_att(ncid, vVarID, '_FillValue', MISSING_R8), &
              'interpolate_centers', 'uReconstructMeridional put_att FillValue '//trim(destroy_file))
call nc_check(nf90_put_att(ncid, vVarID, 'missing_value', MISSING_R8), &
           '   interpolate_centers', 'uReconstructMeridional put_att missing_value '//trim(destroy_file))

! Leave define mode so we can fill the variables.
call nc_check(nf90_enddef(ncid), 'interpolate_centers','enddef '//trim(destroy_file))

! Fill the variables, replacing the first time step
! The start/count arrays are in C-order ... i.e. opposite to the Fortran shape declarations.
! FIXME ... should match the timestep of the parent file ... which hopefully only has one timestep.

call nc_check(nf90_put_var( ncid, uVarID, uhat, start = (/ 1,1,1 /), count=(/ nVertLevels,nCells,1 /)), &
              'interpolate_centers', 'put_var uReconstructZonal '//trim(destroy_file))
call nc_check(nf90_put_var( ncid, vVarID, vhat, start = (/ 1,1,1 /), count=(/ nVertLevels,nCells,1 /)), &
              'interpolate_centers', 'put_var uReconstructMeridional '//trim(destroy_file))

! tidy up
call nc_check(nf90_close(ncid), 'interpolate_centers','close '//trim(destroy_file))

deallocate(lonCell, latCell, uhat, vhat, zGridFace, zGridCenter)

end subroutine interpolate_centers

subroutine interpolate_vertices
end subroutine interpolate_vertices

subroutine interpolate_edges
end subroutine interpolate_edges

end program exhaustion

! <next few lines under version control, do not edit>
! $URL$
! $Id$
! $Revision$
! $Date$
