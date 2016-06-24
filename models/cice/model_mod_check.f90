! DART software - Copyright 2004 - 2013 UCAR. This open source software is
! provided by UCAR, "as is", without charge, subject to all terms of use at
! http://www.image.ucar.edu/DAReS/DART/DART_download
!
! $Id: model_mod_check.f90 8565 2015-09-11 17:16:08Z hkershaw $

program model_mod_check

!----------------------------------------------------------------------
! purpose: test routines
!----------------------------------------------------------------------

use        types_mod, only : r8, metadatalength, MISSING_R8
use    utilities_mod, only : initialize_utilities, finalize_utilities, nc_check, &
                             open_file, close_file, find_namelist_in_file, &
                             check_namelist_read, nmlfileunit, do_nml_file, &
                             do_nml_term, E_MSG, E_ERR, error_handler, file_exist
use     location_mod, only : location_type, set_location, write_location, get_dist, &
                             query_location, LocationDims, get_location, &
                             VERTISLEVEL
use     obs_kind_mod, only : get_raw_obs_kind_name, get_raw_obs_kind_index
use  assim_model_mod, only : netcdf_file_type, aoutput_diagnostics, &
                             init_diag_output, finalize_diag_output

use state_vector_io_mod, only : open_restart_read, open_restart_write, close_restart, &
                             aread_state_restart, awrite_state_restart

use time_manager_mod, only : time_type, set_calendar_type, GREGORIAN, set_time,  &
                             print_date, print_time, operator(-)
use        model_mod, only : static_init_model, get_model_size, get_state_meta_data, &
                             model_interpolate, restart_file_to_sv
use     dart_cice_mod, only : get_cice_restart_filename

use netcdf
use typesizes

implicit none

! version controlled file description for error handling, do not edit
character(len=256), parameter :: source   = &
   "$URL: https://proxy.subversion.ucar.edu/DAReS/DART/branches/rma_cice/models/cice/model_mod_check.f90 $"
character(len=32 ), parameter :: revision = "$Revision: 8565 $"
character(len=128), parameter :: revdate  = "$Date: 2015-09-11 10:16:08 -0700 (Fri, 11 Sep 2015) $"

! standard output string
character(len=256) :: string1, string2

!------------------------------------------------------------------
! The namelist variables
!------------------------------------------------------------------

character(len=256)     :: dart_input_file      = 'dart_ics'
character(len=256)     :: output_file          = 'check_me'
logical                :: advance_time_present = .FALSE.
logical                :: verbose              = .FALSE.
integer                :: test1thru            = -1
real(r8)               :: interp_test_dlon     = 1.0
real(r8)               :: interp_test_dlat     = 1.0
real(r8), dimension(2) :: interp_test_latrange = (/ -90.0,  90.0 /)
real(r8), dimension(2) :: interp_test_lonrange = (/   0.0, 360.0 /)
real(r8), dimension(3) :: loc_of_interest     = -1.0_r8
character(len=metadatalength) :: kind_of_interest = 'ANY'
character(len=metadatalength) :: interp_test_vertcoord = 'VERTISLEVEL'

namelist /model_mod_check_nml/ dart_input_file, output_file, &
                        advance_time_present, test1thru, &
                        loc_of_interest, kind_of_interest, verbose, &
                        interp_test_dlon, interp_test_lonrange, &
                        interp_test_dlat, interp_test_latrange

!----------------------------------------------------------------------

integer :: ios_out, iunit, io, i
integer :: x_size, skip
integer :: mykindindex

type(time_type)       :: model_time, adv_to_time
real(r8), allocatable :: statevector(:)

character(len=metadatalength) :: state_meta(1)
character(len=129) :: cice_input_file  ! set with get_model_analysis_filename() if needed
type(netcdf_file_type) :: ncFileID
type(location_type) :: loc

real(r8) :: interp_val

!----------------------------------------------------------------------
! This portion checks the geometry information.
!----------------------------------------------------------------------

call initialize_utilities(progname='model_mod_check')
call set_calendar_type(GREGORIAN)

write(*,*)
write(*,*)'Reading the namelist to get the input filename.'

call find_namelist_in_file("input.nml", "model_mod_check_nml", iunit)
read(iunit, nml = model_mod_check_nml, iostat = io)
call check_namelist_read(iunit, io, "model_mod_check_nml")

! Record the namelist values used for the run
if (do_nml_file()) write(nmlfileunit, nml=model_mod_check_nml)
if (do_nml_term()) write(     *     , nml=model_mod_check_nml)

loc = set_location(loc_of_interest(1), loc_of_interest(2), loc_of_interest(3), &
                   VERTISLEVEL)

mykindindex = get_raw_obs_kind_index(kind_of_interest)

!----------------------------------------------------------------------
! This harvests all kinds of initialization information
!----------------------------------------------------------------------
if (test1thru > 0) then

   write(*,*)
   write(*,*)
   write(*,*)'------------------ Test #1 - static_init_model ---------------------'

   call static_init_model()

   write(*,*)'------------------ Test #1 complete --------------------------------'

endif

if (test1thru > 1) then

   write(*,*)
   write(*,*)
   write(*,*)'------------------ Test #2 - get_model_size ------------------------'

   x_size = get_model_size()

   write(*,'(''state vector has length'',i10)') x_size
   write(*,*)'------------------ Test #2 complete --------------------------------'

endif

!----------------------------------------------------------------------
! Write a supremely simple restart file. Most of the time, I just use
! this as a starting point for a Matlab function that replaces the
! values with something more complicated.
!----------------------------------------------------------------------

if (test1thru > 2) then

   write(*,*)
   write(*,*)
   write(*,*)'------------------ Test #3 - awrite_state_restart ------------------'

   allocate(statevector(x_size))

   call get_cice_restart_filename(cice_input_file)

   write(*,*)'Reading file from  '//trim(cice_input_file)
   write(*,*)'Writing data into '//trim(output_file)

   call restart_file_to_sv(cice_input_file, statevector, model_time)

   iunit = open_restart_write(output_file,'unformatted')
   call awrite_state_restart(model_time, statevector, iunit)
   call close_restart(iunit)

   call print_date( model_time,'model_mod_check:model date')
   call print_time( model_time,'model_mod_check:model time')

   write(*,*)'------------------ Test #3 complete --------------------------------'

endif

!----------------------------------------------------------------------
! Open a test DART initial conditions file.
! Reads the valid time, the state, and (possibly) a target time.
!----------------------------------------------------------------------

if (test1thru > 3) then

   write(*,*)
   write(*,*)
   write(*,*)'------------------ Test #4 - open_restart_read ---------------------'

   if ( .not. file_exist(dart_input_file) ) then
      write(string1,*) 'cannot open file [', trim(dart_input_file),'] for reading.'
      write(string2,*) "Run 'cice_to_dart' to generate ",trim(dart_input_file)
      call error_handler(E_ERR,'model_mod_check',string1, &
               source,revision,revdate,text2=string2)
   endif

   write(*,*)'Reading ['//trim(dart_input_file)//'] advance_time_present is ', &
             advance_time_present

   iunit = open_restart_read(dart_input_file)
   if ( advance_time_present ) then
      call aread_state_restart(model_time, statevector, iunit, adv_to_time)
   else
      call aread_state_restart(model_time, statevector, iunit)
   endif
   call close_restart(iunit)

   call print_date( model_time,'model_mod_check:model date')
   call print_time( model_time,'model_mod_check:model time')

   write(*,*)'------------------ Test #4 complete --------------------------------'

endif

!----------------------------------------------------------------------
! Output the state vector to a netCDF file ...
! This is the same procedure used by 'perfect_model_obs' & 'filter'
! init_diag_output()
! aoutput_diagnostics()
! finalize_diag_output()
!----------------------------------------------------------------------

if (test1thru > 4) then

   write(*,*)
   write(*,*)
   write(*,*)'------------------ Test #5 - Exercising the netCDF routines --------'

   write(*,*)'Creating '//trim(output_file)//'.nc'

   state_meta(1) = 'restart test'
   ncFileID = init_diag_output(trim(output_file),'just testing a restart', 1, state_meta)

   call aoutput_diagnostics(ncFileID, model_time, statevector, 1)

   call nc_check( finalize_diag_output(ncFileID), 'model_mod_check:main', 'finalize')

   write(*,*)'------------------ Test #5 complete --------------------------------'
endif

!----------------------------------------------------------------------
! Checking get_state_meta_data (and get_state_indices, get_state_kind)
!----------------------------------------------------------------------

if (test1thru > 5) then

   write(*,*)
   write(*,*)
   write(*,*)'------------------ Test #6 check_meta_data() -----------------------'

   skip = 6000000

   do i = 1, x_size, skip
      if ( i > 0 .and. i <= x_size ) call check_meta_data( i )
   enddo

   write(*,*)'------------------ Test #6 complete --------------------------------'

endif

!----------------------------------------------------------------------
! Trying to find the state vector index closest to a given location.
! Checking for valid input is tricky ... we don't know much.
!----------------------------------------------------------------------

if (test1thru > 6) then

   write(*,*)
   write(*,*)
   write(*,*)'------------------ Test #7 - find_closest_gridpoint() --------------'

   if ( loc_of_interest(1) > 0.0_r8 ) then
      write(*,*)'Skipping test because loc_of_interest not fully specified.'
   else
      call find_closest_gridpoint( loc_of_interest )
   endif

   write(*,*)'------------------ Test #7 complete --------------------------------'

endif

!----------------------------------------------------------------------
! Check the interpolation - print initially to STDOUT
!----------------------------------------------------------------------

if (test1thru > 7) then

   write(*,*)
   write(*,*)
   write(*,*)'------------------ Test #8 - single model_interpolate --------------'

   write(*,'('' KIND : '',A)') trim(kind_of_interest)
   write(*,'('' lon/lat/cat '',3(1x,f14.5))')loc_of_interest(1:LocationDims)

   vertcoord = VERTISLEVEL
   loc = set_location(loc_of_interest(1), loc_of_interest(2), loc_of_interest(3), &
                      vertcoord)

   mykindindex = get_raw_obs_kind_index(kind_of_interest)

   call model_interpolate(statevector, loc, mykindindex, interp_val, ios_out)

   if ( ios_out == 0 ) then
      write(*,*)'model_interpolate SUCCESS.'
      write(*,*)'value = ', interp_val
   else
      write(*,*)'model_interpolate WARNING: model_interpolate returned error code ',ios_out
   endif

   write(*,*)'------------------ Test #8 complete --------------------------------'

endif

if (test1thru > 8) then

   write(*,*)
   write(*,*)
   write(*,*)'------------------ Test #9 - Rigorous test of model_interpolate ----'

   ios_out = test_interpolate()

   if ( ios_out == 0 ) then
      write(*,*)'Rigorous model_interpolate SUCCESS.'
   else
      write(*,*)'Rigorous model_interpolate WARNING: model_interpolate had ', ios_out, ' failures.'
   endif

   write(*,*)'------------------ Test #9 complete --------------------------------'

endif

!----------------------------------------------------------------------
! This must be the last few lines of the main program.
!----------------------------------------------------------------------

call finalize_utilities()

!----------------------------------------------------------------------
!----------------------------------------------------------------------

contains


subroutine check_meta_data( iloc )

integer, intent(in) :: iloc
type(location_type) :: loc
integer             :: var_type

! write(*,*)'Checking metadata routines.'

call get_state_meta_data( iloc, loc, var_type)

call write_location(0, loc, fform='formatted', charstring=string1)
write(*,*)' indx ',iloc,' is type ',var_type,' ', &
         trim(get_raw_obs_kind_name(var_type))
write(*,*)' and is at ',trim(string1)
write(*,*)

end subroutine check_meta_data



subroutine find_closest_gridpoint( loc_of_interest )
! Simple exhaustive search to find the indices into the
! state vector of a particular lon/lat/cat. They will
! occur multiple times - once for each state variable.

real(r8), dimension(:), intent(in) :: loc_of_interest

type(location_type) :: loc0, loc1
integer  :: i, var_type, which_cat
real(r8) :: closest, rlon, rlat, rcat, vals(3)
real(r8), allocatable, dimension(:) :: thisdist
character(len=32) :: kind_name
logical :: matched

! Check user input ... if there is no 'catical' ...
if ( (count(loc_of_interest >= 0.0_r8) < 3) .or. &
     (LocationDims < 3 ) ) then
   write(*,*)
   write(*,*)'Interface not fully implemented.'
   return
endif

write(*,'(''Checking for the indices into the state vector that are at'')')
write(*,'(''lon/lat/cat'',3(1x,f14.5))')loc_of_interest(1:LocationDims)

allocate( thisdist(get_model_size()) )
thisdist  = 9999999999.9_r8         ! really far away
matched   = .false.

! Trying to support the ability to specify matching a particular KIND.
! With staggered grids, the closest gridpoint might not be of the kind
! you are interested in. mykindindex = -1 means anything will do.

rlon = loc_of_interest(1)
rlat = loc_of_interest(2)
rcat = loc_of_interest(3)

! Since there can be/will be multiple variables with
! identical distances, we will just cruise once through
! the array and come back to find all the 'identical' values.
do i = 1,get_model_size()

   ! Really inefficient, but grab the 'which_cat' from the
   ! grid and set our target location to have the same.
   ! Then, compute the distance and compare.

   call get_state_meta_data(i, loc1, var_type)

   if ( (var_type == mykindindex) .or. (mykindindex < 0) ) then
      which_cat  = nint( query_location(loc1) )
      loc0        = set_location(rlon, rlat, rcat, which_cat)
      thisdist(i) = get_dist( loc1, loc0, no_vert= .false. )
      matched     = .true.
   endif

enddo

if (.not. matched) then
   write(*,*)'No state vector elements of type '//trim(kind_of_interest)
   return
endif

! Now that we know the distances ... report

closest = minval(thisdist)
if (closest == 9999999999.9_r8) then
   write(*,*)'No closest gridpoint found'
   return
endif


matched = .false.
do i = 1,get_model_size()

   if ( thisdist(i) == closest ) then
      call get_state_meta_data(i, loc1, var_type)
      kind_name = get_raw_obs_kind_name(var_type)
      vals = get_location(loc1)
      write(*,'(''lon/lat/cat'',3(1x,f14.5),'' is index '',i10,'' for '',a)') &
             vals, i, trim(kind_name)
      matched = .true.
   endif

enddo

if ( .not. matched ) then
   write(*,*)'Nothing matched the closest gridpoint'
endif

deallocate( thisdist )

end subroutine find_closest_gridpoint



function test_interpolate()
! function to exercise the model_mod:model_interpolate() function
! This will result in a netCDF file with all salient metadata
! CMB not sure what to do about ncats
! can I just deal with the aggregate quantities?
! writing as if field is 2D but may need to add cat as a dim
integer :: test_interpolate

! Local variables

real(r8), allocatable :: lon(:), lat(:)
real(r8), allocatable :: field(:,:)
integer :: nlon, nlat, 
integer :: ilon, jlat, nfailed, nland
character(len=128) :: ncfilename,txtfilename

character(len=8)      :: crdate      ! needed by F90 DATE_AND_TIME intrinsic
character(len=10)     :: crtime      ! needed by F90 DATE_AND_TIME intrinsic
character(len=5)      :: crzone      ! needed by F90 DATE_AND_TIME intrinsic
integer, dimension(8) :: values      ! needed by F90 DATE_AND_TIME intrinsic

integer :: ncid, nlonDimID, nlatDimID
integer :: VarID, lonVarID, latVarID

test_interpolate = 0   ! normal termination

if ((interp_test_dlon < 0.0_r8) .or. (interp_test_dlat < 0.0_r8)) then
   write(*,*)'Skipping the rigorous interpolation test because one of'
   write(*,*)'interp_test_dlon,interp_test_dlat are < 0.0'
   write(*,*)'interp_test_dlon  = ',interp_test_dlon
   write(*,*)'interp_test_dlat  = ',interp_test_dlat
   return
endif

write( ncfilename,'(a,a)')trim(output_file),'_interptest.nc'
write(txtfilename,'(a,a)')trim(output_file),'_interptest.m'

! round down to avoid exceeding the specified range
nlat  = aint(( interp_test_latrange(2) -  interp_test_latrange(1))/interp_test_dlat) + 1
nlon  = aint(( interp_test_lonrange(2) -  interp_test_lonrange(1))/interp_test_dlon) + 1

iunit = open_file(trim(txtfilename), action='write')
write(iunit,'(''missingvals = '',f12.4,'';'')')MISSING_R8
write(iunit,'(''nlon = '',i8,'';'')')nlon
write(iunit,'(''nlat = '',i8,'';'')')nlat
write(iunit,'(''interptest = [ ... '')')

allocate(lon(nlon), lat(nlat), field(nlon,nlat))
nfailed = 0
nland = 0

do ilon = 1, nlon
   lon(ilon) = interp_test_lonrange(1) + real(ilon-1,r8) * interp_test_dlon
   do jlat = 1, nlat
      lat(jlat) = interp_test_latrange(1) + real(jlat-1,r8) * interp_test_dlat

         loc = set_location(lon(ilon), lat(jlat))

         call model_interpolate(statevector, loc, mykindindex, field(ilon,jlat), ios_out)
         write(iunit,*) field(ilon,jlat)

         ! Interpolation location was over land:
         ! ios_out = 1 . no regular boxes contain wet cell
         ! ios_out = 3 . location over land
         if (ios_out == 1 .or. ios_out == 3) then
            nland = nland + 1
         else if (ios_out /= 0) then
           if (verbose) then
              write(string2,'(''ilon,jlat,lon,lat'',2(1x,i6),2(1x,f14.6))') &
                          ilon,jlat,lon(ilon),lat(jlat)
              write(string1,*) 'interpolation return code was', ios_out
              call error_handler(E_MSG,'test_interpolate',string1,source,revision,revdate,text2=string2)
           endif
           nfailed = nfailed + 1
         endif

   end do
end do

write(iunit,'(''];'')')
write(iunit,'(''datmat = reshape(interptest,nlat,nlon);'')')
write(iunit,'(''datmat(datmat == missingvals) = NaN;'')')
call close_file(iunit)

if (nlon*nlat*ncat - nland == 0) then
   write(*,*)'all grid points in :'
   write(*,*)'   interp_test_dlon  = ',interp_test_dlon
   write(*,*)'   interp_test_dlat  = ',interp_test_dlat
   write(*,*)'are over land.  Modify interp_test range to include points with water.'
else
   write(*,*) 'total interpolations  : ', nlon*nlat - nland
   write(*,*) 'failed interpolations : ', nfailed
endif

! Write out the netCDF file for easy exploration.

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

call nc_check(nf90_def_var(ncid=ncid, name='field', xtype=nf90_double, &
        dimids=(/ nlonDimID, nlatDimID /), varid=VarID), 'test_interpolate', &
                 'field def_var '//trim(ncfilename))
call nc_check(nf90_put_att(ncid, VarID, 'long_name', kind_of_interest), &
           'test_interpolate', 'put_att field long_name '//trim(ncfilename))
call nc_check(nf90_put_att(ncid, VarID, '_FillValue', MISSING_R8), &
           'test_interpolate', 'put_att field FillValue '//trim(ncfilename))
call nc_check(nf90_put_att(ncid, VarID, 'missing_value', MISSING_R8), &
           'test_interpolate', 'put_att field missing_value '//trim(ncfilename))
call nc_check(nf90_put_att(ncid, VarID, 'vertcoord', vertcoord ), &
           'test_interpolate', 'put_att field vertcoord '//trim(ncfilename))

! Leave define mode so we can fill the variables.
call nc_check(nf90_enddef(ncid), &
              'test_interpolate','field enddef '//trim(ncfilename))

! Fill the variables
call nc_check(nf90_put_var(ncid, lonVarID, lon), &
              'test_interpolate','lon put_var '//trim(ncfilename))
call nc_check(nf90_put_var(ncid, latVarID, lat), &
              'test_interpolate','lat put_var '//trim(ncfilename))
call nc_check(nf90_put_var(ncid, VarID, field), &
              'test_interpolate','field put_var '//trim(ncfilename))

! tidy up
call nc_check(nf90_close(ncid), &
             'test_interpolate','close '//trim(ncfilename))

deallocate(lon, lat, field)

test_interpolate = nfailed

end function test_interpolate


end program model_mod_check

! <next few lines under version control, do not edit>
! $URL: https://proxy.subversion.ucar.edu/DAReS/DART/branches/rma_cice/models/cice/model_mod_check.f90 $
! $Id: model_mod_check.f90 8565 2015-09-11 17:16:08Z hkershaw $
! $Revision: 8565 $
! $Date: 2015-09-11 10:16:08 -0700 (Fri, 11 Sep 2015) $
