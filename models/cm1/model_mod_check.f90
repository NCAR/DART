! DART software - Copyright 2004 - 2013 UCAR. This open source software is
! provided by UCAR, "as is", without charge, subject to all terms of use at
! http://www.image.ucar.edu/DAReS/DART/DART_download
!
! $Id$

program model_mod_check

!----------------------------------------------------------------------
! purpose: test routines
!----------------------------------------------------------------------

use             types_mod, only : r8, i8, digits12, metadatalength, missing_r8
use         utilities_mod, only : register_module, find_namelist_in_file, &
                                  check_namelist_read, error_handler, E_MSG, &
                                  open_file, close_file, nc_check
use     mpi_utilities_mod, only : initialize_mpi_utilities, finalize_mpi_utilities
use          location_mod, only : location_type, set_location, write_location, get_dist, &
                                  query_location, LocationDims, get_location
use          obs_kind_mod         ! no 'only' so all kinds available
!use       assim_model_mod, only : init_diag_output, netcdf_file_type, static_init_assim_model
use       assim_model_mod, only : static_init_assim_model
use      time_manager_mod, only : time_type, set_calendar_type, GREGORIAN, &
                                  read_time, get_time, set_time,  &
                                  print_date, get_date, &
                                  print_time, write_time, &
                                  operator(-)
use  ensemble_manager_mod, only : init_ensemble_manager,               &
                                  end_ensemble_manager, ensemble_type,  &
                                  get_my_num_copies, get_ensemble_time, prepare_to_write_to_vars,      &
                                  prepare_to_read_from_vars, &
                                  all_vars_to_all_copies, &
                                  all_copies_to_all_vars, &
                                  copies_in_window, set_num_extra_copies
use            filter_mod, only : filter_set_initial_time, filter_sync_keys_time
use   state_vector_io_mod, only : state_vector_io_init, &
                                  setup_read_write, turn_read_copy_on, turn_write_copy_on,&
                                  filter_read_restart_direct, &
                                  filter_write_restart_direct
use distributed_state_mod, only : create_state_window, free_state_window
use   state_structure_mod, only : get_varid_from_kind, get_num_domains, &
                                  static_init_state_type
use      io_filenames_mod, only : io_filenames_init, get_input_file, set_filenames
use   quality_control_mod, only : set_input_qc, initialize_qc
use             model_mod, only : static_init_model, get_model_size, get_state_meta_data, &
                                  model_interpolate

use typesizes
use netcdf

implicit none

! version controlled file description for error handling, do not edit
character(len=256), parameter :: source   = &
   "$URL$"
character(len=32 ), parameter :: revision = "$Revision$"
character(len=128), parameter :: revdate  = "$Date$"

!------------------------------------------------------------------
! The namelist variables
!------------------------------------------------------------------

character (len=256)           :: input_file  = 'dart_ics.nc'
character (len=256)           :: output_file = 'check_me.nc'
integer(i8)                   :: x_ind   = -1
integer                       :: num_ens = 1
real(r8), dimension(3)        :: loc_of_interest = -1.0_r8
character(len=metadatalength) :: kind_of_interest = 'ANY'
logical                       :: verbose = .FALSE.
integer                       :: run_test = 1
real(r8)               :: interp_test_dy     = 10.0
real(r8)               :: interp_test_dx     = 10.0
real(r8)               :: interp_test_dz    = 10.0
real(r8), dimension(2) :: interp_test_yrange = (/   0.0,  120.0 /)
real(r8), dimension(2) :: interp_test_xrange = (/   0.0,  120.0 /)
real(r8), dimension(2) :: interp_test_zrange = (/  0.0,  100.0 /)

namelist /model_mod_check_nml/ input_file, output_file, &
                        x_ind, num_ens, loc_of_interest, &
                        kind_of_interest, verbose, run_test, &
                        interp_test_dy, interp_test_xrange, &
                        interp_test_dx, interp_test_yrange, &
                        interp_test_dz, interp_test_zrange

!----------------------------------------------------------------------

!> @TODO:
!> ens_size could be in the namelist, and we could use the perturb
!> routines to generate an ensemble for testing.  hardcode to 1 for now.
integer :: ens_size = 1

integer :: iloc
!type(netcdf_file_type) :: StateUnit
type(time_type) :: time1
logical :: read_time_from_file
integer :: init_time_days, init_time_seconds
integer :: iunit, io
integer :: num_failed
integer, allocatable :: ios_out(:)
integer(i8) :: model_size

character(len=512) :: string1, string2

type(time_type) :: model_time

character(len=metadatalength) :: state_meta(1)
type(ensemble_type) :: ens_handle
type(location_type) :: loc

real(r8), allocatable :: interp_vals(:)

integer :: mykindindex

!----------------------------------------------------------------------
! This portion checks the geometry information.
!----------------------------------------------------------------------

write(*,*)
write(*,*) '--------------- RUNNING TEST 1 ----------------'
write(*,*) 'Reading the namelist to get the input filename.'
write(*,*) '-----------------------------------------------'
write(*,*)

call initialize_modules_used()

call set_calendar_type(GREGORIAN)

call find_namelist_in_file("input.nml", "model_mod_check_nml", iunit)
read(iunit, nml = model_mod_check_nml, iostat = io)
call check_namelist_read(iunit, io, "model_mod_check_nml")

! This harvests all kinds of initialization information
call static_init_model()

model_size = get_model_size()
write(*,'(''state vector has length'',i10)') model_size

!----------------------------------------------------------------------
! Write a supremely simple restart file. Most of the time, I just use
! this as a starting point for a Matlab function that replaces the
! values with something more complicated.
!----------------------------------------------------------------------

write(*,*)
write(*,*) '--------------- FINISHED TEST 1 ----------------'
write(*,*)

if ( run_test == 1 ) then 
   call exit(0)
endif

write(*,*)
write(*,*) '--------------- RUNNING TEST 2 ----------------'
write(*,*) ' Read and write trivial restart file.'
write(*,*) '-----------------------------------------------'
write(*,*)

model_time  = set_time(21600, 149446)   ! 06Z 4 March 2010

! Set up the ensemble storage and read in the restart file
call init_ensemble_manager(ens_handle, num_ens, model_size, 1)

! Reading netcdf restart file:
call setup_read_write(1)
call set_filenames(ens_handle, num_ens, "no_inf1", "no_inf2")

!----------------------------------------------------------------------
! Open a test netcdf initial conditions file.
!----------------------------------------------------------------------

call turn_read_copy_on(1)

call filter_set_initial_time(init_time_days, init_time_seconds, time1, read_time_from_file)
call filter_read_restart_direct(ens_handle, time1, 0, read_time_from_file)

call turn_write_copy_on(1)
call filter_write_restart_direct(ens_handle, 0, isprior=.true.)

model_time = time1
call print_date( model_time,'model_mod_check:model date')
call print_time( model_time,'model_mod_check:model time')

!----------------------------------------------------------------------
! Output the state vector to a netCDF file ...
! This is the same procedure used by 'perfect_model_obs' & 'filter'
! init_diag_output()
! aoutput_diagnostics()
! finalize_diag_output()
!
! FIXME: is this still true?
!----------------------------------------------------------------------

write(*,*)
write(*,*) '--------------- FINISHED TEST 2 ----------------'
write(*,*)

if ( run_test == 2 ) then 
   call exit(0)
endif

!! Set up output of truth for state
!state_meta(1) = 'true state'
!if (ens_handle%my_pe == 0) then
!   StateUnit = init_diag_output('True_State', 'true state from control', 1, state_meta)
!endif

!----------------------------------------------------------------------
! Check the meta data
!----------------------------------------------------------------------

write(*,*) '--------------- RUNNING TEST 3 ----------------'
write(*,*) ' Testing get_state_meta_data                   '
write(*,*) '-----------------------------------------------'

if ( x_ind > 0 .and. x_ind <= model_size ) call check_meta_data( x_ind )

write(*,*)
write(*,*) '--------------- FINISHED TEST 3 ----------------'
write(*,*)

if ( run_test == 3 ) then 
   call exit(0)
endif

!----------------------------------------------------------------------
! Check the interpolation - print initially to STDOUT
!----------------------------------------------------------------------

write(*,*)
write(*,*) '--------------- RUNNING TEST 4 ----------------'
write(*,*) ' Testing model_interpolate                     '
write(*,*) '-----------------------------------------------'
write(*,*)

! Create window for forward operators
call set_num_extra_copies(ens_handle, 0)
call create_state_window(ens_handle)

! FIXME: needs to be ensemble_size
loc = set_location(loc_of_interest(1), loc_of_interest(2), loc_of_interest(3))

mykindindex = get_raw_obs_kind_index(kind_of_interest)

allocate(interp_vals(1), ios_out(1))
call model_interpolate(ens_handle, ens_size, loc, mykindindex, interp_vals, ios_out)

if ( all(ios_out(:) == 0) ) then
   write(*,*)'model_interpolate SUCCESS: The interpolated value(s) are ',interp_vals(:)
else
   write(*,*)'model_interpolate ERROR: model_interpolate failed with error code(s) ',ios_out(:)
endif

write(*,*)
write(*,*) '--------------- FINISHED TEST 4 ----------------'
write(*,*)

if ( run_test == 4 ) then 
   call exit(0)
endif

write(*,*)
write(*,*) '--------------- RUNNING TEST 5 ----------------'
write(*,*) ' Testing model_interpolate                     '
write(*,*) '-----------------------------------------------'
write(*,*)

num_failed = test_interpolate()

write(*,*)
write(*,*) '--------------- FINISHED TEST 5 ----------------'
write(*,*)

! finalize model_mod_check
call error_handler(E_MSG,'full_model_mod_check','Finished successfully.',source,revision,revdate)
call finalize_mpi_utilities()

!----------------------------------------------------------------------

contains

!----------------------------------------------------------------------

subroutine check_meta_data( iloc )

integer(i8), intent(in) :: iloc
type(location_type) :: loc
integer             :: var_type
character(len=129)  :: string1

write(*,*)
write(*,*)'Checking metadata routines.'

call get_state_meta_data(ens_handle, iloc, loc, var_type)

call write_location(42, loc, fform='formatted', charstring=string1)
write(*,*)' indx ',iloc,' is type ',var_type,trim(string1)

end subroutine check_meta_data


!----------------------------------------------------------------------

! FIXME: should this all be wrapped into a 'initialize_dart()' routine?

subroutine initialize_modules_used()

! Standard initialization (mpi not needed to use ensemble manager
! since we are enforcing that this run as a single task).
call initialize_mpi_utilities('model_mod_check')

! Initialize modules used that require it
call register_module(source,revision,revdate)

call static_init_assim_model()
call state_vector_io_init()
call static_init_state_type()
!call io_filenames_init()
call initialize_qc()

end subroutine initialize_modules_used
!----------------------------------------------------------------------

function test_interpolate()
! function to exercise the model_mod:model_interpolate() function
! This will result in a netCDF file with all salient metadata 
integer :: test_interpolate

! Local variables

real(r8), allocatable :: lon(:), lat(:), vert(:)
real(r8), allocatable :: field(:,:,:,:)
integer :: nlon, nlat, nvert
integer :: ilon, jlat, kvert, nfailed
character(len=128) :: ncfilename,txtfilename

character(len=8)      :: crdate      ! needed by F90 DATE_AND_TIME intrinsic
character(len=10)     :: crtime      ! needed by F90 DATE_AND_TIME intrinsic
character(len=5)      :: crzone      ! needed by F90 DATE_AND_TIME intrinsic
integer, dimension(8) :: values      ! needed by F90 DATE_AND_TIME intrinsic

integer :: ncid, nlonDimID, nlatDimID, nvertDimID
integer :: VarID, lonVarID, latVarID, vertVarID

test_interpolate = 0   ! normal termination

if ((interp_test_dy < 0.0_r8) .or. (interp_test_dx < 0.0_r8)) then
   write(*,*)'Skipping the rigorous interpolation test because one of'
   write(*,*)'interp_test_dy,interp_test_dx are < 0.0'
   write(*,*)'interp_test_dy  = ',interp_test_dy
   write(*,*)'interp_test_dx  = ',interp_test_dx
   write(*,*)'interp_test_dz = ',interp_test_dz
endif

write( ncfilename,'(a,a)')trim(output_file),'_interptest.nc'
write(txtfilename,'(a,a)')trim(output_file),'_interptest.m'

! round down to avoid exceeding the specified range
nlat  = aint(( interp_test_yrange(2) -  interp_test_yrange(1))/interp_test_dx) + 1
nlon  = aint(( interp_test_xrange(2) -  interp_test_xrange(1))/interp_test_dy) + 1
nvert = aint((interp_test_zrange(2) - interp_test_zrange(1))/interp_test_dz) + 1

iunit = open_file(trim(txtfilename), action='write')
write(iunit,'(''missingvals = '',f12.4,'';'')')MISSING_R8
write(iunit,'(''nlon = '',i8,'';'')')nlon
write(iunit,'(''nlat = '',i8,'';'')')nlat
write(iunit,'(''nvert = '',i8,'';'')')nvert
write(iunit,'(''interptest = [ ... '')')

allocate(lon(nlon), lat(nlat), vert(nvert), field(nlon,nlat,nvert,1))
nfailed = 0

do ilon = 1, nlon
   lon(ilon) = interp_test_xrange(1) + real(ilon-1,r8) * interp_test_dy
   do jlat = 1, nlat
      lat(jlat) = interp_test_yrange(1) + real(jlat-1,r8) * interp_test_dx
      do kvert = 1, nvert
         vert(kvert) = interp_test_zrange(1) + real(kvert-1,r8) * interp_test_dz

         loc = set_location(lon(ilon), lat(jlat), vert(kvert))

         call model_interpolate(ens_handle, ens_size, loc, mykindindex, field(ilon,jlat,kvert,:), ios_out)
         call model_interpolate(ens_handle, ens_size, loc, mykindindex, interp_vals, ios_out)
         ! call model_interpolate(statevector, ens_size, loc, mykindindex, field(ilon,jlat,kvert), ios_out)
         write(iunit,*) field(ilon,jlat,kvert,1)

         if (any(ios_out(:) /= 0)) then
           if (verbose) then
              write(string2,'(''ilon,jlat,kvert,lon,lat,vert'',3(1x,i6),3(1x,f14.6))') &
                          ilon,jlat,kvert,lon(ilon),lat(jlat),vert(kvert)
              write(string1,*) 'interpolation return code was', ios_out
              call error_handler(E_MSG,'test_interpolate',string1,source,revision,revdate,text2=string2)
           endif
           nfailed = nfailed + 1
         endif

      enddo
   end do 
end do
write(iunit,'(''];'')')
write(iunit,'(''datmat = reshape(interptest,nvert,nlat,nlon);'')')
write(iunit,'(''datmat(datmat == missingvals) = NaN;'')')
call close_file(iunit)

write(*,*) 'total interpolations, num failed: ', nlon*nlat*nvert, nfailed

! Write out the netCDF file for easy exploration.

call DATE_AND_TIME(crdate,crtime,crzone,values)
write(string1,'(''YYYY MM DD HH MM SS = '',i4,5(1x,i2.2))') &
                  values(1), values(2), values(3), values(5), values(6), values(7)

call nc_check( nf90_create(path=trim(ncfilename), cmode=NF90_clobber, ncid=ncid), &
                  'test_interpolate', 'open '//trim(ncfilename))
call nc_check( nf90_put_att(ncid, NF90_GLOBAL, 'creation_date' ,trim(string1) ), &
                  'test_interpolate', 'creation put '//trim(ncfilename))
! call nc_check( nf90_put_att(ncid, NF90_GLOBAL, 'parent_file', dart_input_file ), &
!                   'test_interpolate', 'put_att filename '//trim(ncfilename))

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
call nc_check(nf90_put_att(ncid, lonVarID, 'range', interp_test_xrange), &
           'test_interpolate', 'put_att xrange '//trim(ncfilename))
call nc_check(nf90_put_att(ncid, lonVarID, 'cartesian_axis', 'X'),   &
           'test_interpolate', 'lon cartesian_axis '//trim(ncfilename))


call nc_check(nf90_def_var(ncid=ncid, name='lat', xtype=nf90_double, &
        dimids=nlatDimID, varid=latVarID), 'test_interpolate', &
                 'lat def_var '//trim(ncfilename))
call nc_check(nf90_put_att(ncid, latVarID, 'range', interp_test_yrange), &
           'test_interpolate', 'put_att yrange '//trim(ncfilename))
call nc_check(nf90_put_att(ncid, latVarID, 'cartesian_axis', 'Y'),   &
           'test_interpolate', 'lat cartesian_axis '//trim(ncfilename))

call nc_check(nf90_def_var(ncid=ncid, name='vert', xtype=nf90_double, &
        dimids=nvertDimID, varid=vertVarID), 'test_interpolate', &
                 'vert def_var '//trim(ncfilename))
call nc_check(nf90_put_att(ncid, vertVarID, 'range', interp_test_zrange), &
           'test_interpolate', 'put_att zrange '//trim(ncfilename))
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
! call nc_check(nf90_put_att(ncid, VarID, 'vertcoord_string', interp_test_vertcoord ), &
!            'test_interpolate', 'put_att field vertcoord_string '//trim(ncfilename))
! call nc_check(nf90_put_att(ncid, VarID, 'vertcoord', vertcoord ), &
!            'test_interpolate', 'put_att field vertcoord '//trim(ncfilename))

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
call nc_check(nf90_put_var(ncid, VarID, field(:,:,:,1)), &
              'test_interpolate','field put_var '//trim(ncfilename))

! tidy up
call nc_check(nf90_close(ncid), &
             'test_interpolate','close '//trim(ncfilename))

deallocate(lon, lat, vert, field)

test_interpolate = nfailed

end function test_interpolate
!----------------------------------------------------------------------


end program model_mod_check

! <next few lines under version control, do not edit>
! $URL$
! $Id$
! $Revision$
! $Date$
