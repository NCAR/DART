! This code may (or may not) be part of the NOGAPS distribution,
! So it is not protected by the DART copyright agreement.
!
! DART $Id$

program check_me

!----------------------------------------------------------------------
! purpose: test routines ...
!----------------------------------------------------------------------

use        types_mod, only : r8, digits12, metadatalength
use    utilities_mod, only : initialize_utilities, timestamp, nc_check, &
                             open_file, close_file, find_namelist_in_file, &
                             check_namelist_read
use     location_mod, only : location_type, set_location, write_location, get_dist, &
                             query_location, LocationDims, get_location
use     obs_kind_mod, only : get_raw_obs_kind_name
use  assim_model_mod, only : open_restart_read, open_restart_write, close_restart, &
                             aread_state_restart, awrite_state_restart, &
                             netcdf_file_type, aoutput_diagnostics, &
                             init_diag_output, finalize_diag_output
use time_manager_mod, only : time_type, set_calendar_type, GREGORIAN, &
                             read_time, get_time, set_time,  &
                             print_date, get_date, &
                             print_time, write_time, &
                             operator(-)
use        model_mod, only : static_init_model, get_model_size, get_gridsize, &
                             test_interpolate, get_state_meta_data

implicit none

! This is from the original assembla server we used during collaboration.
! character(len=128), parameter :: &
!    source   = "$orgURL: https://svn2.assembla.com/svn/ngdart/check_me.f90 $", &
!    revision = "$orgRevision: 63 $", &
!    revdate  = "$orgDate: 2010-03-05 11:28:41 -0700 (Fri, 05 Mar 2010) $"

! version controlled file description for error handling, do not edit
character(len=256), parameter :: source   = &
   "$URL$"
character(len=32 ), parameter :: revision = "$Revision$"
character(len=128), parameter :: revdate  = "$Date$"

!------------------------------------------------------------------
! The namelist variables
!------------------------------------------------------------------

character (len = 128) :: input_file  = 'dart.ics'
character (len = 128) :: output_file = 'check_me'
logical               :: advance_time_present = .FALSE.
logical               :: verbose              = .FALSE.
integer               :: x_ind = -1
real(r8), dimension(4) :: LocOfInterest = -1.0_r8

namelist /check_me_nml/ input_file, output_file, &
                        advance_time_present, verbose, &
                        x_ind, LocOfInterest

!----------------------------------------------------------------------

integer :: in_unit, out_unit, ios_out, iunit, io, offset
integer :: numlons, numlats, numlevs, x_size
integer :: year, month, day, hour, minute, second
integer :: secs, days

type(time_type)       :: model_time, adv_to_time
real(r8), allocatable :: statevector(:)

character(len=metadatalength) :: state_meta(1)
type(netcdf_file_type) :: ncFileID

!----------------------------------------------------------------------
! This portion checks the geometry information. 
!----------------------------------------------------------------------

call initialize_utilities(progname='check_me', output_flag=verbose)
call set_calendar_type(GREGORIAN)

write(*,*)
write(*,*)'Reading the namelist to get the input filename.'

call find_namelist_in_file("input.nml", "check_me_nml", iunit)
read(iunit, nml = check_me_nml, iostat = io)
call check_namelist_read(iunit, io, "check_me_nml")

write(*,'(''Converting DART file '',A,'' to restart file '',A)') &
     trim(input_file), trim(output_file)

! This harvests all kinds of initialization information
call static_init_model()
x_size = get_model_size()
call get_gridsize(numlons, numlats, numlevs)

write(*,'(''nlons, nlats, nlevs'',3(1x,i10))') numlons,numlats,numlevs
write(*,'(''state vector has length'',i10)') x_size

allocate(statevector(x_size))

!----------------------------------------------------------------------
! Write a supremely simple restart file. Most of the time, I just use
! this as a starting point for a Matlab function that replaces the 
! values with something more complicated.
!----------------------------------------------------------------------

write(*,*)
write(*,*)'Writing a trivial restart file.'

statevector = 1.0_r8;
model_time  = set_time(21600, 149446)   ! 06Z 4 March 2010

iunit = open_restart_write('allones.ics')
call awrite_state_restart(model_time, statevector, iunit)
call close_restart(iunit)

!----------------------------------------------------------------------
! Open the test DART initial conditions file (augmented in Matlab).
! Reads the valid time, the state, and (possibly) a target time.
!----------------------------------------------------------------------

write(*,*)
write(*,*)'Reading the restart file '//trim(input_file)

iunit = open_restart_read(input_file)
if ( advance_time_present ) then
   call aread_state_restart(model_time, statevector, iunit, adv_to_time)
else
   call aread_state_restart(model_time, statevector, iunit)
endif

call close_restart(iunit)
call print_date( model_time,'check_me:model date')
call print_time( model_time,'check_me:model time')

!----------------------------------------------------------------------
! Output the state vector to a netCDF file ...
! This is the same procedure used by 'perfect_model_obs' & 'filter'
! init_diag_output()
! aoutput_diagnostics()
! finalize_diag_output()
!----------------------------------------------------------------------

write(*,*)
write(*,*)'Exercising the netCDF routines.'
write(*,*)'Creating '//trim(output_file)//'.nc'

state_meta(1) = 'restart test'
ncFileID = init_diag_output(trim(output_file),'just testing a restart', 1, state_meta)

call aoutput_diagnostics(ncFileID, model_time, statevector, 1)

call nc_check( finalize_diag_output(ncFileID), 'check_me:main', 'finalize')

!----------------------------------------------------------------------
! Check the interpolation - print initially to STDOUT
!----------------------------------------------------------------------
! Use 500 mb ~ level 30 (in the 42-level version)

write(*,*)
write(*,*)'Testing the interpolation ...'

call test_interpolate(statevector, test_pressure=500.0_r8, &
                      start_lon=142.5_r8)

!----------------------------------------------------------------------
! Checking get_state_meta_data (and get_state_indices, get_state_kind)
! nx = 144; ny=72; nz=42; produce the expected values :
!  U(       1 :  435456)
!  V(  435457 :  870912)
!  T(  870913 : 1306368)
!  Q( 1306369 : 1741824)
! PS( 1741825 : 1752193)    (only 144x72)
!----------------------------------------------------------------------

if ( x_ind > 0 .and. x_ind <= x_size ) call check_meta_data( x_ind )

!----------------------------------------------------------------------
! Trying to find the state vector index closest to a particular ...
! Checking for valid input is tricky ... we don't know much. 
!----------------------------------------------------------------------

if ( LocOfInterest(1) > 0.0_r8 ) call find_closest_gridpoint( LocOfInterest )

!----------------------------------------------------------------------
! When called with 'end', timestamp will call finalize_utilities()
! This must be the last few lines of the main program.
!----------------------------------------------------------------------
call timestamp(string1=source, pos='end')

contains


subroutine check_meta_data( iloc )

integer, intent(in) :: iloc
type(location_type) :: loc
integer             :: var_type
character(len=129)  :: string1

write(*,*)
write(*,*)'Checking metadata routines.'

call get_state_meta_data( iloc, loc, var_type)

call write_location(42, loc, fform='formatted', charstring=string1)
write(*,*)' indx ',iloc,' is type ',var_type,trim(string1)

end subroutine check_meta_data



subroutine find_closest_gridpoint( LocOfInterest )
! Simple exhaustive search to find the indices into the 
! state vector of a particular lon/lat/level. They will 
! occur multiple times - once for each state variable.
real(r8), dimension(:), intent(in) :: LocOfInterest

type(location_type) :: loc0, loc1
integer  :: imclosest
integer  :: i, var_type, which_vert
real(r8) :: closest, rlon, rlat, rlev
real(r8), allocatable, dimension(:) :: thisdist
real(r8), dimension(LocationDims) :: rloc
character(len=32) :: kind_name
logical :: matched

! Check user input ... if there is no 'vertical' ...  
if ( (count(LocOfInterest >= 0.0_r8) < 3) .or. &
     (LocationDims < 3 ) ) then
   write(*,*)
   write(*,*)'Interface not fully implemented.' 
   return
endif

write(*,*)
write(*,'(''Checking for the indices into the state vector that are at'')')
write(*,'(''lon/lat/lev'',3(1x,f10.5))')LocOfInterest(1:LocationDims)

allocate( thisdist(get_model_size()) )
thisdist  = 9999999999.9_r8         ! really far away 
imclosest = get_model_size() + 1    ! out of range
matched   = .false.

rlon = LocOfInterest(1)
rlat = LocOfInterest(2)
rlev = LocOfInterest(3)

! Since there can be/will be multiple variables with
! identical distances, we will just cruise once through 
! the array and come back to find all the 'identical' values.
do i = 1,get_model_size()

   ! Really inefficient, but grab the 'which_vert' from the
   ! grid and set our target location to have the same.
   ! Then, compute the distance and compare.

   call get_state_meta_data(i, loc1, var_type)

   which_vert  = nint( query_location(loc1) )
   loc0        = set_location(rlon, rlat, rlev, which_vert)
   thisdist(i) = get_dist( loc1, loc0, no_vert= .true. )

enddo

closest = minval(thisdist)

! Now that we know the distances ... report 

do i = 1,get_model_size()

   if ( thisdist(i) == closest ) then
      call get_state_meta_data(i, loc1, var_type)
      rloc      = get_location(loc1)
      if (nint(rloc(3)) == nint(rlev)) then
         kind_name = get_raw_obs_kind_name(var_type)
         write(*,'(''lon/lat/lev'',3(1x,f10.5),'' is index '',i10,'' for '',a)') &
             rloc, i, trim(kind_name)
         matched = .true.
      endif
   endif

enddo

if ( .not. matched ) then
   write(*,*)'Nothing matched the vertical.'
endif

deallocate( thisdist )

end subroutine find_closest_gridpoint


end program check_me

! <next few lines under version control, do not edit>
! $URL$
! $Id$
! $Revision$
! $Date$

