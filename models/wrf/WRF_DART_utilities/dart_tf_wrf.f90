! Data Assimilation Research Testbed -- DART
! Copyright 2004, Data Assimilation Initiative, University Corporation for Atmospheric Research
! Licensed under the GPL -- www.gpl.org/licenses/gpl.html
 
PROGRAM dart_tf_wrf

! <next three lines automatically updated by CVS, do not edit>
! $Source$
! $Revision$
! $Date$

use        types_mod, only : r8
use time_manager_mod, only : time_type, write_time, read_time, get_date, set_date, operator(-), &
                             get_time, print_time, set_calendar_type, GREGORIAN, days_per_month
use    utilities_mod, only : get_unit, file_exist, open_file, check_nml_error, close_file, &
                             error_handler, E_ERR, E_MSG, initialize_utilities, &
                             finalize_utilities, register_module, logfileunit
use  wrf_data_module, only : wrf_data, wrf_open_and_alloc, wrf_dealloc, wrf_io, set_wrf_date, &
                             get_wrf_date
use                          netcdf

implicit none

! CVS Generated file description for error handling, do not edit
character(len=128) :: &
source   = "$Source$", &
revision = "$Revision$", &
revdate  = "$Date$"

!-----------------------------------------------------------------------
! Model namelist parameters with default values
!-----------------------------------------------------------------------

logical :: output_state_vector  = .true.  ! output prognostic variables
integer :: num_moist_vars       = 0
integer :: wrf_dt = 300, wrf_dx = 100000
integer :: calendar_type        = GREGORIAN

namelist /model_nml/ output_state_vector, num_moist_vars, wrf_dt, wrf_dx

!-------------------------------------------------------------
! Namelist with default values
! binary_restart_files  == .true.  -> use unformatted file format. 
!                                     Full precision, faster, smaller,
!                                     but not as portable.
! binary_restart_files  == .false.  -> use ascii file format. 
!                                     Portable, but loses precision,
!                                     slower, and larger.

logical  :: binary_restart_files = .false.

namelist /assim_model_nml/ binary_restart_files
!-------------------------------------------------------------

integer :: iunit, dart_unit
logical :: dart_to_wrf, leap

type(wrf_data) :: wrf

real(r8), pointer :: dart(:)
type(time_type)   :: dart_time(2), interval_time
integer           :: number_dart_values, days, seconds, &
                     year, month, day, hour, minute, second
integer           :: ndays, m
integer           :: ndims, idims(2), dimids(2)
integer           :: i, ivtype
character(len=80) :: varname

character(len=19) :: timestring

!----
!  misc stuff

logical, parameter :: debug = .false.
integer            :: mode, io, ierr, var_id

call initialize_utilities
call register_module(source, revision, revdate)
write(logfileunit,*)'STARTING dart_tf_wrf ...'

! Begin by reading the namelist input
if(file_exist('input.nml')) then

   iunit = open_file('input.nml', action = 'read')
   read(iunit, nml = assim_model_nml, iostat = io )
   ierr = check_nml_error(io, 'assim_model_nml')

   rewind(iunit)
   read(iunit, nml = model_nml, iostat = io )
   ierr = check_nml_error(io, 'model_nml')

   call close_file(iunit)

   if ( debug ) then
      ! namelist validation
      write(*, *) 'assim_model_nml read; values are'
      write(*, *) 'binary_restart_files is ', binary_restart_files
      write(*,'(''num_moist_vars = '',i3)')num_moist_vars
   endif
endif

call set_calendar_type(calendar_type)

read(5,*) dart_to_wrf

if ( dart_to_wrf ) then
   call error_handler(E_MSG,'dart_to_wrf', &
       'Converting a dart state vector to a WRF netcdf restart file', &
       source, revision, revdate)
else
   call error_handler(E_MSG,'dart_to_wrf', &
       'Converting a WRF netcdf restart file to a dart state vector', &
       source, revision, revdate)
endif

wrf%n_moist = num_moist_vars

! open wrf data netcdf netCDF file 'wrfinput'
! we get sizes of the WRF geometry and resolution

mode = NF90_NOWRITE                   ! read the netcdf file
if( dart_to_wrf ) mode = NF90_WRITE   ! write to the netcdf file

if(debug) write(6,*) ' wrf_open_and_alloc '
call wrf_open_and_alloc( wrf, 'wrfinput', mode, debug )
if(debug) write(6,*) ' returned from wrf_open_and_alloc '

!---
! allocate space for DART data

if(debug) write(6,*) ' dart_open_and_alloc '
dart_unit = get_unit()
call dart_open_and_alloc( wrf, dart, number_dart_values, &
     dart_unit, dart_to_wrf, binary_restart_files, debug  )
if(debug) write(6,*) ' returned from dart_open_and_alloc '

!----------------------------------------------------------------------
!  get DART data or WRF data

if(debug) write(6,*) ' state input '

if( dart_to_wrf ) then

   call DART_IO( "INPUT ", dart, dart_unit, dart_time, binary_restart_files, debug )
   iunit = get_unit()
   open(unit = iunit, file = 'wrf.info')
   call write_time(iunit, dart_time(1))
   call get_date(dart_time(2), year, month, day, hour, minute, second)
   write (iunit,*) year, month, day, hour, minute, second
   call get_date(dart_time(1), year, month, day, hour, minute, second)
   write (iunit,*) year, month, day, hour, minute, second

   interval_time = dart_time(2) - dart_time(1)
   call get_time(interval_time, seconds, days)

   write (iunit,*) (days *24 *3600) + seconds

   write (iunit,*) wrf_dt
   write (iunit,*) wrf_dx
   write (iunit,*) wrf%we+1
   write (iunit,*) wrf%sn+1
   write (iunit,*) wrf%bt+1
   close(iunit)

else

   call WRF_IO( wrf, "INPUT ", debug )

end if
if(debug) write(6,*) ' returned from state input '

!---
!  translate from DART to WRF, or WRF to DART

if(debug) write(6,*) ' transfer data to_from dart-wrf '

call transfer_dart_wrf ( dart_to_wrf, dart, wrf,    &
     number_dart_values )

if(debug) write(6,*) ' transfer complete '

!---
!  output

if(debug) write(6,*) ' state output '
if( dart_to_wrf ) then
   call get_date(dart_time(2), year, month, day, hour, minute, second)
   call set_wrf_date(timestring, year, month, day, hour, minute, second)
   call check( nf90_inq_varid(wrf%ncid, "Times", var_id) )
   call check( nf90_put_var(wrf%ncid, var_id, timestring) )
   call check( nf90_put_att(wrf%ncid, nf90_global, "START_DATE", timestring) )
   call check( nf90_put_att(wrf%ncid, nf90_global, "JULYR", year) )
   ndays = 0
   leap = (modulo(year,4) == 0)
   if((modulo(year,100).eq.0).and.(modulo(year,400).ne.0))then
      leap=.false.
   endif
   do m = 1, month - 1
      ndays = ndays + days_per_month(m)
      if(leap .and. m == 2) ndays = ndays + 1
   enddo
   ndays = ndays + day
   call check( nf90_put_att(wrf%ncid, nf90_global, "JULDAY", ndays) )
   call WRF_IO( wrf, "OUTPUT", debug )
else
   iunit = get_unit()
   call check( nf90_inq_varid(wrf%ncid, "Times", var_id) )
   call check( nf90_Inquire_Variable(wrf%ncid, var_id, varname, xtype=ivtype, ndims=ndims, dimids=dimids) )
   do i=1,ndims
      call check( nf90_inquire_dimension(wrf%ncid, dimids(i), len=idims(i)) )
      if(debug) write(6,*) ' dimension ',i,idims(i)
   enddo

   call check( nf90_get_var(wrf%ncid, var_id, timestring, start = (/ 1, idims(2) /)) )
   call get_wrf_date(timestring, year, month, day, hour, minute, second)
   dart_time(1) = set_date(year, month, day, hour, minute, second)

   write(6,*) 'Date from wrfinput'
   write(6,*) year, month, day, hour, minute, second

   write(6,*) 'Time written to dart vector file:'
   call print_time(dart_time(1))

   call DART_IO( "OUTPUT", dart, dart_unit, dart_time, binary_restart_files, debug )
end if
if(debug) write(6,*) ' returned from state output '

call check ( nf90_sync(wrf%ncid) )
call check ( nf90_close(wrf%ncid) )

call wrf_dealloc(wrf)
deallocate(dart)

write(logfileunit,*)'FINISHED dart_tf_wrf.'
write(logfileunit,*)

call finalize_utilities ! closes the log file.
 
contains

  ! Internal subroutine - checks error status after each netcdf, prints 
  !                       text message each time an error code is returned. 
  subroutine check(istatus)
    integer, intent ( in) :: istatus 
    if(istatus /= nf90_noerr) call error_handler(E_ERR,'dart_tf_wrf', &
         trim(nf90_strerror(istatus)), source, revision, revdate) 
  end subroutine check

!*****************************************************************************

subroutine dart_open_and_alloc( wrf, dart, n_values, dart_unit, &
     dart_to_wrf, binary_restart_files, debug )

implicit none

integer,           intent(in)  :: dart_unit

logical,           intent(in)  :: dart_to_wrf, binary_restart_files, debug

type(wrf_data),    intent(in)  :: wrf
real(r8),          pointer     :: dart(:)

integer,           intent(out) :: n_values 

! compute number of values in 1D vector

! dry dynamics conponents

n_values = 0

n_values = n_values + (wrf%bt  )*(wrf%sn  )*(wrf%we+1)  ! u field
n_values = n_values + (wrf%bt  )*(wrf%sn+1)*(wrf%we  )  ! v field
n_values = n_values + (wrf%bt+1)*(wrf%sn  )*(wrf%we  )  ! w field
n_values = n_values + (wrf%bt+1)*(wrf%sn  )*(wrf%we  )  ! geopotential field
n_values = n_values + (wrf%bt  )*(wrf%sn  )*(wrf%we  )  ! t field
n_values = n_values + (       1)*(wrf%sn  )*(wrf%we  )  ! dry surface pressure field

! moist variables

if(wrf%n_moist > 0) then
   n_values = n_values + (wrf%bt  )*(wrf%sn  )*(wrf%we  )  ! qv field
endif
if(wrf%n_moist > 1) then
   n_values = n_values + (wrf%bt  )*(wrf%sn  )*(wrf%we  )  ! qc field
endif
if(wrf%n_moist > 2) then
   n_values = n_values + (wrf%bt  )*(wrf%sn  )*(wrf%we  )  ! qr field
endif
if(wrf%n_moist > 3) then
   n_values = n_values + (wrf%bt  )*(wrf%sn  )*(wrf%we  )  ! qi field
endif
if(wrf%n_moist > 4) then
   n_values = n_values + (wrf%bt  )*(wrf%sn  )*(wrf%we  )  ! qs field
endif
if(wrf%n_moist > 5) then
   n_values = n_values + (wrf%bt  )*(wrf%sn  )*(wrf%we  )  ! qg field
endif
if(wrf%n_moist > 6) then
   write(6,*) 'n_moist = ',wrf%n_moist,' is too large.'
   stop
endif

if(debug) write(6,*) ' dart vector length is ',n_values

allocate(dart(n_values))

!  open DART data file

if(dart_to_wrf)  then  !  DART data file should exist, open it
   if ( binary_restart_files ) then
      open( unit=dart_unit,file="dart_wrf_vector",form="unformatted",  &
           status="old",action="read" )
   else
      open( unit=dart_unit,file="dart_wrf_vector",form="formatted",  &
           status="old",action="read" )
   endif
else
   if ( binary_restart_files ) then
      open( unit=dart_unit,file="dart_wrf_vector",form="unformatted",  &
           status="new",action="write" )
   else
      open( unit=dart_unit,file="dart_wrf_vector",form="formatted",  &
           status="new",action="write" )
   endif
end if

end subroutine dart_open_and_alloc

!*****************************************************************************

subroutine dart_io( in_or_out, dart, dart_unit, dart_time, binary_restart_files, debug )

implicit none

character (len=6), intent(in)    :: in_or_out
real(r8),          pointer       :: dart(:)

integer,           intent(in)    :: dart_unit
type(time_type),   intent(inout) :: dart_time(2)
logical,           intent(in)    :: binary_restart_files, debug

if(debug) then
   write(6,*)' in dart_io '
   write(6,*)' mode ',in_or_out
end if

if (in_or_out(1:5) == 'INPUT') then
   if ( binary_restart_files ) then
      dart_time(1) = read_time(dart_unit, "unformatted")
      dart_time(2) = read_time(dart_unit, "unformatted")
   else
      dart_time(1) = read_time(dart_unit)
      dart_time(2) = read_time(dart_unit)
   endif
   if(debug) call print_time(dart_time(1), 'Target time:')
   if(debug) call print_time(dart_time(2), 'Current time:')
   read(dart_unit) dart
else
   rewind(dart_unit)
   if ( binary_restart_files ) then
      call write_time(dart_unit, dart_time(1), "unformatted")
   else
      call write_time(dart_unit, dart_time(1))
   endif
   write(dart_unit) dart
end if

end subroutine dart_io

!*************************************************************************

subroutine transfer_dart_wrf ( dart_to_wrf, dart, wrf, n_values_in)

implicit none

logical :: dart_to_wrf

type(wrf_data)    :: wrf
real(r8), pointer :: dart(:)

integer :: n_values_in

!---

integer           :: in, n_values
character(len=80) :: stringerror

!---


n_values = 0
in = n_values+1
call trans_3d( dart_to_wrf, dart(in:),wrf%u,wrf%we+1,wrf%sn,wrf%bt)
n_values = n_values + (wrf%bt  )*(wrf%sn  )*(wrf%we+1)  ! u field

in = n_values+1
call trans_3d( dart_to_wrf, dart(in:),wrf%v,wrf%we,wrf%sn+1,wrf%bt)
n_values = n_values + (wrf%bt  )*(wrf%sn+1)*(wrf%we  )  ! v field

in = n_values+1
call trans_3d( dart_to_wrf, dart(in:),wrf%w,wrf%we,wrf%sn,wrf%bt+1)
n_values = n_values + (wrf%bt+1)*(wrf%sn  )*(wrf%we  )  ! w field

in = n_values+1
call trans_3d( dart_to_wrf, dart(in:),wrf%ph,wrf%we,wrf%sn,wrf%bt+1)
n_values = n_values + (wrf%bt+1)*(wrf%sn  )*(wrf%we  )  ! geopotential field

in = n_values+1
call trans_3d( dart_to_wrf, dart(in:),wrf%t,wrf%we,wrf%sn,wrf%bt)
n_values = n_values + (wrf%bt  )*(wrf%sn  )*(wrf%we  )  ! t field

in = n_values+1
call trans_2d( dart_to_wrf, dart(in:),wrf%mu,wrf%we,wrf%sn)
n_values = n_values + (       1)*(wrf%sn  )*(wrf%we  )  ! dry surface pressure field

! moist variables

if(wrf%n_moist > 0) then
   in = n_values+1
   call trans_3d( dart_to_wrf, dart(in:),wrf%qv,wrf%we,wrf%sn,wrf%bt)
   n_values = n_values + (wrf%bt  )*(wrf%sn  )*(wrf%we  )  ! qv field
endif
if(wrf%n_moist > 1) then
   in = n_values+1
   call trans_3d( dart_to_wrf, dart(in:),wrf%qc,wrf%we,wrf%sn,wrf%bt)
   n_values = n_values + (wrf%bt  )*(wrf%sn  )*(wrf%we  )  ! qc field
endif
if(wrf%n_moist > 2) then
   in = n_values+1
   call trans_3d( dart_to_wrf, dart(in:),wrf%qr,wrf%we,wrf%sn,wrf%bt)
   n_values = n_values + (wrf%bt  )*(wrf%sn  )*(wrf%we  )  ! qr field
endif
if(wrf%n_moist > 3) then
   in = n_values+1
   call trans_3d( dart_to_wrf, dart(in:),wrf%qi,wrf%we,wrf%sn,wrf%bt)
   n_values = n_values + (wrf%bt  )*(wrf%sn  )*(wrf%we  )  ! qi field
endif
if(wrf%n_moist > 4) then
   in = n_values+1
   call trans_3d( dart_to_wrf, dart(in:),wrf%qs,wrf%we,wrf%sn,wrf%bt)
   n_values = n_values + (wrf%bt  )*(wrf%sn  )*(wrf%we  )  ! qs field
endif
if(wrf%n_moist > 5) then
   in = n_values+1
   call trans_3d( dart_to_wrf, dart(in:),wrf%qg,wrf%we,wrf%sn,wrf%bt)
   n_values = n_values + (wrf%bt  )*(wrf%sn  )*(wrf%we  )  ! qg field
endif
if(wrf%n_moist > 6) then
   write(stringerror,*) 'n_moist = ',wrf%n_moist,' is too large.'
   call error_handler(E_ERR, 'transfer_dart_wrf', &
                             stringerror, source, revision, revdate)
endif

if(n_values /= n_values_in ) then
   write(stringerror,*)' n_values differ in transfer ',n_values, n_values_in
   call error_handler(E_ERR, 'transfer_dart_wrf', &
                             stringerror, source, revision, revdate)
endif

end subroutine transfer_dart_wrf

!*********************************************************************

subroutine trans_2d( one_to_two, a1d, a2d, nx, ny )

implicit none

integer  :: nx,ny
real(r8) :: a1d(:)
real(r8) :: a2d(nx,ny)
logical  :: one_to_two

!---

integer :: i,j,m
character(len=129) :: errstring

i=size(a2d,1)
j=size(a2d,2)
m=size(a1d)

if ( i /= nx .or. &
     j /= ny .or. &
     m < nx*ny) then
   write(errstring,*)'nx, ny, not compatible ',i,j,nx,ny
   call error_handler(E_ERR,'trans_2d',errstring,source,revision,revdate)
endif


if (one_to_two) then

   do j=1,ny
      do i=1,nx
         a2d(i,j) = a1d(i + nx*(j-1))
      enddo
   enddo

else

   do j=1,ny
      do i=1,nx
         a1d(i + nx*(j-1)) = a2d(i,j)
      enddo
   enddo

end if

end subroutine trans_2d

!*********************************************************************

subroutine trans_3d( one_to_three, a1d, a3d, nx, ny, nz )

implicit none

integer  :: nx,ny,nz
real(r8) :: a1d(:)
real(r8) :: a3d(:,:,:)
logical  :: one_to_three

!---

integer :: i,j,k,m
character(len=129) :: errstring

i=size(a3d,1)
j=size(a3d,2)
k=size(a3d,3)
m=size(a1d)

if ( i /= nx .or. &
     j /= ny .or. &
     k /= nz .or. &
     m < nx*ny*nz) then
   write(errstring,*)'nx, ny, nz, not compatible ',i,j,k,nx,ny,nz
   call error_handler(E_ERR,'trans_3d',errstring,source,revision,revdate)
endif

if (one_to_three) then

   do k=1,nz
      do j=1,ny
         do i=1,nx
            a3d(i,j,k) = a1d(i + nx*(j-1) + nx*ny*(k-1) )
         enddo
      enddo
   enddo

else

   do k=1,nz
      do j=1,ny
         do i=1,nx
            a1d(i + nx*(j-1) + nx*ny*(k-1) ) = a3d(i,j,k)
         enddo
      enddo
   enddo

end if

end subroutine trans_3d

END PROGRAM dart_tf_wrf
