PROGRAM dart_tf_wrf

use wrf_data_module

implicit none

integer :: status, unit
logical :: dart_to_wrf
character (len=80) :: name
integer, parameter :: dart_unit = 10

type(wrf_data) :: wrf

real*8, pointer :: dart(:)
integer       :: number_dart_values
integer       :: seconds, days

!----
!  misc stuff

include 'netcdf.inc'
logical, parameter :: debug = .true.
integer :: mode

INTERFACE
subroutine dart_open_and_alloc( wrf, dart, n_values, dart_unit, &
                                dart_to_wrf, debug )

use wrf_data_module
implicit none

logical :: dart_to_wrf
integer :: dart_unit
logical :: debug
type(wrf_data) :: wrf
real*8, pointer :: dart(:)
integer :: n_values 
end subroutine dart_open_and_alloc

subroutine dart_io( in_or_out, dart, dart_unit, n_values, &
                    seconds, days, wrf, debug )

use wrf_data_module
implicit none

logical :: dart_to_wrf
integer :: dart_unit, seconds, days
character (len=6) :: in_or_out
logical :: debug
type(wrf_data) :: wrf
real*8, pointer :: dart(:)
integer :: n_values

end subroutine dart_io

subroutine transfer_dart_wrf ( dart_to_wrf, dart, wrf, n_values_in, debug )

use wrf_data_module
implicit none

logical :: dart_to_wrf
logical :: debug
type(wrf_data) :: wrf
real*8, pointer :: dart(:)
integer :: n_values_in 

end subroutine transfer_dart_wrf

END INTERFACE


!---
!  begin

!dart_to_wrf = .false.  ! this should be an input
!seconds = 60           ! this should be an input
!days = 365             ! this should be an input
read(5,*) dart_to_wrf
write(6,*) ' dart_to_wrf is : ',dart_to_wrf
!read(5,*) seconds
!write(6,*) ' input seconds  : ',seconds
!read(5,*) days
!write(6,*) ' input days     : ',days

wrf%ice_micro = .false.

!  open wrf data file

mode = 0
if( dart_to_wrf ) mode = NF_WRITE

if(debug) write(6,*) ' wrf_open_and_alloc '
  call wrf_open_and_alloc( wrf, mode, debug )
if(debug) write(6,*) ' returned from wrf_open_and_alloc '

!---
! allocate space for DART data

if(debug) write(6,*) ' dart_open_and_alloc '
  call dart_open_and_alloc( wrf, dart, number_dart_values, &
                            dart_unit, dart_to_wrf, debug  )
if(debug) write(6,*) ' returned from dart_open_and_alloc '

!---
!  get DART data or WRF data

if(debug) write(6,*) ' state input '
if( dart_to_wrf ) then
  call DART_IO( "INPUT ", dart, dart_unit, number_dart_values, &
                 seconds, days, wrf, debug )
  unit = dart_unit + 1
  open(unit = unit, file = 'time.dat')
  write (unit,*) seconds
  write (unit,*) days
  close(unit)
else
  call WRF_IO( wrf, "INPUT ", debug )
end if
if(debug) write(6,*) ' returned from state input '

!---
!  translate from DART to WRF, or WRF to DART

if(debug) write(6,*) ' transfer data to_from dart-wrf '

   call transfer_dart_wrf ( dart_to_wrf, dart, wrf,    &
                            number_dart_values, debug )

if(debug) write(6,*) ' transfer complete '

!---
!  output 

if(debug) write(6,*) ' state output '
if( dart_to_wrf ) then
  call WRF_IO( wrf, "OUTPUT", debug )
else
  unit = dart_unit + 1
  open(unit = unit, file = 'time.dat')
  read(unit,*) seconds
  read(unit,*) days
  close(unit)
  call DART_IO( "OUTPUT", dart, dart_unit, number_dart_values, &
       seconds, days, wrf, debug )
end if
if(debug) write(6,*) ' returned from state output '

status = nf_close(wrf%ncid)
if (status /= nf_noerr) call handle_err(4,status)

END 

!******************************************************************************

subroutine wrf_io( wrf, in_or_out, debug )
use wrf_data_module
implicit none
type(wrf_data) wrf
character (len=6) :: in_or_out
integer, dimension(5) :: map, count, start, stride
integer :: k
logical :: debug

map = 1
start = 1
stride = 1

count(1) = 1
count(2) = wrf%bt
count(3) = wrf%sn
count(4) = wrf%we+1
if(debug) write(6,*) ' calling netcdf read for u ', &
           count(1), count(2), count(3), count(4)
call netcdf_read_write_var( wrf%ncid, wrf%u_id, wrf%u,    &
                            start, count, stride, map, in_or_out )

if(debug) write(6,*) ' returned from netcdf read for u '
count(3) = wrf%sn+1
count(4) = wrf%we
if(debug) write(6,*) ' calling netcdf read for v '
call netcdf_read_write_var( wrf%ncid, wrf%v_id, wrf%v,    &
                            start, count, stride, map, in_or_out )
if(debug) write(6,*) ' returned from netcdf read for v '


count(2) = wrf%bt+1
count(3) = wrf%sn
count(4) = wrf%we
if(debug) write(6,*) ' calling netcdf read for w '
call netcdf_read_write_var( wrf%ncid, wrf%w_id, wrf%w,    &
                            start, count, stride, map, in_or_out )
if(debug) write(6,*) ' returned from netcdf read for w '
call netcdf_read_write_var( wrf%ncid, wrf%ph_id, wrf%ph,    &
                            start, count, stride, map, in_or_out )
if( in_or_out == "INPUT ")   &  ! get base state goepot. for full state computation
call netcdf_read_write_var( wrf%ncid, wrf%phb_id, wrf%phb,    &
                            start, count, stride, map, in_or_out )

count(2) = wrf%bt
count(3) = wrf%sn
count(4) = wrf%we
call netcdf_read_write_var( wrf%ncid, wrf%t_id, wrf%t,    &
                            start, count, stride, map, in_or_out )
call netcdf_read_write_var( wrf%ncid, wrf%qv_id, wrf%qv,    &
                            start, count, stride, map, in_or_out )
call netcdf_read_write_var( wrf%ncid, wrf%qc_id, wrf%qc,    &
                            start, count, stride, map, in_or_out )
call netcdf_read_write_var( wrf%ncid, wrf%qr_id, wrf%qr,    &
                            start, count, stride, map, in_or_out )



if(wrf%ice_micro) then
  call netcdf_read_write_var( wrf%ncid, wrf%qi_id, wrf%qi,    &
                              start, count, stride, map, in_or_out )
  call netcdf_read_write_var( wrf%ncid, wrf%qs_id, wrf%qs,    &
                              start, count, stride, map, in_or_out )
  call netcdf_read_write_var( wrf%ncid, wrf%qg_id, wrf%qg,    &
                              start, count, stride, map, in_or_out )
end if

count(2) = wrf%sn
count(3) = wrf%we
count(4) = 1
call netcdf_read_write_var( wrf%ncid, wrf%mu_id, wrf%mu,    &
                            start, count, stride, map, in_or_out )

if(debug) write(6,*) ' in_or_out is ',in_or_out
if( in_or_out(1:5) == "INPUT")   then  ! get base state mu. for full state computation

  if(debug) write(6,*) ' reading mub '
  call netcdf_read_write_var( wrf%ncid, wrf%mub_id, wrf%mub,    &
                              start, count, stride, map, in_or_out )
  if(debug) write(6,*) ' returned from mub read '
end if

if(debug) then
  do k=1,wrf%bt
    write(6,*) ' k, corner vals for u '
    write(6,*) k, wrf%u(1,1,k),wrf%u(wrf%we+1,1,k),  &
               wrf%u(1,wrf%sn,k),wrf%u(wrf%we+1,wrf%sn,k)
  enddo

  write(6,*) ' '

  do k=1,wrf%bt
    write(6,*) ' k, corner vals for v '
    write(6,*) k, wrf%v(1,1,k),wrf%v(wrf%we,1,k),  &
               wrf%v(1,wrf%sn+1,k),wrf%v(wrf%we,wrf%sn+1,k)
  enddo

  write(6,*) ' '

  write(6,*) ' corner vals for mu '
  write(6,*) wrf%mu(1,1),wrf%mu(wrf%we,1),  &
             wrf%mu(1,wrf%sn),wrf%mu(wrf%we,wrf%sn)
  write(6,*) ' corner vals for mub '
  write(6,*) wrf%mub(1,1),wrf%mub(wrf%we,1),  &
             wrf%mub(1,wrf%sn),wrf%mub(wrf%we,wrf%sn)
end if

end subroutine wrf_io

!*************************************************************************

subroutine netcdf_read_write_var( ncid, var_id, var,                    &
                                  start, count, stride, map, in_or_out )
implicit none
include 'netcdf.inc'
integer :: ncid, var_id
real, dimension(*) :: var
character (len=6) :: in_or_out
integer, dimension(5) :: start, count, stride, map
integer :: status


if( in_or_out(1:5) == "INPUT" ) then
  write(6,*) ' call netcdf read ', ncid, var_id
  status = nf_get_var_real(ncid, var_id, var, start, count, stride, map)
  write(6,*) ' returned netcdf read '
else if( in_or_out(1:6) == "OUTPUT" ) then
  status = nf_put_var_real(ncid, var_id, var, start, count, stride, map)
else
  write(6,*) ' unknown IO function for var_id ',var_id, in_or_out
  stop
end if
if (status /= nf_noerr) call handle_err(100,status)

end subroutine netcdf_read_write_var

!**********************************************************************

subroutine handle_err(ifn,status)
implicit none
integer :: ifn, status
character (len=80) :: nf_strerror

write(6,*) ' error for netcdf function ',ifn
write(6,*) ' status code = ',status
write(6,'(a80)') nf_strerror(status)

stop

end subroutine handle_err

!**********************************************************************

subroutine wrf_open_and_alloc( wrf, mode, debug )

use wrf_data_module
implicit none

include 'netcdf.inc'

character (len = 80)      :: path
     integer              :: mode
     integer              :: ncid

integer :: status
character (len=80) :: name
logical :: debug

type(wrf_data) :: wrf


status = nf_open('wrfinput', mode, wrf%ncid) 
if (status /= nf_noerr) call handle_err(1,status)
if(debug) write(6,*) ' wrf%ncid is ',wrf%ncid

! get wrf grid dimensions

status = nf_inq_dimid(wrf%ncid, "bottom_top", wrf%bt_id)
if (status /= nf_noerr) call handle_err(2,status)
status = nf_inq_dim(wrf%ncid, wrf%bt_id, name, wrf%bt)
if (status /= nf_noerr) call handle_err(3,status)

status = nf_inq_dimid(wrf%ncid, "south_north", wrf%sn_id)
if (status /= nf_noerr) call handle_err(2,status)
status = nf_inq_dim(wrf%ncid, wrf%sn_id, name, wrf%sn)
if (status /= nf_noerr) call handle_err(3,status)

status = nf_inq_dimid(wrf%ncid, "west_east", wrf%we_id)
if (status /= nf_noerr) call handle_err(2,status)
status = nf_inq_dim(wrf%ncid, wrf%we_id, name, wrf%we)
if (status /= nf_noerr) call handle_err(3,status)

write(6,*) ' dimensions bt, sn, we are ',wrf%bt,wrf%sn,wrf%we

!---
! get wrf variable ids and allocate space for wrf variables


status = nf_inq_varid(wrf%ncid, "P_TOP", wrf%ptop_id)
if(debug) write(6,*) ' ptop_id = ',wrf%ptop_id
if (status /= nf_noerr) call handle_err(4,status)

status = nf_inq_varid(wrf%ncid, "U", wrf%u_id)
if(debug) write(6,*) ' u_id = ',wrf%u_id
if (status /= nf_noerr) call handle_err(4,status)
allocate(wrf%u(wrf%we+1,wrf%sn,wrf%bt))

status = nf_inq_varid(wrf%ncid, "V", wrf%v_id)
if(debug) write(6,*) ' v_id = ',wrf%v_id
if (status /= nf_noerr) call handle_err(4,status)
allocate(wrf%v(wrf%we,wrf%sn+1,wrf%bt))

status = nf_inq_varid(wrf%ncid, "W", wrf%w_id)
if(debug) write(6,*) ' w_id = ',wrf%w_id
if (status /= nf_noerr) call handle_err(4,status)
allocate(wrf%w(wrf%we,wrf%sn,wrf%bt+1))

status = nf_inq_varid(wrf%ncid, "PH", wrf%ph_id)
if(debug) write(6,*) ' ph_id = ',wrf%ph_id
if (status /= nf_noerr) call handle_err(4,status)
allocate(wrf%ph(wrf%we,wrf%sn,wrf%bt+1))

status = nf_inq_varid(wrf%ncid, "PHB", wrf%phb_id)
if(debug) write(6,*) ' phb_id = ',wrf%phb_id
if (status /= nf_noerr) call handle_err(4,status)
allocate(wrf%phb(wrf%we,wrf%sn,wrf%bt+1))

status = nf_inq_varid(wrf%ncid, "T", wrf%t_id)
if(debug) write(6,*) ' t_id = ',wrf%t_id
if (status /= nf_noerr) call handle_err(4,status)
allocate(wrf%t(wrf%we,wrf%sn,wrf%bt))

status = nf_inq_varid(wrf%ncid, "MU", wrf%mu_id)
if(debug) write(6,*) ' mu_id = ',wrf%mu_id
if (status /= nf_noerr) call handle_err(4,status)
allocate(wrf%mu(wrf%we,wrf%sn))

status = nf_inq_varid(wrf%ncid, "MUB", wrf%mub_id)
if(debug) write(6,*) ' mub_id = ',wrf%mub_id
if (status /= nf_noerr) call handle_err(4,status)
allocate(wrf%mub(wrf%we,wrf%sn))

status = nf_inq_varid(wrf%ncid, "QVAPOR", wrf%qv_id)
if(debug) write(6,*) ' qv_id = ',wrf%qv_id
if (status /= nf_noerr) call handle_err(4,status)
allocate(wrf%qv(wrf%we,wrf%sn,wrf%bt))

status = nf_inq_varid(wrf%ncid, "QCLOUD", wrf%qc_id)
if(debug) write(6,*) ' qc_id = ',wrf%qc_id
if (status /= nf_noerr) call handle_err(4,status)
allocate(wrf%qc(wrf%we,wrf%sn,wrf%bt))

status = nf_inq_varid(wrf%ncid, "QRAIN", wrf%qr_id)
if(debug) write(6,*) ' qr_id = ',wrf%qr_id
if (status /= nf_noerr) call handle_err(4,status)
allocate(wrf%qr(wrf%we,wrf%sn,wrf%bt))

if (wrf%ice_micro) then
  if (status /= nf_noerr) call handle_err(4,status)
  status = nf_inq_varid(wrf%ncid, "QICE", wrf%qi_id)
  if(debug) write(6,*) ' qi_id = ',wrf%qi_id
  allocate(wrf%qi(wrf%we,wrf%sn,wrf%bt))

  if (status /= nf_noerr) call handle_err(4,status)
  status = nf_inq_varid(wrf%ncid, "QSNOW", wrf%qs_id)
  if(debug) write(6,*) ' qs_id = ',wrf%qs_id
  allocate(wrf%qs(wrf%we,wrf%sn,wrf%bt))

  if (status /= nf_noerr) call handle_err(4,status)
  status = nf_inq_varid(wrf%ncid, "QGRAUPEL", wrf%qg_id)
  if(debug) write(6,*) ' qg_id = ',wrf%qg_id
  allocate(wrf%qg(wrf%we,wrf%sn,wrf%bt))

end if

end subroutine wrf_open_and_alloc

!*****************************************************************************

subroutine dart_open_and_alloc( wrf, dart, n_values, dart_unit, &
                                dart_to_wrf, debug )

use wrf_data_module
implicit none

character (len = 80)      :: path
logical :: dart_to_wrf
integer :: dart_unit

integer :: status
character (len=80) :: name
logical :: debug

type(wrf_data) :: wrf
real*8, pointer :: dart(:)

integer :: n_values 

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

n_values = n_values + (wrf%bt  )*(wrf%sn  )*(wrf%we  )  ! qv field
n_values = n_values + (wrf%bt  )*(wrf%sn  )*(wrf%we  )  ! qc field
n_values = n_values + (wrf%bt  )*(wrf%sn  )*(wrf%we  )  ! qr field

if(wrf%ice_micro) then
  n_values = n_values + (wrf%bt  )*(wrf%sn  )*(wrf%we  )  ! qi field
  n_values = n_values + (wrf%bt  )*(wrf%sn  )*(wrf%we  )  ! qs field
  n_values = n_values + (wrf%bt  )*(wrf%sn  )*(wrf%we  )  ! qg field
end if

if(debug) write(6,*) ' dart vector length is ',n_values

allocate(dart(n_values))

!  open DART data file

if(dart_to_wrf)  then  !  DART data file should exist, open it
  open( unit=dart_unit,file="dart_wrf_vector",form="unformatted",  &
        status="old",action="read" )
else
  open( unit=dart_unit,file="dart_wrf_vector",form="unformatted",  &
        status="new",action="write" )
end if

end subroutine dart_open_and_alloc

!*****************************************************************************

subroutine dart_io( in_or_out, dart, dart_unit, n_values, &
                     seconds, days, wrf, debug )

use wrf_data_module
implicit none

character (len = 80)      :: path
logical :: dart_to_wrf
integer :: dart_unit, seconds, days
integer :: seconds_end, days_end

integer :: status
character (len=6) :: in_or_out
logical :: debug

logical, parameter :: test_input = .false.

type(wrf_data) :: wrf
real*8, pointer :: dart(:)

integer :: n_values 

if(debug) then
     write(6,*)' in dart_io '
     write(6,*)' seconds, days = ', seconds, days
     write(6,*)' number of values ', n_values
     write(6,*)' mode ',in_or_out
end if

if (in_or_out(1:5) == 'INPUT') then
  read(dart_unit) seconds, days 
  if(debug) write(6,*) ' seconds and days from dart input ',seconds, days 
  if(.not. test_input) then
     read(dart_unit) seconds_end, days_end
     if(debug) write(6,*) ' seconds_end and days_end from dart input ', &
                             seconds_end, days_end 
  end if
  read(dart_unit) dart
else
  rewind(dart_unit)
  write(dart_unit) seconds, days 
  write(dart_unit) dart
end if

end subroutine dart_io

!*************************************************************************

subroutine transfer_dart_wrf ( dart_to_wrf, dart, wrf, n_values_in, debug )

use wrf_data_module
implicit none

logical :: dart_to_wrf

integer :: status
logical :: debug

type(wrf_data) :: wrf
real*8, pointer :: dart(:)

integer :: n_values_in 

!---

integer :: in, n_values

!---

INTERFACE

subroutine trans_2d( one_to_two, a1d,a2d, nx, ny )

implicit none
integer :: nx,ny
real*8 :: a1d(:)
real :: a2d(nx,ny)
logical ::  one_to_two

end subroutine trans_2d

subroutine trans_3d( one_to_three, a1d,a3d, nx, ny,nz )

implicit none
integer :: nx,ny,nz
real*8 :: a1d(:)
real :: a3d(nx,ny,nz)
logical ::  one_to_three

end subroutine trans_3d

END INTERFACE


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

in = n_values+1
call trans_3d( dart_to_wrf, dart(in:),wrf%qv,wrf%we,wrf%sn,wrf%bt)
n_values = n_values + (wrf%bt  )*(wrf%sn  )*(wrf%we  )  ! qv field

in = n_values+1
call trans_3d( dart_to_wrf, dart(in:),wrf%qc,wrf%we,wrf%sn,wrf%bt)
n_values = n_values + (wrf%bt  )*(wrf%sn  )*(wrf%we  )  ! qc field

in = n_values+1
call trans_3d( dart_to_wrf, dart(in:),wrf%qr,wrf%we,wrf%sn,wrf%bt)
n_values = n_values + (wrf%bt  )*(wrf%sn  )*(wrf%we  )  ! qr field


if(wrf%ice_micro) then

  in = n_values+1
  call trans_3d( dart_to_wrf, dart(in:),wrf%qi,wrf%we,wrf%sn,wrf%bt)
  n_values = n_values + (wrf%bt  )*(wrf%sn  )*(wrf%we  )  ! qi field

  in = n_values+1
  call trans_3d( dart_to_wrf, dart(in:),wrf%qs,wrf%we,wrf%sn,wrf%bt)
  n_values = n_values + (wrf%bt  )*(wrf%sn  )*(wrf%we  )  ! qs field

  in = n_values+1
  call trans_3d( dart_to_wrf, dart(in:),wrf%qg,wrf%we,wrf%sn,wrf%bt)
  n_values = n_values + (wrf%bt  )*(wrf%sn  )*(wrf%we  )  ! qg field

end if

if(n_values /= n_values_in ) then
  write(6,*) ' n_values differ in transfer ',n_values, n_values_in
  write(6,*) ' error exit '
  stop
end if

end subroutine transfer_dart_wrf

!*********************************************************************

subroutine trans_2d( one_to_two, a1d, a2d, nx, ny )

implicit none
integer :: nx,ny
real*8 :: a1d(:)
real :: a2d(nx,ny)
logical ::  one_to_two

!---

integer i,j

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
integer :: nx,ny,nz
real*8 :: a1d(:)
real :: a3d(nx,ny,nz)
logical ::  one_to_three

!---

integer i,j,k

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

