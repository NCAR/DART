! Data Assimilation Research Testbed -- DART
! Copyright 2004, Data Assimilation Initiative, University Corporation for Atmospheric Research
! Licensed under the GPL -- www.gpl.org/licenses/gpl.html
 
PROGRAM dart_tf_wrf

! <next three lines automatically updated by CVS, do not edit>
! $Source$
! $Revision$
! $Date$

use        types_mod, only : r8
use    utilities_mod, only : get_unit, file_exist, open_file, check_nml_error, close_file, &
                             error_handler, E_ERR
use netcdf
use        model_mod, only : netcdf_read_write_var
use wrf_data_module

implicit none

! CVS Generated file description for error handling, do not edit
character(len=128) :: &
source   = "$Source$", &
revision = "$Revision$", &
revdate  = "$Date$"

!-----------------------------------------------------------------------
! model namelist parameters
!-----------------------------------------------------------------------

logical :: output_state_vector = .true.  ! output prognostic variables
integer :: num_moist_vars = 0            ! default value

namelist /model_nml/ output_state_vector, num_moist_vars

integer :: status, iunit, dart_unit
logical :: dart_to_wrf

type(wrf_data) :: wrf

real(r8), pointer :: dart(:)
integer       :: number_dart_values
integer       :: seconds, days

!----
!  misc stuff

include 'netcdf.inc'
logical, parameter :: debug = .false.
integer :: mode, io, ierr

!---
!  begin

! Begin by reading the namelist input                                           
if(file_exist('input.nml')) then

   iunit = open_file('input.nml', action = 'read')
   read(iunit, nml = model_nml, iostat = io )
   ierr = check_nml_error(io, 'model_nml')
   call close_file(iunit)                                                        

   if ( debug ) then
      write(*,'(''num_moist_vars = '',i3)')num_moist_vars
   endif
endif

read(5,*) dart_to_wrf
write(6,*) ' dart_to_wrf is : ',dart_to_wrf

if( (num_moist_vars >= 4).and.(num_moist_vars <= 6) ) then
   wrf%ice_micro = .true.
else
   wrf%ice_micro = .false.
endif

!  open wrf data file

mode = 0
if( dart_to_wrf ) mode = NF_WRITE

if(debug) write(6,*) ' wrf_open_and_alloc '
call wrf_open_and_alloc( wrf, mode, debug )
if(debug) write(6,*) ' returned from wrf_open_and_alloc '

!---
! allocate space for DART data

if(debug) write(6,*) ' dart_open_and_alloc '
dart_unit = get_unit()
call dart_open_and_alloc( wrf, dart, number_dart_values, &
     dart_unit, dart_to_wrf, debug  )
if(debug) write(6,*) ' returned from dart_open_and_alloc '

!---
!  get DART data or WRF data

if(debug) write(6,*) ' state input '
if( dart_to_wrf ) then
   call DART_IO( "INPUT ", dart, dart_unit, number_dart_values, &
        seconds, days, debug )
   iunit = get_unit()
   open(unit = iunit, file = 'time.dat')
   write (iunit,*) seconds
   write (iunit,*) days
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
   call WRF_IO( wrf, "OUTPUT", debug )
else
   iunit = get_unit()
   open(unit = iunit, file = 'time.dat')
   read(iunit,*) seconds
   read(iunit,*) days
   close(iunit)
   call DART_IO( "OUTPUT", dart, dart_unit, number_dart_values, &
        seconds, days, debug )
end if
if(debug) write(6,*) ' returned from state output '

status = nf_close(wrf%ncid)
if (status /= nf_noerr) call error_handler(E_ERR,'main', &
     trim(nf_strerror(status)), source, revision, revdate)

contains

!******************************************************************************

subroutine wrf_io( wrf, in_or_out, debug )

implicit none

type(wrf_data) wrf
character (len=6) :: in_or_out
integer, dimension(5) :: map, count, start, stride
integer :: k
logical :: debug

map = 1
start = 1
stride = 1

count(4) = 1
count(3) = wrf%bt
count(2) = wrf%sn
count(1) = wrf%we+1
if(debug) write(6,*) ' calling netcdf read for u ', &
     count(1), count(2), count(3), count(4)
call netcdf_read_write_var( "U", wrf%ncid, wrf%u_id, wrf%u,    &
     start, count, stride, map, in_or_out, debug, 4 )
if(debug) write(6,*) ' returned from netcdf read for u '

count(2) = wrf%sn+1
count(1) = wrf%we
if(debug) write(6,*) ' calling netcdf read for v '
call netcdf_read_write_var( "V", wrf%ncid, wrf%v_id, wrf%v,    &
     start, count, stride, map, in_or_out, debug, 4 )
if(debug) write(6,*) ' returned from netcdf read for v '


count(3) = wrf%bt+1
count(2) = wrf%sn
count(1) = wrf%we
if(debug) write(6,*) ' calling netcdf read for w '
call netcdf_read_write_var( "W", wrf%ncid, wrf%w_id, wrf%w,    &
     start, count, stride, map, in_or_out, debug, 4 )
if(debug) write(6,*) ' returned from netcdf read for w '
call netcdf_read_write_var( "PH", wrf%ncid, wrf%ph_id, wrf%ph,    &
     start, count, stride, map, in_or_out, debug, 4 )
if( in_or_out == "INPUT ")   &  ! get base state goepot. for full state computation
     call netcdf_read_write_var( "PHB", wrf%ncid, wrf%phb_id, wrf%phb,    &
     start, count, stride, map, in_or_out, debug, 4 )

count(3) = wrf%bt
count(2) = wrf%sn
count(1) = wrf%we
call netcdf_read_write_var( "T", wrf%ncid, wrf%t_id, wrf%t,    &
     start, count, stride, map, in_or_out, debug, 4 )
call netcdf_read_write_var( "QVAPOR", wrf%ncid, wrf%qv_id, wrf%qv,    &
     start, count, stride, map, in_or_out, debug, 4 )
call netcdf_read_write_var( "QCLOUD", wrf%ncid, wrf%qc_id, wrf%qc,    &
     start, count, stride, map, in_or_out, debug, 4 )
call netcdf_read_write_var( "QRAIN", wrf%ncid, wrf%qr_id, wrf%qr,    &
     start, count, stride, map, in_or_out, debug, 4 )



if(wrf%ice_micro) then
   call netcdf_read_write_var( "QICE", wrf%ncid, wrf%qi_id, wrf%qi,    &
        start, count, stride, map, in_or_out, debug, 4 )
   call netcdf_read_write_var( "QSNOW", wrf%ncid, wrf%qs_id, wrf%qs,    &
        start, count, stride, map, in_or_out, debug, 4 )
   call netcdf_read_write_var( "QGRAP", wrf%ncid, wrf%qg_id, wrf%qg,    &
        start, count, stride, map, in_or_out, debug, 4 )
end if

count(3) = 1
count(2) = wrf%sn
count(1) = wrf%we
call netcdf_read_write_var( "MU", wrf%ncid, wrf%mu_id, wrf%mu,    &
     start, count, stride, map, in_or_out, debug, 3 )

if(debug) write(6,*) ' in_or_out is ',in_or_out
if( in_or_out(1:5) == "INPUT")   then  ! get base state mu. for full state computation

   if(debug) write(6,*) ' reading mub '
   call netcdf_read_write_var( "MUB", wrf%ncid, wrf%mub_id, wrf%mub,    &
        start, count, stride, map, in_or_out, debug, 3 )
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

!**********************************************************************

subroutine wrf_open_and_alloc( wrf, mode, debug )

implicit none

include 'netcdf.inc'

integer            :: mode

character (len=80) :: name
logical            :: debug

type(wrf_data)     :: wrf


call check ( nf_open('wrfinput', mode, wrf%ncid) )
if(debug) write(6,*) ' wrf%ncid is ',wrf%ncid

! get wrf grid dimensions

call check ( nf_inq_dimid(wrf%ncid, "bottom_top", wrf%bt_id))
call check ( nf_inq_dim(wrf%ncid, wrf%bt_id, name, wrf%bt))

call check ( nf_inq_dimid(wrf%ncid, "south_north", wrf%sn_id))
call check ( nf_inq_dim(wrf%ncid, wrf%sn_id, name, wrf%sn))

call check ( nf_inq_dimid(wrf%ncid, "west_east", wrf%we_id))
call check ( nf_inq_dim(wrf%ncid, wrf%we_id, name, wrf%we))

write(6,*) ' dimensions bt, sn, we are ',wrf%bt,wrf%sn,wrf%we

!---
! get wrf variable ids and allocate space for wrf variables


call check ( nf_inq_varid(wrf%ncid, "P_TOP", wrf%ptop_id))
if(debug) write(6,*) ' ptop_id = ',wrf%ptop_id

call check ( nf_inq_varid(wrf%ncid, "U", wrf%u_id))
if(debug) write(6,*) ' u_id = ',wrf%u_id
allocate(wrf%u(wrf%we+1,wrf%sn,wrf%bt))

call check ( nf_inq_varid(wrf%ncid, "V", wrf%v_id))
if(debug) write(6,*) ' v_id = ',wrf%v_id
allocate(wrf%v(wrf%we,wrf%sn+1,wrf%bt))

call check ( nf_inq_varid(wrf%ncid, "W", wrf%w_id))
if(debug) write(6,*) ' w_id = ',wrf%w_id
allocate(wrf%w(wrf%we,wrf%sn,wrf%bt+1))

call check ( nf_inq_varid(wrf%ncid, "PH", wrf%ph_id))
if(debug) write(6,*) ' ph_id = ',wrf%ph_id
allocate(wrf%ph(wrf%we,wrf%sn,wrf%bt+1))

call check ( nf_inq_varid(wrf%ncid, "PHB", wrf%phb_id))
if(debug) write(6,*) ' phb_id = ',wrf%phb_id
allocate(wrf%phb(wrf%we,wrf%sn,wrf%bt+1))

call check ( nf_inq_varid(wrf%ncid, "T", wrf%t_id))
if(debug) write(6,*) ' t_id = ',wrf%t_id
allocate(wrf%t(wrf%we,wrf%sn,wrf%bt))

call check ( nf_inq_varid(wrf%ncid, "MU", wrf%mu_id))
if(debug) write(6,*) ' mu_id = ',wrf%mu_id
allocate(wrf%mu(wrf%we,wrf%sn))

call check ( nf_inq_varid(wrf%ncid, "MUB", wrf%mub_id))
if(debug) write(6,*) ' mub_id = ',wrf%mub_id
allocate(wrf%mub(wrf%we,wrf%sn))

call check ( nf_inq_varid(wrf%ncid, "QVAPOR", wrf%qv_id))
if(debug) write(6,*) ' qv_id = ',wrf%qv_id
allocate(wrf%qv(wrf%we,wrf%sn,wrf%bt))

call check ( nf_inq_varid(wrf%ncid, "QCLOUD", wrf%qc_id))
if(debug) write(6,*) ' qc_id = ',wrf%qc_id
allocate(wrf%qc(wrf%we,wrf%sn,wrf%bt))

call check ( nf_inq_varid(wrf%ncid, "QRAIN", wrf%qr_id))
if(debug) write(6,*) ' qr_id = ',wrf%qr_id
allocate(wrf%qr(wrf%we,wrf%sn,wrf%bt))

if (wrf%ice_micro) then

   call check ( nf_inq_varid(wrf%ncid, "QICE", wrf%qi_id))
   if(debug) write(6,*) ' qi_id = ',wrf%qi_id
   allocate(wrf%qi(wrf%we,wrf%sn,wrf%bt))

   call check ( nf_inq_varid(wrf%ncid, "QSNOW", wrf%qs_id))
   if(debug) write(6,*) ' qs_id = ',wrf%qs_id
   allocate(wrf%qs(wrf%we,wrf%sn,wrf%bt))

   call check ( nf_inq_varid(wrf%ncid, "QGRAUPEL", wrf%qg_id))
   if(debug) write(6,*) ' qg_id = ',wrf%qg_id
   allocate(wrf%qg(wrf%we,wrf%sn,wrf%bt))

end if

end subroutine wrf_open_and_alloc


!*****************************************************************************

! Internal subroutine - checks error status after each netcdf, prints 
!                       text message each time an error code is returned. 
subroutine check(istatus)
  integer, intent ( in) :: istatus 
  if(istatus /= nf_noerr) call error_handler(E_ERR,'wrf_open_and_alloc', &
       trim(nf_strerror(istatus)), source, revision, revdate) 
end subroutine check


!*****************************************************************************

subroutine dart_open_and_alloc( wrf, dart, n_values, dart_unit, &
     dart_to_wrf, debug )

implicit none

logical               :: dart_to_wrf
integer               :: dart_unit

logical               :: debug

type(wrf_data)        :: wrf
real(r8), pointer     :: dart(:)

integer               :: n_values 

!-------------------------------------------------------------
! Namelist with default values
! binary_restart_files  == .true.  -> use unformatted file format. 
!                                     Full precision, faster, smaller,
!                                     but not as portable.
! binary_restart_files  == .false.  -> use ascii file format. 
!                                     Portable, but loses precision,
!                                     slower, and larger.

logical  :: binary_restart_files = .false.

namelist /assim_nml/ binary_restart_files
!-------------------------------------------------------------

! Read the namelist input
if(file_exist('input.nml')) then
   iunit = open_file('input.nml', action = 'read')
   ierr = 1
   do while(ierr /= 0)
      read(iunit, nml = assim_nml, iostat = io, end = 11)
      ierr = check_nml_error(io, 'assim_nml')
   enddo
 11 continue
   call close_file(iunit)
endif

! namelist validation
write(*, *) 'assim_nml read; values are'
write(*, *) 'binary_restart_files is ', binary_restart_files

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
   if (binary_restart_files ) then
      open( unit=dart_unit,file="dart_wrf_vector",form="unformatted",  &
           status="old",action="read" )
   else
      open( unit=dart_unit,file="dart_wrf_vector",form="formatted",  &
           status="old",action="read" )
   endif
else
   if (binary_restart_files ) then
      open( unit=dart_unit,file="dart_wrf_vector",form="unformatted",  &
           status="new",action="write" )
   else
      open( unit=dart_unit,file="dart_wrf_vector",form="formatted",  &
           status="new",action="write" )
   endif
end if

end subroutine dart_open_and_alloc

!*****************************************************************************

subroutine dart_io( in_or_out, dart, dart_unit, n_values, &
     seconds, days, debug )

implicit none

integer               :: dart_unit, seconds, days
integer               :: seconds_end, days_end

character (len=6)     :: in_or_out
logical               :: debug

logical, parameter    :: test_input = .false.

real(r8), pointer     :: dart(:)

integer               :: n_values 

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

write(stringerror,*)' n_values differ in transfer ',n_values, n_values_in
if(n_values /= n_values_in ) call error_handler(E_ERR, 'transfer_dart_wrf', &
                             stringerror, source, revision, revdate)

end subroutine transfer_dart_wrf

!*********************************************************************

subroutine trans_2d( one_to_two, a1d, a2d, nx, ny )

implicit none

integer  :: nx,ny
real(r8) :: a1d(:)
real(r8) :: a2d(nx,ny)
logical  :: one_to_two

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

integer  :: nx,ny,nz
real(r8) :: a1d(:)
real(r8) :: a3d(nx,ny,nz)
logical  ::  one_to_three

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

END PROGRAM dart_tf_wrf
