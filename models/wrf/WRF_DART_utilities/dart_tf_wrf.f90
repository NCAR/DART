! Data Assimilation Research Testbed -- DART
! Copyright 2004-2007, Data Assimilation Research Section
! University Corporation for Atmospheric Research
! Licensed under the GPL -- www.gpl.org/licenses/gpl.html
 
PROGRAM dart_tf_wrf

! <next few lines under version control, do not edit>
! $URL$
! $Id$
! $Revision$
! $Date$

use        types_mod, only : r8
use time_manager_mod, only : time_type, write_time, read_time, get_date, set_date, operator(-), &
                             get_time, print_time, set_calendar_type, GREGORIAN, julian_day
use    utilities_mod, only : get_unit, file_exist, open_file, close_file, &
                             error_handler, E_ERR, E_MSG, initialize_utilities, &
                             register_module, logfileunit, nmlfileunit, timestamp, &
                             find_namelist_in_file, check_namelist_read
use  wrf_data_module, only : wrf_data, wrf_open_and_alloc, wrf_dealloc, wrf_io, set_wrf_date, &
                             get_wrf_date
use  assim_model_mod, only : open_restart_read, open_restart_write, aread_state_restart, &
                             awrite_state_restart
use                          netcdf

implicit none

! version controlled file description for error handling, do not edit
character(len=128), parameter :: &
   source   = "$URL$", &
   revision = "$Revision$", &
   revdate  = "$Date$"

type wrf_dom
   type(wrf_data), pointer :: dom(:)
end type wrf_dom

!-----------------------------------------------------------------------
! Model namelist parameters with default values.
!-----------------------------------------------------------------------

logical :: output_state_vector  = .false.  ! state vs. prognostic format
integer :: num_moist_vars       = 3
integer :: num_domains          = 1
integer :: calendar_type        = GREGORIAN
integer :: assimilation_period_seconds = 21600
logical :: surf_obs             = .true.
logical :: soil_data            = .true.
logical :: h_diab               = .false.
logical :: allow_obs_below_vol  = .false.
character(len = 72) :: adv_mod_command = './wrf.exe'
real (kind=r8) :: center_search_half_length = 500000.0_r8
integer :: center_spline_grid_scale = 10
integer :: vert_localization_coord  =  3  ! 1,2,3 == level,pressure,height
!nc -- we are adding these to the model.nml until they appear in the NetCDF files
logical :: polar = .false.
logical :: periodic_x = .false.


namelist /model_nml/ output_state_vector, num_moist_vars, &
                     num_domains, calendar_type, surf_obs, soil_data, h_diab, &
                     adv_mod_command, assimilation_period_seconds, &
                     allow_obs_below_vol, vert_localization_coord, &
                     center_search_half_length, center_spline_grid_scale, &
                     polar, periodic_x

!-------------------------------------------------------------

type(wrf_dom) :: wrf

real(r8), pointer :: dart(:)
type(time_type)   :: dart_time(2)
integer           :: number_dart_values, ndays, &
                     year, month, day, hour, minute, second
integer           :: ndims, idims(2), dimids(2)
integer           :: i, ivtype
character(len=80) :: varname
character(len=19) :: timestring
character(len=1)  :: idom

logical, parameter :: debug = .false.
integer            :: mode, io, var_id, id, iunit, dart_unit

logical :: dart_to_wrf

write(*,*) 'DART to WRF (.true./T) or WRF to DART (.false./F)?'

read(*,*) dart_to_wrf

call initialize_utilities('dart_tf_wrf')
call register_module(source, revision, revdate)

! Begin by reading the namelist input
call find_namelist_in_file("input.nml", "model_nml", iunit)
read(iunit, nml = model_nml, iostat = io)
call check_namelist_read(iunit, io, "model_nml")

! Record the namelist values used for the run ...
write(nmlfileunit, nml=model_nml)
write(     *     , nml=model_nml)

call set_calendar_type(calendar_type)

if ( dart_to_wrf ) then
   call error_handler(E_MSG,'dart_to_wrf', &
       'Converting a dart state vector to a WRF netcdf file', &
       source, revision, revdate)
else
   call error_handler(E_MSG,'dart_to_wrf', &
       'Converting a WRF netcdf file to a dart state vector', &
       source, revision, revdate)
endif

allocate(wrf%dom(num_domains))

wrf%dom(:)%n_moist = num_moist_vars
wrf%dom(:)%surf_obs = surf_obs
wrf%dom(:)%soil_data = soil_data
wrf%dom(:)%h_diab = h_diab

! open wrf data netCDF file 'wrfinput_d0x'
! we get sizes of the WRF geometry and resolution

mode = NF90_NOWRITE                   ! read the netcdf file
if( dart_to_wrf ) mode = NF90_WRITE   ! write to the netcdf file

if(debug) write(*,*) ' wrf_open_and_alloc '
do id=1,num_domains
   write( idom , '(I1)') id
   call wrf_open_and_alloc( wrf%dom(id), 'wrfinput_d0'//idom, mode, debug )
enddo
if(debug) write(*,*) ' returned from wrf_open_and_alloc '

!---
! allocate space for DART data

if(debug) write(*,*) ' dart_open_and_alloc '

call dart_open_and_alloc( wrf, dart, number_dart_values, dart_unit, dart_to_wrf, &
     debug  )
if(debug) write(*,*) ' returned from dart_open_and_alloc '

!----------------------------------------------------------------------
!  get DART or WRF data

if(debug) write(*,*) ' state input '

if( dart_to_wrf ) then

   call aread_state_restart(dart_time(2), dart, dart_unit, dart_time(1))

   iunit = get_unit()
   open(unit = iunit, file = 'wrf.info')
   call write_time(iunit, dart_time(1))
   call write_time(iunit, dart_time(2))
   call get_date(dart_time(2), year, month, day, hour, minute, second)
   write (iunit,FMT='(I4,5I3.2)') year, month, day, hour, minute, second

   write (iunit,*) num_domains
   write (iunit,*) adv_mod_command
   close(iunit)

else

   do id=1,num_domains
      call WRF_IO( wrf%dom(id), "INPUT ", debug )
   enddo

end if
if(debug) write(*,*) ' returned from state input '

!---
!  translate from DART to WRF, or WRF to DART

if(debug) write(*,*) ' transfer data to_from dart-wrf '

call transfer_dart_wrf ( dart_to_wrf, dart, wrf,    &
     number_dart_values )

if(debug) write(*,*) ' transfer complete '

!---
!  output

if(debug) write(*,*) ' state output '
if( dart_to_wrf ) then

!!$   call get_date(dart_time(2), year, month, day, hour, minute, second)
   call set_wrf_date(timestring, year, month, day, hour, minute, second)
   ndays = julian_day(year, month, day)

   do id=1,num_domains
      call check( nf90_inq_varid(wrf%dom(id)%ncid, "Times", var_id) )
      call check( nf90_put_var(wrf%dom(id)%ncid, var_id, timestring) )
      call check( nf90_put_att(wrf%dom(id)%ncid, nf90_global, "START_DATE", timestring) )
      call check( nf90_put_att(wrf%dom(id)%ncid, nf90_global, "JULYR", year) )
      call check( nf90_put_att(wrf%dom(id)%ncid, nf90_global, "JULDAY", ndays) )
      call WRF_IO( wrf%dom(id), "OUTPUT", debug )
   enddo

else
   iunit = get_unit()
   call check( nf90_inq_varid(wrf%dom(1)%ncid, "Times", var_id) )
   call check( nf90_Inquire_Variable(wrf%dom(1)%ncid, var_id, varname, xtype=ivtype, &
        ndims=ndims, dimids=dimids) )
   do i=1,ndims
      call check( nf90_inquire_dimension(wrf%dom(1)%ncid, dimids(i), len=idims(i)) )
      if(debug) write(*,*) ' dimension ',i,idims(i)
   enddo

   call check( nf90_get_var(wrf%dom(1)%ncid, var_id, timestring, start = (/ 1, idims(2) /)) )
   call get_wrf_date(timestring, year, month, day, hour, minute, second)
   dart_time(1) = set_date(year, month, day, hour, minute, second)

   call print_time(dart_time(1),str='Time from wrfinput_d0x:')

   if(file_exist('wrf.info')) then
      open(unit = iunit, file = 'wrf.info')
      dart_time(1) = read_time(iunit)
      close(iunit)
   endif

   call print_time(dart_time(1),str='Time written to dart vector file:')

   call awrite_state_restart(dart_time(1), dart, dart_unit)

end if
if(debug) write(*,*) ' returned from state output '

do id=1,num_domains
   call check ( nf90_sync(wrf%dom(id)%ncid) )
   call check ( nf90_close(wrf%dom(id)%ncid) )
enddo

do id=1,num_domains
   call wrf_dealloc(wrf%dom(id))
enddo
deallocate(dart)

write(logfileunit,*)'FINISHED dart_tf_wrf.'
write(logfileunit,*)

call timestamp(source,revision,revdate,'end') ! That closes the log file, too.
 
contains

  ! Internal subroutine - checks error status after each netcdf, prints 
  !                       text message each time an error code is returned. 
  subroutine check(istatus)
    integer, intent ( in) :: istatus 
    if(istatus /= nf90_noerr) call error_handler(E_ERR,'dart_tf_wrf', &
         trim(nf90_strerror(istatus)), source, revision, revdate) 
  end subroutine check

!*****************************************************************************

subroutine dart_open_and_alloc( wrf, dart, n_values, dart_unit, dart_to_wrf, &
     debug )

implicit none

integer,       intent(out) :: dart_unit

logical,       intent(in)  :: dart_to_wrf, debug

type(wrf_dom), intent(in)  :: wrf
real(r8),      pointer     :: dart(:)

integer,       intent(out) :: n_values 

integer :: id

character(len=80) :: stringerror

! compute number of values in 1D vector

n_values = 0

do id=1,num_domains

! dry dynamics conponents

   n_values = n_values + (wrf%dom(id)%bt  )*(wrf%dom(id)%sn  )*(wrf%dom(id)%we+1)  ! u
   n_values = n_values + (wrf%dom(id)%bt  )*(wrf%dom(id)%sn+1)*(wrf%dom(id)%we  )  ! v
   n_values = n_values + (wrf%dom(id)%bt+1)*(wrf%dom(id)%sn  )*(wrf%dom(id)%we  )  ! w
   n_values = n_values + (wrf%dom(id)%bt+1)*(wrf%dom(id)%sn  )*(wrf%dom(id)%we  )  ! geopotential
   n_values = n_values + (wrf%dom(id)%bt  )*(wrf%dom(id)%sn  )*(wrf%dom(id)%we  )  ! t
   n_values = n_values +                    (wrf%dom(id)%sn  )*(wrf%dom(id)%we  )  ! dry surf. press.
   n_values = n_values +                    (wrf%dom(id)%sn  )*(wrf%dom(id)%we  )  ! skin temperature

! moist variables. Order is qv, qc, qr, qi, qs, qg, qnice.

   if(wrf%dom(id)%n_moist > 7) then
      write(stringerror,*) 'n_moist = ',wrf%dom(id)%n_moist,' is too large.'
      call error_handler(E_ERR, 'dart_open_and_alloc', &
           stringerror, source, revision, revdate)
   else
      n_values = n_values + wrf%dom(id)%n_moist*(wrf%dom(id)%bt)*(wrf%dom(id)%sn)*(wrf%dom(id)%we)
   endif

   if( wrf%dom(id)%surf_obs ) then
      n_values = n_values + 6 * wrf%dom(id)%sn * wrf%dom(id)%we
   endif

   if( wrf%dom(id)%soil_data ) then
      n_values = n_values + 3 * (wrf%dom(id)%sls )*(wrf%dom(id)%sn  )*(wrf%dom(id)%we  )  ! tslb, smois, sh2o
   endif

   if( wrf%dom(id)%h_diab ) then
      n_values = n_values + (wrf%dom(id)%bt  )*(wrf%dom(id)%sn  )*(wrf%dom(id)%we  )
   endif

enddo

if(debug) write(*,*) ' dart vector length is ',n_values

allocate(dart(n_values))

!  open DART data file

if(dart_to_wrf)  then  !  DART data file should exist, open it
   dart_unit = open_restart_read("dart_wrf_vector")
else
   dart_unit = open_restart_write("dart_wrf_vector")
end if

end subroutine dart_open_and_alloc

!*****************************************************************************

subroutine transfer_dart_wrf ( dart_to_wrf, dart, wrf, n_values_in)

implicit none

logical,          intent(in)    :: dart_to_wrf

type(wrf_dom),    intent(inout) :: wrf
real(r8), pointer               :: dart(:)

integer,          intent(in)    :: n_values_in

!---

integer           :: in, n_values,id, end_moist
character(len=80) :: stringerror

!---

n_values = 0

do id=1,num_domains

   in = n_values+1
   call trans_3d( dart_to_wrf, dart(in:),wrf%dom(id)%u,wrf%dom(id)%we+1,wrf%dom(id)%sn,wrf%dom(id)%bt)
   n_values = n_values + (wrf%dom(id)%bt  )*(wrf%dom(id)%sn  )*(wrf%dom(id)%we+1)  ! u

   in = n_values+1
   call trans_3d( dart_to_wrf, dart(in:),wrf%dom(id)%v,wrf%dom(id)%we,wrf%dom(id)%sn+1,wrf%dom(id)%bt)
   n_values = n_values + (wrf%dom(id)%bt  )*(wrf%dom(id)%sn+1)*(wrf%dom(id)%we  )  ! v

   in = n_values+1
   call trans_3d( dart_to_wrf, dart(in:),wrf%dom(id)%w,wrf%dom(id)%we,wrf%dom(id)%sn,wrf%dom(id)%bt+1)
   n_values = n_values + (wrf%dom(id)%bt+1)*(wrf%dom(id)%sn  )*(wrf%dom(id)%we  )  ! w

   in = n_values+1
   call trans_3d( dart_to_wrf, dart(in:),wrf%dom(id)%ph,wrf%dom(id)%we,wrf%dom(id)%sn,wrf%dom(id)%bt+1)
   n_values = n_values + (wrf%dom(id)%bt+1)*(wrf%dom(id)%sn  )*(wrf%dom(id)%we  )  ! geopotential

   in = n_values+1
   call trans_3d( dart_to_wrf, dart(in:),wrf%dom(id)%t,wrf%dom(id)%we,wrf%dom(id)%sn,wrf%dom(id)%bt)
   n_values = n_values + (wrf%dom(id)%bt  )*(wrf%dom(id)%sn  )*(wrf%dom(id)%we  )  ! t

   in = n_values+1
   call trans_2d( dart_to_wrf, dart(in:),wrf%dom(id)%mu,wrf%dom(id)%we,wrf%dom(id)%sn)
   n_values = n_values +                    (wrf%dom(id)%sn  )*(wrf%dom(id)%we  )  ! dry surf. press.

   in = n_values+1
   call trans_2d( dart_to_wrf, dart(in:),wrf%dom(id)%tsk,wrf%dom(id)%we,wrf%dom(id)%sn)
   n_values = n_values +                    (wrf%dom(id)%sn  )*(wrf%dom(id)%we  )  ! skin temperature

! moist variables

   if(wrf%dom(id)%n_moist >= 1) then
      in = n_values+1

!--------Make sure that microphysics variables are not negatives.
      
      end_moist = n_values + (wrf%dom(id)%n_moist)*(wrf%dom(id)%bt)*(wrf%dom(id)%sn)*(wrf%dom(id)%we)
      dart(in:end_moist) = max(0.0_r8,dart(in:end_moist))

      call trans_3d( dart_to_wrf, dart(in:),wrf%dom(id)%qv,wrf%dom(id)%we,wrf%dom(id)%sn,wrf%dom(id)%bt)
      n_values = n_values + (wrf%dom(id)%bt  )*(wrf%dom(id)%sn  )*(wrf%dom(id)%we  )  ! qv
   endif
   if(wrf%dom(id)%n_moist >= 2) then
      in = n_values+1
      call trans_3d( dart_to_wrf, dart(in:),wrf%dom(id)%qc,wrf%dom(id)%we,wrf%dom(id)%sn,wrf%dom(id)%bt)
      n_values = n_values + (wrf%dom(id)%bt  )*(wrf%dom(id)%sn  )*(wrf%dom(id)%we  )  ! qc
   endif
   if(wrf%dom(id)%n_moist >= 3) then
      in = n_values+1
      call trans_3d( dart_to_wrf, dart(in:),wrf%dom(id)%qr,wrf%dom(id)%we,wrf%dom(id)%sn,wrf%dom(id)%bt)
      n_values = n_values + (wrf%dom(id)%bt  )*(wrf%dom(id)%sn  )*(wrf%dom(id)%we  )  ! qr
   endif
   if(wrf%dom(id)%n_moist >= 4) then
      in = n_values+1
      call trans_3d( dart_to_wrf, dart(in:),wrf%dom(id)%qi,wrf%dom(id)%we,wrf%dom(id)%sn,wrf%dom(id)%bt)
      n_values = n_values + (wrf%dom(id)%bt  )*(wrf%dom(id)%sn  )*(wrf%dom(id)%we  )  ! qi
   endif
   if(wrf%dom(id)%n_moist >= 5) then
      in = n_values+1
      call trans_3d( dart_to_wrf, dart(in:),wrf%dom(id)%qs,wrf%dom(id)%we,wrf%dom(id)%sn,wrf%dom(id)%bt)
      n_values = n_values + (wrf%dom(id)%bt  )*(wrf%dom(id)%sn  )*(wrf%dom(id)%we  )  ! qs
   endif
   if(wrf%dom(id)%n_moist >= 6) then
      in = n_values+1
      call trans_3d( dart_to_wrf, dart(in:),wrf%dom(id)%qg,wrf%dom(id)%we,wrf%dom(id)%sn,wrf%dom(id)%bt)
      n_values = n_values + (wrf%dom(id)%bt  )*(wrf%dom(id)%sn  )*(wrf%dom(id)%we  )  ! qg
   endif
   if(wrf%dom(id)%n_moist == 7) then
      in = n_values+1
      call trans_3d( dart_to_wrf, dart(in:),wrf%dom(id)%qnice,wrf%dom(id)%we,wrf%dom(id)%sn,wrf%dom(id)%bt)
      n_values = n_values + (wrf%dom(id)%bt  )*(wrf%dom(id)%sn  )*(wrf%dom(id)%we  )  ! qnice
   endif
   if(wrf%dom(id)%n_moist > 7) then
      write(stringerror,*) 'n_moist = ',wrf%dom(id)%n_moist,' is too large.'
      call error_handler(E_ERR, 'transfer_dart_wrf', &
           stringerror, source, revision, revdate)
   endif

   if( wrf%dom(id)%surf_obs ) then

      in = n_values+1
      call trans_2d( dart_to_wrf, dart(in:),wrf%dom(id)%u10,wrf%dom(id)%we,wrf%dom(id)%sn)
      n_values = n_values + (wrf%dom(id)%sn  )*(wrf%dom(id)%we  )  ! u10

      in = n_values+1
      call trans_2d( dart_to_wrf, dart(in:),wrf%dom(id)%v10,wrf%dom(id)%we,wrf%dom(id)%sn)
      n_values = n_values + (wrf%dom(id)%sn  )*(wrf%dom(id)%we  )  ! v10

      in = n_values+1
      call trans_2d( dart_to_wrf, dart(in:),wrf%dom(id)%t2,wrf%dom(id)%we,wrf%dom(id)%sn)
      n_values = n_values + (wrf%dom(id)%sn  )*(wrf%dom(id)%we  )  ! t2

      in = n_values+1
      call trans_2d( dart_to_wrf, dart(in:),wrf%dom(id)%th2,wrf%dom(id)%we,wrf%dom(id)%sn)
      n_values = n_values + (wrf%dom(id)%sn  )*(wrf%dom(id)%we  )  ! th2

      in = n_values+1
      call trans_2d( dart_to_wrf, dart(in:),wrf%dom(id)%q2,wrf%dom(id)%we,wrf%dom(id)%sn)
      n_values = n_values + (wrf%dom(id)%sn  )*(wrf%dom(id)%we  )  ! q2

      in = n_values+1
      call trans_2d( dart_to_wrf, dart(in:),wrf%dom(id)%ps,wrf%dom(id)%we,wrf%dom(id)%sn)
      n_values = n_values + (wrf%dom(id)%sn  )*(wrf%dom(id)%we  )  ! ps

   endif

   if( wrf%dom(id)%soil_data ) then

     in = n_values+1
     call trans_3d( dart_to_wrf, dart(in:),wrf%dom(id)%tslb,wrf%dom(id)%we,wrf%dom(id)%sn,wrf%dom(id)%sls)
     n_values = n_values + (wrf%dom(id)%sls )*(wrf%dom(id)%sn  )*(wrf%dom(id)%we  )  ! tslb

     in = n_values+1
     call trans_3d( dart_to_wrf, dart(in:),wrf%dom(id)%smois,wrf%dom(id)%we,wrf%dom(id)%sn,wrf%dom(id)%sls)
     n_values = n_values + (wrf%dom(id)%sls )*(wrf%dom(id)%sn  )*(wrf%dom(id)%we  )  ! smois
     
     in = n_values+1
     call trans_3d( dart_to_wrf, dart(in:),wrf%dom(id)%sh2o,wrf%dom(id)%we,wrf%dom(id)%sn,wrf%dom(id)%sls)
     n_values = n_values + (wrf%dom(id)%sls )*(wrf%dom(id)%sn  )*(wrf%dom(id)%we  )  ! sh2o    

   endif

   if( wrf%dom(id)%h_diab ) then

      in = n_values+1
      call trans_3d( dart_to_wrf, dart(in:),wrf%dom(id)%hdiab,wrf%dom(id)%we,wrf%dom(id)%sn,wrf%dom(id)%bt)
      n_values = n_values + (wrf%dom(id)%bt  )*(wrf%dom(id)%sn  )*(wrf%dom(id)%we  )

   endif

enddo

if(n_values /= n_values_in ) then
   write(stringerror,*)' n_values differ in transfer ',n_values, n_values_in
   call error_handler(E_ERR, 'transfer_dart_wrf', &
                             stringerror, source, revision, revdate)
endif

end subroutine transfer_dart_wrf

!*********************************************************************

subroutine trans_2d( one_to_two, a1d, a2d, nx, ny )

implicit none

integer,  intent(in)    :: nx,ny
real(r8), intent(inout) :: a1d(:)
real(r8), intent(inout) :: a2d(nx,ny)
logical,  intent(in)    :: one_to_two

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

integer,  intent(in)    :: nx,ny,nz
real(r8), intent(inout) :: a1d(:)
real(r8), intent(inout) :: a3d(:,:,:)
logical,  intent(in)    :: one_to_three

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
   write(errstring,*)'nx, ny, nz, not compatible ',i,j,k,nx,ny,nz,m
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
