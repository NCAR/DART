! Data Assimilation Research Testbed -- DART
! Copyright 2004, Data Assimilation Initiative, University Corporation for Atmospheric Research
! Licensed under the GPL -- www.gpl.org/licenses/gpl.html

program ensemble_init

! <next three lines automatically updated by CVS, do not edit>
! $Source$
! $Revision$
! $Date$

! Prototype to initialize ensemble for WRF EnKF 
!  
! Reads:   Ensemble size, Ne    from    stdin 
!          Ensemble mean        from    wrfinput_mean, wrfbdy_mean 
!          ith day from climo   from    wrfinput_i, wrfbdy_i 
!                                            (e.g. wrfinput_245) [i < 10^6] 
!  
! Computes: mean of i days from climo, ith deviation from mean
!           IC's, LBC's for ith member =
!                              ensemble mean + scale * ith deviation 
!
! Writes:  ICs  for ith member   to     wrfinput_i
!          LBCs for ith member   to     wrfbdy_i
!          [Note that contents of wrfinput_i, wrfbdy_i are OVERWRITTEN.]
!
! AT PRESENT, ith member inherits "auxilliary variables" (i.e. those
! not read by wrf_io below) from ith day from climo.  To make sure
! that all mems have identical auxilliary variables, could copy 
! wrfinput_mean to wrfinput_mem_i (i=1:Ne) before executing this
! program, and have this program write state variables for ith member
! to wrfinput_mem_i rather than overwriting wrfinput_i.
!

use        types_mod, only : r8
use time_manager_mod, only : time_type, get_date, read_time, set_calendar_type, &
                             GREGORIAN, days_per_month
use  wrf_data_module, only : wrf_data, wrf_bdy_data, &
                             wrf_open_and_alloc, wrfbdy_open_and_alloc, &
                             wrf_dealloc, wrfbdy_dealloc, &
                             wrf_io, wrfbdy_io, &
                             set_wrf_date
use    utilities_mod, only : get_unit, file_exist, open_file, check_nml_error, &
                             close_file, error_handler, E_ERR, initialize_utilities, &
                             finalize_utilities, register_module, logfileunit

use                         netcdf

implicit none

! CVS Generated file description for error handling, do not edit
character(len=128) :: &
source   = "$Source$", &
revision = "$Revision$", &
revdate  = "$Date$"

!-----------------------------------------------------------------------
! Namelist parameters with default values
!-----------------------------------------------------------------------

logical :: output_state_vector = .true.       ! output prognostic variables
integer :: num_moist_vars = 0

namelist /model_nml/ output_state_vector, num_moist_vars
!-----------------------------------------------------------------------

integer :: calendar_type         = GREGORIAN
integer :: iunit, io, ierr, var_id, itime

type(wrf_data)     :: wrf, wrf_mean, wrf_tmp
type(wrf_bdy_data) :: wrf_bdy, wrf_bdy_mean, wrf_bdy_tmp

real(r8)           :: scale       ! each deviation scaled by this amt

character (len=6)  :: imem  
logical            :: debug, leap
integer            :: Ne,                 & ! Ensemble size
                      i

type(time_type)   :: dart_time(2)
integer           :: year, month, day, hour, minute, second
integer           :: ndays, m

character(len=19) :: timestring

read(5,*) Ne       ! Read ensemble size from stdin.
read(5,*) scale    ! Read scaling from stdin.

debug = .false.
!debug = .true.

call initialize_utilities
call register_module(source, revision, revdate)
write(logfileunit,*)'STARTING ensemble_init ...'

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

call set_calendar_type(calendar_type)

wrf%n_moist = num_moist_vars
wrf_mean%n_moist = num_moist_vars
wrf_tmp%n_moist = num_moist_vars
wrf_bdy%n_moist = num_moist_vars
wrf_bdy_mean%n_moist = num_moist_vars
wrf_bdy_tmp%n_moist = num_moist_vars

!-- Allocate arrays for input and bdy data
if(debug) write(6,*) ' wrf_open_and_alloc '
call wrf_open_and_alloc( wrf, 'wrfinput_1', NF90_NOWRITE, debug )
call check ( nf90_close(wrf%ncid) )
if(debug) write(6,*) ' returned from wrf_open_and_alloc '

!-- Read data to be used as ensemble mean (plus open,close netcdf file)
call wrf_open_and_alloc( wrf_mean, 'wrfinput_mean', NF90_NOWRITE, debug )
call wrf_io( wrf_mean, "INPUT ", debug )
call check ( nf90_close(wrf_mean%ncid) )

call wrf_open_and_alloc( wrf_tmp, 'wrfinput_mean', NF90_NOWRITE, debug )
call check ( nf90_close(wrf_tmp%ncid) )

call wrfbdy_open_and_alloc(wrf_bdy, 'wrfbdy_1', NF90_NOWRITE, debug )
call check ( nf90_close(wrf_bdy%ncid) )
call wrfbdy_open_and_alloc(wrf_bdy_mean, 'wrfbdy_mean', NF90_NOWRITE, debug )
call wrfbdy_io( wrf_bdy_mean, "INPUT ", debug )
call check ( nf90_close(wrf_bdy_mean%ncid) )

call wrfbdy_open_and_alloc(wrf_bdy_tmp , 'wrfbdy_mean', NF90_NOWRITE, debug )
call check ( nf90_close(wrf_bdy_tmp%ncid) )

call wrf_add   ( wrf_tmp    , 0.0_r8, wrf_mean    ,  0.0_r8)
call wrfbdy_add( wrf_bdy_tmp, 0.0_r8, wrf_bdy_mean,  0.0_r8)

!-- Compute mean over Ne input and bdy files
do i=1,Ne

   !- Open appropriate files
   write( imem , '(I6)') i
   if(debug) write(6,*) ' OPENING  wrfinput_'//adjustl(trim(imem))
   call check ( nf90_open('wrfinput_'//adjustl(trim(imem)), NF90_NOWRITE, wrf%ncid) )

   if(debug) write(6,*) ' OPENING  wrfbdy_'//adjustl(trim(imem))
   call check ( nf90_open('wrfbdy_'//adjustl(trim(imem)), NF90_NOWRITE, wrf_bdy%ncid) )

   !- Read data
   call wrf_io   ( wrf    , "INPUT ", debug )
   call wrfbdy_io( wrf_bdy, "INPUT ", debug )

   !- Close files
   call check ( nf90_close(wrf%ncid) )
   call check ( nf90_close(wrf_bdy%ncid) )

   !- accumulate sum
   call wrf_add   ( wrf_tmp    , 1.0_r8, wrf    ,  1.0_r8)
   call wrfbdy_add( wrf_bdy_tmp, 1.0_r8, wrf_bdy,  1.0_r8)

enddo

call wrf_add( wrf_tmp    , 1.0_r8/Ne, wrf    ,  0.0_r8)
call wrfbdy_add( wrf_bdy_tmp, 1.0_r8/Ne, wrf_bdy,  0.0_r8)

if (debug) write(6,*) ' --------------------'

iunit = get_unit()
open(unit = iunit, file = 'wrf.info')
dart_time(1) = read_time(iunit)
dart_time(2) = read_time(iunit)
close(iunit)

itime = 1

!-- Compute deviation from mean for each input file; overwrite input file
!   with deviation
do i=1,Ne

   !- Open appropriate files
   write( imem , '(I6)') i
   if(debug) write(6,*) ' OPENING  wrfinput_'//adjustl(trim(imem))
   call check ( nf90_open('wrfinput_'//adjustl(trim(imem)), NF90_NOWRITE, wrf%ncid) ) 
   if(debug) write(6,*) ' OPENING  wrfbdy_'//adjustl(trim(imem))
   call check ( nf90_open('wrfbdy_'//adjustl(trim(imem)), NF90_NOWRITE, wrf_bdy%ncid) ) 

   !- Read data, again
   call wrf_io   ( wrf    , "INPUT ", debug )
   call wrfbdy_io( wrf_bdy, "INPUT ", debug )

   !- Close files
   call check ( nf90_close(wrf%ncid) )
   call check ( nf90_close(wrf_bdy%ncid) )

   !- deviation from mean over input files
   call wrf_add( wrf    , 1.0_r8 , wrf_tmp    , -1.0_r8 )
   call wrfbdy_add( wrf_bdy, 1.0_r8 , wrf_bdy_tmp, -1.0_r8 )

   !- New IC: scaled deviation + chosen ensemble mean 
   call wrf_add( wrf    , scale , wrf_mean   , 1.0_r8 )
   call wrfbdy_add( wrf_bdy, scale , wrf_bdy_mean, 1.0_r8 )

   !- Open same files for writing
   if(debug) write(6,*) ' OPENING  wrfinput_'//adjustl(trim(imem))//' for WRITE'
   call check ( nf90_open('wrfinput_'//adjustl(trim(imem)), NF90_WRITE, wrf%ncid) ) 
   if(debug) write(6,*) ' OPENING  wrfbdy_'//adjustl(trim(imem))//' for WRITE'
   call check ( nf90_open('wrfbdy_'//adjustl(trim(imem)), NF90_WRITE, wrf_bdy%ncid) ) 

   !- Write IC and bdy for ith member
   call wrf_io( wrf    , "OUTPUT", debug )
   if(debug) write(6,*) 'Write boundary conditions'
   call wrfbdy_io( wrf_bdy, "OUTPUT", debug )

   call check( nf90_inq_varid(wrf_bdy%ncid, 'md___thisbdytimee_x_t_d_o_m_a_i_n_m_e_t_a_data_', var_id) )
   call check( nf90_get_var(wrf_bdy%ncid, var_id, timestring, start = (/ 1, itime /)) )
   if(debug) write(6,*) 'Original_thisbdytime = ',timestring
   call check( nf90_Redef(wrf_bdy%ncid) )
   call check( nf90_put_att(wrf_bdy%ncid, NF90_GLOBAL, "Original_thisbdytime", timestring) )
   call check( nf90_enddef(wrf_bdy%ncid) )
   call get_date(dart_time(2), year, month, day, hour, minute, second)
   call set_wrf_date(timestring, year, month, day, hour, minute, second)
   if(debug) write(6,*) 'New thisbdytime = ',timestring
   call check( nf90_put_var(wrf_bdy%ncid, var_id, timestring, start = (/ 1, itime /)) )
   call check( nf90_inq_varid(wrf_bdy%ncid, "Times", var_id) )
   if(debug) write(6,*) 'writing Times = ',timestring
   call check( nf90_put_var(wrf_bdy%ncid, var_id, timestring) )

   call check( nf90_inq_varid(wrf%ncid, "Times", var_id) )
   call check( nf90_put_var(wrf%ncid, var_id, timestring) )

   if(debug) write(6,*) 'writing START_DATE = ',timestring
   call check( nf90_put_att(wrf_bdy%ncid, nf90_global, "START_DATE", timestring) )

   call check( nf90_put_att(wrf%ncid, nf90_global, "START_DATE", timestring) )

   !- Julian year, day, and GMT correspond to the end of the next BD time???

   call get_date(dart_time(1), year, month, day, hour, minute, second)
   call check( nf90_put_att(wrf_bdy%ncid, nf90_global, "JULYR", year) )
   if(debug) write(6,*) 'writing JULYR = ',year
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
   call check( nf90_put_att(wrf_bdy%ncid, nf90_global, "JULDAY", ndays) )
   if(debug) write(6,*) 'writing JULDAY = ',ndays
   call check( nf90_inq_varid(wrf_bdy%ncid, 'md___nextbdytimee_x_t_d_o_m_a_i_n_m_e_t_a_data_', var_id) )
   call check( nf90_get_var(wrf_bdy%ncid, var_id, timestring, start = (/ 1, itime /)) )
   if(debug) write(6,*) 'Original_nextbdytime = ',timestring
   call check( nf90_Redef(wrf_bdy%ncid) )
   call check( nf90_put_att(wrf_bdy%ncid, NF90_GLOBAL, "Original_nextbdytime", timestring) )
   call check( nf90_enddef(wrf_bdy%ncid) )
   call get_date(dart_time(1), year, month, day, hour, minute, second)
   call set_wrf_date(timestring, year, month, day, hour, minute, second)
   if(debug) write(6,*) 'New nextbdytime = ',timestring
   call check( nf90_put_var(wrf_bdy%ncid, var_id, timestring, start = (/ 1, itime /)) )

   !- Close files
   call check ( nf90_close(wrf%ncid) )
   call check ( nf90_close(wrf_bdy%ncid) )

enddo

call wrf_dealloc(wrf)
call wrf_dealloc(wrf_mean)
call wrf_dealloc(wrf_tmp)
call wrfbdy_dealloc(wrf_bdy)
call wrfbdy_dealloc(wrf_bdy_mean)
call wrfbdy_dealloc(wrf_bdy_tmp)

write(logfileunit,*)'FINISHED ensemble_init.'
write(logfileunit,*)

call finalize_utilities ! closes the log file.
 
contains

  ! Internal subroutine - checks error status after each netcdf, prints 
  !                       text message each time an error code is returned. 
  subroutine check(istatus)
    integer, intent (in) :: istatus

    if(istatus /= nf90_noerr) call error_handler(E_ERR, 'ensemble_init', &
       trim(nf90_strerror(istatus)), source, revision, revdate)

  end subroutine check

!---------------------------------------------------------------

 subroutine wrf_add( wrf_a, a, wrf_b, b  )
!
! Does wrf_a = a*wrf_a + b*wrf_b, for each real array in type wrf.
! Other components of wrf_a (integers, logical) are unchanged.
!

 implicit none

 type(wrf_data), intent(inout) :: wrf_a
 type(wrf_data), intent(in   ) :: wrf_b
 real(r8), intent(in)          :: a,b

 !-- Add: wrf_a = a*wrf_a + b*wrf_b, for components 
 !       u,v,w,ph,phb,t,qv,qc,qr, qi,qs,qg, mu,mub of a and b
 wrf_a%u = a * wrf_a%u + b * wrf_b%u
 wrf_a%v = a * wrf_a%v + b * wrf_b%v
 wrf_a%w = a * wrf_a%w + b * wrf_b%w
 wrf_a%ph = a * wrf_a%ph + b * wrf_b%ph
 wrf_a%phb = a * wrf_a%phb + b * wrf_b%phb
 wrf_a%t = a * wrf_a%t + b * wrf_b%t
 wrf_a%qv = a * wrf_a%qv + b * wrf_b%qv
 wrf_a%qc = a * wrf_a%qc + b * wrf_b%qc
 wrf_a%qr = a * wrf_a%qr + b * wrf_b%qr
 wrf_a%qi = a * wrf_a%qi + b * wrf_b%qi
 wrf_a%qs = a * wrf_a%qs + b * wrf_b%qs
 wrf_a%qg = a * wrf_a%qg + b * wrf_b%qg
 wrf_a%mu = a * wrf_a%mu + b * wrf_b%mu
 wrf_a%mub = a * wrf_a%mub + b * wrf_b%mub
 
end subroutine wrf_add
!---------------------------------------------------------------

 subroutine wrfbdy_add( wrfbdy_a, a, wrfbdy_b, b  )
!
! Does wrfbdy_a = a*wrfbdy_a + b*wrfbdy_b, for each real array in type wrfbdy.
! Other components of wrf_a (integers, logical) are unchanged.
!

 implicit none

 type(wrf_bdy_data), intent(inout) :: wrfbdy_a
 type(wrf_bdy_data), intent(in   ) :: wrfbdy_b
 real(r8), intent(in)              :: a,b

 !-- Add: wrfbdy_a = a*wrfbdy_a + b*wrfbdy_b, for components 
 !       u,v,w,ph,t,qv on bndries + their tendencies

 wrfbdy_a%uxs = a * wrfbdy_a%uxs + b * wrfbdy_b%uxs
 wrfbdy_a%uxe = a * wrfbdy_a%uxe + b * wrfbdy_b%uxe
 wrfbdy_a%uys = a * wrfbdy_a%uys + b * wrfbdy_b%uys
 wrfbdy_a%uye = a * wrfbdy_a%uye + b * wrfbdy_b%uye
 wrfbdy_a%utxs = a * wrfbdy_a%utxs + b * wrfbdy_b%utxs
 wrfbdy_a%utxe = a * wrfbdy_a%utxe + b * wrfbdy_b%utxe
 wrfbdy_a%utys = a * wrfbdy_a%utys + b * wrfbdy_b%utys
 wrfbdy_a%utye = a * wrfbdy_a%utye + b * wrfbdy_b%utye

 wrfbdy_a%vxs = a * wrfbdy_a%vxs + b * wrfbdy_b%vxs
 wrfbdy_a%vxe = a * wrfbdy_a%vxe + b * wrfbdy_b%vxe
 wrfbdy_a%vys = a * wrfbdy_a%vys + b * wrfbdy_b%vys
 wrfbdy_a%vye = a * wrfbdy_a%vye + b * wrfbdy_b%vye
 wrfbdy_a%vtxs = a * wrfbdy_a%vtxs + b * wrfbdy_b%vtxs
 wrfbdy_a%vtxe = a * wrfbdy_a%vtxe + b * wrfbdy_b%vtxe
 wrfbdy_a%vtys = a * wrfbdy_a%vtys + b * wrfbdy_b%vtys
 wrfbdy_a%vtye = a * wrfbdy_a%vtye + b * wrfbdy_b%vtye

 wrfbdy_a%phxs = a * wrfbdy_a%phxs + b * wrfbdy_b%phxs
 wrfbdy_a%phxe = a * wrfbdy_a%phxe + b * wrfbdy_b%phxe
 wrfbdy_a%phys = a * wrfbdy_a%phys + b * wrfbdy_b%phys
 wrfbdy_a%phye = a * wrfbdy_a%phye + b * wrfbdy_b%phye
 wrfbdy_a%phtxs = a * wrfbdy_a%phtxs + b * wrfbdy_b%phtxs
 wrfbdy_a%phtxe = a * wrfbdy_a%phtxe + b * wrfbdy_b%phtxe
 wrfbdy_a%phtys = a * wrfbdy_a%phtys + b * wrfbdy_b%phtys
 wrfbdy_a%phtye = a * wrfbdy_a%phtye + b * wrfbdy_b%phtye

 wrfbdy_a%txs = a * wrfbdy_a%txs + b * wrfbdy_b%txs
 wrfbdy_a%txe = a * wrfbdy_a%txe + b * wrfbdy_b%txe
 wrfbdy_a%tys = a * wrfbdy_a%tys + b * wrfbdy_b%tys
 wrfbdy_a%tye = a * wrfbdy_a%tye + b * wrfbdy_b%tye
 wrfbdy_a%ttxs = a * wrfbdy_a%ttxs + b * wrfbdy_b%ttxs
 wrfbdy_a%ttxe = a * wrfbdy_a%ttxe + b * wrfbdy_b%ttxe
 wrfbdy_a%ttys = a * wrfbdy_a%ttys + b * wrfbdy_b%ttys
 wrfbdy_a%ttye = a * wrfbdy_a%ttye + b * wrfbdy_b%ttye

 wrfbdy_a%muxs = a * wrfbdy_a%muxs + b * wrfbdy_b%muxs
 wrfbdy_a%muxe = a * wrfbdy_a%muxe + b * wrfbdy_b%muxe
 wrfbdy_a%muys = a * wrfbdy_a%muys + b * wrfbdy_b%muys
 wrfbdy_a%muye = a * wrfbdy_a%muye + b * wrfbdy_b%muye
 wrfbdy_a%mutxs = a * wrfbdy_a%mutxs + b * wrfbdy_b%mutxs
 wrfbdy_a%mutxe = a * wrfbdy_a%mutxe + b * wrfbdy_b%mutxe
 wrfbdy_a%mutys = a * wrfbdy_a%mutys + b * wrfbdy_b%mutys
 wrfbdy_a%mutye = a * wrfbdy_a%mutye + b * wrfbdy_b%mutye

 wrfbdy_a%qvxs = a * wrfbdy_a%qvxs + b * wrfbdy_b%qvxs
 wrfbdy_a%qvxe = a * wrfbdy_a%qvxe + b * wrfbdy_b%qvxe
 wrfbdy_a%qvys = a * wrfbdy_a%qvys + b * wrfbdy_b%qvys
 wrfbdy_a%qvye = a * wrfbdy_a%qvye + b * wrfbdy_b%qvye
 wrfbdy_a%qvtxs = a * wrfbdy_a%qvtxs + b * wrfbdy_b%qvtxs
 wrfbdy_a%qvtxe = a * wrfbdy_a%qvtxe + b * wrfbdy_b%qvtxe
 wrfbdy_a%qvtys = a * wrfbdy_a%qvtys + b * wrfbdy_b%qvtys
 wrfbdy_a%qvtye = a * wrfbdy_a%qvtye + b * wrfbdy_b%qvtye

end subroutine wrfbdy_add

end program ensemble_init
