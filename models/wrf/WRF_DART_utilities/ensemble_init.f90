! DART software - Copyright UCAR. This open source software is provided
! by UCAR, "as is", without charge, subject to all terms of use at
! http://www.image.ucar.edu/DAReS/DART/DART_download
!
! $Id$

program ensemble_init

! This program is out of date compared to the other utilities in this directory.
! I believe it could call the standard model_mod initialization routine and
! be updated to work with the existing files, but I'm unsure about what exactly
! it does.  So for now, I'm leaving it alone - but it's not going to work
! correctly if used.     nsc 16 feb 2010

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
                             GREGORIAN, julian_day
use  wrf_data_module, only : wrf_data, wrf_bdy_data, &
                             wrf_open_and_alloc, wrfbdy_open_and_alloc, &
                             wrf_dealloc, wrfbdy_dealloc, &
                             wrf_io, wrfbdy_io, &
                             set_wrf_date
use      location_mod, only : location_type, get_location, set_location, & 
                              VERTISHEIGHT
use    utilities_mod, only : get_unit, file_exist, open_file, &
                             close_file, error_handler, E_ERR, E_MSG, initialize_utilities, &
                             register_module, logfileunit, nmlfileunit, &
                             find_namelist_in_file, check_namelist_read, finalize_utilities

use netcdf

implicit none

! version controlled file description for error handling, do not edit
character(len=256), parameter :: source   = &
   "$URL$"
character(len=32 ), parameter :: revision = "$Revision$"
character(len=128), parameter :: revdate  = "$Date$"

! miscellaneous
integer, parameter :: max_state_variables = 100
integer, parameter :: num_state_table_columns = 5
integer, parameter :: num_bounds_table_columns = 4

!-----------------------------------------------------------------------
! Model namelist parameters with default values.
!-----------------------------------------------------------------------

logical :: output_state_vector     = .false.  ! output prognostic variables
logical :: default_state_variables = .true.   ! use default state list?
character(len=129) :: wrf_state_variables(num_state_table_columns,max_state_variables) = 'NULL'
character(len=129) :: wrf_state_bounds(num_bounds_table_columns,max_state_variables) = 'NULL'
integer :: num_moist_vars       = 3
integer :: num_domains          = 1
integer :: calendar_type        = GREGORIAN
integer :: assimilation_period_seconds = 21600
logical :: surf_obs             = .true.
logical :: soil_data            = .true.
logical :: h_diab               = .false.
! Max height a surface obs can be away from the actual model surface
! and still be accepted (in meters)
real (kind=r8) :: sfc_elev_max_diff  = -1.0_r8   ! could be something like 200.0_r8
character(len = 72) :: adv_mod_command = './wrf.exe'
real (kind=r8) :: center_search_half_length = 500000.0_r8
real(r8) :: circulation_pres_level = 80000.0_r8
real(r8) :: circulation_radius     = 108000.0_r8
integer :: center_spline_grid_scale = 10
integer :: vert_localization_coord = VERTISHEIGHT
! Allow (or not) observations above the surface but below the lowest sigma level.
logical :: allow_obs_below_vol = .false.
!nc -- we are adding these to the model.nml until they appear in the NetCDF files
logical :: polar = .false.         ! wrap over the poles
logical :: periodic_x = .false.    ! wrap in longitude or x
logical :: periodic_y = .false.    ! used for single column model, wrap in y
!JPH -- single column model flag
logical :: scm        = .false.    ! using the single column model

! JPH note that soil_data and h_diab are never used and can be removed.
namelist /model_nml/ output_state_vector, num_moist_vars, &
                     num_domains, calendar_type, surf_obs, soil_data, h_diab, &
                     default_state_variables, wrf_state_variables, &
                     wrf_state_bounds, sfc_elev_max_diff, &
                     adv_mod_command, assimilation_period_seconds, &
                     allow_obs_below_vol, vert_localization_coord, &
                     center_search_half_length, center_spline_grid_scale, &
                     circulation_pres_level, circulation_radius, polar, &
                     periodic_x, periodic_y, scm

!-----------------------------------------------------------------------

integer :: iunit, io, var_id, itime

type(wrf_data)     :: wrf, wrf_mean, wrf_tmp
type(wrf_bdy_data) :: wrf_bdy, wrf_bdy_mean, wrf_bdy_tmp

real(r8)           :: scale       ! each deviation scaled by this amt

character (len=6)  :: imem
character (len=1)  :: idom
logical            :: debug = .false.
integer            :: Ne,                 & ! Ensemble size
                      i, id

type(time_type)   :: dart_time(2)
integer           :: year, month, day, hour, minute, second

character(len=19) :: timestring

read(*,*) Ne       ! Read ensemble size from stdin.
read(*,*) scale    ! Read scaling       from stdin.

call initialize_utilities('ensemble_init')
call register_module(source, revision, revdate)

! FIXME: i believe this program should call the wrf static_model_init code
! to initialize the wrf static data, and then it doesn't have to read the
! wrf model_mod namelist.  the num_moist_vars is no longer used, and there's
! an array of fields which are going to be used in the state vector that
! i believe should be used here instead.  see wrf_to_dart and dart_to_wrf
! for now this program should be restructured.

! Begin by reading the namelist input
call find_namelist_in_file("input.nml", "model_nml", iunit)
read(iunit, nml = model_nml, iostat = io)
call check_namelist_read(iunit, io, "model_nml")

write(nmlfileunit, nml=model_nml)
write(     *     , nml=model_nml)

call set_calendar_type(calendar_type)

iunit = get_unit()
open(unit = iunit, file = 'wrf.info')
dart_time(1) = read_time(iunit)
dart_time(2) = read_time(iunit)
close(iunit)

wrf%n_moist = num_moist_vars
wrf_mean%n_moist = num_moist_vars
wrf_tmp%n_moist = num_moist_vars
wrf_bdy%n_moist = num_moist_vars
wrf_bdy_mean%n_moist = num_moist_vars
wrf_bdy_tmp%n_moist = num_moist_vars

wrf%surf_obs = surf_obs
wrf_mean%surf_obs = surf_obs
wrf_tmp%surf_obs = surf_obs

wrf%h_diab = h_diab
wrf_mean%h_diab = h_diab
wrf_tmp%h_diab = h_diab

!-- First, do BC's.

call wrfbdy_open_and_alloc(wrf_bdy, 'wrfbdy_1', NF90_NOWRITE, debug )
call check ( nf90_close(wrf_bdy%ncid) )
call wrfbdy_open_and_alloc(wrf_bdy_mean, 'wrfbdy_mean', NF90_NOWRITE, debug )
call wrfbdy_io( wrf_bdy_mean, "INPUT ", debug )
call check ( nf90_close(wrf_bdy_mean%ncid) )

call wrfbdy_open_and_alloc(wrf_bdy_tmp , 'wrfbdy_mean', NF90_NOWRITE, debug )
call check ( nf90_close(wrf_bdy_tmp%ncid) )

call wrfbdy_add( wrf_bdy_tmp, 0.0_r8, wrf_bdy_mean,  0.0_r8)

!-- Compute mean over Ne input and bdy files
do i=1,Ne

   !- Open appropriate files
   write( imem , '(I6)') i

   if(debug) write(*,*) ' OPENING  wrfbdy_'//adjustl(trim(imem))
   call check ( nf90_open('wrfbdy_'//adjustl(trim(imem)), NF90_NOWRITE, wrf_bdy%ncid) )

   !- Read data
   call wrfbdy_io( wrf_bdy, "INPUT ", debug )

   !- Close files
   call check ( nf90_close(wrf_bdy%ncid) )

   !- accumulate sum
   call wrfbdy_add( wrf_bdy_tmp, 1.0_r8, wrf_bdy,  1.0_r8)

enddo

call wrfbdy_add( wrf_bdy_tmp, 1.0_r8/Ne, wrf_bdy,  0.0_r8)

if (debug) write(*,*) ' --------------------'

itime = 1

!-- Compute deviation from mean for each input file; overwrite input file with deviation
do i=1,Ne

   !- Open appropriate files
   write( imem , '(I6)') i
   if(debug) write(*,*) ' OPENING  wrfbdy_'//adjustl(trim(imem))
   call check ( nf90_open('wrfbdy_'//adjustl(trim(imem)), NF90_NOWRITE, wrf_bdy%ncid) )

   !- Read data, again
   call wrfbdy_io( wrf_bdy, "INPUT ", debug )

   !- Close files
   call check ( nf90_close(wrf_bdy%ncid) )

   !- deviation from mean over input files
   call wrfbdy_add( wrf_bdy, 1.0_r8 , wrf_bdy_tmp, -1.0_r8 )

   !- New IC: scaled deviation + chosen ensemble mean 
   call wrfbdy_add( wrf_bdy, scale , wrf_bdy_mean, 1.0_r8 )

   !- Open same files for writing
   if(debug) write(*,*) ' OPENING  wrfbdy_'//adjustl(trim(imem))//' for WRITE'
   call check ( nf90_open('wrfbdy_'//adjustl(trim(imem)), NF90_WRITE, wrf_bdy%ncid) )

   !- Write bdy for ith member
   if(debug) write(*,*) 'Write boundary conditions'
   call wrfbdy_io( wrf_bdy, "OUTPUT", debug )

   call check( nf90_inq_varid(wrf_bdy%ncid, 'md___thisbdytimee_x_t_d_o_m_a_i_n_m_e_t_a_data_', var_id) )
   call check( nf90_get_var(wrf_bdy%ncid, var_id, timestring, start = (/ 1, itime /)) )
   if(debug) write(*,*) 'Original_thisbdytime = ',timestring
   call check( nf90_Redef(wrf_bdy%ncid) )
   call check( nf90_put_att(wrf_bdy%ncid, NF90_GLOBAL, "Original_thisbdytime", timestring) )
   call check( nf90_enddef(wrf_bdy%ncid) )
   call get_date(dart_time(1), year, month, day, hour, minute, second)
   call set_wrf_date(timestring, year, month, day, hour, minute, second)
   if(debug) write(*,*) 'New thisbdytime = ',timestring
   call check( nf90_put_var(wrf_bdy%ncid, var_id, timestring, start = (/ 1, itime /)) )
   call check( nf90_inq_varid(wrf_bdy%ncid, "Times", var_id) )
   if(debug) write(*,*) 'writing Times = ',timestring
   call check( nf90_put_var(wrf_bdy%ncid, var_id, timestring) )

   if(debug) write(*,*) 'writing START_DATE = ',timestring
   call check( nf90_put_att(wrf_bdy%ncid, nf90_global, "START_DATE", timestring) )

   call check( nf90_inq_varid(wrf_bdy%ncid, 'md___nextbdytimee_x_t_d_o_m_a_i_n_m_e_t_a_data_', var_id) )
   call check( nf90_get_var(wrf_bdy%ncid, var_id, timestring, start = (/ 1, itime /)) )
   if(debug) write(*,*) 'Original_nextbdytime = ',timestring
   call check( nf90_Redef(wrf_bdy%ncid) )
   call check( nf90_put_att(wrf_bdy%ncid, NF90_GLOBAL, "Original_nextbdytime", timestring) )
   call check( nf90_enddef(wrf_bdy%ncid) )
   call get_date(dart_time(2), year, month, day, hour, minute, second)
   call set_wrf_date(timestring, year, month, day, hour, minute, second)
   if(debug) write(*,*) 'New nextbdytime = ',timestring
   call check( nf90_put_var(wrf_bdy%ncid, var_id, timestring, start = (/ 1, itime /)) )

   !- Close files
   call check ( nf90_close(wrf_bdy%ncid) )

enddo

call wrfbdy_dealloc(wrf_bdy)
call wrfbdy_dealloc(wrf_bdy_mean)
call wrfbdy_dealloc(wrf_bdy_tmp)

!-- Now do IC's for all domains.

do id=1,num_domains

   write( idom , '(I1)') id

!-- Allocate arrays for input data

   call wrf_open_and_alloc( wrf, 'wrfinput_d0'//idom//'_1', NF90_WRITE, debug )
   call check ( nf90_close(wrf%ncid) )

!-- Read data to be used as ensemble mean (plus open,close netcdf file)
   call wrf_open_and_alloc( wrf_mean, 'wrfinput_d0'//idom//'_mean', NF90_WRITE, debug )
   call wrf_io( wrf_mean, "INPUT ", debug )
   call check ( nf90_close(wrf_mean%ncid) )

   call wrf_open_and_alloc( wrf_tmp, 'wrfinput_d0'//idom//'_mean', NF90_WRITE, debug )
   call check ( nf90_close(wrf_tmp%ncid) )

   call wrf_add   ( wrf_tmp    , 0.0_r8, wrf_mean    ,  0.0_r8)

!-- Compute mean over Ne input and bdy files
   do i=1,Ne

   !- Open appropriate files
      write( imem , '(I6)') i
      call check ( nf90_open('wrfinput_d0'//idom//'_'//adjustl(trim(imem)), NF90_NOWRITE, wrf%ncid) )

   !- Read data
      call wrf_io   ( wrf    , "INPUT ", debug )

   !- Close files
      call check ( nf90_close(wrf%ncid) )

   !- accumulate sum
      call wrf_add   ( wrf_tmp    , 1.0_r8, wrf    ,  1.0_r8)

   enddo

   call wrf_add( wrf_tmp    , 1.0_r8/Ne, wrf    ,  0.0_r8)

!-- Compute deviation from mean for each input file; overwrite input file with deviation
   do i=1,Ne

   !- Open appropriate files
      write( imem , '(I6)') i
      call check ( nf90_open('wrfinput_d0'//idom//'_'//adjustl(trim(imem)), NF90_NOWRITE, wrf%ncid) )

   !- Read data, again
      call wrf_io   ( wrf    , "INPUT ", debug )

   !- Close files
      call check ( nf90_close(wrf%ncid) )

   !- deviation from mean over input files
      call wrf_add( wrf    , 1.0_r8 , wrf_tmp    , -1.0_r8 )

   !- New IC: scaled deviation + chosen ensemble mean 
      call wrf_add( wrf    , scale , wrf_mean   , 1.0_r8 )

   !- Open same files for writing
      call check ( nf90_open('wrfinput_d0'//idom//'_'//adjustl(trim(imem)), NF90_WRITE, wrf%ncid) )

   !- Write IC for ith member
      call wrf_io( wrf    , "OUTPUT", debug )

      call get_date(dart_time(1), year, month, day, hour, minute, second)
      call set_wrf_date(timestring, year, month, day, hour, minute, second)

      call check( nf90_inq_varid(wrf%ncid, "Times", var_id) )
      call check( nf90_put_var(wrf%ncid, var_id, timestring) )

      call check( nf90_put_att(wrf%ncid, nf90_global, "START_DATE", timestring) )

   !- Close files
      call check ( nf90_close(wrf%ncid) )

   enddo

   call wrf_dealloc(wrf)
   call wrf_dealloc(wrf_mean)
   call wrf_dealloc(wrf_tmp)

enddo

write(logfileunit,*)'FINISHED ensemble_init.'
write(logfileunit,*)

call finalize_utilities('ensemble_init')  ! closes log file.
 
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
 real(r8),       intent(in)    :: a,b

 !-- Add: wrf_a = a*wrf_a + b*wrf_b, for components 
 !       u,v,w,ph,phb,t, qv,qc,qr, qi,qs,qg,qnice, mu,mub, tslb,tsk, of a and b
 wrf_a%u = a * wrf_a%u + b * wrf_b%u
 wrf_a%v = a * wrf_a%v + b * wrf_b%v
 wrf_a%w = a * wrf_a%w + b * wrf_b%w
 wrf_a%ph = a * wrf_a%ph + b * wrf_b%ph
 wrf_a%phb = a * wrf_a%phb + b * wrf_b%phb
 wrf_a%t = a * wrf_a%t + b * wrf_b%t
 wrf_a%mu = a * wrf_a%mu + b * wrf_b%mu
 wrf_a%mub = a * wrf_a%mub + b * wrf_b%mub
 wrf_a%tslb = a * wrf_a%tslb + b * wrf_b%tslb
 wrf_a%tsk = a * wrf_a%tsk + b * wrf_b%tsk
 if(wrf%n_moist > 0) then
    wrf_a%qv = a * wrf_a%qv + b * wrf_b%qv
 endif
 if(wrf%n_moist > 1) then
    wrf_a%qc = a * wrf_a%qc + b * wrf_b%qc
 endif
 if(wrf%n_moist > 2) then
    wrf_a%qr = a * wrf_a%qr + b * wrf_b%qr
 endif
 if(wrf%n_moist > 3) then
    wrf_a%qi = a * wrf_a%qi + b * wrf_b%qi
 endif
 if(wrf%n_moist > 4) then
    wrf_a%qs = a * wrf_a%qs + b * wrf_b%qs
 endif
 if(wrf%n_moist > 5) then
    wrf_a%qg = a * wrf_a%qg + b * wrf_b%qg
 endif
 if(wrf%n_moist > 6) then
    wrf_a%qnice = a * wrf_a%qnice + b * wrf_b%qnice
 endif
 
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
  real(r8),           intent(in)    :: a,b

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

  if(wrfbdy_a%n_moist > 0) then
     wrfbdy_a%qvxs = a * wrfbdy_a%qvxs + b * wrfbdy_b%qvxs
     wrfbdy_a%qvxe = a * wrfbdy_a%qvxe + b * wrfbdy_b%qvxe
     wrfbdy_a%qvys = a * wrfbdy_a%qvys + b * wrfbdy_b%qvys
     wrfbdy_a%qvye = a * wrfbdy_a%qvye + b * wrfbdy_b%qvye
     wrfbdy_a%qvtxs = a * wrfbdy_a%qvtxs + b * wrfbdy_b%qvtxs
     wrfbdy_a%qvtxe = a * wrfbdy_a%qvtxe + b * wrfbdy_b%qvtxe
     wrfbdy_a%qvtys = a * wrfbdy_a%qvtys + b * wrfbdy_b%qvtys
     wrfbdy_a%qvtye = a * wrfbdy_a%qvtye + b * wrfbdy_b%qvtye
  endif

  if(wrfbdy_a%n_moist > 1) then
     wrfbdy_a%qcxs = a * wrfbdy_a%qcxs + b * wrfbdy_b%qcxs
     wrfbdy_a%qcxe = a * wrfbdy_a%qcxe + b * wrfbdy_b%qcxe
     wrfbdy_a%qcys = a * wrfbdy_a%qcys + b * wrfbdy_b%qcys
     wrfbdy_a%qcye = a * wrfbdy_a%qcye + b * wrfbdy_b%qcye
     wrfbdy_a%qctxs = a * wrfbdy_a%qctxs + b * wrfbdy_b%qctxs
     wrfbdy_a%qctxe = a * wrfbdy_a%qctxe + b * wrfbdy_b%qctxe
     wrfbdy_a%qctys = a * wrfbdy_a%qctys + b * wrfbdy_b%qctys
     wrfbdy_a%qctye = a * wrfbdy_a%qctye + b * wrfbdy_b%qctye
  endif
  if(wrfbdy_a%n_moist > 2) then
     wrfbdy_a%qrxs = a * wrfbdy_a%qrxs + b * wrfbdy_b%qrxs
     wrfbdy_a%qrxe = a * wrfbdy_a%qrxe + b * wrfbdy_b%qrxe
     wrfbdy_a%qrys = a * wrfbdy_a%qrys + b * wrfbdy_b%qrys
     wrfbdy_a%qrye = a * wrfbdy_a%qrye + b * wrfbdy_b%qrye
     wrfbdy_a%qrtxs = a * wrfbdy_a%qrtxs + b * wrfbdy_b%qrtxs
     wrfbdy_a%qrtxe = a * wrfbdy_a%qrtxe + b * wrfbdy_b%qrtxe
     wrfbdy_a%qrtys = a * wrfbdy_a%qrtys + b * wrfbdy_b%qrtys
     wrfbdy_a%qrtye = a * wrfbdy_a%qrtye + b * wrfbdy_b%qrtye
  endif
  if(wrfbdy_a%n_moist > 3) then
     wrfbdy_a%qixs = a * wrfbdy_a%qixs + b * wrfbdy_b%qixs
     wrfbdy_a%qixe = a * wrfbdy_a%qixe + b * wrfbdy_b%qixe
     wrfbdy_a%qiys = a * wrfbdy_a%qiys + b * wrfbdy_b%qiys
     wrfbdy_a%qiye = a * wrfbdy_a%qiye + b * wrfbdy_b%qiye
     wrfbdy_a%qitxs = a * wrfbdy_a%qitxs + b * wrfbdy_b%qitxs
     wrfbdy_a%qitxe = a * wrfbdy_a%qitxe + b * wrfbdy_b%qitxe
     wrfbdy_a%qitys = a * wrfbdy_a%qitys + b * wrfbdy_b%qitys
     wrfbdy_a%qitye = a * wrfbdy_a%qitye + b * wrfbdy_b%qitye
  endif
  if(wrfbdy_a%n_moist > 4) then
     wrfbdy_a%qsxs = a * wrfbdy_a%qsxs + b * wrfbdy_b%qsxs
     wrfbdy_a%qsxe = a * wrfbdy_a%qsxe + b * wrfbdy_b%qsxe
     wrfbdy_a%qsys = a * wrfbdy_a%qsys + b * wrfbdy_b%qsys
     wrfbdy_a%qsye = a * wrfbdy_a%qsye + b * wrfbdy_b%qsye
     wrfbdy_a%qstxs = a * wrfbdy_a%qstxs + b * wrfbdy_b%qstxs
     wrfbdy_a%qstxe = a * wrfbdy_a%qstxe + b * wrfbdy_b%qstxe
     wrfbdy_a%qstys = a * wrfbdy_a%qstys + b * wrfbdy_b%qstys
     wrfbdy_a%qstye = a * wrfbdy_a%qstye + b * wrfbdy_b%qstye
  endif
  if(wrfbdy_a%n_moist > 5) then
     wrfbdy_a%qgxs = a * wrfbdy_a%qgxs + b * wrfbdy_b%qgxs
     wrfbdy_a%qgxe = a * wrfbdy_a%qgxe + b * wrfbdy_b%qgxe
     wrfbdy_a%qgys = a * wrfbdy_a%qgys + b * wrfbdy_b%qgys
     wrfbdy_a%qgye = a * wrfbdy_a%qgye + b * wrfbdy_b%qgye
     wrfbdy_a%qgtxs = a * wrfbdy_a%qgtxs + b * wrfbdy_b%qgtxs
     wrfbdy_a%qgtxe = a * wrfbdy_a%qgtxe + b * wrfbdy_b%qgtxe
     wrfbdy_a%qgtys = a * wrfbdy_a%qgtys + b * wrfbdy_b%qgtys
     wrfbdy_a%qgtye = a * wrfbdy_a%qgtye + b * wrfbdy_b%qgtye
  endif
  if(wrfbdy_a%n_moist > 6) then
     wrfbdy_a%qnicexs = a * wrfbdy_a%qnicexs + b * wrfbdy_b%qnicexs
     wrfbdy_a%qnicexe = a * wrfbdy_a%qnicexe + b * wrfbdy_b%qnicexe
     wrfbdy_a%qniceys = a * wrfbdy_a%qniceys + b * wrfbdy_b%qniceys
     wrfbdy_a%qniceye = a * wrfbdy_a%qniceye + b * wrfbdy_b%qniceye
     wrfbdy_a%qnicetxs = a * wrfbdy_a%qnicetxs + b * wrfbdy_b%qnicetxs
     wrfbdy_a%qnicetxe = a * wrfbdy_a%qnicetxe + b * wrfbdy_b%qnicetxe
     wrfbdy_a%qnicetys = a * wrfbdy_a%qnicetys + b * wrfbdy_b%qnicetys
     wrfbdy_a%qnicetye = a * wrfbdy_a%qnicetye + b * wrfbdy_b%qnicetye
  endif
  if(wrfbdy_a%n_moist > 7) then
     write(*,*) 'n_moist = ',wrfbdy_a%n_moist
     call error_handler(E_ERR,'wrfbdy_add', &
          'n_moist is too large.', source, revision, revdate)
     stop
  endif

end subroutine wrfbdy_add

end program ensemble_init

! <next few lines under version control, do not edit>
! $URL$
! $Id$
! $Revision$
! $Date$
