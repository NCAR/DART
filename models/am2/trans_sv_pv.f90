! DART software - Copyright © 2004 - 2010 UCAR. This open source software is
! provided by UCAR, "as is", without charge, subject to all terms of use at
! http://www.image.ucar.edu/DAReS/DART/DART_download
 
program trans_sv_pv

! <next few lines under version control, do not edit>
! $URL$
! $Id$
! $Revision$
! $Date$

!----------------------------------------------------------------------
! purpose: interface between AM2 and DART
!
! method: Read DART state vector ("proprietary" format), but not time(s).
!         Reform state vector back into AM2 fields.
!         Replace those fields on the AM2 initial file with the new values,
!         preserving all other information on the file.
!
! author: Patrick Hofmann 3/7/08
!         based on prog_var_to_vector and vector_to_prog_var by Robert Pincus
! mod:    to read temp_ic (assim_model_state_ic; 2 times) or temp_ud (1 time) and put
!         the fields into the AM2 initial files
!
!----------------------------------------------------------------------

use       types_mod, only : r8
use   utilities_mod, only : get_unit, file_exist, open_file, &
                            initialize_utilities, finalize_utilities, check_namelist_read
use       model_mod, only : model_type, init_model_instance, write_model_init, &
   vector_to_prog_var 
use assim_model_mod, only : assim_model_type, static_init_assim_model, &
   init_assim_model, get_model_size, get_model_state_vector, read_state_restart, &
   open_restart_read, close_restart, get_model_time
use time_manager_mod, only : time_type, read_time, get_time, operator(-)

implicit none

! version controlled file description for error handling, do not edit
character(len=128), parameter :: &
   source   = "$URL$", &
   revision = "$Revision$", &
   revdate  = "$Date$"

integer, external :: iargc

type(assim_model_type) :: x
type(model_type)       :: var
type(time_type)        :: adv_to_time, curr_time, deltat
real(r8), allocatable  :: x_state(:)
integer                :: file_unit, mem_unit, x_size, sec, day, iunit, io, iunit2
character (len = 128)  :: dartSVin, RstFileOut, TrcFileOut
logical                :: do_output = .false.

! Namelist coupler_nml default values
integer :: months = 0, days = 0, hours = 0
integer :: dt_atmos = 1800, dt_ocean = 21600, dt_cpld = 21600
integer, dimension(6) :: current_date = (/2007,1,1,0,0,0/) 
character (len = 8)   :: calendar = "'julian'"
!character (len = 15)  :: current_date = '2006,1,1,0,0,0,'

namelist /coupler_nml/ &
   months, days, hours, dt_atmos, dt_ocean, dt_cpld, current_date, calendar

call initialize_utilities('Trans_sv_pv')

if(file_exist('element1')) do_output = .true.

! Static init assim model calls static_init_model
if (do_output) then
   WRITE(*,'(////A)') '========================================================================='
   PRINT*,'static_init_assim_model in trans_sv_pv'
endif
call static_init_assim_model()
call init_assim_model(x)

! Allocate the instance of the cam model type for storage
call init_model_instance(var)

if(iargc()  == 0) stop "You must specify State Vector and output AM2 files"
call getarg(1, dartSVin)
call getarg(2, RstFileOut)
call getarg(3, TrcFileOut)

! FROM KEVIN'S CAM INTERFACE:
!--------------------------------------------------------------------------------
!if (file_exist( 'temp_ic' )) then
!   file_in = 'temp_ic'
!   file_name = 'caminput.nc'
   ! Get file for DART vector input
!   file_unit = open_restart_read(file_in)
   ! read in target time and state vector from DART and throwing away the time(s)
   ! since those are handled by trans_time.f90
!   call read_state_restart(x, file_unit, adv_to_time)
!else if (file_exist( 'member' )) then
!   mem_unit = open_file ('member')
!   read(mem_unit,'(A)') file_in
!   read(mem_unit,'(A)') file_name
!   PRINT*,' file_in = ',file_in
!   file_unit = open_restart_read(file_in)
   ! read state vector from DART and throw away the time
!   call read_state_restart(x, file_unit)
!endif
!--------------------------------------------------------------------------------

file_unit = open_restart_read(dartSVin)
call read_state_restart(x, file_unit, adv_to_time)
call close_restart(file_unit)

curr_time = get_model_time(x)
deltat = adv_to_time - curr_time
call get_time(deltat,sec,day)

! Read append.nml's coupler_nml values
iunit = 5
iunit2 = 6
open(unit=iunit,file="append.nml",action='read')
open(unit=iunit2,file="newappend.nml",action='write')

read(iunit, nml = coupler_nml, iostat = io)
call check_namelist_read(iunit,io,"coupler_nml")

! Change days and hours to advance
days = day
hours = sec/3600

calendar = "'julian'"

!Write newappend.nml with new days and hours variables
write(iunit2,nml = coupler_nml, iostat = io)
close(iunit)
close(iunit2)

! Get the state part of the assim_model type x
x_size = get_model_size()
allocate(x_state(x_size))
x_state = get_model_state_vector(x)

! decompose vector back into AM2 fields
call vector_to_prog_var (x_state, var)
deallocate (x_state)

! if any of the tracer fields are negative, change them to zero
! also, restrain CF to be between 0 and 1
where(var%tracers < 0) var%tracers = 0
where(var%tracers(:,:,:,3) > 1) var%tracers(:,:,:,3) = 1 

! write fields to the netCDF initial file
call write_model_init(RstFileOut, TrcFileOut, var)

call finalize_utilities()

end program trans_sv_pv
