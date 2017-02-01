! DART software - Copyright UCAR. This open source software is provided
! by UCAR, "as is", without charge, subject to all terms of use at
! http://www.image.ucar.edu/DAReS/DART/DART_download
!
! $Id$

program trans_MITgcm_to_dart

! As the name suggests, this program translates MITgcm restart files
! to a file that dart can use.

use        types_mod, only : r8
use time_manager_mod, only : time_type, write_time, read_time, get_date,  &
                             set_date, operator(-), get_time, print_time, &
                             set_calendar_type, GREGORIAN, julian_day,    &
                             set_time, print_time
use    utilities_mod, only : get_unit, error_handler,  &
                             E_ERR, E_MSG, initialize_utilities,          &
                             finalize_utilities, register_module,         &
                             find_namelist_in_file, check_namelist_read,  &
                             logfileunit, nmlfileunit

implicit none

! version controlled file description for error handling, do not edit
character(len=256), parameter :: source   = &
   "$URL$"
character(len=32 ), parameter :: revision = "$Revision$"
character(len=128), parameter :: revdate  = "$Date$"

!-----------------------------------------------------------------------
! Model namelist parameters with default values
!-----------------------------------------------------------------------

integer  :: model_size        = 539400
integer  :: naz               = 120 
integer  :: nrad              = 31
integer  :: nzed              = 29
integer  :: ntype             = 5
real(r8) :: daz               = 3.00_r8
real(r8) :: drad              = 0.01_r8
real(r8) :: dzed              = 0.005_r8
real(r8) :: inner_rad         = 0.08_r8
real(r8) :: outer_rad         = 0.3_r8
real(r8) :: depth             = 0.14_r8
real(r8) :: delta_t           = 0.1_r8
integer  :: time_step_days    = 0
integer  :: time_step_seconds = 2160

namelist /model_nml/ model_size, naz, nrad, nzed, ntype, daz, drad, dzed, &
   inner_rad, outer_rad, depth, delta_t, time_step_days, time_step_seconds

!-------------------------------------------------------------
! misc local variables
real(r8), allocatable    :: dart(:)
real(r8), allocatable    :: r8seg(:)
type(time_type)          :: dart_time(2)
! type(time_type)          :: dt_temp
integer                  :: icount, j, k, irec
integer                  :: iunit, io, dart_unit
logical, parameter       :: debug = .false.
integer                  :: days, seconds

!-------------------------------------------------------------
! Namelist with default values
! binary_restart_files  == .true.  -> use unformatted file format.
!                                     Full precision, faster, smaller,
!                                     but not as portable.
! binary_restart_files  == .false.  -> use ascii file format.
!                                     Portable, but loses precision,
!                                     slower, and larger.

!QQQ why is this set here and not passed in some way?
logical  :: binary_restart_files = .true.

namelist /assim_model_nml/ binary_restart_files

call initialize_utilities
call register_module(source, revision, revdate)
write(logfileunit,*)'STARTING trans_MITgcm_to_dart ...'

! Read the namelist entry
call find_namelist_in_file("input.nml", "assim_model_nml", iunit)
read(iunit, nml = assim_model_nml, iostat = io)
call check_namelist_read(iunit, io, "assim_model_nml")
! Read the namelist entry
call find_namelist_in_file("input.nml", "model_nml", iunit)
read(iunit, nml = model_nml, iostat = io)
call check_namelist_read(iunit, io, "model_nml")

! Record the namelist values used for the run ...
write(nmlfileunit, nml=assim_model_nml)
write(     *     , nml=assim_model_nml)

write(nmlfileunit, nml=model_nml)
write(     *     , nml=model_nml)

call error_handler(E_MSG,'MITgcm_to_dart',                   &
       'Converting an MITgcm restart file into a dart file', &
       source, revision, revdate)

! allocate space for dart vector
allocate(dart(model_size))

! allocate space for vector used in reading and writing data
allocate(r8seg(naz))

! open dart data file for writing
dart_unit = get_unit()
if ( binary_restart_files ) then
   open( unit=dart_unit,file="dart_vector",form="unformatted",  &
         status="new",action="write" )
else
   open( unit=dart_unit,file="dart_vector",form="formatted",    &
         status="new",action="write" )
endif

! open the pickup file for the MITgcm non-hydrostatic variables
open(unit=3,file='pickup.in.s',status='old',access='direct',recl=naz*8)

! set counter for state
icount = 1

! read in u (uVel in MITgcm-speak)
do k = 1, nzed
   do j = 1, nrad
      irec = (k - 1)*nrad + j
      read(3,rec = irec) r8seg       
      dart(icount:icount + naz - 1) = r8seg(1:naz)
      icount = icount + naz
   end do
end do

! skip u tendencies (gU and gUnm1 in MITgcm-speak)
do k = 1, 2*nzed
   do j = 1, nrad
      irec = (k - 1)*nrad + j + nzed*nrad
      read(3,rec = irec)  r8seg
   end do
end do

! read in v (vVel in MITgcm-speak)
do k = 1, nzed
   do j = 1, nrad
      irec = (k - 1)*nrad + j + nzed*nrad*3
      read(3,rec = irec) r8seg
      dart(icount:icount + naz - 1) = r8seg(1:naz)
      icount = icount + naz
   end do     
end do

! skip v tendencies (gV and gVnm1 in MITgcm-speak)
do k = 1, 2*nzed
   do j = 1, nrad
      irec = (k - 1)*nrad + j + nzed*nrad*4
      read(3,rec = irec)  r8seg
   end do
end do

! read in w and T
do k = 1, 2*nzed
   do j = 1, nrad
      irec = (k - 1)*nrad + j + nzed*nrad*6
      read(3,rec = irec) r8seg
      dart(icount:icount + naz - 1) = r8seg(1:naz)
      icount = icount + naz
   end do
end do

! skip T tendencies (gT and gTnm1 in MITgcm-speak)
do k = 1, 2*nzed
   do j = 1, nrad
      irec = (k - 1)*nrad + j + nzed*nrad*8
      read(3,rec = irec)  r8seg
   end do
end do

! skip salinity (salt in MITgcm-speak)
do k = 1, nzed
   do j = 1, nrad
      irec = (k - 1)*nrad + j + nzed*nrad*10
      read(3,rec = irec) r8seg
   end do
end do

! skip salinity tendencies (gS and gSnm1 in MITgcm-speak)
do k = 1, 2*nzed
   do j = 1, nrad
      irec = (k - 1)*nrad + j + nzed*nrad*11
      read(3,rec = irec)  r8seg
   end do
end do

! skip free surface height (etaN in MITgcm-speak)
do j = 1,nrad
   irec =  j + nzed*nrad*13
   read(3,rec = irec)  r8seg
end do

! close the files
close(unit=3)
close(unit=4)

! now open the pickup file for the hydrostatic variables
open(unit=3,file='pickup_nh.in.s',status='old',access='direct',recl=naz*8)

! read in p (phi_nh in MITgcm-speak)
do k = 1, nzed
   do j = 1, nrad
      irec = (k - 1)*nrad + j 
      read(3,rec = irec) r8seg
      dart(icount:icount + naz - 1) = r8seg(1:naz)
      icount = icount + naz
   end do
end do

! skip the w tendencies (gW and gWnm1 in MITgcm-speak)
do k = 1, 2*nzed
   do j = 1, nrad
      irec = (k - 1)*nrad + j
      read(3,rec = irec)  r8seg
   end do
end do

! close the files
close(unit=3)
close(unit=4)

! write the contents of the second element of dart_time to a file
! for use by trans_MITgcm_to_dart
open(unit=777,file="dart_time.dat",status="old")
read(777,*) days, seconds
dart_time(1) = set_time(seconds,days)
!read(777) dt_temp%seconds
!read(777) dt_temp%days
!dart_time(1)=dt_temp
call print_time(dart_time(1),'dart_time in trans_MITgcm_to_dart')

close(777)

! write the contents of the dart vector into the dart file
call dart_io( "OUTPUT", dart, dart_unit, dart_time, binary_restart_files )

! clean up and exit
write(logfileunit,*)'FINISHED trans_MIGgcm_to_dart'
write(logfileunit,*)      
call finalize_utilities

contains

subroutine dart_io( in_or_out, dart, dart_unit, dart_time, binary_restart_files )

! This subroutine reads and writes information from/to the DART
! file format.  if in_or_out is 'INPUT', then the dart file is
! read into the dart vector, otherwise the contents of the dart
! vector are written to a dart file (specified by dart_unit).
   
implicit none

character (len=6), intent(in)    :: in_or_out
real(r8)                         :: dart(:)

integer,           intent(in)    :: dart_unit
type(time_type),   intent(inout) :: dart_time(2)
logical,           intent(in)    :: binary_restart_files

! This is wild, the ic files have two lines of time at
! the top (current and desired), while the ud files have
! only one line of time (current)
if (in_or_out(1:5) == 'INPUT') then
   if ( binary_restart_files ) then
      dart_time(1) = read_time(dart_unit, "unformatted")
      dart_time(2) = read_time(dart_unit, "unformatted")
   else
      dart_time(1) = read_time(dart_unit)
      dart_time(2) = read_time(dart_unit)
   endif
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



end program trans_MITgcm_to_dart

! <next few lines under version control, do not edit>
! $URL$
! $Id$
! $Revision$
! $Date$
