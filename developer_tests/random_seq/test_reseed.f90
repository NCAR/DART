! DART software - Copyright UCAR. This open source software is provided
! by UCAR, "as is", without charge, subject to all terms of use at
! http://www.image.ucar.edu/DAReS/DART/DART_download
!
! $Id$

! a test of reseeding the random sequence generator, to see
! if different seeds generate overlapping sequences of values.
! one use case in the regular dart code is when running
! perfect_model_obs and generating random draws of gaussian
! error values.  if the program is restarted each time between
! separate model advances, if the same seed is used each time
! then the errors will be the same.  the code now reseeds the
! random number generator with a value based on the state data
! timestamp.  this seed seems to generate random values for
! observation errors which are not repeating.
!
! this program was also helpful in finding a bug in the
! time_manager code with one compiler.  it was generating
! long repeating sequences of values which lead to the
! discovery that the seconds were being rounded off to 0
! when setting a time type and different times would result
! in the same initial seed.

program test_reseed

use        types_mod, only : r8
use    utilities_mod, only : register_module, error_handler, E_MSG, &
                             initialize_utilities, finalize_utilities, &
                             find_namelist_in_file, check_namelist_read, &
                             nmlfileunit, do_nml_file, do_nml_term
use time_manager_mod, only : time_type, operator(+), set_time, generate_seed, &
                             set_calendar_type, print_time, print_date
use   random_seq_mod, only : random_seq_type, init_random_seq, &
                             random_uniform, random_gaussian

implicit none

! version controlled file description for error handling, do not edit
character(len=256), parameter :: source   = &
   "$URL$"
character(len=32 ), parameter :: revision = "$Revision$"
character(len=128), parameter :: revdate  = "$Date$"

type (time_type) :: state_time

integer :: io
integer :: unitnum = 66, iunit
character(len=128) :: fname
integer :: reseed, hours, years

character(len=64) :: calendar_name = 'GREGORIAN'
integer :: nsamples = 1000000
integer :: start_day = 148866
integer :: start_sec = 0

namelist /test_reseed_nml/  calendar_name, nsamples,  &
start_day, start_sec

! --------------------------------------------

call initialize_utilities('test_reseed')
call register_module(source,revision,revdate)

! Read the namelist entry
call find_namelist_in_file("input.nml", "test_reseed_nml", iunit)
read(iunit, nml = test_reseed_nml, iostat = io)
call check_namelist_read(iunit, io, "test_reseed_nml")

! Record the namelist values used for the run ...
if (do_nml_file()) write(nmlfileunit, nml=test_reseed_nml)
if (do_nml_term()) write(     *     , nml=test_reseed_nml)

! test no calendar like low order models use, also test
! the gregorian cal with a date close to current.

call set_calendar_type(calendar_name)
state_time = set_time(start_sec, start_day)

print *, ' '
call print_time(state_time, 'setting start time')
if (calendar_name /= 'NO_CALENDAR') call print_date(state_time, 'is start date')

! -----

call test1(nsamples, .true., state_time, fname = 'uniform_baseline')

call test1(nsamples, .false., state_time, fname = 'gaussian_baseline')

! -----

! these all call test2() with various intervals

call dohourstest(reseed=1, hours=1)

call dohourstest(reseed=10, hours=1)

call dohourstest(reseed=10, hours=6)

call dohourstest(reseed=10, hours=12)

call dohourstest(reseed=100, hours=6)

call dohourstest(reseed=100, hours=12)

call dohourstest(reseed=100, hours=24)

call dohourstest(reseed=1000, hours=24)

call doyearstest(reseed=1, years=1)

call doyearstest(reseed=100, years=10)


! -----

call test3()

! -----


call finalize_utilities()

contains

!-------------------------------------------------------------------------

! output uniform distribution.  seed generator once at the start.

subroutine test1(nreps, unif, start_time, fname)
 integer,               intent(in) :: nreps
 logical,               intent(in) :: unif
 type(time_type),       intent(in) :: start_time
 character(len=*),      intent(in) :: fname

integer :: i, seed, unitnum
type (random_seq_type) :: seq1

open(unit=unitnum, file=fname, action='write')

seed = generate_seed(start_time)
call init_random_seq(seq1, seed)

print *, ' '
if (calendar_name /= 'NO_CALENDAR') then
   call print_date(start_time, 'start time for long seq')
else
   call print_time(start_time, 'start time for long seq')
endif
print *, ' '

do i=1, nreps
   if (unif) then
      write(unitnum, *) random_uniform(seq1)
   else
      write(unitnum, *) random_gaussian(seq1, 0.0_r8, 1.0_r8)
   endif 
end do

close(unitnum)

end subroutine test1

!-------------------------------------------------------------------------

! output uniform distribution.  generates nrep numbers, reseeding
! the generator each 'per' times.

subroutine test2(nreps, per, start_time, delta_time, unif, fname)
 integer,               intent(in) :: nreps, per
 type(time_type),       intent(in) :: start_time, delta_time
 logical,               intent(in) :: unif
 character(len=*),      intent(in) :: fname

integer :: i, j, nextseed
type(time_type) :: t
type(random_seq_type) :: seq1

open(unit=unitnum, file=fname, action='write')

t = start_time
nextseed = generate_seed(t)
call init_random_seq(seq1, nextseed)

print *, ' '
if (calendar_name /= 'NO_CALENDAR') then
   call print_date(t, 'start time for reseed loop')
else
   call print_time(t, 'start time for reseed loop')
endif
call print_time(delta_time, 'delta time for reseed loop')
print *, 'reseeding every ', per, ' times through the loop'

j = 1
do i=1, nreps
   if (unif) then
      write(unitnum, *) random_uniform(seq1)
   else
      write(unitnum, *) random_gaussian(seq1, 0.0_r8, 1.0_r8)
   endif
   if (j == per) then
      t = t + delta_time
      nextseed = generate_seed(t)
      call init_random_seq(seq1, nextseed)
      j = 1
   else
      j = j + 1
   endif
end do

if (calendar_name /= 'NO_CALENDAR') then
   call print_date(t, 'end   time for reseed loop')
else
   call print_time(t, 'end   time for reseed loop')
endif
print *, ' '

close(unitnum)

end subroutine test2

!-------------------------------------------------------------------------

! search to see if reseeding the sequence causes it to repeat
! sequences from other seeds

subroutine test3()

real(r8), allocatable :: history(:)
real(r8) :: next_val
integer :: i, j, k, nextseed
integer, allocatable :: seedhist(:)
type(time_type) :: base_time, state_time, delta_time, delta_time2
type(random_seq_type) :: seq1


base_time = set_time(0, 148000)
state_time = base_time
delta_time = set_time(1)
delta_time2 = set_time(0, 1)

print *, ' '
if (calendar_name /= 'NO_CALENDAR') then
   call print_date(state_time, 'start  time for history loop')
else
   call print_time(state_time, 'start  time for history loop')
endif
call print_time(delta_time,  'delta  time for history loop')
call print_time(delta_time2, 'delta2 time for history loop')
print *, ' '

! generate 10000 consecutive rand nums
allocate(history(10000), seedhist(10000))

state_time = base_time
nextseed = generate_seed(state_time)
call init_random_seq(seq1, nextseed)

do i=1, 10000
   history(i) = random_uniform(seq1)
enddo

! now generate 100 values from various seeds and see if
! they overlap in any of the previously generated vals.
state_time = base_time + delta_time
do i=1, 10000
   nextseed = generate_seed(state_time)
   call init_random_seq(seq1, nextseed)

   do j=1, 100
      next_val = random_uniform(seq1)
      do k=1, 10000-1
         if (history(k) == next_val) then
            print *, 'found match, ', next_val, i, j, k
            next_val = random_uniform(seq1)
            print *, 'next val from ran ', next_val, ' next in seq is ', history(k+1)
            if (history(k+1) == next_val) then
               print *, 'found 2 matching values in a row: ', history(k), history(k+1), i, j, k
            endif
        endif
      enddo
   enddo

   state_time = state_time + delta_time

enddo

deallocate(history)

end subroutine test3

!-------------------------------------------------------------------------

subroutine dohourstest(reseed, hours)
 integer, intent(in) :: reseed, hours

type(time_type) :: delta_time

write(fname, '(A,I4.4,A,I2.2,A)') 'uniform_', reseed, 'dt_', hours, 'h'
delta_time = set_time(hours*3600, 0)

call test2(nsamples, reseed, state_time, delta_time, .true., fname)

write(fname, '(A,I4.4,A,I2.2,A)') 'gaussian_', reseed, 'dt_', hours, 'h'
delta_time = set_time(hours*3600, 0)

call test2(nsamples, reseed, state_time, delta_time, .false., fname)

end subroutine dohourstest

!-------------------------------------------------------------------------

subroutine doyearstest(reseed, years)
 integer, intent(in) :: reseed, years

type(time_type) :: delta_time

write(fname, '(A,I4.4,A,I4.4,A)') 'uniform_', reseed, 'dt_', years, 'y'
delta_time = set_time(21600, nint(years*365.25))

call test2(nsamples, reseed, state_time, delta_time, .true., fname)

write(fname, '(A,I4.4,A,I4.4,A)') 'gaussian_', reseed, 'dt_', years, 'y'
delta_time = set_time(21600, nint(years*265.25))

call test2(nsamples, reseed, state_time, delta_time, .false., fname)

end subroutine doyearstest

!-------------------------------------------------------------------------

end program test_reseed

! <next few lines under version control, do not edit>
! $URL$
! $Id$
! $Revision$
! $Date$
