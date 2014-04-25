! DART software - Copyright 2004 - 2013 UCAR. This open source software is
! provided by UCAR, "as is", without charge, subject to all terms of use at
! http://www.image.ucar.edu/DAReS/DART/DART_download
!
! $Id$

program test_random_gsl

! test the gaussian distribution random number generator routine
! with different means, standard deviations, and iterations.

use      types_mod, only : r8
use  utilities_mod, only : register_module, &
                           initialize_utilities, finalize_utilities
use random_seq_mod, only : random_seq_type, init_random_seq, random_gaussian

implicit none

! version controlled file description for error handling, do not edit
character(len=256), parameter :: source   = &
   "$URL$"
character(len=32 ), parameter :: revision = "$Revision$"
character(len=128), parameter :: revdate  = "$Date$"

type (random_seq_type) :: r
integer :: i, j, n
real(r8) :: r1, sqdist, mean_sqdist
real(r8) :: mean, sd

type t_inputs
 real(r8) :: m
 real(r8) :: stddev
 integer  :: nreps
end type


! to add more tests or change the parameters, specify a test count
! and update the sets of (mean, stddev, repeatcount) here:

integer, parameter :: ntests = 10
type(t_inputs) :: t(ntests) = (/ &
    t_inputs(0.0_r8, 1.0_r8,      1000), &
    t_inputs(0.0_r8, 1.0_r8,   1000000), &
    t_inputs(0.0_r8, 1.0_r8, 100000000), &
    t_inputs(1.0_r8, 5.0_r8, 100000000), &
    t_inputs(1.0_r8, 5.0_r8,     10000), &
    t_inputs(6.0_r8, 1.0_r8, 100000000), &
    t_inputs(6.0_r8, 8.0_r8,     10000), &
    t_inputs(6.0_r8, 8.0_r8,   1000000), &
    t_inputs(6.0_r8, 8.0_r8, 100000000), &
    t_inputs(0.0_r8, 8.0_r8, 100000000)  /)


call initialize_utilities('test_random')
call register_module(source,revision,revdate)

do j=1, ntests

   call init_random_seq(r, 5)

   mean = t(j)%m
   sd = t(j)%stddev
   n = t(j)%nreps

   mean_sqdist = 0.0_r8

   do i = 1, n
      r1 = random_gaussian(r, mean, sd)
      sqdist = (mean-r1)**2
      mean_sqdist = mean_sqdist + sqdist
   end do
   write(*, *) 'input mean, sd, n = ', mean, sd, n
   write(*, *) 'resulting var is ', sqrt(mean_sqdist / n)
   
enddo

call finalize_utilities()

end program test_random_gsl

! <next few lines under version control, do not edit>
! $URL$
! $Id$
! $Revision$
! $Date$
