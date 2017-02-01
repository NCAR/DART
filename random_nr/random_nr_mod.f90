! DART software - Copyright UCAR. This open source software is provided
! by UCAR, "as is", without charge, subject to all terms of use at
! http://www.image.ucar.edu/DAReS/DART/DART_download
!
! $Id$

module random_nr_mod


! THIS MODULE IS DEPRECATED.  USE THE ROUTINES IN random_seq_mod INSTEAD
!
! the old code included these public interfaces and types; the suggested 
! replacements from the random_seq_mod are:
! 
! random_seq_type -> random_seq_type 
!  unchanged: type(random_seq_type) :: r 
! 
! init_ran1(r, seed) -> init_random_seq(r, seed)
!  same usage:  call init_random_seq(r, seed)
!   where seed is optional. old seed needed to be negative or
!   less than 54000 if positive.  the new seed can be
!   any integer value.
! 
! ran1(r) -> random_uniform(r)
!  same usage: val = random_uniform(r)
! 
! gasdev(r) -> random_gaussian(r, mean, stddev)
!  gasdev took only a single arg - the type.  random gaussian takes
!  3 args:  the type, the mean, and the standard_deviation.  to get
!  the equivalent function use: 
!     val = random_gaussian(r, 0.0_r8, 1.0_r8)
!  which gives you a draw from a gaussian distribution centered on
!  0.0, with a standard deviation of 1.0 .  obviously with the 
!  new interface you can directly generate a draw from a gaussian 
!  with a different mean and/or standard deviation.
! 
! the new module also includes two other routines which may be of use:
!  several_random_gaussians(r, mean, standard_deviation, n, rnum)
!  twod_gaussians(r, mean, cov, rnum)
! see the new module html documentation file for more details.
!

end module 

! <next few lines under version control, do not edit>
! $URL$
! $Id$
! $Revision$
! $Date$
