! DART software - Copyright 2004 - 2013 UCAR. This open source software is
! provided by UCAR, "as is", without charge, subject to all terms of use at
! http://www.image.ucar.edu/DAReS/DART/DART_download
!
! $Id$

program test_hist

! generate random values and bin them into histogram bins.
! easy enough to plot with the matlab 'bar()' function.

use      types_mod, only : r4, r8, digits12
use  utilities_mod, only : register_module, error_handler, E_ERR, E_MSG, &
                           initialize_utilities, finalize_utilities
use random_seq_mod, only : random_seq_type, init_random_seq, random_uniform

implicit none

! version controlled file description for error handling, do not edit
character(len=256), parameter :: source   = &
   "$URL$"
character(len=32 ), parameter :: revision = "$Revision$"
character(len=128), parameter :: revdate  = "$Date$"

type (random_seq_type) :: r
real(r8) :: val
integer :: i, bin
integer, parameter :: nbins = 500
integer :: repcount = 100000000
integer :: bincount(nbins)

call initialize_utilities('test_hist')
call register_module(source,revision,revdate)

bincount(:) = 0

call init_random_seq(r, 5)
do i=1, repcount
   val = random_uniform(r)
   bin = int(val*nbins)
   !print *, i, val, bin
   bincount(bin) = bincount(bin)+1
enddo

do i=1, nbins
   !print *, i, bincount(i)
   !print *, (1.0_r8/nbins)*(i-1), bincount(i)
   print *, bincount(i)
enddo

call error_handler(E_MSG, 'test_hist', 'Finished successfully.',&
                   source,revision,revdate)
call finalize_utilities()

end program test_hist

! <next few lines under version control, do not edit>
! $URL$
! $Id$
! $Revision$
! $Date$
