! DART software - Copyright UCAR. This open source software is provided
! by UCAR, "as is", without charge, subject to all terms of use at
! http://www.image.ucar.edu/DAReS/DART/DART_download
!
! $Id$

program test_hist

! generate random values and bin them into histogram bins.
! easy enough to plot with the matlab 'bar()' function.

use      types_mod, only : r4, r8, digits12
use  utilities_mod, only : register_module, error_handler, E_ERR, E_MSG, &
                           initialize_utilities, finalize_utilities,     &
                           open_file, close_file
use random_seq_mod, only : random_seq_type, init_random_seq, random_uniform

implicit none

! version controlled file description for error handling, do not edit
character(len=256), parameter :: source   = &
   "$URL$"
character(len=32 ), parameter :: revision = "$Revision$"
character(len=128), parameter :: revdate  = "$Date$"

type (random_seq_type) :: r
real(r8) :: val
integer :: i, bin, iunit
integer, parameter :: nbins = 500
integer :: repcount = 100000000
integer :: bincount(nbins)

call initialize_utilities('test_hist')
call register_module(source,revision,revdate)

write(*, *) 'creating file "makehist.m" with output in matlab format'

iunit = open_file("makehist.m", 'formatted', 'write')
write(iunit, '(A)') 'bindata = [ ... '
bincount(:) = 0

call init_random_seq(r, 5)
do i=1, repcount
   val = random_uniform(r)
   ! generates a bin number between 0 and nbins-1, 
   ! or possibly equal to nbins w/ roundoff error.  
   ! (tests found 3 cases in 100M samples)
   bin = floor(val*nbins) + 1
   !print *, i, val, bin
   if (bin == nbins+1) bin = nbins
   if (bin < 1 .or. bin > nbins+1) then
      print *, 'error: computed bin = ', bin, ' should be > 1 and  < ', nbins
   endif
   bincount(bin) = bincount(bin)+1
enddo

do i=1, nbins
   !print *, i, bincount(i)
   !print *, (1.0_r8/nbins)*(i-1), bincount(i)
   write(iunit, '(I8,A)') bincount(i), ', ... '
enddo

write(iunit, '(A)') '];'
write(iunit, '(A)') 'bar(bindata);'
write(iunit, '(A)') 'axis([0, 500, 190000, 210000]);'

write(*, *) 'closing "makehist.m" file'

call close_file(iunit)

call finalize_utilities()

end program test_hist

! <next few lines under version control, do not edit>
! $URL$
! $Id$
! $Revision$
! $Date$
