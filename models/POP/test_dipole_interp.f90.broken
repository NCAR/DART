! DART software - Copyright UCAR. This open source software is provided
! by UCAR, "as is", without charge, subject to all terms of use at
! http://www.image.ucar.edu/DAReS/DART/DART_download
!
! DART $Id$

program test_dipole_interp

! I believe that I have done an initial test of the new interpolation
! in the directory /fs/image/home/jla/DART_FIX_POP_INTERP/DART/models/POP/work
! and that it worked successfully. I needed to copy some extra files from
! my old test directories. I ran test_dipole_interp and selected option 2,
! which interpolates from a regular grid to the dipolex3 (I think). I then
! plotted the output with plot_global_grid.m and plot_globe_dif.m. The
! results looked consistent with the results a obtained with the old interp
! version on coral.
! 
! We should definitely do additional clean-up on this testing case, and
! possibly include a x1 test version, plus add a bit more documentation.
! 
! Confirming my test on coral and then testing on blue_fire seems the
! appropriate next steps.

use model_mod, only : test_interpolation

implicit none
integer :: test_casenum

write(*, *) 'There are 3 pairs of tests'
write(*, *) 'Run test 1 followed by test 2'
write(*, *) '    test 3 followed by test 4'
write(*, *) ' or test 5 followed by test 6'
write(*, *)
write(*, *) 'Test Case Menu'
write(*, *)
write(*, *) '1: regular grid ->  dipole grid'
write(*, *) '2: dipole  grid -> regular grid'
write(*, *) '-------------------------------'
write(*, *) '3: regular grid -> regular grid with same grid as dipole in SH'
write(*, *) '4: regular grid with same grid as dipole in SH -> regular grid'
write(*, *) '-------------------------------'
write(*, *) '5: regular grid with same grid as dipole in SH ->  dipole grid'
write(*, *) '6: dipole  grid -> regular grid with same grid as dipole in SH'
write(*, *)
write(*, *) 'Select integer case'
read(*, *) test_casenum

call test_interpolation(test_casenum)

end program test_dipole_interp

! <next few lines under version control, do not edit>
! $URL$
! $Id$
! $Revision$
! $Date$
