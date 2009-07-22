program test_dipole_interp

use model_mod, only : test_interpolation

implicit none
integer :: test_casenum

write(*, *) 'There are 3 pairs of tests'
write(*, *) 'Run test 1 followed by test 2'
write(*, *) '    test 3 followed by test 4'
write(*, *) ' or test 5 followed by test 6'

write(*, *) 'Test Case Menu'
write(*, *) '1: regular grid -> dipole x3'
write(*, *) '2: dipole x3 -> regular grid'
write(*, *) '3: regular grid -> regular grid with same grid as x3 in SH'
write(*, *) '4: regular grid with same grid as x3 in SH -> regular grid'
write(*, *) '5: regular grid with same grid as x3 in SH -> dipole x3'
write(*, *) '6: dipole x3 -> regular grid with same grid as x3 in SH'
write(*, *)
write(*, *) 'Select integer case'
read(*, *) test_casenum
call test_interpolation(test_casenum)

end program test_dipole_interp
