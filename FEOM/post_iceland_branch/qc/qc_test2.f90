! Data Assimilation Research Testbed -- DART
! Copyright 2004, 2005, Data Assimilation Initiative, University Corporation for Atmospheric Research
! Licensed under the GPL -- www.gpl.org/licenses/gpl.html

program qc_test

! <next five lines automatically updated by CVS, do not edit>
! $Source$ 
! $Revision$ 
! $Date$ 
! $Author$ 
! $Name$ 

! Think of dividing area around distribution into 9 boxes defined by
! the two observations. When the observations have the same sign, boxes
! in lower left and upper right summed (2* either one since they are
! symmetric) is two sided probability that something further out would
! have been picked. When x and y have different signs, the right answer
! is the sum of probabilities in upper left and lower right boxes (again
! 2 * either one by symmetry).

use nag_wrap_mod, only : g01haf_wrap

implicit none

! let CVS fill strings ... DO NOT EDIT ...
character(len=128) :: &
   source   = "$Source$", &
   revision = "$Revision$", &
   revdate  = "$Date$"

integer :: i, j, ifail = 0
double precision :: x, y, rho

double precision :: mx, my, inf, p1, p2, p3, p4, p5, p6, p7, p8, p9
double precision :: r1, r2, r3, c1, c2, c3, total

do i = 1, 100
   write(*, *) 'Input x and y and rho '
   read(*, *) x, y, rho
   x = dabs(x)
   mx = -1.0 * x
   y = dabs(y)
   my = -1.0 * y
   inf = 1e10
   p1 = g01haf_wrap(mx, my, rho, ifail)
   p2 = g01haf_wrap(x, my, rho, ifail) - p1
   p3 = g01haf_wrap(inf, my, rho, ifail) - g01haf_wrap(x, my, rho, ifail)
   p4 = g01haf_wrap(mx, y, rho, ifail) - p1
   p5 = g01haf_wrap(x, y, rho, ifail) - g01haf_wrap(x, my, rho, ifail) - p4
   p6 = g01haf_wrap(inf, y, rho, ifail) - g01haf_wrap(x, y, rho, ifail) - p3
   p7 = g01haf_wrap(mx, inf, rho, ifail) - g01haf_wrap(mx, y, rho, ifail)
   p8 = g01haf_wrap(x, inf, rho, ifail) - g01haf_wrap(x, y, rho, ifail) - p7
   p9 = g01haf_wrap(inf, inf, rho, ifail) - g01haf_wrap(inf, y, rho, ifail) -p8 - p7

   write(*, *) 'p1 ', p1
   write(*, *) 'p2 ', p2
   write(*, *) 'p3 ', p3
   write(*, *) 'p4 ', p4
   write(*, *) 'p5 ', p5
   write(*, *) 'p6 ', p6
   write(*, *) 'p7 ', p7
   write(*, *) 'p8 ', p8
   write(*, *) 'p9 ', p9

   r1 = p1 + p2 + p3
   r2 = p4 + p5 + p6
   r3 = p7 + p8 + p9
   c1 = p1 + p4 + p7
   c2 = p2 + p5 + p8
   c3 = p3 + p6 + p9

   write(*, *) 'row1 ', r1
   write(*, *) 'row2 ', r2
   write(*, *) 'row3 ', r3

   write(*, *) 'col1 ', c1
   write(*, *) 'col2 ', c2
   write(*, *) 'col3 ', c3

   total = r1 + r2 + r3
   write(*, *) 'total ', total




   write(*, *) 'qc is ', g01haf_wrap(x, y, rho, ifail)
end do

end program qc_test
