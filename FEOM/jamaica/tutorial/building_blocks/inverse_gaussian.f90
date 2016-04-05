! Data Assimilation Research Testbed -- DART
! Copyright 2004-2007, Data Assimilation Research Section
! University Corporation for Atmospheric Research
! Licensed under the GPL -- www.gpl.org/licenses/gpl.html

program inverse_gaussian

! <next few lines under version control, do not edit>
! $URL$
! $Id$
! $Revision$
! $Date$

implicit none

real :: p, t, num, denom, x

10 write(*, *) 'input value of to be inverted with gaussian cdf'
read(*, *) p

t = sqrt(log(1.0 / p**2))
num = 2.515517 + .802853*t + .010328*t**2
denom = 1.0 + 1.432788*t + .189269*t**2 + .001308*t**3
x = t - num / denom
write(*, *) x

goto 10

end program inverse_gaussian
