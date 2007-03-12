program kurtosis

! Data Assimilation Research Testbed -- DART
! Copyright 2004-2007, Data Assimilation Research Section
! University Corporation for Atmospheric Research
! Licensed under the GPL -- www.gpl.org/licenses/gpl.html
!
! <next few lines under version control, do not edit>
! $URL$
! $Id$
! $Revision$
! $Date$

implicit none

! Look at adjusting ensemble to get a given variance and kurtosis (0 mean for now)

integer, parameter :: n = 10
real :: ens(n), var, kur, power
integer :: i, j, k
real, parameter :: target_var = 1.0, target_kurtosis = 2.0

! Can just make initial distribution uniform over some interval
!do i = 1, n
!   ens(i) = -1.0 + 2.0 * (i - 1.0) / (n - 1.0)
!   write(*, *) i, ens(i)
!end do

! Can also make initial distribution evenly spaced in density
do i = 1, n/2
   call inverse_gaussian(1.0 * i / (n + 1.0), ens(i))
   ens(n + 1 - i) = -1.0 * ens(i)
end do
if(n / 2 * 2 /= n) ens(n/2 + 1) = 0.0

do i = 1, n
   write(*, *) i, ens(i)
end do

! Initial var and kurtosis
call var_kur(ens, n, var, kur)
write(*, *) 'initial var, kurtosis', var, kur


!  Fish for an appropriate power to spread out the outside
! Start out assuming that kurtosis is too small, but this isn't general???
! WANT TO DO A SEARCH EVENTUALLY

! Can do this at the start of the ensemble filter just once and cache it.
do j = 1, 1000
   power = 1.0 + 0.00001 * j
   do i = 1, n / 2 
      k = n + 1 - i
      ens(k) =  ens(k) ** power
      ens(i) = -1.0 * ens(k)
   end do
   call var_kur(ens, n, var, kur)
   write(*, *) j, power, var, kur
   if(kur > target_kurtosis) goto 10
end do

10 continue
! Now normalize the variance
ens = ens * sqrt(target_var / var)
call var_kur(ens, n, var, kur)
write(*, *) 'Final stats ', var, kur
write(*, *) 'final ensemble ', ens



contains

!-------------------------------------------------------
subroutine var_kur(x, n, var, kur)

integer, intent(in) :: n
real, intent(in) :: x(n)
real, intent(out) :: var, kur

var = sum(x(:)**2) / (n - 1.0)
kur = (sum(x(:)**4) / var**2) / n

end subroutine var_kur

!-------------------------------------------------------
subroutine normalize(x, n)

integer, intent(in) :: n
real, intent(inout) :: x(n)

real :: var, kur

call var_kur(x, n, var, kur)
x = x * sqrt(1.0 / var)

end subroutine normalize

!-------------------------------------------------------
subroutine inverse_gaussian(p, x)

real, intent(in) :: p
real, intent(out) :: x

real :: t, num, denom

t = sqrt(log(1.0 / p**2))                                                                
num = 2.515517 + .802853*t + .010328*t**2                                                
denom = 1.0 + 1.432788*t + .189269*t**2 + .001308*t**3                                   
x = -1.0 * (t - num / denom)

end subroutine inverse_gaussian


end program kurtosis

