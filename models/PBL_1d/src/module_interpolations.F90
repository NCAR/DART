MODULE module_interpolations
!
! DART $Id$
!

public spline, linear, seval
private
CONTAINS

  SUBROUTINE spline (n, x, y, b, c, d)
    INTEGER :: n
    REAL :: x(n), y(n), b(n), c(n), d(n)
! 
!   the coefficients b(i), c(i), and d(i), i=1,2,...,n are computed
!  for a cubic interpolating spline
!
!    s(x) = y(i) + b(i)*(x-x(i)) + c(i)*(x-x(i))**2 + d(i)*(x-x(i))**3
!
!    for  x(i) .le. x .le. x(i+1)
!
!  input..
!
!    n = the number of data points or knots (n.ge.2)
!    x = the abscissas of the knots in strictly increasing order
!    y = the ordinates of the knots
!
!  output..
!
!    b, c, d  = arrays of spline coefficients as defined above.
!
!  using  p  to denote differentiation,
!
!    y(i) = s(x(i))
!    b(i) = sp(x(i))
!    c(i) = spp(x(i))/2
!    d(i) = sppp(x(i))/6  (derivative from the right)
!
!  the accompanying function subprogram  seval  can be used
!  to evaluate the spline.
!
!
    INTEGER :: nm1, ib, i
    REAL :: t
!
    nm1 = n-1
    IF ( n .LT. 2 ) RETURN
    IF ( n .LT. 3 ) go to 50
!
!  set up tridiagonal system
!
!  b = diagonal, d = offdiagonal, c = right hand side.
!
    d(1) = x(2) - x(1)
    c(2) = (y(2) - y(1))/d(1)
    DO  i = 2, nm1
       d(i) = x(i+1) - x(i)
       b(i) = 2.*(d(i-1) + d(i))
       c(i+1) = (y(i+1) - y(i))/d(i)
       c(i) = c(i+1) - c(i)
    ENDDO
!
!  end conditions.  third derivatives at  x(1)  and  x(n)
!  obtained from divided differences
!
    b(1) = -d(1)
    b(n) = -d(n-1)
    c(1) = 0.
    c(n) = 0.
    IF ( n .EQ. 3 ) go to 15
    c(1) = c(3)/(x(4)-x(2)) - c(2)/(x(3)-x(1))
    c(n) = c(n-1)/(x(n)-x(n-2)) - c(n-2)/(x(n-1)-x(n-3))
    c(1) = c(1)*d(1)**2/(x(4)-x(1))
    c(n) = -c(n)*d(n-1)**2/(x(n)-x(n-3))
!
!  forward elimination
!
15  DO  i = 2, n
       t = d(i-1)/b(i-1)
       b(i) = b(i) - t*d(i-1)
       c(i) = c(i) - t*c(i-1)
    ENDDO
!
!  back substitution
!
    c(n) = c(n)/b(n)
    DO  ib = 1, nm1
       i = n-ib
       c(i) = (c(i) - d(i)*c(i+1))/b(i)
    ENDDO
!
!  c(i) is now the sigma(i) of the text
!
!  compute polynomial coefficients
!
    b(n) = (y(n) - y(nm1))/d(nm1) + d(nm1)*(c(nm1) + 2.*c(n))
    DO  i = 1, nm1
       b(i) = (y(i+1) - y(i))/d(i) - d(i)*(c(i+1) + 2.*c(i))
       d(i) = (c(i+1) - c(i))/d(i)
       c(i) = 3.*c(i)
    ENDDO
    c(n) = 3.*c(n)
    d(n) = d(n-1)
    RETURN
!
50  b(1) = (y(2)-y(1))/(x(2)-x(1))
    c(1) = 0.
    d(1) = 0.
    b(2) = b(1)
    c(2) = 0.
    d(2) = 0.
    RETURN
  END SUBROUTINE spline
  

  FUNCTION seval(n, u, x, y, b, c, d)
    INTEGER :: n
    REAL :: seval
    REAL  :: u, x(n), y(n), b(n), c(n), d(n)
!
!  this subroutine evaluates the cubic spline function
!
!    seval = y(i) + b(i)*(u-x(i)) + c(i)*(u-x(i))**2 + d(i)*(u-x(i))**3
!
!    where  x(i) .lt. u .lt. x(i+1), using horner's rule
!
!  if  u .lt. x(1) then  i = 1  is used.
!  if  u .ge. x(n) then  i = n  is used.
!
!  input..
!
!    n = the number of data points
!    u = the abscissa at which the spline is to be evaluated
!    x,y = the arrays of data abscissas and ordinates
!    b,c,d = arrays of spline coefficients computed by spline
!
!  if  u  is not in the same interval as the previous call, then a
!  binary search is performed to determine the proper interval.
!
    INTEGER :: i, j, k
    REAL :: dx

    i=1
    IF ( i .GE. n ) i = 1
    IF ( u .LT. x(i) ) go to 10
    IF ( u .LE. x(i+1) ) go to 30
!
!  binary search
!
10  i = 1
    j = n+1
20  k = (i+j)/2
    IF ( u .LT. x(k) ) j = k
    IF ( u .GE. x(k) ) i = k
    IF ( j .GT. i+1 ) go to 20
!
!  evaluate spline
!
30  dx = u - x(i)
    seval = y(i) + dx*(b(i) + dx*(c(i) + dx*d(i)))
    RETURN
  END FUNCTION seval

  REAL FUNCTION linear(x,y,n,xi)

!.....performs linear interpolation for function y given at points x
!.....for points xi, result is output as yi

!.....x are in increasing sequence

    IMPLICIT NONE

!    REAL :: linear
    INTEGER :: i,n
    REAL, intent(in), dimension(n) :: x,y
    REAL :: xi
    INTEGER :: j, nx

    real, dimension(n) :: xtmp,ytmp

    xtmp = x
    ytmp = y

    IF ( n > 1 ) THEN
       i = 1
       DO WHILE ( xtmp(i) < xi .and. i < n ) 
         i = i+1
       ENDDO
       nx = i! + 1
    ELSE
       nx = n
    ENDIF

    DO i=2,nx
       IF (xtmp(i) < xtmp(i-1)) THEN
          IF (ABS(xtmp(i)-xtmp(i-1)) < 1.e-4) THEN
             xtmp(i)=xtmp(i-1)+1.e-4
             ytmp(i)=ytmp(i-1)
          ELSE
             PRINT *,'x is not increasing order',xtmp(i-1),xtmp(i)
             STOP
          ENDIF
       ENDIF
    ENDDO
    
    IF (ABS(xi-xtmp(1)) < 1.e-4) THEN
       linear=ytmp(1)
       RETURN
    ENDIF
    
    IF (xi < xtmp(1)) THEN
       PRINT *,'need to extrapolate, -9999. assigned'
       PRINT *, 'xi less than x(1) ',xi,xtmp(1)
       PRINT*
       PRINT *,'NOTE: This will cause the model to crash and is only ',&
               'appropriate for outputting the analyis'
       PRINT*
!       PRINT *,'Stopping in linear'
       linear = -9999.0
       RETURN
!       STOP
    ENDIF
    
    IF (xi > xtmp(nx)) THEN
!       PRINT *,'need to extrapolate 99999. assigned'
       PRINT *, 'xi greater than x(n) ',xi,xtmp(n)
       PRINT *,'Stopping in linear'
       STOP
    ENDIF
    
    j=1
    DO WHILE (xtmp(j) < xi)
       j=j+1
       IF (ABS(xi-xtmp(j)) < 1.e-4) THEN
          linear=ytmp(j)
          RETURN
       ENDIF
    ENDDO

    linear=ytmp(j)-(xtmp(j)-xi)*(ytmp(j)-ytmp(j-1))/(xtmp(j)-xtmp(j-1))
    RETURN
    
  END FUNCTION linear
  
END MODULE module_interpolations
