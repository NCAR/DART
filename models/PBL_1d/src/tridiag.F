SUBROUTINE tridiag(n,a,b,c,f)

  !! conventions used:
  !! input is lower case
  !! output is upper case

  !! to solve system of linear eqs on tridiagonal matrix n times n
  !! after Peaceman and Rachford, 1955
  !! a,b,c,d - are vectors of order n 
  !! a,b,c - are coefficients on the LHS
  !! d - is initially RHS on the output becomes a solution vector

  IMPLICIT NONE

  INTEGER :: n
  REAL, DIMENSION(n) :: a,b,c,f
  
  INTEGER :: i
  REAL :: p
  REAL, DIMENSION(n) :: q

  c(n)=0.
  q(1)=-c(1)/b(1)
  f(1)=f(1)/b(1)
  
  DO i=2,n
     p=1./(b(i)+a(i)*q(i-1))
     q(i)=-c(i)*p
     f(i)=(f(i)-a(i)*f(i-1))*p
  ENDDO

  DO i=n-1,1,-1
     f(i)=f(i)+q(i)*f(i+1)
  ENDDO
  
END SUBROUTINE tridiag
