c------------------------------------------------------------
      subroutine NR_SOLVE(a, b, n, info)
c------------------------------------------------------------
c
c   Solves a linear system Ax = b by PLU factorization of A
c  
c   Method:
c
c     1. Factor A into three factors: P, a permutation matrix;
c        L, a  unit lower-triangular matrix; U, an upper-triangular
c        matrix
c
c     2. Solve Ly = P^(-1) b for y using forward-substitution
c
c     3. Solve Ux = y for x using back-substitution
c
c   ref: Conte, S. D. and de Boor, C. (1980),
c        Elementary numerical analysis. McGraw-Hill,
c        New York, NY.
c------------------------------------------------------------
c   arguments:
c   
c     n: (input) integer, dimension of linear system
c
c     a: (input/output) real array dimension nxn
c         on input: matrix of coefficients A
c         on exit:  LU factors of A
c
c     b: (input/output) real array dimension n
c         on input: array b
c         on exit: solution x
c
c     info: (output) integer, error flag
c          0 : successful exit
c         -1 : upper triangular factor has one or more zero
c              diagonal entries 
c
c------------------------------------------------------------

      implicit none

      integer :: n, info
      real :: a(n,n), b(n)
      
      real :: d(n), x(n)
      integer :: iflag, ipivot(n)
      
      call FACTOR(a, n, d, ipivot, iflag)

      if (iflag .eq. 0) then
        info = -1
        return
      endif 

      call SUBST(a, ipivot, b, n, x)

      b = x
      info = 0

      return
      end

c------------------------------------------------------------
      subroutine FACTOR(w1, n, d1, ipivot, iflag)
c------------------------------------------------------------
c
c  source : Conte, S. D. and de Boor, C. (1980).
c    Elementary numerical analysis. McGraw-Hill,
c    New York, NY.
c
c  input: 
c    w: nxn array containg the matrix A to be factored
c    n: order of array A
c
c  output:
c    d: real vector of length N to hold row sizes
c    w: array of size nxn containing LU factorization of P*A
c       for some permutation matrix P specified by ipivot
c    ipivot: integer vector of length n indicating that row
c       ipivot(k) was used to eliminate x(k), k=1,...,n
c    iflag: = 1, if an even number of interchanges was carried out
c           =-1, if an odd number of interchanges was carried out
c           = 0, upper triangular factor has one or more zero diaagonal entries
c
c  if iflag .ne. 0 then linear system A*X = B can be solved for X by a call to SUBST
c 
c------------------------------------------------------------
 
      implicit none

      integer :: n, iflag, ipivot(n)
      real ::  d1(n), w1(n,n)
      
      real:: awikod, colmax
      real:: ratio, rowmax, temp

      integer :: istar, istep, i, j, k
      integer :: k1
      integer :: n1, i1 
      
      iflag=1

c... initialize ipivot, d1

      do 20 i=1,n
        ipivot(i)=i
        rowmax=0.
        do 10 j=1,n
          rowmax=max(rowmax,abs(w1(i,j)))
10      continue
        if(rowmax .eq. 0.) then
          iflag=0
          rowmax=1.e0
        endif
      d1(i)=rowmax
20    continue

      if(n.le.1) return

c... factorization
      n1=n-1
      do 70 k=1,n1
        colmax=abs(w1(k,k))/d1(k)
        istar=k
        k1=k+1
        do 30 i=k1,n
          awikod=abs(w1(i,k))/d1(k)
          if(awikod.gt.colmax) then
            colmax=awikod
            istar=i
          endif
30      continue
        if(colmax .eq. 0.) then
          iflag = 0
        else
          if(istar .gt. k) then
c... make K the p[ovot row by interchanging it with the chosen row istar
            iflag = -iflag
            i = ipivot(istar)
            ipivot(istar) = ipivot(k)
            ipivot(k) = i
            temp = d1(istar)
            d1(istar) = d1(k)
            d1(k) = temp
            do j=1, n
              temp=w1(istar,j)
              w1(istar,j)=w1(k,j)
              w1(k,j)=temp
            enddo
          endif

c... eliminate x(k) from rows k+1, ..., n

          do i=k+1, n
            w1(i,k)=w1(i,k)/w1(k,k)
            ratio=w1(i,k)
            do j=k+1, n
               w1(i,j)=w1(i,j)-ratio*w1(k,j)
            enddo 
          enddo

        endif
70    continue
      if(w1(n,n) .eq. 0.) iflag=0

      return

      end

c------------------------------------------------------------
      subroutine SUBST(w, ipivot, b, n, x)
c------------------------------------------------------------
c
c  source : Conte, S. D. and de Boor, C. (1980).
c    Elementary numerical analysis. McGraw-Hill,
c    New York, NY.
c
c  input: 
c    w, ipivot, n : outputs from FACTOR()
c    b : an n-vector giving the rhs of system to b solved
c
c  output:
c    x : an n-vector satisfying A*X = B
c
c------------------------------------------------------------
 
      implicit none

      integer :: n, ipivot(n)
      real :: b(n), w(n,n), x(n)

      integer :: i, j, ip
      real :: sum

      if(n.le.1) then
        x(1)=b(1)/w(1,1)
        return
      endif

      ip=ipivot(1)
      x(1)=b(ip)

      do i=2,n
        sum=0.

        do j=1,i-1
          sum = w(i,j)*x(j) + sum
	enddo

        ip = ipivot(i)
        x(i) = b(ip) - sum
      enddo

      x(n) = x(n)/w(n,n)

      do i = n-1, 1, -1
        sum=0.
        
        do j=i+1, n
          sum = w(i,j)*x(j) + sum
        enddo

        x(i) = (x(i)-sum)/w(i,i)

      enddo

      return

      end

