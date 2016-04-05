      SUBROUTINE GMADDV(A,B,R,IR,AOS)

c  Version 2.0 created by J. Forbes 2/6/91

      COMPLEX A(IR),B(IR),R(IR)

      DO 1 N=1,IR
    1 R(N)=A(N)+AOS*B(N)
      RETURN
      END
