      SUBROUTINE GMADDM(A,B,R,IR,IC,AOS)

c  Version 2.0 created by J. Forbes 2/6/91

      COMPLEX A(IR,IC),B(IR,IC),R(IR,IC)

      DO 1 N=1,IR
      DO 1 M=1,IC
    1 R(N,M)=A(N,M)+AOS*B(N,M)
      RETURN
      END
