      SUBROUTINE ABCRN(A,B,C,K,N,IR,IC,DX)
C     CONVERTS DIFFERENTIAL FORM OF EQ TO FINITE DIFFERENCE FORM

c  Version 2.0 created by J. Forbes 2/3/91

      COMPLEX A(IR,IC),B(IR,IC),C(IR,IC)
      F=1./DX
      IF(K-1) 1,1,2
2     IF(K-N) 3,8,8
1     DO 10 I=1,IR
      DO 10 J=1,IC
      A(I,J)=A(I,J)-F*C(I,J)
10    B(I,J)=F*C(I,J)
      GO TO 6
3     F2=F*F
      F3=F/2.
      F4=-2.*F2
      DO 20 I=1,IR
      DO 20 J=1,IC
      B(I,J)=B(I,J)+F4*C(I,J)
      C(I,J)=F2*C(I,J)+F3*A(I,J)
20    A(I,J)=C(I,J)-F*A(I,J)
      GO TO 6
8     DO 30 I=1,IR
      DO 30 J=1,IC
      B(I,J)=A(I,J)+F*C(I,J)
30    A(I,J)=-F*C(I,J)
6     RETURN
      END
