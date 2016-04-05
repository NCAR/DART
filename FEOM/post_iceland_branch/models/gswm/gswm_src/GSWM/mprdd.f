      SUBROUTINE MPRDD(A,B,R,N,M,L)

c  Version 2.0 created by J. Forbes 2/3/91

      COMPLEX A(N,M),B(N,M),R(N,M)
      DO 1 IL=1,L
      DO 1 IR=1,N
      R(IR,IL)=(0.0,0.0)
      DO 1 IC=1,M
    1 R(IR,IL)=A(IR,IC)*B(IC,IL)+R(IR,IL)
      RETURN
      END
