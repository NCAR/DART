      SUBROUTINE HTVEC

C     ARRAYS OF ALTITUDE (Z) AND STRETCHED VAR. (X) ARE CREATED. 
C     TABLES LATER USED FOR X TO Z CONVERSION BY INTERPOLATION.

      DIMENSION CS(3),DS(2),ZS(2)

      COMMON/HTVECT/XV(501),ZV(501)

      DATA N1,N2,NZMAX/1,501,501/
      DATA CS(1),CS(2),CS(3),DS(1),ZS(1),DS(2),ZS(2)
     1     /4.,.15,.025,.5,1.5,40.,120./

301   FORMAT(/)

      K=0
      DO 22 K=1,20
   22 ZV(K)=.2*(K-1)
      DO 23 KK=1,481
      M=KK+20
      K=KK+3
   23 ZV(M)=FLOAT(K)
      DZ=(N2-N1)/(NZMAX-1)
      ND=(N2-N1)/(NZMAX-1)
      I=0
      DO 20 JZ=N1,N2,ND
      I=I+1
      Z=ZV(I)
      XV(I)=.5*Z*(CS(1)+CS(3))
      DO 20 K=1,2
      ARG=(Z-ZS(K))/DS(K)
      ARG2=ZS(K)/DS(K)
      FACT1=.5*DS(K)*(CS(K+1)-CS(K))
      IF(ARG.GT.10.) GO TO 2
      FACT2=ALOG(COSH(ARG)/COSH(ARG2))
      GO TO 20
2     FACT2=(Z-2.*ZS(K))/DS(K)
20    XV(I)=XV(I)+FACT1*FACT2

 101  FORMAT(19X,58HTHE VECTORS CONTAINING ALTITUDE (Z) AND STRETCHED VA
     1R. (X))
c     write (10, 101)

c     write (10, 301)
 200  format(5(1x,2hZ=,f4.0,1x,2hX=,f6.2))
      nzmax5=nzmax/5
c     do 2005 j=1,nzmax5
c     i=5*(j-1)+1
c     write (10, 200) zv(i),xv(i),zv(i+1),xv(i+1),zv(i+2),xv(i+2),
c    +   zv(i+3),xv(i+3),zv(i+4),xv(i+4)
c2005 continue

      RETURN
      END
