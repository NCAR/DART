
      SUBROUTINE PRNT(X,F,IR,IJ)

c  Modifications for variable dimension F array         ---M.Hagan (3/22/91)

      COMPLEX R(4),F(IR)

      COMMON /MODEi/NZONAL,MOIS,NSS
      COMMON /MODEr/PERIOD,FREQ,FLUX
C     latvect arrays are in colat from the north pole to the south pole
      COMMON/LATVECT/DTHETA,CLATT(91),XLAT(91),SNCLAT(91),CSCLAT(91),
     1    TNCLAT(91),CTNCLAT(91),SIN2I(91),GMLAT(91),DIP(91)

      common/postproc/idecomp,idenpress

100   FORMAT(1X,2HX=,F7.4,1X,2HZ=,F7.3)
200   FORMAT(1X,4HLAT=,F5.1,3X,2HU=,e8.3,1H/,F5.1,2X,2HV=,e8.3,1H/,F5.1,
     1                      2X,2HW=,e8.3,1H/,F5.1,2X,2HT=,e8.3,1H/,F5.1)
400   FORMAT(/)

      CALL CONV(X,Z)
      EXF=EXP(X/2.)
      RADEG = 180./ACOS(-1.)

C     PRINT 100,X,Z
      WRITE(7,100) X,Z

    3 L=0
      DO 1 I=1,IJ
      DO 2 J=1,4
      L=L+1
    2 R(J)=F(L)*EXF

         CALL AMPHZ(R,4)

      VLAT = XLAT(I)*RADEG
C     PRINT 200,VLAT,R
      WRITE(7,200) VLAT,R
    1 CONTINUE
      write(7,400)

      RETURN
      END
