
      SUBROUTINE HOFFI(H,IZ,IN,IJ,C,T,SS)
c  Modifications for 2-hemis; 2-deg grid capabilities
c       Array dimensions set to "upper limit value" (91) --M.Hagan (3/22/91)
      real H(3,5,45),C(91),T(45),S(45),SS(91),S1(3),S2(3)
      DO 1 I=1,45
1     S(I)=H(IZ,IN,I)
      DO 11 I=1,IJ
      COLAT=C(I)
      CALL ATSM(COLAT,T,S,45,1,S1,S2,3)
      CALL ALI(COLAT,S1,S2,FD,3,1.0E-04,IER)
11    SS(I)=FD
      RETURN
      END
