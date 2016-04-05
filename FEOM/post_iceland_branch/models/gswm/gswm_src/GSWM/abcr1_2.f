      SUBROUTINE ABCR1(A,B,C,D,K,N,X,IR,IC,IJ)
C
C   UBC-XTOP correction: missing H added to the denominator of A(IRO,ICO)
C	Note error in Forbes (1982) etc.; see derivation   --M. Hagan (8/25/98)

C This version of abcr1_2 for GSWM98 uses HDIFF modifications found in
c Forbes ~/jforbes/GSWM/PW  --jkh 7/15/98.  Marked by CForbes

c Jeff's spectal16 version w/ call heat16 replaced by call heat3 and
C       call bfunc commented				--M. Hagan (12/15/92)

c This subroutine generates the coefficient arrays for the
c Lindzen-Kuo algorithm
c
c Modified for SR h_diff option; ihdiff switch in driver --M. Hagan (7/15/99)
c
C Corrected error in Polar BC                            --M. Hagan (2/15/02)
C
C Relaced call for forcing with 7yr annual average latent heating with
C   call monthly averaged rates when doing migrating calculations
C                                                        --M. Hagan (3/22/02)
C Adjust latent heating to account for stratiform precipitation 
C   (25% reduction)                                      --M. Hagan (3/26/02)
C
C Stripped down heating and background options for standard migrating
C    calculations for Tomoko Matsuo                      --M. Hagan (3/17/04)

      integer latgrad,polebc
      integer idecomp, idenpress
      integer mx,ny,iz
      COMPLEX A(IR,IC),B(IR,IC),C(IR,IC),D(IR)
      complex XX(4,4),XXX(4),XI
      COMPLEX XMOM(91),YMOM(91)
      COMPLEX E0,E1P,E1M,E2P,E2M
      COMPLEX XIFNST,XINSG,SIG,SIGT,SIGZ,SIGSTR
c     COMPLEX BCA1,BCAIJ,BCB1(4,4),BCBIJ(4,4)
      COMPLEX BCA1(4,4),BCAIJ(4,4),BCB1(4,4),BCBIJ(4,4)

      REAL MU0,K0,MW,H,hdiff

      real sigmat

      real HJ(91)

      COMMON/ATM1/P1(91),P2(91),P3(91),P4(91),P5(91),P6(91),
     +P7(91),P8(91),P9(91),P10(91),P11(91),P12(91),P13(91),P14(91)
      COMMON/ATM2/DP1(91),DP3(91),DP8(91),DP9(91),D2P8(91),P15(91)
      common/interp1/zpt(37,101,3),zpu(37,101,3),zpr(37,101,3),
     +	zpxh(37,101,3),dO2(37,101,3)
      common/interp2/xy(37),y(101),zt(37,101),zu(37,101),zr(37,101),
     +  zxh(37,101),zph(37,101),o2(37,101),sigma
      
      COMMON/LATVECT/DTHETA,CLATT(91),XLAT(91),SNCLAT(91),CSCLAT(91),
     1   TNCLAT(91),CTNCLAT(91),SIN2I(91),GMLAT(91),DIP(91)

      COMMON /MODEi/NZONAL,MOIS,NSS
      COMMON /MODEr/PERIOD,FREQ,FLUX

	common /backatm/latgrad,polebc

      common/postproc/idecomp,idenpress

      DATA OMEGA,W,XI/.72722E-04,6378165.,(0.0,1.0)/

      FN=FLOAT(NZONAL)

C ****** DEFINES BACKGROUND ATMOSPHERE PARAMETERS AT CURRENT X VALUE *****

      CALL CONV(X,Z)

      CALL ATMOS5(Z,IJ,G,GAM,GAM1)
c
c Make a latitude index that represents the equator

      IEQ=(IJ+1)/2
      
C  Initialize "forcing" matrices
       DO 11 I=1,IJ
       HJ(I)=0.0
 11    CONTINUE

c Specify migrating tidal forcing 

c          print *,'Calling SR heat92 from SR abcr1_2 at Z=',z
          CALL heat92(nzonal,mois,z,HJ,ij)

C  Define eddy viscosity
C  Initialize values at all altitudes JKH 8/18/98

       eddyv=10.0
       eddyk=13.6
       deddyv=0.0
       deddyk=0.0

C  Define Newtonian Cooling
      AZ=0.0

      CALL COOL(Z,AZ)
C
C Set minimum threshold for AZ of 1e-14 (value at 110km) jkh 8/10/98
C     because its exponent can otherwise become large...

      if(az.lt.1e-14) az=0.0

C  Define coeff. to incorporate stretching in vertical coordinate
      CALL VARX(Z,F1,F2)
      F12=F1*F1

C  Initialize coefficient matrices
      DO 1 I=1,IR
      D(I)=(0.0,0.0)
      DO 1 J=1,IC
      A(I,J)=(0.0,0.0)
      B(I,J)=(0.0,0.0)
1     C(I,J)=(0.0,0.0)

C     LOWER BOUNDARY CONDITION
      IF(K-N) 3,2,2

2     CONTINUE
      IF(X.LT.0.005) GO TO 98

C     BOUNDARY CONDITION FOR U=V=W=T=0.0 AT XBOT NEQ. 0.0

       DO 97 I=1,IR
       D(I)=(0.0,0.0)
       C(I,I)=(0.0,0.0)
97     A(I,I)=CMPLX(1.0,0.0)
       GO TO 6

C     BOUNDARY CONDITION FOR XBOT = 0.0 (SURFACE)

 98    Continue

       DO 95 I=1,IR
       D(I)=(0.0,0.0)
95     C(I,I)=(1.0,0.0)
       
       DO 96 I=1,IJ

C  Define eddy viscosity

       colat=clatt(i)

          call eddy_gs(0.0,colat,eddyv,eddyk,deddyv,deddyk)

       cs=.0085/(f1*eddyk)-.5
       csp=.00425/(f1*eddyv)-.5

       I1=(I-1)*4+1
       I2=I1+1
       I3=I1+2
       I4=I1+3
       A(I1,I1)=CMPLX(-CSP,0.0)
       A(I2,I2)=CMPLX(-CSP,0.0)
       A(I4,I4)=CMPLX(-CS,0.0)
       A(I3,I3)=(1.0,0.0)
       C(I3,I3)=(0.0,0.0)
96     CONTINUE

       GO TO 6
C
    3 IF(K-1) 4,4,5
C
C     BOUNDARY CONDITION AT XTOP
C
4     DO 8 I=1,IR
      D(I)=(0.0,0.0)
      C(I,I)=(1.0,0.0)
    8 A(I,I)=(.5,0.0)
      DO 9 I=1,IJ
      IRO=(I-1)*4+3
      ICO=(I-1)*4+4
      R=P14(I)
      H=P5(I)

 9    A(IRO,ICO)=-XI*FREQ*R/(F1*G*H)

      GO TO 6

C     COMPUTE COEFFICIENTS OF DIFFERENTIAL EQUATIONS

      
5     DO 66 I=1,IJ

      COLAT=CLATT(I)
      VLAT=XLAT(I)
      XH=P6(I)
      T=P3(I)
      H=P5(I)
      U=P1(I)
      DUDT=DP1(I)
      DUDZ=P2(I)
c *******************************
      RHO=P8(I)
      DRHODT=DP8(I)
      DRHODZ=P9(I)
      D2DZDT=DP9(I)
      D2DTDT=D2P8(I)
      D2DZDZ=P10(I)
      DTDZ=P4(I)
      DTDT=DP3(I)
      D2TDZ=P11(I)
C ************NO Latitudinal Gradients**
      if(latgrad.eq.1)then
         DRHODT=0.0
         D2DZDT=0.0
         D2DTDT=0.0
         DTDT=0.0
      endif

C  Latitude Invariant Dissipation
cc Above Thermosphere, MU0 molecular dissipation; K0 thermal conductivity;
cc Forbes/Gilette

      if(dissinvar.eq.0)then
         MU0=P12(I)
         K0=P13(I)
      elseif(dissinvar.eq.1)then  !equatorial values, 0 blows up
         MU0=p12(ieq)
         K0=p13(ieq)
      endif
      R=P14(I)
      RDR=P15(I)
      MW=P7(I)
      EXF=MU0/RHO

C  Get horizontal diffusion coefficient

      CALL H_DIFF(Z,VLAT,HDIFF)

C  Define eddy viscosity

         call eddy_gs(z,colat,eddyv,eddyk,deddyv,deddyk)

c     pi=acos(-1.)
c     ilt=vlat*180./pi
c     iclt=colat*180./pi
c767	format(1x,f7.2,2i4,f6.1,2e12.3,f6.1,2e12.3)
c     if(z.gt.50..and.z.lt.60.) then
c       print 767,z,ilt,iclt,u,dudt,dudz,t,dtdt,dtdz
c       endif

C  Get ion drag parameters
      GMLT=GMLAT(I)
      XDIP=DIP(I)
      E0=0.0
      E1P=0.0
      E1M=0.0
      E2P=0.0
      E2M=0.0

         CALL IONOS3(Z,VLAT,COLAT,GMLT,XDIP,MW,RHO,E0,E1P,E1M,E2P,E2M)

C  Set effective Rayleigh friction to zero,
C		unless calculating diurnal tide:
      tnc=0.0
      rf=0.0
      rfu=0.0
      rfv=0.0
C  Get additional Newtonian Cooling/Rayleigh Friction
C				damping terms for diurnal:
c Hard-wired GW stress for diurnal migrating tide:

      if (period.eq.1.) CALL FRIC(Z,VLAT,TNC,RF,RFU,RFV) 

C  Define various auxiliary coefficients 
      RT=G*H
      GH=G*H
      RTF1W=RT*F1/W
      GAM1T=GAM1*T
      SLAT=SNCLAT(I)
      STRE=SLAT*W

C  The following changes in definition of SIG from previous versions of this
C  program (I.E., Forbes, JGR, Part I, 1982, or AFGL 1982 Rept, Part I) is
C  introduced to generalize the equations to apply for other than migrating
C  tides.  Now, both zonal wave number and frequency must be specified in 
C  the main program.  In order to accomodate S=0, in contrast to older
C  versions, replace SIG by SIG/FN, except in SIGSTR.  An FN must now appear
C  in the numberator of every wxpression in which SIGSTR appears in the
C  denominator. ---- J. Forbes, 6/30/98.

      SIG=(0.0,0.0)
      SIG=FREQ+U*FN/STRE

C  A complex part to the frequency is introduced to avoid singular points in
C  the numerical solution:
c     SIG=SIG-.02*FREQ*XI/FN

      XINSG=XI/(SIG)
      SIGSTR=SIG*STRE
      XIFNST=-XI*FN/STRE
      CLAT=CSCLAT(I)
      TCOMEGA=2.*OMEGA*CLAT
      TLAT=CTNCLAT(I)
      SIGZ=DUDZ*FN/(STRE*SIG)
      SIGT=(DUDT-CLAT*U/SLAT)*FN/(STRE*SIG)

C** D2DX2 TERMS 

      DO 40 J4=1,4
      DO 40 I4=1,4
   40 XX(I4,J4)=(0.0,0.0)

      XX(1,1)=(EXF+EDDYV)*F12
      XX(2,2)=(EXF+EDDYV)*F12
      XX(3,3)=XINSG*RT*F12*H
      XX(4,4)= ((GAM1/R)*(K0/RHO)+EDDYK)*F12
      DO 41 J4=1,4
      DO 41 I4=1,4
      IRO=(I-1)*4+J4
      ICO=(I-1)*4+I4
   41 C(IRO,ICO)=XX(J4,I4)

C** DDX TERMS 

      DO 50 I4=1,4
      DO 50 J4=1,4
   50 XX(I4,J4)=(0.0,0.0)

      XX(1,1)=EXF*(F12+.66*DTDZ*F1-F2)
     1      +EDDYV*(F12+(DEDDYV+DRHODZ)*F1-F2)
      XX(2,2)=XX(1,1)
      XX(3,3)=XINSG*GH*H*(F12-F2+(DRHODZ-SIGZ)*F1)
      XX(3,4)=R*F1*H
      XX(4,3)=-GAM1T*F1
      XX(4,4)=(GAM1/R)*(K0/RHO)*(F12-2.*F2+RDR*F1+1.33*DTDZ*F1)
      XX(1,3)=RT*F1*FN/SIGSTR
      XX(3,1)=-XX(1,3)*H
      XX(3,2)=XINSG*RTF1W*(DRHODT+TLAT)*H
      XX(2,3)=SIGT*XINSG*RTF1W
      DO 51 I4=1,4
      DO 51 J4=1,4
      IRO=(I-1)*4+J4
      ICO=(I-1)*4+I4
   51 A(IRO,ICO)=XX(J4,I4)

C** D2DXDT TERMS

      XX(2,3)=.5*XINSG*RTF1W/DTHETA
      XX(3,2)=-XX(2,3)*H
      IRR1=(I-1)*4+2
      IRR2=(I-1)*4+3
      ICO1=(I-2)*4+3
      ICO2=ICO1-1
      IF(I-1)63,63,61
   61 A(IRR1,ICO1)=XX(2,3)
      A(IRR2,ICO2)=XX(3,2)
      IF(I-IJ) 63,64,64
   63 ICO1=ICO1+8
      ICO2=ICO2+8
      A(IRR2,ICO2)=-XX(3,2)
      A(IRR1,ICO1)=-XX(2,3)
   64 CONTINUE

C  Set up terms to establish lateral boundary condition for D2DXDT term
C  for v, later on.

      IF(I.GT.1.AND.I.LT.IJ) GO TO 1001

cForbes
      do l=1,4
         do m=1,4
            bca1(l,m)=(0.0,0.0)
            bcaij(l,m)=(0.0,0.0)
         end do
      end do
      
      IF(I.EQ.1) THEN
         BCA1(2,3)=XX(2,3)
         BCA1(3,2)=XX(3,2)
      ELSE
         BCAIJ(2,3)=-XX(2,3)
         bcaij(3,2)=-xx(3,2)
      ENDIF
CForbes

 1001 CONTINUE
C____________________________________________________________________________
C OLD CODE

C  Set up terms to establish lateral boundary condition for D2DXDT term
C  for v, later on.

c	IF(I.GT.1.AND.I.LT.IJ) GO TO 1001

c	IF(I.EQ.1) THEN
c	BCA1=-XX(3,2)
c	ELSE
c	BCAIJ=XX(3,2)
c	ENDIF
c1001	CONTINUE
C OLD CODE
C____________________________________________________________________________

C** D2DT2 TERM

      DO 71 J4=1,4
      DO 71 I4=1,4
   71 XX(J4,I4)=(0.0,0.0)

CForbes
c Define coeff of d2dt2 terms in eqs.
c
      denom=w*w*dtheta*dtheta
      xxx(1)=2.*hdiff/denom
      xxx(2)=xxx(1)-2.*xinsg*rt/denom
      xxx(3)=(0.0,0.0)
      xxx(4)=xxx(1)
c      XXX=-2.*XINSG*RT/(W*W*DTHETA*DTHETA)
CForbes

C** DDT TERMS

      XX(1,2)=RT*FN/(SIGSTR*W)
      XX(2,1)=XX(1,2)
      XX(2,2)=-(DRHODT-SIGT+TLAT)*XINSG*RT/(W*W)
      XX(2,3)=-(XINSG*RT/W)*(.5*F1+DRHODZ)
      XX(2,4)=-R/W
      XX(3,2)=(.5*F1-SIGZ)*XINSG*GH*H/W
      XX(4,2)=-GAM1T/W
      
      DO 79 I4=1,4
      DO 79 J4=1,4
      IRO=(I-1)*4+J4
      ICO1=(I-2)*4+I4
      ICO2=(I)*4+I4
      IF(I-IJ) 76,72,72
   76 B(IRO,ICO2)=.5*XX(J4,I4)/DTHETA
C
C
      IF(I-1) 79,79,72
   72 B(IRO,ICO1)=-.5*XX(J4,I4)/DTHETA
   79 CONTINUE
C
C
CForbes
        IRO1=(I-1)*4+1
        IRO2=IRO1+1
        IRO3=IRO1+2
        IRO4=IRO1+3
        IF(I.LT.IJ) THEN
        B(IRO1,IRO1+4)=B(IRO1,IRO1+4)+.5*XXX(1)
        B(IRO2,IRO2+4)=B(IRO2,IRO2+4)+.5*XXX(2)
        B(IRO3,IRO3+4)=B(IRO3,IRO3+4)+.5*XXX(3)
        B(IRO4,IRO4+4)=B(IRO4,IRO4+4)+.5*XXX(4)
        ENDIF
        IF(I.GT.1) THEN
        B(IRO1,IRO1-4)=B(IRO1,IRO1-4)+.5*XXX(1)
        B(IRO2,IRO2-4)=B(IRO2,IRO2-4)+.5*XXX(2)
        B(IRO3,IRO3-4)=B(IRO3,IRO3-4)+.5*XXX(3)
        B(IRO4,IRO4-4)=B(IRO4,IRO4-4)+.5*XXX(4)
        ENDIF
CForbes

C____________________________________________________________________________
C
C  Set up terms to establish lateral boundary condition.....theta-derivatives
C  of u,v = zero.  Without doing anything, the boundary conditions
C  w=T= zero are implicitly included.
C
	IF(I.GT.1.AND.I.LT.IJ) GO TO 1002
CForbes
	DO L=1,4
	DO M=1,4
	BCB1(M,L)=(0.0,0.0)
	BCBIJ(M,L)=(0.0,0.0)
	END DO
	END DO

	IF(I.EQ.1) THEN
           do l=1,4
              DO M=1,4
                 BCB1(l,m)=-.5*XX(l,m)/DTHETA
              END DO
           end do
           do l=1,4
              bcb1(l,l)=bcb1(l,l)+0.5*xxx(l)
           end do
	ELSE
           do l=1,4
              DO M=1,4
                 BCBIJ(l,m)=.5*XX(l,m)/DTHETA
              END DO
           end do
           do l=1,4
cMHtest       BCBIJ(l,m)=BCBIJ(l,m)+.5*XXX(l)
              BCBIJ(l,l)=BCBIJ(l,l)+.5*XXX(l)
           end do
	ENDIF
CForbes
1002	CONTINUE
C____________________________________________________________________________


C** CONSTANT TERMS

      XX(3,1)=-(.5*F1-SIGZ)*GH*H*FN/SIGSTR
      XX(1,1)=RT*XI*FN*FN/(SIG*STRE*STRE)-E0-RF-RFU-XI*SIG
     1        +EXF*(.25*F12-.5*F2+.33*DTDZ*F1)
     2        +EDDYV*(F12-2.*F2+2.*(DEDDYV+DRHODZ)*F1)/4.
      XX(1,2)=(RT*FN/(SIG*STRE*W))*(DRHODT+TLAT)-(CLAT/STRE)*U-
     1(1./W)*DUDT-TCOMEGA
      XX(1,3)=(GH*FN/SIGSTR)*(.5*F1+DRHODZ)-DUDZ
      XX(1,4)=XIFNST*R
      XX(2,1)=TCOMEGA+2.*(CLAT/STRE)*U-
     +     RT*(CLAT+SIGT*SLAT)*FN/(SIGSTR*STRE)
      XX(2,3)=-XINSG*(RT/W)*(-.5*F1*SIGT+D2DZDT-(DRHODT+SIGT)*DRHODZ)
      XX(2,4)=-(R/W)*DRHODT
      XX(3,2)=(XINSG*GH*H/W)*((TLAT+DRHODT)*.5*F1-TLAT*SIGZ+D2DZDT-
     1   (DRHODZ+SIGZ)*DRHODT)
      XX(3,3)=XINSG*GH*H*(.25*F12-.5*F2+.5*(DRHODZ-SIGZ)*F1+D2DZDZ-
     1   DRHODZ*(DRHODZ-SIGZ))
      XX(3,4)=R*H*(.5*F1+DRHODZ+RDR)
      XX(4,1)=GAM1T*XIFNST
      XX(4,2)=-GAM1T*CLAT/STRE-T*DTDT/W
      XX(4,3)=-GAM1T*.5*F1-T*DTDZ
      XX(4,4)=-XI*SIG-AZ-TNC
     1     +(GAM1/R)*(K0/RHO)*(.25*F12-.5*F2
     2     +.5*F1*(1.33*DTDZ+RDR-F2/F1)+.66*D2TDZ
     3     -.22*DTDZ*DTDZ+.66*RDR*DTDZ-F2*.5*DTDZ/F1)
     4     +EDDYK*(F12-2.*F2+2.*(DEDDYK+DRHODZ)*F1)/4.
CForbes
      XX(2,2)=(XINSG*RT/(W*W))*(-D2DTDT+(DRHODT+SIGT)*DRHODT+
     1     CLAT*SIGT/SLAT+
     1     1./(SLAT*SLAT))-E0*SIN2I(I)-XI*SIG-RF-RFV
     2     +EXF*(.25*F12-.5*F2+.33*DTDZ*F1)
     3     +EDDYV*(F12-2.*F2+2.*(DEDDYV+DRHODZ)*F1)/4.

C Take into account d2dt2 terms
      do l=1,4
         xx(l,l)=xx(l,l)-xxx(l)
      end do
CForbes

      DO 81 I4=1,4
      DO 81 J4=1,4
      IRO=(I-1)*4+J4
      ICO=(I-1)*4+I4
   81 B(IRO,ICO)=XX(J4,I4)

C__________________________________________________________________________

C  The following statements implement a Taylor series expansion near the
C  polar boundaries.

CForbes
c
       if(polebc.eq.1)then
 
         if((i-ij).ge.0)then
           IRR1=(IJ-1)*4+2
           IRR2=IRR1+1
           ICO1=(IJ-2)*4+3
           ICO2=ICO1-1
           A(IRR1,ICO1)=2.*A(IRR1,ICO1)
           A(IRR2,ICO2)=2.*A(IRR2,ICO2)
           ICA1=ICO1+4
           ICA2=ICO2+4
           A(IRR1,ICA1)=A(IRR1,ICA1)-A(IRR1,ICO1)
           A(IRR2,ICA2)=A(IRR2,ICA2)-A(IRR2,ICO2)
 
           DO 92 I4=1,4
             DO J4=1,4
               IRO=(I-1)*4+J4
               ICO1=(I-2)*4+I4
               ICO2=(I-1)*4+I4
               B(IRO,ICO1)=2.*B(IRO,ICO1)
               B(IRO,ICO2)=B(IRO,ICO2)-B(IRO,ICO1)
             END DO
  92       CONTINUE

        IRO1=(I-1)*4+1
        IRO2=IRO1+1
        IRO3=IRO1+2
        IRO4=IRO1+3
        B(IRO1,IRO1-4)=B(IRO1,IRO1-4)-XXX(1) 
        B(IRO2,IRO2-4)=B(IRO2,IRO2-4)-XXX(2) 
        B(IRO3,IRO3-4)=B(IRO3,IRO3-4)-XXX(3) 
        B(IRO4,IRO4-4)=B(IRO4,IRO4-4)-XXX(4) 
        B(IRO1,IRO1)=B(IRO1,IRO1)+XXX(1)*2.
        B(IRO2,IRO2)=B(IRO2,IRO2)+XXX(2)*2.
        B(IRO3,IRO3)=B(IRO3,IRO3)+XXX(3)*2.
        B(IRO4,IRO4)=B(IRO4,IRO4)+XXX(4)*2.
 

         endif

         if((i-1).le.0)then
            IRR1=2
            IRR2=3
            ICO1=7
            ICO2=6
            A(IRR1,ICO1)=2.*A(IRR1,ICO1)
            A(IRR2,ICO2)=2.*A(IRR2,ICO2)
            ICA1=3
            ICA2=2
            A(IRR1,ICA1)=A(IRR1,ICA1)-A(IRR1,ICO1)
            A(IRR2,ICA2)=A(IRR2,ICA2)-A(IRR2,ICO2)
           
            DO 192 I4=1,4
               DO J4=1,4
                  IRO=J4
                  ICO1=4+I4
                  ICO2=I4
                  B(IRO,ICO1)=2.*B(IRO,ICO1)
                  B(IRO,ICO2)=B(IRO,ICO2)-B(IRO,ICO1)
               END DO
  192       continue
        IRO1=1
        IRO2=2
        IRO3=3
        IRO4=4
        B(IRO1,IRO1+4)=B(IRO1,IRO1+4)-XXX(1) 
        B(IRO2,IRO2+4)=B(IRO2,IRO2+4)-XXX(2) 
        B(IRO3,IRO3+4)=B(IRO3,IRO3+4)-XXX(3) 
        B(IRO4,IRO4+4)=B(IRO4,IRO4+4)-XXX(4) 
        B(IRO1,IRO1)=B(IRO1,IRO1)+XXX(1)*2.
        B(IRO2,IRO2)=B(IRO2,IRO2)+XXX(2)*2.
        B(IRO3,IRO3)=B(IRO3,IRO3)+XXX(3)*2.
        B(IRO4,IRO4)=B(IRO4,IRO4)+XXX(4)*2.
         endif  !i-1 le 0
CForbes
C__________________________________________________________________________

C  The following implements the lateral boundary condition that the theta-
C  derivatives of u and v = zero at the poles.  Without doing anything, the
C  boundary condition that w=T= zero at the poles is implicitly included(?).
C  If you want to skip this BC implementation, and use a Taylor expansion
C  instead, set polebc flag to 1 instead of 0
CForbes
           endif                !polebc=1
cHack  elseif(polebc.eq.0)then
        IF(I.EQ.1.or.i.eq.ij) THEN
c All sets of boundary conditions
c for polebc=1:
           NT=4
           B1=0.0
           B2=0.0
           BIJM1=0.0
           BIJ=0.0
           if(polebc.eq.2)then !Taylor expansion
              NT=4
              B1=2.0
              B2=-1.0
              BIJM1=-1.0
              BIJ=2.0
           endif                !polebc=2
        endif                   !i=1 or i=ij
c
        IF(I.EQ.1) THEN
           
           DO I4=1,4
              DO J4=1,NT
                 A(I4,J4)=A(I4,J4)+B1*BCA1(I4,J4)
                 A(I4,J4+4)=A(I4,J4+4)+B2*BCA1(I4,J4)
                 B(I4,J4)=B(I4,J4)+B1*BCB1(I4,J4)
                 B(I4,J4+4)=B(I4,J4+4)+B2*BCB1(I4,J4)
              END DO
           END DO
           
        ELSEif(i.eq.ij)then     ! i=ij
           
           DO I4=1,4
              DO J4=1,NT
                 IRO=(IJ-1)*4+I4
                 ICO=(IJ-1)*4+J4
                 A(IRO,ICO)=A(IRO,ICO)+BIJ*BCAIJ(I4,J4)
                 A(IRO,ICO-4)=A(IRO,ICO-4)+BIJM1*BCAIJ(I4,J4)
                 B(IRO,ICO)=B(IRO,ICO)+BIJ*BCBIJ(I4,J4)
                 B(IRO,ICO-4)=B(IRO,ICO-4)+BIJM1*BCBIJ(I4,J4)
              END DO
           END DO
           
        ENDIF                   ! i=1 or ij


cHack  endif    !polebc=0
cForbes
C___________________________________________________________________________
C OLD CODE
c     elseif(polebc.eq.0)then
c       IF(I.EQ.1) THEN
c          A(3,2)=A(3,2)+BCA1

c          DO I4=1,4
c             DO J4=1,4
c                B(J4,I4)=B(J4,I4)+BCB1(J4,I4)
c             END DO
c          END DO

c       ELSEif (i.eq.ij)then
c          IRO=(IJ-1)*4+3
c          ICO=(IJ-1)*4+2
c          A(IRO,ICO)=A(IRO,ICO)+BCAIJ

c          DO I4=1,4
c             DO J4=1,4
c                IRO=(I-1)*4+J4
c                ICO=(I-1)*4+I4
c                B(IRO,ICO)=B(IRO,ICO)+BCBIJ(J4,I4)
c             END DO
c          END DO
c       ELSE  !in do 66 loop i>1 and i<ij
c          continue
c       ENDIF  !i=1, i=ij
c     endif    !polebc=0
C___________________________________________________________________________
C OLD CODE

      IRO1=(I-1)*4+1
      IRO2=IRO1+1
      IRO3=IRO2+1
      IRO4=IRO3+1
      EXR=EXP(-X/2.)

   70 D(IRO1)=(0.0,0.0)
      D(IRO2)=(0.0,0.0)
      D(IRO3)=(0.0,0.0)

C  Define zonal and meridional momentum source terms
      XMOM(I)=(0.0,0.0)
      YMOM(I)=(0.0,0.0)
      D(IRO1)=EXR*XMOM(I)
      D(IRO2)=EXR*YMOM(I)

600   Continue
      
CTEST print out hhj every 18 deg lat and every alt
c      if(mod(i,6).eq.1)then
C         print*,'hhj abcr',z,clatt(i)*57.3,hhj(i),hhj(i)-heuv(i)
c         print*,'Q abcr',z,clatt(i)*57.3,q(nzonal,2)*rho
c      endif
CTEST end
C
C For Insitu HEAT Forcing or new HEAT92 O3 and H2O forcing,
C	when units are J/m3/s (N.B. J/s=W):
C
CTEST1b heat.f forcing
C      D(IRO4)=-Q(NZONAL,2)*GAM1*EXR/(R*RHO)
c      D(IRO4)=-GAM1*EXR*HHJ(I)/(R*RHO)
       D(IRO4)=-GAM1*EXR*HJ(I)/(R*RHO)
CTEST1b end

C
C For old HEAT3 Forcing due to Insolation Absorption of O3 and H2O:
C
c      D(IRO4)=-1.4*FREQ*EXP(-X/2.)*HJ(I)
CTEST print out a lot of physical parameters at 30,90 colatitude JKH8/27/98
c       if(i.eq.ieq.or.i.eq.10)then
C          dum=sqrt(e0*e0)
C          print 995,z,x,colat*57.3,p13(i),p12(i),p14(i),p7(i),p3(i),
C     +         dum,p5(i),p2(i),dp1(i),p8(i),g
c          print 995,z,colat*57.3,k0,mu0,r,mw,t,h,dudz,dudt,
c     +         rho,hhj(i)
c 995      format(1x,13(1x,e9.3))
c       endif
CTEST print out params end
66    CONTINUE

C Forbes
C  INTRODUCE SCALING SO ALL EQS. ARE OF SIMILAR ORDER, AND ORDER ~ 1

	 DO 48 I=1,IR
	 D(I)=D(I)*1.0E+04
	 DO 48 J=1,IC
	 C(I,J)=C(I,J)*1.0E+04
	 A(I,J)=A(I,J)*1.0E+04
   48    B(I,J)=B(I,J)*1.0E+04
	 DO 49 I=1,IJ
	 IR3=(I-1)*4+3
	 D(IR3)=D(IR3)*1.0E-08
	 DO 49 J=1,IC
	 C(IR3,J)=C(IR3,J)*1.0E-08
	 A(IR3,J)=A(IR3,J)*1.0E-08
   49    B(IR3,J)=B(IR3,J)*1.0E-08               

C Forbes

c     CALL AMPHZ(XMOM,IJ)
c     CALL AMPHZ(YMOM,IJ)
c     CALL AMPHZ(DRAG1,IJ)
c     CALL AMPHZ(DRAG2,IJ)
 1099 FORMAT(1X,2HX=,F6.2,2X,2HZ=,F8.4)
C     PRINT 1099,X,Z
C     WRITE (7,1099) X,Z
 1100 FORMAT(1X,3HEW=,5(E10.3,1H/,F4.1),/,4X,5(E10.3,1H/,F4.1),
     +  /,4X,5(E10.3,1H/,F4.1),/,4X,5(E10.3,1H/,F4.1),
     +  /,4X,5(E10.3,1H/,F4.1),/,4X,5(E10.3,1H/,F4.1),
     +  /,4X,5(E10.3,1H/,F4.1),/,4X,5(E10.3,1H/,F4.1))
C     PRINT 1100,(XMOM(I),I=1,IJ)
 1101 FORMAT(1X,3HNS=,5(E10.3,1H/,F4.1),/,4X,5(E10.3,1H/,F4.1),
     +  /,4X,5(E10.3,1H/,F4.1),/,4X,5(E10.3,1H/,F4.1),
     +  /,4X,5(E10.3,1H/,F4.1),/,4X,5(E10.3,1H/,F4.1),
     +  /,4X,5(E10.3,1H/,F4.1),/,4X,5(E10.3,1H/,F4.1))
C     PRINT 1101,(YMOM(I),I=1,IJ)
 1112 FORMAT(1X,3HXH=,7E10.3,/,4X,7E10.3,/,4X,7E10.3,/,4X,7E10.3,
     +  /,4X,7E10.3,/,4X,7E10.3)
C     PRINT 1112,(P6(I),I=1,IJ)
c     WRITE(10,1112) (P6(I),I=1,IJ)
 1113 FORMAT(1X,3H U=,7E10.3,/,4X,7E10.3,/,4X,7E10.3,/,4X,7E10.3,
     +  /,4X,7E10.3,/,4X,7E10.3)
C     PRINT 1113,(P1(I),I=1,IJ)
c     WRITE(10,1113) (P1(I),I=1,IJ)
 1122 FORMAT(1X,5HXlat=,7E10.3,/,6X,7E10.3,/,6X,7E10.3,/,6X,7E10.3,
     +  /,6X,7E10.3,/,6X,7E10.3)
C     PRINT 1122,(XLAT(I),I=1,IJ)
c     WRITE(10,1122) (XLAT(I),I=1,IJ)
 1102 FORMAT(1X,3HQQ=,7E10.3,/,4X,7E10.3,/,4X,7E10.3,/,4X,7E10.3,
     +  /,4X,7E10.3,/,4X,7E10.3)
c     PRINT 1102,(QQ(I),I=1,IJ)
c     WRITE(10,1102) (QQ(I),I=1,IJ)
C1114 FORMAT(1X,3HHJ=,7E10.3,/,4X,7E10.3,/,4X,7E10.3,/,4X,7E10.3,
C    +  /,4X,7E10.3,/,4X,7E10.3)
C     PRINT 1114,(HJ(I),I=1,IJ)
c     WRITE(10,1114) (HJ(I),I=1,IJ)
 1103 FORMAT(1X,3HE0=,5(F10.2,1H/,F4.1),/,4X,5(F10.2,1H/,F4.1),
     +  /,4X,5(F10.2,1H/,F4.1),/,4X,5(F10.2,1H/,F4.1),
     +  /,4X,5(F10.2,1H/,F4.1),/,4X,5(F10.2,1H/,F4.1),
     +  /,4X,5(F10.2,1H/,F4.1),/,4X,5(F10.2,1H/,F4.1))
C     PRINT 1103,(DRAG0(I),I=1,IJ)
 1104 FORMAT(1X,3HE1=,5(F10.2,1H/,F4.1),/,4X,5(F10.2,1H/,F4.1),
     +  /,4X,5(F10.2,1H/,F4.1),/,4X,5(F10.2,1H/,F4.1),
     +  /,4X,5(F10.2,1H/,F4.1),/,4X,5(F10.2,1H/,F4.1),
     +  /,4X,5(F10.2,1H/,F4.1),/,4X,5(F10.2,1H/,F4.1))
c     PRINT 1104,(DRAG1(I),I=1,IJ)
 1105 FORMAT(1X,3HE2=,5(F10.2,1H/,F4.1),/,4X,5(F10.2,1H/,F4.1),
     +  /,4X,5(F10.2,1H/,F4.1),/,4X,5(F10.2,1H/,F4.1),
     +  /,4X,5(F10.2,1H/,F4.1),/,4X,5(F10.2,1H/,F4.1),
     +  /,4X,5(F10.2,1H/,F4.1),/,4X,5(F10.2,1H/,F4.1))
c     PRINT 1105,(DRAG2(I),I=1,IJ)
 7766 Format(1x,'Eddy Diffusion after Garcia & Solomon, 1985')
 7767 Format(1x,2HX=,f5.2,1x,2hZ=,f6.2,1x,
     1    1x,4hLAT=,f5.2,1x,6HEDDYV=,e10.3,1x,7hDEDDYV=,e10.3)

6     CONTINUE

      RETURN
      END
