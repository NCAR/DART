c  Steady-State Spectral Model for Tides,  CRAY YMP version
c  Eight-order system  including full second-order diffusive terms
c  (Forbes, 1982, JGR, Part I; Forbes et al., 1982, AFGL report, Part I).
c  Updated 2/91-6/91, 2-7/98 jkh


      SUBROUTINE MAIN_2 
c*************************************************************************
C  Updated to include call to SR GEO, which replaces GEOPROC and GEOHT to
C  calculate geopotential height and Hough decomposition or pressure and
C  density perturbations.  JHackney 2/98
c  Updated to include call to SR GEOPROC (calculates geopotential height)
C                                               --M. Hagan (1/95)
c*************************************************************************

      PARAMETER( IR=236, IC=236 )

c     NOTE:  When changing IR,IC must also change recl in all OPEN statements
c            recl = {IR x (IC+1)} x 2 x 8
c		  = 16*(IR*IC + IR)

c     3.0 degree grid, 2 hemispheres 
c     parameter(ir=236, ic=236) 
c     recl=894912
c     5.0 degree grid, 2 hemispheres 
c     parameter(ir=140, ic=140) 
c     recl=315840
c     7.5 degree grid, 2 hemispheres 
c     parameter(ir=92, ic=92) 
c     recl=136896

      COMPLEX ALPH(IR,IC),BET(IR),AUX(IR,IC),FX(IR),ZMX(IC)
      COMPLEX A(IR,IC),B(IR,IC),C(IR,IC),R(IR),FNC(IR),FNC1(IR)

      INTEGER IPVT(IC)
      integer igm, ihd, iyp, ibr, nowind, ihm, latgrad, polebc

      integer heatmodel, dblhme11, diurno3, peako3, skipo3, skipir
      integer call_late,o3conc,iboris, eddymodel,kzzinvar, dissinvar
      integer iondrag, vialdiss, gwstress, idecomp, idenpress
      integer f107,mx,ny,iz

      real hjrm(91),hjim(91)
      real prntmax
      complex gciforce(91)
c
        real euvlat(36),reuv(36,36),ieuv(36,36),euvalt(36)
        real sigmat
        real deuv(36,36,3)

c  The following integer common block contains flags for the dissipation
c  models

      common /boris/iboris,eddymodel,kzzinvar,dissinvar,iondrag,
     +     vialdiss,gwstress

c  The following integer common block contains flags for the heat forcing
c       models

      common /heatforce/heatmodel,dblhme11,diurno3,peako3,skipo3,
     +       skipir,call_late,o3conc,heatnm

      common/lateforce/HJRm,HJIm

      COMMON/MODE/NZONAL,PERIOD,FREQ,MOIS,NSS,FLUX
	common/postproc/idecomp,idenpress,mig
      
      COMMON/SEZUN/SWE(3),NSWE
C interp1 and interp2 arrays are defined fron the north to the south poles:
      common/interp1/zpt(37,101,3),zpu(37,101,3),zpr(37,101,3),
     +	zpxh(37,101,3),dO2(37,101,3)
      common/interp2/xx(37),y(101),zt(37,101),zu(37,101),zr(37,101),
     +  zxh(37,101),zph(37,101),o2(37,101),sigma

      COMMON/HOFFS/H11(91),H12(91),H13(91),H1M2(91),H1M1(91)
     1,H1M4(91),H22(91),H23(91),H24(91),H25(91),H26(91)
     + ,H33(91),H34(91),H35(91),H36(91)

      COMMON/HMODE/I11,I12,I13,I1M2,I1M1,I1M4,I22,I23,I24,I25,I26
     +,I33,I34,I35,I36

c  The following common block will contain control flags for the
c  background windfields.  Put it in the relevant subroutines

	common /backatm/igm,ihd,iyp,ibr,nowind,ihm,latgrad,polebc

      COMMON/gci/gciforce
c
c Common block of integers, arrays, floats, for EUV thermosphere heating

      common/euv/mx,ny,iz,euvlat,euvalt,reuv,ieuv,deuv,sigmat,f107
c

      data ias,inp,iztop,iqbo,ici,icount/1,1,0,0,0,1/

      data I11,I12,I13,I1M2,I1M1,I1M4 /0,0,0,0,0,0/
      data I22,I23,I24,I25,I26 /1,0,0,0,0/
      data I33,I34,I35,I36 /0,0,0,0/
C*********************************************************************
C
C     Logic for ** LUNAR TIDE ** is included in MAIN, ABCR1, and PRNT3, but
C     is commented out in this version

C     COMMON/MOON/OM(91),DOMDT(91),IL                              ** LUNAR **
C     DATA IL/0/                                                   ** LUNAR **
C*********************************************************************

C********Note: Current Printout Ceases at 145km but XTOP=27 ==> Z>245km*******
ccjkh8/98      DATA XBOT,XTOP,DX/0.0,27.,-0.05/
C********This variable and the corresponding logic replaces PRNTDX****
      DATA PRNTDZ/4.0/
C For Rashid:
c     DATA XBOT,XTOP,DX/0.0,10.,-0.05/
c     DATA PRNTDZ/1.0/

      DATA (SWE(I),I=1,3)/-0.39,0.39,0.0/
c
c Set height limits depending on model inputs

      if(heatmodel.eq.6)then    !thermosphere 20.-31.4/-.1 rec.
         xbot=20.0
         xtop=31.0
         dx=-0.1
         prntmax=600.
      else                      !tropos-strato-mlt
         xbot=0.0
         xtop=27.0
         dx=-0.05
         prntmax=400.
      endif
c
      pi=acos(-1.)
      tpi=2.*pi
      FREQ=tpi/(PERIOD*86400.)

      NSWE=3
      IF(MOIS.LT.3) NSWE=1
      IF(MOIS.GT.10) NSWE=1
      IF(MOIS.GE.5.AND.MOIS.LE.8) NSWE=2

      IJ=IC/4

C***************************************************************************

  102 FORMAT(1X,6HMONTH=,I2,1H;,1X,5HFLUX=,F4.0,1H;,1X,7HPERIOD=,F4.1,1H;
     + ,1X,2HS=,I2,1H;,1X,3HDX=,F5.2,1H;,1X,3HIJ=,I2)

  103 FORMAT(1X,6HMONTH=,I2,1H;,1X,5HFLUX=,F4.0)
 
      WRITE(7,102)MOIS,FLUX,PERIOD,NZONAL,DX,IJ
      WRITE(7,103)MOIS,FLUX

C***************************************************************************

      OPEN(UNIT=4,STATUS='SCRATCH',ACCESS='DIRECT',recl=894912)
C  This file is opened when diagnostic output is necessary
C   N.B. Remove corresponding comments in JCL file to retreive output
C        from shavano
C     OPEN(UNIT=8,file='migtide.results',STATUS='new')
C     OPEN(UNIT=10,file='migtide.heat',STATUS='new')

C***************************************************************************

C  ARRAYS OF ALTITUDE (Z) AND STRETCHED VAR. (X) ARE CREATED. 
C  TABLES LATER USED FOR X TO Z CONVERSION BY INTERPOLATION.

      CALL HTVEC

C  SET UP VARIOUS LATITUDE ARRAYS WHICH ARE RE-USED AT EACH HEIGHT

      CALL LATVECS(IJ)

C  SET UP BACKGROUND ATMOSPHERE TABLES FOR LATER INTERPOLATION
cc      print*,'nowind in main_2 before setatmos= ',nowind
      CALL SETATMOS
CTEST stop here to write out background only
c      CLOSE(UNIT=4)
c      return
CTEST background end
C  set up Garcia/Solomon eddy diffusion tables for later interpolation

      if(eddymodel.eq.1)then
         call seteddy_gs(mois)
      endif

C  Now calling setheat to add new forcing.  June 26, 1994 (Julie)
C ....and seth2o (M. Hagan 6/29/94)

C  Migrating Forcing

       if(heatmodel.eq.1.or.heatmodel.eq.4.or.heatmodel.eq.5)then
          if(skipo3.eq.0)then   !run UV forcing
             CALL SETHEAT(MOIS)
          endif

          if(skipir.eq.0)then   !run IR forcing
             CALL SETH2O(MOIS,NZONAL)
          endif
       endif
c
c Thermosphere "EUV" forcing

       if(heatmodel.eq.6)then
          print*,'Setting up EUV forcing'
          call set_euv
       endif

c
c Non-migrating heat forcing: infrared

       if(heatnm.eq.1.or.heatnm.eq.2)then      !run IR forcing
          call seth2onm(mois,nzonal)
       endif

C Include latent heat migrating forcing...
C Get 7-year monthy mean latent heat latitude expansion functions:
C	 for use in SR Latent_heat_tides (latent.f) (M. Hagan 5/96)
c        call_latent flag = 1

      if(heatmodel.eq.2.or.heatmodel.eq.4)then
         if(call_late.eq.0)then  !standard latent heat case: monthly mean
            call getlate(nzonal,period,mois,ij)
         endif
      endif

C Need Hough function (2,2) for Hong and Wang latent heat
C Need Hough function (1,1) for force test
c may also use this call for geo.f? hme11?

      if(heatmodel.eq.3.or.heatmodel.eq.5)then
	call hoff(IJ)
      endif
c
c Non-migrating heat forcing: latent

      if(heatnm.eq.0.or.heatnm.eq.2)then !call latent heating model
         call setgci(nzonal,period,mois,ij)
      endif

c  INTEGRATE FROM XTOP TO XBOT IN STEPS OF DX
c  NUMBER OF STEPS WILL BE
      N=(XBOT-XTOP)/DX+1.001


C  COMPUTE ALPHAs AND BETAs FOR LINDZEN-KUO ALGORITHM
C
c
c Loop runs from top (XTOP) to bottom (XBOT)

      DO 3 I=1,N
         X=DX*FLOAT(I-1)+XTOP
C  SET UP MATRICES
         CALL ABCR1(A,B,C,R,I,N,X,IR,IC,IJ)
C  COMPUTE MATRICES A,B,C AND VECTOR R
         CALL ABCRN(A,B,C,I,N,IR,IC,DX)
C  COMPUTE ALPHA AND BETA MATRICES
         CALL ALPBET2(A,B,C,R,I,N,IR,IC,ALPH,BET,AUX,FX,IPVT,ZMX)
 3    CONTINUE
CTEST put stop here to suppress output after abcr1
C      stop
CTEST end
C  BEGIN BACKWARD SWEEP AND FINAL SOLUTION
C

      DO 4 I=1,IR
    4 FNC1(I)=BET(I)
c
c jkh 8/98 write a comment to the output file instead of blank lines (driver.f)

CTEST300   FORMAT(1H1)
CTEST301   FORMAT(/)
CTEST      WRITE(7,300)
CTEST      WRITE(7,301)

      X=XBOT
      zbot=0.0
      if(idenpress.eq.1.or.idecomp.eq.1) then
         OPEN(UNIT=44,STATUS='SCRATCH',ACCESS='DIRECT',recl=894912)
         write(unit=44,rec=1) (fnc1(il),il=1,ir),xbot,zbot
      endif
c
c write output to unit 7

      CALL PRNT(XBOT,FNC1,IR,IJ)
      nhts=1

      Z1=0.0
      LL=N
      DO 5 I=2,N
      LL=LL-1
      
      READ(UNIT=4,REC=LL) ((ALPH(K,J),K=1,IR),J=1,IC),(BET(K),K=1,IR)

C  COMPUTE FUNCTION
      CALL SOL(ALPH,BET,FNC,FNC1,IR,IC)
      X=X-DX

C  write(7,) RESULTS IN STEPS OF PRNTDZ (specified in DATA statement)

      IPRNT=1
      CALL CONV(X,Z)
      nhts=nhts+1
c
c Geopotential height or hough decomposition printout

      if(idenpress.eq.1.or.idecomp.eq.1) then
         write(unit=44,rec=i) (fnc(il),il=1,ir),x,z
      endif

      TESTZ=Z-Z1
      IF (TESTZ .LT. PRNTDZ) IPRNT=0
c
c Enter the maximum altitude for printing results at the top of this program

      if(z.gt.prntmax) iprnt=0  !150 MLT; 400 Thermo

      IF(IPRNT) 6,6,7
    7 CALL PRNT(X,FNC,IR,IJ)
      CALL CONV(X,Z1)
    6 DO 5 J=1,IR
    5 FNC1(J)=FNC(J)
c
c Geopotential height or hough decomposition printout

      if(idenpress.eq.1.or.idecomp.eq.1) then
         idtemp=2*nhts+ij
c        print *,"In call to Post-Processor nhts=",nhts
         call geo(ir,ij,prntdz,nhts,idtemp)
         CLOSE(UNIT=44)
      endif

      CLOSE(UNIT=4)

      return                    !added 9/11/98 jhackney
      END
