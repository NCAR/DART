
      SUBROUTINE EDDY (ZZ,EDDYV,EDDYK,DEDDYV,DEDDYK)
c
      integer iboris, eddymodel,kzzinvar, dissinvar
      integer iondrag, vialdiss, gwstress, idecomp, idenpress
c
      DIMENSION ZE(31),ED(31),DED(31),SED(5,31),S1(3),S2(3)

c  The following integer common block contains flags for the dissipation
c  models

      common /boris/iboris,eddymodel,kzzinvar,dissinvar,iondrag,
     +     vialdiss,gwstress
c
      DATA IPROF/4/
C
      DATA (SED(1,I),I=1,31)/10.,10.,8.,1.1,1.,1.3,3.2,5.4,7.,8.5,9.5,
     19.8,10.,10.,10.,10.,10.,10.,10.,10.,
     29.5,7.5,4.8,2.5,1.,1.,1.,1.,1.,1.,1./
      DATA (SED(2,I),I=1,31)/10.,10.,8.,1.1,1.,1.3,3.2,6.,10.5,17.4,27.,
     136.5,43.5,48.,50.,50.,50.,50.,50.,50.,
     247.,33.,17.,5.6,1.3,1.,1.,1.,1.,1.,1./
      DATA (SED(3,I),I=1,31)/10.,10.,8.,1.1,1.,1.3,3.2,6.,10.5,18.5,41.,
     1110.,255.,410.,500.,500.,500.,500.,500.,500.,
     2445.,325.,172.,54.,4.5,2.,2.,2.,2.,2.,2./
c     DATA (SED(4,I),I=1,31)/10.,10.,8.,1.1,1.,1.3,3.2,6.,10.5,17.4,27.,
c    1  46.,70.,98.,140.,180.,200.,200.,200.,200.,180.,120.,65.,30.,15.,
c    2  8.,1.,1.,1.,1.,1./
C Trying to replicate dissipation in Vial (1986)	M. Hagan (10/24/96):
      DATA (SED(4,I),I=1,31)/10.,10.,8.,1.1,1.,1.3,3.2,6.,10.5,15.,20.,
     1  25.,30.,35.,50.,65.,100.,200.,200.,200.,200.,200.,100.,60.,40.,
     2  25.,5.,1.,1.,1.,1./
c     DATA (SED(4,I),I=1,31)/10.,10.,8.,1.1,1.,1.3,3.2,6.,10.5,15.,20.,
c    1  25.,30.,35.,50.,65.,100.,250.,300.,300.,300.,250.,100.,60.,40.,
c    2  25.,5.,1.,1.,1.,1./
      DATA (SED(5,I),I=1,31)/10.,10.,8.,1.1,1.,1.3,3.2,6.,10.5,15.,20.,
     1  30.,40.,50.,70.,90.,4*100.,90.,60.,32.,15.,7.5,4.,5*1./
C
      EDDYV=0.0
      EDDYK=0.0
      DEDDYV=0.0
      DEDDYK=0.0
      IF(ZZ.GT.150.) GO TO 10
C
      DO 1 I=1,31
      ED(I)=SED(IPROF,I)
    1 ZE(I)=5.*(I-1)
      H=5000.
      CALL DET3(H,ED,DED,31,IER)
C
      DO 2 I=1,31
    2 DED(I)=DED(I)/ED(I)
C
      CALL ATSM(ZZ,ZE,ED,31,1,S1,S2,3)
      CALL ALI(ZZ,S1,S2,EDDYV,3,1.0E-03,IER)
      CALL ATSM(ZZ,ZE,DED,31,1,S1,S2,3)
      CALL ALI(ZZ,S1,S2,DEDDYV,3,1.0E-06,IER)

      if(vialdiss.eq.0)then
         EDDYK=1.36*EDDYV
         DEDDYK=1.36*DEDDYV

C Comment above and Uncomment the following to TEST********************
C ************* mimic dissipation in Vial (1986); Prantl #=1

      elseif(vialdiss.eq.1)then
         EDDYK=1.*EDDYV
         DEDDYK=1.*DEDDYV
      endif
C ***************************TEST********************
C
   10 RETURN
      END
