C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C######################### COOL_ir.f ########################
C
C  To calculate the IR radiative cooling rate (K/day) due to CO2, 
C   H2O, and O3, for middle atmosphere modeling (8/95).

************************* NOTE :  All units are in MKS *******************

******************************************************************************
*  Dr. Xun Zhu of Johns Hopkins (zhu@gibbs.eps.jhu.edu).                     *
*  The following changes were made:                                          *
*                                                                            *
*  - All the routines in this file were changed to single precision.         *
*                                                                            *
*  - Calculation of air number density was modified to be done only once.    *
*                                                                            *
*  - The air density, AirDen is added to the input parameter list            *
*    AirDen  should have the units of m^-3 .                                 *
*                                                                            *
*  - The input file name is set in msetvar.f.                                *
*                                                                            * 
*  - O3 and H2O are input in the form of volume mixing ratios.               *
*                                                                            *
*  - For now, mixing ratios of carbon dioxide, molecular oxygen and water    *
*    and number density of atomic oxygen are 2D (latitude x pressure).       * 
*                                                                            *
*  - QIR in the routine COOL_IR is changed to the positive sums of the       *
*    contributions from CO2, O3, and H2O, so that the resulting total        *
*    cooling (QIR) is negative.                                              *
*                                                                            *
*  - At the suggestion of X. Zhu, the functions DETIJ, AVT21, AVT32,         *
*    SPT, and DETIJR were incorporated into the body of the main subroutines *
*                                                                            *
*  - The subroutines RADCO2 and RADO3 have been split into two parts:        *
*    CMATCO2 calculates non-LTE CO2 Curtis matrix from zonal mean temperature*
*    CMATO3 calculates non-LTE O3 Curtis matrix from zonal mean temperature  *
*                                                                            *
*    CRATECO2 calculates CO2 cooling rate on model grid                      *
*    CRATEO3 calculates O3 cooling rate on model grid                        *
*                                                                            *
*    Arrays eco2(nzz,nzz,ny) and eo3(nzz,nzz,ny) store non-LTE Curtis        *
*    matrices for CO2 and O3, respectively                                   *
*                                                                            *
*    IF VERICAL MODEL GRID IS CHANGED....                                    *
*    => recalculate curtis matrix (see /aksmithpc2/aksmith/ROSE-t/radiation  *
*    => change parameter k1 in subr. nltec                                   *
*                                                                            * 
*     Modifications done by Rashid Khosravi and Anne Smith                   * 
******************************************************************************

      SUBROUTINE COOL_IR(AirDen, Temp, co2mr, o3int, h2o,
     +                  TCLTOP, QIR, matrx)
    
      use params

      implicit none

      integer ix, iy, k, jj, kk, jsep
      real :: tcltop, div, Z0(nzz), P0(nzz), T0(nzz), GDP(nzz)
      COMMON /CTXZ23/ Z0, P0, T0, GDP   ! P0 read from curtis.dat, in Pascal

      logical matrx                 ! not used

      Real AirDen(nzz, nx, ny)      ! Air density (m^-3)

      Real co2mr(nzz),              ! volume mixing ratio of CO2, specified
     +     Temp(nzz,nx,ny),         ! temperature profile at current time step
     +     o3int(nzz,ny),           ! volume mixing ratios of ozone
     +     h2o(nzz,ny),             !                         water
     +     QIR(nzz,nx,ny)           ! the total cooling rate (K/s)

      Real  QCO2(nzz), QO3(nzz), QH2O(nzz)

      Real  AirDenZ(nzz), TempZ(nzz), OdenZ(nzz),  ! vertical profiles
     +      CO2mrZ(nzz), O3mrZ(nzz), H2OmrZ(nzz)   ! mass mixing ratios

      real :: eeco2(nzz,nzz),eeo3(nzz,nzz),atomico(nzz)

      logical pass1   ! read the Curtis matrices only on first call to QIRMD 
      save  pass1, CO2mrZ, jsep
      data  pass1 /.true./

c  atomic oxygen mixing ratio profile from 2d model
      data atomico / 
     $   .14765E-15, .23172E-15, .34407E-15, .55346E-15, .10230E-14,
     $   .22384E-14, .57503E-14, .28800E-13, .19722E-12, .12716E-11,
     $   .58890E-11, .23797E-10, .73259E-10, .19980E-09, .51435E-09,
     $   .14738E-08, .43645E-08, .13991E-07, .39109E-07, .79436E-07,
     $   .13590E-06, .21770E-06, .33900E-06, .53048E-06, .81931E-06,
     $   .12003E-05, .15826E-05, .18743E-05, .21495E-05, .27104E-05,
     $   .42505E-05, .97450E-05, .29889E-04, .11346E-03, .34680E-03,
     $   .84696E-03, .16714E-02, .29281E-02, .46882E-02, .72132E-02,
     $   .10884E-01, .16792E-01, .25888E-01, .39885E-01, .59933E-01/

      If (pass1) Then  ! read the Curtis matrices
        call inputc
        do k=1,nzz
           CO2mrZ(k) = 1.5225*co2mr(k) ! convert volume to mass m.r. 
        end do
        pass1 = .false.

        JSEP=1     ! 0 => LTE; 1 => nonLTE in z>70 km; 2 => complete nonLTE
        matrx = .true.
      EndIF

      div = 1./float(nx)

      Do iy=1, ny
         do k=1,nzz
            TempZ(k)   = 0.
            AirDenZ(k) = 0.
            Do ix=1, nx
               AirDenZ(k) = AirDenZ(k) + AirDen(k,ix,iy)
               TempZ(k)   = TempZ(k) + Temp(k,ix,iy)
            End Do
            AirDenZ(k) = AirDenZ(k) * div
            TempZ(k)   = TempZ(k) * div
  
            OdenZ(k) = atomico(k) * AirDenZ(k)  ! convert to # density (m^-3)
            O3mrZ(k) = 1.661 * o3int(k, iy)     ! convert to mass m.r.
            H2OmrZ(k) = 0.622*h2o(k,iy)         !          .....
         end do

         CALL cmatCO2(AirDenZ,TempZ,CO2mrZ,OdenZ,eeco2,JSEP)
         CALL cmatO3(AirDenZ,TCLTOP,TempZ,O3mrZ,eeo3,JSEP)

         Do ix=1, nx
            do k=1,nzz
               AirDenZ(k) = AirDen(k,ix,iy)
               TempZ(k)   = Temp(k,ix,iy)
c               OdenZ(k) = 3.e9                 ! convert to # density (m^-3)
            end do

            CALL crateCO2(AirDenZ,TempZ,CO2mrZ,OdenZ,eeco2,QCO2)
            CALL crateO3(AirDenZ,TCLTOP,TempZ,O3mrZ,eeo3,QO3)
            CALL RADH2O(AirDenZ, TempZ, H2OmrZ, QH2O)

            do k=1,nzz
               QIR(k,ix,iy)= (QCO2(k) + QO3(k) + QH2O(k))
     $                         /86400.0 ! total IR cooling rate in K/s
            end do
         End Do
      End Do

      RETURN
      END
      



C  ========================= CO2 matrix : =================================   

      SUBROUTINE cmatCO2(AirDen, T, RR, BOXY, ee, JSEP)

C  To calculates the cooling rate (QCO2) by Curtis matrix CURT for
C  given p, T, [O], and r[CO2] deviated from [CO2]0 around top levels. 
C  CURTZ are the off-line calculated Curtis matrices. RR is mass mixing ratio.
C  The units of [O] and p are  M^-3 and Pascal, respectively.
C  --> quadrature interpolation.

CC  SPT ->
CC  To calculate the temperature dependence of band intensity S
CC  in Houghton's appedices. N=1,2,3 and 4 correspond to 626 
CC  fundamental band, 626 first hot band, isotopes fundamental 
CC  band and total 15 micron band, respectively. Only used in
CC  non-LTE calculations.


      use params

      COMMON /CTXZ23/ Z0(nzz),P0(nzz),T0(nzz),GDP(nzz)
      COMMON /CURTZ2/ R20(nzz),H200(nzz),CURT0(nzz,nzz)
     & ,A10(nzz,nzz),A01(nzz,nzz),A20(nzz,nzz),A02(nzz,nzz),A11(nzz,nzz)

      DIMENSION CURT(nzz,nzz),EE(nzz,nzz),DTJ(nzz)

C  EE(nzz,nzz) is used both for diagonal E matrix and non-LTE Curtix matrix
      DIMENSION T(nzz),RR(nzz),BOXY(nzz),QCO2(nzz)
 
      Dimension AirDen(nzz)

      DATA AA21/1.43/,AA32/1.94/


      DTJ = T-T0
      EE=0.0
      CURT = CURT0
      DO K=1,nzz                 ! Get the Curtis matrix by interpolation
         DO J=1,nzz
            IF(ABS(CURT0(J,K)).gt.0.005) then
               K1=MIN0(j,k)
               K2=MAX0(j,k)
               XX=0.0
               YY=0.0
               DO kk=K1,K2
                  XX = XX + GDP(kk)
                  YY = YY + (T(kk)-T0(kk)) * GDP(kk)
               end do
               DELTA=YY/XX
               CURT(J,K)=CURT0(J,K)+A10(J,K)*DTJ(J)+A01(J,K)*DELTA
               CURT(J,K)=CURT(J,K)+0.5*(A20(J,K)
     &              *DTJ(J)**2+A02(J,K)*DELTA**2)+A11(J,K)*DELTA*DTJ(J)
            end if
         end do
      end do

      IF(JSEP.EQ.0) THEN
         ee = curt
         return
      ENDIF

      DO J=1,nzz
         IF(T(j).LT.200.0) VT1M=2.5E-21
         IF(T(j).GE.200.0) VT1M=2.5E-21*(1.0+0.03*(T(j)-200.))
         AVT21=AirDen(j)*VT1M+BOXY(j)*5.5E-18 ! revised k-rate for [O]
         TRA=T(j)/273.3
         VT2M=1.24E-20*TRA*TRA
         AVT32=AirDen(j)*VT2M+BOXY(j)*5.5E-18 ! revised k-rate for [O]
         TX=T(j)
         IF(TX.GT.320.0) TX=350.0
         IF(TX.LT.130.0) TX=130.0
         X=TX-250.0
         xsq = x*x
         Y1=5.0522-5.4693E-4*X-5.2873E-7*xsq
         Y2=3.8067+7.0636E-3*X-2.1024E-5*xsq+1.964E-8*X*X*X
         Y3=3.2338-5.5612E-4*X-5.3153E-7*xsq
         Y4=5.0961+8.8837E-6*X+3.1662E-8*xsq

         PHIN=((10.0**Y1+10.0**Y3)*AVT21/AA21
     &        +10.0**Y2*AVT32/AA32)/(10.0**Y4)
         EE(J,J)=R20(J)/(RR(J)*2.0*H200(J)*PHIN) ! EJJ=R0/(R*2*H00)
      end do
      CALL NLTEC(CURT,EE,JSEP)

      RETURN
      END



C  ========================= CO2 cooling : =================================   

      SUBROUTINE crateCO2(AirDen, T, RR, BOXY, ee, QCO2)

C  To calculates the cooling rate (QCO2) by Curtis matrix CURT for
C  given p, T, [O], and r[CO2] deviated from [CO2]0 around top levels. 
C  CURTZ are the off-line calculated Curtis matrices. RR is mass mixing ratio.
C  The units of [O] and p are  M^-3 and Pascal, respectively.


      use params

      COMMON /CTXZ23/ Z0(nzz),P0(nzz),T0(nzz),GDP(nzz)
      COMMON /CURTZ2/ R20(nzz),H200(nzz),CURT0(nzz,nzz)
     & ,A10(nzz,nzz),A01(nzz,nzz),A20(nzz,nzz),A02(nzz,nzz),A11(nzz,nzz)

      DIMENSION THE1(nzz),EE(nzz,nzz)

C  EE(nzz,nzz) is non-LTE Curtix matrix
      DIMENSION T(nzz),RR(nzz),BOXY(nzz),QCO2(nzz)
 
      Dimension AirDen(nzz)

      DATA AA21/1.43/,AA32/1.94/

      QCO2=0.0
      THE1=47.612466/(EXP(970.97/T)-1.0)    ! 47.612466=EXP(970.97/250)-1

      DO K=1,nzz
         DO J=1,nzz
            QCO2(J)=QCO2(J)+EE(J,K)*THE1(K)
         end do
      end do

      DO J=1,nzz
         IF(Z0(J).gt.80.0E3) then
            FACX=0.5*(1.0+TANH((Z0(J)-85.0E3)/7.0E3))
            QCO2(J)=QCO2(J)*((RR(J)/R20(J))*FACX+(1.0-FACX))
         end if
      end do

      IF(Z0(nzz).LE.100.01E3) THEN 
         DO K=1,nzz
            XXZ1=(Z0(nzz)+1.0E3-Z0(K))/5.0E3
            QCO2(K)=QCO2(K)*(1.0-EXP(-XXZ1**2))
         end do
      ENDIF

      RETURN
      END
  



      SUBROUTINE INVERT(D,N,M)

C  To invert the matrix D,  M=N+1. Mth row can be any values.

      DIMENSION D(N,M)

      KK=0
      JJ=0
      DO K=1,N
         DO J=1,N
            D(J,M)=0.0
         end do
         D(K,M)=1.0
         JJ=KK+1
         LL=JJ
         KK=KK+1
 20      IF(ABS(D(JJ,KK).le.1.0E-30)) then
            JJ=JJ+1
            IF(JJ.le.N) go to 20
            WRITE(*,98)jj,kk,n
 98         FORMAT('ERROR in input to Curtis matrix',3i4)
            stop
         end if
         IF(LL.ne.JJ) then
            DO MM=1,M
               DTEMP=D(JJ,MM)
               D(LL,MM)=D(JJ,MM)
               D(JJ,MM)=DTEMP
            end do
         end if
         DIV=D(K,K)
         DO LJ=1,M
            D(K,LJ)=D(K,LJ)/DIV
         end do
         DO I=1,N
            IF(I.ne.K) then
               FAC=D(I,K)
               DO LJ=1,M
                  D(I,LJ)=D(I,LJ)-FAC*D(K,LJ)
               end do
            end if
         end do
         DO J=1,N
            D(J,K)=D(J,M)
         end do
      end do

      RETURN
      END


C ======================== Ozone Matrix : ==============================

      SUBROUTINE cmatO3(AirDen, TCLTOP, T, RR, ee, JSEP)

C  To calculate the O3 cooling rate by Curtis matrix CURT.  !  ZM(M)>ZM(1)
c  TCLTOP=effective cloud top temperature.


      use params


      PARAMETER (LBTRP=1)    !!! specified by users

      COMMON /CTXZ23/ Z0(nzz),P0(nzz),T0(nzz),GDP(nzz)
      COMMON /CURTZ3/ R30(nzz),H300(nzz),ZCURT0(nzz,nzz)
     & ,ZA10(nzz,nzz),ZA01(nzz,nzz),ZB01(nzz,nzz),CLB0(nzz,LBTRP)
     & ,A10LB(nzz,LBTRP),A01LB(nzz,LBTRP),B01LB(nzz,LBTRP)

      DIMENSION CURT(nzz,nzz),EE(nzz,nzz),CLB3(nzz),DTJ(nzz)
C  EE(nzz,nzz) is used both for diagonal E matrix and non-LTE Curtix matrix

      Dimension AirDen(nzz)

      DIMENSION T(nzz),RR(nzz),WK1(nzz)


      DTJ = T-T0
      CURT = ZCURT0

      DO K=1,nzz
         DO J=1,nzz             ! Get the Curtis matrix by interpolation
            IF(ABS(ZCURT0(J,K)).gt.0.005) then
               K1=MIN0(j,k)
               K2=MAX0(j,k)
               XX=0.0
               YY=0.0
               XXr=0.0
               YYr=0.0
               DO kk=K1,K2
                  XX=XX+GDP(kk)
                  YY=YY+(T(kk)-T0(kk))*GDP(kk)
                  XXr=XXr+0.2*R30(kk)*GDP(kk)
                  YYr=YYr+(RR(kk)-R30(kk))*GDP(kk)
               end do
               DETIJ=YY/XX
               DETIJR=YYr/XXr
               CURT(J,K)=ZCURT0(J,K)+ZA10(J,K)*DTJ(J)
     &                  +ZA01(J,K)*DETIJ
     &                  +ZB01(J,K)*DETIJR
            end if
         end do
      end do
      DO J=1,nzz
         WK1(J)=(RR(J)/R30(J))*H300(J)
      end do

      DO K=1,nzz
         DO J=1,nzz
            EE(J,K)=CURT(J,K)*WK1(J)
         end do
      end do
      DO K=1,nzz
         DO J=1,nzz
            CURT(J,K)=EE(J,K)
            EE(J,K)=0.0
         end do
      end do

      IF(JSEP.EQ.0) THEN
         ee = curt
         return
      ENDIF

      DO J=1,nzz
         WK1(J)=3.64E-21*AirDen(J) ! WK1=(4.3E-20/11.81)*air density
         EE(J,J)=R30(J)/(RR(J)*2.0*H300(J)*WK1(J)) ! EJJ=R0/(R*2*H00)
      end do
      CALL NLTEC(CURT,EE,JSEP)

      RETURN
      END



C ======================== Ozone cooling : ==============================

      SUBROUTINE crateO3(AirDen, TCLTOP, T, RR, ee, QO3)

C  To calculate the O3 cooling rate by Curtis matrix CURT.  !  ZM(M)>ZM(1)
c  TCLTOP=effective cloud top temperature.


      use params


      PARAMETER (LBTRP=1)    !!! specified by users

      COMMON /CTXZ23/ Z0(nzz),P0(nzz),T0(nzz),GDP(nzz)
      COMMON /CURTZ3/ R30(nzz),H300(nzz),ZCURT0(nzz,nzz)
     & ,ZA10(nzz,nzz),ZA01(nzz,nzz),ZB01(nzz,nzz),CLB0(nzz,LBTRP)
     & ,A10LB(nzz,LBTRP),A01LB(nzz,LBTRP),B01LB(nzz,LBTRP)

      DIMENSION CURT(nzz,nzz),THE1(nzz),EE(nzz,nzz),CLB3(nzz),DTJ(nzz)
C  EE(nzz,nzz) is used both for diagonal E matrix and non-LTE Curtix matrix

      Dimension AirDen(nzz)

      DIMENSION T(nzz),RR(nzz),QO3(nzz),WK1(nzz)

      KLBD=1
      XTHE1=484.0/(EXP(1546.0/TCLTOP)-1.0)

      THE1 = 484.0/(EXP(1546.0/T)-1.0)
      DTJ = T-T0

      DO J=1,nzz
         K1=MIN0(j,klbd)
         K2=MAX0(j,klbd)
         XX=0.0
         YY=0.0
         XXr=0.0
         YYr=0.0
         DO K=K1,K2
            XX=XX+GDP(K)
            YY=YY+(T(K)-T0(K))*GDP(K)
            XXr=XXr+0.2*R30(K)*GDP(K)
            YYr=YYr+(RR(K)-R30(K))*GDP(K)
         end do
         DETIJ=YY/XX
         DETIJR=YYr/XXr
         CLB3(J)=CLB0(J,KLBD)+A10LB(J,KLBD)*DTJ(J)
     &          +A01LB(J,KLBD)*DETIJ
     &          +B01LB(J,KLBD)*DETIJR
      end do

      WK1  = (RR/R30) * H300
      CLB3 = CLB3 * WK1
      QO3  = CLB3 * XTHE1

      DO K=KLBD,nzz
         DO J=1,nzz
            QO3(J)=QO3(J)+EE(J,K)*THE1(K)
         end do
      end do

      RETURN
      END



C ======================== Water vapor cooling : ========================

      SUBROUTINE RADH2O(AirDen, T, RR, QH2O)

C  To calculate the H2O cool-to-space cooling rate.  !  ZM(nzz)>ZM(1)


      use params


      COMMON /CTXZ23/ Z0(nzz),P0(nzz),T0(nzz),GDP(nzz)
      COMMON /GAMWAT/ R10(nzz),GAMR1(nzz),GAMV1(nzz)
     & ,YAR1(nzz),YAV1(nzz),YAR2(nzz),YAV2(nzz)
      DIMENSION GR2(nzz),GV2(nzz),THR(nzz),THV(nzz),PHIV(nzz)
      DIMENSION T(nzz),RR(nzz),QH2O(nzz),WK1(nzz),WK2(nzz)

      Dimension AirDen(nzz)


      DO K=1,nzz
         K1=MIN0(k,nzz)
         K2=MAX0(k,nzz)
         XXr=0.0
         YYr=0.0
         DO kk=K1,K2
            XXr=XXr+0.5*R10(kk)*GDP(kk)
            YYr=YYr+(RR(kk)-R10(kk))*GDP(kk)
         end do
         DETIJR=YYr/XXr
         WK1(K)=DETIJR
         WK2(K)=RR(K)/3.0E-6    ! 3.0E-6 => reference mass mixing ratio
         PHIV(K)=2.1E-21*AirDen(K) ! (4.3E-20/20.5)*air density in m^-3
      end do

      THR = 81.6/(EXP(568.01/T)-1.0)             ! S*B in K/day
      THV = 2730.0/(EXP(2300.8/T)-1.0)           ! S*B in K/day
      GR2 = GAMR1 + YAR1 * WK1 + YAR2 * WK1**2
      GV2 = GAMV1 + YAV1 * WK1 + YAV2 * WK1**2
      QH2O = -WK2 * (THR * GR2
     &     + THV * GV2 / (1.0 + 0.5 * GV2/PHIV))

      RETURN
      END


C  ==================== Input of the Curtis matrices ====================

      SUBROUTINE INPUTC

C --- Input the Curtis matrices for cooling rate calculations


      use params
      use dynam, only : namf40
      use utilities_mod, only : open_file, close_file

      integer :: iunit
      PARAMETER (LBTRP=1)      !!! specified by users

      COMMON /CTXZ23/ Z0(nzz),P0(nzz),T0(nzz),GDP(nzz)
      COMMON /GAMWAT/ R10(nzz),GAMR1(nzz),GAMV1(nzz)
     & ,YAR1(nzz),YAV1(nzz),YAR2(nzz),YAV2(nzz)
      COMMON /CURTZ2/ R20(nzz),H200(nzz),CURT0(nzz,nzz)
     & ,A10(nzz,nzz),A01(nzz,nzz),A20(nzz,nzz),A02(nzz,nzz),A11(nzz,nzz)
      COMMON /CURTZ3/ R30(nzz),H300(nzz),ZCURT0(nzz,nzz)
     & ,ZA10(nzz,nzz),ZA01(nzz,nzz),ZB01(nzz,nzz),CLB0(nzz,LBTRP)
     & ,A10LB(nzz,LBTRP),A01LB(nzz,LBTRP),B01LB(nzz,LBTRP)


  52  FORMAT(8E12.5)
  53  FORMAT(2F12.6,5E12.3)
  55  FORMAT(7E14.5)

!... open radiation file

      iunit = open_file(namf40,form='formatted',action='read')     

c     open(unit=40, file=namf40, form='formatted', status='old')

      DO K=1,nzz
         READ(iunit,'(8E12.5)') Z0(K),T0(K),P0(K),GDP(K),
     &               R20(K),H200(K),R30(K),H300(K)
         READ(iunit,'(2F12.6,5E12.3)') GAMR1(K),GAMV1(K),R10(K),
     &               YAR1(K),YAV1(K),YAR2(K),YAV2(K)
      end do
      READ(iunit,'(7E14.5)') ((CURT0(L1,L2),L1=1,nzz),L2=1,nzz)
      READ(iunit,'(7E14.5)') ((A10(L1,L2),L1=1,nzz),L2=1,nzz)
      READ(iunit,'(7E14.5)') ((A01(L1,L2),L1=1,nzz),L2=1,nzz)
      READ(iunit,'(7E14.5)') ((A20(L1,L2),L1=1,nzz),L2=1,nzz)
      READ(iunit,'(7E14.5)') ((A02(L1,L2),L1=1,nzz),L2=1,nzz)
      READ(iunit,'(7E14.5)') ((A11(L1,L2),L1=1,nzz),L2=1,nzz)
      READ(iunit,'(7E14.5)') ((ZCURT0(L1,L2),L1=1,nzz),L2=1,nzz)
      READ(iunit,'(7E14.5)') ((ZA10(L1,L2),L1=1,nzz),L2=1,nzz)
      READ(iunit,'(7E14.5)') ((ZA01(L1,L2),L1=1,nzz),L2=1,nzz)
      READ(iunit,'(7E14.5)') ((ZB01(L1,L2),L1=1,nzz),L2=1,nzz)
      READ(iunit,'(7E14.5)') ((CLB0(L1,L2),A10LB(L1,L2)
     &             ,A01LB(L1,L2),B01LB(L1,L2),L1=1,nzz),L2=1,LBTRP)

      call close_file(iunit)

      RETURN
      END
  

      


C ==================== get non-LTE Curtis matrix : =======================

      SUBROUTINE NLTEC(CURT,EE,JSEP)

C  To calculate nonLTE curtis matrix from matrices C and E (JGR, 97, 12,790).
C  The output nonLTE curtis matrix is stored in EE.        !  ZM(nzz)>ZM(1)
C  JSEP=0 => LTE; JSEP=1 => nonLTE in top K1 region; JSEP>1 => complete nonLTE


      use params


      PARAMETER (K1=18)    !!! specified by users
      PARAMETER (KMP=nzz+1,K1P=K1+1,K2=nzz-K1)
      DIMENSION CURT(nzz,nzz),EE(nzz,nzz),AINV(nzz,KMP)
      DIMENSION BINV(K1,K1P),WOK(K2,K1)


      IF(JSEP.gt.1) then
         DO K=1,nzz
            DO J=1,nzz
               XX=0.0
               DO LL=1,nzz
                  XX=XX+CURT(J,LL)*EE(LL,K)
               end do
               AINV(J,K)=-XX
            end do
         end do
         DO K=1,nzz
            AINV(K,K)=AINV(K,K)+1.0
         end do
         CALL INVERT(AINV,nzz,KMP)
      end if
      IF(JSEP.eq.1) then
         DO K=1,K1
            DO J=1,K1
               XX=0.0
               DO LL=1,K1
                  XX=XX+CURT(K2+J,K2+LL)*EE(K2+LL,K2+K)
               end do
               BINV(J,K)=-XX
            end do
         end do
         DO K=1,K1
            BINV(K,K)=BINV(K,K)+1.0
         end do
         CALL INVERT(BINV,K1,K1P)

         DO K=1,K2
            DO J=1,nzz
               AINV(J,K)=0.0
            end do
         end do
         DO K=1,K2
            AINV(K,K)=1.0
         end do
         DO K=1,K1
            DO J=1,K1
               AINV(K2+J,K2+K)=BINV(J,K)
            end do
         end do
         DO K=1,K1
            DO J=1,K2
               XX=0.0
               DO LL=1,K1
                  XX=XX+CURT(J,K2+LL)*EE(K2+LL,K2+K)
               end do
               WOK(J,K)=XX
            end do
         end do
         DO K=1,K1
            DO J=1,K2
               XX=0.0
               DO LL=1,K1
                  XX=XX+WOK(J,LL)*BINV(LL,K)
               end do
               AINV(J,K2+K)=XX
            end do
         end do
      end if

      DO K=1,nzz
         DO J=1,nzz
            XX=0.0
            DO LL=1,nzz
               XX=XX+AINV(J,LL)*CURT(LL,K)
            end do
            EE(J,K)=XX
         end do
      end do

      RETURN
      END


