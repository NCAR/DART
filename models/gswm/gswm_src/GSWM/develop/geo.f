c*****************************************************************************
C       Modified 2/98 J. Hackney:  2 subroutines previously called
C       geoht.f and geoproc.f combined into one subroutine for efficiency.
C       Calculation and print options described below are controlled by the
C       0/1 flags idecomp and idenpress.
C
C     This routine is called by main_2 if either idecomp or idenpress=1.
C
C     This subroutine calculates geopotential height between x=0 and 24 (150km)
C	using the surf1 and surfd routines to calculate dvdt and dwdz
C	and prints it out at the altitudes of winds and temps
C					---Maura Hagan (2/93)
C	Modifications to perfor Hough Mode decomposition of geopotential
C		heights --M. Hagan (9/95) geoht.f: idecomp
C       Modifications to calculate and print perturbation pressure and
C               density --M. Hagan (8/95) geoproc.f: idenpress
c*****************************************************************************
C
      subroutine geo(ir,ij,prntdz,nhts,idtemp)

      real zg(nhts),xg(nhts),temp(idtemp)
      real vgr(ij,nhts),wgr(ij,nhts),vgtr(ij,nhts,3),wgtr(ij,nhts,3)
      real vgi(ij,nhts),wgi(ij,nhts),vgti(ij,nhts,3),wgti(ij,nhts,3)
      complex geopot(ij),fnc(ir),xi,u,v,w,t,dwdz,dvdt
      complex ug(ij,nhts),vg(ij,nhts),wg(ij,nhts),tg(ij,nhts)
      complex delrho(ij),delpress(ij), paren

      DIMENSION TR(6,ij),TI(6,ij),XXR(ij),XXI(ij),FR(ij),FI(ij)
      COMPLEX FX(6),GPS(1),GPT(ij)

      COMMON/ATM1/P1(91),P2(91),P3(91),P4(91),P5(91),P6(91),
     +P7(91),P8(91),P9(91),P10(91),P11(91),P12(91),P13(91),P14(91)

      COMMON/ATM2/DP1(91),DP3(91),DP8(91),DP9(91),D2P8(91),P15(91)

      COMMON/LATVECT/DTHETA,CLATT(91),XLAT(91),SNCLAT(91),CSCLAT(91),
     +   TNCLAT(91),CTNCLAT(91),SIN2I(91),GMLAT(91),DIP(91)

      COMMON/MODE/NZONAL,PERIOD,FREQ,MOIS,NSS,FLUX

      COMMON/HOFFS/H11(91),H12(91),H13(91),H1M2(91),H1M1(91)
     +,H1M4(91),H22(91),H23(91),H24(91),H25(91),H26(91)
     + ,H33(91),H34(91),H35(91),H36(91)

      common/postproc/idecomp,idenpress,mig
	
      data wr,xi/6378165.,(0.0,1.0)/

      fn=float(nzonal)
      RADEG = 180./ACOS(-1.)

C Hough Decomposition:  idecomp=1

      If(idecomp.eq.1)then

C Index Value at the equator:

        IEQ=(IJ-1)/2

        IF(NZONAL.lt.2)then
          II11=1
          II12=1
          II13=1
          II1M2=1
          II1M4=1
          II1M1=1
	  if(nzonal.eq.1) print *,"Diurnal Mode Flags are: ",
     +		II11,II12,II13,II1M2,II1M4,II1M1
        elseif(nzonal.eq.2)then
          II22=1
          II23=1
          II24=1
          II25=1
          II26=1
	  print *,"Semidiurnal Mode Flags are: ",
     +		II22,II23,II24,II25,II26
        else    !nzonal>2
          II33=1
          II34=1
          II35=1
          II36=1
        endif

C Tell us where the program is and what it's doing

	print *,"Calling Hoff from Geoht with ij=", IJ

C  COMPUTE HOUGH FUNCTIONS
C  This can be done from outside this routine, in main_2.  Use the idecomp
C  flag to call hoff in main_2

        CALL HOFF(IJ)

C Format Hough Decomposition output

  500   FORMAT(1X,5HGP11=,e8.3,1H/,F4.1,1X,5HGP12=,e8.3,1H/,F4.1,1X,
     +    5HGP13=,e8.3,1H/,F4.1,1X,6HGP1M2=,e8.3,1H/,F4.1,1X,
     +    6HGP1M4=,e8.3,1H/,F4.1,1X,6HGP1M1=,e8.3,1H/,F4.1)
  600   FORMAT(1X,5HGP22=,e8.3,1H/,F4.1,1X,5HGP23=,e8.3,1H/,F4.1,1X,
     +    5HGP24=,e8.3,1H/,F4.1,1X,5HGP25=,e8.3,1H/,F4.1,1X,5HGP26=,
     +    e8.3,1H/,F4.1)
  700    FORMAT(1X,4HT33=,e8.3,1H/,F4.1,1X,4HT34=,e8.3,1H/,F4.1,1X,
     +    4HT35=,e8.3,1H/,F4.1,1X,4HT36=,e8.3,1H/,F4.1)

      endif   !idecomp=1

C Analysis output uses these formats (200 Hough, 201 Den/Prss, 203 Both)

100   FORMAT(1X,2HZ=,F7.3)
200     FORMAT(1X,4HLAT=,F5.1,2X,6HGeoHt=,e8.3,1H/,F4.1,
     +		2X,4HGPS=,e8.3,1H/,F4.1)
201     FORMAT(1X,4HLAT=,F5.1,2X,4HGeo=,e8.3,1H/,F4.1,1x,
     +          4HRho=,e8.3,1H/,F4.1,1x,2HP=,e8.3,1H/,F4.1)
203     FORMAT(1X,4HLAT=,F5.1,2X,6HGeoHt=,e8.3,1H/,F4.1,1x,
     +          4HGPS=,e8.3,1H/,F4.1,1x,4HRho=,e8.3,1H/,F4.1,
     +          1x,2HP=,e8.3,1H/,F4.1)

101   FORMAT(1X,11h97 loop: I=,i3,1x,2HZ=,F7.3)  !not used
102   format(1x,i3,9e11.3)  !not used


C Both analyses use this do loop
     
      do 99 i=1,nhts

	 read(unit=44,rec=i) (fnc(il),il=1,ir),xg(i),zg(i)

	 jj=0
         do 98 j=0,ij-1
	    jj=jj+1
            iu=4*j+1
            iv=4*j+2
            iw=4*j+3
            it=4*j+4
            ug(jj,i)=fnc(iu)
            vg(jj,i)=fnc(iv)
            wg(jj,i)=fnc(iw)
            tg(jj,i)=fnc(it)
            vgr(jj,i)=real(vg(jj,i))
            wgr(jj,i)=real(wg(jj,i))
            vgi(jj,i)=aimag(vg(jj,i))
            wgi(jj,i)=aimag(wg(jj,i))
   98    continue
	
   99 continue

C Both analyses must interpolate data fields

      sigma=1.

      call surf1(ij,nhts,clatt,xg,vgr,ij,vx1,vxm,vy1,vyn,vxy11,
     +		vxym1,vxy1n,vxymn,255,vgtr,temp,sigma,ierr)
      call surf1(ij,nhts,clatt,xg,wgr,ij,wx1,wxm,wy1,wyn,wxy11,
     +		wxym1,wxy1n,wxymn,255,wgtr,temp,sigma,ierr)
      call surf1(ij,nhts,clatt,xg,vgi,ij,vx1,vxm,vy1,vyn,vxy11,
     +		vxym1,vxy1n,vxymn,255,vgti,temp,sigma,ierr)
      call surf1(ij,nhts,clatt,xg,wgi,ij,wx1,wxm,wy1,wyn,wxy11,
     +		wxym1,wxy1n,wxymn,255,wgti,temp,sigma,ierr)

C Both analyses use Z1 for PRNTDZ logic:

      Z1=0.0

      xmod=0.0
      xmodsqr=0.0

      do 97 i=1,nhts
	 x=xg(i)
	 z=zg(i)

	 call varx(z,f1,f2)

	 call atmos5(z,ij,g,gam,gam1)

C Both analyses interpolate

         do 96 j=1,ij
	    u=ug(j,i)
            v=vg(j,i)
            w=wg(j,i)
            t=tg(j,i)

            colat=clatt(j)

            call surfd(colat,x,vr,dvdtr,dvdxr,d2vdtr,d2vdxdtr,d2vdxr,
     +       ij,nhts,clatt,xg,vgr,ij,vgtr,sigma)

            call surfd(colat,x,wwr,dwdtr,dwdxr,d2wdtr,d2wdxdtr,d2wdxr,
     +       ij,nhts,clatt,xg,wgr,ij,wgtr,sigma)

            call surfd(colat,x,vi,dvdti,dvdxi,d2vdti,d2vdxdti,d2vdxi,
     +       ij,nhts,clatt,xg,vgi,ij,vgti,sigma)

	    call surfd(colat,x,wi,dwdti,dwdxi,d2wdti,d2wdxdti,d2wdxi,
     +       ij,nhts,clatt,xg,wgi,ij,wgti,sigma)

            dvdt=cmplx(dvdtr,dvdti)
            dwdz=f1*(cmplx(dwdxr,dwdxi)+.5*w)

            t0=p3(j)
            u0=p1(j)
            R=p14(j)
            rho0=p8(j)
            drhodz=p9(j)
            drhodt=dp8(j)
C TEST ******************* Uncomment for No Mean winds or Latitude Gradients:
c	     u0=0.0
c	     drhodt=0.0
C TEST **********************************************************************
             h=p5(j)
             part2=fn*u0/(wr*snclat(j))

C Ask about this test--obsolete?

             test=.1*freq
             sigg=freq+part2

             paren=(xi*fn*u/(wr*snclat(j))+
     +             dvdt/wr+(ctnclat(j)+drhodt)*v/wr+
     +             drhodz*w+dwdz)
             geopot(j)=h*t/t0+(xi*h/sigg)*paren

C               if(part2.le.test) then
C                  geopot(j)=h*t/t0
C                 print *, 'z=',z,' colat=',colat,' part2=',
C     +                      part2,' test=',test
C                endif

C Hough decomposition (idecomp=1)

             if(idecomp.eq.1) then

C Decompose geopotential results:

               XR=REAL(geopot(j))*exp(x/2.)*snclat(j)
               XXXI=AIMAG(geopot(j))*exp(x/2.)*snclat(j)
               TR(1,j)=XR*(H11(j)*II11+H22(j)*II22+H33(j)*II33)
               TI(1,j)=XXXI*(H11(j)*II11+H22(j)*II22+H33(j)*II33)
               TR(2,j)=XR*(H12(j)*II12+H23(j)*II23+H34(j)*II34)
               TI(2,j)=XXXI*(H12(j)*II12+H23(j)*II23+H34(j)*II34)
               TR(3,j)=XR*(H13(j)*II13+H24(j)*II24+H35(j)*II35)
               TI(3,j)=XXXI*(H13(j)*II13+H24(j)*II24+H35(j)*II35)
               TR(4,j)=XR*(H1M2(j)*II1M2+H25(j)*II25+H36(j)*II36)
               TI(4,j)=XXXI*(H1M2(j)*II1M2+H25(j)*II25+H36(j)*II36)
               TR(5,j)=XR*(H1M4(j)*II1M4+H26(j)*II26)
               TI(5,j)=XXXI*(H1M4(j)*II1M4+H26(j)*II26)
               TR(6,j)=XR*H1M1(j)*II1M1
               TI(6,j)=XXXI*H1M1(j)*II1M1

             endif    !idecomp=1

C Density and Pressure perturbation (idenpress=1)

             if(idenpress.eq.1) then

               delrho(j)=(xi*rho0/sigg)*paren
               delpress(j)=delrho(j)*g*h+t*rho0*R

               delrho(j)=exp(x/2.)*delrho(j)
               delpress(j)=exp(x/2.)*delpress(j)

             endif    !idenpress=1

             xmod=xmod+cabs(geopot(j))
             xmodsqr=xmodsqr+(cabs(geopot(j))**2)
             geopot(j)=exp(x/2.)*geopot(j)
c
   96    continue  !number of latitudes loop

C Hough Mode Decomposition (idecomp=1)

         if(idecomp.eq.1)then

           DO 14 ii=1,6
              DO 13 J=1,IJ
                 XXR(J)=TR(ii,J)
                 XXI(J)=TI(ii,J)

 13	      continue

              DTHE=DTHETA
              CALL QSF(DTHE,XXR,FR,IJ)
              CALL QSF(DTHE,XXI,FI,IJ)

C Integral value at upper limit of latitudinal array:

              FX(ii)=CMPLX(FR(IJ),FI(IJ))

   14      continue

         endif   !idecomp=1

C  write (7,)RESULTS IN STEPS OF PRNTDZ  (Given in DATA statement--main_2.f)
c     
         IPRNT=1
         TESTZ=Z-Z1
         IF (TESTZ .LT. PRNTDZ) IPRNT=0
c  ONLY write(7,) IF Z IS LESS THAN 145 KM
         IF(Z.GT.145.) IPRNT=0

         IF(i.eq.1) go to 7

C Redo this IF test

         IF(IPRNT) 6,6,7
    7    z1=z
         write (7,100) z
c
         if(mig.eq.1)then
            CALL amphz(geopot,IJ)
c
         elseif(mig.eq.0)then
            call amphzgeo(geopot,IJ)
         endif

C Hough Mode Decomposition (idecomp=1)

         if(idecomp.eq.1) then

           DO 1 ii=1,IJ

              IF(NZONAL.lt.2)then 
                GPS(1)=FX(1)*H11(ii)+FX(2)*H12(ii)+FX(3)*H13(ii)
     +          +FX(4)*H1M2(ii)+FX(5)*H1M4(ii)+FX(6)*H1M1(ii)
              elseif(nzonal.eq.2)then
                GPS(1)=FX(1)*H22(ii)+FX(2)*H23(ii)+FX(3)*H24(ii)
     +          +FX(4)*H25(ii)+FX(5)*H26(ii)
              else   !nzonal>2
                GPS(1)=FX(1)*H33(ii)+FX(2)*H34(ii)+FX(3)*H35(ii)+
     +          FX(4)*H36(ii)
              endif

              if(mig.eq.1)then
                 CALL AMPHZ(GPS,1)
              endif

              GPT(ii)=GPS(1)

    1      continue

           if(mig.eq.1)then
              CALL AMPHZ(FX,6)
           endif
           IF(NZONAL.lt.2)then
              write(7,500) (FX(ii),ii=1,6)
           elseif(nzonal.eq.2)then
              write(7,600) (FX(ii),ii=1,5)
           else  !nzonal>2
              write(7,700) (FX(ii),ii=1,4)
           endif

C Output for idecomp alone:

           if(idenpress.eq.0)then
              do 95 j=1,ij
                 VLAT = XLAT(j)*RADEG
                 write (7,200) VLAT,geopot(j),GPT(j)
   95         continue
           endif

         endif  !idecomp

C Density and pressure perturbations (idenpress=1)

         if(idenpress.eq.1)then

            if(mig.eq.1)then
               call amphz(delrho,ij)
               call amphz(delpress,ij)
c
            elseif(mig.eq.0)then
               call amphzgeo(delrho,ij)
               call amphzgeo(delpress,ij)
            endif
c

C Output for idenpress alone:

           if(idecomp.eq.0)then
              do 2001 j=1,ij
                 VLAT = XLAT(j)*RADEG
                 write (7,201) VLAT,geopot(j),delrho(j),delpress(j)
 2001	      continue
           endif

         endif   !idenpress=1

C Output for both calls

         if(idenpress.eq.1.and.idecomp.eq.1)then
           do 2003 j=1,ij
              VLAT = XLAT(j)*RADEG
              write (7,203) VLAT,geopot(j),GPT(j),delrho(j),delpress(j)
 2003      continue
         endif

C Take out the GOTO6 statements

    6    continue

   97 continue  !number of heights loop

  908 format(1x,7hperiod=,f5.3,1x,5hXMOD=,e10.3,1x,8hXMODSQR=,e10.3)
  909 format(1x,5hXMOD=,e10.3,1x,8hXMODSQR=,e10.3)
      print 908, period,xmod,xmodsqr
      write (7,909) xmod,xmodsqr

	return
	end
