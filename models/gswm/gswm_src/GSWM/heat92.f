	subroutine heat92(nzonal,mois,z,hj,ij)
c  This subroutine calculates Ozone+O2 Heating using Strobel Parameterization
C	and obtains H2O Heating from tables read in SR seth2o

      	real hj(91),totja(91,24),xja(91,24),pO3(91),ptO3(91)
      	real pO2(91),ptO2(91)
      	real cja(91,24),totc(91,24),haja(91,24),totha(91,24)
      	real huja(91,24),tothu(91,24),heja(91,24),tothe(91,24)
      	real srcja(91,24),totsrc(91,24),srbja(91,24),totsrb(91,24)
	real H,Hv,dum(25),AA(3),BB(3)
	
	common/interp1/zpt(37,101,3),zpu(37,101,3),zpr(37,101,3),
     +	zpxh(37,101,3),dO2(37,101,3)

      	common/atm1/p1(91),p2(91),p3(91),p4(91),p5(91),p6(91),
     +	p7(91),p8(91),p9(91),p10(91),p11(91),p12(91),p13(91),p14(91)

c  TEMPorary common block to pass gammma-1 to heat92 to convert heating:
c               to units of: deg/day
	common/temp/gam1

C N.B., zph is the array of pressures; zr is log mass density
C defined in SR  SETATMOS on the latitude grid x and altitude grid y

	common/interp2/x(37),y(101),zt(37,101),zu(37,101),zr(37,101),
     +	zxh(37,101),zph(37,101),bO2(37,101),sigma

        common/hth2o/zh2o(41),xlh2o(37),dh2o(37,41),zdh2o(37,41,3),
     +  sigma3,dah2o(37,41),zadh2o(37,41,3),sigma2a

        common/htO3/zO(26),xlO(37),bO3(37,26),btO3(37,26),
     +  dO3(37,26,3),tdO3(37,26,3),sigma2,btO2(37,101),tdO2(37,101,3)

      	common/latvect/dtheta,clatt(91),xlat(91),snclat(91),csclat(91),
     +	tnclat(91),ctnclat(91),sin2i(91),gmlat(91),dip(91)

	data ro,go,poo/6356.766,9.80665,1.01325e+05/

	data fc,sigc/3.7e+5,2.85e-21/
	data fha,sigha/5013.,8.8e-18/
	data b1,b2,xm,xll,xls,sighu/85.,49.,.0127,3055.,2805.,.013/
	data fhe,sighe,sighe2/1200.,4.9e-18,6.6e-24/
C  These values are for the TOTAL heating (Strobel, 1978)
c	data blom2,bsom2,sigl,sigm,sigs/3.43,1.35,2.9e-19,1.54e-18,
c    +  	1.1e-17/
c	data a,b,srbmx/.143,9.64e+8,9.03e-19/
c	data effc,effha,effhe,effhu/1.,1.,1.,1./
C  These values are for the NET heating (Strobel, 1978)
 	data blom2,bsom2,sigl,sigm,sigs/.98,.43,2.9e-19,1.7e-18,
     +  	1.15e-17/
 	data a,b,srbmx/.67,3.44e+9,2.43e-19/
 	data effc,effha,effhe,effhu/.49186,.78828,.74,.80606/

   33	format (1x,'Alt=',f6.1,' Lat=',f4.0,' H=',f7.1,' [O3]=',
     +     e11.3,' [tO3]=',e11.3,' [O2]=',e11.3,' [tO2]=',e11.3)
   44	format (3x,'Chp=',e11.3,' Ha=',e11.3,' Hu=',e11.3,
     +     ' He=',e11.3,' Tot=',e11.3) 
   55	format (2x,'lst=',i3,' X=',f5.0,' zi=',f4.0,' v=',f6.1,
     +     ' Hv=',f7.1,' Ch=',f6.2,' xn3=',e11.3,' xn2=',e11.3)
100     format(1x,i3,6e11.3)
500	format('SR heat92:  Heating (W/m3); Called from SR surfd')
501     format(1x,f5.1,9f10.0,2(/,6x,10f10.0))
503     format(1x,f5.1,9e10.3,2(/,6x,10e10.3))
502     format(1x,f4.0,7e11.3)
111     format(6e11.3)
      
        do 9 i=1,ij
    	pO3(i)=0.0
    	ptO3(i)=0.0
    	pO2(i)=0.0
    	ptO2(i)=0.0
	do 9 ii=1,24
    	totc(i,ii)=0.0
    	totha(i,ii)=0.0
    	tothu(i,ii)=0.0
    	tothe(i,ii)=0.0
    	totsrc(i,ii)=0.0
    	totsrb(i,ii)=0.0
	xja(i,ii)=0.0
	heja(i,ii)=0.0
    	totja(i,ii)=0.0
    9   continue

C calculate o3 heat forcing:

	   if(z.ge.20..and.z.le.126.)then

	      pi=acos(-1.)
	      pio2=pi/2.
	      dtor=pi/180.
	      re=6378.165

C  Calculate solar declination angle as a function of day of year:

	      iyd=15+30*(mois-1)
	      
	      c1=23.5*dtor
	   
	      soldec=atan(tan(c1)*sin(2*pi*(iyd-80.)/365.))


C---------------------------------------------------------------------------
C  Interpolate O2; O3 and tO2; tO3 to values at gridpoints:
C
	      Do 20 i=1,ij
		 colat=clatt(i)
		 H=p5(i)
		 xx=1000.*(re+z)/H
		 plat=(pio2-colat)/dtor
		 
		 call surfd(colat,z,O3,dO3dt,dO3dz,d2O3dt,d2O3dzdt,
     +          d2O3dz,37,26,xlO,zO,bO3,37,dO3,sigma2)

		 call surfd(colat,z,tO3,dtO3dt,dtO3dz,d2tO3dt,
     +          d2tO3dzdt,d2tO3dz,37,26,xlO,zO,btO3,37,tdO3,sigma2)

		 call surfd(colat,z,O2,dO2dt,dO2dz,d2O2dt,d2O2dzdt,
     +          d2O2dz,37,101,x,y,bO2,37,dO2,sigma)

		 call surfd(colat,z,tO2,dtO2dt,dtO2dz,d2tO2dt,
     +          d2tO2dzdt,d2tO2dz,37,101,x,y,btO2,37,tdO2,sigma)

C Heating formulation is in cgs units....need to convert [O3] density

		 pO3(i)=1.e-6*10.**O3
		 ptO3(i)=1.e-6*10.**tO3

		 pO2(i)=10.**O2
		 ptO2(i)=10.**tO2

c	print 33,z,plat,H,pO3(i),ptO3(i),pO2(i),ptO2(i)

		 cscolat=csclat(i)
		 sncolat=snclat(i)

C---------------------------------------------------------------------------
C For every hour of the day, calculate zenith angle and heating:

		 do 20 ihour=1,24

		    tp=(ihour-1)*pi/12.-pi
		    xnu=cscolat*sin(soldec)+cos(soldec)*
     +                  sncolat*cos(tp)
		    zi=acos(xnu)
C
C  Need to calculate Grazing Altitude and corresponding scale height
C     for Zenith Angles larger than PI/2 

		    if(zi.lt.pio2) then
		       v=0.
		       Hv=1.
		    else
		       v=(re+z)*sin(pi-zi)-re
		       if(v.lt.0) then
			  Hv=1.
		       else
			  g=go*(ro/(ro+z))**2.

			  call surfd(colat,v,xh,dxdt,dxdz,d2xdt,
     +                    d2xdzdt,d2xdz,37,101,x,y,zxh,37,zpxh,sigma)

			  call surfd(colat,v,rholog,drhodt,drhodz,
     +                    d2dtdt,d2dzdt,d2dzdz,37,101,x,y,zr,37,
     +                    zpr,sigma)
     

c       print 100, i,z,colat,v,rholog,xh

			  denom=exp(rholog)
			  Hv=poo*exp(-xh)/g/denom
			  
		       endif    !if v lt 0
		    endif	!if zi lt pi02
C
C Calculate O3 Heating due to Chappuis, Hartley, and Huggins bands and
C  Herzberg Continuum and O2 Heating due to Herzberg and Schumann-Runge 
C  Continuum and Schumann-Runge Bands (revised; Strobel, JGR, 1978)
         
C
C  Need O3 and O2 density along the solar ray path (N.B., total gas
c   densities were calculated every 4 km in SR SETHEAT; need to multiply
C   these densities by dz=4e+5 cm.)

		    xn3=4.e+5*ptO3(i)*chp(xx,zi,v,Hv)

		    xn2=4.e+5*ptO2(i)*chp(xx,zi,v,Hv)

		    pzi=zi/dtor         
c	print 55, ihour, xx, pzi, v, Hv, chp(xx,zi,v,Hv), xn3, xn2

		    if (xn3.le.0.) then
		       cja(i,ihour)=0.0
		       haja(i,ihour)=0.0
		       huja(i,ihour)=0.0
		    else
C                
C  N.B. multiply for efficiency factors (net heating rates) to account for
C       chemical energy stored in the molecule
C              
C     CHAPPIUS BANDS    
C                      
		     cja(i,ihour)= effc*pO3(i)*fc*sigc*EXP(-xn3*sigc)
c
C Introduce 30% increase to account for reflected radiation (Lacis + Hansen)
c
		     cja(i,ihour)=cja(i,ihour)+.30*cja(i,ihour)
c
C                                                     
C     HUGGINS BANDS                                  
C                                                   
		       huja(i,ihour)=effhu*pO3(i)*(b1+(b2-b1)*
     +                  EXP(-xn3*sighu*exp(-xm*xll))
     +                  -b2*EXP(-xn3*sighu*exp(-xm*xls)))/xn3/xm
C                                                         
C     HARTLEY REGION                                     
C                                                       
		       haja(i,ihour)=effha*pO3(i)*fha*sigha*
     +                  EXP(-xn3*sigha)

		    endif	!if xn3 le 0

		    if (xn3.le.0..or.xn2.le.0.) then
		       heja(i,ihour)=0.0
		    else
C                                                      
C     HERTZBERG CONTINUUM                             
C                                                    
		       heja(i,ihour)= effhe*fhe*(pO2(i)*sighe2+
     +                  pO3(i)*sighe)*EXP(-xn2*sighe2-xn3*sighe)

		    endif	!if xn3 or xn2 le 0

		    if (xn2.le.0..or.z.lt.60.) then
		       srcja(i,ihour)=0.0
		       srbja(i,ihour)=0.0
		    else
C                                                      
C     SCHUMANN-RUNGE CONTINUUM                             
C                                                    
		     srcja(i,ihour)=pO2(i)*(blom2*exp(-sigl*xn2)+
     +                  (bsom2-blom2)*exp(-sigm*xn2)-bsom2*
     +                  exp(-sigs*xn2))/xn2

C	srcja(i,ihour)=0.0

C
C     SCHUMANN-RUNGE BANDS                             
C                                                    

		       if(xn2.ge.1.e+18) then
			srbja(i,ihour)=pO2(i)*(1./(a*xn2+b*sqrt(xn2)))
		       else
			  srbja(i,ihour)=srbmx*pO2(i)
		       endif

C	srbja(i,ihour)=0.0

		    endif	!if xn2 le 0 or z lt 60
	
C
C Total Ozone + O2 Heating
C
		    xja(i,ihour)=cja(i,ihour)+haja(i,ihour)+
     +              huja(i,ihour)+heja(i,ihour)+srcja(i,ihour)+
     +              srbja(i,ihour)

		    totja(i,ihour)=totja(i,ihour)+xja(i,ihour)
		    totc(i,ihour)=totc(i,ihour)+cja(i,ihour)
		    totha(i,ihour)=totha(i,ihour)+haja(i,ihour)
		    tothu(i,ihour)=tothu(i,ihour)+huja(i,ihour)
		    tothe(i,ihour)=tothe(i,ihour)+heja(i,ihour)
		    totsrc(i,ihour)=totsrc(i,ihour)+srcja(i,ihour)
		    totsrb(i,ihour)=totsrb(i,ihour)+srbja(i,ihour)
                                                
 20		 CONTINUE	!do ihour = 24 hours loop
C                                             
C---------------------------------------------------------------------------

		 do 1 i=1,ij

c Dummy what comes out for comparison with Strobel (1978) (units: deg/day)
c   N.B. erg/cm3/sec*1e+6cm3/m3*3600sec/hour*joule/1e+7erg*gam1/R = deg/hour
c        in an hour=1,24 loop
c  		p14(i)=g*H/T = R
c	fact=(p7(i)/p8(i))*360.*gam1/p14(i)
c Dummy what comes out for comparison with TIMEGCM results (units: erg/g/sec)
c   N.B. erg/cm3/sec*1/rho*1e+3 = erg/g/sec                           
c  		p8(i)=rho(kg/m3)*m3/1e+6cm3*1e+3g/kg=1e-3*rho(g/cm3)
c       do 1 i=1,ij
c 	fact=1e+3/p8(i)
C Convert to J/m3/s for use in abcr1_2.f (M. Hagan 6/28/94)
c   N.B. erg/cm3/sec*1.e+6cm3/m3*joule/1e+7erg= joule/m3/sec

		    fact=.1
		    do 11 ihour=1,24
 11		       dum(ihour)=fact*xja(i,ihour)
		       dum(25)=dum(1)
		       
		       call forit(dum,12,2,AA,BB,ier)

		       hj(i)=sqrt((AA(nzonal+1)**2.)+
     +                       (BB(nzonal+1)**2.))

c   	hj(i)=fact*totja(i)
c   	totc(i)=fact*totc(i)
c   	totha(i)=fact*totha(i)
c   	tothu(i)=fact*tothu(i)
c   	tothe(i)=fact*tothe(i)
c   	totsrc(i)=fact*totsrc(i)
c   	totsrb(i)=fact*totsrb(i)

 1		    continue	! do i=1,ij loop

	   endif		!if z between 20 and 126

C H2O IR Forcing 

	   if(z.le.20)then

C*******************************************************************
C	H2O Heating from tables
C*******************************************************************

	      Do 21 i=1,ij
		 colat=clatt(i)

		 call surfd(colat,z,h2o,dh2odt,dh2odz,d2h2odt,
     +                d2h2odzdt,d2h2odz,37,41,xlh2o,zh2o,dh2o,
     +                37,zdh2o,sigma3)


		 hj(i)=h2o

 21	      continue

	   endif   !if ze le 20 troposphere heating

C TEST PRINT 11/02
c	write(10,500)
c	write(10,501) z,(hj(i),i=1,ij)
c	write(10,503) z,(hj(i),i=1,ij)

	return
	end

      FUNCTION CHP(X,zi,v,Hv)
C                          
      PI=ACOS(-1.)
      PIO2=PI/2.
	dtor=pi/180.
	ziup=100.*dtor
      Re=6378.165
C                      
C For zi.gt.pi/2 Swider Approximation (PSS, 12, 761-782, 1964) is used
c                       
      if(zi.lt.pio2) then
      CHP=CHAPMAN(x,zi)
      else
      CHP=0.
      IF(zi.ge.ziup.or.v.LE.0.0) RETURN 
	Y=1000.*(Re+v)/Hv
      CHP=2.*CHAPMAN(Y,pio2)-CHAPMAN(X,PI-zi)*EXP(Y*(1.-1./SIN(zi)))
      ENDIF
C     
      RETURN
      END  
      FUNCTION CHAPMAN(X,zi)
C                          
C      CALCULATION OF CHAPMAN FUNCTION USING EXPRESSION PROVIDED BY
C          ROBLE (private communication) FROM SR SOLHEAT.F     
C                                                             
C                                                           
        pi=acos(-1.)
      Y=0.5*X*(COS(ZI))**2
      TY=SQRT(Y)
      sy=x*sin(zi)
      IF(TY.GT.8.) GO TO 10
      YERF=(1.0606963+0.55643831*TY)/(1.0619896+1.7245609*TY+TY*TY)
      GO TO 11
  10  YERF=0.56498823/(0.06651874+TY)
  11  CONTINUE
      CHAPMAN=SQRT(0.5*PI*X)*YERF                                                
C                                                        
      RETURN                                            
      END                                              
