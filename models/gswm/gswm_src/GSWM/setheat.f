	subroutine setheat(mois)
C_________________________________________________________________________

C  This subroutine reads in arrays of O3 from CIRA (1991) provided by Keating
C            17 latitudes: -80 to +80 in steps of 10 degrees
C    27 pressure surfaces: .003,.005,.007,.010,.015,.020,.030,.050,.070,.100,
C                          .150,.200,.300,.500,.700,1.00,1.50,2.00,3.00,5.00,
C                          7.00,10.0,15.0,20.0,30.0,50.0 and 70.0 mbar
C
C  The densities are stored in the array dO3(17,27).  Subsequently, the elements
C  of dO3 are converted from units of ppmv to number densities.
C    N.B., [O3] = 2.079E+19*RHO(kg m-3)*O3(ppmv)     
C
C  The array of altitudes, zO3(17,27) is calculated using the output of
C  setatmos.  This array defines the altitudes which correspond to each
C  element of dO3(17,27) 

C_________________________________________________________________________

      dimension zx1(26),zxm(26),zy1(37),zyn(37),temp(89)
      dimension zx12(101),zxm2(101),tempb(239)

      character*10 ochoice1(12),ochoice2(12),ochoice3(12),ochoice4(12)
      character*10 ofile

      integer x1O3(17),d1O3(37,27),z1O3(17,27),lO3(37)

      real xO3(17),zO3(17,27),dO3(17,27),z2(27),d2O3(27),d3O3(17)
      real diO3(17,27)

      real pl(27),s1(4),s2(4) 
      real mzph(17,26),mzr(17,26),m1ph(26),m1r(26)      

      common/interp1/zpt(37,101,3),zpu(37,101,3),zpr(37,101,3),
     +zpxh(37,101,3),dO2(37,101,3)

C N.B., zph is the array of pressures; zr is log mass density
C defined in SR  SETATMOS on the latitude grid x and altitude grid y

      common/interp2/x(37),y(101),zt(37,101),zu(37,101),zr(37,101),
     +zxh(37,101),zph(37,101),O2(37,101),sigma

        common/htO3/alt(26),xlO(37),O3(37,26),tO3(37,26),
     +  zzO3(37,26,3),tzO3(37,26,3),sigma2,tO2(37,101),tzO2(37,101,3)

	data rs,bc,poo,go/8.31432,1.38054e-23,1.01325e+05,9.80665/
	data re,omega/6356766.,.72722e-04/

      data ochoice1/'ciraO3_jan','ciraO3_feb','ciraO3_mar',
     +             'ciraO3_apr','ciraO3_may','ciraO3_jun',
     +             'ciraO3_jul','ciraO3_aug','ciraO3_sep',
     +             'ciraO3_oct','ciraO3_nov','ciraO3_dec'/

c  hard-wired to use CIRA O3

        ofile=ochoice1(mois)
cc
C TEST PRINT 11/02
      print*,'SETHEAT: ofile= ',ofile
c
      pi=acos(-1.)
      pio2=pi/2.
      dtor=pi/180.
      rtod=180./pi
      iyd=15+30*(mois-1)

C  Read in array of latitudes defining CIRA91-O3 grid: -80 to +80 in 10 deg
C  Convert the values to colatitude
C  Write another array in reverse order for print (consistent with SETATMOS)

	OPEN(UNIT=23,file=ofile,form='formatted')

C201	format(8x,17(3x,i3),/)
201	format(8x,17(3x,i3),/)
cc202	format(1x,f6.3,2x,17(1x,f4.2,1x))
202	format(1x,f6.3,1x,17(1x,f5.2))
203	format(/)

100	format(f6.3,17i4)
101	format(/,1x,'CIRA-91 Ozone Density (log*10 m-3) vs pressure',/)
102	format(/,1x,'Altitudes (km) of CIRA-91 Ozone Density Values',/)
111     format(2x,4hlat=,17i4)

300	format(20i4)
301	format(/,1x,'Ozone Density (log*10 m-3) vs altitude',/)
401	format(/,1x,'Ozone Column Density (log*10 m-3) vs altitude',/)
311     format(4hlat=,19i4)

500	format(1x,6e11.3)

	read(23,203)
	read(23,201) (x1O3(i),i=1,17)
C TEST PRINT 11/02
c       print*,'Ozone Density Input:'
c	print '(8x,17(3x,i3),/)', (x1O3(i),i=1,17)


C  Reverse the order of  the array (+80 - -80 - consistent setatmos)
C  *******need colat to be increasing for SURF1/SURFD********
 	do 2 j=1,17
 	il=18-j
 	xO3(il)=float(x1O3(j))
    2  continue

C convert to colatitude:
	do 22 j=1,17
        xO3(j)=pio2-xO3(j)*dtor    
   22  continue

C  Read in pressure level and CIRA91-O3 (ppmv)
C  Reverse the order of  the array (ground up - consistent setatmos)

        j=28    
	do 3 jj=1,27    
        j=j-1  
	read(23,202) pl(j),(diO3(i,j),i=1,17)
C TEST PRINT 11/02
c	print '(1x,f6.3,1x,17(1x,f5.2))', pl(j),(diO3(i,j),i=1,17)
    3  continue

C  Reverse the order of  the array (+80 - -80 - consistent setatmos)
C  *******need colat to be increasing for SURF1/SURFD********
 	do 33 j=1,27
        ii=18
 	do 33 i=1,17
 	ii=ii-1
 	dO3(i,j)=diO3(ii,j)
   33   continue

	CLOSE(UNIT=23)

C  Calculate altitude of pressure levels for each element of dO3 using    
C  background atmosphere (setatmos); store in zO3; construct a  similar
C  grid of mass density to convert O3(ppmv) to [O3]

	ij=0
        do 5 j=1,26  
	ij=ij+1
	do 5 i=3,35,2
	il=(i-1)/2
	mzph(il,ij)=alog(zph(i,j))
	mzr(il,ij)=zr(i,j)
    5  continue
 
        do 7 il=1,17

	do 6 k=1,26 
	alt(k)=float(k-1)*4.
	m1ph(k)=mzph(il,k)
        m1r(k)=exp(mzr(il,k))
    6  continue
 
	do 8 ip=1,27
	xp=alog(pl(ip))
	call atsm(xp,m1ph,alt,26,1,s1,s2,4)
	call ali(xp,s1,s2,zp,4,1.0E-03,ier)
	zO3(il,ip)=zp
	call atsm(xp,m1ph,m1r,26,1,s1,s2,4)
	call ali(xp,s1,s2,xr,4,1.0E-03,ier)

C Convert the units of O3:

	fact=2.079E+19*xr
	dO3(il,ip)=fact*dO3(il,ip)

C *********************************************************************
C Arbitrarily increase the [O3] to test the importance of including
C  diurnal variation in [O3] in heating--user entered preference
C *********************************************************************

	if(diurno3.eq.1)then      !0 is standard

C Suggested by Keating et al. (1990) Figure 9
	   if(zO3(il,ip).le.53.) arg=1.
	   if(zO3(il,ip).gt.53..and.zO3(il,ip).le.70.) arg=.5*(1.+
     +          ((1./85.)*(zO3(il,ip)-53.)**2.))
	   if(zO3(il,ip).gt.70.) arg=.5*4.4
	   dO3(il,ip)=arg*dO3(il,ip)

	elseif(diurno3.eq.2)then  !0 is standard

C Suggested by Bjarnson et al. (1987) Figures 13 and 16
	   arg=1.-.1*exp(-((zO3(il,ip)-50.)/5.)**2)
     +	        +.1*exp(-((zO3(il,ip)-60.)/5.)**2) 
     +	        -.3*exp(-((zO3(il,ip)-75.)/7.5)**2)
     +	        +.4*exp(-((zO3(il,ip)-85.)/5.)**2)
     +	        -.05*exp(-((zO3(il,ip)-95.)/2.5)**2)
	   dO3(il,ip)=arg*dO3(il,ip)

	endif                     !else o3diurn=0 and standard run

C *********************************************************************

    8  continue
 
    7  continue

C Write an array of log[O3]

	do 4 j=1,17
	k=28
	do 4 kk=1,27
	k=k-1
        if(dO3(j,kk).ne.0.0) go to 44
	d1O3(j,k)=dO3(j,kk)
	go to 45
   44   d1O3(j,k)=alog10(dO3(j,kk))*10.
   45   z1O3(j,k)=zO3(j,kk)
    4  continue

C TEST PRINT 11/02
        write(10, 101)
 	write(10,111) (x1O3(i),i=1,17)

	j=28
	do 9 i=1,27
	j=j-1

C TEST PRINT 11/02
        write(10, 100) pl(j),(d1O3(k,i),k=1,17)

    9  continue

C TEST PRINT 11/02
        write(10, 102)
 	write(10,111) (x1O3(i),i=1,17)

	j=28
	do 10 i=1,27
	j=j-1

C TEST PRINT 11/02
 	write(10, 100) pl(j),(z1O3(k,i),k=1,17)

   10  continue

c  Build arrays of [O3] and total column O3 in the 5-deg/4km model grid
c        from 0. to 100.km only

C  preset arrays to 0.
      do 99 j=1,26
      do 99 i=1,37
      O3(i,j)=0.0
      tO3(i,j)=0.0
   99 continue

c  Start by defining the array in the 20-80km range specified by CIRA-91

      do 11 i=1,17
      ii=2*i+1

      do 12 j=1,27
      z2(j)=zO3(i,j)
      d2O3(j)=alog10(dO3(i,j))
   12 continue

      do 13 j=6,21
      z=float(j-1)*4.
      call atsm(z,z2,d2O3,27,1,s1,s2,4)
      call ali(z,s1,s2,xd,4,1.0e-03,ier)
      O3(ii,j)=xd
   13 continue

c  Below 20km make the densities drop off exponentially

      do 113 j=1,5
      z=float(j-1)*4.
      	O320=10.**(O3(ii,6))
	drop=exp(-((z-20.)/5.)**2.)
 	xd=drop*O320
      O3(ii,j)=alog10(xd)
  113 continue


   11 continue


C Write intermediate array of log[O3] in reverse altitude order
C        (consistent setatmos)

	do 4444 j=1,37
	lO3(j)=95-5*j
	xlO(j)=(90.-float(lO3(j)))*dtor
	k=27
	do 4444 kk=1,26
	k=k-1
        if(O3(j,kk).ne.0.0) go to 447
	d1O3(j,k)=O3(j,kk)
	go to 4444
  447   d1O3(j,k)=(O3(j,kk))*10.
 4444  continue

C TEST PRINT 11/02
        write(10, 301)
 	write(10, 311) (lO3(i),i=1,37)

	j=26
	do 449 i=1,26
	j=j-1
        iz=j*4
C TEST PRINT 11/02
        write(10, 300) iz,(d1O3(k,i),k=1,37)
  449  continue

c  Fill in the missing grid points in the 0-80km range

      do 14 j=1,21

      do 15 i=1,17
      ii=2*i+1
      d3O3(i)=O3(ii,j)
   15 continue

C  Take care of -75 to +75 (interpolate)

      do 16 i=1,16
      ii=2*i+2
      xl=xO3(i)+5.*dtor
      xp=xl*rtod
      call atsm(xl,xO3,d3O3,17,1,s1,s2,4)
      call ali(xl,s1,s2,xd,4,1.0e-03,ier)
      O3(ii,j)=xd
   16 continue

C  Extend values at +/-80 to the poles

      O3(1,j)=O3(3,j)
      O3(2,j)=O3(3,j)
      O3(36,j)=O3(35,j)
      O3(37,j)=O3(35,j)

   14 continue

C Set second=0. to turn off 2ndary O3 peak:
C Add a Secondary Peak (ppmv constant with latitude) between 82 and 100 km
C       after Bjarnson et al., JGR, 1987.

	do 55 i=1,37
	   do 55 j=22,26
	      z=float(j-1)*4.
		 second=.8+exp(-((z-88.)/5.)**2.)
	      fact=2.079e+19*exp(zr(i,j))
	      O3(i,j)=(fact*second)+10.**(O3(i,j))
	      O3(i,j)=alog10(O3(i,j))
 55	continue


C Build Array of Total Ozone in a Column at every Latitude
C    NB, the total amount above the corresponding altitude
C    NB, O3 array is log[O3] must convert to [O3] before totalling column

 	do 67 i=1,37
	tO3(i,26)=10.**O3(i,26)
 	do 66 j=1,25
 	tO3(i,26-j)=tO3(i,26-j+1)+(10.**O3(i,26-j))
   66   continue

C Convert column density to log[tO3]

 	do 77 j=1,26
 	tO3(i,j)=alog10(tO3(i,j))
   77   continue
   67   continue

C Write an array of log[O3] in reverse altitude order for printout
C        (consistent setatmos)

	do 544 j=1,37
	lO3(j)=95-5*j
	k=27
	do 544 kk=1,26
	k=k-1
        if(O3(j,kk).ne.0.0) go to 57
	d1O3(j,k)=O3(j,kk)
	go to 544
   57   d1O3(j,k)=(O3(j,kk))*10.
  544  continue

C TEST PRINT 11/02
        write(10, 301)
 	write(10,311) (lO3(i),i=1,37)

	j=26
	do 59 i=1,26
	j=j-1
        iz=j*4
C TEST PRINT 11/02
        write(10, 300) iz,(d1O3(k,i),k=1,37)
   59  continue

C Write an array of log column [O3] in reverse altitude order for printout

	do 444 j=1,37
	lO3(j)=95-5*j
	k=27
	do 444 kk=1,26
	k=k-1
        if(tO3(j,kk).ne.0.0) go to 47
	d1O3(j,k)=tO3(j,kk)
	go to 444
   47   d1O3(j,k)=(tO3(j,kk))*10.
  444  continue

C TEST PRINT 11/02
        write(10, 401)
 	write(10, 311) (lO3(i),i=1,37)

	j=26
	do 49 i=1,26
	j=j-1
        iz=j*4
C TEST PRINT 11/02
        write(10, 300) iz,(d1O3(k,i),k=1,37)
   49  continue

C Build Array of Total O2 in a Column at every Latitude
C    NB, the total amount above the corresponding altitude
C    NB, O2 array is log[O2] must convert to [O2] before totalling column

 	do 68 i=1,37
	if(O2(i,101).gt.0.) then
	tO2(i,101)=10.**O2(i,101)
	else
	tO2(i,101)=0.0
	endif
 	do 86 j=1,100
	if(O2(i,101-j).gt.0.) then
 	tO2(i,101-j)=tO2(i,101-j+1)+(10.**O2(i,101-j))
	else
	tO2(i,101-j)=tO2(i,101-j+1)
	endif
   86   continue

C Convert column density to log[tO2]

 	do 87 j=1,101
	if(tO2(i,j).gt.0.) then
 	tO2(i,j)=alog10(tO2(i,j))
	else
	tO2(i,j)=0.0
	endif
   87   continue
   68   continue

C___________________________________________________________________________
C  Set up arrays for later interpolation in heat92:
C

        sigma2=1.

	call surf1(37,26,xlO,alt,O3,37,zx1,zxm,zy1,zyn,zxy11,zxym1,
     +  zxy1n,zxymn,255,zzO3,temp(1),sigma2,ierr)

	call surf1(37,26,xlO,alt,tO3,37,zx1,zxm,zy1,zyn,zxy11,zxym1,
     +  zxy1n,zxymn,255,tzO3,temp(1),sigma2,ierr)

	call surf1(37,101,x,y,tO2,37,zx12,zxm2,zy1,zyn,zxy11,zxym1,
     +  zxy1n,zxymn,255,tzO2,tempb(1),sigma,ierr)

      return
      end

 
