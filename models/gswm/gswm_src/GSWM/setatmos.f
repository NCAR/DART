	subroutine setatmos
c Modified to include HWM93 wind option		M. Hagan (10/20/98)
c
	real zx1(101),zxm(101),zy1(37),zyn(37),temp(239)
        real zxy11,zxym1,zxy1n,zxymn
c
	integer latgrad,polebc
	integer wback
c
	dimension d(8),tt(2),w(2),sw(25),po(37)
	dimension ih(19)
	dimension ap(1),it(19),iu(19),ifd(19),il(19)
C	dimension zx1(101),zxm(101),zy1(37),zyn(37),temp(239)
	dimension ubr(19)

	character*13  hchoice(12)
	character*13  hfile

      common/interp1/zpt(37,101,3),zpu(37,101,3),zpr(37,101,3),
     +zpxh(37,101,3),do2(37,101,3)
      common/interp2/x(37),y(101),zt(37,101),zu(37,101),zr(37,101),
     +zxh(37,101),zph(37,101),o2(37,101),sigma
      COMMON /MODEi/NZONAL,MOIS,NSS
      COMMON /MODEr/PERIOD,FREQ,FLUX

c  The following common block will contain control flags for the
c  background atmosphere.  Put it in the relevant subroutines

	common /backatm/latgrad,polebc
C
	common /ritecom/wback
C
	data ap/4./
	data (sw(i),i=1,25)/25*1./
	data rs,bc,poo,go/8.31432,1.38054e-23,1.01325e+05,9.80665/
	data re,omega/6356766.,.72722e-04/
C
	data hchoice/'hrdi0-128_jan','hrdi0-128_feb','hrdi0-128_mar',
     +             'hrdi0-128_apr','hrdi0-128_may','hrdi0-128_jun',
     +             'hrdi0-128_jul','hrdi0-128_aug','hrdi0-128_sep',
     +             'hrdi0-128_oct','hrdi0-128_nov','hrdi0-128_dec'/
C
C*************************** OPTIONS ************************************     

C  Hard-wired to
C            overwrite mean winds with interpolated/extrapolated HRDI data
C            0-128 km  (M. Hagan 3/17/04)

C*************************** OPTIONS ************************************     

c Filename with the hrdi information

	hfile=hchoice(mois)
c
	pi=acos(-1.)

	iyd=15+30*(mois-1)

C_________________________________________________________________________

C  This subroutine sets up arrays which ATMOS5 will later use for interpolation
C  and calculation of various partial derivatives.  The arrays set up here are:

C		     zt = zonal mean temp
C		     zu = mean zonal wind
C                    zr = log mass density
C 		     zxh = pressure height

C_________________________________________________________________________
c
c Do loop for all latitude: colat 0 to PI in radians

	do 2 i=1,37
	   colat=float(i-1)*5.*pi/180.
	   x(i)=colat
c
c Do loop for all altitude 0-400km

	   do 2 j=1,101
	      z=float(j-1)*4.
	      y(j)=z
	      
	      deglat=90.-float(i-1)*5.

	      sw(7)=0.
	      sw(8)=0.
	      sw(14)=0.
c
c New findings suggest using SW(10)=0. JKH 4/29/98

	      sw(10)=0.
	      call tselec(sw)
	      call wtselec(sw)

	      call gtd6(iyd,0.,z,deglat,285.,12.,flux,flux,ap,48,d,tt)

C Calculate # density and convert from cm-3 to m-3

	      xn=(d(1)+d(2)+d(3)+d(4)+d(5)+d(7))*1.0e+06
c
c Initialize pressure at ground for each latitude

	      if(j.gt.1) go to 20

	      po(i)=xn*bc*tt(2)
 20	      p=xn*bc*tt(2)
c
c Convert pressure from Pa to mbar for printout only

	      zph(i,j)=p/100.
	      zxh(i,j)=-alog(p/po(i))
c
c Calculate ln mass density and convert from g/cm3 to kg/m3

	      zr(i,j)=alog(d(6)*1.0e+03)
	      zt(i,j)=tt(2)
	      if(d(4).gt.0.) then
		 o2(i,j)=alog10(d(4))
	      else
		 o2(i,j)=0.0
	      endif
c
c End of both latitude and altitude loops

 2	   continue

c jkh 9/4/98 initialize zonal mean wind to zero to 400km

	do j=1,101
	   do i=1,37
	      zu(i,j)=0.0
	   enddo
	enddo

C_________________________________________________________________________

C  Read in HRDI winds between 
C              0 and 128 km (33 altitudes)
C  N.B. These background contain groves winds to 12 km

	      OPEN(UNIT=19,file=hfile,form='formatted',STATUS='old')
	      read(19,903)	!skip top two lines

		 jin=1
		 jend=33
c
 902	      format(4x,19i4)
 903	      format(/)

C HRDI data on GSWM background grid 
C 128km-33rd height:
C delta z=4km, delta lat=10 degrees (interpolated to std. 5 deg. later)
c
c Read in winds according to user's preference
c
	      do 906 j=jin,jend
		 read(19,902) (iu(i),i=1,19)
		 i=1
c
c Overwrite HRDI winds into background
c
		 do k=1,37,2	!latitude every 10 degrees
		    zu(k,j)=float(iu(i))
		    i=i+1
		 end do
 906	      continue
c
c Interpolate to the GSWM grid:  to 128 
c
	      do 904 j=jin,jend
		 do 905 i=2,36,2
		    zu(i,j)=(zu(i+1,j)+zu(i-1,j))*.5
 905		 continue
 904	      continue

	      CLOSE(UNIT=19)
	      
C_________________________________________________________________________


C   Points above 70 degrees; smooth transition to zero at pole

	do 51 j=1,28

	   xlatn=70.*pi/180.
	   do 8 i=1,4
	      xlat=(90.-float(i-1)*5.)*pi/180.
	      zu(i,j)=zu(5,j)*cos(xlat)/cos(xlatn)
 8	      zu(38-i,j)=zu(33,j)*cos(xlat)/cos(xlatn)

 51	   continue


C_________________________________________________________________________

C  Pre-process arrays which are needed later for interpolation.  The following
C  SURF1 subroutine is from "FITPACK", obtained from Dick Valent 
C  (valent@bierstadt.ucar.edu) in NCAR/SCD consulting office.  The routine
C  that later (in ATMOS5) uses the arrays from SURF1 is called SURFD.

	sigma=1.

	call surf1(37,101,x,y,zt,37,zx1,zxm,zy1,zyn,zxy11,zxym1,
     +  zxy1n,zxymn,255,zpt,temp(1),sigma,ierr)

        call surf1(37,101,x,y,zu,37,zx1,zxm,zy1,zyn,zxy11,zxym1,
     +  zxy1n,zxymn,255,zpu,temp(1),sigma,ierr)

	call surf1(37,101,x,y,zxh,37,zx1,zxm,zy1,zyn,zxy11,zxym1,
     +  zxy1n,zxymn,255,zpxh,temp(1),sigma,ierr)

	call surf1(37,101,x,y,zr,37,zx1,zxm,zy1,zyn,zxy11,zxym1,
     +  zxy1n,zxymn,255,zpr,temp(1),sigma,ierr)

	call surf1(37,101,x,y,o2,37,zx1,zxm,zy1,zyn,zxy11,zxym1,
     +  zxy1n,zxymn,255,do2,temp(1),sigma,ierr)

C_________________________________________________________________________

C  Print out array of zonal mean temperature

100	format(20i4)
111     format(4hlat=,19i4)
	ik=0
        do 1 k=1,37,2
	   ik=ik+1
	   il(ik)=(90.01-x(k)*180./pi)
1       continue

101	format(1x,'zonal mean temperature')
c
c Write to the background file if a new month

	if(wback.eq.1)then
 	       write(10, 101)
	       write(10, 111) (il(i),i=1,19)
	endif
c
c Modify loop to print out background up to 400km jkh 8/98

	ix=101
	do 30 i=1,101,2
	   iz=400-(i-1)*4
c	   ix=77
c	   do 30 i=1,77,2
c	      iz=304-(i-1)*4

	   kt=0
	   do 40 k=1,37,2
	      kt=kt+1
	      it(kt)=zt(k,ix)
 40	   continue

	   ix=ix-2
	   if(wback.eq.1)then
	      write(10, 100) iz,(it(k),k=1,19)
	   endif
 30	continue

C  Print out array of zonal mean winds

102	format(1x,'zonal mean winds')

	if(wback.eq.1)then
 	       write(10, 102)
	       write(10, 111) (il(i),i=1,19)
	endif
c
c Modify loop to print out background up to 400km jkh 8/98

	ix=101
	do 31 i=1,101,2
	   iz=400-(i-1)*4
c	   ix=77
c	   do 31 i=1,77,2
c	      iz=304-(i-1)*4

	   kt=0
	   do 41 k=1,37,2
	      kt=kt+1
	      iu(kt)=zu(k,ix)
 41	   continue

	   ix=ix-2
	   if(wback.eq.1)then
	      write(10, 100) iz,(iu(k),k=1,19)
	   endif
 31	continue

 133	format(1x,'pressure (mbar)')
c       write(10, 133)
c	write(10, 111) (il(i),i=1,19)

	ix=26
	do 34 i=1,26  
	   iz=100-(i-1)*4

	   kt=0
	   do 44 k=1,37,2
	      kt=kt+1
	      ih(kt)=zph(k,ix)
 44	   continue

	   ix=ix-1
c	write(10, 100) iz,(ih(k),k=1,19)

34	continue


C  Print out array of normalized intrinsic frequency

103	format(1x,'normalized intrinsic frequency')
c	if(wback.eq.1)then
c 	       write(10, 103)
c	       write(10, 111) (il(i),i=1,19)
c	endif
c
c	ix=77
c	do 32 i=1,77,2
c	iz=304-(i-1)*4
c
c	kt=0
c	do 42 k=1,37,2
c	kt=kt+1
c	sint=sin(x(k))
c	if(k.eq.1) sint=1.
c	if(k.eq.37) sint=1.
c	dfreq=freq+float(nzonal)*zu(k,ix)/(re*sint)
c	dfreq=dfreq/(2.*omega)
c	ifd(kt)=(dfreq*1000.)
c42 	continue
c
c	ix=ix-2
c	if(wback.eq.1)then
c 	write(10, 100) iz,(ifd(k),k=1,19)
c	endif
c32	continue


	return
	end

