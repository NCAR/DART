	subroutine seteddy_gs(mois)

C  This subroutine reads the appropriate Garcia and Solomon eddy
C  diffusion tables, and sets up an array for interpolation later
C  on, by subroutine eddy_gs .
C  Experiment to overwrite values with Khattatov's Kzz. NOTE: If 
C	iboris=1 that other dissipative calls should be commented.
C					M. Hagan (3/96)

	integer idecomp, idenpress
	integer wback

	dimension zx1(13),zxm(13),zy1(30),zyn(30),temp(73)
        dimension rfi(13,30),xi(13),rki(13,14)

	common /ritecom/wback

	common /interpe/zped(13,30,3),xe(13),ye(30),zedd(13,30),sigmae

C________________________________________________________________________________
C  Define horizontal x-array (colatitude, radians) corresponding to
C  Garcia/Solomon eddy coefficient tables.

	pi=acos(-1.)

C  Original created by JMF after FV ====> produces xe(14)=xe(15)
c	dx=7.2*pi/180.

c	xe(1)=0.
c	xe(15)=pi
c	if(mois.gt.6) then
c	do i=3,14
c	xe(i)=float((i-2)*2)*dx
c	end do
c	xe(2)=dx

c	else
c	do i=2,13
c	xe(i)=float((i-2)*2)*dx + dx
c	end do
c	xe(14)=pi
c	endif

C  Modification --M. Hagan (8/94) assumes equally spaced values at
C 	15-deg increments are read in; requires changes to dimensions
C		15====>13
 	dx=15.*pi/180.
 	xe(1)=0.
 	do i=2,13
 	xe(i)=xe(i-1)+dx
 	end do
 
C  __________Read in Garcia/Solomon values for appropriate month ________________

	OPEN(unit=23,file='revised_gs_50.tables',form='formatted',
     +        status='old')

	mi=mois
	if(mi.gt.6) mi=mi-6

	do m=1,mi

	do j=1,26
100     format(4x,i3,13(1x,f4.2))
	read(23,100) ixpo,(rfi(i,j),i=1,13)

	do i=1,13
C rfi arrays need to be inverted for Jan==>Jun, since
C   values read in are for _90 to +90 (No need to invert for months 7-12)
	rfi(i,j)=rfi(i,j)*10.**ixpo
	end do

	end do
	end do

	CLOSE(23)

C  Initialize & construct latitude array for printout

	do i=1,13
	xi(i)=90.-(xe(i)*180./pi)
	do j=1,30
	zedd(i,j)=0.0
	end do
	end do

C  Fill in zedd matrix

        do i=1,13
	L=i
C Inversion of arrays for Jan==>Jun:
        if(mois.le.6) L=14-i

	do j=1,26
	zedd(i,j+4)=rfi(L,j)
        end do

	end do

C  Extrapolate to ground

	do j=1,4
	do i=1,13
	zedd(i,j)=zedd(i,5)*exp(float(j-5)*4./8.)
	end do
	end do

C  Set up y-array (altitudes)

        do j=1,30
	z=float(j-1)*4.
	ye(j)=z

C  Garcia & Solomon (1985) assume momentum deposition/diffusion effects
C  	of breaking diurnal tide at 90km and above (after Lindzen, JGR,
C	9707, 1981) at the equator with exponential decay over 20-deg
C 	latitude and to 5km below --- M. Hagan (8/94)
C Specifically, Kzz=200*exp((z-90)/5)*exp(-lat/20)**2)
C 	Remove diurnal breaking results

C USING SMOOTHED REVISIONS TO TABLES ALREADY CORRECTED
	go to 999
	do i=1,13
	if(z.le.90.) then
	add=200.*exp((z-90.)/5.)*exp(-(xi(i)/20.)**2.)
	try=zedd(i,j)-add
	endif
	if(z.gt.90.) then
	add=200.*exp(-(xi(i)/20.)**2.)
	try=zedd(i,j)-add
	endif
	if(z.gt.70..and.try.lt.50.) then
	try=50.
	endif
	zedd(i,j)=try
	end do	
  999	continue

	end do	

300     format(1x,'Eddy Diffusivity after Garcia & Solomon (1985)',/)
200     format(1x,f4.0,13f5.1)
400     format(5x,13f5.1)
	if(wback.eq.1)then
	   write(10,300) 
	   write(10,400) (xi(i),i=1,13)
	   do j=1,30
	      write(10,200) ye(j),(zedd(i,j),i=1,13)
	   end do
	endif


C_____________________________________________________________________________


C  Pre-process arrays as in setatmos.f  :

  	sigmae=1.0
        call surf1(13,30,xe,ye,zedd,13,zx1,zxm,zy1,zyn,zxy11,zxym1,
     1  zxy1n,zxymn,255,zped,temp(1),sigmae,ierr)
 
        return
        end
