	subroutine seth2o(mois,nzonal)
C_________________________________________________________________________

C  This subroutine reads in arrays of H2O forcing provided by h2osub.f
C	which is based upon the Groves (JATP, 1982) parameterization.
C            35 latitudes: -85 to +85 in steps of 5 degrees
C    41 (irregularly spaced) altitudes: X=0--->10. in steps of dx=1.
C	Units are: J/m3/s; converted from J/kg/s in ~testlat/testh2o.f
C  		arrays used herein are subsequently created in
C			~testlat/plotout/h2onointerp.pro
C
C  The heating rates (J/m3/s) are stored in the array dh2o(37,41).
C  and the corresponding array of altitudes in zh2o(37,41)
C   N.B. Values at the poles are identical to those at +/-85
c
C Modified to include cpability to read in heating rates for additional
C	months obtained by interpolating Groves seasonal averages
C						M. Hagan (6/24/99)(
C_________________________________________________________________________

        dimension zx1(41),zxm(41),zy1(37),zyn(37),temp(119)

	character*13 dchoice(12)
	character*12 schoice(12)
	character*13 dmfile

        real xh2o(41), xlih2o(35), dih2o(35,41)

        common/hth2o/zh2o(41),xlh2o(37),dh2o(37,41),zdh2o(37,41,3),
     +  sigma2,dah2o(37,41),zadh2o(37,41,3),sigma2a


        data dchoice/'h2o_jan.diurn','h2o_feb.diurn','h2o_mar.diurn',
     +             'h2o_apr.diurn','h2o_may.diurn','h2o_jun.diurn',
     +             'h2o_jul.diurn','h2o_aug.diurn','h2o_sep.diurn',
     +             'h2o_oct.diurn','h2o_nov.diurn','h2o_dec.diurn'/

        data schoice/'h2o_jan.semi','h2o_feb.semi','h2o_mar.semi',
     +             'h2o_apr.semi','h2o_may.semi','h2o_jun.semi',
     +             'h2o_jul.semi','h2o_aug.semi','h2o_sep.semi',
     +             'h2o_oct.semi','h2o_nov.semi','h2o_dec.semi'/

       icount=mois
cc       print * ,"In SR SETH2O icount=",icount
cc       print * ," and mois=",mois

       if(nzonal.eq.1) dmfile=dchoice(icount)
       if(nzonal.eq.2) dmfile=schoice(icount)

	pi=acos(-1.)
	pio2=pi/2.
        dtor=pi/180.

C  Read in array of latitudes defining H2O grid: -85 to +85 in 5 deg
C  Convert the values to colatitude
C  Write another array in reverse order for print (consistent with SETATMOS)

	OPEN(UNIT=53,file=dmfile,form='formatted')

201	format(5(5X ,7(f10.5),/))
222	format(2(1x,f5.2))
202	format(5(5X ,7(f10.5),/))
203	format(/)

	read(53,203)
	read(53,201) (xlih2o(i),i=1,35)


C  Reverse the order of  the array (+80 - -80 - consistent setatmos)
C  *******need colat to be increasing for SURF1/SURFD********
 	do 2 j=1,35
 	il=37-j
 	xlh2o(il)=xlih2o(j)
    2  continue
	xlh2o(1)=90.
	xlh2o(37)=-90.
cc	print '(37(1x,f4.0))',(xlh2o(i),i=1,37)

C convert to colatitude:
	do 22 j=1,37
        xlh2o(j)=pio2-xlh2o(j)*dtor    
   22  continue

C  Read in x altitude and corresponding array of h2o forcing at 37 latitudes:

	do 3 j=1,41    
	read(53,222) xh2o(j),zh2o(j)
	read(53,202) (dih2o(i,j),i=1,35)
    3  continue

C  Reverse the order of  the array (+90 - -90 - consistent setatmos)
C  *******need colat to be increasing for SURF1/SURFD********
 	do 34 j=1,41
 	do 33 i=1,35
 	ii=37-i
 	dh2o(ii,j)=dih2o(i,j)
   33   continue
	dh2o(1,j)=dh2o(2,j)
	dh2o(37,j)=dh2o(36,j)
   34   continue
cc	print *,(dh2o(i,5),i=1,37)

	CLOSE(UNIT=53)


C___________________________________________________________________________
C  Set up arrays for later interpolation in heat92:
C

        sigma2=1.

	call surf1(37,41,xlh2o,zh2o,dh2o,37,zx1,zxm,zy1,zyn,zxy11,zxym1,
     +  zxy1n,zxymn,255,zdh2o,temp(1),sigma2,ierr)

      return
      end

 
