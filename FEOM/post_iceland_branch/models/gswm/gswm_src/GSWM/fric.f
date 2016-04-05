	subroutine fric(z,vlat,tnc,rf,rfu,rfv)
      COMMON /MODEi/NZONAL,MOIS,NSS
      COMMON /MODEr/PERIOD,FREQ,FLUX

C  defines newtonian cooling and rayleigh friction terms to
C  parameterize gravity wave stress effects.  Coefficients
C  are defined initially below
C  in units of fraction of wave frequency.  Correct units are implemented
C  before exiting.

C  Logic to include rayleigh friction for non-migrating diurnal tide
C	with maximum valuies at the Equator		M. Hagan 1/11/00

	real relfric(12)
	data relfric/1.,.7,.4,.1,.4,.7,1.,.7,.4,.1,.4,.7/

	pi=acos(-1.)

	deglat=vlat*180./pi

	tnc=0.0
	rf=0.0
	rfu=0.0
	rfv=0.0
C
C Explicitly set GW Stress effects to ZERO outside MLT (M. Hagan 2/98):

        if (z.lt.70..or.z.gt.110.) go to 9

C_______________________  Vial (JGR, 8955-8969, 1986)  ______________________
C_______________________  ***works for Semidiurnal***  ______________________
c
c Units are s-1:
c	 add=2.182e-05
c	 if(z.lt.90.) then add=2.182e-05*exp(-((z-90.)/7.2)**2)
c	 if(z.ge.90..and.z.le.97.) then add=2.182e-05
c	 if(z.gt.97.) then add=2.182e-05*exp(-((z-97.)/7.3)**2)
c	tnc = tnc + add
c	rfu = rfu + add
c	rfv = rfv + add
C Plus molecular diffusion effects (Note: 7-km scale height approx):
c	addm=5.28e-13*exp(z/7.)
c	tnc = tnc + addm
c	rfu = rfu + addm
c	rfv = rfv + addm
c	go to 9
C____________________________  critical line  ___________________________
c
c	if(deglat.ge.0.) add=0.0
c	if(deglat.lt.0.) add=-2.0*sin(vlat)*exp(-((z-30.)/10.)**2)

c       add=0.0
c	tnc=tnc+add
c	rf=rf+add

C____________________________ gravity waves ______________________________

C       in units of /day:

c For migrating diurnal tide:
 	if (nzonal.eq.1) add = 1.5*exp(-((z-90.)/10.)**2.)*sin(2.*vlat)**2.

C       convert units:

 	add = add*(1./(2.*pi))

 	rfu = rfu + add
C RFV=0:
 	rfv = rfv + add
C TEST Add SAO to GW Stress effects--"relfric" defined above (M. Hagan 6/98):

        rfu=relfric(mois)*rfu
        rfv=relfric(mois)*rfv
C_________________________________________________________________________


	tnc=tnc*freq
	rf=rf*freq
	rfu=rfu*freq
	rfv=rfv*freq

   9  continue
	return
	end
