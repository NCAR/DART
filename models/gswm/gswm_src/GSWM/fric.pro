;	subroutine fric(z,vlat,tnc,rf,rfu,rfv)
;	COMMON/MODE/NZONAL,PERIOD,FREQ,MOIS,NSS,FLUX

;  defines newtonian cooling and rayleigh friction terms to
;  parameterize gravity wave stress effects.  Coefficients
;  are defined initially below
;  in units of fraction of wave frequency.  Correct units are implemented
;  before exiting.

;	real relfric(12)
;	data relfric/1.,.7,.4,.1,.4,.7,1.,.7,.4,.1,.4,.7/
	relfric=[1.,.7,.4,.1,.4,.7,1.,.7,.4,.1,.4,.7]
	pi=acos(-1.)
	freq=2*pi/86400.
	mois=3
	rfs=fltarr(59,41)
	rfsn=fltarr(59,41)
	vlat=fltarr(59)
	slat=fltarr(59)
	slatn=fltarr(59)
	alt=fltarr(41)

	for i=1,59 do vlat(i-1)=(3*i-90.)*pi/180.
	for i=1,41 do alt(i-1)=i+69.

	deglat=vlat*180./pi

	for i=1,41 do begin
	z=alt(i-1)
	for j=1,59 do begin
	tnc=0.0
	rf=0.0
	rfu=0.0
	rfun=0.0
	rfv=0.0

; For migrating tide:
 	add = 1.5*exp(-((z-90.)/10.)^2.)*sin(2.*vlat(j-1))^2.
 	slat(j-1) = sin(2.*vlat(j-1))^2.
 	addn = 1.5*exp(-((z-90.)/10.)^2.)*cos(1.3*vlat(j-1))^2.
 	slatn(j-1) = cos(1.3*vlat(j-1))^2.

;       convert units:

 	add = add*(1./(2.*pi))
 	addn = addn*(1./(2.*pi))

 	rfu = rfu + add
 	rfun = rfun + addn
; RFV=0:
 	rfv = rfv + add
; TEST Add SAO to GW Stress effects--"relfric" defined above (M. Hagan 6/98):

        rfu=relfric(mois)*rfu
        rfun=relfric(mois)*rfun
        rfv=relfric(mois)*rfv
	rfs(j-1,i-1)=rfu
	rfsn(j-1,i-1)=rfun
;_________________________________________________________________________


	tnc=tnc*freq
	rf=rf*freq
	rfu=rfu*freq
	rfv=rfv*freq
	endfor
	endfor

	!p.multi=[0,2,2]
 	plot,deglat,slat, $
 	xrange=[-90,90],xticks=6,xminor=3,/xstyle, $
 	xtitle='!5Latitude', ytitle='sin(2.*lat)!u2!n', $
	title='Migrating !7m!5!deff!n Parameterization'

 	plot,deglat,slatn, $
 	xrange=[-90,90],xticks=6,xminor=3,/xstyle, $
 	xtitle='Latitude', ytitle='cos(1.3*lat)!u2!n', $
	title='Non-Migrating !7m!5!deff!n Parameterization'

 	contour,rfs,deglat,alt, $
 	xrange=[-90,90],xticks=6,xminor=3,/xstyle, $
 	yrange=[70,110],yticks=4,yminor=5,/ystyle, $
 	xtitle='Latitude',ytitle='Altitude (km)', $
 	title='GSWM Gravity Wave Drag', $
 	subtitle='Migrating Diurnal Tide'
	xyouts,60.,72.,'s=1'

	contour,rfsn,deglat,alt, $
	xrange=[-90,90],xticks=6,xminor=3,/xstyle, $
	yrange=[70,110],yticks=4,yminor=5,/ystyle, $
	xtitle='Latitude',ytitle='Altitude (km)', $
	title='GSWM Gravity Wave Drag', $
	subtitle='Non-Migrating Diurnal Tide'
	xyouts,60.,72.,'s!9=!51'

	stop
	end
