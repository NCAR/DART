	subroutine h2onm(nzonal,mois,z,hj,ij)
c  This subroutine obtains H2O Heating from tables read in
C	in SR seth2ononmig.f; created from heat92.f (M. Hagan-3/9/95)

      	complex hj(91)
	
        common/hth2o/zh2o(41),xlh2o(37),dh2o(37,41),zdh2o(37,41,3),
     +	sigma2,dah2o(37,41),zadh2o(37,41,3),sigma2a

      	common/latvect/dtheta,clatt(91),xlat(91),snclat(91),csclat(91),
     +	tnclat(91),ctnclat(91),sin2i(91),gmlat(91),dip(91)

501     format(1x,f5.1,4(2x,e11.3,e11.3),/,3x,4(2x,e11.3,e11.3))
      
        do 9 i=1,ij
    	hj(i)=(0.0,0.0)
    9   continue

	if (z.gt.20) go to 22

C*******************************************************************
C	H2O Heating from tables
C*******************************************************************

      Do 21 i=1,ij
	colat=clatt(i)

	call surfd(colat,z,h2o,dh2odt,dh2odz,d2h2odt,d2h2odzdt,
     1       d2h2odz,37,41,xlh2o,zh2o,dh2o,37,zdh2o,sigma2)

	call surfd(colat,z,h2oi,dh2oidt,dh2oidz,d2h2oidt,d2h2oidzdt,
     1       d2h2oidz,37,41,xlh2o,zh2o,dah2o,37,zadh2o,sigma2a)

	hj(i)=cmplx(h2o,h2oi)

   21   continue

c	write(10,501) z,(hj(i),i=1,ij)

   22   continue

	return
	end
