      subroutine h_diff(z,xlat,hdiff)
       real hdiff

      COMMON /MODEi/NZONAL,MOIS,NSS
      COMMON /MODEr/PERIOD,FREQ,FLUX

c  First, the basic hdiff distribution which gives no effect
c  on the s=1 component.  It has relatively small values at
c  low latitudes below 120 km.

	base=1.0e+05 + 9.0e+05*sin(abs(xlat))

	delta=0.0
	if(z.ge.100.and.z.le.150) then
	delta=1.0e+06-base
	end if

	hdiff=base + 0.5*delta*(1.0+sin((z-125.)*3.14159/50.))

c  Then, if other than s=0 or s=1 (any eastward propagating waves
c  and all westward propagating waves with s>1, add more hdiff,
c  increasing with higher wavenumbers.
C  Modified so that eastward waves correspond to negative zonal wave numbers
C					M. Hagan & J. Forbes (1/4/00)

      if(nzonal.lt.0.or.nzonal.gt.1) then
      hdiff = hdiff + abs(float(nzonal))*1.0e+06
C     hdiff = 5.0e+06
      end if

      return 
      end

