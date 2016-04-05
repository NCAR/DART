	subroutine eddy_gs(z,clt,eddyv,eddyk,deddyv,deddyk)

C  This subroutine interpolates the Garcia/Solomon eddy tables set
C  up by subroutine seteddy_gs  .        

        common/interpe/zped(13,30,3),xe(13),ye(30),zedd(13,30),sigmae

	edd=10.
        fact=1.

	if(z.le.116)then
	   ze=z

	   call surfd(clt,ze,edd,dedt,dedz,d2edt,d2edzdt,d2edz,13,30,
     1       xe,ye,zedd,13,zped,sigmae)

	   if(z.gt.70..and.edd.lt.50.) edd=50.

	endif

        eddyv=fact*edd
        deddyv=(1./eddyv)*fact*dedz*.001

	eddyk=1.00*eddyv
	deddyk=1.00*deddyv

	return
	end
