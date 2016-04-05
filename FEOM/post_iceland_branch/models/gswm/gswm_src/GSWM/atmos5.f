	subroutine atmos5(z,ij,g,gam,gam1)

C changes to construct "p" arrays from the north to south pole
C	consistent with changes to LATVECTS	--M. Hagan (12/11/92)

	real koo

C interp1 arrays are defined fron the north to the south poles:
      common/interp1/zpt(37,101,3),zpu(37,101,3),zpr(37,101,3),
     +	zpxh(37,101,3),dO2(37,101,3)
C interp2 arrays are defined fron the north to the south poles:
      common/interp2/x(37),y(101),zt(37,101),zu(37,101),zr(37,101),
     +  zxh(37,101),zph(37,101),o2(37,101),sigma

C latvect arrays are defined fron the north to the south poles:
      common/latvect/dtheta,clatt(91),xlat(91),snclat(91),csclat(91),
     +tnclat(91),ctnclat(91),sin2i(91),gmlat(91),dip(91)

      common/atm1/p1(91),p2(91),p3(91),p4(91),p5(91),p6(91),
     +p7(91),p8(91),p9(91),p10(91),p11(91),p12(91),p13(91),p14(91)
      common/atm2/dp1(91),dp3(91),dp8(91),dp9(91),d2p8(91),p15(91)

	data rs,koo/8.31432,0.015/
	data ro,go,poo/6378.165,9.80665,1.01325e+05/

100	format(1x,i3,6e11.3)

C________________________ Temperature Coefficients _________________________

	do 1 i=1,ij
	colat=clatt(i)
C       write(15,100) i, colat,z

	call surfd(colat,z,t,dtdt,dtdz,d2tdt,d2tdzdt,d2tdz,37,101,
     1       x,y,zt,37,zpt,sigma)
C	write(15,100) i,colat,t,dtdz,dtdt,d2tdz

	p3(i)=t
	p4(i)=.001*dtdz/t
	dp3(i)=dtdt/t
	p11(i)=1.0e-06*d2tdz/t
1       continue

C___________________________ Wind Coefficients ____________________________

	do 2 i=1,ij
	colat=clatt(i)
c       write(15,100) i, colat,z

	call surfd(colat,z,u,dudt,dudz,d2udt,d2udzdt,d2udz,37,101,
     1       x,y,zu,37,zpu,sigma)
c	write(15,100) i,colat,u,dudz,dudt,d2udt,d2udz 

	p1(i)=u
	p2(i)=.001*dudz
	dp1(i)=dudt
2	continue

C______________________________ Pressure Height  ____________________________

        do 3 i=1,ij
	colat=clatt(i)
C       write(15,100) i, colat,z

	call surfd(colat,z,xh,dxdt,dxdz,d2xdt,d2xdzdt,d2xdz,37,101,
     1       x,y,zxh,37,zpxh,sigma)
C	write(15,100) i,colat,xh,dxdz,dxdt,d2xdt,d2xdz

	p6(i)=xh
3	continue

C___________________________ Density Coefficients ____________________________

        do 4 i=1,ij
	colat=clatt(i)
C       write(15,100) i, colat,z

	call surfd(colat,z,rholog,drhodt,drhodz,d2dtdt,d2dzdt,
     1       d2dzdz,37,101,x,y,zr,37,zpr,sigma)
C	write(15,100) i,colat,rholog,drhodt,drhodz,d2dtdt,d2dzdz

	p8(i)=exp(rholog)
	p9(i)=.001*drhodz
	dp8(i)=drhodt
	d2p8(i)=(d2dtdt+drhodt*drhodt)
        p10(i)=(d2dzdz+drhodz*drhodz)*1.0e-06
	dp9(i)=(d2dzdt+drhodz*drhodt)*.001
4	continue

C_____________________________ Other Coefficients ____________________________

        g=go*(ro/(ro+z))**2.
	gam=1.4+.135*(1.+tanh((z-300.)/100.))
	gam1=gam-1.

	do 5 i=1,ij

	arg=-p6(i)
	p5(i)=poo*exp(arg)
	p5(i)=p5(i)/(g*p8(i))
	p14(i)=g*p5(i)/p3(i)
	p7(i)=rs*1000./p14(i)
	p13(i)=koo*(p3(i)**.66)/p7(i)
c
        p12(i)=.266*p13(i)/p14(i)
	p15(i)=-1./p5(i)-p4(i)-p9(i)
5	continue

	return
	end
	
