      SUBROUTINE getlate(nzonal,period,mois,ij)

c----- This subroutine reads in real and imaginary components of 
C	tidal forcing for use in latent_heat_tides
C					M. Hagan (5/2/96)

	character*15 lchoice(2)
	character*15 lfile

      real HJRm(91),HJIm(91)
      real HJR(91,12),HJI(91,12)

      common/lateforce/HJRm,HJIm

	data lchoice/'7yr-12mos.diurn','7yr-12mos.semid'/

      DO 2 k=1,12
      DO 2 I=1,IJ
      HJr(i,k)=0.0
      HJi(i,k)=0.0
2     CONTINUE

100	format(////////////)
c101	format(/)
101	format(///////////)
111	format(8(7(1x,e10.3),/),3(1x,e10.3),/)

	lfile=lchoice(nzonal)
	OPEN(UNIT=103,file=lfile,form='formatted',STATUS='old')
	
      DO 1 nm=1,mois
	if(nm.eq.1) read(103,100)
	if(nm.gt.1) read(103,101)
	read(103,111)(hjr(i,nm),i=1,ij)
	read(103,111)(hji(i,nm),i=1,ij)
1     CONTINUE

	CLOSE(UNIT=103)

      DO 3 i=1,ij
      HJrm(i)=hjr(i,mois)
      HJim(i)=hji(i,mois)
3     CONTINUE
        print*,'In subroutine getlate. Migrating only!!'
 	print *, mois,nzonal,period
 	print *, (hjrm(i), i=1,ij)
 	print *, (hjim(i), i=1,ij)

      RETURN
      END

