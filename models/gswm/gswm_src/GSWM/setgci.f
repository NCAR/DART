	subroutine setgci(nzonal,period,mois,ij)
C_________________________________________________________________________

C  This subroutine reads in complex arrays of tidal heating rates based
C  on Forbes and Zhang's analysis of Global Cloud Imagery Data.
C		latitudes: -90 to +90 in steps of 3 degree
C  for a specific month and zonal wavenumber
C_________________________________________________________________________

	character*18 ochoice(13,12)
	character*18 schoice(13,12)
	character*18 ofile

	real reelgci(91),imaggci(91)
	complex gciforce(91)

      COMMON/gci/gciforce

      data (ochoice(i,1),i=1,13)/'7yr-janW1.diurn_sm',
     +    '7yr-janW2.diurn_sm',
     +	  '7yr-janW3.diurn_sm','7yr-janW4.diurn_sm','7yr-janW5.diurn_sm',
     +	  '7yr-janW6.diurn_sm','7yr-janS0.diurn_sm','7yr-janE1.diurn_sm',
     +	  '7yr-janE2.diurn_sm','7yr-janE3.diurn_sm','7yr-janE4.diurn_sm',
     +	  '7yr-janE5.diurn_sm','7yr-janE6.diurn_sm'/

      data (ochoice(i,2),i=1,13)/'7yr-febW1.diurn_sm',
     +    '7yr-febW2.diurn_sm',
     +	  '7yr-febW3.diurn_sm','7yr-febW4.diurn_sm','7yr-febW5.diurn_sm',
     +	  '7yr-febW6.diurn_sm','7yr-febS0.diurn_sm','7yr-febE1.diurn_sm',
     +	  '7yr-febE2.diurn_sm','7yr-febE3.diurn_sm','7yr-febE4.diurn_sm',
     +	  '7yr-febE5.diurn_sm','7yr-febE6.diurn_sm'/

      data (ochoice(i,3),i=1,13)/'7yr-marW1.diurn_sm',
     +    '7yr-marW2.diurn_sm',
     +	  '7yr-marW3.diurn_sm','7yr-marW4.diurn_sm','7yr-marW5.diurn_sm',
     +	  '7yr-marW6.diurn_sm','7yr-marS0.diurn_sm','7yr-marE1.diurn_sm',
     +	  '7yr-marE2.diurn_sm','7yr-marE3.diurn_sm','7yr-marE4.diurn_sm',
     +	  '7yr-marE5.diurn_sm','7yr-marE6.diurn_sm'/

      data (ochoice(i,4),i=1,13)/'7yr-aprW1.diurn_sm',
     +    '7yr-aprW2.diurn_sm',
     +	  '7yr-aprW3.diurn_sm','7yr-aprW4.diurn_sm','7yr-aprW5.diurn_sm',
     +	  '7yr-aprW6.diurn_sm','7yr-aprS0.diurn_sm','7yr-aprE1.diurn_sm',
     +	  '7yr-aprE2.diurn_sm','7yr-aprE3.diurn_sm','7yr-aprE4.diurn_sm',
     +	  '7yr-aprE5.diurn_sm','7yr-aprE6.diurn_sm'/

      data (ochoice(i,5),i=1,13)/'7yr-mayW1.diurn_sm',
     +    '7yr-mayW2.diurn_sm',
     +	  '7yr-mayW3.diurn_sm','7yr-mayW4.diurn_sm','7yr-mayW5.diurn_sm',
     +	  '7yr-mayW6.diurn_sm','7yr-mayS0.diurn_sm','7yr-mayE1.diurn_sm',
     +	  '7yr-mayE2.diurn_sm','7yr-mayE3.diurn_sm','7yr-mayE4.diurn_sm',
     +	  '7yr-mayE5.diurn_sm','7yr-mayE6.diurn_sm'/

      data (ochoice(i,6),i=1,13)/'7yr-junW1.diurn_sm',
     +    '7yr-junW2.diurn_sm',
     +	  '7yr-junW3.diurn_sm','7yr-junW4.diurn_sm','7yr-junW5.diurn_sm',
     +	  '7yr-junW6.diurn_sm','7yr-junS0.diurn_sm','7yr-junE1.diurn_sm',
     +	  '7yr-junE2.diurn_sm','7yr-junE3.diurn_sm','7yr-junE4.diurn_sm',
     +	  '7yr-junE5.diurn_sm','7yr-junE6.diurn_sm'/

      data (ochoice(i,7),i=1,13)/'7yr-julW1.diurn_sm',
     +    '7yr-julW2.diurn_sm',
     +	  '7yr-julW3.diurn_sm','7yr-julW4.diurn_sm','7yr-julW5.diurn_sm',
     +	  '7yr-julW6.diurn_sm','7yr-julS0.diurn_sm','7yr-julE1.diurn_sm',
     +	  '7yr-julE2.diurn_sm','7yr-julE3.diurn_sm','7yr-julE4.diurn_sm',
     +	  '7yr-julE5.diurn_sm','7yr-julE6.diurn_sm'/

      data (ochoice(i,8),i=1,13)/'7yr-augW1.diurn_sm',
     +    '7yr-augW2.diurn_sm',
     +	  '7yr-augW3.diurn_sm','7yr-augW4.diurn_sm','7yr-augW5.diurn_sm',
     +	  '7yr-augW6.diurn_sm','7yr-augS0.diurn_sm','7yr-augE1.diurn_sm',
     +	  '7yr-augE2.diurn_sm','7yr-augE3.diurn_sm','7yr-augE4.diurn_sm',
     +	  '7yr-augE5.diurn_sm','7yr-augE6.diurn_sm'/

      data (ochoice(i,9),i=1,13)/'7yr-sepW1.diurn_sm',
     +    '7yr-sepW2.diurn_sm',
     +	  '7yr-sepW3.diurn_sm','7yr-sepW4.diurn_sm','7yr-sepW5.diurn_sm',
     +	  '7yr-sepW6.diurn_sm','7yr-sepS0.diurn_sm','7yr-sepE1.diurn_sm',
     +	  '7yr-sepE2.diurn_sm','7yr-sepE3.diurn_sm','7yr-sepE4.diurn_sm',
     +	  '7yr-sepE5.diurn_sm','7yr-sepE6.diurn_sm'/

      data (ochoice(i,10),i=1,13)/'7yr-octW1.diurn_sm',
     +    '7yr-octW2.diurn_sm',
     +	  '7yr-octW3.diurn_sm','7yr-octW4.diurn_sm','7yr-octW5.diurn_sm',
     +	  '7yr-octW6.diurn_sm','7yr-octS0.diurn_sm','7yr-octE1.diurn_sm',
     +	  '7yr-octE2.diurn_sm','7yr-octE3.diurn_sm','7yr-octE4.diurn_sm',
     +	  '7yr-octE5.diurn_sm','7yr-octE6.diurn_sm'/

      data (ochoice(i,11),i=1,13)/'7yr-novW1.diurn_sm',
     +    '7yr-novW2.diurn_sm',
     +	  '7yr-novW3.diurn_sm','7yr-novW4.diurn_sm','7yr-novW5.diurn_sm',
     +	  '7yr-novW6.diurn_sm','7yr-novS0.diurn_sm','7yr-novE1.diurn_sm',
     +	  '7yr-novE2.diurn_sm','7yr-novE3.diurn_sm','7yr-novE4.diurn_sm',
     +	  '7yr-novE5.diurn_sm','7yr-novE6.diurn_sm'/

      data (ochoice(i,12),i=1,13)/'7yr-decW1.diurn_sm',
     +    '7yr-decW2.diurn_sm',
     +	  '7yr-decW3.diurn_sm','7yr-decW4.diurn_sm','7yr-decW5.diurn_sm',
     +	  '7yr-decW6.diurn_sm','7yr-decS0.diurn_sm','7yr-decE1.diurn_sm',
     +	  '7yr-decE2.diurn_sm','7yr-decE3.diurn_sm','7yr-decE4.diurn_sm',
     +	  '7yr-decE5.diurn_sm','7yr-decE6.diurn_sm'/

      data (schoice(i,1),i=1,13)/'7yr-janW1.semid_sm',
     +    '7yr-janW2.semid_sm',
     +	  '7yr-janW3.semid_sm','7yr-janW4.semid_sm','7yr-janW5.semid_sm',
     +	  '7yr-janW6.semid_sm','7yr-janS0.semid_sm','7yr-janE1.semid_sm',
     +	  '7yr-janE2.semid_sm','7yr-janE3.semid_sm','7yr-janE4.semid_sm',
     +	  '7yr-janE5.semid_sm','7yr-janE6.semid_sm'/

      data (schoice(i,2),i=1,13)/'7yr-febW1.semid_sm',
     +    '7yr-febW2.semid_sm',
     +	  '7yr-febW3.semid_sm','7yr-febW4.semid_sm','7yr-febW5.semid_sm',
     +	  '7yr-febW6.semid_sm','7yr-febS0.semid_sm','7yr-febE1.semid_sm',
     +	  '7yr-febE2.semid_sm','7yr-febE3.semid_sm','7yr-febE4.semid_sm',
     +	  '7yr-febE5.semid_sm','7yr-febE6.semid_sm'/

      data (schoice(i,3),i=1,13)/'7yr-marW1.semid_sm',
     +    '7yr-marW2.semid_sm',
     +	  '7yr-marW3.semid_sm','7yr-marW4.semid_sm','7yr-marW5.semid_sm',
     +	  '7yr-marW6.semid_sm','7yr-marS0.semid_sm','7yr-marE1.semid_sm',
     +	  '7yr-marE2.semid_sm','7yr-marE3.semid_sm','7yr-marE4.semid_sm',
     +	  '7yr-marE5.semid_sm','7yr-marE6.semid_sm'/

      data (schoice(i,4),i=1,13)/'7yr-aprW1.semid_sm',
     +    '7yr-aprW2.semid_sm',
     +	  '7yr-aprW3.semid_sm','7yr-aprW4.semid_sm','7yr-aprW5.semid_sm',
     +	  '7yr-aprW6.semid_sm','7yr-aprS0.semid_sm','7yr-aprE1.semid_sm',
     +	  '7yr-aprE2.semid_sm','7yr-aprE3.semid_sm','7yr-aprE4.semid_sm',
     +	  '7yr-aprE5.semid_sm','7yr-aprE6.semid_sm'/

      data (schoice(i,5),i=1,13)/'7yr-mayW1.semid_sm',
     +    '7yr-mayW2.semid_sm',
     +	  '7yr-mayW3.semid_sm','7yr-mayW4.semid_sm','7yr-mayW5.semid_sm',
     +	  '7yr-mayW6.semid_sm','7yr-mayS0.semid_sm','7yr-mayE1.semid_sm',
     +	  '7yr-mayE2.semid_sm','7yr-mayE3.semid_sm','7yr-mayE4.semid_sm',
     +	  '7yr-mayE5.semid_sm','7yr-mayE6.semid_sm'/

      data (schoice(i,6),i=1,13)/'7yr-junW1.semid_sm',
     +    '7yr-junW2.semid_sm',
     +	  '7yr-junW3.semid_sm','7yr-junW4.semid_sm','7yr-junW5.semid_sm',
     +	  '7yr-junW6.semid_sm','7yr-junS0.semid_sm','7yr-junE1.semid_sm',
     +	  '7yr-junE2.semid_sm','7yr-junE3.semid_sm','7yr-junE4.semid_sm',
     +	  '7yr-junE5.semid_sm','7yr-junE6.semid_sm'/

      data (schoice(i,7),i=1,13)/'7yr-julW1.semid_sm',
     +    '7yr-julW2.semid_sm',
     +	  '7yr-julW3.semid_sm','7yr-julW4.semid_sm','7yr-julW5.semid_sm',
     +	  '7yr-julW6.semid_sm','7yr-julS0.semid_sm','7yr-julE1.semid_sm',
     +	  '7yr-julE2.semid_sm','7yr-julE3.semid_sm','7yr-julE4.semid_sm',
     +	  '7yr-julE5.semid_sm','7yr-julE6.semid_sm'/

      data (schoice(i,8),i=1,13)/'7yr-augW1.semid_sm',
     +    '7yr-augW2.semid_sm',
     +	  '7yr-augW3.semid_sm','7yr-augW4.semid_sm','7yr-augW5.semid_sm',
     +	  '7yr-augW6.semid_sm','7yr-augS0.semid_sm','7yr-augE1.semid_sm',
     +	  '7yr-augE2.semid_sm','7yr-augE3.semid_sm','7yr-augE4.semid_sm',
     +	  '7yr-augE5.semid_sm','7yr-augE6.semid_sm'/

      data (schoice(i,9),i=1,13)/'7yr-sepW1.semid_sm',
     +    '7yr-sepW2.semid_sm',
     +	  '7yr-sepW3.semid_sm','7yr-sepW4.semid_sm','7yr-sepW5.semid_sm',
     +	  '7yr-sepW6.semid_sm','7yr-sepS0.semid_sm','7yr-sepE1.semid_sm',
     +	  '7yr-sepE2.semid_sm','7yr-sepE3.semid_sm','7yr-sepE4.semid_sm',
     +	  '7yr-sepE5.semid_sm','7yr-sepE6.semid_sm'/

      data (schoice(i,10),i=1,13)/'7yr-octW1.semid_sm',
     +    '7yr-octW2.semid_sm',
     +	  '7yr-octW3.semid_sm','7yr-octW4.semid_sm','7yr-octW5.semid_sm',
     +	  '7yr-octW6.semid_sm','7yr-octS0.semid_sm','7yr-octE1.semid_sm',
     +	  '7yr-octE2.semid_sm','7yr-octE3.semid_sm','7yr-octE4.semid_sm',
     +	  '7yr-octE5.semid_sm','7yr-octE6.semid_sm'/

      data (schoice(i,11),i=1,13)/'7yr-novW1.semid_sm',
     +    '7yr-novW2.semid_sm',
     +	  '7yr-novW3.semid_sm','7yr-novW4.semid_sm','7yr-novW5.semid_sm',
     +	  '7yr-novW6.semid_sm','7yr-novS0.semid_sm','7yr-novE1.semid_sm',
     +	  '7yr-novE2.semid_sm','7yr-novE3.semid_sm','7yr-novE4.semid_sm',
     +	  '7yr-novE5.semid_sm','7yr-novE6.semid_sm'/

      data (schoice(i,12),i=1,13)/'7yr-decW1.semid_sm',
     +    '7yr-decW2.semid_sm',
     +	  '7yr-decW3.semid_sm','7yr-decW4.semid_sm','7yr-decW5.semid_sm',
     +	  '7yr-decW6.semid_sm','7yr-decS0.semid_sm','7yr-decE1.semid_sm',
     +	  '7yr-decE2.semid_sm','7yr-decE3.semid_sm','7yr-decE4.semid_sm',
     +	  '7yr-decE5.semid_sm','7yr-decE6.semid_sm'/

       if(period.eq.1.) then
       if (nzonal.gt.0) ofile=ochoice(nzonal,mois)
       if (nzonal.le.0) then
       iloc=7+(abs(nzonal))
       ofile=ochoice(iloc,mois)
       endif
       endif

       if(period.eq..5) then
       if (nzonal.gt.0) ofile=schoice(nzonal,mois)
       if (nzonal.le.0) then
       iloc=7+(abs(nzonal))
       ofile=schoice(iloc,mois)
       endif
       endif

          write(6,*) 'Reading latent heating rates from ', ofile

	pi=acos(-1.)
	pio2=pi/2.
        dtor=pi/180.
        rtod=180./pi

C  Read in complex arrays:
	OPEN(UNIT=33,file=ofile,status='old')

201	format(/,8(7(1x,e10.3),/),(3(1x,e10.3)))
203	format(//////////)

	read(33,203)
	read(33,201) (reelgci(i),i=1,ij)
	read(33,201) (imaggci(i),i=1,ij)

	CLOSE(UNIT=33)

	do i=1,ij
	gciforce(i)=cmplx(reelgci(i),imaggci(i))
	enddo

      return
      end

 
