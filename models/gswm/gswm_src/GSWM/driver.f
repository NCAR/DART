	program tide driver
c
c This program opens four files: unit 11 (read) describes the waves to be run.
c  Unit 7 (write)
c contains the program output.  Unit 10 (wrrite) contains the background
c  atmosphere for
c each month.  Unit 14 (read) contains control flags for the run.
c
c Modified for HWM93 option; ihm switch		M. Hagan (10/19/98)
c
c
c Modified for SR h_diff option; ihdiff switch		M. Hagan (7/15/99)
c
c Added idca switch and Redefined mig switch. Modified logic to 
c allow W1 diurnal and W2 semi 7yr monthly DCA forcing  (M. Hagan 7/27/99)
c
C
C Stripped for standard inputs and migrating calculations only (M. Hagan 3/04)
C
	character*8 fname
	character*13 testname
	character*13 backname
	character*12 outname
	character*2 string
	character*72 comment
	character*3 onoff(3)
	character*11 wavtype(3)
	character*9 months(12)
	character*7 windbc(2)
c       
	integer i,j,error,wback
c       
	integer mx,ny,iz
	integer nss,isize
	integer latgrad, polebc
	integer backf(12),mois,nzonal
c       
	real period,flux
	real sigma
c
c Integer common of flag that tells whether to write to background file

	common/ritecom/wback
c
c Integer common of flags for wave characteristics
	
      COMMON /MODEi/NZONAL,MOIS,NSS
      COMMON /MODEr/PERIOD,FREQ,FLUX
c
c Integer common of flags for background atmosphere

	common /backatm/latgrad,polebc

	data backf/0,0,0,0,0,0,0,0,0,0,0,0/
	data onoff/"OFF","ON ","ON "/
	data wavtype/'diurnal    ','semidiurnal','terdiurnal '/
	data months/'January  ','February ','March    ',
     +  'April    ','May      ',
     +  'June     ','July     ','August   ','September',
     +  'October  ','November ','December '/
	data windbc/"Lateral","Taylor "/

	flux=120
	nss=3
	j=1
	error=0
c
C       The input file with the names, months, of waves

	  open(11, file= "waves.inp", status= "old")
c
c Read the number of waves to solve

	  read(11,*)isize	!number of files to run
	  read(11,*)		!blank line
	  read(11,*)		!blank line
	  read(11,*)		!blank line


	print *,isize," Tides will be evaluated."
	print*,' '

	do i= 1, isize
	   wback=1
	   error=0
c
c Read in wave

	   read(11,35)period,nzonal,mois,fname
 35	   format(4x,f3.2,12x,i2,6x,i2,9x,a8)
c
c Error checking

	   if(nzonal.lt.1.or.nzonal.gt.2) error=1
	   if(period.lt.0..or.period.gt.1.) error=1
	   if(mois.lt.0.or.mois.gt.12)error=1
c
c Construct the filenames

	   testname=fname
	   backname=fname
	   outname=fname
	   outname(9:12)=".out"
	   testname(9:13)=".test"
	   backname(9:13)=".back"
c
c
c Open the control file and check for errors

	   open(unit=14, file="GSWM.inp",status="old")
c
c Skip the comments at the top of the input file:  marked with '&'

	   string(1:1)='&'
	   do while(string(1:1).eq.'&')
	      read (14, *) string
	   end do
c Read the input file.  Error message is line 100.

c       
c	   read(14,*)		!---------------
	   read(14,40) comment
	   print*,comment
c
	   read(14,*)		!---------------
	   read(14,10,err=100)latgrad
	   read(14,10,err=100)polebc

c Check if migrating 
	   if(nzonal*period.eq.1.)then
	      print*,i," Migrating Tide Calculation: ",testname
	   endif

c------------------------------
 10	   format(32x,i1)
 20	   format(32x,f3.1)
 30	   format(32x,i2)
 40	   format(A72)
c

c Close the input file
	   
	   close(14)
c
c Error checking

c First, the proper range of values

c       
	   if(latgrad.lt.0.or.latgrad.gt.1)error=1
	   if(polebc.lt.1.or.polebc.gt.2) error=1
c       
c---------------------------------------------------------------------c
c Diagnostic output for debugging

      print*,'============================='
      print*,'Tide Type                       ',wavtype(int(1/period))
      print*,'Tide Period.....................',period
      print*,'Wavenumber......................',nzonal
      print*,'Month...........................',months(mois)
      print*,'-------------------------------------'
      print*,'Latitude Gradients..............',onoff(latgrad+1)
      print*,'Wind BC at Poles................',windbc(polebc)
      print*,'-------------------------------------'
c
c------------------------------
c
c If error

	   if(error.eq.1) then
	      go to 100
	   else
c
c Track number of times one month is requested (for writing background file)

	      backf(mois)=backf(mois)+1

C NOTE: When running multiple wavenumbers only print background out if nzonal=1
CCCC      Only print background once per month, too.
C	  see logical if in SR SETATMOS

c
c Check if .back is already open for the month

	      if(backf(mois).gt.1)then 
		 wback=0
		 print*,'Background file already written for ',months(mois)
	      endif
c       
c Open background file
	      
	      if(wback.eq.1)then
		 open(10, file=backname)
		 print*,'Background file opened for ',months(mois)
	      endif
c       
c Open output file if everything looks OK

	      open (7, file=testname)
	   
c       print*,'Running ',period,',',nzonal,' for ',months(mois)
c
CTEST write comment to output
	      write(7,200) comment
	      write(7,210) testname
	      write(7,*)
 200	      format(A72)
 210	      format(A13)
CTEST write end
c
c Run the progam

	      call main_2
c       
c       print*,'Finished ',period,',',nzonal,' for ',months(mois)
	      print*,' '
	   
	      close(7)		!output
	   
	      if(wback.eq.1)then
		 close(10)	!background
	      endif
	   endif		!if no errors in input
c       
	   go to 1000		!if everything's OK don't print err msg
c       
 100	   write(*,110) i
 110	   format(1x,i2,1x,'Error in format or content of input file.')
	   print*,fname,' not run.'
	   print*,' '
	   
 1000	   continue
	enddo			!all the tides requested

	close(11)		!input file
c

	stop
	end

