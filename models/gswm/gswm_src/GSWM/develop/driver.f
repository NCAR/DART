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
	character*8 fname
	character*13 testname
	character*13 backname
	character*12 outname
	character*2 string
	character*72 comment
	character*3 onoff(3)
	character*11 wavtype(3)
	character*9  months(12)
	character*12 windbc(2)
	character*17 heat_f(7)
	character*20 heat_nmf(3)
	character*8 eddy_m(3)
	character*5 o3file(5)
	character*3 latentin(2)
	character*3 scycle(3)
c       
	integer i,j,error,wback,ipw
c       
	integer mx,ny,iz,f107
	integer nss,isize
	integer igm, ihd, iyp, ibr, nowind, ihm, latgrad, polebc
	integer heatmodel, dblhme11, diurno3, peako3, skipo3, skipir
	integer heatnm,backf(12),mig,mois,nzonal,wback
	integer pwforce,pwheat,pwdelw
	integer call_late, o3conc, iboris, eddymodel,kzzinvar, dissinvar
	integer iondrag, vialdiss, gwstress, idecomp, idenpress
c       
	real period,flux
	real euvlat(36),reuv(36,36),ieuv(36,36),euvalt(36)
	real sigma
	real deuv(36,36,3)
c
c Integer common of flag that tells whether to write to background file

	common/write/wback
c
c Integer common of flags for wave characteristics
	
	common /mode/nzonal,period,freq,mois,nss,flux
c
c Integer common of flags for background atmosphere

	common /backatm/igm,ihd,iyp,ibr,nowind,ihm,latgrad,polebc
c
c Integer common of flags for heat forcing

	common /heatforce/heatmodel,dblhme11,diurno3,peako3,skipo3,
     +     skipir,call_late,o3conc,heatnm,pwforce,pwheat,pwdelw
c
c Integer common of flags for dissipation

	common /boris/iboris,eddymodel,kzzinvar,dissinvar,iondrag,
     +  vialdiss,gwstress
c
c Integer common of flags for postprocessing

	common /postproc/idecomp,idenpress,mig
c
c Common block of integers, arrays, floats, for EUV thermosphere heating

      common/euv/mx,ny,iz,euvlat,euvalt,reuv,ieuv,deuv,sigma,f107
c
	data backf/0,0,0,0,0,0,0,0,0,0,0,0/
	data onoff/"OFF","ON","ON*"/
	data wavtype/'diurnal','semidiurnal','terdiurnal'/
	data months/'January','February','March','April','May',
     +  'June','July','August','September','October','November',
     +  'December'/
	data windbc/'Lateral','Taylor'/
	data heat_f/'Hough11Simple','Heat92','Latent_heat_tides',
     +	    'Honglate','Heat92+Late','Heat92+Hong','ThermosphereEUV'/
	data heat_nmf/'Latent Heat','Infrared','Combined IR + Latent'/
	data eddy_m/"Weak","Eddy_gs","Eddy"/
	data o3file/"None","CIRA","CLAES","MLS","HALOE"/
	data latentin/"MON", "ANN"/
	data scycle /"MIN","MED","MAX"/
c
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
	print*,''

	do i= 1, isize
	   wback=1
	   error=0
c
c Read in wave

	   read(11,35)period,nzonal,mois,fname
 35	   format(4x,f3.2,12x,i2,6x,i2,9x,a8)
	   ipw=0
	   if(period.gt.1.) ipw=1
c
c Error checking

	   if(nzonal.lt.-6.or.nzonal.gt.6) error=1
	   if(period.lt.0.) error=1
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
c Check if migrating or nonmigrating

	   if(nzonal*period.eq.1.)then
	      mig=1
	      print*,i," Migrating Tide Calculation: ",testname
	   else
	      mig=0
	      print*,i,' Non-Migrating Tide Calculation: ',testname
	   endif

	   if (period.gt.1) then
	   mig=-999.
	   print*,i,' Planetary Wave Calculation: ',testname

c------------------------------
c
c Open the control file and check for errors

	   open(unit=14, file="GSWM.inp",status="old")
c
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
	   read(14,10,err=100)nowind
	   read(14,10,err=100)igm
	   read(14,10,err=100)ihd
	   read(14,10,err=100)iyp
	   read(14,10,err=100)ibr
	   read(14,10,err=100)ihm
	   read(14,10,err=100)latgrad
	   read(14,10,err=100)polebc
	   read(14,10,err=100)o3conc
c       
	   read(14,*)		!---------------
c
	   if(mig.eq.1)then
	      read(14,10,err=100)heatmodel
	      read(14,10,err=100)dblhme11
	      read(14,10,err=100)diurno3
	      read(14,10,err=100)peako3
	      read(14,10,err=100)skipo3
	      read(14,10,err=100)skipir
	      read(14,10,err=100)call_late
	      read(14,10,err=100)f107
	      read(14,*)
	      read(14,*)
	      heatnm=-999
	      pwforce=-999

	   elseif(mig.eq.0)then
	      read(14,*)
	      read(14,*)
	      read(14,*)
	      read(14,*)
	      read(14,*)
	      read(14,*)
	      read(14,*)
	      read(14,*)
	      read(14,10,err=100)heatnm
	      heatmodel=-999
	      pwforce=-999

	   elseif(mig.eq.-999)then
	      read(14,10,err=100)pwforce
	      read(14,10,err=100)pwheat
	      read(14,10,err=100)pwdelw
	      heatmodel=-999
	      heatnm=-999

	   endif
c       
	   read(14,*)		!---------------
	   read(14,10,err=100)iboris
	   read(14,10,err=100)eddymodel
c       read(14,10,err=100)kzzinvar
	   read(14,10,err=100)dissinvar
	   read(14,10,err=100)iondrag
	   read(14,10,err=100)vialdiss
	   read(14,10,err=100)gwstress
	   read(14,*)		!---------------
	   read(14,10,err=100)idecomp
	   read(14,10,err=100)idenpress
c
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
	   if(nowind.lt.0.or.nowind.gt.1)error=1
	   if(igm.lt.0.or.igm.gt.2) error=1
	   if(ihd.lt.0.or.ihd.gt.2) error=1
	   if(iyp.lt.0.or.iyp.gt.1) error=1
	   if(ibr.lt.0.or.ibr.gt.2) error=1
	   if(ihm.lt.0.or.ihm.gt.2) error=1
	   if(latgrad.lt.0.or.latgrad.gt.1)error=1
c GSWM98 no horz diffusion: Taylor or Lateral BC only
	   if(polebc.lt.0.or.polebc.gt.1) error=1
	   if(o3conc.lt.0.or.o3conc.gt.4) error=1
c       
	   if(heatmodel.lt.0.or.heatmodel.gt.6)then
	      if(heatnm.eq.-999.and.pwforce.eq.-999)error=1
	   endif
	   if(heatnm.lt.0.or.heatnm.gt.2) then
	      if(heatmodel.eq.-999.and.pwforce.eq.-999)error=1
	   endif
	   if(pwforce.lt.0.or.pwforce.gt.1) then
	      if(heatmodel.eq.-999.and.heatnm.eq.-999)error=1
	   endif
	   if(dblhme11.lt.0.or.dblhme11.gt.1) error=1
	   if(diurno3.lt.0.or.diurno3.gt.2) error=1
	   if(peako3.lt.0.or.peako3.gt.1) error=1
	   if(skipo3.lt.0.or.skipo3.gt.1) error=1
	   if(skipir.lt.0.or.skipir.gt.1) error=1
	   if(call_late.lt.0.or.call_late.gt.1)error=1
c       
	   if(iboris.lt.0.or.iboris.gt.1) error=1
	   if(eddymodel.lt.0.or.eddymodel.gt.2) error=1
	   if(dissinvar.lt.0.or.dissinvar.gt.1) error=1
	   if(iondrag.lt.0.or.iondrag.gt.1) error=1
	   if(vialdiss.lt.0.or.vialdiss.gt.1) error=1
	   if(gwstress.lt.0.or.gwstress.gt.1) error=1
c       
	   if(idecomp.lt.0.or.idecomp.gt.1) error=1
	   if(idenpress.lt.0.or.idenpress.gt.1)error=1
c
c Now specifically incompatible errors between flags------------------c
c
      if(heatmodel.ne.0.and.dblhme11.eq.1)then
         print*,'Warning:  Hough forcing cannot be doubled:  Hough'
         print*,'forcing not specified (HEATING).'
      endif
c
      if(iboris.eq.1)then
         if((iondrag+gwstress).gt.0)then
          print*,'Error:  Disable ionosphere or gravity wave drag'
          print*,'if using Khatattov!  Set to 0.'
          error=1
         elseif(dissinvar.eq.1)then
          print*,'Error:  Disable constant MU0, K0 dissipation if'
          print*,'using Khatattov! Set invariant dissipation to 0.'
          error=1
         elseif(eddymodel.ne.1)then
          print*,'Error: Disable viscous eddy dissipation if using'
          print*,'Khattatov!  Set eddymodel to 1.'
          error=1
         endif
      endif
c
      if(nowind.eq.1.and.igm+ihd+iyp+ibr.gt.0)then
         print*,'Warning: NOWIND =1 and you have requested background'
         print*,'winds.  All background winds will be 0!'
      endif
c
      if(igm.ne.0.and.ihd.eq.2)then
         print*,'Warning: HRDI winds will be used for entire background'
      endif
c
      if(igm.eq.1.and.ihd+iyp+ibr.gt.0)then
         print*,'Warning: Groves winds chosen for background. Will be'
         print*,'overwritten by other winds.'
      endif
c
      if(ihd.eq.2.and.igm+iyp+ibr.gt.0)then
         print*,'Warning: HRDI winds chosen for background. Will be'
         print*,'overwritten by other winds.'
      endif
c
      if(ihm.eq.1.and.igm+ihd+iyp+ibr.gt.0)then
         print*,'Warning: Incompatible wind background. You have'
         print*,'requested HWM93.  All other switches should be 0!'
      endif
c
      if(ihm.eq.2.and.ihd.ne.2)then
         print*,'Warning: Incompatible wind background. You have'
         print*,'requested HWM93+HRDI/Groves.  BOTH switches=2!'
      endif
c
      if(o3conc.gt.4.or.o3conc.lt.0)then
         print*,'Error.  o3conc out of range'
         error=1
      endif
c
      if(o3conc.ne.0.and.heatmodel.ne.1)then
       print*,'Warning:choice of heating model does not use [O3] file!'
      endif
c
      if(o3conc.eq.0.)then
         if(heatmodel.eq.1)then
           print*,'Warning:HEATING model requires an [O3] background'
           print*,'Default HALOE used.'
           o3conc=4
         endif
         if(heatmodel.eq.4)then
           print*,'Warning:HEATING model requires an [O3] background'
           print*,'Default HALOE used.'
           o3conc=4
         endif
         if(heatmodel.eq.5)then
           print*,'Warning:HEATING model requires an [O3] background'
           print*,'Default HALOE used.'
           o3conc=4
         endif
      endif
c
	if(heatmodel.eq.1.or.heatmodel.eq.4.or.heatmodel.eq.5)then
	   if(mod((mois+3),3).ne.1)error=1
	endif
c
      if((skipo3+skipir).gt.0)then
         if(heatmodel.ne.1)then
            print*,'Warning: Ozone or water vapor-forced heating'
            print*,'cannot be suppressed in this heat model!'
         endif
      endif
c
c set the flux

	if(heatmodel.eq.6)then
	   if(f107.eq.1)then
	      flux=70.
	   elseif(f107.eq.2)then
	      flux=130.
	   elseif(f107.eq.3)then
	      flux=200.
	   endif
	else
	   flux=120.
	endif
c
c      if(gwstress.eq.1.and.period.ne.1.)then
c         print*,'Error: Apply GW Stress only for diurnal waves.'
c         error=1
c      endif
c
c      if(vialdiss.eq.1.and.eddymodel.ne.2)then
c         print*,'Warning: Vial eddy dissipation not available.'
c      endif  
c---------------------------------------------------------------------c
c      
c Diagnostic output for debugging

      print*,'============================='
      if(ipw.eq.0.)then
      print*,'Tide Type                       ',wavtype(int(1/period))
      print*,'Tide Period.................',period
      elseif(ipw.eq.1)then
      print*,'Planetary Wave                  '
      print*,'Wave Period.................',period
      endif
      print*,'Wavenumber                   ',nzonal
      print*,'Month...........................',months(mois)
      print*,'-------------------------------------'
      print*,'Background Winds = 0            ',onoff(nowind+1)
      print*,'Groves Winds....................',onoff(igm+1)
      print*,'HRDI Winds                      ',onoff(ihd+1)
      print*,'Portnyagin Winds................',onoff(iyp+1)
      print*,'Randel Winds                    ',onoff(ibr+1)
      print*,'HWM93 Winds                     ',onoff(ihm+1)
      print*,'Latitude Gradients = 0..........',onoff(latgrad+1)
      print*,'Wind BC at Poles:               ',windbc(polebc+1)
      print*,'Background Ozone Concentration..',o3file(o3conc+1)
      print*,'-------------------------------------'
	if(mig.eq.1)then
      print*,'Heat Model                       ',heat_f(heatmodel+1)
      print*,'Double Hough Heating.............',onoff(dblhme11+1)
      print*,'Diurnal O3 Variation             ',onoff(diurno3+1)
      print*,'Second O3 Peak at 92km...........',onoff(peako3+1)
      print*,'Skip O3 Heat Forcing             ',onoff(skipo3+1)
      print*,'Skip H2O Heat Forcing (IR).......',onoff(skipir+1)
      print*,'Call_late                        ',latentin(call_late+1)
      print*,'Solar Cycle......................',scycle(f107)
	elseif(mig.eq.0)then
      print*,'Non-Migrating Heat model         ',heat_nmf(heatnm+1)
	elseif(mig.eq.-999)then
      print*,'Planetary Wave Heat Source         ',onoff(pwheat+1)
      print*,'Planetary Wave Vertical Velocity   ',onoff(pwdelw+1)
	endif
      print*,'-------------------------------------'
      print*,'Khattatov Dissipation            ',onoff(iboris+1)
      print*,'Eddy model.......................',eddy_m(eddymodel+1)
c      print*,'Kzz invar ',onoff(kzzinvar+1)
      print*,'Latitude-constant Dissipation    ',onoff(dissinvar+1)
      print*,'Ion Drag.........................',onoff(iondrag+1)
      print*,'Vial Dissipation (Eddy only)     ',onoff(vialdiss+1)
      print*,'Gravity Wave Stress (24hr only)..',onoff(gwstress+1)
      print*,'-------------------------------------'
      print*,'Hough decomposition              ',onoff(idecomp+1)
      print*,'Den, Press perturbation..........',onoff(idenpress+1)
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
	      print*,''
	   
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
	   print*,''
	   
 1000	   continue
	enddo			!all the tides requested

	close(11)		!input file
c

	stop
	end

