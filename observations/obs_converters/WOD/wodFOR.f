c    This code is not protected by the DART copyright agreement.
c    DART $Id$

      PROGRAM wodFOR
      
c    This program prints out to the screen data from WOD native format
c    ASCII file to the screen. This main program (wodFOR)
c    calls the subroutine WODread (WODread200X if the data file are in WOD05
c    or WOD01 format or the subroutine WODread1998 if the data are in WOD98
c    format). These subroutines do the actual reading of the ASCII format,
c    and load the data into arrays which are passed back to the main program.
c    The main program then work with these data arrays to print out the
c    data on the screen.  
c    
c    It is intended that the subroutine WODread provides an example of how
c    to extract the data and variables from the ASCII format,
c    whereas the main part of the wodFOR program provides an example of how
c    these data can be made accessible/workable as a series of arrays.
c    
c    Comments and suggestions for improving this program would be appreciated.
c    Updates to the World Ocean Data 2005 data and to this program will be posted
c    in the NODC/WOD web site at http://www.nodc.noaa.gov
c
c***********************************************************
c    
c    Missing values used in this dataset = bmiss = -999.99
c    
c***********************************************************
c     
c   Parameters (constants):
c     
c     maxlevel  - maximum number of depth levels, also maximum
c                   number of all types of variables
c     maxcalc   - maximum number of measured and calculated
c                   depth dependent variables
c     kdim      - number of standard depth levels
c     bmiss     - binary missing value marker
c     maxtcode  - maximum number of different taxa variable codes
c     maxtax    - maximum number of taxa sets
c     maxpinf   - maximum number of distinct measured variable
c                 information codes
c
c******************************************************************
      
      parameter (maxlevel=30000, maxcalc=100)
      parameter (kdim=40, bmiss=-999.99)
      parameter (maxtcode=25, maxtax=2000)
      parameter (maxpinf=25)
      
c******************************************************************
c
c   Character Arrays:
c
c     cc        - NODC country code
c     chars     - WOD character data: 1. originators cruise code,
c                                     2. originators station code
c     filename  - file name
c
c*****************************************************************
      
      character*2  cc
      character*15 chars(2)
      character*80 filename
 
c******************************************************************
c
c   Arrays:
c
c     isig()    - number of significant figures in (1) latitude, (2) longitude
c                  and (3) time
c     iprec()   - precision of (1) latitude, (2) longitude, (3) time
c     ip2()     - variable codes for variables in cast
c     ierror()  - whole profile error codes for each variable
c     
c     ipi()     - primary investigators information
c                   1. primary investigators
c                   2. for which variable
c
c     jsig2()   - number of significant figures in each secondary header variable
c     jprec2()  - precision of each secondary header variable
c     sechead() - secondary header variables
c
c     jsigb()   - number of significant figures in each biological variable
c     jprecb()  - precision of each biological variable
c     bio()     - biological data
c
c     depth()   - depth of each measurement
c
c     jtot2()   - number of bytes in each secondary header variable
c     jtotb()   - number of bytes in each biological variable
c
c     msig()    - number of significant figures in each measured variable at
c                  each level of measurement
c     mprec()   - precision of each measured variable at each
c                  level of measurement
c
c     mtot()    - number of digits in each measured variable at
c                  each level of measurement
c
c     temp()    - variable data at each level
c     iderror() - error flags for each variable at each depth level
c     iorigflag()- originators flags for each variable and depth
c
c     isec()    - variable codes for secondary header data
c     ibio()    - variable codes for biological data
c     parminf()  - variable specific information
c     jprecp()   - precision for variable specific information
c     jsigp()    - number of significant figures for variable specific
c                information
c     jtotp()    - number of digits in for variable specific information
c
c     itaxnum() - different taxonomic and biomass variable
c                  codes found in data
c     vtax()    - value of taxonomic variables and biomass variables
c
c     jsigtax() - number of significant figures in taxon values and
c                  biomass variables
c     jprectax()- precision of taxon values and biomass variables
c
c     jtottax() - number of bytes in taxon values 
c     itaxerr() - error codes for taxon data
c     itaxorigerr() - originators error codes for taxon data
c
c     nbothtot()- total number of taxa variables
c     stdz(40) - standard depth levels
c
c*******************************************************************

      integer isig(3), iprec(3), ip2(0:maxlevel), ierror(maxlevel),
     &        ipi(maxlevel,2)
      dimension jsig2(maxlevel),jprec2(maxlevel),sechead(maxlevel)
      dimension jsigb(maxlevel),jprecb(maxlevel),bio(maxlevel)
      dimension depth(maxlevel)
      dimension jtot2(maxlevel),jtotb(maxlevel)
      dimension msig(maxlevel,maxcalc), mprec(maxlevel,maxcalc)
      dimension mtot(maxlevel,maxcalc)
      dimension temp(maxlevel,maxcalc),iderror(maxlevel,0:maxcalc)
      dimension isec(maxlevel),ibio(maxlevel)
      dimension parminf(maxpinf,0:maxcalc),jsigp(maxpinf,0:maxcalc)
      dimension jprecp(maxpinf,0:maxcalc),jtotp(maxpinf,0:maxcalc)
      dimension iorigflag(maxlevel,0:maxcalc)
      dimension itaxnum(maxtcode,maxtax),vtax(0:maxtcode,maxtax)
      dimension jsigtax(maxtcode,maxtax),jprectax(maxtcode,maxtax)
      dimension jtottax(maxtcode,maxtax),itaxerr(maxtcode,maxtax)
      dimension itaxorigerr(maxtcode,maxtax)

      dimension stdz(kdim)

      common /thedata/ depth,temp
      common /flags/ ierror,iderror
      common /oflags/ iorigflag
      common /significant/ msig
      common /precision/ mprec
      common /totfigs/ mtot
      common /second/ jsig2,jprec2,jtot2,isec,sechead
      common /parminfo/ jsigp,jprecp,jtotp,parminf
      common /biology/ jsigb,jprecb,jtotb,ibio,bio
      common /taxon/ jsigtax,jprectax,jtottax,itaxerr,
     &     vtax,itaxnum,nbothtot,itaxorigerr
      
      data stdz/ 0., 10., 20., 30., 50., 75., 100., 125., 150.,
     &     200., 250., 300., 400., 500., 600., 700., 800., 900.,
     &     1000., 1100., 1200., 1300., 1400., 1500., 1750., 2000.,
     &     2500., 3000., 3500., 4000., 4500., 5000., 5500., 6000.,
     &     6500., 7000., 7500., 8000., 8500., 9000./
      
c**************************************************************
c
c     nf is the input file indentification number
c
c**************************************************************

      data nf/11/
      
c**************************************************************
c
c     Get user input file name from which casts will be
c     taken.  Open this file.
c
c**************************************************************

c User can modify the next section to read from a text file listing
c different input data files as a do-loop, for example, as opposed
c a single data input file.

      write(6,*)' '
      write(6,*)'Input File Name:'
      read(5,'(a80)') filename
      write(6,*)' '
      write(6,*)' '      

      open(nf,file=filename,status='old')
      
c**************************************************************
c
c   SUBROUTINE "WODread":  READS IN A SINGLE PROFILE FROM THE ASCII 
c                          FILE AND STORES THE DATA INTO ARRAYS
c   -------------------------------------------------------------------
c
c   Passed Variables:
c     
c     nf      - file identification number for input file
c     jj      - WOD cast number
c     cc      - NODC country code
c     icruise - NODC cruise number
c     iyear   - year of cast
c     month   - month of cast
c     iday    - day of cast
c     time    - time of cast
c     rlat    - latitude of cast
c     rlon    - longitude of cast
c     levels  - number of depth levels of data
c     istdlev - observed (0) or standard (1) levels
c     nvar   - number of variables recorded in cast
c     ip2(i)  - variable codes of variables in cast
c     nsecond - number of secondary header variables
c     nbio    - number of biological variables
c     isig()  - number of significant figures in (1) latitude, (2) longitude
c                and (3) time
c     iprec() - precision of (1) latitude, (2) longitude, (3) time
c     bmiss   - missing value marker
c     ieof    - set to one if end of file has been encountered
c
c   Common/Shared Variables and Arrays (see COMMON area of program):
c
c     depth(x)   - depth in meters (x = depth level)
c     temp(x,y)  - variable data (x = depth level, y = variable ID = ip2(i))
c                ... see also nvar, ip2, istdlev, levels above ...
c     sechead(i) - secondary header data (i = secondary header ID = isec(j))
c     isec(j)    - secondary header ID (j = #sequence (1st, 2nd, 3rd))
c                ... see also nsecond above ...
c     bio(i)     - biology header data (i = biol-header ID = ibio(j))
c     ibio(j)    - biology header ID (j = #sequence (1st, 2nd, 3rd))
c                ... see also nbio above ...
c     nbothtot   - number of taxa set / biomass variables
c     vtax(i,j)  - taxonomic/biomass array, where j = (1..nbothtot)
c                   For each entry (j=1..nbothtot), there are vtax(0,j)
c                   sub-entries.  [Note:  The number of sub-entries is 
c                   variable for each main entry.]  vtax also holds the
c                   value of the sub-entries.
c     itaxnum(i,j)- taxonomic code or sub-code 
c     chars       - WOD character data: 1. originators cruise code,
c                                       2. originators station code
c     npi        - number of PI codes
c     ipi        - Primary Investigator information
c                  1. primary investigator
c                  2. variable investigated
c     

c***************************************************************
      
      iVERSflag = 0             !- default is "WOD-2005"
      ieof = 0                  !- initialize end of file flag
      
      write(6,*)
     *  'Enter number of casts to view (0=view entire file)'
      read(5,*) numcasts
      if ( numcasts .eq. 0 ) numcasts=10000000

      do 50 ij=1,numcasts           !- MAIN CAST LOOP 
         
       chars(1)= '               '
       chars(2)= '               '
 
       
       if (iVERSflag .eq. 0) then
          
c     .   Read in as "WOD-2005" format
          
          call WODread200X(nf,jj,cc,icruise,iyear,month,iday,
     &         time,rlat,rlon,levels,istdlev,nvar,ip2,nsecond,nbio,
     &         isig,iprec,bmiss,ieof,chars,ipi,npi,iVERSflag)
          
c     .   ONLY happens if format rejected (rewind and try as WOD-1999)
          if (iVERSflag .gt. 0) then
             print *,' '
             print *,
     &            ' This data file is not in the WOD-2005 format.',
     &            '    Trying it in the WOD-1998 format.'
             print *,' '
             print *,' '
             rewind(nf) !- rewind file
          endif

       endif

       if (iVERSflag .eq. 1) then 
          
c     .   Read in as "WOD-1998" format          
          
c     .  WODread200X rejected format.  Must be WOD-1998 format.
          
          call WODread1998(nf,jj,cc,icruise,iyear,month,iday,
     &         time,rlat,rlon,levels,istdlev,nvar,ip2,nsecond,nbio,
     &         isig,iprec,bmiss,ieof,chars,ipi,npi)
          
       endif
       
       if ( ieof.gt.0 ) goto 4  !- Exit
       
c***************************************************************
c
c     STANDARD LEVELS OR OBSERVED LEVELS
c     ----------------------------------
c     
c     If this file is on standard levels (istdlev=1), program places the
c     standard depths in the depth() array (otherwise, observed depth values 
c     were read in and stored above by subroutine WODreadxxxx; where xxxx is
c     200X or 1998 depending on the data input format).
c
c***************************************************************

       if (istdlev .eq. 1 .and. ij .eq. 1) then
          
        do 60 i=1,kdim
         depth(i)=stdz(i)
 60     continue
          
       endif
 
c**************************************************************
c     
c     WRITE HEADER INFORMATION TO THE SCREEN
c     --------------------------------------
c
c     cc         - country code (a2)
c     icruise    - WOD cruise identifier (i8)
c     rlat       - latitude (f7.3)
c     rlon       - longitude (f8.3)
c     iyear      - year (i4)
c     month      - month (i2)
c     iday       - day (i2)
c     time       - time (GMT) (f7.2)
c     jj         - WOD cast identifier (i8)
c     levels     - number of depth levels measured (i6)
c     
c**************************************************************
       
 800   format(1x,a2,i8,1x,f7.3,1x,f8.3,1x,i4,1x,i2,1x,i2,
     &      1x,f7.2,1x,i8,1x,i6)

       write(6,*)
     &'----------------------------------------------------------'

       write(6,*) 'Output from ASCII file,  cast# ',ij
       
       write(6,*)
     &'----------------------------------------------------------'

       write(6,*)' '
      
       write(6,*)
     &'CC  cruise Latitde Longitde YYYY MM DD    Time    Cast  #levels'

       write(6,800) cc,icruise,rlat,rlon,iyear,month,iday,
     &      time,jj,levels
       
       write(6,*) ' '
       write(6,*) 'Number of variables in this cast: ',nvar
       write(6,*) ' '
       
c**************************************************************
c
c     WRITE CHARACTER DATA TO THE SCREEN
c     ----------------------------------
c     
c     chars(1)   - Originators cruise identifier
c     chars(2)   - Originators station identifier
c
c**************************************************************

       if ( ( chars(1)(1:1) .ne. ' ' ) .and.
     &     ( chars(1)(1:4) .ne. 'NONE' )) then
        write(6,*) 'Originators Cruise Code: ',chars(1)
       endif

       if ( ( chars(2)(1:1) .ne. ' ' ) .and.
     &     ( chars(2)(1:4) .ne. 'NONE' )) then
        write(6,*) 'Originators Station Code: ',chars(2)
       endif

       write(6,*) ' '

c***************************************************************
c
c     WRITE PRIMARY INVESTIGATOR INFORMATION TO THE SCREEN
c     ----------------------------------------------------
c
c     npi = number of primary investigator entries
c     ipi(1..npi,1)  - PI code
c     ipi(1..npi,2)  - variable for which PI was responsible
c     
c***************************************************************
       
       do 505 n=1,npi

          write(6,'(1x,a21,i5,1x,a20,i3)') 
     &         'Primary Investigator:',ipi(n,1),
     &         ' ... for variable #:',ipi(n,2)
          
 505   continue

       if ( npi .gt. 0 ) write(6,*) ' '


c**************************************************************
c
c     WRITE VARIABLE-CODE (column headings) TO THE SCREEN
c     ----------------------------------------------------
c     
c     nvar          - number of variables (1...nvar)
c     ip2(1..nvar)  - variable code for each variable present
c     
c     Example:  
c       For a cast with just Temperature[1], Oxygen[3], Pressure[25]:
c     
c       The variable sequence is ip2(1)=Temperature, ip2(2)=Oxygen, 
c          ip2(3)=Pressure
c     
c       nvar = 3
c
c       ip2(1) = 1, ip2(2) = 3, ip3(3) = 25
c
c     
c     Note:  If "nvar = 0", biology only cast.
c     
c**************************************************************

c     format(5x,1a,5x,10(3x,i2,11x))

801    format(5x,"z  fo",4x,10(i2,8x,"fo",3x))

       if (nvar .gt. 0) then
          
          write(6,801) (ip2(n),n=1,nvar)
          write(6,*)' '

c**************************************************************
c
c     WRITE DEPTH-DEPENDENT VARIABLE DATA TO THE SCREEN
c     --------------------------------------------------
c
c     Print depth (depth(n)), error flags for depth (iderror(n,0)),
c     each variable (temp(n,1..nvar)), and error flags for each
c     variables (iderror(n,1..nvar))
c
c**************************************************************

 802      format(f7.1,1x,i1,i1,14(f9.3,' (',i1,') ',i1,i1))
          
          do 80 n=1,levels
             write(6,802) depth(n),iderror(n,0),iorigflag(n,0),
     &            (temp(n,ip2(j)),msig(n,ip2(j)),
     &            iderror(n,ip2(j)),iorigflag(n,ip2(j)),j=1,nvar)
 80       continue
          
          write(6,*) ' ' 
          
c***************************************************************
c
c     PRINT ENTIRE-PROFILE ERROR FLAGS
c     ------------------------------------
c
c***************************************************************
          
 8021     format('VarFlag:    ',11x,11(i1,14x))
          
          write(6,8021)(ierror(ip2(j)),j=1,nvar)
          
          write(6,*) ' ' 
          
          
       endif                    !- "if (nvar .gt. 0) then"

c*************************************************************
c
c     WRITE SECONDARY-HEADER INFORMATION TO THE SCREEN
c     ---------------------------------------------
c
c     Print the secondary header code (isec(1..nsecond)) and the value
c     for that secondary header (sechead(isec(n))).
c     
c*************************************************************

 803  format(1x,'Secondary header #',i3,3x,f8.3,' (',i1,')')
 903  format(1x,'Secondary header #',i3,3x,f8.0,' (',i1,')')

      do 85 n = 1,nsecond

       if ( int(sechead(isec(n))) .lt. sechead(isec(n))) then
        write(6,803) isec(n), sechead(isec(n)),jsig2(isec(n))
       else
        write(6,903) isec(n), sechead(isec(n)),jsig2(isec(n))
       endif

85    continue
 
      write(6,*) ' ' 

c*************************************************************
c
c      WRITE VARIABLE SPECIFIC INFORMATION TO THE SCREEN
c      -------------------------------------------------
c
c*************************************************************

 814  format(1x,'Measured Variable #',i3,' Information Code #',i3,
     &       3x,f8.3,' (',i1,')')
 914  format(1x,'Measured Variable #',i3,' Information Code #',i3,
     &       3x,f8.0,' (',i1,')')

      do 86 j0=1,nvar

       j=ip2(j0)

       do 87 i=1,maxpinf

        if ( jtotp(i,j) .gt. 0 ) then

         if ( int(parminf(i,j)) .lt. parminf(i,j)) then
          write(6,814) j, i,parminf(i,j),jsigp(i,j)
         else
          write(6,914) j, i,parminf(i,j),jsigp(i,j)
         endif

         jtotp(i,j)=0

        endif

87     continue

86    continue

c*************************************************************
c
c     WRITE BIOLOGICAL HEADER INFORMATION TO THE SCREEN
c     ------------------------------------------
c     
c     Print the biology header code (ibio(1..nbio)) and the value
c     for that biology header (bio(ibio(n))).
c
c*************************************************************

 804  format(1x,'Biological header #',i3,3x,f13.3,' (',i1,')')

      
      do 90 n = 1,nbio
       write(6,804) ibio(n), bio(ibio(n)),jsigb(ibio(n))
       bio(ibio(n))=bmiss
       jsigb(ibio(n))=0
       ibio(n)=0
90    continue
      nbio=0

      write(6,*) ' ' 
 
c*************************************************************
c     
c     WRITE TAXA SET/BIOMASS VARIABLE INFORMATION TO THE SCREEN
c     ----------------------------------------------------------
c     
c     For each set/variable (1..nbothtot), print the set/variable code 
c     (ivtax = vtax(1,n)) and each member of that set (1..vtax(0,n)), 
c     where the sub-code is itaxnum(n2,n) and the sub-value is vtax(n2,n).
c
c*************************************************************
      
 805  format(5x,'   Code #',i4,3x,f13.3,' (',i1,') ',i1,i1)
      
      do 91 n = 1,nbothtot
         
       intax = vtax(0,n)     
       ivtax = vtax(1,n)

       if ( ivtax .lt. 0. .and. ivtax .gt. -501.) then
        write(6,'(a8,i3,1x,a25,i12," (",i1,")")') 'Taxa-set',n,
     &         ':  Biomass Parameter [1]#',ivtax,jsigtax(1,n)
       else
        write(6,'(a8,i3,1x,a22,4x,i10," (",i1,") ")') 'Taxa-set',n,
     &         ':  Taxonomic Code [1]#',ivtax,jsigtax(1,n)
       endif

       vtax(0,n)=0.
       vtax(1,n)=0.


       do 92 n2 = 2,intax
        write(6,805) itaxnum(n2,n), vtax(n2,n), jsigtax(n2,n),
     &               itaxerr(n2,n),itaxorigerr(n2,n)
        vtax(n2,n)=bmiss
        jsigtax(n2,n)=0
        itaxnum(n2,n)=0
 92    continue
       write(6,*)' '
 91   continue
      nbothtot=0
 
      write(6,*) ' ' 
      


50    continue !- End of MAIN LOOP
      
 4    continue !- EXIT 

      stop
      end

c----------------------------------------------------------------

      SUBROUTINE WODREAD200X(nf,jj,cc,icruise,iyear,month,iday,
     &     time,rlat,rlon,levels,isoor,nvar,ip2,nsecond,nbio,
     &     isig,iprec,bmiss,ieof,chars,ipi,npi,iVERSflag)
      
c     This subroutine reads in the WOD ASCII format and loads it
c     into arrays which are common/shared with the calling program.

c*****************************************************************
c
c   Passed Variables:
c
c     nf       - file identification number for input file
c     jj       - WOD cast number
c     cc       - NODC country code
c     icruise  - NODC cruise number
c     iyear    - year of cast
c     month    - month of cast
c     iday     - day of cast
c     time     - time of cast
c     rlat     - latitude of cast
c     rlon     - longitude of cast
c     levels   - number of depth levels of data
c     isoor    - observed (0) or standard (1) levels
c     nvar     - number of variables recorded in cast
c     ip2      - variable codes of variables in cast
c     nsecond  - number of secondary header variables
c     nbio     - number of biological variables
c     isig     - number of significant figures in (1) latitude, (2) longitude,
c                 and (3) time
c     iprec    - precision of (1) latitude, (2) longitude, (3) time
c     itotfig  - number of digits in (1) latitude, (2) longitude, (3) time
c     bmiss    - missing value marker
c     ieof     - set to one if end of file has been encountered
c     chars    - character data: 1=originators cruise code,
c                                2=originators station code
c     npi      - number of PI codes
c     ipi      - Primary Investigator information
c                  1. primary investigator
c                  2. variable investigated
c
c     iVERSflag  -  set to "1" if data are in WOD-1998 format. 
c                (subroutine exits so 1998 subroutine can be run)
c
c   Common/Shared Variables and Arrays (see COMMON area of program):
c
c     depth(x)   - depth in meters (x = depth level)
c     temp(x,y)  - variable data (x = depth level, y = variable ID = ip2(i))
c                ... see also nvar, ip2, istdlev, levels above ...
c     sechead(i) - secondary header data (i = secondary header ID = isec(j))
c     isec(j)    - secondary header ID (j = #sequence (1st, 2nd, 3rd))
c                ... see also nsecond above ...
c     bio(i)     - biology header data (i = biol-header ID = ibio(j))
c     ibio(j)    - biology header ID (j = #sequence (1st, 2nd, 3rd))
c                ... see also nbio above ...
c     nbothtot   - number of taxa set / biomass variables
c     vtax(i,j)  - taxonomic/biomass array, where j = (1..nbothtot)
c                   For each entry (j=1..nbothtot), there are vtax(0,j)
c                   sub-entries.  [Note:  The number of sub-entries is
c                   variable for each main entry.]  vtax also holds the
c                   value of the sub-entries.
c    itaxnum(i,j)- taxonomic code or sub-code
c    parminf(i,j)- variable specific information
c    origflag(i,j)- originators data flags
c
c***************************************************************


c******************************************************************
c
c   Parameters (constants):
c
c     maxlevel - maximum number of depth levels, also maximum
c                 number of all types of variables
c     maxcalc  - maximum number of measured and calculated
c                 depth dependent variables
c     maxtcode - maximum number of different taxa variable codes
c     maxtax   - maximum number of taxa sets
c     maxpinf - number of distinct variable specific information
c               variables
c
c******************************************************************

      parameter (maxlevel=30000, maxcalc=100)
      parameter (maxtcode=25, maxtax=2000, maxpinf=25)

c******************************************************************
c
c   Character Variables:
c
c     cc       - NODC country code
c     xchar    - dummy character array for reading in each 80
c                 character record
c     aout     - format specifier (used for FORTRAN I/O)
c     ichar    - cast character array
c     
c******************************************************************

      character*2  cc
      character*4  aout
      character*15 chars(2)
      character*80 xchar
      character*1500000 ichar
      
      data aout /'(iX)'/
      
c******************************************************************
c
c    Arrays:
c
c     isig     - number of significant figures in (1) latitude, (2) longitude,
c                 and (3) time
c     iprec    - precision of (1) latitude, (2) longitude, (3) time
c     itotfig  - number of digits in (1) latitude, (2) longitude, (3) time
c     ip2      - variable codes for variables in cast
c     ierror   - whole profile error codes for each variable
c     jsig2    - number of significant figures in each secondary header variable
c     jprec2   - precision of each secondary header variable
c     jtot2    - number of digits in each secondary header variable
c     sechead  - secondary header variables
c     jsigb    - number of significant figures in each biological variable
c     jprecb   - precision of each biological variable
c     jtotb    - number of digits in each biological variable
c     bio      - biological data
c     idsig    - number of significant figures in each depth measurement
c     idprec   - precision of each depth measurement
c     idtot    - number of digits in each depth measurement
c     depth    - depth of each measurement
c     msig     - number of significant figures in each measured variable at
c                 each level of measurement
c     mprec    - precision of each measured variable at each
c                 level of measurement
c     mtot     - number of digits in each measured variable at
c                 each level of measurement
c     temp     - variable data at each level
c     iderror  - error flags for each variable at each depth level
c     iorigflag- originators flags for each variable and depth
c     isec     - variable codes for secondary header data
c     ibio     - variable codes for biological data
c     parminf  - variable specific information
c     jprecp   - precision for variable specific information
c     jsigp    - number of significant figures for variable specific
c                information
c     jtotp    - number of digits in for variable specific information
c     itaxnum  - different taxonomic and biomass variable
c                 codes found in data
c     vtax     - value of taxonomic variables and biomass variables
c     jsigtax  - number of significant figures in taxon values and
c                 biomass variables
c     jprectax - precision of taxon values and biomass variables
c     jtottax  - number of digits in taxon values and biomass
c                 variables
c     itaxerr  - taxon variable error code
c     itaxorigerr - taxon originators variable error code
c     nbothtot - total number of taxa and biomass variables
c     ipi      - Primary investigator informationc
c                 1. primary investigator
c                 2. variable investigated
c
c*******************************************************************

      dimension isig(3), iprec(3), ip2(0:maxlevel), ierror(maxlevel)
      dimension itotfig(3),ipi(maxlevel,2)
      dimension jsig2(maxlevel), jprec2(maxlevel), sechead(maxlevel)
      dimension jsigb(maxlevel), jprecb(maxlevel), bio(maxlevel)
      dimension idsig(maxlevel),idprec(maxlevel), depth(maxlevel)
      dimension jtot2(maxlevel),jtotb(maxlevel),idtot(maxlevel)
      dimension msig(maxlevel,maxcalc), mprec(maxlevel,maxcalc)
      dimension mtot(maxlevel,maxcalc)
      dimension temp(maxlevel,maxcalc),iderror(maxlevel,0:maxcalc)
      dimension isec(maxlevel),ibio(maxlevel)
      dimension parminf(maxpinf,0:maxcalc),jsigp(maxpinf,0:maxcalc)
      dimension jprecp(maxpinf,0:maxcalc),jtotp(maxpinf,0:maxcalc)
      dimension iorigflag(maxlevel,0:maxcalc)
      dimension itaxnum(maxtcode,maxtax),vtax(0:maxtcode,maxtax)
      dimension jsigtax(maxtcode,maxtax),jprectax(maxtcode,maxtax)
      dimension jtottax(maxtcode,maxtax),itaxerr(maxtcode,maxtax)
      dimension itaxorigerr(maxtcode,maxtax)

c*******************************************************************
c     
c   Common Arrays and Variables:
c
c*******************************************************************
      
      common /thedata/ depth,temp
      common /flags/ ierror,iderror
      common /oflags/ iorigflag
      common /significant/ msig
      common /precision/ mprec
      common /totfigs/ mtot
      common /second/ jsig2,jprec2,jtot2,isec,sechead
      common /parminfo/ jsigp,jprecp,jtotp,parminf
      common /biology/ jsigb,jprecb,jtotb,ibio,bio
      common /taxon/ jsigtax,jprectax,jtottax,itaxerr,
     &     vtax,itaxnum,nbothtot,itaxorigerr
      
c******************************************************************
c     
c     Read in the first line of a cast into dummy character
c     variable xchar
c     
c
c     WOD-2005   First byte of each "cast record" is char "A".
c
c     WOD-1998   First byte of each "cast recond" is a number.
c
c******************************************************************

      read(nf,'(a80)',end=500) xchar

      if ( xchar(1:1) .ne. 'B' .and. xchar(1:1) .ne. 'A' ) then

         iVERSflag = 1 !- not WOD-2005 format, must be WOD-1998
         return

      else
         iVERSflag = 0 !- WOD-2005 format
      endif
      
c******************************************************************
c
c     The first seven characters of a cast contain the
c     number of characters which make up the entire cast.  Read
c     this number into nchar
c     
c******************************************************************

      read(xchar(2:2),'(i1)') inc
      write(aout(3:3),'(i1)') inc
      read(xchar(3:inc+2),aout) nchar

c******************************************************************
c
c     Place the first line of the cast into the cast holder
c     character array (ichar)
c
c******************************************************************

      ichar(1:80) = xchar

c******************************************************************
c
c     Calculate the number of full (all 80 characters contain information)
c     lines in this cast.  Subtract one since the first line was
c     already read in.
c
c******************************************************************

      nlines = nchar/80

c*****************************************************************
c
c     Read each line into the dummy variable
c
c*****************************************************************

      do 49 n0 = 2,nlines

       read(nf,'(a80)') xchar

c*****************************************************************
c
c     Place the line into the whole cast array
c
c*****************************************************************

       n = 80*(n0-1)+1
       ichar(n:n+79)=xchar

49    continue

c*****************************************************************
c
c     If there is a last line with partial information, read in
c     this last line and place it into the whole cast array
c
c*****************************************************************

      if ( nlines*80 .lt. nchar .and. nlines .gt. 0) then

       read(nf,'(a80)') xchar

       n = 80*nlines+1
       ichar(n:nchar) = xchar

      endif
       
c*****************************************************************
c
c   Extract header information from the cast array
c
c     jj       - WOD cast number  
c     cc       - NODC country code  
c     icruise  - NODC cruise number
c     iyear    - year of cast
c     month    - month of cast
c     iday     - day of cast
c
c*****************************************************************

      istartc=inc+3
      read(ichar(istartc:istartc),'(i1)') inc
      write(aout(3:3),'(i1)') inc
      read(ichar(istartc+1:istartc+inc),aout) jj
      istartc=istartc+inc+1

      cc = ichar(istartc:istartc+1)
      istartc=istartc+2

      read(ichar(istartc:istartc),'(i1)') inc
      write(aout(3:3),'(i1)') inc
      read(ichar(istartc+1:istartc+inc),aout) icruise
      istartc=istartc+inc+1

      read(ichar(istartc:istartc+3),'(i4)') iyear
      istartc=istartc+4
      read(ichar(istartc:istartc+1),'(i2)') month
      istartc=istartc+2
      read(ichar(istartc:istartc+1),'(i2)') iday
      istartc=istartc+2

c*****************************************************************
c
c   SUBROUTINE "charout":  READS IN AN WOD ASCII FLOATING-POINT
c                          VALUE SEQUENCE (i.e. # sig-figs,
c                          # total figs, precision, value itself).
c                          * THIS WILL BE CALLED TO EXTRACT MOST 
c   Examples:              FLOATING POINT VALUES IN THE WOD ASCII.
c
c     VALUE  Precision    WOD ASCII
c     -----  ---------    ---------
c     5.35       2        332535
c     5.         0        1105
c     15.357     3        55315357
c    (missing)            -
c
c   ---------------------------------------------------------------
c
c  Read in time of cast (time) using CHAROUT subroutine:
c
c     istartc  - position in character array to begin to read
c                 in data
c     isig     - number of digits in data value
c     iprec    - precision of data value
c     ichar    - character array from which to read data
c     time     - data value
c     bmiss    - missing value marker
c
c*****************************************************************

      call charout(istartc,isig(3),iprec(3),itotfig(3),ichar,
     *             time,bmiss)

c*****************************************************************
c
c     Read in latitude (rlat) and longitude (rlon) using CHAROUT:
c     
c        Negative latitude is south.
c        Negative longitude is west.
c     
c*****************************************************************

      call charout(istartc,isig(1),iprec(1),itotfig(3),ichar,
     *             rlat,bmiss)
      call charout(istartc,isig(2),iprec(2),itotfig(3),ichar,
     *             rlon,bmiss)

c*****************************************************************
c     
c     Read in the number of depth levels (levels) using CHAROUT:
c
c*****************************************************************

      read(ichar(istartc:istartc),'(i1)') inc
      write(aout(3:3),'(i1)') inc
      read(ichar(istartc+1:istartc+inc),aout) levels
      istartc=istartc+inc+1

c*****************************************************************
c
c     Read in whether data is on observed levels (isoor=0) or
c     standard levels (isoor=1)
c
c*****************************************************************

      read(ichar(istartc:istartc),'(i1)') isoor
      istartc=istartc+1

c*****************************************************************
c
c     Read in number of variables in cast
c
c*****************************************************************

      read(ichar(istartc:istartc+1),'(i2)') nvar
      istartc=istartc+2

c*****************************************************************
c
c     Read in the variable codes (ip2()), the whole profile
c       error flags (ierror(ip2())), and variable specific
c       information (iorigflag(,ip2()))
c
c*****************************************************************

      do 30 n = 1,nvar

       read(ichar(istartc:istartc),'(i1)') inc
       write(aout(3:3),'(i1)') inc
       read(ichar(istartc+1:istartc+inc),aout) ip2(n)
       istartc=istartc+inc+1

       read(ichar(istartc:istartc),'(i1)') ierror(ip2(n))
       istartc=istartc+1

       read(ichar(istartc:istartc),'(i1)') inc
       write(aout(3:3),'(i1)') inc
       read(ichar(istartc+1:istartc+inc),aout) npinf
       istartc=istartc+inc+1

       do 305 n2=1,npinf

        read(ichar(istartc:istartc),'(i1)') inc
        write(aout(3:3),'(i1)') inc
        read(ichar(istartc+1:istartc+inc),aout) nn
        istartc=istartc+inc+1

        call charout(istartc,jsigp(nn,ip2(n)),jprecp(nn,ip2(n)),
     &  jtotp(nn,ip2(n)),ichar, parminf(nn,ip2(n)),bmiss) 

305    continue

30    continue

c****************************************************************
c
c     Read in number of bytes in character data
c
c****************************************************************

      read(ichar(istartc:istartc),'(i1)') inc
      istartc=istartc+1

      npi=0
      chars(1)(1:4)='NONE'
      chars(2)(1:4)='NONE'

      if ( inc .gt. 0 ) then
       write(aout(3:3),'(i1)') inc
       read(ichar(istartc+1:istartc+inc),aout) inchad
       istartc=istartc+inc

c****************************************************************
c
c    Read in number of character and primary investigator arrays
c
c****************************************************************

      read(ichar(istartc:istartc),'(i1)') ica
      istartc=istartc+1

c****************************************************************
c
c    Read in character and primary investigator data
c      1 - originators cruise code
c      2 - originators station code
c      3 - primary investigators information
c
c****************************************************************

      do 45 nn=1,ica

       read(ichar(istartc:istartc),'(i1)') icn
       istartc=istartc+1

       if ( icn .lt. 3 ) then
        read(ichar(istartc:istartc+1),'(i2)') ns
        istartc=istartc+2
        chars(icn)= '               '
        chars(icn)= ichar(istartc:istartc+ns-1)
        istartc= istartc+ns
       else
        read(ichar(istartc:istartc+1),'(i2)') npi
        istartc=istartc+2
        do 505 n=1,npi
         read(ichar(istartc:istartc),'(i1)') inc
         write(aout(3:3),'(i1)') inc
         read(ichar(istartc+1:istartc+inc),aout) ipi(n,2)
         istartc=istartc+inc+1

         read(ichar(istartc:istartc),'(i1)') inc
         write(aout(3:3),'(i1)') inc
         read(ichar(istartc+1:istartc+inc),aout) ipi(n,1)
         istartc=istartc+inc+1
505     continue
       endif

45    continue

      endif

c****************************************************************
c
c     Read in number of bytes in secondary header variables
c
c****************************************************************

      read(ichar(istartc:istartc),'(i1)') inc
      istartc=istartc+1
      if ( inc .gt. 0 ) then
       write(aout(3:3),'(i1)') inc
       read(ichar(istartc+1:istartc+inc),aout) insec
       istartc=istartc+inc

c****************************************************************
c
c     Read in number of secondary header variables (nsecond)
c
c****************************************************************

       read(ichar(istartc:istartc),'(i1)') inc
       write(aout(3:3),'(i1)') inc
       read(ichar(istartc+1:istartc+inc),aout) nsecond
       istartc=istartc+inc+1

c****************************************************************
c
c     Read in secondary header variables (sechead())
c
c****************************************************************

       do 35 n = 1,nsecond

        read(ichar(istartc:istartc),'(i1)') inc
        write(aout(3:3),'(i1)') inc
        read(ichar(istartc+1:istartc+inc),aout) nn
        istartc=istartc+inc+1

        call charout(istartc,jsig2(nn),jprec2(nn),jtot2(nn),ichar,
     &  sechead(nn),bmiss) 

        isec(n) = nn

35     continue

       endif

c****************************************************************
c
c     Read in number of bytes in biology variables 
c
c****************************************************************

      read(ichar(istartc:istartc),'(i1)') inc
      istartc=istartc+1

      nbio=0
      inbio=0
      if ( inc .gt. 0 ) then
       write(aout(3:3),'(i1)') inc
       read(ichar(istartc+1:istartc+inc),aout) inbio
       istartc=istartc+inc

c****************************************************************
c
c     Read in number of biological variables (nbio)
c
c****************************************************************

      read(ichar(istartc:istartc),'(i1)') inc
      write(aout(3:3),'(i1)') inc
      read(ichar(istartc+1:istartc+inc),aout) nbio
      istartc=istartc+inc+1

c****************************************************************
c
c     Read in biological variables (bio())
c
c****************************************************************

      do 40 n = 1,nbio

       read(ichar(istartc:istartc),'(i1)') inc
       write(aout(3:3),'(i1)') inc
       read(ichar(istartc+1:istartc+inc),aout) nn
       istartc=istartc+inc+1

       call charout(istartc,jsigb(nn),jprecb(nn),jtotb(nn),ichar,
     & bio(nn),bmiss)

       ibio(n) = nn

40    continue

c****************************************************************
c
c     Read in biomass and taxonomic variables
c
c****************************************************************

      read(ichar(istartc:istartc),'(i1)') inc
      write(aout(3:3),'(i1)') inc
      read(ichar(istartc+1:istartc+inc),aout) nbothtot
      istartc=istartc+inc+1

      do 41 n = 1,nbothtot

       itaxtot=0
       read(ichar(istartc:istartc),'(i1)') inc
       write(aout(3:3),'(i1)') inc
       read(ichar(istartc+1:istartc+inc),aout) nn
       istartc=istartc+inc+1

       vtax(0,n)=nn

       do 42 n2 =1,nn

        itaxtot=itaxtot+1

        read(ichar(istartc:istartc),'(i1)') inc
        write(aout(3:3),'(i1)') inc
        read(ichar(istartc+1:istartc+inc),aout) itaxnum(itaxtot,n)
        istartc=istartc+inc+1
        call charout(istartc,jsigtax(itaxtot,n),jprectax(itaxtot,n),
     &   jtottax(itaxtot,n),ichar,vtax(itaxtot,n),bmiss)

        read(ichar(istartc:istartc),'(i1)') itaxerr(itaxtot,n)
        istartc=istartc+1
        read(ichar(istartc:istartc),'(i1)') itaxorigerr(itaxtot,n)
        istartc=istartc+1

42     continue

41    continue
      endif

c****************************************************************
c
c     Read in measured and calculated depth dependent variables
c       along with their individual reading flags
c
c****************************************************************

      do 50 n = 1,levels

       if ( isoor.eq.0 ) then

        call charout(istartc,idsig(n),idprec(n),idtot(n),ichar,
     & depth(n),bmiss)

        read(ichar(istartc:istartc),'(i1)') iderror(n,0)
        istartc=istartc+1
        read(ichar(istartc:istartc),'(i1)') iorigflag(n,0)
        istartc=istartc+1

       endif

       do 55 i = 1,nvar
     
        call charout(istartc,msig(n,ip2(i)),mprec(n,ip2(i)),
     & mtot(n,ip2(i)),ichar,temp(n,ip2(i)),bmiss)

       if ( temp(n,ip2(i)) .gt. bmiss ) then

       read(ichar(istartc:istartc),'(i1)') iderror(n,ip2(i))
       istartc=istartc+1
       read(ichar(istartc:istartc),'(i1)') iorigflag(n,ip2(i))
       istartc=istartc+1

       else
    
        iderror(n,ip2(i))=0
        iorigflag(n,ip2(1))=0
        msig(n,ip2(i))=0
        mprec(n,ip2(i))=0
        mtot(n,ip2(i))=0

       endif

55     continue

50     continue

       return

500    ieof = 1

       return
       end

C-----------------------------------------------

      SUBROUTINE WODREAD1998(nf,jj,cc,icruise,iyear,month,iday,
     &     time,rlat,rlon,levels,isoor,nvar,ip2,nsecond,nbio,
     &     isig,iprec,bmiss,ieof,chars,ipi,npi)
      
c     This subroutine reads in the WOD ASCII format and loads it
c     into arrays which are common/shared with the calling program.

c*****************************************************************
c
c   Passed Variables:
c
c     nf       - file identification number for input file
c     jj       - WOD cast number
c     cc       - NODC country code
c     icruise  - NODC cruise number
c     iyear    - year of cast
c     month    - month of cast
c     iday     - day of cast
c     time     - time of cast
c     rlat     - latitude of cast
c     rlon     - longitude of cast
c     levels   - number of depth levels of data
c     isoor    - observed (0) or standard (1) levels
c     nvar     - number of variables recorded in cast
c     ip2      - variable codes of variables in cast
c     nsecond  - number of secondary header variables
c     nbio     - number of biological variables
c     isig     - number of significant figures in (1) latitude, (2) longitude,
c                 and (3) time
c     iprec    - precision of (1) latitude, (2) longitude, (3) time
c     itotfig  - number of digits in (1) latitude, (2) longitude, (3) time
c     bmiss    - missing value marker
c     ieof     - set to one if end of file has been encountered
c     chars    - character data: 1=originators cruise code,
c                                2=originators station code
c     npi      - number of PI codes
c     ipi      - Primary Investigator information
c                  1. primary investigator
c                  2. variable investigated
c
c   Common/Shared Variables and Arrays (see COMMON area of program):
c
c     depth(x)   - depth in meters (x = depth level)
c     temp(x,y)  - variable data (x = depth level, y = variable ID = ip2(i))
c                ... see also nvar, ip2, istdlev, levels above ...
c     sechead(i) - secondary header data (i = secondary header ID = isec(j))
c     isec(j)    - secondary header ID (j = #sequence (1st, 2nd, 3rd))
c                ... see also nsecond above ...
c     bio(i)     - biology header data (i = biol-header ID = ibio(j))
c     ibio(j)    - biology header ID (j = #sequence(1st, 2nd, 3rd))
c                ... see also nbio above ...
c     nbothtot   - number of taxa set / biomass variables
c     vtax(i,j)  - taxonomic/biomass array, where j = (1..nbothtot)
c                   For each entry (j=1..nbothtot), there are vtax(0,j)
c                   sub-entries.  [Note:  The number of sub-entries is
c                   variable for each main entry.]  vtax also holds the
c                   value of the sub-entries.
c    itaxnum(i,j)- taxonomic code or sub-code
c
c***************************************************************


c******************************************************************
c
c   Parameters (constants):
c
c     maxlevel - maximum number of depth levels, also maximum
c                 number of all types of variables
c     maxcalc  - maximum number of measured and calculated
c                 depth dependent variables
c     maxtcode - maximum number of different taxa variable codes
c     maxtax   - maximum number of taxa sets
c
c******************************************************************

      parameter (maxlevel=30000, maxcalc=100)
      parameter (maxtcode=25, maxtax=2000)

c******************************************************************
c
c   Character Variables:
c
c     cc       - NODC country code
c     xchar    - dummy character array for reading in each 80
c                 character record
c     aout     - format specifier (used for FORTRAN I/O)
c     ichar    - cast character array
c     
c******************************************************************

      character*2  cc
      character*4  aout
      character*15 chars(2)
      character*80 xchar
      character*300000 ichar
      
      data aout /'(iX)'/
      
c******************************************************************
c
c    Arrays:
c
c     isig     - number of significant figures in (1) latitude, (2) longitude,
c                 and (3) time
c     iprec    - precision of (1) latitude, (2) longitude, (3) time
c     itotfig  - number of digits in (1) latitude, (2) longitude, (3) time
c     ip2      - variable codes for variables in cast
c     ierror   - whole profile error codes for each variable
c     jsig2    - number of significant figures in each secondary header variable
c     jprec2   - precision of each secondary header variable
c     jtot2    - number of digits in each secondary header variable
c     sechead  - secondary header variables
c     jsigb    - number of significant figures in each biological variable
c     jprecb   - precision of each biological variable
c     jtotb    - number of digits in each biological variable
c     bio      - biological data
c     idsig    - number of significant figures in each depth measurement
c     idprec   - precision of each depth measurement
c     idtot    - number of digits in each depth measurement
c     depth    - depth of each measurement
c     msig     - number of significant figures in each measured variable at
c                 each level of measurement
c     mprec    - precision of each measured variable at each
c                 level of measurement
c     mtot     - number of digits in each measured variable at
c                 each level of measurement
c     temp     - variable data at each level
c     iderror  - error flags for each variable at each depth level
c     isec     - variable codes for secondary header data
c     ibio     - variable codes for biological data
c     itaxnum  - different taxonomic and biomass variable
c                 codes found in data
c     vtax     - value of taxonomic variables and biomass variables
c     jsigtax  - number of significant figures in taxon values and
c                 biomass variables
c     jprectax - precision of taxon values and biomass variables
c     jtottax  - number of digits in taxon values and biomass
c                 variables
c     itaxerr  - taxon variable error code
c     nbothtot - total number of taxa and biomass variables
c     ipi      - Primary investigator informationc
c                 1. primary investigator
c                 2. variable investigated
c
c*******************************************************************

      dimension isig(3), iprec(3), ip2(0:maxlevel), ierror(maxlevel)
      dimension itotfig(3),ipi(maxlevel,2)
      dimension jsig2(maxlevel), jprec2(maxlevel), sechead(maxlevel)
      dimension jsigb(maxlevel), jprecb(maxlevel), bio(maxlevel)
      dimension idsig(maxlevel),idprec(maxlevel), depth(maxlevel)
      dimension jtot2(maxlevel),jtotb(maxlevel),idtot(maxlevel)
      dimension msig(maxlevel,maxcalc), mprec(maxlevel,maxcalc)
      dimension mtot(maxlevel,maxcalc)
      dimension temp(maxlevel,maxcalc),iderror(maxlevel,0:maxcalc)
      dimension isec(maxlevel),ibio(maxlevel)
      dimension itaxnum(maxtcode,maxtax),vtax(0:maxtcode,maxtax)
      dimension jsigtax(maxtcode,maxtax),jprectax(maxtcode,maxtax)
      dimension jtottax(maxtcode,maxtax),itaxerr(maxtcode,maxtax)
      dimension itaxorigerr(maxtcode,maxtax)

c*******************************************************************
c     
c   Common Arrays and Variables:
c
c*******************************************************************
      
      common /thedata/ depth,temp
      common /flags/ ierror,iderror
      common /significant/ msig
      common /precision/ mprec
      common /totfigs/ mtot
      common /second/ jsig2,jprec2,jtot2,isec,sechead
      common /biology/ jsigb,jprecb,jtotb,ibio,bio
      common /taxon/ jsigtax,jprectax,jtottax,itaxerr,
     &     vtax,itaxnum,nbothtot,itaxorigerr
      
c******************************************************************
c     
c     Read in the first line of a cast into dummy character
c     variable xchar
c     
c******************************************************************

      read(nf,'(a80)',end=500) xchar

c******************************************************************
c
c     The first seven characters of a cast contain the
c     number of characters which make up the entire cast.  Read
c     this number into nchar
c     
c******************************************************************

      read(xchar(1:1),'(i1)') inc
      write(aout(3:3),'(i1)') inc
      read(xchar(2:inc+1),aout) nchar

c******************************************************************
c
c     Place the first line of the cast into the cast holder
c     character array (ichar)
c
c******************************************************************

      ichar(1:80) = xchar

c******************************************************************
c
c     Calculate the number of full (all 80 characters contain information)
c     lines in this cast.  Subtract one since the first line was
c     already read in.
c
c******************************************************************

      nlines = nchar/80

c*****************************************************************
c
c     Read each line into the dummy variable
c
c*****************************************************************

      do 49 n0 = 2,nlines

       read(nf,'(a80)') xchar

c*****************************************************************
c
c     Place the line into the whole cast array
c
c*****************************************************************

       n = 80*(n0-1)+1
       ichar(n:n+79)=xchar

49    continue

c*****************************************************************
c
c     If there is a last line with partial information, read in
c     this last line and place it into the whole cast array
c
c*****************************************************************

      if ( nlines*80 .lt. nchar .and. nlines .gt. 0) then

       read(nf,'(a80)') xchar

       n = 80*nlines+1
       ichar(n:nchar) = xchar

      endif
       
c*****************************************************************
c
c   Extract header information from the cast array
c
c     jj       - WOD cast number  
c     cc       - NODC country code  
c     icruise  - NODC cruise number
c     iyear    - year of cast
c     month    - month of cast
c     iday     - day of cast
c
c*****************************************************************

      istartc=inc+2
      read(ichar(istartc:istartc),'(i1)') inc
      write(aout(3:3),'(i1)') inc
      read(ichar(istartc+1:istartc+inc),aout) jj
      istartc=istartc+inc+1

      cc = ichar(istartc:istartc+1)
      istartc=istartc+2

      read(ichar(istartc:istartc),'(i1)') inc
      write(aout(3:3),'(i1)') inc
      read(ichar(istartc+1:istartc+inc),aout) icruise
      istartc=istartc+inc+1

      read(ichar(istartc:istartc+3),'(i4)') iyear
      istartc=istartc+4
      read(ichar(istartc:istartc+1),'(i2)') month
      istartc=istartc+2
      read(ichar(istartc:istartc+1),'(i2)') iday
      istartc=istartc+2

c*****************************************************************
c
c   SUBROUTINE "charout":  READS IN AN WOD ASCII FLOATING-POINT
c                          VALUE SEQUENCE (i.e. # sig-figs,
c                          # total figs, precision, value itself).
c                          * THIS WILL BE CALLED TO EXTRACT MOST 
c   Examples:              FLOATING POINT VALUES IN THE WOD ASCII.
c
c     VALUE  Precision    WOD ASCII
c     -----  ---------    ---------
c     5.35       2        332535
c     5.         0        1105
c     15.357     3        55315357
c    (missing)            -
c
c   ---------------------------------------------------------------
c
c  Read in time of cast (time) using CHAROUT subroutine:
c
c     istartc  - position in character array to begin to read
c                 in data
c     isig     - number of digits in data value
c     iprec    - precision of data value
c     ichar    - character array from which to read data
c     time     - data value
c   -999.99    - missing value marker (bmiss)
c
c*****************************************************************

      call charout(istartc,isig(3),iprec(3),itotfig(3),ichar,
     *             time,bmiss)

c*****************************************************************
c
c     Read in latitude (rlat) and longitude (rlon) using CHAROUT:
c     
c        Negative latitude is south.
c        Negative longitude is west.
c     
c*****************************************************************

      call charout(istartc,isig(1),iprec(1),itotfig(3),ichar,
     *             rlat,bmiss)
      call charout(istartc,isig(2),iprec(2),itotfig(3),ichar,
     *             rlon,bmiss)

c*****************************************************************
c     
c     Read in the number of depth levels (levels) using CHAROUT:
c
c*****************************************************************

      read(ichar(istartc:istartc),'(i1)') inc
      write(aout(3:3),'(i1)') inc
      read(ichar(istartc+1:istartc+inc),aout) levels
      istartc=istartc+inc+1

c*****************************************************************
c
c     Read in whether data is on observed levels (isoor=0) or
c     standard levels (isoor=1)
c
c*****************************************************************

      read(ichar(istartc:istartc),'(i1)') isoor
      istartc=istartc+1

c*****************************************************************
c
c     Read in number of variables in cast
c
c*****************************************************************

      read(ichar(istartc:istartc+1),'(i2)') nvar
      istartc=istartc+2

c*****************************************************************
c
c     Read in the variable codes (ip2()) and the whole profile
c       error flags (ierror(ip2()))
c
c*****************************************************************

      do 30 n = 1,nvar

       read(ichar(istartc:istartc),'(i1)') inc
       write(aout(3:3),'(i1)') inc
       read(ichar(istartc+1:istartc+inc),aout) ip2(n)
       istartc=istartc+inc+1

       read(ichar(istartc:istartc),'(i1)') ierror(ip2(n))
       istartc=istartc+1

30    continue

c****************************************************************
c
c     Read in number of bytes in character data
c
c****************************************************************

      read(ichar(istartc:istartc),'(i1)') inc
      istartc=istartc+1
      if ( inc .gt. 0 ) then
       write(aout(3:3),'(i1)') inc
       read(ichar(istartc+1:istartc+inc),aout) inchad
       istartc=istartc+inc

c****************************************************************
c
c    Read in number of character and primary investigator arrays
c
c****************************************************************

      npi=0
      chars(1)(1:4)='NONE'
      chars(2)(1:4)='NONE'
      read(ichar(istartc:istartc),'(i1)') ica
      istartc=istartc+1

c****************************************************************
c
c    Read in character and primary investigator data
c      1 - originators cruise code
c      2 - originators station code
c      3 - primary investigators information
c
c****************************************************************

      do 45 nn=1,ica

       read(ichar(istartc:istartc),'(i1)') icn
       istartc=istartc+1

       if ( icn .lt. 3 ) then
        read(ichar(istartc:istartc+1),'(i2)') ns
        istartc=istartc+2
        chars(icn)= '               '
        chars(icn)= ichar(istartc:istartc+ns-1)
        istartc= istartc+ns
       else
        read(ichar(istartc:istartc+1),'(i2)') npi
        istartc=istartc+2
        do 505 n=1,npi
         read(ichar(istartc:istartc),'(i1)') inc
         write(aout(3:3),'(i1)') inc
         read(ichar(istartc+1:istartc+inc),aout) ipi(n,2)
         istartc=istartc+inc+1

         read(ichar(istartc:istartc),'(i1)') inc
         write(aout(3:3),'(i1)') inc
         read(ichar(istartc+1:istartc+inc),aout) ipi(n,1)
         istartc=istartc+inc+1
505     continue
       endif

45    continue

      endif

c****************************************************************
c
c     Read in number of bytes in secondary header variables
c
c****************************************************************

      read(ichar(istartc:istartc),'(i1)') inc
      istartc=istartc+1
      if ( inc .gt. 0 ) then
       write(aout(3:3),'(i1)') inc
       read(ichar(istartc+1:istartc+inc),aout) insec
       istartc=istartc+inc

c****************************************************************
c
c     Read in number of secondary header variables (nsecond)
c
c****************************************************************

       read(ichar(istartc:istartc),'(i1)') inc
       write(aout(3:3),'(i1)') inc
       read(ichar(istartc+1:istartc+inc),aout) nsecond
       istartc=istartc+inc+1

c****************************************************************
c
c     Read in secondary header variables (sechead())
c
c****************************************************************

       do 35 n = 1,nsecond

        read(ichar(istartc:istartc),'(i1)') inc
        write(aout(3:3),'(i1)') inc
        read(ichar(istartc+1:istartc+inc),aout) nn
        istartc=istartc+inc+1

        call charout(istartc,jsig2(nn),jprec2(nn),jtot2(nn),ichar,
     &  sechead(nn),bmiss) 

        isec(n) = nn

35     continue

       endif

c****************************************************************
c
c     Read in number of bytes in biology variables 
c
c****************************************************************

      read(ichar(istartc:istartc),'(i1)') inc
      istartc=istartc+1

      if ( inc .gt. 0 ) then
       write(aout(3:3),'(i1)') inc
       read(ichar(istartc+1:istartc+inc),aout) inbio
       istartc=istartc+inc

c****************************************************************
c
c     Read in number of biological variables (nbio)
c
c****************************************************************

      read(ichar(istartc:istartc),'(i1)') inc
      write(aout(3:3),'(i1)') inc
      read(ichar(istartc+1:istartc+inc),aout) nbio
      istartc=istartc+inc+1

c****************************************************************
c
c     Read in biological variables (bio())
c
c****************************************************************

      do 40 n = 1,nbio

       read(ichar(istartc:istartc),'(i1)') inc
       write(aout(3:3),'(i1)') inc
       read(ichar(istartc+1:istartc+inc),aout) nn
       istartc=istartc+inc+1

       call charout(istartc,jsigb(nn),jprecb(nn),jtotb(nn),ichar,
     & bio(nn),bmiss)

       ibio(n) = nn

40    continue

c****************************************************************
c
c     Read in biomass and taxonomic variables
c
c****************************************************************

      read(ichar(istartc:istartc),'(i1)') inc
      write(aout(3:3),'(i1)') inc
      read(ichar(istartc+1:istartc+inc),aout) nbothtot
      istartc=istartc+inc+1

      do 41 n = 1,nbothtot

       itaxtot=0
       read(ichar(istartc:istartc),'(i1)') inc
       write(aout(3:3),'(i1)') inc
       read(ichar(istartc+1:istartc+inc),aout) nn
       istartc=istartc+inc+1

       vtax(0,n)=nn

       do 42 n2 =1,nn

        itaxtot=itaxtot+1

        read(ichar(istartc:istartc),'(i1)') inc
        write(aout(3:3),'(i1)') inc
        read(ichar(istartc+1:istartc+inc),aout) itaxnum(itaxtot,n)
        istartc=istartc+inc+1
        call charout(istartc,jsigtax(itaxtot,n),jprectax(itaxtot,n),
     &   jtottax(itaxtot,n),ichar,vtax(itaxtot,n),bmiss)

        read(ichar(istartc:istartc),'(i1)') itaxerr(itaxtot,n)
        istartc=istartc+1

42     continue

41    continue
      endif

c****************************************************************
c
c     Read in measured and calculated depth dependent variables
c       along with their individual reading flags
c
c****************************************************************

      do 50 n = 1,levels

       if ( isoor.eq.0 ) then

        call charout(istartc,idsig(n),idprec(n),idtot(n),ichar,
     & depth(n),bmiss)

        read(ichar(istartc:istartc),'(i1)') iderror(n,0)
        istartc=istartc+1

       endif

       do 55 i = 1,nvar
     
        call charout(istartc,msig(n,ip2(i)),mprec(n,ip2(i)),
     & mtot(n,ip2(i)),ichar,temp(n,ip2(i)),bmiss)

       if ( temp(n,ip2(i)) .gt. bmiss ) then

       read(ichar(istartc:istartc),'(i1)') iderror(n,ip2(i))
       istartc=istartc+1

       else
    
        iderror(n,ip2(i))=0

       endif

55     continue

50     continue

       return

500    ieof = 1

       return
       end

C------------------------------------------------------------------
      SUBROUTINE CHAROUT(istartc,jsig,jprec,jtot,ichar,value,bmiss)
      
c     This subroutine reads a single real value from the
c     WOD ASCII format.  This value consists of four
c     components:  # significant figures, # total figures,
c     precision, and the value. 
     
c   Examples:

c     VALUE  Precision    WOD ASCII
c     -----  ---------    ---------
c     5.35       2        332535
c     5.         0        1105
c     15.357     3        55315357
c    (missing)            -           
     
c******************************************************
c     
c   Passed Variables:
c
c     istartc    - starting point to read in data
c     jsig       - number of significant figures in data value
c     jprec      - precision of data value
c     jtot       - number of figures in data value
c     ichar      - character array from which to read data
c     value      - data value
c     bmiss      - missing value marker
c
c*****************************************************

c*****************************************************
c
c   Character Array:
c
c     cwriter    - format statement (FORTRAN I/O)
c
c****************************************************

      character*6 cwriter
      character*(*) ichar
      
      data cwriter /'(fX.X)'/
      
c****************************************************
c     
c     Check if this is a missing value (number of 
c       figures = '-')
c
c****************************************************

      if ( ichar(istartc:istartc) .eq. '-' ) then

       istartc = istartc+1
       value = bmiss
       return

      endif
       
c****************************************************
c
c     Read in number of significant figure, total
c       figures and precision of value
c
c****************************************************

      read(ichar(istartc:istartc),'(i1)') jsig
      read(ichar(istartc+1:istartc+1),'(i1)') jtot
      read(ichar(istartc+2:istartc+2),'(i1)') jprec
      istartc=istartc+3

c****************************************************
c
c     Write these values into a FORTRAN format statement
c
c       e.g. "553" --> '(f5.3)'
c            "332" --> '(f3.2)'
c
c****************************************************

      write(cwriter(3:3),'(i1)') jtot
      write(cwriter(5:5),'(i1)') jprec

c****************************************************
c
c     Read in the data value using thhe FORTRAN 
c       format statement created above (cwriter).
c
c****************************************************

      read(ichar(istartc:istartc+jtot-1),cwriter) value

c****************************************************
c
c     Update the character array position (pointer)
c       and send it back to the calling program.
c
c****************************************************

      istartc=istartc+jtot

      return
      end

c <next few lines under version control, do not edit>
c $URL$
c $Id$
c $Revision$
c $Date$
