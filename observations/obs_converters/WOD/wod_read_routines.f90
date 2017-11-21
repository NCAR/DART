! This code is not protected by the DART copyright agreement.
! DART $Id$

     module WOD_read_routines_mod      

!    calls the subroutine WODread (WODread200X if the data file are in WOD05
!    or WOD01 format or the subroutine WODread1998 if the data are in WOD98
!    format). These subroutines do the actual reading of the ASCII format,
!    and load the data into arrays which are passed back to the main program.

!    It is intended that the subroutine WODread provides an example of how
!    to extract the data and variables from the ASCII format,
!    whereas the main part of the wodFOR program provides an example of how
!    these data can be made accessible/workable as a series of arrays.
!    
!    Comments and suggestions for improving this program would be appreciated.
!    Updates to the World Ocean Data 2005 data and to this program will be posted
!    in the NODC/WOD web site at http://www.nodc.noaa.gov
      
!***********************************************************
!    
!    Missing values used in this dataset = bmiss = -999.99
!    
!***********************************************************
!     
!   Parameters (constants):
!     
!     maxlevel  - maximum number of depth levels, also maximum
!                   number of all types of variables
!     maxcalc   - maximum number of measured and calculated
!                   depth dependent variables
!     kdim      - number of standard depth levels
!     bmiss     - binary missing value marker
!     maxtcode  - maximum number of different taxa variable codes
!     maxtax    - maximum number of taxa sets
!     maxpinf   - maximum number of distinct measured variable
!                 information codes
!
!******************************************************************
      
implicit none
private

      public :: WODreadDART

      integer, public :: maxlevel, maxcalc, kdim, maxtcode, maxtax, maxpinf
      real, public :: bmiss

      parameter (maxlevel=30000, maxcalc=100)
      parameter (kdim=40, bmiss=-999.99)
      parameter (maxtcode=25, maxtax=2000)
      parameter (maxpinf=25)
      
      real,    public :: depth(maxlevel), temp(maxlevel,maxcalc)
      integer, public :: ierror(maxlevel)
      integer, public :: iderror(maxlevel,0:maxcalc)
      integer, public :: iorigflag(maxlevel,0:maxcalc)
      integer, public :: idsig(maxlevel),idprec(maxlevel)
      integer, public :: idtot(maxlevel), itotfig(3)
      integer, public :: msig(maxlevel,maxcalc), mprec(maxlevel,maxcalc)
      integer, public :: mtot(maxlevel,maxcalc)
      integer, public :: jsig2(maxlevel),jprec2(maxlevel)
      real,    public :: sechead(maxlevel)
      integer, public :: jsigp(maxpinf,0:maxcalc)
      real,    public :: parminf(maxpinf,0:maxcalc)
      integer, public :: isec(maxlevel),ibio(maxlevel)
      real,    public :: bio(maxlevel)
      integer, public :: jsigb(maxlevel),jprecb(maxlevel)
      integer, public :: jtot2(maxlevel),jtotb(maxlevel)
      integer, public :: jprecp(maxpinf,0:maxcalc),jtotp(maxpinf,0:maxcalc)
      real,    public :: vtax(0:maxtcode,maxtax)
      integer, public :: itaxnum(maxtcode,maxtax)
      integer, public :: jsigtax(maxtcode,maxtax),jprectax(maxtcode,maxtax)
      integer, public :: jtottax(maxtcode,maxtax),itaxerr(maxtcode,maxtax)
      integer, public :: itaxorigerr(maxtcode,maxtax)

      character*2  cc
      character*4  aout
      character*15 chars(2)
      character*80 xchar
      character*1500000 ichar

      integer :: isig(3), iprec(3)
      integer :: ipi(maxlevel,2)
      
      real :: stdz(kdim)

      ! standard levels.  the actual levels are either in the file,
      ! or there's a flag saying use these ones.
      data stdz/ &
             0.,   10.,   20.,   30.,   50.,   75.,  100.,  125., &
           150.,  200.,  250.,  300.,  400.,  500.,  600.,  700., &
           800.,  900., 1000., 1100., 1200., 1300., 1400., 1500., &
          1750., 2000., 2500., 3000., 3500., 4000., 4500., 5000., &
          5500., 6000., 6500., 7000., 7500., 8000., 8500., 9000./


      data aout /'(iX)'/

      contains

!----------------------------------------------------------------

! original interface:
!
!      SUBROUTINE WODREAD200X(nf,jj,cc,icruise,iyear,month,iday, &
!          time,rlat,rlon,levels,isoor,nvar,ip2,nsecond,nbio, &
!          isig,iprec,bmiss,ieof,chars,ipi,npi,iVERSflag)
      
        subroutine WODreadDART(nf,iyear,month,iday, &              
                          time,rlat,rlon,levels,istdlev,nvar,ip2,nsecond, &
                          bmiss,castid,ieof)

!     This subroutine reads in the WOD ASCII format and loads it
!     into arrays which are common/shared with the calling program.

!*****************************************************************
!
!   Passed Variables:
!
!     nf       - file identification number for input file
!     castid    - WOD cast number
!     !cc       - NODC country code
!     !icruise  - NODC cruise number
!     iyear    - year of cast
!     month    - month of cast
!     iday     - day of cast
!     time     - time of cast
!     rlat     - latitude of cast
!     rlon     - longitude of cast
!     levels   - number of depth levels of data
!     istdlev  - observed (0) or standard (1) levels
!     nvar     - number of variables recorded in cast
!     ip2      - variable codes of variables in cast
!     nsecond  - number of secondary header variables
!     !nbio     - number of biological variables
!     !isig     - number of significant figures in (1) latitude, (2) longitude,
!     !            and (3) time
!     !iprec    - precision of (1) latitude, (2) longitude, (3) time
!     !itotfig  - number of digits in (1) latitude, (2) longitude, (3) time
!     bmiss    - missing value marker
!     ieof     - set to one if end of file has been encountered
!     !chars    - character data: 1=originators cruise code,
!     !                           2=originators station code
!     !npi      - number of PI codes
!     !ipi      - Primary Investigator information
!     !             1. primary investigator
!     !             2. variable investigated
!
!     !iVERSflag  -  set to "1" if data are in WOD-1998 format. 
!     !           (subroutine exits so 1998 subroutine can be run)
!
!   Common/Shared Variables and Arrays (see COMMON area of program):
!
!     depth(x)   - depth in meters (x = depth level)
!     temp(x,y)  - variable data (x = depth level, y = variable ID = ip2(i))
!                ... see also nvar, ip2, istdlev, levels above ...
!     sechead(i) - secondary header data (i = secondary header ID = isec(j))
!     isec(j)    - secondary header ID (j = #sequence (1st, 2nd, 3rd))
!                ... see also nsecond above ...
!     bio(i)     - biology header data (i = biol-header ID = ibio(j))
!     ibio(j)    - biology header ID (j = #sequence (1st, 2nd, 3rd))
!                ... see also nbio above ...
!     nbothtot   - number of taxa set / biomass variables
!     vtax(i,j)  - taxonomic/biomass array, where j = (1..nbothtot)
!                   For each entry (j=1..nbothtot), there are vtax(0,j)
!                   sub-entries.  [Note:  The number of sub-entries is
!                   variable for each main entry.]  vtax also holds the
!                   value of the sub-entries.
!    itaxnum(i,j)- taxonomic code or sub-code
!    parminf(i,j)- variable specific information
!    origflag(i,j)- originators data flags
!
!***************************************************************


!******************************************************************
!
!   Parameters (constants):
!
!     maxlevel - maximum number of depth levels, also maximum
!                 number of all types of variables
!     maxcalc  - maximum number of measured and calculated
!                 depth dependent variables
!     maxtcode - maximum number of different taxa variable codes
!     maxtax   - maximum number of taxa sets
!     maxpinf - number of distinct variable specific information
!               variables
!
!******************************************************************

      !parameter (maxlevel=30000, maxcalc=100)
      !parameter (maxtcode=25, maxtax=2000, maxpinf=25)

!******************************************************************
!
!   Character Variables:
!
!     cc       - NODC country code
!     xchar    - dummy character array for reading in each 80
!                 character record
!     aout     - format specifier (used for FORTRAN I/O)
!     ichar    - cast character array
!     
!******************************************************************
      integer :: nf,iyear,month,iday,levels,istdlev,nvar,ip2(0:maxlevel)
      integer :: nsecond, ieof
      real    :: time,rlat,rlon,bmiss

      
!******************************************************************
!
!    Arrays:
!
!     isig     - number of significant figures in (1) latitude, (2) longitude,
!                 and (3) time
!     iprec    - precision of (1) latitude, (2) longitude, (3) time
!     itotfig  - number of digits in (1) latitude, (2) longitude, (3) time
!     ip2      - variable codes for variables in cast
!     ierror   - whole profile error codes for each variable
!     jsig2    - number of significant figures in each secondary header variable
!     jprec2   - precision of each secondary header variable
!     jtot2    - number of digits in each secondary header variable
!     sechead  - secondary header variables
!     jsigb    - number of significant figures in each biological variable
!     jprecb   - precision of each biological variable
!     jtotb    - number of digits in each biological variable
!     bio      - biological data
!     idsig    - number of significant figures in each depth measurement
!     idprec   - precision of each depth measurement
!     idtot    - number of digits in each depth measurement
!     depth    - depth of each measurement
!     msig     - number of significant figures in each measured variable at
!                 each level of measurement
!     mprec    - precision of each measured variable at each
!                 level of measurement
!     mtot     - number of digits in each measured variable at
!                 each level of measurement
!     temp     - variable data at each level
!     iderror  - error flags for each variable at each depth level
!     iorigflag- originators flags for each variable and depth
!     isec     - variable codes for secondary header data
!     ibio     - variable codes for biological data
!     parminf  - variable specific information
!     jprecp   - precision for variable specific information
!     jsigp    - number of significant figures for variable specific
!                information
!     jtotp    - number of digits in for variable specific information
!     itaxnum  - different taxonomic and biomass variable
!                 codes found in data
!     vtax     - value of taxonomic variables and biomass variables
!     jsigtax  - number of significant figures in taxon values and
!                 biomass variables
!     jprectax - precision of taxon values and biomass variables
!     jtottax  - number of digits in taxon values and biomass
!                 variables
!     itaxerr  - taxon variable error code
!     itaxorigerr - taxon originators variable error code
!     nbothtot - total number of taxa and biomass variables
!     ipi      - Primary investigator informationc
!                 1. primary investigator
!                 2. variable investigated
!
!*******************************************************************

      !dimension isig(3), iprec(3)
      !dimension ipi(maxlevel,2)
      !dimension ip2(0:maxlevel), ierror(maxlevel)
      !dimension itotfig(3),ipi(maxlevel,2)
      !dimension jsig2(maxlevel), jprec2(maxlevel), sechead(maxlevel)
      !dimension jsigb(maxlevel), jprecb(maxlevel), bio(maxlevel)
      !dimension idsig(maxlevel),idprec(maxlevel), depth(maxlevel)
      !dimension jtot2(maxlevel),jtotb(maxlevel),idtot(maxlevel)
      !dimension msig(maxlevel,maxcalc), mprec(maxlevel,maxcalc)
      !dimension mtot(maxlevel,maxcalc)
      !dimension temp(maxlevel,maxcalc),iderror(maxlevel,0:maxcalc)
      !dimension isec(maxlevel),ibio(maxlevel)
      !dimension parminf(maxpinf,0:maxcalc),jsigp(maxpinf,0:maxcalc)
      !dimension jprecp(maxpinf,0:maxcalc),jtotp(maxpinf,0:maxcalc)
      !dimension iorigflag(maxlevel,0:maxcalc)
      !dimension itaxnum(maxtcode,maxtax),vtax(0:maxtcode,maxtax)
      !dimension jsigtax(maxtcode,maxtax),jprectax(maxtcode,maxtax)
      !dimension jtottax(maxtcode,maxtax),itaxerr(maxtcode,maxtax)
      !dimension itaxorigerr(maxtcode,maxtax)

      integer :: n, nn, n0, n2, i, castid
      integer :: nbio, inc, nchar, nlines, istartc, icruise, npinf
      integer :: npi, inchad, ica, icn, ns, insec, inbio, nbothtot, itaxtot

!*******************************************************************
!     
!   Common Arrays and Variables:
!
!*******************************************************************
      
      !common /thedata/ depth,temp
      !common /flags/ ierror,iderror
      !common /oflags/ iorigflag
      !common /significant/ msig
      !common /precision/ mprec
      !common /totfigs/ mtot
      !common /second/ jsig2,jprec2,jtot2,isec,sechead
      !common /parminfo/ jsigp,jprecp,jtotp,parminf
      !common /biology/ jsigb,jprecb,jtotb,ibio,bio
      !common /taxon/ jsigtax,jprectax,jtottax,itaxerr, &
      !    vtax,itaxnum,nbothtot,itaxorigerr
      

!******************************************************************
!     
!     Read in the first line of a cast into dummy character
!     variable xchar
!     
!
!     WOD-2005   First byte of each "cast record" is char "A".
!
!     WOD-1998   First byte of each "cast recond" is a number.
!
!******************************************************************

!print *, 'in library read routine'

      read(nf,'(a80)',end=500) xchar

      if ( xchar(1:1) .ne. 'B' .and. xchar(1:1) .ne. 'A' .and. &
           xchar(1:1) .ne. 'C' ) then

         !iVERSflag = 1 !- not WOD-2005 format, must be WOD-1998
         write(6, *) 'file not in WOD-2005 format; cannot be read'
         stop
         !return

      else
         if ( xchar(1:1) .eq. 'C' ) then
          !iVERSflag=2   !- WOD-2013 format
         else
          !iVERSflag = 0 !- WOD-2005 format
         endif
      endif
      
!print *, 'read record'

!******************************************************************
!
!     The first seven characters of a cast contain the
!     number of characters which make up the entire cast.  Read
!     this number into nchar
!     
!******************************************************************

      read(xchar(2:2),'(i1)') inc
      write(aout(3:3),'(i1)') inc
      read(xchar(3:inc+2),aout) nchar

!******************************************************************
!
!     Place the first line of the cast into the cast holder
!     character array (ichar)
!
!******************************************************************

      ichar(1:80) = xchar

!******************************************************************
!
!     Calculate the number of full (all 80 characters contain information)
!     lines in this cast.  Subtract one since the first line was
!     already read in.
!
!******************************************************************

      nlines = nchar/80

!print *, 'nlines this cast = ', nlines
!*****************************************************************
!
!     Read each line into the dummy variable
!
!*****************************************************************

      do 49 n0 = 2,nlines

       read(nf,'(a80)') xchar

!*****************************************************************
!
!     Place the line into the whole cast array
!
!*****************************************************************

       n = 80*(n0-1)+1
       ichar(n:n+79)=xchar

49    continue

!*****************************************************************
!
!     If there is a last line with partial information, read in
!     this last line and place it into the whole cast array
!
!*****************************************************************

      if ( nlines*80 .lt. nchar .and. nlines .gt. 0) then

       read(nf,'(a80)') xchar

       n = 80*nlines+1
       ichar(n:nchar) = xchar

      endif
       
!*****************************************************************
!
!   Extract header information from the cast array
!
!     castid   - WOD cast number  
!     cc       - NODC country code  
!     icruise  - NODC cruise number
!     iyear    - year of cast
!     month    - month of cast
!     iday     - day of cast
!
!*****************************************************************
!print *, 'all lines in, starting unpack'

      istartc=inc+3
      read(ichar(istartc:istartc),'(i1)') inc
      write(aout(3:3),'(i1)') inc
      read(ichar(istartc+1:istartc+inc),aout) castid
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

!*****************************************************************
!
!   SUBROUTINE "charout":  READS IN AN WOD ASCII FLOATING-POINT
!                          VALUE SEQUENCE (i.e. # sig-figs,
!                          # total figs, precision, value itself).
!                          * THIS WILL BE CALLED TO EXTRACT MOST 
!   Examples:              FLOATING POINT VALUES IN THE WOD ASCII.
!
!     VALUE  Precision    WOD ASCII
!     -----  ---------    ---------
!     5.35       2        332535
!     5.         0        1105
!     15.357     3        55315357
!    (missing)            -
!
!   ---------------------------------------------------------------
!
!  Read in time of cast (time) using CHAROUT subroutine:
!
!     istartc  - position in character array to begin to read
!                 in data
!     isig     - number of digits in data value
!     iprec    - precision of data value
!     ichar    - character array from which to read data
!     time     - data value
!     bmiss    - missing value marker
!
!*****************************************************************

      call charout(istartc,isig(3),iprec(3),itotfig(3),ichar, &
                  time,bmiss)

!print *, '1.'
!*****************************************************************
!
!     Read in latitude (rlat) and longitude (rlon) using CHAROUT:
!     
!        Negative latitude is south.
!        Negative longitude is west.
!     
!*****************************************************************

      call charout(istartc,isig(1),iprec(1),itotfig(3),ichar, &
                  rlat,bmiss)
      call charout(istartc,isig(2),iprec(2),itotfig(3),ichar, &
                  rlon,bmiss)

!print *, '2.'
!*****************************************************************
!     
!     Read in the number of depth levels (levels) using CHAROUT:
!
!*****************************************************************

      read(ichar(istartc:istartc),'(i1)') inc
      write(aout(3:3),'(i1)') inc
      read(ichar(istartc+1:istartc+inc),aout) levels
      istartc=istartc+inc+1

!print *, '3.'
!*****************************************************************
!
!     Read in whether data is on observed levels (istdlev=0) or
!     standard levels (istdlev=1)
!
!*****************************************************************

      read(ichar(istartc:istartc),'(i1)') istdlev
      istartc=istartc+1

!print *, '4.'
!*****************************************************************
!
!     Read in number of variables in cast
!
!*****************************************************************

      read(ichar(istartc:istartc+1),'(i2)') nvar
      istartc=istartc+2

!print *, '5.'
!*****************************************************************
!
!     Read in the variable codes (ip2()), the whole profile
!       error flags (ierror(ip2())), and variable specific
!       information (iorigflag(,ip2()))
!
!*****************************************************************

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

        call charout(istartc,jsigp(nn,ip2(n)),jprecp(nn,ip2(n)), &
       jtotp(nn,ip2(n)),ichar, parminf(nn,ip2(n)),bmiss) 

305    continue

30    continue

!print *, '6.'
!****************************************************************
!
!     Read in number of bytes in character data
!
!****************************************************************

      read(ichar(istartc:istartc),'(i1)') inc
      istartc=istartc+1

      npi=0
      chars(1)(1:4)='NONE'
      chars(2)(1:4)='NONE'

      if ( inc .gt. 0 ) then
       write(aout(3:3),'(i1)') inc
       read(ichar(istartc+1:istartc+inc),aout) inchad
       istartc=istartc+inc

!****************************************************************
!
!    Read in number of character and primary investigator arrays
!
!****************************************************************

      read(ichar(istartc:istartc),'(i1)') ica
      istartc=istartc+1

!****************************************************************
!
!    Read in character and primary investigator data
!      1 - originators cruise code
!      2 - originators station code
!      3 - primary investigators information
!
!****************************************************************

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

!print *, '7.'
!****************************************************************
!
!     Read in number of bytes in secondary header variables
!
!****************************************************************

      read(ichar(istartc:istartc),'(i1)') inc
      istartc=istartc+1
      if ( inc .gt. 0 ) then
       write(aout(3:3),'(i1)') inc
       read(ichar(istartc+1:istartc+inc),aout) insec
       istartc=istartc+inc

!****************************************************************
!
!     Read in number of secondary header variables (nsecond)
!
!****************************************************************

       read(ichar(istartc:istartc),'(i1)') inc
       write(aout(3:3),'(i1)') inc
       read(ichar(istartc+1:istartc+inc),aout) nsecond
       istartc=istartc+inc+1

!****************************************************************
!
!     Read in secondary header variables (sechead())
!
!****************************************************************

       do 35 n = 1,nsecond

        read(ichar(istartc:istartc),'(i1)') inc
        write(aout(3:3),'(i1)') inc
        read(ichar(istartc+1:istartc+inc),aout) nn
        istartc=istartc+inc+1

        call charout(istartc,jsig2(nn),jprec2(nn),jtot2(nn),ichar, &
       sechead(nn),bmiss) 

        isec(n) = nn

35     continue

       endif

!print *, '8.'
!****************************************************************
!
!     Read in number of bytes in biology variables 
!
!****************************************************************

      read(ichar(istartc:istartc),'(i1)') inc
      istartc=istartc+1

      nbio=0
      inbio=0
      if ( inc .gt. 0 ) then
       write(aout(3:3),'(i1)') inc
       read(ichar(istartc+1:istartc+inc),aout) inbio
       istartc=istartc+inc

!****************************************************************
!
!     Read in number of biological variables (nbio)
!
!****************************************************************

      read(ichar(istartc:istartc),'(i1)') inc
      write(aout(3:3),'(i1)') inc
      read(ichar(istartc+1:istartc+inc),aout) nbio
      istartc=istartc+inc+1

!****************************************************************
!
!     Read in biological variables (bio())
!
!****************************************************************

      do 40 n = 1,nbio

       read(ichar(istartc:istartc),'(i1)') inc
       write(aout(3:3),'(i1)') inc
       read(ichar(istartc+1:istartc+inc),aout) nn
       istartc=istartc+inc+1

       call charout(istartc,jsigb(nn),jprecb(nn),jtotb(nn),ichar, &
      bio(nn),bmiss)

       ibio(n) = nn

40    continue

!print *, '9.'
!****************************************************************
!
!     Read in biomass and taxonomic variables
!
!****************************************************************

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
        call charout(istartc,jsigtax(itaxtot,n),jprectax(itaxtot,n), &
        jtottax(itaxtot,n),ichar,vtax(itaxtot,n),bmiss)

        read(ichar(istartc:istartc),'(i1)') itaxerr(itaxtot,n)
        istartc=istartc+1
        read(ichar(istartc:istartc),'(i1)') itaxorigerr(itaxtot,n)
        istartc=istartc+1

42     continue

41    continue
      endif

!print *, '10.'
!****************************************************************
!
!     Read in measured and calculated depth dependent variables
!       along with their individual reading flags
!
!****************************************************************

      do 50 n = 1,levels

       if ( istdlev.eq.0 ) then

        call charout(istartc,idsig(n),idprec(n),idtot(n),ichar, &
                     depth(n),bmiss)

        read(ichar(istartc:istartc),'(i1)') iderror(n,0)
        istartc=istartc+1
        read(ichar(istartc:istartc),'(i1)') iorigflag(n,0)
        istartc=istartc+1

       else
        depth(n) = stdz(n)
       endif

       do 55 i = 1,nvar
     
        call charout(istartc,msig(n,ip2(i)),mprec(n,ip2(i)), &
                     mtot(n,ip2(i)),ichar,temp(n,ip2(i)),bmiss)

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

       ieof = 0
!print *, '11.'
       return

500    ieof = 1

!print *, '12.'
       return
       end subroutine


!------------------------------------------------------------------
      SUBROUTINE CHAROUT(istartc,jsig,jprec,jtot,ichar,value,bmiss)
      
!     This subroutine reads a single real value from the
!     WOD ASCII format.  This value consists of four
!     components:  # significant figures, # total figures,
!     precision, and the value. 
     
!   Examples:

!     VALUE  Precision    WOD ASCII
!     -----  ---------    ---------
!     5.35       2        332535
!     5.         0        1105
!     15.357     3        55315357
!    (missing)            -           
     
!******************************************************
!     
!   Passed Variables:
!
!     istartc    - starting point to read in data
!     jsig       - number of significant figures in data value
!     jprec      - precision of data value
!     jtot       - number of figures in data value
!     ichar      - character array from which to read data
!     value      - data value
!     bmiss      - missing value marker
!
!*****************************************************
     integer :: istartc, jsig, jprec, jtot
     character*500000 :: ichar
     real :: value, bmiss

!*****************************************************
!
!   Character Array:
!
!     cwriter    - format statement (FORTRAN I/O)
!
!****************************************************

      character*6 cwriter
      !character*(*) ichar
      
      data cwriter /'(fX.X)'/
      
!****************************************************
!     
!     Check if this is a missing value (number of 
!       figures = '-')
!
!****************************************************

      if ( ichar(istartc:istartc) .eq. '-' ) then

       istartc = istartc+1
       value = bmiss
       return

      endif
       
!****************************************************
!
!     Read in number of significant figure, total
!       figures and precision of value
!
!****************************************************

      read(ichar(istartc:istartc),'(i1)') jsig
      read(ichar(istartc+1:istartc+1),'(i1)') jtot
      read(ichar(istartc+2:istartc+2),'(i1)') jprec
      istartc=istartc+3

!****************************************************
!
!     Write these values into a FORTRAN format statement
!
!       e.g. "553" --> '(f5.3)'
!            "332" --> '(f3.2)'
!
!****************************************************

      write(cwriter(3:3),'(i1)') jtot
      write(cwriter(5:5),'(i1)') jprec

!****************************************************
!
!     Read in the data value using thhe FORTRAN 
!       format statement created above (cwriter).
!
!****************************************************

      read(ichar(istartc:istartc+jtot-1),cwriter) value

!****************************************************
!
!     Update the character array position (pointer)
!       and send it back to the calling program.
!
!****************************************************

      istartc=istartc+jtot
      return
      end subroutine

      end module WOD_read_routines_mod      

! <next few lines under version control, do not edit>
! $URL$
! $Id$
! $Revision$
! $Date$
