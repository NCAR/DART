c     Data Assimilation Research Testbed -- DART
c     Copyright 2004-2007, Data Assimilation Research Section
c     University Corporation for Atmospheric Research
c     Licensed under the GPL -- www.gpl.org/licenses/gpl.html

      program read_bufr

c     <next few lines under version control, do not edit>
c     $URL$
c     $Id$
c     $Revision$
c     $Date$

c     This is the 4nd version of the BUFR preparation program for the
c     newly revised DART version. Ps and moisture obs are output as well.
c     06/22/2004.  the time is output
c     
c     A bug with the radiosonde T events selection j2 is fixed at 04/06/2005.
c     and surface data are removed from this version.
c     For radiosonde T: T of Pc=1 and 6 is dry T; T of pc=8 is virtual T.
c     For Aircraft T: No moistuerobs and all T are dry T
c     For Surface obs, the T of pc=1 is already virtual T.
c
        REAL*8     R8BFMS
        PARAMETER ( R8BFMS = 10.0E10 )
C                                      "Missing" value for BUFR data
        PARAMETER ( MXR8PM = 10 )
C                                      Maximum number of BUFR parameters
        PARAMETER ( MXR8LV = 400 )
C                                      Maximum number of BUFR levels
        PARAMETER ( MXR8VN = 10 )
C                                      Maximum number of BUFR event sequences
        PARAMETER ( MXR8VT = 6 )
C                                      Maximum number of BUFR variable types
        PARAMETER ( MXSTRL = 80 )
C                                      Maximum size of a string
C
        REAL*8 hdr (MXR8PM), evns(MXR8PM, MXR8LV, MXR8VN, MXR8VT)
C
        COMMON /PREPBC/ hdr, evns, nlev

        PARAMETER (STRLN = 180)
        CHARACTER  outstg*(200), subset*8, inf*19 , outf*20
        CHARACTER  var(MXR8VT)/'P','Q','T','Z','U','V'/
C
        PARAMETER  ( NFILO = 15 )
        INTEGER  iunso ( NFILO )
     +          /   51,   52,   53,   54,   55,
     +              56,   57,   58,   59,   60,
     +              61,   62,   63,   64,   65  /
        CHARACTER*6     filo ( NFILO )
     +          / 'ADPUPA', 'AIRCAR', 'AIRCFT', 'SATWND', 'PROFLR',
     +            'VADWND', 'SATBOG', 'SATEMP', 'ADPSFC', 'SFCSHP',
     +            'SFCBOG', 'SPSSMI', 'SYNDAT', 'ERS1DA', 'GOESND'  /

      dimension tdata(8), udata(8), vdata(8), qdata(8), pdata(8)
      integer wtype, ptype, qtype, ttype
      integer pc_t, pc_q, pc_u, pc_p
      integer tqm, pqm, qqm, uqm, vqm
C
        LOGICAL  found

C----------------------------------
c    for calculation of saturated water vapor for qobs error
      real*8  l0

      eps= .622
      es0= 6.11          ! in mb, qsat pressure at 273.16K
      cl = 4187.
      cpv= 1876.5
      rv = 461.5         ! m2/s2/K
      l0 = 2.501e6
      t0 = 273.16
      fact1 = (cpv - cl) / rv
      fact2 = (l0 + (cl - cpv) * t0) / rv
      fact3 = 1. / t0
      omeps = 1.-eps
c----------------------------------

        inf='prepqm'
        outf='prepqm.out'

        iprof = 0
        lunobs=210
        open(unit=lunobs, file = outf, form='formatted')

        ibufr_unit = 11 
        OPEN  ( UNIT = ibufr_unit, FILE = inf, FORM = 'UNFORMATTED' )
        CALL OPENBF  ( ibufr_unit, 'IN', ibufr_unit )
        CALL DATELEN  ( 10 )
C
C      Get the next station report from the input file.
C
  10    CALL READPB  ( ibufr_unit, subset, idate, ierrpb )
 
        idate00 = idate/100
        hour01 = idate - idate00*100
        if(hour01 .eq. 0.0 ) hour01 = 24
c       print*, 'idate= ', idate, idate00, hour01

        IF ( ierrpb .eq. -1 )  STOP

        iprof = iprof + 1    !! station number
C
C      Set the appropriate output file unit number.
C
        ii = 1
        found = .false.
        DO WHILE  ( ( .not. found ) .and. ( ii .le. NFILO ) )
          IF( subset(1:6) .eq. filo ( ii ) )  THEN
                found = .true.
                iuno = iunso ( ii )

c     OPEN (UNIT = iunso(ii),FILE= 'readpb.out.'//filo(ii) )
      OPEN (UNIT = iuno,FILE= 'adiag.'//filo(ii) )

C      Print the HDR data for this station report. 
cliu
       time0 = hdr(4)
       stat_elv = hdr(5)

c  07/06/2004 for output of exact time
c       print*, 'idate= ', idate, idate00, hour01, time0
         hour01 =  hour01 + time0
c  07/06/2004 for output of exact time

c   modified 04/07/2005
      if(hdr(6) .eq. 120.0 .or. hdr(6) .eq. 130.0 .or. 
     &   hdr(6) .eq. 131.0 .or. hdr(6) .eq. 132.0 .or. 
     &   hdr(6) .eq. 133.0 .or. 
     &   hdr(6) .eq. 220.0 .or. hdr(6) .eq. 221.0 .or. 
     &   hdr(6) .eq. 230.0 .or. hdr(6) .eq. 231.0 .or. 
     &   hdr(6) .eq. 232.0 .or. hdr(6) .eq. 233.0 .or. 
     &   hdr(6) .eq. 242.0 .or. hdr(6) .eq. 243.0 .or. 
     &   hdr(6) .eq. 245.0 .or. hdr(6) .eq. 246.0 .or.
     &   hdr(6) .eq. 252.0 .or. hdr(6) .eq. 253.0 .or. 
     &   hdr(6) .eq. 255.0 ) then   
c
c liu
        WRITE  ( UNIT = iuno, FMT = '("#", 132("-") )' )
        WRITE  ( UNIT = iuno, FMT = '("#", A3,6x,a6,1x,a6,3x,
     +             a5,5x,a4,3x,a5,4x,a4,3x,a4,1x,a3,1x,a3,2x,a6,
     +             6x,a3,6x,a3,4x,a5,1x,a8,2x,a7,6x,a3 )' )
     +           'SID', 'XOB', 'YOB','DHR','ELV','TYP','T29','ITP',
     +           'lev','var','OB','qm', 'pc', 'rc', 'fc', 'an', 'err'

        WRITE  ( UNIT = iuno, FMT = '("#", 132("-") )' )
c liu
       endif
c liu

         ELSE
                ii = ii + 1
         END IF
        END DO

        IF ( ( .not. found ) .and. ( ierrpb .eq. 0 ) )  GO TO 10

C      Print the EVNS data for this station report.
c
       if( hdr(2) .ge. 360.0) hdr(2)=hdr(2) - 360.0
       if( hdr(2) .lt.   0.0) hdr(2)=hdr(2) + 360.0

       d2r = 3.1415926/180.0
       alon = hdr(2)*d2r
       alat = hdr(3)*d2r
       time = hdr(4)

       tdata(2) = alon
       tdata(3) = alat
       tdata(6) = 0
       tdata(7) = iprof + 1000000
       tdata(8) = hour01
          ttype = hdr(6)

       pdata(2) = alon
       pdata(3) = alat
       pdata(6) = 0
       pdata(7) = iprof + 3000000
       pdata(8) = hour01
          ptype = hdr(6)

       qdata(2) = alon
       qdata(3) = alat
       qdata(6) = 0
       qdata(7) = iprof + 5000000
       qdata(8) = hour01
          qtype = hdr(6)

       udata(2) = alon
       udata(3) = alat
       udata(7) = iprof + 2000000
       udata(8) = hour01           !! time

c      U, V are merged to T array and has same type with T in BUFR file
      if(subset(1:6) .eq. 'SATWND' .or. hdr(6) .eq. 221.0 .or. 
     &  hdr(6) .eq. 280.0 .or. hdr(6) .eq. 230.0 .or. 
     &                         hdr(6) .eq. 231.0) then
       wtype=hdr(6)
       else
       wtype=hdr(6)+100  
       endif   
c     !! wind type = t type in other categry.

       vdata(2) = alon
       vdata(3) = alat
       vdata(7) = iprof + 9000000
       vdata(8) = hour01

c   modified 04/07/2005
      if(hdr(6) .eq. 120.0 .or. hdr(6) .eq. 130.0 .or. 
     &   hdr(6) .eq. 131.0 .or. hdr(6) .eq. 132.0 .or. 
     &   hdr(6) .eq. 133.0 .or. 
     &   hdr(6) .eq. 220.0 .or. hdr(6) .eq. 221.0 .or. 
     &   hdr(6) .eq. 230.0 .or. hdr(6) .eq. 231.0 .or. 
     &   hdr(6) .eq. 232.0 .or. hdr(6) .eq. 233.0 .or. 
     &   hdr(6) .eq. 242.0 .or. hdr(6) .eq. 243.0 .or. 
     &   hdr(6) .eq. 245.0 .or. hdr(6) .eq. 246.0 .or.
     &   hdr(6) .eq. 252.0 .or. hdr(6) .eq. 253.0 .or. 
     &   hdr(6) .eq. 255.0 ) then   
c
       DO 200 lv = 1, nlev

c   find out the event before virtural T and Q check step (pc=8.0)
        j2t = 1
        do ievent =1, MXR8VN
        if(evns(3,lv,ievent, 3) .lt. 8.0) then
         j2t = ievent
         go to 201
        endif
        enddo
  201  continue

        j2q = 1
        do ievent =1, MXR8VN
        if(evns(3,lv,ievent, 2) .lt. 8.0) then
         j2q = ievent
         go to 202
        endif
        enddo
  202  continue

         DO kk = 1, MXR8VT

c     for T & Q, the events (pc>=8) are virtual T step. The events (pc=6, 1) are dry T
c     set qm of j2t/j2q event according to 1st event
c
          if(subset(1:6) .eq. 'ADPUPA' .and. var(kk) .eq. 'T') then
           if(j2t.ne.1) evns(7,lv, j2t, kk) = evns(7,lv, 1, kk) 
           jj = j2t
          else if(subset(1:6) .eq. 'ADPUPA' .and. var(kk) .eq. 'Q') then
           if(j2q.ne.1) evns(7,lv, j2q, kk) = evns(7,lv, 1, kk) 
           jj = j2q
          else
           jj = 1                              
c                   ! select 1st event for other
          endif
c
        WRITE( UNIT = outstg, FMT = '( A8, 1X,  
     &         2F7.2, 1X, F7.3, 1X, 2F8.1, 1X, F7.1, 1X, 
     &         F6.1, I4, 1X, A2, 7(1X,F8.1) )' )
     &     (hdr(ii), ii =1, 8),lv, var(kk), (evns(ii,lv,jj,kk), ii=1,7)

               DO mm = 1, 200 
                IF  ( outstg (mm:mm) .eq. '*' ) outstg (mm:mm) = ' '
               END DO

               IF ( outstg (90:154) .ne. ' ') THEN
                WRITE (UNIT = iuno, FMT= '(A150)' ) outstg
               ENDIF
  
         END DO
cliu
c-------------------------------------------
       if(subset(1:6).eq.'ADPUPA') then
c            use the j2t event of pc =6 or 1 if there is T report
         toe = evns ( 7, lv, j2t, 3)
         tob = evns ( 1, lv, j2t, 3) + 273.16
         tqm = evns ( 2, lv, j2t, 3) 
         pc_t = evns (3, lv, j2t, 3)

        else
c            use the 1st event
         toe = evns ( 7, lv, 1, 3)
         tob = evns ( 1, lv, 1, 3) + 273.16
         tqm = evns ( 2, lv, 1, 3) 
         pc_t = evns (3, lv, 1, 3)
       endif

          if(tqm .eq. 0) toe = toe*0.9
          if(tqm .eq. 3) toe = toe*1.2

c-------------------------------------------
c    for moisture obs from ADPUPA, use j2q event 
c     to decimal error for RH (0.2) and mg/kg to g/kg (CAM needs kg/kg)
       if(subset(1:6).eq.'ADPUPA') then
         qoe = evns(7, lv, j2q, 2) * 0.1       
         qob = evns(1, lv, j2q, 2) * 1.0e-3    
         qqm = evns(2, lv, j2q, 2) 
        pc_q = evns(3, lv, j2q, 2)
        else
         qoe = evns(7, lv, 1, 2) * 0.1       
         qob = evns(1, lv, 1, 2) * 1.0e-3    
         qqm = evns(2, lv, 1, 2) 
        pc_q = evns(3, lv, 1, 2)
       endif
          if(qqm .eq. 0) qoe = qoe*0.9
          if(qqm .eq. 3) qoe = qoe*1.2
c
c     For Pressure       ! mb (need to convert to Pa in CAM)
         poe = evns ( 7, lv, 1, 1)     
         pob = evns ( 1, lv, 1, 1)     
         pqm = evns ( 2, lv, 1, 1) 
         pc_p = evns(3, lv, 1, 1)

          if(pqm .eq. 0) poe = poe*0.9
          if(pqm .eq. 3) poe = poe*1.2
         ppb = evns ( 1, lv, 1, 1)            ! mb
         zob = evns ( 1, lv, 1, 4)            ! m

c    for wind
         uoe = evns (7, lv, 1, 5)
         uob = evns (1, lv, 1, 5)
         uqm = evns (2, lv, 1, 5) 
         pc_u = evns(3, lv, 1, 5)
          if(uqm .eq. 0) uoe = uoe*0.9
          if(uqm .eq. 3) uoe = uoe*1.2

         voe = evns ( 7, lv, 1, 6)
         vob = evns ( 1, lv, 1, 6) 
         vqm = evns ( 2, lv, 1, 6) 
          if(vqm .eq. 0) voe = voe*0.9
          if(vqm .eq. 3) voe = voe*1.2

c----------------------------------------------

c   write out temperature observation of ADPUPA, AIRCAR, AIRCFT.
c
      if(subset(1:6) .eq. 'ADPUPA' .or. subset(1:6) .eq. 'AIRCAR' .or.
     &                                  subset(1:6) .eq. 'AIRCFT' )then
        if(toe .lt. 1.e9 .and. tob .lt. 1.e9 ) then
         if(tqm .lt.4 .and. pqm .lt. 4) then

          tdata(1) = toe
          tdata(4) = ppb
          tdata(5) = tob

          write(lunobs, 800) tdata, ttype, tqm,subset(1:6), pc_t
         endif
        endif
       endif
c
c   write out moisture observation from ADPUPA.
c

      if(subset(1:6) .eq. 'ADPUPA' ) then
c
c   compute the error of specific moisture from the fixed RH obs error(20%)
c
      es   = es0 * (tob / t0) ** fact1 * exp( fact2 * (fact3 - 1./tob))
      qsat = eps * es / (pob - omeps * es)
      qoe  = qoe * qsat * 1000.0                 ! to g/kg
c   keep the error > 0.1 g/kg as in NCEP SSI.
      qoe = max(0.1, qoe)

c       print*, 'es= ', es, tob-273.16, qoe

        if(qoe .lt. 1.e9 .and. qob .lt. 1.e9 ) then
         if(qqm .lt.4 .and. pqm .lt. 4) then

        if( tqm .ge. 4) then       ! the T obs is not good for use to calculate qoe
        print*, 'bad = ', qoe
        qoe = 1.0e10
        endif

         qdata(1) = qoe
         qdata(4) = ppb
         qdata(5) = qob

         if(qoe .lt. 9.9) then      !! skip few large qoe obs.
         write(lunobs, 800) qdata, qtype, qqm,subset(1:6), pc_q
         endif

         endif
        endif
       endif

c   write out Ps observation of ADPUPA
c
       if(subset(1:6).eq.'ADPUPA' ) then
        if(poe .lt. 1.e9 .and. pob .lt. 1.e9 ) then
         if(pqm .lt.4 .and. stat_elv .eq. zob) then

         pdata(1) = poe
         pdata(4) = zob
         pdata(5) = pob

         write(lunobs, 800) pdata, ptype, pqm,subset(1:6),pc_p
         endif
        endif
       endif
c
c   write out wind observation of ADPUPA, AIRcraft, SATWND
c
       if(subset(1:6) .eq.'ADPUPA' .or. subset(1:6) .eq. 'AIRCAR' .or.
     &    subset(1:6) .eq.'AIRCFT' .or. subset(1:6) .eq. 'SATWND' )then
        if(uoe .lt. 1.e9 .and. uob .lt. 1.e9 ) then
         if(uqm .lt.4 .and. pqm .lt. 4) then

         udata(1) = uoe
         udata(4) = ppb
         udata(5) = uob
         udata(6) = vob

         vdata(1) = voe
         vdata(4) = ppb
         vdata(5) = vob
         vdata(6) = vob

          write(lunobs, 800) udata, wtype, uqm,subset(1:6),pc_u
          write(lunobs, 800) vdata, wtype, vqm,subset(1:6),pc_u
         endif
        endif
       endif
cliu

  200  continue

       endif

        IF ( ierrpb .eq. 0 )  GO TO 10
c
c 800  format(f4.2, 2f7.3, e12.5, f7.2, f7.2, f9.0, f7.3, f5.0, i3)
  800  format(f4.2,2f7.3,e12.5,f7.2,f7.2,f9.0,f7.3,i4,i2,1x,a6,i2)
        STOP
        END
C-----------------------------------------------------------------------
        SUBROUTINE READPB  ( lunit, subset, idate, iret )
C
C      This subroutine will read and combine the mass and wind subsets
C      of the next station report in the prepbufr file.  It is styled
C      after entry point READNS, and it only requires the prepbufr file
C      to be opened for reading with OPENBF.  The combined station
C      report is returned to the caller in COMMON /PREPBC/.
C      This common area contains the number of levels in the report,
C      a one dimensional array with the header information, and a four
C      dimensional array containing all events from the variables POB,
C      QOB, TOB, ZOB, UOB, and VOB for the report.
C
C      The header array contains the following list of mnemonics:
C
C      SID XOB YOB DHR ELV TYP T29 ITP
C
C      The 4-D array of data, EVNS ( ii, lv, jj, kk ), is indexed
C      as follows:
C
C      "ii" indexes the event data types; these consist of:
C          1) OBservation
C          2) Quality Mark
C          3) Program Code
C          4) Reason Code
C          5) ForeCast value
C          6) ANalysed value
C          7) observation error
C      "lv" indexes the levels of the report
C      "jj" indexes the event stacks
C      "kk" indexes the variable types (p,q,t,z,u,v)
C
C      Note that the structure of this array is identical to one
C      returned from UFBEVN, with an additional (4th) dimension to
C      include the six variable types into the same array.
C
C      The return codes are as follows:
C      iret =  0 - normal return
C           =  1 - the station report within COMMON /PREPBC/ contains the
C                  last available subset from within the prepbufr file
C           = -1 - there are no more subsets available from within the
C                  prepbufr file       
C
        INCLUDE  'readpb.prm'
        CHARACTER*(*)   subset
        CHARACTER*(MXSTRL)      head
     +          / 'SID XOB YOB DHR ELV TYP T29 ITP' /
Cliu
        CHARACTER*(MXSTRL)      ostr ( MXR8VT )
     +          / 'POB PQM PPC PRC PFC PAN POE',
     +            'QOB QQM QPC QRC QFC QAN QOE',
     +            'TOB TQM TPC TRC TFC TAN TOE',
     +            'ZOB ZQM ZPC ZRC ZFC ZAN ZOE',
     +            'UOB WQM WPC WRC UFC UAN WOE',
     +            'VOB WQM WPC WRC VFC VAN WOE'  /
C
        REAL*8          hdr2 ( MXR8PM ),
     +                  evns2 ( MXR8PM, MXR8LV, MXR8VN, MXR8VT )
        REAL*8          r8sid, r8sid2, pob1, pob2
        CHARACTER*8     csid, csid2, subst2
        LOGICAL         match / .true. /
        EQUIVALENCE     ( r8sid, csid ), ( r8sid2, csid2 )
C
        SAVE            match, subst2, idate2
C-----------------------------------------------------------------------

cliu
 1000   continue
cliu

        iret = 0
C
C      If the previous call to this subroutine did not yield matching
C      mass and wind subsets, then READNS is already pointing at an
C      unmatched subset.  Otherwise, call READNS to advance the subset
C      pointer to the next subset.
C
        IF  ( match )  THEN
            CALL READNS  ( lunit, subset, idate, jret )
            IF  ( jret .ne. 0 )  THEN
                iret = -1
                RETURN
            END IF

cliu   select data type
       if(subset .ne. 'ADPUPA' .and. subset .ne. 'AIRCAR' .and.
     &    subset .ne. 'SATWND' .and. subset .ne. 'AIRCFT' ) go to 1000
cliu
        ELSE
            subset = subst2
            idate = idate2
        END IF
C*
C*      Read the HDR and EVNS data for the subset that is currently
C*      being pointed to.
C*
        CALL UFBINT  ( lunit, hdr, MXR8PM, 1, jret, head )
        DO ii = 1, MXR8VT
          CALL UFBEVN ( lunit, evns ( 1, 1, 1, ii ), MXR8PM, MXR8LV,
     +                     MXR8VN, nlev, ostr (ii) )
        END DO
C
C      Now, advance the subset pointer to the following subset and
C      read its HDR data.
C
cliu
 2000   continue
cliu
        CALL READNS  ( lunit, subst2, idate2, jret )
        IF  ( jret .ne. 0 )  THEN
            iret = 1
            RETURN
        END IF

cliu   select data type
       if(subst2.ne. 'ADPUPA' .and. subst2.ne. 'AIRCAR' .and. 
     &    subst2.ne. 'SATWND' .and. subst2.ne. 'AIRCFT' ) go to 2000
c        ! careful about 2000
cliu
        CALL UFBINT  ( lunit, hdr2, MXR8PM, 1, jret, head )
C 
C      Check whether these two subsets have identical SID, YOB, XOB,
C      ELV, and DHR values.  If so, then they are matching mass and
C      wind subsets for a single station report.
C
        match = .true.
C
        IF  ( subset .ne. subst2 )  THEN
            match = .false.
            RETURN
        END IF
C 
        r8sid = hdr (1)
        r8sid2 = hdr2 (1)
        IF  ( csid .ne. csid2 )  THEN
            match = .false.
            RETURN
        END IF
C 
        DO ii = 2, 5
            IF  ( hdr (ii) .ne. hdr2 (ii) )  THEN
                match = .false.
                RETURN
            END IF
        END DO
C
C      Read the EVNS data for the second of the two matching subsets.
C 
        DO ii = 1, MXR8VT
          CALL UFBEVN ( lunit, evns2 ( 1, 1, 1, ii ), MXR8PM, MXR8LV,
     +                     MXR8VN, nlev2, ostr (ii) )
        ENDDO
C
C      Combine the EVNS data for the two matching subsets into a
C      single 4-D array.  Do this by merging the EVNS2 array into the EVNS array.
C
        DO 10 lv2 = 1, nlev2
            DO lv = 1, nlev
                pob1 = evns ( 1, lv, 1, 1 )
                pob2 = evns2 ( 1, lv2, 1, 1 )
                IF  ( pob1 .eq. pob2 )  THEN
C
                  DO kk = 1, MXR8VT
                    DO jj = 1, MXR8VN
                      DO ii = 1, MXR8PM
                        IF  ( evns ( ii, lv, jj, kk ) .eq. R8BFMS ) THEN
                          evns ( ii, lv, jj, kk ) =
     +                          evns2 ( ii, lv2, jj, kk )
                        END IF
                      END DO
                    END DO
                  END DO
                  GO TO 10
                ELSE IF ( (pob2 .gt. pob1) .or. (lv .eq. nlev) ) THEN
C
                  nlev = nlev + 1
                  DO kk = 1, MXR8VT
                    DO jj = 1, MXR8VN
                      DO ii = 1, MXR8PM
                        evns ( ii, nlev, jj, kk ) =
     +                        evns2 ( ii, lv2, jj, kk )
                      END DO
                    END DO
                  END DO
                  GOTO 10
                END IF
            END DO
   10   END DO
C
        RETURN
        END
