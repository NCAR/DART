c    This code is not protected by the DART copyright agreement.
c    DART $Id$

c    Almost identical to the prepbufr.f program - adds 24 hours to obs
c    after midnight and before 3Z of the second day.  intended to be
c    run on 6Z input files only, where it outputs only obs at 3Z exactly.
c
c    Yet another revision of the BUFR preparation program.  It has options
c    to output specific humidity, relative humidity, or dewpoint obs, under
c    namelist control.  New code added by Ryan Torn.  14 Jan 2011
c    (Requires the matching update to the create_real_obs program.)
c
c    This is the 5th version of the BUFR preparation program for the
c    newly revised DART version. Ps and moisture obs are outputed as well.
c     
c    12/2009 added additional time window flavors, specified all of the variable
c    types.   (this is planned to be removed and observation windowing done
c    as a post process.  14jan2011)
c
c    A bug with the radiosonde T events selection j2 is fixed at 04/06/2005.
c    and surface data is included in this version.
c    For radiosonde T: T of Pc=1 and 6 is dry T; T of pc=8 is virtual T.
c    For Aircraft T: No moistuerobs and all T are dry T
c    For Surface obs, the T of pc=1 is already virtual T.
c
c    the READPB() routine currently has a select for obs types based on 
c    name; this should be removed and completely under namelist control.
c    the previously encountered problem was fixed and not related to this.
c    See the prepdecode/docs directory for the key to all the bufr codes.
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

      INTEGER, PARAMETER :: max_otype=300, max_qctype=15 
      REAL, PARAMETER    :: MISSING = -9999.
      INTEGER, PARAMETER :: I_MISSING = -9999
      CHARACTER  outstg*(200), subset*8, inf*19 , outf*20
      CHARACTER  var(MXR8VT)/'P','Q','T','Z','U','V'/
C
      character (len=9),  parameter :: infile  = 'prepqm.in'
      character (len=10), parameter :: outfile = 'prepqm.out'

      PARAMETER  ( NFILO = 15 )
      INTEGER  iunso ( NFILO )
     +          /   51,   52,   53,   54,   55,
     +              56,   57,   58,   59,   60,
     +              61,   62,   63,   64,   65  /
      CHARACTER*6     filo ( NFILO )
C     original list:
     +          / 'ADPUPA', 'AIRCAR', 'AIRCFT', 'SATWND', 'PROFLR',
     +            'VADWND', 'SATBOG', 'SATEMP', 'ADPSFC', 'SFCSHP',
     +            'SFCBOG', 'SPSSMI', 'SYNDAT', 'ERS1DA', 'GOESND'  /
C    updated list if you want to handle quikscat winds.
C     +          / 'ADPUPA', 'AIRCAR', 'AIRCFT', 'SATWND', 'PROFLR',
C     +            'VADWND', 'SATEMP', 'ADPSFC', 'SFCSHP', 'SFCBOG',
C     +            'SPSSMI', 'SYNDAT', 'ERS1DA', 'GOESND', 'QKSWND'/

      dimension tdata(8), udata(8), vdata(8), qdata(8), pdata(8)
      dimension zdata(8)
      integer :: wtype, ptype, qtype, ttype, ztype
c    The pc values are the 'program codes' that tell you what processing
c    was done on this observation.  As of now, these are unused, but could
c    be used for selection or diagnosis.
      integer :: pc_t, pc_q, pc_u, pc_v, pc_p, pc_z
      integer :: tqm, pqm, qqm, uqm, vqm, zqm
      logical :: found, uotype, uqcflag, use_this_data_real, 
     +           use_this_data_int, processed
      logical :: debug = .false.

c    Namelist Parameters
c----------------------------------------------------------------------

      integer :: qctype_use(max_otype)  ! qc values to accept (default all)
      real :: otype_use(max_otype),     ! report types to use (default all)
     +        obs_window       = 1.5,   ! observation time window (+/-hours)
     +        obs_window_upa   = 1.5,   ! sonde time window (+/-hours)
     +        obs_window_air   = 1.5,   ! aircraft obs time window (+/-hours)
     +        obs_window_sfc   = 0.8,   ! surface obs time window (+/-hours)
     +        obs_window_cw    = 1.5,   ! cloud wind obs time window (+/-hours)
     +        land_temp_error  = 2.5,   ! assumed err surface temp. obs (K)
     +        land_wind_error  = 3.5,   ! assumed err surface wind obs (m/s)
     +        land_moist_error = 0.2    ! assumed err surface moist. obs (%)

      namelist /prep_bufr_nml/ obs_window,
     +                         obs_window_upa,
     +                         obs_window_air,
     +                         obs_window_sfc,
     +                         obs_window_cw,
     +                         otype_use,
     +                         qctype_use,
     +                         land_temp_error,
     +                         land_wind_error,
     +                         land_moist_error

c    constants needed to compute saturation water vapor for qobs
c----------------------------------------------------------------------
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

c    read the prep_bufr_nml namelist from input.nml
c----------------------------------------------------------------------
       
      otype_use(:) = MISSING;  qctype_use(:) = I_MISSING
      inml_unit = 15
      open(inml_unit,file='input.nml', status='old',
     &        access='sequential', form='formatted')
      read(inml_unit, nml = prep_bufr_nml)
      close(inml_unit)

      inum_otype=0
      do i = 1, max_otype
        if ( otype_use(i) .eq. MISSING ) exit
        inum_otype = inum_otype + 1
      enddo

      inum_qctype=0
      do i = 1, max_qctype
        if ( qctype_use(i) .eq. MISSING ) exit
        inum_qctype = inum_qctype + 1
      enddo

c     overwrite observation windows if obs_window is defined
c----------------------------------------------------------------------

      if ( obs_window > 0.0 ) then
        print*, "Using the same observation window for all data: ",
     &     obs_window
        obs_window_upa = obs_window
        obs_window_air = obs_window
        obs_window_sfc = obs_window
        obs_window_cw  = obs_window
      endif

c    open the input bufr file and the output text file
c----------------------------------------------------------------------

      iprof = 0
      lunobs=210
      open(unit=lunobs, file = outfile, form='formatted')

      ibufr_unit = 11 
      OPEN  ( UNIT = ibufr_unit, FILE = infile, FORM = 'UNFORMATTED' )
      CALL OPENBF  ( ibufr_unit, 'IN', ibufr_unit )
      CALL DATELEN  ( 10 )

c    read the next station report from the input bufr file
c----------------------------------------------------------------------

10    CALL READPB  ( ibufr_unit, subset, idate, ierrpb )

      if ( debug ) print *, 'next obs type: ', subset(1:6)

      idate00 = idate/100
      hour01 = idate - idate00*100
      if( hour01 .eq. 0.0 ) hour01 = 24
      if ( debug ) print*, 'idate= ', idate, idate00, hour01

      IF ( ierrpb .eq. -1 )  STOP

      iprof = iprof + 1    !! station number

c    set the appropriate output file unit number
c----------------------------------------------------------------------

      ii = 1
      found = .false.
      DO WHILE  ( ( .not. found ) .and. ( ii .le. NFILO ) )
        IF( subset(1:6) .eq. filo ( ii ) )  THEN
          found = .true.
          iuno = iunso ( ii )

          time0 = hdr(4)
          stat_elv = hdr(5)

          if ( debug ) print*, 'idate= ', idate, idate00, hour01, time0
          hour01 =  hour01 + time0

          if ( USE_THIS_DATA_REAL(hdr(6),otype_use,inum_otype) .and. 
     &                    debug ) then  !  print data read in

            WRITE  ( UNIT = iuno, FMT = '("#", 132("-") )' )
            WRITE  ( UNIT = iuno, FMT = '("#", A3,6x,a6,1x,a6,3x,
     &             a5,5x,a4,3x,a5,4x,a4,3x,a4,1x,a3,1x,a3,2x,a6,
     &             6x,a3,6x,a3,4x,a5,1x,a8,2x,a7,6x,a3 )' )
     &             'SID', 'XOB', 'YOB','DHR','ELV','TYP','T29','ITP',
     &             'lev','var','OB','qm', 'pc', 'rc', 'fc', 'an', 'err'

            WRITE  ( UNIT = iuno, FMT = '("#", 132("-") )' )
          
          endif

        ELSE

          ii = ii + 1

        END IF

      END DO

c    check the observation time, skip if outside obs. window
c    the prepbufr.f version tests abs(time0), this one tests hour01
c----------------------------------------------------------------------
      IF ( ( .not. found ) .and. ( ierrpb .eq. 0 ) )  THEN
        if ( debug ) print*, 'record found w/no label match, val was: ',
     &            subset(1:6)
        GO TO 200
      ENDIF
      IF ( subset(1:6).eq.'SATWND' ) then
        IF ( hour01 .gt. obs_window_cw ) then
          if (debug) print*, 'satwnd outside time window, diff was: ',
     &              hour01
          GO TO 200
        ENDIF 
      ELSE IF (subset(1:6).eq.'ADPUPA'.or.
     &         subset(1:6).eq.'SATEMP' ) THEN
        IF ( hour01 .gt. obs_window_upa ) THEN 
         if ( debug ) print*, 
     &              'upper-air outside time window, diff was: ',
     &              abs(time0)
          GO TO 200
        ENDIF
      ELSE IF (subset(1:6).eq.'AIRCFT'.or.subset(1:6).eq.'AIRCAR') THEN
        IF ( hour01 .gt. obs_window_air ) THEN
          if ( debug ) print*,
     &              'aircraft outside time window, diff was: ',
     &              abs(time0)
          GO TO 200
	ENDIF
      ELSE IF (subset(1:6).eq.'ADPSFC'.or.subset(1:6).eq.'SFCSHP') THEN
        IF ( hour01 .gt. obs_window_sfc ) THEN
          if ( debug ) print*,
     &              'surface outside time window, diff was: ',
     &              hour01
          GO TO 200
	ENDIF
      ENDIF

c    place the location data in the appropriate array
c----------------------------------------------------------------------
 
      if( hdr(2) .ge. 360.0) hdr(2) = hdr(2) - 360.0
      if( hdr(2) .lt.   0.0) hdr(2) = hdr(2) + 360.0

      time = hdr(4)

      tdata(2) = hdr(2)
      tdata(3) = hdr(3)
      tdata(6) = 0
      tdata(7) = iprof + 1000000
      tdata(8) = hour01
         ttype = hdr(6)

      pdata(2) = hdr(2)
      pdata(3) = hdr(3)
      pdata(6) = 0
      pdata(7) = iprof + 3000000
      pdata(8) = hour01
         ptype = hdr(6)

      qdata(2) = hdr(2)
      qdata(3) = hdr(3)
      qdata(6) = 0
      qdata(7) = iprof + 5000000
      qdata(8) = hour01
         qtype = hdr(6)

      zdata(2) = hdr(2)
      zdata(3) = hdr(3)
      zdata(6) = 0
      zdata(7) = iprof + 7000000
      zdata(8) = hour01
         ztype = hdr(6)

      udata(2) = hdr(2)
      udata(3) = hdr(3)
      udata(7) = iprof + 2000000
      udata(8) = hour01           !! time

      vdata(2) = hdr(2)
      vdata(3) = hdr(3)
      vdata(7) = iprof + 9000000
      vdata(8) = hour01

c    u and v wind data are merged with T array for some obs. 
c    types.  Account for this and place in appropriate reports
c----------------------------------------------------------------------

      if ( hdr(6) .eq. 120.0 .or. hdr(6) .eq. 130.0 .or.
     &     hdr(6) .eq. 131.0 .or. hdr(6) .eq. 132.0 .or. 
     &     hdr(6) .eq. 133.0 ) then 
        wtype=hdr(6)+100
      else
        wtype=hdr(6)  
      endif   

c    check to see if this observation type is desired 
      if ( .NOT. USE_THIS_DATA_REAL(real(hdr(6)),otype_use,inum_otype) 
     &         ) then
        if (debug) print *, 'this obs type not in use-list, num was: ', 
     &       hdr(6)
        GO TO 200
      endif

      DO lv = 1, nlev  !  loop over all levels in the report

c    find out the event before virtural T and Q check step (pc=8.0)
c----------------------------------------------------------------------

        j2t = 1
        do ievent =1, MXR8VN
          if(evns(3,lv,ievent, 3) .lt. 8.0) then
            j2t = ievent
            exit 
          endif
        enddo

        j2q = 1
        do ievent =1, MXR8VN
          if(evns(3,lv,ievent, 2) .lt. 8.0) then
            j2q = ievent
            exit
          endif
        enddo

        DO kk = 1, MXR8VT

c     for T & Q, the events (pc>=8) are virtual T step. The events (pc=6, 1) are dry T
c     set qm of j2t/j2q event according to 1st event
c
          if((subset(1:6).eq.'ADPUPA' .or. subset(1:6).eq.'ADPSFC' .or.
     &        subset(1:6).eq.'SATEMP' .or. 
     &        subset(1:6).eq.'SFCSHP') .and.   var(kk).eq. 'T') then
            if(j2t.ne.1) evns(7,lv, j2t, kk) = evns(7,lv, 1, kk) 
            jj = j2t
          else if((subset(1:6).eq.'ADPUPA' .or. subset(1:6).eq.'ADPSFC'
     &      .or.   subset(1:6).eq.'SFCSHP') .and.  var(kk).eq.'Q') then
            if(j2q.ne.1) evns(7,lv, j2q, kk) = evns(7,lv, 1, kk) 
            jj = j2q
          else
            jj = 1                              
          endif
c
          WRITE( UNIT = outstg, FMT = '( A8, 1X,  
     &          2F7.2, 1X, F7.3, 1X, 2F8.1, 1X, F7.1, 1X, 
     &          F6.1, I4, 1X, A2, 7(1X,F8.1) )' )
     &     (hdr(ii), ii =1, 8),lv, var(kk), (evns(ii,lv,jj,kk), ii=1,7)

          DO mm = 1, 200 
            IF  ( outstg (mm:mm) .eq. '*' ) outstg (mm:mm) = ' '
          END DO

          if (debug) write(iuno, '(A)') trim(outstg)
        END DO

        processed = .false. 

c    set up the temperature observation data, use the j2t event or 1 if T
c----------------------------------------------------------------------
        if(subset(1:6).eq.'ADPUPA' .or. subset(1:6).eq.'ADPSFC' .or.
     &     subset(1:6).eq.'SATEMP' .or. 
     &     subset(1:6).eq.'SFCSHP' ) then

           toe = evns(7, lv, j2t, 3)
           tob = evns(1, lv, j2t, 3) + 273.16
           tqm = evns(2, lv, j2t, 3) 
          pc_t = evns(3, lv, j2t, 3)

        else  ! use the first event

           toe = evns(7, lv, 1, 3)
           tob = evns(1, lv, 1, 3) + 273.16
           tqm = evns(2, lv, 1, 3) 
          pc_t = evns(3, lv, 1, 3)

        endif
 
        if (pc_t < 0 .or. pc_t > 99) pc_t = 99
        if(tqm .eq. 0) toe = toe*0.9
        if(tqm .eq. 3) toe = toe*1.2

c    set up the moisture observation data, use the j2t event for RH and 
c     compute specific humidity in g/kg
c----------------------------------------------------------------------

        if(subset(1:6).eq.'ADPUPA' .or. subset(1:6).eq.'ADPSFC' .or.
     &     subset(1:6).eq.'SFCSHP' ) then

          qoe  = evns(7, lv, j2q, 2) * 0.1     ! g/kg   
          qob  = evns(1, lv, j2q, 2) * 1.0e-3  ! g/kg
          qqm  = evns(2, lv, j2q, 2) 
          pc_q = evns(3, lv, j2q, 2)

        else

          qoe  = evns(7, lv, 1, 2) * 0.1       ! g/kg    
          qob  = evns(1, lv, 1, 2) * 1.0e-3    ! g/kg
          qqm  = evns(2, lv, 1, 2) 
          pc_q = evns(3, lv, 1, 2)

        endif
        if (pc_q < 0 .or. pc_q > 99) pc_q = 99
        if(qqm .eq. 0) qoe = qoe*0.9
        if(qqm .eq. 3) qoe = qoe*1.2

c    set up the pressure observation data, units are mb
c----------------------------------------------------------------------

        poe  = evns(7, lv, 1, 1)     
        pob  = evns(1, lv, 1, 1)     
        pqm  = evns(2, lv, 1, 1) 
        pc_p = evns(3, lv, 1, 1)

        if (pc_p < 0 .or. pc_p > 99) pc_p = 99
        if(pqm .eq. 0) poe = poe*0.9
        if(pqm .eq. 3) poe = poe*1.2
        ppb = evns(1, lv, 1, 1)            ! mb
        !zob = evns(1, lv, 1, 4)            ! m

c    set up the height observation data, units are m
c----------------------------------------------------------------------

        zoe  = evns(7, lv, 1, 4)     
        zob = evns(1, lv, 1, 4)            ! m
        zqm  = evns(2, lv, 1, 4) 
        pc_z = evns(3, lv, 1, 4)

        if (pc_z < 0 .or. pc_z > 99) pc_z = 99
        if(zqm .eq. 0) zoe = zoe*0.9
        if(zqm .eq. 3) zoe = zoe*1.2

c    set up the wind observation data, units are m/s
c----------------------------------------------------------------------

        uoe  = evns(7, lv, 1, 5)
        uob  = evns(1, lv, 1, 5)
        uqm  = evns(2, lv, 1, 5) 
        pc_u = evns(3, lv, 1, 5)
        if(uqm .eq. 0) uoe = uoe*0.9
        if(uqm .eq. 3) uoe = uoe*1.2

        voe  = evns(7, lv, 1, 6)
        vob  = evns(1, lv, 1, 6) 
        vqm  = evns(2, lv, 1, 6)
        pc_v = evns(3, lv, 1, 6) 

        if (pc_u < 0 .or. pc_u > 99) pc_u = 99
        if (pc_v < 0 .or. pc_v > 99) pc_v = 99
        if(vqm .eq. 0) voe = voe*0.9
        if(vqm .eq. 3) voe = voe*1.2

c----------------------------------------------------------------------
c    this version of the code only selects observations at 3Z, and
c    adds 24 hours to the time.  this is the critical difference
c    between this file and prepbufr.f.  it is intended to be run only
c    on bufr files marked 06Z (which contain data from 3Z to 9Z) to
c    allow obs_seq files from 03:01Z to 03:00Z+1day to be created.

        if( abs(hour01 - 3.0) .gt. 0.0001)  go to 200 
        tdata(8) = hour01 + 24.0
        pdata(8) = hour01 + 24.0
        qdata(8) = hour01 + 24.0
        udata(8) = hour01 + 24.0          !! time
        vdata(8) = hour01 + 24.0


c    write out temperature observation from ADPUPA, AIRCAR, AIRCFT
c----------------------------------------------------------------------
        if (subset(1:6).eq.'ADPUPA' .or. subset(1:6).eq.'AIRCAR' .or.
     &     subset(1:6).eq.'SATEMP' .or. 
     &                                   subset(1:6).eq.'AIRCFT') then
          if (use_this_data_int(tqm,qctype_use,inum_qctype) .and. 
     &        use_this_data_int(pqm,qctype_use,inum_qctype) .and. 
     &        toe .lt. 1.e9 .and. tob .lt. 1.e9                  ) then

            tdata(1) = toe
            tdata(4) = ppb
            tdata(5) = tob

c    in some old files this appears to be out of range
c    and it seems to be unused in converting to an obs_seq so
c    i feel ok setting it to something that will fit in an I2 field.
            if (pc_t < 0 .or. pc_t > 99) pc_t = 99
            write(lunobs, 800) tdata, ttype, tqm, subset(1:6), pc_t
            processed = .true.

          else  
            if ( debug ) print *, 'skip temp, toe,tob,tqm,pqm = ',
     &           toe, tob, tqm, pqm
          endif   ! if use_this_data
        endif

c    write out temperature observation of SFCSHP, ADPSFC
c----------------------------------------------------------------------

        if (subset(1:6).eq.'SFCSHP' .or. subset(1:6).eq.'ADPSFC') then

          if (use_this_data_int(tqm,qctype_use,inum_qctype) .and.
     &        use_this_data_int(pqm,qctype_use,inum_qctype) .and. 
     &        tob .lt. 1.0e9 ) then
 
            if ( ttype .gt. 200 )  ttype = ttype - 100
            if ( ttype .eq. 184 )  ttype = 183
 
            tdata(1) = toe
            if ( toe .ge. 1.e9 ) tdata(1) = land_temp_error
            tdata(4) = zob
            tdata(5) = tob

            if (pc_t < 0 .or. pc_t > 99) pc_t = 99
            write(lunobs, 800) tdata, ttype, tqm, subset(1:6), pc_t
            processed = .true.

           else
             if ( debug ) print *, 'skip temp; tob,tqm,pqm = ',
     &        tob, tqm, pqm
           endif   ! if use_this_data
        endif

c    write out moisture observation from ADPUPA
c----------------------------------------------------------------------

        if ( subset(1:6) .eq. 'ADPUPA' ) then

          if (use_this_data_int(qqm,qctype_use,inum_qctype) .and.
     &        use_this_data_int(pqm,qctype_use,inum_qctype) .and.
     &        qoe .lt. 1.e9 .and. qob .lt. 1.e9             ) then

c           compute the error of specific moisture based on RH obs error

            qoe_rh = qoe

            es   = es0 * (tob / t0) ** fact1 * exp(fact2*(fact3-1./tob))
            qsat = eps * es / (pob - omeps * es)
            qoe  = max(0.1, qoe * qsat * 1000.0) ! to g/kg, set min value

            if ( debug ) print*, 'es= ', es, tob-273.16, qoe

            if( .not. use_this_data_int(tqm,qctype_use,inum_qctype))then
              if ( debug ) print*, 'upa bad = ', qoe  ! the T obs cannot be used for qoe
              qoe = 1.0e10
            endif

            if (qoe .lt. 9.9) then  ! skip large qoe obs.

              if (pc_q < 0 .or. pc_q > 99) pc_q = 99
              processed = .true.

              ! specific humidity
              qdata(1) = qoe
              qdata(4) = ppb
              qdata(5) = qob
              qdata(6) = 0.0
              write(lunobs, 800) qdata, qtype, qqm, subset(1:6), pc_q

              !  write RH data to file
              qdata(1) = qoe_rh
              qdata(5) = qob / (qsat * 1000.0)
              qdata(6) = 1.0
              write(lunobs, 800) qdata, qtype, qqm, subset(1:6), pc_q

              !  write dew point data to file
              eo = (qsat * (1.0-qoe_rh) * pob) / 
     &             (eps + omeps * qsat * (1.0-qoe_rh))
              qdata(1) = tob - 1.0 / (1.0/t0 - (rv/l0) * log(eo/es0))
              eo = (qob * pob * 0.001) / (eps + omeps * qob * 0.001)
              qdata(5) = 1.0 / (1.0/t0 - (rv/l0) * log(eo/es0))
              qdata(6) = 2.0
              write(lunobs, 800) qdata, qtype, qqm, subset(1:6), pc_q

            endif  ! qoe not too large

         else
            if ( debug ) print *, 'skip moist, qoe,qob,qqm,pqm = ',
     &          qoe, qpb, qqm, pqm
         endif   ! if use_this_data

        endif

c    write out moisture observation from SFCSHP ADPSFC
c----------------------------------------------------------------------

        if (subset(1:6).eq.'SFCSHP' .or. subset(1:6).eq.'ADPSFC') then

          if (use_this_data_int(qqm,qctype_use,inum_qctype) .and.
     &        use_this_data_int(pqm,qctype_use,inum_qctype) .and.
     &        qob .lt. 1.0e9                                ) then

           es   = es0 * (tob / t0) ** fact1 * exp(fact2*(fact3-1./tob))
           qsat = eps * es / (pob - omeps * es)
           if (qoe .gt. 1.e9) qoe = land_moist_error
           qoe_rh = qoe
           qoe  = max(0.1, qoe * qsat * 1000.0) ! to g/kg, set min value

           if( .not. use_this_data_int(tqm,qctype_use,inum_qctype)) then
             if ( debug ) print*, 'surface bad = ', qoe  
c             ! the T obs cannot be used for qoe
             qoe = 1.0e10
           endif

           if ( qtype .gt. 200 )  qtype = qtype - 100
           if ( qtype .eq. 184 )  qtype = 183

           if(qoe .lt. 9.9) then   ! skip large qoe obs

             if (pc_q < 0 .or. pc_q > 99) pc_q = 99
             processed = .true.

             ! specific humidity
             qdata(1) = qoe
             qdata(4) = zob
             qdata(5) = qob
             qdata(6) = 0.0
             write(lunobs, 800) qdata, qtype, qqm, subset(1:6), pc_q

             !  write RH data to file
             qdata(1) = qoe_rh
             qdata(5) = qob / (qsat * 1000.0)
             qdata(6) = 1.0
             write(lunobs, 800) qdata, qtype, qqm, subset(1:6), pc_q

             !  write dew point data to file
             eo = (qsat * (1.0-qoe_rh) * pob) / 
     &            (eps + omeps * qsat * (1.0-qoe_rh))
             qdata(1) = tob - 1.0 / (1.0/t0 - (rv/l0) * log(eo/es0))
             eo = (qob * pob * 0.001) / (eps + omeps * qob * 0.001)
             qdata(5) = 1.0 / (1.0/t0 - (rv/l0) * log(eo/es0))
             qdata(6) = 2.0
             write(lunobs, 800) qdata, qtype, qqm, subset(1:6), pc_q

           endif  ! if qoe not too large

          else
               if ( debug ) print *, 'skip moist, qob,qqm,pqm = ',
     &          qob, qqm, pqm
          endif   ! if use_this_data

        endif

c    write out surface pressure observations
c----------------------------------------------------------------------

        if (subset(1:6).eq.'ADPUPA' .or. 
     &      subset(1:6).eq.'SFCSHP' .or. subset(1:6).eq.'ADPSFC') then

          if(use_this_data_int(pqm,qctype_use,inum_qctype) .and. 
     &       poe.lt.1.e9 .and. pob.lt.1.e9 .and. stat_elv.eq.zob) then

            if ( ptype .gt. 200 )  ptype = ptype - 100
            if ( ptype .eq. 184 )  ptype = 183

            pdata(1) = poe
            pdata(4) = zob
            pdata(5) = pob

            if (pc_p < 0 .or. pc_p > 99) pc_p = 99
            write(lunobs, 800) pdata, ptype, pqm, subset(1:6), pc_p
            processed = .true.
    
           else
             if ( debug ) print *, 'skip press, poe,pob,pqm,z,s = ',
     &          poe, pob, pqm, zob, stat_elv
          endif   ! if use_this_data

        endif

c    write out wind observation of ADPUPA, AIRcraft, SATWND
c----------------------------------------------------------------------

        if (subset(1:6).eq.'ADPUPA' .or. subset(1:6).eq.'AIRCAR' .or.
     &      subset(1:6).eq.'AIRCFT' .or. subset(1:6).eq.'SATWND') then

          if(use_this_data_int(uqm,qctype_use,inum_qctype) .and. 
     &       use_this_data_int(pqm,qctype_use,inum_qctype) .and. 
     &       uoe .lt. 1.e9 .and. uob .lt. 1.e9             .and.
     &       voe .lt. 1.e9 .and. vob .lt. 1.e9             ) then

            udata(1) = uoe
            udata(4) = ppb
            udata(5) = uob
            udata(6) = vob
 
            vdata(1) = voe
            vdata(4) = ppb
            vdata(5) = vob
            vdata(6) = uob

            if (pc_u < 0 .or. pc_u > 99) pc_u = 99
            if (pc_v < 0 .or. pc_v > 99) pc_v = 99
            write(lunobs, 800) udata, wtype, uqm, subset(1:6), pc_u
            write(lunobs, 800) vdata, wtype, vqm, subset(1:6), pc_v
            processed = .true.

          else
             if ( debug ) print *, 'skip wind, uoe,uob,uqm,pqm = ',
     &          uoe, uob, uqm, pqm
          endif  ! if use_this_data
  
        endif

c    write out wind observation of SFCSHP ADPSFC
c----------------------------------------------------------------------

        if (subset(1:6).eq.'SFCSHP' .or. subset(1:6).eq.'ADPSFC' ) then

          if (use_this_data_int(uqm,qctype_use,inum_qctype) .and.
     &        use_this_data_int(pqm,qctype_use,inum_qctype) .and. 
     &        uob .lt. 1.e9                                 ) then

            udata(1) = uoe
            if ( uoe .ge. 1.e9 ) udata(1) = land_wind_error
            udata(4) = zob
            udata(5) = uob
            udata(6) = vob

            vdata(1) = voe
            if ( voe .ge. 1.e9 ) vdata(1) = land_wind_error
            vdata(4) = zob
            vdata(5) = vob
            vdata(6) = uob

            if (pc_u < 0 .or. pc_u > 99) pc_u = 99
            if (pc_v < 0 .or. pc_v > 99) pc_v = 99
            write(lunobs, 800) udata, wtype, uqm, subset(1:6), pc_u
            write(lunobs, 800) vdata, wtype, vqm, subset(1:6), pc_v
            processed = .true.

          else
             if ( debug ) print *, 'skip wind, uob,uqm,pqm = ',
     &          uob, uqm, pqm
          endif   ! if use_this_data

        endif

c    write out geopotential height observation of ADPUPA
c----------------------------------------------------------------------

        if (subset(1:6).eq.'ADPUPA') then
 
          if (use_this_data_int(zqm,qctype_use,inum_qctype) .and.
     &        use_this_data_int(pqm,qctype_use,inum_qctype) .and. 
     &        zoe. lt. 1.e9 .and. zob .lt. 1.e9            ) then

            zdata(1) = zoe
            zdata(4) = ppb
            zdata(5) = zob

            if (pc_z < 0 .or. pc_z > 99) pc_z = 99
            write(lunobs, 800) zdata, ztype, zqm, subset(1:6), pc_z
            processed = .true.
    

          else
             if ( debug ) print *, 
     &           'skip geopotential height, zoe,zob,zqm,pqm = ',
     &            zoe, zob, zqm, pqm
          endif   ! if use_this_data

        endif

c----------------------------------------------------------------------
       if (.not. processed) then
         if (debug)print*, 'bot of loop w/o processing, obs was: ',
     &              subset(1:6), ' lv =', lv
       endif

      enddo 
200    continue

      IF ( ierrpb .eq. 0 )  GO TO 10

800   format(f5.2,2f9.4,e12.5,f7.2,f7.2,f9.0,f7.3,i4,i2,1x,a6,i2)
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
        INCLUDE  'prepbufr.prm'
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

Cnsc
 1000   continue
Cnsc

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
     &    subset .ne. 'SATWND' .and. subset .ne. 'AIRCFT' .and. 
     &    subset .ne. 'SATEMP' .and. 
     &    subset .ne. 'ADPSFC' .and. subset .ne. 'SFCSHP') go to 1000
cliu

        ELSE
            subset = subst2
            idate = idate2
        END IF
C
C      Read the HDR and EVNS data for the subset that is currently
C      being pointed to.
C
        CALL UFBINT  ( lunit, hdr, MXR8PM, 1, jret, head )
        DO ii = 1, MXR8VT
          CALL UFBEVN ( lunit, evns ( 1, 1, 1, ii ), MXR8PM, MXR8LV,
     +                     MXR8VN, nlev, ostr (ii) )
        END DO
C
C      Now, advance the subset pointer to the following subset and
C      read its HDR data.
C
Cnsc
 2000   continue
Cnsc
        CALL READNS  ( lunit, subst2, idate2, jret )
        IF  ( jret .ne. 0 )  THEN
            iret = 1
            RETURN
        END IF
C
cliu   select data type
       if(subst2.ne. 'ADPUPA' .and. subst2.ne. 'AIRCAR' .and. 
     &    subst2.ne. 'SATWND' .and. subst2.ne. 'AIRCFT' .and.
     &    subst2.ne. 'SATEMP' .and. 
     &    subst2.ne. 'ADPSFC' .and. subst2.ne. 'SFCSHP') go to 2000
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

C-----------------------------------------------------------------------
        FUNCTION USE_THIS_DATA_INT(ifld, ufld, nfld)
C
C       ifld is the integer value to test, ufld is an array of
C       integers 'nfld' long.  if ifld is found in any of the
C       ufld entries, this function returns true, else false.
C       if nfld is 0, it bypasses the search and always
C       returns true; accepting any value.
C
        integer, intent(in) :: nfld, ifld, ufld(nfld)

        integer             :: nn
        logical             :: use_this_data_int

C       bypass test - always return true
        if (nfld == 0) then
          use_this_data_int = .TRUE.
          return
        endif
        
C       normal search code, look for match in list
        use_this_data_int = .FALSE.
        DO nn = 1, nfld
          IF ( ifld .eq. ufld(nn) ) THEN
            use_this_data_int = .TRUE.
            RETURN
          END IF
        END DO

        RETURN
        END

C-----------------------------------------------------------------------
        FUNCTION USE_THIS_DATA_REAL(rfld, ufld, nfld)
C
C       ifld is the real value to test, ufld is an array of
C       reals 'nfld' long.  if ifld is found in any of the
C       ufld entries, this function returns true, else false.
C       if nfld is 0, it bypasses the search and always
C       returns true; accepting any value.
C

        integer, intent(in) :: nfld
        real, intent(in)    :: rfld, ufld(nfld)

        integer             :: nn
        logical             :: USE_THIS_DATA_REAL

C       bypass test - always return true
        if (nfld == 0) then
          use_this_data_real = .TRUE.
          return
        endif
        
C       normal search code, look for match in list
        USE_THIS_DATA_REAL = .FALSE.
        DO nn = 1, nfld
          IF ( rfld .eq. ufld(nn) ) THEN
            USE_THIS_DATA_REAL = .TRUE.
            RETURN
          END IF
        END DO

        RETURN
        END

C <next few lines under version control, do not edit>
C $URL$
C $Id$
C $Revision$
C $Date$
