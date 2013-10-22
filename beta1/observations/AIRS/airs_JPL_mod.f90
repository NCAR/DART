! This code is not protected by the DART copyright agreement.
! DART $Id$

! adapted from original JPL code, example AIRS readers

module airs_JPL_mod

! the contents of this file are an amalgam of:
!  airs_ret_typ.inc
!  airs_ret_struct.inc
!  airs_ret_rdr.f
! although in several cases they were modified to use
! fortran 90 derived types and free format continuation
! lines, and to not use the unknown syntax 'double precision'.

implicit none
public

      integer   AIRS_RET_GEOXTRACK                      
      parameter(AIRS_RET_GEOXTRACK                      =    30)
      integer   AIRS_RET_GEOTRACK                       
      parameter(AIRS_RET_GEOTRACK                       =    45)
      integer   AIRS_RET_STDPRESSURELEV                 
      parameter(AIRS_RET_STDPRESSURELEV                 =    28)
      integer   AIRS_RET_STDPRESSURELAY                 
      parameter(AIRS_RET_STDPRESSURELAY                 =    28)
      integer   AIRS_RET_AIRSXTRACK                     
      parameter(AIRS_RET_AIRSXTRACK                     =     3)
      integer   AIRS_RET_AIRSTRACK                      
      parameter(AIRS_RET_AIRSTRACK                      =     3)
      integer   AIRS_RET_CLOUD                          
      parameter(AIRS_RET_CLOUD                          =     2)
      integer   AIRS_RET_CHANAMSUA                      
      parameter(AIRS_RET_CHANAMSUA                      =    15)
      integer   AIRS_RET_CHANHSB                        
      parameter(AIRS_RET_CHANHSB                        =     5)
      integer   AIRS_RET_MWHINGESURF                    
      parameter(AIRS_RET_MWHINGESURF                    =     7)
      integer   AIRS_RET_H2OFUNC                        
      parameter(AIRS_RET_H2OFUNC                        =    11)
      integer   AIRS_RET_O3FUNC                         
      parameter(AIRS_RET_O3FUNC                         =     9)
      integer   AIRS_RET_COFUNC                         
      parameter(AIRS_RET_COFUNC                         =     9)
      integer   AIRS_RET_CH4FUNC                        
      parameter(AIRS_RET_CH4FUNC                        =     7)
      integer   AIRS_RET_HINGESURF                      
      parameter(AIRS_RET_HINGESURF                      =   100)
      integer   AIRS_RET_H2OPRESSURELEV                 
      parameter(AIRS_RET_H2OPRESSURELEV                 =    15)
      integer   AIRS_RET_H2OPRESSURELAY                 
      parameter(AIRS_RET_H2OPRESSURELAY                 =    14)

! Record holds an entire granule of airs_ret
type airs_granule_type

! Attributes
        integer*2      NumSO2FOVs
        character*256  processing_level
        character*256  instrument
        character*256  DayNightFlag
        character*256  AutomaticQAFlag
        integer      NumTotalData
        integer      NumProcessData
        integer      NumSpecialData
        integer      NumBadData
        integer      NumMissingData
        integer      NumLandSurface
        integer      NumOceanSurface
        character*256  node_type
        integer      start_year
        integer      start_month
        integer      start_day
        integer      start_hour
        integer      start_minute
        real         start_sec
        integer      start_orbit
        integer      end_orbit
        integer      orbit_path
        integer      start_orbit_row
        integer      end_orbit_row
        integer      granule_number
        integer      num_scansets
        integer      num_scanlines
        real*8       start_Latitude
        real*8       start_Longitude
        real*8       start_Time
        real*8       end_Latitude
        real*8       end_Longitude
        real*8       end_Time
        real           eq_x_longitude
        real*8        eq_x_tai
        integer        orbitgeoqa
        integer*2      num_satgeoqa
        integer*2      num_glintgeoqa
        integer*2      num_moongeoqa
        integer*2      num_ftptgeoqa
        integer*2      num_zengeoqa
        integer*2      num_demgeoqa
        integer*2      num_fpe
        integer*2      LonGranuleCen
        integer*2      LatGranuleCen
        integer*2      LocTimeGranuleCen
        character*256  CO_first_guess
        character*256  CH4_first_guess

! Geolocation fields
        real*8 Latitude(AIRS_RET_GEOXTRACK, &
                                  AIRS_RET_GEOTRACK)
        real*8 Longitude(AIRS_RET_GEOXTRACK, & 
                                   AIRS_RET_GEOTRACK)
        real*8 Time(AIRS_RET_GEOXTRACK,  &
                              AIRS_RET_GEOTRACK)

! Data Fields
        integer*2      RetQAFlag( 30, 45)
        real           satheight( 45)
        real           satroll( 45)
        real           satpitch( 45)
        real           satyaw( 45)
        integer        satgeoqa( 45)
        integer*2      glintgeoqa( 45)
        integer*2      moongeoqa( 45)
        integer        ftptgeoqa( 30, 45)
        integer*2      zengeoqa( 30, 45)
        integer*2      demgeoqa( 30, 45)
        real*8         nadirTAI( 45)
        real*8         sat_lat( 45)
        real*8         sat_lon( 45)
        byte           scan_node_type( 45)
        real           satzen( 30, 45)
        real           satazi( 30, 45)
        real           solzen( 30, 45)
        real           solazi( 30, 45)
        real           glintlat( 45)
        real           glintlon( 45)
        integer*2      sun_glint_distance( 30, 45)
        real           topog( 30, 45)
        real           topog_err( 30, 45)
        real           landFrac( 30, 45)
        real           landFrac_err( 30, 45)
        real           pressStd( 28)
        real           pressH2O( 15)
        real           MWHingeSurfFreqGHz( 7)
        real           latAIRS( 3, 3, 30, 45)
        real           lonAIRS( 3, 3, 30, 45)
        integer*2      Qual_Guess_PSurf( 30, 45)
        real           PSurfStd( 30, 45)
        integer        nSurfStd( 30, 45)
        real           Press_mid_top_bndry( 30, 45)
        integer*2      nStd_mid_top_bndry( 30, 45)
        real           Press_bot_mid_bndry( 30, 45)
        integer*2      nStd_bot_mid_bndry( 30, 45)
        real           PBest( 30, 45)
        real           PGood( 30, 45)
        integer*2      nBestStd( 30, 45)
        integer*2      nGoodStd( 30, 45)
        integer*2      Qual_Temp_Profile_Top( 30, 45)
        integer*2      Qual_Temp_Profile_Mid( 30, 45)
        integer*2      Qual_Temp_Profile_Bot( 30, 45)
        real           TAirStd( 28, 30, 45)
        real           TAirStdErr( 28, 30, 45)
        real           TSurfAir( 30, 45)
        real           TSurfAirErr( 30, 45)
        integer*2      Qual_Surf( 30, 45)
        real           TSurfStd( 30, 45)
        real           TSurfStdErr( 30, 45)
        integer*2      numHingeSurf( 30, 45)
        real           freqEmis( 100, 30, 45)
        real           emisIRStd( 100, 30, 45)
        real           emisIRStdErr( 100, 30, 45)
        integer*2      Qual_MW_Only_Temp_Strat( 30, 45)
        integer*2      Qual_MW_Only_Temp_Tropo( 30, 45)
        real           TAirMWOnlyStd( 28, 30, 45)
        byte           MWSurfClass( 30, 45)
        real           sfcTbMWStd( 7, 30, 45)
        real           EmisMWStd( 7, 30, 45)
        real           EmisMWStdErr( 7, 30, 45)
        integer*2      Qual_MW_Only_H2O( 30, 45)
        real           totH2OMWOnlyStd( 30, 45)
        integer*2      Qual_H2O( 30, 45)
        real           H2OMMRStd( 14, 30, 45)
        real           H2OMMRStdErr( 14, 30, 45)
        real           totH2OStd( 30, 45)
        real           totH2OStdErr( 30, 45)
        real           H2OMMRSat( 14, 30, 45)
        real           H2OMMRSat_liquid( 14, 30, 45)
        integer*2      num_H2O_Func( 30, 45)
        real           H2O_verticality( 11, 30, 45)
        integer*2      Qual_O3( 30, 45)
        real           totO3Std( 30, 45)
        real           totO3StdErr( 30, 45)
        real           O3VMRStd( 28, 30, 45)
        real           O3VMRStdErr( 28, 30, 45)
        integer*2      num_O3_Func( 30, 45)
        real           O3_verticality( 9, 30, 45)
        integer*2      Qual_CO( 30, 45)
        real           CO_total_column( 30, 45)
        integer*2      num_CO_Func( 30, 45)
        integer        CO_trapezoid_layers( 9)
        real           CO_eff_press( 9, 30, 45)
        real           CO_VMR_eff( 9, 30, 45)
        real           CO_VMR_eff_err( 9, 30, 45)
        real           CO_verticality( 9, 30, 45)
        real           CO_dof( 30, 45)
        integer*2      Qual_CH4( 30, 45)
        real           CH4_total_column( 30, 45)
        integer*2      num_CH4_Func( 30, 45)
        integer        CH4_trapezoid_layers( 7)
        real           CH4_eff_press( 7, 30, 45)
        real           CH4_VMR_eff( 7, 30, 45)
        real           CH4_VMR_eff_err( 7, 30, 45)
        real           CH4_verticality( 7, 30, 45)
        real           CH4_dof( 30, 45)
        real           PTropopause( 30, 45)
        real           T_Tropopause( 30, 45)
        real           GP_Tropopause( 30, 45)
        real           GP_Height( 28, 30, 45)
        real           GP_Height_MWOnly( 28, 30, 45)
        real           GP_Surface( 30, 45)
        integer*2      Qual_Cloud_OLR( 30, 45)
        integer        numCloud( 30, 45)
        real           TCldTopStd( 2, 30, 45)
        real           TCldTopStdErr( 2, 30, 45)
        real           PCldTopStd( 2, 30, 45)
        real           PCldTopStdErr( 2, 30, 45)
        real           CldFrcStd( 2, 3, 3, 30, 45)
        real           CldFrcStdErr( 2, 3, 3, 30, 45)
        real           olr( 30, 45)
        real           olr_err( 30, 45)
        integer*2      Qual_clrolr( 30, 45)
        real           clrolr( 30, 45)
        real           clrolr_err( 30, 45)
        integer*2      dust_flag( 3, 3, 30, 45)
        integer*2      spectral_clear_indicator( 3, 3, 30, 45)
        integer*2      num_clear_spectral_indicator( 30, 45)
        real           CC_noise_eff_amp_factor( 30, 45)
        real           CC1_noise_eff_amp_factor( 30, 45)
        real           totCldH2OStd( 30, 45)
        real           totCldH2OStdErr( 30, 45)
        real           CC1_Resid( 30, 45)
        real           CCfinal_Resid( 30, 45)
        real           CCfinal_Noise_Amp( 30, 45)
        real           Tdiff_IR_MW_ret( 30, 45)
        real           Tdiff_IR_4CC1( 30, 45)
        real           TSurfdiff_IR_4CC1( 30, 45)
        real           TSurfdiff_IR_4CC2( 30, 45)
        real           AMSU_Chans_Resid( 30, 45)
        real           TotCld_4_CCfinal( 30, 45)
        real           Surf_Resid_Ratio( 30, 45)
        real           Temp_Resid_Ratio( 30, 45)
        real           Water_Resid_Ratio( 30, 45)
        real           Cloud_Resid_Ratio( 30, 45)
        real           O3_Resid_Ratio( 30, 45)
        real           CO_Resid_Ratio( 30, 45)
        real           CH4_Resid_Ratio( 30, 45)
        real           MWCheck_Resid_Ratio( 30, 45)
        real           O3_dof( 30, 45)
        byte           all_spots_avg( 30, 45)
        byte           MW_ret_used( 30, 45)
        real           Initial_CC_score( 30, 45)
        byte           retrieval_type( 30, 45)
        byte           Startup( 30, 45)
END type airs_granule_type


contains

! This function is autogenerated by the mkezio program to read
! an AIRS swath of type "L2_Standard_atmospheric&surface_product" from file given by the
! file_name argument into a buffer pointed to by the airs_ret_gran
! argument.  The caller owns the buffer.  The entire granule
! is read -- every attribute and field, the whole lat/lon/time
! extent.
!
! Errors opening the file, etc. are fatal and cause STOP.
! Problems reading individual attributes or fields are reported to
! the console but do not interrupt program flow.

      subroutine airs_ret_rdr(file_name, airs_ret_gran)
      character(len=*), intent(in) :: file_name(:)
      type(airs_granule_type) :: airs_ret_gran

      integer statn           ! HDF-EOS status. 0 for success
      integer fid             ! HDF-EOS file ID
      integer swid            ! HDF-EOS swath ID
      integer nchar           ! Number of characters
      character*256 swathname   ! Name of swath
      integer nswath          ! Number of swaths
      integer start(10) /0,0,0,0,0, 0,0,0,0,0/
                                ! start of each dimensions for Swath I/O
                                ! 0 => start with first element
      integer stride(10)/1,1,1,1,1, 1,1,1,1,1/
                                ! stride of each dimensions for Swath I/O
                                ! 1 => use every element
      integer edge(10)        ! size of each dimension for swath I/O
                                ! will be set for each individual read
      integer swopen, swinqswath, swattach
      integer swrdfld, swrdattr
      integer swdetach, swclose

      fid = swopen(file_name, 1)
      if (fid .eq. -1) then
        print *, "Error ", fid, " opening file ", file_name
        stop
      end if

      ! Get name of swath(s)
      nswath = swinqswath(file_name, swathname, nchar)
      if (nswath .ne. 1) then
        print *, "swinqswath found ", nswath, " swaths for file ", &
                 file_name, " Need exactly 1"
        stop
      end if

      ! There's exactly one swath.  Make sure it is the right one.
      if (swathname .ne. &
         'L2_Standard_atmospheric&surface_product') then
        print *, "Error: bad swath name ", swathname, " in file ", &
                 file_name
        print *, "Expected L2_Standard_atmospheric&surface_product"
        stop
      end if

      ! Attach to (open) the one swath.
      swid = swattach(fid, swathname)
      if (swid .eq. -1) then
        print *, "Failed to attach to swath ", swathname, &
                 " in file ", file_name
        stop
      end if

! Attributes
      statn = swrdattr(swid, "NumSO2FOVs", &
                   airs_ret_gran%NumSO2FOVs)
      if (statn .ne. 0) &
        print *, "Error ", statn, " reading attribute ", &
                  "NumSO2FOVs"

      statn = swrdattr(swid, "processing_level", &
                         airs_ret_gran%processing_level)
      if (statn .ne. 0) &
        print *, "Error ", statn, " reading attribute ", &
                  "processing_level"

      statn = swrdattr(swid, "instrument", &
                         airs_ret_gran%instrument)
      if (statn .ne. 0) &
        print *, "Error ", statn, " reading attribute ", &
                  "instrument"

      statn = swrdattr(swid, "DayNightFlag", &
                         airs_ret_gran%DayNightFlag)
      if (statn .ne. 0) &
        print *, "Error ", statn, " reading attribute ", &
                  "DayNightFlag"

      statn = swrdattr(swid, "AutomaticQAFlag", &
                         airs_ret_gran%AutomaticQAFlag)
      if (statn .ne. 0) &
        print *, "Error ", statn, " reading attribute ", &
                  "AutomaticQAFlag"

      statn = swrdattr(swid, "NumTotalData", &
                   airs_ret_gran%NumTotalData)
      if (statn .ne. 0) &
        print *, "Error ", statn, " reading attribute ", &
                  "NumTotalData"

      statn = swrdattr(swid, "NumProcessData", &
                   airs_ret_gran%NumProcessData)
      if (statn .ne. 0) &
        print *, "Error ", statn, " reading attribute ", &
                  "NumProcessData"

      statn = swrdattr(swid, "NumSpecialData", &
                   airs_ret_gran%NumSpecialData)
      if (statn .ne. 0) &
        print *, "Error ", statn, " reading attribute ", &
                  "NumSpecialData"

      statn = swrdattr(swid, "NumBadData", &
                   airs_ret_gran%NumBadData)
      if (statn .ne. 0) &
        print *, "Error ", statn, " reading attribute ", &
                  "NumBadData"

      statn = swrdattr(swid, "NumMissingData", &
                   airs_ret_gran%NumMissingData)
      if (statn .ne. 0) &
        print *, "Error ", statn, " reading attribute ", &
                  "NumMissingData"

      statn = swrdattr(swid, "NumLandSurface", &
                   airs_ret_gran%NumLandSurface)
      if (statn .ne. 0) &
        print *, "Error ", statn, " reading attribute ", &
                  "NumLandSurface"

      statn = swrdattr(swid, "NumOceanSurface", &
                   airs_ret_gran%NumOceanSurface)
      if (statn .ne. 0) &
        print *, "Error ", statn, " reading attribute ", &
                  "NumOceanSurface"

      statn = swrdattr(swid, "node_type", &
                         airs_ret_gran%node_type)
      if (statn .ne. 0) &
        print *, "Error ", statn, " reading attribute ", &
                  "node_type"

      statn = swrdattr(swid, "start_year", &
                   airs_ret_gran%start_year)
      if (statn .ne. 0) &
        print *, "Error ", statn, " reading attribute ", &
                  "start_year"

      statn = swrdattr(swid, "start_month", &
                   airs_ret_gran%start_month)
      if (statn .ne. 0) &
        print *, "Error ", statn, " reading attribute ", &
                  "start_month"

      statn = swrdattr(swid, "start_day", &
                   airs_ret_gran%start_day)
      if (statn .ne. 0) &
        print *, "Error ", statn, " reading attribute ", &
                  "start_day"

      statn = swrdattr(swid, "start_hour", &
                   airs_ret_gran%start_hour)
      if (statn .ne. 0) &
        print *, "Error ", statn, " reading attribute ", &
                  "start_hour"

      statn = swrdattr(swid, "start_minute", &
                   airs_ret_gran%start_minute)
      if (statn .ne. 0) &
        print *, "Error ", statn, " reading attribute ", &
                  "start_minute"

      statn = swrdattr(swid, "start_sec", &
                   airs_ret_gran%start_sec)
      if (statn .ne. 0) &
        print *, "Error ", statn, " reading attribute ", &
                  "start_sec"

      statn = swrdattr(swid, "start_orbit", &
                   airs_ret_gran%start_orbit)
      if (statn .ne. 0) &
        print *, "Error ", statn, " reading attribute ", &
                  "start_orbit"

      statn = swrdattr(swid, "end_orbit", &
                   airs_ret_gran%end_orbit)
      if (statn .ne. 0) &
        print *, "Error ", statn, " reading attribute ", &
                  "end_orbit"

      statn = swrdattr(swid, "orbit_path", &
                   airs_ret_gran%orbit_path)
      if (statn .ne. 0) &
        print *, "Error ", statn, " reading attribute ", &
                  "orbit_path"

      statn = swrdattr(swid, "start_orbit_row", &
                   airs_ret_gran%start_orbit_row)
      if (statn .ne. 0) &
        print *, "Error ", statn, " reading attribute ", &
                  "start_orbit_row"

      statn = swrdattr(swid, "end_orbit_row", &
                   airs_ret_gran%end_orbit_row)
      if (statn .ne. 0) &
        print *, "Error ", statn, " reading attribute ", &
                  "end_orbit_row"

      statn = swrdattr(swid, "granule_number", &
                   airs_ret_gran%granule_number)
      if (statn .ne. 0) &
        print *, "Error ", statn, " reading attribute ", &
                  "granule_number"

      statn = swrdattr(swid, "num_scansets", &
                   airs_ret_gran%num_scansets)
      if (statn .ne. 0) &
        print *, "Error ", statn, " reading attribute ", &
                  "num_scansets"

      statn = swrdattr(swid, "num_scanlines", &
                   airs_ret_gran%num_scanlines)
      if (statn .ne. 0) &
        print *, "Error ", statn, " reading attribute ", &
                  "num_scanlines"

      statn = swrdattr(swid, "start_Latitude", &
                   airs_ret_gran%start_Latitude)
      if (statn .ne. 0) &
        print *, "Error ", statn, " reading attribute ", &
                  "start_Latitude"

      statn = swrdattr(swid, "start_Longitude", &
                   airs_ret_gran%start_Longitude)
      if (statn .ne. 0) &
        print *, "Error ", statn, " reading attribute ", &
                  "start_Longitude"

      statn = swrdattr(swid, "start_Time", &
                   airs_ret_gran%start_Time)
      if (statn .ne. 0) &
        print *, "Error ", statn, " reading attribute ", &
                  "start_Time"

      statn = swrdattr(swid, "end_Latitude", &
                   airs_ret_gran%end_Latitude)
      if (statn .ne. 0) &
        print *, "Error ", statn, " reading attribute ", &
                  "end_Latitude"

      statn = swrdattr(swid, "end_Longitude", &
                   airs_ret_gran%end_Longitude)
      if (statn .ne. 0) &
        print *, "Error ", statn, " reading attribute ", &
                  "end_Longitude"

      statn = swrdattr(swid, "end_Time", &
                   airs_ret_gran%end_Time)
      if (statn .ne. 0) &
        print *, "Error ", statn, " reading attribute ", &
                  "end_Time"

      statn = swrdattr(swid, "eq_x_longitude", &
                   airs_ret_gran%eq_x_longitude)
      if (statn .ne. 0) &
        print *, "Error ", statn, " reading attribute ", &
                  "eq_x_longitude"

      statn = swrdattr(swid, "eq_x_tai", &
                   airs_ret_gran%eq_x_tai)
      if (statn .ne. 0) &
        print *, "Error ", statn, " reading attribute ", &
                  "eq_x_tai"

      statn = swrdattr(swid, "orbitgeoqa", &
                   airs_ret_gran%orbitgeoqa)
      if (statn .ne. 0) &
        print *, "Error ", statn, " reading attribute ", &
                  "orbitgeoqa"

      statn = swrdattr(swid, "num_satgeoqa", &
                   airs_ret_gran%num_satgeoqa)
      if (statn .ne. 0) &
        print *, "Error ", statn, " reading attribute ", &
                  "num_satgeoqa"

      statn = swrdattr(swid, "num_glintgeoqa", &
                   airs_ret_gran%num_glintgeoqa)
      if (statn .ne. 0) &
        print *, "Error ", statn, " reading attribute ", &
                  "num_glintgeoqa"

      statn = swrdattr(swid, "num_moongeoqa", &
                   airs_ret_gran%num_moongeoqa)
      if (statn .ne. 0) &
        print *, "Error ", statn, " reading attribute ", &
                  "num_moongeoqa"

      statn = swrdattr(swid, "num_ftptgeoqa", &
                   airs_ret_gran%num_ftptgeoqa)
      if (statn .ne. 0) &
        print *, "Error ", statn, " reading attribute ", &
                  "num_ftptgeoqa"

      statn = swrdattr(swid, "num_zengeoqa", &
                   airs_ret_gran%num_zengeoqa)
      if (statn .ne. 0) &
        print *, "Error ", statn, " reading attribute ", &
                  "num_zengeoqa"

      statn = swrdattr(swid, "num_demgeoqa", &
                   airs_ret_gran%num_demgeoqa)
      if (statn .ne. 0) &
        print *, "Error ", statn, " reading attribute ", &
                  "num_demgeoqa"

      statn = swrdattr(swid, "num_fpe", &
                   airs_ret_gran%num_fpe)
      if (statn .ne. 0) &
        print *, "Error ", statn, " reading attribute ", &
                  "num_fpe"

      statn = swrdattr(swid, "LonGranuleCen", &
                   airs_ret_gran%LonGranuleCen)
      if (statn .ne. 0) &
        print *, "Error ", statn, " reading attribute ", &
                  "LonGranuleCen"

      statn = swrdattr(swid, "LatGranuleCen", &
                   airs_ret_gran%LatGranuleCen)
      if (statn .ne. 0) &
        print *, "Error ", statn, " reading attribute ", &
                  "LatGranuleCen"

      statn = swrdattr(swid, "LocTimeGranuleCen", &
                   airs_ret_gran%LocTimeGranuleCen)
      if (statn .ne. 0) &
        print *, "Error ", statn, " reading attribute ", &
                  "LocTimeGranuleCen"

      statn = swrdattr(swid, "CO_first_guess", &
                         airs_ret_gran%CO_first_guess)
      if (statn .ne. 0) &
        print *, "Error ", statn, " reading attribute ", &
                  "CO_first_guess"

      statn = swrdattr(swid, "CH4_first_guess", &
                         airs_ret_gran%CH4_first_guess)
      if (statn .ne. 0) &
        print *, "Error ", statn, " reading attribute ", &
                  "CH4_first_guess"


! Geolocation fields
      edge(1) = AIRS_RET_GEOXTRACK
      edge(2) = AIRS_RET_GEOTRACK
      statn = swrdfld(swid, "Latitude", start, stride, edge, &
                      airs_ret_gran%Latitude)
      if (statn .ne. 0)  &
        print *, "Error ", statn, " reading field Latitude"

      statn = swrdfld(swid, "Longitude", start, stride, edge, &
                      airs_ret_gran%Longitude)
      if (statn .ne. 0) &
        print *, "Error ", statn, " reading field Longitude"

      statn = swrdfld(swid, "Time", start, stride, edge, &
                      airs_ret_gran%Time)
      if (statn .ne. 0) &
        print *, "Error ", statn, " reading field Time"


! Data Fields
      edge(2) = 45
      edge(1) = 30
      statn = SWrdfld(swid, "RetQAFlag", &
                   start, stride, edge, &
                   airs_ret_gran%RetQAFlag)
      if (statn .ne. 0) &
        print *, "Error ", statn, " reading field ", &
                  "RetQAFlag"

      edge(1) = 45
      statn = SWrdfld(swid, "satheight", &
                   start, stride, edge, &
                   airs_ret_gran%satheight)
      if (statn .ne. 0) &
        print *, "Error ", statn, " reading field ", &
                  "satheight"

      edge(1) = 45
      statn = SWrdfld(swid, "satroll", &
                   start, stride, edge, &
                   airs_ret_gran%satroll)
      if (statn .ne. 0) &
        print *, "Error ", statn, " reading field ", &
                  "satroll"

      edge(1) = 45
      statn = SWrdfld(swid, "satpitch", &
                   start, stride, edge, &
                   airs_ret_gran%satpitch)
      if (statn .ne. 0) &
        print *, "Error ", statn, " reading field ", &
                  "satpitch"

      edge(1) = 45
      statn = SWrdfld(swid, "satyaw", &
                   start, stride, edge, &
                   airs_ret_gran%satyaw)
      if (statn .ne. 0) &
        print *, "Error ", statn, " reading field ", &
                  "satyaw"

      edge(1) = 45
      statn = SWrdfld(swid, "satgeoqa", &
                   start, stride, edge, &
                   airs_ret_gran%satgeoqa)
      if (statn .ne. 0) &
        print *, "Error ", statn, " reading field ", &
                  "satgeoqa"

      edge(1) = 45
      statn = SWrdfld(swid, "glintgeoqa", &
                   start, stride, edge, &
                   airs_ret_gran%glintgeoqa)
      if (statn .ne. 0) &
        print *, "Error ", statn, " reading field ", &
                  "glintgeoqa"

      edge(1) = 45
      statn = SWrdfld(swid, "moongeoqa", &
                   start, stride, edge, &
                   airs_ret_gran%moongeoqa)
      if (statn .ne. 0) &
        print *, "Error ", statn, " reading field ", &
                  "moongeoqa"

      edge(2) = 45
      edge(1) = 30
      statn = SWrdfld(swid, "ftptgeoqa", &
                   start, stride, edge, &
                   airs_ret_gran%ftptgeoqa)
      if (statn .ne. 0) &
        print *, "Error ", statn, " reading field ", &
                  "ftptgeoqa"

      edge(2) = 45
      edge(1) = 30
      statn = SWrdfld(swid, "zengeoqa", &
                   start, stride, edge, &
                   airs_ret_gran%zengeoqa)
      if (statn .ne. 0) &
        print *, "Error ", statn, " reading field ", &
                  "zengeoqa"

      edge(2) = 45
      edge(1) = 30
      statn = SWrdfld(swid, "demgeoqa", &
                   start, stride, edge, &
                   airs_ret_gran%demgeoqa)
      if (statn .ne. 0) &
        print *, "Error ", statn, " reading field ", &
                  "demgeoqa"

      edge(1) = 45
      statn = SWrdfld(swid, "nadirTAI", &
                   start, stride, edge, &
                   airs_ret_gran%nadirTAI)
      if (statn .ne. 0) &
        print *, "Error ", statn, " reading field ", &
                  "nadirTAI"

      edge(1) = 45
      statn = SWrdfld(swid, "sat_lat", &
                   start, stride, edge, &
                   airs_ret_gran%sat_lat)
      if (statn .ne. 0) &
        print *, "Error ", statn, " reading field ", &
                  "sat_lat"

      edge(1) = 45
      statn = SWrdfld(swid, "sat_lon", &
                   start, stride, edge, &
                   airs_ret_gran%sat_lon)
      if (statn .ne. 0) &
        print *, "Error ", statn, " reading field ", &
                  "sat_lon"

      edge(1) = 45
      statn = SWrdfld(swid, "scan_node_type", &
                   start, stride, edge, &
                   airs_ret_gran%scan_node_type)
      if (statn .ne. 0) &
        print *, "Error ", statn, " reading field ", &
                  "scan_node_type"

      edge(2) = 45
      edge(1) = 30
      statn = SWrdfld(swid, "satzen", &
                   start, stride, edge, &
                   airs_ret_gran%satzen)
      if (statn .ne. 0) &
        print *, "Error ", statn, " reading field ", &
                  "satzen"

      edge(2) = 45
      edge(1) = 30
      statn = SWrdfld(swid, "satazi", &
                   start, stride, edge, &
                   airs_ret_gran%satazi)
      if (statn .ne. 0) &
        print *, "Error ", statn, " reading field ", &
                  "satazi"

      edge(2) = 45
      edge(1) = 30
      statn = SWrdfld(swid, "solzen", &
                   start, stride, edge, &
                   airs_ret_gran%solzen)
      if (statn .ne. 0) &
        print *, "Error ", statn, " reading field ", &
                  "solzen"

      edge(2) = 45
      edge(1) = 30
      statn = SWrdfld(swid, "solazi", &
                   start, stride, edge, &
                   airs_ret_gran%solazi)
      if (statn .ne. 0) &
        print *, "Error ", statn, " reading field ", &
                  "solazi"

      edge(1) = 45
      statn = SWrdfld(swid, "glintlat", &
                   start, stride, edge, &
                   airs_ret_gran%glintlat)
      if (statn .ne. 0) &
        print *, "Error ", statn, " reading field ", &
                  "glintlat"

      edge(1) = 45
      statn = SWrdfld(swid, "glintlon", &
                   start, stride, edge, &
                   airs_ret_gran%glintlon)
      if (statn .ne. 0) &
        print *, "Error ", statn, " reading field ", &
                  "glintlon"

      edge(2) = 45
      edge(1) = 30
      statn = SWrdfld(swid, "sun_glint_distance", &
                   start, stride, edge, &
                   airs_ret_gran%sun_glint_distance)
      if (statn .ne. 0) &
        print *, "Error ", statn, " reading field ", &
                  "sun_glint_distance"

      edge(2) = 45
      edge(1) = 30
      statn = SWrdfld(swid, "topog", &
                   start, stride, edge, &
                   airs_ret_gran%topog)
      if (statn .ne. 0) &
        print *, "Error ", statn, " reading field ", &
                  "topog"

      edge(2) = 45
      edge(1) = 30
      statn = SWrdfld(swid, "topog_err", &
                   start, stride, edge, &
                   airs_ret_gran%topog_err)
      if (statn .ne. 0) &
        print *, "Error ", statn, " reading field ", &
                  "topog_err"

      edge(2) = 45
      edge(1) = 30
      statn = SWrdfld(swid, "landFrac", &
                   start, stride, edge, &
                   airs_ret_gran%landFrac)
      if (statn .ne. 0) &
        print *, "Error ", statn, " reading field ", &
                  "landFrac"

      edge(2) = 45
      edge(1) = 30
      statn = SWrdfld(swid, "landFrac_err", &
                   start, stride, edge, &
                   airs_ret_gran%landFrac_err)
      if (statn .ne. 0) &
        print *, "Error ", statn, " reading field ", &
                  "landFrac_err"

      edge(1) = 28
      statn = SWrdfld(swid, "pressStd", &
                   start, stride, edge, &
                   airs_ret_gran%pressStd)
      if (statn .ne. 0) &
        print *, "Error ", statn, " reading field ", &
                  "pressStd"

      edge(1) = 15
      statn = SWrdfld(swid, "pressH2O", &
                   start, stride, edge, &
                   airs_ret_gran%pressH2O)
      if (statn .ne. 0) &
        print *, "Error ", statn, " reading field ", &
                  "pressH2O"

      edge(1) = 7
      statn = SWrdfld(swid, "MWHingeSurfFreqGHz", &
                   start, stride, edge, &
                   airs_ret_gran%MWHingeSurfFreqGHz)
      if (statn .ne. 0) &
        print *, "Error ", statn, " reading field ", &
                  "MWHingeSurfFreqGHz"

      edge(4) = 45
      edge(3) = 30
      edge(2) = 3
      edge(1) = 3
      statn = SWrdfld(swid, "latAIRS", &
                   start, stride, edge, &
                   airs_ret_gran%latAIRS)
      if (statn .ne. 0) &
        print *, "Error ", statn, " reading field ", &
                  "latAIRS"

      edge(4) = 45
      edge(3) = 30
      edge(2) = 3
      edge(1) = 3
      statn = SWrdfld(swid, "lonAIRS", &
                   start, stride, edge, &
                   airs_ret_gran%lonAIRS)
      if (statn .ne. 0) &
        print *, "Error ", statn, " reading field ", &
                  "lonAIRS"

      edge(2) = 45
      edge(1) = 30
      statn = SWrdfld(swid, "Qual_Guess_PSurf", &
                   start, stride, edge, &
                   airs_ret_gran%Qual_Guess_PSurf)
      if (statn .ne. 0) &
        print *, "Error ", statn, " reading field ", &
                  "Qual_Guess_PSurf"

      edge(2) = 45
      edge(1) = 30
      statn = SWrdfld(swid, "PSurfStd", &
                   start, stride, edge, &
                   airs_ret_gran%PSurfStd)
      if (statn .ne. 0) &
        print *, "Error ", statn, " reading field ", &
                  "PSurfStd"

      edge(2) = 45
      edge(1) = 30
      statn = SWrdfld(swid, "nSurfStd", &
                   start, stride, edge, &
                   airs_ret_gran%nSurfStd)
      if (statn .ne. 0) &
        print *, "Error ", statn, " reading field ", &
                  "nSurfStd"

      edge(2) = 45
      edge(1) = 30
      statn = SWrdfld(swid, "Press_mid_top_bndry", &
                   start, stride, edge, &
                   airs_ret_gran%Press_mid_top_bndry)
      if (statn .ne. 0) &
        print *, "Error ", statn, " reading field ", &
                  "Press_mid_top_bndry"

      edge(2) = 45
      edge(1) = 30
      statn = SWrdfld(swid, "nStd_mid_top_bndry", &
                   start, stride, edge, &
                   airs_ret_gran%nStd_mid_top_bndry)
      if (statn .ne. 0) &
        print *, "Error ", statn, " reading field ", &
                  "nStd_mid_top_bndry"

      edge(2) = 45
      edge(1) = 30
      statn = SWrdfld(swid, "Press_bot_mid_bndry", &
                   start, stride, edge, &
                   airs_ret_gran%Press_bot_mid_bndry)
      if (statn .ne. 0) &
        print *, "Error ", statn, " reading field ", &
                  "Press_bot_mid_bndry"

      edge(2) = 45
      edge(1) = 30
      statn = SWrdfld(swid, "nStd_bot_mid_bndry", &
                   start, stride, edge, &
                   airs_ret_gran%nStd_bot_mid_bndry)
      if (statn .ne. 0) &
        print *, "Error ", statn, " reading field ", &
                  "nStd_bot_mid_bndry"

      edge(2) = 45
      edge(1) = 30
      statn = SWrdfld(swid, "PBest", &
                   start, stride, edge, &
                   airs_ret_gran%PBest)
      if (statn .ne. 0) &
        print *, "Error ", statn, " reading field ", &
                  "PBest"

      edge(2) = 45
      edge(1) = 30
      statn = SWrdfld(swid, "PGood", &
                   start, stride, edge, &
                   airs_ret_gran%PGood)
      if (statn .ne. 0) &
        print *, "Error ", statn, " reading field ", &
                  "PGood"

      edge(2) = 45
      edge(1) = 30
      statn = SWrdfld(swid, "nBestStd", &
                   start, stride, edge, &
                   airs_ret_gran%nBestStd)
      if (statn .ne. 0) &
        print *, "Error ", statn, " reading field ", &
                  "nBestStd"

      edge(2) = 45
      edge(1) = 30
      statn = SWrdfld(swid, "nGoodStd", &
                   start, stride, edge, &
                   airs_ret_gran%nGoodStd)
      if (statn .ne. 0) &
        print *, "Error ", statn, " reading field ", &
                  "nGoodStd"

      edge(2) = 45
      edge(1) = 30
      statn = SWrdfld(swid, "Qual_Temp_Profile_Top", &
                   start, stride, edge, &
                   airs_ret_gran%Qual_Temp_Profile_Top)
      if (statn .ne. 0) &
        print *, "Error ", statn, " reading field ", &
                  "Qual_Temp_Profile_Top"

      edge(2) = 45
      edge(1) = 30
      statn = SWrdfld(swid, "Qual_Temp_Profile_Mid", &
                   start, stride, edge, &
                   airs_ret_gran%Qual_Temp_Profile_Mid)
      if (statn .ne. 0) &
        print *, "Error ", statn, " reading field ", &
                  "Qual_Temp_Profile_Mid"

      edge(2) = 45
      edge(1) = 30
      statn = SWrdfld(swid, "Qual_Temp_Profile_Bot", &
                   start, stride, edge, &
                   airs_ret_gran%Qual_Temp_Profile_Bot)
      if (statn .ne. 0) &
        print *, "Error ", statn, " reading field ", &
                  "Qual_Temp_Profile_Bot"

      edge(3) = 45
      edge(2) = 30
      edge(1) = 28
      statn = SWrdfld(swid, "TAirStd", &
                   start, stride, edge, &
                   airs_ret_gran%TAirStd)
      if (statn .ne. 0) &
        print *, "Error ", statn, " reading field ", &
                  "TAirStd"

      edge(3) = 45
      edge(2) = 30
      edge(1) = 28
      statn = SWrdfld(swid, "TAirStdErr", &
                   start, stride, edge, &
                   airs_ret_gran%TAirStdErr)
      if (statn .ne. 0) &
        print *, "Error ", statn, " reading field ", &
                  "TAirStdErr"

      edge(2) = 45
      edge(1) = 30
      statn = SWrdfld(swid, "TSurfAir", &
                   start, stride, edge, &
                   airs_ret_gran%TSurfAir)
      if (statn .ne. 0) &
        print *, "Error ", statn, " reading field ", &
                  "TSurfAir"

      edge(2) = 45
      edge(1) = 30
      statn = SWrdfld(swid, "TSurfAirErr", &
                   start, stride, edge, &
                   airs_ret_gran%TSurfAirErr)
      if (statn .ne. 0) &
        print *, "Error ", statn, " reading field ", &
                  "TSurfAirErr"

      edge(2) = 45
      edge(1) = 30
      statn = SWrdfld(swid, "Qual_Surf", &
                   start, stride, edge, &
                   airs_ret_gran%Qual_Surf)
      if (statn .ne. 0) &
        print *, "Error ", statn, " reading field ", &
                  "Qual_Surf"

      edge(2) = 45
      edge(1) = 30
      statn = SWrdfld(swid, "TSurfStd", &
                   start, stride, edge, &
                   airs_ret_gran%TSurfStd)
      if (statn .ne. 0) &
        print *, "Error ", statn, " reading field ", &
                  "TSurfStd"

      edge(2) = 45
      edge(1) = 30
      statn = SWrdfld(swid, "TSurfStdErr", &
                   start, stride, edge, &
                   airs_ret_gran%TSurfStdErr)
      if (statn .ne. 0) &
        print *, "Error ", statn, " reading field ", &
                  "TSurfStdErr"

      edge(2) = 45
      edge(1) = 30
      statn = SWrdfld(swid, "numHingeSurf", &
                   start, stride, edge, &
                   airs_ret_gran%numHingeSurf)
      if (statn .ne. 0) &
        print *, "Error ", statn, " reading field ", &
                  "numHingeSurf"

      edge(3) = 45
      edge(2) = 30
      edge(1) = 100
      statn = SWrdfld(swid, "freqEmis", &
                   start, stride, edge, &
                   airs_ret_gran%freqEmis)
      if (statn .ne. 0) &
        print *, "Error ", statn, " reading field ", &
                  "freqEmis"

      edge(3) = 45
      edge(2) = 30
      edge(1) = 100
      statn = SWrdfld(swid, "emisIRStd", &
                   start, stride, edge, &
                   airs_ret_gran%emisIRStd)
      if (statn .ne. 0) &
        print *, "Error ", statn, " reading field ", &
                  "emisIRStd"

      edge(3) = 45
      edge(2) = 30
      edge(1) = 100
      statn = SWrdfld(swid, "emisIRStdErr", &
                   start, stride, edge, &
                   airs_ret_gran%emisIRStdErr)
      if (statn .ne. 0) &
        print *, "Error ", statn, " reading field ", &
                  "emisIRStdErr"

      edge(2) = 45
      edge(1) = 30
      statn = SWrdfld(swid, "Qual_MW_Only_Temp_Strat", &
                   start, stride, edge, &
                   airs_ret_gran%Qual_MW_Only_Temp_Strat)
      if (statn .ne. 0) &
        print *, "Error ", statn, " reading field ", &
                  "Qual_MW_Only_Temp_Strat"

      edge(2) = 45
      edge(1) = 30
      statn = SWrdfld(swid, "Qual_MW_Only_Temp_Tropo", &
                   start, stride, edge, &
                   airs_ret_gran%Qual_MW_Only_Temp_Tropo)
      if (statn .ne. 0) &
        print *, "Error ", statn, " reading field ", &
                  "Qual_MW_Only_Temp_Tropo"

      edge(3) = 45
      edge(2) = 30
      edge(1) = 28
      statn = SWrdfld(swid, "TAirMWOnlyStd", &
                   start, stride, edge, &
                   airs_ret_gran%TAirMWOnlyStd)
      if (statn .ne. 0) &
        print *, "Error ", statn, " reading field ", &
                  "TAirMWOnlyStd"

      edge(2) = 45
      edge(1) = 30
      statn = SWrdfld(swid, "MWSurfClass", &
                   start, stride, edge, &
                   airs_ret_gran%MWSurfClass)
      if (statn .ne. 0) &
        print *, "Error ", statn, " reading field ", &
                  "MWSurfClass"

      edge(3) = 45
      edge(2) = 30
      edge(1) = 7
      statn = SWrdfld(swid, "sfcTbMWStd", &
                   start, stride, edge, &
                   airs_ret_gran%sfcTbMWStd)
      if (statn .ne. 0) &
        print *, "Error ", statn, " reading field ", &
                  "sfcTbMWStd"

      edge(3) = 45
      edge(2) = 30
      edge(1) = 7
      statn = SWrdfld(swid, "EmisMWStd", &
                   start, stride, edge, &
                   airs_ret_gran%EmisMWStd)
      if (statn .ne. 0) &
        print *, "Error ", statn, " reading field ", &
                  "EmisMWStd"

      edge(3) = 45
      edge(2) = 30
      edge(1) = 7
      statn = SWrdfld(swid, "EmisMWStdErr", &
                   start, stride, edge, &
                   airs_ret_gran%EmisMWStdErr)
      if (statn .ne. 0) &
        print *, "Error ", statn, " reading field ", &
                  "EmisMWStdErr"

      edge(2) = 45
      edge(1) = 30
      statn = SWrdfld(swid, "Qual_MW_Only_H2O", &
                   start, stride, edge, &
                   airs_ret_gran%Qual_MW_Only_H2O)
      if (statn .ne. 0) &
        print *, "Error ", statn, " reading field ", &
                  "Qual_MW_Only_H2O"

      edge(2) = 45
      edge(1) = 30
      statn = SWrdfld(swid, "totH2OMWOnlyStd", &
                   start, stride, edge, &
                   airs_ret_gran%totH2OMWOnlyStd)
      if (statn .ne. 0) &
        print *, "Error ", statn, " reading field ", &
                  "totH2OMWOnlyStd"

      edge(2) = 45
      edge(1) = 30
      statn = SWrdfld(swid, "Qual_H2O", &
                   start, stride, edge, &
                   airs_ret_gran%Qual_H2O)
      if (statn .ne. 0) &
        print *, "Error ", statn, " reading field ", &
                  "Qual_H2O"

      edge(3) = 45
      edge(2) = 30
      edge(1) = 14
      statn = SWrdfld(swid, "H2OMMRStd", &
                   start, stride, edge, &
                   airs_ret_gran%H2OMMRStd)
      if (statn .ne. 0) &
        print *, "Error ", statn, " reading field ", &
                  "H2OMMRStd"

      edge(3) = 45
      edge(2) = 30
      edge(1) = 14
      statn = SWrdfld(swid, "H2OMMRStdErr", &
                   start, stride, edge, &
                   airs_ret_gran%H2OMMRStdErr)
      if (statn .ne. 0) &
        print *, "Error ", statn, " reading field ", &
                  "H2OMMRStdErr"

      edge(2) = 45
      edge(1) = 30
      statn = SWrdfld(swid, "totH2OStd", &
                   start, stride, edge, &
                   airs_ret_gran%totH2OStd)
      if (statn .ne. 0) &
        print *, "Error ", statn, " reading field ", &
                  "totH2OStd"

      edge(2) = 45
      edge(1) = 30
      statn = SWrdfld(swid, "totH2OStdErr", &
                   start, stride, edge, &
                   airs_ret_gran%totH2OStdErr)
      if (statn .ne. 0) &
        print *, "Error ", statn, " reading field ", &
                  "totH2OStdErr"

      edge(3) = 45
      edge(2) = 30
      edge(1) = 14
      statn = SWrdfld(swid, "H2OMMRSat", &
                   start, stride, edge, &
                   airs_ret_gran%H2OMMRSat)
      if (statn .ne. 0) &
        print *, "Error ", statn, " reading field ", &
                  "H2OMMRSat"

      edge(3) = 45
      edge(2) = 30
      edge(1) = 14
      statn = SWrdfld(swid, "H2OMMRSat_liquid", &
                   start, stride, edge, &
                   airs_ret_gran%H2OMMRSat_liquid)
      if (statn .ne. 0) &
        print *, "Error ", statn, " reading field ", &
                  "H2OMMRSat_liquid"

      edge(2) = 45
      edge(1) = 30
      statn = SWrdfld(swid, "num_H2O_Func", &
                   start, stride, edge, &
                   airs_ret_gran%num_H2O_Func)
      if (statn .ne. 0) &
        print *, "Error ", statn, " reading field ", &
                  "num_H2O_Func"

      edge(3) = 45
      edge(2) = 30
      edge(1) = 11
      statn = SWrdfld(swid, "H2O_verticality", &
                   start, stride, edge, &
                   airs_ret_gran%H2O_verticality)
      if (statn .ne. 0) &
        print *, "Error ", statn, " reading field ", &
                  "H2O_verticality"

      edge(2) = 45
      edge(1) = 30
      statn = SWrdfld(swid, "Qual_O3", &
                   start, stride, edge, &
                   airs_ret_gran%Qual_O3)
      if (statn .ne. 0) &
        print *, "Error ", statn, " reading field ", &
                  "Qual_O3"

      edge(2) = 45
      edge(1) = 30
      statn = SWrdfld(swid, "totO3Std", &
                   start, stride, edge, &
                   airs_ret_gran%totO3Std)
      if (statn .ne. 0) &
        print *, "Error ", statn, " reading field ", &
                  "totO3Std"

      edge(2) = 45
      edge(1) = 30
      statn = SWrdfld(swid, "totO3StdErr", &
                   start, stride, edge, &
                   airs_ret_gran%totO3StdErr)
      if (statn .ne. 0) &
        print *, "Error ", statn, " reading field ", &
                  "totO3StdErr"

      edge(3) = 45
      edge(2) = 30
      edge(1) = 28
      statn = SWrdfld(swid, "O3VMRStd", &
                   start, stride, edge, &
                   airs_ret_gran%O3VMRStd)
      if (statn .ne. 0) &
        print *, "Error ", statn, " reading field ", &
                  "O3VMRStd"

      edge(3) = 45
      edge(2) = 30
      edge(1) = 28
      statn = SWrdfld(swid, "O3VMRStdErr", &
                   start, stride, edge, &
                   airs_ret_gran%O3VMRStdErr)
      if (statn .ne. 0) &
        print *, "Error ", statn, " reading field ", &
                  "O3VMRStdErr"

      edge(2) = 45
      edge(1) = 30
      statn = SWrdfld(swid, "num_O3_Func", &
                   start, stride, edge, &
                   airs_ret_gran%num_O3_Func)
      if (statn .ne. 0) &
        print *, "Error ", statn, " reading field ", &
                  "num_O3_Func"

      edge(3) = 45
      edge(2) = 30
      edge(1) = 9
      statn = SWrdfld(swid, "O3_verticality", &
                   start, stride, edge, &
                   airs_ret_gran%O3_verticality)
      if (statn .ne. 0) &
        print *, "Error ", statn, " reading field ", &
                  "O3_verticality"

      edge(2) = 45
      edge(1) = 30
      statn = SWrdfld(swid, "Qual_CO", &
                   start, stride, edge, &
                   airs_ret_gran%Qual_CO)
      if (statn .ne. 0) &
        print *, "Error ", statn, " reading field ", &
                  "Qual_CO"

      edge(2) = 45
      edge(1) = 30
      statn = SWrdfld(swid, "CO_total_column", &
                   start, stride, edge, &
                   airs_ret_gran%CO_total_column)
      if (statn .ne. 0) &
        print *, "Error ", statn, " reading field ", &
                  "CO_total_column"

      edge(2) = 45
      edge(1) = 30
      statn = SWrdfld(swid, "num_CO_Func", &
                   start, stride, edge, &
                   airs_ret_gran%num_CO_Func)
      if (statn .ne. 0) &
        print *, "Error ", statn, " reading field ", &
                  "num_CO_Func"

      edge(1) = 9
      statn = SWrdfld(swid, "CO_trapezoid_layers", &
                   start, stride, edge, &
                   airs_ret_gran%CO_trapezoid_layers)
      if (statn .ne. 0) &
        print *, "Error ", statn, " reading field ", &
                  "CO_trapezoid_layers"

      edge(3) = 45
      edge(2) = 30
      edge(1) = 9
      statn = SWrdfld(swid, "CO_eff_press", &
                   start, stride, edge, &
                   airs_ret_gran%CO_eff_press)
      if (statn .ne. 0) &
        print *, "Error ", statn, " reading field ", &
                  "CO_eff_press"

      edge(3) = 45
      edge(2) = 30
      edge(1) = 9
      statn = SWrdfld(swid, "CO_VMR_eff", &
                   start, stride, edge, &
                   airs_ret_gran%CO_VMR_eff)
      if (statn .ne. 0) &
        print *, "Error ", statn, " reading field ", &
                  "CO_VMR_eff"

      edge(3) = 45
      edge(2) = 30
      edge(1) = 9
      statn = SWrdfld(swid, "CO_VMR_eff_err", &
                   start, stride, edge, &
                   airs_ret_gran%CO_VMR_eff_err)
      if (statn .ne. 0) &
        print *, "Error ", statn, " reading field ", &
                  "CO_VMR_eff_err"

      edge(3) = 45
      edge(2) = 30
      edge(1) = 9
      statn = SWrdfld(swid, "CO_verticality", &
                   start, stride, edge, &
                   airs_ret_gran%CO_verticality)
      if (statn .ne. 0) &
        print *, "Error ", statn, " reading field ", &
                  "CO_verticality"

      edge(2) = 45
      edge(1) = 30
      statn = SWrdfld(swid, "CO_dof", &
                   start, stride, edge, &
                   airs_ret_gran%CO_dof)
      if (statn .ne. 0) &
        print *, "Error ", statn, " reading field ", &
                  "CO_dof"

      edge(2) = 45
      edge(1) = 30
      statn = SWrdfld(swid, "Qual_CH4", &
                   start, stride, edge, &
                   airs_ret_gran%Qual_CH4)
      if (statn .ne. 0) &
        print *, "Error ", statn, " reading field ", &
                  "Qual_CH4"

      edge(2) = 45
      edge(1) = 30
      statn = SWrdfld(swid, "CH4_total_column", &
                   start, stride, edge, &
                   airs_ret_gran%CH4_total_column)
      if (statn .ne. 0) &
        print *, "Error ", statn, " reading field ", &
                  "CH4_total_column"

      edge(2) = 45
      edge(1) = 30
      statn = SWrdfld(swid, "num_CH4_Func", &
                   start, stride, edge, &
                   airs_ret_gran%num_CH4_Func)
      if (statn .ne. 0) &
        print *, "Error ", statn, " reading field ", &
                  "num_CH4_Func"

      edge(1) = 7
      statn = SWrdfld(swid, "CH4_trapezoid_layers", &
                   start, stride, edge, &
                   airs_ret_gran%CH4_trapezoid_layers)
      if (statn .ne. 0) &
        print *, "Error ", statn, " reading field ", &
                  "CH4_trapezoid_layers"

      edge(3) = 45
      edge(2) = 30
      edge(1) = 7
      statn = SWrdfld(swid, "CH4_eff_press", &
                   start, stride, edge, &
                   airs_ret_gran%CH4_eff_press)
      if (statn .ne. 0) &
        print *, "Error ", statn, " reading field ", &
                  "CH4_eff_press"

      edge(3) = 45
      edge(2) = 30
      edge(1) = 7
      statn = SWrdfld(swid, "CH4_VMR_eff", &
                   start, stride, edge, &
                   airs_ret_gran%CH4_VMR_eff)
      if (statn .ne. 0) &
        print *, "Error ", statn, " reading field ", &
                  "CH4_VMR_eff"

      edge(3) = 45
      edge(2) = 30
      edge(1) = 7
      statn = SWrdfld(swid, "CH4_VMR_eff_err", &
                   start, stride, edge, &
                   airs_ret_gran%CH4_VMR_eff_err)
      if (statn .ne. 0) &
        print *, "Error ", statn, " reading field ", &
                  "CH4_VMR_eff_err"

      edge(3) = 45
      edge(2) = 30
      edge(1) = 7
      statn = SWrdfld(swid, "CH4_verticality", &
                   start, stride, edge, &
                   airs_ret_gran%CH4_verticality)
      if (statn .ne. 0) &
        print *, "Error ", statn, " reading field ", &
                  "CH4_verticality"

      edge(2) = 45
      edge(1) = 30
      statn = SWrdfld(swid, "CH4_dof", &
                   start, stride, edge, &
                   airs_ret_gran%CH4_dof)
      if (statn .ne. 0) &
        print *, "Error ", statn, " reading field ", &
                  "CH4_dof"

      edge(2) = 45
      edge(1) = 30
      statn = SWrdfld(swid, "PTropopause", &
                   start, stride, edge, &
                   airs_ret_gran%PTropopause)
      if (statn .ne. 0) &
        print *, "Error ", statn, " reading field ", &
                  "PTropopause"

      edge(2) = 45
      edge(1) = 30
      statn = SWrdfld(swid, "T_Tropopause", &
                   start, stride, edge, &
                   airs_ret_gran%T_Tropopause)
      if (statn .ne. 0) &
        print *, "Error ", statn, " reading field ", &
                  "T_Tropopause"

      edge(2) = 45
      edge(1) = 30
      statn = SWrdfld(swid, "GP_Tropopause", &
                   start, stride, edge, &
                   airs_ret_gran%GP_Tropopause)
      if (statn .ne. 0) &
        print *, "Error ", statn, " reading field ", &
                  "GP_Tropopause"

      edge(3) = 45
      edge(2) = 30
      edge(1) = 28
      statn = SWrdfld(swid, "GP_Height", &
                   start, stride, edge, &
                   airs_ret_gran%GP_Height)
      if (statn .ne. 0) &
        print *, "Error ", statn, " reading field ", &
                  "GP_Height"

      edge(3) = 45
      edge(2) = 30
      edge(1) = 28
      statn = SWrdfld(swid, "GP_Height_MWOnly", &
                   start, stride, edge, &
                   airs_ret_gran%GP_Height_MWOnly)
      if (statn .ne. 0) &
        print *, "Error ", statn, " reading field ", &
                  "GP_Height_MWOnly"

      edge(2) = 45
      edge(1) = 30
      statn = SWrdfld(swid, "GP_Surface", &
                   start, stride, edge, &
                   airs_ret_gran%GP_Surface)
      if (statn .ne. 0) &
        print *, "Error ", statn, " reading field ", &
                  "GP_Surface"

      edge(2) = 45
      edge(1) = 30
      statn = SWrdfld(swid, "Qual_Cloud_OLR", &
                   start, stride, edge, &
                   airs_ret_gran%Qual_Cloud_OLR)
      if (statn .ne. 0) &
        print *, "Error ", statn, " reading field ", &
                  "Qual_Cloud_OLR"

      edge(2) = 45
      edge(1) = 30
      statn = SWrdfld(swid, "numCloud", &
                   start, stride, edge, &
                   airs_ret_gran%numCloud)
      if (statn .ne. 0) &
        print *, "Error ", statn, " reading field ", &
                  "numCloud"

      edge(3) = 45
      edge(2) = 30
      edge(1) = 2
      statn = SWrdfld(swid, "TCldTopStd", &
                   start, stride, edge, &
                   airs_ret_gran%TCldTopStd)
      if (statn .ne. 0) &
        print *, "Error ", statn, " reading field ", &
                  "TCldTopStd"

      edge(3) = 45
      edge(2) = 30
      edge(1) = 2
      statn = SWrdfld(swid, "TCldTopStdErr", &
                   start, stride, edge, &
                   airs_ret_gran%TCldTopStdErr)
      if (statn .ne. 0) &
        print *, "Error ", statn, " reading field ", &
                  "TCldTopStdErr"

      edge(3) = 45
      edge(2) = 30
      edge(1) = 2
      statn = SWrdfld(swid, "PCldTopStd", &
                   start, stride, edge, &
                   airs_ret_gran%PCldTopStd)
      if (statn .ne. 0) &
        print *, "Error ", statn, " reading field ", &
                  "PCldTopStd"

      edge(3) = 45
      edge(2) = 30
      edge(1) = 2
      statn = SWrdfld(swid, "PCldTopStdErr", &
                   start, stride, edge, &
                   airs_ret_gran%PCldTopStdErr)
      if (statn .ne. 0) &
        print *, "Error ", statn, " reading field ", &
                  "PCldTopStdErr"

      edge(5) = 45
      edge(4) = 30
      edge(3) = 3
      edge(2) = 3
      edge(1) = 2
      statn = SWrdfld(swid, "CldFrcStd", &
                   start, stride, edge, &
                   airs_ret_gran%CldFrcStd)
      if (statn .ne. 0) &
        print *, "Error ", statn, " reading field ", &
                  "CldFrcStd"

      edge(5) = 45
      edge(4) = 30
      edge(3) = 3
      edge(2) = 3
      edge(1) = 2
      statn = SWrdfld(swid, "CldFrcStdErr", &
                   start, stride, edge, &
                   airs_ret_gran%CldFrcStdErr)
      if (statn .ne. 0) &
        print *, "Error ", statn, " reading field ", &
                  "CldFrcStdErr"

      edge(2) = 45
      edge(1) = 30
      statn = SWrdfld(swid, "olr", &
                   start, stride, edge, &
                   airs_ret_gran%olr)
      if (statn .ne. 0) &
        print *, "Error ", statn, " reading field ", &
                  "olr"

      edge(2) = 45
      edge(1) = 30
      statn = SWrdfld(swid, "olr_err", &
                   start, stride, edge, &
                   airs_ret_gran%olr_err)
      if (statn .ne. 0) &
        print *, "Error ", statn, " reading field ", &
                  "olr_err"

      edge(2) = 45
      edge(1) = 30
      statn = SWrdfld(swid, "Qual_clrolr", &
                   start, stride, edge, &
                   airs_ret_gran%Qual_clrolr)
      if (statn .ne. 0) &
        print *, "Error ", statn, " reading field ", &
                  "Qual_clrolr"

      edge(2) = 45
      edge(1) = 30
      statn = SWrdfld(swid, "clrolr", &
                   start, stride, edge, &
                   airs_ret_gran%clrolr)
      if (statn .ne. 0) &
        print *, "Error ", statn, " reading field ", &
                  "clrolr"

      edge(2) = 45
      edge(1) = 30
      statn = SWrdfld(swid, "clrolr_err", &
                   start, stride, edge, &
                   airs_ret_gran%clrolr_err)
      if (statn .ne. 0) &
        print *, "Error ", statn, " reading field ", &
                  "clrolr_err"

      edge(4) = 45
      edge(3) = 30
      edge(2) = 3
      edge(1) = 3
      statn = SWrdfld(swid, "dust_flag", &
                   start, stride, edge, &
                   airs_ret_gran%dust_flag)
      if (statn .ne. 0) &
        print *, "Error ", statn, " reading field ", &
                  "dust_flag"

      edge(4) = 45
      edge(3) = 30
      edge(2) = 3
      edge(1) = 3
      statn = SWrdfld(swid, "spectral_clear_indicator", &
                   start, stride, edge, &
                   airs_ret_gran%spectral_clear_indicator)
      if (statn .ne. 0) &
        print *, "Error ", statn, " reading field ", &
                  "spectral_clear_indicator"

      edge(2) = 45
      edge(1) = 30
      statn = SWrdfld(swid, "num_clear_spectral_indicator", &
                   start, stride, edge, &
                   airs_ret_gran%num_clear_spectral_indicator)
      if (statn .ne. 0) &
        print *, "Error ", statn, " reading field ", &
                  "num_clear_spectral_indicator"

      edge(2) = 45
      edge(1) = 30
      statn = SWrdfld(swid, "CC_noise_eff_amp_factor", &
                   start, stride, edge, &
                   airs_ret_gran%CC_noise_eff_amp_factor)
      if (statn .ne. 0) &
        print *, "Error ", statn, " reading field ", &
                  "CC_noise_eff_amp_factor"

      edge(2) = 45
      edge(1) = 30
      statn = SWrdfld(swid, "CC1_noise_eff_amp_factor", &
                   start, stride, edge, &
                   airs_ret_gran%CC1_noise_eff_amp_factor)
      if (statn .ne. 0) &
        print *, "Error ", statn, " reading field ", &
                  "CC1_noise_eff_amp_factor"

      edge(2) = 45
      edge(1) = 30
      statn = SWrdfld(swid, "totCldH2OStd", &
                   start, stride, edge, &
                   airs_ret_gran%totCldH2OStd)
      if (statn .ne. 0) &
        print *, "Error ", statn, " reading field ", &
                  "totCldH2OStd"

      edge(2) = 45
      edge(1) = 30
      statn = SWrdfld(swid, "totCldH2OStdErr", &
                   start, stride, edge, &
                   airs_ret_gran%totCldH2OStdErr)
      if (statn .ne. 0) &
        print *, "Error ", statn, " reading field ", &
                  "totCldH2OStdErr"

      edge(2) = 45
      edge(1) = 30
      statn = SWrdfld(swid, "CC1_Resid", &
                   start, stride, edge, &
                   airs_ret_gran%CC1_Resid)
      if (statn .ne. 0) &
        print *, "Error ", statn, " reading field ", &
                  "CC1_Resid"

      edge(2) = 45
      edge(1) = 30
      statn = SWrdfld(swid, "CCfinal_Resid", &
                   start, stride, edge, &
                   airs_ret_gran%CCfinal_Resid)
      if (statn .ne. 0) &
        print *, "Error ", statn, " reading field ", &
                  "CCfinal_Resid"

      edge(2) = 45
      edge(1) = 30
      statn = SWrdfld(swid, "CCfinal_Noise_Amp", &
                   start, stride, edge, &
                   airs_ret_gran%CCfinal_Noise_Amp)
      if (statn .ne. 0) &
        print *, "Error ", statn, " reading field ", &
                  "CCfinal_Noise_Amp"

      edge(2) = 45
      edge(1) = 30
      statn = SWrdfld(swid, "Tdiff_IR_MW_ret", &
                   start, stride, edge, &
                   airs_ret_gran%Tdiff_IR_MW_ret)
      if (statn .ne. 0) &
        print *, "Error ", statn, " reading field ", &
                  "Tdiff_IR_MW_ret"

      edge(2) = 45
      edge(1) = 30
      statn = SWrdfld(swid, "Tdiff_IR_4CC1", &
                   start, stride, edge, &
                   airs_ret_gran%Tdiff_IR_4CC1)
      if (statn .ne. 0) &
        print *, "Error ", statn, " reading field ", &
                  "Tdiff_IR_4CC1"

      edge(2) = 45
      edge(1) = 30
      statn = SWrdfld(swid, "TSurfdiff_IR_4CC1", &
                   start, stride, edge, &
                   airs_ret_gran%TSurfdiff_IR_4CC1)
      if (statn .ne. 0) &
        print *, "Error ", statn, " reading field ", &
                  "TSurfdiff_IR_4CC1"

      edge(2) = 45
      edge(1) = 30
      statn = SWrdfld(swid, "TSurfdiff_IR_4CC2", &
                   start, stride, edge, &
                   airs_ret_gran%TSurfdiff_IR_4CC2)
      if (statn .ne. 0) &
        print *, "Error ", statn, " reading field ", &
                  "TSurfdiff_IR_4CC2"

      edge(2) = 45
      edge(1) = 30
      statn = SWrdfld(swid, "AMSU_Chans_Resid", &
                   start, stride, edge, &
                   airs_ret_gran%AMSU_Chans_Resid)
      if (statn .ne. 0) &
        print *, "Error ", statn, " reading field ", &
                  "AMSU_Chans_Resid"

      edge(2) = 45
      edge(1) = 30
      statn = SWrdfld(swid, "TotCld_4_CCfinal", &
                   start, stride, edge, &
                   airs_ret_gran%TotCld_4_CCfinal)
      if (statn .ne. 0) &
        print *, "Error ", statn, " reading field ", &
                  "TotCld_4_CCfinal"

      edge(2) = 45
      edge(1) = 30
      statn = SWrdfld(swid, "Surf_Resid_Ratio", &
                   start, stride, edge, &
                   airs_ret_gran%Surf_Resid_Ratio)
      if (statn .ne. 0) &
        print *, "Error ", statn, " reading field ", &
                  "Surf_Resid_Ratio"

      edge(2) = 45
      edge(1) = 30
      statn = SWrdfld(swid, "Temp_Resid_Ratio", &
                   start, stride, edge, &
                   airs_ret_gran%Temp_Resid_Ratio)
      if (statn .ne. 0) &
        print *, "Error ", statn, " reading field ", &
                  "Temp_Resid_Ratio"

      edge(2) = 45
      edge(1) = 30
      statn = SWrdfld(swid, "Water_Resid_Ratio", &
                   start, stride, edge, &
                   airs_ret_gran%Water_Resid_Ratio)
      if (statn .ne. 0) &
        print *, "Error ", statn, " reading field ", &
                  "Water_Resid_Ratio"

      edge(2) = 45
      edge(1) = 30
      statn = SWrdfld(swid, "Cloud_Resid_Ratio", &
                   start, stride, edge, &
                   airs_ret_gran%Cloud_Resid_Ratio)
      if (statn .ne. 0) &
        print *, "Error ", statn, " reading field ", &
                  "Cloud_Resid_Ratio"

      edge(2) = 45
      edge(1) = 30
      statn = SWrdfld(swid, "O3_Resid_Ratio", &
                   start, stride, edge, &
                   airs_ret_gran%O3_Resid_Ratio)
      if (statn .ne. 0) &
        print *, "Error ", statn, " reading field ", &
                  "O3_Resid_Ratio"

      edge(2) = 45
      edge(1) = 30
      statn = SWrdfld(swid, "CO_Resid_Ratio", &
                   start, stride, edge, &
                   airs_ret_gran%CO_Resid_Ratio)
      if (statn .ne. 0) &
        print *, "Error ", statn, " reading field ", &
                  "CO_Resid_Ratio"

      edge(2) = 45
      edge(1) = 30
      statn = SWrdfld(swid, "CH4_Resid_Ratio", &
                   start, stride, edge, &
                   airs_ret_gran%CH4_Resid_Ratio)
      if (statn .ne. 0) &
        print *, "Error ", statn, " reading field ", &
                  "CH4_Resid_Ratio"

      edge(2) = 45
      edge(1) = 30
      statn = SWrdfld(swid, "MWCheck_Resid_Ratio", &
                   start, stride, edge, &
                   airs_ret_gran%MWCheck_Resid_Ratio)
      if (statn .ne. 0) &
        print *, "Error ", statn, " reading field ", &
                  "MWCheck_Resid_Ratio"

      edge(2) = 45
      edge(1) = 30
      statn = SWrdfld(swid, "O3_dof", &
                   start, stride, edge, &
                   airs_ret_gran%O3_dof)
      if (statn .ne. 0) &
        print *, "Error ", statn, " reading field ", &
                  "O3_dof"

      edge(2) = 45
      edge(1) = 30
      statn = SWrdfld(swid, "all_spots_avg", &
                   start, stride, edge, &
                   airs_ret_gran%all_spots_avg)
      if (statn .ne. 0) &
        print *, "Error ", statn, " reading field ", &
                  "all_spots_avg"

      edge(2) = 45
      edge(1) = 30
      statn = SWrdfld(swid, "MW_ret_used", &
                   start, stride, edge, &
                   airs_ret_gran%MW_ret_used)
      if (statn .ne. 0) &
        print *, "Error ", statn, " reading field ", &
                  "MW_ret_used"

      edge(2) = 45
      edge(1) = 30
      statn = SWrdfld(swid, "Initial_CC_score", &
                   start, stride, edge, &
                   airs_ret_gran%Initial_CC_score)
      if (statn .ne. 0) &
        print *, "Error ", statn, " reading field ", &
                  "Initial_CC_score"

      edge(2) = 45
      edge(1) = 30
      statn = SWrdfld(swid, "retrieval_type", &
                   start, stride, edge, &
                   airs_ret_gran%retrieval_type)
      if (statn .ne. 0) &
        print *, "Error ", statn, " reading field ", &
                  "retrieval_type"

      edge(2) = 45
      edge(1) = 30
      statn = SWrdfld(swid, "Startup", &
                   start, stride, edge, &
                   airs_ret_gran%Startup)
      if (statn .ne. 0) &
        print *, "Error ", statn, " reading field ", &
                  "Startup"


      ! Final clean-up
      statn = swdetach(swid)
      if (statn .ne. 0) &
        print *, "Error detaching from input file ", file_name
      statn = swclose(fid)
      if (statn .ne. 0) &
        print *, "Error closing input file ", file_name

      return
end subroutine

end module

! <next few lines under version control, do not edit>
! $URL$
! $Id$
! $Revision$
! $Date$
