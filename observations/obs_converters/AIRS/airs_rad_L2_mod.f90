! This code is not protected by the DART copyright agreement.
! DART $Id$

! adapted from original JPL, now GSFC/NASA code, example AIRS readers

module airs_rad_L2_mod

! the contents of this file are an amalgam of:
!  airs_cc_rad_typ.inc
!  airs_cc_rad_struct.inc
!  airs_cc_rad_rdr.f
! from the gsfc example reader program.
! although in several cases they were modified to use
! fortran 90 derived types and free format continuation
! lines, and to not use the unknown syntax 'double precision'.

implicit none
public

integer, parameter :: AIRS_CC_RAD_GEOXTRACK  =   30
integer, parameter :: AIRS_CC_RAD_GEOTRACK   =   45
integer, parameter :: AIRS_CC_RAD_CHANNEL    = 2378
integer, parameter :: AIRS_CC_RAD_AIRSXTRACK =    3
integer, parameter :: AIRS_CC_RAD_AIRSTRACK  =    3
integer, parameter :: AIRS_CC_RAD_MODULE     =   17

! Record holds an entire granule of airs_cc_rad
type airs_granule_type

! Attributes
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
   real           start_sec
   integer      start_orbit
   integer      end_orbit
   integer      orbit_path
   integer      start_orbit_row
   integer      end_orbit_row
   integer      granule_number
   integer      num_scansets
   integer      num_scanlines
   real*8 start_Latitude
   real*8 start_Longitude
   real*8 start_Time
   real*8 end_Latitude
   real*8 end_Longitude
   real*8 end_Time
   real             eq_x_longitude
   real*8 eq_x_tai
   integer*2      LonGranuleCen
   integer*2      LatGranuleCen
   integer*2      LocTimeGranuleCen
   integer*2      num_fpe
   integer        orbitgeoqa
   integer*2      num_satgeoqa
   integer*2      num_glintgeoqa
   integer*2      num_moongeoqa
   integer*2      num_ftptgeoqa
   integer*2      num_zengeoqa
   integer*2      num_demgeoqa
   byte           CalGranSummary
   integer*2      DCR_scan
   character*256  granules_present_L1B

! Geolocation fields
   real*8 Latitude(AIRS_CC_RAD_GEOXTRACK,AIRS_CC_RAD_GEOTRACK)
   real*8 Longitude(AIRS_CC_RAD_GEOXTRACK,AIRS_CC_RAD_GEOTRACK)
   real*8 Time(AIRS_CC_RAD_GEOXTRACK,AIRS_CC_RAD_GEOTRACK)

! Data Fields
   real radiances(AIRS_CC_RAD_CHANNEL,AIRS_CC_RAD_GEOXTRACK,AIRS_CC_RAD_GEOTRACK)
   integer*2 radiances_QC(AIRS_CC_RAD_CHANNEL,AIRS_CC_RAD_GEOXTRACK,AIRS_CC_RAD_GEOTRACK)
   real radiance_err(AIRS_CC_RAD_CHANNEL,AIRS_CC_RAD_GEOXTRACK,AIRS_CC_RAD_GEOTRACK)
   real CldClearParam(AIRS_CC_RAD_AIRSXTRACK,AIRS_CC_RAD_AIRSTRACK,AIRS_CC_RAD_GEOXTRACK,AIRS_CC_RAD_GEOTRACK)
   real nominal_freq(AIRS_CC_RAD_CHANNEL)
   real scanang(AIRS_CC_RAD_GEOXTRACK,AIRS_CC_RAD_GEOTRACK)
   real satheight(AIRS_CC_RAD_GEOTRACK)
   real satroll(AIRS_CC_RAD_GEOTRACK)
   real satpitch(AIRS_CC_RAD_GEOTRACK)
   real satyaw(AIRS_CC_RAD_GEOTRACK)
   real satzen(AIRS_CC_RAD_GEOXTRACK,AIRS_CC_RAD_GEOTRACK)
   real satazi(AIRS_CC_RAD_GEOXTRACK,AIRS_CC_RAD_GEOTRACK)
   real solzen(AIRS_CC_RAD_GEOXTRACK,AIRS_CC_RAD_GEOTRACK)
   real solazi(AIRS_CC_RAD_GEOXTRACK,AIRS_CC_RAD_GEOTRACK)
   real glintlat(AIRS_CC_RAD_GEOTRACK)
   real glintlon(AIRS_CC_RAD_GEOTRACK)
   integer*2 sun_glint_distance(AIRS_CC_RAD_GEOXTRACK,AIRS_CC_RAD_GEOTRACK)
   real*8 nadirTAI(AIRS_CC_RAD_GEOTRACK)
   real*8 sat_lat(AIRS_CC_RAD_GEOTRACK)
   real*8 sat_lon(AIRS_CC_RAD_GEOTRACK)
   byte scan_node_type(AIRS_CC_RAD_GEOTRACK)
   real topog(AIRS_CC_RAD_GEOXTRACK,AIRS_CC_RAD_GEOTRACK)
   real topog_err(AIRS_CC_RAD_GEOXTRACK,AIRS_CC_RAD_GEOTRACK)
   real landFrac(AIRS_CC_RAD_GEOXTRACK,AIRS_CC_RAD_GEOTRACK)
   real landFrac_err(AIRS_CC_RAD_GEOXTRACK,AIRS_CC_RAD_GEOTRACK)
   integer ftptgeoqa(AIRS_CC_RAD_GEOXTRACK,AIRS_CC_RAD_GEOTRACK)
   integer*2 zengeoqa(AIRS_CC_RAD_GEOXTRACK,AIRS_CC_RAD_GEOTRACK)
   integer*2 demgeoqa(AIRS_CC_RAD_GEOXTRACK,AIRS_CC_RAD_GEOTRACK)
   integer satgeoqa(AIRS_CC_RAD_GEOTRACK)
   integer*2 glintgeoqa(AIRS_CC_RAD_GEOTRACK)
   integer*2 moongeoqa(AIRS_CC_RAD_GEOTRACK)
   byte CalFlag(AIRS_CC_RAD_CHANNEL,AIRS_CC_RAD_GEOTRACK)
   byte CalScanSummary(AIRS_CC_RAD_GEOTRACK)
   byte CalChanSummary(AIRS_CC_RAD_CHANNEL)
   byte ExcludedChans(AIRS_CC_RAD_CHANNEL)
   real orbit_phase_deg(AIRS_CC_RAD_GEOTRACK)
   real shift_y0(AIRS_CC_RAD_MODULE,AIRS_CC_RAD_GEOTRACK)
   real scan_freq(AIRS_CC_RAD_CHANNEL,AIRS_CC_RAD_GEOTRACK)
   real Doppler_shift_ppm(AIRS_CC_RAD_GEOXTRACK,AIRS_CC_RAD_GEOTRACK)
   real NeN_L1B(AIRS_CC_RAD_CHANNEL)
   real NeN_L1B_Static(AIRS_CC_RAD_CHANNEL)
   integer*2 dust_flag(AIRS_CC_RAD_GEOXTRACK,AIRS_CC_RAD_GEOTRACK)
   real CC_noise_eff_amp_factor(AIRS_CC_RAD_GEOXTRACK,AIRS_CC_RAD_GEOTRACK)
   real CC1_noise_eff_amp_factor(AIRS_CC_RAD_GEOXTRACK,AIRS_CC_RAD_GEOTRACK)
   real CC1_Resid(AIRS_CC_RAD_GEOXTRACK,AIRS_CC_RAD_GEOTRACK)
   real CCfinal_Resid(AIRS_CC_RAD_GEOXTRACK,AIRS_CC_RAD_GEOTRACK)
   real TotCld_4_CCfinal(AIRS_CC_RAD_GEOXTRACK,AIRS_CC_RAD_GEOTRACK)
   real CCfinal_Noise_Amp(AIRS_CC_RAD_GEOXTRACK,AIRS_CC_RAD_GEOTRACK)
   byte invalid(AIRS_CC_RAD_GEOXTRACK,AIRS_CC_RAD_GEOTRACK)
   byte all_spots_avg(AIRS_CC_RAD_GEOXTRACK,AIRS_CC_RAD_GEOTRACK)
   byte MW_ret_used(AIRS_CC_RAD_GEOXTRACK,AIRS_CC_RAD_GEOTRACK)
   byte bad_clouds(AIRS_CC_RAD_GEOXTRACK,AIRS_CC_RAD_GEOTRACK)
   byte retrieval_type(AIRS_CC_RAD_GEOXTRACK,AIRS_CC_RAD_GEOTRACK)
end type



contains

!   Updated for V6 by Feng Ding in January 2013 

! This function is autogenerated by the mkezio program to read
! an AIRS swath of type "L2_Standard_cloud-cleared_radiance_product" from file given by the
! file_name argument into a buffer pointed to by the airs_cc_rad_gran
! argument.  The caller owns the buffer.  The entire granule
! is read -- every attribute and field, the whole lat/lon/time
! extent.
!
! Errors opening the file, etc. are fatal and cause STOP.
! Problems reading individual attributes or fields are reported to
! the console but do not interrupt program flow.

subroutine airs_cc_rad_rdr(file_name, airs_cc_rad_gran)
      character(len=*) file_name
      TYPE(airs_granule_type) airs_cc_rad_gran

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
        print *, "swinqswath found ", nswath, " swaths for file ",  &
                file_name, " Need exactly 1"
        stop
      end if

      ! There's exactly one swath.  Make sure it is the right one. 
      if (swathname .ne. 'L2_Standard_cloud-cleared_radiance_product') then 
        print *, "Error: bad swath name ", swathname, " in file ",  file_name
        print *, "Expected L2_Standard_cloud-cleared_radiance_product"
        stop
      end if

      ! Attach to (open) the one swath.
      swid = swattach(fid, swathname)
      if (swid .eq. -1) then 
        print *, "Failed to attach to swath ", swathname, &
                " in file ", file_name
        stop
      end if

! Attributes &
      statn = swrdattr(swid, "processing_level", airs_cc_rad_gran%processing_level) 
      if (statn .ne. 0)  &
       print *, "Error ", statn, " reading attribute ", "processing_level"
 
      statn = swrdattr(swid, "instrument", airs_cc_rad_gran%instrument) 
      if (statn .ne. 0)  &
       print *, "Error ", statn, " reading attribute ", "instrument"
 
      statn = swrdattr(swid, "DayNightFlag", airs_cc_rad_gran%DayNightFlag) 
      if (statn .ne. 0)  &
       print *, "Error ", statn, " reading attribute ", "DayNightFlag"
 
      statn = swrdattr(swid, "AutomaticQAFlag", airs_cc_rad_gran%AutomaticQAFlag) 
      if (statn .ne. 0)  &
       print *, "Error ", statn, " reading attribute ", "AutomaticQAFlag"
 
      statn = swrdattr(swid, "NumTotalData", airs_cc_rad_gran%NumTotalData) 
      if (statn .ne. 0)  &
       print *, "Error ", statn, " reading attribute ", "NumTotalData"
 
      statn = swrdattr(swid, "NumProcessData", airs_cc_rad_gran%NumProcessData) 
      if (statn .ne. 0)  &
       print *, "Error ", statn, " reading attribute ", "NumProcessData"

      statn = swrdattr(swid, "NumSpecialData",  airs_cc_rad_gran%NumSpecialData)
      if (statn .ne. 0) &
       print *, "Error ", statn, " reading attribute ", "NumSpecialData"

      statn = swrdattr(swid, "NumBadData", airs_cc_rad_gran%NumBadData)
      if (statn .ne. 0) &
       print *, "Error ", statn, " reading attribute ", "NumBadData"

      statn = swrdattr(swid, "NumMissingData", airs_cc_rad_gran%NumMissingData)
      if (statn .ne. 0) &
       print *, "Error ", statn, " reading attribute ", "NumMissingData"

      statn = swrdattr(swid, "NumLandSurface", airs_cc_rad_gran%NumLandSurface)
      if (statn .ne. 0) &
       print *, "Error ", statn, " reading attribute ", "NumLandSurface"

      statn = swrdattr(swid, "NumOceanSurface", airs_cc_rad_gran%NumOceanSurface)
      if (statn .ne. 0) &
       print *, "Error ", statn, " reading attribute ", "NumOceanSurface"

      statn = swrdattr(swid, "node_type", airs_cc_rad_gran%node_type)
      if (statn .ne. 0) &
       print *, "Error ", statn, " reading attribute ", "node_type"

      statn = swrdattr(swid, "start_year", airs_cc_rad_gran%start_year)
      if (statn .ne. 0) &
       print *, "Error ", statn, " reading attribute ", "start_year"

      statn = swrdattr(swid, "start_month", airs_cc_rad_gran%start_month)
      if (statn .ne. 0) &
       print *, "Error ", statn, " reading attribute ", "start_month"

      statn = swrdattr(swid, "start_day", airs_cc_rad_gran%start_day)
      if (statn .ne. 0) &
       print *, "Error ", statn, " reading attribute ", "start_day"

      statn = swrdattr(swid, "start_hour", airs_cc_rad_gran%start_hour)
      if (statn .ne. 0) &
       print *, "Error ", statn, " reading attribute ", "start_hour"

      statn = swrdattr(swid, "start_minute", airs_cc_rad_gran%start_minute)
      if (statn .ne. 0) &
       print *, "Error ", statn, " reading attribute ", "start_minute"

      statn = swrdattr(swid, "start_sec", airs_cc_rad_gran%start_sec)
      if (statn .ne. 0) &
       print *, "Error ", statn, " reading attribute ", "start_sec"

      statn = swrdattr(swid, "start_orbit", &
                  airs_cc_rad_gran%start_orbit)
      if (statn .ne. 0) &
       print *, "Error ", statn, " reading attribute ", &
                 "start_orbit"

      statn = swrdattr(swid, "end_orbit", &
                  airs_cc_rad_gran%end_orbit)
      if (statn .ne. 0) &
       print *, "Error ", statn, " reading attribute ", &
                 "end_orbit"

      statn = swrdattr(swid, "orbit_path", &
                  airs_cc_rad_gran%orbit_path)
      if (statn .ne. 0) &
       print *, "Error ", statn, " reading attribute ", &
                 "orbit_path"

      statn = swrdattr(swid, "start_orbit_row", &
                  airs_cc_rad_gran%start_orbit_row)
      if (statn .ne. 0) &
       print *, "Error ", statn, " reading attribute ", &
                 "start_orbit_row"

      statn = swrdattr(swid, "end_orbit_row", &
                  airs_cc_rad_gran%end_orbit_row)
      if (statn .ne. 0) &
       print *, "Error ", statn, " reading attribute ", &
                 "end_orbit_row"

      statn = swrdattr(swid, "granule_number", &
                  airs_cc_rad_gran%granule_number)
      if (statn .ne. 0) &
       print *, "Error ", statn, " reading attribute ", &
                 "granule_number"

      statn = swrdattr(swid, "num_scansets", &
                  airs_cc_rad_gran%num_scansets)
      if (statn .ne. 0) &
       print *, "Error ", statn, " reading attribute ", &
                 "num_scansets"

      statn = swrdattr(swid, "num_scanlines", &
                  airs_cc_rad_gran%num_scanlines)
      if (statn .ne. 0) &
       print *, "Error ", statn, " reading attribute ", &
                 "num_scanlines"

      statn = swrdattr(swid, "start_Latitude", &
                  airs_cc_rad_gran%start_Latitude)
      if (statn .ne. 0) &
       print *, "Error ", statn, " reading attribute ", &
                 "start_Latitude"

      statn = swrdattr(swid, "start_Longitude", &
                  airs_cc_rad_gran%start_Longitude)
      if (statn .ne. 0) &
       print *, "Error ", statn, " reading attribute ", &
                 "start_Longitude"

      statn = swrdattr(swid, "start_Time", &
                  airs_cc_rad_gran%start_Time)
      if (statn .ne. 0) &
       print *, "Error ", statn, " reading attribute ", &
                 "start_Time"

      statn = swrdattr(swid, "end_Latitude", &
                  airs_cc_rad_gran%end_Latitude)
      if (statn .ne. 0) &
       print *, "Error ", statn, " reading attribute ", &
                 "end_Latitude"

      statn = swrdattr(swid, "end_Longitude", &
                  airs_cc_rad_gran%end_Longitude)
      if (statn .ne. 0) &
       print *, "Error ", statn, " reading attribute ", &
                 "end_Longitude"

      statn = swrdattr(swid, "end_Time", &
                  airs_cc_rad_gran%end_Time)
      if (statn .ne. 0) &
       print *, "Error ", statn, " reading attribute ", &
                 "end_Time"

      statn = swrdattr(swid, "eq_x_longitude", &
                  airs_cc_rad_gran%eq_x_longitude)
      if (statn .ne. 0) &
       print *, "Error ", statn, " reading attribute ", &
                 "eq_x_longitude"

      statn = swrdattr(swid, "eq_x_tai", &
                  airs_cc_rad_gran%eq_x_tai)
      if (statn .ne. 0) &
       print *, "Error ", statn, " reading attribute ", &
                 "eq_x_tai"

      statn = swrdattr(swid, "LonGranuleCen", &
                  airs_cc_rad_gran%LonGranuleCen)
      if (statn .ne. 0) &
       print *, "Error ", statn, " reading attribute ", &
                 "LonGranuleCen"

      statn = swrdattr(swid, "LatGranuleCen", &
                  airs_cc_rad_gran%LatGranuleCen)
      if (statn .ne. 0) &
       print *, "Error ", statn, " reading attribute ", &
                 "LatGranuleCen"

      statn = swrdattr(swid, "LocTimeGranuleCen", &
                  airs_cc_rad_gran%LocTimeGranuleCen)
      if (statn .ne. 0) &
       print *, "Error ", statn, " reading attribute ", &
                 "LocTimeGranuleCen"

      statn = swrdattr(swid, "num_fpe", &
                  airs_cc_rad_gran%num_fpe)
      if (statn .ne. 0) &
       print *, "Error ", statn, " reading attribute ", &
                 "num_fpe"

      statn = swrdattr(swid, "orbitgeoqa", &
                  airs_cc_rad_gran%orbitgeoqa)
      if (statn .ne. 0) &
       print *, "Error ", statn, " reading attribute ", &
                 "orbitgeoqa"

      statn = swrdattr(swid, "num_satgeoqa", &
                  airs_cc_rad_gran%num_satgeoqa)
      if (statn .ne. 0) &
       print *, "Error ", statn, " reading attribute ", &
                 "num_satgeoqa"

      statn = swrdattr(swid, "num_glintgeoqa", &
                  airs_cc_rad_gran%num_glintgeoqa)
      if (statn .ne. 0) &
       print *, "Error ", statn, " reading attribute ", &
                 "num_glintgeoqa"

      statn = swrdattr(swid, "num_moongeoqa", &
                  airs_cc_rad_gran%num_moongeoqa)
      if (statn .ne. 0) &
       print *, "Error ", statn, " reading attribute ", &
                 "num_moongeoqa"

      statn = swrdattr(swid, "num_ftptgeoqa", &
                  airs_cc_rad_gran%num_ftptgeoqa)
      if (statn .ne. 0) &
       print *, "Error ", statn, " reading attribute ", &
                 "num_ftptgeoqa"

      statn = swrdattr(swid, "num_zengeoqa", &
                  airs_cc_rad_gran%num_zengeoqa)
      if (statn .ne. 0) &
       print *, "Error ", statn, " reading attribute ", &
                 "num_zengeoqa"

      statn = swrdattr(swid, "num_demgeoqa", &
                  airs_cc_rad_gran%num_demgeoqa)
      if (statn .ne. 0) &
       print *, "Error ", statn, " reading attribute ", &
                 "num_demgeoqa"

      statn = swrdattr(swid, "CalGranSummary", &
                  airs_cc_rad_gran%CalGranSummary)
      if (statn .ne. 0) &
       print *, "Error ", statn, " reading attribute ", &
                 "CalGranSummary"

      statn = swrdattr(swid, "DCR_scan", &
                  airs_cc_rad_gran%DCR_scan)
      if (statn .ne. 0) &
       print *, "Error ", statn, " reading attribute ", &
                 "DCR_scan"

      statn = swrdattr(swid, "granules_present_L1B", &
                        airs_cc_rad_gran%granules_present_L1B)
      if (statn .ne. 0) &
       print *, "Error ", statn, " reading attribute ", &
                 "granules_present_L1B"


! Geolocation fields
      edge(1) = AIRS_CC_RAD_GEOXTRACK
      edge(2) = AIRS_CC_RAD_GEOTRACK
      statn = swrdfld(swid, "Latitude", start, stride, edge, &
                     airs_cc_rad_gran%Latitude)
      if (statn .ne. 0)  &
       print *, "Error ", statn, " reading field Latitude"

      statn = swrdfld(swid, "Longitude", start, stride, edge, &
                     airs_cc_rad_gran%Longitude)
      if (statn .ne. 0) &
       print *, "Error ", statn, " reading field Longitude"

      statn = swrdfld(swid, "Time", start, stride, edge, &
                     airs_cc_rad_gran%Time)
      if (statn .ne. 0) &
       print *, "Error ", statn, " reading field Time"


! Data Fields
      edge(3) = 45
      edge(2) = 30
      edge(1) = 2378 
      statn = SWrdfld(swid, "radiances",  &
                  start, stride, edge, &
                  airs_cc_rad_gran%radiances) 
      if (statn .ne. 0)  &
       print *, "Error ", statn, " reading field ", "radiances"

      edge(3) = 45
      edge(2) = 30
      edge(1) = 2378 
      statn = SWrdfld(swid, "radiances_QC",  &
                  start, stride, edge, &
                  airs_cc_rad_gran%radiances_QC) 
      if (statn .ne. 0)  &
       print *, "Error ", statn, " reading field ", "radiances"

      edge(3) = 45
      edge(2) = 30
      edge(1) = 2378 
      statn = SWrdfld(swid, "radiance_err",  &
                  start, stride, edge, &
                  airs_cc_rad_gran%radiance_err) 
      if (statn .ne. 0) &
       print *, "Error ", statn, " reading field ",  "radiance_err"

      edge(4) = 45
      edge(3) = 30
      edge(2) = 3
      edge(1) = 3
      statn = SWrdfld(swid, "CldClearParam", &
                  start, stride, edge, &
                  airs_cc_rad_gran%CldClearParam)
      if (statn .ne. 0) &
       print *, "Error ", statn, " reading field ", "CldClearParam"

      edge(1) = 2378
      statn = SWrdfld(swid, "nominal_freq", &
                  start, stride, edge, &
                  airs_cc_rad_gran%nominal_freq)
      if (statn .ne. 0) &
       print *, "Error ", statn, " reading field ", "nominal_freq"

      edge(2) = 45
      edge(1) = 30
      statn = SWrdfld(swid, "scanang", &
                  start, stride, edge, &
                  airs_cc_rad_gran%scanang)
      if (statn .ne. 0) &
       print *, "Error ", statn, " reading field ", "scanang"

      edge(1) = 45
      statn = SWrdfld(swid, "satheight", &
                  start, stride, edge, &
                  airs_cc_rad_gran%satheight)
      if (statn .ne. 0) &
       print *, "Error ", statn, " reading field ", "satheight"

      edge(1) = 45
      statn = SWrdfld(swid, "satroll", &
                  start, stride, edge, &
                  airs_cc_rad_gran%satroll)
      if (statn .ne. 0) &
       print *, "Error ", statn, " reading field ", &
                 "satroll"

      edge(1) = 45
      statn = SWrdfld(swid, "satpitch", &
                  start, stride, edge, &
                  airs_cc_rad_gran%satpitch)
      if (statn .ne. 0) &
       print *, "Error ", statn, " reading field ", &
                 "satpitch"

      edge(1) = 45
      statn = SWrdfld(swid, "satyaw", &
                  start, stride, edge, &
                  airs_cc_rad_gran%satyaw)
      if (statn .ne. 0) &
       print *, "Error ", statn, " reading field ", &
                 "satyaw"

      edge(2) = 45
      edge(1) = 30
      statn = SWrdfld(swid, "satzen", &
                  start, stride, edge, &
                  airs_cc_rad_gran%satzen)
      if (statn .ne. 0) &
       print *, "Error ", statn, " reading field ", &
                 "satzen"

      edge(2) = 45
      edge(1) = 30
      statn = SWrdfld(swid, "satazi", &
                  start, stride, edge, &
                  airs_cc_rad_gran%satazi)
      if (statn .ne. 0) &
       print *, "Error ", statn, " reading field ", &
                 "satazi"

      edge(2) = 45
      edge(1) = 30
      statn = SWrdfld(swid, "solzen", &
                  start, stride, edge, &
                  airs_cc_rad_gran%solzen)
      if (statn .ne. 0) &
       print *, "Error ", statn, " reading field ", &
                 "solzen"

      edge(2) = 45
      edge(1) = 30
      statn = SWrdfld(swid, "solazi", &
                  start, stride, edge, &
                  airs_cc_rad_gran%solazi)
      if (statn .ne. 0) &
       print *, "Error ", statn, " reading field ", &
                 "solazi"

      edge(1) = 45
      statn = SWrdfld(swid, "glintlat", &
                  start, stride, edge, &
                  airs_cc_rad_gran%glintlat)
      if (statn .ne. 0) &
       print *, "Error ", statn, " reading field ", &
                 "glintlat"

      edge(1) = 45
      statn = SWrdfld(swid, "glintlon", &
                  start, stride, edge, &
                  airs_cc_rad_gran%glintlon)
      if (statn .ne. 0) &
       print *, "Error ", statn, " reading field ", &
                 "glintlon"

      edge(2) = 45
      edge(1) = 30
      statn = SWrdfld(swid, "sun_glint_distance", &
                  start, stride, edge, &
                  airs_cc_rad_gran%sun_glint_distance)
      if (statn .ne. 0) &
       print *, "Error ", statn, " reading field ", &
                 "sun_glint_distance"

      edge(1) = 45
      statn = SWrdfld(swid, "nadirTAI", &
                  start, stride, edge, &
                  airs_cc_rad_gran%nadirTAI)
      if (statn .ne. 0) &
       print *, "Error ", statn, " reading field ", &
                 "nadirTAI"

      edge(1) = 45
      statn = SWrdfld(swid, "sat_lat", &
                  start, stride, edge, &
                  airs_cc_rad_gran%sat_lat)
      if (statn .ne. 0) &
       print *, "Error ", statn, " reading field ", &
                 "sat_lat"

      edge(1) = 45
      statn = SWrdfld(swid, "sat_lon", &
                  start, stride, edge, &
                  airs_cc_rad_gran%sat_lon)
      if (statn .ne. 0) &
       print *, "Error ", statn, " reading field ", &
                 "sat_lon"

      edge(1) = 45
      statn = SWrdfld(swid, "scan_node_type", &
                  start, stride, edge, &
                  airs_cc_rad_gran%scan_node_type)
      if (statn .ne. 0) &
       print *, "Error ", statn, " reading field ", &
                 "scan_node_type"

      edge(2) = 45
      edge(1) = 30
      statn = SWrdfld(swid, "topog", &
                  start, stride, edge, &
                  airs_cc_rad_gran%topog)
      if (statn .ne. 0) &
       print *, "Error ", statn, " reading field ", &
                 "topog"

      edge(2) = 45
      edge(1) = 30
      statn = SWrdfld(swid, "topog_err", &
                  start, stride, edge, &
                  airs_cc_rad_gran%topog_err)
      if (statn .ne. 0) &
       print *, "Error ", statn, " reading field ", &
                 "topog_err"

      edge(2) = 45
      edge(1) = 30
      statn = SWrdfld(swid, "landFrac", &
                  start, stride, edge, &
                  airs_cc_rad_gran%landFrac)
      if (statn .ne. 0) &
       print *, "Error ", statn, " reading field ", &
                 "landFrac"

      edge(2) = 45
      edge(1) = 30
      statn = SWrdfld(swid, "landFrac_err", &
                  start, stride, edge, &
                  airs_cc_rad_gran%landFrac_err)
      if (statn .ne. 0) &
       print *, "Error ", statn, " reading field ", &
                 "landFrac_err"

      edge(2) = 45
      edge(1) = 30
      statn = SWrdfld(swid, "ftptgeoqa", &
                  start, stride, edge, &
                  airs_cc_rad_gran%ftptgeoqa)
      if (statn .ne. 0) &
       print *, "Error ", statn, " reading field ", &
                 "ftptgeoqa"

      edge(2) = 45
      edge(1) = 30
      statn = SWrdfld(swid, "zengeoqa", &
                  start, stride, edge, &
                  airs_cc_rad_gran%zengeoqa)
      if (statn .ne. 0) &
       print *, "Error ", statn, " reading field ", &
                 "zengeoqa"

      edge(2) = 45
      edge(1) = 30
      statn = SWrdfld(swid, "demgeoqa", &
                  start, stride, edge, &
                  airs_cc_rad_gran%demgeoqa)
      if (statn .ne. 0) &
       print *, "Error ", statn, " reading field ", &
                 "demgeoqa"

      edge(1) = 45
      statn = SWrdfld(swid, "satgeoqa", &
                  start, stride, edge, &
                  airs_cc_rad_gran%satgeoqa)
      if (statn .ne. 0) &
       print *, "Error ", statn, " reading field ", &
                 "satgeoqa"

      edge(1) = 45
      statn = SWrdfld(swid, "glintgeoqa", &
                  start, stride, edge, &
                  airs_cc_rad_gran%glintgeoqa)
      if (statn .ne. 0) &
       print *, "Error ", statn, " reading field ", &
                 "glintgeoqa"

      edge(1) = 45
      statn = SWrdfld(swid, "moongeoqa", &
                  start, stride, edge, &
                  airs_cc_rad_gran%moongeoqa)
      if (statn .ne. 0) &
       print *, "Error ", statn, " reading field ", &
                 "moongeoqa"

      edge(2) = 45
      edge(1) = 2378
      statn = SWrdfld(swid, "CalFlag", &
                  start, stride, edge, &
                  airs_cc_rad_gran%CalFlag)
      if (statn .ne. 0) &
       print *, "Error ", statn, " reading field ", &
                 "CalFlag"

      edge(1) = 45
      statn = SWrdfld(swid, "CalScanSummary", &
                  start, stride, edge, &
                  airs_cc_rad_gran%CalScanSummary)
      if (statn .ne. 0) &
       print *, "Error ", statn, " reading field ", &
                 "CalScanSummary"

      edge(1) = 2378
      statn = SWrdfld(swid, "CalChanSummary", &
                  start, stride, edge, &
                  airs_cc_rad_gran%CalChanSummary)
      if (statn .ne. 0) &
       print *, "Error ", statn, " reading field ", &
                 "CalChanSummary"

      edge(1) = 2378
      statn = SWrdfld(swid, "ExcludedChans", &
                  start, stride, edge, &
                  airs_cc_rad_gran%ExcludedChans)
      if (statn .ne. 0) &
       print *, "Error ", statn, " reading field ", &
                 "ExcludedChans"

      edge(1) = 45
      statn = SWrdfld(swid, "orbit_phase_deg", &
                  start, stride, edge, &
                  airs_cc_rad_gran%orbit_phase_deg)
      if (statn .ne. 0) &
       print *, "Error ", statn, " reading field ", &
                 "orbit_phase_deg"

      edge(2) = 45
      edge(1) = 17
      statn = SWrdfld(swid, "shift_y0", &
                  start, stride, edge, &
                  airs_cc_rad_gran%shift_y0)
      if (statn .ne. 0) &
       print *, "Error ", statn, " reading field ", &
                 "shift_y0"

      edge(2) = 45
      edge(1) = 2378
      statn = SWrdfld(swid, "scan_freq", &
                  start, stride, edge, &
                  airs_cc_rad_gran%scan_freq)
      if (statn .ne. 0) &
       print *, "Error ", statn, " reading field ", &
                 "scan_freq"

      edge(2) = 45
      edge(1) = 30
      statn = SWrdfld(swid, "Doppler_shift_ppm", &
                  start, stride, edge, &
                  airs_cc_rad_gran%Doppler_shift_ppm)
      if (statn .ne. 0) &
       print *, "Error ", statn, " reading field ", &
                 "Doppler_shift_ppm"

      edge(1) = 2378
      statn = SWrdfld(swid, "NeN_L1B", &
                  start, stride, edge, &
                  airs_cc_rad_gran%NeN_L1B)
      if (statn .ne. 0) &
       print *, "Error ", statn, " reading field ", &
                 "NeN_L1B"

      edge(1) = 2378
      statn = SWrdfld(swid, "NeN_L1B_Static", &
                  start, stride, edge, &
                  airs_cc_rad_gran%NeN_L1B_Static)
      if (statn .ne. 0) &
       print *, "Error ", statn, " reading field ", &
                 "NeN_L1B_Static"

      edge(2) = 45
      edge(1) = 30
      statn = SWrdfld(swid, "dust_flag", &
                  start, stride, edge, &
                  airs_cc_rad_gran%dust_flag)
      if (statn .ne. 0) &
       print *, "Error ", statn, " reading field ", &
                 "dust_flag"

      edge(2) = 45
      edge(1) = 30
      statn = SWrdfld(swid, "CC_noise_eff_amp_factor", &
                  start, stride, edge, &
                  airs_cc_rad_gran%CC_noise_eff_amp_factor)
      if (statn .ne. 0) &
       print *, "Error ", statn, " reading field ", &
                 "CC_noise_eff_amp_factor"

      edge(2) = 45
      edge(1) = 30
      statn = SWrdfld(swid, "CC1_noise_eff_amp_factor", &
                  start, stride, edge, &
                  airs_cc_rad_gran%CC1_noise_eff_amp_factor)
      if (statn .ne. 0) &
       print *, "Error ", statn, " reading field ", &
                 "CC1_noise_eff_amp_factor"

      edge(2) = 45
      edge(1) = 30
      statn = SWrdfld(swid, "CC1_Resid", &
                  start, stride, edge, &
                  airs_cc_rad_gran%CC1_Resid)
      if (statn .ne. 0) &
       print *, "Error ", statn, " reading field ", &
                 "CC1_Resid"

      edge(2) = 45
      edge(1) = 30
      statn = SWrdfld(swid, "CCfinal_Resid", &
                  start, stride, edge, &
                  airs_cc_rad_gran%CCfinal_Resid)
      if (statn .ne. 0) &
       print *, "Error ", statn, " reading field ", &
                 "CCfinal_Resid"

      edge(2) = 45
      edge(1) = 30
      statn = SWrdfld(swid, "TotCld_4_CCfinal", &
                  start, stride, edge, &
                  airs_cc_rad_gran%TotCld_4_CCfinal)
      if (statn .ne. 0) &
       print *, "Error ", statn, " reading field ", &
                 "TotCld_4_CCfinal"

      edge(2) = 45
      edge(1) = 30
      statn = SWrdfld(swid, "CCfinal_Noise_Amp", &
                  start, stride, edge, &
                  airs_cc_rad_gran%CCfinal_Noise_Amp)
      if (statn .ne. 0) &
       print *, "Error ", statn, " reading field ", &
                 "CCfinal_Noise_Amp"

      edge(2) = 45
      edge(1) = 30
      statn = SWrdfld(swid, "invalid", &
                  start, stride, edge, &
                  airs_cc_rad_gran%invalid)
      if (statn .ne. 0) &
       print *, "Error ", statn, " reading field ", &
                 "invalid"

      edge(2) = 45
      edge(1) = 30
      statn = SWrdfld(swid, "all_spots_avg", &
                  start, stride, edge, &
                  airs_cc_rad_gran%all_spots_avg)
      if (statn .ne. 0) &
       print *, "Error ", statn, " reading field ", &
                 "all_spots_avg"

      edge(2) = 45
      edge(1) = 30
      statn = SWrdfld(swid, "MW_ret_used", &
                  start, stride, edge, &
                  airs_cc_rad_gran%MW_ret_used)
      if (statn .ne. 0) &
       print *, "Error ", statn, " reading field ", &
                 "MW_ret_ued"

      edge(2) = 45
      edge(1) = 30
      statn = SWrdfld(swid, "bad_clouds", &
                  start, stride, edge, &
                  airs_cc_rad_gran%bad_clouds)
      if (statn .ne. 0) &
       print *, "Error ", statn, " reading field ", &
                 "bad_clouds"

      edge(2) = 45
      edge(1) = 30
      statn = SWrdfld(swid, "retrieval_type", &
                  start, stride, edge, &
                  airs_cc_rad_gran%retrieval_type)
      if (statn .ne. 0) &
       print *, "Error ", statn, " reading field ", &
                 "retrieval_type"

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
