! This code is not protected by the DART copyright agreement.
! DART $Id$

! rcs keywords: &RCSfile& 
!               &Revision: 2820 & &Date: 2010-07-26 14:35:27 -0500 (Mon, 26 Jul 2010) &
! SCCS IDENTIFICATION:  %W% %G% 
!                       %U% %P% 
!
!     master list for the instruments identified by NAVDAS
!
      integer, parameter :: mx_nm_instru=99
      integer :: nc(mx_nm_instru)
      integer :: nsat(mx_nm_instru)   !  set to 1 for satellite platforms
      character(len=16) :: c_label(mx_nm_instru) ,c_rec(mx_nm_instru)
      integer :: ival(mx_nm_instru)
      integer :: m
!
      data(nc(m),c_rec(m),nsat(m),c_label(m),ival(m),m=1,mx_nm_instru) &
       /   1,   'SMX',        0,   'sfc land',   400  & ! Land observations (coastal, manual, automated)
       ,   2,   'MTR',        0,   'sfc metar',  400  & ! Land observations (metar)
       ,  10,   'SHX',        0,   'sfc ship',   300  & ! Surface obs from ships, fixed and mobile, drifting buoys
       ,  20,   'MAC',        0,   'milacr',     300  & ! Automated aircraft (Military ACARS)
       , 121,   'MAC',        0,   'milacr_asc', 100  & !	   Military ACARS ascending profile
       ,  21,   'MAC',        0,   'milacr_asd', 100  & !	   Military ACARS ascent (deprofiled)
       , 122,   'MAC',        0,   'milacr_des', 100  & !	   Military ACARS descending profile
       ,  22,   'MAC',        0,   'milacr_ded', 100  & !	   Military ACARS descent (deprofiled)
       ,  23,   'MAC',        0,   'milacr_lvl', 300  & !	   Military ACARS level flight
       ,  25,   'XRX',        0,   'man-airep',  700  & ! Manual AIREP (header XRXX)
       ,  26,   'XRX',        0,   'man-Yairep', 700  & ! Manual AIREP (header YRXX)
       ,  30,   'S0A',        0,   'airep',      300  & ! Aircraft data (aireps)
       , 131,   'S0A',        0,   'airep_asc',  100  & !	   AIREP ascending profile
       ,  31,   'S0A',        0,   'airep_asd',  100  & !	   AIREP ascent (deprofiled)
       , 132,   'S0A',        0,   'airep_des',  100  & !	   AIREP descending profile
       ,  32,   'S0A',        0,   'airep_ded',  100  & !	   AIREP descent (deprofiled)
       ,  33,   'S0A',        0,   'airep_lvl',  300  & !	   AIREP level flight
       ,  34,   'S0A',        0,   'rec-airep',  700  & ! Voice reports from USAF Reconnaissance aircraft
       ,  35,   'AMD',        0,   'amdar',      300  & ! Automated aircraft data (AMDAR)
       , 136,   'AMD',        0,   'amdar_asc',  100  & !	   AMDAR ascending profile
       ,  36,   'AMD',        0,   'amdar_asd',  100  & !	   AMDAR ascent (deprofiled)
       , 137,   'AMD',        0,   'amdar_des',  100  & !	   AMDAR descending profile
       ,  37,   'AMD',        0,   'amdar_ded',  100  & !	   AMDAR descent (deprofiled)
       ,  38,   'AMD',        0,   'amdar_lvl',  300  & !	   AMDAR level flight
       ,  40,   'ACR',        0,   'acars',      300  & ! Automated aircraft (ACARS)
       , 141,   'ACR',        0,   'acars_asc',  100  & !	   ACARS ascending profile
       ,  41,   'ACR',        0,   'acars_asd',  100  & !	   ACARS ascent (deprofiled)
       , 142,   'ACR',        0,   'acars_des',  100  & !	   ACARS descending profile
       ,  42,   'ACR',        0,   'acars_ded',  100  & !	   ACARS descent (deprofiled)
       ,  43,   'ACR',        0,   'acars_lvl',  300  & !	   ACARS level flight
       ,  45,   'MCR',        0,   'mdcrs',      300  & ! Automated aircraft (MDCRS)
       , 146,   'MCR',        0,   'mdcrs_asc',  100  & !	   MDCRS ascending profile
       ,  46,   'MCR',        0,   'mdcrs_asd',  100  & !	   MDCRS ascent (deprofiled)
       , 147,   'MCR',        0,   'mdcrs_des',  100  & !	   MDCRS descending profile
       ,  47,   'MCR',        0,   'mdcrs_ded',  100  & !	   MDCRS descent (deprofiled)
       ,  48,   'MCR',        0,   'mdcrs_lvl',  300  & !	   MDCRS level flight
       ,  50,   'TSX',        1,   'cld wnds1',  400  & ! Satellite-derived wind obs
       ,  51,   'TWX',        1,   'cld wnds2',  400  & ! Satellite-derived wind obs
       ,  63,   'GM5A',       1,   'METEO-5B',   400  & ! AFWA Satellite-derived wind obs
       ,  53,   'GM7A',       1,   'METEO-5A',   400  & ! AFWA Satellite-derived wind obs
       ,  54,   'GM8A',       1,   'METEO-7A',   400  & ! AFWA Satellite-derived wind obs
       ,  52,   'GMT2',       1,   'METEO-5C',   400  & ! SSEC Satellite-derived wind obs (backup)
       ,  59,   'GMT1',       1,   'METEO-7C',   400  & ! SSEC Satellite-derived wind obs (backup)
       ,  55,   'GOSW',       1,   'GOES-W',     400  & ! SSEC Satellite-derived wind obs
       ,  56,   'GOSE',       1,   'GOES-E',     400  & ! SSEC Satellite-derived wind obs
       ,  64,   'G10A',       1,   'GOES-10',    400  & ! AFWA Satellite-derived wind obs
       ,  65,   'G11A',       1,   'GOES-11',    400  & ! AFWA Satellite-derived wind obs
       ,  66,   'G12A',       1,   'GOES-12',    400  & ! AFWA Satellite-derived wind obs
       ,  57,   'GMSA',       1,   'GMSAFW',     400  & ! AFWA Satellite-derived wind obs
       ,  58,   'GMSC',       1,   'GMSC',       400  & ! SSEC Satellite-derived wind obs (backup)
       ,  67,   'G12R',       1,   'G12-RS',     400  & ! SSEC Satellite-derived wind obs (rapid-scan)
       ,  68,   'GMSX',       1,   'GMSX',       400  & ! SSEC Satellite-derived wind obs (backup)
       ,  60,   'ssmi_',      1,   'ssmi ff1',   500  & ! SSM/I wind speed, air-sea EDR
       ,  61,   'SS5',        1,   'ssmi ff2',   500  & ! SSM/I wind speed, air-sea EDR
       ,  70,   'scat_',      1,   'scat winds', 300  & ! scatterometer ocean surface winds
       ,  71,   'qscat',      1,   'qscat wind', 300  & ! QuikSCAT ocean surface winds
       ,  72,   'wsat' ,      1,   'WSat wind',  300  & ! WindSat ocean surface winds
       ,  73,   'ascat' ,     1,   'ascat wind', 300  & ! ASCAT ocean surface winds
       ,  80,   'MTER',       1,   'terra wind', 400  & ! MODIS winds--terra satellite
       ,  81,   'MAQU',       1,   'aqua wind',  400  & ! MODIS winds--aqua satellite
       ,  82,   'AN15',       1,   'n15 wind',   400  & ! AVHRR winds--n15 satellite
       ,  83,   'AN16',       1,   'n16 wind',   400  & ! AVHRR winds--n16 satellite
       ,  84,   'AN17',       1,   'n17 wind',   400  & ! AVHRR winds--n17 satellite
       ,  85,   'AN18',       1,   'n18 wind',   400  & ! AVHRR winds--n18 satellite
       ,  86,   'AMET',       1,   'metop wind', 400  & ! AVHRR winds--metop satellite
       ,  87,   'AN19',       1,   'n19 wind',   400  & ! AVHRR winds--n19 satellite
       ,  90,   'PAB' ,       0,   'Aus synth',  400  & ! Australian sea-level pres synthetic
       , 101,   'S0X',        0,   'raob',         0  & ! Rawinsondes (land,ship,drop,mobil)
       , 110,   'PIB',        0,   'pibal',        0  & ! Pilot balloons (land,ship,mobil)
       , 195,   'analytic' ,  0,   'analytic',     0  & ! synthetic obs derived from analytic conditions
       , 170,   'amsub_ret',  1,   'amsub qprf', 400  & ! AMSU-B retrieved moisture profile
       , 120,   'S0F',        1,   'tovs T',     400  & ! tovs retrieved by nesdis
       , 125,   'BOGUS',      1,   'ssmi prh',   400  & ! bogus ob constructed from 'ssmi TPPW'
       , 126,   'BOGUS',      1,   'WSat prh',   400  & ! bogus ob constructed from 'wnsat TPPW'
       , 180,   'GPS' ,       1,   'gps T',        0  & ! GPS temp *** (added 2/
       , 190,   'tc_syn',     0,   'TC synth',     0  & ! tropical cyclone synthetic observations
       , 191,   'GTO',        0,   'TC synth1',	   0  & ! tropical cyclone synthetic GTO obs
       , 192,   'CTO',        0,   'TC synth2',	   0  & ! tropical cyclone synthetic CTO obs
       , 210,   'atovs_',     1,   'atovs bT',   400  & ! ATOVS brightness temp
       , 220,   'rtovs_',     1,   'rtovs bT',   400  & ! RTOVS brightness temp
       , 230,   'ssmt_',      1,   'ssmt1 bT',   500  & ! SSM/T1 brightness temp
       , 240,   'ssmt2_',     1,   'ssmt2 bT',   500  & ! SSM/T2 brightness temp
       , 250,   'ssmi_',      1,   'ssmi TPPW',  400  & ! ssm/i total precipitable water
       , 251,   'wsatpw',     1,   'WSat TPPW',  400  & ! WindSat total precipitable water
       , 181,   'hirs3_',     1,   'hirs3 bT',   400  & ! HIRS3 brightness temp
       , 183,   'amsub_',     1,   'amsub bT',   400  & ! AMSUB brightness temp
       , 184,   'mhs_',       1,   'mhs bT',     400  & ! MHS   brightness temp
       , 185,   'ssmis_',     1,   'ssmis bT',   400  & ! SSMIS brightness temp
       , 187,   'aqua_',      1,   'aqua bT',    400  & ! AQUA(AIRS/AMSU) brightness temp
       , 188,   'iasi_',      1,   'iasi bT',    400  & ! IASI (METOP) brightness temp
       , 189,   'iasipc',     1,   'iasi PC',    400  & ! IASI principle component
       , 999,   'MSG',        1,   'missing',    999  & ! dummy entry to save space for a future addition
       , 999,   'MSG',        1,   'missing',    999  & ! dummy entry to save space for a future addition
       , 999,   'MSG',        1,   'missing',    999  & ! dummy entry to save space for a future addition
       , 999,   'MSG',        1,   'missing',    999  & ! dummy entry to save space for a future addition
       , 999,   'MSG',        1,   'missing',    999  & ! dummy entry to save space for a future addition
       , 999,   'MSG',        1,   'missing',    999  & ! dummy entry to save space for a future addition
       , 999,   'MSG',        1,   'missing',    999  & ! dummy entry to save space for a future addition
       , 999,   'MSG',        1,   'missing',    999  / ! dummy entry to save space for a future addition

! <next few lines under version control, do not edit>
! $URL$
! $Id$
! $Revision$
! $Date$
