#!/bin/csh -f
#
# DART note: this file started life as:
# /glade/p/cesmdata/cseg/collections/cesm1_1_1/models/ocn/pop2/input_templates/ocn.ecosys.tavg.csh

cat >! $CASEBUILD/pop2conf/ecosys.tavg.nml << EOF
tavg_freq_opt             = 'nday'   'nyear'
tavg_freq                 =  1       1
tavg_file_freq_opt        = 'nmonth' 'nyear'
tavg_file_freq            =  1       1
tavg_start_opt            = 'nstep'  'nstep'
tavg_start                =  0       0
tavg_fmt_in               = 'nc'     'nc'
tavg_fmt_out              = 'nc'     'nc'
ltavg_has_offset_date     = .false.  .false.
tavg_offset_years         =  1       1
tavg_offset_months        =  1       1
tavg_offset_days          =  2       2
ltavg_one_time_header     = .false.  .false.
tavg_stream_filestrings   = 'ecosys.nday1' 'ecosys.nyear1'
EOF

#------------------------------------------------------------------------------------
# For now, set streams manually. You must only set as many streams as are declared
#  in the tavg_nml section. For example, if there are three streams:
#  @ s1 = $my_stream
#  @ s2 = $s1 + 1
#  @ s3 = $s2 + 1
#------------------------------------------------------------------------------------

@ my_stream = $1
if ($my_stream < 1) then
   echo invalid my_stream number  ($my_stream)
   exit 5
endif

@ s1 = 1             # use the base-model stream 1
@ s2 = $my_stream    # use an ecosystem-defined stream
@ s3 = $s2 + 1       # use an ecosystem-defined stream

cat >! $CASEROOT/Buildconf/pop2conf/ecosys_tavg_contents << EOF
$s1  ECOSYS_ATM_PRESS
$s1  ECOSYS_IFRAC
$s1  ECOSYS_XKW
$s1  SCHMIDT_O2
$s1  SCHMIDT_CO2
$s1  IRON_FLUX
$s1  NOx_FLUX
$s1  NHy_FLUX
$s1  PH
$s1  O2SAT
$s1  STF_O2
$s1  CO2STAR
$s1  DCO2STAR
$s1  pCO2SURF
$s1  DpCO2
$s1  FG_CO2
$s1  ATM_CO2
$s1  FvPER_DIC
$s1  FvICE_DIC
$s1  FvPER_ALK
$s1  FvICE_ALK
$s1  PO4
$s1  NO3
$s1  SiO3
$s1  NH4
$s1  Fe
$s1  O2
$s1  O2_ZMIN
$s1  O2_ZMIN_DEPTH
$s1  O2_PRODUCTION
$s1  O2_CONSUMPTION
$s1  AOU
$s1  DIC
$s1  J_DIC
$s1  ALK
$s1  H2CO3
$s1  HCO3
$s1  CO3
$s1  pH_3D
$s1  co3_sat_calc
$s1  zsatcalc
$s1  co3_sat_arag
$s1  zsatarag
$s1  DOC
$s1  DOC_prod
$s1  DOC_remin
$s1  spC
$s1  spChl
$s1  spCaCO3
$s1  diatC
$s1  diatChl
$s1  zooC
$s1  spFe
$s1  diatSi
$s1  diatFe
$s1  diazC
$s1  diazChl
$s1  diazFe
$s1  DON
$s1  DOFe
$s1  DOP
$s1  graze_sp
$s1  graze_diat
$s1  graze_diaz
$s1  sp_agg
$s1  diat_agg
$s1  photoC_sp
$s1  CaCO3_form
$s1  photoC_diat
$s1  photoC_diaz
$s1  photoC_NO3_sp
$s1  photoC_NO3_diat
$s1  photoC_NO3_diaz
$s1  Fe_scavenge
$s1  Fe_scavenge_rate
$s1  diaz_Nfix
$s1  bSi_form
$s1  NITRIF
$s1  DENITRIF
$s1  POC_PROD
$s1  CaCO3_PROD
$s1  SiO2_PROD
$s1  P_iron_PROD
$s1  POC_FLUX_IN
$s1  CaCO3_FLUX_IN
$s1  SiO2_FLUX_IN
$s1  P_iron_FLUX_IN
$s1  PAR_avg
$s1  sp_Fe_lim
$s1  diat_Fe_lim
$s1  diaz_Fe_lim
$s1  sp_N_lim
$s1  diat_N_lim
$s1  sp_PO4_lim
$s1  diat_PO4_lim
$s1  diaz_P_lim
$s1  diat_SiO3_lim
$s1  sp_light_lim
$s1  diat_light_lim
$s1  diaz_light_lim
$s1  DON_prod
$s1  DOFe_prod
$s1  DOP_prod
$s1  sp_loss
$s1  diat_loss
$s1  zoo_loss
$s1  diaz_loss
$s1  Jint_100m_DIC
$s1  Jint_100m_NO3
$s1  Jint_100m_NH4
$s1  Jint_100m_PO4
$s1  Jint_100m_Fe
$s1  Jint_100m_SiO3
$s1  Jint_100m_ALK
$s1  Jint_100m_O2
$s1  Jint_100m_DOC
$s1  tend_zint_100m_DIC
$s1  tend_zint_100m_NO3
$s1  tend_zint_100m_NH4
$s1  tend_zint_100m_PO4
$s1  tend_zint_100m_Fe
$s1  tend_zint_100m_SiO3
$s1  tend_zint_100m_ALK
$s1  tend_zint_100m_O2
$s1  tend_zint_100m_DOC
$s2  photoC_sp_zint
$s2  CaCO3_form_zint
$s2  photoC_diaz_zint
$s2  photoC_diat_zint
$s1  photoC_NO3_sp_zint
$s1  photoC_NO3_diat_zint
$s1  photoC_NO3_diaz_zint
$s2  ECOSYS_IFRAC_2
$s2  ECOSYS_XKW_2
$s2  DpCO2_2
$s2  FG_CO2_2
$s2  STF_O2_2
$s2  spC_zint_100m
$s2  spCaCO3_zint_100m
$s2  diazC_zint_100m
$s2  diatC_zint_100m
$s2  zooC_zint_100m
$s2  spChl_SURF
$s2  diazChl_SURF
$s2  diatChl_SURF
$s3  J_NO3
$s3  J_NH4
$s3  J_PO4
$s3  J_Fe
$s3  J_SiO3
$s3  J_ALK
$s3  UE_O2
$s3  VN_O2
$s3  WT_O2
$s3  KPP_SRC_O2
$s3  DIA_IMPVF_O2
$s3  HDIFE_O2
$s3  HDIFN_O2
$s3  HDIFB_O2
$s3  UE_DOC
$s3  VN_DOC
$s3  WT_DOC
$s3  DIA_IMPVF_DOC
$s3  HDIFE_DOC
$s3  HDIFN_DOC
$s3  HDIFB_DOC
$s3  UE_DIC
$s3  VN_DIC
$s3  WT_DIC
$s3  KPP_SRC_DIC
$s3  DIA_IMPVF_DIC
$s3  HDIFE_DIC
$s3  HDIFN_DIC
$s3  HDIFB_DIC
$s3  UE_Fe
$s3  VN_Fe
$s3  WT_Fe
$s3  KPP_SRC_Fe
$s3  DIA_IMPVF_Fe
$s3  HDIFE_Fe
$s3  HDIFN_Fe
$s3  HDIFB_Fe
EOF

#1  dust_FLUX_IN
#1   DON_remin
#1   DOFe_remin
#1   DOP_remin
#1   photoFe_diaz
#1   photoFe_diat
#1   photoFe_sp
#1  Jint_PO4
#1  Jint_NO3
#1  Jint_SiO3
#1  Jint_NH4
#1  Jint_Fe
#1  Jint_O2
#1  Jint_DIC
#1  Jint_ALK
#1  Jint_DOC
#1  Jint_spC
#1  Jint_spChl
#1  Jint_spCaCO3
#1  Jint_diatC
#1  Jint_diatChl
#1  Jint_zooC

