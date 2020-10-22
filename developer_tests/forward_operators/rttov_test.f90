! DART software - Copyright UCAR. This open source software is provided
! by UCAR, "as is", without charge, subject to all terms of use at
! http://www.image.ucar.edu/DAReS/DART/DART_download

!> test of rttov interfaces

program rttov_test

use         types_mod, only : r8

use     utilities_mod, only : register_module, initialize_utilities, &
                              error_handler, E_ERR, E_MSG, &
                              finalize_utilities 

use      location_mod, only : location_type

use      obs_kind_mod, only : get_index_for_type_of_obs

! The obs_def_rttov_mod.f90 has a bunch of includes for the RTTOV modules.

use obs_def_rttov_mod, only : atmos_profile_type, &
                              trace_gas_profile_type, &
                              cloud_profile_type, &
                              aerosol_profile_type, &
                              rttov_sensor_type, &
                              visir_metadata_type, &
                              mw_metadata_type, &
                              get_rttov_sensor, &
                              sensor_runtime_setup, &
                              do_forward_model, &
                              read_sensor_db_file, &
                              atmos_profile_setup, &
                              trace_gas_profile_setup, &
                              aerosol_profile_setup, &
                              cloud_profile_setup

use       rttov_types, only : rttov_options, rttov_options_scatt

   implicit none

   integer, parameter :: ens_size  = 3
   integer, parameter :: numlevels = 54

   type(location_type)          :: location
   type(atmos_profile_type)     :: atmos
   type(trace_gas_profile_type) :: trace_gas
   type(cloud_profile_type)     :: clouds
   type(aerosol_profile_type)   :: aerosols

   integer, parameter :: test_ir_plat_id   = 9
   integer, parameter :: test_ir_sat_id    = 2
   integer, parameter :: test_ir_sensor_id = 11

   integer, parameter :: test_mw_plat_id   = 9
   integer, parameter :: test_mw_sat_id    = 2
   integer, parameter :: test_mw_sensor_id = 3

   integer :: test_ir_instrument(3)
   integer :: test_mw_instrument(3)
   integer :: sensor_instrument(3)

   character(len=512) :: rttov_sensor_db_file = 'rttov_test_me.csv'
   type(rttov_sensor_type), pointer :: sensor_ir
   type(rttov_sensor_type), pointer :: sensor_mw

   logical :: use_q2m              = .true.
   logical :: use_uv10m            = .true.
   logical :: use_wfetch           = .true.
   logical :: use_water_type       = .true.
   logical :: use_salinity         = .true.
   logical :: supply_foam_fraction = .true.
   logical :: use_sfc_snow_frac    = .true.

   logical :: ozone_data  = .false.
   logical :: co2_data    = .false.
   logical :: n2o_data    = .false.
   logical :: ch4_data    = .false.
   logical :: co_data     = .false.
   logical :: so2_data    = .false.

   logical :: add_aerosl  = .false.
   integer :: aerosl_type = 1  ! OPAC
   integer :: flavor

   logical :: cfrac_data          = .true.
   logical :: clw_data            = .true.
   integer :: clw_scheme          = 1
   logical :: rain_data           = .false.
   logical :: ciw_data            = .true.
   integer :: ice_scheme          = 1
   integer :: idg_scheme          = 1
   logical :: use_icede           = .false.
   logical :: snow_data           = .false.
   logical :: graupel_data        = .false.
   logical :: hail_data           = .false.
   logical :: w_data              = .false.
   logical :: htfrtc_simple_cloud = .false.

   logical :: first_lvl_is_sfc  = .false.
   logical :: mw_clear_sky_only = .false.
   logical :: do_lambertian     = .false.
   logical :: use_totalice      = .false.
   logical :: use_zeeman        = .false.

   real(r8) :: radiances(ens_size)
   integer  :: error_status(ens_size)

   type(visir_metadata_type), pointer :: visir_mds(:)
   type(mw_metadata_type),    pointer :: mw_mds(:)
   type(visir_metadata_type), pointer :: visir_md
   type(mw_metadata_type),    pointer :: mw_md

   ! Profiles of t, q, p, and clw from RTTOV
   ! Original location: rttov_test/profile-datasets/standard54lev_clw/{001 - 003}

   ! K 
   real(r8), parameter :: t_data(*) = (/ &
         178.1742, 166.1031, 200.9644, &
         188.2529, 176.2058, 213.3252, &
         205.1321, 196.4627, 224.7168, &
         221.334,  217.0163, 235.2652, &
         237.077,  237.2116, 245.0305, &
         252.3965, 253.1156, 253.9712, &
         261.5169, 264.7174, 261.5057, &
         267.6492, 272.1215, 265.6042, &
         269.7398, 275.4912, 262.2031, &
         263.9926, 270.8047, 252.5451, &
         257.2797, 262.9886, 243.4565, &
         251.061,  255.8046, 235.2337, &
         245.3975, 249.3036, 227.8755, &
         240.2603, 243.4317, 221.3266, &
         235.5384, 238.1507, 218.3437, &
         231.2017, 233.7377, 216.5764, &
         227.2464, 229.8147, 215.4607, &
         223.4022, 226.9817, 215.2631, &
         219.8331, 224.819,  215.2,    &
         216.4919, 223.0744, 215.2,    &
         212.2492, 221.3649, 215.2,    &
         207.1632, 219.7497, 215.2,    &
         202.3874, 218.143,  215.4279, &
         198.0815, 216.785,  216.0201, &
         195.2696, 215.7,    216.5855, &
         198.3492, 215.7,    217.1217, &
         204.5386, 215.7,    217.6306, &
         210.5458, 215.716,  218.1151, &
         216.794,  216.2984, 218.577,  &
         222.5243, 222.0199, 219.0153, &
         228.1523, 227.5903, 219.4315, &
         233.7824, 233.0306, 221.2409, &
         239.2519, 238.2408, 225.7942, &
         244.4527, 243.2519, 230.2292, &
         249.4698, 248.1434, 234.4881, &
         254.3639, 252.8513, 238.5656, &
         259.0342, 257.3262, 242.4896, &
         263.4607, 261.5477, 246.2457, &
         267.7589, 265.3643, 249.8067, &
         271.8353, 268.9933, 253.2246, &
         275.7137, 272.4373, 256.4401, &
         279.3804, 275.7036, 259.4958, &
         282.8149, 278.7564, 262.0787, &
         285.1086, 281.6015, 263.6424, &
         286.9082, 284.2282, 265.0843, &
         289.0168, 286.3178, 266.4233, &
         291.3232, 288.0253, 267.6449, &
         293.4089, 289.5693, 268.7503, &
         295.3068, 290.968,  269.7518, &
         296.9916, 292.2075, 270.6383, &
         298.4603, 293.288,  271.4111, &
         299.7152, 294.2112, 272.0715, &
         300.7582, 294.9785, 272.6202, &
         301.5907, 295.591,  273.0583/)

   ! moisture MR (g/kg)
   real(r8), parameter :: q_data(*) = (/ &
        0.000878,  0.000881,  0.000874, &
        0.001458,  0.001357,  0.001380, &
        0.002198,  0.001843,  0.001836, &
        0.002879,  0.002278,  0.002247, &
        0.003382,  0.002680,  0.002618, &
        0.003716,  0.003022,  0.002869, &
        0.003732,  0.003246,  0.003028, &
        0.003732,  0.003369,  0.003078, &
        0.003684,  0.003421,  0.003110, &
        0.003527,  0.003395,  0.003110, &
        0.003348,  0.003282,  0.003080, &
        0.003135,  0.003155,  0.003047, &
        0.002940,  0.003100,  0.003016, &
        0.002763,  0.003061,  0.002989, &
        0.002600,  0.003002,  0.002964, &
        0.002436,  0.002924,  0.002941, &
        0.002251,  0.002808,  0.002919, &
        0.002099,  0.002700,  0.002899, &
        0.001999,  0.002583,  0.002860, &
        0.001791,  0.002418,  0.002824, &
        0.001685,  0.002221,  0.002801, &
        0.001621,  0.002095,  0.002799, &
        0.001625,  0.002002,  0.002799, &
        0.001727,  0.001960,  0.002799, &
        0.001817,  0.001999,  0.002847, &
        0.001991,  0.002066,  0.002914, &
        0.002663,  0.002338,  0.002977, &
        0.003953,  0.003408,  0.003089, &
        0.006087,  0.005998,  0.003579, &
        0.016127,  0.017736,  0.005300, &
        0.037237,  0.051778,  0.011865, &
        0.084644,  0.120799,  0.022572, &
        0.165268,  0.201043,  0.035125, &
        0.282768,  0.291349,  0.059046, &
        0.447403,  0.400520,  0.102923, &
        0.672443,  0.567801,  0.169275, &
        0.955907,  0.756625,  0.282357, &
        1.293410,  0.963454,  0.399943, &
        1.782302,  1.245328,  0.517178, &
        2.229163,  1.674571,  0.678396, &
        2.620337,  2.237942,  0.856932, &
        3.659562,  2.921221,  1.112063, &
        4.967239,  3.601064,  1.342873, &
        6.751134,  4.607604,  1.536489, &
        8.590907,  5.597516,  1.714924, &
        9.945896,  6.583699,  1.873173, &
       10.905132,  7.529462,  2.016856, &
       11.769977,  8.382210,  2.148540, &
       12.920812,  9.302317,  2.300865, &
       13.997381, 10.129011,  2.435632, &
       14.932796, 10.847873,  2.553067, &
       15.729741, 11.460740,  2.653370, &
       16.390504, 11.969164,  2.736708, &
       16.916947, 12.374422,  2.803219/)

   ! hPa 
   real(r8), parameter :: p_data(*) = (/ &
         0.005,    0.005,    0.005,    &
         0.0131,   0.0131,   0.0131,   &
         0.0304,   0.0304,   0.0304,   &
         0.0644,   0.0644,   0.0644,   &
         0.1263,   0.1263,   0.1263,   &
         0.2324,   0.2324,   0.2324,   &
         0.4052,   0.4052,   0.4052,   &
         0.6749,   0.6749,   0.6749,   &
         1.0801,   1.0801,   1.0801,   &
         1.6691,   1.6691,   1.6691,   &
         2.5011,   2.5011,   2.5011,   &
         3.6462,   3.6462,   3.6462,   &
         5.1864,   5.1864,   5.1864,   &
         7.215,    7.215,    7.215,    &
         9.8368,   9.8368,   9.8368,   &
         13.1672,  13.1672,  13.1672,  &
         17.3308,  17.3308,  17.3308,  &
         22.4601,  22.4601,  22.4601,  &
         28.6937,  28.6937,  28.6937,  &
         36.1735,  36.1735,  36.1735,  &
         45.043,   45.043,   45.043,   &
         55.4433,  55.4433,  55.4433,  &
         67.5109,  67.5109,  67.5109,  &
         81.3744,  81.3744,  81.3744,  &
         97.1505,  97.1505,  97.1505,  &
         114.9415, 114.9415, 114.9415, &
         134.8318, 134.8318, 134.8318, &
         156.8846, 156.8846, 156.8846, &
         181.1394, 181.1394, 181.1394, &
         207.6092, 207.6092, 207.6092, &
         236.2784, 236.2784, 236.2784, &
         267.1012, 267.1012, 267.1012, &
         300.,     300.,     300.,      &
         334.8648, 334.8648, 334.8648, &
         371.5529, 371.5529, 371.5529, &
         409.8893, 409.8893, 409.8893, &
         449.6677, 449.6677, 449.6677, &
         490.6516, 490.6516, 490.6516, &
         532.5769, 532.5769, 532.5769, &
         575.1538, 575.1538, 575.1538, &
         618.0706, 618.0706, 618.0706, &
         660.9965, 660.9965, 660.9965, &
         703.5863, 703.5863, 703.5863, &
         745.4841, 745.4841, 745.4841, &
         786.3278, 786.3278, 786.3278, &
         825.7546, 825.7546, 825.7546, &
         863.4047, 863.4047, 863.4047, &
         898.9275, 898.9275, 898.9275, &
         931.9853, 931.9853, 931.9853, &
         962.2587, 962.2587, 962.2587, &
         989.451,  989.451,  989.451,  &
         1013.292, 1013.292, 1013.292, &
         1033.544, 1033.544, 1033.544, &
         1050.,    1050.,    1050./)

   ! cloud liquid water (kg/kg) 
   real(r8), parameter :: clw_prof_data(*) = (/ &
         0.00000E+00, 0.00000E+00, 0.00000E+00, &
         0.00000E+00, 0.00000E+00, 0.00000E+00, &
         0.00000E+00, 0.00000E+00, 0.00000E+00, &
         0.00000E+00, 0.00000E+00, 0.00000E+00, &
         0.00000E+00, 0.00000E+00, 0.00000E+00, &
         0.00000E+00, 0.00000E+00, 0.00000E+00, &
         0.00000E+00, 0.00000E+00, 0.00000E+00, &
         0.00000E+00, 0.00000E+00, 0.00000E+00, &
         0.00000E+00, 0.00000E+00, 0.00000E+00, &
         0.00000E+00, 0.00000E+00, 0.00000E+00, &
         0.00000E+00, 0.00000E+00, 0.00000E+00, &
         0.00000E+00, 0.00000E+00, 0.00000E+00, &
         0.00000E+00, 0.00000E+00, 0.00000E+00, &
         0.00000E+00, 0.00000E+00, 0.00000E+00, &
         0.00000E+00, 0.00000E+00, 0.00000E+00, &
         0.00000E+00, 0.00000E+00, 0.00000E+00, &
         0.00000E+00, 0.00000E+00, 0.00000E+00, &
         0.00000E+00, 0.00000E+00, 0.00000E+00, &
         0.00000E+00, 0.00000E+00, 0.00000E+00, &
         0.00000E+00, 0.00000E+00, 0.00000E+00, &
         6.02780E-05, 1.95948E-05, 8.94624E-05, &
         1.00000E-04, 1.00000E-04, 7.95915E-05, &
         0.00000E+00, 1.86917E-05, 0.00000E+00, &
         0.00000E+00, 0.00000E+00, 0.00000E+00, &
         0.00000E+00, 0.00000E+00, 0.00000E+00, &
         0.00000E+00, 0.00000E+00, 0.00000E+00, &
         0.00000E+00, 0.00000E+00, 0.00000E+00, &
         0.00000E+00, 0.00000E+00, 0.00000E+00, &
         0.00000E+00, 0.00000E+00, 7.54063E-05, &
         0.00000E+00, 0.00000E+00, 3.69476E-05, &
         0.00000E+00, 0.00000E+00, 0.00000E+00, &
         0.00000E+00, 0.00000E+00, 0.00000E+00, &
         0.00000E+00, 0.00000E+00, 0.00000E+00, &
         0.00000E+00, 0.00000E+00, 0.00000E+00, &
         0.00000E+00, 0.00000E+00, 0.00000E+00, &
         0.00000E+00, 0.00000E+00, 0.00000E+00, &
         0.00000E+00, 0.00000E+00, 0.00000E+00, &
         0.00000E+00, 0.00000E+00, 0.00000E+00, &
         0.00000E+00, 0.00000E+00, 1.77796E-06, &
         2.29149E-05, 2.98885E-05, 5.87432E-05, &
         8.08015E-05, 8.72882E-05, 1.00000E-04, &
         1.00000E-04, 1.41726E-04, 1.00000E-04, &
         1.00000E-04, 1.92606E-04, 1.00000E-04, &
         1.00000E-04, 1.19949E-04, 1.00000E-04, &
         1.00000E-04, 3.23937E-05, 1.00000E-04, &
         7.80533E-05, 0.00000E+00, 6.50497E-05, &
         3.96129E-05, 0.00000E+00, 3.01453E-05, &
         4.85138E-06, 0.00000E+00, 0.00000E+00, &
         0.00000E+00, 0.00000E+00, 0.00000E+00, &
         0.00000E+00, 0.00000E+00, 0.00000E+00, &
         0.00000E+00, 0.00000E+00, 0.00000E+00, &
         0.00000E+00, 0.00000E+00, 0.00000E+00, &
         0.00000E+00, 0.00000E+00, 0.00000E+00, &
         0.00000E+00, 0.00000E+00, 0.00000E+00 /)

   real(r8), parameter :: cfrac_prof_data(*) = (/ &
         1.0d0, 1.0d0, 1.0d0, &
         1.0d0, 1.0d0, 1.0d0, &
         1.0d0, 1.0d0, 1.0d0, &
         1.0d0, 1.0d0, 1.0d0, &
         1.0d0, 1.0d0, 1.0d0, &
         1.0d0, 1.0d0, 1.0d0, &
         1.0d0, 1.0d0, 1.0d0, &
         1.0d0, 1.0d0, 1.0d0, &
         1.0d0, 1.0d0, 1.0d0, &
         1.0d0, 1.0d0, 1.0d0, &
         1.0d0, 1.0d0, 1.0d0, &
         1.0d0, 1.0d0, 1.0d0, &
         1.0d0, 1.0d0, 1.0d0, &
         1.0d0, 1.0d0, 1.0d0, &
         1.0d0, 1.0d0, 1.0d0, &
         1.0d0, 1.0d0, 1.0d0, &
         1.0d0, 1.0d0, 1.0d0, &
         1.0d0, 1.0d0, 1.0d0, &
         1.0d0, 1.0d0, 1.0d0, &
         1.0d0, 1.0d0, 1.0d0, &
         1.0d0, 1.0d0, 1.0d0, &
         1.0d0, 1.0d0, 1.0d0, &
         1.0d0, 1.0d0, 1.0d0, &
         1.0d0, 1.0d0, 1.0d0, &
         1.0d0, 1.0d0, 1.0d0, &
         1.0d0, 1.0d0, 1.0d0, &
         1.0d0, 1.0d0, 1.0d0, &
         1.0d0, 1.0d0, 1.0d0, &
         1.0d0, 1.0d0, 1.0d0, &
         1.0d0, 1.0d0, 1.0d0, &
         1.0d0, 1.0d0, 1.0d0, &
         1.0d0, 1.0d0, 1.0d0, &
         1.0d0, 1.0d0, 1.0d0, &
         1.0d0, 1.0d0, 1.0d0, &
         1.0d0, 1.0d0, 1.0d0, &
         1.0d0, 1.0d0, 1.0d0, &
         1.0d0, 1.0d0, 1.0d0, &
         1.0d0, 1.0d0, 1.0d0, &
         1.0d0, 1.0d0, 1.0d0, &
         1.0d0, 1.0d0, 1.0d0, &
         1.0d0, 1.0d0, 1.0d0, &
         1.0d0, 1.0d0, 1.0d0, &
         1.0d0, 1.0d0, 1.0d0, &
         1.0d0, 1.0d0, 1.0d0, &
         1.0d0, 1.0d0, 1.0d0, &
         1.0d0, 1.0d0, 1.0d0, &
         1.0d0, 1.0d0, 1.0d0, &
         1.0d0, 1.0d0, 1.0d0, &
         1.0d0, 1.0d0, 1.0d0, &
         1.0d0, 1.0d0, 1.0d0, &
         1.0d0, 1.0d0, 1.0d0, &
         1.0d0, 1.0d0, 1.0d0, &
         1.0d0, 1.0d0, 1.0d0, &
         1.0d0, 1.0d0, 1.0d0/)

   integer, parameter :: surftype_data(*)  = (/ &
      1, 1, 1/)
   integer, parameter :: watertype_data(*) = (/ &
      0, 1, 1/)

   real(r8), parameter :: sfc_t_data(*) = (/ &
      300.2, 294.8, 272.2/)

   real(r8), parameter :: salinity_data(*) = (/ &
      35.0, 33.0, 37.0/)

   real(r8), parameter :: foam_fraction_data(*) = (/ &
      0.0, 0.3, 1.0/)

   real(r8), parameter :: snow_fraction_data(*) = (/ &
      0.0, 0.0, 0.0/)

   real(r8), parameter :: specularity_data(*) = (/ &
      0.0, 0.3, 1.0/)

   real(r8), parameter :: fastem_data(*) = (/ &
      3.0,  3.0,  3.0, &
      5.0,  5.0,  5.0, &
     15.0, 15.0, 15.0, &
      0.1,  0.1,  0.1, &
      0.3,  0.3,  0.3/)

   real(r8), parameter :: s2m_t_data(*) = (/ &
    298.7, 295.2,  272.2 /)

   ! 2m moisture MR (g/kg)
   real(r8), parameter :: s2m_q_data(*) = (/ &
    15.128594151295280, 11.453320846477883, 2.7345502236699266/)

   real(r8), parameter :: sfc_p_data(*) = (/ &
   1077.0,  1063.0,  1100.0 /)

   real(r8), parameter :: s10m_u_data(*) = (/ &
      5.0,  0.0,  -3.0 /)
   real(r8), parameter :: s10m_v_data(*) = (/ &
      2.0,  7.0,  -1.0 /)

   real(r8), parameter :: wfetch_data(*) = (/ &
      100000.0, 50000.0, 200000.0 /)

   real(r8), parameter :: zenangle_data(*) = (/ &
      45.0, 50.0, 45.0 /)

   real(r8), parameter :: azangle_data(*) = (/ &
      0.0, 185.0, 170.0 /)

   real(r8), parameter :: sunzenangle_data(*) = (/ &
      45.0, 185.0, 35.0/)

   real(r8), parameter :: sunazangle_data(*) = (/ &
      179.0, 179.0, 179.0/)

   real(r8), parameter :: lat_data(*) = (/ &
      15.0, 45.0, 45.0/)

   real(r8), parameter :: lon_data(*) = (/ &
      300.0, 270.0, 180.0/)

   real(r8), parameter :: elev_data(*) = (/ &
      0.0, 0.0, 0.0/)

   real(r8), parameter :: mag_field_data(*) = (/ &
      0.2, 0.5, 0.7/)

   real(r8), parameter :: cosbk_data(*) = (/ &
      0.0, 0.0, 0.5/)

   integer :: channel
   integer :: i, j, ind1d

   type(rttov_options),       pointer :: opts        => null() ! Options for RTTOV-DIRECT 
   type(rttov_options_scatt), pointer :: opts_scatt  => null() ! Options for RTTOV-SCATT

   ! version controlled file description for error handling, do not edit
   character(len=*), parameter :: source   = 'rttov_test.f90'
   character(len=*), parameter :: revision = ''
   character(len=*), parameter :: revdate  = ''

   character(len=512) :: string1, string2, string3

   call initialize_utilities('rttov_test')
   call register_module(source,revision,revdate) 

   allocate(visir_mds(3))
   allocate(mw_mds(3))

   ! setup: create a dummy rttov_sensor_db file
   open(unit=142, file=rttov_sensor_db_file)

   write(142,'(A)')'<DART Obs Type>,<Platform ID>,<Satellite ID>,<Sensor ID>,' // &
      '<Sensor type>,<Coef file name>,<Channels to use list (optional)>'
   write(142,'(A,I0,A,I0,A,I0,A)')'EOS_2_AMSUA_TB,',test_mw_plat_id,',',test_mw_sat_id,',',&
      test_mw_sensor_id,',mw,rtcoef_eos_2_amsua.dat'
   write(142,'(A,I0,A,I0,A,I0,A)')'EOS_2_AIRS_RADIANCE,',test_ir_plat_id,',',test_ir_sat_id,',',&
      test_ir_sensor_id,',ir,rtcoef_eos_2_airs.H5,168,900'

   close(142)

   ! test 1: reading of newly created sensor file
   call read_sensor_db_file(rttov_sensor_db_file)

   ! test 2.1: check that the IR instrument id matches what is expected
   test_ir_instrument(1) = test_ir_plat_id
   test_ir_instrument(2) = test_ir_sat_id
   test_ir_instrument(3) = test_ir_sensor_id

   sensor_ir => get_rttov_sensor(test_ir_instrument)

   sensor_instrument(1) = sensor_ir % platform_id
   sensor_instrument(2) = sensor_ir % satellite_id
   sensor_instrument(3) = sensor_ir % sensor_id

   if (any(sensor_instrument /= test_ir_instrument)) then
      write(string1,*) 'ERROR: returned sensor_instrument (IR) ',sensor_instrument,&
         ' differed from expected:',test_ir_instrument
      call error_handler(E_ERR,'rttov_test',string1,&
         source, revision, revdate)
   end if

   if (size(sensor_ir % channels) /= 2) then
      write(string1,*) 'ERROR: returned sensor_instrument (IR) ',sensor_instrument,&
         ' channel size ',size(sensor_ir % channels),'differed from expected:',2
      
      call error_handler(E_ERR,'rttov_test',string1,&
         source, revision, revdate)
   end if

   ! test 2.2: check that the MW instrument id matches what is expected
   test_mw_instrument(1) = test_mw_plat_id
   test_mw_instrument(2) = test_mw_sat_id
   test_mw_instrument(3) = test_mw_sensor_id

   sensor_mw => get_rttov_sensor(test_mw_instrument)

   sensor_instrument(1) = sensor_mw % platform_id
   sensor_instrument(2) = sensor_mw % satellite_id
   sensor_instrument(3) = sensor_mw % sensor_id

   if (any(sensor_instrument /= test_mw_instrument)) then
      write(string1,*) 'ERROR: returned sensor_instrument (MW) ',sensor_instrument,&
         ' differed from expected:',test_mw_instrument

      call error_handler(E_ERR,'rttov_test',string1,&
         source, revision, revdate)
   end if

   ! test 3: exercise the IR forward (direct) operator

   ! test 3.1: setup atmos, trace gas, aerosols, and clouds
   call atmos_profile_setup(atmos, ens_size, numlevels, use_q2m, & 
      use_uv10m, use_wfetch, use_water_type, use_salinity,       &
      supply_foam_fraction, use_sfc_snow_frac)

   call trace_gas_profile_setup(trace_gas, ens_size, numlevels,  &
      ozone_data, co2_data, n2o_data, ch4_data, co_data, so2_data)

   if (add_aerosl) then
      call aerosol_profile_setup(aerosols, ens_size, numlevels,  &
         aerosl_type)
   end if

   call cloud_profile_setup(clouds, ens_size, numlevels,     &
      cfrac_data, clw_data, clw_scheme, rain_data, ciw_data, &
      ice_scheme, use_icede, snow_data, graupel_data,        &
      hail_data, w_data, htfrtc_simple_cloud)   

   ! test 3.2: fill in the profiles for the IR values
   do j=1,numlevels
      do i=1,ens_size
         ind1d = (j-1)*ens_size+i
         atmos%temperature(i,j) = t_data(ind1d)
         atmos%moisture(i,j)    = q_data(ind1d)/1000.d0
         atmos%pressure(i,j)    = 100.d0*p_data(ind1d) ! hPa -> Pa
         clouds%clw(i,j)        = clw_prof_data(ind1d)
         clouds%cfrac(i,j)      = cfrac_prof_data(ind1d)
      end do
   end do

   do i=1,ens_size
      atmos%surftype(i)      = surftype_data(i)
      atmos%water_type(i)    = watertype_data(i)
      atmos%skin_temp(i)     = sfc_t_data(i)
      atmos%sfc_salinity(i)  = salinity_data(i)
      atmos%sfc_foam_frac(i) = foam_fraction_data(i)
      atmos%sfc_snow_frac(i) = snow_fraction_data(i)
      atmos%s2m_t(i)         = s2m_t_data(i)
      atmos%s2m_q(i)         = s2m_q_data(i)/1000.d0
      atmos%sfc_p(i)         = 100.d0*sfc_p_data(i)  ! hPa -> Pa
      atmos%s10m_u(i)        = s10m_u_data(i)
      atmos%s10m_v(i)        = s10m_v_data(i)
      atmos%wfetch(i)        = wfetch_data(i)
      atmos%sfc_elev(i)      = elev_data(i)
   end do

   do i=1,ens_size
      visir_mds(i)%sat_az      = azangle_data(i)
      visir_mds(i)%sat_ze      = zenangle_data(i)
      visir_mds(i)%sun_az      = sunazangle_data(i)
      visir_mds(i)%sun_ze      = sunzenangle_data(i)
      visir_mds(i)%platform_id = sensor_ir % platform_id
      visir_mds(i)%sat_id      = sensor_ir % satellite_id
      visir_mds(i)%sensor_id   = sensor_ir % sensor_id
      visir_mds(i)%specularity = specularity_data(i)
   end do

   ! test 3.3: specify IR options
   allocate(opts)

   ! these options have descriptions above, except for where defaults are set outside of user-control
   opts % interpolation % addinterp        = .true.      ! Allow interpolation of input profile
   opts % interpolation % reg_limit_extrap = .true.      ! Extrapolate beyond top of model intelligently 
   opts % interpolation % interp_mode      = 1    
   opts % interpolation % lgradp           = .false.     ! Do not allow TL/AD/K of user pressure levels
   opts % interpolation % spacetop         = .true.      ! Treat model top as space boundary (not recommended to be false)

   opts % config % do_checkinput    = .true.
   opts % config % apply_reg_limits = .true.
   opts % config % verbose          = .true.
   opts % config % fix_hgpl         = .false.

   opts % rt_all % switchrad              = .false.                 ! switch to BT for AD/K if microwave
   opts % rt_all % do_lambertian          = do_lambertian
   opts % rt_all % lambertian_fixed_angle = .true.
   opts % rt_all % rad_down_lin_tau       = .false.
   opts % rt_all % use_q2m                = use_q2m
   opts % rt_all % addrefrac              = .false.
   opts % rt_all % plane_parallel         = .false.

   ! test 3.4: setup the IR sensor runtime
   call sensor_runtime_setup(sensor_ir,             &
                             ens_size=ens_size,     &
                             nlevs=numlevels,       &
                             opts=opts)

   visir_md => visir_mds(1)
   mw_md    => null()
   channel = 168

   flavor = get_index_for_type_of_obs('MSG_4_SEVIRI_RADIANCE')

   ! test 3.5: call the IR forward operator
   call do_forward_model(ens_size, numlevels, flavor, location, &
      atmos, trace_gas, clouds, aerosols, sensor_ir, channel, &
      first_lvl_is_sfc, mw_clear_sky_only, clw_scheme,        &
      ice_scheme, idg_scheme, aerosl_type, do_lambertian,     &
      use_totalice, use_zeeman, radiances, error_status,      &
      visir_md=visir_md, mw_md=mw_md)

   write(string1,*) 'radiances:',radiances

   call error_handler(E_MSG,'rttov_test',string1,&
      source, revision, revdate)

   ! test 4: exercise the MW scatt forward operator

   ! test 4.1: fill in the profiles for the MW values
   do i=1,ens_size
      ind1d = 0
      mw_mds(i)%fastem_p1 = fastem_data(ind1d+i)
      ind1d = 3
      mw_mds(i)%fastem_p2 = fastem_data(ind1d+i)
      ind1d = 6
      mw_mds(i)%fastem_p3 = fastem_data(ind1d+i)
      ind1d = 9
      mw_mds(i)%fastem_p4 = fastem_data(ind1d+i)
      ind1d = 12
      mw_mds(i)%fastem_p5 = fastem_data(ind1d+i)

      mw_mds(i)%sat_az      = azangle_data(i)
      mw_mds(i)%sat_ze      = zenangle_data(i)
      mw_mds(i)%platform_id = sensor_mw % platform_id
      mw_mds(i)%sat_id      = sensor_mw % satellite_id
      mw_mds(i)%sensor_id   = sensor_mw % sensor_id

      mw_mds(i)%mag_field = mag_field_data(i)
      mw_mds(i)%cosbk     = cosbk_data(i)
   end do

   ! test 4.2: set the opts_scatt properties
   opts % rt_all % switchrad            = .true.                  ! switch to BT for AD/K if microwave
   opts % rt_mw % apply_band_correction = .true.
   opts % rt_mw % clw_data              = .false.
   opts % rt_mw % clw_calc_on_coef_lev  = .false.                 ! calculate CLW optical depths on input levels
   opts % rt_mw % clw_scheme            = clw_scheme
   opts % rt_mw % clw_cloud_top         = 322.d0 ! hPa
   opts % rt_mw % fastem_version        = 6
   opts % rt_mw % supply_foam_fraction  = supply_foam_fraction

   allocate(opts_scatt)
   opts_scatt % config % do_checkinput    = .true.
   opts_scatt % config % apply_reg_limits = .true.
   opts_scatt % config % verbose          = .true.
   opts_scatt % config % fix_hgpl         = .false.

   opts_scatt % rad_down_lin_tau      = .false.
   opts_scatt % apply_band_correction = .true.
   opts_scatt % interp_mode           = 1
   opts_scatt % lgradp                = .false.               ! Do not allow TL/AD/K of user pressure levels
   opts_scatt % reg_limit_extrap      = .true.                ! intelligently extend beyond the model top
   opts_scatt % fastem_version        = 6
   opts_scatt % supply_foam_fraction  = supply_foam_fraction
   opts_scatt % use_q2m               = use_q2m
   opts_scatt % lradiance             = .false.               ! Calculate Brightness Temperatures 
   opts_scatt % lusercfrac            = .false.               ! Have RTTOV-SCATT calculate the effective cloud fraction 
   opts_scatt % cc_threshold          = -1

   ! test 4.1: setup the mw sensor
   call sensor_runtime_setup(sensor_mw,             &
                             ens_size=ens_size,     &
                             nlevs=numlevels,       &
                             opts=opts,             &
                             opts_scatt=opts_scatt)

   visir_md => null()
   mw_md    => mw_mds(1)
   channel  = 1

   ! test 4.4: call the MW scatt forward operator
   call do_forward_model(ens_size, numlevels, flavor, location,       &
      atmos, trace_gas, clouds, aerosols, sensor_mw, channel, &
      first_lvl_is_sfc, mw_clear_sky_only, clw_scheme,        &
      ice_scheme, idg_scheme, aerosl_type, do_lambertian,     &
      use_totalice, use_zeeman, radiances, error_status,      &
      visir_md=visir_md, mw_md=mw_md)

   write(string1,*) 'TB:',radiances

   call error_handler(E_MSG,'rttov_test',string1,&
      source, revision, revdate)

   call finalize_utilities('rttov_test')

end program rttov_test
