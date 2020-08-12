#!/bin/csh -f
#
# DART software - Copyright UCAR. This open source software is provided
# by UCAR, "as is", without charge, subject to all terms of use at
# http://www.image.ucar.edu/DAReS/DART/DART_download
 
#===========================================================================
# Convert a bunch of Ameriflux tower files
#===========================================================================

cat << EOF >! baseinput.nml
&preprocess_nml
    input_obs_kind_mod_file = '../../../assimilation_code/modules/observations/DEFAULT_obs_kind_mod.F90',
   output_obs_kind_mod_file = '../../../assimilation_code/modules/observations/obs_kind_mod.f90',
     input_obs_def_mod_file = '../../../observations/forward_operators/DEFAULT_obs_def_mod.F90',
    output_obs_def_mod_file = '../../../observations/forward_operators/obs_def_mod.f90',
   input_files              = '../../../observations/forward_operators/obs_def_tower_mod.f90'
   /

&obs_kind_nml
   /

&location_nml
   /

&utilities_nml
   module_details = .false.,
   termlevel      = 2
   /

&obs_sequence_nml
   write_binary_obs_sequence = .false.
   /

EOF

#===========================================================================
# Bartlett Experimental Forest [New Hampshire]
# lon/lat = -71.288077,44.0646397
# station height = 272 m
# instrument height =  24.5 m
# Time Zone = Eastern = -5 UTC
# Deciduous Broad-leaf forest

foreach YEAR ( 2004 2005 2006 )

   cat << EOF >! stationdata.txt

&level4_to_obs_nml
   text_input_file = '../data/USBar${YEAR}_L4_h.txt',
   obs_out_file    = '/glade/p/image/Observations/FluxTower/obs_seq.USBar.${YEAR}',
   year            = ${YEAR},
   timezoneoffset  = -5,
   latitude        =  44.0646397,
   longitude       = -71.288077,
   elevation       = 272,
   flux_height     = 24.5,
   maxgoodqc       = 3,
   verbose         = .TRUE.
   /

EOF

   cat baseinput.nml stationdata.txt >! input.nml
   ./level4_to_obs || exit 1

end


#===========================================================================
# Barrow Arctic Science Consortium [Alaska]
# lon/lat = -156.62588,71.322525
# station height = 1 m
# instrument height =  1.9 m
# Time Zone = Alaskan = -9 UTC
# Open shrubland

foreach YEAR ( 1999      2001 )

   cat << EOF >! stationdata.txt

&level4_to_obs_nml
   text_input_file = '../data/USBrw${YEAR}_L4_h.txt',
   obs_out_file    = '/glade/p/image/Observations/FluxTower/obs_seq.USBrw.${YEAR}',
   year            = ${YEAR},
   timezoneoffset  = -9,
   latitude        =  71.322525,
   longitude       = -156.62588,
   elevation       = 1,
   flux_height     = 1.9,
   maxgoodqc       = 3,
   verbose         = .TRUE.
   /
EOF

   cat baseinput.nml stationdata.txt >! input.nml
   ./level4_to_obs || exit 2

end

#===========================================================================
# Harvard Forest [Massachusetts] Harvard Forest
# lon/lat = -72.171478,42.5377556
# station height = 340 m
# instrument height =  29 m
# Time Zone = Eastern = -5 UTC
# Mixed forest

foreach YEAR (           1992 1993 1994 1995 1996 1997 1998 1999 \
               2000 2001 2002 2003 2004 2005 2006                )

   cat << EOF >! stationdata.txt

&level4_to_obs_nml
   text_input_file = '../data/USHa1${YEAR}_L4_h.txt',
   obs_out_file    = '/glade/p/image/Observations/FluxTower/obs_seq.USHa1.${YEAR}',
   year            = ${YEAR},
   timezoneoffset  = -5,
   latitude        =  42.5377556,
   longitude       = -72.171478,
   elevation       = 340,
   flux_height     = 29,
   maxgoodqc       = 3,
   verbose         = .TRUE.
   /
EOF

   cat baseinput.nml stationdata.txt >! input.nml
   ./level4_to_obs || exit 3

end

#===========================================================================
# US-NR1 Niwot Ridge [Colorado] Niwot Ridge Mountain Research Station
# lon/lat = -105.5464,40.0328778
# station height = 3050 m
# instrument height =  21.5 m
# Time Zone = Mountain = -7 UTC
# Evergreen Needle-leaf forest

foreach YEAR (                                         1998 1999 \
               2000 2001 2002 2003 2004 2005 2006 2007           )

   cat << EOF >! stationdata.txt

&level4_to_obs_nml
   text_input_file = '../data/USNR1${YEAR}_L4_h.txt',
   obs_out_file    = '/glade/p/image/Observations/FluxTower/obs_seq.USNR1.${YEAR}',
   year            = ${YEAR},
   timezoneoffset  = -7,
   latitude        =  40.0328778,
   longitude       = -105.5464,
   elevation       = 3050,
   flux_height     = 21.5,
   maxgoodqc       = 3,
   verbose         = .TRUE.
   /
EOF

   cat baseinput.nml stationdata.txt >! input.nml
   ./level4_to_obs || exit 4

end

#===========================================================================
# US-SP3 Ordway Swisher Biological [Florida] (Donaldson)
# lon/lat = -82.163283, 29.7547667
# station height = 36 m
# instrument height =  24.3 m
# Time Zone = Eastern = -5 UTC
# Evergreen Broad-leaf forest

foreach YEAR ( 1999 2000 2001 2002 2003 2004 )

   cat << EOF >! stationdata.txt

&level4_to_obs_nml
   text_input_file = '../data/USSP3${YEAR}_L4_h.txt',
   obs_out_file    = '/glade/p/image/Observations/FluxTower/obs_seq.USSP3.${YEAR}',
   year            = ${YEAR},
   timezoneoffset  = -5,
   latitude        =  29.7547667,
   longitude       = -82.163283,
   elevation       = 36,
   flux_height     = 24.3,
   maxgoodqc       = 3,
   verbose         = .TRUE.
   /
EOF

   cat baseinput.nml stationdata.txt >! input.nml
   ./level4_to_obs || exit 5

end

#===========================================================================
# US-SRM Santa Rita Exp. Range [Arizona] (Santa Rita Mesquite Savanna)
# lon/lat = -110.86611, 31.82143
# station height = 1116 m
# instrument height =  7.82 m
# Time Zone = Mountain = -7 UTC
# Open Shrubland

foreach YEAR ( 2004 2005 2006 )

   cat << EOF >! stationdata.txt

&level4_to_obs_nml
   text_input_file = '../data/USSRM${YEAR}_L4_h.txt',
   obs_out_file    = '/glade/p/image/Observations/FluxTower/obs_seq.USSRM.${YEAR}',
   year            = ${YEAR},
   timezoneoffset  = -7,
   latitude        =  31.82143,
   longitude       = -110.86611,
   elevation       = 1116,
   flux_height     = 7.82,
   maxgoodqc       = 3,
   verbose         = .TRUE.
   /
EOF

   cat baseinput.nml stationdata.txt >! input.nml
   ./level4_to_obs || exit 6

end

#===========================================================================
# US-WCr Treehaven [Wisconsin] (Willow Creek)
# lon/lat = -90.07985917, 45.80592667
# station height = 515 m
# instrument height =  29.6 m
# Time Zone = Central = -6 UTC
# Deciduous Broad-leaf forest

foreach YEAR ( 1999 2000 2001 2002 2003 2004 2005 2006 )

   cat << EOF >! stationdata.txt

&level4_to_obs_nml
   text_input_file = '../data/USWCr${YEAR}_L4_h.txt',
   obs_out_file    = '/glade/p/image/Observations/FluxTower/obs_seq.USWcr.${YEAR}',
   year            = ${YEAR},
   timezoneoffset  = -6,
   latitude        =  45.80592667,
   longitude       = -90.07985917,
   elevation       = 515,
   flux_height     = 29.6,
   maxgoodqc       = 3,
   verbose         = .TRUE.
   /
EOF

   cat baseinput.nml stationdata.txt >! input.nml
   ./level4_to_obs || exit 7

end

#===========================================================================
# US-Wrc Wind River Crane Site [Washington] (Wind River Experimental Forest)
# lon/lat = -121.9519111, 45.82048889
# station height = 371 m
# instrument height =  70 m
# Time Zone = Pacific = -8 UTC
# Evergreen Needle-leaf forest

foreach YEAR ( 1998 1999 2000 2001 2002      2004 2005 2006 )

   cat << EOF >! stationdata.txt

&level4_to_obs_nml
   text_input_file = '../data/USWrc${YEAR}_L4_h.txt',
   obs_out_file    = '/glade/p/image/Observations/FluxTower/obs_seq.USWrc.${YEAR}',
   year            = ${YEAR},
   timezoneoffset  = -8,
   latitude        =  45.82048889,
   longitude       = -121.9519111,
   elevation       = 371,
   flux_height     = 70,
   maxgoodqc       = 3,
   verbose         = .TRUE.
   /
EOF

   cat baseinput.nml stationdata.txt >! input.nml
   ./level4_to_obs || exit 8

end

exit 0


