#!/bin/csh
#
# DART software - Copyright UCAR. This open source software is provided
# by UCAR, "as is", without charge, subject to all terms of use at
# http://www.image.ucar.edu/DAReS/DART/DART_download
#
################################################################################
#
# madis_archive.csh
#
# core mechanics from Ryan Torn, SUNY Albany
# adapted by G. Romine to get a full obs set from MADIS archive
#
# calls the get_madis routine for input - EDIT for your acct info,
# give this script an analysis date in yyyymmddhh format as an argument
# An appropriate driver script would call this routine with 
# dates for each needed analysis window.
#
# observations are windowed by type as determined by the vars below
#
# wrf preprocessor is run to domain check obs, etc... with all madis
# obs processed through the auxiliary stream inputs. 
#
# Before running, make sure you have built all of the needed observation
# converters, configured as you want, and with obs errors set as you deem 
# appropriate. You must also create a template directory with input.nml 
# and wrfinput_d01 valid for the domain you plan to use
#
# The output should be a single obs_sequence file for each analysis time 
# placed in outputdir/datea
#
################################################################################

   set datea            = ${1}

   #  System specific commands
   setenv   REMOVE 'rm -rf'
   setenv   COPY 'cp -pfr'
   setenv   MOVE 'mv -f'
   setenv   LINK 'ln -fs'
   setenv   WGET /contrib/wget/bin/wget

   # Paths
   set BASE_DIR         = /path/to/new/project
   set OUTPUT_DIR       = ${BASE_DIR}/output
   set OBSPROC_DIR      = ${BASE_DIR}/obsproc
   set DART_DIR         = /path/to/DART
   set TEMPLATE_DIR     = ${BASE_DIR}/template

   # MADIS rawinsonde conversion parameters
   set sigwnd  = .false.
   set sigtmp  = .false.

   # MADIS cloud track winds parameters
   set include_ir  = "T"
   set include_vis = "T"
   set include_wv  = "T"

   # MADIS obs type windows - ***IN MINUTES***  EDIT as you see fit!!!!!
   set METAR_WINDOW    = 30
   set MARINE_WINDOW   = 30
   set ACARS_WINDOW    = 60
   set RAWIN_WINDOW    = 60
   set SATCL_WINDOW    = 60
   ################################################################################
   # Should not need to edit below                                                #
   ################################################################################

   if ( ! -d ${OUTPUT_DIR}/${datea} ) mkdir -p ${OUTPUT_DIR}/${datea}

   mkdir -p $OBSPROC_DIR
   cd $OBSPROC_DIR
   ${REMOVE} obs_seq.*
   ${COPY} ${TEMPLATE_DIR}/input.nml input.nml
   ${LINK} ${TEMPLATE_DIR}/wrfinput_d01 wrfinput_d01

   set yyyy     = `echo $datea | cut -b1-4`
   set mm       = `echo $datea | cut -b5-6`
   set dd       = `echo $datea | cut -b7-8`
   set hh       = `echo $datea | cut -b9-10`
   set dstring  = `echo $datea 0 -w | ${DART_DIR}/models/wrf/work/advance_time`
   set gdate    = `echo $datea 0 -g | ${DART_DIR}/models/wrf/work/advance_time`

   #  Copy data from MADIS servers
     foreach fhr ( -1 0 )

        set datef = `echo $datea $fhr | ${DART_DIR}/models/wrf/work/advance_time`
# ACARS
        ${BASE_DIR}/get_madis $datef acars acars_input.nc

        if ( -e acars_input.nc.gz ) then
           gunzip -f acars_input.nc.gz
           ${DART_DIR}/observations/MADIS/work/convert_madis_acars
           ${REMOVE} acars_input.nc
           if ( -e obs_seq.acars ) then
              set gdate1 = (`echo $datea -${ACARS_WINDOW}m -g | ${DART_DIR}/models/wrf/work/advance_time`)
              set gdate2 = (`echo $datea  ${ACARS_WINDOW}m -g | ${DART_DIR}/models/wrf/work/advance_time`)
              echo 'acars window ' $gdate1 $gdate2
              ${REMOVE} script.sed input.nml
              cat >! script.sed << EOF
              /filename_seq/c\
              filename_seq = 'obs_seq.acars',
              /filename_seq_list/c\
              filename_seq_list = '',
              /filename_out/c\
              filename_out = 'obs_seq.acars_output'
              /first_obs_days/c\
              first_obs_days = ${gdate1[1]},
              /first_obs_seconds/c\
              first_obs_seconds = ${gdate1[2]},
              /last_obs_days/c\
              last_obs_days = ${gdate2[1]},
              /last_obs_seconds/c\
              last_obs_seconds = ${gdate2[2]},
              /edit_copies/c\
              edit_copies = .false.,
EOF
              sed -f script.sed ${TEMPLATE_DIR}/input.nml >! input.nml
              ${DART_DIR}/models/wrf/work/obs_sequence_tool
              ${MOVE} obs_seq.acars_output obs_seq.acars
           endif

        endif
# METAR
        ${BASE_DIR}/get_madis $datef metar metar_input.nc

        if ( -e metar_input.nc.gz ) then
           gunzip -f metar_input.nc.gz
           ${DART_DIR}/observations/MADIS/work/convert_madis_metar
           ${REMOVE} metar_input.nc
           if ( -e obs_seq.metar ) then
              set gdate1 = (`echo $datea -${METAR_WINDOW}m -g | ${DART_DIR}/models/wrf/work/advance_time`)
              set gdate2 = (`echo $datea  ${METAR_WINDOW}m -g | ${DART_DIR}/models/wrf/work/advance_time`)
              echo 'metar window ' $gdate1 $gdate2
              ${REMOVE} script.sed input.nml
              cat >! script.sed << EOF
              /filename_seq/c\
              filename_seq = 'obs_seq.metar',
              /filename_seq_list/c\
              filename_seq_list = '',
              /filename_out/c\
              filename_out = 'obs_seq.metar_output'
              /first_obs_days/c\
              first_obs_days = ${gdate1[1]},
              /first_obs_seconds/c\
              first_obs_seconds = ${gdate1[2]},
              /last_obs_days/c\
              last_obs_days = ${gdate2[1]},
              /last_obs_seconds/c\
              last_obs_seconds = ${gdate2[2]},
              /edit_copies/c\
              edit_copies = .false.,
EOF
              sed -f script.sed ${TEMPLATE_DIR}/input.nml >! input.nml
              ${DART_DIR}/models/wrf/work/obs_sequence_tool
              ${MOVE} obs_seq.metar_output obs_seq.metar
           endif

        endif
# RAOB
        ${BASE_DIR}/get_madis $datef raob rawin_input.nc

        if ( -e rawin_input.nc.gz ) then
           gunzip -f rawin_input.nc.gz
           echo $sigwnd $sigtmp | ${DART_DIR}/observations/MADIS/work/convert_madis_rawin
           ${REMOVE} rawin_input.nc
           if ( -e obs_seq.rawin ) then
              set gdate1 = (`echo $datea -${RAWIN_WINDOW}m -g | ${DART_DIR}/models/wrf/work/advance_time`)
              set gdate2 = (`echo $datea  ${RAWIN_WINDOW}m -g | ${DART_DIR}/models/wrf/work/advance_time`)
              echo 'rawin window ' $gdate1 $gdate2
              ${REMOVE} script.sed input.nml
              cat >! script.sed << EOF
              /filename_seq/c\
              filename_seq = 'obs_seq.rawin',
              /filename_seq_list/c\
              filename_seq_list = '',
              /filename_out/c\
              filename_out = 'obs_seq.rawin_output'
              /first_obs_days/c\
              first_obs_days = ${gdate1[1]},
              /first_obs_seconds/c\
              first_obs_seconds = ${gdate1[2]},
              /last_obs_days/c\
              last_obs_days = ${gdate2[1]},
              /last_obs_seconds/c\
              last_obs_seconds = ${gdate2[2]},
              /edit_copies/c\
              edit_copies = .false.,
EOF
              sed -f script.sed ${TEMPLATE_DIR}/input.nml >! input.nml
              ${DART_DIR}/models/wrf/work/obs_sequence_tool
              ${MOVE} obs_seq.rawin_output obs_seq.rawin
           endif

        endif

# MARITIME
        ${BASE_DIR}/get_madis $datef maritime marine_input.nc

        if ( -e marine_input.nc.gz ) then
           gunzip -f marine_input.nc.gz
           ${DART_DIR}/observations/MADIS/work/convert_madis_marine
           ${REMOVE} marine_input.nc

           if ( -e obs_seq.marine ) then
              set gdate1 = (`echo $datea -${MARINE_WINDOW}m -g | ${DART_DIR}/models/wrf/work/advance_time`)
              set gdate2 = (`echo $datea  ${MARINE_WINDOW}m -g | ${DART_DIR}/models/wrf/work/advance_time`)
              echo 'marine window ' $gdate1 $gdate2
              ${REMOVE} script.sed input.nml
              cat >! script.sed << EOF
              /filename_seq/c\
              filename_seq = 'obs_seq.marine',
              /filename_seq_list/c\
              filename_seq_list = '',
              /filename_out/c\
              filename_out = 'obs_seq.marine_output'
              /first_obs_days/c\
              first_obs_days = ${gdate1[1]},
              /first_obs_seconds/c\
              first_obs_seconds = ${gdate1[2]},
              /last_obs_days/c\
              last_obs_days = ${gdate2[1]},
              /last_obs_seconds/c\
              last_obs_seconds = ${gdate2[2]},
              /edit_copies/c\
              edit_copies = .false.,
EOF
              sed -f script.sed ${TEMPLATE_DIR}/input.nml >! input.nml
              ${DART_DIR}/models/wrf/work/obs_sequence_tool
              ${MOVE} obs_seq.marine_output obs_seq.marine
           endif
        endif

# SAT CLD
        ${BASE_DIR}/get_madis $datef HDW satwnd_input.nc

        if ( -e satwnd_input.nc.gz ) then
           gunzip -f satwnd_input.nc.gz
           echo $include_ir $include_vis $include_wv | ${DART_DIR}/observations/MADIS/work/convert_madis_satwnd
           ${REMOVE} satwnd_input.nc
           if ( -e obs_seq.satwnd ) then
              set gdate1 = (`echo $datea -${SATCL_WINDOW}m -g | ${DART_DIR}/models/wrf/work/advance_time`)
              set gdate2 = (`echo $datea  ${SATCL_WINDOW}m -g | ${DART_DIR}/models/wrf/work/advance_time`)
              echo 'satwnd window ' $gdate1 $gdate2
              ${REMOVE} script.sed input.nml
              cat >! script.sed << EOF
              /filename_seq/c\
              filename_seq = 'obs_seq.satwnd',
              /filename_seq_list/c\
              filename_seq_list = '',
              /filename_out/c\
              filename_out = 'obs_seq.satwnd_output'
              /first_obs_days/c\
              first_obs_days = ${gdate1[1]},
              /first_obs_seconds/c\
              first_obs_seconds = ${gdate1[2]},
              /last_obs_days/c\
              last_obs_days = ${gdate2[1]},
              /last_obs_seconds/c\
              last_obs_seconds = ${gdate2[2]},
              /edit_copies/c\
              edit_copies = .false.,
EOF
              sed -f script.sed ${TEMPLATE_DIR}/input.nml >! input.nml
              ${DART_DIR}/models/wrf/work/obs_sequence_tool
              ${REMOVE} obs_seq.satwnd
              ${MOVE} obs_seq.satwnd_output obs_seq.old
           endif
        endif

     end

   # If you do not use the wrf/dart obs preprocessor, use the obs_sequence_tool
   # to paste the various obs_seq files together into a single obs_seq file.

   # run wrf/dart obs preprocessor here.  does superobs, trim by wrf grid location,
   # inflate ob errors along boundaries, etc.
   echo ${gdate[1]}, ${gdate[2]} | ${DART_DIR}/models/wrf/work/wrf_dart_obs_preprocess

   ${MOVE} obs_seq.new ${OUTPUT_DIR}/${datea}/obs_seq.out

exit 0


