#!/bin/csh
#
# DART software - Copyright UCAR. This open source software is provided
# by UCAR, "as is", without charge, subject to all terms of use at
# http://www.image.ucar.edu/DAReS/DART/DART_download

# diagnostics_obs.csh - shell script that computes observation
#                       specific diagnostics.
#
# $1 - analysis date
# $2 - parameter file
#
#     created Aug. 2009 Ryan Torn, U. Albany

  set datea     = ${1}
  set paramfile = ${2}
  source $paramfile

  cd $OBS_DIAG_DIR
  ${COPY} ${RUN_DIR}/input.nml input.nml
  set gdate  = (`echo $datea 0 -g | ${DART_DIR}/models/wrf/work/advance_time`)
  set yyyy2  = `echo $datea | cut -b1-4`
  set mm2    = `echo $datea | cut -b5-6`
  set dd2    = `echo $datea | cut -b7-8`
  set hh2    = `echo $datea | cut -b9-10`

  #  Determine appropriate dates for observation diagnostics
  @ nhours = $OBS_VERIF_DAYS * 24
  set datef = `echo $datea -${nhours} | ${DART_DIR}/models/wrf/work/advance_time`
  set yyyy1 = `echo $datef | cut -b1-4`
  set mm1   = `echo $datef | cut -b5-6`
  set dd1   = `echo $datef | cut -b7-8`
  set hh1   = `echo $datef | cut -b9-10`

  @ half_bin  = $ASSIM_INT_HOURS / 2
  set datefbs = `echo $datef -${half_bin} | ${DART_DIR}/models/wrf/work/advance_time`
  set fbs_yyyy1 = `echo $datefbs | cut -b1-4`
  set fbs_mm1   = `echo $datefbs | cut -b5-6`
  set fbs_dd1   = `echo $datefbs | cut -b7-8`
  set fbs_hh1   = `echo $datefbs | cut -b9-10`

  set datefbe = `echo $datef ${half_bin} | ${DART_DIR}/models/wrf/work/advance_time`
  set fbe_yyyy1 = `echo $datefbe | cut -b1-4`
  set fbe_mm1   = `echo $datefbe | cut -b5-6`
  set fbe_dd1   = `echo $datefbe | cut -b7-8`
  set fbe_hh1   = `echo $datefbe | cut -b9-10`

  set datelbe = `echo $datea ${half_bin} | ${DART_DIR}/models/wrf/work/advance_time`
  set lbe_yyyy1 = `echo $datelbe | cut -b1-4`
  set lbe_mm1   = `echo $datelbe | cut -b5-6`
  set lbe_dd1   = `echo $datelbe | cut -b7-8`
  set lbe_hh1   = `echo $datelbe | cut -b9-10`

  while ( $datef <= $datea )

     if ( -e ${OUTPUT_DIR}/${datef}/obs_seq.final )  ${LINK} ${OUTPUT_DIR}/${datef}/obs_seq.final obs_seq.final_${datef}
     set datef = `echo $datef $ASSIM_INT_HOURS | ${DART_DIR}/models/wrf/work/advance_time`

  end
  readlink -f obs_seq.final_* >! flist

  cat >! script.sed << EOF
  /obs_sequence_name/c\
  obs_sequence_name = '',
  /obs_sequence_list/c\
  obs_sequence_list = 'flist',
  /first_bin_center/c\
  first_bin_center =  ${yyyy1}, ${mm1}, ${dd1}, ${hh1}, 0, 0,
  /last_bin_center/c\
  last_bin_center  =  ${yyyy2}, ${mm2}, ${dd2}, ${hh2}, 0, 0,
  /filename_seq /c\
  filename_seq = 'obs_seq.final',
  /filename_seq_list/c\
  filename_seq_list = '',
  /filename_out/c\
  filename_out = 'obs_seq.final_reduced',
  /first_obs_days/c\
  first_obs_days = -1,
  /first_obs_seconds/c\
  first_obs_seconds = -1,
  /last_obs_days/c\
  last_obs_days = -1,
  /last_obs_seconds/c\
  last_obs_seconds = -1,
  /edit_copies/c\
  edit_copies        = .true.,
  /new_copy_index/c\
  new_copy_index     = 1, 2, 3, 4, 5,
  /first_bin_start/c\
  first_bin_start    = ${fbs_yyyy1}, ${fbs_mm1}, ${fbs_dd1}, ${fbs_hh1}, 0, 0,
  /first_bin_end/c\
  first_bin_end      = ${fbe_yyyy1}, ${fbe_mm1}, ${fbe_dd1}, ${fbe_hh1}, 0, 0,
  /last_bin_end/c\
  last_bin_end       = ${lbe_yyyy1}, ${lbe_mm1}, ${lbe_dd1}, ${lbe_hh1}, 0, 0,
EOF


  sed -f script.sed ${RUN_DIR}/input.nml >! input.nml

  # create the state-space diagnostic summary

  ${DART_DIR}/models/wrf/work/obs_diag || exit 1
  ${MOVE} obs_diag_output.nc ${OUTPUT_DIR}/${datea}/.
  ${MOVE} `ls -1 observation_locations.*.dat | tail -1` ${OUTPUT_DIR}/${datea}/observation_locations.dat

  # create a netCDF file with the original observation data (may not have some of the unusual metadata)

  ${DART_DIR}/models/wrf/work/obs_seq_to_netcdf
  ${MOVE} obs_epoch* ${OUTPUT_DIR}/${datea}/
  ${REMOVE} *.txt obs_seq.final_* flist observation_locations.*.dat

  # prune the obs_seq.final and store result  keeps first 5 copies? why not set num_output_obs = 0
  # is it the time subsetting that is of interest?

  ${LINK} ${OUTPUT_DIR}/${datea}/obs_seq.final .
  ${DART_DIR}/models/wrf/work/obs_sequence_tool
  ${MOVE} obs_seq.final_reduced ${OUTPUT_DIR}/${datea}/.
  ${REMOVE} obs_seq.final

  # process the mean analysis increment

  cd ${OUTPUT_DIR}/${datea}
  ${COPY} ${SHELL_SCRIPTS_DIR}/mean_increment.ncl .
  echo "ncl ${OUTPUT_DIR}/${datea}/mean_increment.ncl" >! nclrun.out
  chmod +x nclrun.out
  ./nclrun.out

  touch ${OUTPUT_DIR}/${datea}/obs_diags_done

exit 0

