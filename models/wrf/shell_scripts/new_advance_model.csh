#!/bin/csh
#
# DART software - Copyright UCAR. This open source software is provided
# by UCAR, "as is", without charge, subject to all terms of use at
# http://www.image.ucar.edu/DAReS/DART/DART_download
 
# Shell script to run the WRF model from DART input.
# where the model advance is executed as a separate process.
#
# This script performs the following:
# 1.  Creates a temporary directory to run a WRF realization (see options)
# 2.  Copies or links the necessary files into the temporary directory
# 3.  Converts DART state vectors to wrf input
# 4.  Updates LBCs (optionally draws perturbations from WRF-Var random covariances)
# 5.  Writes a WRF namelist from a template
# 6.  Runs WRF
# 7.  Checks for incomplete runs
# 8.  Converts wrf output to DART state vectors

# NOTES:
# 1.  This version executes da_wrfvar.exe in serial (no mpirun)
# 2.  If the ensemble mean assim_model_state_ic_mean is present in the 
# $CENTRALDIR, it is converted to a WRF netCDF format.
# It is then used in update_wrf_bc to calculate the deviation from the mean.
# This deviation from the mean is then added at the end of the interval to
# calculate new boundary tendencies. The magnitude of the perturbation added
# at the end of the interval is controlled by infl. The purpose is to increase
# time correlation at the lateral boundaries.

# POSSIBLE TODO
# 1.  modularization?  
# 2.  error checking all over

# this next section could be somewhere else in some docs.  If so
# reference is needed ("for more information about how to run this, see ...")

#-----snip--------------

#-------------------------------------------------------
# Dependencies (user responsibility)
#-------------------------------------------------------
# REQUIRED:
# 1. advance_time (from DART), located in your $CENTRALDIR
# 2. one of either (da_wrfvar.exe and pert_wrf_bc) or update_wrf_bc if you 
# want to run real-data cases with specified LBCs.  Elaborated below.
# 3. directory $CENTRALDIR/WRF_RUN containing all the WRF run-time files
# (typically files with data for the physics: LANDUSE.TBL, RRTM_DATA, etc
# but also anything else you want to link into the wrf-run directory.  If
# using WRF-Var then be.dat should be in there too.
# 4. wrf.exe, located in your $CENTRALDIR 
# 5. A wrfinput_d01 file in your $CENTRALDIR.  
# 6. namelist.input in your $CENTRALDIR for use as a template.  This file 
# should include the WRF-Var namelists if you are using WRF-Var (v3.1++ required).
#
# OPTIONAL:
# ####EITHER 1 or 2 is required for specified LBC runs
# 1.  da_wrfvar.exe (version 3.1 or later) and pert_wrf_bc in your $CENTRALDIR/WRF_RUN.
# In this case you also need be.dat (the be.dat.cv3 file from the WRF-Var 
# distribution) in your $CENTRALDIR/WRF_RUN, and WRF-Var namelists in 
# your $CENTRALDIR/namelist.input
# 2.  update_wrf_bc in your $CENTRALDIR for using pre-existing LBC files. Pre-existing LBC files should live in $CENTRALDIR/WRF

# File naming conventions:
#  mean wrfinput - wrfinput_d0X_${gday}_${gsec}_mean
#  mean wrfbdy - wrfbdy_d01_${gday}_${gsec}_mean
#  wrfbdy members - wrfbdy_d01_${gday}_${gsec}_${member}

#-----snip--------------

echo "new_advance_model.csh is running in `pwd`"

# Arguments are the process number of caller, the number of state copies
# belonging to that process, and the name of the filter_control_file for
# that process
set process = $1
set num_domains = $2
set control_file = $3
set num_states = 1      # forcing option of only one model advance per execution
set paramfile = $4      # Need this to load modules/environment
source $paramfile
# MULTIPLE DOMAINS - pass along the # of domains here?  We just default a value of 1 for the second variable, process is the ensemble member #


# Setting to vals > 0 saves wrfout files,
# will save all member output files <= to this value
set save_ensemble_member = 0
set delete_temp_dir = False

# set this to true if you want to maintain complete individual wrfinput/output
# for each member (to carry through non-updated fields)
set individual_members = true

# next line ensures that the last cycle leaves everything in the temp dirs
if ( $individual_members == true ) set delete_temp_dir = false

set  myname = $0
set  CENTRALDIR = `pwd`
echo $CENTRALDIR
set  WRFOUTDIR  = ${CENTRALDIR}/WRFOUT
set  REMOVE = '/bin/rm -rf'
set    COPY = '/bin/cp -p'
set    MOVE = '/bin/mv -f'
set      LN = '/bin/ln -sf'
unalias cd
unalias ls

setenv TARGET_CPU_LIST -1

# if process 0 go ahead and check for dependencies here
if ( $process == 0 ) then

   if ( ! -x ${CENTRALDIR}/advance_time ) then
     echo ABORT\: advance_model.csh could not find required executable dependency ${CENTRALDIR}/advance_time
     exit 1
   endif

   if ( ! -d WRF_RUN ) then
      echo ABORT\: advance_model.csh could not find required data directory ${CENTRALDIR}/WRF_RUN, which contains all the WRF run-time input files
      exit 1
   endif

   if ( ! -x ${CENTRALDIR}/WRF_RUN/da_wrfvar.exe ) then
      echo
      echo WARNING\: advance_model.csh could not find optional executable dependency ${CENTRALDIR}/WRF_RUN/da_wrfvar.exe
      echo
      if ( ! -x update_wrf_bc ) then
         # if the boundary conditions are specified, we need update_wrf_bc.  otherwise, it's ok if it isn't found.
         set SPEC_BC = `grep specified ${CENTRALDIR}/namelist.input | grep true | wc -l`
         if ( $SPEC_BC > 0 ) then
            echo ABORT\: advance_model.csh could not find required executable dependency ${CENTRALDIR}/update_wrf_bc
            exit 1
         endif
      endif

   else

     echo
     echo WARNING\: da_wrfvar.exe found, using it to update LBCs on the fly
     echo
     if ( ! -x ${CENTRALDIR}/pert_wrf_bc ) then
        echo ABORT\: advance_model.csh could not find required executable dependency ${CENTRALDIR}/pert_wrf_bc
        exit 1
     endif
     if ( ! -r ${CENTRALDIR}/WRF_RUN/be.dat ) then
        echo ABORT\: advance_model.csh could not find required readable dependency ${CENTRALDIR}/WRF_RUN/be.dat
        exit 1
     endif
     if ( ! -e ${CENTRALDIR}/bc_pert_scale ) then
        echo WARNING\:  using default VAR covariance scales
     endif

   endif

endif # process 0 dependency checking

# set this flag here for all processes so we don't have to keep checking
if ( -x ${CENTRALDIR}/WRF_RUN/da_wrfvar.exe ) then
   set USE_WRFVAR = 1
   echo "use wrfvar set"
else
   set USE_WRFVAR = 0
endif
 
# set this flag here if the radar additive noise script is found
if ( -e ${CENTRALDIR}/add_noise.csh ) then
   set USE_NOISE = 1
else
   set USE_NOISE = 0
endif
if ( -e ${CENTRALDIR}/replace_wrf_fields ) then
   set USE_REPLACE = 1
else
   set USE_REPLACE = 0
endif

# give the filesystem time to collect itself
sleep 5

# Each parallel task may need to advance more than one ensemble member.
# This control file has the actual ensemble number, the input filename,
# and the output filename for each advance.  Be prepared to loop and
# do the rest of the script more than once.
set USE_WRFVAR = 1
set state_copy = 1
set ensemble_member_line = 1
set input_file_line = 2    # MUTIPLE DOMAINS - this doesn't work for multiple domains, need something more general, maybe just a header name
set output_file_line = 3

# MULTIPLE DOMAINS - need a way to tell this shell script if there are multiple wrf domains in the analysis

while($state_copy <= $num_states)     # MULTIPLE DOMAINS - we don't expect advance model to run more than one member anymore. Reuse num_states for # domains?

   set ensemble_member = `head -n $ensemble_member_line ${CENTRALDIR}/${control_file} | tail -n 1`
   set input_file      = `head -n $input_file_line      ${CENTRALDIR}/${control_file} | tail -n 1`
   set output_file     = `head -n $output_file_line     ${CENTRALDIR}/${control_file} | tail -n 1`

   set infl = 0.0

   #  create a new temp directory for each member unless requested to keep and it exists already
   set temp_dir = "advance_temp${ensemble_member}"
   cd $temp_dir

#   if ( ( -d $temp_dir ) & ( $individual_members == "true" ) ) then
#
#      cd $temp_dir
##      set rmlist = ( `ls | grep -v wrfinput_d0.` )
#      set rmlist = ( `ls | egrep -v 'wrfinput_d0.|wrf.info'` )
#      echo $rmlist
#      ${REMOVE} $rmlist
#
#   else
#
#      ${REMOVE} $temp_dir >& /dev/null
#      mkdir -p $temp_dir
#      cd $temp_dir
#
#      ${COPY} ${CENTRALDIR}/wrfinput_d0? .
#
#   endif

   # link WRF-runtime files (required) and be.dat (if using WRF-Var)
   ${LN} ${CENTRALDIR}/WRF_RUN/*       .

   # Copy DART namelist if necessary
   if ( ! -e input.nml) then
      ${COPY} ${CENTRALDIR}/input.nml       .
   endif
   # append LSM data from previous cycle
   if ( -e ${CENTRALDIR}/append_lsm_data ) then
      ${LN} ${CENTRALDIR}/LSM/lsm_data_${ensemble_member}.nc lsm_data.nc
      ${CENTRALDIR}/append_lsm_data
      ${REMOVE} lsm_data.nc
   endif

   # nfile is required when using MPICH to run wrf.exe
   # nfile is machine specific.  Not needed on all platforms

   hostname >! nfile
   hostname >>! nfile

   #  Add number of domains information

# MULTIPLE_DOMAINS - need a more general instrument here
   if ( -e ${CENTRALDIR}/moving_domain_info ) then

      set MY_NUM_DOMAINS = `head -n 1 ${CENTRALDIR}/moving_domain_info | tail -n 1`
      ${MOVE} input.nml input.nml--
#      sed /num_domains/c\ "   num_domains = ${MY_NUM_DOMAINS}," input.nml-- >! input.nml
      cat >! script.sed << EOF
      /num_domains/c\
      num_domains = ${MY_NUM_DOMAINS},
EOF
      sed -f script.sed input.nml-- >! input.nml
      ${REMOVE} input.nml--

   endif

# DMODS - we don't have this option right now, and don't need to convert a file
#
#   # if a mean state ic file exists convert it to a wrfinput_mean netcdf file
#   if ( -e ${CENTRALDIR}/assim_model_state_ic_mean ) then
#      ${LN} ${CENTRALDIR}/assim_model_state_ic_mean dart_wrf_vector
#      ${CENTRALDIR}/dart_to_wrf >&! out.dart_to_wrf_mean
#      ${COPY} wrfinput_d01 wrfinput_mean
#      ${REMOVE} wrf.info dart_wrf_vector
#   endif
#
#   # ICs for this wrf run; Convert DART file to wrfinput netcdf file
#   ${MOVE} ${CENTRALDIR}/${input_file} dart_wrf_vector 
#   ${CENTRALDIR}/dart_to_wrf >&! out.dart_to_wrf
#   ${REMOVE} dart_wrf_vector

  set stuff_vars = $increment_vars_a

# may want multiple lists here, e.g. do we want w from the analysis?

  set stuff_str = ''   # these are variables we want to cycle
  set i = 1
  set num_vars = $#stuff_vars
  while ( $i < $num_vars )
     set stuff_str = `echo ${stuff_str}$stuff_vars[$i],`
     @ i ++
  end
  set stuff_str = `echo ${stuff_str}$stuff_vars[$num_vars]`
  echo "stuff var " ${stuff_str}

  set dn = 1
  while ( $dn <= $num_domains )

     set dchar = `echo $dn + 100 | bc | cut -b2-3` 
     set dchar = `printf %02d $dn` 
     set icnum = `echo $ensemble_member + 10000 | bc | cut -b2-5`
     set icnum = `printf %04d $ensemble_member`

     set this_file = filter_restart_d${dchar}.${icnum}

     if ( -e ../${this_file} ) then
        ncks -A -v ${stuff_str} ../${this_file} wrfinput_d${dchar}
     else
        echo "WARNING: ../${this_file} does not exist ..."
        echo "WARNING: this is expected for the first cycle ONLY!"
     endif

     @ dn ++  #
  end  
   #  Move and remove unnecessary domains    MULTIPLE DOMAINS - this problably needs to be removed to avoid confusion
   if ( -e ${CENTRALDIR}/moving_domain_info ) then

      set REMOVE_STRING = `cat ${CENTRALDIR}/remove_domain_info`
      if ( ${#REMOVE_STRING} > 0 )  ${REMOVE} ${REMOVE_STRING}

      set n = 1
      set NUMBER_FILE_MOVE = `cat ${CENTRALDIR}/rename_domain_info | wc -l`
      while ( $n <= $NUMBER_FILE_MOVE )
         ${MOVE} `head -n ${n} ${CENTRALDIR}/rename_domain_info | tail -n 1`
         @ n ++
      end

   endif

   # The program dart_to_wrf has created the file wrf.info.
   # Time information is extracted from wrf.info.
   # (bc in the following few lines is the calculator program,
   # not boundary conditions.)

# DMODS - note the wrf.info file was pre-generated, not from dart_to_wrf
   set secday = `head -n 1 wrf.info`
   set targsecs = $secday[1]
   set targdays = $secday[2]
   set targkey = `echo "$targdays * 86400 + $targsecs" | bc`

   set secday = `head -n 2 wrf.info | tail -n 1`
   set wrfsecs = $secday[1]
   set wrfdays = $secday[2]
   set wrfkey = `echo "$wrfdays * 86400 + $wrfsecs" | bc`

   echo "wrf.info is read" 
   echo $USE_WRFVAR
   # Find all BC's file available and sort them with "keys".
   # NOTE: this needs a fix for the idealized wrf case in which there are no
   # boundary files (also same for global wrf).  right now some of the
   # commands below give errors, which are ok to ignore in the idealized case
   # but it is not good form to generate spurious error messages.

   # check if LBCs are "specified" (in which case wrfbdy files are req'd)
   # and we need to set up a key list to manage target times
   set SPEC_BC = `grep specified ${CENTRALDIR}/namelist.input | grep true | wc -l`

   if ( $SPEC_BC > 0 ) then

      if ( $USE_WRFVAR ) then
         set bdyfiles = `ls ${CENTRALDIR}/WRF/wrfbdy_d01_*_mean`
      else
         set bdyfiles = `ls ${CENTRALDIR}/WRF/wrfbdy_d01_*_${ensemble_member} | grep -v mean`
      endif
      echo $bdyfiles
      set keylist = ()
      foreach f ( $bdyfiles )
	 set day = `echo $f | awk -F_ '{print $(NF-2)}'`
	 set sec = `echo $f | awk -F_ '{print $(NF-1)}'`
         set key = `echo "$day * 86400 + $sec" | bc`
         set keylist = ( $keylist $key )
      end

      set keys = `echo $keylist | sort`
   else  #  idealized WRF with non-specified BCs

      set keys = ( $targkey )

   endif

   set cal_date    = `head -n 3 wrf.info | tail -n 1`
   set START_YEAR  = $cal_date[1]
   set START_MONTH = $cal_date[2]
   set START_DAY   = $cal_date[3]
   set START_HOUR  = $cal_date[4]
   set START_MIN   = $cal_date[5]
   set START_SEC   = $cal_date[6]

   set START_STRING = ${START_YEAR}-${START_MONTH}-${START_DAY}_${START_HOUR}:${START_MIN}:${START_SEC}
   set datea = ${cal_date[1]}${cal_date[2]}${cal_date[3]}${cal_date[4]}

   set MY_NUM_DOMAINS    = `head -n 4 wrf.info | tail -n 1`
   set ADV_MOD_COMMAND   = `head -n 5 wrf.info | tail -n 1`

   #  Code for dealing with TC nests
   if ( -e ${CENTRALDIR}/fixed_domain_info ) then

      set MY_NUM_DOMAINS = `head -n 1 ${CENTRALDIR}/fixed_domain_info | tail -n 1`

   else if ( -e ${CENTRALDIR}/moving_domain_info ) then

      set MY_NUM_DOMAINS = `head -n 2 ${CENTRALDIR}/moving_domain_info | tail -n 1`
      ${MOVE} input.nml input.nml--
      cat >! script.sed << EOF
      /num_domains/c\
      num_domains = ${MY_NUM_DOMAINS},
EOF

#      sed /num_domains/c\ "   num_domains = ${MY_NUM_DOMAINS}," input.nml-- >! input.nml
      sed -f script.sed input.nml-- >! input.nml
      ${REMOVE} input.nml--

   endif 

   # Find the next BC's file available.

   @ ifile = 1
   while ( $keys[${ifile}] <= $wrfkey )
      if ( $ifile < $#bdyfiles ) then
         @ ifile ++
      else
         echo No boundary file available to move beyond
         echo $START_STRING
         exit 1
      endif
   end

   # radar additive noise option.  if shell script is available
   # in the centraldir, it will be called here.
   if ( $USE_NOISE ) then
      ${CENTRALDIR}/add_noise.csh $wrfsecs $wrfdays $state_copy $ensemble_member $temp_dir $CENTRALDIR
   endif

## run the replace_wrf_fields utility to update the static fields
   if ( $USE_REPLACE ) then
      echo ../wrfinput_d01 wrfinput_d01 | ${CENTRALDIR}/replace_wrf_fields
   endif


   ###############################################################
   # Advance the model with new BC until target time is reached. #
   ###############################################################

   while ( $wrfkey < $targkey )

      set iday = `echo "$keys[$ifile] / 86400" | bc`
      set isec = `echo "$keys[$ifile] % 86400" | bc`

      # Copy the boundary condition file to the temp directory if needed.
      if ( $SPEC_BC > 0 ) then

         if ( $USE_WRFVAR ) then
            ${COPY} ${CENTRALDIR}/WRF/wrfbdy_d01_${iday}_${isec}_mean               wrfbdy_d01
         else
            ${COPY} ${CENTRALDIR}/WRF/wrfbdy_d01_${iday}_${isec}_${ensemble_member} wrfbdy_d01
         endif

      endif

      if ( $targkey > $keys[$ifile] ) then
         set INTERVAL_SS = `echo "$keys[$ifile] - $wrfkey" | bc`
      else
         set INTERVAL_SS = `echo "$targkey - $wrfkey" | bc`
      endif

      set INTERVAL_MIN = `expr $INTERVAL_SS \/ 60`
#GSR      set INTERVAL_MIN = 60

      set END_STRING = `echo ${START_STRING} ${INTERVAL_SS}s -w | ${CENTRALDIR}/advance_time`
      set END_YEAR  = `echo $END_STRING | cut -c1-4`
      set END_MONTH = `echo $END_STRING | cut -c6-7`
      set END_DAY   = `echo $END_STRING | cut -c9-10`
      set END_HOUR  = `echo $END_STRING | cut -c12-13`
      set END_MIN   = `echo $END_STRING | cut -c15-16`
      set END_SEC   = `echo $END_STRING | cut -c18-19`
#      set datef = ${END_STRING[1]}${END_STRING[2]}${END_STRING[3]}${END_STRING[4]}

      # Update boundary conditions.
      # WARNING: da_wrfvar.exe will only work correctly if running WRF V3.1 or later!
      # If it is found in the central dir, use it to regnerate perturbed boundary files
      # Otherwise, do the original call to update_wrf_bc
      if ( $USE_WRFVAR ) then

         #  Set the covariance perturbation scales using file or default values
         if ( -e ${CENTRALDIR}/bc_pert_scale ) then
            set pscale = `head -n 1 ${CENTRALDIR}/bc_pert_scale | tail -n 1`
            set hscale = `head -n 2 ${CENTRALDIR}/bc_pert_scale | tail -n 1`
            set vscale = `head -n 3 ${CENTRALDIR}/bc_pert_scale | tail -n 1`
         else
            set pscale = 0.25
            set hscale = 1.0
            set vscale = 1.5
         endif
         @ iseed2 = $ensemble_member * 10

         ${REMOVE} script.sed
         cat >! script.sed << EOF
            /analysis_date/c\
            analysis_date = \'${END_STRING}.0000\',
            /as1/c\
            as1 = ${pscale}, ${hscale}, ${vscale},
            /as2/c\
            as2 = ${pscale}, ${hscale}, ${vscale},
            /as3/c\
            as3 = ${pscale}, ${hscale}, ${vscale},
            /as4/c\
            as4 = ${pscale}, ${hscale}, ${vscale},
            /as5/c\
            as5 = ${pscale}, ${hscale}, ${vscale},
            /seed_array1/c\
            seed_array1 = 1${END_MONTH}${END_DAY}${END_HOUR},
            /seed_array2/c\
            seed_array2 = $iseed2,
            /start_year/c\
            start_year = ${END_YEAR},
            /start_month/c\
            start_month = ${END_MONTH},
            /start_day/c\
            start_day = ${END_DAY},
            /start_hour/c\
            start_hour = ${END_HOUR},
            /start_minute/c\
            start_minute = ${END_MIN},
            /start_second/c\
            start_second = ${END_SEC},
            /end_year/c\
            end_year = ${END_YEAR},
            /end_month/c\
            end_month = ${END_MONTH},
            /end_day/c\
            end_day = ${END_DAY},
            /end_hour/c\
            end_hour = ${END_HOUR},
            /end_minute/c\
            end_minute = ${END_MIN},
            /end_second/c\
            end_second = ${END_SEC},
            /max_dom/c\
            max_dom = 1,
EOF
# The EOF on the line above MUST REMAIN in column 1.

         sed -f script.sed ${CENTRALDIR}/namelist.input >! namelist.input

         ${LN} ${CENTRALDIR}/WRF/wrfinput_d01_${targdays}_${targsecs}_mean ./fg
################################
## instead of running wrfda, just add static pertubations from the pert bank
#  note the static perturbation path is defined in the ncl script
#         setenv TARGET_CPU_RANGE        "-1"	
#         mpiexec_mpt dplace -s 1  ${CENTRALDIR}/WRF_RUN/da_wrfvar.exe >>&! out.wrfvar
          cp fg wrfvar_output
          cp ${CENTRALDIR}/add_bank_perts.ncl .
          set cmd3 = "ncl 'MEM_NUM=${ensemble_member}' 'PERTS_DIR="\""${PERTS_DIR}"\""' ${CENTRALDIR}/advance_temp${ensemble_member}/add_bank_perts.ncl"
          ${REMOVE} nclrun3.out
          cat >! nclrun3.out << EOF
          $cmd3
EOF
          chmod +x nclrun3.out
          ./nclrun3.out >& add_perts.out
          if ( -z add_perts.err ) then
            echo "Perts added to member ${ensemble_member}"
          else
             echo "Error! Non-zero status returned from add_bank_perts.ncl. Check ${RUN_DIR}/advance_temp${ensemble_member}/add_perts.err."
             cat add_perts.err
             exit
          endif
################################
         cp namelist.input namelist.input.3dvar
	 if ( -e rsl.out.0000 ) cat rsl.out.0000 >> out.wrfvar

         ${MOVE} wrfvar_output wrfinput_next
         ${LN} wrfinput_d01 wrfinput_this
         ${LN} wrfbdy_d01 wrfbdy_this

         # if wrfinput_mean file found, rename it
         if ( -e wrfinput_mean ) then
            ${MOVE} wrfinput_mean   wrfinput_this_mean
            ${MOVE} fg              wrfinput_next_mean
         endif

         ${CENTRALDIR}/pert_wrf_bc >&! out.pert_wrf_bc
         ${REMOVE} wrfinput_this wrfinput_next wrfbdy_this
         if ( -e wrfinput_this_mean ) ${REMOVE} wrfinput_this_mean wrfinput_next_mean

      else  # Update boundary conditions from existing wrfbdy files

         echo $infl | ${CENTRALDIR}/update_wrf_bc >&! out.update_wrf_bc

      endif

      ${REMOVE} script.sed namelist.input
      cat >! script.sed << EOF
         /run_hours/c\
         run_hours                  = 0,
         /run_minutes/c\
         run_minutes                = 0,
         /run_seconds/c\
         run_seconds                = ${INTERVAL_SS},
         /start_year/c\
         start_year                 = ${MY_NUM_DOMAINS}*${START_YEAR},
         /start_month/c\
         start_month                = ${MY_NUM_DOMAINS}*${START_MONTH},
         /start_day/c\
         start_day                  = ${MY_NUM_DOMAINS}*${START_DAY},
         /start_hour/c\
         start_hour                 = ${MY_NUM_DOMAINS}*${START_HOUR},
         /start_minute/c\
         start_minute               = ${MY_NUM_DOMAINS}*${START_MIN},
         /start_second/c\
         start_second               = ${MY_NUM_DOMAINS}*${START_SEC},
         /end_year/c\
         end_year                   = ${MY_NUM_DOMAINS}*${END_YEAR},
         /end_month/c\
         end_month                  = ${MY_NUM_DOMAINS}*${END_MONTH},
         /end_day/c\
         end_day                    = ${MY_NUM_DOMAINS}*${END_DAY},
         /end_hour/c\
         end_hour                   = ${MY_NUM_DOMAINS}*${END_HOUR},
         /end_minute/c\
         end_minute                 = ${MY_NUM_DOMAINS}*${END_MIN},
         /end_second/c\
         end_second                 = ${MY_NUM_DOMAINS}*${END_SEC},
         /history_interval/c\
         history_interval           = ${MY_NUM_DOMAINS}*${INTERVAL_MIN},
         /frames_per_outfile/c\
         frames_per_outfile         = ${MY_NUM_DOMAINS}*1,
         /max_dom/c\
         max_dom                    = ${MY_NUM_DOMAINS},
EOF
# The EOF on the line above MUST REMAIN in column 1.

      if ( -e ${CENTRALDIR}/fixed_domain_info ) then

         set nx_string      = `head -n 2 ${CENTRALDIR}/fixed_domain_info | tail -n 1`
         set ny_string      = `head -n 3 ${CENTRALDIR}/fixed_domain_info | tail -n 1`
         set i_start_str    = `head -n 4 ${CENTRALDIR}/fixed_domain_info | tail -n 1`
         set j_start_str    = `head -n 5 ${CENTRALDIR}/fixed_domain_info | tail -n 1`

         cat >> script.sed << EOF
         /e_we/c\
         e_we                                = ${nx_string},
         /e_sn/c\
         e_sn                                = ${ny_string},
         /i_parent_start/c\
         i_parent_start                      = ${i_start_str},
         /j_parent_start/c\
         j_parent_start                      = ${j_start_str},
EOF

      else if ( -e ${CENTRALDIR}/moving_domain_info ) then

         set nx_string      = `head -n 3  ${CENTRALDIR}/moving_domain_info | tail -n 1`
         set ny_string      = `head -n 4  ${CENTRALDIR}/moving_domain_info | tail -n 1`
         set i_start_str    = `head -n 5  ${CENTRALDIR}/moving_domain_info | tail -n 1`
         set j_start_str    = `head -n 6  ${CENTRALDIR}/moving_domain_info | tail -n 1`
         set input_file     = `head -n 7  ${CENTRALDIR}/moving_domain_info | tail -n 1`
         set num_move_str   = `head -n 8  ${CENTRALDIR}/moving_domain_info | tail -n 1`
         set id_move_str    = `head -n 9  ${CENTRALDIR}/moving_domain_info | tail -n 1`
	 set move_time_str  = `head -n 10 ${CENTRALDIR}/moving_domain_info | tail -n 1`
	 set x_move_string  = `head -n 11 ${CENTRALDIR}/moving_domain_info | tail -n 1`
	 set y_move_string  = `head -n 12 ${CENTRALDIR}/moving_domain_info | tail -n 1`

         cat >> script.sed << EOF
         /e_we/c\
         e_we                                = ${nx_string},
         /e_sn/c\
         e_sn                                = ${ny_string},
         /i_parent_start/c\
         i_parent_start                      = ${i_start_str},
         /j_parent_start/c\
         j_parent_start                      = ${j_start_str},
	 /input_from_file/c\
	 input_from_file                     = ${input_file},
	 /num_moves/c\
         num_moves                           = ${num_move_str},
         /move_id/c\
         move_id                             = ${id_move_str}
         /move_interval/c\
         move_interval                       = ${move_time_str}
         /move_cd_x/c\
         move_cd_x                           = ${x_move_string}
         /move_cd_y/c\
         move_cd_y                           = ${y_move_string}
EOF

      endif

      if ( $ensemble_member <= 1 ) then
         echo "     /auxhist1_interval/c\"            >> script.sed
         echo "     auxhist1_interval      = 0, 3, 3" >> script.sed
      endif

      sed -f script.sed ${CENTRALDIR}/namelist.input >! namelist.input

      #-------------------------------------------------------------
      #
      # HERE IS A GOOD PLACE TO GRAB FIELDS FROM OTHER SOURCES
      # AND STUFF THEM INTO YOUR wrfinput_d0? FILES
      #
      #------------------------------------------------------------

      # clean out any old rsl files
      if ( -e rsl.out.integration )  ${REMOVE} rsl.*

      # run WRF here
      setenv MPI_SHEPHERD FALSE
      ${ADV_MOD_COMMAND} >>&! rsl.out.integration

      if ( -e rsl.out.0000 ) cat rsl.out.0000 >> rsl.out.integration
      ${COPY} rsl.out.integration ${WRFOUTDIR}/wrf.out_${targdays}_${targsecs}_${ensemble_member}
#      sleep 1

      set SUCCESS = `grep "wrf: SUCCESS COMPLETE WRF" rsl.out.integration | cat | wc -l`
      if ($SUCCESS == 0) then
         echo $ensemble_member >>! ${CENTRALDIR}/blown_${targdays}_${targsecs}.out
         echo "Model failure! Check file " ${CENTRALDIR}/blown_${targdays}_${targsecs}.out
         echo "for a list of failed ensemble_members, and check here for the individual output files:"
         echo " ${CENTRALDIR}/wrf.out_${targdays}_${targsecs}_${ensemble_member}  "
         exit -1
      endif

      if ( -e ${CENTRALDIR}/append_precip_to_diag ) then
         set dn = 1  ;  which ncks >& /dev/null
         while ( $dn <= $num_domains )
            ncks -h -O -F -v RAINC,RAINNC wrfout_d0${dn}_${END_STRING} wrf_precip.nc
            ${MOVE} wrf_precip.nc ${CENTRALDIR}/wrf_precip_d0${dn}_${END_STRING}_${ensemble_member}
            @ dn ++
         end
      endif

#      if ( $ensemble_member <= 1 ) then
#         set dn = 2
#         while ( $dn <= $num_domains )
#            ${COPY} psfc_data_d0${dn}.nc wrfout.nc
#            ${CENTRALDIR}/domain_psfc_tend
#            ${MOVE} psfc_stats.nc ${CENTRALDIR}/psfc_stat_d0${dn}.nc
#            @ dn ++
#         end
#      endif

# zip up the wrfin file
      set dn = 1
      while ( $dn <= $num_domains )
         ${MOVE} wrfinput_d0${dn} wrfinput_d0${dn}_${ensemble_member} 
         gzip wrfinput_d0${dn}_${ensemble_member} &
         @ dn ++
      end

# forecast date
      set dn = 1
      while ( $dn <= $num_domains )
         if ( $ensemble_member <= $save_ensemble_member ) ${COPY} wrfout_d0${dn}_${END_STRING} ${WRFOUTDIR}/wrfout_d0${dn}_${END_STRING}_${ensemble_member}
# if the wrfinput file zip operation is finished, wrfinput_d0${dn}_$ensemble_member should no 
# longer be in the directory
# test for this, and wait if the zip operation is not yet finished
          while ( -e wrfinput_d0${dn}_${ensemble_member} )
             sleep 3
             touch ${CENTRALDIR}/HAD_TO_WAIT
          end
          ${MOVE} wrfinput_d0${dn}_${ensemble_member}.gz ../WRFIN/wrfinput_d0${dn}_${ensemble_member}.gz
          ${MOVE} wrfout_d0${dn}_${END_STRING} wrfinput_d0${dn}
          @ dn ++
      end

      ${REMOVE} wrfout*

      set START_YEAR  = $END_YEAR
      set START_MONTH = $END_MONTH
      set START_DAY   = $END_DAY
      set START_HOUR  = $END_HOUR
      set START_MIN   = $END_MIN
      set START_SEC   = $END_SEC
      set wrfkey      = $keys[$ifile]
      @ ifile ++

   end

   ##############################################
   # At this point, the target time is reached. #
   ##############################################
   # withdraw LSM data to use in next cycle   This is remnant from the Lanai days, we now pull soil state
   # together with everything else
   if ( -e ${CENTRALDIR}/fixed_domain_info )  set MY_NUM_DOMAINS = 1
   if ( -e ${CENTRALDIR}/append_lsm_data ) then
      set dn = 1
      while ( $dn <= $num_domains )
         ncks -h -F -A -a -v TSLB,SMOIS,SH2O,TSK wrfinput_d0${dn} lsm_data.nc
         ncrename -h -v TSLB,TSLB_d0${dn} -v SMOIS,SMOIS_d0${dn} -v SH2O,SH2O_d0${dn} -v TSK,TSK_d0${dn} \
                  -d west_east,west_east_d0${dn} -d south_north,south_north_d0${dn} \
                  -d soil_layers_stag,soil_layers_stag_d0${dn} lsm_data.nc
	 @ dn ++
      end
      ${REMOVE} ${CENTRALDIR}/LSM/lsm_data_${ensemble_member}.nc
      ${MOVE} lsm_data.nc ${CENTRALDIR}/LSM/lsm_data_${ensemble_member}.nc
   endif

   if ( -e ${CENTRALDIR}/fixed_domain_info || -e ${CENTRALDIR}/moving_domain_info ) then
      ln -sf ${CENTRALDIR}/wrfinput_d01 wrfinput_d01_base
      ${CENTRALDIR}/recalc_wrf_base >&! out.recalc_wrf_base
   endif
#   extract the cycle variables
#   # create new input to DART (taken from "wrfinput")
#   ${CENTRALDIR}/wrf_to_dart >&! out.wrf_to_dart
#   ${MOVE} dart_wrf_vector ${CENTRALDIR}/${output_file}
  set num_vars = $#extract_vars_a
  set extract_str_a = ''
  set i = 1
  while ( $i < $num_vars )
     set extract_str_a = `echo ${extract_str_a}$extract_vars_a[$i],`
     @ i ++
  end
  set extract_str_a = `echo ${extract_str_a}$extract_vars_a[$num_vars]`
  echo ${extract_str_a}

  set num_vars = $#extract_vars_b
  set extract_str_b = ''
  set i = 1
  while ( $i < $num_vars )
     set extract_str_b = `echo ${extract_str_b}$extract_vars_b[$i],`
     @ i ++
  end
  set extract_str_b = `echo ${extract_str_b}$extract_vars_b[$num_vars]`
  echo ${extract_str_b}


# MULTIPLE DOMAINS - loop through wrf files that are present
  set dn = 1
  while ( $dn <= $num_domains )
     set dchar = `echo $dn + 100 | bc | cut -b2-3`
     set icnum = `echo $ensemble_member + 10000 | bc | cut -b2-5`
     set outfile =  prior_d${dchar}.${icnum}
     if ( $dn == 1) then
       ncks -O -v ${extract_str_a} wrfinput_d${dchar} ../$outfile   # MULTIPLE DOMAINS - output file is incomplete filename?
     else
       ncks -O -v ${extract_str_b} wrfinput_d${dchar} ../$outfile   # MULTIPLE DOMAINS - output file is incomplete filename?
     endif
     @ dn ++
     echo "should have made $outfile"
  end
# MULTIPLE DOMAINS - may need to remove below to avoid confusion
   if ( -e ${CENTRALDIR}/moving_domain_info && $ensemble_member == 1 ) then
      set dn = 2
      while ( $dn <= $num_domains )
         ${COPY} wrfinput_d0${dn} ${CENTRALDIR}/wrfinput_d0${dn}_new
         @ dn ++
      end
   endif


   touch ${CENTRALDIR}/done_member_$ensemble_member

   cd $CENTRALDIR

   #  delete the temp directory for each member if desired
   if ( $delete_temp_dir == true )  ${REMOVE} ${temp_dir}
   echo "Ensemble Member $ensemble_member completed"

   # and now repeat the entire process for any other ensemble member that
   # needs to be advanced by this task.
   # don't expect this to ever be run 
   @ state_copy ++
   @ ensemble_member_line += 3
   @ input_file_line += 3
   @ output_file_line += 3

end

# Remove the filter_control file to signal completion
# Is there a need for any sleeps to avoid trouble on completing moves here?
      ${REMOVE} $control_file
      if ($SUCCESS == 1) then
         echo " done_member_$ensemble_member"
      endif
exit 0

