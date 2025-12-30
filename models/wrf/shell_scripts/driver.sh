#!/bin/bash
#
# DART software - Copyright UCAR. This open source software is provided
# by UCAR, "as is", without charge, subject to all terms of use at
# http://www.image.ucar.edu/DAReS/DART/DART_download

#   driver.sh - script that performs assimilation cycling
#
#   run script as: ./driver.sh 2024051906 param.sh >& run.out &
#
#      This provides an input argument of the first
#      analysis time in yyyymmddhh format.
########################################################################
# Set the correct values here
paramfile="$(readlink -f "${2}")" # Get absolute path for param.sh from command line arg
datefnl=2024051912 # target date   YYYYMMDDHH  # set this appropriately #%%%#
########################################################################

source "$paramfile"
uname -a
cd "${RUN_DIR}"

#  First determine the appropriate analysis date

if (( $# > 0 )); then
  datea="${1}" # starting date
  export restore=1   # set the restore variable
  echo 'starting a restore'
else
  echo "please enter a date: yyyymmddhh"
  exit 1
fi

touch "${RUN_DIR}/cycle_started_${datea}"

while true; do

   if [[ ! -d "${OUTPUT_DIR}/${datea}" && "${restore}" == "1" ]]; then
      ${REMOVE} "${RUN_DIR}/ABORT_RETRO"
      echo 'exiting because output directory does not exist and this is a restore'
      exit 0
   fi

   datep="$(echo "${datea} -${ASSIM_INT_HOURS}"   | "${DART_DIR}/models/wrf/work/advance_time")"
   read -r -a gdate  < <(echo "${datea} 0 -g"                  | "${DART_DIR}/models/wrf/work/advance_time")
   read -r -a gdatef < <(echo "${datea} ${ASSIM_INT_HOURS} -g" | "${DART_DIR}/models/wrf/work/advance_time")
   wdate="$(echo "${datea} 0 -w"                  | "${DART_DIR}/models/wrf/work/advance_time")"
   hh="${datea:8:2}"

   echo 'ready to check inputs'
   domains="${NUM_DOMAINS}"   # from the param file

   #  Check to make sure all input data exists
   dn=1
   while (( dn <= domains )); do
      dchar="$(echo "${dn} + 100" | bc | cut -b2-3)"

      for infile in \
         "wrfinput_d${dchar}_${gdate[0]}_${gdate[1]}_mean" \
         "wrfinput_d${dchar}_${gdatef[0]}_${gdatef[1]}_mean" \
         "wrfbdy_d01_${gdatef[0]}_${gdatef[1]}_mean" \
         "obs_seq.out"
      do
         if [[ ! -e "${OUTPUT_DIR}/${datea}/${infile}" ]]; then
            echo  "${OUTPUT_DIR}/${datea}/${infile} is missing!  Stopping the system"
            touch ABORT_RETRO
            exit 2
         fi
      done

      (( dn++ ))
   done        # loop through domains

   #  Clear the advance_temp directory, write in new template file, and
   # overwrite variables with the compact prior netcdf files
   #
   #  NOTE that multiple domains might be present, but only looking for domain 1

   if [[ "${SUPER_PLATFORM}" == "LSF queuing system" ]]; then
      ic_queue="caldera"
      logfile="${RUN_DIR}/ic_gen.log"
      sub_command=( bsub -q "${ic_queue}" -W 00:05 -o "${logfile}" -n 1 -P "${COMPUTER_CHARGE_ACCOUNT}" )
   elif [[ "${SUPER_PLATFORM}" == "derecho" ]]; then
      ic_queue="main"
      # NOTE: qsub flags here match the original intent; users may adjust for their environment.
      sub_command=( qsub -l "select=1:ncpus=128:mpiprocs=128:mem=5GB" -l "walltime=00:03:00" -q "${ic_queue}" -A "${COMPUTER_CHARGE_ACCOUNT}" -j oe -k eod -N icgen )
   else
      echo "Unknown SUPER_PLATFORM='${SUPER_PLATFORM}'"
      exit 2
   fi

   echo "this platform is $SUPER_PLATFORM and the job submission command is ${sub_command[*]}"

   dn=1
   while (( dn <= domains )); do
      dchar="$(echo "${dn} + 100" | bc | cut -b2-3)"
      n=1
      while (( n <= NUM_ENS )); do
         ensstring="$(echo "${n} + 10000" | bc | cut -b2-5)"
         if [[ -e "${OUTPUT_DIR}/${datep}/PRIORS/prior_d${dchar}.${ensstring}" ]]; then

            if (( dn == 1 )) && [[ -d "${RUN_DIR}/advance_temp${n}" ]]; then
               ${REMOVE} "${RUN_DIR}/advance_temp${n}"
            fi

            mkdir -p "${RUN_DIR}/advance_temp${n}"
            ${LINK} "${OUTPUT_DIR}/${datea}/wrfinput_d${dchar}_${gdate[0]}_${gdate[1]}_mean" \
                    "${RUN_DIR}/advance_temp${n}/wrfinput_d${dchar}"
         else
            echo "${OUTPUT_DIR}/${datep}/PRIORS/prior_d${dchar}.${ensstring} is missing! Stopping the system"
            touch ABORT_RETRO
            exit 3
         fi
         (( n++ ))
      done  # loop through ensemble members
      (( dn++ ))
   done   # loop through domains

   # Fire off a bunch of small jobs to create the initial conditions for the short model forecast.
   # the prep_ic.sh script creates a file "${RUN_DIR}/ic_d${dchar}_${n}_ready" to signal a
   # successful completion.
   # NOTE : Submit commands here are system specific and work for this tutorial, users may want/need to change
   #        for their system and/or production.

   n=1
   while (( n <= NUM_ENS )); do
      if [[ "${SUPER_PLATFORM}" == "derecho" ]]; then   # can't pass along arguments in the same way
         "${sub_command[@]}" -v "mem_num=${n},date=${datep},domain=${domains},paramf=${paramfile}" "${SHELL_SCRIPTS_DIR}/prep_ic.sh"
      else
         # LSF: pass args directly
         "${sub_command[@]}" "${SHELL_SCRIPTS_DIR}/prep_ic.sh" "${n}" "${datep}" "${dn}" "${paramfile}"
      fi
      (( n++ ))
   done  # loop through ensemble members

   # If any of the queued jobs has not completed in 5 minutes, run them manually
   # cleanup any failed stuffs
   # NOTE : No automated cleanup for queued jobs. User may want to add system specific monitoring.
   dn=1
   while (( dn <= domains )); do
      dchar="$(echo "${dn} + 100" | bc | cut -b2-3)"
      n=1
      loop=1
      while (( n <= NUM_ENS )); do
         if [[ -e "${RUN_DIR}/ic_d${dchar}_${n}_ready" ]]; then
            ${REMOVE} "${RUN_DIR}/ic_d${dchar}_${n}_ready"
            (( n++ ))
            loop=1
         else
            echo "waiting for ic member $n in domain $dn"
            sleep 5
            (( loop++ ))
            if (( loop > 60 )); then    # wait 5 minutes for the ic file to be ready, else run manually
               echo "gave up on ic member $n - redo"
               "${SHELL_SCRIPTS_DIR}/prep_ic.sh" "${n}" "${datep}" "${dn}" "${paramfile}"
               # If manual execution of script, shouldn't queued job be killed?
            fi
         fi
      done
      (( dn++ ))
   done   # loop through domains

   mkdir -p "${OUTPUT_DIR}/${datea}/logs"
   ${MOVE} icgen.o* "${OUTPUT_DIR}/${datea}/logs/" 2>/dev/null || true

   #  Get wrfinput source information
   dn=1
   while (( dn <= domains )); do
      dchar="$(echo "${dn} + 100" | bc | cut -b2-3)"
      ${COPY} "${OUTPUT_DIR}/${datea}/wrfinput_d${dchar}_${gdate[0]}_${gdate[1]}_mean" "wrfinput_d${dchar}"
      (( dn++ ))
   done

   # Copy the inflation files from the previous time and for all  domains
   # The ADAPTIVE_INFLATION variable is set in scripts/param.sh and should
   # be consistent with DART's input.nml inflation setting (inf_flavor)

   if [[ "${ADAPTIVE_INFLATION}" == "1" ]]; then
      # Create the home for inflation and future state space diagnostic files
      mkdir -p "${RUN_DIR}/Inflation_input" "${RUN_DIR}/Output"

      dn=1
      while (( dn <= domains )); do
         dchar="$(echo "${dn} + 100" | bc | cut -b2-3)"

         if [[ -e "${OUTPUT_DIR}/${datep}/Inflation_input/input_priorinf_mean_d${dchar}.nc" ]]; then

            ${LINK} "${OUTPUT_DIR}/${datep}/Inflation_input/input_priorinf_mean_d${dchar}.nc" "${RUN_DIR}/."
            ${LINK} "${OUTPUT_DIR}/${datep}/Inflation_input/input_postinf_mean_d${dchar}.nc"  "${RUN_DIR}/."

            ${LINK} "${OUTPUT_DIR}/${datep}/Inflation_input/input_priorinf_sd_d${dchar}.nc"   "${RUN_DIR}/."
            ${LINK} "${OUTPUT_DIR}/${datep}/Inflation_input/input_postinf_sd_d${dchar}.nc"   "${RUN_DIR}/."

         else

            echo "${OUTPUT_DIR}/${datep}/Inflation_input/input_priorinf_mean_d${dchar}.nc file does not exist. Stopping"
            echo "If first assimilation cycle make sure fill_inflation_restart was used to generate mean and sd inflation files"
            touch ABORT_RETRO
            exit 3

         fi
         (( dn++ ))
      done # Loop through domains

   fi   # ADAPTIVE_INFLATION file check

   ${LINK} "${OUTPUT_DIR}/${datea}/obs_seq.out" .
   ${REMOVE} "${RUN_DIR}/WRF"
   ${REMOVE} "${RUN_DIR}/prev_cycle_done"
   ${LINK} "${OUTPUT_DIR}/${datea}" "${RUN_DIR}/WRF"

   #  Run filter to generate the analysis
   ${REMOVE} script.sed 2>/dev/null || true

   assimilate_job="${RUN_DIR}/assimilate.sh"

   if [[ "${SUPER_PLATFORM}" == "LSF queuing system" ]]; then

     cat > "${assimilate_job}" <<EOF
#!/bin/bash
#==================================================================
#BSUB -J assimilate_${datea}
#BSUB -o assimilate_${datea}.%J.log
#BSUB -P ${COMPUTER_CHARGE_ACCOUNT}
#BSUB -W ${FILTER_TIME}
#BSUB -q ${FILTER_QUEUE}
#BSUB -n ${FILTER_CORES}
#BSUB -x
#BSUB -R "span[ptile=${FILTER_PTILE}]"
#==================================================================

${SHELL_SCRIPTS_DIR}/assimilate.sh "${datea}" "${paramfile}"
EOF

     chmod +x "${assimilate_job}"

     if [[ -n "${reservation:-}" ]]; then
       bsub -U "$(/contrib/lsf/get_my_rsvid)" < "${assimilate_job}"
     else
       bsub < "${assimilate_job}"
     fi

   elif [[ "${SUPER_PLATFORM}" == "derecho" ]]; then

     cat > "${assimilate_job}" <<EOF
#!/bin/bash
#=================================================================
#PBS -N assimilate_${datea}
#PBS -j oe
#PBS -A ${COMPUTER_CHARGE_ACCOUNT}
#PBS -l walltime=${FILTER_TIME}
#PBS -q ${FILTER_QUEUE}
#PBS -l job_priority=${FILTER_PRIORITY}
#PBS -m ae
#PBS -M ${EMAIL}
#PBS -k eod
#PBS -l select=${FILTER_NODES}:ncpus=${FILTER_PROCS}:mpiprocs=${FILTER_MPI}
#=================================================================

${SHELL_SCRIPTS_DIR}/assimilate.sh "${datea}" "${paramfile}"
EOF

     chmod +x "${assimilate_job}"
     qsub "${assimilate_job}"
   fi


   cd "${RUN_DIR}"  

   filter_thresh_min="$(echo "${FILTER_TIME}" | cut -b3-4)"
   filter_thresh_hr="$(echo "${FILTER_TIME}" | cut -b1-1)"
   filter_thresh=$(( 10#${filter_thresh_min} * 60 + 10#${filter_thresh_hr} * 3600 ))

   while [[ ! -e filter_done ]]; do

      # Check the timing.  If it took longer than the time allocated, abort.
      if [[ -e filter_started ]]; then

         start_time="$(head -1 filter_started)"
         end_time="$(date +%s)"
         total_time=$(( end_time - start_time ))

         if (( total_time > filter_thresh )); then

            # If the job needs to be aborted ... we need to qdel the hanging job

            echo "Time exceeded the maximum allowable time.  Exiting."
            touch ABORT_RETRO
            ${REMOVE} filter_started
            exit 5

         fi

      fi
      sleep 10

   done

   echo "filter is done, cleaning up"

   ${MOVE} icgen.o* "${OUTPUT_DIR}/${datea}/logs/" 2>/dev/null || true
   ${REMOVE} "${RUN_DIR}/filter_started" \
             "${RUN_DIR}/filter_done" \
             "${RUN_DIR}/obs_seq.out" \
             "${RUN_DIR}/postassim_priorinf"* \
             "${RUN_DIR}/preassim_priorinf"* 2>/dev/null || true

   if [[ -e assimilate.sh ]]; then
      ${REMOVE} "${RUN_DIR}/assimilate.sh"
   fi

   echo "Listing contents of rundir before archiving at $(date)"
   ls -l *.nc blown* dart_log* filter_* input.nml obs_seq* Output/inf_ic* 2>/dev/null || true
   mkdir -p "${OUTPUT_DIR}/${datea}/Inflation_input" "${OUTPUT_DIR}/${datea}/WRFIN" "${OUTPUT_DIR}/${datea}/PRIORS" "${OUTPUT_DIR}/${datea}/logs"

   extract_str=""
   if declare -p increment_vars_a &>/dev/null; then
      for v in "${increment_vars_a[@]}"; do
         if [[ -z "${extract_str}" ]]; then
            extract_str="${v}"
         else
            extract_str="${extract_str},${v}"
         fi
      done
   fi

   # Create an analysis increment file that has valid static data.
   # First, create the difference of a subset of variables
   # Second, create a netCDF file with just the static data
   # Third, append the static data onto the difference.
   dn=1
   while (( dn <= domains )); do
        dchar="$(echo "${dn} + 100" | bc | cut -b2-3)"
        ncdiff -F -O -v "${extract_str}" "postassim_mean_d${dchar}.nc" "preassim_mean_d${dchar}.nc" "analysis_increment_d${dchar}.nc"
        ncks  -F -O -x -v "${extract_str}" "postassim_mean_d${dchar}.nc" "static_data_d${dchar}.nc"
        ncks  -A "static_data_d${dchar}.nc" "analysis_increment_d${dchar}.nc"

        # Move diagnostic and obs_seq.final data to storage directories
        #
        if (( dn == 1 )) && [[ -e obs_seq.final ]]; then
             ${MOVE} obs_seq.final "${OUTPUT_DIR}/${datea}/."
             if [[ $? -ne 0 ]]; then
                 echo "failed moving ${RUN_DIR}/obs_seq.final"
                 touch BOMBED
             fi
        else
              if (( dn == 1 )); then
                 echo "${OUTPUT_DIR}/obs_seq.final does not exist and should."
                 echo "Stopping driver.sh"
                 exit 0
                 
              fi
        fi

        for FILE in \
          "postassim_mean_d${dchar}.nc" "preassim_mean_d${dchar}.nc" \
          "postassim_sd_d${dchar}.nc"   "preassim_sd_d${dchar}.nc" \
          "analysis_increment_d${dchar}.nc" \
          "output_mean_d${dchar}.nc" "output_sd_d${dchar}.nc"
        do
          if [[ -e "${FILE}" && -s "${FILE}" ]]; then
             ${MOVE} "${FILE}" "${OUTPUT_DIR}/${datea}/."
             if [[ $? -ne 0 ]]; then
                echo "failed moving ${RUN_DIR}/${FILE}"
                touch BOMBED
             fi
          else
             echo "${OUTPUT_DIR}/${FILE} does not exist and should."
             ls -l
             touch BOMBED
          fi
        done

      (( dn++ ))
   done # loop through domains

   echo "past the analysis file moves"

   # Move inflation files to storage directories
   # The output inflation file is used as the input for the next cycle,
   # so rename the current inflation output for the next cycle input.
   cd "${RUN_DIR}"

   if [[ "${ADAPTIVE_INFLATION}" == "1" ]]; then

      dn=1
      while (( dn <= domains )); do
         dchar="$(echo "${dn} + 100" | bc | cut -b2-3)"
         old_file=( "input_postinf_mean_d${dchar}.nc"  "input_postinf_sd_d${dchar}.nc"  "input_priorinf_mean_d${dchar}.nc"  "input_priorinf_sd_d${dchar}.nc" )
         new_file=( "output_postinf_mean_d${dchar}.nc" "output_postinf_sd_d${dchar}.nc" "output_priorinf_mean_d${dchar}.nc" "output_priorinf_sd_d${dchar}.nc" )

         nfiles="${#new_file[@]}"
         i=0
         while (( i < nfiles )); do
            if [[ -e "${new_file[$i]}" && -s "${new_file[$i]}" ]]; then
               ${MOVE} "${new_file[$i]}" "${OUTPUT_DIR}/${datea}/Inflation_input/${old_file[$i]}"
               if [[ $? -ne 0 ]]; then
                  echo "failed moving ${RUN_DIR}/Output/${new_file[$i]}"
                  touch BOMBED
               fi
            fi
            (( i++ ))
         done
         (( dn++ ))
      done  # loop through domains

      echo "past the inflation file moves"
   fi   # adaptive_inflation file moves

   # Submit jobs to integrate ensemble members to next analysis time ...
   # BEFORE calculating the observation-space diagnostics for the existing cycle.

   echo "ready to integrate ensemble members"

   # Removing old start_member and done_member diagnostics
   if [[ -e "${RUN_DIR}/start_member_1" ]]; then
      ${REMOVE} "${RUN_DIR}/start_member_"* "${RUN_DIR}/done_member_"* 2>/dev/null || true
   fi

   n=1
   while (( n <= NUM_ENS )); do

     if [[ "${SUPER_PLATFORM}" == "LSF queuing system" ]]; then

       jobfile="assim_advance_mem${n}.sh"
       cat > "${jobfile}" <<EOF
#!/bin/bash
#==================================================================
#BSUB -J assim_advance_${n}
#BSUB -o assim_advance_${n}.%J.log
#BSUB -P ${COMPUTER_CHARGE_ACCOUNT}
#BSUB -W ${ADVANCE_TIME}
#BSUB -q ${ADVANCE_QUEUE}
#BSUB -n ${ADVANCE_CORES}
#BSUB -x
#BSUB -R "span[ptile=${ADVANCE_PTILE}]"
#==================================================================

set -uo pipefail

# Run the actual advance script (no scheduler headers inside it)
"${SHELL_SCRIPTS_DIR}/assim_advance.sh" "${datea}" "${n}" "${paramfile}"
EOF
       chmod +x "${jobfile}"

       if [[ -n "${reservation:-}" ]]; then
         echo "MEMBER ${n} USING RESERVATION," "$(/contrib/lsf/get_my_rsvid)"
         bsub -U "$(/contrib/lsf/get_my_rsvid)" < "${jobfile}"
       else
         bsub < "${jobfile}"
       fi

     elif [[ "${SUPER_PLATFORM}" == "derecho" ]]; then

       jobfile="assim_advance_mem${n}.sh"

       cat > "${jobfile}" <<EOF
#!/bin/bash
#=================================================================
#PBS -N assim_advance_${n}
#PBS -j oe
#PBS -A ${COMPUTER_CHARGE_ACCOUNT}
#PBS -l walltime=${ADVANCE_TIME}
#PBS -q ${ADVANCE_QUEUE}
#PBS -l job_priority=${ADVANCE_PRIORITY}
#PBS -m a
#PBS -M ${EMAIL}
#PBS -k eod
#PBS -l select=${ADVANCE_NODES}:ncpus=${ADVANCE_PROCS}:mpiprocs=${ADVANCE_MPI}
#=================================================================

set -uo pipefail

cd "${RUN_DIR}"

# Run the actual advance script (no scheduler headers inside it)
"${SHELL_SCRIPTS_DIR}/assim_advance.sh" "${datea}" "${n}" "${paramfile}"
EOF
       chmod +x "${jobfile}"

       qsub "${jobfile}"
	
     fi

     (( n++ ))
   done



   #  Compute Diagnostic Quantities (in the background)

   if [[ -e obs_diag.log ]]; then
      ${REMOVE} obs_diag.log
   fi
   "${SHELL_SCRIPTS_DIR}/diagnostics_obs.sh" "${datea}" "${paramfile}" >& "${RUN_DIR}/obs_diag.log" &

   # ---------------------------------------------------------------------------
   #  Check to see if all of the ensemble members have advanced
   # ---------------------------------------------------------------------------

   advance_thresh_min="$(echo "${ADVANCE_TIME}" | cut -b3-4)"
   advance_thresh_hr="$(echo "${ADVANCE_TIME}" | cut -b1-1)"
   advance_thresh=$(( 10#${advance_thresh_min} * 60 + 10#${advance_thresh_hr} * 3600 ))

   n=1
   while (( n <= NUM_ENS )); do

      ensstring="$(echo "${n} + 10000" | bc | cut -b2-5)"
      keep_trying=true
      max_retry=1

      while [[ "${keep_trying}" == "true" ]]; do

         #  Wait for the script to start
         while [[ ! -e "${RUN_DIR}/start_member_${n}" ]]; do

            if [[ "${SUPER_PLATFORM}" == "LSF queuing system" ]]; then

               if [[ "$(bjobs -w | grep -c "assim_advance_${n}" || true)" -eq 0 ]]; then
                  echo "assim_advance_${n} is missing from the queue"
                  if [[ -n "${reservation:-}" ]]; then
                     echo "MEMBER ${n} USING RESERVATION," "$(/contrib/lsf/get_my_rsvid)"
                     bsub -U "$(/contrib/lsf/get_my_rsvid)" < "assim_advance_mem${n}.sh"
                  else
                     bsub < "assim_advance_mem${n}.sh"
                  fi
               fi

            elif [[ "${SUPER_PLATFORM}" == "derecho" ]]; then

               if [[ "$(qstat -wa | grep -c "assim_advance_${n}" || true)" -eq 0 ]]; then

                  echo "Warning, detected that assim_advance_${n} is missing from the queue"
                  echo "If this warning leads to  missing output from ensemble ${n}"
                  echo "consider enabling the qsub command within keep_trying while statement in driver.sh"

                  #qsub assim_advance_mem${n}.sh
               fi

            fi
            sleep 5

         done

         start_time="$(head -1 "start_member_${n}")"
         echo "Member $n has started.  Start time $start_time"

         #  Wait for the output file
         while true; do

            current_time="$(date +%s)"
            length_time=$(( current_time - start_time ))

            if [[ -e "${RUN_DIR}/done_member_${n}" ]]; then

               #  If the output file already exists, move on
               keep_trying=false
               break

            elif (( length_time > advance_thresh )); then

               #  If WRF member has failed 2 resubmission attempts, immediately stop driver.sh
               if (( max_retry > 2 )); then

                  echo "Stopping the driver.sh script! The WRF ensemble member ${n}"
                  echo "has exceeded the maximum resubmission attempts (2) without completing."
                  echo "This typically means the WRF integration has failed."
                  echo "Check your BASE_DIR/rundir/advance_temp${n} directory and locate"
                  echo "the WRF rsl.out.0000 or rsl.error.0000 log files for further information."
                  echo "If applicable, check the DART analysis_increment.nc from previous assimilation step"
                  exit 7

               fi

               # The WRF job did not complete. Resubmit to queue
               ${REMOVE} "start_member_${n}"
               echo "Did not find the member done file, WRF run did not complete"
               echo "Attempting resubmission $max_retry"
               (( max_retry++ ))

               if [[ "${SUPER_PLATFORM}" == "LSF queuing system" ]]; then

                  if [[ -n "${reservation:-}" ]]; then
                     echo "MEMBER ${n} USING RESERVATION," "$(/contrib/lsf/get_my_rsvid)"
                     bsub -U "$(/contrib/lsf/get_my_rsvid)" < "assim_advance_mem${n}.sh"
                  else
                     bsub < "assim_advance_mem${n}.sh"
                  fi

               elif [[ "${SUPER_PLATFORM}" == "derecho" ]]; then

                  qsub "assim_advance_mem${n}.sh"
                  sleep 5
               fi
               break

            fi
            sleep 15    # this might need to be longer, though I moved the done flag lower in the
                        # advance_model.sh to hopefully avoid the file moves below failing
         done

      done

      #  Move output data to correct location
      dn=1
      while (( dn <= domains )); do
          dchar="$(echo "${dn} + 100" | bc | cut -b2-3)"
          echo "moving ${n} ${ensstring} for domain ${dn}"

          if (( dn == 1 )); then
             ${MOVE} "${RUN_DIR}/assim_advance_${n}.o"*              "${OUTPUT_DIR}/${datea}/logs/." 2>/dev/null || true
             ${MOVE} "WRFOUT/wrf.out_${gdatef[0]}_${gdatef[1]}_${n}" "${OUTPUT_DIR}/${datea}/logs/." 2>/dev/null || true
             ${REMOVE} "start_member_${n}" "done_member_${n}"
             if [[ -e "assim_advance_mem${n}.sh" ]]; then
                ${REMOVE} "assim_advance_mem${n}.sh"
             fi
             pert="$(cat "${RUN_DIR}/advance_temp${n}/mem${n}_pert_bank_num" 2>/dev/null || true)"
             echo "Member $n uses perturbation bank ensemble member $pert" >>  "${OUTPUT_DIR}/${datea}/pert_bank_members.txt"
          fi

          ${MOVE} "WRFIN/wrfinput_d${dchar}_${n}.gz"              "${OUTPUT_DIR}/${datea}/WRFIN/."
          ${MOVE} "${RUN_DIR}/prior_d${dchar}.${ensstring}"       "${OUTPUT_DIR}/${datea}/PRIORS/."
          ${REMOVE} "filter_restart_d${dchar}.${ensstring}"

         (( dn++ ))
      done #loop through domains

      (( n++ ))
   done

   # ---------------------------------------------------------------------------
   # All ensemble members should have advanced by now.
   # ---------------------------------------------------------------------------

   if [[ -e obs_prep.log ]]; then
      ${REMOVE} obs_prep.log
   fi

      #  Clean everything up and finish

      #  Move DART-specific data to storage directory
      ${COPY} input.nml "${OUTPUT_DIR}/${datea}/."
      ${MOVE} "${RUN_DIR}/dart_log.out" "${RUN_DIR}/dart_log.nml" "${RUN_DIR}/"*.log "${OUTPUT_DIR}/${datea}/logs/." 2>/dev/null || true

      #  Remove temporary files from both the run directory and old storage directories
      ${REMOVE} "${OUTPUT_DIR}/${datep}/wrfinput_d"*"_mean" "${RUN_DIR}/wrfinput_d"* "${RUN_DIR}/WRF" 2>/dev/null || true

      #  Prep data for archive
      cd "${OUTPUT_DIR}/${datea}"
      gzip -f wrfinput_d*_"${gdate[0]}"_"${gdate[1]}"_mean wrfinput_d*_"${gdatef[0]}"_"${gdatef[1]}"_mean wrfbdy_d*_mean
      tar -cvf retro.tar obs_seq.out wrfin*.gz wrfbdy_d*.gz
      tar -rvf dart_data.tar obs_seq.out obs_seq.final wrfinput_d*.gz wrfbdy_d*.gz \
                            Inflation_input/* logs/* input.nml
      ${REMOVE} wrfinput_d*_"${gdate[0]}"_"${gdate[1]}"_mean.gz wrfbdy_d*.gz
      gunzip -f wrfinput_d*_"${gdatef[0]}"_"${gdatef[1]}"_mean.gz

      cd "${RUN_DIR}"
      ${MOVE} "${RUN_DIR}/assim"*".o"*            "${OUTPUT_DIR}/${datea}/logs/." 2>/dev/null || true
      ${REMOVE} "${RUN_DIR}/input_priorinf_"* 2>/dev/null || true
      ${REMOVE} "${RUN_DIR}/static_data"* 2>/dev/null || true
      touch prev_cycle_done
      touch "${RUN_DIR}/cycle_finished_${datea}"
      if [[ -e "cycle_started_${datea}" ]]; then
         rm -f "${RUN_DIR}/cycle_started_${datea}"
      fi

      # If doing a reanalysis, increment the time if not done.  Otherwise, let the script exit
      if [[ "${restore}" == "1" ]]; then
         if [[ "${datea}" == "${datefnl}" ]]; then
            echo "Reached the final date "
            echo "Script exiting normally"
            exit 0
         fi
         datea="$(echo "${datea} ${ASSIM_INT_HOURS}" | "${DART_DIR}/models/wrf/work/advance_time")"
      else
         echo "Script exiting normally cycle ${datea}"
         exit 0
      fi

done

