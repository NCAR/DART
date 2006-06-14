#!/bin/csh
#
# Data Assimilation Research Testbed -- DART
# Copyright 2004, Data Assimilation Initiative, University Corporation for Atmospheric Research
# Licensed under the GPL -- www.gpl.org/licenses/gpl.html

#set echo

set queue = regular
set project_number = 86850054

set exp = rad
set exp_dir = jan_rad8

set month = 01
set chunk1 = 1
set chunkn = 31

set seconds = 0
set days = 146827  

set restart_in_file_name = filter_ics

set script_dir = /ptmp/$user/work

# -------------------------------------------------------------------------
# Have an overall outer loop for the days
#  can this loop divided into several scripts and these
#  scripts submitted at the same time but run sequentially?
# thanks alot.
#---------------------------
set i = $chunk1
while($i <= $chunkn)

   if($i < 10) then
      set day = 0${i}
   else
      set day = ${i}
   endif
   set job_step = step_${day}
   set script_filename = ${exp}_${job_step}.lsf.csh

   rm -f ${script_filename}

cat > ${script_filename} << EOF_HEAD
#### LSF options for BSUB
### -J      job name    (master script job.csh presumes filter.xxxx.log)
### -o      output listing filename 
### -P      account number
### -q      queue
### -n      number of tasks (processors)
### -x      exclusive use of node
### -R "span[ptile=(num procs you want on each node)]"
#
#BSUB -J ${exp}_${job_step}
#BSUB -W 06:00                          # run time limit
EOF_HEAD

   if( ${i} > ${chunk1}) then
      @ im1 = $i - 1
      if($im1 < 10) then
         set pre_job_step = step_0${im1}
      else
         set pre_job_step = step_${im1}
      endif
   
cat >> ${script_filename} << EOF_DEPENDENCY
#BSUB -w done(${exp}_${pre_job_step})                   # dependecy
EOF_DEPENDENCY
   endif

cat >> ${script_filename} << EOF_TAIL
#BSUB -o ${exp}_${job_step}.log
#BSUB -P ${project_number}
#BSUB -q ${queue}
#BSUB -B
#BSUB -N

   echo ' '
   echo starting iteration ${day}
   echo ' '

   cd ${script_dir}

   set work_dir = /ptmp/$user/${exp_dir}/${exp}/${month}_${day}

   if(!-d \${work_dir}) then
     mkdir -p \${work_dir}
   endif

   cd \${work_dir}

   ln -sf  /ptmp/hliu/WRF_bdy_2003jan   WRF
   ln -sf  namelist.input_50km   namelist.input

# Prepare observation file.
# for excess

   ln -sf /home/lightning/hliu/DART3/gps/obseq_rad/allexc/obs_seq2003${month}${day} obs_seq.out

# Prepare the input.nml

   rm -f script.sed
   echo "s/my_days/${days}/" > script.sed
   echo "s/my_seconds/${seconds}/" >> script.sed
   echo "s/my_restart_in/${restart_in_file_name}/" >> script.sed

#  change here -----------------------------------------

   sed -f script.sed ${script_dir}/input.nml.template.radonly > input.nml

#  change here -----------------------------------------

# Run filter

   ${script_dir}/filter >& fff.${day}

# Move the netcdf files to an output directory
#  mv -v  Prior_Diag.nc Posterior_Diag.nc obs_seq.final ../${exp_dir}/${exp}/${month}_${day}
   mv input.nml run.input.nml
   cp ${exp}_${job_step}.log .

# Move along to next iteration
   echo ' '
   echo ending iteration ${day}
   echo ' '
EOF_TAIL

   set restart_in_file_name = filter_restart

   echo bsub < ${script_filename}

   @ i ++
   @ days ++
end

