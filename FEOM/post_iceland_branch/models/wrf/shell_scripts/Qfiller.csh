#!/bin/csh
#
# Data Assimilation Research Testbed -- DART
# Copyright 2004-2006, Data Assimilation Research Section 
# University Corporation for Atmospheric Research
# Licensed under the GPL -- www.gpl.org/licenses/gpl.html
#
# <next three lines automatically updated by CVS, do not edit>
# $Id$
# $Source$
# $Name$
#
# DESCRIPTION:
#
# This script is designed to create a series of scripts that will be 
# submitted to a queue. The batch job directives in each of the scripts
# will ensure that the next script will start to execute only upon
# successful completion of the previous script.
#
# The intent here is that an experiment usually consists of looping over
# a series of observation sequence files. Each observation sequence file
# will be 'processed' by one script. Each script will have to perform a
# bit of work 'setting the table', etc. THIS script actually BUILDS the
# other scripts. 

set myname = $0
set REMOVE = 'rm -rf'
set   COPY = 'cp -p'
set   MOVE = 'mv -f'

# Create a clean temporary directory and go there
${REMOVE} $temp_dir
mkdir -p  $temp_dir
cd        $temp_dir


set queue = regular
set project_number = 86850054

set exp = rad
set exp_dir = test
mkdir -p ../${exp_dir}/${exp}

set month = 01
set chunk1 = 1
set chunkn = 3

set seconds = 0
set days = 146827  

set restart_in_file_name = filter_ics

#set script_dir = /ptmp/$user/work
setenv script_dir  `pwd`
echo $script_dir 

rm *.log *.err *.out

# ----------------------------------------
# Have an overall outer loop for the days
#-----------------------------------------
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
#!/bin/csh
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
#BSUB -x                                # exclusive use of node
#BSUB -n 1                              # sum of number of tasks
#BSUB -R "span[ptile=1]"                # number of processes per node
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
#BSUB -o ${exp}_${job_step}.out
#BSUB -e ${exp}_${job_step}.err
#BSUB -P ${project_number}
#BSUB -q ${queue}
#BSUB -B
#BSUB -N

set echo

   echo ' '
   echo starting iteration ${day}
   echo ' '

   cd ${script_dir}

rm -f  *.nc  
rm -rf /ptmp/hliu/wrfrun*

rm -f assim_model_state_ud*
rm -f assim_model_state_ic*

rm -f WRF namelist.input 


   set out_dir = /ptmp/$user/${exp_dir}/${exp}/${month}_${day}

   if(! -d \${out_dir}) then
     mkdir -p \${out_dir}
   endif

#  ln -s   /ptmp/hliu/WRF_bdy_2003jan   WRF
#  ln -sf  namelist.input_50km   namelist.input

   ln -s   /ptmp/hliu/WRF_bdy_200km   WRF
   ln -sf  namelist.input_200km   namelist.input

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
   mv -v  Prior_Diag.nc Posterior_Diag.nc obs_seq.final \${out_dir}/.
   mv input.nml \${out_dir}/run.input.nml
#  cp ${exp}_${job_step}.log .

# Move along to next iteration
   echo ' '
   echo ending iteration ${day}
   echo ' '

EOF_TAIL

   set restart_in_file_name = filter_restart

   bsub < ${script_filename}

   @ i ++
   @ days ++
end

