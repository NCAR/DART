#!/bin/csh
# Shell script to do repeated segment integrations
# BUT, don't know how to do MATLAB from script yet

set exp = Guam_test

# WARNING filter_restart_31 won't be copied to filter_ics; already done?

set output_dir = 01
set chunk1 = 1
set chunkn = 2

set save_freq = 2
set mod_save = 1

set input = input_

if (-d ${exp}) then
   echo ${exp} already exists
else
   mkdir ${exp}
   cp namelistin ${exp}/namelistin
endif

# Have an overall outer loop
set i = $chunk1
while($i <= $chunkn)
   echo ' '
   echo ' '
   echo starting iteration $i
   echo ' '
   echo ' '

# remove clminput files from previous experiment
   if (-e clminput_1.nc && $i == $chunk1) rm clminput_*.nc

   cp $input$i.nml input.nml

   if ($i == 1) then
      cp filter_restart_31 filter_ics
      csh ./init_advance_model.csh
   else if ($i == $chunk1) then
      @ j = $i - 1
      cp ${exp}/${output_dir}_$j/filter_restart filter_ics
      cp ${exp}/${output_dir}_$j/CLM/clminput*.nc .
      csh ./init_advance_model.csh
   else
      mv filter_restart filter_ics
   endif


# Run filter with async=2 in input.nml
   ./filter >&! filter_$i.out

# Move the netcdf files to an output directory
   mkdir ${exp}/${output_dir}_$i
   mkdir ${exp}/${output_dir}_$i/CLM
   mv Prior_Diag.nc Posterior_Diag.nc    ${exp}/${output_dir}_$i
   mv obs_seq.final                      ${exp}/${output_dir}_$i
   mv input.nml                          ${exp}/${output_dir}_$i
   mv cam_out_temp1                      ${exp}/${output_dir}_$i
#   cp obs_seq_jan$i.out                  ${exp}/${output_dir}_$i/obs_seq.out

# Copy the filter restart to start files
   cp filter_restart ${exp}/${output_dir}_$i
   cp clminput_*.nc ${exp}/${output_dir}_$i/CLM
   if ($i % $save_freq != $mod_save) then
      @ j = $i - 1
#      if (-e ${exp}/${output_dir}_$i/filter_restart)    \
#         rm ${exp}/${output_dir}_$j/filter_restart 
#      if (-e ${exp}/${output_dir}_$i/CLM/clminput_20.nc) \
#         rm ${exp}/${output_dir}_$j/CLM/clminput*.nc 
      if (-e ${exp}/${output_dir}_$i/CLM/clminput_20.nc) then
         cd ${exp}/${output_dir}_$j
         gzip -r CLM
         tar cf clminput.gz.tar CLM
         rm -rf CLM
         cd ../..
      else
         echo 'NO clminput.80;  ABORTING' 
         exit 1
      endif       
   endif

# Move along to next iteration
   echo ' '
   echo ending iteration $i
   echo ' '
   @ i++

end
