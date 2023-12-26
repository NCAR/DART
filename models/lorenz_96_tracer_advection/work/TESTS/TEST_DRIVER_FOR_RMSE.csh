# Script to do bit wise reproducibility testing for QCEFF filters applied
# to the L96 tracer model.

# Change to work directory (one level up)
cd ..

# Compile without mpi
quickbuild.sh nompi

# Create a single step obs_sequenc
create_obs_sequence < TESTS/create_obs_sequence_input

# Generate the 1000 timestep obs_seq.in file
create_fixed_network_seq < TESTS/create_fixed_network_seq_in

# List of CSV files
set csv_file_list = ('all_bnrhf_qceff_table.csv' 'neg_qceff_table.csv' 'TESTS\/one_below_qceff_table.csv')
@ low_csv_file  = 1
@ high_csv_file = 3

# Appropriate model namelist settings for the positive_tracer and bounded_above
set positive_tracer = ('.true.' '.false.' '.false.')
set bounded_above = ('.false.' '.false.' '.true.')

@ csv_file_index = $low_csv_file
while($csv_file_index <= $high_csv_file)
   echo $csv_file_index
   set csv_file = $csv_file_list[$csv_file_index]
   echo $csv_file

   # Set the values of the tracer params
   set POSITIVE_TRACER_VAL = $positive_tracer[$csv_file_index]
   set BOUNDED_ABOVE_VAL = $bounded_above[$csv_file_index]

   # The list  of ensemble sizes to test
   set ens_size_list = (6 10 15 20 25 30 40 50 60 70 80 90 100 110 120 130 140 150 160 170 180 240 320 640)
   @ low_ens_size = 1
   @ high_ens_size = 24 

   # Set up input.nml to do the initial perfect_model_obs run
   cp TESTS/TEST_BASE_INPUT.NML input.nml

   vi input.nml << HERE
:1,\$s/T_QCEFF_TABLE_FILENAME/$csv_file/
:1,\$s/T_READ_INPUT_STATE_FROM_FILE/.false./
:1,\$s/T_POST_INF_FLAVOR/5/
:1,\$s/T_POSITIVE_TRACER/$POSITIVE_TRACER_VAL/
:1,\$s/T_BOUNDED_ABOVE_IS_ONE/$BOUNDED_ABOVE_VAL/
:wq
HERE
   perfect_model_obs

   # Do the next perfect_model iteration and the filter
   cp perfect_output.nc perfect_input.nc

   cp TESTS/TEST_BASE_INPUT.NML input.nml
   vi input.nml << HERE
:1,\$s/T_QCEFF_TABLE_FILENAME/$csv_file/
:1,\$s/T_READ_INPUT_STATE_FROM_FILE/.true./
:1,\$s/T_POST_INF_FLAVOR/5/
:1,\$s/T_POSITIVE_TRACER/$POSITIVE_TRACER_VAL/
:1,\$s/T_BOUNDED_ABOVE_IS_ONE/$BOUNDED_ABOVE_VAL/
:wq
HERE

   perfect_model_obs

   # Loop over a range of ensemble sizes for filter
   @ ens_size_index = $low_ens_size
   while($ens_size_index <= $high_ens_size)
      # Get the ensemble size from the table
      set ens_size = $ens_size_list[$ens_size_index]

      cp TESTS/TEST_BASE_INPUT.NML input.nml
      vi input.nml << HERE
:1,\$s/T_QCEFF_TABLE_FILENAME/$csv_file/
:1,\$s/T_READ_INPUT_STATE_FROM_FILE/.true./
:1,\$s/T_FILTER_INPUT/perfect_input.nc/
:1,\$s/T_ENS_SIZE/$ens_size/
:1,\$s/T_PERTURB_FROM_SINGLE_INSTANCE/.true./
:1,\$s/T_POST_INF_FLAVOR/5/
:1,\$s/T_POSITIVE_TRACER/$POSITIVE_TRACER_VAL/
:1,\$s/T_BOUNDED_ABOVE_IS_ONE/$BOUNDED_ABOVE_VAL/
:wq
HERE

      filter

      # Generate matlab rmse time series plots
      matlab21a -nodesktop -nosplash -nodisplay << HERE
plot_total_err

figure(1)
print -dpng fig1.png
figure(2)
print -dpng fig2.png
HERE

      # Move the files to the test directory
      mv fig1.png TESTS/fig1-$csv_file_index-$ens_size.png
      mv fig2.png TESTS/fig2-$csv_file_index-$ens_size.png

      @ ens_size_index ++
   end
   @ csv_file_index ++

end

