# Script to do bit wise reproducibility testing for QCEFF filters applied
# to the L96 tracer model.

rm test_output

# Change to work directory (one level up)
cd ..

# Compile with mpi
quickbuild.sh

# Create a single step obs_sequenc
create_obs_sequence < TESTS/create_obs_sequence_input

# Generate the 1000 timestep obs_seq.in file
create_fixed_network_seq < TESTS/create_fixed_network_seq_in

cp TESTS/one_below_qceff_table.csv .

# List of CSV files
set csv_file_list = ('all_eakf_qceff_table.csv' 'all_bnrhf_qceff_table.csv' 'one_below_qceff_table.csv' 'one_below_qceff_table.csv')
@ low_csv_file  = 1
@ high_csv_file = 4

# Appropriate model namelist settings for the positive_tracer and bounded_above
set positive_tracer = ('.true.' '.true.' '.false.' '.false.')
set bounded_above = ('.false.' '.false.' '.true.' '.true.')
set post_inf_flavor = (0 0 0 5)

@ csv_file_index = $low_csv_file
while($csv_file_index <= $high_csv_file)
   set csv_file = $csv_file_list[$csv_file_index]

   # Set the values of the tracer params
   set POSITIVE_TRACER_VAL = $positive_tracer[$csv_file_index]
   set BOUNDED_ABOVE_VAL = $bounded_above[$csv_file_index]
   set POST_INF_FLAVOR_VAL = $post_inf_flavor[$csv_file_index]

   echo $csv_file post_inf_flavor is $POST_INF_FLAVOR_VAL >> TESTS/test_output

   # The list  of ensemble sizes to test
   set ens_size_list = (12 17 20 33 40 57 80 111 160)
   @ low_ens_size = 1
   @ high_ens_size = 9

   # The list of PE counts to test
   set pe_count_list = (1 2 3 7)
   @ low_pe_count = 1
   @ high_pe_count = 4

   # Set up input.nml to do the initial perfect_model_obs run
   cp TESTS/TEST_BASE_INPUT.NML input.nml

   vi input.nml << HERE
:1,\$s/T_QCEFF_TABLE_FILENAME/$csv_file/
:1,\$s/T_READ_INPUT_STATE_FROM_FILE/.false./
:1,\$s/T_POST_INF_FLAVOR/$POST_INF_FLAVOR_VAL/
:1,\$s/T_POSITIVE_TRACER/$POSITIVE_TRACER_VAL/
:1,\$s/T_BOUNDED_ABOVE_IS_ONE/$BOUNDED_ABOVE_VAL/
:wq
HERE
   perfect_model_obs

   # Do the next perfect_model iteration and the filter
   cp perfect_output.nc perfect_input.nc

   # Do ensemble size of 160 so that subsequent ICs with mpirun will have same ICs
   cp TESTS/TEST_BASE_INPUT.NML input.nml
   vi input.nml << HERE
:1,\$s/T_QCEFF_TABLE_FILENAME/$csv_file/
:1,\$s/T_READ_INPUT_STATE_FROM_FILE/.true./
:1,\$s/T_FILTER_INPUT/perfect_input.nc/
:1,\$s/T_ENS_SIZE/160/
:1,\$s/T_POST_INF_FLAVOR/$POST_INF_FLAVOR_VAL/
:1,\$s/T_PERTURB_FROM_SINGLE_INSTANCE/.true./
:1,\$s/T_POSITIVE_TRACER/$POSITIVE_TRACER_VAL/
:1,\$s/T_BOUNDED_ABOVE_IS_ONE/$BOUNDED_ABOVE_VAL/
:wq
HERE

   perfect_model_obs
   filter

   # Do  the next cycle of DA that has non-random filter ensemble ICs
   cp perfect_output.nc perfect_input.nc
   cp filter_output.nc filter_input.nc

   perfect_model_obs

   # Loop over a range of ensemble sizes
   @ ens_size_index = $low_ens_size
   while($ens_size_index <= $high_ens_size)
      # Get the ensemble size from the table
      set ens_size = $ens_size_list[$ens_size_index]

      # Loop over a range of PE counts
      @ pe_count_index = $low_pe_count
      while($pe_count_index <= $high_pe_count)
         set pe_count = $pe_count_list[$pe_count_index]

         cp TESTS/TEST_BASE_INPUT.NML input.nml
         vi input.nml << HERE
:1,\$s/T_QCEFF_TABLE_FILENAME/$csv_file/
:1,\$s/T_READ_INPUT_STATE_FROM_FILE/.true./
:1,\$s/T_FILTER_INPUT/filter_input.nc/
:1,\$s/T_ENS_SIZE/$ens_size/
:1,\$s/T_PERTURB_FROM_SINGLE_INSTANCE/.false./
:1,\$s/T_POST_INF_FLAVOR/$POST_INF_FLAVOR_VAL/
:1,\$s/T_POSITIVE_TRACER/$POSITIVE_TRACER_VAL/
:1,\$s/T_BOUNDED_ABOVE_IS_ONE/$BOUNDED_ABOVE_VAL/
:wq
HERE

         # Make sure there is no file around in case filter fails
         rm filter_output.nc

         mpirun --oversubscribe -np $pe_count filter
   
         echo -n 'ens_size = ' $ens_size, 'pes = ' $pe_count '  ' >> TESTS/test_output
         rm one_var_temp.nc
         ncrcat -d location,1,1 filter_output.nc one_var_temp.nc
         ncks -V -C -v state_variable_mean one_var_temp.nc | tail -3 | head -1 >> TESTS/test_output
         rm one_var_temp.nc
      
         @ pe_count_index ++
      end

      @ ens_size_index ++
   end
   @ csv_file_index ++

end

# Compare to the base line file

echo _________________
echo

cd TESTS
diff -q test_output BASELINE_OUTPUT

if ($status == 0) then
   echo TEST PASSED: test_output is the same as BASELINE_OUTPUT
else
   echo TEST FAILED: test_output differs from BASELINE_OUTPUT
endif
