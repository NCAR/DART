#!/bin/bash

ensemble_projectname=Default_DART_Ensemble
#output_main_projectname=Test_Posterior_CovInf_All_Observe
#output_main_projectname=No_DA_Reference_All_Observe
output_main_projectname=Moha_Prior_Reference_Sparse
ensemble_directory=/Users/xinguan/DATA
obs_repository=/Users/xinguan/DATA/obs_repository/L96
n_ensemble_members=10
is_regular_obs=1                         # 0 = No, 1 = Yes
is_regular_obs_2=1                       # 0 = No, 1 = Yes
is_irregular_obs=0                       # 0 = No, 1 = Yes
is_random_obs=0                          # 0 = No, 1 = Yes | Note: not fully random, consistent in time
                                         #                 |     : for fully random obs, generate manually, then set is_full_random_input=1
n_tasks=2
is_clean_output_directory=1              # 0 = No, 1 = Yes

is_perfect_model_obs=1                   # 0 = No, 1 = Yes
is_generate_obs_input=1                  # 0 = No, 1 = Yes
is_run_filter=1                          # 0 = No, 1 = Yes
is_save_obs_file=1                       # 0 = No, 1 = Yes
obs_save_filepath_header=obs_seq_sparse
is_grab_obs_file=1                       # 0 = No, 1 = Yes
obs_grab_filepath_header=obs_seq_sparse
is_full_random_input=0                   # 0 = No, 1 = Yes | REQUIRES is_generate_obs_input=0

prior_inflation_type=5                   # 0 = No Inflation, 5 = Moha Inflation, 6 = My Inflation
posterior_inflation_type=0               # 0 = No Inflation, 5 = Moha Inflation, 6 = My Inflation
initial_prior_inflation=1.00
initial_posterior_inflation=1.00
inflation_increment_prior=0.002
inflation_increment_posterior=0.002
experiment_count=10
is_tuning=0                              # 0 = False, 1 = True

###########################
#### End of User Input ####
###########################

for i_experiment in $(seq 1 1 ${experiment_count}); do
	output_projectname="${output_main_projectname}_${i_experiment}"
	if [ ${is_tuning} -eq 0 ]; then
		obs_save_filepath="${obs_save_filepath_header}_${i_experiment}.out"
		obs_grab_filepath="${obs_grab_filepath_header}_${i_experiment}.out"
		printf -v pmo_file_index '%05d' ${i_experiment}
    else
		obs_save_filepath="${obs_save_filepath_header}_tuning.out"
		obs_grab_filepath="${obs_grab_filepath_header}_tuning.out"
		printf -v pmo_file_index '%05d' 1
		if [ ${i_experiment} -ne 1 ]; then
			is_generate_obs_input=0
			is_grab_obs_file=1
			is_save_obs_file=1
		fi
		initial_prior_inflation=$(echo "scale=10; ${initial_prior_inflation}+${inflation_increment_prior}" | bc)
		initial_posterior_inflation=$(echo "scale=10; ${initial_posterior_inflation}+${inflation_increment_posterior}" | bc)

		echo "Prior Inflation of ${initial_prior_inflation}"
		echo "Posterior Inflation of ${initial_posterior_inflation}"
	fi

	if [ ${is_clean_output_directory} -eq 1 ]; then
		if [ -d ${ensemble_directory}/${output_projectname} ]; then
			rm ${ensemble_directory}/${output_projectname}/*.nc
		fi
	fi

	cp runtime_create_input.sh.template runtime_create_input.sh
	sed -i .temp -e "s@PRIOR_INFLATION_TYPE@${prior_inflation_type}@"               runtime_create_input.sh
	sed -i .temp -e "s@POSTERIOR_INFLATION_TYPE@${posterior_inflation_type}@"       runtime_create_input.sh
	sed -i .temp -e "s@INITIAL_PRIOR_INFLATION@${initial_prior_inflation}@"         runtime_create_input.sh
	sed -i .temp -e "s@INITIAL_POSTERIOR_INFLATION@${initial_posterior_inflation}@" runtime_create_input.sh
	sed -i .temp -e "s@PMO_FILE_INDEX@${pmo_file_index}@"                           runtime_create_input.sh
	sh runtime_create_input.sh


	if [ ${is_perfect_model_obs} -eq 1 ]; then 
    	if [ ${is_generate_obs_input} -eq 1 ]; then
			num_obs_types=$((is_regular_obs+is_regular_obs_2+is_irregular_obs+is_random_obs))
			if [ ${num_obs_types} -eq 0 ]; then
				echo "No obs given, run cannot occur"
				exit 1
			fi
			if [ ${is_regular_obs} -eq 1 ]; then
				sh generate_create_obs_sequence_input_regular.sh
				./create_obs_sequence < COS_input.txt
			fi
			if [ ${is_regular_obs_2} -eq 1 ]; then
				sh generate_create_obs_sequence_input_regular_2.sh
				./create_obs_sequence < COS_input.txt
			fi
	        if [ ${is_irregular_obs} -eq 1 ]; then
			    sh generate_create_obs_sequence_input_irregular.sh
			    ./create_obs_sequence < COS_input.txt
	    	fi
			if [ ${is_random_obs} -eq 1 ]; then
				python generate_COS_input_random.py
				./create_obs_sequence < COS_input.txt
			fi

			if [ ${num_obs_types} -gt 1 ]; then
            	if [ -f ObsSeqFileList.txt ]; then
			    	rm ObsSeqFileList.txt
		  		fi
				if [ ${is_regular_obs} -eq 1 ]; then
			        sh generate_create_fixed_network_seq.sh 0
			        ./create_fixed_network_seq < CFNS_input.txt
				fi
				if [ ${is_regular_obs_2} -eq 1 ]; then
			        sh generate_create_fixed_network_seq.sh 6
			        ./create_fixed_network_seq < CFNS_input.txt
				fi
				if [ ${is_irregular_obs} -eq 1 ]; then
			        sh generate_create_fixed_network_seq.sh 1
			        ./create_fixed_network_seq < CFNS_input.txt
				fi
				if [ ${is_random_obs} -eq 1 ]; then
					sh generate_create_fixed_network_seq.sh 4
					./create_fixed_network_seq < CFNS_input.txt
				fi

			    ./obs_sequence_tool

		    elif [ ${is_regular_obs} -eq 1 ]; then
			    sh generate_create_fixed_network_seq.sh 2
		    	./create_fixed_network_seq < CFNS_input.txt
        	elif [ ${is_irregular_obs} -eq 1 ]; then
	        	sh generate_create_fixed_network_seq.sh 3
	        	./create_fixed_network_seq < CFNS_input.txt
		    elif [ ${is_random_obs} -eq 1 ]; then
				sh generate_create_fixed_network_seq.sh 5
        	    ./create_fixed_network_seq < CFNS_input.txt
        	else
				echo "CANNOT PARSE OBS TYPE"
				exit 1
	        fi
		elif [ ${is_full_random_input} -eq 1 ]; then
			python generate_full_random_obs_in.py
	    fi

    	mpirun -n ${n_tasks} ./perfect_model_obs

	elif [ ${is_grab_obs_file} -eq 1 ]; then
		mv ${obs_repository}/${obs_grab_filepath} ./obs_seq.out
	fi

	if [ ${is_run_filter} -eq 1 ]; then
        sh generate_filter_file_list.sh ${ensemble_directory}/${ensemble_projectname} ${ensemble_directory}/${output_projectname} ${n_ensemble_members}   
		mpirun -n ${n_tasks} ./filter
	fi

	if [ $? -eq 0 ]; then
		for thisoutputfile in analysis.nc forecast.nc postassim.nc preassim.nc true_state.nc; do
        	if [ -f ${thisoutputfile} ]; then
	        	mv ${thisoutputfile} ${ensemble_directory}/${output_projectname}
	    	fi
    	done
	fi

	if [ ${is_save_obs_file} -eq 1 ]; then
    	mv obs_seq.out ${obs_repository}/${obs_save_filepath}
	fi

	if [ -f obs_seq.final ]; then
		mv obs_seq.final obs_seq_${output_projectname}.final
	fi
done
