#!/bin/bash

ensemble_projectname=Default_DART_Ensemble
output_projectname=Debugging
#output_projectname=Full_Obs_No_DA_Reference
ensemble_directory=/Users/xinguan/DATA
obs_repository=/Users/xinguan/DATA/obs_repository/L96
n_ensemble_members=10
is_regular_obs=1                         # 0 = No, 1 = Yes
is_irregular_obs=0                       # 0 = No, 1 = Yes
is_random_obs=0                          # 0 = No, 1 = Yes | Note: not fully random, consistent in time
                                         #                 |     : for fully random obs, generate manually, then set is_full_random_input=1
n_tasks=4

is_perfect_model_obs=1                   # 0 = No, 1 = Yes
is_generate_obs_input=0                  # 0 = No, 1 = Yes
is_run_filter=1                          # 0 = No, 1 = Yes
is_save_obs_file=1                       # 0 = No, 1 = Yes
obs_save_filepath=obs_seq_dense_1.out
is_grab_obs_file=1                       # 0 = No, 1 = Yes
obs_grab_filepath=obs_seq_dense_1.out
is_full_random_input=0                   # 0 = No, 1 = Yes | REQUIRES is_generate_obs_input=0
is_create_inflation_initial_file=0       # 0 = No, 1 = Yes

###########################
#### End of User Input ####
###########################

sh create_input.sh

if [ ${is_create_inflation_initial_file} -eq 1 ]; then
    ./fill_inflation_restart
fi

if [ ${is_perfect_model_obs} -eq 1 ]; then 
    if [ ${is_generate_obs_input} -eq 1 ]; then

        if [ ${is_regular_obs} -eq 0 ]; then
	        if [ ${is_irregular_obs} -eq 0 ]; then
				if [ ${is_random_obs} -eq 0 ]; then
		            echo "No obs given, run cannot occur"
		            exit 1
				else
					python generate_COS_input_random.py
					./create_obs_sequence < COS_input.txt
				fi
	        else
		        sh generate_create_obs_sequence_input_irregular.sh
		        ./create_obs_sequence < COS_input.txt
				if [ ${is_random_obs} -eq 1 ]; then
					python generate_COS_input_random.py
					./create_obs_sequence < COS_input.txt
				fi
	        fi
        else
	        sh generate_create_obs_sequence_input_regular.sh
	        ./create_obs_sequence < COS_input.txt
	        if [ ${is_irregular_obs} -eq 1 ]; then
		        sh generate_create_obs_sequence_input_irregular.sh
		        ./create_obs_sequence < COS_input.txt
	        fi
			if [ ${is_random_obs} -eq 1 ]; then
				python generate_COS_input_random.py
				./create_obs_sequence < COS_input.txt
			fi
        fi

		num_obs_types=$((is_regular_obs + is_irregular_obs + is_random_obs))
		if [ ${num_obs_types} -gt 1 ]; then
            if [ -f ObsSeqFileList.txt ]; then
			    rm ObsSeqFileList.txt
		    fi

			if [ ${is_regular_obs} -eq 1 ]; then
		        sh generate_create_fixed_network_seq.sh 0
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

