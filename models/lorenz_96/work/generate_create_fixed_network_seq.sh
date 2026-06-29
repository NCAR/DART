#!/bin/bash

regularity=$1                # 0 = Regular, 1 = Irregular, 2 = Agnostic Regular,
                             # 3 = Agnostic Irregular, 4 = Random, 5 = Agnostic Random, 6 = Regular 2

obs_template_basename=set_def
output_file_basename=obs_seq
obs_frequency_type=1         # 1 = Regular, 2 = Irregular

# regular obs
n_obs_in_sequence=3000
initial_time=(0 0)           # days, seconds
obs_interval=(1 0)

# regular obs 2
n_obs_in_sequence_2=3000
initial_time_2=(0 43200)
obs_interval_2=(1 0)

# irregular obs
max_obs_number=200000
n_obs_cycle_irreg=10000
is_gen_irreg_times=1         # 0 = No, 1 = Yes
#irreg_time_days=(0 0 0 1 2 3 3 3 4 5 6 6 6 7 8)
#irreg_time_seconds=(7200 21600 43200 0 0 7200 21600 43200 0 0 7200 21600 43200 0 0)

# random obs
max_random_obs_number=10000
random_obs_end_day=75

#########################
### End of User Input ###
#########################

if [ ${is_gen_irreg_times} -eq 1 ]; then
	irreg_time_days=()
	irreg_time_seconds=()
    pattern_counter=0
    cycle_counter=0
    day_counter=0
	while [ ${cycle_counter} -lt ${n_obs_cycle_irreg} ]; do
        if [ ${pattern_counter} -eq 2 ]; then
            pattern_counter=0
			cycle_counter=$((cycle_counter+4))
            for i_seconds in 0 21600 43200 64800; do
				irreg_time_seconds+=(${i_seconds})
				irreg_time_days+=(${day_counter})
	        done
			day_counter=$((day_counter+1))
		elif [ ${pattern_counter} -gt 2 ]; then
			echo "ERROR SETTING UP IRREGULAR OBSERVATION TIMES"
			exit 1
		else
			pattern_counter=$((pattern_counter+1))
			day_counter=$((day_counter+1))
		fi
	done
fi


if [ ${regularity} -eq 0 ]; then
	obs_template_name=${obs_template_basename}_regular.out
    obs_frequency_type=1
	output_file_name=${output_file_basename}_regular.in
	echo ${output_file_name} >> ObsSeqFileList.txt
elif [ ${regularity} -eq 6 ]; then
	obs_template_name=${obs_template_basename}_regular_2.out
    obs_frequency_type=1
	output_file_name=${output_file_basename}_regular_2.in
	echo ${output_file_name} >> ObsSeqFileList.txt
	initial_time=(${initial_time_2[0]} ${initial_time_2[1]})
elif [ ${regularity} -eq 1 ]; then
	obs_template_name=${obs_template_basename}_irregular.out
    obs_frequency_type=2          
    output_file_name=${output_file_basename}_irregular.in
	echo ${output_file_name} >> ObsSeqFileList.txt
elif [ ${regularity} -eq 2 ]; then
	obs_template_name=${obs_template_basename}_regular.out
    obs_frequency_type=1          
    output_file_name=${output_file_basename}.in
elif [ ${regularity} -eq 3 ]; then
	obs_template_name=${obs_template_basename}_irregular.out
    obs_frequency_type=2          
    output_file_name=${output_file_basename}.in
elif [ ${regularity} -eq 4 ]; then
	obs_template_name=${obs_template_basename}_random.out
	output_file_name=${output_file_basename}_random.in
	echo ${output_file_name} >> ObsSeqFileList.txt
	python generate_CFNS_input_random.py ${obs_template_name} ${output_file_name} ${max_random_obs_number} ${random_obs_end_day}
	exit $?
elif [ ${regularity} -eq 5 ]; then
	obs_template_name=${obs_template_basename}_random.out                                                
    output_file_name=${output_file_basename}_random.in
    python generate_CFNS_input_random.py ${obs_template_name} ${output_file_name} ${max_random_obs_number} ${random_obs_end_day}
    exit $?	
fi

echo ${obs_template_name} > CFNS_input.txt
echo ${obs_frequency_type} >> CFNS_input.txt
if [ ${obs_frequency_type} -eq 1 ]; then
    echo ${n_obs_in_sequence} >> CFNS_input.txt
    echo "${initial_time[0]} ${initial_time[1]}" >> CFNS_input.txt
    echo "${obs_interval[0]} ${obs_interval[1]}" >> CFNS_input.txt
    echo ${output_file_name} >> CFNS_input.txt
else
	echo ${max_obs_number} >> CFNS_input.txt
	irreg_obs_counter=0
	while [ ${irreg_obs_counter} -lt ${n_obs_cycle_irreg} ]; do
	    echo "${irreg_time_days[${irreg_obs_counter}]} ${irreg_time_seconds[${irreg_obs_counter}]}" >> CFNS_input.txt
		irreg_obs_counter=$((irreg_obs_counter+1))
    done
	echo "-1 0" >> CFNS_input.txt
    echo ${output_file_name} >> CFNS_input.txt
fi
