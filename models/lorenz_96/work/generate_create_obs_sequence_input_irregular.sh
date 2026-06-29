#!/bin/bash

n_obs=5
max_num_obs=10000
copies_of_data=0                                   # 0 to set just a definition
num_qc_values_per_qc_field=0
variable_name=RAW_STATE_VARIABLE
#variable_location=(0.1 0.3 0.5 0.7 0.9)           # between 0 and 1
is_gen_variable_location=1                         # 0 = No, 1 = Yes
time=(0 0)                                         # days seconds
variance=1
file_name=set_def_irregular.out

#########################
### End of User Input ###
#########################

if [ ${is_gen_variable_location} -eq 1 ]; then
    variable_location=()
	n_obs=0
	for i_var_loc in $(seq 0.0 0.025 0.975); do
        variable_location+=(${i_var_loc})
		n_obs=$((n_obs+1))
	done
fi

echo ${max_num_obs} > COS_input.txt
echo ${copies_of_data} >> COS_input.txt
echo ${num_qc_values_per_qc_field} >> COS_input.txt
this_obs_counter=0
while [ ${this_obs_counter} -lt ${n_obs} ]; do
    echo 1 >> COS_input.txt
    echo ${variable_name} >> COS_input.txt
    echo ${variable_location[${this_obs_counter}]} >> COS_input.txt
    echo "${time[0]} ${time[1]}" >> COS_input.txt
	echo ${variance} >> COS_input.txt
	this_obs_counter=$((this_obs_counter+1))
done
echo -1 >> COS_input.txt
echo ${file_name} >> COS_input.txt
