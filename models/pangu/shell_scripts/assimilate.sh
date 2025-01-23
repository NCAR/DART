#!/bin/bash

#PBS  -N dart_pangu_test
#PBS  -A xxxxxxx
#PBS  -q casper
#PBS  -l select=1:ncpus=1:mem=16G
#PBS  -l walltime=00:40:00
#PBS  -o assim.out
#PBS  -j oe

# to run: qsub ./assimilate.sh
# this will cause memory issue: ./assimilate.sh >& assim.out &


####### begin user defined section ######
num_instances=3
num_cycles=1
start_datetime="2024-01-20-12"
output_dir="./era5_data"   
obs_dir='/glade/work/chennuo/code/cice_old/DART/observations/obs_converters/NCEP/prep_bufr/data/20240120/prepout'
filter_script="./test_dart.csh"

# modify the virtual environment activation command (e.g. conda/venv/mamba)
export activate='/glade/u/apps/opt/conda/bin/activate'
export deactivate='conda deactivate'
####### end user defined section ######

old_date=$start_datetime
source ${activate} pgdart

####### modify the ens_size in input.nml 
ex input.nml <<ex_end
g;ens_size ;s;= .*;= ${num_instances};
g;num_output_state_members ;s;= .*;= ${num_instances};
# g;num_output_obs_members ;s;= .*;= ${num_instances};
wq
ex_end

####### modify input and output file list for filter
output_file="filter_input_list.txt" 
rm $output_file
for ((i=1; i<=num_instances; i++)); do
    printf "pginput_%04d.nc\n" $i >> $output_file
done

output_file="filter_output_list.txt"
rm $output_file
for ((i=1; i<=num_instances; i++)); do
    printf "pgoutput_%04d.nc\n" $i >> $output_file
done

start_time=$(date +%s)
######### begin cycling loop

for ((i=1; i<=num_cycles; i++)); do 
    
    echo "------ STARTING CYCLE: $i -------"
    echo $old_date
    
    echo "------ BEGIN -----"
    echo "------ RUN PANGU WEATHER-----"
    ######### run the pangu model #########
    ### TODO: either run the ensemble members parallelly or run them on GPU
    ######### input: $output_dir/input_surface_000x.postassim-yyyy-mm-dd-hh.npy
    ######### input: $output_dir/input_upper_000x.postassim-yyyy-mm-dd-hh.npy
    ######### output: $output_dir/output_surface_000x.forecast-yyyy-mm-dd-hh+06.npy
    ######### output: $output_dir/output_upper_000x.forecast-yyyy-mm-dd-hh+06.npy
    
    ######### change the model onnx file location in interface_cpu/gpu.py
    inst=1
    while (( inst <= num_instances ))
    do
    	inst_string=$(printf "_%04d" $inst)
    	python inference_cpu.py $old_date $inst_string $output_dir
     
    	let inst++
    done
    echo "------ RUN PANGU WEATHER-----"
    echo "------ DONE -----"
    
    
    ######### add time increment
    date_string="${old_date%-*} ${old_date##*-}"
    formatted_date=$(date -d "$date_string + 6 hours" "+%Y-%m-%d %H:%M:%S")
    new_date="${formatted_date:0:10}-${formatted_date:11:2}"
    new_date_short="${new_date//-/}"

    ln -sf $obs_dir/obs_seq$new_date_short obs_seq.out
    
    echo "------ BEGIN -----"
    echo "------ CONVERT PANGU NPY INTO NETCDF-----"
    ######### convert pangu model output to netcdf file #########
    ######### input: $output_dir/output_surface_000x.forecast-yyyy-mm-dd-hh.npy
    ######### input: $output_dir/output_upper_000x.forecast-yyyy-mm-dd-hh.npy
    ######### output: pginput_000x.nc
    
    inst=1
    while (( inst <= num_instances ))
    do
    	inst_string=$(printf "_%04d" $inst)
    	python convert_pgout_to_nc.py $new_date $inst_string $output_dir
     
    	let inst++
    done
    
    echo "------ CONVERT PANGU NPY INTO NETCDF-----"
    echo "------ DONE -----"
    
    
    echo "------ BEGIN -----"
    echo "------ DART FILTER -----"
    ######### run ./filter  #########
    ######### >TODO: any place in the shell script requires the fill_inflation?
    ######### requires input.nml, sampling_error_correction_table.nc
    ######### new obs_seq.out already linked
    
    job_id=$(qsub ${filter_script} | cut -d '.' -f 1)
    
    ######### this is run ./filter for cycling?  #########
    ######### waiting 30 seconds for job to get into the quque
    sleep 30
    ######### Wait for the job to finish
    while qstat | grep -q $job_id
        do
            sleep 30
        done
    
    echo "------ DART FILTER -----"
    echo "------ DONE -----"
    
    
    ######### convert pgoutput_00xx.nc into npy files that pangu model takes in
    # input: pgoutput_00xx.nc
    # output: data/input_surface_000x.postassim-yyyy-mm-dd-hh.npy
    # output: data/input_upper_000x.postassim-yyyy-mm-dd-hh.npy
    
    echo "------ BEGIN -----"
    echo "------ CONVERT pgoutput_000x INTO NPY -----"
    
    inst=1
    while (( inst <= num_instances ))
    do
    	inst_string=$(printf "_%04d" $inst)
       # !!!!!!!  IMPORTANT WARNING !!!!!
       # for time step one, use ensemble version convert_dartout_to_npy_mem.py, or modify the pgout
    	python convert_dartout_to_npy.py $new_date $inst_string $output_dir
     
    	let inst++
    done
    
    echo "------ CONVERT pgoutput_000x INTO NPY-----"
    echo "------ DONE -----"
    
    ######### rename the dart output_mean and output_sd
    mv output_mean.nc $output_dir/output_mean_$new_date_short.nc
    mv output_sd.nc $output_dir/output_sd_$new_date_short.nc
    mv obs_seq.final $output_dir/obs_seq.final_$new_date_short
    
    
    # repeat and cycle
    
    old_date=$new_date

done
######### end cycling loop

${deactivate}

end_time=$(date +%s)
runtime=$((end_time - start_time))
runtime_minutes=$(echo "scale=2; $runtime / 60" | bc)

echo "The script took $runtime_minutes minutes to run."
