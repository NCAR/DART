#!/bin/bash

#############################
### Basic Run Information ###
#############################

ensemble_size=10

################
### Run Type ###
################

is_qceff=0              # 0 = No, 1 = Yes
is_perfect_model_obs=1  # 0 = No, 1 = Yes

##########################
### Run Specifications ###
##########################

# Run Logistics #
n_tasks_per_node=4
single_precision_output=1              # 0 = No, 1 = Yes
vector_read_piecewise=1                # 0 = No, 1 = Yes
regression_factor_type=1               # 1 = via sampling theory, and size
                                       # 2 = use archived time mean (exists for Lorenz_96)
                				       # 3 = use bgrid archived file (? available)
regression_factor_input_file=time_mean_reg
use_sampling_error_correction=0        # 0 = No, 1 = Yes
assimilation_period_seconds=30
parallel_mode=0                      # 0 = All calls are serial, model advanced by subroutine call
                                     # 2 = Filter calls in parallel, model advanced by shell script
				                     # 4 = Filter and NWP calls in parallel, model advanced by parallel job script
input_filelist=lorenz_input_file_list.txt
output_filelist=lorenz_output_file_list.txt

# Model #
number_of_locations=40                 # Also number of state variables (sans inflation variables)
lorenz_forcing=8.00
lorenz_timestep=0.05
real_equivalent_timestep_days=0
real_equivalent_timestep_seconds=3600

# Perfect Model Obs #
is_perfect_model_obs_from_file=1               # 0 = No, 1 = Yes
is_write_perfect_model_member_output_file=0    # 0 = No, 1 = Yes
perfect_model_file_input=/Users/xinguan/DATA/Synthetic_Truth/l96_mem00001.nc      # Run filepath that determines synthetic truth
perfect_model_file_output=pmo_output.nc    # Synthetic truth file (may have changed from input)
perfect_model_obs_input=obs_seq.in             # Describes location, inaccuracy, and other metadata of observations
perfect_model_obs_output=obs_seq.out           # Output observations
# To be continued ...

# QCEFF #
is_qceff_already_made=1       # 0 = No, 1 = Yes
qceff_table_name=
# qceff_params=() # On the fly generation of QCEFF table currently unavailable

# Filter #
model_advance_shell_script=./advance_model.csh
stages_to_write=(forecast preassim postassim analysis output) # input, forecast, preassim, postassim, analysis, output
                                                              # Note that "output" controls printing of posterior members
is_output_members=1                # 0 = No, 1 = Yes
n_output_state_members=${ensemble_size}
n_output_obs_members=${ensemble_size}
filter_obs_input=obs_seq.out
filter_obs_output=obs_seq.final

# Inflation | All inflation parameters have two inputs--key in as length 2 array--the first for the prior, the second for the posterior
inflation_type=(6 6) # 0=No inflation, 2=Spatial and Temporally varying Gaussian inflation
                                         # 3=Temporally varying Gaussian inflation, 4=Relaxation to Prior Spread
                                         # 5=Spatially and Temporally varying Gamma inflation
										 # 6=Covariance_Only_Inflation :D
get_inflation_from_restart_file=(0 0)    # 0=No, 1=Yes
initial_inflation=(1.1 1.1)          # While this is the inflation factor if not getting inflation from restart file
                                         # If inflation_type=4 (RTPS), the second entry is the alpha parameter
                                         # (alpha=0.3 means 30% prior spread, 70% posterior spread)
initial_inflation_std=(0.6 0.3)
inflation_lower_bound=(1.0 1.0)
inflation_upper_bound=(10000.0 10000.0)
inflation_std_lower_bound=(0.3 0.1)
inflation_std_max_step=(1.05 1.05)       # Only used for flavour=5, where std can change.
                                         # Denotes maximum (proportional) change in std with each step.
                                         # Value of 1 means std can not increase
inflation_damping=(1.0 1.0)
write_own_initial_inflation=(0 0)        # 0=No, 1=Yes # For generating your own initial inflation file
own_initial_inflation_mean=(1.03 1.08)
own_initial_inflation_std=(0.6 0.6)
own_initial_inflation_reference_state_file=\'reference_lorenz_member.nc\' # THIS ONE DOESN'T WORK --> NEED SOMETHING WITH TIME DIMENSION

# Localisation #
localisation_half_radius=0.2
adaptive_localisation_n_obs_threshold=-1      # -1 to turn off adaptive localisation

# Quality Control #
quality_control_grade_threshold=3.0   # Acceptable "grade" of observation. Grades are given by observation source
quality_control_std_threshold=-1      # Acceptable number of sigma away from prior mean an observation can be
use_special_outlier_function=0        # 0 = No, 1 = Yes | Programme then looks for failed_outlier() to run if Yes

# Obs Sequence Tool #
# Mostly used to combine obs_seq.in for use in perfect_model_obs #
is_incompatible_metadata=0                # 0 = No, 1 = Yes
new_metadata_name=observation           # for obs_sequence_tool metadata.
obs_seq_filelist=ObsSeqFileList.txt
obs_seq_tool_output_file=obs_seq.in
observation_days=(0 300)                    # Start, End
observation_seconds=(0 0)                 # Start, End

#########################
### END OF USER INPUT ###
#########################

cp input.nml.template input.nml

# Settling Filepathing and Basic Toggles
sed -i -e "s@INPUTFILES@${input_filelist}@" input.nml
sed -i -e "s@OUTPUTFILES@${output_filelist}@" input.nml
sed -i -e "s@PARALLELMODE@${parallel_mode}@g" input.nml
sed -i -e "s@N_OUTPUT_STATE@${n_output_state_members}@" input.nml
sed -i -e "s@ENSEMBLE_SIZE@${ensemble_size}@" input.nml
sed -i -e "s@PMO_OBS_IN@${perfect_model_obs_input}@" input.nml
sed -i -e "s@PMO_OBS_OUT@${perfect_model_obs_output}@" input.nml
sed -i -e "s@PMO_FILE_IN@${perfect_model_file_input}@" input.nml
sed -i -e "s@PMO_FILE_OUT@${perfect_model_file_output}@" input.nml
sed -i -e "s@N_OUTPUT_OBS@${n_output_obs_members}@" input.nml
sed -i -e "s@FILTER_OBS_IN@${filter_obs_input}@" input.nml
sed -i -e "s@FILTER_OBS_OUT@${filter_obs_output}@" input.nml
sed -i -e "s@OWNINITINFINPUT_FILES@${own_initial_inflation_reference_state_file}@" input.nml
sed -i -e "s@NTASKSPERNODE@${n_tasks_per_node}@" input.nml
sed -i -e "s@REGRESSIONFACTORTYPE@${regression_factor_type}@" input.nml
sed -i -e "s@REGRESSIONFACTORINPUTFILE@${regression_factor_input_file}@" input.nml
sed -i -e "s@MODELMODASSIMSECS@${assimilation_period_seconds}@" input.nml
sed -i -e "s@NUMLORENZLOC@${number_of_locations}@" input.nml
sed -i -e "s@LORENZF@${lorenz_forcing}@" input.nml
sed -i -e "s@LORENZDT@${lorenz_timestep}@" input.nml
sed -i -e "s@RWEQDTDAYS@${real_equivalent_timestep_days}@" input.nml
sed -i -e "s@RWEQDTSECS@${real_equivalent_timestep_seconds}@" input.nml
sed -i -e "s@ASSIMCUTOFF@${localisation_half_radius}@" input.nml
sed -i -e "s@ADAPLOCTHRESH@${adaptive_localisation_n_obs_threshold}@" input.nml
sed -i -e "s@QCGRADE@${quality_control_grade_threshold}@" input.nml
sed -i -e "s@QCSTD@${quality_control_std_threshold}@" input.nml
sed -i -e "s@OBSSEQFILELIST@${obs_seq_filelist}@" input.nml
sed -i -e "s@OBSSEQOUTPUTFILE@${obs_seq_tool_output_file}@" input.nml
#sed -i -e "s@@@" input.nml
#sed -i -e "s@@@" input.nml
#sed -i -e "s@@@" input.nml
#sed -i -e "s@@@" input.nml
#sed -i -e "s@@@" input.nml
#sed -i -e "s@@@" input.nml

# Settling Boolean Replacements
if [ ${is_perfect_model_obs_from_file} -eq 0 ]; then
	sed -i -e "s@BOOL_PMO_IN_FILE@.false.@" input.nml
else
	sed -i -e "s@BOOL_PMO_IN_FILE@.true.@" input.nml
fi

if [ ${is_write_perfect_model_member_output_file} -eq 0 ]; then
	sed -i -e "s@BOOL_PMO_OUT_FILE@.false.@" input.nml
else
	sed -i -e "s@BOOL_PMO_OUT_FILE@.true.@" input.nml
fi

if [ ${is_output_members} -eq 0 ]; then
	sed -i -e "s@BOOL_OUTPUT_MEMBERS@.false.@" input.nml
else
	sed -i -e "s@BOOL_OUTPUT_MEMBERS@.true.@" input.nml
fi

if [ ${single_precision_output} -eq 1 ]; then
    sed -i -e "s/SQUISHMYOUTPUT/.true./" input.nml
else
    sed -i -e "s/SQUISHMYOUTPUT/.false./" input.nml
fi

if [ ${use_sampling_error_correction} -eq 1 ]; then
    sed -i -e "s/SAMPLINGERRORCORRECTION/.true./" input.nml
else
    sed -i -e "s/SAMPLINGERRORCORRECTION/.false./" input.nml
fi

if [ ${use_special_outlier_function} -eq 1 ]; then
    sed -i -e "s/BOOLQCSPECIALCODE/.true./" input.nml
else
    sed -i -e "s/BOOLQCSPECIALCODE/.false./" input.nml
fi

if [ ${is_incompatible_metadata} -eq 1 ]; then
    sed -i -e "s/BOOLOBSSEQMETADATA/.true./" input.nml
    sed -i -e "s/OBSSEQOBSMETANAME/${new_metadata_name}/" input.nml
else
    sed -i -e "s/BOOLOBSSEQMETADATA/.false./" input.nml
fi

#if [  -eq 0 ]; then
#	sed -i -e "s@@.false.@" input.nml
#else
#	sed -i -e "s@@.true.@" input.nml
#fi

# Settling Looped Toggles
for inflation_index in 0 1; do
    sed -i -e "s/INFFLAVOUR${inflation_index}/${inflation_type[${inflation_index}]}/" input.nml
    sed -i -e "s/INITIALINFLATION${inflation_index}/${initial_inflation[${inflation_index}]}/" input.nml
    sed -i -e "s/INFLATIONLOWERBOUND${inflation_index}/${inflation_lower_bound[${inflation_index}]}/" input.nml
    sed -i -e "s/INFLATIONUPPERBOUND${inflation_index}/${inflation_upper_bound[${inflation_index}]}/" input.nml
    sed -i -e "s/INITIALINFLATIONSTD${inflation_index}/${initial_inflation_std[${inflation_index}]}/" input.nml
    sed -i -e "s/INFLATIONSTDLOWERBOUND${inflation_index}/${inflation_std_lower_bound[${inflation_index}]}/" input.nml
    sed -i -e "s/INFLATIONSTDMAXSTEP${inflation_index}/${inflation_std_max_step[${inflation_index}]}/" input.nml
    sed -i -e "s/INFLATIONDAMPING${inflation_index}/${inflation_damping[${inflation_index}]}/" input.nml
    sed -i -e "s/OWNINITINFMEAN${inflation_index}/${own_initial_inflation_mean[${inflation_index}]}/" input.nml
    sed -i -e "s/OWNINITINFSTD${inflation_index}/${own_initial_inflation_std[${inflation_index}]}/" input.nml
    if [ ${get_inflation_from_restart_file[${inflation_index}]} -eq 1 ]; then
        sed -i -e "s/BOOLINFINITRESTART${inflation_index}/.true./" input.nml
		if [ ${inflation_type[${inflation_index}]} -eq 5 ]; then
            sed -i -e "s/BOOLINFSTDINITRESTART${inflation_index}/.true./" input.nml
		else
            sed -i -e "s/BOOLINFSTDINITRESTART${inflation_index}/.false./" input.nml
		fi
    else
        sed -i -e "s/BOOLINFINITRESTART${inflation_index}/.false./" input.nml
        sed -i -e "s/BOOLINFSTDINITRESTART${inflation_index}/.false./" input.nml
    fi
    if [ ${write_own_initial_inflation[${inflation_index}]} -eq 1 ]; then
        sed -i -e "s/BOOLINFRESTARTWRITE${inflation_index}/.true./" input.nml
    else
        sed -i -e "s/BOOLINFRESTARTWRITE${inflation_index}/.false./" input.nml
    fi
done

StageWriteCounter=${#stages_to_write[@]}
while [ ${StageWriteCounter} -gt 0 ]; do
	StageWriteCounter=$((StageWriteCounter-1))
	if [ ${StageWriteCounter} -eq 0 ]; then
		sed -i -e "s@STAGES_TO_WRITE@\"${stages_to_write[${StageWriteCounter}]}\"@" input.nml
	else
		sed -i -e "s@STAGES_TO_WRITE@STAGES_TO_WRITE, \"${stages_to_write[${StageWriteCounter}]}\"@" input.nml
	fi
done

# Settling Extra Arguments
if [ ${parallel_mode} -eq 2 ]; then
	sed -i -e "s@ASYNCADDON@\n   adv_ens_command              = \"${model_advance_shell_script}\"@g" input.nml
else
	sed -i -e "s@ASYNCADDON@@g" input.nml
fi

# Settling QCEFF
if [ ${is_qceff} -eq 1 ]; then
	if [ ${is_qceff_already_made} -eq 1 ]; then
		sed -i -e "s@QCEFF_TABLE_NAME@${qceff_table_name}@" input.nml
	else
		echo 'On the fly generation of QCEFF table currently unavailable. Turning off QCEFF'
		sed -i -e "s@QCEFF_TABLE_NAME@@" input.nml
	fi
else
	sed -i -e "s@QCEFF_TABLE_NAME@@" input.nml
fi

# Settling Obs Sequence Tool
for obs_seq_tool_index in 0 1; do
    sed -i -e "s@OBSSEQDAYS${obs_seq_tool_index}@${observation_days[${obs_seq_tool_index}]}@" input.nml
    sed -i -e "s@OBSSEQSECONDS${obs_seq_tool_index}@${observation_seconds[${obs_seq_tool_index}]}@" input.nml
done

