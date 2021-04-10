# Tests

## General
Test scripts for different platforms/machines are found in `local/`. In each subdirectory, scripts call pytests supplying the variables

```
pytest  -s  -vv  --pdb  --exp_yaml $yaml_file  --exp_answer_file $answer_file  --use_existing_build 
```

where `--use_existing_build` is optional and skips dart and wrf_hydro compiling if they already exist (to speed up test development).


## Cheyenne
Currently (12/20) just setting up tests to run on chyenne. The single experiment test is run by 

```
local/cheyenne/suite_1.sh
```

It 

  1. Uses the florence_cutout domain (found on cheyenne).
  1. Sets up the experiment (dir structures, generates the initial ensemble, and creates some obs_seq files) defined by `experiment_test_files/florence_cutout/cheyenne.yaml` and `experiment_test_files/florence_cutout/constructor.py`.
  1. Runs the filter experiment (3 members for 3 hours)
  1. Checks the stats of the member posteriors at hour 3, see `experiment_test_files/florence_cutout/answers.yaml`
  

