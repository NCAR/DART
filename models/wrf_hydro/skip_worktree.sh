#!/bin/bash

array=( \
        work/path_names_filter 
        work/path_names_model_mod_check
        work/path_names_perfect_model_obs
        ../../assimilation_code/modules/utilities/mpi_utilities_mod.f90
        ../../assimilation_code/modules/utilities/null_mpi_utilities_mod.f90
        ../noah/work/path_names_filter
        ../noah/work/path_names_model_mod_check
        ../noah/work/path_names_perfect_model_obs )

for ff in ${array[@]}; do
    echo $ff
    git checkout -- $ff
    git update-index --skip-worktree $ff || exit 1
done


exit 0
