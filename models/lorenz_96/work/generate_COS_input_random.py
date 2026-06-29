#!/bin/python

# This script should generate the input textfile for ./create_obs_sequence as COS_input.txt

import numpy as np

n_max_obs=10000
copies_of_data=0
qc_fields=0

RNG = np.random.default_rng()
n_obs = RNG.integers(10,101)
obs_locations = [RNG.random() for i_obs in range(n_obs)]
var_name="RAW_STATE_VARIABLE"
output_name="set_def_random.out"

with open("COS_input.txt", "w") as my_file:
    my_file.write("{}\n".format(n_max_obs))
    my_file.write("{}\n".format(copies_of_data))
    my_file.write("{}\n".format(qc_fields))
    for i_obs in range(n_obs):
        my_file.write("1\n")
        my_file.write("{}\n".format(var_name))
        my_file.write("{}\n".format(obs_locations[i_obs]))
        my_file.write("0 0\n") # obs time
        my_file.write("{}\n".format(RNG.uniform(0.5,3))) # variance
    my_file.write("-1\n")
    my_file.write("{}".format(output_name))



