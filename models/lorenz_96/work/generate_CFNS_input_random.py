#!/bin/python

# This script should generate the input textfile for ./create_fixed_network_seq called CFNS_intput.txt

import argparse
import numpy as np

parser = argparse.ArgumentParser()
parser.add_argument("obs_template_name", type=str)
parser.add_argument("output_file_name", type=str)
parser.add_argument("max_obs_number", type=int)
parser.add_argument("end_day", type=int)
args = parser.parse_args()

end_day = args.end_day
max_obs_number = args.max_obs_number
output_file_name = args.output_file_name
obs_template_name = args.obs_template_name
obs_frequency_type = 2

RNG = np.random.default_rng()
n_obs = RNG.integers(10, int(min(201,max_obs_number/10)))
observation_times = [(RNG.integers(end_day), RNG.integers(86400)) for i_obs in range(n_obs)]

with open("CFNS_input.txt", "w") as my_file:
    my_file.write("{}\n".format(obs_template_name))
    my_file.write("{}\n".format(obs_frequency_type))
    my_file.write("{}\n".format(max_obs_number))
    for i_obs in range(n_obs):
        my_file.write("{} {}\n".format(observation_times[i_obs][0], observation_times[i_obs][1]))
    my_file.write("-1 0\n{}".format(output_file_name))

        






