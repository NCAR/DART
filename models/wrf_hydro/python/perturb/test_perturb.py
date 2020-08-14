#!/usr/bin/python

#######################################################
# Calling
# python test_perturb.py chrtout_file n_ens_members out_dir
#
# example usage:
# python \
#    test_perturb.py \
#    /Users/jamesmcc/Downloads/04233300_DART/DART.201306-201309.channelOnly/FORCING.CHRTOUT/201308090000.CHRTOUT_DOMAIN1 \
#    5 \
#    /Users/jamesmcc/

import numpy as np
import os
from perturb_channel_only_forcing import *
import sys

#chrtout_file = '/Users/jamesmcc/Downloads/04233300_DART/DART.201306-201309.channelOnly/FORCING.CHRTOUT/201308090000.CHRTOUT_DOMAIN1'
#n_ens_members = 5
chrtout_file  = sys.argv[1]
n_ens_members = int(sys.argv[2])
out_dir       = sys.argv[3]

#######################################################
# Noise model for qSfcLatRunoff
# 0) Additive noise,
# 1) Zero-mean,
# 2) frac: Standard deviation is a fraction of the value,
# 3) min: Truncated below at min,
# 4) size: number of samples,
# 5) Closure takes a single value argument.

def close_trunc_gauss_sd_pct_value(frac, min, size):
    def the_closure(x):
        return np.maximum(x+np.random.normal(0.0, frac*x[0], size), min)
    return the_closure


trunc_gauss_sd_pct_value = close_trunc_gauss_sd_pct_value(.2, 0, n_ens_members)

# Noise model to be applied station-wise


def close_noise_by_station(variable_str, noise_func) :
    def the_closure(dataset):
        return dataset[variable_str].groupby('station').apply(noise_func)
    return the_closure


qsfclat_noise_by_station = close_noise_by_station('qSfcLatRunoff', trunc_gauss_sd_pct_value)
noise_model = qsfclat_noise_by_station

#######################################################
# Noise model for qBucket
## 0) Additive noise,
## 1) Zero-mean,
## 2) Standard deviation is a fixed at .02
## 3) min: Truncated below at min,
## 4) size: number of samples,
## 5) Closure takes a single value argument.
def close_trunc_gauss(min, size):
    def the_closure(x):
        return np.maximum(x+np.random.normal(0.0, .2, size), min)
    return the_closure

trunc_gauss = close_trunc_gauss(0, n_ens_members)

qbucket_noise_by_station = close_noise_by_station('qBucket', trunc_gauss)


#######################################################
# Make some noise!
newfiles = perturb_channel_only_forcing(chrtout_file, n_ens_members,
                                        qsfclat_noise_by_station,
                                        qbucket_noise_by_station,
                                        out_dir)
                                        
