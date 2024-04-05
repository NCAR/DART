#!/usr/bin/env python
# coding: utf-8

# When given a NetCDF file from a deterministic run of WRF-Hydro, this script
# creates daily obs_sequence files, denoted by YYYYMMDD strings appended to
# the end of the obs_seq string.

# This script needs one third-party module: netcdf4-python.
# On CISL resources (Cheyenne and Casper) please run these two commands first:
# > module load python
# > ncar_pylib

# Then run this script:
# > python pmo_osse.py

# IMPORT STANDARD LIBRARIES

from __future__ import print_function
from __future__ import division
import os
import time
import sys
import datetime
from math import pi

# IMPORT THIRD PARTY MODULE(S)
import netCDF4

# PRINT SCRIPT INFORMATION
# Print script data in case output is redirected from standard output to file.
this_script = sys.argv[0]
this_script = this_script[0:-3]
# This prints the name of the script
print('Running', this_script)
# This prints the last time the script was modified
print('Which was last modified at',
      time.ctime(os.path.getmtime(os.path.realpath(__file__))))
print('\n')

# CONSTANTS
# Get the name of the user running the script so we can output to their work
# directory.
user = os.environ['USER']
# Frequency of output
ntimes_per_day = 24
# Following the observation error procedure from
# create_identity_streamflow_obs.f90 we define a mininum error bound and a
# observational error fraction (40%) and pick whichever is larger.
min_err                = 0.1
max_err                = 1000000.0
obs_fraction_for_error = 0.4
deg2rad                = pi/180.0

# STRINGS AND PATHS
# Change the following strings to match the experiment names and locations of
# the files.
domain_name = 'drb'
pmo_name    = 'osse_id2'
input_path  = '/glade/work/gharamti/DART/DART_development/models/wrf_hydro/pmo/drb_mem0/'
output_path = '/glade/work/' + user + '/wrfhydro_dart/'

# Check to see if the output path exists.
obs_seq_path = output_path + domain_name + '/obs_seqs/' + pmo_name + '/'
if not os.path.exists(obs_seq_path):
    # If not, create the path.
    print('Making directory for obs_seq files:', obs_seq_path)
    os.makedirs(obs_seq_path)

# Perfect Model Output
# This file, when created by DART'S PMO, is typically called perfect_output.nc
pmo_path = input_path + 'pmo_drb_truth_des_gages.nc'
pmo_all_path = input_path + 'pmo_drb_truth_all_gages.nc'
# Route Link File
rl_path  = input_path + 'RouteLink.nc'

# PMO
# Get the necessary data from the PMO file
pmo_pointer = netCDF4.Dataset(pmo_path)
ntimes = len(pmo_pointer.dimensions['time'])
print('ntimes:', ntimes)

pmo_all_pointer = netCDF4.Dataset(pmo_all_path)

# TIMES
# Use the units string on the time variable to assign integers for year, month,
# day, etc
pmo_time     = pmo_pointer.variables['time']
start_year   = int(pmo_time.units[12:16])
start_month  = int(pmo_time.units[17:19])
start_day    = int(pmo_time.units[20:22])
start_hour   = int(pmo_time.units[23:25])
start_minute = 0
start_second = 0

print('Start date:', start_year, start_month, start_day, start_hour)

# Create an integration start_time using the integers from the time units
# string.
integration_start_time = datetime.datetime(start_year, start_month, start_day,
                                           start_hour, start_minute,
                                           start_second)

# If obs_seq files should only be output after a certain day, change it here.
# The default behavior is that the output_start_day is the same as the
# deterministic run start day.
output_start_day = datetime.datetime(start_year, start_month, start_day)
# output_start_day = datetime.datetime(2018, 9, 7)

# This loop starts the observation sequence loop only after the
# output_start_day specified by the user.
start_time = 0
nday = -1
ndays = int(ntimes/ntimes_per_day)

for iday in range(0, ndays):
    if output_start_day > integration_start_time+datetime.timedelta(days=iday):
        nday = nday + 1
        start_time = start_time + ntimes_per_day

# DART uses time since 1601-01-01 00:00:00
overall_start_time = datetime.datetime(1601, 1, 1, 0, 0, 0)
# Get the time since DART start time at the beginning of the file
file_start_time_delta = integration_start_time-overall_start_time
print('DART start time:', file_start_time_delta.seconds, file_start_time_delta.days)

# Get feature information from the perfect obs file
nfeatures     = len(pmo_all_pointer.dimensions['feature_id'])
nfeatures_des = len(pmo_pointer.dimensions['feature_id']) 
print('RL gauges:', nfeatures, ', desired ones:', nfeatures_des)

nobs     = ntimes*nfeatures
nobs_day = ntimes_per_day*nfeatures

pmo_reach_id     = pmo_all_pointer.variables['feature_ids']
pmo_time         = pmo_all_pointer.variables['time']
pmo_streamflow   = pmo_all_pointer.variables['streamflow']

pmo_des_reach_id = pmo_pointer.variables['feature_ids'] 

# ROUTELINK
# Get the necessary data from the Route Link file.
rl_pointer = netCDF4.Dataset(rl_path)
# Get the variables from the Route Link file.
lat  = rl_pointer.variables['lat']
lon  = rl_pointer.variables['lon']
alt  = rl_pointer.variables['alt']
link = rl_pointer.variables['link']

# GAGE LISTS
# Build lists with the following:
# 1. Index of the link with the desired gage
ilinks = []
# 2. Latitude of link
lats = []
# 3. Longitude of link
lons = []
# 4. Altitude of link
alts = []

# Loop through the reach ids to build the lists
print('\n')
print('Looping through the links in the Route Link file to get the location '
      'data for each desired gage.')
print('Thank you', user, 'for your patience.')
print('\n')

gg = 0 
for ipmo_reach, this_pmo_reach in enumerate(pmo_reach_id):
    for ilink, this_link in enumerate(link):
        if this_pmo_reach == this_link:
            gg = gg + 1
            print('Gauge no:', gg)
            print('Feature ID:', this_pmo_reach, 'and Link Index:', ilink+1) 
            print('Location: lat', lat[ilink], ', lon', lon[ilink], ', alt', alt[ilink], '\n')
            ilinks.append(ilink)
            this_lat = lat[ilink]*deg2rad
            lats.append(this_lat)
            this_lon = lon[ilink]
            if this_lon < 0.0:
                this_lon = this_lon+360.0
            this_lon = this_lon*deg2rad
            lons.append(this_lon)
            alts.append(alt[ilink])

# OBS SEQUENCE LOOP
# Loop through the times in the PMO file to build the obs_seq files

for itime in range(start_time, ntimes):
    # The commented line assumes that we want to make the observation time half
    # an observation period ahead of when it actually occurs so that the window
    # is centered on when the observation was taken. Do we want this?
    # If so, uncomment the next line and comment the following one.
    # this_time = file_start_time_delta+datetime.timedelta(hours=itime) - \
    #             datetime.timedelta(hours=0.5)
    this_time = file_start_time_delta+datetime.timedelta(hours=itime)

    # If the time index modulo the number of times per day equals zero, then
    # that means we're at the start of a day and we should thus open a new
    # obs_seq file designated with the proper YYYYMMDD string and write the
    # appropriate header information to it.
    if itime % ntimes_per_day == 0:
        # Get the YYYYMMDD string of this new day so that we can name the
        # obs_seq file appropriately.
        nday = nday+1

        file_start_day = integration_start_time+datetime.timedelta(days=nday)
        time_string = str(file_start_day.year) + \
            str(file_start_day.month).zfill(2) + \
            str(file_start_day.day).zfill(2)

        # Append 'obs_seq.' with the YYYYMMDD string
        obs_seq_file = obs_seq_path + 'obs_seq.' + time_string
        print('Writing observation sequences for day', str(nday).zfill(2),
              'to:', obs_seq_file)

        # Create the file pointer for writing
        obs_seq = open(obs_seq_file, 'w')

        # Write the header strings 'obs_sequence', 'obs_kind_definitions',
        # etc to the file
        print(' obs_sequence', file=obs_seq)
        print('obs_kind_definitions', file=obs_seq)
        # There aren't any obs_kind_definitions because this file only contains
        # identity obs.
        print('           0', file=obs_seq)
        print('num_copies:            1  num_qc:            1', file=obs_seq)
        print('num_obs:            ', nobs_day, ' max_num_obs:  ', nobs_day,
              file=obs_seq)
        print(' observation', file=obs_seq)
        print('QC VALUE', file=obs_seq)
        print('first:', '1'.rjust(8), 'last:', str(nobs_day).rjust(8),
              file=obs_seq)
        # Reset the obs counter so the OBS line is correct and the linked list
        # strings are correct as well.
        iobs = -1

    for ifeature in range(0, nfeatures):
        # Now we're looping through the actual gages and writing them, their
        # linked list, state_vector strings and observation error variance
        # to the obs_sequence files.
        iobs = iobs + 1
        this_obs = pmo_streamflow[itime, ifeature]

        # The observation error standard deviation is specified as the larger
        # of the observation magnitude times error fraction or the minimum
        # error threshold.
        if any(pmo_des_reach_id == pmo_reach_id[ifeature]):
           obs_err = max(this_obs*obs_fraction_for_error, min_err)
        else: 
           obs_err = max_err

        # Square it to get the variance
        obs_var = obs_err*obs_err

        # The observations are 1 indexed, but python loops are 0 indexed so
        # we add 1 to the observation index before writing to the file
        print('OBS', str(iobs+1).rjust(10), file=obs_seq)
        # Write the value of the observation
        print(this_obs, file=obs_seq)
        # write the QC value, 0
        print('1.000000000000000E+000', file=obs_seq)
        
        # The linked list line has three configurations
        if iobs == 0:
            # If it's the first observation, the first integer is -1
            print('-1'.rjust(4), str(iobs+2).rjust(4), '-1'.rjust(4),
                  file=obs_seq)
        elif iobs == nobs_day-1:
            # If it's the last observation the second integer is -1
            print(str(iobs).rjust(4), '-1'.rjust(4), '-1'.rjust(4),
                  file=obs_seq)
        else:
            # If it's any other observation the first integer is the obs
            # number minus 1 (assuming the observations are 1 indexed) and
            # the second integer is obs number plus one (again assuming
            # observations are 1 indexed).
            print(str(iobs).rjust(4), str(iobs+2).rjust(4), '-1'.rjust(4),
                  file=obs_seq)
        # Then we have the obdef section of the observation.
        print('obdef', file=obs_seq)
        # This is a 3-D observation with....
        print('loc3d', file=obs_seq)
        # ...latitude, longitude, altitude and the -1 denotes that the
        # vertical coordinate is a surface value, VERTISSURFACE.
        print(str(lons[ifeature]).rjust(12), str(lats[ifeature]).rjust(12),
              str(alts[ifeature]).rjust(12), '3'.rjust(12), file=obs_seq)
        print('kind', file=obs_seq)
        # Since these are identity observations they're the negative of the
        # position within the state vector.
        print('        -'+str(ilinks[ifeature]+1), file=obs_seq)
        # The time of the observation is days and seconds since 1601-01-01
        # 00:00:00
        print(this_time.seconds, this_time.days, file=obs_seq)
        # Finally, write the observation error variance
        print(obs_var, file=obs_seq)

    # If the next time modulo the number of times per day equals zero, for
    # example if we just wrote the observations for 11PM (23 hours), then we
    # just wrote the last appropriate observations to this file and we need
    # to close the file.
    if itime+1 % ntimes_per_day == 0:
        obs_seq.close()


