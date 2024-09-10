#!/usr/bin/env python3
import numpy as np
import netCDF4 as nc
import matplotlib.pyplot as plt

#################################################################
####################### SCRIPT PARAMETERS #######################
#################################################################

# input and output datafiles
input_file_path = 'bats_bottle.txt'
output_txt_path = 'bats_climatology.txt'
output_nc_path  = 'bats_climatology.nc'

first_data_line = 61        # file is read starting at this line (so that header information is ignored)
month_columns   = [18, 19]  # first and last columns in the data file containing month values
depth_columns   = [64, 69]  # first and last columns in the data file containing depth values
num_obs_types   = 6         # number of variables being read from the datafile
max_num_data    = 16660     # upper bound on the number of data points to be processed

# first and last columns in the data file containing values for each variable
obs_val_columns = [[113, 119],  # O2
                   [137, 143],  # CO2
                   [145, 151],  # Alk
                   [153, 159],  # NO31
                   [170, 176],  # PO41
                   [178, 184]]  # Si1

# layer interface depths for data in prog_z.nc output file from MARBL, run "ncdump -v "z_i" prog_z.nc" to generate.
marbl_depths = [0, 5, 15, 25, 40, 62.5, 87.5, 112.5, 137.5, 175, 225, 275, 350, 450, 
                550, 650, 750, 850, 950, 1050, 1150, 1250, 1350, 1450, 1625, 1875, 2250, 
                2750, 3250, 3750, 4250, 4750, 5250, 5750, 6250]

# set to 'True' in order to produce a histogram of BATS observation vs depth for each month
plot_sample_hist = False

# human-readable names of the variables being read from the file
obs_names  = ["O2", "CO2", "Alk", "NO31", "PO41", "Si1"]

#################################################################
####################### MAIN PROGRAM ############################
#################################################################

clim_means    = np.zeros((num_obs_types, 12, len(marbl_depths)))
sq_means      = np.zeros((num_obs_types, 12, len(marbl_depths)))
sample_counts = np.zeros((num_obs_types, 12, len(marbl_depths)), dtype=int)

with open(input_file_path, "r") as bats_data:
    for line in bats_data.readlines()[first_data_line-1:]:
    
        # extracting the depth and month where this observation was taken
        depth_val   = float(line[depth_columns[0]-1:depth_columns[1]])
        month_index = int(line[month_columns[0]-1:month_columns[1]]) - 1

        # assigning this observation to one of the MARBL layers
        depth_index = len(marbl_depths) - 1
    
        while((depth_val < marbl_depths[depth_index]) and (depth_index > 0)):
            depth_index -= 1
    
        # adding this observation into the climatology
        for obs_type_index in range(num_obs_types):
            obs_val = float(line[obs_val_columns[obs_type_index][0]-1:obs_val_columns[obs_type_index][1]])
        
            if(abs(obs_val + 999.0) > 1e-8):    # this filters out 'missing' values
                prev_mean    = clim_means[obs_type_index, month_index, depth_index]
                prev_sqmean  = sq_means[obs_type_index, month_index, depth_index]
                prev_samples = sample_counts[obs_type_index, month_index, depth_index]

                new_mean   = (prev_samples*prev_mean + obs_val)/(prev_samples + 1)
                new_sqmean = (prev_samples*prev_sqmean + obs_val**2)/(prev_samples + 1)

                clim_means[obs_type_index, month_index, depth_index] = new_mean 
                sq_means[obs_type_index, month_index, depth_index]   = new_sqmean
            
                sample_counts[obs_type_index, month_index, depth_index] += 1

                print("month =",month_index,", type =",obs_type_index,", depth =",marbl_depths[depth_index],", value =",obs_val)


# setting up a CSV file to write the climatology into; this will be read by the DART data converter.
output_txt = open(output_txt_path, "w")

# setting up a netCDF file to write the climatology into; this will be useful for creating MARBL-DART diagnostics.
output_nc = nc.Dataset(output_nc_path, "w")

output_nc.createDimension("Month")
output_nc.createDimension("Layer")
output_nc.createVariable("Depth", "double", ("Layer",))
output_nc["Depth"][:] = marbl_depths

for obs_type_index in range(num_obs_types):
    output_nc.createVariable(obs_names[obs_type_index]+"_value", "double", ("Month", "Layer",))
    output_nc.createVariable(obs_names[obs_type_index]+"_error_sd", "double", ("Month", "Layer",))
    output_nc.createVariable(obs_names[obs_type_index]+"_samples", "int", ("Month", "Layer",))

for month_index in range(12):
    for obs_type_index in range(num_obs_types):
        for depth_index in range(len(marbl_depths)):
            samples = sample_counts[obs_type_index, month_index, depth_index]
            output_nc[obs_names[obs_type_index]+"_samples"][month_index, depth_index] = samples

            if(samples > 0):
                depth   = marbl_depths[depth_index]
                mean    = clim_means[obs_type_index, month_index, depth_index]
                sq_mean = sq_means[obs_type_index, month_index, depth_index]

                variability_err = (sq_mean - mean**2)/samples
                obs_error       = (.1*mean)**2
                uncertainty     = np.sqrt(variability_err + obs_error)

                output_txt.write(str(month_index+1)+","+str(obs_type_index+1)+","+str(depth)+","+str(mean)+","+str(uncertainty)+"\n")

                output_nc[obs_names[obs_type_index]+"_value"][month_index, depth_index] = mean
                output_nc[obs_names[obs_type_index]+"_error_sd"][month_index, depth_index] = uncertainty

output_txt.close()
output_nc.close()

if(plot_sample_hist):
    month_names = ["January", "February", "March", "April", "May", "June", "July", "August", "September", "October", "November", "December"]
    month_plot  = 12*[None]

    fig = plt.figure(figsize = (6, 70))

    for month_index in range(12):
        month_plot[month_index] = fig.add_subplot(12, 1, month_index + 1)
        month_plot[month_index].set_xlabel("Depth (meters)")
        month_plot[month_index].set_xscale("log")
        month_plot[month_index].set_ylabel("Number of Data Points ("+month_names[month_index]+")")
        month_plot[month_index].set_yscale("log")
        
        (month_index == 0) and month_plot[month_index].set_title("Depth Distribution of BATS Data")

        for obs_type_index in range(num_obs_types):
            month_plot[month_index].plot(marbl_depths, sample_counts[obs_type_index, month_index, :], marker = "o", markerfacecolor = "none", label = obs_names[obs_type_index])
        
        month_plot[month_index].legend()

    plt.legend()
    plt.savefig("bats_depth_hist.pdf")
    plt.close(fig)
