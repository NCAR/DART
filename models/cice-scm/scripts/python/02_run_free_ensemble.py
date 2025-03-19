#################################################################
## IMPORT NECESSARY LIBRARIES                                ##
###############################################################
# This tool MUST be run with an environment in which f90nml exists
import f90nml 
import glob
import os
import sys
import xarray as xr
import numpy as np

###############################################################
## SET MANUAL INPUTS HERE                                    ##
###############################################################

# set the free case name and the spinup it is based on
case_name = sys.argv[1]
user = sys.argv[2]
spinup_case = sys.argv[3]

# set the relevant directory names 
# where is this project repo located on your system?
project_dir = '/glade/work/'+user+'/Projects/cice-scm-da/'
# where is the Icepack installation located on your system?
icepack_dir = '/glade/work/'+user+'/Icepack/'
# where is the scratch directory Icepack will write output to?
scratch_dir = '/glade/derecho/scratch/'+user+'/'

# set machine name and compiler
machine = 'derecho'
compiler = 'intel'

# choose a simulation length and find the ensemble size and spinup length
spinup_length = len(glob.glob(scratch_dir + '/ICEPACK_RUNS/'+spinup_case+'/mem0001/restart/*.nc'))
ensemble_size = len(glob.glob(scratch_dir + '/ICEPACK_RUNS/'+spinup_case+'/mem*'))
simulation_years = 5

# choose a lateral flux option
flux = sys.argv[4]

# choose a location (this is read assuming that the spinup case is named "spinup_LOCATION")
location = spinup_case[7:]

if location != 'default':
    # comd = 'cp '+project_dir+'/data/forcings/'+location+'/JRA55/ATM_FORCING_*.txt '+project_dir+'/data/forcings/JRA55/'
    # os.system(comd)
    # comd = 'cp '+project_dir+'/data/forcings/'+location+'/ISPOL_2004/OCN_FORCING_*.txt '+project_dir+'/data/forcings/ISPOL_2004/'
    # os.system(comd)

    locs = {'Barents': [1.309, 0.698132], 
            'CoastalCanada': [1.41372, 6.24828], 
            'SibChuk': [1.3183906501, 3.0446500856], 
            'CentralArctic': [1.53589, 0]}
# else:
#     comd = 'cp '+project_dir+'/data/forcings/OCN_PERTS/Original/OCN_FORCING_*.txt '+project_dir+'/data/forcings/ISPOL_2004/'
#     os.system(comd)
#     locs = {'default': [0, 0]}

# This code comes equipped with the ability to perturb param-
# eters in the sea ice code. Your choice of perturbations should
# be the same as those used in the spinup ensemble. To do so,
# you must specify the same parameters supplied to the 
# 01_spinup_ensemble.py script. 

# The code currently supports the following perturbation options:
# R_snw: Snow grain radius tuning parameter (unitless)
# ksno: Snow thermal conductivity (W/m/K)
# dragio: Ice-ocean drag coefficient (unitless)
# hi_ssl: Ice surface scattering layer thickness (m)
# hs_ssl: Snow surface scattering layer thickness (m)
# rsnw_mlt: Snow melt rate (kg/m^2/s)
# rhoi: Ice density (kg/m^3)
# rhos: Snow density (kg/m^3)
# Cf: ratio of ridging work to PE change in ridging (unitless)
# atm: atmospheric forcing 
# ocn: oceanic forcing

###############################################################
## BEGIN LAUNCH- DO NOT EDIT BELOW THIS LINE                 ##
###############################################################
#-----------------------------------------#
# 0. Read parameter inputs
#-----------------------------------------#
if len(sys.argv) > 5:
    perturb = [sys.argv[i] for i in range(5, len(sys.argv))]
else:
    perturb = []

#-----------------------------------------#
# 1. Set parameter details
#-----------------------------------------#
# Defaults
# dR_snw = -2.0
dR_snw = 1.5
dksno = 0.3
dCf = 17
ddragio = 0.00536
dhi_ssl = 0.05
dhs_ssl = 0.04
drsnw_mlt = 1500.0
drhoi = 917.0
drhos = 330.0
ddt_mlt = 1.0

# Perturbed
parameters = xr.open_dataset(project_dir+ '/data/forcings/ICE_PERTS/parameters_30_cice5.nc')
zeros = np.zeros(ensemble_size)

if 'R_snw' in perturb:
    R_snw = list(parameters.R_snw.values)
else:
    R_snw = list(zeros+dR_snw)

if 'ksno' in perturb:
    ksno = list(parameters.ksno.values)
else:
    ksno = list(zeros+dksno)

if 'dragio' in perturb:
    dragio = list(parameters.dragio.values)
else:
    dragio = list(zeros+ddragio)

if 'hi_ssl' in perturb:
    hi_ssl = list(parameters.hi_ssl.values)
else:
    hi_ssl = list(zeros+dhi_ssl)

if 'hs_ssl' in perturb:
    hs_ssl = list(parameters.hs_ssl.values)
else:
    hs_ssl = list(zeros+dhs_ssl)

if 'rsnw_mlt' in perturb:
    rsnw_mlt = list(parameters.rsnw_melt.values)
else:
    rsnw_mlt = list(zeros+drsnw_mlt)

if 'rhoi' in perturb:
    rhoi = list(parameters.rhoi.values)
else:
    rhoi = list(zeros+drhoi)

if 'rhos' in perturb:
    rhos = list(parameters.rhos.values)
else:
    rhos = list(zeros+drhos)

if 'Cf' in perturb:
    Cf = list(parameters.Cf.values)
else:
    Cf = list(zeros+dCf)

if 'dt_mlt' in perturb:
    print('dt_mlt is not in pre-calculated parameter perturbations!')
    dt_mlt = list(zeros+ddt_mlt)
else:
    dt_mlt = list(zeros+ddt_mlt)

#-----------------------------------------#
# 2. Set up an Icepack case
#-----------------------------------------#
# go to the Icepack directory
os.chdir(icepack_dir)

# set up a new case 
comd = './icepack.setup -c '+case_name+' -m '+machine+' -e '+compiler
os.system(comd)

#-----------------------------------------#
# 3. Build the case 
#-----------------------------------------#
# go to the case directory
os.chdir(case_name)

# build the case 
comd = './icepack.build'
os.system(comd)

# check that the case was built correctly
storage_dir = scratch_dir + '/ICEPACK_RUNS/' + case_name
if os.path.exists(storage_dir) is False:
    print('Model did not build correctly! Please rebuild model.')
else:
    print('Model with following perturbations has been built:', perturb)

#-----------------------------------------#
# 4. Run the case
#-----------------------------------------#
simulation_length = 8760*simulation_years
year_init = 2011
year_end = year_init + simulation_years
final_year = '{0}'.format('%04d' % year_end)

#-------------------------------------#
# # 5. Begin cycling through ensemble   
#-------------------------------------#
mem = 1
while mem <= ensemble_size:
    inst_string ='{0}'.format('%04d' % mem) 
    print('Running member '+inst_string+'...')

    #create history and restart directories for the run
    os.chdir(storage_dir)
    os.makedirs('mem' + inst_string + '/history/')
    os.makedirs('mem' + inst_string + '/restart/')

    # link the model executable for each member to the main one built for the case
    os.symlink(storage_dir+'/icepack','mem' + inst_string+'/icepack')
    
    # begin working on an individual ensemble member
    os.chdir(storage_dir+'/mem' + inst_string)
    
    # handle restarts
    runtype_flag = True
    restart_flag = True
    
    ## if you used the cycling spinup option, your restart file will be named something like '2012-01-01.year10'
    # restart_file = scratch_dir + '/ICEPACK_RUNS/'+spinup_case+'/mem'+inst_string+'/restart/iced.2012-01-01-00000.year'+str(spinup_length)+'.nc'
    ## standard restart file option
    restart_file = scratch_dir + '/ICEPACK_RUNS/'+spinup_case+'/mem'+inst_string+'/restart/iced.2011-01-01-00000.nc'
    
    # read namelist template
    # namelist = f90nml.read(project_dir + '/data/templates/ICEPACK_input.nml.template_noflux')
    namelist = f90nml.read(project_dir + '/data/templates/ICEPACK_input.nml.template_JRA55_flux')

    # set case settings 
    namelist['setup_nml']['year_init'] = year_init
    namelist['setup_nml']['npt'] = simulation_length
    namelist['setup_nml']['restart'] = restart_flag
    namelist['setup_nml']['runtype_startup'] = runtype_flag
    namelist['setup_nml']['dumpfreq'] = 'd'
    if restart_flag is not True:
        namelist['setup_nml']['ice_ic'] = 'default'
    else:
        namelist['setup_nml']['restart_dir'] = './restart/'
        namelist['setup_nml']['ice_ic'] = restart_file
    if location != 'default':
        namelist['setup_nml']['input_lat'] = locs[location][0]
        namelist['setup_nml']['input_lon'] = locs[location][1]
    else:
        print('Using default location for forcing files.')

    # change namelist parameters of note
    namelist['thermo_nml']['ksno'] = ksno[mem-1]
    namelist['thermo_nml']['rhoi'] = rhoi[mem-1]
    namelist['shortwave_nml']['r_snw'] = R_snw[mem-1]
    namelist['shortwave_nml']['hi_ssl'] = hi_ssl[mem-1]
    namelist['shortwave_nml']['hs_ssl'] = hs_ssl[mem-1]
    namelist['shortwave_nml']['rsnw_mlt'] = rsnw_mlt[mem-1]
    namelist['shortwave_nml']['dt_mlt'] = dt_mlt[mem-1]
    namelist['snow_nml']['rhos']= rhos[mem-1]
    namelist['dynamics_nml']['Cf'] = Cf[mem-1]
    namelist['dynamics_nml']['dragio'] = dragio[mem-1]

    # set namelist forcing options
    namelist['forcing_nml']['data_dir'] = project_dir + '/data/forcings/'+location+'/free/'
    namelist['forcing_nml']['ycycle'] = simulation_years
    if 'atm' in perturb:
        namelist['forcing_nml']['atm_data_file'] = 'ATM_FORCING_'+inst_string+'.txt'
    else:
        namelist['forcing_nml']['atm_data_file'] = 'ATM_FORCING_0001.txt'
    if 'ocn' in perturb:
        namelist['forcing_nml']['ocn_data_file'] = 'OCN_FORCING_PERT_'+inst_string+'.txt'
    else:
        namelist['forcing_nml']['ocn_data_file'] = 'OCN_FORCING_'+inst_string+'.txt'
    
    namelist['forcing_nml']['lateral_flux_type'] = flux
    # write namelist to needed file
    namelist.write('icepack_in',force=True)

    # run the member instance and dump output 
    comd = './icepack > icepack.out'
    os.system(comd)

    # check the output file for successful model completion
    check_finished = 'ICEPACK COMPLETED SUCCESSFULLY'
    txt = open('icepack.out').readlines()
    if check_finished not in txt:
        AssertionError('Icepack did not run correctly! Process stopped.')
    else:
        print('Icepack ran successfully!')
    
    # advance to the next ensemble member 
    mem += 1

check_restarts = glob.glob(storage_dir + '/mem*/restart/iced.'+final_year+'-01-01-00000.nc')
# print(check_restarts)
if len(check_restarts) != ensemble_size:
    AssertionError('Free run did not complete as expected. Some restarts are missing. Processed stopped.')
else:   
    print('Free run complete. Submitting postprocessing...')

output_path = '/glade/work/'+user+'/Projects/cice-scm-da/data/processed/ensemble/'+case_name+'/'
if os.path.exists(output_path) == False:
    os.mkdir(output_path)

files = sorted(glob.glob('/glade/derecho/scratch/'+user+'/ICEPACK_RUNS/'+case_name+'/mem*/history/*.nc'))
DS = []
for file in files:
    ds = xr.open_dataset(file).isel({'ni':2}).drop(['ntrcr','ni','trcr','trcrn'])
    DS.append(ds)

ens_ds = xr.concat(DS, dim='member')

ens_ds = ens_ds.resample(time = '1D').mean()
ens_ds['hi'] = ens_ds.vice/ens_ds.aice

ens_ds.to_netcdf(output_path+'/postprocessed_ensemble.nc')
ens_ds.mean(dim='member').to_netcdf(output_path+'/postprocessed_ensemble_mean.nc')
ens_ds.std(dim='member', ddof=1).to_netcdf(output_path+'/postprocessed_ensemble_std.nc')

# check if files have been written
if len(glob.glob(output_path+'/postprocessed_ens*.nc')) == 3:
    print('Postprocessing complete! PROCESSED FINISHED.')
else:
    print('Postprocessing failed!')