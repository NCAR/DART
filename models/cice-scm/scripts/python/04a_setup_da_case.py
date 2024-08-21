###############################################################
## IMPORT NECESSARY LIBRARIES                                ##
###############################################################
# This tool MUST be run with an environment in which f90nml exists
import f90nml 
import glob
import os
import xarray as xr
import numpy as np
import sys

from datetime import datetime

###############################################################
## SET MANUAL INPUTS HERE                                    ##
###############################################################

# set the free cases name and the spinup it is based on
case_name = sys.argv[1] # SIT_f1_NORM_test
user = sys.argv[2]
spinup_case = sys.argv[3] # spinup_test
flux = sys.argv[4]

# set the relevant directory names
# where is this project repo located on your system?
project_dir = '/glade/work/'+user+'/Projects/cice-scm-da/'
# where is the Icepack installation located on your system?
icepack_dir = '/glade/work/'+user+'/Icepack/'
# where is the DART installation located on your system?
dart_dir = '/glade/work/'+user+'/dart_manhattan/'
# where is the scratch directory Icepack will write output to?
scratch_dir = '/glade/derecho/scratch/'+user+'/'

# set the machine name and compiler
machine = 'derecho'
compiler = 'intel'

# set assimilation times
first_assim_time = datetime(2011, 1, 2)
model_init_time = datetime(2011, 1, 1)

###############################################################
## BEGIN LAUNCH- DO NOT EDIT BELOW THIS LINE                 ##
###############################################################

# determine the ensemble size
ensemble_size = len(glob.glob(scratch_dir + '/ICEPACK_RUNS/'+spinup_case+'/mem*'))

# determine the number of days between model initation and assimilation
days_to_assim = first_assim_time - model_init_time
assim_date_str = '{0}-{1}-{2}'.format('%04d'%first_assim_time.year, '%02d'%first_assim_time.month, '%02d'%first_assim_time.day)

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
dR_snw = -2.0
dksno = 0.3
dCf = 17
ddragio = 0.00536
dhi_ssl = 0.05
dhs_ssl = 0.04
drsnw_mlt = 1500.0
drhoi = 917.0
drhos = 330.0

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
    # AssertionError('Model did not build correctly! Please rebuild model.')
else:
    print('Model with following perturbations has been built:', perturb)

#-----------------------------------------#
# 4. Run the case
#-----------------------------------------#
simulation_length = 24 * days_to_assim.days
year_init = model_init_time.year

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
    restart_file = scratch_dir + '/ICEPACK_RUNS/'+spinup_case+'/mem'+inst_string+'/restart/iced.2012-01-01-00000.year10.nc'
    
    # read namelist template
    namelist = f90nml.read(project_dir + '/data/templates/ICEPACK_input.nml.template_1.3.1_flux')

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

    # change namelist parameters of note
    namelist['thermo_nml']['ksno'] = ksno[mem-1]
    namelist['thermo_nml']['rhoi'] = rhoi[mem-1]
    namelist['shortwave_nml']['r_snw'] = R_snw[mem-1]
    namelist['shortwave_nml']['hi_ssl'] = hi_ssl[mem-1]
    namelist['shortwave_nml']['hs_ssl'] = hs_ssl[mem-1]
    namelist['shortwave_nml']['rsnw_mlt'] = rsnw_mlt[mem-1]
    namelist['snow_nml']['rhos']= rhos[mem-1]
    namelist['dynamics_nml']['Cf'] = Cf[mem-1]
    namelist['dynamics_nml']['dragio'] = dragio[mem-1]

    # set namelist forcing options
    namelist['forcing_nml']['data_dir'] = project_dir + '/data/forcings/'
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
    # check_finished = 'ICEPACK_COMPLETED SUCCESSFULLY'
    # txt = open('icepack.out').readlines()
    # if check_finished not in txt:
    #     print('Icepack did not run correctly! Process stopped.')
    # else:
    #     print('Icepack ran successfully!')
    
    # advance to the next ensemble member 
    mem += 1

check_restarts = glob.glob(storage_dir + '/mem*/restart/iced.'+assim_date_str+'-00000.nc')
# print(check_restarts)
if len(check_restarts) != ensemble_size:
    print('Free run did not complete as expected. Some restarts are missing. Processed stopped.')
else:   
    print('PROCESS COMPLETE. Please check member directories.')
