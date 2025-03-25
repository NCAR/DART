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
import shutil

from datetime import datetime, timedelta
from time import sleep, perf_counter
import pandas as pd
###############################################################
## SET MANUAL INPUTS HERE                                    ##
###############################################################

# set the assimilation branch case name
case = sys.argv[1] # SIT_f1_NORM_test
user = sys.argv[2]

# set the source free case name
branch_case = sys.argv[3] 

# set the assimilation dates
initialization_date = datetime(int(sys.argv[4]),int(sys.argv[5]),int(sys.argv[6]))

# choose number of years to forecast
no_years = int(sys.argv[7])

###############################################################
## DEFINE HELPER FUNCTIONS                                   ##
###############################################################

def construct_case(case):

    # set the relevant directory names
    # where is the Icepack installation located on your system?
    icepack_dir = '/glade/work/'+user+'/Icepack/'
    # where is the scratch directory Icepack will write output to?
    scratch_dir = '/glade/derecho/scratch/'+user+'/'

    #-----------------------------------------#
    # 2. Set up an Icepack case
    #-----------------------------------------#
    # go to the Icepack directory
    os.chdir(icepack_dir)

    # set up a new case
    comd = './icepack.setup -c '+case+' -m derecho -e intel'
    os.system(comd)

    #-----------------------------------------#
    # 3. Build the case 
    #-----------------------------------------#
    # go to the case directory
    os.chdir(case)

    # build the case
    comd = './icepack.build'
    os.system(comd)

    # check that the case built correctly
    storage_dir = scratch_dir + '/ICEPACK_RUNS/' + case
    if os.path.exists(storage_dir) is False:
        print('Model did not build correctly! Please rebuild model.')
        # AssertionError('Model did not build correctly! Please rebuild model.')
    else:
        print('Model has been built!')


def forecast(case, branch_case, assim_date, forecast_years):

    case_dir  = '/glade/derecho/scratch/'+user+'/ICEPACK_RUNS/'+case+'/'
    branch_dir = '/glade/derecho/scratch/'+user+'/ICEPACK_RUNS/'+branch_case+'/'


    # Get date information
    year = assim_date.year
    mon = assim_date.month
    day = assim_date.day
    date_str = '{0}'.format('%04d'%year) + '-{0}'.format('%02d'%mon) + '-{0}'.format('%02d'%day)
    
    # Set up directory
    ensemble_size = len(glob.glob('/glade/derecho/scratch/'+user+'/ICEPACK_RUNS/'+branch_case+'/mem*'))

    # Move filter output to case directory
    if os.path.exists(case_dir + '/forecasts/'+date_str+'/') is False:
        os.makedirs(case_dir + '/forecasts/'+date_str+'/')
   
    mem = 1
    while mem <= ensemble_size:
        inst_string = '{0}'.format('%04d'%mem) 
        os.chdir(case_dir)
        os.makedirs('mem' + inst_string+'/history/')
        os.makedirs('mem' + inst_string+'/restart/')
                    
        # link the model executable for each member to the main one built for the case
        os.symlink(case_dir+'/icepack','mem' + inst_string+'/icepack')
    
        # Move to the case directory member folder to run the model forward 
        os.chdir(case_dir+'/mem'+inst_string+'/')

        # Copy namelist from branch_case
        shutil.copy(branch_dir +'/mem'+inst_string+'/icepack_in', 'icepack_in')

        # Edit the namelist file in order to run a forecast with icepack
        icepack_nml = f90nml.read('icepack_in')
        icepack_nml['setup_nml']['ice_ic'] = branch_dir + '/mem'+inst_string+'/restart/iced.'+date_str+'-00000.nc'
        icepack_nml['setup_nml']['npt'] = 24*365*forecast_years
        icepack_nml['setup_nml']['runtype_startup'] = False
        icepack_nml.write('icepack_in', force=True)

        # Run Icepack
        comd = './icepack > icepack.out'
        os.system(comd)
        sleep(2)


        # Clean up 
        shutil.copy('history/icepack.h.{0}{1}{2}'.format('%04d'%year, '%02d'%mon, '%02d'%day)+'.nc', case_dir+'/forecasts/'+date_str+'/icepack.h.'+date_str+'_'+inst_string+'.nc')
        shutil.move('ice_diag.full_ITD', case_dir+'/forecasts/'+date_str+'/ice_diag.full_ITD_'+inst_string)
        os.remove('icepack.out')
        comd = 'rm ice_diag.*'
        os.system(comd)

        # Advance to the next member
        mem += 1 

    print('Done forecasting '+case+' for '+date_str+'...')
    # end = perf_counter()

    # print(f'Execution time was {end - start:0.4f} seconds')
    return

###############################################################
## PERFORM CYCLING- DO NOT EDIT BELOW THIS LINE              ##
###############################################################

start = perf_counter()

construct_case(case)
forecast(case, branch_case, initialization_date, no_years)

end = perf_counter()
print(f'Execution time was {(end - start)/60:0.4f} minutes')
    
