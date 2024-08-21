###############################################################
## IMPORT NECESSARY LIBRARIES                                ##
###############################################################
import xarray as xr
import glob
import os
import sys
from time import perf_counter

###############################################################
## SET MANUAL INPUTS HERE                                    ##
###############################################################

case = sys.argv[1]
truth = int(sys.argv[2])

if sys.argv[3] == 'all':
    stages = ['preassim','pre_filter','post_filter','postprocessed','analysis','forecast']
else:
    stages = [sys.argv[i] for i in range(3,len(sys.argv))]

user = 'mollyw'
output_path = '/glade/work/'+user+'/Projects/cice-scm-da/data/processed/'+case+'/'
if not os.path.exists(output_path):
    os.makedirs(output_path)

###############################################################
## PERFORM POSTPROCESSING-- DO NOT EDIT BELOW                ##
###############################################################

start = perf_counter()

print('Working on case '+case+'...')
print('Postprocessing stages: '+str(stages))

for stage in stages:
    if stage in ['input','preassim','analysis','output']:
        base_path = '/glade/derecho/scratch/'+user+'/ICEPACK_RUNS/'+case+'/output_files/'
    elif stage in ['pre_filter','post_filter','postprocessed']:
        base_path = '/glade/derecho/scratch/'+user+'/ICEPACK_RUNS/'+case+'/analyses/'
    elif stage in'forecast':
        base_path = '/glade/derecho/scratch/'+user+'/ICEPACK_RUNS/'+case+'/forecasts/'
    else:
        print('stage '+stage+' not recognized')
        break
    
    if stage in ['preassim','analysis','pre_filter','postprocessed','post_filter','forecast']:
        # Cycle over ensemble members and concatenate in time
        if 'output_files' in base_path:
            for i in range(1,30):
                os.system('ncecat '+base_path+'????-??-??/'+stage+'_member_00'+f"{i:02d}"+'.nc '+output_path+stage+'_member_00'+f"{i:02d}"+'.nc')
                os.system('ncrename -d record,time '+output_path+stage+'_member_00'+f"{i:02d}"+'.nc')
        elif 'analyses' in base_path:
            for i in range(1,31):
                if i == truth:
                    if stage == 'postprocessed':
                        os.system('ncecat '+base_path+'????-??-??/truth_data_00'+f"{i:02d}"+'.nc '+output_path+'truth_data_00'+f"{i:02d}"+'.nc')
                        os.system('ncrename -d record,time '+output_path+'truth_data_00'+f"{i:02d}"+'.nc')
                    else:
                        continue
                else:
                    os.system('ncecat '+base_path+'????-??-??/'+stage+'_restart_00'+f"{i:02d}"+'.nc '+output_path+stage+'_member_00'+f"{i:02d}"+'.nc')
                    os.system('ncrename -d record,time '+output_path+stage+'_member_00'+f"{i:02d}"+'.nc')
        elif 'forecast' in base_path:
            for i in range(1,31):
                os.system('ncecat '+base_path+'????-??-??/icepack.h.????-??-??_00'+f"{i:02d}"+'.nc '+output_path+stage+'_member_00'+f"{i:02d}"+'.nc')
                # take an average over the time dimension
                os.system('ncwa -O -C -a time '+output_path+stage+'_member_00'+f"{i:02d}"+'.nc '+output_path+stage+'_member_00'+f"{i:02d}"+'.nc')
                # drop the time dimension (which should now have length 1) and rename the record dimension to time
                os.system('ncks -O -C -x -v time '+output_path+stage+'_member_00'+f"{i:02d}"+'.nc '+output_path+stage+'_member_00'+f"{i:02d}"+'.nc')
                os.system('ncrename -d record,time '+output_path+stage+'_member_00'+f"{i:02d}"+'.nc')
        else:
            print('Assimilation stage not recognized! Exiting.')
            break
        
        # Using python, concatenate all the ensemble member files 
        DS = []
        files = glob.glob(output_path+stage+'_member_00*.nc')
        for file in files:
            ds = xr.open_dataset(file)
            DS.append(ds)
        ens_ds = xr.concat(DS,dim='ensemble_member')
        ens_ds.to_netcdf(output_path+stage+'_ens.nc')

    # if the *_ens.nc file exists, delete the individual member files using os.remove()
    if os.path.exists(output_path+stage+'_ens.nc'):
        for f in glob.glob(output_path+stage+'_member_00*.nc'):
            os.remove(f)

    # Concatenate all the mean and sd files across time if they exist
    if os.path.exists(base_path+'2011-01-02/'+stage+'_mean.nc'):
        os.system('ncecat '+base_path+'????-??-??/'+stage+'_mean.nc '+output_path+stage+'_mean.nc')
        os.system('ncecat '+base_path+'????-??-??/'+stage+'_sd.nc '+output_path+stage+'_sd.nc')
        os.system('ncrename -d record,time '+output_path+stage+'_mean.nc')
        os.system('ncrename -d record,time '+output_path+stage+'_sd.nc')


    # Concatenate all the prior inflation files across time, if they exist
    if os.path.exists(base_path+'2011-01-02/'+stage+'_priorinf_mean.nc'):
        os.system('ncecat '+base_path+'????-??-??/'+stage+'_priorinf_mean.nc '+output_path+stage+'_priorinf_mean.nc')
        os.system('ncecat '+base_path+'????-??-??/'+stage+'_priorinf_sd.nc '+output_path+stage+'_priorinf_sd.nc')
        os.system('ncrename -d record,time '+output_path+stage+'_priorinf_mean.nc')
        os.system('ncrename -d record,time '+output_path+stage+'_priorinf_sd.nc')

    # if the stage is forecast, pre_filter, postprocessed, or post_filter, calculated the mean and sd across the ensemble_member dimension using nco
    if stage in ['forecast','pre_filter','postprocessed','post_filter']:
        os.system('ncwa -a ensemble_member '+output_path+stage+'_ens.nc '+output_path+stage+'_mean.nc')
        os.system('ncbo -O '+output_path+stage+'_ens.nc '+output_path+stage+'_mean.nc '+output_path+stage+'_devs.nc')
        os.system('ncwa -O -y rmssdn -a ensemble_member '+output_path+stage+'_devs.nc '+output_path+stage+'_sd.nc')
        # remove the temporary files
        os.remove(output_path+stage+'_devs.nc')

end = perf_counter()
print('Postprocessing completed! Files are located in '+output_path)
print(f'Execution time was {(end - start)/60:0.4f} minutes')
   