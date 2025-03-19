import xarray as xr
import numpy as np
import os
import sys
import pandas as pd
import glob

def convert_times(input_times):
    new_times = np.zeros((input_times.shape[0]), dtype=np.int64)
    for t in range(input_times.shape[0]):
        new_times[t] = int(pd.to_datetime(str(input_times[t])).strftime('%Y%m%d%H%M'))
    return new_times

year1 = sys.argv[1]
year2 = sys.argv[2]
purpose = sys.argv[3]
ens_size = int(sys.argv[4])

LatsOfInterest = [75.53822, 81, 88, 75]
LonsOfInterest = [174.4456, 358, 0, 40]
PointsOfInterest = ['SibChuk', 'CoastalCanada' 'CentralArctic', 'Barents']
JRA55_forcing_path = '/glade/derecho/scratch/mollyw/C-ICESAT-2/ATMOSPHERE_FORCING/'

mems = np.arange(1, ens_size + 1)
temperature_data = xr.open_dataset(JRA55_forcing_path + 'mem01_JRA.v1.5_t_10_' + year1 + '.nc')

lon = temperature_data.longitude.values
lat = temperature_data.latitude.values
temperature_data.close()

jra_times = pd.date_range(year1+'-01-01', str(int(year2)+1)+'-01-01', freq='3H')
jra_times = jra_times[~((jra_times.month == 2) & (jra_times.day == 29))]
jra_times = jra_times[:-1]
ip_times = pd.date_range(year1+"-01-01", str(int(year2)+1)+"-01-01", freq='1H')
ip_times = ip_times[~((ip_times.month == 2) & (ip_times.day == 29))]
ip_times = ip_times[:-1]

jra_times = convert_times(jra_times)
ip_times = convert_times(ip_times)


for m in range(mems.shape[0]):
    label = str(m + 1).zfill(4)        
    print("working on mem{0:02d}".format(m + 1) + "...")

    t_10_files = [JRA55_forcing_path + '/mem{0:02d}_JRA.v1.5_t_10_'.format(m + 1) + str(yr) + '.nc' for yr in range(int(year1), int(year2)+1)]
    u_10_files = [JRA55_forcing_path + '/mem{0:02d}_JRA.v1.5_u_10_'.format(m + 1) + str(yr) + '.nc' for yr in range(int(year1), int(year2)+1)]
    v_10_files = [JRA55_forcing_path + '/mem{0:02d}_JRA.v1.5_v_10_'.format(m + 1) + str(yr) + '.nc' for yr in range(int(year1), int(year2)+1)]
    q_10_files = [JRA55_forcing_path + '/mem{0:02d}_JRA.v1.5_q_10_'.format(m + 1) + str(yr) + '.nc' for yr in range(int(year1), int(year2)+1)]
    swdn_files = [JRA55_forcing_path + '/mem{0:02d}_JRA.v1.5_swdn_'.format(m + 1) + str(yr) + '.nc' for yr in range(int(year1), int(year2)+1)]
    lwdn_files = [JRA55_forcing_path + '/mem{0:02d}_JRA.v1.5_lwdn_'.format(m + 1) + str(yr) + '.nc' for yr in range(int(year1), int(year2)+1)]
    prec_files = [JRA55_forcing_path + '/mem{0:02d}_JRA.v1.5_prec_'.format(m + 1) + str(yr) + '.nc' for yr in range(int(year1), int(year2)+1)]

    print('Concatenating files...')
    t_10 = xr.open_mfdataset(t_10_files).load()
    u_10 = xr.open_mfdataset(u_10_files).load()
    v_10 = xr.open_mfdataset(v_10_files).load()
    q_10 = xr.open_mfdataset(q_10_files).load()
    swdn = xr.open_mfdataset(swdn_files).load()
    lwdn = xr.open_mfdataset(lwdn_files).load()
    prec = xr.open_mfdataset(prec_files).load()
        
    for p in range(len(LatsOfInterest)):
        lat_pt = LatsOfInterest[p]
        lon_pt = LonsOfInterest[p]
        a_j = abs(lon - lon_pt)
        a_i = abs(lat - lat_pt)
        j = np.unravel_index(a_j.argmin(), a_j.shape)
        i = np.unravel_index(a_i.argmin(), a_i.shape)

        # Make a data directory for ocean forcing data
        PoI_path = '/glade/work/mollyw/Projects/cice-scm-da/data/forcings/' + PointsOfInterest[p]+'/'+purpose+'/JRA55/'
        if not os.path.exists(PoI_path):
            os.makedirs(PoI_path)

        file_path = os.path.join(PoI_path, 'ATM_FORCING_{0}.txt'.format(label))
        if os.path.exists(file_path):
            print("File {0} already exists. Skipping...".format(file_path))
            continue
        
        print("Extracting point data for ", PointsOfInterest[p])
        t_10_pt = list(t_10.t_10.values[:, i, j][:, 0])
        u_10_pt = list(u_10.u_10.values[:, i, j][:, 0])
        v_10_pt = list(v_10.v_10.values[:, i, j][:, 0])
        q_10_pt = list(q_10.q_10.values[:, i, j][:, 0])
        swdn_pt = list(swdn.swdn.values[:, i, j][:, 0])
        lwdn_pt = list(lwdn.lwdn.values[:, i, j][:, 0])
        prec_pt = list(prec.prec.values[:, i, j][:, 0])
        
        print('Interpolating data to hourly timesteps...')

        interp_t10 = np.interp(ip_times, jra_times, t_10_pt)
        interp_u10 = np.interp(ip_times, jra_times, u_10_pt)
        interp_v10 = np.interp(ip_times, jra_times, v_10_pt)
        interp_q10 = np.interp(ip_times, jra_times, q_10_pt)
        interp_swdn = np.interp(ip_times, jra_times, swdn_pt)
        interp_lwdn = np.interp(ip_times, jra_times, lwdn_pt)
        interp_prec = np.interp(ip_times, jra_times, prec_pt)

        with open(file_path, 'w') as new_file:
            new_file.writelines('#DSWSFC      DLWSFC    WNDU10    WNDV10    TEMP2M    SPECHUM    PRECIP\n')
            new_file.writelines('#W/m**2      W/m**2    m/s       m/s       K         Kg/Kg      kg/m**2/s\n')
            for l in range(interp_t10.shape[0]):
                new_file.writelines('{0:10.5f} {1:10.5f} {2:10.5f} {3:10.5f} {4:10.5f} {5:10.5f} {6:10.8f}\n'.format(
                    interp_swdn[l], interp_lwdn[l], interp_u10[l], interp_v10[l], interp_t10[l], interp_q10[l], interp_prec[l]))
                
    del t_10, q_10, u_10, v_10, swdn, lwdn, prec

print('Done creating ATM_FORCING files for {0}'.format(year1) + ' through {}'.format(year2)+' for all points of interest!')