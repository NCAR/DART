import xarray as xr
import numpy as np
import os
import sys
import pandas as pd

def convert_times(input_times):
    new_times = np.zeros((input_times.shape[0]), dtype=np.int64)
    for t in range(input_times.shape[0]):
        new_times[t] = int(pd.to_datetime(str(input_times[t])).strftime('%Y%m%d%H%M'))
    return new_times

year1 = sys.argv[1]
year2 = sys.argv[2]
purpose = sys.argv[3]
ens_size = int(sys.argv[4])

LatsOfInterest = [88] #[75.53822, 81, 88, 75]
LonsOfInterest = [0] #[174.4456, 358, 0, 40]
PointsOfInterest = ['CentralArctic']#['Sib_Chuk', 'CoastalCanada', 'CentralArctic', 'Barents']
JRA55_forcing_path = '/glade/derecho/scratch/mollyw/C-ICESAT-2/ATMOSPHERE_FORCING/'

mems = np.arange(1, ens_size + 1)
temperature_data = xr.open_dataset(JRA55_forcing_path + 'mem01_JRA.v1.5_t_10_' + year1 + '.nc')

lon = temperature_data.longitude.values
lat = temperature_data.latitude.values
temperature_data.close()

jra_times = pd.date_range(year1+'-01-01', str(int(year1)+1)+'-01-01', freq='3H')
jra_times = jra_times[~((jra_times.month == 2) & (jra_times.day == 29))]
jra_times = jra_times[:-1]
ip_times = pd.date_range(year1+"-01-01", str(int(year1)+1)+"-01-01", freq='1H')
ip_times = ip_times[~((ip_times.month == 2) & (ip_times.day == 29))]
ip_times = ip_times[:-1]

jra_times = convert_times(jra_times)
ip_times = convert_times(ip_times)

# print((int(year2)-int(year1) + 1 )*365*8)
# print(len(jra_times))
# print((int(year2)-int(year1) + 1)*365*24)
# print(len(ip_times))

for p in range(len(LatsOfInterest)):
    lat_pt = LatsOfInterest[p]
    lon_pt = LonsOfInterest[p]
    a_j = abs(lon - lon_pt)
    a_i = abs(lat - lat_pt)
    j = np.unravel_index(a_j.argmin(), a_j.shape)
    i = np.unravel_index(a_i.argmin(), a_i.shape)

    # Make a data directory for ocean forcing data
    PoI_path = '/glade/work/mollyw/Projects/cice-scm/data/forcings/' + PointsOfInterest[p]+'/'+purpose+'/JRA55/'
    if not os.path.exists(PoI_path):
        os.makedirs(PoI_path)

    for m in range(mems.shape[0]):
        label = str(m + 1).zfill(4)
        file_path = os.path.join(PoI_path, 'ATM_FORCING_{0}.txt'.format(label))
        if os.path.exists(file_path):
            print("File {0} already exists. Skipping...".format(file_path))
            continue
        print("working on mem{0:02d}".format(m + 1) + "...")
        
        t_10 = []; u_10 = []; v_10 = []; q_10 = []; swdn = []; lwdn = []; prec=[]

        for yr in range(int(year1), int(year2)+1):
            t_10.append(xr.open_dataset(JRA55_forcing_path + '/mem{0:02d}_JRA.v1.5_t_10_'.format(m + 1) + str(yr) + '.nc').t_10.values[:, i, j][:, 0])
            u_10.append(xr.open_dataset(JRA55_forcing_path + '/mem{0:02d}_JRA.v1.5_u_10_'.format(m + 1) + str(yr) + '.nc').u_10.values[:, i, j][:, 0])
            v_10.append(xr.open_dataset(JRA55_forcing_path + '/mem{0:02d}_JRA.v1.5_v_10_'.format(m + 1) + str(yr) + '.nc').v_10.values[:, i, j][:, 0])
            q_10.append(xr.open_dataset(JRA55_forcing_path + '/mem{0:02d}_JRA.v1.5_q_10_'.format(m + 1) + str(yr) + '.nc').q_10.values[:, i, j][:, 0])
            swdn.append(xr.open_dataset(JRA55_forcing_path + '/mem{0:02d}_JRA.v1.5_swdn_'.format(m + 1) + str(yr) + '.nc').swdn.values[:, i, j][:, 0])
            lwdn.append(xr.open_dataset(JRA55_forcing_path + '/mem{0:02d}_JRA.v1.5_lwdn_'.format(m + 1) + str(yr) + '.nc').lwdn.values[:, i, j][:, 0])
            prec.append(xr.open_dataset(JRA55_forcing_path + '/mem{0:02d}_JRA.v1.5_prec_'.format(m + 1) + str(yr) + '.nc').prec.values[:, i, j][:, 0])

        interp_t10 = np.interp(ip_times, jra_times, t_10)
        interp_u10 = np.interp(ip_times, jra_times, u_10)
        interp_v10 = np.interp(ip_times, jra_times, v_10)
        interp_q10 = np.interp(ip_times, jra_times, q_10) 
        interp_swdn = np.interp(ip_times, jra_times, swdn)
        interp_lwdn = np.interp(ip_times, jra_times, lwdn)
        interp_prec = np.interp(ip_times, jra_times, prec)

        with open(file_path, 'w') as new_file:
            new_file.writelines('#DSWSFC      DLWSFC    WNDU10    WNDV10    TEMP2M    SPECHUM    PRECIP\n')
            new_file.writelines('#W/m**2      W/m**2    m/s       m/s       K         Kg/Kg      kg/m**2/s\n')
            for l in range(interp_t10.shape[0]):
                new_file.writelines('{0:10.5f} {1:10.5f} {2:10.5f} {3:10.5f} {4:10.5f} {5:10.5f} {6:10.8f}\n'.format(
                    interp_swdn[l], interp_lwdn[l], interp_u10[l], interp_v10[l], interp_t10[l], interp_q10[l], interp_prec[l]))
        
print('Done creating ATM_FORCING files for {0}'.format(year) + ' for all points of interest!')

# a = abs(lon-lon_pt) + abs(lat-lat_pt)
# i,j = np.unravel_index(a.argmin(),a.shape)
# # i = 345
# # j = 118
# print(lat[i,j])
# print(lon[i,j])
# T = data.T.values[:,i,j]
# S = data.S.values[:,i,j]
# hblt = data.hblt.values[:,i,j]
# U = data.U.values[:,i,j]
# V = data.V.values[:,i,j]
# dhdx = data.dhdx.values[:,i,j]
# dhdy = data.dhdy.values[:,i,j]
# qdp = data.qdp.values[:,i,j]
# data.close()

# #hblt = np.array([6.1535325050354,6.48408842086792,11.9306182861328,21.7429275512695,38.2438545227051,122.865844726562,75.3419570922852,55.5070915222168,39.6110305786133,27.3990707397461,9.60380172729492,6.30194234848023])
# #T = np.array([-1.74758059,-1.75429083,-1.75728937,-1.76056446,-1.75164911,-1.73909664,-1.38452551,-0.48927524,-0.20830296,-1.37489022,-1.72806866,-1.73854388])
# if pert_ocean:
#   for m in range(0,mems.shape[0]):
#     np.random.seed(m*2)
#     hold = scipy.stats.truncnorm.rvs(a=-10,b=10,loc=0,scale=0.1,size=(6,12))
#     hold_T = T + hold[0]
#     hold_S = S + hold[1]
#     hold_hblt = hblt + hold[2]
#     hold_U = U + hold[3]
#     hold_V = V + hold[4]
#     hold_qdp = qdp + hold[5]
#     label = str(m+1).zfill(4)
#     print('Working on => {0}'.format(label))
#     new_file = open('NEW_OCN_FORCING_PERT_{0}.txt'.format(label),'w')
#     new_file.writelines('{0:10.13f} {1:10.13f} {2:10.13f}\n    {3:10.13f} {4:10.13f} {5:10.13f}\n    {6:10.13f} {7:10.13f} {8:10.13f}\n    {9:10.13f} {10:10.13f} {11:10.13f}\n'.format(hold_T[0],hold_T[1],hold_T[2],hold_T[3],hold_T[4],hold_T[5],hold_T[6],hold_T[7],hold_T[8],hold_T[9],hold_T[10],hold_T[11]))
#     new_file.writelines('{0:10.13f} {1:10.13f} {2:10.13f}\n    {3:10.13f} {4:10.13f} {5:10.13f}\n    {6:10.13f} {7:10.13f} {8:10.13f}\n    {9:10.13f} {10:10.13f} {11:10.13f}\n'.format(hold_S[0],hold_S[1],hold_S[2],hold_S[3],hold_S[4],hold_S[5],hold_S[6],hold_S[7],hold_S[8],hold_S[9],hold_S[10],hold_S[11]))
#     new_file.writelines('{0:10.13f} {1:10.13f} {2:10.13f}\n    {3:10.13f} {4:10.13f} {5:10.13f}\n    {6:10.13f} {7:10.13f} {8:10.13f}\n    {9:10.13f} {10:10.13f} {11:10.13f}\n'.format(hold_hblt[0],hold_hblt[1],hold_hblt[2],hold_hblt[3],hold_hblt[4],hold_hblt[5],hold_hblt[6],hold_hblt[7],hold_hblt[8],hold_hblt[9],hold_hblt[10],hold_hblt[11]))
#     new_file.writelines('{0:10.13f} {1:10.13f} {2:10.13f}\n    {3:10.13f} {4:10.13f} {5:10.13f}\n    {6:10.13f} {7:10.13f} {8:10.13f}\n    {9:10.13f} {10:10.13f} {11:10.13f}\n'.format(hold_U[0],hold_U[1],hold_U[2],hold_U[3],hold_U[4],hold_U[5],hold_U[6],hold_U[7],hold_U[8],hold_U[9],hold_U[10],hold_U[11]))
#     new_file.writelines('{0:10.13f} {1:10.13f} {2:10.13f}\n    {3:10.13f} {4:10.13f} {5:10.13f}\n    {6:10.13f} {7:10.13f} {8:10.13f}\n    {9:10.13f} {10:10.13f} {11:10.13f}\n'.format(hold_V[0],hold_V[1],hold_V[2],hold_V[3],hold_V[4],hold_V[5],hold_V[6],hold_V[7],hold_V[8],hold_V[9],hold_V[10],hold_V[11]))
#     new_file.writelines('{0:10.13f} {1:10.13f} {2:10.13f}\n    {3:10.13f} {4:10.13f} {5:10.13f}\n    {6:10.13f} {7:10.13f} {8:10.13f}\n    {9:10.13f} {10:10.13f} {11:10.13f}\n'.format(dhdx[0],dhdx[1],dhdx[2],dhdx[3],dhdx[4],dhdx[5],dhdx[6],dhdx[7],dhdx[8],dhdx[9],dhdx[10],dhdx[11]))
#     new_file.writelines('{0:10.13f} {1:10.13f} {2:10.13f}\n    {3:10.13f} {4:10.13f} {5:10.13f}\n    {6:10.13f} {7:10.13f} {8:10.13f}\n    {9:10.13f} {10:10.13f} {11:10.13f}\n'.format(dhdy[0],dhdy[1],dhdy[2],dhdy[3],dhdy[4],dhdy[5],dhdy[6],dhdy[7],dhdy[8],dhdy[9],dhdy[10],dhdy[11]))
#     new_file.writelines('{0:10.13f} {1:10.13f} {2:10.13f}\n    {3:10.13f} {4:10.13f} {5:10.13f}\n    {6:10.13f} {7:10.13f} {8:10.13f}\n    {9:10.13f} {10:10.13f} {11:10.13f}\n'.format(hold_qdp[0],hold_qdp[1],hold_qdp[2],hold_qdp[3],hold_qdp[4],hold_qdp[5],hold_qdp[6],hold_qdp[7],hold_qdp[8],hold_qdp[9],hold_qdp[10],hold_qdp[11]))
#     new_file.close()
# else:  
#   for m in range(0,mems.shape[0]):
#     label = str(m+1).zfill(4)
#     print('Working on => {0}'.format(label))
#     new_file = open('NEW_OCN_FORCING_{0}.txt'.format(label),'w')
#     new_file.writelines('{0:10.13f} {1:10.13f} {2:10.13f}\n    {3:10.13f} {4:10.13f} {5:10.13f}\n    {6:10.13f} {7:10.13f} {8:10.13f}\n    {9:10.13f} {10:10.13f} {11:10.13f}\n'.format(T[0],T[1],T[2],T[3],T[4],T[5],T[6],T[7],T[8],T[9],T[10],T[11]))
#     new_file.writelines('{0:10.13f} {1:10.13f} {2:10.13f}\n    {3:10.13f} {4:10.13f} {5:10.13f}\n    {6:10.13f} {7:10.13f} {8:10.13f}\n    {9:10.13f} {10:10.13f} {11:10.13f}\n'.format(S[0],S[1],S[2],S[3],S[4],S[5],S[6],S[7],S[8],S[9],S[10],S[11]))
#     new_file.writelines('{0:10.13f} {1:10.13f} {2:10.13f}\n    {3:10.13f} {4:10.13f} {5:10.13f}\n    {6:10.13f} {7:10.13f} {8:10.13f}\n    {9:10.13f} {10:10.13f} {11:10.13f}\n'.format(hblt[0],hblt[1],hblt[2],hblt[3],hblt[4],hblt[5],hblt[6],hblt[7],hblt[8],hblt[9],hblt[10],hblt[11]))
#     new_file.writelines('{0:10.13f} {1:10.13f} {2:10.13f}\n    {3:10.13f} {4:10.13f} {5:10.13f}\n    {6:10.13f} {7:10.13f} {8:10.13f}\n    {9:10.13f} {10:10.13f} {11:10.13f}\n'.format(U[0],U[1],U[2],U[3],U[4],U[5],U[6],U[7],U[8],U[9],U[10],U[11]))
#     new_file.writelines('{0:10.13f} {1:10.13f} {2:10.13f}\n    {3:10.13f} {4:10.13f} {5:10.13f}\n    {6:10.13f} {7:10.13f} {8:10.13f}\n    {9:10.13f} {10:10.13f} {11:10.13f}\n'.format(V[0],V[1],V[2],V[3],V[4],V[5],V[6],V[7],V[8],V[9],V[10],V[11]))
#     new_file.writelines('{0:10.13f} {1:10.13f} {2:10.13f}\n    {3:10.13f} {4:10.13f} {5:10.13f}\n    {6:10.13f} {7:10.13f} {8:10.13f}\n    {9:10.13f} {10:10.13f} {11:10.13f}\n'.format(dhdx[0],dhdx[1],dhdx[2],dhdx[3],dhdx[4],dhdx[5],dhdx[6],dhdx[7],dhdx[8],dhdx[9],dhdx[10],dhdx[11]))
#     new_file.writelines('{0:10.13f} {1:10.13f} {2:10.13f}\n    {3:10.13f} {4:10.13f} {5:10.13f}\n    {6:10.13f} {7:10.13f} {8:10.13f}\n    {9:10.13f} {10:10.13f} {11:10.13f}\n'.format(dhdy[0],dhdy[1],dhdy[2],dhdy[3],dhdy[4],dhdy[5],dhdy[6],dhdy[7],dhdy[8],dhdy[9],dhdy[10],dhdy[11]))
#     new_file.writelines('{0:10.13f} {1:10.13f} {2:10.13f}\n    {3:10.13f} {4:10.13f} {5:10.13f}\n    {6:10.13f} {7:10.13f} {8:10.13f}\n    {9:10.13f} {10:10.13f} {11:10.13f}\n'.format(qdp[0],qdp[1],qdp[2],qdp[3],qdp[4],qdp[5],qdp[6],qdp[7],qdp[8],qdp[9],qdp[10],qdp[11]))
#     new_file.close() 

