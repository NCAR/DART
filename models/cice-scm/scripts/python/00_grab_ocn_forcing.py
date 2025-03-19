import xarray as xr
import numpy as np
import os
# import scipy.stats

# pert_ocean = False

#lat_pt = 88.3935
#lon_pt = 80.51758

# lat_pt = 75.53822
# lon_pt = 174.4456
LatsOfInterest = [75.53822, 81, 88, 75]
LonsOfInterest = [174.4456, 358, 0, 40]
PointsOfInterest = ['Sib_Chuk', 'CoastalCanada', 'CentralArctic', 'Barents']
mems = np.arange(1,31)
data = xr.open_dataset('/glade/work/jrotondo/icepack_da/initialization_data/forcings/ISPOL_2004/pop_frc.b.e21.BW1850.f09_g17.CMIP6-piControl.001.190514.nc')
lon = data.xc.values
lat = data.yc.values

for p in range(0, len(LatsOfInterest)):
    lat_pt = LatsOfInterest[p]
    lon_pt = LonsOfInterest[p]
    a = abs(lon-lon_pt) + abs(lat-lat_pt)
    i,j = np.unravel_index(a.argmin(),a.shape)

    # make a data directory for ocean forcing data
    PoI_path = '/glade/work/jrotondo/icepack_da/initialization_data/forcings/ISPOL_2004/'+PointsOfInterest[p]
    if not os.path.exists(PoI_path):
        os.makedirs(PoI_path)

    T = data.T.values[:,i,j]
    S = data.S.values[:,i,j]
    hblt = data.hblt.values[:,i,j]
    U = data.U.values[:,i,j]
    V = data.V.values[:,i,j]
    dhdx = data.dhdx.values[:,i,j]
    dhdy = data.dhdy.values[:,i,j]
    qdp = data.qdp.values[:,i,j]

    for m in range(0,mems.shape[0]):
        label = str(m+1).zfill(4)
        print('Working on => ' + str(lat[i,j]) + 'N, '+str(lon[i,j])+'E, member {0}'.format(label))
        new_file = open(PoI_path + '/OCN_FORCING_{0}.txt'.format(label),'w')
        new_file.writelines('{0:10.13f} {1:10.13f} {2:10.13f}\n    {3:10.13f} {4:10.13f} {5:10.13f}\n    {6:10.13f} {7:10.13f} {8:10.13f}\n    {9:10.13f} {10:10.13f} {11:10.13f}\n'.format(T[0],T[1],T[2],T[3],T[4],T[5],T[6],T[7],T[8],T[9],T[10],T[11]))
        new_file.writelines('{0:10.13f} {1:10.13f} {2:10.13f}\n    {3:10.13f} {4:10.13f} {5:10.13f}\n    {6:10.13f} {7:10.13f} {8:10.13f}\n    {9:10.13f} {10:10.13f} {11:10.13f}\n'.format(S[0],S[1],S[2],S[3],S[4],S[5],S[6],S[7],S[8],S[9],S[10],S[11]))
        new_file.writelines('{0:10.13f} {1:10.13f} {2:10.13f}\n    {3:10.13f} {4:10.13f} {5:10.13f}\n    {6:10.13f} {7:10.13f} {8:10.13f}\n    {9:10.13f} {10:10.13f} {11:10.13f}\n'.format(hblt[0],hblt[1],hblt[2],hblt[3],hblt[4],hblt[5],hblt[6],hblt[7],hblt[8],hblt[9],hblt[10],hblt[11]))
        new_file.writelines('{0:10.13f} {1:10.13f} {2:10.13f}\n    {3:10.13f} {4:10.13f} {5:10.13f}\n    {6:10.13f} {7:10.13f} {8:10.13f}\n    {9:10.13f} {10:10.13f} {11:10.13f}\n'.format(U[0],U[1],U[2],U[3],U[4],U[5],U[6],U[7],U[8],U[9],U[10],U[11]))
        new_file.writelines('{0:10.13f} {1:10.13f} {2:10.13f}\n    {3:10.13f} {4:10.13f} {5:10.13f}\n    {6:10.13f} {7:10.13f} {8:10.13f}\n    {9:10.13f} {10:10.13f} {11:10.13f}\n'.format(V[0],V[1],V[2],V[3],V[4],V[5],V[6],V[7],V[8],V[9],V[10],V[11]))
        new_file.writelines('{0:10.13f} {1:10.13f} {2:10.13f}\n    {3:10.13f} {4:10.13f} {5:10.13f}\n    {6:10.13f} {7:10.13f} {8:10.13f}\n    {9:10.13f} {10:10.13f} {11:10.13f}\n'.format(dhdx[0],dhdx[1],dhdx[2],dhdx[3],dhdx[4],dhdx[5],dhdx[6],dhdx[7],dhdx[8],dhdx[9],dhdx[10],dhdx[11]))
        new_file.writelines('{0:10.13f} {1:10.13f} {2:10.13f}\n    {3:10.13f} {4:10.13f} {5:10.13f}\n    {6:10.13f} {7:10.13f} {8:10.13f}\n    {9:10.13f} {10:10.13f} {11:10.13f}\n'.format(dhdy[0],dhdy[1],dhdy[2],dhdy[3],dhdy[4],dhdy[5],dhdy[6],dhdy[7],dhdy[8],dhdy[9],dhdy[10],dhdy[11]))
        new_file.writelines('{0:10.13f} {1:10.13f} {2:10.13f}\n    {3:10.13f} {4:10.13f} {5:10.13f}\n    {6:10.13f} {7:10.13f} {8:10.13f}\n    {9:10.13f} {10:10.13f} {11:10.13f}\n'.format(qdp[0],qdp[1],qdp[2],qdp[3],qdp[4],qdp[5],qdp[6],qdp[7],qdp[8],qdp[9],qdp[10],qdp[11]))
        new_file.close() 


data.close()




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

