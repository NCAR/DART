
#imports (standard ones first, then the ones I wrote)
import numpy as np #for data structures
from scipy.io import netcdf #to read DART data
from scipy.interpolate import griddata #for interpolation
from scipy import interpolate
import subprocess #for moving files
import sys
import cPickle as pickle
import datetime, time #for keeping track of GITM files
import math #for pi
import matplotlib
matplotlib.use('Agg') #so as not to need a display or X-server
import matplotlib.pyplot as plt #for plotting
from matplotlib.patches import Polygon
from mpl_toolkits.mplot3d import Axes3D #for 3d plotting

import read_gitm_bin #module I wrote for reading gitm binary files
import read_gps_txt #module I wrote for reading gps text files
import plot_tec_lib #module I wrote for plotting tec

def tecify(time, Lon, Lat, Alt, Eds): 
    """
    compute VTEC from electron density (Eds)
    
    INPUTS:
    - time[ti] - np.array of "datetime" times
    - Lon[n] - np.array of lons
    - Lat[t] - np.array of lats
    - Alt[a] - np.array of alts
    - Eds[n,t,a,ti] - np.array over [lon,lat,alt,time]
    
    OUTPUTS:
    - Vtec_r[n,t,ti] - np.array over [lon,lat,time]
    """
    dh = np.diff(Alt) #Alt(i+1)-Alt(i)
    H = np.diag(np.ones(len(Alt), float)) + np.diag(np.ones(len(Alt)-1, float),k=1)
    H = H[0:-1,:] #np.dot(H,I) will result in sum of neighbors = [I(2)+I(1), I(3)+I(2), ..., I(50)+I(49)]
    Vtec_r=np.zeros( (len(Lon), len(Lat), len(time)), float )
    for ti in range(len(time)): #I just ran it for 1 day
        Vtec = np.zeros( (len(Lon),len(Lat)), float )
        for n in range(len(Lon)):
            for t in range(len(Lat)):
                Vtec_r[n,t,ti] = np.dot(dh, np.dot(H,Eds[n,t,:,ti]))/2.0 * 10.0**-16 #tec=(A2-A1)*(f2+f1)/2 - trapezoidal rule, 10^-16 is the TECU
    return Vtec_r



#contains - I use this keyword as separation of function defenitions from actual calls. this is opposite of Fortran, but oh well
fn_gitm = sys.argv[1] #where the truth  files are
fn_gps = sys.argv[2] #where the gps observation text file is
sim = sys.argv[3] in ['true', 'True', '.true.', 'T', 't', '1', True] #was this gps file simulated (True) or real (False)? eg sim = 'True' for simulated. Affects what gets comparred to what below (search sim) and what gets plotted

#read GITM binary files
timeT1 = datetime.datetime(2002, 12, 1, 0, 5, 0) #when did you start you truth simulation? and note that GITM produces first file AFTER one full step (ie if you started at 0 minutes and asked it to write every 5 minutes, the first file will show up at 5 min)
timeT = np.array([ timeT1 + datetime.timedelta(minutes=m) for m in range(0, 30, 5)]) #assumes you ran gitm for 30 minutes with step of 5 minutes
print 'Reading GITM binary files... '
#to save time in subsequent reads, I only read truth data once and pickle it (into fn+'.p') and use the pickled file everytime afterwards
try: 
    f = open( fn_gitm+'.p', 'rb' ) 
    print 'pickled file found, reading it'
    (LonT, LatT, AltT, EdsT) = pickle.load(f)
    f.close()
except:
    print 'pickled file not found, making one'
    (LonT, LatT, AltT, EdsT) = read_gitm_bin.read(timeT, fn_gitm)
    f = open( fn_gitm+'.p', 'wb' ) 
    pickle.dump((LonT, LatT, AltT, EdsT) , f)
    f.close()

print 'Done w/ GITM binary files.'
LonT=LonT*180/math.pi #convert to degrees
LatT=LatT*180/math.pi #convert to degrees
VtecT = tecify(timeT, LonT, LatT, AltT, EdsT)

#read DART netcdf files
print 'Reading DART netcdf files... '
f = netcdf.netcdf_file('../work/preassim.nc', 'r')
#f.dimensions
#for d in dir(f): print d
#for d in f.variables: print d
timeE1 = f.variables['time'].data #timeE.units 'days since 1601-01-01 00:00:00'
#convert dart time to unix epoch time (seconds since 1970 1 1)
# to do that, add the original dart base of 1601 and subtract the unix base of 1970
tdiff = (datetime.datetime(1970,1,1)-datetime.datetime(1601,1,1)).total_seconds()
timeE = timeE1*86400 - tdiff 
timeE = np.array([ datetime.datetime.utcfromtimestamp(t) for t in timeE ]) #convert unix epoch to datetime
#tt = (timeT - datetime.datetime(1970,1,1)).total_seconds() #convert timeT to unix epoch time

LonE = f.variables['LON'].data
LatE = f.variables['LAT'].data
AltE = f.variables['ALT'].data
f.variables['ie_IDensityS'].dimensions
type(f.variables['ie_IDensityS'].dimensions)
#f.variables['ie_IDensityS'].T.dimensions
EdsEr = f.variables['ie_IDensityS'].data[:,0,:,:,:].T #np.array over [time,copy,alt,lat,lon], where copy=0 means ensemble mean, 1 ens. var, possibly ens. members themself (if you change num_output_state_members in input.nml to nonzero values) and then -2 is inflation_value at every step (see section 9 in DART tutorial and inf_initial in input.nml) and -1 is inflation_sdandard_deviation (see inf_sd_initial in input.nml)
#Transpose is needed because dart writes it in [time,alt,lat,lon], whereas I want [lon,lat,alt,time], so only a reversal (=transpose) is needed
EdssdEr = f.variables['ie_IDensityS'].data[:,1,:,:,:].T #standard deviation
f107Er = f.variables['f107'].data[:,0,0] #[time,copy,wavelength] - mean, just 1 wavelength for now
f107sdEr = f.variables['f107'].data[:,1,0] #standard deviation
f.close()

f = netcdf.netcdf_file('../work/analysis.nc', 'r') 
EdsEo = f.variables['ie_IDensityS'].data[:,0,:,:,:].T #mean
EdssdEo = f.variables['ie_IDensityS'].data[:,1,:,:,:].T #standard deviation
f107Eo = f.variables['f107'].data[:,0,0] #mean
f107sdEo = f.variables['f107'].data[:,1,0] #sd
f.close()

f = netcdf.netcdf_file('../work/obs_diag_output.nc', 'r') 
#for d in dir(f): print d
#for d in f.variables: print d
#for d in range(21):  print d, ''.join(f.variables['CopyMetaData'][d]) 
obs_time = f.variables['time'].data
obs_time = obs_time*86400 - (datetime.datetime(1970,1,1)-datetime.datetime(1601,1,1)).total_seconds()
obs_time = np.array([ datetime.datetime.utcfromtimestamp(t) for t in obs_time ]) #convert unix e
obs_vposs_r = f.variables['GND_GPS_VTEC_guess'][:,0,3,0]
obs_vused_r = f.variables['GND_GPS_VTEC_guess'][:,1,3,0]
obs_vrmse_r = f.variables['GND_GPS_VTEC_guess'][:,6,3,0]
obs_vbias_r = f.variables['GND_GPS_VTEC_guess'][:,7,3,0]
obs_vspread_r = f.variables['GND_GPS_VTEC_guess'][:,8,3,0]
obs_vtotspread_r = f.variables['GND_GPS_VTEC_guess'][:,9,3,0]
obs_vbadqc_r = f.variables['GND_GPS_VTEC_guess'][:,10,3,0]
obs_vtruth_r = f.variables['GND_GPS_VTEC_guess'][:,11,3,0]
obs_vensm_r = f.variables['GND_GPS_VTEC_guess'][:,12,3,0]

obs_vposs_o = f.variables['GND_GPS_VTEC_analy'][:,0,3,0]
obs_vused_o = f.variables['GND_GPS_VTEC_analy'][:,1,3,0]
obs_vrmse_o = f.variables['GND_GPS_VTEC_analy'][:,6,3,0]
obs_vbias_o = f.variables['GND_GPS_VTEC_analy'][:,7,3,0]
obs_vspread_o = f.variables['GND_GPS_VTEC_analy'][:,8,3,0]
obs_vtotspread_o = f.variables['GND_GPS_VTEC_analy'][:,9,3,0]
obs_vbadqc_o = f.variables['GND_GPS_VTEC_analy'][:,10,3,0]
obs_vtruth_o = f.variables['GND_GPS_VTEC_analy'][:,11,3,0]
obs_vensm_o = f.variables['GND_GPS_VTEC_analy'][:,12,3,0]
f.close()
print 'Done w/ DART netcdf files.'



VtecE = tecify(timeE, LonE, LatE, AltE, EdsEo)


#load gps vtec data
timeD, LonD, LatD, VtecD, VtecsdD = read_gps_txt.read(fn_gps)


if sim:  
    timeD1= timeD #if sim data, no need to shift time
else:
    timeD1 = timeD + datetime.timedelta(seconds=150) #timeD is shifted by 2.5 minutes so it matches timeT (only for plotting purposes)

tcommon = sorted(set(timeT) & set(timeE) & set(timeD1)) #plot only at common times

#sequence of 2dplots of only data
#execfile('plot_data_2d.py')

#sequence of 3dplots of data with eakf estimates and with truth
#execfile('plot_gps_3d.py')

#fit a surface to gps data at every timestep
VtecDi = np.zeros( (len(LonE), len(LatE), len(tcommon) ), float)
for ti in range(len(tcommon)): 
    x = LonD[timeD1==tcommon[ti]]
    y = LatD[timeD1==tcommon[ti]]
    z = VtecD[timeD1==tcommon[ti]]
    VtecDi[:,:,ti] = griddata((x, y), z, np.meshgrid(LonE, LatE), method='linear', fill_value=0).T #BEWARE of fill_value=0!!!!!, Transpose is needed because meshgrid creates Lon that varies inside a row instead of column!
    
#more legit interpolation way (including time), but doesn't work yet
#VtecDi = np.zeros( (len(LonE), len(LatE), len(tcommon) ), float)
#loni, lati = np.meshgrid(LonE, LatE)
#for ti in range(len(timeE)): 
#    lon = LonD[abs(timeD-timeE[ti]) < datetime.timedelta(minutes=30)]
#    lat = LatD[abs(timeD-timeE[ti]) < datetime.timedelta(minutes=30)]
#    time = time.mktime(timeD[abs(timeD-timeE[ti]) < datetime.timedelta(minutes=30)].timetuple())
#    vtec = VtecD[abs(timeD-timeE[ti]) < datetime.timedelta(minutes=30)]
#    timei = np.empty_like(loni, object)
#    timei[:] = timeE[ti]
#    griddata((lon, lat, time), vtec, (loni, lati, timei), method='linear') # fill_value=0
#    VtecDi[:,:,ti] = griddata((lon, lat, time), vtec, (loni, lati, timei), method='linear') # fill_value=0

#sequence of 3dplots
#execfile('plot_gps_3d_interp.py')

#ensemble - data
VtecED = np.zeros( (len(LonE), len(LatE), len(tcommon) ), float)
for ti in range(len(tcommon)): 
    if sim: #simulated data, want to compare est with truth  
        VtecED[:,:,ti] = abs(VtecE[:,:,np.where(timeE==list(tcommon)[ti])[0][0]] - VtecT[:,:,np.where(timeT==list(tcommon)[ti])[0][0]])
    else: #real data, want to compare est with real
        VtecED[:,:,ti] = abs(VtecE[:,:,np.where(timeE==list(tcommon)[ti])[0][0]] - VtecDi[:,:,ti])

#interpolate the diff so it looks nicer on plots
res = 1.0 #lat,lon resolution to which interpolate the data
loni = np.arange(min(LonE)+res,max(LonE),res) #what lons do you want to interpolate it to
lati = np.arange(min(LatE)+res,max(LatE),res) #what lons do you want to interpolate it to
VtecEDi = np.zeros( (len(loni), len(lati), len(tcommon) ), float)
for ti in range(len(tcommon)): 
    f = interpolate.RectBivariateSpline(LonE,LatE,VtecED[:,:,ti])
    for t in range(len(lati)):
        for n in range(len(loni)):
            VtecEDi[n,t,ti] = f( loni[n], lati[t] ) 

#memory was an issue, so here is how to check:
# for i in dir():
#     try: 
#         print (i, eval(i).nbytes )
#     except: 
#         print (i, sys.getsizeof(eval(i)) )

#sequence of 3dplots
#execfile('plot_diff_2d.py')
#execfile('plot_diff_2d_f107.py')

lona=13 #index of Lon closest to Ann Arbor
lata=13 #index of Lat closest to Ann Arbor
alta=34 #index of Alt I chose
print 'LonT[lona], LatT[lata], AltT[alta], EdsT[lona,lata,alta,0], EdsEr[lona,lata,alta,0], EdsEo[lona,lata,alta,0]', LonT[lona], LatT[lata], AltT[alta], EdsT[lona,lata,alta,0], EdsEr[lona,lata,alta,0], EdsEo[lona,lata,alta,0]
ind = obs_vrmse_r>0

if sim: #simulated data, want to compare est with truth  
#f107 vs time plot
    plot_tec_lib.time_plot(timeT,[148.0]*len(timeT),'$F_{10.7}$ measured', timeE,f107Er,f107Eo,f107sdEr,f107sdEo,'$F_{10.7}$ estimated','$F_{10.7}$ estimated $\pm$ SD', (0,24),(100,300),'$F_{10.7}$  [SFU]', 'f107_00')
    
#Eds[13,13,34,:] vs time plot
    plot_tec_lib.time_plot(timeT, EdsT[lona,lata,alta,:], '$D_{e,aa}$ measured', timeE, EdsEr[lona,lata,alta,:], EdsEo[lona,lata,alta,:], EdssdEr[lona,lata,alta,:], EdssdEo[lona,lata,alta,:], '$D_{e,aa}$ estimated', '$D_{e,aa}$ estimated $\pm$ SD', (0,24), (None,None), '$D_{e,AA}$  [$e$ $m^{-3}$]', 'daa00')
#obs_diag rmse vs time plot
    
    plot_tec_lib.time_plot(obs_time[ind], obs_vtruth_r[ind],'ignore', obs_time[ind], obs_vrmse_r[ind], obs_vrmse_o[ind], obs_vspread_r[ind], obs_vspread_o[ind], 'Average VTEC error', 'Average VTEC spread ($\pm$ SD)', (0,24),(0,10),'Average VTEC error  [TECU]', 'obs_rmse00')
    
else:
#f107 vs time plot
    plot_tec_lib.time_plot(timeT,[148.0]*len(timeT),'ignore', timeE,f107Er,f107Eo,f107sdEr,f107sdEo,'$F_{10.7}$ estimated','$F_{10.7}$ estimated $\pm$ SD', (0,24),(100,300),'$F_{10.7}$  [SFU]', 'f107_00')
    
#Eds[13,13,34,:] vs time plot
    plot_tec_lib.time_plot(timeT, EdsT[lona,lata,alta,:], 'ignore', timeE, EdsEr[lona,lata,alta,:], EdsEo[lona,lata,alta,:], EdssdEr[lona,lata,alta,:], EdssdEo[lona,lata,alta,:], '$D_{e,aa}$ estimated', '$D_{e,aa}$ estimated $\pm$ SD', (0,24), (None,None), '$D_{e,AA}$  [$e$ $m^{-3}$]', 'daa00')
    
#obs_diag rmse vs time plot
    plot_tec_lib.time_plot(obs_time[ind], obs_vtruth_r[ind],'ignore', obs_time[ind], obs_vrmse_r[ind], obs_vrmse_o[ind], obs_vspread_r[ind], obs_vspread_o[ind], 'Average VTEC error', 'Average VTEC spread ($\pm$ SD)', (0,24),(0,16),'Average VTEC error  [TECU]', 'obs_rmse00')


#Eds[13,13,:,:] vs time, alt subplots
#plot_tec_lib.time_plot_alt(timeT, EdsT[lona,lata,:,:], AltT, '$D_{e,aa}$ measured', timeE, EdsEr[lona,lata,:,:], EdsEo[lona,lata,:,:], EdssdEr[lona,lata,:,:], EdssdEo[lona,lata,:,:], '$D_{e,aa}$ estimated', '$D_{e,aa}$ estimated $\pm$ SD', (0,24), (1e8,1e12), '$D_{e,AA}$  [$e$ $m^{-3}$]', 'daa_alt00')

#error vs time plot
#plot_tec_lib.time_plot(timeT,[1e32.0]*len(timeT),'ignore', np.array(list(tcommon)),VtecED[:,:,ti],f107Eo,f107sdEr,f107sdEo,'$F_{10.7}$ estimated','$F_{10.7}$ estimated $\pm$ SD', (0,24),(100,300),'$F_{10.7}$  [SFU]', 'f107_00')


#obs_diag mean vs time plot
ind = obs_vrmse_r>0
plot_tec_lib.time_plot(obs_time[ind], obs_vtruth_r[ind],'VT', obs_time[ind], obs_vensm_r[ind], obs_vensm_o[ind], obs_vspread_r[ind], obs_vspread_o[ind], 'ensm', 'spread', (None,None),(None,None),'Vensm', 'obs_ens00')


#write interpolated Truth data into a dart obs file
#execfile('write_gps_txt.py')

# DART $Id$
# from Alexey Morozov

# <next few lines under version control, do not edit>
# $URL$
# $Revision$
# $Date$
