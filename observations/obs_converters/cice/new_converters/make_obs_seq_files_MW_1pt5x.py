import pandas as pd
from datetime import datetime,timedelta
import xarray as xr
import numpy as np
import glob
from itertools import chain
import os
from datetime import timedelta
from os import path
import sys
# from scipy.interpolate import CubicSpline
# from mpl_toolkits.basemap import Basemap

def lon_lat_to_cartesian(lon, lat, R = 1):
    """
    calculates lon, lat coordinates of a point on a sphere with
    radius R
    """
    lon_r = np.radians(lon)
    lat_r = np.radians(lat)

    x =  R * np.cos(lat_r) * np.cos(lon_r)
    y = R * np.cos(lat_r) * np.sin(lon_r)
    z = R * np.sin(lat_r)
    return x,y,z

#data = xr.open_dataset('seaice_obs_errs.nc')
#sic_vals = data.sic_vals.values
#sic_err = data.sic_err.values
#sit_vals = data.sit_vals.values
#sit_err = data.sit_err.values
#max_sit_val = sit_vals[-1]
#max_sit_err = sit_err[-1]
#data.close()
#sic_err_fun = CubicSpline(sic_vals,sic_err)
#sit_err_fun = CubicSpline(sit_vals,sit_err)

datelist = pd.date_range(start=datetime(2019,1,2),end=datetime(2020,1,1)).to_pydatetime()
ens_truth = 11

print('random truth is: '+str(ens_truth))

print('Opening tracks file..')
track_data = xr.open_dataset('simulated_tracks_2019.nc')
track_yearmonday = track_data.year_month_day.values
print(track_yearmonday)
track_hours = np.array([int(str(d)[8:-2]) for d in track_data.dates.values])
print(track_hours)
#print('Opening seaice data...')
#seaice_data = xr.open_dataset('../spinupCase-CAM-2013-ENS-PARAMS-RERUN/spinupCase-CAM-2013-ENS-PARAMS-RERUN/run/output/0025/history_timeseries.nc')
#seaice_yearmonday = np.array([int(d.strftime('%Y%m%d')) for d in seaice_data.time.values])
#seaice_lat = seaice_data.TLAT.values
#seaice_lon = seaice_data.TLON.values


for t in range(0,datelist.shape[0]):
# for t in range(0,1):
  count = 0
  print('--------------------------------')
  total_lat = []
  total_lon = []
  print(datelist[t])
  curr_date = int(datelist[t].strftime('%Y%m%d'))
  check_date = int((datelist[t]-timedelta(days=1)).strftime('%Y%m%d'))
  year = (datelist[t]).strftime('%Y')
  mon = (datelist[t]).strftime('%m')
  day_old = (datelist[t]-timedelta(days=1)).strftime('%d')
  day = (datelist[t]).strftime('%d')
  
  year_old = (datelist[t]-timedelta(days=1)).strftime('%Y')
  mon_old = (datelist[t]-timedelta(days=1)).strftime('%m')
  day_old = (datelist[t]-timedelta(days=1)).strftime('%d')

  diff = datelist[t] - datetime(1601,1,1,0)
  comd = 'cat > test.sed << EOF\ns#{days}#'+str(diff.days)+'#g\ns#{seconds}#'+str(diff.seconds)+'#g\nEOF'
  os.system(comd)
  comd = 'sed -f test.sed input.nml.template > input.nml'
  os.system(comd)
  os.system('rm test.sed')
  #sys.exit()
  ind_track = np.where(track_yearmonday==check_date)[0]
  print(ind_track)
  
  if ind_track.shape[0] == 0:
    print('Dates could not be found in of the two files...')
    sys.exit()

  # print(track_hours)
  track_lat = track_data.latitude.values[ind_track[0]]
  track_lon = track_data.longitude.values[ind_track[0]]
  track_hrs = track_hours[ind_track]
  # print(track_hrs)
  print('Number of sat. points for that date:'+str(track_lat.shape[0]))
  ind = np.where(track_lat>30)[0]
  track_lat = track_lat[ind]
  track_lon = track_lon[ind]
  # track_hrs = track_hrs[ind]
  print('Number of sat. points after filtering for NH:'+str(track_lat.shape[0]))
  unique_hrs = np.unique(track_hrs)
  # This section of code cycles through the different hours in the observations. We don't have
  # this in our files, the hour written is 0000

  # for every unique hour 
  for p in range(0,unique_hrs.shape[0]):
    hold_lat = []
    hold_lon = []
    hold_conc = []
    hold_thic = []
    print('Working on hour:'+str(unique_hrs[p]))
    ind2 = np.where(track_hrs==unique_hrs[p])[0]
    hold_trk_lat = track_lat #[ind2]
    print(hold_trk_lat)
    hold_trk_lon = track_lon #[ind2]
    label = f'{unique_hrs[p]*3600:05}'
    save_date = 'obs_seq.in_'+str(check_date)+'_'+label
    #restart = '/glade/scratch/mollyw/free30/run/free30.cice_0024.r.{0}-{1}-{2}-{3}.nc'.format(year_old,mon_old,day_old, label)
    path = '/glade/campaign/univ/uwas0070/mwieringa/assimilation/cases/free30_temp/run/truth/mem_00'+str(ens_truth)+'/free30.cice_00'+str(ens_truth)+'.h.{0}-{1}-{2}.nc'.format(year_old,mon_old,day_old)
    print(path)
    if len(glob.glob(path)) == 0:
      print('There is no seaice data for '+str(curr_date)+' at hour '+str(label))
      continue
    seaice_data = xr.open_dataset(path)    
    seaice_lat = seaice_data.TLAT.values
    seaice_lon = seaice_data.TLON.values
    seaice_conc = seaice_data.aice.values[0,:,:]
    seaice_thick = seaice_data.hi.values[0,:,:]
    mask = seaice_data.tmask.values[:,:]
    seaice_data.close()
    for pp in range(0,hold_trk_lat.shape[0]):
      a = abs(seaice_lat-hold_trk_lat[pp])+abs(seaice_lon-hold_trk_lon[pp])
      #a = np.sqrt((seaice_lat-hold_trk_lat[pp])**2+(seaice_lon-hold_trk_lon[pp])**2)
      i,j = np.unravel_index(np.nanargmin(a),a.shape)
      #print('--------------------------')
      #print('Track points:'+str(track_lat[pp])+','+str(track_lon[pp]))
      #print('Seaice points:'+str(seaice_lat[i,j])+','+str(seaice_lon[i,j]))
      #print('--------------------------')
      
      if hold_trk_lat[pp] > 90. or hold_trk_lat[pp]<-90:
        print('STOP')
        sys.exit()
      if np.isnan(seaice_conc[i,j]) or mask[i,j]== 0:# or seaice_conc[i,j] == 0 or seaice_conc[i,j]<0.01:
        continue
      #elif seaice_conc[i,j] == 0.0 and mean_conc[i,j] == 0.0:
      #  continue
      #elif seaice_conc[i,j] >= 0.01 and mean_conc[i,j] == 0.0:
      #  continue
      elif hold_trk_lat[pp] < 45.:
        continue
      else:
        hold_lat.append(hold_trk_lat[pp])
        hold_lon.append(hold_trk_lon[pp])
        hold_conc.append(seaice_conc[i,j])
        hold_thic.append(seaice_thick[i,j])
    total_lat.append(hold_lat)
    total_lon.append(hold_lon)
    hold_lat = np.array(hold_lat)
    hold_lon = np.array(hold_lon)
    if hold_lat.shape[0] == 0 or hold_lon.shape[0] == 0:
      print('No obs over Arctic sea ice so skip time....')
      continue
    hold_conc = np.array(hold_conc)
    hold_thic = np.array(hold_thic)
    if unique_hrs[p] < 10:
      hour_label = '0'+str(unique_hrs[p])
    else:
      hour_label = str(unique_hrs[p])
    print('Number of sat. points after filtering over seaice:'+str(hold_lat.shape[0])) 
    num_obs = hold_lat.shape[0]*3
    comd = 'rm input.txt'
    os.system(comd)
    new_file = open('input.txt','w')
    new_file.writelines(str(num_obs)+'\n')
    #new_file.writelines('2\n')
    new_file.writelines('0\n')
    new_file.writelines('0\n')
    for l in range(0,hold_lat.shape[0]):
      count += 1
      new_file.writelines('0\n')
      new_file.writelines('12\n')
      new_file.writelines('-1\n')
      new_file.writelines('0\n')
      new_file.writelines(str(hold_lon[l])+'\n')
      new_file.writelines(str(hold_lat[l])+'\n')
      new_file.writelines(year+' '+mon+' '+day+' '+str(00)+' 00 00\n')
      #new_file.writelines('0.05\n')
      #error = sic_err_fun(hold_conc[l])
      if hold_conc[l] == 0:
       error = 0.075
      else:
       error = hold_conc[l]*0.225
      new_file.writelines(str(error**2)+'\n')

      count += 1
      new_file.writelines('1\n')
      new_file.writelines('11\n')
      new_file.writelines('-1\n')
      new_file.writelines('0\n') 
      if hold_lon[l] == np.degrees(5.628433346061153) and hold_lat[l] == np.degrees(1.102595383759034):
       print('STOP')
       sys.exit()   
      new_file.writelines(str(hold_lon[l])+'\n')
      new_file.writelines(str(hold_lat[l])+'\n')
      new_file.writelines(year+' '+mon+' '+day+' '+str(00)+' 00 00\n')
      new_file.writelines('0.005625\n')
      
      count += 1
      new_file.writelines('1\n')
      new_file.writelines('15\n')
      new_file.writelines('-1\n')
      new_file.writelines('0\n') 
      if hold_lon[l] == np.degrees(5.628433346061153) and hold_lat[l] == np.degrees(1.102595383759034):
       print('STOP')
       sys.exit()   
      new_file.writelines(str(hold_lon[l])+'\n')
      new_file.writelines(str(hold_lat[l])+'\n')
      new_file.writelines(year+' '+mon+' '+day+' '+str(00)+' 00 00\n')
      ##if hold_thic[l] > max_sit_val:
      ##  error = max_sit_err
      ##elif hold_thic[l] == 0.0:
      ##  error = sit_err_fun(hold_thic[l])
      ##else:
      ##  error = sit_err_fun(hold_thic[l]/hold_conc[l])
      # error = hold_thic[l]*0.1 
      new_file.writelines('0.0225\n')
      # new_file.writelines('0.00\n') 
    #sys.exit()  
    # new_file.writelines('-1\n')
    new_file.writelines(save_date)
    new_file.close()
    #comd = 'rm cice.r.nc'
    #os.system(comd)
    #comd = 'ln -s '+restart+' cice.r.nc'
    #print(comd)
    #os.system(comd)
    comd = './create_obs_sequence < input.txt > output.create_obs_sequence'
    print(comd)
    os.system(comd)
    txt = open('output.create_obs_sequence').readlines()
    if ' Finished ... at YYYY MM DD HH MM SS = \n' not in txt:
      print('create obs sequence did not finish correctly')
      sys.exit()
    else:
      os.system('rm output.create_obs_sequence')
    comd = 'rm obs_seq.in'
    os.system(comd)
    comd = 'ln -s '+save_date+' obs_seq.in'
    os.system(comd)
    comd = 'rm input_file.nc'
    os.system(comd)
    comd = 'ln -s '+path+' input_file.nc'
    os.system(comd)
    comd = './perfect_model_obs > output.perfect_model_obs_'+str(check_date)+'_'+label
    print(comd)
    os.system(comd)
    txt = open('output.perfect_model_obs_'+str(check_date)+'_'+label).readlines()
    if ' Finished ... at YYYY MM DD HH MM SS = \n' not in txt:
      print('perfect mode obs did not finish correctly')
      sys.exit()
    else:
      print('delete file')
      os.system('rm output.perfect_model_obs_'+str(check_date)+'_'+label)
    comd = 'mv obs_seq.out obs_seq.'+str(check_date)+'_'+label
    os.system(comd)
    #sys.exit()
  comd = 'ls obs_seq.'+str(check_date)+'_* > obs_seq_flist'
  os.system(comd)
  comd = './obs_sequence_tool > output.obs_sequence_tool'
  print(comd)
  os.system(comd)
  txt = open('output.obs_sequence_tool').readlines()
  if ' Finished ... at YYYY MM DD HH MM SS = \n' not in txt:
    print('obs sequence tool did not finish correctly')
    sys.exit()
  else:
    print('delete file')
    os.system('rm output.obs_sequence_tool')
  comd = 'mv obs_seq.out obs_seq.'+str(curr_date)[0:4]+'-'+str(curr_date)[4:6]+'-'+str(curr_date)[6:]+'-00000'
  os.system(comd)
  if os.path.exists('obs_seq.'+str(curr_date)[0:4]+'-'+str(curr_date)[4:6]+'-'+str(curr_date)[6:]+'-00000'):
    comd = 'mv obs_seq.in_'+str(check_date)+'_* obs_seq.'+str(check_date)+'_* raw_files/'
    os.system(comd)
    comd = 'mv obs_seq.'+str(curr_date)[0:4]+'-'+str(curr_date)[4:6]+'-'+str(curr_date)[6:]+'-00000 ./OBS/multivariate/mem_00'+str(ens_truth)+'/stdx1pt5/'
    os.system(comd)
    os.system('rm obs_seq_flist')
  else:
    print('Obs file was not created correctly....STOP')
    sys.exit() 
  # debug = True
  # if debug:
  #   import matplotlib.pyplot as plt 
  #   from mpl_toolkits.basemap import Basemap
  #   fig,ax = plt.subplots(1,1,figsize=(8,8))#,constrained_layout=True)
  #   m = Basemap(projection='npstere',boundinglat=55.,lon_0=270.,resolution='l')
  #   m.drawstates(color='#444444',linewidth=1.00)
  #   m.drawcoastlines(color='#444444',linewidth=1.0)
  #   m.drawcountries(color='#444444',linewidth=1.00)
  #   m.fillcontinents(color='grey')
  #   x,y = m(seaice_lon,seaice_lat)
  #   m.contourf(x,y,seaice_conc,cmap=plt.get_cmap('gray'))
  #   xpt,ypt = m(list(chain.from_iterable(total_lon)),list(chain.from_iterable(total_lat)))
  #   m.plot(xpt,ypt,'ok')
  #   plt.show()


