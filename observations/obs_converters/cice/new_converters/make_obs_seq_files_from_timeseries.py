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

# define function to fine the lat, lon coordinates of points on a sphere
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


# set a list of observation dates to convert 
datelist = pd.date_range(start=datetime(2019,4,1),end=datetime(2020,1,1)).to_pydatetime()

# set the name of the experiment to pull from
exp = 'free_dynamic_repeat'

# set the ensemble truth member tag value
ens_truth = 24
print('random truth is: ' + str(ens_truth))

# open the track file (observation locations) and extract date information
print('Opening tracks file..')
track_data = xr.open_dataset('simulated_tracks_2019.nc')
track_yearmonday = track_data.year_month_day.values
print(track_yearmonday)
track_hours = np.array([int(str(d)[8:-2]) for d in track_data.dates.values])
print(track_hours)

# open the sea ice timeseries data file
path = glob.glob('/glade/scratch/mollyw/archive/'+exp+'/ice/proc/tseries/day_1/'+exp+'.cice_00'+str(ens_truth)+'.h1.aice.*.nc')
print(path[0])
sic_data = xr.open_dataset(path[0])

# begin cycling through the prescribed dates
for t in range(0,datelist.shape[0]):
  count = 0
  print('--------------------------------')
  total_lat = []
  total_lon = []

  # print the date for sanity
  print(datelist[t])

  # find the current date (the date of the assimilation, associated with the restart),
  # and the previous date (associated with the history files that provide SIC errors)
  rest_date = int(datelist[t].strftime('%Y%m%d'))
  hist_date = int((datelist[t]-timedelta(days=1)).strftime('%Y%m%d'))

  # evaluate the model times to be a comparable format
  model_times = [int(sic_data.time.values[i].strftime('%Y%m%d')) for i in range(0,len(sic_data.time))]
  model_times = np.array(model_times) 

  # isolate the year, mon, and day of the restart date
  year = (datelist[t]).strftime('%Y')
  mon = (datelist[t]).strftime('%m')
  day = (datelist[t]).strftime('%d')

  # isolate the year, mon, and day of the history date
  year_old = (datelist[t]-timedelta(days=1)).strftime('%Y')
  mon_old = (datelist[t]-timedelta(days=1)).strftime('%m')
  day_old = (datelist[t]-timedelta(days=1)).strftime('%d')

  # deal with the calendar index by determining the no. of days since 1601-01-01
  diff = datelist[t] - datetime(1601,1,1,0)

  # replace index where necessary in the input.nml file 
  comd = 'cat > test.sed << EOF\ns#{days}#'+str(diff.days)+'#g\ns#{seconds}#'+str(diff.seconds)+'#g\nEOF'
  os.system(comd)
  comd = 'sed -f test.sed input.nml.template > input.nml'
  os.system(comd)
  os.system('rm test.sed')

  # get indices for where the track data that matches the desired date
  ind_track = np.where(track_yearmonday==rest_date)[0]
  print(ind_track)

  # get indices for where the model sic data matches the desired date
  time_ind = np.where(model_times==hist_date)[0]
  print(time_ind)
  
  # raise an exception if there are no places where the track data matches the desired date 
  if ind_track.shape[0] == 0 or time_ind.shape[0] == 0:
    print('Dates could not be found in one of the two files...')
    sys.exit()

  # grab the lat, lon, and hour of the tracks in question
  track_lat = track_data.latitude.values[ind_track[0]]
  track_lon = track_data.longitude.values[ind_track[0]]
  track_hrs = track_hours[ind_track]

  # write out the no. of points available (filtered for NH as well)
  print('Number of sat. points for that date:'+str(track_lat.shape[0]))
  ind = np.where(track_lat>30)[0]
  track_lat = track_lat[ind]
  track_lon = track_lon[ind]
  print('Number of sat. points after filtering for NH:'+str(track_lat.shape[0]))
  unique_hrs = np.unique(track_hrs)
  # This section of code cycles through the different hours in the observations. We don't have
  # this in our files, the hour written is 0000

  # grab the appropriate sea ice model data
  seaice_lat = sic_data.TLAT.values
  seaice_lon = sic_data.TLON.values 
  seaice_conc = sic_data.aice.values[time_ind,:,:].squeeze()
  mask = sic_data.tmask.values[:,:]
  sic_data.close()

  # for every unique hour...  
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

    save_date = 'obs_seq.in_'+str(rest_date)+'_'+label
    restart = '/glade/scratch/mollyw/'+exp+'/run/'+exp+'.cice_0024.r.{0}-{1}-{2}-{3}.nc'.format(year,mon,day,label)
    if len(glob.glob(path[0])) == 0:
      print('There is no seaice data for '+str(rest_date)+' at hour '+str(label))
      continue
    
    for pp in range(0,hold_trk_lat.shape[0]):
      a = abs(seaice_lat-hold_trk_lat[pp])+abs(seaice_lon-hold_trk_lon[pp])
      i,j = np.unravel_index(np.nanargmin(a),a.shape)

      np.unravel_index
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
        # hold_thic.append(seaice_thick[i,j])
    total_lat.append(hold_lat)
    total_lon.append(hold_lon)
    hold_lat = np.array(hold_lat)
    hold_lon = np.array(hold_lon)
    if hold_lat.shape[0] == 0 or hold_lon.shape[0] == 0:
      print('No obs over Arctic sea ice so skip time....')
      continue
    hold_conc = np.array(hold_conc)
    # hold_thic = np.array(hold_thic)
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
       error = 0.05
      else:
       error = hold_conc[l]*0.15
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
      new_file.writelines('0.0025\n')
      
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
      new_file.writelines('0.01\n')
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
    comd = 'ln -s '+restart+' input_file.nc'
    os.system(comd)
    comd = './perfect_model_obs > output.perfect_model_obs_'+str(hist_date)+'_'+label
    print(comd)
    os.system(comd)
    txt = open('output.perfect_model_obs_'+str(hist_date)+'_'+label).readlines()
    if ' Finished ... at YYYY MM DD HH MM SS = \n' not in txt:
      print('perfect model obs did not finish correctly')
      sys.exit()
    else:
      print('delete file')
      os.system('rm output.perfect_model_obs_'+str(hist_date)+'_'+label)
    comd = 'mv obs_seq.out obs_seq.'+str(hist_date)+'_'+label
    os.system(comd)
    #sys.exit()
  comd = 'ls obs_seq.'+str(hist_date)+'_* > obs_seq_flist'
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
  comd = 'mv obs_seq.out obs_seq.'+str(rest_date)[0:4]+'-'+str(rest_date)[4:6]+'-'+str(rest_date)[6:]+'-00000'
  os.system(comd)
  if os.path.exists('obs_seq.'+str(rest_date)[0:4]+'-'+str(rest_date)[4:6]+'-'+str(rest_date)[6:]+'-00000'):
    comd = 'mv obs_seq.in_'+str(hist_date)+'_* obs_seq.'+str(hist_date)+'_* raw_files/'
    os.system(comd)
    comd = 'mv obs_seq.'+str(rest_date)[0:4]+'-'+str(rest_date)[4:6]+'-'+str(rest_date)[6:]+'-00000 ./OBS/validate/multivariate/mem_00'+str(ens_truth)+'/std1/'
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


