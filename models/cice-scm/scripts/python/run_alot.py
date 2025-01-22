import pandas as pd
from datetime import datetime,timedelta
import xarray as xr
import numpy as np
import glob
from itertools import chain
import os
from datetime import timedelta
from os import path
from scipy.interpolate import CubicSpline

def remove_leapyears(datelist):
  new_array = []
  for d in datelist:
    if int(d.strftime('%m')) == 2 and int(d.strftime('%d')) == 29:
      continue
    else:
      new_array.append(d)
  return np.array(new_array)

#datelist = pd.date_range(start=datetime(2011,1,2),end=datetime(2011,1,31)).to_pydatetime()
#datelist = pd.date_range(start=datetime(2011,4,1),end=datetime(2011,5,31)).to_pydatetime()
datelist = pd.date_range(start=datetime(2011,1,2),end=datetime(2011,3,31)).to_pydatetime()
datelist = pd.date_range(start=datetime(2011,1,2),end=datetime(2011,1,12)).to_pydatetime()
datelist = remove_leapyears(datelist)

for t in range(0,datelist.shape[0]):
#for t in range(0,365):
  count = 0
  print('--------------------------------')
  print(datelist[t])
  curr_date = (datelist[t].strftime('%Y%m%d'))
  
  comd = './CICE-SCM-DART.csh '+curr_date+'00'
  print(comd)
  os.system(comd)



