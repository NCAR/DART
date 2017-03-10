#!/usr/bin/env python
#
# DART software - Copyright UCAR. This open source software is provided
# by UCAR, "as is", without charge, subject to all terms of use at
# http://www.image.ucar.edu/DAReS/DART/DART_download
#
# DART $Id$
#
# CREDIT: This script was donated to DART by Fred Castruccio Thanks Fred!

import numpy as np
import os
import sys
import subprocess
import shutil


obs_in_dir  = '/glade/p/work/fredc/OBS/SLA_alongtrack'
obs_out_dir = '/glade/p/work/fredc/Observations/AVISO/data'


ylst = [2008]
mlst = [1]
if (len(mlst) == 0):
  mlst = [1,2,3,4,5,6,7,8,9,10,11,12]

plst = ['j1','en','g2']

try:
  os.stat(obs_out_dir)
except:
  os.mkdir(obs_out_dir) 

for year in ylst:
  daysinmonth = [31,28,31,30,31,30,31,31,30,31,30,31]
  if (year%4 == 0):
    daysinmonth = [31,29,31,30,31,30,31,31,30,31,30,31]
  for month in mlst:
    for day in range(1,daysinmonth[month-1]+1):
      for platform in plst:
        try:
          fobs = '{}/dt_global_{}_sla_vfec_{:04d}{:02d}{:02d}_20140106.nc'.format(obs_in_dir,platform,year,month,day)
          obs_seq_out = '{}/obs_seq.{}.{:04d}{:02d}{:02d}'.format(obs_out_dir,platform,year,month,day)
          cmd = ['convert_aviso', fobs]
          print 'convert_aviso {}'.format(os.path.basename(fobs))
          pipe = subprocess.Popen(cmd)
          stdout, stderr = pipe.communicate()
          shutil.move('obs_seq.aviso',obs_seq_out)
        except OSError as err:
          print('WARNING',err.errno,err.strerror)

exit(0)    

# <next few lines under version control, do not edit>
# $URL$
# $Revision$
# $Date$

