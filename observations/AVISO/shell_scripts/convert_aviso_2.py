#!/usr/bin/env python
#
# DART software - Copyright UCAR. This open source software is provided
# by UCAR, "as is", without charge, subject to all terms of use at
# http://www.image.ucar.edu/DAReS/DART/DART_download
#
# DART $Id$
#
# CREDIT: This script was donated to DART by Fred Castruccio Thanks Fred!
# Romain Escudier of Rutgers added some logic to avoid recreating files
# that already exist, and the ability to convert a particular year. 

import numpy as np
import os
import sys
import subprocess
import shutil
import glob

year = '2011'

obs_in_dir  = '/glade/p/cesm/omwg_dev/AVISO/SLA/along_track'

# Pick an output directory ...
obs_out_dir = '/glade/p/work/${USER}/Observations/AVISO/data'
obs_out_dir = '/glade/p/work/${USER}/Data/ObsData/AVISO/ATerr003/SATS/'+year
obs_out_dir = '/glade/scratch/${USER}/Observations/data'

try:
  os.stat(obs_out_dir)
except:
  os.mkdir(obs_out_dir) 

# convert_aviso can take two command-line arguments.
# The second argument is a file containing the list of observation error standard deviations
# for each observation. If this file is not present, the namelist value is used. 
# Fred is the only one to use the second argument - I am not sure how he created the file.

# If the output file already exists, it is not created again.
# If the output file creation failed, the script continues ...

flst = glob.glob(obs_in_dir + '/*/'+year+'/*')

for fobs in flst:

  try:
    platform = fobs.split('_')[-5]
    year  = np.int(fobs.split('_')[-2][:4])
    month = np.int(fobs.split('_')[-2][4:6])
    day   = np.int(fobs.split('_')[-2][6:])
    obs_seq_out = '{}/obs_seq.{}.{:04d}{:02d}{:02d}'.format(obs_out_dir,platform,year,month,day)
    if not os.path.isfile(obs_seq_out):
       cmd = ['../work/convert_aviso', fobs]
       print 'convert_aviso {}'.format(os.path.basename(fobs))
       pipe = subprocess.Popen(cmd)
       stdout, stderr = pipe.communicate()
       if os.path.isfile('obs_seq.aviso'):
          shutil.move('obs_seq.aviso',obs_seq_out)
  except OSError as err:
    print('WARNING',err.errno,err.strerror)

exit(0)    

# <next few lines under version control, do not edit>
# $URL$
# $Revision$
# $Date$

