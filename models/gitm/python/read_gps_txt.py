import numpy as np #for data structures
import datetime, time #for keeping track of GITM files

import subprocess


def read(fn):
    f = open(fn, 'r'); 
    n_lines = int(subprocess.check_output('wc -l < ' + fn, shell=True))-1 #number of nonheader lines in the file, -1 because don't need to preallocate for header
    timeD = np.empty(n_lines, datetime.datetime ) #preallocation will be slightly too big, as I throw away the missing values
    LonD = np.zeros( n_lines, float ) #D stands for data
    LatD = np.zeros( n_lines, float )
    VtecD = np.zeros( n_lines, float ) 
    VtecsdD = np.zeros( n_lines, float ) #sd stands for standard deviation (uncertainty, error in the measurement)

    line = f.readline().split() #read the header
    lat_ind = line.index('GDLAT') #index of the latitude column - in some files this will be 6, in some 9
    i = 0 #counter
    print 'Reading GPS text file... '
    for line in f:
        if 'missing' not in line: #remove any lines with missing values
            line = line.split()
            d = np.array(map(int, line[0:6])) #date part, 6 is noninclusive
            timeD[i] = datetime.datetime(d[0], d[1], d[2], d[3], d[4], d[5])
            LatD[i], LonD[i], VtecD[i], VtecsdD[i] = np.array(map(float, line[lat_ind:])) #data part
            i += 1
    f.close()
    print 'Done w/ GPS text file.'

#remove empty entries (empty entries appear at the end of the arrays and are due to throwing away missing values)
    timeD = timeD[0:i]
    LonD = LonD[0:i]
    LatD = LatD[0:i]
    VtecD = VtecD[0:i]
    VtecsdD = VtecsdD[0:i]
    
    LonD[LonD < 0] = LonD[LonD < 0] + 360 #make 0 <= lon <=360 
    
    return timeD, LonD, LatD, VtecD, VtecsdD

# DART $Id$
# from Alexey Morozov

# <next few lines under version control, do not edit>
# $URL$
# $Revision$
# $Date$
