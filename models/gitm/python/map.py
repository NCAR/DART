#plotting the map without "basemap". If you want just the shores - uncomment line 25 and comment out 27-31
import numpy as np
import matplotlib.pyplot as plt

f = open('map.txt', 'r'); elev = np.array(map(int, f.read().split())).reshape((360, 180)).T; f.close()
C = plt.contour(np.linspace(0,359,360), np.linspace(-90,89,180), elev, 0, colors='black') #just a contour at 0 elevation

# DART $Id$
# from Alexey Morozov
#
# <next few lines under version control, do not edit>
# $URL$
# $Revision$
# $Date$
