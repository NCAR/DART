#!/usr/bin/env python
#
# This code is not protected by the DART copyright agreement.
# DART $Id$

__author__    = "Rahul Mahajan"
__email__     = "rahul.mahajan@nasa.gov"
__copyright__ = "Copyright 2012, NASA / GSFC / GMAO"

#############################################################
# Generate an observation network and sequence for DART
#############################################################

import os, sys
import numpy                 as     np
from   matplotlib            import pyplot
from   mpl_toolkits.basemap  import Basemap

# User options - start
kmax                 = 128/2       # should match with spectral_mod.f90
lmax                 = 64/2        # should match with spectral_mod.f90
channel_center       = 45.0        # should match with spectral_mod.f90
channel_width        = 40.0        # should match with spectral_mod.f90
netsize              = 10          # number of observations on each level
errvar               = 0.0625      # error variance for each observation
nassim               = 400         # number of observation times in sequence
init_time_days       = 145731      # initial time (day) in sequence [01/01/2000 00:00:00]
init_time_seconds    = 0           # initial time (seconds) in sequence
period_time_days     = 0           # period of observations in sequence (days)
period_time_seconds  = 21600       # period of observations in sequence (seconds)
verbose              = True        # write information to stdout
clean                = True        # clean-up
# User options - end

# insure the same sequence of random numbers EVERY TIME
np.random.seed(0)

# create grid
if (verbose): print 'creating grid ...'
lons = np.zeros(2*kmax)
lats = np.zeros(2*lmax)
levs = np.zeros(2,dtype=int)
for i in range(0,2*kmax): lons[i] = 360.0*i/(2*kmax)
for i in range(0,2*lmax): lats[i] = channel_center + channel_width * (float(i)/(2*lmax) - 0.5)
for i in range(0,2):      levs[i] = i+1

# get observation index, convert to corresponding state, x, y index
if (verbose): print 'creating observation network ...'
state_ind = np.zeros(2 * netsize, dtype=int)
lon_ind   = np.zeros(2 * netsize, dtype=int)
lat_ind   = np.zeros(2 * netsize, dtype=int)
lev_ind   = np.zeros(2 * netsize, dtype=int)
for lev in levs:
    obs_ind = np.sort(np.random.random_integers(0,high=2*kmax*2*lmax,size=netsize))
    index = np.arange(0,netsize,1) + netsize*(lev-1)
    state_ind[index] = obs_ind + (2*kmax*2*lmax) * (lev-1)
    [lon_ind[index], lat_ind[index]] = np.unravel_index(obs_ind,[2*kmax,2*lmax],order='F')
    lev_ind[index] = lev-1

# draw a map of the observation network
if (verbose): print 'plotting observation network ...'
fig = pyplot.figure()
pyplot.title('Surface / Tropopause Observation Network')
map = Basemap(projection='cyl', \
              llcrnrlon=0,urcrnrlon=360,llcrnrlat=-85,urcrnrlat=85, \
              resolution='c')
map.drawcoastlines(color='0.8')
map.fillcontinents(color='0.8')
map.drawmeridians(np.arange(0, 360+45, 45),labels=[0,0,0,1],color='1.0',dashes=[0,1])
map.drawparallels(np.arange(-90,90,30),    labels=[0,1,0,0],color='1.0',dashes=[0,1])

# plot the first location for legend purposes
x,y = map(lons[lon_ind[0]],lats[lat_ind[0]])
pyplot.plot(x,y,'ro',markersize=5,markeredgecolor='r',label='Surface (%d)'    % netsize)
x,y = map(lons[lon_ind[-1]],lats[lat_ind[-1]])
pyplot.plot(x,y,'co',markersize=5,markeredgecolor='c',label='Tropopause (%d)' % netsize)
pyplot.legend(loc=0)

# first plot grid
lon, lat = np.meshgrid(lons, lats)
x,y = map(lon,lat)
pyplot.plot(x,y,'k.',markersize=0.5)

# second surface observation locations
x,y = map(lons[lon_ind[:netsize]],lats[lat_ind[:netsize]])
pyplot.plot(x,y,'ro',markersize=5,markeredgecolor='r')

# third tropopause observation locations
x,y = map(lons[lon_ind[netsize:]],lats[lat_ind[netsize:]])
pyplot.plot(x,y,'co',markersize=5,markeredgecolor='c')

pyplot.savefig('obs_network.png',dpi=100,orientation='portrait', format='png')
pyplot.savefig('obs_network.eps',dpi=300,orientation='landscape',format='eps')

os.system('display obs_network.png &')

# save the locations of the observations to file for later use
if (verbose): print 'saving observation network to disk ...'
fh = open('obs_network.txt','w')
fh.write('%d\n' % (2 * netsize))
for i in range(0,netsize): fh.write('%f %f %f\n' % (lats[lat_ind[i]], lons[lon_ind[i]], levs[lev_ind[i]]) )
for i in range(0,netsize): fh.write('%f %f %f\n' % (lats[lat_ind[i+netsize]], lons[lon_ind[i+netsize]], levs[lev_ind[i+netsize]]) )
fh.close()

# create a input file for create_obs_sequence in DART
if (verbose): print 'preparing for create_obs_sequence ...'
fh = open('create_obs_sequence.INPUT','w')

# Input upper bound on number of observations in sequence
fh.write('%d\n' % (len(state_ind)))

# Input number of copies of data (0 for just a definition)
fh.write('%d\n' % (0))

# Input number of quality control values per field (0 or greater)
fh.write('%d\n' % (0))

for i in range(0,len(state_ind)):

    # input a -1 if there are no more obs
    fh.write('%d\n' % (i+100000001))

    # Input -1 * state variable index for identity observations
    # OR input the name of the observation kind from table below:
    # OR input the integer index, BUT see documentation...
    #fh.write('%d\n' % (-1 * state_ind[i]))
    fh.write('%s\n' % ('POTENTIAL_TEMPERATURE'))
    #fh.write('%d\n' % (1))

    # Vertical coordinate options
    #      -2  --> vertical coordinate undefined
    #      -1  --> surface
    #       1  --> model level
    #       2  --> pressure
    #       3  --> height
    #       4  --> scale height
    fh.write('%d\n' % (1))
    # Vertical coordinate model level
    fh.write('%f\n' % levs[lev_ind[i]])
    # Input longitude: value 0 to 360.0 or a negative number for
    # Uniformly distributed random location in the horizontal
    fh.write('%f\n' % lons[lon_ind[i]])
    # Input latitude: value -90.0 to 90.0
    fh.write('%f\n' % lats[lat_ind[i]])

    # input time in days and seconds (as integers)
    fh.write('%d %d\n' % (0, 0))

    # Input error variance for this observation definition
    fh.write('%f\n' % (errvar))

#Input filename for sequence (  set_def.out   usually works well)
fh.write('%s\n' % ('set_def.out'))

fh.close()

if (verbose): print 'launching create_obs_sequence ...'
cmd = './create_obs_sequence < create_obs_sequence.INPUT'
os.system(cmd)

# create a input file for create_fixed_network_sequence in DART
if (verbose): print 'preparing for create_fixed_network_seq ...'
fh = open('create_fixed_network_seq.INPUT','w')

#Input filename for sequence (  set_def.out   usually works well)
fh.write('%s\n' % ('set_def.out'))

# To input a regularly repeating time sequence enter 1
# To enter an irregular list of times enter 2
fh.write('%d\n' % (1))

# Input number of observation times in sequence
fh.write('%d\n' % (nassim))

# Input initial time in sequence
# input time in days and seconds (as integers)
fh.write('%d %d\n' % (init_time_days, init_time_seconds))

# Input period of obs in sequence in days and seconds
fh.write('%d %d\n' % (period_time_days, period_time_seconds))

# What is output file name for sequence (  obs_seq.in   is recommended )
fh.write('%s\n' % ('obs_seq.in'))

fh.close()

if (verbose): print 'launching create_fixed_network_seq ...'
cmd = './create_fixed_network_seq < create_fixed_network_seq.INPUT'
os.system(cmd)

# cleaning up
if ( clean ):
    print 'cleaning up ...'
    os.system('rm -f set_def.out create_obs_sequence.INPUT create_fixed_network_seq.INPUT')

if (verbose): print 'all done ...'
sys.exit(0)

# <next few lines under version control, do not edit>
# $URL$
# $Revision$
# $Date$
