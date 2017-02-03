import struct #for converting binary fortran to numbers, to read GITM data
import numpy as np #for data structures

def read(time, prefix):
    """
    returns a tuple of lon, lat, alt, electron_dens
    
    INPUTS:
    - time[ti] - np.array of "datetime" times
    - prefix - path to the files as well as their type. example: '../work/advance_temp_e21/data/3DALL_t'
    
    OUTPUTS:
    - LonT[n] - np.array of lons
    - LatT[t] - np.array of lats
    - AltT[a] - np.array of alts
    - EdsT[n,t,a,ti] - np.array over [lon,lat,alt,time]
    """
    for t in time: 
        filename = prefix + str(t.year-2000).zfill(2) + str(t.month).zfill(2) + str(t.day).zfill(2) + '_'+ str(t.hour).zfill(2) + str(t.minute).zfill(2) + str(t.second).zfill(2) + '.bin'
        #print filename

        f = open(filename, mode='rb')
    
        nv=struct.unpack("i", f.read(4))[0]/8 #how many numbers I'm about to read
        v=struct.unpack("d"*nv, f.read(nv*8)) #read the data
        nv=struct.unpack("i", f.read(4))[0]/8 #how many numbers did I just read?
        
        nv=struct.unpack("i", f.read(4))[0]/4 #how many numbers I'm about to read
        nLonsTotal, nLatsTotal, nAltsTotal=struct.unpack("i"*nv, f.read(nv*4)) #read the data
        nv=struct.unpack("i", f.read(4))[0]/4 #how many numbers did I just read?
        
        nv=struct.unpack("i", f.read(4))[0]/4 #how many numbers I'm about to read
        nVars=struct.unpack("i"*nv, f.read(nv*4))[0] #read the data
        nv=struct.unpack("i", f.read(4))[0]/4 #how many numbers did I just read?
        
        Variables=[]
        for iVar in range(nVars):
            nv=struct.unpack("i", f.read(4))[0] #how many numbers I'm about to read
            Variables.append( f.read(nv) )
            if Variables[iVar].find('Longitude') > 0 : iVarn=iVar
            if Variables[iVar].find('Latitude') > 0 : iVart=iVar
            if Variables[iVar].find('Altitude') > 0 : iVara=iVar
            if Variables[iVar].find('[e-') > 0 : iVarEDS=iVar
            nv=struct.unpack("i", f.read(4))[0] #how many numbers did I just read?
        
        nv=struct.unpack("i", f.read(4))[0]/4 #how many numbers I'm about to read
        iYear, iMonth, iDay, iHour, iMinute, iSecond, iMilli = struct.unpack("i"*nv, f.read(nv*4)) #read the data
        nv=struct.unpack("i", f.read(4))[0]/4 #how many numbers did I just read?
        
        if t==time[0]: #the first time, preallocate
            AllData = np.zeros( (nLonsTotal*nLatsTotal*nAltsTotal, nVars), float)
            EdsT = np.zeros( (nLonsTotal-4, nLatsTotal-4, nAltsTotal-4, len(time)), float)
        
        for iVar in range(nVars):
            nv=struct.unpack("i", f.read(4))[0]/8 #how many numbers I'm about to read
            AllData[:,iVar] = np.array(struct.unpack("d"*nv, f.read(nv*8)))
            nv=struct.unpack("i", f.read(4))[0]/8 #how many numbers did I just read?
        
        f.close()
        
#strip off the ghost cells
        EdsT[:,:,:,np.where(time==t)[0][0]] = AllData[:,iVarEDS].reshape( (nLonsTotal, nLatsTotal, nAltsTotal), order='F')[2:nLonsTotal-2, 2:nLatsTotal-2, 2:nAltsTotal-2]
        

    LonT = AllData[:,iVarn].reshape( (nLonsTotal, nLatsTotal, nAltsTotal), order='F')[2:nLonsTotal-2, 0, 0]
    LatT = AllData[:,iVart].reshape( (nLonsTotal, nLatsTotal, nAltsTotal), order='F')[0, 2:nLatsTotal-2, 0]
    AltT = AllData[:,iVara].reshape( (nLonsTotal, nLatsTotal, nAltsTotal), order='F')[0, 0, 2:nAltsTotal-2]
    
    return (LonT, LatT, AltT, EdsT)

# DART $Id$
# from Alexey Morozov

# <next few lines under version control, do not edit>
# $URL$
# $Revision$
# $Date$
