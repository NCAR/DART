#interpolate truth to write "higher res" data above US
unc = 0.05 #standard deviation uncertainty in TEC in percent (0.05 implies 5% uncertainty)
res = 2.0 #lat,lon resolution to which interpolate the data
loni = np.arange(240.,281.,res) #what lons do you want to interpolate it to
lati = np.arange(30.,51.,res) #what lons do you want to interpolate it to
LonTsi = np.zeros( (len(loni)*len(lati)*len(timeT) ), float) #longitude truth selectively interpolated
LatTsi = np.zeros( (len(loni)*len(lati)*len(timeT) ), float)
timeTsi = np.zeros( (len(loni)*len(lati)*len(timeT) ), object)
VtecTsi = np.zeros( (len(loni)*len(lati)*len(timeT) ), float)
string = 'YEAR MONTH DAY HOUR MINUTE SECOND IGNORE_I IGNORE_I IGNORE_I GDLAT GDLON TEC DTEC\n' #header
for ti in range(len(timeT)): 
    f = interpolate.RectBivariateSpline(LonT,LatT,VtecT[:,:,ti])
    for t in range(len(lati)):
        for n in range(len(loni)):
            #print n + t*len(loni) + ti*len(loni)*len(lati)
            LonTsi[n + t*len(loni) + ti*len(loni)*len(lati)] = loni[n]
            LatTsi[n + t*len(loni) + ti*len(loni)*len(lati)] = lati[t]
            timeTsi[n + t*len(loni) + ti*len(loni)*len(lati)] = timeT[ti]
            VtecTsi[n + t*len(loni) + ti*len(loni)*len(lati)] = f( loni[n], lati[t] ) 
            l = loni[n]
            if l>180.: l=l-360.
            string = string + str(timeT[ti].year) + ' ' + str(timeT[ti].month).zfill(2) + ' ' + str(timeT[ti].day).zfill(2) + ' ' + str(timeT[ti].hour).zfill(2) + ' ' + str(timeT[ti].minute).zfill(2) + ' ' + str(timeT[ti].second).zfill(2) + ' 1 1 1 ' + ('%7.1f' % lati[t]) + ' ' +  ('%7.1f' % l) + ' ' + ('%7.2f' % f( loni[n], lati[t] )[0][0]) + ' ' + ('%7.2f' % (f( loni[n], lati[t] )*unc)[0][0]) + '\n'

#write the first part of the string (US data), as it gets too long if keep waiting
fn = 'gps021201_GITM.txt'
file_gps = open(fn, 'w')
file_gps.write(string)


#add the data to the string from EUROPE

loni = np.arange(0.,41.,res) #what lons do you want to interpolate it to
#lati = np.arange(-10.,11.,res) #Africa
lati = np.arange(40.,61.,res) #what lons do you want to interpolate it to
LonTsi = np.zeros( (len(loni)*len(lati)*len(timeT) ), float) #longitude truth selectively interpolated
LatTsi = np.zeros( (len(loni)*len(lati)*len(timeT) ), float)
timeTsi = np.zeros( (len(loni)*len(lati)*len(timeT) ), object)
VtecTsi = np.zeros( (len(loni)*len(lati)*len(timeT) ), float)
string = ''
for ti in range(len(timeT)): 
    f = interpolate.RectBivariateSpline(LonT,LatT,VtecT[:,:,ti])
    for t in range(len(lati)):
        for n in range(len(loni)):
            #print n + t*len(loni) + ti*len(loni)*len(lati)
            LonTsi[n + t*len(loni) + ti*len(loni)*len(lati)] = loni[n]
            LatTsi[n + t*len(loni) + ti*len(loni)*len(lati)] = lati[t]
            timeTsi[n + t*len(loni) + ti*len(loni)*len(lati)] = timeT[ti]
            VtecTsi[n + t*len(loni) + ti*len(loni)*len(lati)] = f( loni[n], lati[t] ) 
            l = loni[n]
            if l>180.: l=l-360.
            string = string + str(timeT[ti].year) + ' ' + str(timeT[ti].month).zfill(2) + ' ' + str(timeT[ti].day).zfill(2) + ' ' + str(timeT[ti].hour).zfill(2) + ' ' + str(timeT[ti].minute).zfill(2) + ' ' + str(timeT[ti].second).zfill(2) + ' 1 1 1 ' + ('%7.1f' % lati[t]) + ' ' +  ('%7.1f' % l) + ' ' + ('%7.2f' % f( loni[n], lati[t] )[0][0]) + ' ' + ('%7.2f' % (f( loni[n], lati[t] )*unc)[0][0]) + '\n'


#write the Europe data
file_gps.write(string)
file_gps.close()

di = '/home/morozova/HEAD/DART13/observations/gnd_gps_vtec/work/' #where to copy the file
print subprocess.check_output('cp ' + fn + ' ' + di , shell=True) #copy the file and print status if needed
command = 'sed -i\'.tmp\' \'s/text_input_file.*=.*/text_input_file = "' + fn + '",/\' ' + di + 'input.nml'
print command
print subprocess.check_output(command , shell=True) #change the name in input.nml
command = 'cd ' + di + '; ./text_to_obs'
print command
print subprocess.check_output(command, shell=True) #convert .txt into dart obs file
print subprocess.check_output('pwd', shell=True) 
fn2 = 'obs_seq_01to02_tec_9x9_useu.out'
command = 'cp ' + di + 'obs_seq.out /home/morozova/HEAD/DART13/models/gitm/work/' + fn2
print command
print subprocess.check_output(command, shell=True) #convert .txt into dart obs file

# DART $Id$
# from Alexey Morozov

# <next few lines under version control, do not edit>
# $URL$
# $Revision$
# $Date$
