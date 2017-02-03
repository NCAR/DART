#plot difference along with f107
al = 0.3 #transparency
#vi,va = 0,max(VtecED[:,:,-1].flatten()) #specify zlim and colorlimits
vi, va = 0., 100.
levels = np.arange(vi,va)
f = open('map.txt', 'r'); elev = np.array(map(int, f.read().split())).reshape((360, 180)).T; f.close() #load the Earth map
for ti in range(len(tcommon)): #sorted(tcommon): #timeE[1:3]:#
    fig = plt.figure(1)
    fig.clf() #in case it was already open
    ax = plt.subplot2grid((3,1), (0,0), rowspan=2)
#map
    N,T= np.meshgrid(np.linspace(0,359,360),np.linspace(-90,89,180))
    ma = ax.contour(N, T, elev, 0, colors='black') #the Earth map plotted on the bottom of the z-axis
#diff
    XE,YE = np.meshgrid(loni, lati) 
    contED = ax.contourf(XE, YE, VtecEDi[:,:,ti].T, levels, alpha=al, cmap=plt.cm.jet, extend='both') #cmap=plt
#    XE,YE = np.meshgrid(LonE, LatE) 
#    contED = ax.contourf(XE, YE, VtecED[:,:,ti].T, 50, alpha=al, vmin=vi, vmax=va, cmap=plt.cm.hot) #cmap=plt.cm.hot
#    cbar=plt.colorbar(contED)    
#    cbar.ax.set_ylabel('$|\Delta VTEC|$')
    #cbar.orientation='horizontal'
    
#labels, etc
    #ax.set_title(str(ti))
    ax.grid(True)
    ax.set_xticks(np.linspace(0,360,7) );
    ax.set_yticks(np.linspace(-90,90,7) );
#    ax.set_xlabel('Longitude [deg]');
    ax.set_ylabel('Latitude [deg]');
    ax.text(0, .01, str(tcommon[ti]), fontsize=10,  transform=ax.transAxes)
    
#f107 plot
    ax = plt.subplot2grid((3,1), (2, 0))
    fT, = ax.plot([0,24],[148,148], '--k', label='$F_{10.7}$ as measured by NOAA')
    fi, = ax.plot([0,24],[130,130], 'g', label='Initial $F_{10.7}$ estimate')
    tE = np.array([(t.day-1)*24 + t.hour + float(t.minute)/60.0 for t in timeE])
    fEo, = ax.plot(tE,f107Eo, label='EAKF ensemble mean')
    x = np.float32(np.concatenate(( tE, tE[::-1]))) 
    y = np.concatenate(( f107Eo+f107sdEo, (f107Eo-f107sdEo)[::-1]))
    p = np.zeros((len(x),2), np.float32)
    p[:,0] = x
    p[:,1] = y
    po = Polygon(p , True, alpha=0.3, color='b', label='EAKF ensemble mean +-SD')
    fsdEo = ax.add_patch(po)
#labels, etc
    plt.rcParams['legend.fontsize'] = 10
    h, l = ax.get_legend_handles_labels()
    ax.legend(h, l, loc=4)
    fdot = ax.plot([(t.day-1)*24 + t.hour + float(t.minute)/60.0 for t in timeE[timeE==tcommon[ti]]],f107Eo[timeE==tcommon[ti]],'o',color='red')
#labels, etc
#ax.set_title(str(ti))
    ax.set_xticks(np.linspace(0,24,24/3+1) );
#ax.set_yticks(np.linspace(-90,90,7) );
    ax.set_xlim(0,24)
    ax.set_ylim(100,300)
    ax.set_xlabel('Time [hr]');
    ax.set_ylabel('$F_{10.7}$ [SFU]');
    ax.grid(True)
    #    fig.canvas.draw()    
    plt.savefig('vtec' + str(ti).zfill(2) + '_diff_f.png', bbox_inches=0)
    if ti==0:   plt.savefig('vtec00_diff_f.eps', bbox_inches=0)
    print ti

print subprocess.check_output('convert *_diff_f.png vtec00_diff_f.gif', shell=True) #convert to animated gif 

# DART $Id$
# from Alexey Morozov

# <next few lines under version control, do not edit>
# $URL$
# $Revision$
# $Date$
