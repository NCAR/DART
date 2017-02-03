#plot difference
al = 0.3 #transparency
#vi,va = 0,max(VtecED[:,:,-1].flatten()) #specify zlim and colorlimits
if sim:
    vi, va = 0., 10.
    levels = np.linspace(vi,va,51)
else:
    vi, va = 0., 50.
    levels = np.linspace(vi,va,51)

f = open('map.txt', 'r'); elev = np.array(map(int, f.read().split())).reshape((360, 180)).T; f.close() #load the Earth map
for ti in range(len(tcommon)): #sorted(tcommon): #timeE[1:3]:#
    fig = plt.figure(1)
    fig.set_size_inches((8,4))
    fig.clf() #in case it was already open
    ax = fig.add_subplot(111)
#map
    N,T= np.meshgrid(np.linspace(0,359,360),np.linspace(-90,89,180))
    ma = ax.contour(N, T, elev, 0, colors='black') #the Earth map plotted on the bottom of the z-axis
#diff
    XE,YE = np.meshgrid(loni, lati) 
    contED = ax.contourf(XE, YE, VtecEDi[:,:,ti].T, levels, alpha=al, cmap=plt.cm.jet, extend='both') #cmap=plt
#    XE,YE = np.meshgrid(LonE, LatE) 
#    contED = ax.contourf(XE, YE, VtecED[:,:,ti].T, 50, alpha=al, vmin=vi, vmax=va, cmap=plt.cm.hot) #cmap=plt.cm.hot
    cbar=plt.colorbar(contED)    
    if sim:  
        cbar.ax.set_ylabel('$|V_{GITM} - V_{EAKF}|$')
    else:
        cbar.ax.set_ylabel('$|V_{GPS} - V_{EAKF}|$')
    #cbar.orientation='horizontal'
    
#labels, etc
    #ax.set_title(str(ti))
    plt.grid(True)
    ax.set_xticks(np.linspace(0,360,7) );
    ax.set_yticks(np.linspace(-90,90,7) );
    ax.set_xlabel('Longitude [deg]');
    ax.set_ylabel('Latitude [deg]');
    ax.text(0, .01, str(tcommon[ti]), fontsize=10,  transform=ax.transAxes)
    #    fig.canvas.draw()    
    plt.savefig('diff' + str(ti).zfill(2) + '.png', bbox_inches=0)
    plt.savefig('diff' + str(ti).zfill(2) + '.eps', bbox_inches=0)
    #if ti==0:   plt.savefig('diff0.eps', bbox_inches=0)
    print ti
    #plt.ion()

print subprocess.check_output('convert diff*.png diff00.gif', shell=True) #convert to animated gif 

# DART $Id$
# from Alexey Morozov

# <next few lines under version control, do not edit>
# $URL$
# $Revision$
# $Date$
