#plot difference
al = 0.5 #transparency
#vi,va = 0,max(VtecED[:,:,-1].flatten()) #specify zlim and colorlimits
vi, va = 0., 20.
#levels = np.arange(vi,va)
f = open('map.txt', 'r'); elev = np.array(map(int, f.read().split())).reshape((360, 180)).T; f.close() #load the Earth map
for ti in range(len(tcommon)): #sorted(tcommon): #timeE[1:3]:#
    fig = plt.figure(1)
    fig.set_size_inches((8,4))
    fig.clf() #in case it was already open
    ax = fig.add_subplot(111)
#map
    N,T= np.meshgrid(np.linspace(0,359,360),np.linspace(-90,89,180))
    ma = ax.contour(N, T, elev, 0, colors='black') #the Earth map plotted on the bottom of the z-axis
#data
    x = LonD[timeD1==tcommon[ti]]
    y = LatD[timeD1==tcommon[ti]]
    z = VtecD[timeD1==tcommon[ti]]
    scatD = ax.scatter(x, y, c=z, cmap=plt.cm.jet, s=25, vmin=vi, vmax=va, alpha=al, lw=0) 
    cbar = plt.colorbar(scatD)
    cbar.ax.set_ylabel('TEC  [TECU]')
    #cbar.orientation='horizontal'
    
#labels, etc
    #ax.set_title(str(ti))
    plt.grid(True)
    ax.set_xlim(0,360)
    ax.set_ylim(-90,90)
    ax.set_xticks(np.linspace(0,360,7) );
    ax.set_yticks(np.linspace(-90,90,7) );
    ax.set_xlabel('Longitude [deg]');
    ax.set_ylabel('Latitude [deg]');
    ax.text(0, .01, str(tcommon[ti]), fontsize=10,  transform=ax.transAxes)
    #    fig.canvas.draw()    
    plt.savefig('data' + str(ti).zfill(2) + '.png', bbox_inches=0)
    if ti==0:   plt.savefig('data00.eps', bbox_inches=0)
    print ti
    #plt.ion()

print subprocess.check_output('convert data*.png data00.gif', shell=True) #convert to animated gif 

# DART $Id$
# from Alexey Morozov

# <next few lines under version control, do not edit>
# $URL$
# $Revision$
# $Date$
