#3dplot
al = 0.3 #transparency
vi,va = 0,max(VtecE[:,:,-1].flatten()) #specify zlim and colorlimits
f = open('map.txt', 'r'); elev = np.array(map(int, f.read().split())).reshape((360, 180)).T; f.close() #load the Earth map
tii=0 #counter
for ti in tcommon: #timeE[1:3]:#
    print tii
    fig = plt.figure(1)
    fig.clf() #in case it was already open
    ax = fig.add_subplot(111, projection='3d')
    
#map
    N,T= np.meshgrid(np.linspace(0,359,360),np.linspace(-90,89,180))
    ma = ax.contour(N, T, elev, 0, zs=0, zdir='z', colors='black') #the Earth map plotted on the bottom of the z-axis
    
#GITM truth
    XT,YT = np.meshgrid(LonT, LatT) 
    surfT = ax.plot_surface(XT, YT, VtecT[:,:,np.where(timeT==ti)[0][0]].T, rstride=1, cstride=1, alpha=al, color='blue', linewidth=0, vmin=vi, vmax=va) #cmap=plt.cm.hot
    #plt.colorbar(surfT)    
    
#DART ensemble
    XE,YE = np.meshgrid(LonE, LatE) 
    surfE = ax.plot_surface(XE, YE, VtecE[:,:,np.where(timeE==ti)[0][0]].T, rstride=1, cstride=1, alpha=al, color='red', linewidth=0, vmin=vi, vmax=va)
    
#GPS data
    XE,YE = np.meshgrid(LonE, LatE) 
    surfD = ax.plot_surface(XE, YE, VtecDi[:,:,tii].T, rstride=1, cstride=1, alpha=al, color='green', linewidth=0, vmin=vi, vmax=va)
    
#legend
    p1 = plt.Rectangle((0, 0), 1, 1, fc='b') #normal legend doesn't work for mplot3d
    p2 = plt.Rectangle((0, 0), 1, 1, fc='r')
    p3 = plt.Rectangle((0, 0), 1 ,1, fc='g')
    plt.rcParams['legend.fontsize'] = 10
    ax.legend([p1, p2, p3], ['no-DA GITM', 'GITM with EAKF', 'GPS VTEC data'], loc=1);
    
#labels, etc
    #ax.set_title(str(ti))
    ax.set_zlim(vi,va);
    ax.set_xticks(np.linspace(0,360,7) );
    ax.set_yticks(np.linspace(-90,90,7) );
    ax.set_xlabel('Longitude [deg]');
    ax.set_ylabel('Latitude [deg]');
    ax.set_zlabel('VTEC [TECU]');
    ax.text(0, 0, 0, str(ti), fontsize=10,  transform=ax.transAxes)
    
    #    fig.canvas.draw()    
    plt.savefig('interp' + str(tii).zfill(2) + '.png', bbox_inches=0)
    if tii==0:   plt.savefig('interp00.eps', bbox_inches=0)
    tii += 1
    #plt.ion()

print subprocess.check_output('convert interp*.png interp00.gif', shell=True) #convert to animated gif 

# DART $Id$
# from Alexey Morozov

# <next few lines under version control, do not edit>
# $URL$
# $Revision$
# $Date$
