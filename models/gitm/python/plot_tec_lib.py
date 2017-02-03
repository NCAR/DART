import numpy as np #for data structures
import matplotlib
matplotlib.use('Agg') #so as not to need a display or X-server
import matplotlib.pyplot as plt #for plotting
from matplotlib.patches import Polygon
from mpl_toolkits.mplot3d import Axes3D #for 3d plotting


def time_plot(timeT,VarT,lT, timeEr,VarEr,VarEo,VarsdEr,VarsdEo,lE,lsdE, xli,yli,yla, fn):
    """
    plot some quantity coming from data and estimated on the same plot as functions of time.
    
    INPUTS:
    - timeT[t] - truth time: np.array of datetime objects signifying the time instances when data (or truth) were taken
    - VarT[t] - truth variable: np.array over [time] of data
    - lt - label to be used for the true quantity in legend

    - timeEr[t] - ensemble time (r for prior, but hopefully the prior time is equal to posterior): np.array of datetime objects for ensemble members
    - VarEr[t] - ensemble (or estimated) variable prior: np.array over [time] of values prior to assimilation
    - VarEo[t] - ensemble (or estimated) variable posterior: np.array over [time] of values after the assimilation
    - VarsdEr[t] - ensemble (or estimated) variable uncertainty (standard deviation, ensemble spread) prior: np.array over [time] of prior sd's
    - VarsdEo[t] - ensemble (or estimated) variable standard deviation posterior: np.array over [time] of posterior sd's
    - lE - label to be used for the estimated quantity in legend

    - xli - xlim
    - yli - ylim
    - yla - ylabel

    - fn - filename without extension, where this plot will be saved as .png and as .eps


    OUTPUTS:
    - fn + '.png' - the plot saved as a png
    - fn + '.eps' - the plot saved as an eps

    EXAMPLES:
    t=np.array([datetime.datetime(2002,12,1), datetime.datetime(2002,12,2)])
    plots_tec_alex.time_plot(t,np.array([1,1]),'lT', t,np.array([2,2]),np.array([2.1,2.1]),np.array([.1,.1]),np.array([.2,.2]),'lE','lsdE', (-1,25),(0,3),'yla', 'hi')
    """
    VarE=np.zeros(len(VarEr)*2,float)
    VarE[0:-1:2]=VarEr
    VarE[1::2]=VarEo
    VarsdE=np.zeros(len(VarsdEr)*2,float)
    VarsdE[0:-1:2]=VarsdEr
    VarsdE[1::2]=VarsdEo
    timeE=np.zeros(len(timeEr)*2,object)
    timeE[0:-1:2]=timeEr
    timeE[1::2]=timeEr
    tE = np.array([(t.day-1)*24 + t.hour + float(t.minute)/60.0 for t in timeE])
    tT = np.array([(t.day-1)*24 + t.hour + float(t.minute)/60.0 for t in timeT])
#f107 plot
    fig = plt.figure(1)
    fig.clf() #in case it was already open
    fig.set_size_inches((8,4))
    ax = fig.add_subplot(111)
    x = np.float32(np.concatenate(( tE, tE[::-1]))) 
    y = np.concatenate(( VarE+VarsdE, (VarE-VarsdE)[::-1]))
    p = np.zeros((len(x),2), np.float32)
    p[:,0] = x
    p[:,1] = y
    po = Polygon(p, True, alpha=0.3, color=(0.1, 0.6, 1), label=lsdE) #True appears as an answer to "Close the polygon?" 
    fsdEo = ax.add_patch(po)
    
    if lT!='ignore':
        fT, = ax.plot(tT, VarT, '--k', label=lT)
    
    fE, = ax.plot(tE, VarE,   'b', label=lE)

#labels, etc
    plt.rcParams['legend.fontsize'] = 10
    h, l = ax.get_legend_handles_labels()
    ax.legend(h, l, loc=1)
#ax.set_title(str(ti))
#ax.set_xticks(np.linspace(0,360,7) );
#ax.set_yticks(np.linspace(-90,90,7) );
    ax.grid(True)
    ax.set_xticks(np.linspace(0,24,24/3+1) );
    ax.set_xlim(xli)
    ax.set_ylim(yli)
    ax.set_xlabel('Hours since 00UT 12/1/2002');
    ax.set_ylabel(yla);
#ax.text(0, 0, str(tcommon[ti]), fontsize=10,  transform=ax.transAxes)
    fig.canvas.draw()    
    plt.savefig(fn + '.png', bbox_inches=0)
    plt.savefig(fn + '.eps', bbox_inches=0)
#plt.ion()

def time_plot_alt(timeT,VarTa,AltT,lT, timeEr,VarEra,VarEoa,VarsdEra,VarsdEoa,lE,lsdE, xli,yli,yla, fn):
    """
    same thing as time_plot, but N row-subplots for altitude, or really any other dimension
    
    INPUTS:
    - timeT[t] - truth time: np.array of datetime objects signifying the time instances when data (or truth) were taken
    - VarTa[a,t] - truth variable: np.array over [alt,time] of data
    - AltT[a] - the coordinate axis over which we are iterating [alt]
    - lt - label to be used for the true quantity in legend

    - timeEr[t] - ensemble time (r for prior, but hopefully the prior time is equal to posterior): np.array of datetime objects for ensemble members
    - VarEr[a,t] - ensemble (or estimated) variable prior: np.array over [alt,time] of values prior to assimilation
    - VarEo[a,t] - ensemble (or estimated) variable posterior: np.array over [alt,time] of values after the assimilation
    - VarsdEr[a,t] - ensemble (or estimated) variable uncertainty (standard deviation, ensemble spread) prior: np.array over [time] of prior sd's
    - VarsdEo[a,t] - ensemble (or estimated) variable standard deviation posterior: np.array over [time] of posterior sd's
    - lE - label to be used for the estimated quantity in legend

    - xli - xlim
    - yli - ylim
    - yla - ylabel

    - fn - filename without extension, where this plot will be saved as .png and as .eps


    OUTPUTS:
    - fn + '.png' - the plot saved as a png
    - fn + '.eps' - the plot saved as an eps

    EXAMPLES:
    plots_tec_alex.time_plot_alt(timeT, EdsT[lona,lata,:,:], AltT, '$D_{e,aa}$ measured', timeE, EdsEr[lona,lata,:,:], EdsEo[lona,lata,:,:], EdssdEr[lona,lata,:,:], EdssdEo[lona,lata,:,:], '$D_{e,aa}$ estimated', '$D_{e,aa}$ estimated $\pm$ SD', (0,24), (1e8,1e12), '$D_{e,AA}$  [$e$ $m^{-3}$]', 'vtec_daa_alt')
    """
    fig = plt.figure(1)
    fig.clf() #in case it was already open
    fig.set_size_inches((8,50))
    for iAlt in range(len(VarTa[:,0])-1,0,-1):
        #print iAlt
        VarT = VarTa[iAlt,:]
        VarEr = VarEra[iAlt,:]
        VarEo = VarEoa[iAlt,:]
        VarsdEr = VarsdEra[iAlt,:]
        VarsdEo = VarsdEoa[iAlt,:]
        VarE=np.zeros(len(VarEr)*2,float)
        VarE[0:-1:2]=VarEr
        VarE[1::2]=VarEo
        VarsdE=np.zeros(len(VarsdEr)*2,float)
        VarsdE[0:-1:2]=VarsdEr
        VarsdE[1::2]=VarsdEo
        timeE=np.zeros(len(timeEr)*2,object)
        timeE[0:-1:2]=timeEr
        timeE[1::2]=timeEr
        tE = np.array([(t.day-1)*24 + t.hour + float(t.minute)/60.0 for t in timeE])
        tT = np.array([(t.day-1)*24 + t.hour + float(t.minute)/60.0 for t in timeT])
        
        ax = fig.add_subplot(len(VarTa[:,0]),1,len(VarTa[:,0])-iAlt)
        #ensemble standard deviation patch
        x = np.float32(np.concatenate(( tE, tE[::-1]))) 
        y = np.concatenate(( VarE+VarsdE, (VarE-VarsdE)[::-1]))
        p = np.zeros((len(x),2), np.float32)
        p[:,0] = x
        p[:,1] = y
        po = Polygon(p, True, alpha=0.3, color=(0.1, 0.6, 1), label=lsdE) #True appears as an answer to "Close the polygon?" 
        fsdEo = ax.add_patch(po)
        #truth line
        fT, = ax.plot(tT, VarT, '--k', label=lT)
        #ensemble mean line
        fE, = ax.plot(tE, VarE,   'b', label=lE)

#labels, etc
        plt.rcParams['legend.fontsize'] = 10
        h, l = ax.get_legend_handles_labels()        
#ax.set_title(str(ti))
#ax.set_xticks(np.linspace(0,360,7) );
#ax.set_yticks(np.linspace(-90,90,7) );
        ax.set_yscale('log')
        ax.text(.01, .02, 'Alt = ' + ('%3.0f' % (AltT[iAlt]/1000)), fontsize=10,  transform=ax.transAxes)  #+ ('%02i' % iAlt ) +','
        ax.grid(True)
        ax.set_xticks(np.linspace(0,24,24/3+1) );
        ax.set_xlim(xli)
        ax.set_ylim(yli)
        #ax.set_ylabel(yla);
#ax.text(0, 0, str(tcommon[ti]), fontsize=10,  transform=ax.transAxes)
        #fig.canvas.draw()    
        fig.subplots_adjust(hspace = .001)
        ax.set_xticklabels(())
    
    ax.legend(h, l, loc=4)
    ax.set_xticklabels(np.linspace(0,24,24/3+1) );
    ax.set_xlabel('Hours since 00UT 12/1/2002');    
    plt.savefig(fn + '.png', bbox_inches=0)
    plt.savefig(fn + '.eps', bbox_inches=0)
#plt.ion()

# DART $Id$
# from Alexey Morozov

# <next few lines under version control, do not edit>
# $URL$
# $Revision$
# $Date$
