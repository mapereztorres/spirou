import numpy as np

from matplotlib.offsetbox import OffsetImage, AnnotationBbox
import matplotlib.patches as mpatches
#matplotlib.rc_file_defaults()
plt.style.use(['bmh','SPIworkflow/spi.mplstyle'])

## All diagnostic plots together 

fig, (ax1, ax2,ax3,ax4) = plt.subplots(4, 1, sharex=True)
fig.subplots_adjust(hspace=0)
ax1.plot(x, v_orb/1e5*np.ones(len(x)), color='r', linestyle='dotted')
ax1.plot(x, v_alf/1e5*np.ones(len(x)), color='g', linestyle='dashdot')
ax1.plot(x, v_sw/1e5*np.ones(len(x)), color='b', linestyle='dashed')
ax1.plot(x, v_rel/1e5*np.ones(len(x)), color='k', linestyle='solid')
ax1.plot(x, v_sound/1e5*np.ones(len(x)), color='orange', linestyle=(0,(1,3.5)))
red_patch = mpatches.Patch(color='red', label=r'v_${\rm orb}$')
black_patch = mpatches.Patch(color='black', label=r'v_${\rm rel}$')
green_patch = mpatches.Patch(color='green', label=r'v_${\rm alf}$')
blue_patch = mpatches.Patch(color='blue', label=r'v_${\rm sw}$')
orange_patch = mpatches.Patch(color='orange', label=r'v_${\rm sound}$')
ax1.set_xscale('log')
ax1.set_yscale('log')          
ax1.legend(handles=[blue_patch,red_patch,black_patch,orange_patch,green_patch],loc='upper left',fontsize=20,facecolor='white',edgecolor='white', framealpha=0)


ax2.plot(x, np.abs(B_r)*np.ones(len(x)), color='r', linestyle='dotted')
ax2.plot(x, B_phi*np.ones(len(x)), color='g', linestyle='dashdot')
ax2.plot(x, B_sw*np.ones(len(x)), color='b', linestyle='solid')
ax2.plot(x, B_sw*np.ones(len(x))*np.sqrt(geom_f), color='k', linestyle='dashed')
red_patch = mpatches.Patch(color='red', label=r'B$_{\rm r}$')
green_patch = mpatches.Patch(color='green', label=r'B$_{\rm phi}$')
blue_patch = mpatches.Patch(color='blue', label=r'B$_{\rm tot}$')
black_patch = mpatches.Patch(color='black', label=r'B$_{\rm perp}$')
ax2.set_xscale('log')
ax2.set_yscale('log')            
ax2.legend(handles=[red_patch,green_patch,blue_patch,black_patch],loc='upper right',fontsize=20,facecolor='white',edgecolor='white', framealpha=0)

ax3.plot(x, M_A*np.ones(len(x)), color='k', lw=lw)
ax3.set_xscale('log')
ax3.set_yscale('log')


ax4.plot(x, P_B_sw*np.ones(len(x)), color='b', linestyle='dashdot')
ax4.plot(x, P_dyn_sw*np.ones(len(x)), color='r', linestyle='solid')
ax4.plot(x, P_th_sw*np.ones(len(x)), color='g', linestyle=(0,(1,3.5)))
# LPM - CHECK IF USEFUL
#ax4.plot(x, P_sw*np.ones(len(x)), color='k', ls=(0, (1, 2)))
#ax4.plot(x, P_B_planet*np.ones(len(x)), color='magenta', ls=(0, (1, 2)))

black_patch = mpatches.Patch(color='black', label=r'P$_{\rm sw}$') #r'Primary T$_{\rm eff}$')
red_patch = mpatches.Patch(color='red', label=r'P$_{\rm dyn_{\rm sw}}$')
green_patch = mpatches.Patch(color='green', label=r'P$_{\rm th_{\rm sw}}$')
blue_patch = mpatches.Patch(color='blue', label=r'P$_{\rm B_{\rm sw}}$')
magenta_patch = mpatches.Patch(color='magenta', label=r'P$_{\rm B_planet}$')             
ax4.set_xscale('log')
ax4.set_yscale('log')  
ax4.legend(handles=[blue_patch,red_patch,green_patch],loc='upper right',fontsize=20,facecolor='white',edgecolor='white', framealpha=0)

ax1.set_ylabel(r"v $[\rm km$ $s^{-1}] $")
ax2.set_ylabel(r"$B_{\rm sw}$ $[G]$")
ax3.set_ylabel(r"$M_A$")
ax4.set_ylabel(r"P $[\rm erg$ $\rm cm^{-3}] $")
        
ax1.set_facecolor("white")
ax2.set_facecolor("white")
ax3.set_facecolor("white")
ax4.set_facecolor("white")

ax1.axvline(x = xnom, ls='--', color='k', lw=2)
ax2.axvline(x = xnom, ls='--', color='k', lw=2)
ax3.axvline(x = xnom, ls='--', color='k', lw=2)
ax4.axvline(x = xnom, ls='--', color='k', lw=2)
        
if STUDY == "D_ORB":
    ax4.set_xlim([2,x[-1]])
if M_A[-1]>1:
    ax3.axhline(y = 1, ls='-.', color='grey', lw=2)   
ax4.set_xlabel(xlabel,fontsize=20)

fig.set_figwidth(8)
fig.set_figheight(20)
#diagnostic_string = "-Bplanet" + str(B_planet_arr[loc_pl]) + "G" +'-'+'T_corona'+str(T_corona/1e6)+'MK'+'-'+'SPI_at_'+str(R_ff_in/R_star)+'R_star'+'.pdf' 
diagnostic_string = "{:.1f}".format(B_star) + "G" + "-Bplanet" + str(B_planet_arr[loc_pl]) + "G" + '-'+"{:.1e}".format(BETA_EFF_MIN)+'-'+"{:.1e}".format(BETA_EFF_MAX)+'-'+'T_corona'+str(T_corona/1e6)+'MK'+'SPI_at_'+str(R_ff_in/R_star)+'R_star'

#if Bfield_geom_arr[ind] == 'open_parker_spiral':
    #out_diagnos =  FOLDER + '/' + STUDY + "_" + str(Exoplanet.replace(" ", "_")) + "-diagnostic" + "-Open-Bstar" + diagnostic_string 
#    geometry = "-Open-spiral-Bstar" 
#elif Bfield_geom_arr[ind]== 'closed_dipole':
    #out_diagnos =  FOLDER + '/' + STUDY + "_" + str(Exoplanet.replace(" ", "_")) + "-diagnostic" + "-Closed-Bstar" + diagnostic_string 
#    geometry = "-Closed-dipole-Bstar"
#else:
#    geometry = "-Closed-PFSS-Bstar"
out_diagnos =  FOLDER + '/' +"diagnostic-" + STUDY + "_" + str(Exoplanet.replace(" ", "_")) +  geometry + diagnostic_string +'.pdf' 
plt.savefig(out_diagnos,bbox_inches='tight')

