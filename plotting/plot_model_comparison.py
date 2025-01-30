import numpy as np

from matplotlib.offsetbox import OffsetImage, AnnotationBbox
import matplotlib.patches as mpatches
from matplotlib.lines import Line2D
#matplotlib.rc_file_defaults()
plt.style.use(['bmh','SPIworkflow/spi.mplstyle'])

##comparison plots for the fluxes
alfven_wing_parker=pd.read_csv(FOLDER+'/CSV/'+'Flux_'+STUDY + "_" + str(Exoplanet.replace(" ", "_")) + '-open-parker-spiral-Bstar'+"{:.1f}".format(B_star)+'G-Bplanet' + str(B_planet_arr[loc_pl]) + 'G' + '-'+"{:.1e}".format(BETA_EFF_MIN)+'-'+"{:.1e}".format(BETA_EFF_MAX)+'-'+'T_corona'+str(T_corona/1e6)+'MK'+'SPI_at_'+str(R_ff_in/R_star)+'R_star'+'_alfven_wing_model.csv')
reconnection_parker=pd.read_csv(FOLDER+'/CSV/'+'Flux_'+STUDY + "_" + str(Exoplanet.replace(" ", "_")) + '-open-parker-spiral-Bstar'+"{:.1f}".format(B_star)+'G-Bplanet' + str(B_planet_arr[loc_pl]) + 'G' + '-'+"{:.1e}".format(BETA_EFF_MIN)+'-'+"{:.1e}".format(BETA_EFF_MAX)+'-'+'T_corona'+str(T_corona/1e6)+'MK'+'SPI_at_'+str(R_ff_in/R_star)+'R_star'+'_reconnection_model.csv')
alfven_wing_dipole=pd.read_csv(FOLDER+'/CSV/'+'Flux_'+STUDY + "_" + str(Exoplanet.replace(" ", "_")) + '-closed-dipole-Bstar'+"{:.1f}".format(B_star)+'G-Bplanet' + str(B_planet_arr[loc_pl]) + 'G' + '-'+"{:.1e}".format(BETA_EFF_MIN)+'-'+"{:.1e}".format(BETA_EFF_MAX)+'-'+'T_corona'+str(T_corona/1e6)+'MK'+'SPI_at_'+str(R_ff_in/R_star)+'R_star'+'_alfven_wing_model.csv')
reconnection_dipole=pd.read_csv(FOLDER+'/CSV/'+'Flux_'+STUDY + "_" + str(Exoplanet.replace(" ", "_")) + '-closed-dipole-Bstar'+"{:.1f}".format(B_star)+'G-Bplanet' + str(B_planet_arr[loc_pl]) + 'G' + '-'+"{:.1e}".format(BETA_EFF_MIN)+'-'+"{:.1e}".format(BETA_EFF_MAX)+'-'+'T_corona'+str(T_corona/1e6)+'MK'+'SPI_at_'+str(R_ff_in/R_star)+'R_star'+'_reconnection_model.csv')
#alfven_wing_hybrid=pd.read_csv(FOLDER+'/'+STUDY + "_" + str(Exoplanet.replace(" ", "_")) + '-hybrid-Bstar'+"{:.1f}".format(B_star)+'G-Bplanet' + str(B_planet_arr[loc_pl]) + 'G' + '-'+"{:.1e}".format(BETA_EFF_MIN)+'-'+"{:.1e}".format(BETA_EFF_MAX)+'-'+'T_corona'+str(T_corona/1e6)+'MK'+'SPI_at_'+str(R_ff_in/R_star)+'R_star'+'_alfven_wing_model.csv')
#reconnection_hybrid=pd.read_csv(FOLDER+'/'+STUDY + "_" + str(Exoplanet.replace(" ", "_")) + '-hybrid-Bstar'+"{:.1f}".format(B_star)+'G-Bplanet' + str(B_planet_arr[loc_pl]) + 'G' + '-'+"{:.1e}".format(BETA_EFF_MIN)+'-'+"{:.1e}".format(BETA_EFF_MAX)+'-'+'T_corona'+str(T_corona/1e6)+'MK'+'SPI_at_'+str(R_ff_in/R_star)+'R_star'+'_reconnection_model.csv')
alfven_wing_pfss_parker=pd.read_csv(FOLDER+'/CSV/'+'Flux_'+STUDY + "_" + str(Exoplanet.replace(" ", "_")) + '-pfss-parker-Bstar'+"{:.1f}".format(B_star)+'G-Bplanet' + str(B_planet_arr[loc_pl]) + 'G' + '-'+"{:.1e}".format(BETA_EFF_MIN)+'-'+"{:.1e}".format(BETA_EFF_MAX)+'-'+'T_corona'+str(T_corona/1e6)+'MK'+'SPI_at_'+str(R_ff_in/R_star)+'R_star'+'_alfven_wing_model.csv')
reconnection_pfss_parker=pd.read_csv(FOLDER+'/CSV/'+'Flux_'+STUDY + "_" + str(Exoplanet.replace(" ", "_")) + '-pfss-parker-Bstar'+"{:.1f}".format(B_star)+'G-Bplanet' + str(B_planet_arr[loc_pl]) + 'G' + '-'+"{:.1e}".format(BETA_EFF_MIN)+'-'+"{:.1e}".format(BETA_EFF_MAX)+'-'+'T_corona'+str(T_corona/1e6)+'MK'+'SPI_at_'+str(R_ff_in/R_star)+'R_star'+'_reconnection_model.csv')
plt.figure(figsize=(8, 7.5))
lw=2
ax2 = plt.subplot2grid((1, 1), (0, 0), rowspan=1, colspan=1)
ax2.set_facecolor("white")    
ax2.set_ylim([YLIMLOW * 10**-2, YLIMHIGH * 10**3])     
ax2.set_yscale('log')
ax2.plot(alfven_wing_parker[STUDY], alfven_wing_parker['FLUX'], color='orange')
ax2.plot(reconnection_parker[STUDY], reconnection_parker['FLUX'], color='blue')
ax2.plot(alfven_wing_dipole[STUDY], alfven_wing_dipole['FLUX'], color='orange', linestyle='dotted')
ax2.plot(reconnection_dipole[STUDY], reconnection_dipole['FLUX'], color='blue', linestyle='dotted')
#ax2.plot(alfven_wing_hybrid[STUDY], alfven_wing_hybrid['FLUX'], color='orange', linestyle='dashed')
#ax2.plot(reconnection_hybrid[STUDY], reconnection_hybrid['FLUX'], color='blue', linestyle='dashed')
ax2.plot(alfven_wing_pfss_parker[STUDY], alfven_wing_pfss_parker['FLUX'], color='orange', linestyle='dashed')
ax2.plot(reconnection_pfss_parker[STUDY], reconnection_pfss_parker['FLUX'], color='blue', linestyle='dashed')
legend_elements = [
Line2D([0], [0], color='orange', lw=lw, label='Alfvén wing (Open Parker Spiral)'),
Line2D([0], [0], color='blue', lw=lw, label='Reconnection (Open Parker Spiral)'),
Line2D([0], [0], color='orange', linestyle='dotted', lw=lw, label='Alfvén wing (Dipole)'),
Line2D([0], [0], color='blue', linestyle='dotted', lw=lw, label='Reconnection (Dipole)'),
#Line2D([0], [0], color='orange', linestyle='dashed', lw=2, label='Alfvén wing (Hybrid)'),
#Line2D([0], [0], color='blue', linestyle='dashed', lw=2, label='Reconnection (Hybrid)'),
Line2D([0], [0], color='orange', linestyle='dashed', lw=lw, label='Alfvén wing (PFSS-Parker)'),
Line2D([0], [0], color='blue', linestyle='dashed', lw=lw, label='Reconnection (PFSS-Parker)'),
]
ax2.set_xlabel(xlabel,fontsize=20)
ax2.set_ylabel(r"Flux density [mJy]")
ax2.legend(handles=legend_elements, loc='lower left', fontsize=12, facecolor='white', edgecolor='white', framealpha=0)
ax2.set_xlim(left=2)     

plt.savefig(FOLDER+'/'+'Flux'+'_model_comparison-'+ STUDY + "_" + str(Exoplanet.replace(" ", "_")) + '-Bstar'+"{:.1f}".format(B_star)+'G-Bplanet' + str(B_planet_arr[loc_pl]) + 'G' + '-'+"{:.1e}".format(BETA_EFF_MIN)+'-'+"{:.1e}".format(BETA_EFF_MAX)+'-'+'T_corona'+str(T_corona/1e6)+'MK'+'SPI_at_'+str(R_ff_in/R_star)+'R_star'+'.pdf')


##comparison plots for the stellar wind magnetic field
FOLDER + '/CSV/' +"diagnostic-" + STUDY + "_" + str(Exoplanet.replace(" ", "_")) +  geometry + diagnostic_string +'_B_sw.csv'
diagnostic_string = "{:.1f}".format(B_star) + "G" + "-Bplanet" + str(B_planet_arr[loc_pl]) + "G" + '-'+"{:.1e}".format(BETA_EFF_MIN)+'-'+"{:.1e}".format(BETA_EFF_MAX)+'-'+'T_corona'+str(T_corona/1e6)+'MK'+'SPI_at_'+str(R_ff_in/R_star)+'R_star'

bsw_parker=pd.read_csv(FOLDER + '/CSV/' +'diagnostic-'+ STUDY + "_" + str(Exoplanet.replace(" ", "_")) + '-open-parker-spiral-Bstar'+"{:.1f}".format(B_star) + "G" + "-Bplanet" + str(B_planet_arr[loc_pl]) + "G" + '-'+"{:.1e}".format(BETA_EFF_MIN)+'-'+"{:.1e}".format(BETA_EFF_MAX)+'-'+'T_corona'+str(T_corona/1e6)+'MK'+'SPI_at_'+str(R_ff_in/R_star)+'R_star'+'_B_sw.csv')

bsw_dipole=pd.read_csv(FOLDER + '/CSV/' +'diagnostic-'+ STUDY + "_" + str(Exoplanet.replace(" ", "_")) + '-closed-dipole-Bstar'+"{:.1f}".format(B_star) + "G" + "-Bplanet" + str(B_planet_arr[loc_pl]) + "G" + '-'+"{:.1e}".format(BETA_EFF_MIN)+'-'+"{:.1e}".format(BETA_EFF_MAX)+'-'+'T_corona'+str(T_corona/1e6)+'MK'+'SPI_at_'+str(R_ff_in/R_star)+'R_star'+'_B_sw.csv')

#bsw_hybrid=pd.read_csv(FOLDER + '/' +'diagnostic-'+ STUDY + "_" + str(Exoplanet.replace(" ", "_")) + '-hybrid-Bstar'+"{:.1f}".format(B_star) + "G" + "-Bplanet" + str(B_planet_arr[loc_pl]) + "G" + '-'+"{:.1e}".format(BETA_EFF_MIN)+'-'+"{:.1e}".format(BETA_EFF_MAX)+'-'+'T_corona'+str(T_corona/1e6)+'MK'+'SPI_at_'+str(R_ff_in/R_star)+'R_star'+'_B_sw.csv')

bsw_pffs_parker=pd.read_csv(FOLDER + '/CSV/' +'diagnostic-'+ STUDY + "_" + str(Exoplanet.replace(" ", "_")) + '-pfss-parker-Bstar'+"{:.1f}".format(B_star) + "G" + "-Bplanet" + str(B_planet_arr[loc_pl]) + "G" + '-'+"{:.1e}".format(BETA_EFF_MIN)+'-'+"{:.1e}".format(BETA_EFF_MAX)+'-'+'T_corona'+str(T_corona/1e6)+'MK'+'SPI_at_'+str(R_ff_in/R_star)+'R_star'+'_B_sw.csv')

plt.figure(figsize=(8, 7.5))
lw=2
ax2 = plt.subplot2grid((1, 1), (0, 0), rowspan=1, colspan=1)
ax2.set_facecolor("white")  
ax2.set_xscale('log')
ax2.set_yscale('log')
ax2.plot(bsw_parker[STUDY], bsw_parker['Bsw'], color='blue', linestyle='dashed')
ax2.plot(bsw_dipole[STUDY], bsw_dipole['Bsw'], color='orange', linestyle='dotted')
#ax2.plot(bsw_hybrid[STUDY], bsw_hybrid['Bsw'], color='blue', linestyle='dashed')
ax2.plot(bsw_pffs_parker[STUDY], bsw_pffs_parker['Bsw'], color='black')
legend_elements = [
Line2D([0], [0], color='blue', linestyle='dashed', lw=lw, label='Open Parker Spiral'),
Line2D([0], [0], color='orange', linestyle='dotted', lw=lw, label='Dipole'),
#Line2D([0], [0], color='blue', linestyle='dashed', lw=2, label='Hybrid'),
Line2D([0], [0], color='black', lw=lw, label='PFSS-Parker'),
]

ax2.legend(handles=legend_elements, loc='lower left', fontsize=12, facecolor='white', edgecolor='white', framealpha=0)
ax2.set_xlim(left=2)     
ax2.set_xlabel(xlabel,fontsize=20)
ax2.set_ylabel(r"$B_{\rm sw}$ $[G]$")
plt.savefig(FOLDER + '/' +'B_sw_'+'model_comparison-'+ STUDY + "_" + str(Exoplanet.replace(" ", "_")) + '-hybrid-Bstar'+"{:.1f}".format(B_star) + "G" + "-Bplanet" + str(B_planet_arr[loc_pl]) + "G" + '-'+"{:.1e}".format(BETA_EFF_MIN)+'-'+"{:.1e}".format(BETA_EFF_MAX)+'-'+'T_corona'+str(T_corona/1e6)+'MK'+'SPI_at_'+str(R_ff_in/R_star)+'R_star'+'.pdf')




##comparison plots for M_A
bsw_parker=pd.read_csv(FOLDER + '/CSV/' +'diagnostic-'+ STUDY + "_" + str(Exoplanet.replace(" ", "_")) + '-open-parker-spiral-Bstar'+"{:.1f}".format(B_star) + "G" + "-Bplanet" + str(B_planet_arr[loc_pl]) + "G" + '-'+"{:.1e}".format(BETA_EFF_MIN)+'-'+"{:.1e}".format(BETA_EFF_MAX)+'-'+'T_corona'+str(T_corona/1e6)+'MK'+'SPI_at_'+str(R_ff_in/R_star)+'R_star'+'_M_A.csv')

bsw_dipole=pd.read_csv(FOLDER + '/CSV/' +'diagnostic-'+ STUDY + "_" + str(Exoplanet.replace(" ", "_")) + '-closed-dipole-Bstar'+"{:.1f}".format(B_star) + "G" + "-Bplanet" + str(B_planet_arr[loc_pl]) + "G" + '-'+"{:.1e}".format(BETA_EFF_MIN)+'-'+"{:.1e}".format(BETA_EFF_MAX)+'-'+'T_corona'+str(T_corona/1e6)+'MK'+'SPI_at_'+str(R_ff_in/R_star)+'R_star'+'_M_A.csv')

#bsw_hybrid=pd.read_csv(FOLDER + '/' +'diagnostic-'+ STUDY + "_" + str(Exoplanet.replace(" ", "_")) + '-hybrid-Bstar'+"{:.1f}".format(B_star) + "G" + "-Bplanet" + str(B_planet_arr[loc_pl]) + "G" + '-'+"{:.1e}".format(BETA_EFF_MIN)+'-'+"{:.1e}".format(BETA_EFF_MAX)+'-'+'T_corona'+str(T_corona/1e6)+'MK'+'SPI_at_'+str(R_ff_in/R_star)+'R_star'+'_M_A.csv')

bsw_pffs_parker=pd.read_csv(FOLDER + '/CSV/' +'diagnostic-'+ STUDY + "_" + str(Exoplanet.replace(" ", "_")) + '-pfss-parker-Bstar'+"{:.1f}".format(B_star) + "G" + "-Bplanet" + str(B_planet_arr[loc_pl]) + "G" + '-'+"{:.1e}".format(BETA_EFF_MIN)+'-'+"{:.1e}".format(BETA_EFF_MAX)+'-'+'T_corona'+str(T_corona/1e6)+'MK'+'SPI_at_'+str(R_ff_in/R_star)+'R_star'+'_M_A.csv')


plt.figure(figsize=(8, 7.5))
lw=2
ax2 = plt.subplot2grid((1, 1), (0, 0), rowspan=1, colspan=1)
ax2.set_facecolor("white")  
ax2.set_xscale('log')
ax2.set_yscale('log')
ax2.plot(bsw_parker[STUDY], bsw_parker['M_A'], color='blue')
ax2.plot(bsw_dipole[STUDY], bsw_dipole['M_A'], color='orange', linestyle='dotted')
#ax2.plot(bsw_hybrid[STUDY], bsw_hybrid['M_A'], color='black', linestyle='dashed')
ax2.plot(bsw_pffs_parker[STUDY], bsw_pffs_parker['M_A'], color='black', linestyle='dashed')
legend_elements = [
Line2D([0], [0], color='blue', lw=lw, label='Open Parker Spiral'),
Line2D([0], [0], color='orange', linestyle='dotted', lw=lw, label='Dipole'),
#Line2D([0], [0], color='black', linestyle='dashed', lw=2, label='Hybrid'),
Line2D([0], [0], color='black', linestyle='dashed', lw=lw, label='PFSS-Parker'),
]

ax2.legend(handles=legend_elements, loc='upper left', fontsize=12, facecolor='white', edgecolor='white', framealpha=0)
ax2.set_xlim(left=2)     
ax2.set_xlabel(xlabel,fontsize=20)
ax2.set_ylabel(r"$M_A$")
plt.savefig(FOLDER + '/'+'M_A_'+'model_comparison-'+ STUDY + "_" + str(Exoplanet.replace(" ", "_")) + '-hybrid-Bstar'+"{:.1f}".format(B_star) + "G" + "-Bplanet" + str(B_planet_arr[loc_pl]) + "G" + '-'+"{:.1e}".format(BETA_EFF_MIN)+'-'+"{:.1e}".format(BETA_EFF_MAX)+'-'+'T_corona'+str(T_corona/1e6)+'MK'+'SPI_at_'+str(R_ff_in/R_star)+'R_star'+'.pdf')






#print(FOLDER + '/' +'diagnostic-'+ STUDY + "_" + str(Exoplanet.replace(" ", "_")) + '-pfss-parker-Bstar'+"{:.1f}".format(B_star) + "G" + "-Bplanet" + str(B_planet_arr[loc_pl]) + "G" + '-'+"{:.1e}".format(BETA_EFF_MIN)+'-'+"{:.1e}".format(BETA_EFF_MAX)+'-'+'T_corona'+str(T_corona/1e6)+'MK'+'SPI_at_'+str(R_ff_in/R_star)+'R_star'+'_B_sw.csv')
'''

'''
