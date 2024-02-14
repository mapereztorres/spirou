# # Star-planet Interaction (model)
# ## Sub-Alfvenic flux for both a Parker spiral magnetic field configuration and a closed dipolar field configuration


# import statements
import os
import shutil
import numpy as np
import pandas as pd
import matplotlib
import matplotlib.pyplot as plt
from matplotlib.offsetbox import OffsetImage, AnnotationBbox
from scipy.special import lambertw

#matplotlib.rc_file_defaults()
plt.style.use(['bmh','SPIworkflow/spi.mplstyle'])

from output import OutputWriter

### Getting parameters to predict SPI radio emission
#
# Observing parameters: Observin frequency, assumed rms noise, Delta_nu
# Stellar parameters: T_corona, n_base_corona, B_field, Isothermality of plasma
# Geometry of the sub-Alfvénic interaction: alpha, theta_M
#
from SPIworkflow.__init__ import *

# Import useful constants and functions to be used in the code
from SPIworkflow.constants import *
import SPIworkflow.SPIutils as spi

from SPIworkflow.load_data import get_spi_data, create_data_tables, load_target, define_lists


# Create output directory for the results 
# Return df_planets and df_no_planets
# Create CARMENES tables for targets 
# with planets only, and with no planets, unless those already exist
# 
outdir, df_planets, df_no_noplanets = create_data_tables()

# Read in the input data to estimate radio emission from SPI
#from SPIworkflow.import_table import *
#source_data = './INPUT/SPI-sources_planets_MASTER.csv'
#source_data = './INPUT/SPI-sources-sample5.csv'
#
if selection_criteria == False:
     data = get_spi_data(infile_data=source_data)
else:
     data = get_spi_data(infile_data=source_data,distance_max=15, p_orb_max = 10, bfield_min=100,bfield_max=1000.0, dec_min=-90)


# Setting the stellar magnetic field geometry and the value of the 
# intensity of the planetary magnetic field
# 
# Stellar magnetic field geometry
# The convention is that Bfield_geom_arr = 1 => open Parker spiral geometry; 
#                        Bfield_geom_arr = 0 - closed dipolar geometry
#Bfield_geom_arr = [0]

# magnetized_pl_arr is like a False/True array, magfield_planet is the modulus of the magnetic field
#magnetized_pl_arr= [0]
#Bfield_pl = 0.5
planet_name_list,star_name_list, d_star_list, mass_star_list, radius_star_list, p_rot_list,  bfield_star_list, a_list,  p_orb_list,eccentricity_list,q_list,Q_list, mass_planet_list, radius_planet_list,T_cor_list, m_dot_list, nu_pl_list, nu_cycl_list, rho_pl_list, B_pl_list, B_sw_list, v_alf_list, M_A_list, Flux_r_S_ZL_min_list, Flux_r_S_ZL_max_list, P_Bpl_list, P_dyn_list, P_th_list, P_Bsw_list, Rmp_list=define_lists()
### Select data in the array to run the code on
print(source_data)
print(data)

### Table initialization
# 
# Create column for M_star_dot to fill it with values
data['M_star_dot(M_sun_dot)']=''
# If bfield_star(gauss) is missing, set it to np.nan
data['bfield_star(gauss)'].replace('', np.nan, inplace=True)
# If p_rot is missing, set it to np.nan
#data['p_rot(days)'].replace('', np.nan, inplace=True)
# Remove targets without p_rot
data.dropna(subset=['p_rot(days)'], inplace=True)
# Do not use stars with P_rot smaller than 10 days
data = data[data['p_rot(days)'] > 10.0]
data['radius_planet(r_earth)'].replace('', np.nan, inplace=True)
data.reset_index(inplace=True) # to prevent funny jumps in the indices



if COMPUTE_ALL:
    planet_array = range(round(len(data)/1))
else:
    planet_array = compute_planets
    
for indi in planet_array:
    starname,d, R_star, M_star, P_rot_star, B_star, Exoplanet, Mp, Rp, r_orb, P_orb,eccentricity, q, Q = load_target(data, indi)

    # Fill B_star column if empty. Uses original units from table
    if pd.isna(B_star):
        data['bfield_star(gauss)'][indi] = spi.B_starmass(star_mass=data['mass_star(m_sun)'][indi],Prot=data['p_rot(days)'][indi])
        # Eventually include uncertainties in B _star
        data['e_bfield_star'][indi]='TBD'
    B_star = data['bfield_star(gauss)'][indi]            # Stellar surface magnetic field

    data['M_star_dot(M_sun_dot)'][indi] = spi.Mdot_star(R_star=data['radius_star(r_sun)'][indi],
        M_star=data['mass_star(m_sun)'][indi], Prot_star=data['p_rot(days)'][indi])/M_sun_dot
    M_star_dot = data['M_star_dot(M_sun_dot)'][indi]     # Stellar mass loss rate in solar units 
    print('M_star_dot :',M_star_dot)

    # Electron gyrofrequency and ECM bandwidth 
    gyrofreq = e*B_star/(2*np.pi * m_e * c) # in cgs units
    Delta_nu_cycl = 0.5 * gyrofreq # Hz - width of ECMI emission  assumed to be  (0.5 * gyrofreq), 

    # Common properties for star and planet
    # 
    M_star_msun = M_star / M_sun # Stellar mass in units of solar mass
    Omega_star = 2.0*np.pi / P_rot_star # Angular rotation velocity of the star

    d_orb_max = r_orb/R_star  + 10 # Max. orbital distance, in units of R_star
    print('R_star = ', R_star)
    print('r_orb = ', r_orb)
    print('P_orb = ', P_orb)
    Nsteps = int(2*d_orb_max)

    #d_orb = np.linspace(1.002, 10, Nsteps) * R_star # Array of (orbital) distances to the star
    if sweep=="RAD":
      d_orb = np.linspace(1.02, d_orb_max, Nsteps) * R_star # Array of (orbital) distances to the star, in cm 
    else:
      d_orb=[r_orb]
    #d_orb = np.linspace(1.02, 210, Nsteps) * R_star # Array of (orbital) distances to the star
    #print(len(d_orb))
    v_orb = (G * M_star/d_orb)**0.5 # Orbital (Keplerian) speed of planet as f(distance to star), in cm/s
    v_corot = d_orb * Omega_star # Corotation speed (cm/s)

    #Omega_planet = np.ones(len(d_orb)) * Omega_earth # array of angular speeds of the planet, in  s^(-1)
    Omega_planet =  v_orb / d_orb # Angular speed of the planet, in s^(-1). NOte that it's an array

    # Get wind composition, from the fraction of protons
    X_e, mu, m_av = spi.wind_composition(X_p)

    # Compute stellar wind velocity at each value of d_orb
    v_sound, r_sonic, v_sw = spi.v_stellar_wind(d_orb, M_star, T_corona, m_av)

    v_sw_base = v_sw[0]    # Stellar wind velocity at the closest distance to the star
     

    # Plasma number density at base of the corona
    n_base_corona = spi.n_wind(M_star_dot, R_star, v_sw_base, m_av) 

    # Maximum plasma frequency at the base of the corona. If the ECM
    # freq is less than the plasma frequency, the emission is
    # completely absorbed 
    nu_plasma_corona = spi.plasma_freq(n_base_corona) # in Hz

    #print("V_sound = {0:.3f} km/s; V_sw at the base = {1:.3f} km/s".format(vsound/1e5, v_sw_base/1e5))    
    
    # Eq. 23 of Turnpenney+18 - Second term of RHS 
    # The vector v_rel = v_sw - v_orb (Eq. 5 in Saur+13, and see also Fig. 1 in Turnpenney+18)
    # 
    v_rel = np.sqrt(v_orb**2 + v_sw**2) # Relative speed between stellar wind and obstacle
    v_rel_angle = np.arctan(v_orb/v_sw) # Angle between radial vector and relative velocity
    
    # n_sw_planet - Number density of the wind at orbital distance to the planet. 
    # If the stellar plasma is assumed to be isothermal, then 
    # the density falls down as ~ R^(-2) * v_sw^(-1).
    # Alternatively, we fix the density at the distance of the planet from the host star.
    if isothermal:
        #n_sw_planet = n_sw_base / (d_orb/R_star)**2 / (v_sw/v_sw_base) # Plasma density at distance (R/R_star)
        n_sw_planet = spi.n_wind(M_star_dot, d_orb, v_sw, m_av) # Plasma number density at distance (R/R_star)
    else:
        n_sw_planet = np.ones(len(d_orb)) * 1e4  # fixed                 

    rho_sw_planet = m_av * n_sw_planet #wind density at the distance to the planet, in g * cm^(-3)

    for ind in range(len(Bfield_geom_arr)):
        for ind1 in range(len(magnetized_pl_arr)):
            if Bfield_geom_arr[ind]:
                print("\nOpen Parker magnetic field geometry")
            else:
                print("\nClosed dipolar magnetic field geometry")
           
            # Magnetic field geometry
            # open_field - defines the geometry of the magnetic field
            open_field = Bfield_geom_arr[ind]
            B_r, B_phi, B_sw, B_ang, theta, geom_f = spi.get_bfield_comps(open_field, B_star, d_orb, R_star, v_corot, v_sw, v_rel_angle)
            
            # Compute Alfvén parameters in the stellar wind at a distance d_orb 
            v_alf, M_A, v_alf_r, M_A_radial = spi.get_alfven(rho_sw_planet, B_sw, B_r, v_rel, v_sw)

            # defines whether planet is unmagnetized (magnetized_pl_arr[ind1] = 0), or magnetized (magnetized_pl_arr[ind1] = 1)
            if magnetized_pl_arr[ind1]: # magnetized planet
                print('Magnetized planet\n')
                if B_planet_law == 'Sano':
                    # Planetary magnetic field, using Sano's (1993) scaling law, in units of B_earth 
                    # Assumes a tidally locked planet, i.e., the rotation period of the
                    # planet equals its orbital one. 
                    # WARNING: For small rotation periods, the inferred magnetic field
                    # is too large to be reliable at all.
                    r_core, rho_core, magn_moment_planet, B_planet = spi.bfield_sano(M_planet = Mp / M_earth, 
                                               R_planet = Rp / R_earth, 
                                               Omega_rot_planet = Omega_planet / Omega_earth)  
                    print('B_planet after applying Sano scaling law:', B_planet)
                    B_planet *= bfield_earth  # B_planet, in Tesla
                else: 
                    B_planet = B_planet_default  # B_planet, in Tesla
                B_planet    *=  Tesla2Gauss #  B_planet, in Gauss 
                B_planet_arr = np.ones(len(d_orb)) * B_planet  
                
            else:  # unmagnetized planet
                print('Unmagnetized planet\n')
                B_planet_arr  = np.zeros(len(d_orb)) # unmagnetized planet
            print('format of B_planet_arr', type(B_planet_arr)) 
            #
            # Effective radius of the obstacle
            # Case 1. À la Saur+2013. 
            R_planet_eff_Saur = spi.get_Rmp_Saur(Rp, theta_M, B_planet_arr, B_sw)


            # Case 2. À la Zarka (2007), Turnpenney+2018, etc.
            #
            # Compute radius of magnetopause, Rmp as balance of wind and planet's
            # pressures
            
            # Planet pressure - only the magnetic component is considered
            P_B_planet  = spi.get_P_B_planet(B_planet_arr) 

            # Stellar wind pressure
            P_dyn_sw = spi.get_P_dyn_sw(n_sw_planet, mu, v_rel) 
            P_th_sw  = spi.get_P_th_sw(n_sw_planet, T_corona)
            P_B_sw   = spi.get_P_B_sw(B_sw)

            # Radius of magnetopause, in cm
            Rmp = spi.get_Rmp(P_B_planet, P_dyn_sw, P_th_sw, P_B_sw) * Rp

            # The effective radius is the radius of the magnetopause
            R_planet_eff = Rmp

            R_planet_eff[ R_planet_eff < Rp] = Rp # R_planet_eff cannot be smaller than Rp

            # Compute min and max speed of electrons emitting via ECM, in units of the speed of light 
            beta_min = spi.beta_keV(Ekin_min);  beta_max = spi.beta_keV(Ekin_max)
            
            # Beam solid angle covered by the ECM emission, in sterradians
            if which_beam_solid_angle == 'Jupiter-Io':
                bsa_Omega = 1.6  # Value obtained from the DAM emission from Jupiter-Io 
                Omega_min = bsa_Omega; Omega_max = bsa_Omega
            else:
                # Get the minimum and maximum values of the beam solid angle for the cone of emission 
                Omega_min, Omega_max = spi.beam_solid_angle(beta_min, beta_max)
            
            # Get Poynting flux using Eq. 55 in Saur+2013 (S_poynt) and Eq. 8 in Zarka
            # 2007 (S_poyn_ZL), in cgs units. They coincide, except for a factor 2. 
            S_poynt, S_poynt_ZL = spi.get_S_poynt(R_planet_eff, B_sw, v_alf, v_rel, M_A, alpha, geom_f)

            # Get fluxes at Earth, in cgs units for both Saur+ (Flux_r_S...) and
            # Zarka/Lanza (Flux_r_S_ZL...)
            Flux_r_S_min, Flux_r_S_max, Flux_r_S_ZL_min, Flux_r_S_ZL_max = spi.get_Flux(Omega_min, Omega_max, 
                                                          Delta_nu_cycl, d, S_poynt, S_poynt_ZL)

            """
            Moving parts of plotting outside the loop
            """
            # Find out the position of the planet in the distance array
            d_diff = np.abs((d_orb-r_orb)/R_star)
            loc_pl = np.where(d_diff == d_diff.min())

            """
            ###########################################################################
            ####                  PLOTTING                                         ####
            ###########################################################################
            
            ### Plot received flux density as a function of distance from the star
            ###
            #matplotlib.rc_file_defaults()
            #plt.style.use(['bmh', '/home/torres/Dropbox/python/styles/paper.mplstyle'])

            lw = 3

            # Kepler's third law, with d_orb_mark in units of R_star, 
            # so that period_mark is in days.
            #
            #period_mark = np.array([1, 10, 20, 40, 80, 100, 120, 140, 160,])
            period_mark = np.array([1, 10, 30, 60, 100, 200, 500, 1000, 2000])
            d_orb_mark = (period_mark/yr)**(2/3) * M_star_msun**(1/3) * (au/R_star)

            # Plot only Flux density vs. orbital distance, or also log(M_A) vs. d_orb 
            multiple = 1
            if multiple:
                plt.figure(figsize=(8,11))
                ax1 = plt.subplot2grid((3,1),(0,0),rowspan=1,colspan=1)
                ax2 = plt.subplot2grid((3,1),(1,0),rowspan=2,colspan=1)
            else:
                plt.figure(figsize=(8,7.5))
                ax2 = plt.subplot2grid((1,1),(0,0),rowspan=1,colspan=1)


            x   = d_orb/R_star # (distance array, in units of R_star)

            y_min = np.log10(Flux_r_S_min) # minimum flux (array), Saur/Turnpenney model
            y_max = np.log10(Flux_r_S_max) # maximum flux (array)
            y_min_ZL = np.log10(Flux_r_S_ZL_min) # minimum flux (array), Zarka/Lanza model
            y_max_ZL = np.log10(Flux_r_S_ZL_max) # maximum flux (array)

            ax1.plot(x, np.log10(M_A), color='k', lw=lw)
            #ax2.plot(x, y_min, lw=lw, color='orange', lw=lw, label="Saur/Turnpenney model")
            #ax2.plot(x, y_max, lw=lw, color='orange')
            
            # Fill with color for ZL and ST models
            #
            ax2.fill_between(x, y_min_ZL, y_max_ZL, color="blue", alpha=0.7, label="Zarka/Lanza model")
            ax2.fill_between(x, y_min, y_max,color="orange", alpha=0.7, label="Saur/Turnpenney model")


            #ax2.plot(dvec/R_star,np.log10(Flux_r_S), lw=lw, label="Saur/Turnpenney model")
            #ax2.plot(dvec/R_star,np.log10(Flux_r_S_ZL), lw=lw,label = "Zarka/Lanza model")
            #
            #ax2.fill_between([np.amin(dvec)/R_star,np.amax(dvec)/R_star], \
            #                 [np.log10(Poynt_min),np.log10(Poynt_min)], \
            #                 [np.log10(Poynt_max),np.log10(Poynt_max)],color="orange",alpha=0.4)

            ax11 = ax1.twiny()
            ax11.set_xticks(d_orb_mark)
            ax11.set_xticklabels(period_mark)
            ax11.tick_params(top=False,which='minor')
            ax1.tick_params(top=False,which='both')
            ax1.set_xticklabels([])
            #ax2.tick_params(labeltop=False, labelright=True)

            # Axis limits
            draw_all_xlim = True
            if draw_all_xlim:
                xmin = np.amin(d_orb)/R_star
                xmax = np.amax(d_orb)/R_star
            else:
                xmin = 2.5
                xmax = 25

            ax11.set_xlim([xmin, xmax])
            ax1.set_xlim([xmin, xmax])
            ax2.set_xlim([xmin, xmax])

            #ax1.set_xscale('log')
            #ax2.set_xscale('log')
            #VALORES ORIGINALES
            ymin_ax1=-3
            ymin=-3
            ymax=3
            ax1.set_ylim([ymin_ax1,0])
            ax2.set_ylim([ymin,ymax])
            #VALORES TEST
            #ax1.set_ylim([-50,20])
            #ax2.set_ylim([-50,80])

            # Draw vertical line at average position of planet
            ax1.axvline(x=r_orb/R_star, ls='--', color='k', lw=2)
            ax2.axvline(x=r_orb/R_star, ls='--', color='k', lw=2)
            
            #
            ax11.set_xlabel("Orbital period [days]")
            ax2.set_xlabel(r"Distance / Stellar radius")
            ax1.set_ylabel(r"${\rm log} (M_A)$")
            ax2.set_ylabel(r"${\rm log}$ (Flux density [mJy])")

            #Draw also Alfven radial Mach number
            draw_M_A_radial = 0
            if draw_M_A_radial:
                ax12 = ax1.twinx() #instantiate 2nd axes that shares the same x-axis
                ax12.set_ylim([-3,0])
                color = 'tab:blue'            
                ax12.set_ylabel('')
                ax12.plot(x, np.log10(M_A_radial), color=color, lw=lw)
                ax12.set_ylabel(r"${\rm log} (M_{A, \rm radial})$", color=color)
                ax12.tick_params(axis='y', labelcolor=color)
            
            # draw 3*rms upper limit
            draw_rms = 1
            if draw_rms:
                ax2.axhline(y=np.log10(3*rms), ls='-.', color='grey', lw=2)
                xpos = d_orb_max/6
                ypos = np.log10(4*rms)
                ax2.text(x=xpos,y=ypos,s=r"3$\times$RMS",fontsize='small')
            
            ax1.plot([xmin, xmax],[22,22],'k--') # This line needs to be modified to take into account different
            ax2.legend(loc=1)

            # draw a little Earth at the planet position for visualization purposes
            draw_earth = 1
            if draw_earth:
                paths = ['./pics/earth.png']
                x = [r_orb/R_star]
                y = [np.log10(3*rms)]
                for x0, y0, path in zip(x, y, paths):
                    ab_earth = AnnotationBbox(spi.getImage(path), (x0, y0), frameon=False)
                    ax2.add_artist(ab_earth)            
    
            #Print out relevant input and output parameters, including the expected flux received at Earth 
            # from the SPI at the position of the planet
            # To this end, first find out the position of the planet in the distance array
            d_diff = np.abs((d_orb-r_orb)/R_star)
            loc_pl = np.where(d_diff == d_diff.min())

            B_planet_loc = round(float(B_planet_arr[loc_pl]/(bfield_earth*Tesla2Gauss)), 2) # Planetary magnetic field, in units of Bfield_earth
            #B_planet_loc = int(B_planet_arr[loc_pl]) # Planetary magnetic field, in Gauss
            
            #ax2.text(x=60,y=1.8,s=r"$B_\ast = \,{\rm G}$")
            #ax2.text(x=60,y=1.4,s=r"$B_{\rm pl} = 1 \,{\rm G}$")
            #
            #ax2.text(x=12,y=-1.2,s=r"$B_{\rm cycl}$   = " + str(B_star) + " G ",fontsize='small')
            #ax2.text(x=12,y=-1.6,s=r"$B_{\rm planet}$ = " + str(B_planet_arr) + " G ", fontsize='small')
            #ax2.text(x=12,y=-2.0,s=r"$n_{\rm corona}$ = " + str(n_sw_base/1e7) + "x10$^7$ cm$^{-3}$ ", fontsize='small')
            #ax2.text(x=3,y=0.1+np.log10(3*rms),s=r"Requested 3$\sigma$", fontsize='x-small')
            xpos =xmax*0.8
            ypos_offset = (ymax-ymin)/8
            ypos = (ymax-ymin)/4 + ypos_offset
            d_ypos = (ymax-ymin)/12
            ax2.text(x=xpos,y=ypos,s=r"$B_\star$    = " + str(f'{B_star:.1f}') + " G ",fontsize='small')
            ax2.text(x=xpos,y=ypos-d_ypos,s=r"$B_{\rm planet}$ = " + str(f'{B_planet_loc:.1f}') + r"$B_{\rm Earth}$", fontsize='small')
            #ax2.text(x=xpos,y=ypos-2*d_ypos, s=r"$n_{\rm corona}$ = " + str(n_sw_base/1e7) + "x10$^7$ cm$^{-3}$ ", fontsize='small')
            ax2.text(x=xpos,y=ypos-2*d_ypos, s=r"$\dot{M}$ = " + str(f'{M_star_dot:.1f}') + "$\dot{M}_\odot$", fontsize='small')
        
            # Create OUTPUT folder if it doesn't exist
            FOLDER = 'OUTPUT/'+str(Exoplanet.replace(" ", "_"))
            if not(os.path.isdir(FOLDER)):
                os.system('mkdir OUTPUT/'+str(Exoplanet.replace(" ", "_")))
            else:
                print(FOLDER+' already exists.')
            common_string = str(B_star)+"G"+"-Bplanet"+str(B_planet_arr[loc_pl])+"G"+'-'+str(eps_min*100)+'-'+str(eps_max*100)+'percent'             
            if open_field:
                # ax1.text(x=0, y=1, s= Exoplanet + " - Open field")
                #outfile = Exoplanet + "-Open-Bstar"+ common_string
                outfile = str(Exoplanet.replace(" ", "_")) + '/' + str(Exoplanet.replace(" ", "_")) + "-Open-Bstar" + common_string 
            else:
                # ax1.text(x=0, y=1, s= Exoplanet + " - Closed field")
                #outfile = Exoplanet + "-Closed-Bstar"+ common_string
                outfile = str(Exoplanet.replace(" ", "_")) + '/' + str(Exoplanet.replace(" ", "_")) + "-Closed-Bstar" + common_string 
            # Variable to send output to files (plotout= True), or show them in
            # the terminal (plotout = False) 
            if plotout:
                plt.tight_layout()
                outfilePDF = os.path.join(outdir, outfile+".pdf")
                plt.savefig(outfilePDF)
                outfilePNG = os.path.join(outdir, outfile+".png")
                plt.savefig(outfilePNG)
                plt.close()
            else:
                plt.tight_layout()
                plt.show()

            ###########################################################
            ################### Send OUTPUT to external text file/s
            ###########################################################

            outfileTXT = os.path.join(outdir, outfile+'.txt')
            out_to_file = OutputWriter(outfileTXT)
            out_to_file.write_parameters(T_corona, M_star_dot, mu, d, R_star, M_star, P_rot_star, B_star, 
                Exoplanet, Rp, Mp, r_orb, P_orb, loc_pl, n_base_corona, nu_plasma_corona, gyrofreq,
                Flux_r_S_min, Flux_r_S_max, rho_sw_planet, n_sw_planet, v_sw_base, Flux_r_S_ZL_min,
                Flux_r_S_ZL_max, v_sw, v_rel, v_alf, M_A, B_sw, Rmp, R_planet_eff)

            # Print out the expected flux received at Earth from the SPI at the position of the planet

            # First, find out the position of the planet in the distance array
            d_diff = np.abs((d_orb-r_orb)/R_star)
            location_pl = np.where(d_diff == d_diff.min())

            #print('S_poynt: {0} erg/s\n'.format(S_poynt))
            #print('Flux_ST: {0} mJy\n'.format(Flux_r_S_max))    
            #print('Flux_ZL: {0} mJy\n'.format(Flux_r_S_ZL_max))    

            # Print out minimum and maximum values of flux density at the planet location
            #print("B_star =", B_star, " G; " "B_planet = ", B_planet_arr, " G ")
            #print("\nPrint out Poynting Flux at the planet location")
            #print("Saur/Turnpenney (erg/s): ", S_poynt[location_pl])
            #print("\nPrint out Poynting Flux at the first cell")
            #print("Saur/Turnpenney (erg/s): ", S_poynt[0])

            print("\nPrint out minimum and maximum values of flux density at the planet location")
            print('B_planet_loc = {0:.3f} G'.format(B_planet_loc * bfield_earth*Tesla2Gauss))
            #print('Rmp / Rp = {0:.3f} '.format(Rmp[loc_pl]/Rp))
            print("Saur/Turnpenney (mJy): ", Flux_r_S_min[location_pl], Flux_r_S_max[location_pl])
            print("Zarka/Lanza: (mJy)", Flux_r_S_ZL_min[location_pl], Flux_r_S_ZL_max[location_pl])
            print(f"Done with planet {Exoplanet}")
            #print("\nPrint out minimum and maximum values of flux density at the first cell")
            #print("Saur/Turnpenney (mJy): ", Flux_r_S_min[0], Flux_r_S_max[0])
            #print("Zarka/Lanza: (mJy)", Flux_r_S_ZL_min[0], Flux_r_S_ZL_max[0])
            """


            #### TEMPORARY TABLE
            ####################
            
            #output_table=pd.copy(data)
            #output=output[['planet_name','star_name','d_star(pc)','mass_star(m_sun)','radius_star(r_sun)',
            print(Bfield_geom_arr[ind],magnetized_pl_arr[ind1])
            if ((Bfield_geom_arr[ind]==0) and (magnetized_pl_arr[ind1]==True)):            
                planet_name_list.append(Exoplanet)
                star_name_list.append(starname)
                d_star_list.append("{:.2f}".format(d/pc))
                mass_star_list.append("{:.3f}".format(M_star/M_sun))
                radius_star_list.append("{:.3f}".format(R_star/R_sun))
                p_rot_list.append("{:.3f}".format(P_rot_star/day))
                bfield_star_list.append("{:.2f}".format(B_star))
                a_list.append("{0:.3f}".format(r_orb/au))
                p_orb_list.append("{0:.3f}".format(P_orb))
                eccentricity_list.append("{0:.3f}".format(eccentricity))
                q_list.append("{0:.3f}".format(q/au))
                Q_list.append("{0:.3f}".format(Q/au))    
                mass_planet_list.append("{:.2f}".format(Mp/M_earth))
                radius_planet_list.append("{:.2f}".format(Rp/R_earth))
                T_cor_list.append("{:.2e}".format(T_corona))
                m_dot_list.append("{:.2e}".format(M_star_dot))
                nu_pl_list.append("{:.2f}".format(nu_plasma_corona/1e6))
                nu_cycl_list.append("{:.2f}".format(gyrofreq/1e6))
                rho_pl_list.append("{:.2f}".format(rho_sw_planet[loc_pl][0]))
                B_pl_list.append("{:.2e}".format(B_planet[loc_pl][0]))
                B_sw_list.append("{:.2e}".format(B_sw[loc_pl][0]))
                v_alf_list.append("{:.2e}".format(v_alf[loc_pl][0]))
                M_A_list.append("{:.2e}".format(M_A[loc_pl][0]))
                Flux_r_S_ZL_min_list.append("{:.2e}".format(Flux_r_S_ZL_min[loc_pl][0]))          
                Flux_r_S_ZL_max_list.append("{:.2e}".format(Flux_r_S_ZL_max[loc_pl][0]))
                P_Bpl_list.append("{:.2e}".format(P_B_planet[loc_pl][0]))
                P_dyn_list.append("{:.2e}".format(P_dyn_sw[loc_pl][0]))
                P_th_list.append("{:.2e}".format(P_th_sw[loc_pl][0]))
                P_Bsw_list.append("{:.2e}".format(P_B_sw[loc_pl][0]))
                Rmp_list.append("{:.2e}".format(Rmp[loc_pl][0]/Rp))
            
# dictionary of lists 
parameters = {'planet_name': planet_name_list, 'star_name': star_name_list, 'd_star(pc)': d_star_list,
              'mass_star(m_sun)':mass_star_list, 'radius_star(r_sun)': radius_star_list, 
              'p_rot(days)': p_rot_list, 'bfield_star(gauss)': bfield_star_list, 'a(au)':a_list,
              'p_orb(days)':p_orb_list,'eccentricity':eccentricity_list,'q':q_list,'Q':Q_list,
              'mass_planet(m_earth)':mass_planet_list, 'radius_planet(r_earth)': radius_planet_list,
              'T_cor(K)':T_cor_list, 'mass_loss_rate(solar)':m_dot_list, 'nu_pl':nu_pl_list,
              'nu_cycl':nu_cycl_list,'rho_pl':rho_pl_list,'B_pl':B_pl_list,'B_sw':B_sw_list,'v_alf':v_alf_list,
              'M_A':M_A_list,'Flux_r_S_ZL_min':Flux_r_S_ZL_min_list, 'Flux_r_S_ZL_max':Flux_r_S_ZL_max_list,
              'P_Bpl':P_Bpl_list,'P_dyn':P_dyn_list,'P_th':P_th_list,'P_Bsw':P_Bsw_list,'Rmp':Rmp_list
}   
output = pd.DataFrame(parameters)

            
# Generate table with useful SPI parameters to generate various plots
output.to_csv('OUTPUT/out_table.csv')

# Generate plots from data in out_table.csv
os.system('python plot_out_table.py')
