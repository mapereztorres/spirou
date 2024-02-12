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
            
from SPIworkflow.load_data import get_spi_data, create_data_tables, load_target


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
#Bfield_geom_arr = [0, 1]
Bfield_geom_arr = [0]

# B_planet_arr is like a False/True array, magfield_planet is the modulus of the magnetic field
B_planet_arr= [0]
Bfield_pl = 0.5

### Select data in the array to run the code on
print(source_data)
print(data)
#star_array = range(len(data))
#star_array = [0, 1, 2]
star_array = [0, 1]

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
data['radius_planet(r_earth)'].replace('', np.nan, inplace=True)
data.reset_index(inplace=True) # to prevent funny jumps in the indices

for indi in star_array:
    d, R_star, M_star, P_rot_star, B_star, Exoplanet, Mp, Rp, r_orb, P_orb = load_target(data, indi)

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


    # Common properties for star and planet
    # 
    M_star_msun = M_star / M_sun # Stellar mass in units of solar mass
    Omega_star = 2.0*np.pi / P_rot_star # Angular rotation velocity of the star

    d_orb_max = r_orb/R_star  + 10 # Max. orbital distance, in units of R_star
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

    for ind in Bfield_geom_arr:
        for ind1 in B_planet_arr:
            # Magnetic field geometry
            # open_field - defines the geometry of the magnetic field
            open_field = Bfield_geom_arr[ind]
            B_r, B_phi, B_sw, B_ang, theta, geom_f = spi.get_bfield_comps(open_field, B_star, d_orb, R_star, v_corot, v_sw, v_rel_angle)
            
            # Compute Alfvén parameters in the stellar wind at a distance d_orb 
            v_alf, M_A, v_alf_r, M_A_radial = spi.get_alfven(rho_sw_planet, B_sw, B_r, v_rel, v_sw)

            # defines whether planet is unmagnetized (B_planet_arr[ind1] = 0), or magnetized (B_planet_arr[ind1] = 1)
            if B_planet_arr[ind1]: # magnetized planet
                # Planetary magnetic field, using Sano's (1993) scaling law, in units of B_earth 
                # This is a simple Sano(1993) scaling law dependence, assuming a tidally locked planet, 
                # core_radius equal to the Earth radius, and core density equal to that of the Earth.
                B_planet = spi.bfield_sano(M_planet = Mp/M_earth, R_planet =
                            Rp/R_earth, Omega_rot_planet =
                            Omega_planet/Omega_earth) 
                
                B_planet *= bfield_earth * Tesla2Gauss # in Gauss 
                # For now, force B_planet = Bfield_pl
                B_planet  = np.ones(len(d_orb)) * Bfield_pl 
            else:  # unmagnetized planet
                B_planet  = np.zeros(len(d_orb)) # unmagnetized planet
            
            #
            # Effective radius of the obstacle

            # Case 1. À la Saur+2013. 
            #
            # Effective radius of the Alfvén wing, in units of R_p (R_obst in Eq. 57 of Saur+2013, A&A)
            # It depends on the orientation, theta_M, of the intrinsic planetary
            # magnetic field (B_planet) wrt the external magnetic field (B_sw).
            #
            R_planet_eff = Rp * np.sqrt(3*np.cos(theta_M/2)) * (B_planet/B_sw)**(1./3.) # in cm
            R_planet_eff[ R_planet_eff < Rp] = Rp # R_planet_eff cannot be smaller than Rplanet    

            # Case 2. À la Zarka (2007), Turnpenney+2018, etc.
            #
            # Compute radius of magnetopause, Rmp as balance of wind and planet's
            # pressures
            
            # Planet pressure - only the magnetic component is considered
            P_B_planet  = spi.get_P_B_planet(B_planet) 

            # Stellar wind pressure
            P_dyn_sw = spi.get_P_dyn_sw(n_sw_planet, mu, v_rel) 
            P_th_sw  = spi.get_P_th_sw(n_sw_planet, T_corona)
            P_B_sw   = spi.get_P_B_sw(B_sw)

            # Radius of magnetopause, in cm
            Rmp = spi.get_Rmp(P_B_planet, P_dyn_sw, P_th_sw, P_B_sw) * Rp

            # The effective radius is the radius of the magnetopause
            R_planet_eff = Rmp
            R_planet_eff[ R_planet_eff < Rp] = Rp # R_planet_eff cannot be smaller than Rp


            # Total Poynting flux, as in Saur+2013 - Eq. 55 (page 7 of 20)
            # Applies if  M_A is small (<< 1)
            # Note that for the geometric factor, we follow Turnpenney's definition, so 
            # the factor is sin^2(theta), not cos^2(theta)
            # Saur says that the power is "per hemisphere", as Zarka below
            #
            # Total Poynting flux (S_mks), in mks units [kg * m * s^(-2) * A^(-2)]
            mu_0 = 4*np.pi*1e-7 # magnetic permeability in vacuum, in mks units
            # Poynting flux, in mks units
            S_poynt_mks = 2 * np.pi * (R_planet_eff/1e2)**2 * (alpha*M_A)**2  \
                            * (v_alf/1e2) * (B_sw/1e4)**2 / mu_0 * geom_f
            S_poynt = S_poynt_mks * 1e7 # Total Poynting flux, in cgs units (erg/s) 
            
            # Total Poynting flux, as in Lanza 2009 (Eq. 8) and Zarka 2007 (Eq. 9) 
            # They have a value which is different by a factor 2 * M_A * alpha^2
            # In addition, they include a geometric factor of order 1/2.
            #
            #ZL_factor = 0.5
            #S_poynt_ZL = S_poynt * ZL_factor / (2 * M_A * alpha**2 * geom_f)

            # Total Poynting flux, as in Zarka 2007 (Eq. 8), but using the
            # Alfvén conductance as defined in Neubaur.
            # Eq. 8 in Zarka2007 explicityly states that it is the power "per hemisphere",
            # In this sense, this is the same as in the expresion by Saur, 
            # so there seems to be a factor of two discrepancy, 
            # if taken into account that v_rel = v_alf * M_A. 
            #
            S_poynt_ZL_mks = 1./ np.sqrt(1 + 1/M_A**2) *  (v_rel/1e2) \
                            * (B_sw/1e4)**2 * geom_f / mu_0 * np.pi*(R_planet_eff/1e2)**2 
            S_poynt_ZL     = S_poynt_ZL_mks * 1e7  # in cgs units
            
            # Beam solid angle covered by the ECM emission
            # It depends on the kinetic energy of electrons (as beta is determined from them), in keV
            #
            #Ekin_min = 10 ; Ekin_max = 511        
            Ekin_min = 20 ; Ekin_max = 200
            beta_min = spi.beta_keV(Ekin_min) 
            beta_max = spi.beta_keV(Ekin_max)
            
            # BEAM SOLID ANGLE
            # We assume an emission cone with half-opening angle theta and angular width d_theta.
            # theta and d_theta are related to the speed of the electrons as theta \approx d_theta \approx v/c = beta
            # 
            # The solid angle, Omega, of a cone with half-opening angle, theta, 
            # is Omega = 2*np.pi * (1. - np.cos(theta))
            # Therefore, a cone sheet with inner half-opening angle, theta, and angular width d_theta
            # has a solid angle Omega = 2*np.pi* ( np.cos (theta - d_theta/2) - np.cos(theta +d_theta/2)) 
            
            #Range of beam solid angles (Omega) for the emitter depends on beta_min and beta_max
            Omega_1 = 2*np.pi * (np.cos(np.arccos(beta_min) - beta_min/2) - np.cos(np.arccos(beta_min) + beta_min/2)) 
            Omega_2 = 2*np.pi * (np.cos(np.arccos(beta_max) - beta_max/2) - np.cos(np.arccos(beta_max) + beta_max/2)) 
            
            # beam solid angle of the emitter / 4*pi 
            # Check why I (MPT) wrote " /4*pi" in the previous line 
            Omega_min = min(Omega_1, Omega_2)
            Omega_max = max(Omega_1, Omega_2)

            
            #print("beta_min of electrons (E_k={0:.1f} keV) = {1:.3f}; \n beta_max of electrons (E_k={2:.1f} keV) = {3:.3f}".\
            #      format(Ekin_min, beta_min, Ekin_max, beta_max))
            #print("Half-opening angle of emission cone: theta_min = {0:.1f} deg; theta_max = {1:.1f} deg".\
            #       format(np.arccos(beta_max)*180./np.pi, np.arccos(beta_min)*180/np.pi))
            #print("(Angular) width of emission cone: theta_min = {0:.1f} deg; theta_max = {1:.1f} deg".\
            #       format(beta_min*180./np.pi, beta_max*180/np.pi))
            #print("solid angles: Min = {0:.3f}; Max = {1:.3f}".format(Omega_min, Omega_max))
            #print(" ")

            # in-band radio power received from one whole hemisphere
            #power  = 2*np.pi * flux * mJy * d**2 * Delta_nu_obs 

            # The range of allowed powers, considering the beamed solid angle
            # and the possible total bandwidth
            # 
            gyrofreq = e*B_star/(2*np.pi * m_e * c)
            Delta_nu_cycl = 0.5 * gyrofreq # width of ECMI emission = (0.5 * gyrofreq)
        
            #Fix flux_min = flux_max = flux
            #flux_min = flux_max = flux 

            #
            #power_min = power/(2*np.pi) * (flux_min/flux) * Omega_min * Delta_nu_cycl/Delta_nu_obs
            #power_max = power/(2*np.pi) * (flux_max/flux) * Omega_max * Delta_nu_cycl/Delta_nu_obs
            #
            # Range of values for the star-ward Poynting flux
            #Poynt_min = power_min / eps_max 
            #Poynt_max = power_max / eps_min 

            #print("In-band power (from 1 whole hemisphere) =  {0:.2e}".format(power))
            #print("power (from 1 sr)  =  {0:.2e}".format(power/2/np.pi * Delta_nu_cycl/Delta_nu_obs))
            #print("power_min =  {0:.2e}; power_max = {1:.2e}".format(power_min, power_max))
            #print("Poynting fluxes: Min = {0:.2e}; Max = {1:.2e}".format(Poynt_min, Poynt_max))
            # Flux density received at Earth (from the theoretically expected Poynting flux)
            # 
            
            #beam solid angle of the ECMI emission, in sterradians. There is no need to fix this value, as this is
            #given by the population of electrons. If we don't fix it, we can constrain other parameters, e.g., 
            # eps, the efficienty factor in converting energy into Poynting flux.
            bsa_Omega = 1.6 

            # For simplicity, take Omega_min = Omega_max = Omega
            #Omega_min = Omega_max = Omega
            Omega_min = bsa_Omega
            Omega_max = bsa_Omega
            
            dilution_factor_min = eps_min / (Omega_max * d**2 * Delta_nu_cycl) 
            dilution_factor_max = eps_max / (Omega_min * d**2 * Delta_nu_cycl)

            # Min and Max expected flux density to be received for Saur-Turnpenney model, in erg/s/Hz/cm2
            Flux_r_S_min = S_poynt * dilution_factor_min
            Flux_r_S_min *= 1e26 # Flux density, in mJy
            Flux_r_S_max = S_poynt * dilution_factor_max
            Flux_r_S_max *= 1e26 # Flux density, in mJy

            # Min and Max expècted flux density to be received for Zarka-Lanza model, in erg/s/Hz/cm2
            Flux_r_S_ZL_min = S_poynt_ZL * dilution_factor_min
            Flux_r_S_ZL_min *= 1e26 # Flux density, in mJy
            Flux_r_S_ZL_max = S_poynt_ZL * dilution_factor_max
            Flux_r_S_ZL_max *= 1e26 # Flux density, in mJy
            

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

            #ax1.set_xscale('log')
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

            B_planet_loc = round(float(B_planet[loc_pl]/(bfield_earth*Tesla2Gauss)), 2) # Planetary magnetic field, in units of Bfield_earth
            #B_planet_loc = int(B_planet[loc_pl]) # Planetary magnetic field, in Gauss
            
            #ax2.text(x=60,y=1.8,s=r"$B_\ast = \,{\rm G}$")
            #ax2.text(x=60,y=1.4,s=r"$B_{\rm pl} = 1 \,{\rm G}$")
            #
            #ax2.text(x=12,y=-1.2,s=r"$B_{\rm cycl}$   = " + str(B_star) + " G ",fontsize='small')
            #ax2.text(x=12,y=-1.6,s=r"$B_{\rm planet}$ = " + str(B_planet) + " G ", fontsize='small')
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
        
            # save all plots in a specific folder for each planet  
            print('mkdir OUTPUT/'+str(Exoplanet.replace(" ", "_")))
            try:
                os.system('mkdir OUTPUT/'+str(Exoplanet.replace(" ", "_")))
            except:
                pass 
            common_string = str(B_star)+"G"+"-Bplanet"+str(B_planet[loc_pl])+"G"+'-'+str(eps_min*100)+'-'+str(eps_max*100)+'percent'             
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

            # Send output to external file/s

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
            #print("B_star =", B_star, " G; " "B_planet = ", B_planet, " G ")
            #print("\nPrint out Poynting Flux at the planet location")
            #print("Saur/Turnpenney (erg/s): ", S_poynt[location_pl])
            #print("\nPrint out Poynting Flux at the first cell")
            #print("Saur/Turnpenney (erg/s): ", S_poynt[0])

            print("\nPrint out minimum and maximum values of flux density at the planet location")
            print('B_planet = {0:.3f} G'.format(B_planet_loc * bfield_earth*Tesla2Gauss))
            #print('Rmp / Rp = {0:.3f} '.format(Rmp[loc_pl]/Rp))
            print("Saur/Turnpenney (mJy): ", Flux_r_S_min[location_pl], Flux_r_S_max[location_pl])
            print("Zarka/Lanza: (mJy)", Flux_r_S_ZL_min[location_pl], Flux_r_S_ZL_max[location_pl])
            print(f"Done with planet {Exoplanet}")
            #print("\nPrint out minimum and maximum values of flux density at the first cell")
            #print("Saur/Turnpenney (mJy): ", Flux_r_S_min[0], Flux_r_S_max[0])
            #print("Zarka/Lanza: (mJy)", Flux_r_S_ZL_min[0], Flux_r_S_ZL_max[0])


