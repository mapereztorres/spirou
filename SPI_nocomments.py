import os
import shutil
import numpy as np
import pandas as pd
import matplotlib
import matplotlib.pyplot as plt
from matplotlib.offsetbox import OffsetImage, AnnotationBbox
from scipy.special import lambertw

matplotlib.rc_file_defaults()

from SPIworkflow.__init__ import *

from SPIworkflow.constants import *
import SPIworkflow.SPIutils as spi
            
from SPIworkflow.data import get_spi_data, create_data_tables



outdir, df_planets, df_no_noplanets = create_data_tables()

if selection_criteria == False:
     data = get_spi_data(infile_data=source_data)
else:
     data = get_spi_data(infile_data=source_data,distance_max=15, p_orb_max = 10, bfield_min=100,bfield_max=1000.0, dec_min=-90)


mu = (0.5*m_p + 0.5*m_e)/(m_p + m_e) # mean "molecular weight"
m_av = mu * m_p  # average mass density

vsound = np.sqrt(k_B * T_corona / m_av) 


Bfield_geom_arr = [0, 1]


Bp0_arr= [0, 1]

star_array = range(len(data))
star_array = [8]

for indi in star_array:
    d      = data['d_star(pc)'][indi] * pc               # Distance to stellar system , in  cm
    R_star = data['radius_star(r_sun)'][indi] * R_sun    # Stellar radius in cm
    M_star = data['mass_star(m_sun)'][indi] * M_sun      # Stellar mass in g,
    P_rot_star = float(data['p_rot(days)'][indi]) * day  # Rotation period  of star, in sec
    B_star = data['bfield_star(gauss)'][indi]            # Stellar surface magnetic field
  
    # Planet - 
    Exoplanet = data['planet_name'][indi]
    Rp = data['radius_planet(r_earth)'][indi]*R_earth # Planetary radius, in cm
    Mp = float(data['mass_planet(m_earth)'][indi])*M_earth # Planetary mass, in grams
    r_orb  = data['a(au)'][indi]*au    # orbital distance, in cm
    P_orb = data['p_orb(days)'][indi] #orbital period of planet, in days
    
    for ind in Bfield_geom_arr:
        for ind1 in Bp0_arr:

            Omega_star = 2.0*np.pi / P_rot_star # Angular rotation velocity of the star
            M_star_msun = M_star / M_sun # Stellar mass in units of solar mass
        
            r_sonic =  G * M_star / (2 * vsound**2) # Radius of sonic point

            d_orb_max = r_orb/R_star  + 10 # Max. orbital distance, in units of R_star
            Nsteps = int(2*d_orb_max)
            d_orb = np.linspace(1.02, d_orb_max, Nsteps) * R_star # Array of (orbital) distances to the star, in cm 
            v_orb = (G * M_star/d_orb)**0.5 # Orbital (Keplerian) speed of planet as f(distance to star), in cm/s
            v_corot = d_orb * Omega_star # Corotation speed (cm/s)

            Omega_planet =  v_orb / d_orb # Angular speed of the planet, in s^(-1). NOte that it's an array

            D_r = (d_orb/r_sonic)**(-4) * np.exp(4*(1 - r_sonic/d_orb) - 1)
            v_sw2 = np.zeros(len(d_orb), dtype=complex)
            v_sw  = np.zeros(len(d_orb))

            for i in range(len(d_orb)):
                if (d_orb[i]/r_sonic) >= 1.0:
                    v_sw2[i] = -vsound**2 * lambertw(-D_r[i], k=-1)
                else: 
                    v_sw2[i] = -vsound**2 * lambertw(-D_r[i], k=0)
                v_sw[i]  = np.sqrt(v_sw2[i].real)

            v_sw_base = v_sw[0]    # Stellar wind velocity at the closest distance to the star
            
            n_base_corona = spi.n_wind(M_star_dot, R_star, v_sw_base, mu) 

            nu_plasma_corona = spi.plasma_freq(n_base_corona) # in Hz

            v_rel = np.sqrt(v_orb**2 + v_sw**2) # Relative speed between stellar wind and obstacle
            v_rel_angle = np.arctan(v_orb/v_sw) # Angle between radial vector and relative velocity
            
            if isothermal:
                n_sw_planet = spi.n_wind(M_star_dot, d_orb, v_sw, mu) # Plasma number density at distance (R/R_star)
            else:
                n_sw_planet = 1e4  # fixed                 
                
            open_field = Bfield_geom_arr[ind]
            
            if open_field: 
                B_r = B_star * (d_orb/R_star)**(-2) # Stellar wind B-field at (R/R_star), Eqn 20 Turnpenney 2018
                B_phi = B_r * v_corot/v_sw # Azimuthal field (Eqn 21 Turnpenney 2018)
                B_sw = np.sqrt(B_r**2 + B_phi**2) # Total stellar wind B-field at planet orbital distance
                B_ang = np.arctan(B_phi/B_r) # Angle the B-field makes with the radial direction
                theta = np.absolute(B_ang - v_rel_angle) 
                geom_f = (np.sin(theta))**2 # Geometric factor in efficiency 
            else:
                phi = 0. # azimuth, measured from the North magnetic pole of the star (in degrees)
                phi *= np.pi/180. # to radians

                B_r   = -2 * B_star * (d_orb/R_star)**(-3) * np.cos(phi) # Radial component of the dipole magnetic field of the stellar wind as f(distance to star)
                B_phi = - B_star * (d_orb/R_star)**(-3) * np.sin(phi) # Azimuthal component of the dipole magnetic field 
                B_sw = np.sqrt(B_r**2 + B_phi**2) # Total dipole magnetic field 

                geom_f = 1.0 # Geometric factor. 1 for closed dipole configuration, different for the open field configuration

            rho_sw = m_av * n_sw_planet #wind density, in g * cm^(-3)
            v_alf = B_sw / np.sqrt(4.0 * np.pi * rho_sw) * 1./np.sqrt(1 + (B_sw**2/(4.*np.pi * rho_sw * c**2)))
            M_A   = v_rel/v_alf # Alfven mach number
            
            mu_0_cgs = 1.0 # magnetic permeability in vacuum, in cgs units
            v_alf_r = B_r / np.sqrt(mu_0_cgs * rho_sw) # in cm/s
            M_A_radial = np.abs(v_sw / v_alf_r)

            Bp0 = Bp0_arr[ind1]
            print(Bp0)
            if Bp0:
                Bp = spi.bfield_sano(M_planet = Mp/M_earth, R_planet =
                        Rp/R_earth, Omega_rot_planet =
                        Omega_planet/Omega_earth) 
                
                Bp *= bfield_earth * Tesla2Gauss # in Gauss 
                Bp = np.ones(len(d_orb))  # For now, force Bp=1.0 Gauss
            else:
                Bp = np.zeros(len(d_orb)) # unmagnetized planet
            
            Rp_eff = Rp * np.sqrt(3*np.cos(theta_M/2)) * (Bp/B_sw)**(1./3.) # in cm
            Rp_eff[ Rp_eff < Rp] = Rp # Rp_eff cannot be smaller than Rplanet    

            P_Bp     = spi.get_P_Bp(Bp) 
            P_dyn_sw = spi.get_P_dyn_sw(n_sw_planet, mu, v_rel) 
            P_th_sw  = spi.get_P_th_sw(n_sw_planet, T_corona)
            P_B_sw   = spi.get_P_B_sw(B_sw)
            Rmp = spi.get_Rmp(P_Bp, P_dyn_sw, P_th_sw, P_B_sw) * Rp
            Rp_eff = Rmp
            Rp_eff[ Rp_eff < Rp] = Rp # Rp_eff cannot be smaller than Rp
            mu_0 = 4*np.pi*1e-7 # magnetic permeability in vacuum, in mks units
            # Poynting flux, in mks units
            S_poynt_mks = 2 * np.pi * (Rp_eff/1e2)**2 * (alpha*M_A)**2  \
                            * (v_alf/1e2) * (B_sw/1e4)**2 / mu_0 * geom_f
            S_poynt = S_poynt_mks * 1e7 # Total Poynting flux, in cgs units (erg/s) 
            S_poynt_ZL_mks = 1./ np.sqrt(1 + 1/M_A**2) *  (v_rel/1e2) \
                            * (B_sw/1e4)**2 * geom_f / mu_0 * np.pi*(Rp_eff/1e2)**2 
            S_poynt_ZL     = S_poynt_ZL_mks * 1e7  # in cgs units
            Ekin_min = 20 ; Ekin_max = 200
            beta_min = spi.beta_keV(Ekin_min) 
            beta_max = spi.beta_keV(Ekin_max)
            
            Omega_1 = 2*np.pi * (np.cos(np.arccos(beta_min) - beta_min/2) - np.cos(np.arccos(beta_min) + beta_min/2)) 
            Omega_2 = 2*np.pi * (np.cos(np.arccos(beta_max) - beta_max/2) - np.cos(np.arccos(beta_max) + beta_max/2)) 
            
            Omega_min = min(Omega_1, Omega_2)
            Omega_max = max(Omega_1, Omega_2)
            gyrofreq = e*B_star/(2*np.pi * m_e * c)
            Delta_nu_cycl = 0.5 * gyrofreq # width of ECMI emission = (0.5 * gyrofreq)
            bsa_Omega = 1.6 

            Omega_min = bsa_Omega
            Omega_max = bsa_Omega
            
            dilution_factor_min = eps_min / (Omega_max * d**2 * Delta_nu_cycl) 
            dilution_factor_max = eps_max / (Omega_min * d**2 * Delta_nu_cycl)

            Flux_r_S_min = S_poynt * dilution_factor_min
            Flux_r_S_min *= 1e26 # Flux density, in mJy
            Flux_r_S_max = S_poynt * dilution_factor_max
            Flux_r_S_max *= 1e26 # Flux density, in mJy

            Flux_r_S_ZL_min = S_poynt_ZL * dilution_factor_min
            Flux_r_S_ZL_min *= 1e26 # Flux density, in mJy
            Flux_r_S_ZL_max = S_poynt_ZL * dilution_factor_max
            Flux_r_S_ZL_max *= 1e26 # Flux density, in mJy
            
            lw = 3
            period_mark = np.array([1, 10, 30, 60, 100, 200, 500, 1000, 2000])
            d_orb_mark = (period_mark/yr)**(2/3) * M_star_msun**(1/3) * (au/R_star)
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

            B_pl_loc = round(float(Bp[loc_pl]/(bfield_earth*Tesla2Gauss)), 2) # Planetary magnetic field, in units of Bfield_earth
            #B_pl_loc = int(Bp[loc_pl]) # Planetary magnetic field, in Gauss
            
            #ax2.text(x=60,y=1.8,s=r"$B_\ast = \,{\rm G}$")
            #ax2.text(x=60,y=1.4,s=r"$B_{\rm pl} = 1 \,{\rm G}$")
            #
            #ax2.text(x=12,y=-1.2,s=r"$B_{\rm cycl}$   = " + str(B_star) + " G ",fontsize='small')
            #ax2.text(x=12,y=-1.6,s=r"$B_{\rm planet}$ = " + str(Bp) + " G ", fontsize='small')
            #ax2.text(x=12,y=-2.0,s=r"$n_{\rm corona}$ = " + str(n_sw_base/1e7) + "x10$^7$ cm$^{-3}$ ", fontsize='small')
            #ax2.text(x=3,y=0.1+np.log10(3*rms),s=r"Requested 3$\sigma$", fontsize='x-small')
            xpos =xmax*0.8
            ypos_offset = (ymax-ymin)/8
            ypos = (ymax-ymin)/4 + ypos_offset
            d_ypos = (ymax-ymin)/12
            ax2.text(x=xpos,y=ypos,s=r"$B_\star$    = " + str(B_star) + " G ",fontsize='small')
            ax2.text(x=xpos,y=ypos-d_ypos,s=r"$B_{\rm planet}$ = " + str(B_pl_loc) + r"$B_{\rm Earth}$", fontsize='small')
            #ax2.text(x=xpos,y=ypos-2*d_ypos, s=r"$n_{\rm corona}$ = " + str(n_sw_base/1e7) + "x10$^7$ cm$^{-3}$ ", fontsize='small')
            ax2.text(x=xpos,y=ypos-2*d_ypos, s=r"$\dot{M}$ = " + str(M_star_dot) + "$M_\odot$", fontsize='small')
        
             
            common_string = str(B_star)+"G"+"-Bplanet"+str(Bp[loc_pl])+"G"
            if open_field:
                ax1.text(x=0, y=1, s= Exoplanet + " - Open field")
                outfile = Exoplanet + "-Open-Bstar"+ common_string
            else:
                ax1.text(x=0, y=1, s= Exoplanet + " - Closed field")
                outfile = Exoplanet + "-Closed-Bstar"+ common_string
            
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

            
            outfileTXT = os.path.join(outdir, outfile+'.txt')
            with open(outfileTXT, 'w') as f:
                f.write('# INPUT PARAMETERS:               ########\n') 
                f.write('#                                 ########\n') 
                f.write('# GENERIC (fixed for all targets) ########\n') 
                f.write('#                                 ########\n') 
                f.write('T_corona                 = {0:.0e} K\n'.format(T_corona))
                f.write('Stellar wind mass loss rate = {0:.3f} Sun mass loss rate\n'.format(M_star_dot))
                f.write('mean molecular weight    = {0:.2f}\n'.format(mu))
                #f.write('mu_0 [mks units] = {0:.0e}\n \n'.format(mu_0))
                f.write('#                                 ########\n') 
                f.write('# STAR:                           ########\n') 
                f.write('#                                 ########\n') 
                #f.write('B_star = {0:.0f} G; B_planet = {1:.0f} G\n'.format(B_star, Bp[loc_pl]))
                f.write('Distance to star         = {0:.2f} pc\n'.format(d/pc))
                f.write('Stellar radius           = {0:.3f} Rsun\n'.format(R_star/R_sun))
                f.write('Stellar mass             = {0:.3f} Msun\n'.format(M_star/M_sun))
                f.write('Stellar rotation period  = {0:.3f} days\n'.format(P_rot_star/day))
                f.write('Stellar magnetic field   = {0:.2f} Gauss\n'.format(B_star))
                f.write('#                                 ########\n') 
                f.write('# PLANET:                         ########\n') 
                f.write('#                                 ########\n') 
                f.write('Exoplanet name           = {0}\n'.format(Exoplanet)) 
                f.write('Planetary radius         = {0:.3f} R_Earth\n'.format(Rp/R_earth))
                f.write('Planetary mass           = {0:.3f} M_Earth\n'.format(Mp/M_earth))
                f.write('Semimajor axis of orbit  = {0:.3f} au\n'.format(r_orb/au))
                f.write('Orbital period           = {0:.3f} days\n'.format(P_orb))
                #f.write('Stellar wind particle density at the base = {0:.0e} cm-3\n'.format(n_sw_base))
                f.write('#                                 ########\n') 
                f.write('#                                 ########\n') 
                f.write('# OUTPUT PARAMETERS:              ########\n')
                f.write('#                                 ########\n') 
                f.write('n_base_corona = {0:.3e} cm^-3\n'.format(n_base_corona))
                f.write('nu_plasma_corona = {0:.2e} MHz\n'.format(nu_plasma_corona/1e6))
                f.write('ECMI freq (fundamental) = {0:.0f} MHz\n'.format(gyrofreq/1e6))
                f.write('Flux_ST: ({0}, {1}) mJy\n'.format(Flux_r_S_min[loc_pl], Flux_r_S_max[loc_pl]))
                #f.write('Flux_ZL: ({0}, {1}) mJy\n'.format(Flux_r_S_ZL_min[loc_pl], Flux_r_S_ZL_max[loc_pl]))
                f.write('rho_sw at r_orb: {0} \n'.format(rho_sw[loc_pl]))
                f.write('n_sw_planet at r_orb: {0} \n'.format(n_sw_planet[loc_pl]))
                f.write('n_sw_planet at base: {0} \n'.format(n_sw_planet[0]))
                f.write('v_sw at base of the wind: {0} \n'.format(v_sw_base))
                f.write('Flux_ZL: ({0}, {1}) mJy\n'.format(Flux_r_S_ZL_min[loc_pl], Flux_r_S_ZL_max[loc_pl]))
                f.write('v_sw at the base of the wind: {0} \n'.format(v_sw_base))
                f.write('v_sw at r_orb: {0} \n'.format(v_sw[loc_pl]))
                f.write('v_sw(r_orb)/v_sw_base: {0} \n'.format(v_sw[loc_pl]/v_sw_base))
                f.write('v_rel at r_orb: {0} \n'.format(v_rel[loc_pl]))
                f.write('v_alf at r_orb: {0} \n'.format(v_alf[loc_pl]))
                f.write('M_A at r_orb: {0} \n'.format(M_A[loc_pl]))
                #f.write('v_rel/v_alf at r_orb: {0} \n'.format(v_rel[loc_pl]/v_alf[loc_pl]))
                f.write('B_sw at r_orb: {0} \n'.format(B_sw[loc_pl]))
                f.write('Rmp at r_orb: {0} \n'.format(Rmp[loc_pl]/Rp))
                f.write('Rp_eff at r_orb: {0} \n'.format(Rp_eff[loc_pl]/Rp))

            

            
            # Print out the expected flux received at Earth from the SPI at the position of the planet

            # First, find out the position of the planet in the distance array
            d_diff = np.abs((d_orb-r_orb)/R_star)
            location_pl = np.where(d_diff == d_diff.min())

            #print('S_poynt: {0} erg/s\n'.format(S_poynt))
            #print('Flux_ST: {0} mJy\n'.format(Flux_r_S_max))    
            #print('Flux_ZL: {0} mJy\n'.format(Flux_r_S_ZL_max))    

            # Print out minimum and maximum values of flux density at the planet location
            #print("B_star =", B_star, " G; " "B_planet = ", Bp, " G ")
            #print("\nPrint out Poynting Flux at the planet location")
            #print("Saur/Turnpenney (erg/s): ", S_poynt[location_pl])
            #print("\nPrint out Poynting Flux at the first cell")
            #print("Saur/Turnpenney (erg/s): ", S_poynt[0])

            print("\nPrint out minimum and maximum values of flux density at the planet location")
            print('B_planet = {0:.3f} G'.format(B_pl_loc * bfield_earth*Tesla2Gauss))
            #print('Rmp / Rp = {0:.3f} '.format(Rmp[loc_pl]/Rp))
            print("Saur/Turnpenney (mJy): ", Flux_r_S_min[location_pl], Flux_r_S_max[location_pl])
            print("Zarka/Lanza: (mJy)", Flux_r_S_ZL_min[location_pl], Flux_r_S_ZL_max[location_pl])
            print(f"Done with planet {Exoplanet}")
            #print("\nPrint out minimum and maximum values of flux density at the first cell")
            #print("Saur/Turnpenney (mJy): ", Flux_r_S_min[0], Flux_r_S_max[0])
            #print("Zarka/Lanza: (mJy)", Flux_r_S_ZL_min[0], Flux_r_S_ZL_max[0])


