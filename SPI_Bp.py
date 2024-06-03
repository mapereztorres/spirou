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
os.system('clear')
matplotlib.rc_file_defaults()
#plt.style.use(['bmh','spi.mplstyle'])
plt.style.use(['bmh','SPIworkflow/spi.mplstyle'])
### Getting parameters to predict SPI radio emission
#
# Observing parameters: Observin frequency, assumed rms noise, Delta_nu
# Stellar parameters: T_corona, n_base_corona, B_field, Isothermality of plasma
# Geometry of the sub-Alfvénic interaction: alpha, theta_M
#
from SPIworkflow.data_selection import *

from SPIworkflow.__init__ import *

# Import useful constants and functions to be used in the code
from SPIworkflow.constants import *
import SPIworkflow.SPIutils as spi
print('M_star_dot', M_star_dot)
print('10**M_star_dot :',10**M_star_dot)         
from SPIworkflow.data import get_spi_data, create_data_tables

# In the future, use function n_wind (under SPIworkflow)
# ## Instead of using particle density, we'll use M_dot 
# ### Currently, this function isn't used for the calculations
#
# n_w = n_wind(M_star_dot, r0=R_sun, v_r0=1.0)


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

print(source_data)
print(data)
# Assume fully ionized, purely hydrogen plasma (=> 50% protons, 50% electrons)
mu = (0.5*m_p + 0.5*m_e)/(m_p + m_e) # mean "molecular weight"
m_av = mu * m_p  # average mass density

# Isothermal Parker wind (assumed)
# Isothermal sound speed, in cm/s- Depends only on the Temperature of the stellar corona
vsound = np.sqrt(k_B * T_corona / m_av) 


# Setting the stellar magnetic field geometry and the value of the 
# intensity of the planetary magnetic field
# 
# Stellar magnetic field geometry
# The convention is that Bfield_geom_arr = 1 => open Parker spiral geometry; 
#                        Bfield_geom_arr = 0 - closed dipolar geometry
#Bfield_geom_arr = [0, 1]
Bfield_geom_arr = [0, 1]

# Magnetic field of the planet
#  Bp0_arr = 0 - unmagnetized planet, i.e., B_planet = 0 G; 
#  Bp0_arr = 1 - magnetized planet, i.e., B_planet > 0 G; 

# Magnetic field of the planet
#Bp0_arr= [0, 1]

### Select data in the array to run the code on
#
#star_array = range(len(data))
#star_array = [0, 1, 2]
#star_array = [21]
#star_array = np.linspace(0,22,23)
#star_array = [15]
print(star_array)
M_star_dot_set=np.copy(M_star_dot)
#M_star_dot = 1.4020249747030762 #10**(0.147)#0.1#10.23 # For GJ 436 (R_star=0.464 R_sun), this results in n_sw_base ~ 1e7 cm^-3 
print('M_star_dot :',M_star_dot_set)
#M_star_dot=np.copy(M_star_dot_set)
print('10**M_star_dot :',10**M_star_dot)
print('log(M_star_dot) :',np.log10(M_star_dot))
print(data.loc[[star_array[0]]])
#print(data['d_star(pc)','radius_star(r_sun)','mass_star(m_sun)','p_rot(days)','bfield_star(gauss)'])
for indi in star_array:
	print(data.loc[[star_array[0]]])
	#if indi==2:
	#	rms=0.022
	#if indi==18:
	#	rms=0.023	
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
	print('d :',d/pc)
	print('R_star :',R_star/R_sun)
	print('M_star :',M_star/M_sun)
	print('P_rot_star :',P_rot_star/day)
	print('B_star :', B_star)
	print('Exoplanet :',Exoplanet)
	print('Rp :', Rp/R_earth)
	print('Mp :',Mp/M_earth)
	print('r_orb :',r_orb/au)
	print('P_orb :',P_orb)
	for ind in Bfield_geom_arr:
	#for ind1 in Bp0_arr:
		npoints=40000
		Bmin=0
		Bmax=4
		Bp0_arr=np.linspace(Bmin, Bmax, npoints, endpoint=True)
		print(r_orb/R_star)
		# Common properties for star and planet
		# 
		Omega_star = 2.0*np.pi / P_rot_star # Angular rotation velocity of the star
		M_star_msun = M_star / M_sun # Stellar mass in units of solar mass

		r_sonic =  G * M_star / (2 * vsound**2) # Radius of sonic point

		# MPT: The two lines below ("Nsteps" and "d_orb" should be outside this loop, 
		# but for some reason to be understood, the first time results in v_sw[0] = 0, which in turn results
		# in a division by zero. This causes exiting the program. So, for now, we keep those lines inside the loop, which 
		# prevents the zeroth value. Maybe the line v_sw = np.zeros(len(d_orb)) is causing the trouble? Check it
		#Nsteps = 10
		#d_orb_max = 2*r_orb/R_star # Max. orbital distance, in units of R_star
		d_orb_max = r_orb/R_star  + 10 # Max. orbital distance, in units of R_star
		Nsteps = int(2*d_orb_max)
		#d_orb = np.linspace(1.002, 10, Nsteps) * R_star # Array of (orbital) distances to the star
		#d_orb = np.linspace(1.02, d_orb_max, Nsteps) * R_star # Array of (orbital) distances to the star, in cm 
		d_orb=[r_orb]
		#d_orb = np.linspace(1.02, 210, Nsteps) * R_star # Array of (orbital) distances to the star
		#print(len(d_orb))
		v_orb = (G * M_star/d_orb)**0.5 # Orbital (Keplerian) speed of planet as f(distance to star), in cm/s
		v_corot = r_orb * Omega_star # Corotation speed (cm/s)

		#Omega_planet = np.ones(len(d_orb)) * Omega_earth # array of angular speeds of the planet, in  s^(-1)
		Omega_planet =  v_orb / d_orb # Angular speed of the planet, in s^(-1). NOte that it's an array
		print('d_orb',d_orb)
		
		#
		# The wind speed is computed as in Turnpenney+18
		#
		# using the Lambert W function
		# D_r as in Eq. 18 of Turnpenney+2018, which is taken from eq. 15 in Cranmer 2004
		D_r = (d_orb/r_sonic)**(-4) * np.exp(4*(1 - r_sonic/d_orb) - 1)
		#print(D_r)
		v_sw2 = np.zeros(len(d_orb), dtype=complex)
		v_sw  = np.zeros(len(d_orb))
		
		for i in range(len(d_orb)):
			#print(d_orb,r_orb)
			if (d_orb[i]/r_sonic) >= 1.0:
				v_sw2[i] = -vsound**2 * lambertw(-D_r[i], k=-1)
			else: 
				v_sw2[i] = -vsound**2 * lambertw(-D_r[i], k=0)
			# The actual speeed is the real part of v_sw2[i]
			v_sw[i]  = np.sqrt(v_sw2[i].real)
		print("v_sw: ",v_sw)
		print('r_sonic: ',r_sonic)
		print(1.02* R_star/r_sonic)
		#print((r_orb/r_sonic)**(-4) * np.exp(4*(1 - r_sonic/r_orb) - 1))
		#print(-vsound**2 * lambertw(-(r_orb/r_sonic)**(-4) * np.exp(4*(1 - r_sonic/r_orb) - 1), k=-1))
		#print(-vsound**2 * lambertw(-D_r, k=-1))
		#print(np.sqrt(-vsound**2 * lambertw(-D_r, k=-1).real))
		#v_sw_base = v_sw[0]    # Stellar wind velocity at the closest distance to the star
		D_base = (1.02* R_star/r_sonic)**(-4) * np.exp(4*(1 - r_sonic/(1.02* R_star)) - 1)
		print('D_base: ',D_base)
		print('lambertw(-D_base, k=0) :',lambertw(-D_base, k=0))
		v_sw_base = np.sqrt( -vsound**2 * lambertw(-D_base, k=0) )    # Stellar wind velocity at the closest distance to the star
		print("v_sw_base: ",v_sw_base)
		#M_star_dot_all=[0.1,1,5,10,50]

		#M_star_dot_all=[0.147]#[0.1]
		M_star_dot_all=[M_star_dot]
		print("M_star_dot_all: ",M_star_dot_all)
		#for M_star_dot in M_star_dot_all:
		# Plasma number density at base of the corona
		n_base_corona=[]
		print('##############################################')
		print('before spi.n_wind')
		print('##############################################')
		print('M_star_dot_elem :',M_star_dot_all)
		print('R_star :',R_star)
		print('v_sw_base :',v_sw_base)
		print('mu :', mu)
		print('##############################################')
		for M_star_dot_elem in M_star_dot_all:
			basewind=spi.n_wind(M_star_dot_elem, R_star, v_sw_base, mu)
			n_base_corona.append(basewind)

			
			print('n_base_corona :',n_base_corona)#print(n_base_corona)
			# Maximum plasma frequency at the base of the corona. If the ECM
			# freq is less than the plasma frequency, the emission is
			# completely absorbed 
			nu_plasma_corona = spi.plasma_freq(n_base_corona) # in Hz
		print('n_base_corona :',n_base_corona)
		print('plasma_freq_corona :',9*np.sqrt(n_base_corona))
		
		print('##############################################')
		print('after spi.n_wind')
		print('##############################################')
		print('M_star_dot_elem :',M_star_dot_all)
		print('R_star :',R_star)
		print('v_sw_base :',v_sw_base)
		print('mu :', mu)
		print('##############################################')
		print(nu_plasma_corona)
		#print(nu_plasma_corona)
		#M_star_dot_all=M_star_dot
		M_star_dot_all=np.copy(M_star_dot)
		#print("V_sound = {0:.3f} km/s; V_sw at the base = {1:.3f} km/s".format(vsound/1e5, v_sw_base/1e5))    

		# Eq. 23 of Turnpenney+18 - Second term of RHS 
		# The vector v_rel = v_sw - v_orb (Eq. 5 in Saur+13, and see also Fig. 1 in Turnpenney+18)
		# 
		v_rel = np.sqrt(v_orb**2 + v_sw**2) # Relative speed between stellar wind and obstacle
		v_rel_angle = np.arctan(v_orb/v_sw) # Angle between radial vector and relative velocity

		# n_sw_planet - Number density of the wind at orbital distance to the planet. 
		# 
		# If the stellar plasma is assumed to be isothermal, then 
		# the density falls down as ~ R^(-2) * v_sw^(-1).
		#
		# Alternatively, we fix the density at the distance of the planet from the host star.
		#
		print('##############################################')
		print('before spi.n_wind')
		print('##############################################')
		print('M_star_dot_all :',M_star_dot_all)
		print('r_orb :',r_orb)
		print('v_sw :',v_sw)
		print('mu :', mu)
		print('##############################################')
		if isothermal:
		#n_sw_planet = n_sw_base / (d_orb/R_star)**2 / (v_sw/v_sw_base) # Plasma density at distance (R/R_star)
			n_sw_planet = spi.n_wind(M_star_dot_all, r_orb, v_sw, mu) # Plasma number density at distance (R/R_star)

		else:
			n_sw_planet = 1e4  # fixed                 
		print("n_sw_planet: ",n_sw_planet)
		print('plasma_freq_planet :',9*np.sqrt(n_sw_planet))
		print('##############################################')
		print('after spi.n_wind')
		print('##############################################')
		print('M_star_dot_all :',M_star_dot_all)
		print('r_orb :',r_orb)
		print('v_sw :',v_sw)
		print('mu :', mu)
		print('##############################################')
		# Magnetic field geometry
		# open_field - defines the geometry of the magnetic field
		open_field = Bfield_geom_arr[ind]
		print("M_star_dot_all: ",M_star_dot_all)
		if open_field: 
			# Open Parker Spiral - Falls down with distances as R^(-2) rather than R^(-3) as in the dipole case
			B_r = B_star * (r_orb/R_star)**(-2) # Stellar wind B-field at (R/R_star), Eqn 20 Turnpenney 2018
			B_phi = B_r * v_corot/v_sw # Azimuthal field (Eqn 21 Turnpenney 2018)
			B_sw = np.sqrt(B_r**2 + B_phi**2) # Total stellar wind B-field at planet orbital distance

			# Eq. 23 of Turnpenney 2018 -  First term of RHS
			B_ang = np.arctan(B_phi/B_r) # Angle the B-field makes with the radial direction

			# Angle between the stellar wind magnetic field and the impinging plasma velocity
			# Eqn 23 in Turnpenney 2018. It's also Eq. 13 in Zarka 2007
			theta = np.absolute(B_ang - v_rel_angle) 

			# theta is the angle between the B_sw (the insterstellar magnetic field), and the
			# incident stellar wind velocity.  See Fig. 1 in Turnpenney+2018
			#
			geom_f = (np.sin(theta))**2 # Geometric factor in efficiency 
		else:
			# Closed, dipolar configuration - It falls with distance as R^(-3)
			# B_star - magnetic field at the magnetic equator on the stellar surface
			# 
			phi = 0. # azimuth, measured from the North magnetic pole of the star (in degrees)
			phi *= np.pi/180. # to radians

			B_r   = -2 * B_star * (r_orb/R_star)**(-3) * np.cos(phi) # Radial component of the dipole magnetic field of the stellar wind as f(distance to star)
			B_phi = - B_star * (r_orb/R_star)**(-3) * np.sin(phi) # Azimuthal component of the dipole magnetic field 
			B_sw = np.sqrt(B_r**2 + B_phi**2) # Total dipole magnetic field 

			geom_f = 1.0 # Geometric factor. 1 for closed dipole configuration, different for the open field configuration

		# Alfven speed and Mach Number
		rho_sw = m_av * n_sw_planet #wind density, in g * cm^(-3)
		print("rho_sw: ",rho_sw)
		#v_alf = 2.18e11 * B_sw / np.sqrt(n_sw_planet) # Alfven speed at the distance of the planet, in cm/s
		#v_alf = B_sw / np.sqrt(4.0 * np.pi * rho_sw) 
		# Relativistically corrected Alfvén speed, as v_alf must be less than the speed of light
		v_alf = B_sw / np.sqrt(4.0 * np.pi * rho_sw) * 1./np.sqrt(1 + (B_sw**2/(4.*np.pi * rho_sw * c**2)))
		print('B_sw :',B_sw)
		print('B_r :',B_r)
		print('B_phi :',B_phi)
		print('v_corot :',v_corot)
		print(r_orb, Omega_star)
		print('P_rot_star :',P_rot_star)
		print('v_sw :',v_sw)
		print('np.sqrt(4.0 * np.pi * rho_sw) :', np.sqrt(4.0 * np.pi * rho_sw))
		print('1./np.sqrt(1 + (B_sw**2/(4.*np.pi * rho_sw * c**2))) :',1./np.sqrt(1 + (B_sw**2/(4.*np.pi * rho_sw * c**2))))
		print("v_alf: ",v_alf)
		M_A   = v_rel/v_alf # Alfven mach number
		M_A=M_A*np.ones(npoints)
		print("M_A: ",M_A)
		#Radial Alfvén speed
		mu_0_cgs = 1.0 # magnetic permeability in vacuum, in cgs units
		v_alf_r = B_r / np.sqrt(mu_0_cgs * rho_sw) # in cm/s
		M_A_radial = np.abs(v_sw / v_alf_r)
		
		#print('Relative speed = {0:.1e} km/s \n', format(v_rel/1e5))
		#print('Alfvén speed   = {0:.1e} km/s \n', format(v_alf/1e5))
		#print('Alfvén Mach number = {0:.3f} \n', format(M_A))
		#print('v_alf_rel/v_alf  = {0:.1e} km/s \n', format(v_alf_rel/v_alf))

		#print('Wind speed   = {0:.1e} km/s \n', format(v_sw/1e5))
		#print('Alfvén radial speed = {0:.1e} km/s \n', format(v_alf_r/1e5))
		#print('Alfvén radial Mach number = {0:.3f} \n', format(M_A_radial))

		# defines whether planet is unmagnetized (Bp0=0), or magnetized (Bp0 = 1)
		Bp=Bp0_arr
		'''
		for Bp0 in Bp0_arr:
			Bp.append(Bp0) 
			print(Bp0)
			if Bp0:
				# Planetary magnetic field, using Sano's (1993) scaling law, in units of B_earth 
				# This is a simple Sano(1993) scaling law dependence, assuming a tidally locked planet, 
				# core_radius equal to the Earth radius, and core density equal to that of the Earth.
				Bp = spi.bfield_sano(M_planet = Mp/M_earth, R_planet =
					Rp/R_earth, Omega_rot_planet =
					Omega_planet/Omega_earth) 

				Bp *= bfield_earth * Tesla2Gauss # in Gauss 
				Bp = np.ones(len(d_orb))  # For now, force Bp=1.0 Gauss
				# For now, force Bp=1.0 Gauss
			else:
				Bp = np.zeros(len(d_orb)) # unmagnetized planet
		'''
		print('Bp: ',Bp)
		print(len(Bp))
		#print('Bp: {0:.2e} Gauss'.format(Bp))

		#
		# Effective radius of the obstacle

		# Case 1. À la Saur+2013. 
		#
		# Effective radius of the Alfvén wing, in units of R_p (R_obst in Eq. 57 of Saur+2013, A&A)
		# It depends on the orientation, theta_M, of the intrinsic planetary magnetic field (Bp) 
		# wrt the external magnetic field (B_sw).
		#
		Rp_eff = Rp * np.sqrt(3*np.cos(theta_M/2)) * (Bp/B_sw)**(1./3.) # in cm
		print('Rp: ',Rp)
		print('Rp_eff: ',Rp_eff)
		Rp_eff[ Rp_eff < Rp] = Rp # Rp_eff cannot be smaller than Rplanet    
		print('Rp_eff :')
		print(Rp_eff)
		
		Rlarger=np.copy(Rp_eff)
		Rlarger[ Rlarger == Rp] = 0
		# Case 2. À la Zarka (2007), Turnpenney+2018, etc.
		#
		# Compute radius of magnetopause, Rmp as balance of wind and planet's
		# pressures

		# Planet pressure - only the magnetic component is considered
		P_Bp     = spi.get_P_Bp(Bp) 

		# Stellar wind pressure
		P_dyn_sw = spi.get_P_dyn_sw(n_sw_planet, mu, v_rel) 
		P_th_sw  = spi.get_P_th_sw(n_sw_planet, T_corona)
		P_B_sw   = spi.get_P_B_sw(B_sw)

		# Radius of magnetopause, in cm
		Rmp = spi.get_Rmp(P_Bp, P_dyn_sw, P_th_sw, P_B_sw) * Rp
		print('Rmp: ',Rmp)
		# The effective radius is the radius of the magnetopause
		Rp_eff = Rmp
		Rp_eff[ Rp_eff < Rp] = Rp # Rp_eff cannot be smaller than Rp
		print('Rp_eff: ',Rp_eff)

		# Total Poynting flux, as in Saur+2013 - Eq. 55 (page 7 of 20)
		# Applies if  M_A is small (<< 1)
		# Note that for the geometric factor, we follow Turnpenney's definition, so 
		# the factor is sin^2(theta), not cos^2(theta)
		# Saur says that the power is "per hemisphere", as Zarka below
		#
		# Total Poynting flux (S_mks), in mks units [kg * m * s^(-2) * A^(-2)]
		mu_0 = 4*np.pi*1e-7 # magnetic permeability in vacuum, in mks units
		# Poynting flux, in mks units
		S_poynt_mks = 2 * np.pi * (Rp_eff/1e2)**2 * (alpha*M_A)**2  \
			    * (v_alf/1e2) * (B_sw/1e4)**2 / mu_0 * geom_f
		S_poynt = S_poynt_mks * 1e7 # Total Poynting flux, in cgs units (erg/s) 
		print("S_poynt: ",S_poynt)
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
			    * (B_sw/1e4)**2 * geom_f / mu_0 * np.pi*(Rp_eff/1e2)**2 
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
		print("Flux_r_S_min: ",Flux_r_S_min)
		print("Flux_r_S_max: ",Flux_r_S_max)
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
			plt.figure(figsize=(11,11))
			#ax1 = plt.subplot2grid((3,1),(0,0),rowspan=1,colspan=1)
			ax2 = plt.subplot2grid((1,1),(0,0),rowspan=2,colspan=1)
		else:
			plt.figure(figsize=(8,7.5))
			ax2 = plt.subplot2grid((1,1),(0,0),rowspan=1,colspan=1)

		plt.yticks(fontsize=20)
		plt.xticks(fontsize=20)
		#x   = d_orb/R_star # (distance array, in units of R_star)
		Bp0_arr=np.linspace(Bmin, Bmax, npoints, endpoint=True)
		#x=np.log10(Bp0_arr)
		x=Bp0_arr
		print(x)
		y_min = np.log10(Flux_r_S_min) # minimum flux (array), Saur/Turnpenney model
		y_max = np.log10(Flux_r_S_max) # maximum flux (array)
		y_min_ZL = np.log10(Flux_r_S_ZL_min) # minimum flux (array), Zarka/Lanza model
		y_max_ZL = np.log10(Flux_r_S_ZL_max) # maximum flux (array)
		#y_min = Flux_r_S_min # minimum flux (array), Saur/Turnpenney model
		#y_max = Flux_r_S_max # maximum flux (array)
		#y_min_ZL = Flux_r_S_ZL_min # minimum flux (array), Zarka/Lanza model
		#y_max_ZL = Flux_r_S_ZL_max # maximum flux (array)
		print("x: ",x)
		print(np.log10(M_A))
		#ax1.plot(x, np.log10(M_A), color='k', lw=lw)
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

		#ax11 = ax1.twiny()
		#ax11.set_xticks(d_orb_mark)
		#ax11.set_xticklabels(period_mark)
		#ax11.tick_params(top=False,which='minor')
		#ax1.tick_params(top=False,which='both')
		#ax1.set_xticklabels([])
		#ax2.tick_params(labeltop=False, labelright=True)

		# Axis limits
		draw_all_xlim = True
		if draw_all_xlim:
			#xmin = np.amin(d_orb)/R_star
			#xmax = np.amax(d_orb)/R_star
			xmin = Bmin
			xmax= Bmax
		else:
			xmin = 2.5
			xmax = 25

		#ax11.set_xlim([xmin, xmax])
		#ax1.set_xlim([xmin, xmax])
		ax2.set_xlim([xmin, xmax])
		print('min :',y_min_ZL,y_min)
		print(min(min(y_min_ZL),min(y_min)))
		#print(max(y_max_ZL,y_max))
		ax2.set_ylim([min(min(y_min_ZL),min(y_min)), max(max(y_max_ZL),max(y_max))])
		#ax1.set_xscale('log')
		#ax2.set_xscale('log')
		#VALORES ORIGINALES
		#ymin_ax1=-3
		#ymin=-3
		#ymax=3
		#ymin=Flux_r_S_min
		'''
		ymin=np.log10(np.min(Flux_r_S_min))
		ymax=np.log10(np.max(Flux_r_S_max))
		print(ymin,ymax)
		#ax1.set_ylim([ymin_ax1,0])
		
		ax1.set_ylim([-2,0])
		ax2.set_ylim([ymin-2,ymax])
		'''
		#ymin_ax1=min(np.log10(M_A))#-3
		#ymin=-3
		#ymax=3
		#ymin=Flux_r_S_min
		print(int(np.min(y_min)*100)/100)
		print(int(np.min(y_min_ZL)*100)/100)
		min_array=[int(np.min(y_min)*100)/100,int(np.min(y_min_ZL)*100)/100]
		print(np.min(min_array))
		#print(np.min(int(np.min(y_min)*100)/100,int(np.min(y_min_ZL)*100)/100))
		ymin=np.min(min_array)

		print(int(np.min(y_max)*100)/100)
		print(int(np.min(y_max_ZL)*100)/100)
		max_array=[int(np.max(y_max)*100)/100,int(np.max(y_max_ZL)*100)/100]
		print(np.max(max_array))
		#print(np.min(int(np.min(y_min)*100)/100,int(np.min(y_min_ZL)*100)/100))
		ymax=np.max(max_array)			
		
		ymin=-3
		ymax=1
		#ymax=np.max(np.log10(np.max(Flux_r_S_max)),np.log10(np.min(y_max_ZL)))
		#ymin=ymin+2
		ymax=ymax+2
		ymin=-2.5
		ymax=0.5
		print(ymin,ymax)
		#ax1.set_ylim([ymin_ax1,0])#ax1.set_ylim([ymin_ax1,0])
		#ax2.set_ylim([ymin,ymax])		
		
		
		#VALORES TEST
		#ax1.set_ylim([-50,20])
		#ax2.set_ylim([-50,80])

		# Draw vertical line at average position of planet
		#ax1.axvline(x=r_orb/R_star, ls='--', color='k', lw=2)
		#ax2.axvline(x=r_orb/R_star, ls='--', color='k', lw=2)
		#ax2.axhline(y=np.log10(3*0.01), color='r', linestyle='--')
		ax2.axhline(y=np.log10(3*0.03), color='black', linestyle='--')
		#ax2.axhline(y=np.log10(0.12), color='black', linestyle='--')
		#ax2.axhline(y=np.log10(0.1), color='black', linestyle='--')
		xpos=(xmax*2)/3
		#ypos=np.log10(3*0.03)
		#ax2.text(x=xpos,y=ypos,s=r"3$\times 30 uJy",fontsize='medium',color='black')
		#ypos=np.log10(3*0.01)
		#ax2.text(x=xpos,y=ypos,s=r"3$\times 10 uJy",fontsize='medium', color='r')
		#
		#ax11.set_xlabel("Orbital period [days]")
		#ax2.set_xlabel(r"log(Planet magnetic field (G))")
		#ax2.set_xlabel(r"Planet magnetic field [G]",fontsize=30)
		#ax1.set_ylabel(r"${\rm log} (M_A)$")
		#ax2.set_ylabel(r"${\rm log}$ (Flux density [mJy])",fontsize=30)
		xpos=(xmax+xmin)/2
		
		ypos=np.log10(3*rms)+0.5-0.45#+0.15
		#ax2.text(x=xpos,y=ypos,s=r"$100 \mu$ Jy",fontsize='medium',color='black')
		#ax2.text(x=xpos,y=ypos,s=r"$3~\sigma$",fontsize='medium',color='black')
		#ax2.text(x=xpos,y=ypos,s="3 * max rms",fontsize='medium',color='black')
		#Draw also Alfven radial Mach number
		'''
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
		'''
		#ax1.plot([xmin, xmax],[22,22],'k--') # This line needs to be modified to take into account different
		#ax2.legend(loc=1)

		# draw a little Earth at the planet position for visualization purposes
		'''
		draw_earth = 0
		if draw_earth:
			paths = ['./pics/earth.png']
			x = [r_orb/R_star]
			y = [np.log10(3*rms)]
			for x0, y0, path in zip(x, y, paths):
			    ab_earth = AnnotationBbox(spi.getImage(path), (x0, y0), frameon=False)
			    ax2.add_artist(ab_earth)            
		'''
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
		#ax2.text(x=xpos,y=ypos,s=r"$B_\star$    = " + str(B_star) + " G ",fontsize='small')
		#ax2.text(x=xpos,y=ypos-d_ypos,s=r"$B_{\rm planet}$ = " + str(B_pl_loc) + r"$B_{\rm Earth}$", fontsize='small')
		#ax2.text(x=xpos,y=ypos-2*d_ypos, s=r"$n_{\rm corona}$ = " + str(n_sw_base/1e7) + "x10$^7$ cm$^{-3}$ ", fontsize='small')
		#ax2.text(x=xpos,y=ypos-2*d_ypos, s=r"$\dot{M}$ = " + str(M_star_dot) + "$M_\odot$", fontsize='small')
		#ax2.text(x=0.25,y=0.7,s='Saur/Turnpenny',color='orange',fontsize='x-large')
		#ax2.text(x=0.25,y=0.6,s='Zarka/Lanza',color='blue',fontsize='x-large')
		#ax2.text(x=(0.2*xmax+5*xmin)/2,y=(2.8*ymax+ymin)/2,s='Saur/Turnpenney',color='orange',fontsize='x-small')
		#ax2.text(x=(0.2*xmax+5*xmin)/2,y=(2.6*ymax+ymin)/2,s='Zarka/Lanza',color='blue',fontsize='x-small')
		ax2.set_facecolor("white")
		ax2.set_facecolor("white")
		#ax2.xaxis.set_tick_params(labeltop=True)
		ax2.xaxis.set_tick_params(labeltop=True,width=4,length=10,top=True)
		ax2.tick_params(axis='x', which='minor',top=True)
		ax2.tick_params(axis='x', labelsize=25,pad=10)
		#ax2.set_xticklabels(fontsize=20)
		
		ax2.yaxis.set_tick_params(width=4,length=10)
		ax2.tick_params(axis='y', labelsize=25,pad=10)
		#ax1.xaxis.set_tick_params(labeltop=True)
		
		ax2.set_xlabel(r"Planet magnetic field [G]",fontsize=30)
		ax2.set_ylabel(r"${\rm log}$ (Flux density [mJy])",fontsize=30)
		#ax3 = ax2.secondary_xaxis('top')	
		#ax2.set_ylim([-1,1])	
		#ax2.set_ylim([-4.5,0])
		ax2.set_ylim([-1.5,1.5])
		
		if str(Exoplanet.replace(" ", "_"))=='GJ_486_b':
			ax2.axvline(x=0.9, color='gray', linestyle='-')
			ax2.set_ylim([-2,1.2])
		print('mkdir OUTPUT/'+str(Exoplanet.replace(" ", "_")))
		try:
			os.system('mkdir OUTPUT/'+str(Exoplanet.replace(" ", "_")))
		except:
			pass
		common_string = str(B_star)+"G"+"-Bpl"+str(Bp[loc_pl])+"G"+' '+str(eps_min*100)+'-'+str(eps_max*100)+'perc'
		if open_field:
			#ax1.text(x=0, y=1, s= Exoplanet + " - Open field")
			outfile = str(Exoplanet.replace(" ", "_"))+'/'+ Exoplanet + "-Open-Bstar"+ common_string+'_B_pl_'+'T_corona'+"{:.2e}".format(T_corona)
		else:
			#ax1.text(x=0, y=1, s= Exoplanet + " - Closed field")
			outfile = str(Exoplanet.replace(" ", "_"))+'/'+ Exoplanet + "-Closed-Bstar"+ common_string+'_B_pl_'+'T_corona'+"{:.2e}".format(T_corona)

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
		M_star_dot_2=M_star_dot_set
		M_star_dot_2=10**M_star_dot_2
		print(M_star_dot_2)
		outfileTXT = os.path.join(outdir, outfile+'.txt')
		with open(outfileTXT, 'w') as f:
			f.write('# INPUT PARAMETERS:               ########\n')
			f.write('min {0:.3f}:\n'.format(min(min(y_min_ZL),min(y_min))))
			f.write('max {0:.3f}:\n'.format(max(max(y_max_ZL),max(y_max)))) 
			f.write('#                                 ########\n') 
			f.write('# GENERIC (fixed for all targets) ########\n') 
			f.write('#                                 ########\n') 
			f.write('T_corona                 = {0:.0e} K\n'.format(T_corona))
			f.write('Stellar wind mass loss rate = {0:.3f} Sun mass loss rate\n'.format(M_star_dot_2))
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
			#f.write('Semimajor axis of orbit  = {0:.3f} au\n'.format(r_orb/au))
			f.write('Semimajor axis of orbit  = {0:.3f} au = {1:.3f} R\u209B\u209C \n'.format(r_orb/au,r_orb/(R_star)))
			f.write('Orbital period           = {0:.3f} days\n'.format(P_orb))
			#f.write('Stellar wind particle density at the base = {0:.0e} cm-3\n'.format(n_sw_base))
			f.write('#                                 ########\n') 
			f.write('#                                 ########\n') 
			f.write('# OUTPUT PARAMETERS:              ########\n')
			f.write('#                                 ########\n') 
			n_base_corona=np.real(n_base_corona[0])
			nu_plasma_corona=np.real(nu_plasma_corona[0])
			print(n_base_corona, type(n_base_corona))
			print(nu_plasma_corona/1e6,type(nu_plasma_corona/1e6))
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
			#f.write('B_sw at r_orb: {0} \n'.format(B_sw[loc_pl]))
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


