# The expressions used here are from the following links:
#https://scholarlypublications.universiteitleiden.nl/access/item%3A3485847/view
#https://ar5iv.labs.arxiv.org/html/0803.4507
from SPIworkflow.__init__ import *
import scipy.integrate as integral

import SPIworkflow.SPIutils as spi
import numpy as np

def get_gaunt(T,nu):
    """
    OUTPUT - g:  float - Gaunt factor
    INPUT  - T:  float - Temperature of the wind,  in Kelvin
           - nu: float - Observing frequency, in Hz
    The Gaunt factor  is computed using Eq. X in fulanito et al. YYYY. 
    It also takes into account the ionization state (Z - read from __init__.py)
    """
    g = 10.6 + 1.9 * np.log10(T) - 1.26 * np.log10(Z * nu)
    return g

def ff_absorption(M_star, nu, T, m_av, X_p, mdot, R_ff_in, R_ff_out, NSTEPS_FF,R_star):
    '''
    mdot: stellar mass-loss rate, in units of the Sun mass-loss rate
    nu: frequency of observation, in Hz
    T: temperature of the corona, which would be constant for an isothermal wind, in K
    m_av - average particle mass of the stellar wind (in gr)    
    X_p  - Fraction of protons
    R_ff_in: Distance from stellar center where SPI emission takes place, in cm (>= 1.0)
    R_ff_out: Distance from stellar center where free-free absorption becomes negligible, in cm
    '''
    
    dist_absorption = np.linspace(R_ff_in, R_ff_out, NSTEPS_FF)
    #dist_absorption = np.logspace(np.log10(R_ff_in),np.log10(R_ff_out),NSTEPS_FF)
    v_sound, r_sonic, v_sw = spi.v_stellar_wind(dist_absorption, M_star, T, m_av)
    n_sw = spi.n_wind(mdot, dist_absorption, v_sw, m_av) 
    n_p  = n_sw * X_p
    n_e  = n_sw * (1 - X_p)
    g = get_gaunt(T, nu)
    knu=3.692*10**8*(1-np.exp((-h*nu)/(k_B*T)))*Z**2*g*T**(-1/2)*nu**(-3)
    alphanu = knu * n_e * n_p
    taunu=integral.trapezoid(alphanu,dist_absorption)
    return v_sw,n_sw, knu,alphanu,taunu
    

'''  
nu_ecm = B_star * 2.8e6 # cyclotron freq, in Hz

# Generate an array of distances ranging from where the SPI emission takes place above the stellar surface to where the free-free absorption would become irrelevant
Nsteps=10000
altitude = 1.5 #altitude over stellar surface where SPI takes place, in R/R_star units
outerlimit = r_orb*500 #limit for integration of free-free absorption, in cm
#*1.495978707*10**13
dist_absorption = np.linspace(data['radius_star(r_sun)'][indi]*6.957*10**10*altitude, outerlimit, Nsteps)

X_e, mu, m_av = spi.wind_composition(X_p = 0.5)
v_sound, r_sonic, v_sw = spi.v_stellar_wind(dist_absorption, M_star, T_corona, m_av)
# Plasma number density at distance (R/R_star)
g = get_gaunt(T_corona,nu_ecm)

#for mdot in M_star_dot_arr:
n_sw = spi.n_wind(mdot, dist_absorption, v_sw, m_av) 
nu_plasma = spi.plasma_freq(np.max(n_sw))
knu,alphanu,taunu = get_taunu(nu_ecm,T_corona,g,dist_absorption,n_sw)
#absorption factor
absor = np.exp(-taunu)
'''


