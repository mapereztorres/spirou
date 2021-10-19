import numpy as np
import pandas as pd
import os
import shutil

import matplotlib
matplotlib.rc_file_defaults()

import matplotlib.pyplot as plt
plt.style.use(['bmh','/home/torres/Dropbox/python/styles/paper.mplstyle'])

# Import useful constants to be used in the code
from SPIworkflow.constants import *

### Setting up the relevant parameters to predict SPI radio emission
from SPIworkflow.__init__ import *

### READ IN THE DATA 
data = pd.read_csv("./INPUT/SPI-sources.csv")
#print(data[89:90])

# Sound speed, in cm/s - Depends only on the Temperature of the stellar corona
vsound = np.sqrt(5/3 * k_B * T_corona/m_h) 
#print("vsound = {0:.2f} km/s".format(vsound/1e5))


#####################################
### START OF FUNCTION DEFINITIONS ###
#####################################

def Kepler_r(M_star, P_orb):
    """Computes the orbital radius of a planet orbiting a star, given the orbital
    period.

    INPUT:  M_star - Star mass, in units of M_sol
            P_orb  - Orbital period, in seconds 
    """
    r_orb = (G * M_star)**(1/3) * (P_orb/2/np.pi)**(2/3)
    return r_orb

def n_wind(M_star_dot=3e-14, d=7e10, v_sw=25.6e5):
    """ Computes the particle density of the stellar wind at some distance d from the
        center of the star.
        OUTPUT: n_sw - particle density of the stellar wind at distance d, in #/cm^3
        INPUT:  M_star_dot - stellar mass-loss rate, M_sun/yr
                d          - distance from the center of the star, in cm
                v_sw       - Speed of stellar wind at distance d, in cm/s
    """
    M_star_dot *= M_sun/yr2sec  # mass-loss rate, in grams/sec
    m_av  =  1.92e-24 # average particle mass of the solar wind, in grams
    rho = M_star_dot/(4*np.pi * d**2 * v_sw*1e5) # Density of the stellar wind, in gr/cm3
    n_sw = rho / m_av
    return n_sw

#####################################
###  END OF FUNCTION DEFINITIONS  ###
#####################################

M_sun_dot = 2e-14 # Mass-loss rate of Sun, in M_sol/yr
# Proxima Cen
M_star_dot = 0.035*M_sun_dot
n_w = n_wind(M_star_dot, d=0.145*R_sun, v_sw=18.1)
#print("Particle density is n_sw = {0:.2e} #/cm^3".format(n_w))

