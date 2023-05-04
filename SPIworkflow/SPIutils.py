import numpy as np
import matplotlib
import matplotlib.pyplot as plt
from matplotlib.offsetbox import OffsetImage, AnnotationBbox

from SPIworkflow.__init__ import *
from SPIworkflow.constants import *


# FUNCTION DEFINITIONS 
def getImage(path):
    return OffsetImage(plt.imread(path, format="jpg"), zoom=.02)

def Kepler_r(M_star, P_orb):
    """Computes the orbital radius of a planet (=semimajor axis, in units of au)
       orbiting a star of mass M_star (in solar masses) and orbital period, P_orb (in days)
    """
    r_orb = (G * M_star*M_sun)**(1/3) * (P_orb*day/2/np.pi)**(2/3)/au
    return r_orb

def beta_keV(E_kin=511):
    """ Computes the velocity of an electron with kinetic energy E_k, in units of the speed of light
        OUTPUT: beta (= particle speed / speed of light)
        INPUT: E_kin - kinetic energy of the particle, in keV
    """ 
    E_0 = 511.0 # Rest mass of the electron, in keV
    beta = np.sqrt(1 - (1 + E_kin/E_0)**(-2)) 
    return beta

def bfield_sano(M_planet = 1.0, R_planet = 1.0, Omega_rot_planet = 1.0):
    """ Computes the surface magnetic field strength of a planet using the Sano scaling law.
        (From Sano, J. Geomag. Geolectr., 45, 65-77, 1993), 
        OUTPUT: B_planet: Surface magnetic field of the planet, in units of B_Earth 
        INPUT : M_planet: Planet radius, in units of the Earth mass
                R_planet: in units of the Earth radius
                Omega_rot_planet: Rotational speed of the planet. It is assumed that the 
                        planets are tidally locked, hence P_rot_planet = P_orb_planet. 
                        In units of the rotational speed of Earth (radians/sec) 
                rho_core - Planet outer core density, in g/cm^3. For now, fixed to the
                value of the density of the Earth's outer core. 
    """
    # Scaling law for the planet's core radius, from Curtis & Ness (1986)
    r_core = M_planet**0.44 # in units of r_core_Earth
    Magn_moment_planet = Omega_rot_planet * r_core**(7/2) # Magnetic moment, in units of Earth magn. moment
    B_planet  = Magn_moment_planet / R_planet**3 # in units of B_earth
 
    return B_planet

def Rp_eff_func(Rp, theta_M, Bp, B_tot):
    """ Computes the effective obstacle radius for the planet. 
        If the planet is magnetized, it is normally larger than the planet radius (Rp).
        If unmagnetized, then the effective radius is made equal to Rp.
        OUTPUT: Rp_eff - Effective planet radius, in cm 
        INPUT : Rp - Planet radius (cm)
                theta_M - angle  (radians), of the intrinsic planetary magnetic field (Bp) 
                          wrt the external magnetic field (B_tot).
                Bp -  planet magnetic field (G)
                B_tot - Magnetic field of the stellar wind at planet position (G)
    """
    if (Bp > 0.0):
        Rp_eff = Rp * np.sqrt(3*np.cos(theta_M/2)) * (Bp/B_tot)**(1./3.)
        Rp_eff[Rp_eff<Rp]=Rp # Rp_eff cannot be smaller than Rplanet
    else:
        Rp_eff = Rp
    
    return Rp_eff

def Lrad_leto(B_star=1.0, R_star=1.0, P_rot=1.0):
    """ Returns the radio luminosity of an early-type magnetic star, 
        using the scaling law of Leto et al. (2021, MNRAS)
        OUTPUT: Lnu_rad, in erg/s/Hz
        INPUT:  B_star, in Gauss
                R_star, in R_sun
                P_rot , in days
    """
    alpha = 13.6; beta = 1.94 # from Leto+2021
    nu_ecm = B_star * 2.8e6 # cyclotron freq, in Hz 
    Lnu_rad = 10**alpha * (B_star/1e3)**beta * R_star**(2*beta) * P_rot**(-beta) # erg/s/Hz
    Ltot_rad = Lnu_rad * nu_ecm
    return Ltot_rad
    
def n_wind(M_star_dot=1.0, d=7e10, v_sw=25.6e5, mu=0.5):
    """ Computes the particle density of the stellar wind at some distance d from the
        center of the star.
        OUTPUT: n_sw - particle density of the stellar wind at distance d, in #/cm^3
        INPUT:  M_star_dot - stellar mass-loss rate, in units of the Sun mass-loss rate
                d          - distance from the center of the star, in cm 
                v_sw       - Speed of stellar wind at distance d, in cm/s
                mu         - mean molecular weight in the wind
    """
    M_sun_dot = 2e-14 # Sun mass-loss rate, equal to 2e-14 gr/cm^3
    M_star_dot *= M_sun_dot * M_sun/yr2sec  # mass-loss rate, in grams/sec
    m_av  =  mu * m_p # average particle mass of the solar wind, in grams
    rho = M_star_dot/(4*np.pi * d**2 * v_sw) # Density of the stellar wind, in gr/cm3
    n_sw = rho / m_av
    return n_sw

def plasma_freq(n_e = 1.0):
    """ Returns the plasma frequency.
        Any wave emitting at a frequency less than the plasma frequency will
        not propate in this medium.

        OUTPUT: nu_plasma (plasma frequency), in Hz
        INPUT : n_e - electron density, in cm^(-3)
    """
    nu_plasma = 9.0e3 * np.sqrt(n_e)
    return nu_plasma
