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
    """Computes the orbital radius of a planet (=semimajor axis)
       orbiting a star of mass M_star (in solar masses) and orbital period, P_orb (in days)
    """
    r_orb = (G * M_star)**(1/3) * (P_orb/2/np.pi)**(2/3)
    return r_orb

def beta_keV(E_kin=511):
    """ Computes the velocity of an electron with kinetic energy E_k, in units of the speed of light
        OUTPUT: beta (= particle speed / speed of light)
        INPUT: E_kin - kinetic energy of the particle, in keV
    """ 
    E_0 = 511.0 # Rest mass of the electron, in keV
    beta = np.sqrt(1 - (1 + E_kin/E_0)**(-2)) 
    return beta

def bfield_Sano(v_orb=1.0, r_orb=1.0, M_planet = M_earth, rho_core=rho_core_earth):
    """ Computes the surface magnetic field strength of a planet using the Sano scaling law.
        OUTPUT: B_planet, in Gauss
        INPUT : v_orb, M_planet, R_planet, rho_core
                v_orb - orbital speed, in cm/s
                M_planet - Planet mass, in g
                rho_core - Planet outer core density, in g/cm^3. For now, fixed to the
                value of the oute core density of the Earth.
    """
    Omega_earth = 2*np.pi / 86400.  # Angular speed of Earth, in [s^(-1)] = 2*pi radians/day
    Omega_planet = v_orb / r_orb  # Angular speed of planet, in [s^(-1)]. Assumes the planet is tidally locked
    r_core = r_core_earth * (M_planet/M_earth)**0.44 # Scaling law from Curtis & NEss (1986)
    scaling_law = (rho_core / rho_core_earth )**(1/2) * (r_core / r_core_earth)**(7/2) * (Omega_planet / Omega_earth)
    magnetic_moment_planet = magnetic_moment_earth * scaling_law
    B_planet  = 2 * magnetic_moment_planet / r_core**3 # in Tesla.  
    B_planet *= 1e-4 # In Gauss
 
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
    

