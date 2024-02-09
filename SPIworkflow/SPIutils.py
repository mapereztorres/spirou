import numpy as np
import matplotlib
import matplotlib.pyplot as plt
from matplotlib.offsetbox import OffsetImage, AnnotationBbox

from SPIworkflow.__init__ import *
from SPIworkflow.constants import *


# FUNCTION DEFINITIONS 

def get_bfield_comps(open_field, B_star, d_orb, R_star, v_corot, v_sw, v_rel_angle):
    """Computes the radial, azimuthal components of the stellar wind magnetic field for
    an open (Parker spiral) and a closed (dipolar) magnetic field topologies of the host
    star.
    OUTPUT:
    INPUT: 
    ...
    """
    if open_field: 
        # Open Parker Spiral - Falls down with distances as R^(-2) rather than R^(-3) as in the dipole case
        B_r = B_star * (d_orb/R_star)**(-2) # Stellar wind B-field at (R/R_star), Eqn 20 Turnpenney 2018
        B_phi = B_r * v_corot/v_sw # Azimuthal field (Eqn 21 Turnpenney 2018)
        B_sw = np.sqrt(B_r**2 + B_phi**2) # Total stellar wind B-field at planet orbital distance
    else:
        # Closed, dipolar configuration - It falls with distance as R^(-3)
        # B_star - magnetic field at the magnetic equator on the stellar surface
        # 
        phi = 0. # azimuth, measured from the North magnetic pole of the star (in degrees)
        phi *= np.pi/180. # to radians
        B_r   = -2 * B_star * (d_orb/R_star)**(-3) * np.cos(phi) # Radial component of the dipole magnetic field of the stellar wind as f(distance to star)
        B_phi = - B_star * (d_orb/R_star)**(-3) * np.sin(phi) # Azimuthal component of the dipole magnetic field 
        B_sw = np.sqrt(B_r**2 + B_phi**2) # Total dipole magnetic field 

    # Eq. 23 of Turnpenney 2018 -  First term of RHS
    B_ang = np.arctan(B_phi/B_r) # Angle the B-field makes with the radial direction

    # Angle between the stellar wind magnetic field and the impinging plasma velocity
    # Eqn 23 in Turnpenney 2018. It's also Eq. 13 in Zarka 2007
    theta = np.absolute(B_ang - v_rel_angle) 

    geom_f = 1.0 # Geometric factor. 1 for closed dipole configuration, different for the open field configuration
    if open_field:
        # theta is the angle between the B_sw (the insterstellar magnetic field), and the
        # incident stellar wind velocity.  See Fig. 1 in Turnpenney+2018
        #
        geom_f = (np.sin(theta))**2 # Geometric factor in efficiency 
    return B_r, B_phi, B_sw, B_ang, theta, geom_f

def getImage(path):
    return OffsetImage(plt.imread(path, format="jpg"), zoom=.02)

def Kepler_r(M_star, P_orb):
    """Computes the orbital radius of a planet (=semimajor axis, in units of au)
       orbiting a star of mass M_star (in solar masses) and orbital period, P_orb (in days)
    """
    r_orb = (G * M_star*M_sun)**(1/3) * (P_orb*day/2/np.pi)**(2/3)/au
    return r_orb

def Kepler_P(M_star, r_orb):
    """Computes the orbital period (P_orb) of a planet (in days) of a planet
       orbiting a star of mass M_star (in solar masses) and semi-major axis (r_orb), in au.
    """
    
    P_orb = 2*np.pi * (G * M_star*M_sun)**(-1/2) * (r_orb*au)**(3/2)/day
    return P_orb

def beta_keV(E_kin=511):
    """ Computes the velocity of an electron with kinetic energy E_k, in units of the speed of light
        OUTPUT: beta (= particle speed / speed of light)
        INPUT: E_kin - kinetic energy of the particle, in keV
    """ 
    E_0 = 511.0 # Rest mass of the electron, in keV
    beta = np.sqrt(1 - (1 + E_kin/E_0)**(-2)) 
    return beta

def Rp_Zeng(Mp = 1.0):
    """ Computes the planet radius using the mass-radius relationship in 
        in Zeng et al. (PNAS, 2019): 
        https://www.pnas.org/content/116/20/9723
        R [R_earth] = f * (M/M_earth)^(1/3.7)  (for rocky planets, f=1.)
        
        OUTPUT: Rp - planet radius, in units of Earth radius 
        INPUT:  Mp - planet mass, in units of Earth mass
        We assume f = 1 (Earth-like rocky planet)
    """
    f = 1.0
    Rp = f * Mp**(1/3.7)
    return Rp

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

def R_planet_eff_func(Rp, theta_M, B_planet, B_tot):
    """ Computes the effective obstacle radius for the planet. 
        If the planet is magnetized, it is normally larger than the planet radius
        (R_planet).
        If unmagnetized, then the effective radius is made equal to Rp.
        OUTPUT: R_planet_eff - Effective planet radius, in cm 
        INPUT : Rp - Planet radius (cm)
                theta_M - angle  (radians), of the intrinsic planetary magnetic field
                (B_planet) 
                          wrt the external magnetic field (B_tot).
                B_planet -  planet magnetic field (G)
                B_tot    - Magnetic field of the stellar wind at planet position (G)
    """
    if (B_planet > 0.0):
        R_planet_eff = Rp * np.sqrt(3*np.cos(theta_M/2)) * (Bp/B_tot)**(1./3.)
        R_planet_eff[R_planet_eff < Rp] = Rp # R_planet_eff cannot be smaller than R_planet
    else:
        R_planet_eff = Rp
    
    return R_planet_eff

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
    M_sun_dot = 2e-14 # Sun mass-loss rate, in Msun/yr
    M_star_dot *= M_sun_dot * M_sun/yr2sec  # mass-loss rate, in grams/sec
    m_av  =  mu * m_p # average particle mass of the solar wind, in grams
    rho = M_star_dot/(4*np.pi * d**2 * v_sw) # Density of the stellar wind, in gr/cm3
    n_sw = rho / m_av
    return n_sw

def plasma_freq(n_e = 1.0):
    """ Returns the plasma frequency.
        Any wave emitting at a frequency less than the plasma frequency will
        not propate in this medium.

        OUTPUT: 
              nu_plasma - plasma frequency, in Hz

        INPUT : n_e - electron density, in cm^(-3)
    """
    nu_plasma = 9.0e3 * np.sqrt(n_e)
    return nu_plasma

def get_P_dyn_sw(n_sw=1e7, mu=0.5, v_rel=1e7):
    """Computes the dynamical pressure of the stellar wind at the orbital distance 
               of the planet, P_dyn_sw, in erg/cm3

       OUTPUT: P_dyn_sw, in erg/cm3

       INPUT : n_sw  - number density of the stellar wind at the planet's
                       orbital position, in  # cm^(-3)
               mu    - mean molecular weight 
               v_rel - Relative velocity between the stellar wind to the
                       orbital velocity of the planet, in cm/s
    """
    rho_sw  = mu * m_p * n_sw
    P_dyn_sw = rho_sw * v_rel**2

    return P_dyn_sw

def get_P_th_sw(n_sw=1e7, T_e=2e6):
    """Computes the thermal pressure of the stellar wind at the orbital distance 
                of the planet, P_th_sw, in erg/cm3

       OUTPUT: P_th_sw, in erg/cm3

       INPUT : n_sw - number density of the stellar wind at the planet's
                      orbital position, in  # cm^(-3)
               k_B  - Boltzmann's constant
               T_e  - Temperature of the stellar wind ath the planet's
                      position. 
    """
    P_th_sw = n_sw * k_B * T_e

    return P_th_sw

def get_P_B_sw(B_sw=1.0):
    """Computes the  pressure of the stellar wind at the orbital distance 
               of the planet, P_dyn_sw, in erg/cm3

       OUTPUT: 
              P_B_sw, in erg/cm3

       INPUT : 
              B_sw  - Stellar wind magnetic field, in Gauss
    """

    P_B_sw = B_sw**2 / (8 * np.pi)

    return P_B_sw

def get_P_B_planet(B_planet = 1.0):
    """Computes the magnetic pressure of planet at the equator, P_Bp, in erg/cm3

       OUTPUT: 
              P_B_planet - magnetic pressure of planet at the equator, in erg/cm3

       INPUT: 
              B_planet - planetary magnetic field (dipole) at the poles, in Gauss

    """
    P_B_planet = (B_planet / 2)**2 / (8*np.pi)

    return P_B_planet

def R_magnetopause(B_planet = 2.0, P_dyn_sw = 1.0, P_th_sw=1.0, B_sw = 1.0):
    """Returns the radius of the planetary magnetopause, in units of Rp

    It's calculated as the radius that balances the pressure between the
    stellar wind (ram, thermal, and magnetic) and the planet (magnetic)
    pressure. Note that, in the equation below, the ram and thermal components
    of the planetary pressure are neglected.

    OUTPUT: R_mp - Radius of the planetary magnetopause, in cm

    INPUT : B_planet - planetary magnetic field (dipole) at the poles, in Gauss
            P_dyn_sw - dynamic pressure of the stellar wind at the orbital distance 
                       of the planet, in erg/cm3
            P_th_sw  - thermal pressure of the stellar wind at the orbital distance 
                       of the planet, in erg/cm3
            B_sw     - Stellar wind magnetic field, in Gauss
    """
    # Pressure of the planetary magnetic field at the equator 
    P_B_planet_equator = (B_planet/2)**2 / (8*np.pi)

    # Pressure of the stellar wind magnetic field (at the equator)
    P_B_sw = B_sw**2 / (8 * np.pi)
   
    # Total pressure of the stellar wind 
    P_sw = P_dyn_sw + P_th_sw + P_B_sw
    
    # Expression in Donati & Vidotto (2017).
    #
    # Similar to others in, e.g., Zarka (2007), Turnpenney+2018, etc., but in
    # units of Rp
    Rmp = 2**(1./3.) * (P_B_planet_equator / P_sw)**(1./6) 

    return Rmp

def get_Rmp(P_Bp=1.0, P_dyn_sw=1.0, P_th_sw=1.0, P_B_sw=1.0):
    """Computes the radius of the planetary magnetopause, in units of planet radii

    It's calculated as the radius that balances the pressure between the
    stellar wind (ram, thermal, and magnetic) and the planet (magnetic)
    pressure. Note that, in the equation below, the ram and thermal components
    of the planetary pressure are neglected.

    OUTPUT: Rmp - Radius of the planetary magnetopause, in units of Rp

    INPUT : P_Bp     - pressure of the planetary magnetic field (dipole) at the equator
            P_dyn_sw - dynamic pressure of the stellar wind at the orbital distance 
                       of the planet, in erg/cm3
            P_th_sw  - thermal pressure of the stellar wind at the orbital distance 
                       of the planet, in erg/cm3
            B_sw     - Stellar wind magnetic field, in Gauss
    """
    # Total pressure of the planet (only magnetic pressure is considered)
    P_planet = P_Bp

    # Total pressure of the stellar wind 
    P_sw = P_dyn_sw + P_th_sw + P_B_sw
    
    # Expression in Donati & Vidotto (2017).
    # Similar to others in, e.g., Zarka (2007), Turnpenney+2018, etc.
    Rmp = 2**(1./3.) * (P_planet / P_sw)**(1./6) 

    return Rmp
    
def B_starmass(star_mass,Prot):
  """Calculation of the Stellar magnetic field at the surface as a function of the stellar mass and its rotation period
  The value of the magnetic field depends on the Rossby number (Ro), as defined in https://www.aanda.org/articles/aa/pdf/2022/06/aa43251-22.pdf.
  There are two relations, depending no whether the star rotates fast (Ro<0.13) or slow (Ro>0.13)
  
  The Rossby number in turn is calculated using the stellar mass, as defined in equation 6 of https://academic.oup.com/mnras/article/479/2/2351/5045253.
  
  
  OUTPUT: B_mass     -  Stellar magnetic field at the surface (in G)
  
  
  INPUT : star_mass  -  Mass of the star (solar units)
          Prot       -  Rotation period of the star (in days)
  
  
  
  """
  #https://www.aanda.org/articles/aa/pdf/2022/06/aa43251-22.pdf (table )
  #https://academic.oup.com/mnras/article/479/2/2351/5045253  (eq 6)
  C=2.33
  D=-1.50
  E=0.31
  alpha=-1.26
  alpha_fast=-0.11
  tau_mass=10**(C+D*star_mass+E*star_mass**2)
  Ro_mass=Prot/tau_mass
  if Ro_mass>0.13:
    print('Slow rotator')
    B_mass= 199*Ro_mass**alpha
  elif Ro_mass<0.13:
    print('Fast rotator')
    B_mass= 2050*Ro_mass**alpha_fast    
  else:
    print('Uncertain regime: Ro=',Ro_mass)
  return B_mass

def B_color(starname,star_mass,Prot):
  #in development
  from astroquery.simbad import Simbad
  Simbad.add_votable_fields('flux(V)')
  Simbad.add_votable_fields('flux(K)')
  #https://www.aanda.org/articles/aa/pdf/2022/06/aa43251-22.pdf (table )
  #https://academic.oup.com/mnras/article/479/2/2351/5045253  (eq 5)
  A=0.64
  B=0.25
  alpha=-1.26
  #eachstar=data['star_name'][indi]
  result_table = Simbad.query_object(starname)
  print(result_table['MAIN_ID'][0],result_table['FLUX_V'][0],result_table['FLUX_K'][0])
  #pdn['B_color'][ind]=result_table['FLUX_B'][0]
  V_color=result_table['FLUX_V'][0]
  K_color=result_table['FLUX_K'][0]  
  VK=V_color-K_color
  logtau_color=A+B*VK
  tau_color=10**(A+B*VK)
  Ro_color=Prot/tau_color
  B_color= 199*Ro_color**alpha

  return B_color


def Mdot_star(R_star, M_star, Prot_star):
    """
    OUTPUT: Stellar mass-loss rate, in Msolar/yr
    INPUT: Radius of star (in units of Rsun)
           M_star    (in units of Msun)
           Prot_star (days)
           
           It uses Eq. 7 in Johnstone and Güdel (2015, A&A, 577, A27) = J+2015
           to estimate the mass-loss rate of a low-mass main sequence star. 
           M_dot_sun_fit and Omega_sun_fit are taken as defined in J+2015. 
           The exponents a and b below are best fits obtained by J+2015, using the 
           fixed values of M_dot_sun_fit (1.4e-14 Msolar/yr) and Omega_sun_fit (2.67e-6
           rad/sec). 
           NOTE that M_dot_sun_fit is taken to be 1.4e-14 Msolar/yr in J+2015, 
           instead of the usual value of 2E-14 Msolar/yr.

    """
    #R_sun = 7e10 #Sun radius, in cm 
    #M_sun = 2e33 #Sun mass, in g
    R_sun = 1; M_sun = 1
    a = 4./3
    b = -10./3 
    M_dot_sun_fit = 1.4e-14 #Sun mass-loss rate, in g/s
    Omega_sun_fit = 2.67e-6 # Sun angular speed, in rad/s - equals to Prot_sun = 27.2367 days
    Omega_star = 2*np.pi / (Prot_star*86400)
    Mdot = M_dot_sun_fit * (R_star/R_sun)**2 * (Omega_star/Omega_sun_fit)**a * (M_star/M_sun)**b
    return Mdot
