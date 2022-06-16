import numpy as np

###########################
# Useful, general constants
###########################

c   = 2.99792458e10             #  speed of light, in cm/s
m_e = 9.1093897e-28             # electron mass, in g
m_p = 1.6726219e-24             # proton mass, in g
e   = 4.803204712570263e-10     # electron charge (esu)
k_B = 1.380658e-16              # Boltzmann's constant (erg/K)
m_h = 1.673534e-24              # Hydrogen mass, in g

pc = pc2cm   = 3.0857e18        # 1 pc, in cm
au      = 1.49597871e+13        # astonomical unit, in cm
mJy     = 1.0e-26               # 1 mJy, in erg/s/cm2/Hz

G     = 6.67408e-8              # Gravitational constant (cgs)
M_sun = Msun = 1.9891e33        # Mass of Sun, in g
R_sun = Rsun = 6.9599e10        # Radius of the Sun, in cm
M_earth = 5.9742e+27            # in g
R_earth = 6.378136e+8           # in cm
M_jup = M_jupiter  = 1.8987e+30 # in g
R_jup = R_jupiter  = 7.1492e+9  # in cm

yr = 365.24                     # yr, in days
day = 86400.                    # day,in seconds
yr2sec = yr*day                 #yr, in seconds

Omega_earth = 2*np.pi / 86400.  # Angular speed of Earth, in [s^(-1)] = 2*pi radians/day

Tesla2Gauss = 1e4  # 1 Tesla, in Gauss
Gauss2Tesla = 1e-4 # 1 Gauss, in Teslas

###############################
# Specific constants for SPI.py
###############################

bfield_earth = 3e-4              # Earth surface magnetic field at the pole, in Teslas
r_core_earth = 3.486e6           # Earth outer core radius, in [meters], from Sano 1993
rho_core_earth = 1.0e4           # Outer core density of Earth, in kg * m^(-3)

