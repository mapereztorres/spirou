import numpy as np
from SPIworkflow.constants import *


#####################################################
# ARRAY OF PLANETS TO COMPUTE RADIO EMISSION FROM SPI
#####################################################
# If COMPUTE_ALL = True, calcualte SPI radio emission for all targets in table.
COMPUTE_ALL = False
# If COMPUTE_ALL = False, then calculate the SPI radio emission for planets in array
# which_planets
which_planets = [20, 21]

###################################
### INPUT TABLE
###################################
# If EXOPLANETS = True, then we use SPI-targets.csv, which contains (confirmed) planets
EXOPLANETS = True

# If true, it reads a table with stellar systems hosting planets, so 
# the code expects info from the planets.
if EXOPLANETS == True:
    source_data = './INPUT/SPI-targets.csv'
else:
    source_data = './INPUT/SPI-NO-planets.csv'

selection_criteria = False


#####################################################
# ARRAY OF PLANETS TO COMPUTE RADIO EMISSION FROM SPI
#####################################################
#STUDY = "D_ORB"
STUDY = "M_DOT"
#STUDY = "B_PL"

#######################################################################
#  STUDY CASES
#  STUDY == 'D_ORB' - Predicted flux as a function of orbital distance
#  STUDY == 'M_DOT' - Predicted flux as a function of star mass-loss rate
#  STUDY == 'B_PL'  - Predicted flux as a function of planetary magnetic field
#######################################################################
#  STUDY = 'M_DOT' SETUP
#
# M_DOT_MIN, M_DOT_MAX: Minimum and maximum mass-loss rates to carry out the study of M_dot
# In units of M_dot_sun 
# M_DOT_STRETCH: Number of points per dex in the STUDY of M_DOT
M_DOT_STRETCH = 100
M_DOT_MIN = 1e-1
M_DOT_MAX = 1e+1

#  STUDY = 'B_PL' SETUP
#
# B_PL_MIN, B_PL_MAX: Minimum and maximum planetary magnetic field to carry out 
# the study of "B_PL". Values in Gauss
#
STEP = 0.05
B_PL_MIN = 0
B_PL_MAX = 20


### SETTING UP VALUES TO PREDICT SPI RADIO EMISSION

#####################################
# MAGNETIC FIELD SETUP
#####################################
# Setting the stellar magnetic field geometry and the value of the 
# intensity of the planetary magnetic field
# 
# Stellar magnetic field geometry
# The convention is that Bfield_geom_arr = 0 - closed dipolar geometry
#                        Bfield_geom_arr = 1 => open Parker spiral geometry; 
Bfield_geom_arr = [0]

# magnetized_pl_arr is a [False,True] array
# False: Unmagnetized planet 
# True : Magnetized planet
# 
#magnetized_pl_arr = [False, True]
magnetized_pl_arr = [True]

# Computation of planetary magnetic field
# B_pl_law = 'Sano' => Uses Sano's scaling law (Sano 1993)
# B_pl_law = 'None' => Doesn't use any scaling law. Uses B_planet_default instead.
B_planet_law = 'Sano'
#B_planet_law = 'None'

# Default planetary magnetic field, in Tesla
B_PLANET_DEFAULT = bfield_earth 


###################################
# OBSERVED PARAMETERS
###################################

# Observing frequency, in  Hz
#freq_obs = 400e6
# Assumed rms noise figure, in mJy
rms  = 0.030

# Representative bandwidth of the ECMI emission, in Hz
#Delta_nu_obs = freq_obs/2 

# Observed flux density at some given frequency, in mJy  (This should be modified in the
# future)
#flux = 1.0 

###################################
# STELLAR PARAMETERS
###################################

# Temperature of the corona,  
# base density at the corona 
# magnetic field at the surface of the star
# "isothermality" of stellar plasma 

# Base density using the empirical law from Peres+2004 (ApJ)
T_corona = 2.0e6 #A standard value (from soft X-ray observations of a number of M-dwarf stars)

# Density at the base of the corona nbase = 4.3e6*(T_corona/1e6)**4.2  (for
# Solar Type stars)
#n_sw_base = 1.0e8  # For a Sun-like star
# n_sw_base = 1.0e7  # Appropriate for an M-dwarf star?

# Mass-loss rate of the wind, in units of the Sun mass-loss-rate (2e-14 M_sol/year)
#M_star_dot = 0.23 # For GJ 436 (R_star=0.464 R_sun), this results in n_sw_base ~ 1e7 cm^-3 

# Stellar magnetic field at the surface of the star pole, in gauss 
# Now it's read from an external file
#B_star = 307.

#####################################
# STELLAR WIND
#####################################

# Is the stellar wind plasma assumed to be isothermal?
isothermal = True 

# Fully ionized, purely hydrogen stellar wind plasma (=> 50% protons, 50% electrons)
X_p = 0.5 # fraction of protons

#####################################
# SUB-ALFVENIC INTERACTION parameters
#####################################
 
# alpha - relative strength of the sub-Alfvénic interaction. 
# We assume that the planet has a highly conductive atmosphere, i.e., alpha = 1
alpha = 1

# Efficiency factor to convert Poynting flux into ECM radio emission.
#eps_min = 0.01; eps_max = 0.11
eps_min = 0.002; eps_max = 0.02

# theta_M - Angle between planetary field and stellar field (rad) in the planet rest frame
theta_M = np.pi/2

# Minimum and maximum kinetic energy of electrons emitting via ECM, in keV
# Electrons are mildly relativistic, so 20 to 200 keV energies should be good values, 
# but a choice from 10 to 511 keV (=m_e * c^2) should be also  OK. 
Ekin_min = 20 
Ekin_max = 200

##
which_beam_solid_angle = 'Jupiter-Io'
#which_beam_solid_angle = 'Computed'




#####################################
# PLOTTING AND WRITING SETUP
#####################################
## LINEWIDTH 
LW = 3
# PLOT ALSO Alfvén Mach Number, M_A?
PLOT_M_A = True
# PLOTOUT = True => plot graphs to external files
# PLOTOUT = False => plot in terminal
PLOTOUT = True 
# DRAW_RMS?
DRAW_RMS = True
# DRAW Little Earth?
DRAW_EARTH = True
