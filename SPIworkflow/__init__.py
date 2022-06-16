import numpy as np

### SETTING UP VALUES TO PREDICT SPI RADIO EMISSION

###################################
# OBSERVED PARAMETERS
###################################

# Observing frequency, in  Hz
#freq_obs = 400e6
# Assumed rms noise figure, in mJy
rms  = 0.022

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
n_sw_base = 1.0e7  # Appropriate for an M-dwarf star?

# Stellar magnetic field at the surface of the star pole, in gauss 
# Now it's read from an external file
#B_star = 307.

#####################################
# STELLAR WIND
#####################################

# Is the stellar plasma assumed to be isothermal?
isothermal = True 

#####################################
# SUB-ALFVENIC INTERACTION parameters
#####################################
 
# alpha - relative strength of the sub-Alfvénic interaction. 
# We assume that the planet has a highly conductive atmosphere, i.e., alpha = 1
alpha = 1

# Efficiency factor to convert Poynting flux into ECM radio emission.
eps_min = 0.01; eps_max = 0.1

# theta_M - Angle between planetary field and stellar field (rad) in the planet rest frame
theta_M = np.pi/2


#####################################
# PLOTTING AND WRITING SETUP
#####################################

# plotout = True => plot graphs to external files
# plotout = False => plot in terminal
plotout = True 
 
