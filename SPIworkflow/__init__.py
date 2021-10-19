### SETTING UP VALUES TO PREDICT SPI RADIO EMISSION

# OBSERVED PARAMETERS

# Observing frequency, in  Hz
freq_obs = 400e6
# Assumed rms noise figure, in mJy
rms  = 0.013 

# Representative bandwidth of the ECMI emission, in Hz
Delta_nu_obs = freq_obs/2 

# Observed flux density at some given frequency, in mJy  (This should be modified in the
# future)
flux = 1.0 

# STELLAR PARAMETERS
# Temperature of the corona,  base density at the corona, 

# Base density using the empirical law from Peres+2004 (ApJ)
#
T_corona = 2.0e6 #A standard value (from soft X-ray observations of a number of M-dwarf stars)

# Density at the base of the corona nbase = 4.3e6*(T_corona/1e6)**4.2 
n_sw_base = 1e7

# Stellar magnetic field at the surface of the pole, in gauss 
B_star = 53.

# Is the stellar plasma assumed to be isothermal?
isothermal = True 

#
# alpha - relative strength of the sub-Alfvénic interaction We assume that the planet 
# has a highly conductive atmosphere => alpha = 1
#
# theta_M - Angle between planetary field and stellar field (rad) in the planet rest frame
alpha = 1
theta_M = np.pi/2

