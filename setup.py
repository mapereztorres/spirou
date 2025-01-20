import numpy as np
from SPIworkflow.constants import *

#####################################################
# ARRAY OF PLANETS TO COMPUTE RADIO EMISSION FROM SPI
#####################################################
# Choose whether to carry out the computations for a 
# table containing multiple targets (INPUT_TABLE = True), or
# for a single target (INPUT_TABLE = False). 

# Uncomment the line that applies
INPUT_TABLE = True
#INPUT_TABLE = False

#######################################################################
#  STUDY CASES
#  STUDY == 'D_ORB' - Predicted flux as a function of orbital distance
#  STUDY == 'M_DOT' - Predicted flux as a function of star mass-loss rate
#  STUDY == 'B_PL'  - Predicted flux as a function of planetary magnetic field
#  
#######################################################################

# Uncomment the line that applies
STUDY = "D_ORB"
#STUDY = "M_DOT"
#STUDY = "B_PL"

# STUDY = "D_ORB" SETUP
#D_ORB_LIM = np.nan
D_ORB_LIM = 3000

#  STUDY = 'M_DOT' SETUP
#
# M_DOT_MIN, M_DOT_MAX: Minimum and maximum mass-loss rates to carry out the study of M_dot
# In units of M_dot_sun 
# M_DOT_STRETCH: Number of points per dex in the STUDY of M_DOT
M_DOT_STRETCH = 50
M_DOT_MIN = 1e-1
M_DOT_MAX = 1e+2

#  STUDY = 'B_PL' SETUP
#
# B_PL_MIN, B_PL_MAX: Minimum and maximum planetary magnetic field to carry out 
# the study of "B_PL". Values in Gauss
#
STEP = 0.05
B_PL_MIN = 0
B_PL_MAX = 4 

##
# Distance (from the centre of the star) where SPI emission takes place (in units of R_star)
R_SPI = 1.0

####################################################
# Stellar wind FREE ABSORPTION of SPI radio emission
####################################################

### Consider free-free absorption (True => Yes; False => No)
freefree = False

# Ionization state (Z = 1 - fully ionized hydrogen)
Z = 1 

# Distance, measured from R_SPI, in units of R_star, where free-free absorption becomes negligible. 
R_ff_OBSERVER = 10000

### NSTEPS_FF: number of points for the distance array
#NSTEPS_FF = 1000000
NSTEPS_FF = 10000
#####################################################
### SETTING UP VALUES TO PREDICT SPI RADIO EMISSION
#####################################################

#####################################
# STELLAR MAGNETIC FIELD SETUP
#####################################
# Setting the stellar magnetic field geometry and the value of the 
# intensity of the planetary magnetic field
# 
# Stellar magnetic field geometry
# The convention is that Bfield_geom_arr = 0 - closed dipolar geometry
#                        Bfield_geom_arr = 1 => open Parker spiral geometry; 
#Bfield_geom_arr = [0,1]
Bfield_geom_arr = [0] 

# MAGN_OBLIQ - magnetic obliquity. Angle betw the magnetic and rotation axes of the star
# (in degrees). Fixed to zero for simplicity. 
MAGN_OBLIQ = 0.
 
#####################################
# PLANET MAGNETIC FIELD SETUP
#####################################

# magnetized_pl_arr is a [False,True] array
# False: Unmagnetized planet 
# True : Magnetized planet
# 
#magnetized_pl_arr = [False, True]
magnetized_pl_arr = [True]

# Default planetary magnetic field, in Tesla
B_PLANET_DEFAULT = bfield_earth

# Setting the stellar magnetic field geometry and the value of the 
# Computation of planetary magnetic field 
# B_pl_law = 'Sano' => Uses Sano's scaling law (Sano 1993)
# B_pl_law = 'None' => Doesn't use any scaling law. Uses B_PLANET_DEFAULT instead.
#B_planet_law = 'Sano'
B_planet_law = 'None'

# K_MAGNETOPAUSE - factor by which the magnetopause currents enhance
# the magnetospheric magnetic field at the magnetopause, which is a value
# between 2 and 3. ADD REFERENCE
K_MAGNETOPAUSE = 2.0

###################################
# OBSERVED PARAMETERS
###################################

# Observing frequency, in  Hz
#freq_obs = 400e6
# Assumed rms noise figure, in mJy
RMS = 0.015
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

# Default value: 2.0e6 K. Solar value
T_CORONA_DEF = 2.0e6 




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
 
# ALPHA_SPI - relative strength of the sub-Alfvénic interaction. 
# We assume that the planet has a highly conductive atmosphere, i.e., ALPHA_SPI = 1
ALPHA_SPI = 1

# Efficiency factor in converting Poynting flux into ECM radio emission.
#
# We use the naming convention in Zarka (2018) and Zarka (2024), in the 
# Handbook of Exoplanets, i.e., "beta" (his Eq. 1 in Zarka 2024), 
# which we rename as "beta_eff". 
# Zarka (2024) shows that beta_eff is in the range from 1e-4 (BETA_EFF_MIN) to 1e-2
# (BETA_EFF_MAX)
# 
BETA_EFF_MIN = 1e-4; BETA_EFF_MAX = 1e-2 

# Fraction of Poynting flux going to dissipated power (Eq. 2 in Zarka 2024)
# Zarka (2024) quotes a value of 0.2 +/- 0.1
EPSILON = 0.2 

# theta_M - Angle between planetary field and stellar field (rad) in the planet rest frame
THETA_M = np.pi/2

# Minimum and maximum kinetic energy of electrons emitting via ECM, in keV
# Electrons are mildly relativistic, so 20 to 200 keV energies should be good values, 
# but a choice from 10 to 511 keV (=m_e * c^2) should be also  OK. 
EKIN_MIN = 1 
EKIN_MAX = 20

## Beam solid angle of the ECM radio emission
## 
## Choose whether to compute the beam solid angle (COMPUTE_BSA = True), 
## or use fixed values (COMPUTE_BSA = False)
COMPUTE_BSA = False

## The standard avlue for OMEGA_MIN for SPI is
## taken from the Io-Jupiter interaction.
## OMEGA_JUPITER_IO = 0.16 sterradians (DAM emission from a single flux tube)

OMEGA_MIN = OMEGA_JUPITER_IO
# OMEGA_MAX due to star-planet interaction is at most a few times OMEGA_MIN
# Note: Many papers have wrongly assumed OMEGA_MAX = 1.6 (sterr), which is the 
# beam solid angle of the whole Jupiter auroral oval. 

# OMEGA_MAX = OMEGA_MIN
OMEGA_MAX = 3*OMEGA_JUPITER_IO

#####################################
# PLOTTING AND WRITING SETUP
#####################################
#Plot graph for M_A
PLOT_M_A = False
## LINEWIDTH 
LW = 3
# PLOTOUT = True => plot graphs to external files
# PLOTOUT = False => plot in terminal
PLOTOUT = True 
# DRAW_RMS?
DRAW_RMS = True
# DRAW Little Earth?
DRAW_EARTH = True
LIMS_MA = True
LIM_MA_LOW = 1e-2
LIM_MA_HIGH = 1e0
FLUX_LIMS = True
FLUX_LOW  = 3*RMS * 1e-1
FLUX_HIGH = 3*RMS * 1e2
YLIMLOW   = 1e-3
YLIMHIGH  = 1e2


### WELCOME TO SPIROU
print('###########################################################')
print('###                                                     ###')
print('### SPIROU code @Pérez-Torres & Peña-Moñino 2025        ###')
print('###                                                     ###')
print('### If you use fully or partial parts of this code,     ###')
print('### please acknowledge our work by citing the paper     ###') 
print('### XXXX and indicating the url address for the code:   ###') 
print('### https://github.com/mapereztorres/spirou             ###') 
print('###                                                     ###')
print('###########################################################\n')

###
print('You have chosen the following setup parameters for SPIROU:\n')

###
if INPUT_TABLE == False:
    print(f'RUNNING THE CODE FOR A SINGLE OBJECT: INPUT_TABLE = {INPUT_TABLE}\n')
else:
    print(f'RUNNING THE CODE FOR MULTIPLE OBJECTS: INPUT_TABLE = {INPUT_TABLE}\n')

###
if STUDY == "D_ORB":
    print('CARRYING OUT A STUDY OF RADIO EMISSION VS ORBITAL SEPARATION: STUDY == D_ORB\n')
elif STUDY == "M_DOT":
    print('CARRYING OUT A STUDY OF RADIO EMISSION vs STELLAR MASS LOSS RATE: STUDY == M_DOT\n ')
elif STUDY == "B_PL":
    print('CARRYING OUT A STUDY OF RADIO EMISSION VS PLANETARY MAGNETIC FIELD: STUDY == B_PL\n')

###

if len(Bfield_geom_arr)==2:
    print('RUNNING FOR BOTH AN OPEN PARKER SPIRAL AND \n A (CLOSED) DIPOLAR MAGNETIC FIELD GEOMETRIES\n')
else:
    if Bfield_geom_arr[0] == 0:
        print('RUNNING FOR A CLOSED DIPOLAR MAGNETIC FIELD GEOMETRY\n')    
    else:
        print('RUNNING FOR AN OPEN PARKER SPIRAL MAGNETIC FIELD GEOMETRY\n')    
        
###    
if len(magnetized_pl_arr)==2:    
    print('RUNNING FOR BOTH A MAGNETIZED AND A NON-MAGNETIZED PLANET\n')
else:
    if magnetized_pl_arr[0] == FALSE:
        print('RUNNING FOR A NON-MAGNETIZED PLANET\n')    
    else:
        print('RUNNING FOR A MAGNETIZED PLANET\n')    
        
if (B_planet_law) == 'Sano':
    print('THE PLANETARY MAGNETIC FIELD IS ESTIMATED FOLLOWING SANO\'s SCALING LAW\n')
else:
    print('THE PLANETARY MAGNETIC FIELD IS ASSUMED TO BE THAT OF THE EARTH\n')


###
if freefree == True:
    print('FREE FREE ABSORPTION EFFECTS ARE CONSIDERED\n')
else:
    print('FREE FREE ABSORPTION EFFECTS ARE IGNORED\n')    
    
###
print(f'THE EFFICIENCY FACTOR IN CONVERTING POYNTING FLUX \nTO RADIO EMISSION RANGES FROM {BETA_EFF_MIN:.1e} UP TO {BETA_EFF_MAX:.1e}\n')

###
if COMPUTE_BSA == False: 
    print(f'THE BEAM SOLID ANGLE OF THE AURORAL EMISSION RANGES \nFROM {OMEGA_MIN:.2f} STERRADIANS UP TO {OMEGA_MAX:.2f} STERRADIANS\n')
else:
    print(f'THE BEAM SOLID ANGLE OF THE AURORAL EMISSION RANGES WILL BE COMPUTED BASED\nON THE VALUES OF THE MINIMUM AND MAXIMUM KINETIC ENERGY OF THE ELECTRONS\nEMITTING VIA THE ECM MECHANISM: FROM {EKIN_MIN:.0f} KEV UP TO {EKIN_MAX:.0f} KEV\n')

### 
print(f'You are assuming an observing rms value of {RMS*1e3:.1f} μJy\n')    

print('###########################################################\n')
