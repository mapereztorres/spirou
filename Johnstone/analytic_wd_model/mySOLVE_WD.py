
import math
import matplotlib
import matplotlib.pyplot as plt
import numpy
from scipy import optimize
import pickle

matplotlib.rcParams['axes.titlesize'] = 20
matplotlib.rcParams['axes.labelsize'] = 20
matplotlib.rcParams['xtick.labelsize'] = 20
matplotlib.rcParams['ytick.labelsize'] = 20
matplotlib.rcParams['lines.linewidth'] = 3
matplotlib.rcParams['savefig.bbox'] = 'tight'
matplotlib.rcParams['axes.linewidth'] = 2.0
matplotlib.rcParams['xtick.major.width'] = 2.0
matplotlib.rcParams['xtick.minor.width'] = 2.0
matplotlib.rcParams['ytick.major.width'] = 2.0
matplotlib.rcParams['ytick.minor.width'] = 2.0
matplotlib.rcParams['ytick.major.pad'] = 6
matplotlib.rcParams['ytick.minor.pad'] = 6
matplotlib.rcParams['xtick.major.pad'] = 6
matplotlib.rcParams['xtick.minor.pad'] = 6


# =======================================================================
# =======================================================================

def get_nominator_denominator(X):
  
  R = X[0]
  Vr = X[1]
  
  Ar = Fb * Vr**0.5 / ( R * ( 4.0*3.142*Fm )**0.5 ) 
  
  Vphi = R * Omega * ( (L*Vr**2.0)/(R**2.0 * Omega) - Ar**2.0 ) / (Vr**2.0 - Ar**2.0)
  
  Aphi = Ar * (Vphi - R*Omega) / Vr
  
  nominator = (Vr**2.0 - Ar**2.0) * (2.0*cs**2.0 + Vphi**2.0 - G*Mstar/R) + 2.0*Vr*Vphi*Ar*Aphi
  denominator = (Vr**2.0 - Ar**2.0) * (Vr**2.0 - cs**2.0) - Vr**2.0 * Aphi**2.0
    
  return [nominator , denominator]


# =======================================================================
# =======================================================================

def get_slow(rA,vA,L):
  global nClick , R0click , V0click
  # starting points for solution

  R0 = rA/35.0
  V0 = vA/2.0
  
  if (nClick < nClickMax):
    R0click,V0click = guess_slow()
  R0 = R0click * rA
  V0 = V0click * vA
  nClick += 1
  
  slow = optimize.root(get_nominator_denominator, [R0,V0],method='lm',tol=10.0**(-25))
  
  return slow.x

# =======================================================================
# =======================================================================

def get_fast(rA,vA,L):
  
  # starting points for solution
  R0 = rA*1.2
  V0 = vA*1.1
  

  fast = optimize.root(get_nominator_denominator, [R0,V0],method='lm',tol=10.0**(-25))
  
  return fast.x
  
# =======================================================================
# =======================================================================

def get_slow_fast(V):
  
  

  # get Alfven point
  vA = cs / math.sqrt(V)
  rA = Fb  / math.sqrt(Fm*4.0*3.142*vA)
  
  
  # get L (angular momentum)
  L = Omega * rA**2.0

  
  # get slow and fast points
  Rslow,Vslow = get_slow(rA,vA,L)
  Rfast,Vfast = get_fast(rA,vA,L)
  
  
  return Rslow,Vslow,Rfast,Vfast



# =======================================================================
# =======================================================================

def get_LHS_minus_RHS(Vr):
  
    
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
  U = G * Mstar / (rA * vA**2.0)
  V = cs**2.0 / vA**2.0 
  W = (rA * Omega)**2.0 / vA**2.0

  x = R / rA
  u = Vr / vA
  LHS = W + 0.5 * u**2.0 + 0.5 * W * x**2.0 * (1-u)**2.0 / (1-x**2.0*u)**2.0 
  
  RHS = V * math.log(x**2.0) + V * math.log(u) + U/x + W * x**2.0 * (1.0 - u) / (1.0 - x**2.0 * u) + const
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 

    
  return LHS - RHS
  

# =======================================================================
# =======================================================================

def find_solution(Vmin,Vmax):
  
  # first test if fa and fb have different signs 
  
  bam = get_LHS_minus_RHS(Vmin)*get_LHS_minus_RHS(Vmax)
  
  
  if (bam < 0):
    Vsolution = optimize.brentq(get_LHS_minus_RHS,Vmin,Vmax)
  else:
    # make small grid and find minimum
    nV = 1000
    Vgrid = numpy.linspace(Vmin,Vmax,nV)
    funcgrid = numpy.zeros(nV)
    for i in range(0,nV):
      funcgrid[i] = get_LHS_minus_RHS(Vgrid[i])
    # get min
    Vsolution = Vgrid[numpy.where(funcgrid == numpy.amin(funcgrid))]
    if (len(Vsolution) > 1):
      Vsolution = Vsolution[0]    

    Vsolution = -1
  
  return Vsolution
    
# =======================================================================
# =======================================================================

def plot_zeros():
  
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  # make grid of Rs and Vs

  RsMin = Rstar*1.1
  RsMax = Rfast*1.3
  VsMin = vA/1000
  VsMax = Vfast*1.3

  nRs = 500
  nVs = 500

  RsAll = numpy.linspace(RsMin,RsMax,nRs)
  VsAll = numpy.linspace(VsMin,VsMax,nVs)

  grid1 = numpy.zeros((nRs,nVs))
  grid2 = numpy.zeros((nRs,nVs))


  # scroll through possibilities
  for iR in range(0,nRs):
    for iV in range(0,nVs):
      
      
      Ar = Fb * VsAll[iV]**0.5 / ( RsAll[iR] * ( 4.0*3.142*Fm )**0.5 ) 
      
      Vphi = RsAll[iR] * Omega * ( (L*VsAll[iV]**2.0)/(RsAll[iR]**2.0 * Omega) - Ar**2.0 ) / (VsAll[iV]**2.0 - Ar**2.0)
      
      Aphi = Ar * (Vphi - RsAll[iR]*Omega) / VsAll[iV]
      
      grid1[iV][iR] = (VsAll[iV]**2.0 - Ar**2.0) * (2.0*cs**2.0 + Vphi**2.0 - G*Mstar/RsAll[iR]) + 2.0*VsAll[iV]*Vphi*Ar*Aphi
      grid2[iV][iR] = (VsAll[iV]**2.0 - Ar**2.0) * (VsAll[iV]**2.0 - cs**2.0) - VsAll[iV]**2.0 * Aphi**2.0



  bam1 = numpy.where(grid1 > 0)
  bam2 = numpy.where(grid1 < 0)
  grid1[bam1] = 10
  grid1[bam2] = -10

  bam3 = numpy.where(grid2 > 0)
  bam4 = numpy.where(grid2 < 0)
  grid2[bam3] = 10
  grid2[bam4] = -10


  X,Y = numpy.meshgrid(RsAll/rA,VsAll/vA)

  plt.figure()

  levels = numpy.array([ -numpy.amax(grid1) , 0.0 , numpy.amax(grid1) ])
  plt.contour(X,Y,grid1,levels)

  levels = numpy.array([ - numpy.amin(grid2) , 0.0 , numpy.amax(grid2)])
  plt.contour(X,Y,grid2,levels)

  plt.plot([Rslow/rA],[Vslow/vA],'x',ms=20,mew=4)
  plt.plot([Rfast/rA],[Vfast/vA],'x',ms=20,mew=4)

  plt.grid()
  plt.show()
  plt.close()




# =======================================================================
# =======================================================================


def on_click(event):
  global rClick , vClick
  rClick, vClick = event.xdata, event.ydata
  plt.close()

def guess_slow():
  
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  # make grid of Rs and Vs

  RsMin = Rstar*1.05
  RsMax = rA
  VsMin = vA/1000
  VsMax = vA

  nRs = 100
  nVs = 100

  RsAll = numpy.linspace(RsMin,RsMax,nRs)
  VsAll = numpy.linspace(VsMin,VsMax,nVs)

  grid1 = numpy.zeros((nRs,nVs))
  grid2 = numpy.zeros((nRs,nVs))


  # scroll through possibilities
  for iR in range(0,nRs):
    for iV in range(0,nVs):
      
      
      Ar = Fb * VsAll[iV]**0.5 / ( RsAll[iR] * ( 4.0*3.142*Fm )**0.5 ) 
      Vphi = RsAll[iR] * Omega * ( (L*VsAll[iV]**2.0)/(RsAll[iR]**2.0 * Omega) - Ar**2.0 ) / (VsAll[iV]**2.0 - Ar**2.0)
      
      Aphi = Ar * (Vphi - RsAll[iR]*Omega) / VsAll[iV]
      
      grid1[iV][iR] = (VsAll[iV]**2.0 - Ar**2.0) * (2.0*cs**2.0 + Vphi**2.0 - G*Mstar/RsAll[iR]) + 2.0*VsAll[iV]*Vphi*Ar*Aphi
      grid2[iV][iR] = (VsAll[iV]**2.0 - Ar**2.0) * (VsAll[iV]**2.0 - cs**2.0) - VsAll[iV]**2.0 * Aphi**2.0



  bam1 = numpy.where(grid1 > 0)
  bam2 = numpy.where(grid1 < 0)
  grid1[bam1] = 10
  grid1[bam2] = -10

  bam3 = numpy.where(grid2 > 0)
  bam4 = numpy.where(grid2 < 0)
  grid2[bam3] = 10
  grid2[bam4] = -10


  X,Y = numpy.meshgrid(RsAll/rA,VsAll/vA)




  plt.figure()

  levels = numpy.array([ -numpy.amax(grid1) , 0.0 , numpy.amax(grid1) ])
  plt.contour(X,Y,grid1,levels)

  levels = numpy.array([ - numpy.amin(grid2) , 0.0 , numpy.amax(grid2)])
  plt.contour(X,Y,grid2,levels)
  
  plt.grid()


  bam = plt.connect('button_press_event', on_click)
  
  

  plt.show()
  
  
  return (rClick,vClick)


    
# =======================================================================
# =======================================================================

def get_isothermal_Parker(Rmin,Rmax,nGrid):
  
  # get the isothermal parker wind
  
  # get sonic point radius
  Rc = G * Mstar / (2.0 * cs**2.0)
  
  # make grid of R and Vr
  RGrid = numpy.linspace(Rmin,Rmax,nGrid)
  VrGrid = numpy.zeros(nGrid)
  
  for iR in range(0,nGrid):
    
    Vr0 = cs
    dVr = cs/10
    if (RGrid[iR] < Rc):
      dVr = - dVr  
    
    LHS = Vr0 * math.exp(-Vr0**2.0 / (2.0 * cs**2.0) ) 
    RHS = cs * (Rc/RGrid[iR])**2.0 * math.exp(-2.0 * Rc/RGrid[iR] + 3.0/2.0)
    
    Vr = Vr0
    sensitivity = 1.0e-8
    while (abs(LHS/RHS-1.0) > sensitivity):
      
      # save old LHS 
      LHSold = LHS
      
      # update Vr
      Vr = Vr + dVr
      
      # calculate new LHS 
      LHS = Vr * math.exp(-Vr**2.0 / (2.0 * cs**2.0) ) 
      
      # see if crossed solution
      if ( (LHS/RHS-1.0) * (LHSold/RHS-1.0) < 0):
	dVr = -dVr/2.0
      
    VrGrid[iR] = Vr

  return (RGrid,VrGrid)
  
# =======================================================================
# =======================================================================


  
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

# SOME SETUP

# some constants
Msun = 1.99e33 # g
Rsun = 6.955e10 # cm
year = 24.0 * 60.0 * 60.0 * 365.2425 # s
OmegaSun = 2.67e-6 # rad/s
kB = 1.38e-16 # erg/K
gamma = 1.0 # 5.0/3.0
mp = 1.67e-24 # g
AU = 1.496e13 # cm
G = 6.67e-8 # cgs


# ----------------------------------------------------------------------------------------------------
# HERE YOU SETUP THE SIMULATION PARAMETERS 

# physical values of the system
Mstar = 0.5			# stellar mass in Msun
Rstar = 0.5			# stellar equatorial radius in Rsun
Omega = 5.0			# stellar rotation rate in solar value (assumed 2.67e-6 rad/s)
Mdot  = 1.0e-14			# mass loss rate in Msun/year
Twind = 2.0			# wind temperature (MPT: in 1e6 K???)
Bstar = 100.0			# magnetic field strength in G

# values for output of results
filename = "WDanalytic"		# output filename
Rmin = Rstar*1.1		# radius of inner boundary
Rmax = AU*1.0			# radius of outer boundary
nR = 200			# number of grid points in between boundaries 

# some solver values
V0 = 0.005			# starting value of V
dV = 0.005
sensitivity = 1.0e-7		# desired sensitivity (make smaller for more accurate result, but make larger if solver gets stuck)
nClickMax = 10 			# number of times to ask for user input on interations

# ----------------------------------------------------------------------------------------------------







# convert input units to CGS
Mstar = Mstar * Msun
Rstar = Rstar * Rsun
Omega = Omega * OmegaSun
Mdot = Mdot * Msun / year
Twind = Twind * 1.0e6


# need as global variables for clicking on slow solution. Only vary nClickMax!!! 
rClick = 1 # keep constant
vClick = 1 # keep constant
nClick = 0 # keep constant
R0click = 1 # keep constant
V0click = 1 # keep constant





# ----------------------------------------------------------------------------------------------------
# IN THIS SECTION, THE CODE WILL SOLVE THE WD MODEL AS DESCRIBED IN THE APPENDIX OF THE PAPER

# get fluxes
Fm = Mdot / (4.0 * 3.142)
Fb = Bstar * Rstar**2.0

# get sound speed
cs = math.sqrt( gamma * kB * Twind / (0.6 * mp) )



# fit V

V = V0

# get Alfven point
vA = cs / math.sqrt(V)
rA = Fb  / math.sqrt(Fm*4.0*3.142*vA)


# get L (angular momentum)
L = Omega * rA**2.0

Rslow,Vslow,Rfast,Vfast = get_slow_fast(V)


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
U = G * Mstar / (rA * vA**2.0)
V = cs**2.0 / vA**2.0 
W = (rA * Omega)**2.0 / vA**2.0

xSLOW = Rslow / rA
uSLOW = Vslow / vA
constSLOW = W + 0.5 * uSLOW**2.0 + 0.5 * W * xSLOW**2.0 * (1-uSLOW)**2.0 / (1-xSLOW**2.0*uSLOW)**2.0 - V * math.log(xSLOW**2.0) - V * math.log(uSLOW) - U/xSLOW - W * xSLOW**2.0 * (1.0 - uSLOW) / (1.0 - xSLOW**2.0 * uSLOW)

xFAST = Rfast / rA
uFAST = Vfast / vA
constFAST = W + 0.5 * uFAST**2.0 + 0.5 * W * xFAST**2.0 * (1-uFAST)**2.0 / (1-xFAST**2.0*uFAST)**2.0 - V * math.log(xFAST**2.0) - V * math.log(uFAST) - U/xFAST - W * xFAST**2.0 * (1.0 - uFAST) / (1.0 - xFAST**2.0 * uFAST)
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 



while (abs(constSLOW/constFAST-1) > sensitivity):
  
  # save previous values
  Vlast = V
  ratiolast = constSLOW / constFAST
  
  # update V
  V = V + dV
  
  # get Alfven point
  vA = cs / math.sqrt(V)
  rA = Fb  / math.sqrt(Fm*4.0*3.142*vA)
  
  
  # get L (angular momentum)
  L = Omega * rA**2.0

  Rslow,Vslow,Rfast,Vfast = get_slow_fast(V)
  
  
  
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
  U = G * Mstar / (rA * vA**2.0)
  V = cs**2.0 / vA**2.0 
  W = (rA * Omega)**2.0 / vA**2.0

  xSLOW = Rslow / rA
  uSLOW = Vslow / vA
  constSLOW = W + 0.5 * uSLOW**2.0 + 0.5 * W * xSLOW**2.0 * (1-uSLOW)**2.0 / (1-xSLOW**2.0*uSLOW)**2.0 - V * math.log(xSLOW**2.0) - V * math.log(uSLOW) - U/xSLOW - W * xSLOW**2.0 * (1.0 - uSLOW) / (1.0 - xSLOW**2.0 * uSLOW)

  xFAST = Rfast / rA
  uFAST = Vfast / vA
  constFAST = W + 0.5 * uFAST**2.0 + 0.5 * W * xFAST**2.0 * (1-uFAST)**2.0 / (1-xFAST**2.0*uFAST)**2.0 - V * math.log(xFAST**2.0) - V * math.log(uFAST) - U/xFAST - W * xFAST**2.0 * (1.0 - uFAST) / (1.0 - xFAST**2.0 * uFAST)
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 

    
  ratio = constSLOW / constFAST

  print 'Fitting V: ' , V , dV , constSLOW , constFAST , ratio
  
  # test for crossing solution
  if ((ratiolast-1.0)*(ratio-1.0) < 0.0):
    dV = -dV/2.0
  


Rslow,Vslow,Rfast,Vfast = get_slow_fast(V)



# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
# get the final solution parameters
U = G * Mstar / (rA * vA**2.0)
V = cs**2.0 / vA**2.0 
W = (rA * Omega)**2.0 / vA**2.0

xSLOW = Rslow / rA
uSLOW = Vslow / vA
constSLOW = W + 0.5 * uSLOW**2.0 + 0.5 * W * xSLOW**2.0 * (1-uSLOW)**2.0 / (1-xSLOW**2.0*uSLOW)**2.0 - V * math.log(xSLOW**2.0) - V * math.log(uSLOW) - U/xSLOW - W * xSLOW**2.0 * (1.0 - uSLOW) / (1.0 - xSLOW**2.0 * uSLOW)

xFAST = Rfast / rA
uFAST = Vfast / vA
constFAST = W + 0.5 * uFAST**2.0 + 0.5 * W * xFAST**2.0 * (1-uFAST)**2.0 / (1-xFAST**2.0*uFAST)**2.0 - V * math.log(xFAST**2.0) - V * math.log(uFAST) - U/xFAST - W * xFAST**2.0 * (1.0 - uFAST) / (1.0 - xFAST**2.0 * uFAST)


plot_zeros()

const = 0.5 * ( constSLOW + constFAST ) 

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 



# ----------------------------------------------------------------------------------------------------
# Now we have the solution, we should calculate the radial structures of relevant parameters
# These are:-
#	density
#	Vr
#	Vphi
#	Br
#	Bphi


# setup arrays for holding grid radial and Vr values
Rgrid = numpy.logspace(math.log10(Rmin),math.log10(Rmax),nR)
Vrgrid = numpy.zeros(nR)

# also add array for failed grid points. These are simply removed from the final solution
FAILEDgrid = numpy.zeros(nR) # = 1 if no solution found at that point



# find isothermal parker wind solution
Riso , Vriso = get_isothermal_Parker(Rmin,Rmax,nR)

print Vriso
print Riso
# scroll through grid points and get radial outflow speed
for iR in range(0,nR):
  
  R = Rgrid[iR]
  
  if (abs(R/Rslow-1.0) < 0.001): # at slow point
    Vrgrid[iR] = Vslow
    #print 'hernknniose' , Vrgrid[iR] , Rslow
   
  elif (abs(R/rA-1.0) < 0.001): # at Alfven point
    Vrgrid[iR] = rA
    
  elif (abs(R/Rfast-1.0) < 0.001): # at fast point
    Vrgrid[iR] = Vfast 
    
  elif (R < Rslow): # between star and slow point
    Vrgrid[iR] = find_solution(0.01,Vslow)
    #print 'here' , Vrgrid[iR]
               
  elif (R < rA): # between slow point and Alfven point
    Vmin = numpy.amax([Vrgrid[iR-1],Vslow])
    Vrgrid[iR] = find_solution(Vmin,vA)
     
  elif (R < Rfast): # between Alfven point and fast point
    Vrgrid[iR] = find_solution(Vrgrid[iR-1],Vfast)
     
  else: # outside of fast point
    Vrgrid[iR] = find_solution(Vrgrid[iR-1],Vfast*5.0)
  
  if (Vrgrid[iR] == -1):
    
    FAILEDgrid[iR] = 1
    Vrgrid[iR] = Vrgrid[iR-1]
    


# get rid of failed grid points
Rgrid = Rgrid[numpy.where(FAILEDgrid != 1)]
Vrgrid = Vrgrid[numpy.where(FAILEDgrid != 1)]


# get radial component of velocity
Argrid = Fb * Vrgrid**0.5 / ( Rgrid * ( 4.0*3.142*Fm )**0.5 ) 
Vphigrid = Rgrid * Omega * ((L * Vrgrid**2.0)/(Rgrid**2.0 * Omega) - Argrid**2.0) / (Vrgrid**2.0 - Argrid**2.0)



# get Michel velocity
Vm = ( Rstar**4.0 * Bstar**2.0 * Omega**2.0 / Mdot )**(1.0/3.0)
#print Vm


# ----------------------------------------------------------------------------------------------------
# Now output the results into a plot, a pickle file, and a ascii text file

# output plot

plt.figure()


plt.plot(Rgrid/AU,Vrgrid/1.0e5,'b',label='Weber-Davis')
plt.plot(Rgrid/AU,Vphigrid/1.0e5,'b--',label='Weber-Davis Vphi')

plt.plot(Riso/AU,Vriso/1.0e5,'k:',label='Parker solution')


plt.axvline(x=Rslow/AU, color='r', linestyle=':',label='Slow radius')
plt.axvline(x=rA/AU, color='g', linestyle=':',label='Alfven radius')
plt.axvline(x=Rfast/AU, color='y', linestyle=':',label='Fast radius')

plt.axhline(y=Vm/1.0e5, color='k', linestyle=':')


plt.grid()

#plt.xlim([Rstar*1.1/AU,1.0])
plt.xlim([Rstar*1.1/AU,0.2])
plt.legend(loc=4)

plt.xlabel("Distance from star (AU)")
plt.ylabel("Radial outflow speed (km s$^{-1}$)")

plt.savefig(filename + ".png")
plt.savefig(filename + ".pdf")



# save data into pickle file
with open(filename+'.pickle', 'w') as f:  
    pickle.dump([Rgrid,Vrgrid,Vphigrid,Rslow,Vslow,rA,vA,Rfast,Vfast,Riso,Vriso,Vm], f)


# use this code to read in the results from the pickle file
#with open(filename+'.pickle') as f:  
#    Rgrid,Vrgrid,Vphigrid,Rslow,Vslow,rA,vA,Rfast,Vfast,Riso,Vriso,Vm = pickle.load(f)






# output to ascii file
filename = filename + '.dat'
outfile = open(filename,'w')

# write header with variables and units
outfile.write('Alfven radius = ' + str(rA/AU) + ' AU \n')
outfile.write('Outflow speed at Alfven point = ' + str(vA/1.0e5) + ' km/s \n')
outfile.write('Slow point radius = ' + str(Rslow/AU) + ' AU \n')
outfile.write('Outflow speed at slow point = ' + str(Vslow/1.0e5) + ' km/s \n')
outfile.write('Fast point radius = ' + str(Rfast/AU) + ' AU \n')
outfile.write('Outflow speed at fast point = ' + str(Vfast/1.0e5) + ' km/s \n')
header = "R (cm)    density (g/cm^3)    Vr (cm/s)    Vphi (cm/s)     Br (G)     Bphi (G) \n"
outfile.write(header)

for i in range(0,len(Rgrid)):
  
  # get values in desired units
  R = Rgrid[i] # in cm
  Vr = Vrgrid[i] # in cm/s
  density = Mdot / (4.0 * 3.142 * R**2.0 * Vr) # in g/cm^3
  Vphi = Vphigrid[i]
  Br = Fb / R**2.0
  Bphi = Br * ( (Vphi - R * Omega) / (Vr) )
  
  line = str(R) + '    ' + str(density) + '    ' + str(Vr) + '    ' + str(Vphi) + '    ' + str(Br) + '    ' + str(Bphi) + ' \n'
  outfile.write(line)


outfile.close()




