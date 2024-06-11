import numpy as np

def B_shulyak(P_rot):
    """ Estimates the magnetic field (in Gauss), in case measurements are not available
        based on a fit to the data in Shulyak+2019
        Fit by Luis Peña : log10(B) = m*log10(P(d))+n, where m= -0.2964946637065301 +-
        0.05044615665184572 and n= 3.5429917300425564 +- 0.025728405404311927
      
      INPUT: P_rot - float, in days 
      OUTPUT: Magnetic field estimate, in Gauss
    """
    m = -0.29649; n =  3.54299
    B_field = 10**(m * np.log10(P_rot) + n)
    return B_field
