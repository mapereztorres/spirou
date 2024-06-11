import pandas as pd
from astropy import units as u
from astropy.coordinates import SkyCoord

def ang_separation(filename='./INPUT/SPI-sources.csv', ind1=0, ind2=33):
    """Compute the angular separation, in degrees, for two sources in a table

    Parameters
    ----------
    filename : string (optional)
    ind1     : index of 1st source (optional)
    ind2     : index of 2nd source (optional)

    Returns
    -------
    sep : angular separation between those two sources, in degrees

    Test:

    ang_separation() computes, by default, the 
    angular separation between Barnard's star and GJ 1214 - 10 deg 35 arcmin
    #ra1 = data['RA(h)'][0]; dec1 = data['DEC(deg)'][0]
    #ra2 = data['RA(h)'][33]; dec2 = data['DEC(deg)'][33]; 
    """
    data = pd.read_csv(filename)
    ra1 = data['RA(h)'][ind1]; dec1 = data['DEC(deg)'][ind1]
    ra2 = data['RA(h)'][ind2]; dec2 = data['DEC(deg)'][ind2]; 

    c1 = SkyCoord(ra1, dec1, unit=(u.hourangle, u.deg))
    c2 = SkyCoord(ra2, dec2, unit=(u.hourangle, u.deg))

    sep = c1.separation(c2)
    return sep
