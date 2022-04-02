def n_wind(M_star_dot=3e-14, r0=7e10, v_r0=25.6e5):
    """ Computes the particle density of the stellar wind at some distance d from the
        center of the star.
        OUTPUT: n_sw - particle density of the stellar wind at distance d, in #/cm^3
        INPUT:  M_star_dot - stellar mass-loss rate, in M_sun/yr
                r0         - radial distance from the center of the star, in cm
                v_r0       - Speed of stellar wind at the reference radial distance r0, in cm/s
    """
    M_star_dot *= M_sun/yr2sec  # mass-loss rate, in grams/sec
    m_av  =  1.92e-24 # average particle mass of the solar wind, in grams
    rho = M_star_dot/(4*np.pi * r0**2 * v_r0*1e5) # Density of the stellar wind, in gr/cm3
    n_sw = rho / m_av
    return n_sw

M_sun_dot = 2e-14 # Mass-loss rate of Sun, in M_sol/yr
# Proxima Cen
#M_star_dot = 0.035*M_sun_dot
M_star_dot = M_sun_dot
#n_w = n_wind(M_star_dot, d=0.145*R_sun, v_sw=18.1)
n_w = n_wind(M_star_dot, r0=R_sun, v_r0=1.0)
print("Particle density is n_sw = {0:.2e} #/cm^3".format(n_w))

