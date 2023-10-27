from SPIworkflow.constants import *

# class to handle output to file
class OutputWriter:
    def __init__(self, outfileTXT):
        self.outfileTXT = outfileTXT
    
    def write_parameters(self, T_corona, M_star_dot, mu, d, R_star, M_star, P_rot_star, B_star,
        Exoplanet, Rp, Mp, r_orb, P_orb, 
        loc_pl, n_base_corona, nu_plasma_corona, gyrofreq,
        Flux_r_S_min, Flux_r_S_max, rho_sw, n_sw_planet, v_sw_base, Flux_r_S_ZL_min,
        Flux_r_S_ZL_max, v_sw, v_rel, v_alf, M_A, B_sw, Rmp, Rp_eff):

        with open(self.outfileTXT, 'w') as f:
             f.write('# INPUT PARAMETERS:               ########\n') 
             f.write('#                                 ########\n') 
             f.write('# GENERIC (fixed for all targets) ########\n') 
             f.write('#                                 ########\n') 
             f.write('T_corona                 = {0:.0e} K\n'.format(T_corona))
             f.write('Stellar wind mass loss rate = {0:.3f} Sun mass loss rate\n'.format(M_star_dot))
             f.write('mean molecular weight    = {0:.2f}\n'.format(mu))
             #f.write('mu_0 [mks units] = {0:.0e}\n \n'.format(mu_0))
             f.write('#                                 ########\n') 
             f.write('# STAR:                           ########\n') 
             f.write('#                                 ########\n') 
             #f.write('B_star = {0:.0f} G; B_planet = {1:.0f} G\n'.format(B_star, Bp[loc_pl]))
             f.write('Distance to star         = {0:.2f} pc\n'.format(d/pc))
             f.write('Stellar radius           = {0:.3f} Rsun\n'.format(R_star/R_sun))
             f.write('Stellar mass             = {0:.3f} Msun\n'.format(M_star/M_sun))
             f.write('Stellar rotation period  = {0:.3f} days\n'.format(P_rot_star/day))
             f.write('Stellar magnetic field   = {0:.2f} Gauss\n'.format(B_star))
             f.write('#                                 ########\n') 
             f.write('# PLANET:                         ########\n') 
             f.write('#                                 ########\n') 
             f.write('Exoplanet name           = {0}\n'.format(Exoplanet)) 
             f.write('Planetary radius         = {0:.3f} R_Earth\n'.format(Rp/R_earth))
             f.write('Planetary mass           = {0:.3f} M_Earth\n'.format(Mp/M_earth))
             f.write('Semimajor axis of orbit  = {0:.3f} au\n'.format(r_orb/au))
             f.write('Orbital period           = {0:.3f} days\n'.format(P_orb))
             #f.write('Stellar wind particle density at the base = {0:.0e} cm-3\n'.format(n_sw_base))
             f.write('#                                 ########\n') 
             f.write('#                                 ########\n') 
             f.write('# OUTPUT PARAMETERS:              ########\n')
             f.write('#                                 ########\n') 
             f.write('n_base_corona = {0:.3e} cm^-3\n'.format(n_base_corona))
             f.write('nu_plasma_corona = {0:.2e} MHz\n'.format(nu_plasma_corona/1e6))
             f.write('ECMI freq (fundamental) = {0:.0f} MHz\n'.format(gyrofreq/1e6))
             f.write('Flux_ST: ({0}, {1}) mJy\n'.format(Flux_r_S_min[loc_pl], Flux_r_S_max[loc_pl]))
             #f.write('Flux_ZL: ({0}, {1}) mJy\n'.format(Flux_r_S_ZL_min[loc_pl], Flux_r_S_ZL_max[loc_pl]))
             f.write('rho_sw at r_orb: {0} \n'.format(rho_sw[loc_pl]))
             f.write('n_sw_planet at r_orb: {0} \n'.format(n_sw_planet[loc_pl]))
             f.write('n_sw_planet at base: {0} \n'.format(n_sw_planet[0]))
             f.write('v_sw at base of the wind: {0} \n'.format(v_sw_base))
             f.write('Flux_ZL: ({0}, {1}) mJy\n'.format(Flux_r_S_ZL_min[loc_pl], Flux_r_S_ZL_max[loc_pl]))
             f.write('v_sw at the base of the wind: {0} \n'.format(v_sw_base))
             f.write('v_sw at r_orb: {0} \n'.format(v_sw[loc_pl]))
             f.write('v_sw(r_orb)/v_sw_base: {0} \n'.format(v_sw[loc_pl]/v_sw_base))
             f.write('v_rel at r_orb: {0} \n'.format(v_rel[loc_pl]))
             f.write('v_alf at r_orb: {0} \n'.format(v_alf[loc_pl]))
             f.write('M_A at r_orb: {0} \n'.format(M_A[loc_pl]))
             #f.write('v_rel/v_alf at r_orb: {0} \n'.format(v_rel[loc_pl]/v_alf[loc_pl]))
             f.write('B_sw at r_orb: {0} \n'.format(B_sw[loc_pl]))
             f.write('Rmp at r_orb: {0} \n'.format(Rmp[loc_pl]/Rp))
             f.write('Rp_eff at r_orb: {0} \n'.format(Rp_eff[loc_pl]/Rp))

