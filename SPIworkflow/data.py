import os
import pandas as pd


def get_spi_data(infile_data='./INPUT/SPI-targets.csv',
        out_latex=True, outfile='latex_table.tex', 
        distance_min = 0.1, distance_max=1000.0, 
        p_orb_min=0.1, p_orb_max=1000.0, 
        bfield_min=10., bfield_max=10000., 
        dec_min = -90.0, dec_max=90.): 
    """
    Read in the input data to estimate radio emission from SPI

    Returns
    -------
    data : pandas.Dataframe
        The SPI useful data, after having applied different masks

    latex_table.tex: Latex table 
        Latex table of the data pandas.Dataframe (optional)

    Parameters
    ----------
    infile_data : string (optional)
        location to read in the data

    out_latex : Boolean (optional)
        If True, create table

    outfile   : string (optional)
        Name of the output latex table

    distance_min, distance_max : float (optional)
        Distance cuts (min and max), in pc

    p_orb_min, p_orb_max : float (optional)
        Orbital period cuts, in days

    bfield_min, bfield_max : float (optional)
        Magnetic field cuts, in Gauss

    dec_min, dec_max : float (optional)
        Declination cuts, in degrees
    """
    #df=df.rename(columns={"star_radius(R_Sun)": "radius_star(r_sun)", "Planet_radius(R_Earth)": "radius_planet(r_earth)"
    #                  , "P_orb(days)": "p_orb(days)", "M_star": "mass_star(m_sun)"
    #                   , "semi_major_axis(AU)": "a(au)", "RA(deg)": "ra(deg)"
    #                   , "DEC(deg)": "dec(deg)", "Planet_mass(M_Earth)": "mass_planet(m_earth)"
    #                   , "starname": "star_name", "P_rot(days)": "p_rot(days)"
    #                   , "<B>": "bfield_star(gauss)", "sptype": "SP_TYPE"
    #                   , "distance(pc)": "d_star(pc)"
    #                  })

    # Read in data file using pandas
    df = pd.read_csv(infile_data)

    # Create new pandas.Dataframe with a subset of rows (easy to display in,
    # e.g., jupyter notebook.
    if infile_data =='./INPUT/SPI-sources_planets_MASTER.csv':
        df2 = df[["star_name", "SP_TYPE", "d_star(pc)", "mass_star(m_sun)",
            "radius_star(r_sun)", "p_rot(days)", "bfield_star(gauss)",
            "planet_name", "p_orb(days)", "a(au)", "radius_planet(r_earth)",
            "mass_planet(m_earth)","ra(hms)","dec(dms)", "ra(deg)", "dec(deg)"]]

    elif infile_data =='./INPUT/my-SPI-sources.csv':
        df2 = df[["planet_name", "star_name", 
            "ra(hms)","dec(dms)", "ra(deg)", "dec(deg)", 
            "d_star(pc)", "radius_star(r_sun)", "mass_star(m_sun)", "p_rot(days)", "bfield_star(gauss)",
            "radius_planet(r_earth)", "mass_planet(m_earth)", "p_orb(days)", "a(au)" ]]
    else: 
        df2 = df[["planet_name", "star_name",  
            "d_star(pc)", "radius_star(r_sun)", "mass_star(m_sun)", "p_rot(days)", "bfield_star(gauss)", 
            "e_bfield_star", "radius_planet(r_earth)", "mass_planet(m_earth)", "p_orb(days)", "a(au)" ]]

    #
    #df2.to_csv(r'./INPUT/SPI-sources.csv', index=True, header=True)

    # Define masks
    #
    #mask_Rpl = df2['radius_planet(r_earth)'] < 1.5
    #p_rot_mask = df2['p_rot(days)'] < 200.0
    distance_mask_min =  df2['d_star(pc)'] > distance_min
    distance_mask_max =  df2['d_star(pc)'] < distance_max 
    p_orb_mask_min = df2['p_orb(days)'] > p_orb_min
    p_orb_mask_max = df2['p_orb(days)'] < p_orb_max
    bfield_mask_min = df2['bfield_star(gauss)'] > bfield_min
    bfield_mask_max = df2['bfield_star(gauss)'] < bfield_max
    #
    #MPT@20230116 -  The masks in declination give a TypeError when reading
    # './INPUT/SPI-sources.csv' 
    # (but not when reading the file './INPUT/SPI-sources_planets_MASTER.csv'    
    # TypeError: '>' not supported between instances of 'str' and 'float'
    # Needs to be understood and fixed.
    # It may have to do with "NaN" values, which may conflict
    #
    #declination_mask_min = df2['dec(deg)'] > dec_min
    #declination_mask_max = df2['dec(deg)'] < dec_max

    # Apply masks
    #data = df2[distance_mask_min & distance_mask_max & p_orb_mask_min &
    #        p_orb_mask_max & bfield_mask_min & bfield_mask_max &
    #        declination_mask_min & declination_mask_max]
    data = df2[distance_mask_min & distance_mask_max & p_orb_mask_min &
            p_orb_mask_max & bfield_mask_min & bfield_mask_max]


    # Create column with cyclotron freq, in GHz
    # It gives a SetWarningCopy message, but it's seems to work fine.
    #
    #data['freq_cycl(ghz)'] = data.loc[:,('bfield_star(gauss)')]*2.8e-3
    data['freq_cycl(ghz)'] = data['bfield_star(gauss)']*2.8e-3

    data.reset_index(inplace=True)


    # Generate also latex_table? out_latex==True => Return. It will truncate 
    if out_latex==True:
        with open('./OUTPUT/latex_table.tex', 'w') as tf:
            tf.write(data.to_latex(index=False))
        #print(data.to_latex(index=False))

    return data

def create_data_tables(infile_data='./INPUT/SPI-targets.csv'):
    """
    Read in an input table with stars and, optional, planets, and
    creates two tables: one including targets hosting planets, and another
    one with targets that do not host planets.

    Returns
    -------
    df_planets : pandas.Dataframe
         It contains the targets with planets 

    df_no_planets : pandas.Dataframe
         It contains the targets without planets 

    outfile_planets: File, for the moment, it's fixed:  'StarTable-CARMENES_only_planets.csv')

    outfile_no_planets: File, for the moment, it's fixed:  'StarTable-CARMENES_no_planets.csv')

    Parameters
    ----------
    infile_data : string (optional)
        location to read in the data

    """

    # Create OUTPUT directory for keeping output from code.
    outdir = 'OUTPUT'
    try:
        os.mkdir(outdir)
        print('Directory', outdir, 'created.\n')
    except FileExistsError:
        print('Directory', outdir, 'already exists.\n')

    df = pd.read_csv(infile_data)
    # Copy Dataframes to generate two tables: one for sources with planets,
    # and a second one for sources not hosting planets
    df_planets=df.copy()
    df_planets=df_planets.dropna(subset=['planet_name'], inplace=False)

    df_no_planets = df.copy()
    df_no_planets = df_no_planets[df_no_planets['planet_name'].isnull()]

    # Create CARMENES tables for targets with planets only, and with no planets,
    # unless they already exist
    # 
    outfile_planets = os.path.join(outdir, 'StarTable-CARMENES_only_planets.csv')
    outfile_no_planets = os.path.join(outdir, 'StarTable-CARMENES_no_planets.csv')

    if os.path.exists(outfile_planets): 
        print(f'File {outfile_planets} already exists. I will not overwrite it.\n')
    else:
        print(f'Creating table file: {outfile_planets}.\n')
        df_planets.to_csv(outfile_planets)

    if os.path.exists(outfile_no_planets): 
        print(f'File {outfile_no_planets} already exists. I will not overwrite it.\n')
    else:
        print(f'Creating table file: {outfile_no_planets}.\n')
        df_no_planets.to_csv(outfile_no_planets)

    return outdir, df_planets, df_no_planets 
