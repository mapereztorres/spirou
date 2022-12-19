import os
import pandas as pd


def get_spi_data(infile_data='./INPUT/SPI-sources_planets_MASTER.csv',
        out_latex=True, outfile='latex_table.tex', 
        distance_min = 0.1, distance_max=15.0, 
        p_orb_min=0.1, p_orb_max=10.0, 
        bfield_min=10., bfield_max=180., 
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

    # Read in data file
    df = pd.read_csv(infile_data)

    # Create new pandas.Dataframe with a subset of rows (easy to display in,
    # e.g., jupyter notebook
    df2 = df[["star_name", "SP_TYPE", "d_star(pc)", "mass_star(m_sun)",
        "radius_star(r_sun)", "p_rot(days)", "bfield_star(gauss)",
        "planet_name", "p_orb(days)", "a(au)", "radius_planet(r_earth)",
        "mass_planet(m_earth)","ra(hms)","dec(dms)", "ra(deg)", "dec(deg)"]]

    # Apply masks and generate new pd Dataframe

    #mask_Rpl = df2['radius_planet(r_earth)'] < 1.5
    #p_rot_mask = df2['p_rot(days)'] < 200.0
    distance_mask_min =  df2['d_star(pc)'] > distance_min
    distance_mask_max =  df2['d_star(pc)'] < distance_max 
    p_orb_mask_min = df2['p_orb(days)'] > p_orb_min
    p_orb_mask_max = df2['p_orb(days)'] < p_orb_max
    bfield_mask_min = df2['bfield_star(gauss)'] > bfield_min
    bfield_mask_max = df2['bfield_star(gauss)'] < bfield_max
    declination_mask_min = df2['dec(deg)'] > dec_min
    declination_mask_max = df2['dec(deg)'] < dec_max
    df3 = df2[distance_mask_min & distance_mask_max & p_orb_mask_min &
            p_orb_mask_max & bfield_mask_min & bfield_mask_max &
            declination_mask_min & declination_mask_max]
    df3.reset_index(inplace=True)

    # Create pandas Dataframe to be returned, with subset of useful rows
    data = df3[["star_name", "ra(hms)", "dec(dms)", "SP_TYPE", "d_star(pc)", "p_rot(days)", "bfield_star(gauss)",
          "planet_name", "p_orb(days)", "radius_planet(r_earth)", "mass_planet(m_earth)"]]

    # Generate also latex_table? out_latex==True => Return. It will truncate 
    if out_latex==True:
        with open('./OUTPUT/latex_table.tex', 'w') as tf:
            tf.write(data.to_latex(index=False))
        #print(data.to_latex(index=False))

    return data

def create_data_tables(infile_data='./INPUT/SPI-sources_planets_MASTER.csv',
        outdir='OUTPUT'):
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

    return df_planets, df_no_planets 
