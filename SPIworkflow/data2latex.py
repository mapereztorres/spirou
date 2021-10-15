import pandas as pd

def get_latex_table(data='./INPUT/SPI-sources.csv'):
    df = pd.DataFrame(data, columns = ['planet name', 'RA (h:m:s)', 'DEC (d:m:s)', 'd_pc', 'P_rot', 'planet_radius (R_Earth)', 'orbital_period(days)'])
    for indi in [0, 5, 31, 33, 89, 90, 91]:
        with open('sources.tex', 'w') as tf:
            tf.write(df.to_latex(index=False))
