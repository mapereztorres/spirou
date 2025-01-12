![**S**tar-**P**lanet **I**nteraction and **R**adio **O**bservations **U**nited](pics/spirou-logo.png)


# SPIROU

**S**tar-**P**lanet **I**nteraction and **R**adio **O**bservations **U**nited

`SPIROU` is a public Python code to predict the radio emission from
Sub-Alfvénic star-planet interaction (SPI), via the electron-cyclotron maser
instability mechanism.  The main script is `spirou.py`.

which computes the
expected flux density due to this interaction. 

The code can compute the expected flux density from SPI (1) as a function of the orbital
distance (i.e., the semi-major axis) to the planet,  (2) as a function of the planetary magnetic field, 
or (3) as a function of the stellar wind mass loss rate.  For cases 2 and 3, 

If the semi-major axis is
fixed, tthe code can compute the expected flux density as 

The current version of the code assumes an isothermal Parker wind, and takes
into account the free-free extinction of the radio emission within the stellar
wind.  The code uses two different geometries of the stellar wind magnetic field: a
closed dipolar one, or an open Parker spiral. The effective radius of the
planet, i.e., the magnetopause distance, is obtained as an equilibrium of
pressures between the wind and the planet. 





## Files

```
├── spirou.yml     - yml file to create a python environment to run the code
├── spirou.py      - The main file
```


First, create the conda spirou environment:
 
```
conda env create -f spirou.yml 
```

This will create the conda environment `spirou`. 
To activate this environment, use 

```
$ conda activate spirou
```

To deactivate an active environment, use 

```
$ conda deactivate
```

In case you need to update the evnironment, run 

```
conda activate spirou 
conda env update --name spirou --file spirou.yml --prune
```

## Directories

```
├── INPUT     - Contains input files (data tables to be fed into the code)
├── OUTPUT    - Contains output files (created at first run of code)
├── pics      - Folder containing pictures 
├── README.md - This file
├── output.py - module containing function for output table 
├── SPIworkflow - Folder containing the python package 
    ──  __init__.py  - Initialition file for the code 
    ── constants.py  - File containing useful constants 
    ── load_data.py  - File containing auxiliary functions to read in data to feed modelling
    ── SPIutils.py   - File containing many useful functions for spirou.py
```

SPIworkflow is actually a Python package in itself. It requires at least the
existence of the file ``__init__.py``. Everything in that file becomes
available once the package is imported. 

