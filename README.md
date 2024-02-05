# SPI

Code to estimate radio emission from star-planet (Sub-Alfvénic) interaction
via the electron-cyclotron mechanism.  Initially, it was a Jupyter notebook
that grew too much, hence it is now split into several pieces of python code.


## Files

```
├── spi.yml     - yml file to create a python environment to run the code
├── SPI.ipynb   - Jupyter notebook for the code - outdated
├── SPI.py      - Current version, in continuous development
├── SPI.py.orig - original version of SPI.py.orig - to be removed eventually
├── SPI_mdot.py      - Current version, in continuous development. Calculates emission as a function of the stellar Ṁ.
├── SPI_Bp.py      - Current version, in continuous development. Calculates emission as a function of the magnetic field of the planet.
```


'''
First, create the conda spi environment:

`conda env create -f spi.yml' 

This will create the conda environment `spi'. 

To activate this environment, use `$ conda activate spi'

To deactivate an active environment, use `$ conda deactivate' 

In case you need to update the evnironment, run 

`conda activate spi
conda env update --name spi --file spi.yml --prune'
'''

## Directories

```
.
├── INPUT     - Contains input files (data tables to be fed into the code)
├── OUTPUT    - Contains output files 
├── OUTPUT.old
├── OUTPUT.veryold
├── pics     - Folder containing pictures 
├── README.md - This file
├── specific  - Folder containing Jupyter nb and code specific for some sources 
├── SPIworkflow - Folder containing the python package 
    ──  __init__.py  - necessary for the folder to become a Python package
    ── constants.py  - useful constants for the code
    ── data.py       - Auxiliary function to read in data to feed modelling
    ── SPIutils.py   - useful functions for SPI.py
    ── <Other_files> - currently not used by the code
└── testdir  - Folder with code to test functions, subroutines, etc.
```

The main piece of code is SPI.py.

SPIworkflow is actually a Python package in itself. It requires at least
the existence of the file ``__init__.py``. Eveything in that file becomes
available once the package is imported. 



## Limitations 

The current code uses the approximation of an isothermal Parker wind for the
stellar wind of all stars. Therefore, stellar rotation and magnetic field
effects on the wind are neglected.  


## List of things to be developed/modified/etc. 

Obtain the stellar wind number density at the base of the corona,
``n_sw_base``, from an assumed/measured/estimated mass loss rate for the star.
Currently, ``n_sw_base`` is read from __init__.py.

Generate plots of Flux density as a function of the planetary magnetic field, for a given
orbital distance (and the rest of the parameters being fixed).

Generate plots of Flux density as a function of the density at the corona of
the star, for a given orbital distance (and the rest of the parameters being
fixed).


