# SPI

Code to estimate radio emission from star-planet (Sub-Alfvénic) interaction
via the electron-cyclotron mechanism.  Initially, it was a Jupyter notebook
that grew too much, hence it is now split into several pieces of python code.


## Files

```
├── spi.yml     - yml file to create a python environment to run the code
├── SPI.py      - Current version, in continuous development
├── SPI_mdot.py      - Current version, in continuous development. Calculates emission as a function of the stellar Ṁ. - Not used anymore
├── SPI_Bp.py      - Current version, in continuous development. Calculates emission as a function of the magnetic field of the planet. - Not used anymore
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
├── pics      - Folder containing pictures 
├── README.md - This file
├── specific  - Folder containing Jupyter nb and code specific for some sources 
├── output.py - module containing function for output (table and, eventually, graphs)
├── SPIworkflow - Folder containing the python package 
    ──  __init__.py  - File needed for the folder to become a Python package
    ── constants.py  - File containing useful constants for the code
    ── load_data.py  - File containing auxiliary functions to read in data to feed modelling
    ── SPIutils.py   - File containing useful functions for SPI.py
└── testdir          - Folder with code to test functions, subroutines, etc.
```

The main piece of code is SPI.py.

SPIworkflow is actually a Python package in itself. It requires at least
the existence of the file ``__init__.py``. Everything in that file becomes
available once the package is imported. 


## Limitations 

The current code uses the approximation of an isothermal Parker wind for the
stellar wind of all stars. Therefore, stellar rotation and magnetic field
effects on the wind are neglected.  


## List of things to be developed/modified/etc. 

