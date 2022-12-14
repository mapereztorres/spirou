# JupyterSPI

Code to estimate radio emission from star-planet (Sub-Alfvénic) interaction
via the electron-cyclotron mechanism.  Initially, it was a Jupyter notebook
that grew too much, hence it is now split into several pieces of python code.

## Files

```
├── SPI.ipynb - Jupyter notebook for the code - outdated
├── SPI.py    - Current version, in continuous development
├── SPI.py.orig - original version of SPI.py.orig - to be removed eventually
├── SPIutils.py - Utility functions for SPI.py, so SPI.py is not too large 
```

## Directories

```
.
├── INPUT     - Contains input files (data tables to be fed into the code)
├── OUTPUT    - Contains output files 
├── OUTPUT.old
├── OUTPUT.veryold
├── pics     - pictures to use
├── README.md - This file
├── specific  - Jupyter nb and code specific for some sources 
├── SPIworkflow -  Folder containing __init__.py and constants.py 
└── testdir  - Folder with code to test functions, subroutines, etc.
```

The main piece of code is SPI.py.
SPIutils.py contains functions used by SPI.py, to keep the main program not
too large. 
Definitions and constants used in the code are under SPIworkflow, in files
``__init__.py`` and ``constants.py``. 



## Limitations 

The current code uses the approximation of an isothermal Parker wind for the
stellar wind of all stars. Therefore, stellar rotation and magnetic field
effects on the wind are neglected.  


## List of things to be developed/modified/etc. 

Obtain the stellar wind number density at the base of the corona,
``n_sw_base``, from an assumed/measured/estimated mass loss rate for the star.
Currently, ``n_sw_base`` is read from __init__.py.
