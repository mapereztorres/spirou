# SPIROU

**S**tar-**P**lanet **I**nteraction and **R**adio **O**bservations **U**nited

Code to predict the radio emission from star-planet (Sub-Alfvénic) interaction
via the electron-cyclotron mechanism.  


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
├── OUTPUT    - Contains output files (it's created at the first run of the
code)
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

