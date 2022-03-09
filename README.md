# XtraCrysPy
A crystallographic visualizer for Python.

Features:
- Plot molecular systems with minimal setup
- Plot systems from DFT input- or output- files
- Interact with the system to caluclate distances, angles, or receive atomic information
- Display relaxation steps from QE relax output files
- Display reciprocal space features
- Plot iso-surfaces (in development)
- Band structure and DoS plots from QE or PAOFLOW data
- Electron Transport plots from PAOFLOW data

Fury currently has a bug with the cylinder routine, used for displaying bonds.
Clone the forked fury reposity, from my GitHub, and check out the glyph\_orientation branch
(https://github.com/Sassafrass6/fury/tree/glyph\_orientation)

Requirements:
- python 3.8
- numpy 1.19
- matplotlib 3.5
- [fury 0.8](https://github.com/fury-gl/fury)
  
Installation (only use --user if you do not have permission to install python packages):  
- `python setup.py install`  
- `python setup.py install --user`  
  
Usage:

  See examples/CrystalView for crystal and Brillouin zone plotting:

    From scratch:
      python caffeine.py
      python example1.py

    From inputfile:
      python example2.py
      python example3.py

    Relax QE outputfile:
      python example4.py


  See examples/PlotTools for QE and PAOFLOW plotting functions:
    Close windows to advance plot script

    python example01_qe.py
    python example01_pao.py
