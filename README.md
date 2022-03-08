# XtraCrysPy
A crystallographic visualizer for Python.

Features:
- Plot molecular systems with minimal setup
- Plot systems from DFT input- or output- files
- Interact with the system to caluclate distances, angles, or receive atomic information
- Display relaxation steps from QE relax output files
- Display reciprocal space features
- Plot iso-surfaces (in development)

Fury currently has a bug with the cylinder routine, used for displaying bonds.
Clone the forked fury reposity, from my GitHub, and check out the glyph\_orientation branch
(https://github.com/Sassafrass6/fury/tree/glyph_orientation)

Requirements:  
- python 3.7  
- numpy 1.17  
- [fury 0.8](https://github.com/fury-gl/fury)
  
Installation (only use --user if you do not have permission to install python packages):  
- `python setup.py install`  
- `python setup.py install --user`  
  
Usage:  
  See examples
