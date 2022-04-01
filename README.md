# XtraCrysPy
A Python tool for visualizing atomic systems and properties of condensed matter.

Features:
- Plot molecular systems with minimal setup
- Plot systems from DFT input- or output- files
- Interact with the system to caluclate distances, angles, or receive atomic information
- Display relaxation steps from QE relax output files
- Display reciprocal space features
- Plot real space iso-surfaces (charge density, etc)
- Plot reciprocal iso-surfaces (Fermi surface, spin-texture, etc)
- Band structure and DoS plots from QE or PAOFLOW data
- Electron Transport plots from PAOFLOW data

Supported File Types:
- Quantum ESPRESSO input or output
- POSCAR
- CIF (in progress)
- XSF
- BXSF (coming soon)

Fury currently has a bug with the cylinder routine, used for displaying bonds.
Clone the forked fury reposity, from my GitHub, and check out the glyph\_orientation branch
(https://github.com/Sassafrass6/fury/tree/glyph_orientation)

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
      python silicon_fcc.py

    From inputfile (with main.py script):
      python main.py <inputfile>
      python main.py <inputfile> <bond_length>

    From inputfile (programmatic):
      python SnTe_2D.py
      python example_qe_in.py
      python example_poscar.py

    Relax QE outputfile:
      python example_relax.py

    Charge Density:
      python example_charge_density.py

    Reciprocal Space:
      python example_BZ.py

    Fermi Surfaces:
      python example_fermi_surface1.py
      python example_fermi_surface2.py
      python example_fermi_surface_clip.py
      python example_colored_fermi_surface.py
      python example_textured_fermi_surface_clip.py

  See examples/PlotTools for QE and PAOFLOW plotting functions:
    Close windows to advance plot script

    python example01_qe.py
    python example01_pao.py
