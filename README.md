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
- CIF
- XSF
- CP2K (coming soon)
- BXSF (coming soon)
- Any ASE supported format

Requirements:
- python 3.8
- numpy 1.19
- [fury 0.8](https://github.com/fury-gl/fury)
- [ase](https://wiki.fysik.dtu.dk/ase/)
  
Installation (only use --user if you do not have permission to install python packages):  
- `python setup.py install`  
- `python setup.py install --user`  
  
Usage:

  Display an atomic model from inputfile, with the examples/main.py script:
    python main.py <inputfile>

  Control Inputs:
    'u' : Toggle UI visibility
    'a' : Toggle Axis visibility
    'b' : Toggle Boundary visiblity
    '>' : Step forward in relax or MD
    '<' : Step backward in relax or MD
    CTRL + ('>' or '<') : Step 5% through the relaxation or MD steps
    SHIFT + 'c' : Reset camera to default position
    SHIFT + 's' : Take snapshot
    CTRL + 'w' : Exit
    Arrow Keys : Rotate model
    SHIFT + Arrow Keys : Translate camera
    CTRL + Arrow Keys (+ SHIFT) : Rotate (or Translate) in smaller steps


See examples directory for crystal and Brillouin zone plotting:

  From scratch:
    python caffeine.py
    python silicon_fcc.py

  From inputfile (with main.py script):
    python main.py <inputfile>
    python main.py <inputfile> <bond_length>

  From inputfile (programmatic):
    python SnTe_2D.py
    python example_cif.py
    python example_qe_in.py
    python example_poscar.py

  Relax QE outputfile:
    python example_relax.py

  LAMMPS MD trajectory:
    python example_lammps_md.py

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

