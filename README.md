# XtraCrysPy
## A Python tool for visualizing atomic systems and spatial properties of condensed matter.

![Alt text](https://github.com/Sassafrass6/XtraCrysPy/blob/master/examples/img/Caffeine.png?raw=true) ![Alt text](https://github.com/Sassafrass6/XtraCrysPy/blob/master/examples/img/Si_charge_density.png?raw=true)

![Alt text](https://github.com/Sassafrass6/XtraCrysPy/blob/master/examples/img/Colored_surface.png?raw=true) ![Alt text](https://github.com/Sassafrass6/XtraCrysPy/blob/master/examples/img/Textured_surface.png?raw=true)

## Features:
- Plot molecular systems with minimal setup
- Plot systems from DFT or MD input/output files
- Interact with the system to caluclate distances, angles, or receive atomic information
- Display relaxation or MD steps from QE relax output, raw coordinates, LAMMPS trajectory files, or CP2K (xyz)
- Animate trajectories, with control over frames per second and steps per frame.
- Display reciprocal space features
- Plot real space iso-surfaces (charge density, etc)
- Plot reciprocal iso-surfaces (Fermi surface, spin-texture, etc)

## Supported File Types:
- Quantum ESPRESSO input or output
- POSCAR
- CIF
- XSF
- CP2K
- BXSF (coming soon)
- Any ASE supported format

## Requirements:
- python 3.8
- numpy 1.19
- [fury 0.8](https://github.com/fury-gl/fury)
- [ase](https://wiki.fysik.dtu.dk/ase/)

## Installation:

The [Anaconda](https://www.anaconda.com/) Python distribution is recommended.

### Install with pip:
-  pip install xtracryspy
### Install from source:
-  python setup.py install
-  python setup.py install --user

(only use --user if you do not have permission to install python packages)

## Control Inputs:
- 'u' : Toggle UI visibility
- 'a' : Toggle Axis visibility
- 'b' : Toggle Boundary visiblity
- 'c' : Toggle constraint of atoms within cell
- 'k' : Toggle the reciprocal lattice vector visibility
- 's' : Togle the selection type panel visibility
- 'n' : Toggle the cell repition panel visibility
- '>' : Step forward in relax or MD
- '<' : Step backward in relax or MD
- CTRL + ('>' or '<') : Step 5% through the relaxation or MD steps
- SPACE : Pause or resume an active animation
- SHIFT + 'o' : Report the camera position and orientation
- SHIFT + 'c' : Reset camera to default position
- SHIFT + 's' : Take snapshot
- CTRL + 'w' : Exit
- Arrow Keys : Rotate model
- SHIFT + Arrow Keys : Translate camera
- CTRL + Arrow Keys (+ SHIFT) : Rotate (or Translate) in smaller steps

## Display an atomic model from inputfile, with the examples/main.py script:
- python main.py <input_file>
- python main.py <input_file> <bond_distance>

## Examples:
Available in the examples directory

From scratch:
- python CNO3.py
- python caffeine.py
- python silicon_fcc.py
- python example_motion.py
- python example_textured_isosurface.py

From inputfile (with main.py script):
- python main.py <input_file>
- python main.py <input_file> <bond_distance>

From inputfile (programmatic):
- python SnTe_2D.py
- python example_cif.py
- python example_cp2k.py
- python example_qe_in.py
- python example_poscar.py

Relax QE outputfile:
- python example_relax.py

LAMMPS MD trajectory:
- python example_lammps_md.py

Charge Density:
- python example_charge_density.py

Reciprocal Space:
- python example_BZ.py

Isosurfaces:
- python example_isosurface1.py
- python example_isosurface2.py
- python example_isosurface_clip.py
- python example_colored_isosurface.py
- python example_textured_isosurface.py
- python example_textured_isosurface_clip.py
