from XtraCrysPy import Atomic,Model
import numpy as np

if '__main__' == __name__:

  # Atomic positions
  atoms = np.array([[2.235, 0.000, 0.000],
                    [3.300, 2.003, 0.000],
                    [0.000, 1.252, 0.000],
                    [0.990, 1.973, 0.000],
                    [2.269, 1.268, 0.000]])
  # Label for each atom and colors
  labels = ['O']*3 + ['N'] + ['C']

  # Bond distance in Angstrom
  bonds = 1.5

  # Visual information: species, lattice, atoms, bonds, units (lattice is specified in angstrom)
  a_info = {'species':labels,
            'abc':atoms,
            'bonds':bonds,
            'aunit':'angstrom'}

  xcp = Atomic.Atomic(params=a_info, atom_res=(6,10))
  xcp.start_crystal_view()

