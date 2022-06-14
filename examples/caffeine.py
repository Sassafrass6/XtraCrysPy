from XtraCrysPy import Atomic,Model
import numpy as np

if '__main__' == __name__:

  # Cell parameters (Angstrom)
  cubic_vecs = np.array([[7.33,0,0],[0.,6.36,0],[0,0,1.80]])

  # Atomic basis (crystal coordinates)
  atoms = np.array([[0.49,0.91,0.5],[0.0,0.43,0.5],[0.29,0.3,0.5],[0.73,0.52,0.5],[0.24,0.67,0.5],[0.62,0.2,0.5],[0.54,0.54,0.5],[0.48,0.34,0.5],[0.43,0.73,0.5],[0.17,0.46,0.5],[0.77,0.31,0.5],[0.23,0.08,0.5],[0.86,0.69,0.5],[0.11,0.85,0.5],[0.91,0.25,0.5],[0.28,0.0,0.0],[0.08,0.07,0.5],[0.28,0.0,1.0],[1.0,0.63,0.5],[0.84,0.79,0.0],[0.84,0.79,1.0],[0.18,1.0,0.5],[0.03,0.83,0.99],[0.03,0.83,0.01]])

  # Label for each atom and colors
  labels = ['O']*2 + ['N']*4 + ['C']*8 + ['H']*10

  # Bond distance in Angstrom
  bonds = 1.5

  # Visual information: species, lattice, atoms, bonds, units (lattice is specified in angstrom)
  a_info = {'species':labels, 'lattice':cubic_vecs, 'abc':atoms,
            'bonds':bonds, 'lunit':'angstrom', 'aunit':'crystal'}

  model = Model.Model(params=a_info)
  xcp = Atomic.Atomic(model=model, boundary=False)
  xcp.start_crystal_view()

