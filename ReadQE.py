import numpy as np
import re

def read_qe_file ( self, fname, ftype='scf.in' ):
  '''
  Read necessary information from QE file. (Currently only input files are readable)

  Arguments:
    fname (str): Quantum Espresso inputfile name
    ftype (str): Type of QE file
  '''
  if ftype == 'scf.in':
    with open(fname) as f:
      strip_int = lambda s : int(re.search(r'\d+',s).group())
      strip_float = lambda s : float(re.search(r'\d+.\d+',s).group())
      natoms = 0
      self.spec = []
      self.atoms = []
      self.cell_param= 6*[0]

      line = f.readline()
      while line != '':
        for l in line.split(','):
          if 'ibrav' in l:
            self.ibrav = strip_int(l)
          elif 'celldm' in l:
            dm = strip_int(l) - 1
            self.cell_param[dm] = strip_float(l)
          elif 'nat' in l:
            self.natoms = natoms = strip_int(l)
          elif 'ATOMIC_POSITIONS' in l:
            ls = l.split()
            self.coord_type = ls[1][1:-1] if len(ls)>1 else 'alat'
            while natoms > 0:
              ls = f.readline().split()
              if len(ls) == 4:
                natoms -= 1
                self.spec.append(ls[0])
                self.atoms.append([float(v) for v in ls[1:]])
        line = f.readline()
      self.cell_param[1] *= self.cell_param[0]
      self.cell_param[2] *= self.cell_param[0]

def crystal_conversion ( atoms, coords, coord_type ):
  '''
  Convert positions of 'atoms' from 'coord_type' to crystal coordinates

  pos (list or ndarray): List of atoms, each with 3 positions to convert to crystal
  coord (list or ndarray): List of 3 3d lattice vectors
  coord_type (str): String - 'angstrom'
  '''

  if coord_type == 'angstrom':
    conv =  0.529177210
    return np.array([[a[i]/(conv*np.linalg.norm(coords[i])) for i in range(3)] for a in atoms])
  else:
    raise ValueError('coord_type %s not implemented'%coord_type)


def qe_lattice ( ibrav, cell_param ):
  '''
  Compute the lattice vectors, given cell parameters and ibrav number

  Arguments:
    ibrav (int): Quantum Espresso ibrav number
    cell_param (list): The 6 Quantum Espresso cell parameters
  '''

  if ibrav == 1:
    A = cell_param[0]
    coords = A * np.array([[1,0,0],[0,1,0],[0,0,0]])
  elif ibrav == 2:
    A = cell_param[0]
    coords = A/2 * np.array([[-1,0,1],[0,1,1],[-1,1,0]])
  elif ibrav == 3:
    A = cell_param[0]
    coords = A/2 * np.array([[1,1,1],[-1,1,1],[-1,-1,1]])
  elif ibrav == -3:
    A = cell_param[0]
    coords = A/2 * np.array([[-1,1,1],[1,-1,1],[1,1,-1]])
  elif ibrav == 4:
    A,C = cell_param[0],cell_param[2]
    coords = A * np.array([[1,0,0],[-.5,-np.sqrt(3)/2,0],[0,0,C/A]])
  elif ibrav == 8:
    A,B,C = cell_param[0],cell_param[1],cell_param[2]
    coords = np.array([[A,0,0],[0,B,0],[0,0,C]])
  else:
    raise ValueError('Lattice not defined for ibrav %d'%ibrav)

  return coords
