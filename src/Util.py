import numpy as np

def read_qe_file ( self, fname, ftype='scf.in' ):
  '''
  Read necessary information from QE file. (Currently only input files are readable)

  Arguments:
    fname (str): Quantum Espresso inputfile name
    ftype (str): Type of QE file
  '''
  import re
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

def read_bxsf ( fname ):
  '''
  '''
  bands = None
  b_vec = None
  with open(fname) as f:
    l = f.readline()
    while l != '':
      if 'BANDGRID_3D_BANDS' in l:
        l = f.readline()
        nbnd = int(l.split()[0])
        nx,ny,nz = (int(v) for v in f.readline().split()); f.readline()
        b_vec = np.array([[float(v) for v in f.readline().split()] for _ in range(3)])
        bands = np.zeros((nx,ny,nz,nbnd), dtype=float)
        for n in range(nbnd):
          f.readline()
          for i in range(nx):
            for j in range(ny):
              ls = f.readline().split()
              for k in range(nz):
                bands[i,j,k,n] = float(ls[k])
        break
      l = f.readline()
  return b_vec,bands

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
    cell_param (list): The 6 Quantum Espresso cell parameters (A,B,C,Cos(AB),Cos(AC),Cos(BC))
  '''
  if ibrav == 1:
    A = cell_param[0]
    coords = A * np.array([[1,0,0],[0,1,0],[0,0,0]])
  elif ibrav == 2:
    A = cell_param[0]
    coords = A/2 * np.array([[-1,0,1],[0,1,1],[-1,1,0]])
  elif ibrav == (3 or -3):
    A = cell_param[0]
    sgn = np.sign(ibrav)
    coords = A/2 * np.array([[sgn,1,1],[-sgn,sgn,1],[-sgn,-sgn,sgn]])
  elif ibrav == 4:
    A,C = cell_param[0],cell_param[2]
    coords = A * np.array([[1,0,0],[-.5,-np.sqrt(3)/2,0],[0,0,C/A]])
  elif ibrav == (5 or -5):
    A,CG = cell_param[0],cell_param[3]
    tx,ty,tz = np.sqrt((1-c)/2),np.sqrt((1-c)/6),np.sqrt((1+2*c)/3)
    if ibrav == 5:
      coords = A * np.array([[tx,-ty,tz],[0,2*ty,tz],[-tx,-ty,tz]])
    else:
      u = tz + 2*np.sqrt(2)*ty
      v = tz - np.sqrt(2)*ty
      coords = A/np.sqrt(3) * np.array([[u,v,v],[v,u,v],[v,v,u]])
  elif ibrav == 6:
    A,C = cell_param[0],cell_param[2]
    coords = np.array([[A,0,0],[0,A,0],[0,0,C]])
  elif ibrav == 7:
    A,C = cell_param[0],cell_param[2]
    coords = A/2 * np.array([1,-1,C/A],[1,1,C/A],[-1,-1,C/A])
  elif ibrav == 8:
    A,B,C = cell_param[0],cell_param[1],cell_param[2]
    coords = np.array([[A,0,0],[0,B,0],[0,0,C]])
  elif ibrav == (9 or -9):
    A,B,C = cell_param[0],cell_param[1],cell_param[2]
    sng = np.sign(ibrav)
    coords = .5 * np.array([[A,sng*B,0],[-sgn*A,B,0],[0,0,2*C]])
  elif ibrav == 91:
    A,B,C = cell_param[0],cell_param[1],cell_param[2]
    coords = .5 * np.array([[2*A,0,0],[0,B,-C],[0,B,C]])
  elif ibrav == 10:
    A,B,C = cell_param[0],cell_param[1],cell_param[2]
    coords = .5 * np.array([[A,0,C],[A,B,0],[0,B,C]])
  elif ibrav == 11:
    A,B,C = cell_param[0],cell_param[1],cell_param[2]
    coords = .5 * np.array([[A,B,C],[-A,B,C],[-A,-B,C]])
  elif ibrav == (12 or -12):
    A,B,C = cell_param[0],cell_param[1],cell_param[2]
    if ibrav == 12:
      G = np.arccos(cell_param[3])
      coords = np.array([A,0,0],[B*np.cos(G),B*np.sin(G),0],[0,0,C])
    else:
      G = np.arccos(cell_param[4])
      coords = np.array([A,0,0],[0,B,0],[C*np.cos(G),0,C*np.sin(G)])
  elif ibrav == (13 or -13):
    A,B,C = cell_param[0],cell_param[1],cell_param[2]
    if ibrav == 13:
      G = np.arccos(cell_param[3])
      coords = np.array([A/2,0,-C/2],[B*np.cos(G),B*np.sin(G),0],[A/2,0,C/2])
    else:
      G = np.arccos(cell_param[4])
      coords = np.array([A/2,-B/2,0],[A/2,B/2,0],[C*np.cos(G),0,C*np.sin(G)])
  elif ibrav == 14:
    A,B,C,AB,AC,BC = tuple([cell_param[i] for i in range(6)])
    G0,G1,G2 = np.arccos(AB),np.arccos(AC),np.arccos(BC)
    v1,v2 = [A,0,0],[B*np.cos(G0),B*np.sin(G0),0]
    v3 = [C*np.cos(G1),C*(np.cos(G2)-np.cos(G1)*np.cos(G0))/np.sin(G0),C*np.sqrt(1+2*np.cos(G2)*np.cos(G1)*np.cos(G0)-np.cos(G2)**2-np.cos(G1)**2-np.cos(G0)**2)/np.sin(G0)]
    coords = np.array([v1,v2,v3])
  else:
    raise ValueError('Lattice not defined for ibrav %d'%ibrav)

  return coords

def bravais_boundaries ( b_vec ):
    '''
    Create the Brillouin Zone boundary in the vpython window

    Arguments:
      b_vec (list or ndarray): List of 3 3-d vectors, representing the reciprocal lattice vectors

    Returns:
      (list): List of tuples, each containing two 3-d points between which to draw an edge
    '''
    from numpy.linalg import det,norm,solve

    indices = [[0,0,1],[0,0,-1],[0,1,0],[0,-1,0],[0,1,1],[0,-1,-1],[1,0,0],[-1,0,0],[1,0,1],[-1,0,-1],[1,1,0],[-1,-1,0],[1,1,1],[-1,-1,-1],[1,1,-1],[-1,-1,1],[1,-1,1],[-1,1,-1],[-1,1,1],[1,-1,-1],[-1,1,0],[1,-1,0],[1,0,-1],[-1,0,1],[0,1,-1],[0,-1,1]]
    G = [np.sum([[i[0]],[i[1]],[i[2]]]*b_vec,axis=0) for i in indices]
    incl_G = np.ones(len(G))

    # Determine which G vectors lie on BZ boundary
    #   If G/2 is closer to Gamma than to another reciprocal lattice vector
    for i,g1 in enumerate(G):
      half_G = g1/2
      for j,g2 in enumerate(G[:len(G)//2]):
        if i != j and norm(half_G) >= norm(half_G-g2):
          incl_G[i] = 0
    planes = np.array([g for i,g in enumerate(G) if incl_G[i]])

    corners = []
    # Search combinations of planes and save intersections as corners
    for i,p1 in enumerate(planes):
      for j,p2 in enumerate(planes[i+1:]):
        for k,p3 in enumerate(planes[j+1:]):
          M = np.array([p1,p2,p3])
          sqr = lambda v : np.sum(v*v)
          magG = .5 * np.array([sqr(p1),sqr(p2),sqr(p3)])
          if not np.isclose(det(M),0.):
            c = solve(M,magG)
            corners.append(c)

    # Set near zero vlues to zero explicitly
    for i,c in enumerate(corners):
      for j,v in enumerate(c):
        if np.isclose(v,0):
          corners[i][j] = 0.

    # Select corners closer to Gamma than to another reciprocal lattice vector
    # Then, eliminate any repeated corners
    incl_C = np.ones(len(corners))
    for i,c in enumerate(corners):
      for g in G:
        if norm(c) > norm(c-g)+1e-10:
          incl_C[i] = 0

    # Eliminate duplicate corners
    corners = np.unique(np.array([t for i,t in enumerate(corners) if incl_C[i]]), axis=0)

    pairs = []
    # Compute pairs of points which share two planes
    for ci,c1 in enumerate(corners):
      for c2 in corners[ci+1:]:
        for pi,p1 in enumerate(planes):
          for p2 in planes[pi+1:]:
            dists = [np.abs(np.sum(p*c)-np.sum(p*p)/2) for p in [p1,p2] for c in [c1,c2]]
            if all(np.isclose(d,0.) for d in dists):
              pairs.append((c1,c2))

    return pairs
