import numpy as np

def read_scf_file ( self, fname ):
  '''
  Read necessary information from QE file. (Currently only input files are readable)

  Arguments:
    fname (str): Quantum Espresso inputfile name
  '''
  import re
  with open(fname) as f:
    strip_int = lambda s : int(re.search(r'\d+',s).group())
    strip_float = lambda s : float(re.search(r'\d+.\d+',s).group())
    natoms = 0
    self.atoms = []
    self.basis_labels = []
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
          self.coord_type = ls[1] if len(ls)>1 else 'alat'
          while natoms > 0:
            ls = f.readline().split()
            if len(ls) > 0:
              natoms -= 1
              self.basis_labels.append(ls[0])
              self.atoms.append([float(v) for v in ls[1:4]])
      line = f.readline()
    self.atoms = np.array(self.atoms)
    self.cell_param[1] *= self.cell_param[0]
    self.cell_param[2] *= self.cell_param[0]

def read_relax_file ( self, fname ):
  '''
  Read relax inputfile

  Arguments:
    fname (str): Quantum Espresso inputfile name
  '''
  import re
  self.relax_poss = []
  self.cell_param = 6*[0]
  self.relax_lattices = []
  self.coord_type = 'relax'
  strip_int = lambda s : int(re.search(r'\d+',s).group())
  strip_float = lambda s : float(re.search(r'\d+.\d+',s).group())
  with open(fname) as f:
    l = f.readline()
    while l != '':
      if 'atoms/cell' in l:
        self.natoms = int(l.split()[-1])
      elif 'bravais' in l:
        self.ibrav = int(l.split()[3])
      elif 'celldm' in l:
        ls = l.split()
        for i in range(0,len(ls),2):
          ip = strip_int(ls[i])-1
          self.cell_param[ip] = strip_float(ls[i+1])
      elif 'CELL_PARAMETERS' in l:
        alat = strip_float(l)
        avec = np.zeros((3,3), dtype=float)
        for i in range(3):
          ls = f.readline().split()
          avec[i] = np.array([alat*float(s) for s in ls])
        self.relax_lattices.append(avec)
      if 'ATOMIC_POSITIONS' in l:
        if self.natoms is None:
          raise ValueError('natoms not found in relax file. Are you sure that a QE relax output file was provided?')
        ls = l.split()
        if len(ls) > 1:
          self.coord_type = ls[1]
        self.basis_labels = []
        apos = np.zeros((self.natoms,3), dtype=float)
        for i in range(self.natoms):
          ls = f.readline().split()
          self.basis_labels.append(ls[0])
          apos[i,:] = np.array([float(v) for v in ls[1:]])
        self.relax_poss.append(apos)
      l = f.readline()
  self.cell_param[1] *= self.cell_param[0]
  self.cell_param[2] *= self.cell_param[0]

def read_cube_file ( self, fname ):
  '''
  Read GAUSSIAN cube files

  Arguments:
    fname (str): .cube formatted filename
  '''
  ignore_lines = 2
  n = []
  species = set()
  self.lattice = np.zeros((3, 3))
  self.atoms = []
  with open(fname) as f:
    for _ in range(ignore_lines):
      f.readline()
    # read number of atoms and position of origin
    l = f.readline().split()
    self.natoms = int(l[0])
    self.origin = [float(v) for v in l[1:]]
    # read size and transform matrix
    for i in range(3):
      l = f.readline().split()
      n.append(abs(int(l[0])))
      # format implies possibility of mixed inputs, however.. just assume this isn't the case
      self.coord_type = "BOHR" if n[0] > 0 else "ANGSTROM"
      self.lattice[i] = [float(v) for v in l[1:]]
    # read atomic info
    for i in range(self.natoms):
      l = f.readline().split()
      # ignore charge for now..
      species.add(int(l[0]))
      self.atoms.append([float(v) for v in l[2:]])
    
    # sort out species things (TODO add if needed)
    #self.nspec = len(species)
    #self.spec = dict()
    
    self.volume_data = np.zeros(tuple(n))
    i = j = k = 0
    # only works if python > 3.8!
    while l := f.readline().split():
      for val in l:
        try:
          self.volume_data[i, j, k] = float(val)
        except IndexError:
          print("number of data points do not match up with grid dimensions")
          # debug this if it's encountered, for now just return
          return
        k += 1
        if k >= n[2]:
          k = 0
          j += 1
          if j >= n[1]:
            j = 0
            i += 1

def read_bxsf ( fname ):
  '''
  Read BXSF file, originally formatted for XCrysDen
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
    tx,ty,tz = np.sqrt((1-CG)/2),np.sqrt((1-CG)/6),np.sqrt((1+2*CG)/3)
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
    coords = .5 * np.array([A,-A,C],[A,A,C],[-A,-A,C])
  elif ibrav == 8:
    A,B,C = cell_param[0],cell_param[1],cell_param[2]
    coords = np.array([[A,0,0],[0,B,0],[0,0,C]])
  elif ibrav == (9 or -9):
    A,B,C = cell_param[0],cell_param[1],cell_param[2]
    sgn = np.sign(ibrav)
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

def constrain_atoms_to_unit_cell ( lattice, atoms ):
  '''
  Force atoms to lie inside of the unit cell.

  Arguments:
    lattice (ndarray): 3x3 array of lattice vectors
    atoms (list): List of 3 component vectors representing the atomic basis

  Returns:
    (list): Updated list of atomic positions
  '''

  def inside_cell ( apos ):
    ''' Helper function returning True if the position is interior to the cell '''
    for i,l in enumerate(lattice):
      tvec = np.linalg.inv(lattice).T @ apos
      for v in tvec:
        if v >= 1 or v < 0:
          return False
    return True

  for i,a in enumerate(atoms):
    while not inside_cell(atoms[i]):
      nlat = np.linalg.inv(lattice).T @ atoms[i]
      for j,n in enumerate(nlat):
        if n < 0:
          atoms[i] += lattice[j]
        elif n >= 1:
          atoms[i] -= lattice[j]

  return atoms

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
  G = [i@b_vec for i in indices]
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
  sqr = lambda v : np.sum(v*v)
  # Search combinations of planes and save intersections as corners
  for i,p1 in enumerate(planes):
    for j,p2 in enumerate(planes[i+1:]):
      for k,p3 in enumerate(planes[j+1:]):
        M = np.array([p1,p2,p3])
        magG = .5 * np.array([sqr(p1),sqr(p2),sqr(p3)])
        if not np.isclose(det(M),0.):
          c = solve(M,magG)
          corners.append(c)

  corners = np.array(corners)
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
          incl = True
          for p in [(p1,c1), (p1,c2), (p2,c1), (p2,c2)]:
            d = np.sum(p[0]*(p[1] - p[0]/2))
            if not np.isclose(d, 0.):
              incl = False
              break
          if incl:
            pairs.append((c1,c2))

  # Adjust planes to represent midpoints between reciprocal lattice points
  planes = np.array(planes)/2

  return planes, np.array(pairs)

def blender_xml ( fnpath, natoms, species, labels, atomic_poss, bonds, frame, cameras ):
  with open(fnpath, 'w') as f:
    f.write('<?xml version="1.0"?>\n')
    f.write('<XCP_DATA>\n')

    f.write('\t<SCENE>\n')
    if cameras is not None:
      for i,c in enumerate(cameras.values()):
        f.write('\t\t<CAMERA id="%d">\n'%i)
        f.write('\t\t\t<POSITION>%f %f %f</POSITION>'%tuple(c['position']))
        f.write('\t\t\t<ROTATION>%f %f %f</ROTATION>'%tuple(c['rotation']))
        f.write('\t\t</CAMERA>\n')
    f.write('\t</SCENE>\n')

    f.write('\t<SYSTEM>\n')
    for k,v in species.items():
      f.write('\t\t<SPECIES id="%d" label="%s">\n'%(species[k]['id'], k))
      f.write('\t\t\t<RADIUS>%f</RADIUS>\n'%species[k]['radius'])
      f.write('\t\t\t<COLOR>%f %f %f %f</COLOR>\n'%species[k]['color'])
      f.write('\t\t</SPECIES>\n')

    for i,a in enumerate(atomic_poss):
      key = labels[i%natoms]
      f.write('\t\t<ATOM id="%d" species="%d">'%(i,species[key]['id']))
      f.write('%f %f %f</ATOM>\n'%tuple(a))

    DEFAULT_BOND_TYPE = 1
    for i,b in enumerate(bonds):
      tup = (i, DEFAULT_BOND_TYPE, b[0], b[1])
      f.write('\t\t<BOND id="%d" type="%d" A="%d" B="%d"></BOND>\n'%tup)
    f.write('\t</SYSTEM>\n')

    vertices,edges = frame[0],frame[1]
    f.write('\t<FRAME>\n')
    for i,v in enumerate(vertices):
      f.write('\t\t<VERTEX id="%d">%f %f %f</VERTEX>\n'%((i,)+tuple(v)))
    for i,e in enumerate(edges):
      f.write('\t\t<EDGE id="%d" A="%d" B="%d"></EDGE>\n'%(i, e[0], e[1]))
    f.write('\t</FRAME>\n')

    f.write('</XCP_DATA>\n')
