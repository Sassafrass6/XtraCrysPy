import vpython as vp
import numpy as np

class XtraCrysPy:

  def __init__ ( self, inputfile=None, fmt="SCF", lattice=None, basis=None, basis_labels=None, origin=[0,0,0], species=None, bonds=None ):
    '''
    Initialize the XtraCrysPy object, creating a canvas and computing the corresponding lattice

    Arguments:
      inputfile (str): Filename of a quantum espresso inputfile
      relax (bool): True if a QE relax output file is the inputfile
      lattice (list): List of 3 3d vectors, representing the lattice parameters
      basis (list): List of N 3d vectors, representing the positions of each of N atoms
      origin (list): List of 3 points, representing the x,y,z position of the grid's center
      species (list): List of N species strings corresponding to atoms in 'basis' (e.g. ['Si','Si'])
      bonds (dict or dict): Maximum bond distance or dictionary of form {'Sp1_Sp2':max_bond_dist}
    '''

    self.spec = {}           # Dicitonary of atomic species
    self.nspec = 0           # Number of species
    self.natoms = 0          # Number of atoms
    self.view = None         # VPython View
    self.ibrav = None        # Bravais lattice ID, following QE indexing
    self.atoms = None        # Atomic positions (basis)
    self.bonds = None        # Dictionary of bond distances
    self.volume_data = None  # numpy array that describes volume data
    self.format = fmt        # format of the file
    self.cameras = None      # New feature which could provide camera locations to Blender
    self.lattice = None      # Unit vectors (lattice)
    self.origin = origin     # Origin of the figure. Default [0,0,0]
    self.coord_type = None   # Units (angstrom, bohr, alat, crystal, manual)
    self.relax_poss = None   # Atomic positions for each step of relaxation
    self.cell_param = None   # Average of norms for lattice?
    self.relax_poss = None   # Sets of coordinates by relaxation step
    self.relax_index = None  # Index of the current relax step being modeled
    self.relax_steps = None  # Number of relaxation steps, if applicable
    self.basis_labels = None # Species label for each atom in the basis

    self.WHITE = (1,1,1,1)
    self.BLACK = (0,0,0,1)
    self.DEFAULT_RADIUS = 1
    self.BOHR_TO_ANGSTROM = .52917720

    relax = (format == "RELAX")
    if inputfile is None:
      self.coord_type = 'manual'
      if lattice is None or basis is None or species is None:
        print('Lattice and Basis not defined. Only \'plot_bxsf\' will function.')
      else:
        self.atoms = np.array(basis)
        self.spec = species
        self.natoms = len(basis)
        self.lattice = np.array(lattice)
        self.basis_labels = basis_labels
        self.cell_param = [np.abs(np.mean([np.linalg.norm(v) for v in self.lattice]))]
    else:
      # Read coords from file
      from .Util import qe_lattice
      if self.format == "SCF":
        from .Util import read_relax_file
        read_relax_file(self, inputfile)
      elif self.format == "RELAX":
        from .Util import read_scf_file
        read_scf_file(self, inputfile)
      elif self.format == "CUBE":
        from .Util import read_cube_file
        read_cube_file(self, inputfile)
      else:
        raise ValueError("format has not been implemented yet")
      self.lattice = qe_lattice(self.ibrav, self.cell_param)

      if relax:
        # Setup model to reflect the first relax position
        self.relax_steps = len(self.relax_lattices)
        if self.relax and self.relax_steps > 0:
          self.relax_index = 0
          self.lattice = self.relax_lattices[0]
        else:
          raise ValueError('No relax steps were found.')

        if 'angstrom' in self.coord_type:
          for i in range(len(self.relax_poss)):
            self.relax_poss[i] /= self.BOHR_TO_ANGSTROM
          else:
            self.relax_poss[i][:,:] = np.array([self.atomic_position(a) for a in self.relax_poss[i]])
        self.atoms = self.relax_poss[0]

      else:
        if 'angstrom' in self.coord_type:
          self.atoms = np.array(self.atoms) / self.BOHR_TO_ANGSTROM
        elif 'alat' in self.coord_type or 'crystal' in self.coord_type:
          self.atoms =  np.array([self.atomic_position(a) for a in self.atoms])
          print('Crystal & Alat coordinate types require more testing')

    if species is not None:
      if self.basis_labels is None:
        raise ValueError('Must specify basis_labels if the species dictionary is provided.')
      if len(species) < len(list(set(self.basis_labels))):
        raise ValueError('Must specify every label which appears in basis_labels')
      for bl in self.basis_labels:
        if bl not in species:
          raise ValueError('Must provide entry for %s in species dictionary.'%bl)
          
      self.spec = species
      self.nspec = len(species)
      for i,tup in enumerate(self.spec.items()):
        k,v = tup
        self.spec[k]['id'] = i
        if 'color' not in self.spec[k]:
          self.spec[k]['color'] = self.WHITE
        if 'radius' not in self.spec[k]:
          self.spec[k]['radius'] = self.DEFAULT_RADIUS

    else:
      if self.basis_labels is None:
        self.nspec = self.natoms
        self.spec = {i:{'id':i, 'radius':self.DEFAULT_RADIUS, 'color':self.WHITE} for i in range(self.natoms)}
      else:
        unique_spec = list(set(self.basis_labels))
        for u in unique_spec:
          self.spec[u] = {'id':self.nspec, 'radius':self.DEFAULT_RADIUS, 'color':self.WHITE}
          self.nspec += 1

    # Sort the bond distances
    if bonds is not None:
      self.bonds = {'%s_%s'%tuple(sorted(k.split('_'))):bonds[k] for k in bonds.keys()}

    # Tranlsate atoms lying outside the unti cell back inside
    if self.lattice is not None and self.atoms is not None and not self.ibrav == 1:
      from .Util import constrain_atoms_to_unit_cell
      self.atoms = constrain_atoms_to_unit_cell(self.lattice, self.atoms)

  def reciprocal_lattice ( self, lattice=None ):
    '''
    Construct reciprocal lattice vectors

    lattice (list or ndarray): 3 3d vectors of the real space lattice
    '''
    if lattice is None:
      if self.lattice is None:
        raise ValueError('No lattice to convert')
      lattice = self.lattice
    return 2*np.pi*np.linalg.inv(lattice).T

  def atomic_position ( self, v, lattice=None ):
    '''
    Calculate the atom's position within the unit cell

    Arguments:
      v (list or ndarray): List of 3 values, each to weight one of the 3 lattice vectors
      lattice (list or ndarray): 3 3d vectors of the real space lattice
    '''
    if lattice is None:
      if self.lattice is None:
        raise ValueError('No lattice to convert')
      lattice = self.lattice
    return v @ lattice

  def get_boundary_positions ( self, lattice=None, nx=1, ny=1, nz=1 ):
    '''
    Compute the boundary endpoints for a lattice, repeating nx,ny,nz times in the respective direction.

    Arguments:
      lattice (ndarray): Array of 3 3-vectors, representing the lattice. None will take self.lattice
      nx (int): Number of cells to draw in the x direction
      ny (int): Number of cells to draw in the y direction
      nz (int): Number of cells to draw in the z direction
    '''
    if lattice is None:
      if self.lattice is None:
        raise ValueError('No lattice to convert')
      lattice = self.lattice

    edges = []
    vertex_ind = 0
    vertices = [[0,0,0]]
    vertex_map = {(0,0,0):0}
    for ix in range(nx+1):
      for iy in range(ny+1):
        for iz in range(nz+1):
          orig_ind = vertex_map[(ix,iy,iz)]
          orig = self.atomic_position([ix,iy,iz], lattice)
          if ix < nx:
            vertex_ind += 1
            vertices += [orig+lattice[0]]
            edges += [[orig_ind, vertex_ind]]
            vertex_map[(ix+1,iy,iz)] = vertex_ind
          if iy < ny:
            vertex_ind += 1
            vertices += [orig+lattice[1]]
            edges += [[orig_ind, vertex_ind]]
            vertex_map[(ix,iy+1,iz)] = vertex_ind
          if iz < nz:
            vertex_ind += 1
            vertices += [orig+lattice[2]]
            edges += [[orig_ind, vertex_ind]]
            vertex_map[(ix,iy,iz+1)] = vertex_ind
    vertices = np.array(vertices) + self.origin
    return vertices, edges

  def get_atomic_positions ( self, lattice=None, atoms=None, nx=1, ny=1, nz=1 ):
    '''
    Compute the coordinates of each atomic position for a lattice and basis, repeating nx,ny,nz times in the respective direction.

    Arguments:
      lattice (ndarray): Array of 3 3-vectors, representing the lattice. None will take self.lattice
      atoms (ndarray): Array of N 3-vectors, representing the basis. None will take self.atoms
      nx (int): Number of cells to draw in the x direction
      ny (int): Number of cells to draw in the y direction
      nz (int): Number of cells to draw in the z direction
    '''
    if lattice is None:
      if self.lattice is None:
        raise ValueError('No lattice to convert')
      lattice = self.lattice
    if atoms is None:
      if self.atoms is None:
        raise ValueError('No lattice to convert')
      atoms = self.atoms

    positions = []
    for ix in range(nx):
      for iy in range(ny):
        for iz in range(nz):
          c_pos = self.origin + self.atomic_position([ix,iy,iz], lattice)
          for i,a in enumerate(atoms):
            positions += [c_pos + a]
    return positions

  
  def get_bond_pairs ( self, atoms ):
    '''
    '''
    bonds = []
    if self.bonds is None or len(atoms) < 2:
      return bonds

    na = self.natoms
    for i,k in enumerate(self.bonds.keys()):
      sA,sB = tuple(k.split('_'))
      for m in range(len(atoms)):
        for n in range(m+1, len(atoms)):
          l1,l2 = self.basis_labels[m%na],self.basis_labels[n%na]
          if (l1 == sA and l2 == sB) or (l1 == sB and l2 == sA):
            if np.linalg.norm(atoms[m] - atoms[n]) <= self.bonds[k]:
              bonds += [[m, n]]
    return bonds
    

  def get_BZ_corners ( self, rlat=None ):
    '''
    Compute the corner positions of the Brillouin Zone

    Arguments:
      rlat (list or ndarray): 3 3d vectors representing the reciprocal lattice. If 'None' vectors will be computed from 'self.lattice'
    '''
    if rlat is None:
      if self.rlattice is None:
        if self.lattice is None:
          raise ValueError('No lattice to convert')
        self.rlattice = self.reciprocal_lattice(self.lattice)
      rlat = self.rlattice
    self.view.draw_BZ_points(points, color, rlat)

  def get_color_list ( self, nx, ny, nz ):
    colors = []
    for ix in range(nx):
      for iy in range(ny):
        for iz in range(nz):
          for i in range(len(atoms)):
            colors += [self.spec[self.basis_labels[i]]['color']]
    return colors

  def write_blender_xml ( self, directory='./', fname='xcpy_system.xml', frame=True, nx=1, ny=1, nz=1 ):
    from os.path import join as opjoin
    from .Util import blender_xml

    atom_positions = self.get_atomic_positions(nx=nx, ny=ny, nz=nz)
    frame = self.get_boundary_positions(nx=nx, ny=ny, nz=nz)
    bonds = self.get_bond_pairs(atom_positions)

    blender_xml(opjoin(directory,fname), self.natoms, self.spec, self.basis_labels, atom_positions, bonds, frame, self.cameras)

  def start_cryspy_view ( self, title='', w_width=1000, w_height=750, perspective=False, f_color=(1,1,1), bg_color=(0,0,0), nx=1, ny=1, nz=1 ):
    from .View import View
    frame = self.get_boundary_positions(nx=nx, ny=ny, nz=nz)
    atoms = self.atoms if not self.relax else self.relax_poss
    bonds = self.get_bond_pairs(atoms)
    print(bonds)
    model = [self.lattice, atoms, self.spec, self.basis_labels, bonds]
    self.view = View(title,w_width,w_height,self.origin,model,perspective,frame,f_color,bg_color,nx,ny,nz)
