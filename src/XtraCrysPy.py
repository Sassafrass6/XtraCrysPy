import vpython as vp
import numpy as np

class XtraCrysPy:

  def __init__ ( self, inputfile=None, relax=False, lattice=None, basis=None, basis_labels=None, origin=[0,0,0], species=None, bonds=None ):
    '''
    Initialize the XtraCrysPy object, creating a canvas and computing the corresponding lattice

    Arguments:
      inputfile (str): Filename of a quantum espresso inputfile
      relax (bool): True if a QE relax output file is the inputfile
      lattice (list): List of 3 3d vectors, representing the lattice parameters
      basis (list): List of N 3d vectors, representing the positions of each of N atoms
      origin (list): List of 3 points, representing the x,y,z position of the grid's center
      species (list): List of N species strings corresponding to atoms in 'basis' (e.g. ['Si','Si'])
      bonds (dict): Maximum bond distance or dictionary of form {'Sp1_Sp2':max_bond_dist}
    '''

    self.spec = {}           # Dicitonary of atomic species
    self.nspec = 0           # Number of species
    self.natoms = 0          # Number of atoms
    self.view = None         # VPython View
    self.ibrav = None        # Bravais lattice ID, following QE indexing
    self.atoms = None        # Atomic positions (basis)
    self.bonds = None        # Dictionary of bond distances
    self.relax = relax       # Boolean. True if there are relaxation steps
    self.cameras = None      # New feature which could provide camera locations to Blender
    self.lattice = None      # Unit vectors (lattice)
    self.origin = origin     # Origin of the figure. Default [0,0,0]
    self.coord_type = None   # Units (angstrom, bohr, alat, crystal, manual)
    self.relax_poss = None   # Atomic positions for each step of relaxation
    self.relax_poss = None   # Sets of coordinates by relaxation step
    self.basis_labels = None # Species label for each atom in the basis

    self.WHITE = (1,1,1,1)
    self.BLACK = (0,0,0,1)
    self.DEFAULT_RADIUS = 1
    self.BOHR_TO_ANGSTROM = .52917720

    if inputfile is None:
      self.coord_type = 'manual'
      if lattice is None or basis is None:
        print('Lattice and Basis not defined. Only \'plot_bxsf\' will function.')
      else:
        self.atoms = np.array(basis)
        self.spec = species
        self.bonds = bonds
        self.natoms = len(basis)
        self.lattice = np.array(lattice)
        self.basis_labels = basis_labels
    else:
      # Read coords from file
      from .Util import qe_lattice
      if relax:
        from .Util import read_relax_file
        data = read_relax_file(inputfile)
        self.ibrav,self.relax_poss = data[0],data[1]
        self.relax_lattices,self.basis_labels = data[2],data[3]
        self.cell_param,self.coord_type = data[4],data[5]
        self.natoms = self.relax_poss.shape[1]
      else:
        from .Util import read_scf_file
        data = read_scf_file(inputfile)
        self.ibrav,self.atoms = data[0],data[1]
        self.basis_labels,celldm,cell_param = data[2],data[3],data[4]
        self.coord_type = data[5]
        self.natoms = self.atoms.shape[0]

      if relax:
        # Setup model to reflect the first relax position
        relax_steps = len(self.relax_lattices)
        if self.relax and relax_steps > 0:
          self.relax_index = 0
          self.lattice = self.relax_lattices[0]
        else:
          raise ValueError('No relax steps were found.')

        if 'angstrom' in self.coord_type:
          self.relax_poss /= self.BOHR_TO_ANGSTROM
        else:
          for i in range(relax_steps):
            self.relax_poss[i,:,:] = np.array([self.atomic_position(a,lattice=self.relax_lattices[i]) for a in self.relax_poss[i]])
        self.atoms = self.relax_poss[0]

      else:
        if self.ibrav == 0:
          self.lattice = cell_param
        else:
          self.lattice = qe_lattice(self.ibrav, celldm)
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
        self.basis_labels = range(self.natoms)
        self.spec = {i:{'id':i, 'radius':self.DEFAULT_RADIUS, 'color':self.WHITE} for i in self.basis_labels}
      else:
        self.spec = {}
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

  
  def draw_BZ ( self, rlat=None ):
    '''
    Draw the Brillouin Zone

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
    from .model import get_bond_pairs
    from .Util import blender_xml

    atom_positions = self.get_atomic_positions(nx=nx, ny=ny, nz=nz)
    frame = self.get_boundary_positions(nx=nx, ny=ny, nz=nz)
    bonds = get_bond_pairs(self.natoms, atom_positions, self.bonds, self.basis_labels)

    blender_xml(opjoin(directory,fname), self.natoms, self.spec, self.basis_labels, atom_positions, bonds, frame, self.cameras)

  def start_cryspy_view ( self, title='', w_width=1000, w_height=750, perspective=False, f_color=None, bg_color=(0,0,0), nx=1, ny=1, nz=1 ):
    '''
    Arguments:
      f_color (tuple): (R,G,B) color for the unit cell Frame. None defaults to white with no inital display.
    '''
    from .View import View
    atoms = self.atoms if not self.relax else self.relax_poss
    lattice = self.lattice if not self.relax else self.relax_lattices
    model = [lattice, atoms, self.spec, self.basis_labels, self.bonds]
    self.view = View(title,w_width,w_height,self.origin,model,perspective,f_color,bg_color,nx,ny,nz)

  def plot_spin_texture ( self, fermi_fname, spin_fname, colors=None, e_up=1, e_dw=-1, title='', w_width=1000, w_height=750, f_color=(1,1,1), bg_color=(0,0,0)):
    '''
    Plots the spin texture read from the format output by PAOFLOW
    Arguments:
      fermi_fname (str): Name for Fermi surface file, representing energy eigenvalues in BZ for such band
      spin_fname (str): Name for spin texture file, representing spin direction in BZ for such band
      e_up (float): BZ points with energies higher than e_up wont be plotted
      e_dw (float): BZ points with energies lower than e_dw wont be plotted 
    '''
    from scipy.fftpack import fftshift
    from .View import View

    self.recip_space = True
    fdata,sdata = np.load(fermi_fname),np.load(spin_fname)
    eig,spins = fdata['nameband'],sdata['spinband']

    eig = fftshift(eig,axes=(0,1,2))
    spins = np.moveaxis(np.real(spins[:,:,:,:]),spins.ndim-1,0)

    colscale = np.mean(eig)
    _,nx,ny,nz = spins.shape
    vp_shift = .5 * np.array([1,1,1])
    points,directions = [],[]
    c_append = False
    if colors is None:
      colors = []
      c_append = True
    rlat = self.reciprocal_lattice(self.lattice)
    for x in range(nx):
      for y in range(ny):
        for z in range(nz):
          bv = eig[x,y,z]
          if bv > e_dw and bv < e_up:
            points.append(([x/nx,y/ny,z/nz]-vp_shift) @ rlat)
            directions.append(spins[:,x,y,z])
            if c_append:
              colors.append([bv/colscale,.5,.1])

    self.view = View(title,w_width,w_height,self.origin,None,False,f_color,bg_color,1,1,1)
    self.view.draw_arrows(rlat, points, directions, colors, .01)


  def plot_bxsf ( self, fname, iso=[0], bands=[0], colors=[[0,1,0]], normals=True, title='', w_width=1000, w_height=750, perspective=False, f_color=(1,1,1), bg_color=(0,0,0), write_obj=False ):
    '''
    Create the Brillouin Zone boundary and bsxf points at iso value for each band in the vpython window

    Arguemnts:
      fname (str): Name of bxsf file
      iso (list): List of floats corresponding to the isosurface values for each respective band in 'bands'
      bands (list): List of integers representing the index of the band to plot
      colors (list): List of 3-d RGB color vectors for each band. If ignored, each band will be green
      normals (bool): True adds normals to triangle vertices, improving surface visibility
      write_obj (bool): True will write the triange vertex, face, and normal data to an obj file
    '''
    from .View import View
    from .Util import read_bxsf

    if len(iso) != len(bands):
      raise ValueError("Each band in 'bands' should have one corresponding isosurface in 'iso'")
    if len(iso) != len(colors) and len(colors) != 1:
      raise ValueError("Specify 1 color to plot all bands in the same color, or specify 1 color for each band.")

    self.recip_space = True
    b_vec,data = read_bxsf(fname)

    if np.max(bands) >= data.shape[-1]:
      raise ValueError("'bands' contains an index too large for the dataset")
    self.view = View(title,w_width,w_height,self.origin,None,perspective,f_color,bg_color,1,1,1)

    self.view.draw_surface(False, b_vec, data, iso, bands, colors, normals, write_obj)

  def plot_xsf ( self, fname, iso=0, color=[0,1,0], normals=True, title='', w_width=1000, w_height=750, perspective=False, f_color=(1,1,1), bg_color=(0,0,0), write_obj=False ):
    '''
    Create the Brillouin Zone boundary and bsxf points at iso value for each band in the vpython window

    Arguemnts:
      fname (str): Name of bxsf file
      iso (list): List of floats corresponding to the isosurface values for each respective band in 'bands'
      bands (list): List of integers representing the index of the band to plot
      colors (list): List of 3-d RGB color vectors for each band. If ignored, each band will be green
      normals (bool): True adds normals to triangle vertices, improving surface visibility
      write_obj (bool): True will write the triange vertex, face, and normal data to an obj file
    '''
    from .View import View
    from .Util import read_xsf

    pvec,pcoord,data = read_xsf(fname)

    self.view = View(title,w_width,w_height,self.origin,None,perspective,f_color,bg_color,1,1,1)

    self.view.draw_surface(True, pvec, data, [iso], [0], [color], normals, write_obj)
