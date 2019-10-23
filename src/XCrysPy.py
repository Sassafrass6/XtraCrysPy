import vpython as vp
import numpy as np

class XCrysPy:

  def __init__ ( self, qe_fname=None, lattice=None, basis=None, draw_cell=True, origin=[0,0,0], species=None, spec_col=None, perspective=True, w_width=1200, w_height=750, bg_col=(0,0,0), bnd_col=(1,1,1), nx=1, ny=1, nz=1, coord_axes=False, boundary=True, bond_dists=0. ):
    '''
    Initialize the CrysPy object, creating a canvas and computing the corresponding lattice

    Arguments:
      qe_fname (str): Filename of a quantum espresso inputfile
      lattice (list): List of 3 3d vectors, representing the lattice parameters
      basis (list): List of N 3d vectors, representing the positions of each of N atoms
      draw_cell (bool): Draw the cell on creation of XCrysPy object
      origin (list): List of 3 points, representing the x,y,z position of the grid's center
      species (list): List of N species strings corresponding to atoms in 'basis' (e.g. ['Si','Si'])
      sepc_col (dict): Dictionary linking atom identfiers (e.g. 'Si' or 'He') to color tuples (R,G,B)
      perspective (bool): True for perspective view (with fov set for the camera)
      w_width (int): Vpython window width
      w_height (int): Vpython window height
      bg_col (tuple): Background color (R,G,B)
      bnd_col (tuple): Brillouin Zone boundary color (R,G,B)
      nx,ny,nz (int): Number of cells to draw in the x,y,z direction
      coord_axes (bool): Draw the cartesian axes for reference
      boundary (bool): True to draw the boundary of the lattice
      bond_dists (int or dict): Maximum bond distance or dictionary of form {'species':max_bond_dist}
    '''
    from .View import View
    from .Util import qe_lattice,read_scf_file,read_relax_file

    self.natoms = 0
    self.spec = None
    self.atoms = None
    self.lattice = None
    self.relax_index = 0
    self.coord_type = None
    self.eval_dist = False
    self.eval_angle = False
    self.relax_coords = None
    self.recip_space = False

    if qe_fname is None:
      self.coord_type = 'manual'
      if (lattice or basis or species) is None:
        print('Lattice and Basis not defined. Only \'plot_bxsf\' will function.')
      else:
        self.atoms = basis
        self.spec = species
        self.natoms = len(basis)
        self.lattice = np.array(lattice)
        self.cell_param = [np.abs(np.mean([np.linalg.norm(v) for v in self.lattice]))]
    else:
      if 'relax' in qe_fname:
        read_relax_file(self,qe_fname)
      else:
        read_scf_file(self,qe_fname)
      self.lattice = qe_lattice(self.ibrav, self.cell_param)

      if 'angstrom' in self.coord_type:
        conv = .529177210 # Reciprocal of (Bohr to Anstrom)
        self.atoms = np.array(self.atoms) / conv
      elif 'crystal' in self.coord_type:
        self.atoms = [self.atomic_position(a) for a in self.atoms]
        print('Crystal coordinate types require more testing')
      elif 'relax' in self.coord_type:
        if len(self.relax_lattices) > 0:
          self.lattice = self.relax_lattices[0]
        for i in range(len(self.relax_poss)):
          self.relax_poss[i][:,:] = np.array([self.atomic_position(a) for a in self.relax_poss[i]])
        self.atoms = self.relax_poss[0]

    if species is None and self.spec is None:
      self.spec = [i for i in range(self.natoms)]
    if spec_col is None:
      self.spec_col = {v:(1,1,1) for v in self.spec}
    else:
      self.spec_col = spec_col

    self.atom_radius,self.bond_radius = .7,.07
    title = 'CrysPy' if qe_fname is None else qe_fname

    self.canvas = vp.canvas(title=title+'\n', width=w_width, height=w_height, background=self.vector(bg_col))
    self.canvas.bind('click', self.click)

    anch = self.canvas.title_anchor
    self.disp_menu = vp.menu(choices=['Atoms', 'Bonds'], pos=anch, bind=self.disp_menu_change)
    self.sel_menu = vp.menu(choices=['Select', 'Distance', 'Angle'], pos=self.canvas.title_anchor, bind=self.sel_menu_change)
    self.sel_menu_text = vp.wtext()

    sel_nums = [str(i+1) for i in range(6)]
    self.sel_nx = vp.menu(choices=sel_nums, pos=anch, bind=self.sel_nx_cells, selected=str(nx))
    self.sel_ny = vp.menu(choices=sel_nums, pos=anch, bind=self.sel_ny_cells, selected=str(ny))
    self.sel_nz = vp.menu(choices=sel_nums, pos=anch, bind=self.sel_nz_cells, selected=str(nz))

    text = 'Draw Cell Boundaries'
    self.sel_bounary = vp.checkbox(text=text, pos=anch, bind=self.toggle_boundary, checked=boundary)

    text = 'Perspective View'
    self.sel_fov = vp.checkbox(text=text, pos=anch, bind=self.toggle_fov, checked=perspective)

    self.view = View(self.canvas,origin,perspective,bnd_col,nx,ny,nz,coord_axes,boundary,bond_dists)

    if draw_cell and self.lattice is not None:
      if 'relax' in self.coord_type:
        self.draw_relax()
      else:
        self.draw_cell(self.lattice, self.atoms)

  def disp_menu_change ( self, m ):
    if m.selected == 'Atoms':
      self.atom_radius = .7
      self.bond_radius = .07
    if m.selected == 'Bonds':
      self.atom_radius = .25
      self.bond_radius = .15
    self.draw_cell(self.lattice, self.atoms)

  def sel_menu_change ( self, m ):
    '''
    React to the atom selection menu selction changing

    Arguments:
      m (vpython widget): Menu widget
    '''
    if m.selected == 'Select':
      self.atom_select()
    elif m.selected == 'Distance':
      self.atom_dist()
    elif m.selected == 'Angle':
      self.atom_angle()
    else:
      raise ValueError('No such selection.')

  def toggle_boundary ( self, m ):
    self.view.boundary = not self.view.boundary
    self.draw_cell(self.lattice, self.atoms)

  def toggle_fov ( self, m ):
    self.view.canvas.fov = self.view.oFOV if m.checked else .01
    self.draw_cell(self.lattice, self.atoms)

  def sel_num_cells ( self, m, ind):
    self.view.cell_dim[ind] = int(m.selected)
    self.view.reset_selection()
    self.draw_cell(self.lattice, self.atoms)

  def sel_nx_cells ( self, m ):
    self.sel_num_cells(m, 0)

  def sel_ny_cells ( self, m ):
    self.sel_num_cells(m, 1)

  def sel_nz_cells ( self, m ):
    self.sel_num_cells(m, 2)

  def atom_select ( self ):
    '''
    Allow the selection of atoms
    '''
    self.view.reset_selection()
    self.eval_dist = False
    self.eval_angle = False

  def atom_dist ( self ):
    '''
    Allow the calculation of distance between atoms upon successive selection of two atoms
    '''
    self.view.reset_selection()
    self.eval_dist = True
    self.eval_angle = False

  def atom_angle ( self ):
    '''
    Allow the calculation of angle between atoms upon successive selection of three atoms
    '''
    self.view.reset_selection()
    self.eval_dist = False
    self.eval_angle = True

  def click ( self ):
    '''
    Handle mouse click events. Used to select atoms for distance calculation
    '''
    new_atom = self.canvas.mouse.pick
    if self.eval_dist:
      v = self.view.distance_selection(new_atom)
    elif self.eval_angle:
      v = self.view.angle_selection(new_atom)
    else:
      v = self.view.single_selection(new_atom)
    if v is not None:
      self.sel_menu_text.text = v

  def relax_step ( self, sgn ):
    '''
    Make a step forward or backward in relaxation

    Arguments:
      sgn (int): +1 => forward, -1 => backward
    '''

    # Ensure that the new index still indexes a valid relaxation step
    top = sgn > 0 and self.relax_index < len(self.relax_poss)-1
    bot = sgn < 0 and self.relax_index > 0

    if top or bot:
      self.relax_index += sgn
      self.relax_text.text = str(self.relax_index)
      self.atoms = self.relax_poss[self.relax_index]
      if len(self.relax_lattices) > 0:
        self.lattice = self.relax_lattices[self.relax_index]

      self.draw_cell(self.lattice, self.atoms)

  def relax_step_forward ( self ):
    self.relax_step(1)

  def relax_step_backward ( self ):
    self.relax_step(-1)

  def vector ( self, a ):
    '''
    Return list a as a vpython vector

    Aruments:
      a (list or ndarray): List or ndarray with 3 components (x,y,z)
    '''
    return vp.vector(a[0],a[1],a[2])

  def reciprocal_lattice ( self, lattice ):
    '''
    Construct reciprocal lattice vectors

    lattice (list or ndarray): 3 3d vectors of the real space lattice
    '''
    return 2*np.pi*np.linalg.inv(lattice).T

  def atomic_position ( self, v, lattice=None ):
    '''
    Calculate the atom's position within the unit cell

    Arguments:
      v (list or ndarray): List of 3 values, each to weight one of the 3 lattice vectors
      lattice (list or ndarray): 3 3d vectors of the real space lattice
    '''
    if lattice is None:
      lattice = self.lattice
    return np.sum(v*lattice, axis=1)

  def draw_cell ( self, lattice, atoms ):
    '''
    Draw the cell by calling the View's draw_cell method with the current cell configurations
    '''
    if lattice is not None and atoms is not None and not self.recip_space:
      self.view.draw_cell(lattice, atoms, self.spec, self.spec_col, self.atom_radius, self.bond_radius)

  def draw_relax ( self ):
    '''
    Start a relaxation visualization. Saves argument information for future drawings
    Creates Forward & Backward buttons & draws the first cell
    '''
    self.relax_backward = vp.button(text='<-', bind=self.relax_step_backward)
    self.relax_forward = vp.button(text='->', bind=self.relax_step_forward)
    self.relax_text = vp.wtext(text='0')

    self.draw_cell(self.relax_lattices[0], self.relax_poss[0])

  def draw_BZ_points ( self, points, color=None, rlat=None ):
    '''
    Draw points inside of the BZ

    Arguments:
      points (list or ndarray): List of points (x,y,z) to plot in the BZ, each between -Pi/2 and Pi/2
      color (tuple): Tuple of (R,G,B) with each between 0 & 1
      rlat (list or ndarray): 3 3d vectors representing the reciprocal lattice. If 'None' vectors will be computed from 'self.lattice'
    '''
    if rlat is None:
      rlat = self.rlattice = self.reciprocal_lattice(self.lattice)
    self.view.draw_BZ_points(points, color, rlat)

  def plot_spin_texture ( self, fermi_fname, spin_fname, e_up=1, e_dw=-1 ):
    '''
    Plots the spin texture read from the format output by PAOFLOW

    Arguments:
      fermi_fname (str): Name for Fermi surface file, representing energy eigenvalues in BZ for such band
      spin_fname (str): Name for spin texture file, representing spin direction in BZ for such band
      e_up (float): BZ points with energies higher than e_up wont be plotted
      e_dw (float): BZ points with energies lower than e_dw wont be plotted 
    '''
    from scipy.fftpack import fftshift

    self.recip_space = True
    fdata,sdata = np.load(fermi_fname),np.load(spin_fname)
    eig,spins = fdata['nameband'],sdata['spinband']

    eig = fftshift(eig,axes=(0,1,2))
    spins = np.moveaxis(np.real(spins[:,:,:,:]),spins.ndim-1,0)

    self.view.draw_arrows(spins, eig, self.reciprocal_lattice(self.lattice), e_up, e_dw)

  def plot_bxsf ( self, fname, iso=[0], bands=[0], colors=[[0,1,0]] ):
    '''
    Create the Brillouin Zone boundary and bsxf points between 'fermiup' and 'fermidw' in the vpython window

    Arguemnts:
      fname (str): Name of bxsf file
      iso (list): List of floats corresponding to the isosurface values for each respective band in 'bands'
      bands (list): List of integers representing the index of the band to plot
      colors (list): List of 3-d RGB color vectors for each band. If ignored, each band will be green
    '''
    from .Util import read_bxsf

    if len(iso) != len(bands):
      raise ValueError("Each band in 'bands' should have one corresponding isosurface in 'iso'")
    if len(iso) != len(colors) and len(colors) != 1:
      raise ValueError("Specify 1 color to plot all bands in the same color, or specify 1 color for each band.")

    self.recip_space = True
    b_vec,data = read_bxsf(fname)

    if np.max(bands) > data.shape[-1]:
      raise ValueError("'nbnd' too large to plot all bands in 'bands'")

    self.view.draw_bxsf(b_vec, data, iso, bands, colors)
