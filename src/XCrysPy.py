import vpython as vp
import numpy as np

class XCrysPy:

  def __init__ ( self, qe_fname=None, lattice=None, basis=None, origin=[0,0,0], species=None, spec_col=None, perspective=True, w_width=1200, w_height=750, bg_col=(0,0,0), bnd_col=(1,1,1) ):
    '''
    Initialize the CrysPy object, creating a canvas and computing the corresponding lattice

    Arguments:
      qe_fname (str): Filename of a quantum espresso inputfile
      lattice (list): List of 3 3d vectors, representing the lattice parameters
      basis (list): List of N 3d vectors, representing the positions of each of N atoms
      origin (list): List of 3 points, representing the x,y,z position of the grid's center
      species (list): List of N species strings corresponding to atoms in 'basis' (e.g. ['Si','Si'])
      sepc_col (dict): Dictionary linking atom identfiers (e.g. 'Si' or 'He') to color tuples (R,G,B)
      w_width (int): Vpython window width
      w_height (int): Vpython window height
      bg_col (tuple): Background color (R,G,B)
      bnd_col (tuple): Brillouin Zone boundary color (R,G,B)
    '''
    from .Util import qe_lattice,read_qe_file,crystal_conversion
    self.canvas = vp.canvas(title='CrysPy', width=w_width, height=w_height, background=self.vector(bg_col))
    self.dist_text = 'Distance (Disabled)'
    self.angle_text = 'Angle (Disabled)'
    self.dist_otext = '0.000000 angstroms'
    self.angle_otext = '0.000000 degrees'
    self.dist_button = vp.button(text=self.dist_text, bind=lambda:self.toggle_dist())
    self.dist_obutton = vp.button(text=self.dist_otext, bind=lambda:0)
    self.angle_button = vp.button(text = self.angle_text, bind=lambda:self.toggle_angle())
    self.angle_obutton = vp.button(text = self.angle_otext, bind=lambda:0)
    self.canvas.bind('click', self.click)

    if not perspective:
      self.canvas.fov = .01
    self.canvas.forward = vp.vector(0,+1,0)
    self.canvas.up = vp.vector(0,0,1)
    self.orient_lights()

    self.spec = None
    self.bonds = None
    self.vAtoms = None
    self.arrows = None
    self.BZ_corners = None
    self.coord_axes = None
    self.eval_dist = False
    self.eval_angle = False
    self.selected_atoms = []
    self.selected_bonds = []
    self.selected_colors = []
    self.origin = np.array(origin)

    self.bnd_col = self.vector(bnd_col) if bnd_col is not None else vp.vector(1,1,1)

    if qe_fname is None:
      if (lattice or basis or species) is None:
        print('Lattice and Basis not defined. Only \'plot_bxsf\' will function.')
        return
      self.atoms = basis
      self.spec = species
      self.natoms = len(basis)
      self.lattice = np.array(lattice)
    else:
      read_qe_file(self,qe_fname)
      self.lattice = qe_lattice(self.ibrav, self.cell_param)
      if self.coord_type == 'angstrom':
        self.atoms = crystal_conversion(self.atoms, self.lattice, self.coord_type)

    # Construct reciprocal lattice vectors
    self.rlattice = 2*np.pi*np.linalg.inv(self.lattice).T

    if species is None and self.spec is None:
      self.spec = [i for i in range(self.natoms)]
    if spec_col is None:
      self.specD = {v:(1,1,1) for v in self.spec}
    else:
      self.specD = spec_col

  def calc_dist ( self, a, b ):
      '''
      Calculate the distance between two atoms

      Arguments:
        a (Atom or vpython.sphere): First atom
        b (Atom or vpython.sphere): Second atom
      '''
      dist = b.pos - a.pos
      text = '%f angstroms'%dist.mag
      print('Distance = %s'%text)
      self.dist_obutton.text = text

  def calc_angle ( self, atoms ):
    '''
    Calculate the angle between three atoms. Angle between vectors created from indices 0-1 and 2-1

    Arguments:
      atoms (list): list of 3 vpython spheres
    '''
    v1 = atoms[0].pos - atoms[1].pos
    v2 = atoms[2].pos - atoms[1].pos
    angle = np.arccos(v1.dot(v2)/(v1.mag*v2.mag))
    text = '%f degrees'%np.degrees(angle)
    print('Angle = %s (%f radians)'%(text,angle))
    self.angle_obutton.text = text

  def toggle_dist ( self ):
    '''
    Toggle the calculation of distance between atoms upon successive selection of two atoms
    '''
    self.reset_selection()
    if self.eval_angle:
      self.toggle_angle()
    self.eval_dist = not self.eval_dist
    self.dist_button.text = 'Distance (%s)'%('Enabled' if self.eval_dist else 'Disabled')

  def toggle_angle ( self ):
    '''
    Toggle the calculation of angle between atoms upon successive selection of two atoms
    '''
    self.reset_selection()
    if self.eval_dist:
      self.toggle_dist()
    self.eval_angle = not self.eval_angle
    self.angle_button.text = 'Angle (%s)'%('Enabled' if self.eval_angle else 'Disabled')

  def select_atom ( self, atom ):
    '''
    Select and highlight atom

    Arguments:
      atom (vpython.sphere): Sphere object to select
    '''
    select_col = vp.vector(0,1,1)
    col = atom.color
    if len(self.selected_atoms) > 0:
      bpos1 = self.selected_atoms[-1].pos
      bpos2 = atom.pos
      self.selected_bonds.append(vp.curve({'pos':bpos1,'color':select_col},{'pos':bpos2,'color':select_col}))
    self.selected_atoms.append(atom)
    self.selected_colors.append(vp.vector(col.x, col.y, col.z))
    atom.color = select_col

  def pop_prev_atom ( self ):
    '''
    Deselect the most recently selected atom
    '''
    old_atom = self.selected_atoms.pop()
    old_atom.color = self.selected_colors.pop()
    if len(self.selected_bonds) > 0:
      bond = self.selected_bonds.pop()
      bond.visible = False
      del bond
    return old_atom

  def reset_selection ( self ):
    '''
    Unselect all selected atoms
    '''
    if len(self.selected_atoms) != 0:
      for i,a in enumerate(self.selected_atoms):
        a.color = self.selected_colors[i]
      for b in self.selected_bonds:
        b.visible = False
        del b
      self.selected_atoms = []
      self.selected_bonds = []
      self.selected_colors = []

  def click ( self ):
    '''
    Handle mouse click events. Used to select atoms for distance calculation
    '''
    new_atom = self.canvas.mouse.pick
    if self.eval_dist:
      if isinstance(new_atom, vp.sphere):
        if len(self.selected_atoms) == 2:
          self.reset_selection()
          self.select_atom(new_atom)
        elif len(self.selected_atoms) == 1:
          if new_atom in self.selected_atoms:
            self.pop_prev_atom()
          else:
            self.select_atom(new_atom)
            self.calc_dist(self.selected_atoms[0], self.selected_atoms[1])
        else:
          self.select_atom(new_atom)
    elif self.eval_angle:
      if isinstance(new_atom, vp.sphere):
        if len(self.selected_atoms) == 3:
          self.reset_selection()
          self.select_atom(new_atom)
        elif len(self.selected_atoms) == 2:
          if new_atom in self.selected_atoms:
            self.pop_prev_atom()
          else:
            self.select_atom(new_atom)
            self.calc_angle(self.selected_atoms)
        elif len(self.selected_atoms) == 1 and self.selected_atoms[0] == new_atom:
          self.pop_prev_atom()
        else:
          self.select_atom(new_atom)
    else:
        if isinstance(new_atom, vp.sphere):
          if len(self.selected_atoms) == 0:
            self.select_atom(new_atom)
          else:
            old_atom = self.pop_prev_atom()
            if old_atom != new_atom:
              self.select_atom(new_atom)
        else:
          if len(self.selected_atoms) > 0:
            self.reset_selection()

  def figure_cell_params ( self, ibrav ):
    pass

  def orient_lights ( self, direction=None ):
    '''
    Orient the lights to face the "front" of the sample
    '''
    if direction is None:
      for l in self.canvas.lights:
        cp = l.direction
        l.direction = vp.vector(cp.x,-cp.y,cp.z)

  def vector ( self, a ):
    '''
    Return list a as a vpython vector

    Aruments:
      a (list or ndarray): List or ndarray with 3 components (x,y,z)
    '''
    return vp.vector(a[0],a[1],a[2])

  def get_cell_pos ( self, ix, iy, iz ):
    '''
    Return the cartesian origin of the unit cell at index [ix,iy,iz]

    Arguments:
      ix (int): x index
      iy (int): y index
      iz (int): z index
    '''
    return self.vector(np.sum([v*self.lattice[i] for i,v in enumerate([ix,iy,iz])],axis=0))

  def atomic_position ( self, v, lattice=None ):
    '''
    Calculate the atom's position within the unit cell

    Arguments:
      v (list or ndarray): List of 3 values, each to weight one of the 3 lattice vectors
    '''
    if lattice is None:
      lattice = self.lattice
    return np.sum(v*lattice, axis=1)

  def clear_canvas ( self ):
    '''
    Resets the canvas to its default state
    '''
    self.reset_selection()
    self.dist_obutton.text = self.dist_otext

    def vpobject_destructor ( obj ):
      if obj is not None:
        for o in obj:
          o.visible = False
          del o
      return None

    if self.vAtoms is not None:
      self.vAtoms = vpobject_destructor([a.vpy_sph for a in self.vAtoms])
    self.bonds = vpobject_destructor(self.bonds)
    self.BZ_corners = vpobject_destructor(self.BZ_corners)
    self.coord_axes = vpobject_destructor(self.coord_axes)

  def draw_coord_axis ( self, offset=[-10,0,0], length=1. ):
    '''
    Draw the Coordinate Axes

    Arguments:
      length (int): Length of the coordiniate axis arrows
    '''
    apos = np.array(offset)
    vpa = lambda a : vp.arrow(pos=self.vector(apos),axis=self.vector(length*np.array(a)))
    self.coord_axes = [vpa([2,0,0]),vpa([0,2,0]),vpa([0,0,2])]

  def draw_bonds ( self, dist=1. ):
    '''
    Draw the bonds if atoms are closer than 'dist' together.
    To specifiy bonds between two types of atoms pass a dictionary into dist with key/value
    pairs '%s_%s':val, where %s are the respective species and val is the maximum bond distance.

    Arguments:
      dits (float or dict): A float describes global bond distance, while a dictionary can specify distances for various pairs of atoms
    '''
    if self.vAtoms is None or len(self.vAtoms) == 1:
      raise ValueError('Not enough Atoms have been drawn')

    if not isinstance(dist, dict):
      self.bonds = [vp.curve({'pos':a.pos,'color':a.col},{'pos':b.pos,'color':b.col}) for i,a in enumerate(self.vAtoms) for j,b in enumerate(self.vAtoms) if i!=j and (a.pos-b.pos).mag<=dist]
    else:
      self.bonds = []
      dist = {'%s_%s'%tuple(sorted(k.split('_'))):dist[k] for k in dist.keys()}
      for i,a in enumerate(self.vAtoms):
        for j,b in enumerate(self.vAtoms):
          if i != j:
            key = '%s_%s'%tuple(sorted([a.species,b.species]))
            if key in dist:
              if (a.pos-b.pos).mag <= dist[key]:
                self.bonds.append(vp.curve({'pos':a.pos,'color':a.col},{'pos':b.pos,'color':b.col}))

  def draw_cell ( self, nx=1, ny=1, nz=1, boundary=False ):
    '''
    Create the cell simulation in the vpython environment

    Arguments:
      nx (int): Number of cells to draw in x direction
      ny (int): Number of cells to draw in y direction
      nz (int): Number of cells to draw in z direction
      boundary (bool): Draw cell boundaries
    '''
    from .Atom import Atom
    if boundary:
      lines = []
      lat = self.lattice
      for ix in range(nx):
        for iy in range(ny):
          for iz in range(nz):
            corner = np.sum(lat, axis=0)
            alines = [[3*[0],a] for a in lat]
            alines += [[p,corner] for p in [lat[i]+lat[j] for i in range(3) for j in range(i+1,3)]]
            alines += [[lat[k],lat[i]+lat[j]] for i in range(3) for j in range(i+1,3) for k in [i,j]]
            orig = self.get_cell_pos(ix,iy,iz)
            lines += [[orig+self.vector(al[0]),orig+self.vector(al[1])] for al in alines]
      lines = [vp.curve(l[0],l[1],color=self.bnd_col) for l in lines]

    self.vAtoms = []
    for ix in range(nx):
      for iy in range(ny):
        for iz in range(nz):
          c_pos = self.vector(self.origin) + self.get_cell_pos(ix,iy,iz)
          for i,a in enumerate(self.atoms):
            a_pos = c_pos + self.vector(self.atomic_position(a))
            color = self.vector(self.specD[self.spec[i]])
            self.vAtoms.append(Atom(a_pos,col=color,species=self.spec[i]))

    self.canvas.center = self.vector(np.mean([[v.pos.x,v.pos.y,v.pos.z] for v in self.vAtoms], axis=0))
    self.draw_coord_axis()

  def draw_BZ_boundary ( self, b_vec=None ):
    '''
    Draw the Brillouin Zone boundary in the vpython window

    Arguments:
      b_vec (list or ndarray): 3 3-d vectors representing the reciprocal lattice vectors
    '''
    from .Util import bravais_boundaries
    if b_vec is None:
      b_vec = self.rlattice

    self.BZ_curves = []
    self.BZ_planes,corners = bravais_boundaries(b_vec)
    for c in corners:
      self.BZ_curves.append(vp.curve(self.vector(c[0]), self.vector(c[1]), color=self.bnd_col))

  def plot_spin_texture ( self, fermi_fname, spin_fname, e_up=10, e_dw=-10 ):
    '''
    Plots the spin texture read from the format output by PAOFLOW

    Arguments:
      fermi_fname (str): Name for Fermi surface file, representing energy eigenvalues in BZ for such band
      spin_fname (str): Name for spin texture file, representing spin direction in BZ for such band
      e_up (float): BZ points with energies higher than e_up wont be plotted
      e_dw (float): BZ points with energies lower than e_dw wont be plotted 
    '''
    from scipy.fftpack import fftshift

    fdata, sdata = np.load(fermi_fname), np.load(spin_fname)
    eig, spin = fdata['nameband'], sdata['spinband']
    nx,ny,nz,_ = spin.shape

    eig = fftshift(eig,axes=(0,1,2))
    S = np.moveaxis(np.real(spin[:,:,:,:]),spin.ndim-1,0)

    self.draw_BZ_boundary(self.rlattice)

    self.arrows = []
    colscale = np.mean(eig)
    vp_shift = .5 * np.array([1,1,0])
    for x in range(nx):
      for y in range(ny):
        for z in range(nz):
          bv = eig[x,y,z]
          if bv > e_dw and bv < e_up:
            point = self.atomic_position([x/nx,y/ny,z/nz]-vp_shift, lattice=self.rlattice)
            self.arrows.append(vp.arrow(pos=self.vector(point),axis=self.vector(S[:,x,y,z]), length=.01, color=self.vector([bv/colscale,.5,.1])))


  def plot_bxsf ( self, fname, iso=[0], bands=[0], colors=[[0,1,0]] ):
    '''
    Create the Brillouin Zone boundary and bsxf points between 'fermiup' and 'fermidw' in the vpython window

    Arguemnts:
      fname (str): Name of bxsf file
      iso (list): List of floats corresponding to the isosurface values for each respective band in 'bands'
      bands (list): List of integers representing the index of the band to plot
      colors (list): List of 3-d RGB color vectors for each band. If ignored, each band will be green
    '''
    from numpy.linalg import det,norm,solve
    from scipy.fftpack import fftshift
    from .Util import read_bxsf

    if len(iso) != len(bands):
      raise ValueError("Each band in 'bands' should have one corresponding isosurface in 'iso'")

    b_vec,data = read_bxsf(fname)
    nx,ny,nz,nbnd = data.shape

    if np.max(bands) > nbnd:
      raise ValueError("'nbnd' too large to plot all bands in 'bands'")

    draw_flag = vp.text(text='Drawing...', align='center', color=vp.vector(1,0,0))
    draw_flag.up = self.canvas.up

    self.draw_BZ_boundary(b_vec)

    poss = []
    data = fftshift(data,axes=(0,1,2))
    vp_shift = .5 * (np.array([1,1,1]) - [1/nx,1/ny,1/nz])

    def inside_BZ ( pnt ):
      # Helper function to test whether points lie within the BZ
      pmag = np.linalg.norm(pnt)
      for p in self.BZ_planes:
        if pmag > np.linalg.norm(p):
          return False
      return True

    # Test points between x & x+1 to see whether they house an 'iso' point
    for i,b in enumerate(bands):
      bpos = []
      for x in range(nx-1):
        for y in range(ny-1):
          for z in range(nz-1):
            dA = data[x,y,z][b]
            dBL = [data[x+v[0],y+v[1],z+v[2]][b] for v in [(1,0,0),(0,1,0),(0,0,1)]]
            point = [x/nx,y/ny,z/nz] - vp_shift
            point = np.sum([[v] for v in point]*b_vec,axis=0)
            if dA == iso[i]:
              if inside_BZ(point):
                bpos.append(self.vector(point))
            for dB in dBL:
              if (dA > iso[i] and dB < iso[i]) or (dA < iso[i] and dB > iso[i]):
                if inside_BZ(point):
                  bpos.append(self.vector(point))
      if len(bpos) > 0:
        poss.append(bpos)

    self.vAtoms = []
    static_color = len(colors) < len(bands)
    srad = np.min([.5/n for n in (nx,ny,nz)])
    for i,p in enumerate(poss):
      col = colors[0] if static_color else colors[i]
      self.vAtoms.append(vp.points(pos=p,color=self.vector(col)))

    self.draw_coord_axis(length=.1*np.linalg.norm(b_vec[0]),offset=[-1,0,0])

    draw_flag.visible = False
    del draw_flag
