from ReadQE import read_qe_file,qe_lattice,crystal_conversion
from Atom import Atom
import vpython as vp
import numpy as np

class CrysPy:

  def __init__ ( self, qe_fname=None, lattice=None, basis=None, origin=[0,0,0], species=None, spec_col=None, w_width=1000, w_height=700 ):
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
    '''

    self.canvas = vp.canvas(title='CrysPy', width=w_width, height=w_height, background=vp.color.black)
    self.dist_button = vp.button(text='Distance (Disabled)', bind=lambda:self.toggle_dist())
    self.text_button = vp.button(text='0.000000 angstroms', bind=lambda:0)
    self.canvas.bind('click', self.click)

    self.canvas.forward = vp.vector(0,+1,0)
    self.canvas.up = vp.vector(0,0,1)
    self.orient_lights()

    self.spec = None
    if qe_fname is None:
      self.atoms = basis
      self.natoms = len(basis)
      self.lattice = np.array(lattice)
    else:
      read_qe_file(self,qe_fname)
      self.lattice = qe_lattice(self.ibrav, self.cell_param)
      if self.coord_type == 'angstrom':
        self.atoms = crystal_conversion(self.atoms, self.lattice, self.coord_type)
      

    if spec_col is None:
      self.spec = [i for i in range(self.natoms)]
      self.specD = {v:(1,1,1) for v in self.spec}
    else:
      if self.spec is None:
        self.spec = species
      self.specD = spec_col
      assert(len(self.spec) == self.natoms)

    self.vAtoms = None
    self.eval_dist = False
    self.selected_atoms = []
    self.selected_colors = []
    self.origin = np.array(origin)

    # Construct Recip Lattice
    self.rlattice = []

  def calc_dist ( self, a, b ):
      '''
      Calculate the distance between two atoms

      Arguments:
        a (Atom or vpython.sphere): First atom
        b (Atom or vpython.sphere): Second atom
      '''
      dist = b.pos - a.pos
      text = '%f angstroms'%dist.mag
      print(text)
      self.text_button.text = text

  def toggle_dist ( self ):
    '''
    Toggle the calculation of distance between atoms upon successive selection of two atoms
    '''
    self.reset_selection()
    self.eval_dist = not self.eval_dist
    self.dist_button.text = 'Distance (%s)'%('Enabled' if self.eval_dist else 'Disabled')

  def select_atom ( self, atom ):
    '''
    Select and highlight atom

    Arguments:
      atom (vpython.sphere): Sphere object to select
    '''
    select_col = vp.vector(0,1,1)
    col = atom.color
    self.selected_atoms.append(atom)
    self.selected_colors.append(vp.vector(col.x, col.y, col.z))
    atom.color = select_col

  def reset_selection ( self ):
    '''
    Unselect all selected atoms
    '''
    if len(self.selected_atoms) != 0:
      for i,a in enumerate(self.selected_atoms):
        a.color = self.selected_colors[i]
      self.selected_atoms = []
      self.selected_colors = []

  def click ( self ):
    '''
    Handle mouse click events. Used to select atoms for distance calculation
    '''
    new_atom = self.canvas.mouse.pick
    is_sphere = lambda v : isinstance(v, vp.sphere)
    if self.eval_dist:
      if is_sphere(new_atom):
        if len(self.selected_atoms) == 2:
          self.reset_selection()
          self.select_atom(new_atom)
        elif len(self.selected_atoms) == 1:
          if new_atom == self.selected_atoms[0]:
            old_atom = self.selected_atoms.pop()
            old_atom.color = self.selected_colors.pop()
          else:
            self.select_atom(new_atom)
            self.calc_dist(self.selected_atoms[0], self.selected_atoms[1])
        else:
          self.select_atom(new_atom)
    else:
        if len(self.selected_atoms) > 0:
          self.reset_selection()
        if is_sphere(new_atom):
          self.select_atom(new_atom)

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

  def draw_brillouin_zone ( self ):
    planes = []
    poss = []
    planes.append(vp.box(length=.01,height=4,width=4,pos=vp.vector(0,0,0),axis=vp.vector(1,0,0)))
    quit()

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
      for i,a in enumerate(self.vAtoms):
        for j,b in enumerate(self.vAtoms):
          if i != j:
            key = None
            tk = '%s_%s'%(a.species,b.species)
            if tk in dist:
              key = tk
            tk = '%s_%s'%(b.species,a.species)
            if tk in dist:
              key = tk
            if key is not None:
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
      lines = [vp.curve(l[0],l[1]) for l in lines]

    self.vAtoms = []
    for ix in range(nx):
      for iy in range(ny):
        for iz in range(nz):
          c_pos = self.vector(self.origin) + self.get_cell_pos(ix,iy,iz)
          for i,a in enumerate(self.atoms):
            a_pos = c_pos + self.vector(self.atomic_position(a))
            color = self.vector(self.specD[self.spec[i]])
            self.vAtoms.append(Atom(a_pos,col=color,species=self.spec[i]))

    apos = [np.min([v.pos.x for v in self.vAtoms])-5, 0, 0]
    vpa = lambda a : vp.arrow(pos=self.vector(apos),axis=self.vector(a))
    coord = [vpa([2,0,0]),vpa([0,2,0]),vpa([0,0,2])]
    self.canvas.center = self.vector(np.mean([[v.pos.x,v.pos.y,v.pos.z] for v in self.vAtoms], axis=0))

  def atomic_position ( self, v ):
    '''
    Calculate the atom's position within the unit cell

    Arguments:
      v (list or ndarray): List of 3 values, each to weight one of the 3 lattice vectors
    '''
    return np.sum(v*self.lattice, axis=1)
