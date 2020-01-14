import vpython as vp
import numpy as np

class View:

  def __init__ ( self, canvas, origin, perspective, bnd_col, nx, ny, nz, coord_axes, boundary, bond_dists, bond_thickness, atom_radii):
    '''
    Initialize the CrysPy object, creating a canvas and computing the corresponding lattice

    Arguments:
      canvas (vpython.scene): The vpython window
      origin (list): List of 3 points, representing the x,y,z position of the grid's center
      perspective (bool): Draw the view in perspective mode if True
      bg_col (tuple): Background color (R,G,B)
      bnd_col (tuple): Brillouin Zone or Cell boundary color (R,G,B)
      nx (int): Number of cells to draw in x direction
      ny (int): Number of cells to draw in y direction
      nz (int): Number of cells to draw in z direction
      coord_axes (bool): Draw the coordinate system
      boundary (bool): Draw cell boundaries
      bond_thickness (float): Thickness of the drawn bonds
      atom_radii: (float or dict): Float to specify atom radius, or dict of form {'Species':radius}
    '''

    self.canvas = canvas
    self.boundary = boundary
    self.cell_dim = [nx,ny,nz]
    self.coord_axes = coord_axes
    self.bond_dists = bond_dists
    self.atom_radii = atom_radii
    self.origin = np.array(origin)
    self.perspective = perspective
    self.bnd_thck = bond_thickness
    self.bnd_col = self.vector(bnd_col)

    self.oFOV = self.canvas.fov
    if not perspective:
      self.canvas.fov = .01
    self.canvas.forward = vp.vector(0,+1,0)
    self.canvas.up = vp.vector(0,0,1)
    self.orient_lights()

    self.bonds = None
    self.arrows = None
    self.vAtoms = None
    self.BZ_planes = None
    self.coord_axis = None
    self.bound_curve = None
    self.bond_radius = None

    self.selected_atoms = []
    self.selected_bonds = []
    self.selected_colors = []

  
  def orient_lights ( self, direction=None ):
    '''
    Orient the lights to face the "front" of the sample

    Arguments:
      direction (list or ndarray): (x,y,z) direction for the light
    '''
    if direction is None:
      for l in self.canvas.lights:
        cp = l.direction
        l.direction = vp.vector(cp.x,-cp.y,cp.z)
    else:
      l.direction = self.vector(direction)

  def vector ( self, a ):
    '''
    Return list a as a vpython vector

    Aruments:
      a (list or ndarray): List or ndarray with 3 components (x,y,z)
    '''
    return vp.vector(a[0],a[1],a[2])
  
  def calc_dist ( self, a, b ):
      '''
      Calculate the distance between two atoms

      Arguments:
        a (Atom or vpython.sphere): First atom
        b (Atom or vpython.sphere): Second atom
      '''
      dist = b.pos - a.pos
      text = '%f angstroms'%(dist.mag * .529177)
      print('Distance = %s'%text)
      return text

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
    return text

  def get_cell_pos ( self, lat, ix, iy, iz ):
    '''
    Return the cartesian origin of the unit cell at index [ix,iy,iz]

    Arguments:
      lat (list or ndarray): Lattice vectors
      ix (int): x index
      iy (int): y index
      iz (int): z index
    '''
    return self.vector(np.dot([ix,iy,iz],lat))

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

  def single_selection ( self, atom ):
    '''
    Select a single atom at a time

    Arguments:
      atom (vpython object): The selected object
    '''
    if isinstance(atom, vp.sphere):
      if len(self.selected_atoms) == 0:
        self.select_atom(atom)
      else:
        old_atom = self.pop_prev_atom()
        if old_atom != atom:
          self.select_atom(atom)
    else:
      if len(self.selected_atoms) > 0:
        self.reset_selection()

  def angle_selection ( self, atom ):
    '''
    Select up to three atoms (1 2 & 3). Angle between 1_2 and 2_3 calculated after third is selected

    Arguments:
      atom (vpython object): The selected object
    '''
    if isinstance(atom, vp.sphere):
      if len(self.selected_atoms) == 3:
        self.reset_selection()
        self.select_atom(atom)
      elif len(self.selected_atoms) == 2:
        if atom in self.selected_atoms:
          self.pop_prev_atom()
        else:
          self.select_atom(atom)
          return self.calc_angle(self.selected_atoms)
      elif len(self.selected_atoms) == 1 and self.selected_atoms[0] == atom:
        self.pop_prev_atom()
      else:
        self.select_atom(atom)

  def distance_selection ( self, atom ):
    '''
    Select a up to two atoms (1 & 2). Distance between 1_2 calculated after second is selected

    Arguments:
      atom (vpython object): The selected object
    '''
    if isinstance(atom, vp.sphere):
      if len(self.selected_atoms) == 2:
        self.reset_selection()
        self.select_atom(atom)
      elif len(self.selected_atoms) == 1:
        if atom in self.selected_atoms:
          self.pop_prev_atom()
        else:
          self.select_atom(atom)
          return self.calc_dist(self.selected_atoms[0], self.selected_atoms[1])
      else:
        self.select_atom(atom)

  def select_atom ( self, atom ):
    '''
    Highlight atom and draw connections to previously selected atoms

    Arguments:
      atom (vpython.sphere): Sphere object to select
    '''
    select_col = vp.vector(0,1,1)
    col = atom.color
    if len(self.selected_atoms) > 0:
      dSel = lambda p : {'pos':p,'color':select_col,'radius':self.bond_radius}
      bpos1 = self.selected_atoms[-1].pos
      bpos2 = atom.pos
      self.selected_bonds.append(vp.curve(dSel(bpos1),dSel(bpos2)))
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

  def clear_canvas ( self ):
    '''
    Resets the canvas to its default state
    '''
    self.canvas.center = self.vector(self.origin)
    self.reset_selection()

    def vpobject_destructor ( obj ):
      if obj is not None:
        for o in obj:
          o.visible = False
          del o
      return None

    if self.vAtoms is not None:
      self.vAtoms = vpobject_destructor([a.vpy_sph for a in self.vAtoms])
    self.bonds = vpobject_destructor(self.bonds)
    self.coord_axis = vpobject_destructor(self.coord_axis)
    self.bound_curve = vpobject_destructor(self.bound_curve)

  def draw_coord_axis ( self, offset=[-10,0,0], length=1. ):
    '''
    Draw the Coordinate Axes

    Arguments:
      length (int): Length of the coordiniate axis arrows
    '''
    apos = np.array(offset)
    vpa = lambda a : vp.arrow(pos=self.vector(apos),axis=self.vector(length*np.array(a)))
    self.coord_axis = [vpa([2,0,0]),vpa([0,2,0]),vpa([0,0,2])]

  def draw_bonds ( self ):
    '''
    Draw the bonds if atoms are closer than 'dist' together.
    To specifiy bonds between two types of atoms pass a dictionary into dist with key/value
    pairs '%s_%s':val, where %s are the respective species and val is the maximum bond distance.

    Arguments:
      dits (float or dict): A float describes global bond distance, while a dictionary can specify distances for various pairs of atoms
    '''
    dist = self.bond_dists
    dList = lambda a : {'pos':a.pos, 'color':a.col, 'radius':self.bond_radius}
    if not isinstance(dist, dict):
      self.bonds = [vp.curve(dList(a),dList(b)) for i,a in enumerate(self.vAtoms) for j,b in enumerate(self.vAtoms) if i!=j and (a.pos-b.pos).mag<=dist]
    else:
      self.bonds = []
      dist = {'%s_%s'%tuple(sorted(k.split('_'))):dist[k] for k in dist.keys()}
      for i,a in enumerate(self.vAtoms):
        for j,b in enumerate(self.vAtoms):
          if i != j:
            key = '%s_%s'%tuple(sorted([a.species,b.species]))
            if key in dist:
              if (a.pos-b.pos).mag <= dist[key]:
                self.bonds.append(vp.curve(dList(a),dList(b)))

  def draw_cell ( self, lat, atoms, species, spec_col, atom_radius, bond_radius ):
    '''
    Create the cell simulation in the vpython environment

    Arguments:
      lat (list or ndarray): Lattice vectors
      atoms (list or ndarray): Atoms to draw on lattice
    '''
    from .Atom import Atom

    self.clear_canvas()
    nx,ny,nz = self.cell_dim
    self.bond_radius = bond_radius if self.bnd_thck is None else self.bnd_thck

    atom_radius = atom_radius if self.atom_radii is None else self.atom_radii
    if not isinstance(atom_radius, dict):
      rad = atom_radius
      atom_radius = {s:rad for s in species}

    if self.boundary:
      self.bound_curve = []
      for ix in range(nx):
        for iy in range(ny):
          for iz in range(nz):
            corner = np.sum(lat, axis=0)
            alines = [[3*[0],a] for a in lat]
            alines += [[p,corner] for p in [lat[i]+lat[j] for i in range(3) for j in range(i+1,3)]]
            alines += [[lat[k],lat[i]+lat[j]] for i in range(3) for j in range(i+1,3) for k in [i,j]]
            orig = self.get_cell_pos(lat,ix,iy,iz)
            self.bound_curve += [vp.curve(orig+self.vector(al[0]),orig+self.vector(al[1]),color=self.bnd_col) for al in alines]

    self.vAtoms = []
    for ix in range(nx):
      for iy in range(ny):
        for iz in range(nz):
          c_pos = self.vector(self.origin) + self.get_cell_pos(lat,ix,iy,iz)
          for i,a in enumerate(atoms):
            a_pos = c_pos + self.vector(a)
            color = self.vector(spec_col[species[i]])
            self.vAtoms.append(Atom(a_pos,col=color,species=species[i],radius=atom_radius[species[i]]))

    if len(self.vAtoms) > 1:
      self.draw_bonds()

    self.canvas.center = self.vector(np.mean([[v.pos.x,v.pos.y,v.pos.z] for v in self.vAtoms], axis=0))
    if self.coord_axes:
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

    self.bound_curve = []
    self.BZ_planes,corners = bravais_boundaries(b_vec)
    for c in corners:
      self.bound_curve.append(vp.curve(self.vector(c[0]), self.vector(c[1]), color=self.bnd_col))

  def draw_BZ_points ( self, points, color, rlat ):
    '''
    Draw points inside of the BZ

    Arguments:
      points (list or ndarray): List of points (x,y,z) to plot in the BZ, each between -Pi/2 and Pi/2
      color (tuple): Tuple of (R,G,B) with each between 0 & 1
      rlat (list or ndarray): 3 3d vectors representing the reciprocal lattice
    '''
    self.clear_canvas()

    if self.BZ_planes is None:
      self.draw_BZ_boundary(b_vec=rlat)

    col = vp.vector(1,1,1) if color is None else self.vector(color)
    points = [self.vector(p) for p in points]
    if len(points) > 0:
      self.points_BZ = vp.points(pos=points, color=col)
    

  def draw_arrows ( self, rlat, points, directions, colors, size ):
    '''
    Draw the arrows for spin texture plots

    Arugments:
      rlat (list or ndarray): 3 3d vectors representing the reciprocal lattice
      points (list or ndarray): List of 3d vectors with positions of points (x,y,z) to plot
      directions (list or ndarray): List of 3d vectors with directions of arrows for each point.
      colors (list or ndarray): List of 3d vectors representing colors (R,G,B) for each arrow
      rlat (list or ndarray): 3 3d vectors representing the reciprocal lattice
      size (float): Size of the arrow. Default is .05
    '''
    self.clear_canvas()

    draw_flag = vp.text(text='Drawing...', align='center', color=vp.vector(1,0,0))
    draw_flag.up = self.canvas.up

    self.draw_BZ_boundary(rlat)

    self.arrows = []
    col = self.vector([.8,.8,0])
    for i,p in enumerate(points):
      pnt = self.vector(p)
      adir = self.vector(directions[i])
      if colors is not None:
        col = self.vector(colors[i])
      self.arrows.append(vp.arrow(pos=pnt,axis=adir, length=size, color=col))

    draw_flag.visible = False
    del draw_flag

  def draw_bxsf ( self, rlat, data, iso, bands, colors ):
    '''
    Create the Brillouin Zone boundary and bsxf points between 'fermiup' and 'fermidw' in the vpython window

    Arguemnts:
      rlat (list or ndarray): 3 3d vectors representing the reciprocal lattice
      data (ndarray): eigenvalue data for each band and point in the BZ to plot
      iso (list): List of floats corresponding to the isosurface values for each respective band in 'bands'
      bands (list): List of integers representing the index of the band to plot
      colors (list): List of 3-d RGB color vectors for each band. If ignored, each band will be green
    '''
    from numpy.linalg import det,norm,solve
    from scipy.fftpack import fftshift

    draw_flag = vp.text(text='Drawing...', align='center', color=vp.vector(1,0,0))
    draw_flag.up = self.canvas.up

    self.draw_BZ_boundary(rlat)

    poss = []
    nx,ny,nz,nbnd = data.shape
    data = fftshift(data, axes=(0,1,2))
    vp_shift = .5 * (np.array([1,1,1]) - [1/nx,1/ny,1/nz])

    # Helper function to test whether points lie within the BZ
    def inside_BZ ( pnt ):
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
            point = np.dot(point,rlat)
            if dA == iso[i]:
              if inside_BZ(point):
                bpos.append(self.vector(point))
            for dB in dBL:
              if (dA > iso[i] and dB < iso[i]) or (dA < iso[i] and dB > iso[i]):
                if inside_BZ(point):
                  bpos.append(self.vector(point))
      if len(bpos) > 0:
        poss.append(bpos)

    # Plot the points from each band
    self.vAtoms = []
    srad = np.min([.5/n for n in (nx,ny,nz)])
    for i,p in enumerate(poss):
      col = colors[0] if len(colors)<2 else colors[i]
      self.vAtoms.append(vp.points(pos=p,color=self.vector(col)))

    if self.coord_axes:
      self.draw_coord_axis(length=.1*np.linalg.norm(rlat[0]),offset=[-1,0,0])

    draw_flag.visible = False
    del draw_flag
