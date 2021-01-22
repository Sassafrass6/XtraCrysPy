import vpython as vp
import numpy as np

class View:

  def __init__ ( self, cname, w_width, w_height, origin, model, perspective, bnd_col, bg_col, nx, ny, nz ):
    '''
    Initialize the CrysPy object, creating a canvas and computing the corresponding lattice

    Arguments:
      canvas (vpython.scene): The vpython window
      cname (str): None or the String name of the canvas
      w_width (int): Width of the window
      w_height (int): Height of the window
      origin (ndarray): 3-vector array, the system origin.
      model (list): Model components: [lattice, atoms, species, atom labels, bonds]
      perspective (bool): Draw the view in perspective mode if True
      bg_col (tuple): Background color (R,G,B). None defaults to white
      bnd_col (tuple): Brillouin Zone or Cell boundary color (R,G,B)
      nx (int): Number of cells to draw in x direction
      ny (int): Number of cells to draw in y direction
      nz (int): Number of cells to draw in z direction
    '''

    self.model = model       # Model, containing lattice, species, and bond information
    self.relax_poss = None   # Atomic positions for each relax step
    self.relax_index = None  # Index of the current relax step being modeled

    if model != None:
      self.bonds = model[4]
      self.species = model[2]
      self.basis_labels = model[3]
      if len(self.model[1].shape) == 3:
        self.relax_index = 0
        self.relax_poss = self.model[1]
        self.atoms = self.relax_poss[0]
        self.relax_lattices = self.model[0]
        self.lattice = self.relax_lattices[0]
      else:
        self.atoms = model[1]
        self.lattice = model[0]
      self.natoms = len(self.atoms)
    
    self.cell_dim = [nx,ny,nz]
    self.origin = np.array(origin)
    self.perspective = perspective
    self.bnd_thck = None
    self.bnd_col = bnd_col

    self.unit_anch = None
    self.unit_menu = None
    self.eval_dist = False
    self.eval_angle = False
    self.recip_space = False

    self.boundary = True
    if self.bnd_col is None:
      self.boundary = False
      self.bnd_col = (1,1,1)
    self.bnd_col = self.vector(self.bnd_col)

    self.arrows = None
    self.vAtoms = None
    self.vBonds = None
    self.BZ_planes = None
    self.coord_axes = None
    self.bound_curve = None
    self.bond_radius = None

    self.selected_atoms = []
    self.selected_bonds = []
    self.selected_colors = []

    self.setup_canvas(cname, w_width, w_height, bg_col, nx, ny, nz, perspective)

    self.oFOV = self.canvas.fov
    if not perspective:
      self.canvas.fov = .01
    self.canvas.forward = vp.vector(0,+1,0)
    self.canvas.up = vp.vector(0,0,1)
    self.orient_lights()

    if model != None:
      self.draw_cell()

  def setup_canvas ( self, cname, w_width, w_height, bg_col, nx, ny, nz, perspective):
    '''
    Create the Canvas object and initialize the caption text and buttons.

    Arguments:
      cname (str): Filename of a quantum espresso inputfile
      draw_cell (bool): Draw the cell on creation of XtraCrysPy object
      w_width (int): Vpython window width
      w_height (int): Vpython window height
      bg_col (tuple): Background color (R,G,B)
      bnd_col (tuple): Brillouin Zone boundary color (R,G,B)
      nx,ny,nz (int): Number of cells to draw in the x,y,z direction
      boundary (bool): True to draw the boundary of the lattice
      perspective (bool): Flag to set the FOV as perspective mode
    '''
    self.atom_radius,self.bond_radius = .7,.07
    title = '\tXtraCrysPy' if cname=='' else cname

    self.canvas = vp.canvas(title=title+'\n', width=w_width, height=w_height, background=self.vector(bg_col))
    self.canvas.bind('click', self.click)

    anch = self.canvas.caption_anchor
    self.canvas.caption = '\nOptions:\t\t\t\tTools:\n'

    text = 'Cell Boundaries'
    self.sel_bounary = vp.checkbox(text=text, pos=anch, bind=self.toggle_boundary, checked=self.boundary)

    self.canvas.append_to_caption('     \t')
    self.disp_menu = vp.menu(choices=['Atom Primary','Bond Primary'], pos=anch, bind=self.sel_disp_menu)
    self.sel_menu = vp.menu(choices=['Select Atom','Distance','Angle'], pos=anch, bind=self.sel_menu_change)
    self.unit_menu = vp.menu(choices=['Units','Angstrom','Bohr','Degree','Radian'], pos=anch, bind=lambda m:None, selected='Units')

    if self.relax_poss is not None:
      self.canvas.append_to_caption('  \t')
      self.relax_backward = vp.button(text='<-', bind=self.relax_step_backward)
      self.relax_forward = vp.button(text='->', bind=self.relax_step_forward)
      self.relax_text = vp.wtext(text='Step: 0')

    text = 'Perspective View'
    self.canvas.append_to_caption('\n')
    self.sel_fov = vp.checkbox(text=text, pos=anch, bind=self.toggle_fov, checked=perspective)
    self.canvas.append_to_caption('\t\t')
    self.sel_menu_text = vp.wtext()

    self.canvas.append_to_caption('\n\nNumber of Cells:\n   Nx      Ny      Nz\n')

    sel_nums = [str(i+1) for i in range(6)]
    self.sel_nx = vp.menu(choices=sel_nums, pos=[2,0], bind=self.sel_nx_cells, selected=str(nx))
    self.sel_ny = vp.menu(choices=sel_nums, pos=anch, bind=self.sel_ny_cells, selected=str(ny))
    self.sel_nz = vp.menu(choices=sel_nums, pos=anch, bind=self.sel_nz_cells, selected=str(nz))

  def sel_disp_menu ( self, m ):
    if m.selected == 'Atom Primary':
      self.atom_radius = .7
      self.bond_radius = .07
    if m.selected == 'Bond Primary':
      self.atom_radius = .25
      self.bond_radius = .25
    self.draw_cell()

  def sel_menu_change ( self, m ):
    '''
    React to the atom selection menu selction changing

    Arguments:
      m (vpython widget): Menu widget
    '''
    if 'Select' in m.selected:
      self.atom_select()
    elif 'Distance' in m.selected:
      self.atom_dist()
    elif 'Angle' in m.selected:
      self.atom_angle()
    else:
      raise ValueError('No such selection.')

  def toggle_boundary ( self, m ):
    self.boundary = not self.boundary
    self.draw_cell()

  def toggle_fov ( self, m ):
    self.canvas.fov = self.oFOV if m.checked else .01
    self.draw_cell()

  def sel_num_cells ( self, m, ind):
    self.cell_dim[ind] = int(m.selected)
    self.reset_selection()
    self.draw_cell()

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
    self.reset_selection()
    self.eval_dist = False
    self.eval_angle = False

  def atom_dist ( self ):
    '''
    Allow the calculation of distance between atoms upon successive selection of two atoms
    '''
    self.reset_selection()
    self.eval_dist = True
    self.eval_angle = False
    unit = self.unit_menu.selected
    if unit != 'Angstrom' and unit != 'Bohr':
      self.unit_menu.selected = 'Bohr'

  def atom_angle ( self ):
    '''
    Allow the calculation of angle between atoms upon successive selection of three atoms
    '''
    self.reset_selection()
    self.eval_dist = False
    self.eval_angle = True
    unit = self.unit_menu.selected
    if unit != 'Degree' and unit != 'Radian':
      self.unit_menu.selected = 'Degree'

  def click ( self ):
    '''
    Handle mouse click events. Used to select atoms for distance calculation
    '''
    new_atom = self.canvas.mouse.pick
    if self.eval_dist:
      v = self.distance_selection(new_atom)
    elif self.eval_angle:
      v = self.angle_selection(new_atom)
    else:
      v = self.single_selection(new_atom)
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
      self.relax_text.text = 'Step: %d'%self.relax_index
      self.atoms = self.relax_poss[self.relax_index]
      if len(self.relax_lattices) > 0:
        self.lattice = self.relax_lattices[self.relax_index]

      self.draw_cell()

  def relax_step_forward ( self ):
    self.relax_step(1)

  def relax_step_backward ( self ):
    self.relax_step(-1)

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
    return vp.vector(a[0], a[1], a[2])
  
  def calc_dist ( self, a, b ):
      '''
      Calculate the distance between two atoms

      Arguments:
        a (Atom or vpython.sphere): First atom
        b (Atom or vpython.sphere): Second atom
      '''
      dist = b.pos - a.pos
      dist = dist.mag
      if self.unit_menu.selected == 'Angstrom':
        dist *= .529177
        text = '%f angstroms'%dist
      else:
        text = '%f a.u.'%dist
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
    angle = np.arccos((v1.dot(v2))/(v1.mag*v2.mag))
    if self.unit_menu.selected == 'Degree':
      text = '%f degrees'%np.degrees(angle)
    else:
      text = '%f radians'%angle
    print('Angle is %s'%text)
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
    return np.array([ix,iy,iz]) @ lat

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
    self.vBonds = vpobject_destructor(self.vBonds)
    self.coord_axes = vpobject_destructor(self.coord_axes)
    self.bound_curve = vpobject_destructor(self.bound_curve)

  def draw_coord_axes ( self, offset=[-10,0,0], length=1. ):
    '''
    Draw the Coordinate Axes

    Arguments:
      length (int): Length of the coordiniate axis arrows
    '''
    apos = np.array(offset)
    vpa = lambda a : vp.arrow(pos=self.vector(apos),axis=self.vector(length*np.array(a)))
    self.coord_axes = [vpa([2,0,0]),vpa([0,2,0]),vpa([0,0,2])]

  def draw_bonds ( self, atoms ):
    '''
    Draw the bonds if atoms are closer than 'dist' together.
    To specifiy bonds between two types of atoms pass a dictionary into dist with key/value
    pairs '%s_%s':val, where %s are the respective species and val is the maximum bond distance.

    Arguments:
      dits (float or dict): A float describes global bond distance, while a dictionary can specify distances for various pairs of atoms
    '''
    from .model import get_bond_pairs

    self.vBonds = []
    dList = lambda a : {'pos':a.pos, 'color':a.col, 'radius':self.bond_radius}
    for pair in get_bond_pairs(self.natoms, atoms, self.bonds, self.basis_labels):
      a1 = self.vAtoms[pair[0]]
      a2 = self.vAtoms[pair[1]]
      self.vBonds.append(vp.curve(dList(a1),dList(a2)))

  def draw_cell ( self ):
    '''
    Create the cell simulation in the vpython environment

    Arguments:
      lat (list or ndarray): Lattice vectors
      atoms (list or ndarray): Atoms to draw on lattice
    '''
    from .Atom import Atom

    self.clear_canvas()
    nx,ny,nz = self.cell_dim

    lat = self.lattice
    if self.boundary:
      self.bound_curve = []
      for ix in range(nx):
        for iy in range(ny):
          for iz in range(nz):
            corner = np.sum(lat, axis=0)
            alines = [[3*[0],a] for a in lat]
            alines += [[p,corner] for p in [lat[i]+lat[j] for i in range(3) for j in range(i+1,3)]]
            alines += [[lat[k],lat[i]+lat[j]] for i in range(3) for j in range(i+1,3) for k in [i,j]]
            orig = self.vector(self.get_cell_pos(lat,ix,iy,iz))
            self.bound_curve += [vp.curve(orig+self.vector(al[0]),orig+self.vector(al[1]),color=self.bnd_col) for al in alines]

    atoms = self.atoms
    species = self.species
    labels = self.basis_labels

    patoms = []
    self.vAtoms = []
    for ix in range(nx):
      for iy in range(ny):
        for iz in range(nz):
          c_pos = self.origin + self.get_cell_pos(lat,ix,iy,iz)
          for i,a in enumerate(atoms):
            spec = species[labels[i]]
            a_pos = c_pos + a
            patoms.append(a_pos)
            self.vAtoms.append(Atom(self.vector(a_pos), col=self.vector(spec['color']), radius=spec['radius']))

    if len(self.vAtoms) > 1:
      self.draw_bonds(patoms)

    self.canvas.center = self.vector(np.mean([[v.pos.x,v.pos.y,v.pos.z] for v in self.vAtoms], axis=0))
    if not self.coord_axes is None:
      self.draw_coord_axs()

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

  def draw_surface ( self, r_space, lat, data, iso, bands, colors, normals, write_obj, ):
    '''
    Draw a surface.
    If in the reciprocal space, create the Brillouin Zone boundary and bsxf points between 'fermiup' and 'fermidw' in the vpython window

    Arguemnts:
      r_space (bool): True for real space, false for reciprocal.
      lat (list or ndarray): 3 3d vectors representing the lattice or reciprocal lattice
      data (ndarray): eigenvalue data for each band and point in the BZ to plot
      iso (list): List of floats corresponding to the isosurface values for each respective band in 'bands'
      bands (list): List of integers representing the index of the band to plot
      colors (list): List of 3-d RGB color vectors for each band. If ignored, each band will be green
      normals (bool): True adds normals to triangle vertices, improving surface visibility
      write_obj (bool): True will write the triange vertex, face, and normal data to an obj file
    '''
    from .MarchingCubes import marching_cubes
    from numpy.linalg import det,norm,solve
    from scipy.fftpack import fftshift

    draw_flag = vp.text(text='Drawing...', align='center', color=vp.vector(1,0,0))
    draw_flag.up = self.canvas.up

    if r_space:
      if self.boundary:
        self.bound_curve = []
        nx=ny=nz=1
        for ix in range(nx):
          for iy in range(ny):
            for iz in range(nz):
              corner = np.sum(lat, axis=0)
              alines = [[3*[0],a] for a in lat]
              alines += [[p,corner] for p in [lat[i]+lat[j] for i in range(3) for j in range(i+1,3)]]
              alines += [[lat[k],lat[i]+lat[j]] for i in range(3) for j in range(i+1,3) for k in [i,j]]
              orig = self.get_cell_pos(lat,ix,iy,iz)
              for i in range(3):
                orig -= 0.5*lat[i]
              orig = self.vector(orig)
              self.bound_curve += [vp.curve(orig+self.vector(al[0]),orig+self.vector(al[1]),color=self.bnd_col) for al in alines]
    else:
      self.draw_BZ_boundary(lat)

    poss = []
    nx,ny,nz,nbnd = data.shape

    self.vAtoms = []
    for i,b in enumerate(bands):
      self.vAtoms += marching_cubes(r_space, data[:,:,:,b], iso[i], lat, self.BZ_planes, colors[i], write_obj=write_obj)

    if not self.coord_axes is None:
      self.draw_coord_axes(length=.1*np.linalg.norm(lat[0]),offset=[-1,0,0])

    draw_flag.visible = False
    del draw_flag
