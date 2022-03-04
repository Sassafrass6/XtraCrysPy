import numpy as np

class XtraCrysPy:

  def __init__ ( self, size=(1024, 1024), structure=None ):
    '''

    Arguments:
    '''
    from fury import window

    self.aposs = None
    self.atoms = None
    self.bonds = None
    self.natoms = None
    self.picker = None
    self.sel_inds = []
    self.sel_cols = []
    self.stype = 'info'
    self.scolor = np.array((0,210,210))
    
    self.scene = window.Scene()
    self.scene.background((1,1,1))
    self.smanager = window.ShowManager(self.scene, size=size, order_transparent=True)
    self.smanager.initialize()

    self.axes = None


  def set_atom_color ( self, mem, index, nvert, col ):
    for i in range(index*nvert, (index+1)*nvert):
      mem[i] = np.array(col)


  def push_atom ( self, mem, index, nvert ):
    self.sel_inds.append(index)
    self.sel_cols.append(mem[index*nvert].copy())
    self.set_atom_color(mem, index, nvert, self.scolor)


  def pop_atom ( self, mem, index, nvert ):
    aind = self.sel_inds.index(index)
    self.sel_inds.pop(aind)
    col = self.sel_cols.pop(aind)
    self.set_atom_color(mem, index, nvert, col)


  def left_click ( self, obj, event ):
    from fury.utils import colors_from_actor,update_actor,vertices_from_actor

    pos = self.picker.event_position(self.smanager.iren)
    pinfo = self.picker.pick(pos, self.smanager.scene)

    vertices = vertices_from_actor(obj)
    colors = colors_from_actor(obj, 'colors')

    nvert = int(vertices.shape[0]/self.natoms)
    index = int(np.floor(pinfo['vertex']/nvert))

    if self.stype == 'info':
      if index in self.sel_inds:
        self.pop_atom(colors, index, nvert)
      else:
        self.push_atom(colors, index, nvert)
    elif self.stype == 'distance':
      if index in self.sel_inds:
        self.pop_atom(colors, index, nvert)
      elif len(self.sel_inds) == 0:
        self.push_atom(colors, index, nvert)
      elif len(self.sel_inds) == 1:
        self.push_atom(colors, index, nvert)
        ## Draw Line
        ## Compute Distance
    update_actor(obj)

  def start_crystal_view ( self ):
    from fury import pick

    self.picker = pick.PickingManager()
    self.atoms.AddObserver('LeftButtonPressEvent', self.left_click, 1)
    self.smanager.start()


  def render_atomic_model ( self, model, nsc=(1,1,1), bond_type='stick' ):
    '''
    '''
    from fury import actor
    import numpy as np

    ainfo,binfo = model.lattice_atoms_bonds(*nsc, bond_type)
 
    self.aposs = ainfo[0]
    self.natoms = ainfo[0].shape[0]
    if self.natoms > 0:
      self.atoms = actor.sphere(centers=ainfo[0], colors=ainfo[1], radii=ainfo[2])
      self.scene.add(self.atoms)

    self.bonds = []
    for i in range(binfo[0].shape[0]):
      tbond = actor.cylinder(binfo[0][i], binfo[1][i], binfo[2][i], radius=binfo[3][i], heights=binfo[4][i], resolution=20)
      self.scene.add(tbond)
      self.bonds.append(tbond)


  def show_axes ( self, show=True ):
    '''
    '''
    return
    from fury import actor

    self.axes = actor.axes()
    self.scene.add(self.axes)
    self.scene.reset_camera()
