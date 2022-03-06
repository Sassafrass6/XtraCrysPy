from fury.utils import colors_from_actor,update_actor,vertices_from_actor
import numpy as np

class XtraCrysPy:

  def __init__ ( self, size=(1024, 1024), axes=True, perspective=False ):
    '''
    Arguments:
    '''
    from fury import ui,window

    self.wsize = size
    self.picker = None
    self.sel_inds = []
    self.sel_cols = []
    self.stype = 'info'
    self.scolor = np.array((0,210,210))

    self.scene = window.Scene()
    self.smanager = window.ShowManager(self.scene, size=size, order_transparent=True)
    self.smanager.initialize()

    self.axes = axes
    if axes:
      self.axis_eles = [None]*3
      self.axis_lines = [None]*3
      self.axis_width = ax_wid = 200
      self.ax_panel = ui.Panel2D(size=(ax_wid,ax_wid), color=(0,0,0))
      self.ax_panel.center = (ax_wid/2, ax_wid/2)
      self.scene.add(self.ax_panel)

      axis_vecs = np.array([[1,0,0],[0,1,0],[0,0,1]])
      centers = self.axis_endpoint_positions(axis_vecs)
      for i in range(3):
        self.axis_eles[i] = ui.Disk2D(outer_radius=8, center=centers[i], color=axis_vecs[i])
        self.axis_lines[i] = ui.Line2D((ax_wid/2,ax_wid/2), centers[i], color=axis_vecs[i])
        self.axis_lines[i].width = 5

        self.scene.add(self.axis_eles[i])
        self.scene.add(self.axis_lines[i])


  def axis_endpoint_positions ( self, vecs ):
    import numpy as np
    centers = np.empty((3,3), dtype=float)
    for i in range(3):
      wid = self.axis_width / 2
      centers[i,:] = np.array([wid,wid,0]) + 80*vecs[i]
    return centers[:,:2]


  def left_click ( self, obj, event ):
    pass


  def update_axes ( self, caller, event ):

    camera = self.scene.GetActiveCamera()

    y = np.array(camera.GetViewUp())
    y /= np.linalg.norm(y)
    z = np.array(camera.GetFocalPoint()) - np.array(camera.GetPosition())
    z /= np.linalg.norm(z)
    x = np.cross(y, z)
    x /= np.linalg.norm(x)

    rot = []
    unit = np.eye(3)
    cam_vecs = np.array([-x,y,z])
    for i,u in enumerate(unit[:2]):
      v = np.cross(cam_vecs[i], u)
      r = np.array([[0,-v[2],v[1]],[v[2],0,-v[0]],[-v[1],v[0],0]])
      rot.append(np.eye(3) + r + r@r/(1+cam_vecs[i].dot(u)))
      cam_vecs = rot[-1] @ cam_vecs.T

    axis_vecs = np.eye(3) @ np.linalg.inv(rot[0]) @ np.linalg.inv(rot[1])

    for i,v in enumerate(axis_vecs):
        axis_vecs[i] /= np.linalg.norm(v)
    centers = self.axis_endpoint_positions(axis_vecs)
    axis_order = np.argsort([v[2] for v in axis_vecs])
    for i,o in enumerate(axis_order):
      color = np.zeros(3); color[o] = 1
      self.axis_lines[i].color = color
      self.axis_lines[i].p2 = centers[o]
      self.axis_eles[i].color = color
      self.axis_eles[i].center = centers[o]


  def start_crystal_view ( self ):
    from fury import pick

    self.picker = pick.PickingManager()
    self.atoms.AddObserver('LeftButtonPressEvent', self.left_click, 1)
    if self.axes:
      self.smanager.iren.AddObserver('InteractionEvent', self.update_axes)
    self.smanager.start()


class XCP_Atoms ( XtraCrysPy ):


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


  def selection_logic ( self, colors, index, nvert ):

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


  def left_click ( self, obj, event ):

    pos = self.picker.event_position(self.smanager.iren)
    pinfo = self.picker.pick(pos, self.smanager.scene)

    vertices = vertices_from_actor(obj)
    colors = colors_from_actor(obj, 'colors')

    nvert = int(vertices.shape[0]/self.natoms)
    index = int(np.floor(pinfo['vertex']/nvert))

    self.selection_logic(colors, index, nvert)

    update_actor(obj)


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

