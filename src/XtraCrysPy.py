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
    self.sel_bnds = []
    self.units = 'bohr'
    self.runits = 'degree'
    self.sel_type = 'distance'
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


  def angle ( self, ai0, ai1, ai2 ):
    from numpy.linalg import norm

    v1 = self.aposs[ai0] - self.aposs[ai1]
    v2 = self.aposs[ai2] - self.aposs[ai1]
    conv = 1 if self.runits == 'radian' else 180/np.pi
    return conv * np.arccos((v1.dot(v2))/(norm(v1)*norm(v2)))


  def distance ( self, ai1, ai2 ):
    dist = np.sqrt(np.sum((self.aposs[ai2]-self.aposs[ai1])**2))
    if self.units == 'angstrom':
      from .conversion import BOHR_ANG
      dist *= BOHR_ANG
    else:
      pass # Bohr
    return dist


  def set_atom_color ( self, mem, index, nvert, col ):
    for i in range(index*nvert, (index+1)*nvert):
      mem[i] = np.array(col)


  def pop_sbond ( self ):
    tbond = self.sel_bnds.pop()
    self.scene.rm(tbond)


  def push_sbond ( self, ai1, ai2 ):
    from fury import actor

    conn = self.aposs[ai1] - self.aposs[ai2]
    dist = np.linalg.norm(conn)
    cent = (self.aposs[ai1] + self.aposs[ai2]) / 2

    if self.bond_type == 'stick':
      brad = 0.01 + 1 / (4 * dist)
    elif self.bond_type == 'primary':
      prad = np.min([r for k,r in self.model.radii.items()])
      brad = 0.01 + 1.5 * prad / (2 * dist)
    tbond = actor.cylinder([cent], [conn], [self.scolor/255], radius=brad, heights=dist, resolution=20)

    self.scene.add(tbond)
    self.sel_bnds.append(tbond)


  def pop_atom ( self, mem, index, nvert ):
    aind = self.sel_inds.index(index)
    self.sel_inds.pop(aind)
    col = self.sel_cols.pop(aind)
    self.set_atom_color(mem, index, nvert, col)


  def push_atom ( self, mem, index, nvert ):
    self.sel_inds.append(index)
    self.sel_cols.append(mem[index*nvert].copy())
    self.set_atom_color(mem, index, nvert, self.scolor)


  def selection_logic ( self, colors, index, nvert ):

    if index in self.sel_inds:
      nsi = len(self.sel_inds)
      sind = self.sel_inds.index(index)
      npop = nsi - sind
      if nsi == 1 or sind == 0:
        npop -= 1
        self.pop_atom(colors, self.sel_inds[-1], nvert)
      for _ in range(npop):
        self.pop_sbond()
        self.pop_atom(colors, self.sel_inds[-1], nvert)

    else:
      if self.sel_type == 'info':
        self.push_atom(colors, index, nvert)

      elif self.sel_type == 'distance':
        if len(self.sel_inds) == 0:
          self.push_atom(colors, index, nvert)
        elif len(self.sel_inds) == 1:
          self.push_atom(colors, index, nvert)
          self.push_sbond(*self.sel_inds)
          dist = self.distance(*self.sel_inds)
          print('Distance between atoms {} and {}:'.format(*self.sel_inds))
          print('{} {}\n'.format(self.distance(*self.sel_inds),self.units))

      elif self.sel_type == 'angle':
        if len(self.sel_inds) == 0:
          self.push_atom(colors, index, nvert)
        elif len(self.sel_inds) in [1,2]:
          self.push_atom(colors, index, nvert)
          self.push_sbond(*self.sel_inds[-2:])
          if len(self.sel_inds) == 3:
            ## Do the stuff
            print('Angle between atoms {}, {}, and {}:'.format(*self.sel_inds))
            print('{} {}\n'.format(self.angle(*self.sel_inds),self.runits))

      elif self.sel_type == 'chain':
        if len(self.sel_inds) == 0:
          self.push_atom(colors, index, nvert)
        else:
          self.push_atom(colors, index, nvert)
          self.push_sbond(*self.sel_inds[-2:])


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

    self.bond_type = bond_type

    if model.units != self.units:
      self.units = model.units

    ainfo,binfo = model.lattice_atoms_bonds(*nsc, bond_type)
 
    self.model = model
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

